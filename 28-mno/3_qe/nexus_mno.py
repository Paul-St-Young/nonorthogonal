#!/usr/bin/env python
import numpy as np
from nexus import obj

def apply_machine_settings(machine):
  if machine != 'quartz' and (not machine.startswith('ws')):
    raise NotImplementedError('cannot handle machine=%s yet'%machine)
  # end if
  from nexus import settings,Job
  settings(
    runs       = 'default',
    machine    = machine,
    pseudo_dir = './pseudo',
  )
  if machine == 'quartz':
    scf_job = Job(nodes=1,cores=4,hours=6)
    p2q_job = Job(nodes=1,serial=True,minutes=15)
    dmc_job = Job(nodes=1,cores=36,hours=24,app='/g/g91/yang41/soft/master_qmcpack/build/bin/qmcpack_comp')
  else: # workstation, defaults should do
    scf_job = Job()
    p2q_job = Job()
    dmc_job = Job()
  # end if

  jobs = {'scf':scf_job,'p2q':p2q_job,'dmc':dmc_job}
  return jobs
# end def apply_machine_settings

def get_structure(fname):
  from nexus import Structure
  struct = Structure()
  struct.read(fname)
  return struct
# end def

def gamma_scf_input(func,exx,ecut,scf_job,system):
  kgrid  = [1,1,1]
  kshift = [0,0,0]

  myid = '%s-%d' % ( func,round(exx*100.) )

  pwscf_inputs = obj(
    identifier = myid + '-scf',
    path       = myid + '/scf',
    job        = scf_job,
    system     = system,
    conv_thr   = 1e-8,
    input_type = 'scf',
    input_dft    = func,
    exx_fraction = exx,
    pseudos = ['Mn.BFD.upf','O.BFD.upf'],
    ecut    = ecut,
    kgrid   = kgrid,
    kshift  = kshift
  )
  return pwscf_inputs
# end def gamma_scf_input

def inputs_to_scan_hse_exx(exx_list,scf_job,system):
  njob = len(exx_list)
  print('writing %d inputs'%njob)

  # setting defaults
  func   = 'hse'
  ecut   = 80

  inputs = []
  for exx in exx_list:
    pwscf_inputs = gamma_scf_input(func,exx,ecut,scf_job,system.copy())
    inputs.append(pwscf_inputs.copy())
  # end for

  return inputs
# end def

def p2q_input_from_scf(scf,p2q_job):
  p2q_input = obj(
    identifier = scf.identifier.replace('-scf','-p2q'),
    path = scf.path,
    outdir = scf.input.control.outdir,
    write_psir = False,
    job = p2q_job,
    dependencies = (scf,'orbitals')
  )
  return p2q_input
# end def

def hf_jastrows():
  rcut    = 2.958196
  mn_j1coeff_text = '-0.6780760628 -1.145683226 -0.7354441945 -0.5265647833 -0.2481322844 -0.1165250299 -0.04045642633 -0.02247358863'
  o_j1coeff_text  = '-0.8618402525 -0.3445886344 -0.2663570607 -0.2221275873 -0.0973739713 -0.1114785209 -0.01759614592 -0.0463741941'
  j2uu_text = '0.3106257352 0.2245131638 0.1647228383 0.1118443597 0.07590016152 0.04741092948 0.02459199409 0.01175994094'
  j2ud_text = '0.4393978199 0.2999445804 0.2106775191 0.1424467661 0.09722690905 0.05934016809 0.03218846169 0.01604788512'

  mn_j1coeff = map(float,mn_j1coeff_text.split())
  o_j1coeff = map(float,o_j1coeff_text.split())
  j2uu    = map(float,j2uu_text.split())
  j2ud    = map(float,j2ud_text.split())

  from nexus import generate_jastrow1,generate_jastrow2
  j1mn = generate_jastrow1(function='bspline',size=len(mn_j1coeff),coeff=[mn_j1coeff],elements=['Mn'])
  j1o = generate_jastrow1(function='bspline',size=len(o_j1coeff),coeff=[o_j1coeff],elements=['O'])
  j2 = generate_jastrow2(function='bspline',size=len(j2uu),coeff=[j2uu,j2ud],init=None) 
  return j1mn,j1o,j2
# end def

def gamma_dmc_input(p2q,dmc_job,system):
  from nexus import vmc,dmc

  myid = p2q.identifier.replace('-p2q','')

  # dmc time steps to try
  tss = [0.01,0.005]
  # vmc correlation time in a.u.
  correlation_time = 1.0

  jas = hf_jastrows()
  vmc_input = obj(
    warmupsteps =  40,
    blocks      = 200,
    steps       =  10,
    substeps    =   3,
    timestep    = 1.0,
    walkers     = 144,
    samples     =4320, # dmc walkers
    checkpoint  = 0
  )
  dmc_input = obj(
    warmupsteps = 40,
    blocks      = 1000,
    steps       = 40,   # will be overwritten
    timestep    = 0.02, # will be overwritten
    checkpoint  = 100
  )
  calcs = [vmc(**vmc_input)]
  for ts in tss:
    dmc_input.steps   = int(round(correlation_time/ts))
    dmc_input.timestep = ts
    calcs.append(dmc(**dmc_input))
  # end for ts
  mysystem = system.copy()
  mysystem.structure.kpoints = np.array([[0.,0.,0.]]) # system must have kpoints to assign twistnums
  dmc_inputs  = obj(
    identifier  = myid + '-dmc',
    path        = myid + '/dmc',
    job         = dmc_job,
    input_type  = 'basic',
    system      = mysystem,
    bconds      = 'ppp',    # periodic in xyz directions
    calculations = calcs,
    twistnum    = 0,
    estimators   = [],
    jastrows     = jas,
    pseudos      = ['Mn.BFD.xml','O.BFD.xml'],
    dependencies = [(p2q,'orbitals')]
  )
  return dmc_inputs
# end def gamma_dmc_input

def fake_structure():
  struct = get_structure('mno.xsf')
  struct.change_units('B')
  from nexus import Structure
  fake_struct = Structure(axes=struct.axes,elem=struct.elem,pos=struct.pos,units='A')
  fake_struct.write('ang_mno.xsf')
# end def

if __name__ == '__main__':
  jobs = apply_machine_settings('quartz')

  #fake_structure()
  struct = get_structure('mno.xsf')
  from nexus import generate_physical_system
  system = generate_physical_system(structure=struct,Mn=15,O=6)

  from nexus import generate_pwscf, generate_pw2qmcpack, generate_qmcpack
  scf_sims = []
  p2q_sims = []
  dmc_sims = []

  ## one test calculation
  #myinput = gamma_scf_input('hse',0.3,80,jobs['scf'],system)
  #scf = generate_pwscf(**myinput)
  #scf_sims.append(scf)

  # scan exchange fraction
  exx2scan   = np.linspace(0,1,10)
  hse_inputs = inputs_to_scan_hse_exx(exx2scan,jobs['scf'],system)
  for myinput in hse_inputs:
    scf = generate_pwscf(**myinput)
    p2q_inputs = p2q_input_from_scf(scf,jobs['p2q'])
    p2q = generate_pw2qmcpack(**p2q_inputs)
    dmc_inputs = gamma_dmc_input(p2q,jobs['dmc'],system)
    sdmc = generate_qmcpack(**dmc_inputs)
    scf_sims.append(scf)
    p2q_sims.append(p2q)
    dmc_sims.append(sdmc)
  # end for

  from nexus import ProjectManager
  pm = ProjectManager()
  pm.add_simulations(scf_sims+p2q_sims+dmc_sims)
  pm.run_project()

# end __main__
