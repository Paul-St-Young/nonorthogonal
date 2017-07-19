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
    pseudo_dir = './pseudo'
  )
  if machine == 'quartz':
    scf_job = Job(nodes=1,cores=1,hours=4)
    p2q_job = Job(nodes=1,serial=True,minutes=5)
    dmc_job = Job(nodes=20,hours=24,app='/g/g91/yang41/soft/master_qmcpack/build/bin/qmcpack_comp')
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
    input_type = 'scf',
    input_dft    = func,
    exx_fraction = exx,
    pseudos = ['C.BFD.upf'],
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
  ecut   = 160

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
  rcut    = 4.810456
  j1coeff_text = '-0.2881867991 -0.2164434839 -0.1447128309 -0.06841238183 -0.02040495028 -0.0009045054505 0.001882471792 -0.0001330702941'
  j2uu_text = '0.4072085225 0.2869348431 0.1938006249 0.123350439 0.07387213081 0.04090015102 0.02000131129 0.008510794954'
  j2ud_text = '0.617109107 0.3988359485 0.249481554 0.1494738224 0.0855004812 0.04638206743 0.02267617002 0.009360860175'

  j1coeff = map(float,j1coeff_text.split())
  j2uu    = map(float,j2uu_text.split())
  j2ud    = map(float,j2ud_text.split())

  from nexus import generate_jastrow1,generate_jastrow2
  j1 = generate_jastrow1(function='bspline',size=len(j1coeff),coeff=[j1coeff],elements=['C'])
  j2 = generate_jastrow2(function='bspline',size=len(j2uu),coeff=[j2uu,j2ud],init=None) 
  return j1,j2
# end def

def gamma_dmc_input(p2q,dmc_job,system):
  from nexus import vmc,dmc

  myid = p2q.identifier.replace('-p2q','')

  # dmc time steps to try
  tss = [0.025,0.01]
  # vmc correlation time in a.u.
  correlation_time = 0.2

  j1,j2 = hf_jastrows()
  vmc_input = obj(
    warmupsteps =  40,
    blocks      = 200,
    steps       =  10,
    substeps    =   3,
    timestep    = 1.0,
    samples     =5760 # dmc walkers
  )
  dmc_input = obj(
    warmupsteps = 40,
    blocks      = 100,
    steps       = 40,   # will be overwritten
    timestep    = 0.02  # will be overwritten
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
    jastrows     = [j1,j2],
    pseudos      = ['C.BFD.xml'],
    dependencies = [(p2q,'orbitals')]
  )
  return dmc_inputs
# end def gamma_dmc_input

if __name__ == '__main__':
  jobs = apply_machine_settings('quartz')
  struct = get_structure('c2-16.xsf')
  from nexus import generate_physical_system
  system = generate_physical_system(structure=struct,C=4)

  from nexus import generate_pwscf, generate_pw2qmcpack, generate_qmcpack
  scf_sims = []
  p2q_sims = []
  dmc_sims = []
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
  #myinput = gamma_scf_input('pbe0',1.0,160,jobs['scf'],system)
  #scf = generate_pwscf(**myinput)
  #scf_sims.append(scf)

  from nexus import ProjectManager
  pm = ProjectManager()
  pm.add_simulations(scf_sims+p2q_sims+dmc_sims)
  pm.run_project()

# end __main__