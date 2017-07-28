#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
from nexus import obj

def apply_machine_settings(machine):
  if machine != 'quartz' and (not machine.startswith('ws')):
    raise NotImplementedError('cannot handle machine=%s yet'%machine)
  # end if
  from nexus import settings,Job
  settings(
    runs       = 'scan',
    machine    = machine,
    pseudo_dir = './pseudo',
    skip_submit= True
  )
  if machine == 'quartz':
    pbe_job = Job(nodes=1,cores=4,minutes=30,queue='pdebug')
    hse_job = Job(nodes=1,cores=36,hours=1,account='qmchhp')
    lda_plus_u_job = Job(nodes=1,cores=4,minutes=30,queue='pdebug')
    p2q_job = Job(nodes=1,serial=True,minutes=15,queue='pdebug')
    opt_job = Job(nodes=16,hours=2,app='/g/g91/yang41/soft/master_qmcpack/build/bin/qmcpack_comp',account='qmchhp')
    dmc_job = Job(nodes=16,hours=4,account='qmchhp',app='/g/g91/yang41/soft/master_qmcpack/build/bin/qmcpack_comp')
  else: # workstation, defaults should do
    pbe_job = Job()
    hse_job = Job()
    lda_plus_u_job= Job()
    p2q_job = Job()
    opt_job = Job()
    dmc_job = Job()
  # end if

  jobs = {'pbe':pbe_job,'hse':hse_job,'ldau':lda_plus_u_job,'p2q':p2q_job,'opt':opt_job,'dmc':dmc_job}
  return jobs
# end def apply_machine_settings

def get_structure(fname):
  # read xsf file using nexus
  from nexus import Structure
  struct = Structure()
  struct.read(fname)

  # extract structure specifications
  axes = struct.axes
  elem = struct.elem
  pos  = struct.pos
  units= struct.units

  # setup AFM structure 
  assert len(elem) == 4 # !!!! hard code for formula unit cell for now
  assert elem[0] == 'Mn'
  elem[0] = 'Mn1'
  assert elem[1] == 'Mn'
  elem[1] = 'Mn2'
  afm_struct = Structure(axes=axes,elem=elem,pos=pos,units=units)

  return afm_struct
# end def get_structure

# =================== DFT inputs ===================

def default_scf_input(func,exx,ecut,scf_job,system,hubbard_u=None,kgrid=[1,1,1],kshift=[0,0,0]):

  myid = '%s-%d-ecut%d' % ( func,round(exx*100.),ecut )

  pwscf_inputs = obj(
    identifier = myid + '-scf',
    path       = myid + '/scf',
    job        = scf_job,
    system     = system,
    conv_thr   = 1e-8,
    input_type = 'scf',
    input_dft    = func,
    exx_fraction = exx,
    start_mag    = {'Mn1':-0.5,'Mn2':0.5,'O':0},
    pseudos = ['Mn.BFD.upf','Mn.BFD.upf','O.BFD.upf'],
    ecut    = ecut,
    kgrid   = kgrid,
    kshift  = kshift,
  )

  if hubbard_u is not None:
    pwscf_inputs['hubbard_u'] = hubbard_u
  # end if

  return pwscf_inputs
# end def default_scf_input

def inputs_to_scan_hse_exx(exx_list,scf_job,system):
  # setting defaults
  func   = 'hse'
  ecut   = 320

  inputs = []
  for exx in exx_list:
    pwscf_inputs = default_scf_input(func,exx,ecut,scf_job,system.copy())
    inputs.append(pwscf_inputs.copy())
  # end for
  return inputs
# end def

def append_exx_scan(scf_sims,jobs,system):
  subdir     = 'exx'
  exx2scan   = [0.0,0.25,0.5,0.75,1.0]
  hse_inputs = inputs_to_scan_hse_exx(exx2scan,jobs['hse'],system)
  for myinput in hse_inputs:
    myinput.path = os.path.join(subdir,myinput.path)
    scf = generate_pwscf(**myinput)
    scf_sims.append(scf)
  # end for
# end def

def inputs_to_scan_lda_plus_u(u_list,scf_job,system):
  # setting defaults
  func   = 'lda'
  ecut   = 320
  exx    = 0.

  inputs = []
  for myu in u_list:
    pwscf_inputs = default_scf_input(func,exx,ecut,scf_job,system.copy(),hubbard_u={'Mn1':myu,'Mn2':myu,'O':0})
    mypath = pwscf_inputs.path
    subdir = os.path.dirname(mypath)
    mydir  = os.path.basename(mypath)
    pwscf_inputs.path = os.path.join( subdir, 'u%3.2f'%myu, mydir)
    inputs.append(pwscf_inputs.copy())
  # end for
  return inputs
# end def

def append_u_scan(scf_sims,jobs,system):
  subdir     = 'ldau'
  u2scan   = [0.5,5.0,10.0]
  hse_inputs = inputs_to_scan_lda_plus_u(u2scan,jobs['ldau'],system)
  for myinput in hse_inputs:
    myinput.path = os.path.join(subdir,myinput.path)
    scf = generate_pwscf(**myinput)
    scf_sims.append(scf)
  # end for
# end def

# =================== QMC inputs ===================

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

  # get coefficients 
  mn1_j1coeff_text = '-2.916291331 -2.79855945 -2.558033332 -2.226421313 -1.843365663 -1.411587226 -0.893282313 -0.4764249498'
  mn2_j1coeff_text = '-1.187103165 -1.064759167 -0.8450658369 -0.5410089754 -0.247669876 -0.02845378027 0.08280421186 0.11118505'
  o_j1coeff_text  = '-0.8078107508 -0.7385166615 -0.658887247 -0.5487006968 -0.4312140372 -0.3270768785 -0.2014092658 -0.08796479085'
  j2uu_text = '0.3072025819 0.2417595805 0.1812174375 0.1331879365 0.09271434913 0.0604246262 0.03300067526 0.01541665883'
  j2ud_text = '0.4420184348 0.3137089166 0.2241094917 0.1622985713 0.113953872 0.07556330908 0.04267892588 0.02095134822'

  mn1_j1coeff = map(float,mn1_j1coeff_text.split())
  mn2_j1coeff = map(float,mn2_j1coeff_text.split())
  o_j1coeff = map(float,o_j1coeff_text.split())

  j1_coeffs = [mn1_j1coeff,mn2_j1coeff,o_j1coeff]
  j1_nsize = len(j1_coeffs[0])
  for j1_coeff in j1_coeffs:
    assert len(j1_coeff) == j1_nsize
  # end for

  j2uu    = map(float,j2uu_text.split())
  j2ud    = map(float,j2ud_text.split())

  # use coefficients: j1_coeffs, j2uu, j2ud
  from nexus import generate_jastrow1,generate_jastrow2
  j1 = generate_jastrow1(function='bspline',size=j1_nsize
    ,coeff=j1_coeffs,elements=['Mn1','Mn2','O'])

  j2 = generate_jastrow2(function='bspline',size=len(j2uu),coeff=[j2uu,j2ud],init=None)
  return j1,j2
# end def hf_jastrows

def gamma_opt_input(p2q,opt_job,system):
  from nexus import loop, linear

  myid = p2q.identifier.replace('-p2q','-opt')
  nscf_dir = os.path.basename(p2q.path)
  mypath = p2q.path.replace(nscf_dir,'opt')

  linopt = obj(
    energy = 0.95,
    reweightedvariance = 0.05,
    unreweightedvariance = 0.0,
    warmupsteps =  40,
    blocks      = 200,
    steps       =  10,
    substeps    =   3,
    timestep    = 1.0,
    walkers     = 16,
    samples     = 144000,
    checkpoint  = 0
  )
  calcs = [loop(max=5,qmc=linear(**linopt))]

  init_jas = hf_jastrows()

  mysystem = system.copy()
  mysystem.structure.kpoints = np.array([[0.,0.,0.]]) # system must have kpoints to assign twistnums
  opt_inputs  = obj(
    identifier  = myid,
    path        = mypath,
    job         = opt_job,
    input_type  = 'basic',
    system      = mysystem,
    bconds      = 'ppp',    # periodic in xyz directions
    calculations = calcs,
    twistnum    = 0,
    estimators   = [],
    jastrows     = init_jas,
    pseudos      = ['Mn.BFD.xml','O.BFD.xml'],
    dependencies = [(p2q,'orbitals')]
  )
  return opt_inputs
# end def gamma_opt_input

def gamma_dmc_input(p2q,opt,dmc_job,system):
  from nexus import vmc,dmc

  myid = p2q.identifier.replace('-p2q','-dmc')
  nscf_dir = os.path.basename(p2q.path)
  mypath = p2q.path.replace(nscf_dir,'dmc')

  # dmc time steps to try
  tss = [0.01,0.005]
  # vmc correlation time in a.u.
  correlation_time = 1.0

  vmc_input = obj(
    warmupsteps =  40,
    blocks      = 200,
    steps       =  10,
    substeps    =   3,
    timestep    = 1.0,
    walkers     = 16,
    samples     =4320, # dmc walkers
    checkpoint  = 0
  )
  dmc_input = obj(
    warmupsteps = 40,
    blocks      = 40,
    steps       = 40,   # will be overwritten
    timestep    = 0.02, # will be overwritten
    checkpoint  = 0
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
    identifier  = myid,
    path        = mypath,
    job         = dmc_job,
    input_type  = 'basic',
    system      = mysystem,
    bconds      = 'ppp',    # periodic in xyz directions
    calculations = calcs,
    twistnum    = 0,
    estimators   = [],
    jastrows     = [],
    pseudos      = ['Mn.BFD.xml','O.BFD.xml'],
    dependencies = [(p2q,'orbitals'),(opt,'jastrow')]
  )
  return dmc_inputs
# end def gamma_dmc_input

if __name__ == '__main__':
  xsf_file = '../struct/mno.xsf'
  jobs = apply_machine_settings('quartz')

  struct = get_structure(xsf_file)
  from nexus import generate_physical_system
  system = generate_physical_system(structure=struct,Mn1=15,Mn2=15,O=6)

  # run DFT
  from nexus import generate_pwscf
  scf_sims = []
  append_exx_scan(scf_sims,jobs,system)
  append_u_scan(scf_sims,jobs,system)

  # take each DFT wavefunction and run DMC
  from nexus import generate_pw2qmcpack, generate_qmcpack
  p2q_sims = []
  opt_sims = []
  dmc_sims = []
  for scf in scf_sims:
    p2q_inputs = p2q_input_from_scf(scf,jobs['p2q'])
    p2q = generate_pw2qmcpack(**p2q_inputs)
    opt_inputs = gamma_opt_input(p2q,jobs['opt'],system)
    opt = generate_qmcpack(**opt_inputs)
    dmc_inputs = gamma_dmc_input(p2q,opt,jobs['dmc'],system)
    sdmc = generate_qmcpack(**dmc_inputs)
    p2q_sims.append(p2q)
    opt_sims.append(opt)
    dmc_sims.append(sdmc)
  # end for scf
  
  from nexus import ProjectManager
  pm = ProjectManager()
  pm.add_simulations(scf_sims+p2q_sims+opt_sims+dmc_sims)
  pm.run_project()

# end __main__
