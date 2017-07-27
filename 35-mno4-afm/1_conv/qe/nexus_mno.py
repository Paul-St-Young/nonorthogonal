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
    runs       = 'convergence',
    machine    = machine,
    pseudo_dir = './pseudo',
  )
  if machine == 'quartz':
    pbe_job = Job(nodes=1,cores=4,minutes=30)
    hse_job = Job(nodes=1,cores=4,hours=6)
    lda_plus_u_job = Job(nodes=1,cores=4,hours=6)
  else: # workstation, defaults should do
    pbe_job = Job()
    hse_job = Job()
    lda_plus_u_job= Job()
  # end if

  jobs = {'pbe':pbe_job,'hse':hse_job,'ldau':lda_plus_u_job}
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
    disk_io = 'none'
  )

  if hubbard_u is not None:
    pwscf_inputs['hubbard_u'] = hubbard_u
  # end if

  return pwscf_inputs
# end def default_scf_input

def inputs_to_scan_ecut(ecut_list,scf_job,system):
  # setting defaults
  func = 'pbe'
  exx  = 0.0

  inputs = []
  for ecut in ecut_list:
    pwscf_inputs = default_scf_input(func,exx,ecut,scf_job,system.copy())
    inputs.append(pwscf_inputs.copy())
  # end for
  return inputs
# end def

def append_ecut_scan(scf_sims,jobs,system):
  subdir = 'ecut'
  inputs = inputs_to_scan_ecut([20,40,80,160,320,640],jobs['pbe'],system)
  for myinput in inputs:
    myinput.path = os.path.join(subdir,myinput.path)
    scf = generate_pwscf(**myinput)
    scf_sims.append(scf)
  # end for
# end def

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

if __name__ == '__main__':
  xsf_file = '../../struct/mno.xsf'
  jobs = apply_machine_settings('quartz')

  struct = get_structure(xsf_file)
  from nexus import generate_physical_system
  system = generate_physical_system(structure=struct,Mn1=15,Mn2=15,O=6)

  from nexus import generate_pwscf, generate_pw2qmcpack, generate_qmcpack
  scf_sims = []
  append_ecut_scan(scf_sims,jobs,system)
  append_exx_scan(scf_sims,jobs,system)
  append_u_scan(scf_sims,jobs,system)

  from nexus import ProjectManager
  pm = ProjectManager()
  pm.add_simulations(scf_sims)
  pm.run_project()

  # analysis
  data = []
  for scf in scf_sims:
    sa = scf.load_analyzer_image()
    data.append( sa.to_dict() )
  # end for
  df = pd.DataFrame(data)
  df.to_json( 'mno4_qe_ecut_exx.json' )

# end __main__
