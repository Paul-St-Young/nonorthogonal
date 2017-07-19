#!/usr/bin/env python
import numpy as np
import pandas as pd
from nexus import obj

def apply_machine_settings(machine,run_name):

  if machine != 'quartz' and (not machine.startswith('ws')):
    raise NotImplementedError('cannot handle machine=%s yet'%machine)
  # end if

  if run_name not in ['bfd','hgh']:
    raise NotImplementedError('please add an \"if block\" for pseudopotentials %s'%run_name)
  # end if

  from nexus import settings,Job
  settings(
    runs       = run_name,
    machine    = machine,
    pseudo_dir = './pseudo'
  )

  if machine == 'quartz':
    scf_job = Job(nodes=1,cores=4,minutes=30)
    p2q_job = Job(nodes=1,serial=True,minutes=5)
    dmc_job = Job(nodes=1,cores=36,minutes=30,app='/g/g91/yang41/soft/master_qmcpack/build/bin/qmcpack_comp')
  else: # workstation, defaults should do
    scf_job = Job()
    p2q_job = Job(serial=True)
    dmc_job = Job()
  # end if

  if run_name == 'bfd':
    dft_pseudos = ['Mn.BFD.upf','O.BFD.upf']
    qmc_pseudos = ['Mn.BFD.xml','O.BFD.upf']
  elif run_name == 'hgh':
    dft_pseudos = ['Mn.pbe-sp-hgh.UPF','O.pbe-hgh.UPF']
    qmc_pseudos = None
  # end if

  jobs = {'scf':scf_job,'p2q':p2q_job,'dmc':dmc_job}
  pseudos = {'dft':dft_pseudos,'qmc':qmc_pseudos}
  return jobs, pseudos
# end def apply_machine_settings

def get_structure(fname):
  from nexus import Structure
  struct = Structure()
  struct.read(fname)
  return struct
# end def

def gamma_scf_input(func,exx,ecut,scf_job,system,dft_pseudos):
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
    pseudos = dft_pseudos,
    ecut    = ecut,
    kgrid   = kgrid,
    kshift  = kshift
  )
  return pwscf_inputs
# end def gamma_scf_input

if __name__ == '__main__':

  run_id = 'jtk_p'
  myname = 'bfd'

  jobs, pseudos = apply_machine_settings('quartz',myname)
  struct = get_structure('mno.xsf')
  axes = struct.axes
  elem = struct.elem
  pos  = struct.pos
  units= struct.units

  from nexus import generate_physical_system, Structure
  from nexus import generate_pwscf, generate_pw2qmcpack, generate_qmcpack
  scf_sims = []
  p2q_sims = []
  dmc_sims = []

  ecut = 320
  strains = [0.95,0.98,1.0,1.02,1.05]
  for strain in strains:
    new_struct = Structure(axes=axes*strain,elem=elem,pos=pos*strain,units=units)
    system = generate_physical_system(structure=new_struct.copy(),Mn=15,O=6)
    myinput = gamma_scf_input('pbe',strain,ecut,jobs['scf'],system,pseudos['dft'])
    scf = generate_pwscf(**myinput)
    scf_sims.append(scf)
  # end for strain

  from nexus import ProjectManager
  pm = ProjectManager()
  pm.add_simulations(scf_sims+p2q_sims+dmc_sims)
  pm.run_project()

  # analysis
  data = []
  for scf in scf_sims:
    sa = scf.load_analyzer_image()
    data.append( sa.to_dict() )
  # end for
  df = pd.DataFrame(data)

  df.to_json( '%s-%s.json' % (run_id,myname) )

# end __main__
