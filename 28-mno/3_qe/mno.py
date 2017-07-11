#!/usr/bin/env python

def apply_machine_settings():
  from nexus import settings
  settings(
    runs       = 'default',
    machine    = 'ws4',
    pseudo_dir = './pseudo'
  )
# end def apply_machine_settings

def mno_structure(fname):
  from nexus import Structure
  struct = Structure()
  struct.read(fname)
  return struct
# end def

if __name__ == '__main__':
  apply_machine_settings()
  struct = mno_structure('../1_cell/mno.xsf')
  from nexus import generate_physical_system
  system = generate_physical_system(structure=struct)

  ecut = 40
  func = 'pbe'
  exx  = 0.0
  myid = '%s-%d' % ( func,round(exx*100.) )
  from nexus import obj, Job
  pwscf_inputs = obj(
    identifier = myid + '-scf',
    path       = myid + '/scf',
    job        = Job(),
    input_type = 'scf',
    input_dft = 'pbe',
    system  = system.copy(),
    pseudos = ['Mn.pbe-sp-van.UPF','O.pbe-van_bm.UPF'],
    ecut    = ecut,
    kgrid   = [1,1,1],
    kshift  = [0,0,0]  # gamma point
  )

  from nexus import generate_pwscf
  scf = generate_pwscf(**pwscf_inputs)

  from nexus import ProjectManager
  pm = ProjectManager()
  pm.add_simulations(scf)
  pm.run_project()
