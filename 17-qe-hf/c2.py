#!/usr/bin/env python
import numpy as np
from generic import obj
from nexus import settings,Job
from nexus import Structure,generate_physical_system
from nexus import generate_pwscf, generate_pw2qmcpack
from nexus import generate_qmcpack,loop,linear,vmc,dmc
from nexus import ProjectManager

if __name__ == '__main__':
  # machine settings
  settings(
    runs    = 'gamma',
    machine = 'ws16',
    pseudo_dir = './pseudo'
  )
  scf_job  = Job(app='pw.x',serial=True)
  p2q_job  = Job(app='pw2qmcpack.x',serial=True)
  opt_job  = Job(app='qmcpack')

  # define crystal structure
  alat = 3.6
  axes = (np.ones((3,3))-np.eye(3))*alat/2.
  elem = ['C','C']
  bravais_basis = [[0,0,0],[0.25,0.25,0.25]]
  pos = np.dot(bravais_basis,axes)

  structure = Structure(
    axes  = axes,
    elem  = elem,
    pos   = pos,
    units = 'A',
  )

  # simulation parameters
  ecut   = 20
  #kgrid0 = [2,2,2]
  #kgrid  = [2,2,2]
  kgrid0 = [1,1,1]
  kshift0= [0,0,0]
  kgrid  = [1,1,1]
  kshift = [0,0,0]

  system = generate_physical_system(
    structure = structure,
    C     = 4., # specify number of valence electrons (valency)
    kgrid = [2,2,2]
  )

  sims = []
  for ecut in [20,40,80,160]:
      subdir = 'ecut%d'%ecut

      # generate simulation objects
      scf_inputs = obj(
        input_type = 'scf',
        input_dft  = 'pbe0',
        exx_fraction = 1.0,
        path       = subdir+'/scf',
        identifier = 'scf',
        system    = system.copy(),
        pseudos   = ['C.BFD.upf'],
        kgrid     = kgrid0,
        kshift    = kshift0,
        job       = scf_job,
        ecut      = ecut
      )
      scf = generate_pwscf(**scf_inputs)
      sims.append(scf)

      p2q = generate_pw2qmcpack(
        path   = scf.path,
        outdir = scf.input.control.outdir,
        identifier   = 'p2q',
        write_psir   = False,
        job          = p2q_job,
        dependencies = (scf,'orbitals')
      )
      sims.append(p2q)

      linopt = obj(
        energy               = 0.95,
        unreweightedvariance = 0.0,
        reweightedvariance   = 0.05,
        timestep             = 1.0,
        substeps             = 5,
        warmupsteps          = 10,
        blocks               = 200,
        samples              = 8192,
        usedrift             = True
      )

      opt = generate_qmcpack(
        path = subdir+'/opt',
        identifier = 'opt',
        job  = opt_job,
        input_type = 'basic',
        twistnum   = 0,
        system     = system.copy(),
        bconds     = 'ppp',
        jastrows   = [('J1','cusp',6),('J2','size',8)],
        pseudos    = ['C.BFD.xml'],
        calculations = [loop(max=8,qmc=linear(**linopt))],
        dependencies = (p2q,'orbitals')
      )
      sims.append(opt)

  # run simulations
  pm = ProjectManager()
  pm.add_simulations(sims)
  pm.run_project()

# end __main__
