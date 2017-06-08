#!/usr/bin/env python
import numpy as np
from nexus import settings, Job, obj
from nexus import PwscfInput
from nexus import generate_pwscf,generate_pw2qmcpack,generate_qmcpack

if __name__ == '__main__':

    settings(
        runs       = 'default',
        pseudo_dir = '/home/yyang173/Desktop/phases/pseudo',
        machine    = 'ws16',
        generate_only = 0,
        status_only = 0,
        sleep = 3
    )

    #scf_job = Job(cores=4,app='pw.x',app_options='-nk 4')
    scf_job = Job(serial=True,app='pw.x')
    p2q_job = Job(serial=True
      ,app='/home/yyang173/soft/qmcpack-espresso-5.3.0/bin/pw2qmcpack.x')

    pi = PwscfInput('./scf.in')
    system = pi.return_system(H=1)
    system.structure.inversion_symmetrize_kpoints()

    scf_input = obj(
        identifier = 'scf',
        path       = 'scf',
        input_type = 'scf',
        pseudos    = ['H.coulomb-ae.UPF'],
        occupations= 'smearing',
        smearing   = 'f-d',
        degauss    = 1e-4,#0.0019, # 300 K, close enough to zero
        ecut       = 500, #20,
        input_dft  = 'lda',
        system     = system,
        job        = scf_job,
        kgrid      = (1,1,1),
        kshift     = (0,0,0)
    )
    scf = generate_pwscf(**scf_input)

    p2q_input = obj(
        identifier = 'p2q',
        path       = scf_input.path,
        dependencies = (scf,'orbitals'),
        job = p2q_job
    )
    p2q = generate_pw2qmcpack(**p2q_input)

    from nexus import ProjectManager
    pm = ProjectManager()
    pm.add_simulations(scf,p2q)
    pm.run_project()
