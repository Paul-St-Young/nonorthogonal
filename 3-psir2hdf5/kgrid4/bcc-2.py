#!/usr/bin/env python
import numpy as np
from nexus import settings, Job, obj, PwscfInput, generate_physical_system
from nexus import generate_pwscf,generate_pw2qmcpack,generate_qmcpack
from nexus import vmc,dmc

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
    dmc_job = Job(app='qmcpack_comp')
    sims = []

    # !!!! return_system pseudizes the ions
    #pi = PwscfInput('./scf.in')
    #system = pi.return_system(H=1)
    #system.structure.inversion_symmetrize_kpoints()

    lat_basis = np.array([[0,0,0],[0.5,0.5,0.5]])
    natom = len(lat_basis)

    alat   = 3.77945227
    axes   = alat*np.eye(3) # cubic cell
    system = generate_physical_system(
        axes   = axes,
        pos    = np.dot(lat_basis,axes),
        elem   = natom*['H'],
        units  = 'B',
        net_charge = 0,
        net_spin   = 0,
        kgrid      = (4,4,4), # system must have kpoints to assign twistnums
        kshift     = (0,0,0)
    )

    scf_input = obj(
        identifier = 'scf',
        path       = 'scf',
        input_type = 'scf',
        pseudos    = ['H.coulomb-ae.UPF'],
        occupations= 'smearing',
        smearing   = 'f-d',
        degauss    = 0.0019, # 300 K, close enough to zero
        ecut       = 500,
        input_dft  = 'lda',
        system     = system,
        job        = scf_job,
        kgrid      = (4,4,4),
        kshift     = (0,0,0)
    )
    scf = generate_pwscf(**scf_input)
    sims.append(scf)

    p2q_input = obj(
        identifier = 'p2q',
        path       = scf_input.path,
        dependencies = (scf,'orbitals'),
        job = p2q_job
    )
    p2q = generate_pw2qmcpack(**p2q_input)
    sims.append(p2q)

    # dmc inputs
    vmc_input = obj(
        warmupsteps =  40,
        blocks      =  20,
        steps       =  10,
        timestep    = 1.0,
        samples     =  64 # dmc walkers
    )
    dmc_input = obj(
        blocks      = 40,
        steps       = 10,
        timestep    = 0.01
    )
    sdmc = generate_qmcpack(
        twistnum    = 0,
        identifier  = 'dmc',
        path        = 'dmc',
        job         = dmc_job,
        input_type  = 'basic',
        system      = system,
        bconds      = 'ppp',
        jastrows    = [],
        calculations= [vmc(**vmc_input),dmc(**dmc_input)],
        dependencies= [(p2q,'orbitals')]
    )
    sims.append(sdmc)

    from nexus import ProjectManager
    pm = ProjectManager()
    pm.add_simulations(sims)
    pm.run_project()
