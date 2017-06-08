#! /usr/bin/env python
from generic import obj
from nexus import settings,Job
from nexus import generate_physical_system
from nexus import generate_gamess, generate_convert4qmc
from nexus import generate_cusp_correction
from nexus import generate_qmcpack,loop,linear,vmc,dmc

if __name__ == '__main__':
    settings(
        runs          = 'mcscf', 
        status_only   = 0,
        generate_only = 0,
        sleep         = 3,
        machine       = 'ws16',
        ericfmt       = '/opt/intel17/apps/gamess_2016R1/auxdata/ericfmt.dat'
        )
    gms_job = Job(app='rungms',serial=True)


    sims = []

    atom = generate_physical_system(
        type       = 'dimer',
        dimer      = ['B','H'],
        separation = 2.362239,
        units      = 'B',
        net_charge = 0,
        net_spin   = 0
        )

    rhf_inputs = obj(
        identifier = 'rhf',
        system     = atom,
        job        = gms_job,
        scftyp     = 'rhf', 
        runtyp     = 'energy', 
        icharg     = 0, 
        ispher     = 1, 
        mult       = 1, 
        units      = 'bohr',
        mwords     = 200,
        guess      = 'huckel', 
        group      = 'c1', 
        gbasis     = 'cct'
    )
    rhf = generate_gamess(
        path = 'mcscf',
        **rhf_inputs
    )
    sims.append(rhf)

    mcscf_inputs = obj(
        identifier = 'mcscf',
        system     = atom,
        job        = gms_job,
        scftyp     = 'mcscf', 
        runtyp     = 'energy', 
        icharg     = 0, 
        ispher     = 1, 
        mult       = 1, 
        units      = 'bohr',
        mwords     = 200,
        guess      = 'moread', 
        norb       = 69,
        group      = 'c1', 
        ncore      = 0, 
        nels       = 6, 
        nact       = 9,  
        prttol     = 0.001, 
        mcscf      = obj(maxit=2000),
        fullnr     = True,
        memddi     = 100,
        gbasis     = 'cct'
    )
    mcscf = generate_gamess(
        path = 'mcscf',
        dependencies = (rhf,'orbitals'),
        **mcscf_inputs
    )
    sims.append(mcscf)

    mcscf_rerun = generate_gamess(
        path = 'mcscf_rerun',
        **mcscf_inputs
        )
    mcscf_rerun.depends(mcscf,'orbitals')
    sims.append(mcscf_rerun)

    c4q = generate_convert4qmc(
        identifier = 'c4q',
        path = 'c4q',
        job  = Job(app='/home/yyang173/soft/qmcpack/build/bin/convert4qmc',serial=True),
        dependencies = (mcscf_rerun,'orbitals')
    )
    sims.append(c4q)

    cc = generate_cusp_correction(
      identifier = 'cuspcorr',
      path       = 'cusp',
      job        = Job(serial=True,app="qmcpack"),
      system     = atom,
      dependencies = [(c4q,'orbitals')]#,(c4q,'particles')]
    )

    # optimization inputs
    linopt = obj(
	energy               = 0.5,
	unreweightedvariance = 0.0,
	reweightedvariance   = 0.5,
	timestep             = 0.5,
	substeps             = 10,
	warmupsteps          = 50,
	blocks               = 200,
	samples              = 16384,
	usedrift             = True
    )
    opt = generate_qmcpack(
	identifier  = 'opt',
	path        = 'opt',
	job         = Job(app="qmcpack"),
	input_type  = 'basic',
	system      = atom,
	bconds      = 'nnn',    # non-periodic in xyz directions
	jastrows    = [('J1','rcut',3.5),
		       ('J2','rcut',4.5)],
	calculations = [loop(max=3,qmc=linear(**linopt))],
	dependencies = [(c4q,'orbitals'),(cc,'cuspcorr')]
    )
    sims.append(opt)
    
    # dmc inputs
    vmc_input = obj(
	warmupsteps =  40,
	blocks      = 200,
	steps       =  10,
	substeps    =   3,
	timestep    =  .1,
	samples     =  64 # dmc walkers
    )
    dmc_input = obj(
	warmupsteps = 40,
	blocks      = 40,
	steps       = 100,
	timestep    = 0.01
    )
    sdmc = generate_qmcpack(
	identifier  = 'dmc',
	path        = 'dmc',
	job         = Job(app='qmcpack'),
	input_type  = 'basic',
	system      = atom,
	bconds      = 'nnn',    # non-periodic in xyz directions
	jastrows    = [('J1','rcut',3.5), # rcut must be provided for an open system
		       ('J2','rcut',4.5)],
	calculations = [vmc(**vmc_input),dmc(**dmc_input)],
	dependencies = [(c4q,'orbitals'),(c4q,'particles'),(cc,'cuspcorr'),(opt,'jastrow')]
    )
    sims.append(sdmc)

    from nexus import ProjectManager
    pm = ProjectManager()
    pm.add_simulations(sims)
    pm.run_project()

# end __main__
