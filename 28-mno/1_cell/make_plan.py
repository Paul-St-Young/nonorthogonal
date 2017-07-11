import pickle
from JobEnsemble import JobEnsemble
import Recipes
from PySCFRunner import LocalPySCFRunner
from QWalkRunner import LocalQWalkRunner

if __name__ == '__main__':

  # pseudopotential library
  bfd_lib_xml = '/home/yang41/soft/autogenv2/BFD_Library.xml'

  # read structure
  mycif = open('mno.cif','r').read()

  # make a list of jobs to run
  pyscf_opts = {'cif':mycif,'method':'RHF','kpts':[1,1,1],'gs':[4,4,4],'bfd_library':bfd_lib_xml,'spin':0}
  sim = Recipes.PySCFQWalk('rhf-mno',
    pyscf_opts    = pyscf_opts,
    variance_opts = {'iterations':5},
    energy_opts   = {'total_nstep':8192},
    dmc_opts      = {'timesteps':[0.03]},
    pyscfrunner   = LocalPySCFRunner(),
    qwalkrunner   = LocalQWalkRunner()
  )

  # pickle the plan
  plan = JobEnsemble([sim])
  with open('plan.pickle','wb') as fout:
    pickle.dump(plan,fout)

# end __main__
