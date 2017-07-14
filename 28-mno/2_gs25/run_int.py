from qmctools import integrals_from_chkfile
from pyscf.pbc import df
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nproc = comm.Get_size()

#integrals_from_chkfile.eri_to_h5("fcidump", rank, nproc, "bfd.h5") #RHF
integrals_from_chkfile.eri_to_h5("fcidump", rank, nproc, "bfd.h5",orthoAO=True,wfnName="wfn.dat",wfnPHF="phfrun.det")

comm.Barrier()

if rank==0:
  integrals_from_chkfile.combine_eri_h5("fcidump", nproc)
