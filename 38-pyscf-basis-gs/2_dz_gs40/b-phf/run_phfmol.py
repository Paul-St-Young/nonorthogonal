#!/usr/bin/env python
import os

def phfmol_inp_template():
  ftext = '''&ham_data
  nbas   = {nbas:d},
  fdet   = 'phfrun.det',
  lhdf   = .true.,
  fints  = '{fints:s}' /
!
&nel_data
  nup    = {nup:d},
  ndn    = {ndn:d} /
!
&mtd_data
  ndet_max = {ndet:d},
  wfntyp = '{wfntyp:s}',
  method = 'hf' /
!
&opt_data
  iprint = 3, ! print hamiltonian and overlap matrices
  igstyp = {igstyp:d},
  igsnmx = {igsnmx:d},
  igsmix = {igsmix:d},
  gtol   = {gtol:f},
  maxit  = {maxit:d} /
!
&prj_data
  imult  = {imult:d},
  ngrdb  = {ngrdb:d} /
!
&stt_data
  lnewdt = {lnewdt:s} /
!'''
  return ftext
# end def phfmol_inp_template

def phfmol_submission_template():
  ftext =  '''#!/bin/bash
#SBATCH --job-name gen_dets
#SBATCH -N 4
#SBATCH -t 30
#SBATCH -p pdebug

BIN=~/soft/phfmol/phfmol.x
export OMP_NUM_THREADS=36

cp ../phfrun.det .
srun -N 4 $BIN
  '''
  return ftext
# end def phfmol_submission_template


def get_integrals_from_chkfile(chkfile,
  fcidump_prefix = 'fcidump',
  orthoAO = True, 
  wfnName = 'wfn.dat',
  wfnPHF  = 'phfrun.det'):
  from mpi4py import MPI
  from my_qmctools import new_integrals_from_chkfile

  # get rank and nproc
  comm  = MPI.COMM_WORLD
  rank  = comm.Get_rank()
  nproc = comm.Get_size()

  # get integrals using many MPI groups
  new_integrals_from_chkfile.eri_to_h5(fcidump_prefix, rank, nproc, chkfile
    ,orthoAO=orthoAO ,wfnName=wfnName ,wfnPHF=wfnPHF)

  # collect integrals from the MPI groups
  comm.Barrier()
  if rank == 0:
    new_integrals_from_chkfile.combine_eri_h5(fcidump_prefix, nproc)
  # end if

  return nproc
# end def

if __name__ == '__main__':
  # using PySCF chkfile, construct many non-orthogonal determinants
  scf_dir = '../a-uhf'
  chkfile = os.path.join(scf_dir,'gs40.h5')
  fockh5  = os.path.join(scf_dir,'fock.h5')
  ndet    = 10
  wfntyp  = 'uhf'
  igstyp  = 0 # new determinant guess type - 0=read file phfrun.det
  igsnmx  = 3 # number of orbitals to mix
  igsmix  = 3 # strength of mixing - 3=medium
  gtol    = 1e-4
  maxit   = 1000
  imult   = 1  # spin multiplicity, default = 2*(nup-ndn)+1
  ngrdb   = 16 # number of grid points for beta integration
  lnewdt  = '.true.' # add a new determinant to expansion

  import os
  import h5py
  import subprocess as sp
  from datetime import datetime
  from pyscf.pbc import scf

  if not os.path.isfile(chkfile):
    raise RuntimeError('chkfile %s does not exist. Run PySCF first.')
  # end if

  # construct mf from chkfile
  cell,scf_rec = scf.chkfile.load_scf(chkfile)
  mf = scf.RHF(cell)
  mf.__dict__.update(scf_rec)

  # get integrals from chkfile containing hcore and fock
  fci_prefix = 'fcidump'
  dump_fname = '%s.h5'%fci_prefix
  if not os.path.isfile(dump_fname):
    print('dumping integrals to file ...')
    print(datetime.now())
    nproc = get_integrals_from_chkfile(fockh5,fcidump_prefix=fci_prefix)
    print(datetime.now())
    print('used %d processor(s), dump done.' % nproc)
  # end if

  inp_name = 'phfrun.inp'
  sub_name = 'dets.sbatch'
  phf_dir  = 'dets'
  if not os.path.isdir(phf_dir):
    sp.check_call(['mkdir',phf_dir])

    # write phfmol input
    nup,ndn = mf.cell.nelec
    ftext = phfmol_inp_template()
    inp_content = ftext.format(
      nbas  = mf.mo_coeff.shape[-1],
      fints = os.path.join('../',dump_fname), # !!!! hard code one dir up
      nup   = nup,
      ndn   = ndn,
      ndet  = ndet,
      wfntyp= wfntyp,
      igstyp= igstyp,
      igsnmx= igsnmx,
      igsmix= igsmix,
      gtol  = gtol,
      maxit = maxit,
      imult = imult,
      ngrdb = ngrdb,
      lnewdt= lnewdt
    )

    # write phfmol submission file
    sub_content = phfmol_submission_template()

    # send inputs
    phf_inp_loc = os.path.join(phf_dir,inp_name)
    phf_sub_loc = os.path.join(phf_dir,sub_name)
    with open(phf_inp_loc,'w') as f:
      f.write(inp_content)
    # end with
    with open(phf_sub_loc,'w') as f:
      f.write(sub_content)
    # end with
  # end if

# end __main__
