#!/usr/bin/env python
import numpy as np
import pandas as pd
import h5py

def build_cell(gs,verbose=3):
  from pyscf.pbc import gto
  #alat = 2.34301125 # angstrom !!!! wrong
  alat = 4.427649569275 # angstrom
  pseudo = {'Mn':'bfd','O':'bfd'}
  from pyscf.gto.basis import parse
  basis={
  'Mn':parse('''
  Mn s
    23.6450680 -0.014659 
    13.4620290 0.223661 
    8.21262900 -0.564535 
    1.85994300 0.554600 
    0.88975400 0.534536 
    0.42256000 0.142766 
   Mn p
    17.8033170 0.003171 
    8.83499500 -0.080060 
    4.16134000 0.185735 
    2.12919600 0.375367 
    1.06371300 0.388624 
    0.51014100 0.192450 
    0.21130200 0.026146 
   Mn d
    9.25373600 0.076285 
    4.36301600 0.247722 
    1.83482100 0.330050 
    0.75611800 0.334518 
    0.29609200 0.263545 
   Mn f
    1.0998003483 1.00000 
   Mn s
  0.2 1.0
   Mn s
  0.6000000000000001 1.0
   Mn p
  0.2 1.0
   Mn p
  0.6000000000000001 1.0
   Mn d
  0.2 1.0
   Mn d
  0.6000000000000001 1.0
  '''),
  'O':parse('''
  O s
    0.268022 0.304848 
    0.573098 0.453752 
    1.225429 0.295926 
    2.620277 0.019567 
    5.602818 -0.128627 
    11.980245 0.012024 
    25.616801 0.000407 
    54.775216 -0.000076 
   O p
    0.333673 0.255999 
    0.666627 0.281879 
    1.331816 0.242835 
    2.660761 0.161134 
    5.315785 0.082308 
    10.620108 0.039899 
    21.217318 0.004679 
   O d
    0.669340 1.000000 
   O f
    1.423104 1.000000 
   O s
  0.2 1.0
   O s
  0.6000000000000001 1.0
   O p
  0.2 1.0
   O p
  0.6000000000000001 1.0
  '''),
  }

  axes_text = '''
   {alat:12.10f}  {half:12.10f}  {half:12.10f}
   {half:12.10f}  {alat:12.10f}  {half:12.10f}
   {half:12.10f}  {half:12.10f}  {alat:12.10f}
  '''.format(alat=alat,half=alat/2.)
  atom_text = '''
   Mn      0.00000      0.00000      0.00000
   Mn      {alat:12.10f} {alat:12.10f} {alat:12.10f}
   O       {half:12.10f} {half:12.10f} {half:12.10f}
   O       {half3:12.10f} {half3:12.10f} {half3:10.6}
  '''.format(alat=alat,half=alat/2.,half3=3.*alat/2.)

  cell = gto.M(
    verbose = verbose,
    gs      = gs,
    a       = axes_text,
    atom    = atom_text,
    basis   = basis,
    pseudo  = pseudo,
    unit    = 'angstrom'
  )
  return cell
# end def build_cell

def run_pyscf(gs,verbose=3,chkfile_name='bfd.h5'):
  import os
  from pyscf.pbc import scf
  #from mpi4pyscf.pbc import df as mpidf

  cell = build_cell(gs,verbose)
  kpt = cell.get_abs_kpts([.0,.0,.0])  # gamma point calculation
  mf = scf.UHF(cell,exxdiv=None)       # exxdiv=None means no Ewald sum
  #mf = scf.newton(mf)                  # use Newton converger
  #mf.with_df = mpidf.FFTDF(cell,kpt)   # use mpi for FFT

  # throw out linearly dependent basis functions
  #def eig(h, s, eps=1e-10):
  #  d, t = np.linalg.eigh(s)
  #  x = t[:,d>eps] / np.sqrt(d[d>eps])
  #  xhx = reduce(np.dot, (x.T.conj(), h, x))
  #  e, c = np.linalg.eigh(xhx)
  #  c = np.dot(x, c)
  #  return e, c
  ## end def eig
  #mf.eig = eig
  
  # tweak converger
  mf.conv_tol       = 1e-7
  mf.direct_scf_tol = 1e-7
  mf.max_cycle      = 50
  mf.diis_start_cycle = 1

  # run or load
  mf.chkfile = chkfile_name # set restart file
  if os.path.isfile(chkfile_name):
    from pyscf import lib
    mf.__dict__.update(lib.chkfile.load(chkfile_name,'scf'))
  else:
    mf.kernel()
  return mf
# end def run_pyscf

def modify_1rdm(dm,mn_dorb_indices,ndorb=5):
  """ impose anti-ferromagnetic on UHF density matrix
  Inputs:
    dm: 3D numpy array of shape (nspin,nao,nao) - expect nspin=2 for up and down 
    mn_dorb_indices: a list of int, one for each Mn atom - index of the first d orbital
    ndorb: int, number of d orbitals - should be 5 right? 
  Output:
    new_dm: 3D numpy array, copied from dm, then modified """

  nspin,nao,nao1 = dm.shape
  assert nspin == 2
  new_dm = dm.copy()
  for ispin in [0,1]:

    # AFM order
    first_up = True
    if ispin == 0:
      first_up = False

    set_up = first_up
    for ibegin in mn_dorb_indices:
      # set all d electrons
      iend   = ibegin + ndorb
      if set_up:
        new_dm[ispin,ibegin:iend,ibegin:iend] = 1.0 * np.eye(ndorb) 
      else:
        new_dm[ispin,ibegin:iend,ibegin:iend] = 0.0 * np.eye(ndorb) 
      # end if
      set_up = not set_up 
    # end for ibegin

  # end for ispin
  return new_dm
# end def modify_1rdm

def write_uhf_multideterminant_spo(det_list,ndet,nfill,mf):
  import pandas as pd
  from pyscf_orbital_routines import multideterminant_orbitals

  # ndet0: number of available determinants
  ndet0,nspin,nmo,nmo1 = det_list.shape

  # split into up and down determinants to reuse RHF routines
  det_list_up = np.zeros([ndet,nmo,nmo],dtype=complex)
  det_list_dn = np.zeros([ndet,nmo,nmo],dtype=complex)
  for ispin in range(nspin): # up and down
    for idet in range(ndet):
      if ispin == 0:
        det_list_up[idet,:,:] = det_list[idet,ispin,:,:]
      elif ispin == 1:
        det_list_dn[idet,:,:] = det_list[idet,ispin,:,:]
      else:
        raise RuntimeError('cannot handle ispin=%d'%ispin)
      # end if
    # end for idet
  # end for ispin

  gvecs_up, eig_df_up = multideterminant_orbitals(det_list_up,nfill,mf.cell,mf.mo_coeff[0],mf.kpt)
  gvecs_dn, eig_df_dn = multideterminant_orbitals(det_list_dn,nfill,mf.cell,mf.mo_coeff[1],mf.kpt)
  assert np.allclose(gvecs_up,gvecs_dn)
  eig_df = pd.concat([eig_df_up,eig_df_dn])
  return gvecs_up,eig_df
# end def write_uhf_multideterminant_spo

if __name__ == '__main__':

  from datetime import datetime

  skip_scf   = True
  verbosity  = 3
  grid_shape = [25]*3
  ndet       = 10

  # run pyscf to obtain mean-field object
  print('running HF...')
  print(datetime.now())
  mf = run_pyscf(grid_shape,verbose=verbosity)
  print(datetime.now())
  print(' HF done')

  # make sure the SCF cycled is converged to the AFM state
  if not skip_scf:
    if not mf.converged:

      #from pyscf.scf.hf import init_guess_by_minao
      #init_dm = init_guess_by_minao(mf.cell)
      #np.savetxt('minao.dat',init_dm)

      # get current 1RDM
      dm = mf.make_rdm1()
      np.savetxt('1rdm_up.dat',dm[0])
      np.savetxt('1rdm_dn.dat',dm[1])
      # enforce anti-ferromagnetic order
      new_dm  = modify_1rdm(dm,[12,46])
      np.savetxt('new_1rdm_up.dat',new_dm[0])
      np.savetxt('new_1rdm_dn.dat',new_dm[1])
      # rerun scf
      print('rerunning HF...')
      print(datetime.now())
      mf.kernel(new_dm)
      print(' HF done')
      print(datetime.now())
    # end if

    ## check state!
    #labels = mf.cell.ao_labels(fmt=None)
    #with open('basis.txt','w') as fout:
    #  fout.write( ['%s  %s\n' % (label[2],label[3]) for label in labels] )
    ## end with
    print("Mulliken population:")
    mp  = mf.mulliken_pop()

    # save Hcore and Fock matrices for fcidump.h5
    hcore = mf.get_hcore()
    fock = (hcore + mf.get_veff())
    if mf.chkfile:
      with h5py.File(mf.chkfile) as fh5:
        fh5['scf/hcore'] = hcore
        fh5['scf/fock'] = fock
      # end with
    # end if not mf.converged
  # end if skip_scf
  #assert mf.converged

  # ==== run run_int.py to generate fcidump.h5 ==== #
  # ====  then run phfmol.x to generate determinants ==== #
  # ====  finally run read_det_list.py to write 'det_list.dat' ==== #

  # save Kohn-Sham eigensystem to dataframe
  import os
  if not (os.path.isfile('gvectors.dat') and os.path.isfile('eigensystem.json')):
    nspin,nao,nmo = mf.mo_coeff.shape
    nfill = mf.cell.tot_electrons()
    det_list = np.loadtxt('det_list.dat').view(complex).reshape([ndet,nspin,nmo,nmo])
    gvecs,df = write_uhf_multideterminant_spo(det_list,ndet,nfill,mf)
    np.savetxt('gvectors.dat',gvecs)
    df.reset_index().to_json('eigensystem.json')
  else:
    gvecs = np.loadtxt('gvectors.dat')
    df = pd.read_json('eigensystem.json').set_index(['ikpt','ispin','istate'])
  # end if

  # write QMCPACK wavefunction file
  from pyscf_orbital_routines import generate_pwscf_h5
  generate_pwscf_h5(mf.cell,gvecs,eig_df,pseudized_charge={'Mn':10,'O':2})

# end __main__
