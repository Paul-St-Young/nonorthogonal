#!/usr/bin/env python
import os
import numpy as np

def check_step3(ndet,det_dir):
  import os
  for idet in range(ndet):
    fout = os.path.join(det_dir,'phfrun.out.%d' % idet)
    if not os.path.isfile(fout):
      break
    # end if
  # end for idet
  return idet
# end def

def step5_generate_qmcpack_wavefunction_file(detlist_fname,wf_h5_fname,nfill,mf):

  from pyscf_orbital_routines import save_multideterminant_orbitals, generate_pwscf_h5
  detlist = np.loadtxt(detlist_fname).view(complex)

  # original
  gvecs,eig_df = save_multideterminant_orbitals(detlist,nfill,mf)
  generate_pwscf_h5(mf.cell,gvecs,eig_df,pseudized_charge={'C':2},h5_fname=wf_h5_fname)

# end def step5_generate_qmcpack_wavefunction_file

def step5_1_single_det_h5s(detlist_fname,nfill,mf):
  # split
  ndet = len(detlist)
  for idet in range(ndet):
    dets = np.array([detlist[idet]])
    gvecs,eig_df = save_multideterminant_orbitals(dets,nfill,mf)
    generate_pwscf_h5(mf.cell,gvecs,eig_df,pseudized_charge={'C':2},h5_fname='det%d.h5'%(idet))
  # end for ifile
# end def step5_1_single_det_h5s

def dump_molecular_orbitals(mf):
  from pyscf_orbital_routines import get_pyscf_psir
  norb = 16 # !!!! hard-code write out the first four orbitals
  for iorb in range(norb):
    moR = get_pyscf_psir(mf.mo_coeff[:,iorb],mf.cell)
    np.savetxt('moR%d.dat'%iorb,moR.flatten().view(float))
  # end for iorb
# end def dump_molecular_orbitals

if __name__ == '__main__':
  # import scripts for each step
  # =============================
  import sys
  sys.path.insert(0,'./scripts')
  from step1_run_pyscf import step1_run_pyscf
  from step2_dump_integral_table import step2_dump_integral_table
  from step3_generate_determinants import step3_generate_determinants
  from step4_read_determinants import step4_read_determinants
  from step6_write_qmcpack_input import step6_write_qmcpack_input

  # define parameters
  # =============================
  ndet = 50               # number of determinants in expansion
  det_dir = 'gen_dets'    # folder to store determinants
  nfill = 16              # number of occupied orbitals in each determinant
  grid_shape = (16,16,16) # shape of real-space grid

  # execute steps
  # =============================

  # each step generates a local file. If the required file exists, then skip the corresponding step
  check_points = {
    1:'vdz.h5',
    2:'fcidump.dat',
    # step 3 is complicated (run phfmol), check using a function
    4:'det_list.dat',
    5:'dets.h5',
    6:'msd.xml'
  }

  with_ewald = -44.01179058
  no_ewald   = -37.3387542793
  ewald_correction = no_ewald - with_ewald

  from datetime import datetime

  print('performing Hatree-Fock calculation')
  print(datetime.now())
  mf = step1_run_pyscf(grid_shape)          # generate pyscf checkpoint file: vdz.h5
  print(mf.e_tot)
  print(datetime.now())

  #dump_molecular_orbitals(mf)

  if not os.path.isfile( check_points[2] ):
    step2_dump_integral_table(mf)             # generate integral table: fcidump.dat
  # end if

  if not os.path.isdir( det_dir ):
    step3_generate_determinants(ndet,det_dir) # generate determinants using phfmol 
  # end if

  # make sure determinants are generated
  ndone = check_step3(ndet,det_dir)
  if ndone != ndet-1:
    raise RuntimeError('Still running step 3: %d/%d determinants done' % (ndone+1,ndet)) 
  # end if

  if not os.path.isfile( check_points[4] ) or not os.path.isdir('coeff_list'):
    step4_read_determinants(ndet,det_dir)     # read determinants and expansion coefficients: detlist.dat, coeff_list
    assert os.path.isfile( check_points[4] )
    assert os.path.isdir('coeff_list')
  # end if

  if not os.path.isfile( check_points[5] ):
    try:
      mf
    except:
      mf = step1_run_pyscf(grid_shape)          # generate pyscf checkpoint file: vdz.h5
    # end try
    print('Fourier transform to generate wavefunction file')
    print(datetime.now())
    step5_generate_qmcpack_wavefunction_file(check_points[4],check_points[5],nfill,mf) # read determinants and create hdf5 file containing all orbitals: dets.h5
    #step5_1_single_det_h5s(check_points[4],nfill,mf)
    print(datetime.now())
  # end if


  if not os.path.isfile( check_points[6] ):
    try:
      cell = mf.cell
    except:
      from step1_run_pyscf import build_carbon_cell
      cell = build_carbon_cell(grid_shape,verbose=3)
    # end try
    nup,ndn = cell.nelec
    step6_write_qmcpack_input(check_points[6],cell,check_points[5],nup,ndn,proj_id='c2') # write msd.xml
  # end if

# end __main__
