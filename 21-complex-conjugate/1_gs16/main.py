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
  if os.path.isfile(wf_h5_fname):
    return # already done

  from pyscf_orbital_routines import save_multideterminant_orbitals, generate_pwscf_h5
  detlist = np.loadtxt(detlist_fname).view(complex)
  gvecs,eig_df = save_multideterminant_orbitals(detlist,nfill,mf)
  generate_pwscf_h5(mf.cell,gvecs,eig_df,pseudized_charge={'C':2},h5_fname=wf_h5_fname)

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
  ndet = 5             # number of determinants in expansion
  det_dir = 'gen_dets' # folder to store determinants
  nfill = 4            # number of occupied orbitals in each determinant
  grid_shape = (16,16,16) # shape of real-space grid

  # execute steps
  # =============================
  from datetime import datetime
  print(datetime.now())
  print('performing Hatree-Fock calculation')
  mf = step1_run_pyscf(grid_shape)          # generate pyscf checkpoint file: vdz.h5
  print(mf.e_tot)
  print(datetime.now())
  step2_dump_integral_table(mf)             # generate integral table: fcidump.dat
  step3_generate_determinants(ndet,det_dir) # generate determinants using phfmol 
  # make sure determinants are generated
  ndone = check_step3(ndet,det_dir)
  if ndone != ndet-1:
    raise RuntimeError('Still running step 3: %d/%d determinants done' % (ndone+1,ndet)) 
  # end if

  step4_read_determinants(ndet,det_dir)     # read determinants and expansion coefficients: detlist.dat, ci_coeffs.dat
  assert os.path.isfile('det_list.dat')

  print(datetime.now())
  print('Fourier transform to generate wavefunction file')
  step5_generate_qmcpack_wavefunction_file('det_list.dat','dets.h5',nfill,mf) # read determinants and create hdf5 file containing all orbitals: dets.h5
  print(datetime.now())

  #from step1_run_pyscf import build_carbon_cell
  #cell = build_carbon_cell(grid_shape,verbose=3)
  step6_write_qmcpack_input('msd.xml',mf.cell,'dets.h5',4,4) # write msd.xml


# end __main__
