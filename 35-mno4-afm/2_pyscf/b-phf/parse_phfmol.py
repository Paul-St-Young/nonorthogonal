#!/usr/bin/env python
import numpy as np
from mmap import mmap

def single_det_energies(phfrun_out,ndet):
  from phfmol_parsing_routines import all_idx_with_label, read_phfrun_det
  with open(phfrun_out,'r+') as f:
    mm = mmap(f.fileno(),0)
  # end with

  # get hamiltonian and overlap matrices in the determinant basis
  hmmt = read_phfrun_det(mm,'hmmt matrix',ndet) 
  ovmt = read_phfrun_det(mm,'ovmt matrix',ndet) 

  # imaginary parts should be zero
  zvec = np.zeros(ndet)
  assert np.allclose( np.diag(hmmt.imag), zvec )
  assert np.allclose( np.diag(ovmt.imag), zvec )
  
  diag = np.diag(hmmt.real/ovmt.real)
  return diag
# end def single_det_energies

def get_diagonal(ham_fname,save_diff=True):
  # get single determinant energies from old phfmol output
  from mmap import mmap
  with open(ham_fname,'r+') as f:
    mm = mmap(f.fileno(),0)
  # end with
  from phfmol_parsing_routines import read_next_det

  idx0,hmmt0 = read_next_det(mm,prefix="Hamiltonian")
  # read the second set of hmmt and ovmt
  idx,hmmt = read_next_det(mm,prefix="Hamiltonian")
  idx,ovmt = read_next_det(mm,prefix="Overlap")

  nx = int(np.ceil(np.sqrt(len(hmmt))))

  hmat = hmmt.reshape(nx,nx,order='F')
  omat = ovmt.reshape(nx,nx,order='F')

  diag = hmat.diagonal()/omat.diagonal()
  diff = diag-diag[0]
  assert np.allclose(diff.imag,0.0)
  if save_diff:
    np.savetxt('diag.dat',diff.real,fmt='%6.4f')
  # end if
  return diag
# end def get_diagonal

if __name__ == '__main__':
  # create det_list.dat by parsing phfmol input and outputs
  import os
  from mmap import mmap
  from phfmol_parsing_routines import get_val, parse_determinants

  det_dir     = 'dets'
  phf_infile  = os.path.join(det_dir,'phfrun.inp')        # phf input file, needed to read ndet & nmo
  phf_outfile = os.path.join(det_dir,'phfrun.out')        # phf output file, needed to read hmmt & ovmt
  dets_fname  = os.path.join(det_dir,'determinants1.det') # file containing the determinants

  # read nbas & ndet from phf input
  with open(phf_infile,'r+') as f:
    mm = mmap(f.fileno(),0)
  # end with
  nmo  = int( get_val('nbas',mm).strip(',') )
  ndet = int( get_val('ndet_max',mm).strip(',') )

  # read determinant and create det_list.dat
  det_list_fname = 'det_list.dat'
  if not os.path.isfile(det_list_fname):
    import numpy as np
    det_list = parse_determinants(dets_fname,nmo)
    np.savetxt(det_list_fname,det_list.view(float))
  # end if

  # read hamiltonian and overlap matrices in determinant basis
  diag = single_det_energies(phf_outfile,ndet)
  diff = diag - diag[0]
  np.savetxt('diag.dat',diff)

# end __main__
