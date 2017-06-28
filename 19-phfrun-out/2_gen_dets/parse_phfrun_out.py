#!/usr/bin/env python
import numpy as np

def read_phfrun_det_part(mm,mo_header_idx,nbas,real_or_imag):
  """ read either the real or the complex part of the determinant printed in phfrun.out
  Inputs: 
    mm: mmap of phfrun.out
    mo_header_idx: memory index of the beginning of ther determinant
    nbas: number of basis functions, which determine the shape of the determinant (nbas,nbas)
    real_or_imag: must be either 'real' or 'imag' 
  Outputs: 
    idet: index of the determinant read from the mo_header line
    mydet: either the real or the imaginary part of the determinant 
  """

  if (real_or_imag != 'real') and (real_or_imag != 'imag'):
    raise InputError('real_or_imag must be one of "real" or "imag"')
  # end if

  mydet = np.zeros([nbas,nbas])

  mm.seek(mo_header_idx)
  mo_header = mm.readline()
  ridx = mo_header.index(real_or_imag)
  assert mo_header[ridx:ridx+len(real_or_imag)] == real_or_imag
  didx = mo_header.index('det')
  tokens = mo_header[didx:].split('=')
  new_tokens = tokens[1].split(']')
  idet = int(new_tokens[0])

  col_line = mm.readline()
  col_idx  = np.array(col_line.split(),dtype=int) -1 # -1 to start index at 0

  nblock = int( np.ceil(nbas/len(col_idx)) )
  first_block = True
  for iblock in range(nblock):
    if first_block: # already read col_idx
      first_block = False
    else:
      col_line = mm.readline()
      col_idx = np.array(col_line.split(),dtype=int) -1 # -1 to start index at 0
    # end if
    for ibas in range(nbas):
      row_line   = mm.readline()
      row_tokens = row_line.split()
      irow = int(row_tokens[0]) -1 # -1 to start index at 0
      row = np.array(row_tokens[1:],dtype=float)
      mydet[irow,col_idx] = row.copy()
    # end ibas
  # end iblock

  return idet,mydet
# end def read_phfrun_det_part

def read_next_phfrun_det(mm,cur_idx,nbas,mo_coeff_line = 'MO coefficients'):
  mm.seek(cur_idx)
  # find new determinant
  idx1 = mm.find(mo_coeff_line)
  if idx1 == -1:
    raise RuntimeError('no new determinant found, is lnewdt=.true.?')
  # end if

  # 1. read real part
  idet,det_real = read_phfrun_det_part(mm,idx1,nbas,'real')
  
  # 2. read imaginary part
  idx2 = mm.find(mo_coeff_line) # find next MO coefficient line
  idetm,det_imag = read_phfrun_det_part(mm,idx2,nbas,'imag')

  # ensure real and imaginary parts are read for the same determinant
  assert idetm == idet

  det = det_real + 1j*det_imag
  return idet, det
# end def

def parse_phfrun_out(phfrun_out,require_idet
  ,success_line = 'The minimization terminated without errors.'):
  from mmap import mmap
  with open(phfrun_out,'r+') as f:
    mm = mmap(f.fileno(),0)
  # end with

  # read number of basis to allocate memory for determinant
  idx = mm.find('nbas')
  mm.seek(idx)
  nbas_line = mm.readline()
  nbas = int(nbas_line.split('=')[-1])

  # check that phfmol.x succeeded
  idx = mm.find(success_line)
  assert idx!=-1
  mm.seek(idx)

  idet,det = read_next_phfrun_det(mm,idx,nbas)
  assert idet == require_idet

  return det
# end def parse_phfrun_out

def first_two_phf_dets(phfrun_out):
  """ read the first two determinants in phfrun_out, check if they are equal """
  from mmap import mmap
  with open(phfrun_out,'r+') as f:
    mm = mmap(f.fileno(),0)
  # end with

  # read number of basis to allocate memory for determinant
  idx = mm.find('nbas')
  mm.seek(idx)
  nbas_line = mm.readline()
  nbas = int(nbas_line.split('=')[-1])

  # read 2 determinants from output
  idet1, det1 = read_next_phfrun_det(mm,idx,nbas)
  idet2, det2 = read_next_phfrun_det(mm,mm.tell(),nbas)

  return det1,det2
# end def first_two_phf_dets

def compare_with_manual(ndet=20,nbas=48):
  """ compare formatted output with manual output """

  # 1. read last determinant from phfrun_out
  mydetlist = fake_detlist(ndet)

  # 2. read last determinant from dets_out
  from read_dets import get_multi_determinants
  ci_coeff,detlist = get_multi_determinants('dets.%d'%ndet)

  # 3. compare
  for idet in range(ndet):
    det = mydetlist[idet].reshape(nbas,nbas)
    det_manual = detlist[idet].reshape(nbas,nbas)
    print np.linalg.norm(det-det_manual.T)
# end def

def fake_detlist(ndet):
  detlist = []
  for idet in range(ndet):
    phfout = 'phfrun.out.%d' % idet
    det = parse_phfrun_out(phfout,idet+1)
    detlist.append(det.flatten())
  # end for idet
  return np.array(detlist)
# end def fake_detlist

def get_diagonal(ham_fname,save_diff=True):
  from mmap import mmap
  from read_dets import read_next_det
  with open(ham_fname,'r+') as f:
    mm = mmap(f.fileno(),0)
  # end with

  idx0,hmmt0 = read_next_det(mm,prefix="Hamiltonian")
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

  ## this test demonstrates that determinants are not edited in final_diag
  #for idet in range(20):
  #  phfout = 'phfrun.out.%d' % idet
  #  det1,det2 = first_two_phf_dets(phfout)
  #  assert np.allclose(det1,det2)
  ## end for

  # this test shows that phfrun.out and dets.* contain the same determinants
  #compare_with_manual()

  mydetlist = fake_detlist(21)
  np.savetxt('out_detlist.dat',mydetlist.view(float))

  newlist = np.loadtxt('out_detlist.dat').view(complex)
  assert np.allclose(mydetlist,newlist)

  imag = np.linalg.norm(mydetlist.imag,axis=1)
  np.savetxt('imag.dat',imag)

  diag = get_diagonal('mat.out.20',save_diff=True)

# end __main__
