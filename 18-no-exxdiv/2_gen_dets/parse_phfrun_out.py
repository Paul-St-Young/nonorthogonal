#!/usr/bin/env python
import numpy as np

def read_phfrun_det(mm,mo_header_idx,nbas,real_or_imag):

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
# end def read_phfrun_det

def parse_phfrun_out(phfrun_out,require_idet
 ,success_line = 'The minimization terminated without errors.'
 , mo_coeff_line = 'MO coefficients'):
  from mmap import mmap
  with open(phfrun_out,'r+') as f:
    mm = mmap(f.fileno(),0)
  # end with

  # read number of basis to allocate memory for determinant
  idx = mm.find('nbas')
  mm.seek(idx)
  nbas_line = mm.readline()
  nbas = int(nbas_line.split('=')[-1])
  det_imag = np.zeros([nbas,nbas])

  # check that phfmol.x succeeded
  idx = mm.find(success_line)
  assert idx!=-1
  mm.seek(idx)

  # find new determinant
  idx = mm.find(mo_coeff_line)
  if idx == -1:
    raise RuntimeError('no new determinant found, is lnewdt=.true.?')
  # end if
  mm.seek(idx)

  # 1. read real part
  idet,det_real = read_phfrun_det(mm,idx,nbas,'real')
  assert idet == require_idet
  
  # 2. read imaginary part
  idx = mm.find(mo_coeff_line) # find next MO coefficient line
  idet,det_imag = read_phfrun_det(mm,idx,nbas,'imag')
  assert idet == require_idet

  det = det_real + 1j*det_imag
  return det
# end def parse_phfrun_out

def compare_with_dets():
  phfrun_out = 'phfrun.out.20'
  
  det = parse_phfrun_out(phfrun_out,21)
  #np.savetxt('mydet.21',det.view(float))

  from read_dets import save_det
  idet = 20
  save_det(idet)
  det_path = 'coeff%d/detlist_conj0.dat' % idet
  detlist = np.loadtxt(det_path).view(complex)
  assert len(detlist)==21
  det_ref = detlist[-1].reshape(det.shape).T

  print( np.linalg.norm(det-det_ref) )

  import matplotlib.pyplot as plt
  fig,ax_arr = plt.subplots(1,2)
  ax0 = ax_arr[0]
  ax1 = ax_arr[1]
  ax0.matshow(np.absolute(det_ref))
  ax1.matshow(np.absolute(det))
  plt.show()
# end def compare_with_dets

def fake_detlist(ndet):
  detlist = []
  for idet in range(ndet):
    phfout = 'phfrun.out.%d' % idet
    det = parse_phfrun_out(phfout,idet+1)
    detlist.append(det.flatten())
  # end for idet
  return np.array(detlist)
# end def fake_detlist

if __name__ == '__main__':
  #compare_with_dets()
  mydetlist = fake_detlist(21)
  np.savetxt('out_detlist.dat',mydetlist.view(float))

  newlist = np.loadtxt('out_detlist.dat').view(complex)
  assert np.allclose(mydetlist,newlist)

# end __main__
