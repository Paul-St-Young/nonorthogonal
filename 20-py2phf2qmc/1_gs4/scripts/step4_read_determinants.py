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

def fake_detlist(ndet,det_dir):
  import os
  detlist = []
  for idet in range(ndet):
    phfout = os.path.join(det_dir,'phfrun.out.%d' % idet)
    det = parse_phfrun_out(phfout,idet+1)
    detlist.append(det.flatten())
  # end for idet
  return np.array(detlist)
# end def fake_detlist

def str2comp(text):
  tokens = text.strip(' ()\n').split(',')
  real,imag = map(float,tokens)
  return real+1j*imag
# end def str2comp

def build_ci_coefficients(ndet,det_dir):
  import os
  from mmap import mmap
  coeff_list = []
  for idet in range(ndet):
    phfout = os.path.join(det_dir,'phfrun.out.%d' % idet)
    with open(phfout,'r+') as f:
      mm  = mmap(f.fileno(),0)
    # end with

    idx1 = mm.find('expansion coefficients')
    mm.seek(idx1)
    mm.readline()
    # get the second output of expansion coefficients
    idx2 = mm.find('expansion coefficients')
    mm.seek(idx2)
    mm.readline()

    ci_coeff = []
    for k in range(idet+1):
      line = mm.readline().strip()
      if line.startswith('(') and line.endswith(')'):
        pair = line.split(',')
        if len(pair) == 2:
          ci_coeff.append(str2comp(line))
        elif len(pair) == 3:
          for num in line.split():
            ci_coeff.append(str2comp(num))
          # end for
        else:
          raise RuntimeError('cannot handle %d numbers on a line'%len(pair))
        # end if
      # end if
    # end for k 
    next_line = mm.readline().strip()
    one_more_coeff = next_line.startswith('(') and next_line.endswith(')')
    assert not one_more_coeff
    assert len(ci_coeff) == idet+1
    coeff_list.append(ci_coeff)
  # end for idet
  return coeff_list
# end def build_ci_coefficients

def read_next_det(mm,prefix='Determinant',max_nao=1024):
    """ mm: mmap object
        max_nao: maximum number of atominc orbital in the determinant """

    idx = mm.find(prefix)
    if idx == -1:
      return -1
    mm.seek(idx)

    det_data = []
    det_line = mm.readline()
    det_idx = int( det_line.split(':')[-1] )
    for i in range(max_nao*max_nao):
      cur_idx = mm.tell()
      line = mm.readline().strip()
      if line.startswith('(') and line.endswith(')'):
        pair = line.split(',')
        if len(pair) == 2:
          det_data.append(str2comp(line))
        elif len(pair) == 3:
          for num in line.split():
            det_data.append(str2comp(num))
          # end for
        else:
          raise RuntimeError('cannot handle %d numbers on a line'%len(pair))
        # end if
      else:
        mm.seek(cur_idx)
        break
      # end if
    # end for i

    return det_idx, np.array(det_data,dtype=complex)
# end def read_next_det

def get_diagonal(ham_fname,save_diff=True):
  from mmap import mmap
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

def step4_read_determinants(ndet,det_dir,detlist_fname='det_list.dat',coeff_list_dir='coeff_list'):
  import os
  already_done = os.path.isfile(detlist_fname) and os.path.isdir(coeff_list_dir)
  if not already_done:
    mydetlist = fake_detlist(ndet,det_dir)
    np.savetxt(detlist_fname,mydetlist.view(float))
    coeff_list = build_ci_coefficients(ndet,det_dir)
    if not os.path.isdir(coeff_list_dir):
      os.system('mkdir '+coeff_list_dir)
    # end if
    for k in range(len(coeff_list)):
      ci_file = os.path.join(coeff_list_dir,'ci%d.dat'%k)
      np.savetxt(ci_file,coeff_list[k])
    # end for k
  # end if

  diag = get_diagonal(os.path.join(det_dir,'mat.out.%d'%(ndet-1)),save_diff=True)

  ## to read the determinants:
  #newlist = np.loadtxt('detlist.dat').view(complex)
  #assert np.allclose(mydetlist,newlist)

