#!/usr/bin/env python
import numpy as np
import sys
sys.path.insert(0,'../../utils')
from parsing import str2comp

def read_next_det(mm,prefix='Determinant',max_nao=1024):
    """ Read the next Determinant in file pointed to by memory map mm
    Inputs:
      mm: mmap object
      max_nao: int, maximum number of atominc orbital in the determinant
    Outputs:
      det_idx: int, index of determinant
      det_data: 1D numpy array of size (nmo*nmo,), stores flattened determinant
    Effect:
      mm will point to the end of the parsed determinant
    """

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
      if i >= (max_nao*max_nao-1):
        raise RuntimeError('may need to increase max_nao')
      # emd if
    # end for i

    return det_idx, np.array(det_data,dtype=complex)
# end def read_next_det

def read_header(mm,header_loc=0):
  def val(line,sep='=',pos=-1):
    return line.split(sep)[pos].strip()
  def get_val(name,mm,sep='=',pos=-1):
    idx = mm.find(name)
    if idx == -1:
        raise RuntimeError('%s not found'%name)
    # end if
    mm.seek(idx)
    line = mm.readline()
    return val(line,sep,pos)

  mm.seek(header_loc) # default is to rewind file

  # make sure this file is what I expect
  assert mm.find('FullMO') != -1
  assert mm.find('CMajor') != -1
  orb_type = get_val('TYPE',mm)
  assert orb_type == 'rotated'
  mm.seek(header_loc)
  det_format = get_val('FORMAT',mm)
  assert det_format == 'PHF'
  mm.seek(header_loc)
  uhf = int(get_val('UHF',mm))

  # get number of determinants
  mm.seek(header_loc)
  nci = int(get_val('NCI',mm))

  # rewind
  mm.seek(header_loc)

  return uhf,nci
# end def read_header

def parse_determinants(fname,nmo):
  """ read a list of determinants in the MO basis from file fname
  Inputs:
    fname: str, name of file that contains all the determinants 
    nmo: int, number of MOs, each determinant should be (nmo,nmo) 
  Outputs:
    det_list: 3D numpy array of shape (nspin,ndet,nmo,nmo), for RHF nspin=1, for UHF nspin=2
  """

  from mmap import mmap
  with open(fname,'r+') as f:
    mm = mmap(f.fileno(),0)
  # end with

  uhf,ndet = read_header(mm)
  nspin = 1 # assume RHF
  if uhf: # if UHF, then two determinants per iteration
    nspin = 2 
  # end if uhf

  det_list = np.zeros([ndet,nspin,nmo,nmo],dtype=complex)
  for idet in range(ndet):
    det_idx, mat_vec = read_next_det(mm)
    my_nmo_sq = len(mat_vec) 
    if uhf:
      my_nmo_sq /= 2
    # end if
    my_nmo = int(round( np.sqrt(my_nmo_sq) ))
    if my_nmo != nmo:
      raise RuntimeError('wrong number of MOs, nmo=%d given in argument, but %d in determinant' % (nmo,my_nmo))
    # end if
    
    my_det = mat_vec.reshape(nspin,nmo,nmo,order='F') # phfmol outputs in column major, checked in read_header
    det_list[idet,:,:,:] = my_det.copy() 
  # end for

  return det_list
# end def parse_determinants

if __name__ == '__main__':
  nmo = 116 # number of basis functions
  fname = 'hdf5/determinants1.det'
  det_list = parse_determinants(fname,nmo)
  np.savetxt('det_list.dat',det_list.flatten().view(float))
