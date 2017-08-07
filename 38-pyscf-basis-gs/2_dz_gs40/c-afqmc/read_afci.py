#!/usr/bin/env python
import numpy as np
from phfmol_parsing_routines import str2comp

def read_afci(fname,max_nci=1024):
  # parse AFQMC output such as the following:
  #0 old: (0.772346192025,-0.171864042663) new: (-0.762126952533,-0.212646176028)
  #1 old: (-0.151697172987,-0.0952639270359) new: (0.156543176654,-0.0870717041184)
  #2 old: (0.159814060355,-0.08041839808) new: (-0.155316880244,-0.0887936841126)
  #3 old: (0.0357742843283,0.0914831384828) new: (-0.0405831015317,0.0894537804596)
  #4 old: (-0.0349066361869,0.0872319007228) new: (0.0302238530109,0.0889628936522)
  #5 old: (0.099621817391,0.00529909313764) new: (-0.0997626462114,0)

  old_ci_list = []
  new_ci_list = []

  with open(fname,'r+') as f:
    for line in f:
      if 'old' in line and 'new' in line:
        itext,old,ci_old,new,ci_new = line.split()
        assert old == 'old:'
        assert new == 'new:'
        old_ci_list.append(str2comp(ci_old))
        new_ci_list.append(str2comp(ci_new))
      # end if
    # end for
  # end with

  assert len(old_ci_list) == len(new_ci_list)
  if len(old_ci_list) >= max_nci:
    raise RuntimeError('may need to increase max_nci')
  # end if

  return np.array(old_ci_list),np.array(new_ci_list)
# end def read_afci

if __name__ == '__main__':

  dets2read = [25,50,100,150,200]#[5,10,25,50,100] # !!!! hard code
  for ndet in dets2read:
    det_fname = 'get%d/out%d' % (ndet,ndet)
    ci_old,ci_new = read_afci(det_fname)
    np.savetxt('afcoeff%d.dat'%ndet,ci_new.view(float))
    
    plot = False
    if plot:
      import matplotlib.pyplot as plt
      fig,ax = plt.subplots(1,1)
      ax.set_xlabel('old CI coefficient')
      ax.set_ylabel('new CI coefficient')
      ax.plot(ci_old.real,ci_new.real,marker='x',mew=1,ls='',label='real part')
      ax.plot(ci_old.imag,ci_new.imag,marker='+',mew=1,ls='',label='imaginary part')
      ax.set_xlim(-0.22,0.22)
      ax.set_ylim(-0.22,0.22)
      ax.legend()
      plt.show()
    # end if plot
  # end for ndet

# end __main__
