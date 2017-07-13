#!/usr/bin/env python
import numpy as np

if __name__ == '__main__':
  import matplotlib.pyplot as plt
  fig,ax = plt.subplots(1,1)
  
  for folder in ['2_vdz','3_vtz','4_vqz']:
    label = folder.split('_')[-1]
    fname = '%s/gen_dets/ke.dat' % folder
    data  = np.loadtxt(fname)
    ax.plot(data[:,0],data[:,1],'.',label=label)
  # end for folder

  ax.legend()
  plt.show()
