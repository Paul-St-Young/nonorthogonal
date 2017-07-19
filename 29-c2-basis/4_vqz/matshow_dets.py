#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
  nbas = 100
  det_list = np.loadtxt('det_list.dat').view(complex)
  ndet = len(det_list)
  dets = det_list.reshape(ndet,nbas,nbas)

  ndet2show = 20
  nfill = 4
  fig,ax_arr = plt.subplots(1,ndet2show)
  for idet in range(ndet2show):
    ax = ax_arr[idet]

    if idet != 0:
      ax.get_yaxis().set_visible(False)

    det = dets[idet,:,:nfill]
    img = ax.matshow( np.absolute(det), cmap=plt.get_cmap('bone') )

    ax.set_xticks([])
    ax.set_xlabel('Det. %d'%idet,rotation=45)
    ax.set_ylabel('Hatree-Fock Orbital Index',fontsize=16)

  # end for idet
  cbaxis = fig.add_axes([0.925,0.1,0.01,0.8])
  fig.colorbar(img,cax=cbaxis)
  fig.savefig('c2-qz-%ddets.png'%ndet,dpi=480)
  plt.show()

# end __main__
