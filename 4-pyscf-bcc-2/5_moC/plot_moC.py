#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    moC0 = np.loadtxt('../2_orb/moC.dat')
    nao  = len(moC0)
    moC_arr = np.loadtxt('moC.dat')
    nk   = len(moC_arr)

    ikpt = 0
    imo  = 0

    # kpoint selection
    moC = moC_arr[ikpt].reshape(nao,nao)

    # mo selection
    ref = moC0[:,imo]
    kmo = moC[:,imo]
    myx = range(len(ref))

    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.plot(myx,ref,'o')
    ax.plot(myx,kmo)

    plt.show()

# end __main__
