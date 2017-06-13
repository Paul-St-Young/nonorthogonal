#!/usr/bin/env python
import numpy as np

if __name__ == '__main__':

    aoR = np.loadtxt('../4_aoR/aoR.dat')
    moC_data = np.loadtxt('../5_moC/moC.dat')
    moC = moC_data.reshape(8,10,10)

    # LCAO using matrix multiplication
    data = np.zeros([8,729,10])
    for ikpt in range(8):
        dot_moRs = np.dot(aoR,moC[ikpt])
        data[ikpt] = dot_moRs.copy()
        for imo in range(5):
            dot_moR  = dot_moRs[:,imo]

            # check LCAO
            chk_moR = np.zeros(aoR.shape[0])
            for iao in range(aoR.shape[1]):
                chk_moR += aoR[:,iao]*moC[ikpt,iao,imo]
            # end for 

            assert np.allclose(dot_moR,chk_moR)
    np.savetxt('moR.dat',data.reshape(8,7290))

# end __main__
