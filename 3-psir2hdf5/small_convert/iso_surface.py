#!/usr/bin/env python
import numpy as np
import h5py
import matplotlib.pyplot as plt
from pwscf_h5 import PwscfH5

def psig_to_psir(psig,gvec,nx):
    assert len(psig) == len(gvec)
    kgrid = np.zeros([nx,nx,nx])
    for ig in range(len(gvec)):
        kgrid[tuple(gvec[ig])] = psig[ig]
    # end for
    psir = np.fft.fftshift( np.fft.ifftn(kgrid) ) * nx**3.
    return psir
# end def

if __name__ == '__main__':
    # 3D image
    from skimage import measure
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    # read orbital
    wf = PwscfH5()
    wf.read('pwscf.pwscf.h5')
    gvec = wf.get('gvectors')
    psig_arr = wf.psig()
    psir_arr = wf.psir()

    psig = psig_arr[:,0] + 1j*psig_arr[:,1]
    mypsir = psig_to_psir(psig,gvec,12)
    psir = psir_arr[:,:,:,0] + 1j*psir_arr[:,:,:,1]

    nx,ny,nz = psir.shape
    #val = (psir.conj()*psir).real # min: -0.0304, max: 0.03562, then squre
    val = (mypsir.conj()*mypsir).real
    print val.min(),val.max()
    verts, faces, normals, values = measure.marching_cubes_lewiner(
        val ,1.0)

    fig = plt.figure()
    ax  = fig.add_subplot(111,projection='3d')

    mesh = Poly3DCollection(verts[faces])
    mesh.set_edgecolor('k')
    ax.add_collection3d(mesh)

    ax.set_xlim(0,nx)
    ax.set_ylim(0,ny)
    ax.set_zlim(0,nz)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()

# end __main__
