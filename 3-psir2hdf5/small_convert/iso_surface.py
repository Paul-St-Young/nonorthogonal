#!/usr/bin/env python
import numpy as np
import h5py
import matplotlib.pyplot as plt
from pwscf_h5 import PwscfH5

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
    psir = psir_arr[:,:,:,0] + 1j*psir_arr[:,:,:,1]

    nx,ny,nz = psir.shape
    verts, faces, normals, values = measure.marching_cubes_lewiner(
        (psir.conj()*psir).real,0.0012) # min: -0.0304, max: 0.03562

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
