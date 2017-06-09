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

def isosurf(ax,vol,level=None):
    """ draw iso surface of volumetric data on matplotlib axis at given level
    Inputs:
      ax: matplotlib axis with projection='3d' 
      vol: 3D volumetric data as a numpy array (nx,ny,nz) 
      level: value of iso surface
    Output:
      None
    Effect:
      draw on ax """

    nx,ny,nz = vol.shape
    lmin,lmax = vol.mean(),vol.max()

    if level is None: # set level to average if none given
        level = 0.5*(lmin*lmax)
    else: # check isosurface level
        if level<lmin or level>lmax:
            raise RuntimeError('level must be >%f and < %f'%(lmin,lmax))
        # end if
    # end if

    # make marching cubes
    verts, faces, normals, values = measure.marching_cubes_lewiner(
        vol, level)

    # plot surface
    mesh = Poly3DCollection(verts[faces])
    mesh.set_edgecolor('k')
    ax.add_collection3d(mesh)
    ax.set_xlim(0,nx)
    ax.set_ylim(0,ny)
    ax.set_zlim(0,nz)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
# end def isosurf

if __name__ == '__main__':
    # 3D image
    from skimage import measure
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection

    # read orbital
    loc = {'ikpt':0,'ispin':0,'istate':0}

    wf = PwscfH5()
    wf.read('pwscf.pwscf.h5')
    gvec = wf.get('gvectors')
    psig_arr = wf.psig(**loc)
    psir_arr = wf.psir(**loc)

    psig = psig_arr[:,0] + 1j*psig_arr[:,1]
    mypsir = psig_to_psir(psig,gvec,12)
    psir = np.fft.fftshift( psir_arr[:,:,:,0] + 1j*psir_arr[:,:,:,1] )

    nx,ny,nz = psir.shape
    valr = (psir.conj()*psir).real # min: -0.0304, max: 0.03562, then squre
    val  = (mypsir.conj()*mypsir).real

    fig = plt.figure()
    ax  = fig.add_subplot(121,projection='3d',aspect=1)
    ax1 = fig.add_subplot(122,projection='3d',aspect=1)

    isosurf(ax,val,1.0)
    isosurf(ax1,valr,1e-3)

    plt.show()

# end __main__
