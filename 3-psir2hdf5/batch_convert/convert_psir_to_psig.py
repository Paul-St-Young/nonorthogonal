#!/usr/bin/env python
import numpy as np
import h5py
import matplotlib.pyplot as plt
from pwscf_h5 import PwscfH5
from skimage import measure
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Yubo "Paul" Yang, June 9 2017
#  convert between psi_g and psi_r in pwscf.h5 file from pw2qmcpack
#  note: under write_psir, change tmp_evc to evc and use npw instead of sym_npw
# !!!! only tested on cubic grid with nx points on each dimension

def psig_to_psir(psig,gvec,psir_shape):
    assert len(psig) == len(gvec)
    kgrid = np.zeros(psir_shape,dtype=complex)
    for ig in range(len(gvec)):
        kgrid[tuple(gvec[ig])] = psig[ig]
    # end for
    psir = np.fft.fftshift( np.fft.ifftn(kgrid) ) * np.prod(psir_shape)
    return psir
# end def

def psir_to_psig(psir,gvec):
    kgrid = np.fft.fftn( np.fft.fftshift(psir) )
    psig  = np.zeros(len(gvec),dtype=complex)
    for ig in range(len(gvec)):
        psig[ig] = kgrid[tuple(gvec[ig])]
    # end for
    return psig/np.prod(psir.shape)
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
        level = 0.5*(lmin+lmax)
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

def check_psir_psig(loc,wf,gvec):
    psig_arr = wf.psig(**loc)
    psir_arr = wf.psir(**loc)

    # original psi_r and psi_g
    psig = psig_arr[:,0] + 1j*psig_arr[:,1]
    psir = np.fft.fftshift( psir_arr[:,:,:,0] + 1j*psir_arr[:,:,:,1] )

    # convert psir to psig
    alat   = 3.77945227
    psig_from_psir = psir_to_psig(psir,gvec)*alat**3.
    assert np.allclose( psig_from_psir,psig )

    # convert psig to psir
    psir_from_psig = psig_to_psir(psig,gvec,psir.shape) / alat**3.
    assert np.allclose(psir_from_psig,psir)

    return psir,psir_from_psig
# def

if __name__ == '__main__':
    # 3D image

    wf = PwscfH5()
    wf.read('pwscf.pwscf.h5')
    gvec = wf.get('gvectors')

    # read orbital
    for ikpt in range(8):
        for istate in range(4):
            loc = {'ikpt':ikpt,'ispin':0,'istate':istate}

            pair = check_psir_psig(loc,wf,gvec)

            #"""
            if (ikpt==1) and (istate==0):
                # check visually

                val  = np.absolute(pair[0])
                valr = np.absolute(pair[1])
                print val.min(),val.max()
                print valr.min(),valr.max()
               
                fig = plt.figure()
                ax  = fig.add_subplot(121,projection='3d',aspect=1)
                ax1 = fig.add_subplot(122,projection='3d',aspect=1)
                isosurf(ax,val,0.02)
                isosurf(ax1,valr,0.02)
                plt.show()
                break
            # end if
            #"""

# end __main__
