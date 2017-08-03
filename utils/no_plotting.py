def plot_styles():
  color_map_basis = {
  'double-zeta':'k',
  'triple-zeta':'g',
  'quadruple-zeta':'b'
  }

  ls_map_wf = {
    'Slater':':',
    'Slater-Jastrow':'-'
  }

  marker_map_method = {
    'phfmol':'o',
    'VMC':'^',
    'DMC':'s'
  }

  return color_map_basis, ls_map_wf, marker_map_method
# end def plot_styles

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
    from skimage import measure
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    nx,ny,nz = vol.shape
    lmin,lmax = vol.min(),vol.max()

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
    return mesh
# end def isosurf
