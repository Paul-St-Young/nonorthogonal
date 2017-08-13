
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

def isosurf(ax,vol,level_frac=0.5):
    """ draw iso surface of volumetric data on matplotlib axis at given level
    Inputs:
      ax: matplotlib axis with projection='3d' 
      vol: 3D volumetric data as a numpy array (nx,ny,nz) 
      level_frac: float 0.0->1.0 isosurface value as a fraction between min and max
    Output:
      mesh: Poly3DCollection object
    Effect:
      draw on ax """
    from skimage import measure
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    nx,ny,nz = vol.shape
    lmin,lmax = vol.min(),vol.max()

    level = lmin + level_frac*(lmax+lmin)
    if level<lmin or level>lmax:
        raise RuntimeError('level must be >%f and < %f'%(lmin,lmax))
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

def draw_cell(ax,axes,pos,atom_color='b',draw_super=True):
  atoms = []
  dots  = ax.plot(pos[:,0],pos[:,1],pos[:,2],'o',c=atom_color,ms=10)
  atoms.append(dots)
  if draw_super:
    import numpy as np
    from itertools import product
    for ix,iy,iz in product(range(2),repeat=3):
      if ix==iy==iz==0:
        continue
      shift = (np.array([ix,iy,iz])*axes).sum(axis=0)
      spos  = (shift.reshape(-1,1,3) + pos).reshape(-1,3)
      dots  = ax.plot(spos[:,0],spos[:,1],spos[:,2],'o',c='gray',ms=10,alpha=0.8)
      atoms.append(dots)

  # show simulation cell
  cell = []
  for idim in range(3):
    line = ax.plot([0,axes[idim,0]],[0,axes[idim,1]],[0,axes[idim,2]],c='k',lw=2)
    cell.append(line)

  return atoms,cell
# end def
