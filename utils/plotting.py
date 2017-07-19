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
