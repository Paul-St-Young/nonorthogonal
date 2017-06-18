#!/usr/bin/env python
import numpy as np
from lxml import etree
from input_xml import InputXml

def save_vps(fxml,fname_fmt='vps-n{n:d}-{l:s}.dat'):
  inp  = InputXml()
  doc  = etree.parse(fxml)
  root = doc.getroot()
  pots = root.findall('.//vps')

  for pot in pots:
    vps_attribs = inp.node2dict(pot)
    myn = int(vps_attribs['principal-n'])
    myl = vps_attribs['l']
    fname = fname_fmt.format(n=myn,l=myl)

    fr_node  = pot.find('.//radfunc')
    entry = inp.radial_function(fr_node)
    if entry['type'] != 'linear':
      raise NotImplementedError('only implemented linear grid')
    # end if
    if entry['units'] != 'bohr':
      raise NotImplementedError('hard-coded for bohr')
    # end if
    myx = np.linspace(entry['ri'],entry['rf'],entry['npts'])
    myy = entry['rval']
    data = np.array([myx,myy]).T
    np.savetxt(fname,data)
  # end for pot
# end def save_vps

if __name__ == '__main__':

  #save_vps('C.BFD.xml')

  import matplotlib.pyplot as plt
  fig,ax = plt.subplots(1,1)

  for myl in ['s','p']:
    fname = 'vps-n0-%s.dat'%myl
    data  = np.loadtxt(fname)
    myx,myy = data.T
    ax.plot(myx,myy,label='%s'%myl)
  # end for
  ax.legend(loc=0)
  plt.show()
# end __main__
