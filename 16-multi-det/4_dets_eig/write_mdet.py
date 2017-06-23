#!/usr/bin/env python
import numpy as np
from lxml import etree

def occupy(istart,nfill,ntot):
  """ return strings like 11110000, 00001111, which represent occupation of single-particle states """
  occ_arr = ['0'] * ntot
  occ_arr[istart:(istart+nfill)] = ['1'] * nfill
  text = ''.join(occ_arr)
  return text
# end def occupy

def multideterminant_from_ci(ci_coeff,nfill,nstate):
  ndet = len(ci_coeff)
  node = etree.Element('multideterminant',{'optimize':'no','spo_up':'spo-up','spo_dn':'spo-dn'})
  detlist = etree.Element('detlist',{'size':str(ndet),'type':'DETS','nca':'0','ncb':'0','nea':str(nfill),'neb':str(nfill),'cutoff':'1e-16','nstates':str(nstate)})
  node.append(detlist)
  for idet in range(ndet):
    coeff_text = '(%f,%f)' % (ci_coeff[idet].real,ci_coeff[idet].imag)
    alpha = occupy(idet*nfill,nfill,nstate)
    beta = occupy(idet*nfill,nfill,nstate)
    det = etree.Element('ci',{'id':'CIcoeff_%d'%idet,'coeff':coeff_text,'alpha':alpha,'beta':beta})
    detlist.append(det)
  # end for idet
  return node
# end def multideterminant_from_ci

def sposet_builder(h5_file,nstate):
  sponode = etree.Element('sposet_builder',{'type':'bspline','href':h5_file,'tilematrix':'1 0 0 0 1 0 0 0 1','twistnum':'0','source':'ion0','version':'0.10','meshfactor':'1.0','precision':'double'})
  sponode.append(etree.Element('sposet',{'type':'bspline','name':'spo-up','size':str(nstate),'spindataset':'0'}))
  sponode.append(etree.Element('sposet',{'type':'bspline','name':'spo-dn','size':str(nstate),'spindataset':'0'}))
  return sponode
# end def sposet_builder

def main():
  nfill = 4 # !!!! hard-code 4 orbitals per determinant
  ci_coeff = np.loadtxt('../3_parse_dets/ci_coeff.dat').view(complex)
  h5_file = '../4_dets_eig/pyscf2qmcpack.h5'

  ndet = len(ci_coeff)
  nstate = ndet*nfill
  multidet_node = multideterminant_from_ci(ci_coeff,nfill,nstate)
  spo_node = sposet_builder(h5_file,nstate)

  wf = etree.Element('wavefunction',{'name':'psi0','target':'e'})
  wf.append(spo_node)
  wf.append(multidet_node)
  doc = etree.ElementTree(wf)
  doc.write('mdet.xml',pretty_print=True)

if __name__ == '__main__':
  main()
# end __main__
