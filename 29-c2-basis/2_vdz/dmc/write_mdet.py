#!/usr/bin/env python
import numpy as np
from lxml import etree

def read_fortran_array(fname):
  ci_list = []
  with open(fname,'r') as f:
    for line in f:
      tokens = line.strip(' ()\n').split(',')
      real,imag = map(float,tokens)
      ci_list.append(real+1j*imag)
    # end for
  # end with
  ci_arr = np.array(ci_list)
  return ci_arr
# end def

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
  sponode = etree.Element('sposet_builder',{'type':'bspline','href':h5_file,'tilematrix':'1 0 0 0 1 0 0 0 1','twistnum':'0','source':'ion0','version':'0.10','meshfactor':'1.0','fftgrid':'65 65 65','precision':'double'})
  sponode.append(etree.Element('sposet',{'type':'bspline','name':'spo-up','size':str(nstate),'spindataset':'0'}))
  sponode.append(etree.Element('sposet',{'type':'bspline','name':'spo-dn','size':str(nstate),'spindataset':'0'}))
  return sponode
# end def sposet_builder

def save_mdet(fname,ci_coeff,h5_file,nfill):
  ndet = len(ci_coeff)
  nstate = ndet*nfill

  # build <sposet_builder>
  spo_node = sposet_builder(h5_file,nstate)

  # build <determinantset>
  multidet_node = multideterminant_from_ci(ci_coeff,nfill,nstate)
  detset_node = etree.Element('determinantset')
  detset_node.append(multidet_node)

  wf = etree.Element('wavefunction',{'name':'psi0','target':'e'})
  wf.append(spo_node)
  wf.append(detset_node)
  doc = etree.ElementTree(wf)
  doc.write(fname,pretty_print=True)
# end def save_mdet
