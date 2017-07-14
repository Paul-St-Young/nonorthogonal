#!/usr/bin/env python
# collect database containing columns:
#  ndet energy_mean energy_error wf method basis
import os
import numpy as np
import pandas as pd

def read_phfmol_ke(fname,ewald):
  phfmol = np.loadtxt(fname)

  mydf = pd.DataFrame()
  mydf['ndet']         = phfmol[:,0].astype(int) + 1
  mydf['energy_mean']  = phfmol[:,1].astype(float) + ewald
  mydf['energy_error'] = 0.0
  mydf['method']       = 'phfmol'
  mydf['wf']           = 'Slater'
  return mydf
# end def read_phfmol_ke

def parse_qmcas_output(fname):
  # parse outputs such as:
  #  2_vdz/mdets/detsci49_af/c2  series 0  -10.452004 +/- 0.006424    1.6   1.601253 +/- 0.068660    1.1   0.1532

  def basis_ndet_from_myid(myid):
    tokens    = myid.split('/')
    basis_str = tokens[-4]
    det_str   = tokens[-2]
    basis = basis_str.split('_')[-1]
    ndet = det_str.replace('detsci','').replace('_af','')
    basis_name_map = {
      'vdz':'double-zeta',
      'vtz':'triple-zeta',
      'vqz':'quadruple-zeta'
    }
    return basis_name_map[basis],int(ndet) + 1
  # end def basis_ndet_from_myid

  data = {'ndet':[],'energy_mean':[],'energy_error':[],'method':[],'iqmc':[],'basis':[]}
  with open(fname) as f:
    for line in f:
      tokens = line.split()
      if len(tokens) == 12:
        # determine basis,ndet
        myid   = tokens[0]
        basis,ndet = basis_ndet_from_myid(myid)

        # determine method
        iqmc = int(tokens[2])
        if iqmc == 0:
          method = 'VMC'
        else:
          method = 'DMC'
        # end if

        # read energy
        energy_mean  = float(tokens[3])
        energy_error = float(tokens[5])

        data['ndet'].append(ndet)
        data['basis'].append(basis)
        data['method'].append(method)
        data['iqmc'].append(iqmc)
        data['energy_mean'].append(energy_mean)
        data['energy_error'].append(energy_error)
      # end if
    # end for 
  # end with
  mydf = pd.DataFrame(data)
  return mydf
# end def parse_qmcas_output

def collect_all():
  ewald = -2.695783 # ewald correction

  folder_map = {
    'double-zeta':'2_vdz',
    'triple-zeta':'3_vtz',
    'quadruple-zeta':'4_vqz'
  }

  data = []

  # collect phfmol results
  phfmol_dir = 'gen_dets'
  phfmol_efile = 'ke.dat'
  for basis in ['double-zeta','triple-zeta','quadruple-zeta']:
    subdir = folder_map[basis]

    fname = os.path.join(subdir,phfmol_dir,phfmol_efile)
    mydf  = read_phfmol_ke(fname,ewald=ewald)
    mydf['basis'] = basis
    data.append(mydf)
  # end for basis

  wf_fname_map = {
    'Slater':'s_vmc.dat',
    'Slater-Jastrow':'sj_qmc.dat'
  }
  # collect Slater
  for wf in ['Slater','Slater-Jastrow']:
    fname = wf_fname_map[wf]
    qdf = parse_qmcas_output(fname)
    qdf['wf'] = wf
    data.append(qdf)
  # end for wf

  df_all = pd.concat(data).reset_index(drop=True)
  return df_all
# end def collect_all

if __name__ == '__main__':

  df = collect_all()
  df.to_json('energy.json')

