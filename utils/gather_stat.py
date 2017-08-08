#!/usr/bin/env python
import h5py
import pandas as pd

# ==== Level 1 ==== 
class StatFile:
  def __init__(self):
    pass
  def read(self,fname):
    self.fp = h5py.File(fname)
  # end def
  def fill_value_lists(self):
    self.value_list = []
    self.value_sq_list = []
    self.fp.visit( self.find_value )
  # end def
  def find_value(self,name):
    if 'value' in name:
      if name.endswith('squared'):
        self.value_sq_list.append(name)
      elif name.endswith('/value'):
        self.value_list.append(name)
      else:
        raise RuntimeError('unrecogonized name %s'%name)
      # end if
    # end if
  # end def find_value

# end class StatFile

def path_value_dataframe(fname,nequil):
  # find all ['value','value_squared'], store in stath5.value_list & value_sq_list
  stath5 = StatFile()
  stath5.read(fname)
  stath5.fill_value_lists()

  # construct value dataframe
  data = []
  for iname in range(len(stath5.value_list)):
    path = stath5.value_list[iname]
    val  = stath5.fp[path].value[nequil:].mean(axis=0)
    entry = {'h5path':path,'value':val}
    data.append(entry)
  # end for iname
  df = pd.DataFrame(data)

  # TODO: add value_squared to df
  #  name.replace('_squared','')

  return df
# end def path_value_dataframe
# =================

if __name__ == '__main__':
  """
gofr_e_0_0: [u'cutoff', u'delta', u'value', u'value_squared']
sk: [u'kpoints', u'value', u'value_squared']
sksp: [u'd', u'u'] -> [u'value', u'value_squared']
mloc: [u'value', u'value_squared']
spin_density: [u'd', u'u'] -> [u'value', u'value_squared']
  """
  nequil = 1
  folder = '.'

  # step 1: reel'em in
  import subprocess as sp
  ftext = sp.check_output(['find',folder,'-path','*.stat.h5'])
  fname_list = ftext.split('\n')[:-1]
  df_list = []
  for fname in fname_list:
    df = path_value_dataframe(fname,nequil)
    df['fname'] = fname
    df_list.append(df)
  # end for

# end __main__
