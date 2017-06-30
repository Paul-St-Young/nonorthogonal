#!/usr/bin/env python
import os
import numpy as np

def set_ci_coeff(fname,ci_coeff,comp=False):
  """ set the CI coefficients of a multideterminant input 'fname' """
  ndet = len(ci_coeff)

  from input_xml import InputXml
  inp = InputXml()
  inp.read(fname)

  detlist_node = inp.find('.//detlist')
  nsize = int( detlist_node.get('size') )
  assert nsize == ndet

  ci_nodes = detlist_node.findall('.//ci')
  for idet in range(ndet):
    ci_node = ci_nodes[idet]
    ci_id   = ci_node.get('id')

    try:
      idx = int( ci_id.split('_')[-1] )
      assert idx == idet
    except:
      raise RuntimeError('failed to check @id of ci, uncomment this line to skip the check')
    # end try

    ci = ci_coeff[idet]
    ci_text = str(ci)
    if comp:
      ci_text = '(%6.4f,%6.4f)' % (ci.real,ci.imag)
    ci_node.set( 'coeff',ci_text )
  # end for idet

  inp.write(fname)
  del inp
# end def set_ci_coeff

def change_submission_name(sub_file,name):
  import os
  os.system( 'sed -i "s/my_epic_name/%s/" %s' % (name,sub_file) )
# end def


if __name__ == '__main__':

  template_dir = 'noqc_det34'
  input_name   = 'vmc.xml'
  sub_fname    = 'qmc.msub'
  ndet = 5

  from itertools import combinations
  for (idet,jdet) in combinations(range(5),2):
    run_dir = 'det%d%d'%(idet,jdet)
    ci_coeff = np.zeros(ndet,dtype=complex)
    ci_coeff[idet] = 1.0
    ci_coeff[jdet] = 1.0
    if not os.path.isdir(run_dir):
      os.system('cp -r %s %s' % (template_dir,run_dir) )
      fname = os.path.join(run_dir,input_name)
      set_ci_coeff(fname,ci_coeff,comp=True)
      change_submission_name( os.path.join(run_dir,sub_fname), run_dir )
      # !!!! hack submission in here
      os.system('cd %s; sbatch %s' % (run_dir,sub_fname))
    # end if
  # end for 

# end __main__
