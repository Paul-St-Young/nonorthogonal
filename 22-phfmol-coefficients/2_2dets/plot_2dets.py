#!/usr/bin/env python
import numpy as np
def get_xlabel(results): 
  xlabels = []
  for i in range(len(results)):
    xlabel = '%d%d' % (results[i,0],results[i,1])
    xlabels.append(xlabel)
  # end for i
  return xlabels
# end def get_xlabel
def plot_result(ax,results,label,fmt='x-'):
  myx = range(len(results))
  line = ax.errorbar(myx,results[:,2],yerr=results[:,3],fmt=fmt,label=label)
  return line
# end def
def plot_phfmol(ax,results,label):
  myx = range(len(results))
  line = ax.plot(myx,results[:,2],c='k',marker='o',ms=8,label=label)
  return line
# end def 

if __name__ == '__main__':
  # 9 columns: idet jdet E Ee corr  V Ve corr  ratio
  real_results = np.loadtxt('real_code/results.dat')
  comp_results = np.loadtxt('comp_code/results.dat')
  phfmol_results = np.loadtxt('gen_dets/paire.dat')
  npair = len(real_results)

  import matplotlib.pyplot as plt
  fig,ax = plt.subplots(1,1)
  ax.set_xlim(-0.5,npair-0.5)
  ax.set_xticks(range(npair))
  ax.set_xlabel('determinant pair',fontsize=16)
  ax.set_ylabel('total energy (ha)',fontsize=16)

  rline = plot_result(ax,real_results,label='RealType')
  cline = plot_result(ax,comp_results,label='ValueType')
  mline = plot_phfmol(ax,phfmol_results,label='phfmol')
  ax.legend(loc='lower right')

  # set labels
  xlabels = get_xlabel(real_results)
  ax.set_xticklabels(xlabels)

  clabels = get_xlabel(comp_results)
  plabels = get_xlabel(phfmol_results)
  for i in range(len(xlabels)):
    assert xlabels[i] == clabels[i]
    assert xlabels[i] == plabels[i]
  # end for i
  
  fig.tight_layout()
  #fig.savefig('det_pairs.eps')
  plt.show()
# end __main__
