#!/usr/bin/env python
import os
import subprocess as sp

if __name__ == '__main__':
  subfix = 'scf'
  cmd = 'find . -name *-%s.sbatch.in' % subfix
  proc = sp.Popen(cmd,shell=True,stdout=sp.PIPE,stderr=sp.PIPE)
  out,err = proc.communicate()

  for fname in out.split('\n')[:-1]:
    sub_dir  = os.path.dirname(fname)
    sub_file = os.path.basename(fname)

    sub_cmd = 'cd %s;sbatch %s'%(sub_dir,sub_file)
    sp.check_call(sub_cmd,shell=True)
  # end for
# end __main__
