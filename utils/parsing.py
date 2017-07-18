import numpy as np
import pandas as pd

def str2comp(line):
  tokens = line.strip(' ()\n').split(',')
  real,imag = map(float,tokens)
  return real+1j*imag
# end def str2comp

def read_fortran_array(fname):
  ci_list = []
  with open(fname,'r') as f:
    for line in f:
      ci_list.append( str2comp(line) )
    # end for
  # end with
  ci_arr = np.array(ci_list)
  return ci_arr
# end def

def parse_qmcas_output(fname):
  # parse outputs such as:
  #  ./detsci49_af/c2  series 0  -10.452004 +/- 0.006424    1.6   1.601253 +/- 0.068660    1.1   0.1532

  data = {'myid':[],'energy_mean':[],'energy_error':[],'method':[],'iqmc':[]}
  with open(fname) as f:
    for line in f:
      tokens = line.split()
      if len(tokens) == 12:
        # determine basis,ndet
        myid   = tokens[0]

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

        # organize data in table
        label_and_data = zip(
          ['myid','energy_mean','energy_error','method','iqmc'],
          [myid,energy_mean,energy_error,method,iqmc]
        )
        for key,val in label_and_data:
          data[key].append(val)
        # end for
      # end if
    # end for 
  # end with
  mydf = pd.DataFrame(data)
  return mydf
# end def parse_qmcas_output
