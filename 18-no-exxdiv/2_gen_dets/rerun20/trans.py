#!/usr/bin/env python

if __name__ == '__main__':
  import numpy as np

  ci_list = []
  with open('ham20.out','r') as f:
    for line in f:
      tokens = line.strip(' ()\n').split(',')
      real_imag = map(float,tokens)
      for num in real_imag:
        ci_list.append(num)
      # end for
    # end for

  # end with

  np.savetxt('ci_coeff.dat',ci_list)
# end __main__
