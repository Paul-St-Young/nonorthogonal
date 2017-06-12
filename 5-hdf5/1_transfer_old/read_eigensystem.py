#!/usr/bin/env python
import os
import h5py
import numpy as np
import pandas as pd
from pwscf_h5 import PwscfH5

if __name__ == '__main__':

    ref_fname = './ref/pwscf.pwscf.h5'
    ref = PwscfH5()
    ref.read(ref_fname)
    gvec = ref.get('gvectors')
    df = ref.eigensystem()

    #new  = h5py.File('test.h5','w')
    #ref.create_electrons_group(new,gvec,df)
