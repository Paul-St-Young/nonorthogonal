#!/usr/bin/env python
import os
import numpy as np
import pandas as pd

if __name__ == '__main__':

    df0 = pd.read_json('mydf.json')
    df  = df0.set_index(['ikpt','ispin','istate'],drop=True)

    ikpt   = 0
    ispin  = 0
    istate = 0

    psig_arr = np.array( df.loc[(ikpt,ispin,istate),'evector'] )
    psig = psig_arr[:,0]+1j*psig_arr[:,1]
    print psig.shape
