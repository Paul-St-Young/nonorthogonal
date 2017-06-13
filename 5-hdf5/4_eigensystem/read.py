#!/usr/bin/env python
import os
import numpy as np
import pandas as pd

if __name__ == '__main__':

    df0 = pd.read_json('mydf.json')
    df  = df0.set_index(['ikpt','ispin','istate'],drop=True)
    print df.tail()
