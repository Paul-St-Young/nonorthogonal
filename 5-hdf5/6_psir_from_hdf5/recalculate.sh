#!/bin/bash

# rebuild eigensystem dataframe
cd ../4_eigensystem
python eig_from_KRKS.py

# rebuild pwscf.pwscf.h5
cd ../5_pyscf_hdf5
python pyscf2qmcpack.py

cd ../6_psir_from_hdf5
python psir.py
