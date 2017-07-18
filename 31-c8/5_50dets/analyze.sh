#!/bin/bash

qmca -q ev --sac mdets/detsci*/*scalar.dat > s_vmc.dat
qmca -q ev --sac */dmc/detsci*/*scalar.dat > sj_qmc.dat
