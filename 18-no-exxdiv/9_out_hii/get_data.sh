#!/bin/bash
qmca -q ev --sac det*/*scalar.dat > results.dat
#tail -n +2 results.dat | awk '{print $4}' > data.dat
#awk '{print $1+10.257441}' data.dat > diff.dat
