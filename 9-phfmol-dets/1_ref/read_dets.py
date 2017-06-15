#!/usr/bin/env python
from one_det import read_next_det
from mmap import mmap

if __name__ == '__main__':
    det_fname = 'next/dets.2'
    with open(det_fname,'r+') as f:
        mm = mmap(f.fileno(),0)
    # end with

    dets = []
    for idet in range(3):
        idx,det = read_next_det(mm)
        print idx
        dets.append(det)

    print [det.shape for det in dets]
