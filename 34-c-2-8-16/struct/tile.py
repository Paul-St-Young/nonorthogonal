#!/usr/bin/env python

if __name__ == '__main__':
  from nexus import Structure
  struct = Structure()
  struct.read('c2.xsf')

  for volfac in [4.,8.]:
    Topt,ropt = struct.opt_tilematrix(volfac=volfac)
    new_struct = struct.tile(Topt)
    natom = new_struct.pos.shape[0]
    new_struct.write('c%d.xsf'%natom)
  # end for volfac
# end __main__
