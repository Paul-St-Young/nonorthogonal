import sys
import pyscf
import numpy
from pyscf.pbc import gto,scf
from pyscf.pbc.scf import KRHF as RHF
from pyscf.pbc.scf import KUHF as UHF
from pyscf.pbc.dft import KRKS as RKS
from pyscf.pbc.dft import KUKS as UKS
from pyscf2qwalk import print_qwalk
basis={
'Mn':pyscf.gto.basis.parse('''
Mn s
  23.6450680 -0.014659 
  13.4620290 0.223661 
  8.21262900 -0.564535 
  1.85994300 0.554600 
  0.88975400 0.534536 
  0.42256000 0.142766 
 Mn p
  17.8033170 0.003171 
  8.83499500 -0.080060 
  4.16134000 0.185735 
  2.12919600 0.375367 
  1.06371300 0.388624 
  0.51014100 0.192450 
  0.21130200 0.026146 
 Mn d
  9.25373600 0.076285 
  4.36301600 0.247722 
  1.83482100 0.330050 
  0.75611800 0.334518 
  0.29609200 0.263545 
 Mn f
  1.0998003483 1.00000 
 Mn s
0.2 1.0
 Mn s
0.6000000000000001 1.0
 Mn p
0.2 1.0
 Mn p
0.6000000000000001 1.0
 Mn d
0.2 1.0
 Mn d
0.6000000000000001 1.0
'''),
'O':pyscf.gto.basis.parse('''
O s
  0.268022 0.304848 
  0.573098 0.453752 
  1.225429 0.295926 
  2.620277 0.019567 
  5.602818 -0.128627 
  11.980245 0.012024 
  25.616801 0.000407 
  54.775216 -0.000076 
 O p
  0.333673 0.255999 
  0.666627 0.281879 
  1.331816 0.242835 
  2.660761 0.161134 
  5.315785 0.082308 
  10.620108 0.039899 
  21.217318 0.004679 
 O d
  0.669340 1.000000 
 O f
  1.423104 1.000000 
 O s
0.2 1.0
 O s
0.6000000000000001 1.0
 O p
0.2 1.0
 O p
0.6000000000000001 1.0
'''),
}
mol=gto.M(verbose=4,
gs=[4, 4, 4],
atom='''Mn 0.0 0.0 0.0
O -3.270036660486693 -2.0024802438743468 1.377207824448945e-07
''',
a=''' -2.9975315493471464 0.0 0.9037894892875293 -2.180022986975804 -1.3349859392589383 -1.8075805107124707 -1.3625144246044611 -2.6699718785178765 0.9037894892875293 ''',
basis=basis,
spin=0,
ecp='bfd')
mol.charge=0
kpts=mol.make_kpts([1, 1, 1])
m=RHF(mol,kpts)
m.max_cycle=50
m.direct_scf_tol=1e-07
m.chkfile='pyscf_driver.py.chkfile'
m.conv_tol=1e-07
m.diis_start_cycle=1
init_dm=pyscf.scf.rhf.init_guess_by_minao(mol)
dm_kpts= [init_dm for k in range(len(kpts))]
m.xc="pbe,pbe"
print('E(HF) =',m.kernel(dm_kpts))
print_qwalk(mol,m)
print ("All_done")