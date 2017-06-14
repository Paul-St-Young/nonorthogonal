#! /usr/bin/env python

from generic import obj
from nexus import settings,Job,run_project
from nexus import generate_physical_system
from nexus import generate_gamess,FormattedGroup
from nexus import generate_convert4qmc
from nexus import generate_qmcpack,vmc
from nexus import generate_cusp_correction
from debug import *

data = FormattedGroup('''BH..RATZ..MCSCF
C1
B 5.0 0.0000 0.0000 0.0000
S   14
  1  33360.2170000              0.0000583        
  2   4972.0952000              0.0004565        
  3   1125.6417000              0.0024089        
  4    316.4913600              0.0101631        
  5    102.0072600              0.0360955        
  6     36.2958730              0.1059810        
  7     13.9714100              0.2395654        
  8      5.7411560              0.3737726        
  9      2.4942680              0.3081465        
 10      1.1142020              0.0776928        
 11      0.4215490              0.0014797        
 12      0.1696330              0.0004991        
 13      0.0685350             -0.0000262        
 14      0.0239870              0.0000236        
S   14
  1  33360.2170000             -0.0000121        
  2   4972.0952000             -0.0000950        
  3   1125.6417000             -0.0005028        
  4    316.4913600             -0.0021264        
  5    102.0072600             -0.0076744        
  6     36.2958730             -0.0232179        
  7     13.9714100             -0.0569468        
  8      5.7411560             -0.1062577        
  9      2.4942680             -0.1457694        
 10      1.1142020             -0.0430285        
 11      0.4215490              0.3719441        
 12      0.1696330              0.5652090        
 13      0.0685350              0.1921269        
 14      0.0239870              0.0219128        
S   14
  1  33360.2170000              0.0000121        
  2   4972.0952000              0.0000945        
  3   1125.6417000              0.0004997        
  4    316.4913600              0.0021192        
  5    102.0072600              0.0076519        
  6     36.2958730              0.0234045        
  7     13.9714100              0.0581156        
  8      5.7411560              0.1116739        
  9      2.4942680              0.1705623        
 10      1.1142020              0.0172676        
 11      0.4215490             -0.9989044        
 12      0.1696330             -0.2265722        
 13      0.0685350              1.0294550        
 14      0.0239870              0.2213786        
S   14
  1  33360.2170000             -0.0000180        
  2   4972.0952000             -0.0001292        
  3   1125.6417000             -0.0007771        
  4    316.4913600             -0.0027820        
  5    102.0072600             -0.0125600        
  6     36.2958730             -0.0301270        
  7     13.9714100             -0.1113291        
  8      5.7411560             -0.1047497        
  9      2.4942680             -0.4598841        
 10      1.1142020              0.6740179        
 11      0.4215490              1.5042519        
 12      0.1696330             -2.4760430        
 13      0.0685350              0.5776670        
 14      0.0239870              0.7429993        
S   14
  1  33360.2170000              0.0000215        
  2   4972.0952000              0.0001376        
  3   1125.6417000              0.0009649        
  4    316.4913600              0.0028031        
  5    102.0072600              0.0165450        
  6     36.2958730              0.0306156        
  7     13.9714100              0.1680322        
  8      5.7411560              0.0264152        
  9      2.4942680              0.8205322        
 10      1.1142020             -2.1545330        
 11      0.4215490              1.0217582        
 12      0.1696330              1.7061016        
 13      0.0685350             -3.0837640        
 14      0.0239870              1.8083378        
P   9
  1     55.0000000              0.0010186        
  2     13.3661010              0.0073958        
  3      4.1353940              0.0322514        
  4      1.4812560              0.0983464        
  5      0.6021360              0.2318024        
  6      0.2556360              0.3412436        
  7      0.1111510              0.2932983        
  8      0.0476480              0.1695951        
  9      0.0166770              0.1507569        
P   9
  1     55.0000000             -0.0007776        
  2     13.3661010             -0.0057768        
  3      4.1353940             -0.0226134        
  4      1.4812560             -0.0696834        
  5      0.6021360             -0.2551146        
  6      0.2556360             -0.3607786        
  7      0.1111510              0.0612395        
  8      0.0476480              0.2959938        
  9      0.0166770              0.6629278        
P   9
  1     55.0000000              0.0007590        
  2     13.3661010              0.0059151        
  3      4.1353940              0.0136009        
  4      1.4812560              0.0320628        
  5      0.6021360              0.6295261        
  6      0.2556360              0.1732649        
  7      0.1111510             -0.8917853        
  8      0.0476480             -0.3535550        
  9      0.0166770              0.8469041        
P   9
  1     55.0000000             -0.0032210        
  2     13.3661010             -0.0325794        
  3      4.1353940             -0.1313091        
  4      1.4812560             -0.6809573        
  5      0.6021360             -0.3491983        
  6      0.2556360              1.3937462        
  7      0.1111510             -0.3529796        
  8      0.0476480             -0.7747354        
  9      0.0166770              0.5768738        
D   4
  1      1.2000000              0.0803458        
  2      0.4200000              0.4219358        
  3      0.1470000              0.4639167        
  4      0.0514500              0.2855734        
D   4
  1      1.2000000             -0.1224628        
  2      0.4200000             -0.5604733        
  3      0.1470000             -0.0894260        
  4      0.0514500              0.9015938        
D   4
  1      1.2000000              0.5130462        
  2      0.4200000              0.5097938        
  3      0.1470000             -1.3189990        
  4      0.0514500              0.8663855        
F   3
  1      0.8500000              0.1825219        
  2      0.3400000              0.3946453        
  3      0.1360000              0.6104375        
F   3
  1      0.8500000             -0.4773313        
  2      0.3400000             -0.6578114        
  3      0.1360000              0.9513672        

H  1.0   0.00000  0.00000  2.362239
S   8
  1    188.6144500              0.0009639        
  2     28.2765960              0.0074920        
  3      6.4248300              0.0375954        
  4      1.8150410              0.1433950        
  5      0.5910630              0.3486363        
  6      0.2121490              0.4382974        
  7      0.0798910              0.1651066        
  8      0.0279620              0.0210229        
S   8
  1    188.6144500             -0.0013119        
  2     28.2765960             -0.0103451        
  3      6.4248300             -0.0504953        
  4      1.8150410             -0.2073855        
  5      0.5910630             -0.4350885        
  6      0.2121490             -0.0247297        
  7      0.0798910              0.3225260        
  8      0.0279620              0.7072754        
S   8
  1    188.6144500              0.0024224        
  2     28.2765960              0.0203382        
  3      6.4248300              0.0896394        
  4      1.8150410              0.4422907        
  5      0.5910630              0.5757144        
  6      0.2121490             -0.9802890        
  7      0.0798910             -0.6721538        
  8      0.0279620              1.1417685        
S   8
  1    188.6144500             -0.01157010       
  2     28.2765960             -0.08371540       
  3      6.4248300             -0.44516630       
  4      1.8150410             -1.14627100       
  5      0.5910630              2.50318710       
  6      0.2121490             -1.58284930       
  7      0.0798910              0.03096569       
  8      0.0279620              0.30862864       
P   4
  1      2.3050000              0.1127902        
  2      0.8067500              0.4185075        
  3      0.2823620              0.4700077        
  4      0.0988270              0.1826260        
P   4
  1      2.3050000             -0.2108688        
  2      0.8067500             -0.5943796        
  3      0.2823620              0.0896889        
  4      0.0988270              0.8611634        
P   4
  1      2.3050000              0.7599501        
  2      0.8067500              0.1646159        
  3      0.2823620             -1.3710140        
  4      0.0988270              1.0593155        
D   3
  1      1.8190000              0.2705134        
  2      0.7276000              0.5510125        
  3      0.2910400              0.3310866        
D   3
  1      1.8190000             -0.7938035        
  2      0.7276000             -0.0914252        
  3      0.2910400              0.8620033  

''' ) 

settings(
    runs          = 'mcscf', 
    status_only   = 0,
    generate_only = 0,
    sleep         = 3,
    machine       = 'ws16',
    ericfmt       = '/opt/intel17/apps/gamess_2016R1/auxdata/ericfmt.dat'
    )
gms_job = Job(app='rungms',serial=True)


sims = []

atom = generate_physical_system(
    type       = 'dimer',
    dimer      = ['B','H'],
    separation = 2.362239,
    units      = 'B',
    net_charge = 0,
    net_spin   = 0
    )

rhf_inputs = obj(
    identifier = 'rhf',
    system     = atom,
    job        = gms_job,
    scftyp     = 'rhf', 
    runtyp     = 'energy', 
    icharg     = 0, 
    ispher     = 1, 
    mult       = 1, 
    units      = 'bohr',
    mwords     = 200,
    guess      = 'huckel', 
    group      = 'c1', 
    gbasis     = 'cct'
    #data       = data
)
rhf = generate_gamess(
    path = 'mcscf',
    **rhf_inputs
)
sims.append(rhf)

mcscf_inputs = obj(
    identifier = 'mcscf',
    system     = atom,
    job        = gms_job,
    scftyp     = 'mcscf', 
    runtyp     = 'energy', 
    icharg     = 0, 
    ispher     = 1, 
    mult       = 1, 
    units      = 'bohr',
    mwords     = 200,
    guess      = 'moread', 
    norb       = 69,
    group      = 'c1', 
    ncore      = 0, 
    nels       = 6, 
    nact       = 9,  
    prttol     = 0.001, 
    mcscf      = obj(maxit=2000),
    #data       = data,
    fullnr     = True,
    memddi     = 100,
    gbasis     = 'cct'
    )
mcscf = generate_gamess(
    path = 'mcscf',
    dependencies = (rhf,'orbitals'),
    **mcscf_inputs
)
sims.append(mcscf)

mcscf_rerun = generate_gamess(
    path = 'mcscf_rerun',
    **mcscf_inputs
    )
mcscf_rerun.depends(mcscf,'orbitals')
sims.append(mcscf_rerun)

c4q = generate_convert4qmc(
    identifier = 'c4q',
    path = 'c4q',
    job  = Job(app='/home/yyang173/soft/qmcpack/build/bin/convert4qmc',serial=True),
    dependencies = (mcscf,'orbitals')
)
sims.append(c4q)

# ======= main ======
if __name__=='__main__':
    from nexus import ProjectManager
    pm = ProjectManager()
    pm.add_simulations(sims)
    pm.run_project()
# end __main__

'''
analyze = 0
if analyze:

    ma = mcscf.load_analyzer_image()
    mr = mcscf_rerun.load_analyzer_image()
    qa = vmc_scf.load_analyzer_image()
    le_scf = qa.qmc[0].scalars.LocalEnergy
    qa = vmc_cc.load_analyzer_image()
    le_cc = qa.qmc[0].scalars.LocalEnergy

    print 
    print 'Energies'
    print '  mcscf      :',ma.energy_components.total
    print '  mcscf_rerun:',mr.energy_components.total
    print '  qmcpack scf: {0} +/- {1}'.format(le_scf.mean,le_scf.error)
    print '  qmcpack cc : {0} +/- {1}'.format(le_cc.mean,le_cc.error)

#end if
'''
