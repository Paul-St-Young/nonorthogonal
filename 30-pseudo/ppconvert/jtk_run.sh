#!/bin/bash

pot=Mn.BFD
ref="1s(2)2p(6)3d(3)2s(2)"
loc=1 # 0 for s, 1 for p, 2 for d

ppconvert --local_channel $loc  --gamess_pot ${pot}.gamess --s_ref $ref --p_ref $ref --d_ref $ref --f_ref $ref --log_grid --upf ${pot}.upf>&convert_upf.out

ppconvert --local_channel $loc  --gamess_pot ${pot}.gamess --s_ref $ref --p_ref $ref --d_ref $ref --f_ref $ref --xml ${pot}.xml>&convert_xml.out
