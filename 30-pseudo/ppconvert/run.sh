#!/bin/bash

#ref="1s(2)2s(2)2p(6)3d(5)"
ref="1s(2)2s(2)2p(6)"
./ppconvert --gamess_pot Mn.BFD.gamess --log_grid --s_ref $ref --p_ref $ref --d_ref $ref --upf Mn.BFD.upf --xml Mn.BFD.xml

ref="1s(2)"
./ppconvert --gamess_pot O.BFD.gamess --log_grid --s_ref $ref --p_ref $ref --d_ref $ref --upf O.BFD.upf --xml O.BFD.xml
