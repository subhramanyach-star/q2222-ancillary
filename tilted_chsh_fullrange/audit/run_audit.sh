#!/bin/bash
# Audit spot re-runs (fresh processes, archived certificates). Times in HKT.
cd "$(dirname "$0")/../fullrange"
A=../audit   # same relative path from fullrange/ and smalltheta/
md5sum vgen.py c_Wa.json c_V1h.json c_V3h.json c_V3l.json c_Wa_negpin_a.json c_V3l_negmid_0.json c_Wa_neg.json > $A/md5_before.txt
(cd ../smalltheta; md5sum sgen.py tgen3.py swap_data.json c_V1s.json c_V1s_ext.json c_V3t.json c_V3m.json hats_V3_V1.json c_V1s_nega.json c_V3t_negep.json >> $A/md5_before.txt)
nice python3 vgen.py c_Wa.json --strips 12 13 > $A/a_Wa12.log 2>&1
nice python3 vgen.py c_Wa.json --strips 23 24 > $A/a_Wa23.log 2>&1
nice python3 vgen.py c_V1h.json --strips 7 8 > $A/a_V1h7.log 2>&1
nice python3 vgen.py c_V3h.json --strips 23 24 > $A/a_V3h23.log 2>&1
nice python3 vgen.py c_V3l.json --strips 0 1 > $A/a_V3l0.log 2>&1
nice python3 vgen.py c_Wa_negpin_a.json > $A/n_Wa_negpin_a.log 2>&1
nice python3 vgen.py c_V3l_negmid_0.json > $A/n_V3l_negmid_0.log 2>&1
cd ../smalltheta
nice python3 sgen.py c_V1s.json > $A/a_V1s.log 2>&1
nice python3 sgen.py c_V1s_ext.json --tstrips 3 4 > $A/a_V1s_ext3.log 2>&1
nice python3 sgen.py c_V1s_nega.json > $A/n_V1s_nega.log 2>&1
nice python3 tgen3.py c_V3t.json --tstrips 0 1 > $A/a_V3t0.log 2>&1
nice python3 tgen3.py c_V3t_negep.json --tstrips 0 1 > $A/n_V3t_negep.log 2>&1
echo DONE $(date) > $A/done
