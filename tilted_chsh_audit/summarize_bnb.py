import json, sys
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 50
S_LO = Fr(int(mp.floor(mp.sin(mp.pi/4) * 10**30)), 10**30)
S_HI = Fr(int(mp.ceil(mp.sin(3*mp.pi/8) * 10**30)), 10**30)
for V in sys.argv[1:]:
    recs = {}
    for line in open('bnb_%s.jsonl' % V):
        d = json.loads(line); recs[d['strip']] = d
    n = 4096
    idx = sorted(recs)
    allok = all(recs[i]['ok'] for i in idx)
    full = idx == list(range(n))
    # exact tiling of [S_LO, S_HI]
    tile = full and Fr(recs[0]['s_lo']) == S_LO and Fr(recs[n-1]['s_hi']) == S_HI and \
        all(Fr(recs[i]['s_hi']) == Fr(recs[i+1]['s_lo']) for i in range(n-1)) if full else None
    tests = {}
    for d in recs.values():
        for k, v in d['tests'].items(): tests[k] = tests.get(k, 0) + v
    print(V, 'strips done %d/%d' % (len(idx), n), 'all ok:', allok, 'exact tiling of S_Theta:', tile,
          'sub-strips:', sum(d['substrips'] for d in recs.values()), 'c-leaf boxes:', sum(d['cboxes'] for d in recs.values()),
          'tests:', tests, 'cpu s (sum of strips): %.0f' % sum(d['seconds'] for d in recs.values()),
          'theta/pi range: [%.6f, %.6f]' % (recs[idx[0]]['theta_over_pi'][0], recs[idx[-1]]['theta_over_pi'][1]))
    bad = [i for i in idx if not recs[i]['ok']]
    if bad: print('  FAILED strips:', bad[:20])
