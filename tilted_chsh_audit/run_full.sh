#!/bin/bash
# Full independent B&B run for V1 then V3 (checkpointed: re-running skips finished strips).
cd "$(dirname "$0")"
for V in V1 V3; do
  echo "start $V $(date '+%F %T %Z')"
  s=$(date +%s)
  python3 -c "
import resource, runpy, sys
sys.argv = ['bnb.py', '$V', '--nstrips', '4096', '--minw', '0.0001', '--out', 'bnb_$V.jsonl']
try:
    runpy.run_path('bnb.py', run_name='__main__')
finally:
    print('maxrss_MB', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, file=sys.stderr)
" > bnb_$V.log 2> bnb_$V.err
  echo "end $V $(date '+%F %T %Z') wall $(( $(date +%s) - s )) s"
done
