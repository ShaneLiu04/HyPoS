#!/bin/bash
# Collect medians from ar004 baseline runs. Usage: bash collect_ar004.sh <rootdir>
# Layout: <rootdir>/<case>_r{1,2,3}/performance_report.json
set -u
ROOT="${1:-/root/ar004-baseline}"
echo "case | iter_ms_med | iters_med | final_res_med | raw_iter_ms"
for d in "$ROOT"/*/; do
  name=$(basename "$d")
  case_=${name%_r[0-9]}
  r=${name##*_r}
  [ "$case_" = "$name" ] && continue
  [ "$r" -ge 1 ] 2>/dev/null && [ "$r" -le 3 ] 2>/dev/null || continue
  f="$d/performance_report.json"
  [ -f "$f" ] || continue
  v=$(grep -o '"iter_time_ms": *[0-9.eE-]*' "$f" | grep -o '[0-9.eE-]*$')
  echo "MS $case_ $v"
  v=$(grep -o '"iterations": *[0-9]*' "$f" | grep -o '[0-9]*$')
  echo "IT $case_ $v"
  v=$(grep -o '"final_residual": *[0-9.eE-]*' "$f" | grep -o '[0-9.eE-]*$')
  echo "FR $case_ $v"
done | awk '
  /^MS / { ms[$2] = ms[$2] " " $3; cases[$2]=1 }
  /^IT / { it[$2] = it[$2] " " $3 }
  /^FR / { fr[$2] = fr[$2] " " $3 }
  END {
    # portable (mawk-safe): insertion sort helper for <=3 values
    for (c in ms) {
      split(ms[c], a, " ")
      na = 0
      for (k in a) { na++; tmp[na] = a[k] + 0 }
      for (x = 2; x <= na; x++) { v = tmp[x]; y = x - 1; while (y >= 1 && tmp[y] > v) { tmp[y+1] = tmp[y]; y-- } tmp[y+1] = v }
      m = tmp[int((na+1)/2)]
      split(it[c], b, " ")
      nb = 0
      for (k in b) { nb++; tmi[nb] = b[k] + 0 }
      for (x = 2; x <= nb; x++) { v = tmi[x]; y = x - 1; while (y >= 1 && tmi[y] > v) { tmi[y+1] = tmi[y]; y-- } tmi[y+1] = v }
      mi = tmi[int((nb+1)/2)]
      split(fr[c], e, " ")
      nf = 0
      for (k in e) { nf++; tmf[nf] = e[k] + 0 }
      for (x = 2; x <= nf; x++) { v = tmf[x]; y = x - 1; while (y >= 1 && tmf[y] > v) { tmf[y+1] = tmf[y]; y-- } tmf[y+1] = v }
      mf = tmf[int((nf+1)/2)]
      printf "%s | %.4f | %d | %.6e |%s\n", c, m, mi, mf, ms[c]
    }
  }'
