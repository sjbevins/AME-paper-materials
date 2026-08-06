#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
PYTHON=${PYTHON:-python3}

run_and_compare() {
  local script="$1"
  local expected="$2"
  local actual="${script%.py}.actual.txt"
  echo "==> $script"
  "$PYTHON" "$script" > "$actual"
  diff -u --strip-trailing-cr "$expected" "$actual"
  rm -f "$actual"
}

run_and_compare check_consistency.py check_consistency.expected.txt
run_and_compare verify_ame_12_5.py verify_ame_12_5.expected.txt
run_and_compare verify_ame_18_11.py verify_ame_18_11.expected.txt
run_and_compare verify_ame_18_13.py verify_ame_18_13.expected.txt

if "$PYTHON" -c 'import numpy, numba' >/dev/null 2>&1; then
  run_and_compare verify_17_11_numba.py verify_17_11_numba.expected.txt
else
  echo "SKIP q=11 all-deletion census (NumPy/Numba unavailable)"
fi

if command -v g++ >/dev/null 2>&1 && echo '#include <boost/multiprecision/cpp_int.hpp>' | g++ -x c++ -std=c++17 -fsyntax-only - >/dev/null 2>&1; then
  echo "==> corrected template census"
  g++ -O3 -std=c++17 template_screen_generic.cpp -o .template_screen_generic
  ./.template_screen_generic 11 > template_q11.actual.txt 2> template_q11.actual.stderr
  ./.template_screen_generic 13 > template_q13.actual.txt 2> template_q13.actual.stderr
  diff -u --strip-trailing-cr template_q11_corrected.expected.txt template_q11.actual.txt
  diff -u --strip-trailing-cr template_q13_corrected.expected.txt template_q13.actual.txt
  rm -f .template_screen_generic template_q11.actual.txt template_q13.actual.txt template_q11.actual.stderr template_q13.actual.stderr
else
  echo "SKIP template census (g++ or Boost headers unavailable)"
fi

echo 'ALL PACKAGED VERIFICATION CHECKS PASSED'
