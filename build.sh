#!/usr/bin/env bash
# Compile Overburden with `juliac --trim=safe` into build/bin/overburden, and
# smoke test the executable on the BAM files in test/data.
# Requires Julia >= 1.13. Usage: ./build.sh
# Set JULIAC_ENV to reuse an environment with JuliaC installed.
# The BAM files were created from the SAM files with
# `samtools view --no-PG -b -o X.bam X.sam`.
set -euo pipefail
cd "$(dirname "$0")"

JULIAC_ENV="${JULIAC_ENV:-$(mktemp -d)}"
julia --project="$JULIAC_ENV" -e 'using Pkg; Pkg.add("JuliaC")'
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# --bundle fails if the target directory already exists.
rm -rf build
# Note: If verification fails, juliac prints "Verifier error" lines, and exits non-zero.
julia --project="$JULIAC_ENV" -e 'using JuliaC; JuliaC._main_cli(ARGS)' -- \
    --output-exe overburden --trim=safe --bundle build .

EXE=./build/bin/overburden
fail() {
    echo "SMOKE TEST FAILED: $1" >&2
    exit 1
}

out="$("$EXE" test/data/sorted.bam)" || fail "exit code on sorted input"
[[ "$out" == "selected templates: 3" ]] || fail "unexpected output on sorted input: $out"

status=0
err="$("$EXE" test/data/unsorted.bam 2>&1)" || status=$?
[[ $status -eq 1 ]] || fail "exit code $status on unsorted input"
[[ "$err" == *"not coordinate sorted"* ]] || fail "unexpected error on unsorted input: $err"

status=0
"$EXE" test/data/missing.bam 2>/dev/null || status=$?
[[ $status -eq 1 ]] || fail "exit code $status on missing input"

status=0
"$EXE" 2>/dev/null || status=$?
[[ $status -eq 2 ]] || fail "exit code $status with no arguments"

echo "Smoke test OK"
