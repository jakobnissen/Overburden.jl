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

OUT="$(mktemp -d)"
trap 'rm -rf "$OUT"' EXIT

err="$("$EXE" -d 1 test/data/sorted.bam "$OUT/out.bam" 2>&1)" || fail "exit code on sorted input"
[[ "$err" == "Wrote 4 records of 3 templates" ]] || fail "unexpected output on sorted input: $err"
# Check the output is a valid BAM file, by reading it back
"$EXE" -d 1 "$OUT/out.bam" "$OUT/out2.bam" 2>/dev/null || fail "exit code on reading the output"

status=0
"$EXE" -d 1 test/data/sorted.bam "$OUT/out.bam" 2>/dev/null || status=$?
[[ $status -eq 1 ]] || fail "exit code $status on existing output"

status=0
err="$("$EXE" test/data/unsorted.bam "$OUT/unsorted.bam" 2>&1)" || status=$?
[[ $status -eq 1 ]] || fail "exit code $status on unsorted input"
[[ "$err" == *"not coordinate sorted"* ]] || fail "unexpected error on unsorted input: $err"
[[ ! -e "$OUT/unsorted.bam" ]] || fail "output left behind after unsorted input"

status=0
"$EXE" test/data/missing.bam "$OUT/missing.bam" 2>/dev/null || status=$?
[[ $status -eq 1 ]] || fail "exit code $status on missing input"

for args in "" "test/data/sorted.bam" "-d 0 test/data/sorted.bam $OUT/out.bam"; do
    status=0
    # shellcheck disable=SC2086
    "$EXE" $args 2>/dev/null || status=$?
    [[ $status -eq 2 ]] || fail "exit code $status with arguments '$args'"
done

out="$("$EXE" --help)" || fail "exit code with --help"
[[ "$out" == *"Usage: overburden"* ]] || fail "unexpected output with --help: $out"
out="$("$EXE" --version)" || fail "exit code with --version"
[[ "$out" == "overburden "* ]] || fail "unexpected output with --version: $out"

echo "Smoke test OK"
