#!/bin/bash
# Benchmark script: Compare Python SWORD2 vs Rust SWORD2
# Tests correctness (output comparison) and performance (runtime)
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUST_BIN="$SCRIPT_DIR/sword2-rs/target/release/sword2"
PYTHON_BIN="python3 $SCRIPT_DIR/SWORD2.py"
TEST_DIR="/home/user/test_results"
SIZES=(50 150 300)

echo "=============================================="
echo "  SWORD2: Python vs Rust Benchmark"
echo "=============================================="
echo ""

# Generate test PDB files of different sizes
echo "[1/4] Generating test PDB files..."
python3 -c "
import numpy as np
import sys

def make_test_pdb(n_residues, filename):
    lines = []
    atom_idx = 1
    atoms = [('N', (-0.5, 0.0, 0.0)), ('CA', (0.0, 0.0, 0.0)), ('C', (1.0, 0.0, 0.0)), ('O', (1.5, 1.0, 0.0))]
    aa = 'ALA'
    for i in range(n_residues):
        angle = i * 100 * np.pi / 180
        z = i * 1.5
        x_off = 2.3 * np.cos(angle)
        y_off = 2.3 * np.sin(angle)
        for aname, (dx, dy, dz) in atoms:
            x = x_off + dx
            y = y_off + dy
            z_coord = z + dz
            lines.append(f'ATOM  {atom_idx:5d}  {aname:<3s} {aa} A{i+1:4d}    {x:8.3f}{y:8.3f}{z_coord:8.3f}  1.00  0.00           {aname[0]:>2s}')
            atom_idx += 1
    lines.append('END')
    with open(filename, 'w') as f:
        f.write('\n'.join(lines) + '\n')

for n in [50, 150, 300]:
    fn = '$TEST_DIR/test_{}res.pdb'.format(n)
    make_test_pdb(n, fn)
    print(f'  Created {fn} ({n} residues)')
"
echo ""

# Run benchmarks
echo "[2/4] Running benchmarks (no energies, no plots)..."
echo ""
printf "%-10s | %-15s | %-15s | %-10s\n" "Size" "Python (s)" "Rust (s)" "Speedup"
printf "%-10s-+-%-15s-+-%-15s-+-%-10s\n" "----------" "---------------" "---------------" "----------"

for SIZE in "${SIZES[@]}"; do
    PDB_FILE="$TEST_DIR/test_${SIZE}res.pdb"
    PY_OUT="$TEST_DIR/bench_python_${SIZE}"
    RS_OUT="$TEST_DIR/bench_rust_${SIZE}"

    # Clean previous results
    rm -rf "$PY_OUT" "$RS_OUT"

    # Run Python version
    PY_START=$(date +%s%N)
    $PYTHON_BIN -i "$PDB_FILE" -c A -o "$PY_OUT" -e -l 2>/dev/null || true
    PY_END=$(date +%s%N)
    PY_TIME=$(echo "scale=2; ($PY_END - $PY_START) / 1000000000" | bc)

    # Run Rust version
    RS_START=$(date +%s%N)
    $RUST_BIN -i "$PDB_FILE" -c A -o "$RS_OUT" -e -l --base-dir "$SCRIPT_DIR" 2>/dev/null || true
    RS_END=$(date +%s%N)
    RS_TIME=$(echo "scale=2; ($RS_END - $RS_START) / 1000000000" | bc)

    # Calculate speedup
    if [ "$(echo "$RS_TIME > 0" | bc)" -eq 1 ]; then
        SPEEDUP=$(echo "scale=2; $PY_TIME / $RS_TIME" | bc)
    else
        SPEEDUP="N/A"
    fi

    printf "%-10s | %-15s | %-15s | %-10s\n" "${SIZE} res" "${PY_TIME}s" "${RS_TIME}s" "${SPEEDUP}x"
done
echo ""

# Verify output correctness
echo "[3/4] Verifying output correctness..."
echo ""
ALL_PASS=true

for SIZE in "${SIZES[@]}"; do
    PY_DIR=$(ls -d "$TEST_DIR/bench_python_${SIZE}/test_${SIZE}res_A"* 2>/dev/null | head -1)
    RS_DIR="$TEST_DIR/bench_rust_${SIZE}/test_${SIZE}res_A"

    if [ -z "$PY_DIR" ] || [ ! -d "$RS_DIR" ]; then
        echo "  [${SIZE} res] SKIP - output directories not found"
        continue
    fi

    # Compare SWORD2_summary.txt
    if diff -q "$PY_DIR/SWORD2_summary.txt" "$RS_DIR/SWORD2_summary.txt" >/dev/null 2>&1; then
        echo "  [${SIZE} res] SWORD2_summary.txt: PASS"
    else
        echo "  [${SIZE} res] SWORD2_summary.txt: FAIL"
        diff "$PY_DIR/SWORD2_summary.txt" "$RS_DIR/SWORD2_summary.txt" || true
        ALL_PASS=false
    fi

    # Compare JSON (content, not formatting)
    JSON_MATCH=$(python3 -c "
import json, sys
try:
    with open('$PY_DIR/SWORD2_summary.json') as f: py = json.load(f)
    with open('$RS_DIR/SWORD2_summary.json') as f: rs = json.load(f)
    print('PASS' if py == rs else 'FAIL')
except Exception as e:
    print(f'ERROR: {e}')
" 2>&1)
    echo "  [${SIZE} res] SWORD2_summary.json: $JSON_MATCH"
    if [ "$JSON_MATCH" != "PASS" ]; then ALL_PASS=false; fi

    # Compare PEELING_summary.txt
    if diff -q "$PY_DIR/PEELING_summary.txt" "$RS_DIR/PEELING_summary.txt" >/dev/null 2>&1; then
        echo "  [${SIZE} res] PEELING_summary.txt: PASS"
    else
        echo "  [${SIZE} res] PEELING_summary.txt: FAIL"
        diff "$PY_DIR/PEELING_summary.txt" "$RS_DIR/PEELING_summary.txt" || true
        ALL_PASS=false
    fi

    # Compare mapping file
    if diff -q "$PY_DIR/mapping_auth_resnums.txt" "$RS_DIR/mapping_auth_resnums.txt" >/dev/null 2>&1; then
        echo "  [${SIZE} res] mapping_auth_resnums.txt: PASS"
    else
        echo "  [${SIZE} res] mapping_auth_resnums.txt: FAIL"
        ALL_PASS=false
    fi

    echo ""
done

# Run with energies on 150-res for comparison
echo "[4/4] Running benchmark WITH energies (150 residues)..."
echo ""

PDB_FILE="$TEST_DIR/test_150res.pdb"
PY_OUT="$TEST_DIR/bench_energy_python_150"
RS_OUT="$TEST_DIR/bench_energy_rust_150"
rm -rf "$PY_OUT" "$RS_OUT"

# Python with energies
PY_START=$(date +%s%N)
$PYTHON_BIN -i "$PDB_FILE" -c A -o "$PY_OUT" -l 2>/dev/null || true
PY_END=$(date +%s%N)
PY_TIME=$(echo "scale=2; ($PY_END - $PY_START) / 1000000000" | bc)

# Rust with energies
RS_START=$(date +%s%N)
$RUST_BIN -i "$PDB_FILE" -c A -o "$RS_OUT" -l --base-dir "$SCRIPT_DIR" 2>/dev/null || true
RS_END=$(date +%s%N)
RS_TIME=$(echo "scale=2; ($RS_END - $RS_START) / 1000000000" | bc)

if [ "$(echo "$RS_TIME > 0" | bc)" -eq 1 ]; then
    SPEEDUP=$(echo "scale=2; $PY_TIME / $RS_TIME" | bc)
else
    SPEEDUP="N/A"
fi

printf "%-20s | %-15s | %-15s | %-10s\n" "Task" "Python (s)" "Rust (s)" "Speedup"
printf "%-20s-+-%-15s-+-%-15s-+-%-10s\n" "--------------------" "---------------" "---------------" "----------"
printf "%-20s | %-15s | %-15s | %-10s\n" "150 res + energies" "${PY_TIME}s" "${RS_TIME}s" "${SPEEDUP}x"
echo ""

# Verify energy outputs match
PY_DIR_E=$(ls -d "$TEST_DIR/bench_energy_python_150/test_150res_A"* 2>/dev/null | head -1)
RS_DIR_E="$TEST_DIR/bench_energy_rust_150/test_150res_A"

if [ -n "$PY_DIR_E" ] && [ -d "$RS_DIR_E" ]; then
    if diff -q "$PY_DIR_E/SWORD2_summary.txt" "$RS_DIR_E/SWORD2_summary.txt" >/dev/null 2>&1; then
        echo "  Energy SWORD2_summary.txt: PASS"
    else
        echo "  Energy SWORD2_summary.txt: FAIL"
        diff "$PY_DIR_E/SWORD2_summary.txt" "$RS_DIR_E/SWORD2_summary.txt" || true
    fi

    if diff -q "$PY_DIR_E/PEELING_summary.txt" "$RS_DIR_E/PEELING_summary.txt" >/dev/null 2>&1; then
        echo "  Energy PEELING_summary.txt: PASS"
    else
        echo "  Energy PEELING_summary.txt: FAIL"
        diff "$PY_DIR_E/PEELING_summary.txt" "$RS_DIR_E/PEELING_summary.txt" || true
    fi
fi

echo ""
echo "=============================================="
if $ALL_PASS; then
    echo "  ALL CORRECTNESS CHECKS PASSED"
else
    echo "  SOME CHECKS FAILED - see above"
fi
echo "=============================================="
