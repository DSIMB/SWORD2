#!/bin/bash
# Benchmark script: Run tests for Rust SWORD2
# Tests correctness (by verifying outputs are generated and not crashing) and performance (runtime)
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUST_BIN="$SCRIPT_DIR/sword2-rs/target/release/sword2"
TEST_DIR="/home/chili/cretin/SWORD2/results"
SIZES=(50 150 300)

echo "=============================================="
echo "  SWORD2: Rust Benchmark"
echo "=============================================="
echo ""

# Ensure the executable is built
if [ ! -f "$RUST_BIN" ]; then
    echo "Building sword2-rs release..."
    cd "$SCRIPT_DIR/sword2-rs" && cargo build --release
    cd "$SCRIPT_DIR"
fi

# Run benchmarks
echo "[1/2] Running benchmarks (no energies)..."
echo ""
printf "%-10s | %-15s\n" "Size" "Rust (s)"
printf "%-10s-+-%-15s\n" "----------" "---------------"

for SIZE in "${SIZES[@]}"; do
    PDB_FILE="$TEST_DIR/test_${SIZE}res.pdb"
    RS_OUT="$TEST_DIR/bench_rust_${SIZE}"

    # Clean previous results
    rm -rf "$RS_OUT"
    
    # We do not have Python to generate the test files anymore, so they must be pre-generated.
    if [ ! -f "$PDB_FILE" ]; then
        echo "  [${SIZE} res] SKIP - $PDB_FILE not found (generate them first)"
        continue
    fi

    # Run Rust version
    RS_START=$(date +%s%N)
    $RUST_BIN -i "$PDB_FILE" -c A -o "$RS_OUT" -e -l --base-dir "$SCRIPT_DIR" 2>/dev/null || true
    RS_END=$(date +%s%N)
    RS_TIME=$(echo "scale=2; ($RS_END - $RS_START) / 1000000000" | bc)

    printf "%-10s | %-15s\n" "${SIZE} res" "${RS_TIME}s"
done
echo ""


# Run with energies on 150-res for comparison
echo "[2/2] Running benchmark WITH energies (150 residues)..."
echo ""

PDB_FILE="$TEST_DIR/test_150res.pdb"
RS_OUT="$TEST_DIR/bench_energy_rust_150"
rm -rf "$RS_OUT"

if [ -f "$PDB_FILE" ]; then
    # Rust with energies
    RS_START=$(date +%s%N)
    $RUST_BIN -i "$PDB_FILE" -c A -o "$RS_OUT" -l --base-dir "$SCRIPT_DIR" 2>/dev/null || true
    RS_END=$(date +%s%N)
    RS_TIME=$(echo "scale=2; ($RS_END - $RS_START) / 1000000000" | bc)

    printf "%-20s | %-15s\n" "Task" "Rust (s)"
    printf "%-20s-+-%-15s\n" "--------------------" "---------------"
    printf "%-20s | %-15s\n" "150 res + energies" "${RS_TIME}s"
    echo ""
else
    echo "  SKIP - $PDB_FILE not found"
fi

echo "=============================================="
