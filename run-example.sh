#!/usr/bin/env bash
set -eu
pkg=$(cd "$(dirname "$0")" && pwd)
run=${1:-carbon-run}
mkdir "$run"
mkdir -p "$run/export/atom"
cp "$pkg"/example/BSR/* "$run/export/atom/"
ln -s "$pkg/BIN" "$run/deps"
cd "$run/export/atom"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

../../deps/bsr_Splinewaves nlsp=2 > waves.stdout
../../deps/bsr_H_jk ih1=1 ih2=2 > H.stdout
# Establish the angular bank in the order used for the dipole files.
../../deps/mult3 cfg.002 cfg.001 E1 mult_bnk_002001 > seed.stdout
../../deps/bsr_d_jk bin=../../deps klsp1=1 klsp2=2 ipar0=-1 > dipoles.stdout
../../deps/bsr_Splinedata > grid.stdout

for file in H D00 D001 Splinedata Splinewaves; do
    cmp "$file" "$pkg/example/RMT/$file"
    echo "$file matches the supplied RMT file."
done
