#!/usr/bin/env bash
set -eu
src=$(cd "$(dirname "$0")" && pwd)
: "${BSR_DIR:?Set BSR_DIR to the BSR build containing bin, lib and mod}"
bsr=$(cd "$BSR_DIR" && pwd)
fc=${FC:-gfortran}
out="$src/BIN"

for program in mult3 bsr_dmat3; do
    test -x "$bsr/bin/$program" || { echo "Missing $bsr/bin/$program" >&2; exit 1; }
done
mkdir "$out"
com="$bsr/lib/libzcom.a"
conf="$bsr/lib/libzconfLS.a"
build() {
    name=$1
    shift
    "$fc" -O2 -ffree-line-length-none -I"$bsr/mod" "$src/$name.f90" "$@" -o "$out/$name"
}
build bsr_Splinedata "$bsr/lib/libBS.a" "$conf" "$com"
build bsr_Splinewaves "$com"
build bsr_H_jk "$com"
build bsr_d_jk "$conf" "$com"
ln -s "$bsr/bin/mult3" "$out/mult3"
ln -s "$bsr/bin/bsr_dmat3" "$out/bsr_dmat3"
