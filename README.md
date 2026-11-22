# BSR–RMT interface utilities

This repository contains the four utilities described in Table 1 of *Multiconfigurational and Spin–Orbit Extensions of the BSR–RMT Interface*. They take a completed set of inner-region BSR data and produce the complete set of atomic input files in the format required by RMT.

- **bsr_Splinedata.f90** writes `Splinedata` using the radial grid in the BSR file `knot.dat`.
- **bsr_Splinewaves.f90** writes `Splinewaves` from the channel coefficients found in the BSR eigenvector files `rsol.nnn`.
- **bsr_H_jk.f90** writes `H` from the BSR Hamiltonian files `h.nnn`.
- **bsr_d_jk.f90** writes the dipole index `D00` and one dipole file (`D001`, `D002`, etc.) for each dipole-allowed pair of symmetries.

## Carbon example

The example is the carbon calculation with spin–orbit coupling included, described in Section 4.2.1 of the paper. There are two symmetry blocks (J = 0, even, and J = 1, odd) which produce one dipole file, `D001`. Figures 2 and 3 show the dipole response and spectrum obtained with these atomic inputs. The corresponding curves are labeled BSR-SO.

The completed BSR files can be found in `example/BSR/`. The five RMT input files to recreate are in `example/RMT/`: `H`, `D00`, `D001`, `Splinedata`, and `Splinewaves`.

## Reproducing the example

To regenerate the RMT input files, the Fano–Racah BSR build used in the paper is required. The matching BSR source files and build instructions are available from the corresponding author. Once the utilities are built, `Splinedata`, `Splinewaves`, and `H` can be produced directly from the supplied BSR data. To prepare the dipole files for RMT, `bsr_d_jk` also calls the BSR programs `mult3` and `bsr_dmat3`.

After BSR is built and `BSR_DIR` is set according to its build instructions, run the following commands from the repository directory. The supplied build script was tested with GNU Fortran.

```bash
bash build.sh
bash run-example.sh
```

The build script creates the utility executables in `BIN/`. The run script generates the RMT input files in `carbon-run/export/atom/` and compares them with the supplied files in `example/RMT/`.

Author credits and references are in `NOTICE.txt`.
