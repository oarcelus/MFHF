# MFHF - Mean-field Hartree-Fock solver

MFHF is a Fortran research code for solving low-energy Hubbard-type models in the self-consistent mean-field Hartree-Fock approximation. It starts from a localized/Wannier-basis Hamiltonian and local Coulomb interactions prepared from a prior DFT calculation (Quantum-Espresso + Wannier90 workflow). The solver produces the electronic structure, density matrices, and Hartree-Fock potentials used by the companion HFPP post-processing code.

## Features

- Self-consistent Hartree-Fock solution of a low-energy Hubbard model.
- Wannier-basis hopping Hamiltonian and a local four-index Coulomb interaction matrix.
- Collinear calculations and a spin-orbit branch controlled by `spinorb` in `hf.in`.
- Explicit Coulomb matrices or spherical interaction parameterization.
- Optional density-of-states and band-structure output.

The program entry point is `main.f90`. Input/model handling is in `parameters.f90`; the Hartree-Fock Hamiltonian, density, energy, and plotting routines are in `hamiltonian.f90`, `density.f90`, `energy.f90`, and `plot.f90`.

## Requirements and build

The supplied `Makefile` is configured for Intel Fortran and Intel MKL:

```sh
make
```

It invokes `ifort -mkl` and uses LAPACK diagonalization. The recipe does not specify `-o hf.x`, so the generated executable name is compiler- and platform-dependent. Run the executable actually produced by the compiler, or add an explicit output name in a local build configuration.

## Inputs and running

The executable reads fixed filenames from its current working directory. Start from a copy of `Examples/NaFePO4/` or create a separate calculation directory containing the required inputs.

| File | Role |
| --- | --- |
| `hf.in` | Positional calculation controls |
| `hf_struct.in` | Lattice vectors and site positions |
| `hf_kpt.in` | k-point mesh/list |
| `hf_ham.in` | Wannier-basis hopping Hamiltonian |
| `hf_mag.in` | Initial magnetic/density-matrix data |
| `hf_coulomb.in` | Local Coulomb matrix; required when `lsph = .false.` |
| `hf_soc.in` | Spin-orbit input when required by the selected SOC mode |
| `hf_path.in` | Band path when `bandplot = .true.` |

After the ignored heading line, `hf.in` gives the following values in order: `spinorb`, `soctype`, `bondcoord`, `nelec`, `nscf`, `nener`, `etemp`, `ediff`, `mix`, `dosplot`, `smear`, `bandplot`, and `lsph`. The provided NaFePO4 example is the best reference for the expected fixed-format inputs.

Run the built executable from that calculation directory. Output files are written with replacement semantics, so preserve results or use a new directory before re-running.

## Outputs

MFHF always writes:

- `hf.out` - self-consistency log, energy, occupations, and magnetization.
- `eig.out` and `wav.out` - eigenvalues/Fermi energy and eigenvectors.
- `bond.out` - real-space bond information.
- `pot.out` and `dens.out` - converged Hartree-Fock potential and density matrix.

These files form the main interface to HFPP. Depending on the requested calculation, MFHF can also write `dos.out`, `eigpath.out`, `wavpath.out`, `bands.1.out`, `bands.2.out`, or `bands.out`. Density-of-states output is not implemented for `spinorb = .true.` in the current code.

## Scientific background

1. I. V. Solovyev, ["Combining DFT and many-body methods to understand correlated materials"](https://iopscience.iop.org/article/10.1088/0953-8984/20/29/293201), *Journal of Physics: Condensed Matter* **20**, 293201 (2008). https://doi.org/10.1088/0953-8984/20/29/293201
2. I. V. Solovyev, ["Self-consistent linear response for the spin-orbit interaction related properties"](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.024417), *Physical Review B* **90**, 024417 (2014). https://doi.org/10.1103/PhysRevB.90.024417

## License

This project is distributed under the [MIT License](LICENSE).
