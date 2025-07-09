# haddock-restraints

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13362093.svg)](https://doi.org/10.5281/zenodo.13362093)
[![Crates.io Version](https://img.shields.io/crates/v/haddock-restraints)](https://crates.io/crates/haddock-restraints)
[![Crates.io License](https://img.shields.io/crates/l/haddock-restraints)](https://crates.io/crates/haddock-restraints)
[![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F-green)](https://fair-software.eu)

[![tests](https://github.com/haddocking/haddock-restraints/actions/workflows/test.yml/badge.svg)](https://github.com/haddocking/haddock-restraints/actions/workflows/test.yml)
![docs.rs](https://img.shields.io/docsrs/haddock-restraints)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/cc008f968e394457ae63650cccfd27da)](https://app.codacy.com/gh/haddocking/haddock-restraints/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![Crates.io Total Downloads](https://img.shields.io/crates/d/haddock-restraints)](https://crates.io/crates/haddock-restraints)

## Usage

Go to [wenmr.science.uu.nl/haddock-restraints](https://wenmr.science.uu.nl/haddock-restraints) for a graphical user interface

Check [bonvinlab.org/haddock-restraints](https://bonvinlab.org/haddock-restraints) for a **user guide** on how to use the code as a command-line application.

Visit [docs.rs/haddock-restraints](https://docs.rs/haddock-restraints) for the **developer documentation**, and how to use it as a library in your code.

See [haddocking/haddock-restraints-wasm](https://github.com/haddocking/haddock-restraints-wasm) for the **web assembly** bindings.


## Commands

- [`tbl`: Generates a TBL file](https://www.bonvinlab.org/haddock-restraints/tbl.html)
- [`ti`: Generate true-interface restraints from a PDB file](https://www.bonvinlab.org/haddock-restraints/ti.html)
- [`restraint`: Generate Unambiguous restraints to keep molecule (including heteroatoms) together during docking](https://www.bonvinlab.org/haddock-restraints/restraint.html)
- [`interface`: List residues in the interface](https://www.bonvinlab.org/haddock-restraints/interface.html)
- [`z`: Generate Z-restraints for a protein](https://www.bonvinlab.org/haddock-restraints/z.html)
- [`unambig-ti`: Generate unambiguous true-interface restraints from a PDB file](https://www.bonvinlab.org/haddock-restraints/unambig-ti.html)

## Install

- Download the [latest binary in the release page](https://github.com/haddocking/haddock-restraints/releases/latest)

OR

- Install it with [`cargo`](https://www.rust-lang.org/tools/install)

  ```bash
  cargo install haddock-restraints
  ```

## Execute

```bash
$ haddock-restraints -h
Generate restraints to be used in HADDOCK

Usage: haddock-restraints <COMMAND>

Commands:
  tbl         Generate TBL file from input file
  ti          Generate true-interface restraints from a PDB file
  unambig-ti  Generate unambiguous true-interface restraints from a PDB file
  restraint   Generate unambiguous restraints to keep molecules together during docking
  interface   List residues in the interface
  z           Generate Z-restraints for a protein
  help        Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```

## Planned features

- [x] Generate `.tbl` files from an input file (tbl command)
- [x] Define passive residues based on surface accessibility (tbl command - `surface_as_passive`)
- [x] Define passive residues around active ones (tbl command - `passive_from_active`)
- [x] Support for N interactors; 2-body, 3-body, 4-body, etc (tbl command)
- [x] Support for multiple interaction sites in the same interactor (tbl command)
- [x] Generate _true-interface_ restraints for benchmarking (ti command)
- [x] Create unambiguous restraints to keep molecules together during docking (restraint command)
- [x] Filter out buried residues (tbl command)
- [x] List residues in the interface (interface command)
- [x] Add Z-restraints to keep molecules aligned in the Z-axis (z command)
- [x] Specify atom subsets
- [ ] Template based restraints
- [ ] ~Generate random-restraints~ done via CNS

## Troubleshooting

### `Unable to find libclang`

```bash
sudo apt-get install libclang-dev
```
