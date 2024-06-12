# haddock-restraints

[![Crates.io Version](https://img.shields.io/crates/v/haddock-restraints)](https://crates.io/crates/haddock-restraints)
[![Crates.io Total Downloads](https://img.shields.io/crates/d/haddock-restraints)](https://crates.io/crates/haddock-restraints)
[![Crates.io License](https://img.shields.io/crates/l/haddock-restraints)](https://crates.io/crates/haddock-restraints)

[![tests](https://github.com/haddocking/haddock-restraints/actions/workflows/test.yml/badge.svg)](https://github.com/haddocking/haddock-restraints/actions/workflows/test.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/cc008f968e394457ae63650cccfd27da)](https://app.codacy.com/gh/haddocking/haddock-restraints/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)

A standalone command-line application to generate restraints to be used in HADDOCK.

## Commands

- [`tbl`: Generate TBL file from input file](#tbl-command)
- [`ti`: Generate true-interface restraints from a PDB file](#ti-command)
- [`restraint`: Generate Unambiguous restraints to keep molecules together during docking](#restraint-command)
- [`interface`: List residues in the interface](#interface-command)

## Planned features

- [x] Generate `.tbl` files from an input file
- [x] Define passive residues based on surface accessibility (`surface_as_passive`)
- [x] Define passive residues around active ones (`passive_from_active`)
- [x] Support for N interactors; 2-body, 3-body, 4-body, etc
- [x] Support for multiple interaction sites in the same interactor
- [x] Generate _true-interface_ restraints for benchmarking
- [x] Create unambiguous restraints to keep molecules together during docking
- [x] Filter out buried residues
- [x] List residues in the interface
- [ ] Template based restraints
- [ ] Specify atom subsets
- [ ] ~Generate random-restraints~ done via CNS

## Usage

### Install

> [Install Rust](https://www.rust-lang.org/tools/install) if you don't have it yet - later we will provide the execution binaries 🤓

```bash
cargo install haddock-restraints
```

### Execute

```bash
$ haddock-restraints -h
Generate restraints to be used in HADDOCK

Usage: haddock-restraints <COMMAND>

Commands:
  tbl        Generate TBL file from input file
  ti         Generate true-interface restraints from a PDB file
  restraint  Generate Unambiguous restraints to keep molecules together during docking
  interface  List residues in the interface
  help       Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```

### `ti` command

```bash
$ ./haddock-restraints ti -h
Generate true-interface restraints from a PDB file

Usage: haddock-restraints ti <INPUT> <CUTOFF>

Arguments:
  <INPUT>   PDB file
  <CUTOFF>  Cutoff distance for interface residues

Options:
  -h, --help  Print help
```

Example:

```bash
./haddock-restraints ti examples/2oob.pdb 5.0 > ti.tbl
```

### `tbl` command

```bash
$ ./haddock-restraints tbl -h
Generate TBL file from input file

Usage: haddock-restraints tbl <INPUT>

Arguments:
  <INPUT>  Input file

Options:
  -h, --help  Print help
```

Example:

```bash
./haddock-restraints tbl examples/restraints.json > ambig.tbl
```

Check the [examples](https://github.com/rvhonorato/haddock-restraints/tree/main/examples) folder for examples of restraint files.

The mandatory fields are:

- `id`: an integer that identifies the interactor
- `chain`: the chain of the interactor
- `active`: a list of residues that are active in the interaction
- `passive`: a list of residues that are passive in the interaction
- `target`: a list of integers that identifies the interactors that the current interactor interacts with

Optional fields are:

- `structure`: the PDB file that contains the structure of the interactor
- `passive_from_active`: if true, the passive residues are defined based on the active residues (_requires structure_)
- `surface_as_passive`: if true, the passive residues are defined based on the surface accessibility of the residues (_requires structure_)
- `filter_buried`: if true, the buried residues are filtered out (_requires structure_)
- `filter_buried_cutoff`: the cutoff to consider a residue as buried, default = 0.7 (_requires structure_)

```json
[
  {
    "id": 1,
    "chain": "A",
    "active": [
      934,
      939
    ],
    "passive": [],
    "structure": "2oob.pdb",
    "target": [
      2
    ],
    "passive_from_active": true,
    "filter_buried": true
  },
  {
    "id": 2,
    "chain": "B",
    "active": [
      68
    ],
    "passive": [],
    "target": [
      1
    ]
  },
  {
    "id": 3,
    "chain": "B",
    "active": [],
    "passive": [],
    "target": [
      1
    ],
    "structure": "2oob.pdb",
    "surface_as_passive": true
  }
]
```

### `restraint` command

```bash
$ ./haddock-restraints restraint -h
Generate Unambiguous restraints to keep molecules together during docking

Usage: haddock-restraints restraint <INPUT>

Arguments:
  <INPUT>  PDB file

Options:
  -h, --help  Print help
```

Example:

```bash
./haddock-restraints restraint examples/2oob_w_gaps.pdb > unambiguous.tbl
```

### `interface` command

```bash
$ ./haddock-restraints interface -h
List residues in the interface

Usage: haddock-restraints interface <INPUT> <CUTOFF>

Arguments:
  <INPUT>   PDB file
  <CUTOFF>  Cutoff distance for interface residues

Options:
  -h, --help  Print help
```

Example:

```bash
./haddock-restraints interface examples/2oob.pdb 5.0

Chain A: [931, 933, 934, 936, 937, 938, 940, 941, 946, 950]
Chain B: [6, 8, 42, 44, 45, 46, 47, 48, 49, 66, 68, 69, 70]
```

***

## Troubleshooting

### `/usr/bin/ld: cannot find -lc++: No such file or directory`

```bash
sudo apt-get install libc++-dev libc++abi-dev
```

### `Unable to find libclang`

```bash
sudo apt-get install libclang-dev
```
