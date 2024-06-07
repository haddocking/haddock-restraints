# haddock-restraints

![Crates.io Version](https://img.shields.io/crates/v/haddock-restraints)
![Crates.io Total Downloads](https://img.shields.io/crates/d/haddock-restraints)
![Crates.io License](https://img.shields.io/crates/l/haddock-restraints)

A standalone command-line application to generate restraints to be used in HADDOCK.

## Planned features

- [x] Generate `.tbl` files from an input file
- [x] Define passive residues based on surface accessibility (`surface_as_passive`)
- [x] Define passive residues around active ones (`passive_from_active`)
- [x] Support for N interactors; 2-body, 3-body, 4-body, etc
- [x] Support for multiple interaction sites in the same interactor
- [X] Generate _true-interface_ restraints for benchmarking
- [x] Create unambiguous restraints to keep molecules together during docking
- [ ] Specify atom subsets
- [ ] Filter out active residues that are not accessible
- [ ] ~Generate random-restraints~ done via CNS

## Usage

### Install

```bash
cargo install haddock-restraints
```

### Execute

```bash
$ haddock-restraints -h
Generate restraints to be used in HADDOCK

Usage: haddock-restraints <COMMAND>

Commands:
  tbl   Generate TBL file from input file
  ti    Generate true-interface restraints from a PDB file
  help  Print this message or the help of the given subcommand(s)

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

```json
[
  {
    "id": 1,
    "chain": "A",
    "active": [
      934
    ],
    "passive": [],
    "structure": "2oob.pdb",
    "target": [
      2
    ],
    "passive_from_active": true
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

## Troubleshooting

### `/usr/bin/ld: cannot find -lc++: No such file or directory`

```bash
sudo apt-get install libc++-dev libc++abi-dev
```

### `Unable to find libclang`

```bash
sudo apt-get install libclang-dev
```
