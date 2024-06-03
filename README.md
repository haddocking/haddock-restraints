# haddock-restraints

A standalone command-line application to generate restraints to be used in HADDOCK.

## Planned features

- [x] Generate `.tbl` files from an input file
- [x] Define passive residues based on surface accessibility (`surface_as_passive`)
- [x] Define passive residues around active ones (`passive_from_active`)
- [x] Support for N interactors; 2-body, 3-body, 4-body, etc
- [x] Support for multiple interaction sites in the same interactor
- [ ] Filter out active residues that are not accessible
- [X] Generate _true-interface_ restraints for benchmarking
- [ ] Generate random-restraints
- [ ] Create unambiguous restraints to keep molecules together during docking

## Usage

### Clone

```bash
git clone https://github.com/rvhonorato/haddock-restraints.git && cd haddock-restraints
```

### Build

```bash
cargo build --release
cp target/release/haddock-restraints .
```

### Execute

```bash
$ ./haddock-restraints -h
A tool to process different input files

Usage: haddock-restraints <COMMAND>

Commands:
  tbl   Generate TBL file from input file
  ti    Generate true-interface restraints from a PDB file
  help  Print this message or the help of the given subcommand(s)

Options:
  -h, --help  Print help
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
$ ./haddock-restraints ti examples/2oob.pdb 5.0 > ti.tbl
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
Check the `examples/` folder for examples of restraint files.

The mandatory fields are:
- `id`: an integer that identifies the interactor
- `chain`: the chain of the interactor
- `active`: a list of residues that are active in the interaction
- `passive`: a list of residues that are passive in the interaction
- `target`: a list of integers that identifies the interactors that the current interactor interacts with

Optional fields are:
- `structure`: the PDB file that contains the structure of the interactor
- `passive_from_active`: if true, the passive residues are defined based on the active residues (*requires structure*)
- `surface_as_passive`: if true, the passive residues are defined based on the surface accessibility of the residues (*requires structure*)

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

###  `/usr/bin/ld: cannot find -lc++: No such file or directory`

```bash
sudo apt-get install libc++-dev libc++abi-dev
```

### `Unable to find libclang`

```bash
sudo apt-get install libclang-dev
```
