# haddock-restraints

🚧 🚧 🚧

A standalone command-line application to generate restraints to be used in HADDOCK.

🚧 🚧 🚧

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
./haddock-restraints examples/restraints.json
```

### Input

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


## Planned features

- [x] Generate `.tbl` files from an input file
- Filter out residues that are not accessible
- [x] Define passive residues based on surface accessibility
- Generate _true-interface_ restraints for benchmarking
- [x] Support for N interactors; 2-body, 3-body, 4-body, etc
- [x] Support for multiple interaction sites in the same interactor
- Generate random-restraints

## Troubleshooting

###  `/usr/bin/ld: cannot find -lc++: No such file or directory`

```bash
sudo apt-get install libc++-dev libc++abi-dev
```

### `Unable to find libclang`

```bash
sudo apt-get install libclang-dev
```