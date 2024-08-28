# `tbl` configuration file

Restraints are very powerful tools to guide the docking process, to use its full potential, in `haddock-restraints` we use a configuration file to define the restraints. This file is a JSON file that contains the information needed to generate the restraints.

To use the full potential of ambiguous restraints we have introduced the concept of _interactors_, which are arbitrary groups of residues that can be used to define the restraints.

`haddock-restraints` supports an unlimited number of interactors, and each interactor can have an unlimited number of residues. This allows you to define complex restraints that can be used to guide the docking process.

- [Mandatory fields](#mandatory-fields)
- [Optional fields](#optional-fields)
- [Reference guide](#reference-guide) with detailed explanation of the fields in the configuration file

***

## Mandatory fields

- `id`: an integer that identifies the interactor
- `chain`: the chain of the interactor
- `active`: a list of residues that are active in the interaction
- `passive`: a list of residues that are passive in the interaction
- `target`: a list of integers that identifies the interactors that the current interactor interacts with

A minimal configuration file would look like this:

```json
[
  {
    "id": 1,
    "chain": "A",
    "active": [950],
    "passive": [],
    "target": [2]
  },
  {
    "id": 2,
    "chain": "B",
    "active": [41, 42, 43, 44, 45],
    "passive": [],
    "target": [1]
  }
]
```

## Optional fields

- `structure`: the PDB file that contains the structure of the interactor
- `passive_from_active`: if true, the passive residues are defined based on the active residues (_requires structure_)
- `surface_as_passive`: if true, the passive residues are defined based on the surface accessibility of the residues (_requires structure_)
- `filter_buried`: if true, the buried residues are filtered out (_requires structure_)
- `filter_buried_cutoff`: the cutoff to consider a residue as buried, default = 0.7 (_requires structure_)
- `target_distance`: the distance to consider two residues as interacting, default = 2.0
- `lower_margin`: the lower margin to consider two residues as interacting, default = 2.0
- `upper_margin`: the upper margin to consider two residues as interacting, default = 0.0

Note: Check [this paper](https://doi.org/10.1038/s41596-018-0018-5) for a deeper explanation about the target distance and the margins

A configuration file with optional fields would look like this:

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

## Reference guide

- `id`
  - type: integer
  - description: an integer that identifies the interactor

- `chain`
  - type: string
  - description: the chain of the interactor

- `active`
  - type: list of integers
  - description: a list of residues that are active in the interaction

- `passive`
  - type: list of integers
  - description: a list of residues that are passive in the interaction

- `target`
  - type: list of integers
  - description: a list of integers that identifies the interactors that the current interactor interacts with

Optional fields are:

- `structure`
  - type: string
  - description: the PDB file that contains the structure of the interactor. If using relative paths, they should be relative to the configuration file

- `passive_from_active`
  - type: boolean
  - description: define passive residues are defined based on the active residues
  - requires: `structure`

- `surface_as_passive`
  - type: boolean
  - description: define passive residues are defined based on the surface accessibility of the residues
  - requires: `structure`

- `filter_buried`
  - type: boolean
  - description: filter out buried residues
  - requires: `structure`

- `filter_buried_cutoff`
  - type: float
  - description: the cutoff to consider a residue as buried, default = 0.7
  - requires: `structure`
