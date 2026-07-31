# Visualizing restraints with PyMOL (`--pml`)

The `tbl`, `ti`, `unambig-ti`, and `restraint` subcommands all accept a `--pml` option that, in
addition to generating the restraints file, writes out a PyMOL script (`.pml`) to visualize the
restraints network.

The script colors **active** residues red, **passive** residues green, and draws dashed lines
between the residues involved in each restraint, making it easy to inspect whether the network
makes sense before starting a docking run.

## Usage

Pass `--pml` with an output filename to any of the supported subcommands. For example:

```bash
haddock-restraints tbl path/to/config.json --pml network.pml > restraints.tbl
```

Then load the structure and the script together in PyMOL:

```bash
pymol path/to/structure.pdb network.pml
```

or, from within an already open PyMOL session:

```
load path/to/structure.pdb
@network.pml
```
