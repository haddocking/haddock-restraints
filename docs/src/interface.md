# List residues in the interface

This command is a helpful when you want to know which residues are in the interface of a protein-protein complex. It calculates the distance between the residues of two chains and if the distance is less than a given cutoff, the residues are considered to be in the interface.

> The result of this might be useful to generate restraints using the [`tbl` command](./tbl.md) - to define active and/or passive residues.

## Usage

To run the `interface` subcommand, you just need to provide the path to the PDB file and the cutoff distance. For example:

```bash
haddock-restraints interface path/to/complex.pdb 5.0
```
