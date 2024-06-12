# Generate _true-interface_ restraints from a PDB file

This is a very specific type of restraint, that is used to restrain the interface of a protein-protein complex most commonly used to benchmark the efficiency of a docking workflow or protocol.

What this command does is to calculate the distance between the residues of two chains and if the distance is less than a given cutoff, the residues are considered to be in the interface.

Based on this, `haddock-restraints` fills in the `active` and `passive` fields and provides you with a `.tbl` file that can be used to restrain the interface of the protein-protein complex.

## Usage

To run the `ti` subcommand, you just need to provide the path to the PDB file and the cutoff distance. For example:

```bash
haddock-restraints ti path/to/complex.pdb 5.0 > ti.tbl
```
