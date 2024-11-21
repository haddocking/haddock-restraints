# Generate _unambiguous true-interface_ restraints from a PDB file

This is a more strict version of the true-interface restraints that assigns
specific pairs to each residue in the interface.

For each atom in each of the chains, the algorithm defines the pairs based on
atomic interactions between any atom of a given residue in the interface and
any other atom in an interacting chain. The restraints are later defined
considering the CA-CA distances of such pairs; thus allowing for
side-chain flexibility.

## Usage

To run the `unambig-ti` subcommand, you just need to provide the path to the
PDB file and the cutoff distance. For example:

```bash
haddock-restraints unambig-ti path/to/complex.pdb 5.0 > unambig.tbl
```
