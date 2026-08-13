# Generate a PyMOL visualization from a `.tbl` file

`haddock-restraints` can generate a PyMOL (`.pml`) visualization directly from an existing
`.tbl` restraints file and the PDB(s) it refers to — without needing a `config.json`.

This is the standalone equivalent of the [`--pml`](./pml.md) option available on `tbl`, `ti`,
`unambig-ti`, and `restraint`, for when you already have a `.tbl` file (handed to you, or from a
previous run) instead of the inputs that generated it.

## Usage

```bash
haddock-restraints tbl2pml restraints.tbl complex.pdb --output network.pml
```

Multiple PDBs can be passed if the restraints span more than one structure file; each gets its
own `load` line in the generated script:

```bash
haddock-restraints tbl2pml restraints.tbl chainA.pdb chainB.pdb --output network.pml
```

```bash
pymol network.pml
```

Like the rest of `--pml`, active residues are colored red, passive residues green, with dashed
lines drawn between restrained residue pairs.

> **Note**: `tbl2pml` checks that the PDB(s) you pass exist, but does not validate that the
> `segid` values in the `.tbl` match chain IDs present in them — a mismatch simply renders
> nothing for that residue in PyMOL, the same way an unresolved selection would in any other
> `.pml` script.
