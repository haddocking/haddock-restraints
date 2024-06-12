# Generate unambiguous restraints to keep molecules together during docking

In some scenarios, you might be docking proteins that contain **structural gaps**, but you still want to keep the molecules together during the docking process. This is where the `restraint` subcommand comes in handy.

It will generate **unambiguous restraints** to keep the molecules together during docking.

## Usage

To run the `restraint` subcommand, you just need to provide the path to the PDB file. For example:

```bash
haddock-restraints restraint path/to/complex.pdb > unambig.tbl
```
