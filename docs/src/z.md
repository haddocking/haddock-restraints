# Generate `z` restraints to keep the molecule aligned to the Z-axis

This subcommand generates restraints to keep the molecule aligned to the Z-axis. This is useful when you want to keep the molecule in a specific orientation during docking.

As input you need to pass at least one selection of residues, this will be used to define a plane perpendicular to the Z-axis. The restraints will be generated to keep the molecule aligned to this plane.

## Usage

To run the `z` subcommand, you need to provide the path to the PDB file, the selection of residues, the output file which will contain the shape beads, the number of beads to generate, and the distance between the beads. For example:

```bash
./haddock-restraints z \
  --residues 19,83,145,119,167 \
  path/to/the/input.pdb \
  path/to/the/output/shape.pdb \
  20 \ # spacing between the beads in angstrom - 20A
  6 \  # grid size in dimension - 6x6
  > z_restraints.tbl
```

Further, follow the HADDOCK documentation **\<pending\>** to use the generated restraints in the docking protocol.
