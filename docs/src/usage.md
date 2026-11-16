# Usage

`haddock-restraints` is a command-line tool that is organized with subcommands.
Each subcommand has its own set of options and arguments. The main subcommands are:

- [`tbl`: Generate .tbl file from a configuration file](./tbl.md)

- [`ti`: Generate _true-interface_ restraints from a PDB file](./ti.md)

- [`restraint`: Generate _unambiguous restraints_ to keep molecules together during docking](./restraint.md)

- [`interface`: List residues in the interface](./interface.md)

- [`z`: Generate restraints to keep the molecule aligned to the Z-axis](./z.md)

- [`unambig-ti`: Generate unambiguous true-interface restraints](./unambig-ti.md)

- [`tbl2pml`: Generate a PyMOL visualization from a `.tbl` file and its PDB(s)](./tbl2pml.md)

## Structures

- Only the first model of a multi-model PDB file (ensemble) is used.
- The surface accessibility (`filter_buried`, `surface_as_passive`) is computed as in HADDOCK3:
  - with the NACCESS atom radii of HADDOCK3, which include the modified residues it supports (e.g. NEP, SEP, PTR), and the element radius for other atoms;
  - hydrogens and `HETATM` records are ignored, also when the file has no element column;
  - a residue is on the surface when the relative accessibility of its main chain or of its side chain is at least 0.15 (15%).
