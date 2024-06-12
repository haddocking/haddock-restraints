# Generate `.tbl` file from a configuration file

`haddock-restraints` can generate a `.tbl` file from a configuration file. The configuration file is a JSON file that contains the information needed to generate the restraints.

To use this file you first need to get familiar with a few concepts:

- [Active/Passive residues](./active_passive.md)

- [Configuration file](./configuration_file.md)

## Usage

To run the `tbl` subcommand, you just need to provide the path to the configuration file. For example:

```bash
haddock-restraints tbl path/to/config.json > restraints.tbl
```
