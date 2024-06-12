# Welcome to the `haddock-restraints` guide

![image](./banner_home-mini.jpg)

`haddock-restraints` is a command-line tool that helps you to generate restraints for the HADDOCK docking software.

The generation of restraints is a crucial step in the preparation of a docking experiment. The restraints define the spatial constraints that guide the docking software in the search for the best conformation of the complex.

[HADDOCK](https://www.bonvinlab.org/software/#haddock) uses a `.tbl` format to define the restraints. If you are using the [web service](https://wenmr.science.uu.nl), these are defined automatically via the web interface. However, if you are running HADDOCK locally, you need to provide these restraints yourself.

This is where `haddock-restraints` comes in. It helps you to generate the restraints in the `.tbl` format, based on the input structures you provide.

See the [INSTALLATION](./installation.md) section to get started and proceed to the [USAGE](./usage.md) section to learn how to use the tool.

## Getting help

If you encounter any issues or have any questions, please open an issue on the [GitHub repository](https://github.com/haddocking/haddock-restraints), contact us at _bonvinlab.support@uu.nl_ or join the [BioExcel forum](https://ask.bioexcel.eu) and post your question there.
