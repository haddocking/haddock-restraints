# Installation

`haddock-restraints` is a Rust command-line tool. To install it, you need to have the Rust toolchain installed on your system. You can install Rust by following the instructions on the [official website](https://www.rust-lang.org/tools/install).

It is also published in [Crates.io](https://crates.io/crates/haddock-restraints), so you can install it using the following command using Cargo, the Rust package manager:

```bash
cargo install haddock-restraints
```

Now it should be available in your system. You can check if it is installed by running:

```bash
haddock-restraints -h
```

## From source

If you prefer to install it from source, you can clone the repository and build it using Cargo:

```bash
git clone https://github.com/haddocking/haddock-restraints
cd haddock-restraints
cargo build --release
```

***

Go ahead and proceed to the [USAGE](./usage.md) section to learn how to use the tool.
