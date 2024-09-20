# Installation

`haddock-restraints` is a Rust command-line tool. You can either download the pre-compiled binary, download from the Rust package manager or build it from source.

## Download

Get a pre-compiled binary for your system at the link below:

<div style="">
  <style>
    .feedback-button {
      display: inline-block;
      padding: 10px 20px;
      background-color: #337ab7;
      color: white !important;
      text-decoration: none;
      border-radius: 5px;
      font-weight: bold;
      transition: background-color 0.3s ease;
    }
    .feedback-button:hover {
      background-color: #23527c;
    }
  </style>

<a href="https://github.com/haddocking/haddock-restraints/releases/tag/v0.5.0" class="feedback-button" target="_blank">Download</a>

</div>

## Build from Cargo

To install the tool using Cargo, you need to have the Rust toolchain installed on your system. You can install Rust by following the instructions on the [official website](https://www.rust-lang.org/tools/install).

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

---

Go ahead and proceed to the [USAGE](./usage.md) section to learn how to use the tool.
