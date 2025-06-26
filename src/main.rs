use clap::{Parser, Subcommand};

mod cli;

#[derive(Parser)]
#[command(name = "haddock-restraints")]
#[command(about = "Generate restraints to be used in HADDOCK", long_about = None, version)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    #[command(about = "Generate TBL file from input file")]
    Tbl {
        #[arg(help = "Input file")]
        input: String,

        #[arg(
            long,
            help = "PyMol Script (.pml) output file",
            value_name = "output.pml"
        )]
        pml: Option<String>,
    },
    #[command(about = "Generate true-interface restraints from a PDB file")]
    Ti {
        #[arg(help = "PDB file")]
        input: String,
        #[arg(help = "Cutoff distance for interface residues")]
        cutoff: f64,
        #[arg(
            long,
            help = "PyMol Script (.pml) output file",
            value_name = "output.pml"
        )]
        pml: Option<String>,
    },
    #[command(about = "Generate unambiguous true-interface restraints from a PDB file")]
    UnambigTi {
        #[arg(help = "PDB file")]
        input: String,
        #[arg(help = "Cutoff distance for interface residues")]
        cutoff: f64,
        #[arg(
            long,
            help = "PyMol Script (.pml) output file",
            value_name = "output.pml"
        )]
        pml: Option<String>,
    },
    #[command(about = "Generate unambiguous restraints to keep molecules together during docking")]
    Restraint {
        #[arg(help = "PDB file")]
        input: String,
        #[arg(
            long,
            help = "PyMol Script (.pml) output file",
            value_name = "output.pml"
        )]
        pml: Option<String>,
    },
    #[command(about = "List residues in the interface")]
    Interface {
        #[arg(help = "PDB file")]
        input: String,
        #[arg(help = "Cutoff distance for interface residues")]
        cutoff: f64,
    },
    #[command(about = "Generate Z-restraints for a protein")]
    Z {
        #[arg(required = true, help = "Input file")]
        input: String,
        #[arg(
            required = true,
            help = "Filename of the output shape beads to be used in Z-restraints"
        )]
        output: String,
        #[arg(long, required = true, help = "Group of residue indexes (can be specified multiple times)", value_parser = parse_residues, number_of_values = 1)]
        residues: Vec<Vec<isize>>,
        #[arg(
            required = true,
            help = "Spacing between two beads in Angstrom",
            default_value = "2.0"
        )]
        grid_spacing: f64,
        #[arg(required = true, help = "Size in xy dimension", default_value = "10")]
        grid_size: usize,
    },
}

// Parse a comma-separated list of residues
fn parse_residues(arg: &str) -> Result<Vec<isize>, String> {
    arg.split(',')
        .map(|s| {
            s.trim()
                .parse::<isize>()
                .map_err(|_| format!("Invalid number: {}", s))
        })
        .collect()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    match &cli.command {
        Commands::Tbl { input, pml } => {
            cli::handle_gen_tbl(input, pml);
        }
        Commands::Ti { input, cutoff, pml } => {
            cli::handle_ti(input, cutoff, pml);
        }
        Commands::UnambigTi { input, cutoff, pml } => {
            cli::handle_unambig_ti(input, cutoff, pml);
        }
        Commands::Restraint { input, pml } => {
            cli::handle_restraint_bodies(input, pml);
        }
        Commands::Interface { input, cutoff } => {
            cli::handle_list_interface(input, cutoff);
        }
        Commands::Z {
            input,
            output,
            residues,
            grid_size,
            grid_spacing,
        } => {
            cli::handle_generate_z_restraints(input, output, residues, grid_size, grid_spacing);
        }
    }

    Ok(())
}
