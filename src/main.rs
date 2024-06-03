mod air;
mod input;
mod interactor;
mod sasa;
mod structure;
use air::Air;
use interactor::Interactor;
use std::error::Error;

use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(name = "haddock-restraints")]
#[command(about = "A tool to process different input files", long_about = None)]
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
    },
    #[command(about = "Generate true-interface restraints from a PDB file")]
    Ti {
        #[arg(help = "PDB file")]
        input: String,
        #[arg(help = "Cutoff distance for interface residues")]
        cutoff: f64,
    },
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let cli = Cli::parse();

    match &cli.command {
        Commands::Tbl { input } => {
            gen_tbl(input);
        }
        Commands::Ti { input, cutoff } => {
            let _ = true_interface(input, cutoff);
        }
    }

    Ok(())
}

fn gen_tbl(input_file: &str) {
    println!("Reading input file: {}", input_file);
    let mut interactors = input::read_json_file(input_file).unwrap();

    let wd = std::path::Path::new(input_file)
        .parent()
        .unwrap()
        .to_str()
        .unwrap();

    interactors.iter_mut().for_each(|interactor| {
        if !interactor.structure().is_empty() {
            interactor.set_structure(format!("{}/{}", wd, interactor.structure()).as_str());

            if interactor.passive_from_active() {
                interactor.set_passive_from_active();
            }

            if interactor.surface_as_passive() {
                interactor.set_surface_as_passive();
            }
        }
    });

    let air = Air::new(interactors);

    let _tbl = air.gen_tbl().unwrap();
    println!("{}", _tbl);
}

fn true_interface(input_pdb: &str, cutoff: &f64) -> Result<(), Box<dyn Error>> {
    println!("Reading PDB file: {}", input_pdb);

    let pdb = match pdbtbx::open_pdb(input_pdb, pdbtbx::StrictnessLevel::Loose) {
        Ok((pdb, _warnings)) => pdb,
        Err(e) => {
            panic!("Error opening PDB file: {:?}", e);
        }
    };

    let true_interface = structure::get_true_interface(&pdb, *cutoff);
    let chains_in_contact = structure::get_chains_in_contact(&pdb, *cutoff);

    // Sort the true_interface by chain id
    let mut true_interface: Vec<_> = true_interface.iter().collect();
    true_interface.sort_by(|a, b| a.0.cmp(b.0));

    let mut interactors: Vec<Interactor> = Vec::new();
    for (index, (chain_id, residues)) in true_interface.iter().enumerate() {
        let mut interactor = Interactor::new(index as u16);
        interactor.set_chain(chain_id);
        interactor.set_active(residues.iter().map(|&residue| residue as i16).collect());

        // Assign the targets
        for (chain_i, chain_j) in chains_in_contact.iter() {
            let target_chain = if chain_i == *chain_id {
                chain_j
            } else {
                chain_i
            };
            if let Some(target_index) = pdb.chains().position(|chain| chain.id() == *target_chain) {
                interactor.add_target(target_index as u16);
            }
        }

        interactors.push(interactor);
    }
    interactors.iter().for_each(|interactor| {
        println!("{:?}", interactor);
    });

    // Make the restraints
    let air = Air::new(interactors);
    let tbl = air.gen_tbl().unwrap();

    println!("{}", tbl);

    Ok(())
}
