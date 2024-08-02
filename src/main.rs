mod air;
mod input;
mod interactor;
mod sasa;
mod structure;
use air::Air;
use interactor::Interactor;
use std::{error::Error, vec};

use clap::{Parser, Subcommand};

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
    },
    #[command(about = "Generate true-interface restraints from a PDB file")]
    Ti {
        #[arg(help = "PDB file")]
        input: String,
        #[arg(help = "Cutoff distance for interface residues")]
        cutoff: f64,
    },
    #[command(about = "Generate Unambiguous restraints to keep molecules together during docking")]
    Restraint {
        #[arg(help = "PDB file")]
        input: String,
    },
    #[command(about = "List residues in the interface")]
    Interface {
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
        Commands::Restraint { input } => {
            let _ = restraint_bodies(input);
        }
        Commands::Interface { input, cutoff } => {
            let _ = list_interface(input, cutoff);
        }
    }

    Ok(())
}

fn gen_tbl(input_file: &str) {
    let mut interactors = input::read_json_file(input_file).unwrap();

    interactors.iter_mut().for_each(|interactor| {
        if !interactor.structure().is_empty() {
            if interactor.passive_from_active() {
                interactor.set_passive_from_active();
            }

            if interactor.surface_as_passive() {
                interactor.set_surface_as_passive();
            }

            if interactor.filter_buried() {
                interactor.remove_buried_residues();
            }
        }
    });

    let air = Air::new(interactors);

    let tbl = air.gen_tbl().unwrap();
    println!("{}", tbl);
}

fn true_interface(input_pdb: &str, cutoff: &f64) -> Result<(), Box<dyn Error>> {
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

    // Make the restraints
    let air = Air::new(interactors);
    let tbl = air.gen_tbl().unwrap();

    println!("{}", tbl);

    Ok(())
}

fn restraint_bodies(input_file: &str) -> Result<(), Box<dyn Error>> {
    // Read PDB file
    let pdb = match pdbtbx::open_pdb(input_file, pdbtbx::StrictnessLevel::Loose) {
        Ok((pdb, _warnings)) => pdb,
        Err(e) => {
            panic!("Error opening PDB file: {:?}", e);
        }
    };

    // Find in-contiguous chains
    let gaps = structure::find_structural_gaps(&pdb);

    // Create the interactors
    let mut interactors: Vec<Interactor> = Vec::new();
    let mut counter = 0;
    gaps.iter().for_each(|g| {
        // let (chain, res_i, res_j) = g;
        let mut interactor_i = Interactor::new(counter);
        counter += 1;
        let mut interactor_j = Interactor::new(counter);
        interactor_j.add_target(counter - 1);
        interactor_i.add_target(counter);
        counter += 1;

        interactor_i.set_chain(g.chain.as_str());
        interactor_i.set_active(vec![g.res_i as i16]);
        interactor_i.set_active_atoms(vec![g.atom_i.clone()]);
        interactor_i.set_passive_atoms(vec![g.atom_j.clone()]);
        interactor_i.set_target_distance(g.distance);
        interactor_i.set_lower_margin(0.0);
        interactor_i.set_upper_margin(0.0);

        interactor_j.set_chain(g.chain.as_str());
        interactor_j.set_passive(vec![g.res_j as i16]);

        interactors.push(interactor_i);
        interactors.push(interactor_j);
    });

    let air = Air::new(interactors);
    let tbl = air.gen_tbl().unwrap();

    println!("{}", tbl);

    Ok(())
}

fn list_interface(input_file: &str, cutoff: &f64) -> Result<(), Box<dyn Error>> {
    let pdb = match pdbtbx::open_pdb(input_file, pdbtbx::StrictnessLevel::Loose) {
        Ok((pdb, _warnings)) => pdb,
        Err(e) => {
            panic!("Error opening PDB file: {:?}", e);
        }
    };

    let true_interface = structure::get_true_interface(&pdb, *cutoff);

    for (chain_id, residues) in true_interface.iter() {
        let mut sorted_res = residues.iter().collect::<Vec<_>>();
        sorted_res.sort();

        println!("Chain {}: {:?}", chain_id, sorted_res);
    }

    Ok(())
}
