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

/// Generates a table of Ambiguous Interaction Restraints (AIRs) based on input data.
///
/// This function reads interactor data from a JSON file, processes the interactors
/// according to specified flags, and generates a table of AIRs.
///
/// # Arguments
///
/// * `input_file` - A string slice that holds the path to the input JSON file.
///
/// # Functionality
///
/// 1. Reads interactor data from the specified JSON file.
/// 2. For each interactor with a non-empty structure:
///    a. If `passive_from_active` is set, derives passive residues from active ones.
///    b. If `surface_as_passive` is set, treats surface residues as passive.
///    c. If `filter_buried` is set, removes buried residues from consideration.
/// 3. Creates an `Air` instance from the processed interactors.
/// 4. Generates a table of AIRs.
/// 5. Prints the generated table to stdout.
///
/// # Panics
///
/// This function will panic if:
/// - The JSON file cannot be read or parsed.
/// - The `Air` instance fails to generate the table.
///
/// # Notes
///
/// - The function assumes that the JSON file is properly formatted and contains valid interactor data.
/// - The output is printed directly to stdout and is not returned or saved to a file.
/// - Error handling is minimal; most errors will result in a panic.
///
/// # Dependencies
///
/// This function relies on several other modules and functions:
/// - `input::read_json_file` for reading the JSON input.
/// - `Air` struct and its methods for processing the interactors and generating the table.
/// - Various methods on the `Interactor` struct for processing individual interactors.
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

/// Analyzes the true interface of a protein structure and generates Ambiguous Interaction Restraints (AIRs).
///
/// This function reads a PDB file, identifies the true interface between chains based on a distance cutoff,
/// creates interactors for each chain involved in the interface, and generates AIRs.
///
/// # Arguments
///
/// * `input_pdb` - A string slice that holds the path to the input PDB file.
/// * `cutoff` - A reference to a f64 value specifying the distance cutoff (in Angstroms) for determining interfaces.
///
/// # Returns
///
/// A `Result<(), Box<dyn Error>>` which is Ok(()) if the function completes successfully, or an Error if something goes wrong.
///
/// # Functionality
///
/// 1. Opens and parses the PDB file.
/// 2. Identifies the true interface residues and chains in contact using the specified cutoff.
/// 3. Creates `Interactor` instances for each chain involved in the interface.
/// 4. Assigns active residues and target chains to each interactor.
/// 5. Generates AIRs using the created interactors.
/// 6. Prints the generated AIR table to stdout.
///
/// # Panics
///
/// This function will panic if:
/// - The PDB file cannot be opened or parsed.
///
/// # Notes
///
/// - The function uses a loose strictness level when parsing the PDB file.
/// - The true interface is determined based on the provided distance cutoff.
/// - Interactors are created only for chains that are part of the identified interface.
/// - The output is printed directly to stdout and is not returned or saved to a file.
///
/// # Dependencies
///
/// This function relies on several other modules and functions:
/// - `pdbtbx` for opening and parsing PDB files.
/// - `structure::get_true_interface` and `structure::get_chains_in_contact` for interface analysis.
/// - `Interactor` struct for representing protein chains and their interactions.
/// - `Air` struct for generating the AIR table.
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

/// Generates distance restraints between non-contiguous bodies in a protein structure.
///
/// This function analyzes a PDB file to identify separate bodies within the protein structure,
/// calculates gaps between these bodies, and generates Ambiguous Interaction Restraints (AIRs)
/// to maintain the relative positions of these bodies.
///
/// # Arguments
///
/// * `input_file` - A string slice that holds the path to the input PDB file.
///
/// # Returns
///
/// A `Result<(), Box<dyn Error>>` which is Ok(()) if the function completes successfully,
/// or an Error if something goes wrong.
///
/// # Functionality
///
/// 1. Opens and parses the PDB file using a loose strictness level.
/// 2. Identifies separate bodies within the protein structure.
/// 3. Calculates gaps between these bodies.
/// 4. For each gap:
///    a. Creates two interactors, one for each side of the gap.
///    b. Sets up the interactors with appropriate chain, residue, and atom information.
///    c. Configures distance restraints based on the gap distance.
/// 5. Generates AIRs using the created interactors.
/// 6. Prints the generated AIR table to stdout.
///
/// # Panics
///
/// This function will panic if:
/// - The PDB file cannot be opened or parsed.
///
/// # Notes
///
/// - The function uses a loose strictness level when parsing the PDB file.
/// - Interactors are created in pairs, one for each side of a gap between bodies.
/// - The distance restraints are set to exactly match the gap distance (lower and upper margins are set to 0.0).
/// - The output is printed directly to stdout and is not returned or saved to a file.
///
/// # Dependencies
///
/// This function relies on several other modules and functions:
/// - `pdbtbx` for opening and parsing PDB files.
/// - `structure::find_bodies` for identifying separate bodies in the protein structure.
/// - `structure::create_iter_body_gaps` for calculating gaps between bodies.
/// - `Interactor` struct for representing parts of the protein involved in restraints.
/// - `Air` struct for generating the AIR table.
fn restraint_bodies(input_file: &str) -> Result<(), Box<dyn Error>> {
    // Read PDB file
    let pdb = match pdbtbx::open_pdb(input_file, pdbtbx::StrictnessLevel::Loose) {
        Ok((pdb, _warnings)) => pdb,
        Err(e) => {
            panic!("Error opening PDB file: {:?}", e);
        }
    };

    // Find in-contiguous chains
    let bodies = structure::find_bodies(&pdb);
    let gaps = structure::create_iter_body_gaps(&bodies);

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

/// Lists the interface residues for each chain in a protein structure.
///
/// This function analyzes a PDB file to identify the interface residues between chains
/// based on a specified distance cutoff, and prints the results.
///
/// # Arguments
///
/// * `input_file` - A string slice that holds the path to the input PDB file.
/// * `cutoff` - A reference to a f64 value specifying the distance cutoff (in Angstroms)
///              for determining interface residues.
///
/// # Returns
///
/// A `Result<(), Box<dyn Error>>` which is Ok(()) if the function completes successfully,
/// or an Error if something goes wrong.
///
/// # Functionality
///
/// 1. Opens and parses the PDB file using a loose strictness level.
/// 2. Identifies the interface residues using the specified cutoff distance.
/// 3. For each chain involved in the interface:
///    a. Sorts the interface residues.
///    b. Prints the chain ID and the sorted list of interface residue numbers.
///
/// # Panics
///
/// This function will panic if:
/// - The PDB file cannot be opened or parsed.
///
/// # Output Format
///
/// The function prints to stdout in the following format:
/// ```
/// Chain A: [1, 2, 3, 4, 5]
/// Chain B: [10, 11, 12, 13]
/// ```
/// Where each line represents a chain involved in the interface, followed by a sorted list
/// of its interface residue numbers.
///
/// # Notes
///
/// - The function uses a loose strictness level when parsing the PDB file.
/// - Interface residues are determined based on the provided distance cutoff.
/// - The output is printed directly to stdout and is not returned or saved to a file.
///
/// # Dependencies
///
/// This function relies on several other modules and functions:
/// - `pdbtbx` for opening and parsing PDB files.
/// - `structure::get_true_interface` for identifying interface residues.
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
