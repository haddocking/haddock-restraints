mod air;
mod input;
mod interactor;
mod sasa;
mod structure;
mod utils;
use air::Air;
use core::panic;
use interactor::Interactor;
use std::collections::HashMap;
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
            gen_tbl(input, pml);
        }
        Commands::Ti { input, cutoff, pml } => {
            let _ = true_interface(input, cutoff, pml);
        }
        Commands::UnambigTi { input, cutoff, pml } => {
            let _ = unambig_ti(input, cutoff, pml);
        }
        Commands::Restraint { input, pml } => {
            let _ = restraint_bodies(input, pml);
        }
        Commands::Interface { input, cutoff } => {
            let _ = list_interface(input, cutoff);
        }
        Commands::Z {
            input,
            output,
            residues,
            grid_size,
            grid_spacing,
        } => {
            let _ = generate_z_restraints(input, output, residues, grid_size, grid_spacing);
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
fn gen_tbl(input_file: &str, pml: &Option<String>) {
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

    if let Some(output_f) = pml {
        air.gen_pml(output_f)
    };
}

/// Generates Unambiguous Topological Interactions (TIs) from a protein structure.
///
/// This function reads a PDB file, identifies the closest residue pairs based on a specified distance cutoff,
/// and creates unambiguous interactors for each residue pair.
///
/// # Arguments
///
/// * input_file - A string slice that holds the path to the input PDB file.
/// * cutoff - A reference to a f64 value specifying the distance cutoff (in Angstroms) for determining interactions.
///
/// # Returns
///
/// A Result<String, Box<dyn Error>> containing the generated TBL (Topological Restraints List) if successful.
///
/// # Functionality
///
/// 1. Loads the PDB file using the structure::load_pdb function.
/// 2. Finds the closest residue pairs within the specified distance cutoff.
/// 3. Creates Interactor instances for each residue pair.
/// 4. Assigns chains, active/passive residues, and atoms to the interactors.
/// 5. Generates the Topological Restraints List (TBL) using the created interactors.
/// 6. Prints the generated TBL to stdout.
///
/// # Panics
///
/// This function will panic if:
/// - The PDB file cannot be opened or parsed.
fn unambig_ti(
    input_file: &str,
    cutoff: &f64,
    pml: &Option<String>,
) -> Result<String, Box<dyn Error>> {
    let pdb = match structure::load_pdb(input_file) {
        Ok(pdb) => pdb,
        Err(e) => {
            panic!("Error opening PDB file: {:?}", e);
        }
    };
    let pairs = structure::get_closest_residue_pairs(&pdb, *cutoff);

    let mut interactors: Vec<Interactor> = Vec::new();
    let mut counter = 0;
    pairs.iter().for_each(|g| {
        let mut interactor_i = Interactor::new(counter);
        counter += 1;
        let mut interactor_j = Interactor::new(counter);
        interactor_j.add_target(counter - 1);
        interactor_i.add_target(counter);
        counter += 1;

        interactor_i.set_chain(g.chain_i.as_str());
        interactor_i.set_active(vec![g.res_i as i16]);
        interactor_i.set_active_atoms(vec![g.atom_i.clone()]);
        interactor_i.set_target_distance(g.distance);

        interactor_j.set_chain(g.chain_j.as_str());
        interactor_j.set_passive(vec![g.res_j as i16]);
        interactor_j.set_passive_atoms(vec![g.atom_j.clone()]);

        interactors.push(interactor_i);
        interactors.push(interactor_j);
    });

    // Make the restraints
    let air = Air::new(interactors);
    let tbl = air.gen_tbl().unwrap();

    println!("{}", tbl);

    if let Some(output_f) = pml {
        air.gen_pml(output_f)
    };

    Ok(tbl)
}

/// Analyzes the true interface of a protein structure and generates Ambiguous Interaction Restraints (AIRs).
///
/// This function reads a PDB file, identifies the true interface between chains based on a distance cutoff,
/// creates interactors for each chain involved in the interface, and generates AIRs.
///
/// # Arguments
///
/// * `input_file` - A string slice that holds the path to the input PDB file.
/// * `cutoff` - A reference to a f64 value specifying the distance cutoff (in Angstroms) for determining interfaces.
///
/// # Returns
///
/// A `Result<String, Box<dyn Error>>` which is Ok(String) if the function completes successfully, or an Error if something goes wrong.
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
fn true_interface(
    input_file: &str,
    cutoff: &f64,
    pml: &Option<String>,
) -> Result<String, Box<dyn Error>> {
    // Read PDB file
    let pdb = match structure::load_pdb(input_file) {
        Ok(pdb) => pdb,
        Err(e) => {
            panic!("Error opening PDB file: {:?}", e);
        }
    };

    let true_interface = structure::get_true_interface(&pdb, *cutoff);
    let chains_in_contact = structure::get_chains_in_contact(&pdb, *cutoff);

    // Sort the true_interface by chain id
    let mut true_interface: Vec<_> = true_interface.iter().collect();
    true_interface.sort_by(|a, b| a.0.cmp(b.0));

    // NOTE: Here the IDs of the interactors are their position in the PDB file; this is so that
    // we can handle any order of chains.
    let mut interactors: Vec<Interactor> = Vec::new();
    for (chain_id, residues) in true_interface.iter() {
        // Get what is the position of this chain in the PDB, this will be its ID
        let target_id = pdb
            .chains()
            .position(|chain| chain.id() == *chain_id)
            .unwrap();

        let mut interactor = Interactor::new(target_id as u16);
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

    if let Some(output_f) = pml {
        air.gen_pml(output_f)
    };

    Ok(tbl)
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
fn restraint_bodies(input_file: &str, pml: &Option<String>) -> Result<String, Box<dyn Error>> {
    // Read PDB file
    let pdb = match structure::load_pdb(input_file) {
        Ok(pdb) => pdb,
        Err(e) => {
            panic!("Error opening PDB file: {:?}", e);
        }
    };

    // Find in-contiguous chains
    let bodies = structure::find_bodies(&pdb);
    let mut gaps = structure::create_iter_body_gaps(&bodies);

    // Find same-chain ligands
    let ligand_gaps = structure::gaps_around_ligand(&pdb);

    // NOTE: One restraint per atom of the ligand might be too much, apply some filtering
    // - no duplicated atoms in the protein should be used
    // - select every-other pair
    let filtered_gaps = structure::filter_unique_by_atom_j(ligand_gaps);
    let every_other_gaps: Vec<structure::Gap> = filtered_gaps
        .into_iter()
        .enumerate()
        .filter_map(|(idx, g)| {
            if idx % 2 == 0 {
                // select even
                Some(g)
            } else {
                None
            }
        })
        .collect();

    gaps.extend(every_other_gaps);

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
        interactor_i.set_target_distance(g.distance);
        interactor_i.set_lower_margin(0.0);
        interactor_i.set_upper_margin(0.0);

        interactor_j.set_chain(g.chain.as_str());
        interactor_j.set_passive(vec![g.res_j as i16]);
        interactor_j.set_passive_atoms(vec![g.atom_j.clone()]);

        interactors.push(interactor_i);
        interactors.push(interactor_j);
    });

    let air = Air::new(interactors);
    let tbl = air.gen_tbl().unwrap();

    println!("{}", tbl);

    if let Some(output_f) = pml {
        air.gen_pml(output_f)
    };
    Ok(tbl)
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
    let pdb = match structure::load_pdb(input_file) {
        Ok(pdb) => pdb,
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

fn generate_z_restraints(
    input_file: &str,
    output_file: &str,
    selections: &[Vec<isize>],
    grid_size: &usize,
    grid_spacing: &f64,
) -> Result<(), Box<dyn Error>> {
    let pdb = match structure::load_pdb(input_file) {
        Ok(pdb) => pdb,
        Err(e) => {
            panic!("Error opening PDB file: {:?}", e);
        }
    };

    // // DEVELOPMENT, move the pdb to the origin --------------------------------------------------
    // let mut debug_pdb = pdb.clone();
    // structure::move_to_origin(&mut debug_pdb);
    // let output_path = Path::new("input.pdb");
    // let file = File::create(output_path)?;
    // pdbtbx::save_pdb_raw(
    //     &debug_pdb,
    //     BufWriter::new(file),
    //     pdbtbx::StrictnessLevel::Strict,
    // );
    // // -----------------------------------------------------------------------------------------

    let atoms1: Vec<pdbtbx::Atom>;
    let atoms2: Vec<pdbtbx::Atom>;

    let mut restraints: HashMap<usize, Vec<structure::Bead>> = HashMap::new();

    if selections.len() >= 2 {
        (atoms1, atoms2) = structure::find_furthest_selections(selections, &pdb);
    } else {
        atoms1 = structure::get_atoms_from_resnumbers(&pdb, &selections[0]);
        atoms2 = vec![
            pdbtbx::Atom::new(false, 1, "CA", 0.0, 0.0, 0.0, 1.0, 0.0, "C", 0)
                .expect("Failed to create atom"),
        ];
    }

    let center1 = structure::calculate_geometric_center(&atoms1);
    let center2 = structure::calculate_geometric_center(&atoms2);

    // Project endpoints onto global Z-axis and center at origin
    let min_z = center1.z.min(center2.z);
    let max_z = center1.z.max(center2.z);
    let half_length = (max_z - min_z) / 2.0;

    // Generate grids at both ends, perpendicular to global Z-axis
    let grid_beads1 = structure::generate_grid_beads(-half_length, *grid_size, *grid_spacing);
    let grid_beads2 = structure::generate_grid_beads(half_length, *grid_size, *grid_spacing);

    restraints.insert(0, grid_beads1.clone());
    restraints.insert(1, grid_beads2.clone());

    let mut all_beads = Vec::new();
    all_beads.extend(grid_beads1);
    all_beads.extend(grid_beads2);

    // It can be that `selections` contains more than 2 selections, if that's the case, we need to place more grids in between
    if selections.len() > 2 {
        for (i, selection) in selections.iter().enumerate().skip(2) {
            let atoms = structure::get_atoms_from_resnumbers(&pdb, selection);
            let center = structure::calculate_geometric_center(&atoms);
            let grid_beads = structure::generate_grid_beads(center.z, *grid_size, *grid_spacing);
            restraints.insert(i, grid_beads.clone());
            all_beads.extend(grid_beads);
        }
    }

    // Write the beads to a PDB file
    structure::write_beads_pdb(&all_beads, output_file)?;

    let mut interactors: Vec<Interactor> = Vec::new();
    let mut counter = 0;
    let restraint_distance = ((grid_spacing / 2.0) - 2.0).max(2.0);
    selections
        .iter()
        .enumerate()
        .for_each(|(index, selection)| {
            let beads = restraints.get(&(index)).unwrap();
            let z = beads[0].position.z;

            let comparison_operator = if z >= 0.0 { "ge" } else { "le" };

            selection.iter().for_each(|resnum| {
                let mut interactor_i = Interactor::new(counter);
                counter += 1;
                let mut interactor_j = Interactor::new(counter);
                interactor_j.add_target(counter - 1);
                interactor_i.add_target(counter);
                counter += 1;

                interactor_i.set_chain("A");
                interactor_i.set_active(vec![*resnum as i16]);
                interactor_i.set_active_atoms(vec!["CA".to_string()]);
                interactor_i.set_passive_atoms(vec!["SHA".to_string()]);
                interactor_i.set_target_distance(restraint_distance);
                interactor_i.set_lower_margin(restraint_distance);
                interactor_i.set_upper_margin(0.0);

                interactor_j.set_chain("S");
                interactor_j
                    .set_wildcard(format!("and attr z {} {:.3}", comparison_operator, z).as_str());

                interactors.push(interactor_i);
                interactors.push(interactor_j);
            });
        });

    let air = Air::new(interactors);
    let tbl = air.gen_tbl().unwrap();

    println!("{}", tbl);

    Ok(())
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_true_interface() {
        let expected_tbl = r#"assign ( resid 933 and segid A )
       (
        ( resid 46 and segid B )
     or
        ( resid 47 and segid B )
       ) 2.0 2.0 0.0

assign ( resid 950 and segid A )
       (
        ( resid 46 and segid B )
     or
        ( resid 47 and segid B )
       ) 2.0 2.0 0.0

assign ( resid 46 and segid B )
       (
        ( resid 933 and segid A )
     or
        ( resid 950 and segid A )
       ) 2.0 2.0 0.0

assign ( resid 47 and segid B )
       (
        ( resid 933 and segid A )
     or
        ( resid 950 and segid A )
       ) 2.0 2.0 0.0

"#;

        let opt: Option<String> = None;
        match true_interface("tests/data/complex.pdb", &3.0, &opt) {
            Ok(tbl) => assert_eq!(tbl, expected_tbl),
            Err(_e) => (),
        };
    }
    #[test]
    fn test_true_interface_ba() {
        let expected_tbl = r#"assign ( resid 933 and segid A )
       (
        ( resid 46 and segid B )
     or
        ( resid 47 and segid B )
       ) 2.0 2.0 0.0

assign ( resid 950 and segid A )
       (
        ( resid 46 and segid B )
     or
        ( resid 47 and segid B )
       ) 2.0 2.0 0.0

assign ( resid 46 and segid B )
       (
        ( resid 933 and segid A )
     or
        ( resid 950 and segid A )
       ) 2.0 2.0 0.0

assign ( resid 47 and segid B )
       (
        ( resid 933 and segid A )
     or
        ( resid 950 and segid A )
       ) 2.0 2.0 0.0

"#;

        let opt: Option<String> = None;
        match true_interface("tests/data/complex_BA.pdb", &3.0, &opt) {
            Ok(tbl) => assert_eq!(tbl, expected_tbl),
            Err(_e) => (),
        };
    }

    #[test]
    fn test_restraint_bodies() {
        let expected_tbl = r"assign ( resid 1 and segid A and name CA ) ( resid 4 and segid A and name CA ) 10.2 0.0 0.0

assign ( resid 2 and segid A and name CA ) ( resid 8 and segid A and name CA ) 14.2 0.0 0.0

";

        let opt: Option<String> = None;
        match restraint_bodies("tests/data/gaps.pdb", &opt) {
            Ok(tbl) => assert_eq!(tbl, expected_tbl),
            Err(_e) => (),
        }
    }
    #[test]
    fn test_unambigti() {
        let opt: Option<String> = None;
        let expected_tbl = "assign ( resid 2 and segid A and name CA ) ( resid 10 and segid B and name CA ) 9.1 2.0 0.0\n\n";

        match unambig_ti("tests/data/two_res.pdb", &5.0, &opt) {
            Ok(tbl) => assert_eq!(tbl, expected_tbl),
            Err(_e) => (),
        }
    }
}
