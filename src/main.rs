mod air;
mod input;
mod interactor;
mod sasa;
mod structure;
use air::Air;
use core::panic;
use interactor::Interactor;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
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

    Z {
        #[arg(help = "Input file")]
        input: String,
        #[arg(long, required = true, help = "Group of residue indexes (can be specified multiple times)", value_parser = parse_residues, number_of_values = 1)]
        residues: Vec<Vec<isize>>,
        #[arg(help = "Spacing between two beads in Angstrom", default_value = "2.0")]
        grid_spacing: f64,
        #[arg(help = "Size in xy dimension", default_value = "10")]
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
        Commands::Z {
            input,
            residues,
            grid_size,
            grid_spacing,
        } => {
            let _ = generate_z_restraints(input, residues, grid_size, grid_spacing);
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

fn generate_z_restraints(
    input_file: &str,
    selections: &[Vec<isize>],
    grid_size: &usize,
    grid_spacing: &f64,
) -> Result<(), Box<dyn Error>> {
    let pdb = match pdbtbx::open_pdb(input_file, pdbtbx::StrictnessLevel::Loose) {
        Ok((pdb, _warnings)) => pdb,
        Err(e) => {
            panic!("Error opening PDB file: {:?}", e);
        }
    };

    // DEVELOPMENT, move the pdb to the origin --------------------------------------------------
    let mut debug_pdb = pdb.clone();
    structure::move_to_origin(&mut debug_pdb);
    let output_path = Path::new("input.pdb");
    let file = File::create(output_path)?;
    pdbtbx::save_pdb_raw(
        &debug_pdb,
        BufWriter::new(file),
        pdbtbx::StrictnessLevel::Strict,
    );
    // -----------------------------------------------------------------------------------------

    let atoms1: Vec<pdbtbx::Atom>;
    let atoms2: Vec<pdbtbx::Atom>;

    let mut restraints: HashMap<isize, Vec<structure::Bead>> = HashMap::new();

    if selections.len() >= 2 {
        (atoms1, atoms2) = structure::find_furthest_selections(selections, &pdb);
    } else {
        atoms1 = structure::get_atoms_from_resnumbers(&pdb, &selections[0]);
        atoms2 = vec![pdbtbx::Atom::new(
            false, // hetero
            1,     // serial_number
            "CA",  // atom_name
            0.0,   // x
            0.0,   // y
            0.0,   // z
            1.0,   // occupancy
            0.0,   // b_factor
            "C",   // element
            0,     // charge
        )
        .expect("Failed to create atom")];
    }

    // let z_axis = structure::calculate_z_axis(&atoms1, &atoms2);
    let center1 = structure::calculate_geometric_center(&atoms1);
    let center2 = structure::calculate_geometric_center(&atoms2);

    // Project endpoints onto global Z-axis and center at origin
    let min_z = center1.z.min(center2.z);
    let max_z = center1.z.max(center2.z);
    // let mid_z = (min_z + max_z) / 2.0;
    let half_length = (max_z - min_z) / 2.0;

    // Project endpoints onto global Z-axis
    // let start_z = center1.z.min(center2.z);
    // let end_z = center1.z.max(center2.z);
    // let start = nalgebra::Vector3::new(0.0, 0.0, start_z);
    // let end = nalgebra::Vector3::new(0.0, 0.0, end_z);

    // // Generate beads along the global Z-axis
    // let num_axis_beads = 10; // Number of beads along the axis
    // let axis_beads = structure::generate_axis_beads(-half_length, half_length, num_axis_beads);

    // Generate grids at both ends, perpendicular to global Z-axis
    let grid_beads1 = structure::generate_grid_beads(-half_length, *grid_size, *grid_spacing);
    let grid_beads2 = structure::generate_grid_beads(half_length, *grid_size, *grid_spacing);

    restraints.insert(0, grid_beads1.clone());
    restraints.insert(1, grid_beads2.clone());

    // let atoms3 = structure::get_atoms_from_resnumbers(&pdb, &selections[2]);
    // let center_atoms3 = structure::calculate_geometric_center(&atoms3);
    // let grid_beads3 = structure::generate_grid_beads(center_atoms3.z, *grid_size, *grid_spacing);

    // Combine all beads
    // let mut all_beads = axis_beads;
    let mut all_beads = Vec::new();
    all_beads.extend(grid_beads1);
    all_beads.extend(grid_beads2);
    // all_beads.extend(grid_beads3);

    structure::write_beads_pdb(&all_beads, "z_beads.pdb")?;

    let mut interactors: Vec<Interactor> = Vec::new();
    let mut counter = 0;
    selections
        .iter()
        .enumerate()
        .for_each(|(index, selection)| {
            let beads = restraints.get(&(index as isize)).unwrap();
            // Get the z coordinates of the first bead
            let z = beads[0].position.z;
            //

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
                interactor_i.set_lower_margin(0.0);
                interactor_i.set_upper_margin(0.0);
                interactor_i.set_target_distance(2.0);

                interactor_j.set_chain("S");
                interactor_j.set_wildcard(format!("attr z gt {:.2}", z).as_str());

                interactors.push(interactor_i);
                interactors.push(interactor_j);
            });
        });

    // interactors.iter().for_each(|interactor| {
    //     println!("Interactor: {:?}", interactor);
    // });

    let air = Air::new(interactors);
    let tbl = air.gen_tbl().unwrap();

    println!("{}", tbl);

    Ok(())
}

// sele s1, resid 19+83+145+167 and input
// sele s2, resid 98+101+126+129 and input
// sele s3, resid 23+62+87+111+116+153+163 and input
