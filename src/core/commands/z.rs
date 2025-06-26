use crate::*;
use std::{collections::HashMap, error::Error};

// TODO: Docstring
pub fn generate_z_restraints(
    pdb: pdbtbx::PDB,
    output_file: &str,
    selections: &[Vec<isize>],
    grid_size: &usize,
    grid_spacing: &f64,
) -> Result<(), Box<dyn Error>> {
    // // DEVELOPMENT, move the pdb to the origin --------------------------------------------------
    // let mut debug_pdb = pdb.clone();
    // move_to_origin(&mut debug_pdb);
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

    let mut restraints: HashMap<usize, Vec<Bead>> = HashMap::new();

    if selections.len() >= 2 {
        (atoms1, atoms2) = find_furthest_selections(selections, &pdb);
    } else {
        atoms1 = get_atoms_from_resnumbers(&pdb, &selections[0]);
        atoms2 = vec![
            pdbtbx::Atom::new(false, 1, "CA", 0.0, 0.0, 0.0, 1.0, 0.0, "C", 0)
                .expect("Failed to create atom"),
        ];
    }

    let center1 = calculate_geometric_center(&atoms1);
    let center2 = calculate_geometric_center(&atoms2);

    // Project endpoints onto global Z-axis and center at origin
    let min_z = center1.z.min(center2.z);
    let max_z = center1.z.max(center2.z);
    let half_length = (max_z - min_z) / 2.0;

    // Generate grids at both ends, perpendicular to global Z-axis
    let grid_beads1 = generate_grid_beads(-half_length, *grid_size, *grid_spacing);
    let grid_beads2 = generate_grid_beads(half_length, *grid_size, *grid_spacing);

    restraints.insert(0, grid_beads1.clone());
    restraints.insert(1, grid_beads2.clone());

    let mut all_beads = Vec::new();
    all_beads.extend(grid_beads1);
    all_beads.extend(grid_beads2);

    // It can be that `selections` contains more than 2 selections, if that's the case, we need to place more grids in between
    if selections.len() > 2 {
        for (i, selection) in selections.iter().enumerate().skip(2) {
            let atoms = get_atoms_from_resnumbers(&pdb, selection);
            let center = calculate_geometric_center(&atoms);
            let grid_beads = generate_grid_beads(center.z, *grid_size, *grid_spacing);
            restraints.insert(i, grid_beads.clone());
            all_beads.extend(grid_beads);
        }
    }

    // Write the beads to a PDB file
    write_beads_pdb(&all_beads, output_file)?;

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
mod tests {}
