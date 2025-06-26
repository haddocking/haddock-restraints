use std::error::Error;

use crate::*;

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
pub fn restraint_bodies(input_file: &str, pml: &Option<String>) -> Result<String, Box<dyn Error>> {
    // Read PDB file
    let pdb = match load_pdb(input_file) {
        Ok(pdb) => pdb,
        Err(e) => {
            panic!("Error opening PDB file: {:?}", e);
        }
    };

    // Find in-contiguous chains
    let bodies = find_bodies(&pdb);
    let mut gaps = create_iter_body_gaps(&bodies);

    // Find same-chain ligands
    let ligand_gaps = gaps_around_ligand(&pdb);

    // NOTE: One restraint per atom of the ligand might be too much, apply some filtering
    // - no duplicated atoms in the protein should be used
    // - select every-other pair
    let filtered_gaps = filter_unique_by_atom_j(ligand_gaps);
    let every_other_gaps: Vec<Gap> = filtered_gaps
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

#[cfg(test)]
mod tests {

    use super::*;
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
}
