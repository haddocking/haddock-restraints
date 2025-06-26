use crate::*;
use std::error::Error;

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
pub fn true_interface(
    input_file: &str,
    cutoff: &f64,
    pml: &Option<String>,
) -> Result<String, Box<dyn Error>> {
    // Read PDB file
    let pdb = match load_pdb(input_file) {
        Ok(pdb) => pdb,
        Err(e) => {
            panic!("Error opening PDB file: {:?}", e);
        }
    };

    let true_interface = get_true_interface(&pdb, *cutoff);
    let chains_in_contact = get_chains_in_contact(&pdb, *cutoff);

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
pub fn unambig_ti(
    input_file: &str,
    cutoff: &f64,
    pml: &Option<String>,
) -> Result<String, Box<dyn Error>> {
    let pdb = match load_pdb(input_file) {
        Ok(pdb) => pdb,
        Err(e) => {
            panic!("Error opening PDB file: {:?}", e);
        }
    };
    let pairs = get_closest_residue_pairs(&pdb, *cutoff);

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

/// Lists the interface residues for each chain in a protein structure.
///
/// This function analyzes a PDB file to identify the interface residues between chains
/// based on a specified distance cutoff, and prints the results.
///
/// # Arguments
///
/// * `input_file` - A string slice that holds the path to the input PDB file.
/// * `cutoff` - A reference to a f64 value specifying the distance cutoff (in Angstroms)
///   for determining interface residues.
///
/// # Returns
///
/// A `Result<(), Box<dyn Error>>` which is Ok(()) if the function completes successfully,
/// or an Error if something goes wrong.
///
pub fn list_interface(input_file: &str, cutoff: &f64) -> Result<(), Box<dyn Error>> {
    let pdb = match load_pdb(input_file) {
        Ok(pdb) => pdb,
        Err(e) => {
            panic!("Error opening PDB file: {:?}", e);
        }
    };

    let true_interface = get_true_interface(&pdb, *cutoff);

    for (chain_id, residues) in true_interface.iter() {
        let mut sorted_res = residues.iter().collect::<Vec<_>>();
        sorted_res.sort();

        println!("Chain {}: {:?}", chain_id, sorted_res);
    }

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
    fn test_unambigti() {
        let opt: Option<String> = None;
        let expected_tbl = "assign ( resid 2 and segid A and name CA ) ( resid 10 and segid B and name CA ) 9.1 2.0 0.0\n\n";

        match unambig_ti("tests/data/two_res.pdb", &5.0, &opt) {
            Ok(tbl) => assert_eq!(tbl, expected_tbl),
            Err(_e) => (),
        }
    }
}
