use std::collections::{HashMap, HashSet};

use kd_tree::KdTree;
use pdbtbx::Residue;

pub fn neighbor_search(
    pdb: pdbtbx::PDB,
    target_residues_numbers: Vec<&pdbtbx::Residue>,
    radius: f64,
) -> Vec<isize> {
    let points: Vec<[f64; 3]> = pdb
        .atoms()
        .map(|atom| {
            let x = atom.x();
            let y = atom.y();
            let z = atom.z();
            [x, y, z]
        })
        .collect();

    let kdtree = KdTree::build_by_ordered_float(points.clone());

    // Find the coordinates of each target residue
    // let mut neighbors: Vec<&Atom> = Vec::new();
    let mut result: HashSet<isize> = HashSet::new();
    for residue in &target_residues_numbers {
        let query_vec: Vec<[f64; 3]> = residue
            .atoms()
            .map(|atom| {
                let x = atom.x();
                let y = atom.y();
                let z = atom.z();
                [x, y, z]
            })
            .collect();

        for query in query_vec {
            let found_vec = kdtree.within_radius(&query, radius);
            for found in found_vec {
                // Loop over the whole PDB and find an atom that matches the found point
                // FIXME: Optimize this
                for residue in pdb.residues() {
                    for atom in residue.atoms() {
                        let x = atom.x();
                        let y = atom.y();
                        let z = atom.z();
                        if x == found[0] && y == found[1] && z == found[2] {
                            result.insert(residue.serial_number());
                        }
                    }
                }
            }
        }
    }
    // Remove the target residues from the result
    for residue in target_residues_numbers {
        result.remove(&residue.serial_number());
    }

    // Sort the result
    let mut result: Vec<isize> = result.into_iter().collect();
    result.sort();

    result
}

pub fn get_residues(pdb: &pdbtbx::PDB, target: Vec<isize>) -> Vec<&Residue> {
    let mut result: Vec<&Residue> = Vec::new();
    for residue in pdb.residues() {
        if target.contains(&residue.serial_number()) {
            result.push(residue);
        }
    }
    result
}

pub fn get_true_interface(pdb: &pdbtbx::PDB, cutoff: f64) -> HashMap<String, HashSet<isize>> {
    let mut true_interface: HashMap<String, HashSet<isize>> = HashMap::new();
    for chain_i in pdb.chains() {
        for chain_j in pdb.chains() {
            if chain_i.id() == chain_j.id() {
                continue;
            }
            for res_i in chain_i.residues() {
                for res_j in chain_j.residues() {
                    // Check if any of the atoms are touching each other
                    for atom_i in res_i.atoms() {
                        for atom_j in res_j.atoms() {
                            let distance = atom_i.distance(atom_j);
                            if distance < cutoff {
                                true_interface
                                    .entry(chain_i.id().to_string())
                                    .or_default()
                                    .insert(res_i.serial_number());
                                true_interface
                                    .entry(chain_j.id().to_string())
                                    .or_default()
                                    .insert(res_j.serial_number());
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    true_interface
}

pub fn get_chains_in_contact(pdb: &pdbtbx::PDB, cutoff: f64) -> HashSet<(String, String)> {
    let mut chains_in_contact: HashSet<(String, String)> = HashSet::new();

    for chain_i in pdb.chains() {
        for chain_j in pdb.chains() {
            if chain_i.id() == chain_j.id() {
                continue;
            }
            for res_i in chain_i.residues() {
                for res_j in chain_j.residues() {
                    // Check if any of the atoms are touching each other
                    for atom_i in res_i.atoms() {
                        for atom_j in res_j.atoms() {
                            let distance = atom_i.distance(atom_j);
                            if distance < cutoff {
                                chains_in_contact
                                    .insert((chain_i.id().to_string(), chain_j.id().to_string()));
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    chains_in_contact
}
