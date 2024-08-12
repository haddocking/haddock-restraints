use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use kd_tree::KdTree;
use nalgebra::Vector3;
use pdbtbx::Residue;
use pdbtbx::{Atom, PDB};
use rand::rngs::StdRng;
use rand::seq::SliceRandom;
use rand::SeedableRng;
use std::fs::File;
use std::io::Write;

pub struct Bead {
    position: Vector3<f64>,
}

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

#[derive(Debug)]
pub struct Gap {
    pub chain: String,
    pub atom_i: String,
    pub atom_j: String,
    pub res_i: isize,
    pub res_j: isize,
    pub distance: f64,
}

pub fn find_bodies(pdb: &pdbtbx::PDB) -> HashMap<isize, Vec<(isize, &str, &pdbtbx::Atom)>> {
    // Check if the distance of a given atom to its next one is higher than 4A

    // Get only the `CA` atoms
    let mut ca_atoms: Vec<(&str, isize, &pdbtbx::Atom)> = Vec::new();
    pdb.chains().for_each(|chain| {
        chain.residues().for_each(|residue| {
            residue.atoms().for_each(|atom| {
                if atom.name() == "CA" {
                    ca_atoms.push((chain.id(), residue.serial_number(), atom));
                }
            });
        });
    });
    let mut bodies: HashMap<isize, Vec<(isize, &str, &pdbtbx::Atom)>> = HashMap::new();
    let mut body_id = 0;
    for (i, j) in ca_atoms.iter().zip(ca_atoms.iter().skip(1)) {
        let (chain_i, res_i, atom_i) = i;
        let (chain_j, _res_j, atom_j) = j;
        if chain_i != chain_j {
            continue;
        }
        let distance = atom_i.distance(atom_j);
        if distance > 4.0 {
            body_id += 1;
        }
        bodies
            .entry(body_id)
            .or_default()
            .push((*res_i, chain_i, atom_i));
    }

    bodies
}

pub fn create_iter_body_gaps(
    bodies: &HashMap<isize, Vec<(isize, &str, &pdbtbx::Atom)>>,
) -> Vec<Gap> {
    let mut rng = StdRng::seed_from_u64(42);
    let mut pairs: Vec<Gap> = Vec::new();
    let body_ids: Vec<isize> = bodies.keys().cloned().collect();
    for (i, &body_id1) in body_ids.iter().enumerate() {
        for &body_id2 in body_ids[i + 1..].iter() {
            if let (Some(atoms1), Some(atoms2)) = (bodies.get(&body_id1), bodies.get(&body_id2)) {
                if atoms1.len() >= 2 && atoms2.len() >= 2 {
                    let selected1: Vec<_> = atoms1.choose_multiple(&mut rng, 2).cloned().collect();
                    let selected2: Vec<_> = atoms2.choose_multiple(&mut rng, 2).cloned().collect();

                    for i in 0..2 {
                        pairs.push(Gap {
                            chain: selected1[i].1.to_string(),
                            atom_i: selected1[i].2.name().to_string(),
                            atom_j: selected2[i].2.name().to_string(),
                            res_i: selected1[i].0,
                            res_j: selected2[i].0,
                            distance: selected1[i].2.distance(selected2[i].2),
                        });
                    }
                }
            }
        }
    }

    pairs
}

// pub fn calculate_z_axis(selection1: &[pdbtbx::Atom], selection2: &[pdbtbx::Atom]) -> Vector3<f64> {
//     let center1 = calculate_geometric_center(selection1);
//     let center2 = calculate_geometric_center(selection2);
//     (center2 - center1).normalize()
// }

pub fn calculate_geometric_center(atoms: &[pdbtbx::Atom]) -> Vector3<f64> {
    let mut center = Vector3::zeros();
    for atom in atoms {
        center += Vector3::new(atom.x(), atom.y(), atom.z());
    }
    center / (atoms.len() as f64)
}

// pub fn generate_axis_beads(start_z: f64, end_z: f64, num_beads: usize) -> Vec<Bead> {
//     let mut beads = Vec::with_capacity(num_beads);
//     let length = end_z - start_z;

//     for i in 0..num_beads {
//         let t = (i as f64) / ((num_beads - 1) as f64);
//         let z = start_z + t * length;
//         let position = Vector3::new(0.0, 0.0, z);
//         beads.push(Bead { position });
//     }
//     beads
// }

pub fn generate_grid_beads(center_z: f64, grid_size: usize, spacing: f64) -> Vec<Bead> {
    let mut beads = Vec::with_capacity(grid_size * grid_size);

    let grid_offset = (grid_size as f64 - 1.0) * spacing / 2.0;

    for i in 0..grid_size {
        for j in 0..grid_size {
            let x = (i as f64) * spacing - grid_offset;
            let y = (j as f64) * spacing - grid_offset;

            let final_pos = Vector3::new(x, y, center_z);
            beads.push(Bead {
                position: final_pos,
            });
        }
    }

    beads
}

pub fn write_beads_pdb(beads: &[Bead], filename: &str) -> std::io::Result<()> {
    let mut file = File::create(filename)?;

    for (i, bead) in beads.iter().enumerate() {
        writeln!(
            file,
            "{:<6}{:>5} {:<4}{:1}{:<3} {:1}{:>4}{:1}   {:>8.3}{:>8.3}{:>8.3}{:>6.2}{:>6.2}          {:>2}{:>2}",
            "ATOM",
            i + 1,
            "SHA",
            " ",
            "SHA",
            "S",
            i + 1,
            " ",
            bead.position.x,
            bead.position.y,
            bead.position.z,
            1.00,
            1.00,
            "H",
            ""
        )?;
    }
    writeln!(file, "END")?;

    Ok(())
}

pub fn find_furthest_selections(selections: &[Vec<isize>], pdb: &PDB) -> (Vec<Atom>, Vec<Atom>) {
    let sele: HashMap<usize, (Vec<Atom>, Vector3<f64>)> = selections
        .iter()
        .enumerate()
        .map(|(i, sel)| {
            let atoms = get_atoms_from_resnumbers(pdb, sel);
            let center = calculate_geometric_center(&atoms);
            (i, (atoms, center))
        })
        .collect();

    // Find the two selections that are the furthest apart
    let mut max_distance = 0.0;
    let mut atoms1 = Vec::new();
    let mut atoms2 = Vec::new();
    for (i, j) in (0..selections.len()).tuple_combinations() {
        let distance = (sele[&i].1 - sele[&j].1).norm();
        if distance > max_distance {
            max_distance = distance;
            atoms1 = sele[&i].0.clone();
            atoms2 = sele[&j].0.clone();
        }
    }

    (atoms1, atoms2)
}

pub fn get_atoms_from_resnumbers(pdb: &PDB, selection: &[isize]) -> Vec<Atom> {
    pdb.residues()
        .filter(|res| selection.contains(&res.serial_number()))
        .flat_map(|res| res.atoms())
        .cloned() // This clones each &Atom into an owned Atom
        .collect()
}

pub fn move_to_origin(pdb: &mut PDB) {
    let center = calculate_geometric_center(&pdb.atoms().cloned().collect::<Vec<_>>());
    for atom in pdb.atoms_mut() {
        let x = atom.x() - center.x;
        let y = atom.y() - center.y;
        let z = atom.z() - center.z;
        let _ = atom.set_pos((x, y, z));
    }
}

#[cfg(test)]
mod tests {

    use std::env;

    use super::*;

    #[test]
    fn test_get_true_interface() {
        let pdb_path = env::current_dir().unwrap().join("tests/data/complex.pdb");

        let pdb = pdbtbx::open_pdb(pdb_path.to_str().unwrap(), pdbtbx::StrictnessLevel::Loose)
            .unwrap()
            .0;

        let observed_true_interface = get_true_interface(&pdb, 5.0);

        let mut expected_true_interface: HashMap<String, HashSet<isize>> = HashMap::new();
        let chain_a = "A".to_string();
        let res_a: HashSet<isize> = vec![934, 936, 933, 946, 950, 938, 940, 941, 937, 931]
            .into_iter()
            .collect();
        let chain_b = "B";
        let res_b: HashSet<isize> = vec![66, 48, 68, 49, 46, 45, 44, 47, 69, 6, 70, 8, 42]
            .into_iter()
            .collect();

        expected_true_interface.insert(chain_a, res_a);
        expected_true_interface.insert(chain_b.to_string(), res_b);

        assert_eq!(observed_true_interface, expected_true_interface);
    }
}
