use core::f64;
use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use kd_tree::KdTree;
use nalgebra::Vector3;
use pdbtbx::{Atom, Element, PDB, ReadOptions};
use pdbtbx::{PDBError, Residue};
use rand::SeedableRng;
use rand::prelude::IndexedRandom;
use rand::rngs::StdRng;
use std::fs::File;
use std::io::BufReader;
use std::io::Cursor;
use std::io::Write;

#[derive(Clone)]
pub struct Bead {
    pub position: Vector3<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Pair {
    pub chain_i: String,
    pub res_i: isize,
    pub atom_i: String,
    pub chain_j: String,
    pub res_j: isize,
    pub atom_j: String,
    pub distance: f64,
}

/// The `Finder` struct contains two hashmaps to quickly retrieve the chain ID and residue number associated with a given atom by its serial number. These lookups are initialized from the provided PDB file during instantiation.
///
struct Finder {
    chain_lookup: HashMap<usize, String>,
    residue_lookup: HashMap<usize, isize>,
}

impl Finder {
    /// Creates a new `Finder` instance.
    ///
    /// This method initializes the `Finder` struct by extracting the chain and residue information from the provided `PDB` reference. It constructs two hashmaps: `chain_lookup` (for mapping atom serial numbers to chain IDs) and `residue_lookup` (for mapping atom serial numbers to residue serial numbers).
    ///
    fn new(pdb: &PDB) -> Self {
        let chain_lookup = pdb
            .chains()
            .flat_map(|chain| {
                let chain_id = chain.id().to_string();
                chain.residues().flat_map(move |residue| {
                    residue.atoms().map({
                        let v = chain_id.clone();
                        move |atom| (atom.serial_number(), v.clone())
                    })
                })
            })
            .collect();

        let residue_lookup = pdb
            .chains()
            .flat_map(|chain| {
                chain.residues().flat_map(|residue| {
                    let res_num = residue.serial_number();
                    residue
                        .atoms()
                        .map(move |atom| (atom.serial_number(), res_num))
                })
            })
            .collect();

        Finder {
            chain_lookup,
            residue_lookup,
        }
    }

    /// Finds the chain ID associated with a given atom.
    ///
    /// This method looks up the chain ID for an atom based on its serial number. If no chain is found, the default value "A" is returned.
    ///
    fn find_chain_id(&self, atom: &Atom) -> String {
        self.chain_lookup
            .get(&atom.serial_number())
            .unwrap_or(&"A".to_string())
            .clone()
    }

    /// Finds the residue serial number associated with a given atom.
    ///
    /// This method looks up the residue serial number for an atom based on its serial number. If no residue number is found, the default value `0` is returned.
    ///
    fn find_residue_number(&self, atom: &Atom) -> isize {
        *self.residue_lookup.get(&atom.serial_number()).unwrap_or(&0)
    }
}

/// Performs a neighbor search for residues within a specified radius of target residues.
///
/// This function uses a K-d tree to efficiently find residues that are within a given radius
/// of any atom in the target residues.
///
/// # Arguments
///
/// * `pdb` - A `pdbtbx::PDB` structure representing the entire protein.
/// * `target_residues_numbers` - A vector of references to `pdbtbx::Residue` objects representing the target residues.
/// * `radius` - A `f64` value specifying the search radius in Angstroms.
///
/// # Returns
///
/// A `Vec<isize>` containing the sorted serial numbers of residues that are within the specified
/// radius of any atom in the target residues, excluding the target residues themselves.
///
/// # Algorithm
///
/// 1. Constructs a K-d tree from all atom coordinates in the PDB structure.
/// 2. For each atom in the target residues:
///    a. Performs a radius search in the K-d tree.
///    b. Identifies the residues corresponding to the atoms found in the search.
/// 3. Removes the target residues from the result set.
/// 4. Sorts the resulting residue serial numbers.
///
/// # Notes
///
/// - The function uses the `KdTree` data structure for efficient spatial searching.
/// - The current implementation may have performance bottlenecks in the atom-to-residue mapping step.
/// - The function assumes that atom coordinates are unique for identification purposes.
///
/// # TODO
///
/// - Optimize the atom-to-residue mapping step for better performance.
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

/// Retrieves specific residues from a PDB structure based on their serial numbers.
///
/// # Arguments
///
/// * `pdb` - A reference to a `pdbtbx::PDB` structure representing the entire protein.
/// * `target` - A vector of `isize` values representing the serial numbers of the residues to retrieve.
///
/// # Returns
///
/// A `Vec<&Residue>` containing references to the residues whose serial numbers match those in the `target` vector.
///
/// # Notes
///
/// - This function performs a linear search through all residues in the PDB structure.
/// - The order of residues in the result vector matches the order they appear in the PDB structure,
///   not necessarily the order of serial numbers in the `target` vector.
/// - If a serial number in `target` doesn't match any residue in the PDB, it is simply ignored.
pub fn get_residues(pdb: &pdbtbx::PDB, target: Vec<isize>) -> Vec<&Residue> {
    let mut result: Vec<&Residue> = Vec::new();
    for residue in pdb.residues() {
        if target.contains(&residue.serial_number()) {
            result.push(residue);
        }
    }
    result
}

/// Identifies the true interface residues between chains in a PDB structure.
///
/// This function determines interface residues by finding pairs of residues from different chains
/// that have atoms within a specified cutoff distance of each other.
///
/// # Arguments
///
/// * `pdb` - A reference to a `pdbtbx::PDB` structure representing the entire protein.
/// * `cutoff` - A `f64` value specifying the maximum distance (in Angstroms) between atoms
///   for residues to be considered part of the interface.
///
/// # Returns
///
/// A `HashMap<String, HashSet<isize>>` where:
/// - The keys are chain identifiers (as strings).
/// - The values are `HashSet`s of residue serial numbers (as `isize`) that are part of the interface for that chain.
///
/// # Algorithm
///
/// 1. Iterates over all pairs of chains in the PDB structure.
/// 2. For each pair of chains, compares all residues between the two chains.
/// 3. For each pair of residues, checks if any pair of atoms (one from each residue) is within the cutoff distance.
/// 4. If atoms are within the cutoff, both residues are added to the interface set for their respective chains.
///
/// # Notes
///
/// - This function performs an exhaustive search, which may be computationally expensive for large structures.
/// - The cutoff is applied to atom-atom distances, not residue-residue distances.
/// - A residue is considered part of the interface if any of its atoms are within the cutoff of any atom in a residue from another chain.
/// - The function does not distinguish between different types of atoms or residues.
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

/// Identifies pairs of chains in a PDB structure that are in contact with each other.
///
/// This function determines which chains are in contact by finding pairs of residues from different chains
/// that have atoms within a specified cutoff distance of each other.
///
/// # Arguments
///
/// * `pdb` - A reference to a `pdbtbx::PDB` structure representing the entire protein.
/// * `cutoff` - A `f64` value specifying the maximum distance (in Angstroms) between atoms
///   for chains to be considered in contact.
///
/// # Returns
///
/// A `HashSet<(String, String)>` where each tuple represents a pair of chain identifiers that are in contact.
/// The order of chain identifiers in each tuple is arbitrary.
///
/// # Algorithm
///
/// 1. Iterates over all pairs of chains in the PDB structure.
/// 2. For each pair of chains, compares all residues between the two chains.
/// 3. For each pair of residues, checks if any pair of atoms (one from each residue) is within the cutoff distance.
/// 4. If atoms are within the cutoff, the pair of chain identifiers is added to the result set.
///
/// # Notes
///
/// - This function performs an exhaustive search, which may be computationally expensive for large structures.
/// - The cutoff is applied to atom-atom distances, not residue-residue or chain-chain distances.
/// - A pair of chains is considered in contact if any atom from one chain is within the cutoff distance of any atom from the other chain.
/// - The function does not distinguish between different types of atoms or residues.
/// - Each pair of chains appears only once in the result set, regardless of the number of contacts between them.
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

/// Represents two residues in a protein chain.
///
/// A `Gap` is defined by two residues that are sequential in the primary structure
/// but have a distance between their atoms that exceeds a certain threshold in the
/// tertiary structure.
#[derive(Debug)]
pub struct Gap {
    /// The chain identifier where the gap is located.
    pub chain: String,

    /// The name of the atom in the first residue used for distance calculation.
    pub atom_i: String,

    /// The name of the atom in the second residue used for distance calculation.
    pub atom_j: String,

    /// The sequence number of the first residue.
    pub res_i: isize,

    /// The sequence number of the second residue.
    pub res_j: isize,

    /// The distance between the specified atoms of the two residues.
    pub distance: f64,
}

/// Filters a vector of `Gap` elements, retaining only the unique elements based on the `atom_j` value.
///
/// This function removes duplicate `Gap` elements where the `atom_j` value is the same across multiple entries,
/// keeping only the first occurrence of each unique `atom_j`. This is useful when you only need to retain
/// one `Gap` entry for each unique `atom_j` regardless of the `atom_i` or other fields.
///
pub fn filter_unique_by_atom_j(gaps: Vec<Gap>) -> Vec<Gap> {
    let mut seen = HashSet::new();
    let mut unique: Vec<Gap> = Vec::new();

    gaps.into_iter().for_each(|g| {
        if !seen.contains(&g.atom_j) {
            seen.insert(g.atom_j.clone());
            unique.push(g)
        }
    });

    unique
}

/// Identifies separate bodies within a protein structure based on CA atom distances.
///
/// This function segments the protein into separate "bodies" by analyzing the distances
/// between consecutive CA (alpha carbon) atoms in each chain. A new body is defined when
/// the distance between consecutive CA atoms exceeds 4 Angstroms.
///
/// # Arguments
///
/// * `pdb` - A reference to a `pdbtbx::PDB` structure representing the entire protein.
///
/// # Returns
///
/// A `HashMap<isize, Vec<(isize, &str, &pdbtbx::Atom)>>` where:
/// - The key is a body identifier (starting from 0).
/// - The value is a vector of tuples, each containing:
///   - The residue serial number (`isize`)
///   - The chain identifier (`&str`)
///   - A reference to the CA atom (`&pdbtbx::Atom`)
///
/// # Algorithm
///
/// 1. Extracts all CA atoms from the PDB structure, along with their chain and residue information.
/// 2. Iterates through consecutive CA atoms within each chain.
/// 3. If the distance between consecutive CA atoms exceeds 4 Angstroms, a new body is started.
/// 4. Groups CA atoms into bodies based on this distance criterion.
///
/// # Notes
///
/// - This function only considers CA atoms for body identification.
/// - The 4 Angstrom cutoff is hard-coded and may not be suitable for all use cases.
/// - Bodies are numbered sequentially starting from 0.
/// - Gaps between chains always result in a new body, regardless of spatial proximity.
/// - The function assumes that CA atoms within a chain are listed in sequence order.
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

/// Creates a list of gaps between different bodies in a protein structure.
///
/// This function analyzes the spatial relationships between separate bodies of a protein
/// structure and generates a list of `Gap` instances representing the distances between
/// randomly selected pairs of atoms from different bodies.
///
/// # Arguments
///
/// * `bodies` - A reference to a `HashMap<isize, Vec<(isize, &str, &pdbtbx::Atom)>>` representing
///   the bodies of the protein structure, as returned by the `find_bodies` function.
///
/// # Returns
///
/// A `Vec<Gap>` containing `Gap` instances that represent the distances between pairs of
/// atoms from different bodies.
///
/// # Algorithm
///
/// 1. Initializes a random number generator with a fixed seed (42) for reproducibility.
/// 2. Iterates over all unique pairs of bodies.
/// 3. For each pair of bodies with at least 2 atoms each:
///    a. Randomly selects 2 atoms from each body.
///    b. Creates 2 `Gap` instances, one for each pair of selected atoms.
/// 4. Collects all created `Gap` instances into a vector.
///
/// # Notes
///
/// - The function uses a fixed random seed (42) for reproducibility. This means that
///   for the same input, it will always produce the same output.
/// - Only bodies with at least 2 atoms are considered for gap creation.
/// - The function creates 2 gaps for each pair of bodies, using different randomly selected atoms.
/// - The `Gap` instances contain information about the chain, atom names, residue numbers,
///   and the distance between the selected atoms.
/// - This function assumes that the input `bodies` HashMap is structured as expected,
///   with each value being a vector of tuples (residue number, chain ID, atom reference).
///
/// # Safety
///
/// This function uses `unwrap()` on the results of `choose_multiple()`. It assumes that
/// the bodies contain at least 2 atoms each, as checked earlier in the function.
pub fn create_iter_body_gaps(
    bodies: &HashMap<isize, Vec<(isize, &str, &pdbtbx::Atom)>>,
) -> Vec<Gap> {
    let mut rng = StdRng::seed_from_u64(42);
    let mut pairs: Vec<Gap> = Vec::new();
    let mut body_ids: Vec<isize> = bodies.keys().cloned().collect();
    body_ids.sort();
    for (i, &body_id1) in body_ids.iter().enumerate() {
        for &body_id2 in body_ids[i + 1..].iter() {
            if let (Some(atoms1), Some(atoms2)) = (bodies.get(&body_id1), bodies.get(&body_id2))
                && atoms1.len() >= 2 && atoms2.len() >= 2 {
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

/// Loads a PDB structure from a file path
///
/// # Arguments
/// * `input_pdb` - Path to the PDB file
///
/// # Returns
/// Result containing either the parsed PDB or a list of errors
pub fn load_pdb(input_pdb: &str) -> Result<PDB, Vec<PDBError>> {
    std::fs::read_to_string(input_pdb)
        .map_err(|e| {
            vec![PDBError::new(
                pdbtbx::ErrorLevel::GeneralWarning,
                "File read error",
                format!("Failed to read file {}: {}", input_pdb, e),
                pdbtbx::Context::None,
            )]
        })
        .and_then(|content| process_pdb_contents(&content))
}

/// Loads a PDB structure from string content
///
/// # Arguments
/// * `content` - PDB file contents as a string
///
/// # Returns
/// Result containing either the parsed PDB or a list of errors
pub fn load_pdb_from_content(content: &str) -> Result<PDB, Vec<PDBError>> {
    process_pdb_contents(content)
}

/// Processes raw PDB content by removing remarks and padding lines
///
/// # Arguments
/// * `content` - Raw PDB content
///
/// # Returns
/// Processed PDB content ready for parsing
pub fn process_pdb_contents(pdb_content: &str) -> Result<PDB, Vec<PDBError>> {
    let processed_content = process_pdb_string(pdb_content);

    let mut opts = ReadOptions::new();
    opts.set_format(pdbtbx::Format::Pdb)
        .set_level(pdbtbx::StrictnessLevel::Loose);

    let cursor = Cursor::new(processed_content.into_bytes());
    let reader = BufReader::new(cursor);

    match opts.read_raw(reader) {
        Ok((pdb, _)) => Ok(pdb),
        Err(e) => Err(e),
    }
}

/// Processes raw PDB content by removing remarks and padding lines
///
/// # Arguments
/// * `content` - Raw PDB content
///
/// # Returns
/// Processed PDB content ready for parsing
fn process_pdb_string(content: &str) -> String {
    let mut output = String::with_capacity(content.len());

    // Process lines with single allocation
    for line in content.lines() {
        if !line.starts_with("REMARK") {
            output.push_str(line);
            if line.len() < 80 {
                output.push_str(&" ".repeat(80 - line.len()));
            }
            output.push('\n');
        }
    }

    output
}

/// Finds the closest residue pairs between different protein chains within a specified distance cutoff.
///
/// This function analyzes inter-chain interactions by identifying the closest atom pairs between residues
/// from different chains, using the alpha carbon (CA) atoms for distance calculation.
///
/// # Arguments
///
/// * pdb - A reference to a parsed PDB structure.
/// * cutoff - A floating-point value specifying the maximum distance (in Angstroms) for considering residue pairs.
///
/// # Returns
///
/// A Vec<Pair> containing the closest residue pairs within the specified cutoff, sorted by distance.
///
/// # Functionality
///
/// 1. Iterates through all unique chain pairs in the PDB structure.
/// 2. For each chain pair, finds the closest atoms between residues.
/// 3. Selects pairs where the inter-residue distance is less than the specified cutoff.
/// 4. Uses alpha carbon (CA) atoms for distance calculation and pair representation.
/// 5. Sorts the resulting pairs by their inter-alpha carbon distance.
///
/// # Notes
///
/// - Only inter-chain interactions are considered.
/// - The distance is calculated based on the closest atom pairs between residues.
/// - The resulting pairs are sorted from shortest to longest distance.
pub fn get_closest_residue_pairs(pdb: &pdbtbx::PDB, cutoff: f64) -> Vec<Pair> {
    let mut closest_pairs = Vec::new();

    let chains: Vec<_> = pdb.chains().collect();

    for i in 0..chains.len() {
        for j in (i + 1)..chains.len() {
            let chain_i = &chains[i];
            let chain_j = &chains[j];

            for res_i in chain_i.residues() {
                let mut min_distance = f64::MAX;
                let mut closest_pair = None;

                for res_j in chain_j.residues() {
                    let mut atom_dist = f64::MAX;
                    for atom_i in res_i.atoms() {
                        for atom_j in res_j.atoms() {
                            let distance = atom_i.distance(atom_j);
                            if distance < atom_dist {
                                atom_dist = distance;
                            }
                        }
                    }
                    if atom_dist < min_distance {
                        min_distance = atom_dist;
                        closest_pair = Some((res_j, res_i));
                    }
                }

                if min_distance < cutoff
                    && let Some((res_j, res_i)) = closest_pair {
                        let ca_i = res_i.atoms().find(|atom| atom.name() == "CA").unwrap();
                        let ca_j = res_j.atoms().find(|atom| atom.name() == "CA").unwrap();
                        let ca_ca_dist = ca_i.distance(ca_j);
                        closest_pairs.push(Pair {
                            chain_i: chain_i.id().to_string(),
                            res_i: res_i.serial_number(),
                            atom_i: ca_i.name().to_string(),
                            chain_j: chain_j.id().to_string(),
                            res_j: res_j.serial_number(),
                            atom_j: ca_j.name().to_string(),
                            distance: ca_ca_dist,
                        });
                    }
            }
        }
    }

    closest_pairs.sort_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap());

    closest_pairs
}

/// Finds pairs of atoms between a ligand and nearby protein residues within a certain distance.
///
/// This function processes the given PDB structure to identify heteroatoms (ligand atoms) and
/// their neighboring atoms within a specified distance threshold (5.0 Å). It then finds the closest
/// protein atom for each ligand atom and returns a vector of `Pair` structures containing information
/// about the ligand atom and its closest protein atom. The result is a list of ligand-protein pairs
/// and their corresponding distance, residue, and chain information.
///
pub fn gaps_around_ligand(pdb: &PDB) -> Vec<Gap> {
    let ligand_atoms: Vec<&Atom> = pdb
        .atoms()
        .filter(|a| a.hetero() && *a.element().unwrap() != Element::H)
        .collect();

    let ligand_residues: Vec<&Residue> = pdb
        .residues()
        .filter(|r| r.atoms().any(|a| a.hetero()))
        .collect();

    let neighbors = neighbor_search(pdb.clone(), ligand_residues, 5.0);
    let protein_atoms = get_atoms_from_resnumbers(pdb, &neighbors);

    let mut result: Vec<Gap> = Vec::new();

    let finder = Finder::new(pdb);

    for ligand_atom in ligand_atoms {
        let closest_protein_atom = protein_atoms
            .iter()
            .filter(|a| *a.element().unwrap() != Element::H)
            .map(|a| (a, ligand_atom.distance(a)))
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .unwrap();

        let ligand_chain = finder.find_chain_id(ligand_atom);
        let protein_chain = finder.find_chain_id(closest_protein_atom.0);

        if ligand_chain == protein_chain {
            result.push(Gap {
                chain: ligand_chain,
                res_i: finder.find_residue_number(ligand_atom),
                atom_i: ligand_atom.name().to_string(),
                res_j: finder.find_residue_number(closest_protein_atom.0),
                atom_j: closest_protein_atom.0.name().to_string(),
                distance: closest_protein_atom.1,
            })
        }
    }

    result
}

#[cfg(test)]
mod tests {

    use std::env;

    use tempfile::NamedTempFile;

    use super::*;

    #[test]
    fn test_get_true_interface() {
        let pdb = load_pdb("tests/data/complex.pdb").unwrap();

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

    #[test]
    fn test_load_pdb_short_lines() {
        let pdb_path = env::current_dir()
            .unwrap()
            .join("tests/data/leu_short_lines.pdb");
        let pdb = load_pdb(pdb_path.to_str().unwrap()).unwrap();
        assert_eq!(pdb.atoms().count(), 8);
    }

    #[test]
    fn test_load_pdb_normal_lines() {
        let pdb_path = env::current_dir()
            .unwrap()
            .join("tests/data/leu_normal_lines.pdb");
        let pdb = load_pdb(pdb_path.to_str().unwrap()).unwrap();
        assert_eq!(pdb.atoms().count(), 8);
    }
    // Updated test helper function that works with older Rust versions
    fn create_test_pdb(remarks: usize, atoms: usize) -> String {
        let mut pdb = String::new();

        // Add remarks
        for i in 0..remarks {
            pdb.push_str(&format!("REMARK {} Test remark\n", i));
        }

        // Add atoms (ATOM records)
        for i in 0..atoms {
            pdb.push_str(&format!(
                "ATOM  {:>5}  N   ALA A{:>4}    {:>8}{:>8}{:>8}  1.00  0.00           N\n",
                i + 1,
                i + 1,
                format!("{:.3}", i as f32),
                format!("{:.3}", i as f32),
                format!("{:.3}", i as f32)
            ));
        }

        pdb
    }

    #[test]
    fn test_process_pdb_string_removes_remarks() {
        let input = "REMARK 1 Test\nATOM      1  N   ALA A   1       0.000   0.000   0.000\n";
        let processed = process_pdb_string(input);

        assert!(!processed.contains("REMARK"));
        assert!(processed.contains("ATOM"));
    }

    #[test]
    fn test_process_pdb_string_pads_lines() {
        let short_line = "ATOM      1  N   ALA A   1       0.000   0.000   0.000";
        let processed = process_pdb_string(short_line);

        // Check line is padded to 80 chars (including newline)
        assert_eq!(processed.lines().next().unwrap().len(), 80);
        assert_eq!(processed.chars().count(), 81); // 80 + newline
    }

    #[test]
    fn test_process_pdb_string_preserves_long_lines() {
        let long_line =
            "ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N  EXTRA";
        let processed = process_pdb_string(long_line);

        // Line should remain unchanged (except newline)
        assert_eq!(processed.trim(), long_line);
    }

    #[test]
    fn test_process_pdb_contents_valid_pdb() {
        let pdb_content = create_test_pdb(3, 5);
        let result = process_pdb_contents(&pdb_content);

        assert!(result.is_ok());
        let pdb = result.unwrap();
        assert_eq!(pdb.atoms().count(), 5);
    }

    #[test]
    fn test_process_pdb_contents_empty_input() {
        let result = process_pdb_contents("");
        assert!(result.is_err());

        if let Err(errors) = result {
            assert_eq!(errors.len(), 1);
        }
    }

    #[test]
    fn test_process_pdb_contents_invalid_pdb() {
        let invalid_content = "This is not a PDB file\n";
        let result = process_pdb_contents(invalid_content);

        assert!(result.is_err());
        if let Err(errors) = result {
            assert!(!errors.is_empty());
        }
    }

    #[test]
    fn test_process_pdb_string_empty_input() {
        let processed = process_pdb_string("");
        assert_eq!(processed, "");
    }

    #[test]
    fn test_process_pdb_string_mixed_content() {
        let input = "\
REMARK 1 Test
ATOM      1  N   ALA A   1       0.000   0.000   0.000
REMARK 2 Another
ATOM      2  N   ALA A   2       1.000   1.000   1.000
";
        let processed = process_pdb_string(input);
        let lines: Vec<&str> = processed.lines().collect();

        assert_eq!(lines.len(), 2); // Only ATOM lines should remain
        assert!(lines[0].starts_with("ATOM"));
        assert!(lines[1].starts_with("ATOM"));
    }

    #[test]
    fn test_process_pdb_contents_with_real_file() {
        // Create a temporary PDB file
        let mut file = NamedTempFile::new().unwrap();
        let pdb_content = create_test_pdb(2, 3);
        write!(file, "{}", pdb_content).unwrap();

        // Test with file content
        let file_content = std::fs::read_to_string(file.path()).unwrap();
        let result = process_pdb_contents(&file_content);

        assert!(result.is_ok());
        let pdb = result.unwrap();
        assert_eq!(pdb.atoms().count(), 3);
    }

    #[test]
    fn test_get_closest_residue_pairs() {
        let pdb = load_pdb("tests/data/two_res.pdb").unwrap();
        let observed_pairs = get_closest_residue_pairs(&pdb, 5.0);
        let expected_pairs = &[Pair {
            chain_i: "A".to_string(),
            chain_j: "B".to_string(),
            atom_i: "CA".to_string(),
            atom_j: "CA".to_string(),
            res_i: 2,
            res_j: 10,
            distance: 9.1,
        }];
        assert_eq!(observed_pairs.len(), expected_pairs.len());
        let pair = &observed_pairs[0];
        assert_eq!(pair.chain_i, expected_pairs[0].chain_i);
        assert_eq!(pair.chain_j, expected_pairs[0].chain_j);
        assert_eq!(pair.atom_i, expected_pairs[0].atom_i);
        assert_eq!(pair.atom_j, expected_pairs[0].atom_j);
        assert_eq!(pair.res_i, expected_pairs[0].res_i);
        assert_eq!(pair.res_j, expected_pairs[0].res_j);
        assert!((pair.distance - 9.1).abs() < 0.1);
    }
    #[test]
    fn test_gaps_around_ligand() {
        let pdb = load_pdb("tests/data/prot_ligand.pdb").unwrap();
        let observed = gaps_around_ligand(&pdb);
        let expected = &[
            Gap {
                chain: "E".to_string(),
                res_i: 351,
                atom_i: "N1".to_string(),
                res_j: 123,
                atom_j: "N".to_string(),
                distance: 3.0400192433601467,
            },
            Gap {
                chain: "E".to_string(),
                res_i: 351,
                atom_i: "C2".to_string(),
                res_j: 327,
                atom_j: "CE2".to_string(),
                distance: 3.79,
            },
        ];

        // Compare chain, residue, and atom names
        assert_eq!(observed[0].chain, expected[0].chain);
        assert_eq!(observed[0].res_i, expected[0].res_i);
        assert_eq!(observed[0].atom_i, expected[0].atom_i);
        assert_eq!(observed[0].res_j, expected[0].res_j);
        assert_eq!(observed[0].atom_j, expected[0].atom_j);

        assert_eq!(observed[1].chain, expected[1].chain);
        assert_eq!(observed[1].res_i, expected[1].res_i);
        assert_eq!(observed[1].atom_i, expected[1].atom_i);
        assert_eq!(observed[1].res_j, expected[1].res_j);
        assert_eq!(observed[1].atom_j, expected[1].atom_j);

        // TODO: Find a smart way to check the distances
        //
        // assert_eq!(observed[0].distance, expected[0].distance);
        // assert_eq!(observed[1].distance, expected[1].distance);
    }

    #[test]
    fn test_finder() {
        let pdb = load_pdb("tests/data/complex.pdb").unwrap();
        let finder = Finder::new(&pdb);

        assert_eq!(finder.chain_lookup.keys().len(), 924);
        assert_eq!(finder.residue_lookup.keys().len(), 924)
    }

    #[test]
    fn test_finder_find_chain_id() {
        let pdb = load_pdb("tests/data/complex.pdb").unwrap();
        let finder = Finder::new(&pdb);

        let found = finder.find_chain_id(pdb.atom(1).expect(""));
        assert_eq!(found, "A");

        let found = finder.find_chain_id(pdb.atom(700).expect(""));
        assert_eq!(found, "B");
    }

    #[test]
    fn test_finder_find_residue_number() {
        let pdb = load_pdb("tests/data/complex.pdb").unwrap();
        let finder = Finder::new(&pdb);

        let found = finder.find_residue_number(pdb.atom(1).unwrap());

        assert_eq!(found, 929)
    }

    #[test]
    fn test_filter_unique_by_atom_j() {
        let gaps = vec![
            Gap {
                chain: "A".to_string(),
                atom_i: "C".to_string(),
                atom_j: "B".to_string(),
                res_i: 1,
                res_j: 2,
                distance: 3.5,
            },
            Gap {
                chain: "A".to_string(),
                atom_i: "A".to_string(),
                atom_j: "B".to_string(),
                res_i: 3,
                res_j: 4,
                distance: 2.8,
            },
            Gap {
                chain: "B".to_string(),
                atom_i: "C".to_string(),
                atom_j: "D".to_string(),
                res_i: 5,
                res_j: 6,
                distance: 4.1,
            },
            Gap {
                chain: "A".to_string(),
                atom_i: "E".to_string(),
                atom_j: "B".to_string(),
                res_i: 7,
                res_j: 8,
                distance: 1.5,
            },
        ];

        let filtered = filter_unique_by_atom_j(gaps);

        // The function should remove duplicates based on `atom_j` and keep only the first occurrence.
        assert_eq!(filtered.len(), 2); // "B" and "D" should remain
        assert_eq!(filtered[0].atom_j, "B");
        assert_eq!(filtered[1].atom_j, "D");

        // Verify the correct `atom_j` is retained
        assert!(
            filtered
                .iter()
                .all(|gap| gap.atom_j != "B" || gap.atom_i == "C")
        );
    }

    #[test]
    fn test_no_duplicates() {
        let gaps = vec![
            Gap {
                chain: "A".to_string(),
                atom_i: "C".to_string(),
                atom_j: "B".to_string(),
                res_i: 1,
                res_j: 2,
                distance: 3.5,
            },
            Gap {
                chain: "A".to_string(),
                atom_i: "A".to_string(),
                atom_j: "D".to_string(),
                res_i: 3,
                res_j: 4,
                distance: 2.8,
            },
        ];

        let filtered = filter_unique_by_atom_j(gaps);
        assert_eq!(filtered.len(), 2); // No duplicates, so the length should be the same
    }

    #[test]
    fn test_empty_input() {
        let gaps: Vec<Gap> = Vec::new();
        let filtered = filter_unique_by_atom_j(gaps);
        assert!(filtered.is_empty()); // Should return an empty list
    }

    #[test]
    fn test_single_element() {
        let gaps = vec![Gap {
            chain: "A".to_string(),
            atom_i: "C".to_string(),
            atom_j: "B".to_string(),
            res_i: 1,
            res_j: 2,
            distance: 3.5,
        }];

        let filtered = filter_unique_by_atom_j(gaps);
        assert_eq!(filtered.len(), 1); // Single element, should remain
        assert_eq!(filtered[0].atom_j, "B");
    }
}
