use crate::sasa;
use crate::structure;
use core::panic;
use serde::Deserialize;
use std::collections::HashSet;

/// Represents an interactor in a molecular system.
///
/// This struct contains information about a specific interactor, including
/// its identification, residues, atoms, target interactions, and various
/// parameters for interaction calculations.
#[derive(Deserialize, Debug, Clone)]
pub struct Interactor {
    /// Unique identifier for the interactor.
    id: u16,

    /// Chain identifier for the interactor.
    chain: String,

    /// Set of active residue numbers.
    active: HashSet<i16>,

    /// Optional list of active atom names.
    active_atoms: Option<Vec<String>>,

    /// Set of passive residue numbers.
    pub passive: HashSet<i16>,

    /// Optional list of passive atom names.
    passive_atoms: Option<Vec<String>>,
    target: HashSet<u16>,

    /// Optional target distance for interactions.
    target_distance: Option<f64>,

    /// Optional lower margin for distance calculations.
    lower_margin: Option<f64>,

    /// Optional upper margin for distance calculations.
    upper_margin: Option<f64>,

    /// Optional path to the structure file.
    structure: Option<String>,

    /// Optional flag to determine if passive residues should be derived from active ones.
    passive_from_active: Option<bool>,

    /// Optional flag to treat surface residues as passive.
    surface_as_passive: Option<bool>,

    /// Optional flag to filter buried residues.
    filter_buried: Option<bool>,

    /// Optional cutoff value for buried residue filtering.
    filter_buried_cutoff: Option<f64>,
}

#[allow(clippy::too_many_arguments)]
impl Interactor {
    /// Creates a new `Interactor` instance with default values.
    ///
    /// This method initializes a new `Interactor` with the given ID and default values for all other fields.
    /// It's marked with `#[allow(clippy::too_many_arguments)]` to suppress warnings about the number of fields,
    /// even though this constructor doesn't actually take multiple arguments.
    ///
    /// # Arguments
    ///
    /// * `id` - A `u16` that specifies the unique identifier for the new `Interactor`.
    ///
    /// # Returns
    ///
    /// A new `Interactor` instance with the specified ID and default values for all other fields.
    ///
    pub fn new(id: u16) -> Self {
        Interactor {
            id,
            chain: String::new(),
            active: HashSet::new(),
            passive: HashSet::new(),
            target: HashSet::new(),
            structure: None,
            passive_from_active: None,
            surface_as_passive: None,
            filter_buried: None,
            filter_buried_cutoff: None,
            active_atoms: None,
            passive_atoms: None,
            wildcard: None,
            target_distance: None,
            lower_margin: None,
            upper_margin: None,
        }
    }

    /// Checks if the Interactor is in a valid state.
    ///
    /// This method performs two validity checks:
    /// 1. Ensures that the target set is not empty.
    /// 2. Verifies that there's no overlap between active and passive residues.
    ///
    /// # Returns
    ///
    /// - `Ok(true)` if the Interactor is valid.
    /// - `Err(&str)` with an error message if any validity check fails.
    ///
    pub fn is_valid(&self) -> Result<bool, &str> {
        if self.target.is_empty() {
            return Err("Target residues are empty");
        }
        if self.active.intersection(&self.passive).next().is_some() {
            return Err("Active/Passive selections overlap");
        }
        Ok(true)
    }

    /// Sets passive residues based on the active residues and their neighboring residues.
    ///
    /// This method performs the following steps:
    /// 1. Opens the PDB file specified in the `structure` field.
    /// 2. Retrieves the residues corresponding to the active residues.
    /// 3. Performs a neighbor search to find residues within 5.0 Å of the active residues.
    /// 4. Adds these neighboring residues to the passive set.
    ///
    /// # Panics
    ///
    /// This method will panic if:
    /// - The `structure` field is `None`.
    /// - There's an error opening or parsing the PDB file.
    ///
    /// # Side Effects
    ///
    /// This method modifies the `passive` set of the `Interactor`, adding new residues based on
    /// the neighbor search results.
    ///
    /// # Dependencies
    ///
    /// This method relies on external functions from the `pdbtbx` and `structure` modules:
    /// - `pdbtbx::open_pdb`
    /// - `structure::get_residues`
    /// - `structure::neighbor_search`
    ///
    pub fn set_passive_from_active(&mut self) {
        match pdbtbx::open_pdb(
            self.structure.clone().unwrap(),
            pdbtbx::StrictnessLevel::Loose,
        ) {
            Ok((pdb, _warnings)) => {
                let residues = structure::get_residues(
                    &pdb,
                    self.active.iter().map(|x| *x as isize).collect(),
                );

                let neighbors = structure::neighbor_search(pdb.clone(), residues, 5.0);

                // Add these neighbors to the passive set
                neighbors.iter().for_each(|x| {
                    self.passive.insert(*x as i16);
                });
            }
            Err(e) => {
                panic!("Error opening PDB file: {:?}", e);
            }
        }
    }

    /// Sets surface residues as passive based on their solvent accessible surface area (SASA).
    ///
    /// This method performs the following steps:
    /// 1. Opens the PDB file specified in the `structure` field.
    /// 2. Calculates the SASA for all residues in the structure.
    /// 3. Identifies surface residues (those with relative SASA > 0.7) on the same chain as the interactor.
    /// 4. Adds these surface residues to the passive set.
    ///
    /// # Panics
    ///
    /// This method will panic if:
    /// - The `structure` field is `None`.
    /// - There's an error opening or parsing the PDB file.
    ///
    /// # Side Effects
    ///
    /// This method modifies the `passive` set of the `Interactor`, adding new residues based on
    /// the SASA calculation results.
    ///
    /// # Dependencies
    ///
    /// This method relies on external functions from the `pdbtbx` and `sasa` modules:
    /// - `pdbtbx::open_pdb`
    /// - `sasa::calculate_sasa`
    ///
    /// # Note
    ///
    /// The threshold for considering a residue as "surface" is set to 0.7 relative SASA.
    /// This value may need to be adjusted based on specific requirements.
    pub fn set_surface_as_passive(&mut self) {
        match pdbtbx::open_pdb(
            self.structure.clone().unwrap(),
            pdbtbx::StrictnessLevel::Loose,
        ) {
            Ok((pdb, _warnings)) => {
                let sasa = sasa::calculate_sasa(pdb.clone());

                // Add these neighbors to the passive set
                sasa.iter().for_each(|r| {
                    // If the `rel_sasa_total` is more than 0.7 then add it to the passive set
                    if r.rel_sasa_total > 0.7 && r.chain == self.chain {
                        self.passive.insert(r.residue.serial_number() as i16);
                    }
                });
            }
            Err(e) => {
                panic!("Error opening PDB file: {:?}", e);
            }
        }
    }

    /// Removes buried residues from both active and passive sets based on solvent accessible surface area (SASA).
    ///
    /// This method performs the following steps:
    /// 1. Opens the PDB file specified in the `structure` field.
    /// 2. Calculates the SASA for all residues in the structure.
    /// 3. Identifies buried residues (those with relative SASA below a certain cutoff) on the same chain as the interactor.
    /// 4. Removes these buried residues from both the active and passive sets.
    ///
    /// # Panics
    ///
    /// This method will panic if:
    /// - The `structure` field is `None`.
    /// - There's an error opening or parsing the PDB file.
    ///
    /// # Side Effects
    ///
    /// This method modifies both the `active` and `passive` sets of the `Interactor`,
    /// removing residues based on the SASA calculation results.
    ///
    /// # Dependencies
    ///
    /// This method relies on external functions from the `pdbtbx` and `sasa` modules:
    /// - `pdbtbx::open_pdb`
    /// - `sasa::calculate_sasa`
    ///
    /// # Note
    ///
    /// The default threshold for considering a residue as "buried" is set to 0.7 relative SASA.
    /// This can be customized by setting the `filter_buried_cutoff` field of the `Interactor`.
    pub fn remove_buried_residues(&mut self) {
        match pdbtbx::open_pdb(
            self.structure.clone().unwrap(),
            pdbtbx::StrictnessLevel::Loose,
        ) {
            Ok((pdb, _warnings)) => {
                let sasa = sasa::calculate_sasa(pdb.clone());

                let sasa_cutoff = self.filter_buried_cutoff.unwrap_or(0.7);

                sasa.iter().for_each(|r| {
                    // If the `rel_sasa_total` is more than 0.7 then add it to the passive set
                    if r.rel_sasa_total < sasa_cutoff && r.chain == self.chain {
                        // This residue is not accessible, remove it from the passive and active sets
                        self.passive.remove(&(r.residue.serial_number() as i16));
                        self.active.remove(&(r.residue.serial_number() as i16));
                    }
                });
            }
            Err(e) => {
                panic!("Error opening PDB file: {:?}", e);
            }
        }
    }

    /// Returns the unique identifier of the Interactor.
    ///
    /// # Returns
    ///
    /// A `u16` representing the ID of the Interactor.
    pub fn id(&self) -> u16 {
        self.id
    }

    /// Returns the chain identifier of the Interactor.
    ///
    /// # Returns
    ///
    /// A string slice (`&str`) representing the chain of the Interactor.
    pub fn chain(&self) -> &str {
        &self.chain
    }

    /// Returns a reference to the set of active residues.
    ///
    /// # Returns
    ///
    /// A reference to a `HashSet<i16>` containing the active residue numbers.
    pub fn active(&self) -> &HashSet<i16> {
        &self.active
    }

    /// Returns a reference to the set of passive residues.
    ///
    /// # Returns
    ///
    /// A reference to a `HashSet<i16>` containing the passive residue numbers.
    pub fn passive(&self) -> &HashSet<i16> {
        &self.passive
    }

    pub fn wildcard(&self) -> &str {
        match &self.wildcard {
            Some(wildcard) => wildcard,
            None => "",
        }
    }

    pub fn wildcard(&self) -> &str {
        match &self.wildcard {
            Some(wildcard) => wildcard,
            None => "",
        }
    }

    /// Returns a reference to the set of target interactor IDs.
    ///
    /// # Returns
    ///
    /// A reference to a `HashSet<u16>` containing the IDs of target interactors.
    pub fn target(&self) -> &HashSet<u16> {
        &self.target
    }

    /// Returns the structure file path of the Interactor.
    ///
    /// # Returns
    ///
    /// A string slice (`&str`) representing the structure file path, or an empty string if not set.
    pub fn structure(&self) -> &str {
        match &self.structure {
            Some(structure) => structure,
            None => "",
        }
    }

    /// Sets the structure file path for the Interactor.
    ///
    /// # Arguments
    ///
    /// * `structure` - A string slice containing the path to the structure file.
    pub fn set_structure(&mut self, structure: &str) {
        self.structure = Some(structure.to_string());
    }

    /// Sets the chain identifier for the Interactor.
    ///
    /// # Arguments
    ///
    /// * `chain` - A string slice containing the chain identifier.
    pub fn set_chain(&mut self, chain: &str) {
        self.chain = chain.to_string();
    }

    /// Sets the active residues for the Interactor.
    ///
    /// # Arguments
    ///
    /// * `active` - A vector of `i16` values representing the active residue numbers.
    pub fn set_active(&mut self, active: Vec<i16>) {
        self.active = active.into_iter().collect();
    }

    /// Sets the passive residues for the Interactor.
    ///
    /// # Arguments
    ///
    /// * `passive` - A vector of `i16` values representing the passive residue numbers.
    pub fn set_passive(&mut self, passive: Vec<i16>) {
        self.passive = passive.into_iter().collect();
    }

    pub fn set_wildcard(&mut self, wildcard: &str) {
        self.wildcard = Some(wildcard.to_string());
    }

    pub fn set_wildcard(&mut self, wildcard: &str) {
        self.wildcard = Some(wildcard.to_string());
    }

    /// Sets the target distance for the Interactor.
    ///
    /// # Arguments
    ///
    /// * `distance` - A `f64` value representing the target distance.
    pub fn set_target_distance(&mut self, distance: f64) {
        self.target_distance = Some(distance);
    }

    /// Sets the lower margin for the Interactor.
    ///
    /// # Arguments
    ///
    /// * `margin` - A `f64` value representing the lower margin.
    pub fn set_lower_margin(&mut self, margin: f64) {
        self.lower_margin = Some(margin);
    }

    /// Sets the upper margin for the Interactor.
    ///
    /// # Arguments
    ///
    /// * `margin` - A `f64` value representing the upper margin.
    pub fn set_upper_margin(&mut self, margin: f64) {
        self.upper_margin = Some(margin);
    }

    /// Returns whether passive residues should be derived from active residues.
    ///
    /// # Returns
    ///
    /// A `bool` indicating if passive residues should be derived from active ones.
    pub fn passive_from_active(&self) -> bool {
        self.passive_from_active.unwrap_or(false)
    }

    /// Returns whether surface residues should be treated as passive.
    ///
    /// # Returns
    ///
    /// A `bool` indicating if surface residues should be treated as passive.
    pub fn surface_as_passive(&self) -> bool {
        self.surface_as_passive.unwrap_or(false)
    }

    /// Returns whether buried residues should be filtered.
    ///
    /// # Returns
    ///
    /// A `bool` indicating if buried residues should be filtered.
    pub fn filter_buried(&self) -> bool {
        self.filter_buried.unwrap_or(false)
    }

    /// Adds a target interactor ID.
    ///
    /// # Arguments
    ///
    /// * `target` - A `u16` value representing the ID of the target interactor to add.
    pub fn add_target(&mut self, target: u16) {
        self.target.insert(target);
    }

    /// Sets the active atoms for the Interactor.
    ///
    /// # Arguments
    ///
    /// * `atoms` - A vector of `String`s representing the active atom names.
    pub fn set_active_atoms(&mut self, atoms: Vec<String>) {
        self.active_atoms = Some(atoms);
    }

    /// Sets the passive atoms for the Interactor.
    ///
    /// # Arguments
    ///
    /// * `atoms` - A vector of `String`s representing the passive atom names.
    pub fn set_passive_atoms(&mut self, atoms: Vec<String>) {
        self.passive_atoms = Some(atoms);
    }

    /// Creates a block of restraints for the Interactor.
    ///
    /// This method generates a string representation of restraints for the Interactor,
    /// based on its active residues and the provided target residues.
    ///
    /// # Arguments
    ///
    /// * `target_res` - A vector of tuples, each containing a chain identifier (&str)
    ///                  and a residue number (&i16) for the target residues.
    ///
    /// # Returns
    ///
    /// A `String` containing the formatted block of restraints.
    ///
    /// # Details
    ///
    /// The method performs the following steps:
    /// 1. Sorts the active residues and target residues.
    /// 2. Determines if multiline formatting is needed (for multiple target residues).
    /// 3. For each active residue:
    ///    - Formats atom strings for active and passive atoms.
    ///    - Creates assignment strings for each target residue.
    ///    - Adds distance constraints.
    ///
    /// The resulting block includes all necessary restraints formatted according
    /// to the Interactor's properties and the provided target residues.
    pub fn create_block(&self, passive_res: Vec<PassiveResidues>) -> String {
        let mut block = String::new();
        let mut _active: Vec<i16> = self.active().iter().cloned().collect();
        _active.sort();

        // Sort the target residues by residue number
        let mut passive_res: Vec<PassiveResidues> = passive_res.clone();
        passive_res.sort_by(|a, b| a.res_number.cmp(&b.res_number));

        // Check if need to use multiline separation
        let multiline = passive_res.len() > 1;

        for resnum in _active {
            let atom_str = format_atom_string(&self.active_atoms);

            let mut assign_str = format!(
                "assign ( resid {} and segid {}{} {})",
                resnum,
                self.chain(),
                atom_str,
                &self.wildcard()
            );

            if multiline {
                assign_str += "\n       (\n";
            }

            block.push_str(assign_str.as_str());

            // panic!("Target res: {:?}", target_res);

            let res_lines: Vec<String> = passive_res
                .iter()
                .enumerate()
                .map(|(index, res)| {
                    let passive_atom_str = format_atom_string(&self.passive_atoms);

                    let mut res_line = String::new();
                    if multiline {
                        res_line.push_str(
                            format!(
                                "        ( {} segid {}{} {})\n",
                                res.res_number
                                    .map_or(String::new(), |num| format!("resid {} and", num)),
                                res.chain_id,
                                passive_atom_str,
                                res.wildcard
                            )
                            .as_str(),
                        );
                    } else {
                        res_line.push_str(
                            format!(
                                " ( {} segid {}{} {})",
                                res.res_number
                                    .map_or(String::new(), |num| format!("resid {} and", num)),
                                res.chain_id,
                                passive_atom_str,
                                res.wildcard
                            )
                            .as_str(),
                        );
                    }

                    if index != passive_res.len() - 1 {
                        res_line.push_str("     or\n");
                    }
                    res_line
                })
                .collect();

            block.push_str(&res_lines.join(""));

            let distance_string = format_distance_string(
                &self.target_distance,
                &self.lower_margin,
                &self.upper_margin,
            );
            if multiline {
                block.push_str(format!("       ) {}\n\n", distance_string).as_str());
            } else {
                block.push_str(format!(" {}\n\n", distance_string).as_str())
            }
        }
        block
    }
}

#[derive(Debug, Clone)]
pub struct PassiveResidues<'a> {
    chain_id: &'a str,
    res_number: Option<i16>,
    wildcard: &'a str,
}

/// Collects residue numbers from a vector of Interactors.
///
/// This function gathers both active and passive residue numbers from each Interactor,
/// along with their corresponding chain identifiers.
///
/// # Arguments
///
/// * `interactors` - A vector of references to Interactor objects.
///
/// # Returns
///
/// A vector of tuples, where each tuple contains:
/// - A string slice representing the chain identifier
/// - A reference to an i16 representing the residue number
///
pub fn collect_residues(interactors: Vec<&Interactor>) -> Vec<PassiveResidues> {
    let mut resnums = Vec::new();
    for interactor in interactors {
        let active = interactor.active().iter().map(|&x| PassiveResidues {
            chain_id: interactor.chain(),
            res_number: Some(x),
            wildcard: interactor.wildcard(),
        });

        let passive = interactor.passive().iter().map(|&x| PassiveResidues {
            chain_id: interactor.chain(),
            res_number: Some(x),
            wildcard: interactor.wildcard(),
        });

        resnums.extend(active);
        resnums.extend(passive);

        // If both active and passive are empty, add a single ResidueIdentifier with None as res_number
        if interactor.active().is_empty() && interactor.passive().is_empty() {
            resnums.push(PassiveResidues {
                chain_id: interactor.chain(),
                res_number: None,
                wildcard: interactor.wildcard(),
            });
        }
    }
    resnums
}

/// Formats a distance string based on target, lower, and upper bounds.
///
/// This function creates a formatted string representing distance constraints.
/// If any of the input values are None, default values are used.
///
/// # Arguments
///
/// * `target` - An Option<f64> representing the target distance.
/// * `lower` - An Option<f64> representing the lower bound of the distance.
/// * `upper` - An Option<f64> representing the upper bound of the distance.
///
/// # Returns
///
/// A String containing the formatted distance values, with one decimal place precision.
///
pub fn format_distance_string(
    target: &Option<f64>,
    lower: &Option<f64>,
    upper: &Option<f64>,
) -> String {
    let target = match target {
        Some(target) => target,
        None => &2.0,
    };

    let lower = match lower {
        Some(lower) => lower,
        None => &2.0,
    };

    let upper = match upper {
        Some(upper) => upper,
        None => &0.0,
    };

    format!("{:.1} {:.1} {:.1}", target, lower, upper)
}

/// Formats a string representing atom names for use in constraints.
///
/// This function takes an optional vector of atom names and formats them
/// into a string suitable for use in constraint definitions.
///
/// # Arguments
///
/// * `atoms` - An Option<Vec<String>> containing atom names.
///
/// # Returns
///
/// A String containing the formatted atom names, or an empty string if no atoms are provided.
///
pub fn format_atom_string(atoms: &Option<Vec<String>>) -> String {
    match atoms {
        Some(atoms) => {
            let atoms: Vec<String> = atoms.iter().map(|x| format!(" and name {}", x)).collect();
            atoms.join("")
        }
        None => "".to_string(),
    }
}

#[cfg(test)]
mod tests {

    use crate::interactor::{Interactor, PassiveResidues};

    #[test]
    fn test_create_block_multiline() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_chain("A");

        let observed = interactor.create_block(vec![
            PassiveResidues {
                chain_id: "B",
                res_number: Some(2),
                wildcard: "",
            },
            PassiveResidues {
                chain_id: "B",
                res_number: Some(3),
                wildcard: "",
            },
        ]);

        let block = "assign ( resid 1 and segid A )\n       (\n        ( resid 2 and segid B )\n     or\n        ( resid 3 and segid B )\n       ) 2.0 2.0 0.0\n\n";

        assert_eq!(observed, block);
    }

    #[test]
    fn test_create_block_oneline() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_chain("A");

        let observed = interactor.create_block(vec![PassiveResidues {
            chain_id: "B",
            res_number: Some(2),
            wildcard: "",
        }]);

        let block = "assign ( resid 1 and segid A ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n";

        assert_eq!(observed, block);
    }

    #[test]
    fn test_create_block_active_atoms() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_chain("A");
        interactor.set_active_atoms(vec!["CA".to_string()]);

        let observed = interactor.create_block(vec![PassiveResidues {
            chain_id: "B",
            res_number: Some(2),
            wildcard: "",
        }]);

        let block =
            "assign ( resid 1 and segid A and name CA ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n";

        assert_eq!(observed, block);
    }

    #[test]
    fn test_create_block_passive_atoms() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_chain("A");
        interactor.set_passive_atoms(vec!["CA".to_string()]);

        let observed = interactor.create_block(vec![PassiveResidues {
            chain_id: "B",
            res_number: Some(2),
            wildcard: "",
        }]);

        let block =
            "assign ( resid 1 and segid A ) ( resid 2 and segid B and name CA ) 2.0 2.0 0.0\n\n";

        assert_eq!(observed, block);
    }

    #[test]
    fn test_create_block_active_passive_atoms() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_chain("A");
        interactor.set_active_atoms(vec!["CA".to_string()]);
        interactor.set_passive_atoms(vec!["CB".to_string()]);

        let observed = interactor.create_block(vec![PassiveResidues {
            chain_id: "B",
            res_number: Some(2),
            wildcard: "",
        }]);

        let block =
            "assign ( resid 1 and segid A and name CA ) ( resid 2 and segid B and name CB ) 2.0 2.0 0.0\n\n";

        assert_eq!(observed, block);
    }

    #[test]
    fn test_create_multiline_block_active_passive_atoms() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_chain("A");
        interactor.set_active_atoms(vec!["CA".to_string()]);
        interactor.set_passive_atoms(vec!["CB".to_string()]);

        let observed = interactor.create_block(vec![
            PassiveResidues {
                chain_id: "B",
                res_number: Some(2),
                wildcard: "",
            },
            PassiveResidues {
                chain_id: "B",
                res_number: Some(3),
                wildcard: "",
            },
        ]);

        let block = "assign ( resid 1 and segid A and name CA )\n       (\n        ( resid 2 and segid B and name CB )\n     or\n        ( resid 3 and segid B and name CB )\n       ) 2.0 2.0 0.0\n\n";

        assert_eq!(observed, block);
    }

    #[test]
    fn test_create_block_with_distance() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_chain("A");
        interactor.set_target_distance(5.0);
        interactor.set_lower_margin(0.0);

        let observed = interactor.create_block(vec![PassiveResidues {
            chain_id: "B",
            res_number: Some(2),
            wildcard: "",
        }]);

        let block = "assign ( resid 1 and segid A ) ( resid 2 and segid B ) 5.0 0.0 0.0\n\n";

        assert_eq!(observed, block);
    }
}
