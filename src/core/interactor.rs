use crate::core::sasa;
use crate::core::structure;
use crate::load_pdb;
use pdbtbx::PDB;
use pdbtbx::PDBError;
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

    /// Set of target interactor IDs.
    target: HashSet<u16>,

    /// Optional target distance for interactions.
    target_distance: Option<f64>,

    /// Optional lower margin for distance calculations.
    lower_margin: Option<f64>,

    /// Optional upper margin for distance calculations.
    upper_margin: Option<f64>,

    /// Optional path to the structure file.
    structure: Option<String>,

    /// Optional PDB object.
    pdb: Option<PDB>,

    /// Optional flag to determine if passive residues should be derived from active ones.
    passive_from_active: Option<bool>,

    /// Optional radius to define the neighbor search radius
    passive_from_active_radius: Option<f64>,

    /// Optional flag to treat surface residues as passive.
    surface_as_passive: Option<bool>,

    /// Optional flag to filter buried residues.
    filter_buried: Option<bool>,

    /// Optional cutoff value for buried residue filtering.
    filter_buried_cutoff: Option<f64>,

    /// Optional wildcard value for the interactor.
    wildcard: Option<String>,
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
            pdb: None,
            passive_from_active: None,
            passive_from_active_radius: None,
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
    /// This method relies on external functions from the `structure` module:
    /// - `structure::get_residues`
    /// - `structure::neighbor_search`
    /// - `structure::load_pdb`
    ///
    pub fn set_passive_from_active(&mut self) {
        if let Some(pdb) = &self.pdb {
            let residues =
                structure::get_residues(pdb, self.active.iter().map(|x| *x as isize).collect());

            let search_cutoff = self.passive_from_active_radius.unwrap_or(6.5);
            let neighbors = structure::neighbor_search(pdb.clone(), residues, search_cutoff);

            // Add these neighbors to the passive set
            neighbors.iter().for_each(|x| {
                self.passive.insert(*x as i16);
            });
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
    /// This method relies on external functions from the `sasa` and `structure` modules:
    /// - `sasa::calculate_sasa`
    /// - `structure::load_pdb`
    ///
    /// # Note
    ///
    /// The threshold for considering a residue as "surface" is set to 0.7 relative SASA.
    /// This value may need to be adjusted based on specific requirements.
    pub fn set_surface_as_passive(&mut self) {
        if let Some(pdb) = &self.pdb {
            let sasa = sasa::calculate_sasa(pdb.clone());

            // Add these neighbors to the passive set
            sasa.iter().for_each(|r| {
                // If the `rel_sasa_total` is more than 0.7 then add it to the passive set
                if r.rel_sasa_total > 0.7 && r.chain == self.chain {
                    self.passive.insert(r.residue.serial_number() as i16);
                }
            });
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
    /// This method relies on external functions from the `sasa` and `structure` module:
    /// - `sasa::calculate_sasa`
    /// - `structure::load_pdb`
    ///
    /// # Note
    ///
    /// The default threshold for considering a residue as "buried" is set to 0.7 relative SASA.
    /// This can be customized by setting the `filter_buried_cutoff` field of the `Interactor`.
    pub fn remove_buried_residues(&mut self) {
        if let Some(pdb) = &self.pdb {
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

    /// Returns a reference to a set of active atoms strings.
    ///
    /// # Returns
    ///
    /// A reference to a `Option<Vec<String>>` containing active atom names.
    pub fn active_atoms(&self) -> &Option<Vec<String>> {
        &self.active_atoms
    }

    /// Returns a reference to the set of passive residues.
    ///
    /// # Returns
    ///
    /// A reference to a `HashSet<i16>` containing the passive residue numbers.
    pub fn passive(&self) -> &HashSet<i16> {
        &self.passive
    }

    /// Returns a reference to a set of passive atoms strings.
    ///
    /// # Returns
    ///
    /// A reference to a `Option<Vec<String>>` containing passive atom names.
    pub fn passive_atoms(&self) -> &Option<Vec<String>> {
        &self.passive_atoms
    }

    /// Returns the wildcard string associated with this Interactor.
    ///
    /// # Returns
    ///
    /// - If a wildcard is set, returns a string slice (`&str`) containing the wildcard value.
    /// - If no wildcard is set, returns an empty string slice.
    ///
    /// # Notes
    ///
    /// - This method provides read-only access to the wildcard value.
    /// - The wildcard is typically used to represent any residue or atom in certain contexts.
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

    /// Loads a PDB structure from the given path and stores it.
    ///
    /// # Returns
    ///
    /// `Ok(())` on success, or `Vec<PDBError>` if loading failed.
    pub fn load_structure(&mut self, structure_path: &str) -> Result<(), Vec<PDBError>> {
        match load_pdb(structure_path) {
            Ok(pdb) => {
                self.structure = Some(structure_path.to_string());
                self.pdb = Some(pdb);
                Ok(())
            }
            Err(e) => Err(e),
        }
    }

    /// Returns a reference to the stored PDB structure, if any.
    ///
    /// # Returns
    ///
    /// Reference to an `Option<PDB>` containing the structure.
    pub fn pdb(&self) -> &Option<PDB> {
        &self.pdb
    }

    /// Stores a PDB structure in the Interactor.
    ///
    /// # Arguments
    ///
    /// * `pdb` - The PDB structure to store
    pub fn set_pdb(&mut self, pdb: PDB) {
        self.pdb = Some(pdb)
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

    /// Sets the wildcard string for this Interactor.
    ///
    /// This method allows you to set or update the wildcard value associated with the Interactor.
    ///
    /// # Arguments
    ///
    /// * `wildcard` - A string slice (`&str`) that specifies the new wildcard value.
    ///
    /// # Notes
    ///
    /// - This method will overwrite any previously set wildcard value.
    /// - The wildcard is stored as an owned `String`, so the input `&str` is cloned.
    /// - An empty string is a valid wildcard value, though its interpretation may depend on the context.
    /// - The wildcard is typically used to represent any residue or atom in certain contexts.
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

    /// Sets whether passive residues should be derived from active residues.
    ///
    /// # Arguments
    ///
    /// * `cutoff` - A `f64` value representing the cutoff to filter buried residues.
    pub fn set_filter_buried_cutoff(&mut self, cutoff: f64) {
        self.filter_buried_cutoff = Some(cutoff);
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
    ///   and a residue number (&i16) for the target residues.
    ///
    /// # Returns
    ///
    /// A `String` containing the formatted block of restraints.
    ///
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
            // Create the `assign` statement
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

            // Loop over the passive residues
            let res_lines: Vec<String> = passive_res
                .iter()
                .enumerate()
                .map(|(index, res)| {
                    let atom_str = format_atom_string(res.atom_str);

                    let mut res_line = String::new();
                    if multiline {
                        res_line.push_str(
                            format!(
                                "        ( {} segid {}{} {})\n",
                                res.res_number
                                    .map_or(String::new(), |num| format!("resid {} and", num)),
                                res.chain_id,
                                atom_str,
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
                                atom_str,
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

    pub fn make_pml_string(&self, passive_res: Vec<PassiveResidues>) -> String {
        let mut pml = String::new();
        let mut _active: Vec<i16> = self.active().iter().cloned().collect();
        _active.sort();

        let mut passive_res: Vec<PassiveResidues> = passive_res.clone();
        passive_res.sort_by(|a, b| a.res_number.cmp(&b.res_number));

        for resnum in _active {
            let identifier = format!("{}-{}", resnum, self.chain);
            let active_sel = format!("resi {} and name CA and chain {}", resnum, self.chain);

            for passive_resnum in &passive_res {
                let passive_sel = format!(
                    "resi {} and name CA and chain {}",
                    passive_resnum.res_number.unwrap(),
                    passive_resnum.chain_id
                );

                pml.push_str(
                    format!(
                        "distance {}, ({}), ({})\n",
                        identifier, active_sel, passive_sel
                    )
                    .as_str(),
                )
            }
        }

        pml
    }
}

#[derive(Debug, Clone)]
pub struct PassiveResidues<'a> {
    pub chain_id: &'a str,
    pub res_number: Option<i16>,
    wildcard: &'a str,
    // TODO: ADD THE ATOM ATOM NAMES HERE, THEY SHOULD BE USED WHEN GENERATING THE BLOCK
    atom_str: &'a Option<Vec<String>>,
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
            atom_str: interactor.active_atoms(),
        });

        let passive = interactor.passive().iter().map(|&x| PassiveResidues {
            chain_id: interactor.chain(),
            res_number: Some(x),
            wildcard: interactor.wildcard(),
            atom_str: interactor.passive_atoms(),
        });

        resnums.extend(active);
        resnums.extend(passive);

        // If both active and passive are empty, add a single ResidueIdentifier with None as res_number
        if interactor.active().is_empty() && interactor.passive().is_empty() {
            resnums.push(PassiveResidues {
                chain_id: interactor.chain(),
                res_number: None,
                wildcard: interactor.wildcard(),
                atom_str: &None,
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
        Some(atoms) if atoms.len() > 1 => {
            let atoms: String = atoms
                .iter()
                .map(|x| {
                    if x.contains("-") || x.contains("+") {
                        format!(r#""{}""#, x)
                    } else {
                        format!("{}", x)
                    }
                })
                .collect::<Vec<String>>()
                .join(" or ");

            format!(" and name ({})", atoms)
        }
        Some(atoms) if atoms.len() == 1 => {
            if atoms[0].contains("-") || atoms[0].contains("+") {
                format!(r#" and name "{}""#, atoms[0])
            } else {
                format!(" and name {}", atoms[0])
            }
        },
        _ => "".to_string(),
    }
}

#[cfg(test)]
mod tests {

    use crate::core::interactor::{Interactor, PassiveResidues, format_atom_string};

    #[test]
    fn test_format_atom_string() {
        let atom_str = format_atom_string(&Some(vec!["O".to_string()]));
        let expected_atom_str = " and name O".to_string();
        assert_eq!(atom_str, expected_atom_str)
    }

    #[test]
    fn test_format_atom_string_multiple() {
        let atom_str = format_atom_string(&Some(vec!["O".to_string(), "CA".to_string()]));
        let expected_atom_str = " and name (O or CA)".to_string();
        assert_eq!(atom_str, expected_atom_str)
    }

    #[test]
    fn test_format_atom_string_special_chars() {
        let atom_str = format_atom_string(&Some(vec!["ZN+2".to_string()]));
        let expected_atom_str = " and name \"ZN+2\"".to_string();
        assert_eq!(atom_str, expected_atom_str)
    }

    #[test]
    fn test_format_atom_string_multiple_special_chars() {
        let atom_str = format_atom_string(&Some(vec!["ZN+2".to_string(), "FE-3".to_string()]));
        let expected_atom_str = " and name (\"ZN+2\" or \"FE-3\")".to_string();
        assert_eq!(atom_str, expected_atom_str)
    }

    #[test]
    fn test_format_atom_string_multiple_hybrid_chars() {
        let atom_str = format_atom_string(&Some(vec!["ZN+2".to_string(), "CA".to_string()]));
        let expected_atom_str = " and name (\"ZN+2\" or CA)".to_string();
        assert_eq!(atom_str, expected_atom_str)
    }

    #[test]
    fn test_valid_interactor() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_passive(vec![2]);
        interactor.add_target(2);

        assert_eq!(interactor.is_valid(), Ok(true));
    }

    #[test]
    fn test_invalid_interactor_empty() {
        let interactor = Interactor::new(1);

        assert_eq!(interactor.is_valid(), Err("Target residues are empty"));
    }

    #[test]
    fn test_invalid_interactor_overlap() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_passive(vec![1]);
        interactor.add_target(2);

        assert_eq!(
            interactor.is_valid(),
            Err("Active/Passive selections overlap")
        );
    }

    #[test]
    fn test_set_passive_from_active() {
        let mut interactor = Interactor::new(1);
        interactor.load_structure("tests/data/complex.pdb").unwrap();
        interactor.set_active(vec![1]);
        interactor.passive_from_active_radius = Some(5.0);
        interactor.set_passive_from_active();

        let expected_passive = [16, 15, 18, 3, 19, 61, 56, 17, 2, 62, 63];

        assert_eq!(
            interactor.passive(),
            &expected_passive.iter().cloned().collect()
        );
    }

    #[test]
    fn test_set_surface_as_passive() {
        let mut interactor = Interactor::new(1);
        interactor.load_structure("tests/data/complex.pdb").unwrap();
        interactor.set_chain("A");
        interactor.set_surface_as_passive();

        let expected_passive = [
            938, 965, 953, 944, 933, 958, 966, 972, 931, 936, 961, 929, 943, 954, 932, 945, 942,
            957, 955, 947, 940, 941, 937, 964, 970, 930, 969, 968, 950, 952, 959, 971, 967, 956,
            946, 960, 962, 935, 948, 951, 934, 939,
        ];

        assert_eq!(
            interactor.passive(),
            &expected_passive.iter().cloned().collect()
        );
    }

    #[test]
    fn test_remove_buried_active_residues() {
        let mut interactor = Interactor::new(1);

        interactor.load_structure("tests/data/complex.pdb").unwrap();
        interactor.set_chain("A");
        interactor.filter_buried = Some(true);
        interactor.filter_buried_cutoff = Some(0.7);
        interactor.set_active(vec![949, 931]);
        interactor.remove_buried_residues();

        let expected_active = [931];

        assert_eq!(
            interactor.active(),
            &expected_active.iter().cloned().collect()
        );
    }

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
                atom_str: &None,
            },
            PassiveResidues {
                chain_id: "B",
                res_number: Some(3),
                wildcard: "",
                atom_str: &None,
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
            atom_str: &None,
        }]);

        let block = "assign ( resid 1 and segid A ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n";

        assert_eq!(observed, block);
    }

    #[test]
    fn test_create_block_oneline_atom_subset() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_chain("A");
        interactor.set_active_atoms(vec!["CA".to_string(), "CB".to_string()]);

        let observed = interactor.create_block(vec![PassiveResidues {
            chain_id: "B",
            res_number: Some(2),
            wildcard: "",
            atom_str: &None,
        }]);

        let block = "assign ( resid 1 and segid A and (name CA or name CB) ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n";

        assert_eq!(observed, block);
    }

    #[test]
    fn test_create_block_multiline_atom_subset() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_chain("A");
        interactor.set_active_atoms(vec!["CA".to_string(), "CB".to_string()]);
        interactor.set_passive_atoms(vec!["CA".to_string(), "CB".to_string()]);
        let observed = interactor.create_block(vec![
            PassiveResidues {
                chain_id: "B",
                res_number: Some(2),
                wildcard: "",
                atom_str: &None,
            },
            PassiveResidues {
                chain_id: "B",
                res_number: Some(3),
                wildcard: "",
                atom_str: &None,
            },
        ]);

        let block = "assign ( resid 1 and segid A and (name CA or name CB) )\n       (\n        ( resid 2 and segid B )\n     or\n        ( resid 3 and segid B )\n       ) 2.0 2.0 0.0\n\n";

        assert_eq!(observed, block);
    }

    #[test]
    fn test_create_block_multiline_atom_subset_passive() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_chain("A");
        interactor.set_active_atoms(vec!["CA".to_string(), "CB".to_string()]);
        let observed = interactor.create_block(vec![
            PassiveResidues {
                chain_id: "B",
                res_number: Some(2),
                wildcard: "",
                atom_str: &Some(vec!["N".to_string(), "C".to_string()]),
            },
            PassiveResidues {
                chain_id: "B",
                res_number: Some(3),
                wildcard: "",
                atom_str: &None,
            },
        ]);

        let block = "assign ( resid 1 and segid A and (name CA or name CB) )\n       (\n        ( resid 2 and segid B and (name N or name C) )\n     or\n        ( resid 3 and segid B )\n       ) 2.0 2.0 0.0\n\n";

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
            atom_str: &None,
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

        let observed = interactor.create_block(vec![PassiveResidues {
            chain_id: "B",
            res_number: Some(2),
            wildcard: "",
            atom_str: &Some(vec!["CA".to_string()]),
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

        let observed = interactor.create_block(vec![PassiveResidues {
            chain_id: "B",
            res_number: Some(2),
            wildcard: "",
            atom_str: &None,
        }]);

        let block =
            "assign ( resid 1 and segid A and name CA ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n";

        assert_eq!(observed, block);
    }

    #[test]
    fn test_create_multiline_block_active_passive_atoms() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_chain("A");
        interactor.set_active_atoms(vec!["CA".to_string()]);

        let observed = interactor.create_block(vec![
            PassiveResidues {
                chain_id: "B",
                res_number: Some(2),
                wildcard: "",
                atom_str: &Some(vec!["CB".to_string()]),
            },
            PassiveResidues {
                chain_id: "B",
                res_number: Some(3),
                wildcard: "",
                atom_str: &Some(vec!["N".to_string()]),
            },
        ]);

        let block = "assign ( resid 1 and segid A and name CA )\n       (\n        ( resid 2 and segid B and name CB )\n     or\n        ( resid 3 and segid B and name N )\n       ) 2.0 2.0 0.0\n\n";

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
            atom_str: &None,
        }]);

        let block = "assign ( resid 1 and segid A ) ( resid 2 and segid B ) 5.0 0.0 0.0\n\n";

        assert_eq!(observed, block);
    }

    #[test]
    fn test_create_block_with_wildcard() {
        let mut interactor = Interactor::new(1);
        interactor.set_active(vec![1]);
        interactor.set_chain("A");
        interactor.set_wildcard("and attr z gt 42.00 ");

        let observed = interactor.create_block(vec![PassiveResidues {
            chain_id: "B",
            res_number: Some(2),
            wildcard: "",
            atom_str: &None,
        }]);

        let block = "assign ( resid 1 and segid A and attr z gt 42.00 ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n";

        assert_eq!(observed, block);
    }
}
