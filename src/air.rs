use std::collections::HashSet;

use crate::{interactor, utils};

/// Represents the Air (Ambiguous Interaction Restraints) structure.
///
/// This struct holds a collection of interactors and provides methods
/// to find interaction partners and generate tables.
pub struct Air(Vec<interactor::Interactor>);

impl Air {
    /// Creates a new `Air` instance with the given interactors.
    ///
    /// # Arguments
    ///
    /// * `interactors` - A vector of `Interactor` objects.
    ///
    /// # Returns
    ///
    /// A new `Air` instance.
    pub fn new(interactors: Vec<interactor::Interactor>) -> Self {
        Air(interactors)
    }

    /// Finds potential interaction partners for a given interactor.
    ///
    /// # Arguments
    ///
    /// * `interactor` - The interactor to find partners for.
    ///
    /// # Returns
    ///
    /// A vector of references to `Interactor` objects that are potential partners.
    pub fn find_partners(
        &self,
        interactor: &interactor::Interactor,
    ) -> Vec<&interactor::Interactor> {
        self.0
            .iter()
            .filter(|x| x.target().contains(&interactor.id()) && x.id() != interactor.id())
            .collect()
    }

    /// Generates a table representation of the AIR data.
    ///
    /// This method validates all interactors and creates a string representation
    /// of the AIR data, including interaction blocks for each valid interactor
    /// with partners.
    ///
    /// # Returns
    ///
    /// * `Ok(String)` - A string containing the generated table.
    /// * `Err(&str)` - An error message if any interactor is invalid.
    pub fn gen_tbl(&self) -> Result<String, &str> {
        let mut tbl = String::new();

        for interactor in self.0.iter() {
            match interactor.is_valid() {
                Ok(_) => (),
                Err(err) => {
                    eprintln!("## Interactor {} is not valid ##", interactor.id());
                    return Err(err);
                }
            }

            let partners = self.find_partners(interactor);

            if partners.is_empty() {
                continue;
            }

            // let header = append_header(i);
            // tbl.push_str(&header);

            let target_res = interactor::collect_residues(partners);
            let block = interactor.create_block(target_res);
            tbl.push_str(&block);
        }
        Ok(tbl)
    }

    pub fn gen_pml(&self, output_f: &str) {
        let mut pml = String::new();

        // General settings
        pml.push_str("set label_size, 0\n");
        pml.push_str("set dash_gap, 0\n");
        pml.push_str("set dash_color, yellow\n");

        let mut active: HashSet<(i16, &str)> = HashSet::new();
        let mut passive: HashSet<(i16, &str)> = HashSet::new();

        for interactor in self.0.iter() {
            let partners = self.find_partners(interactor);

            if partners.is_empty() {
                continue;
            }

            interactor.active().iter().for_each(|r| {
                active.insert((*r, interactor.chain()));
            });

            let target_res = interactor::collect_residues(partners);

            target_res.iter().for_each(|r| {
                let resnum = r.res_number.unwrap_or(0);
                passive.insert((resnum, r.chain_id));
            });

            let block = interactor.make_pml_string(target_res);
            pml.push_str(&block);
        }

        pml.push_str("color white\n");
        passive.iter().for_each(|(resnum, chain)| {
            pml.push_str(format!("color green, (resi {} and chain {})\n", resnum, chain).as_str())
        });
        active.iter().for_each(|(resnum, chain)| {
            pml.push_str(format!("color red, (resi {} and chain {})\n", resnum, chain).as_str())
        });

        utils::write_string_to_file(&pml, output_f).expect("Could not write pml")
    }
}

/// Generates a header string for AIR restraints.
///
/// This function is currently not used in the implementation (prefixed with _).
///
/// # Arguments
///
/// * `index` - The index of the selection, used in the header.
///
/// # Returns
///
/// A formatted string representing the header.
fn _append_header(index: usize) -> String {
    format!(
        "!========================================!\n\
         ! HADDOCK AIR restraints for selection {} !\n\
         !========================================!\n",
        index + 1
    )
}
