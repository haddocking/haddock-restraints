use crate::*;

/// Generates a table of Ambiguous Interaction Restraints (AIRs) based on input data.
///
/// This function reads interactor data from a JSON file, processes the interactors
/// according to specified flags, and generates a table of AIRs.
///
/// # Arguments
///
/// * `input_file` - A string slice that holds the path to the input JSON file.
///
/// # Functionality
///
/// 1. Reads interactor data from the specified JSON file.
/// 2. For each interactor with a non-empty structure:
///    a. If `passive_from_active` is set, derives passive residues from active ones.
///    b. If `surface_as_passive` is set, treats surface residues as passive.
///    c. If `filter_buried` is set, removes buried residues from consideration.
/// 3. Creates an `Air` instance from the processed interactors.
/// 4. Generates a table of AIRs.
/// 5. Prints the generated table to stdout.
///
pub fn gen_tbl(input_file: &str, pml: &Option<String>) {
    let mut interactors = read_json_file(input_file).unwrap();

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

    if let Some(output_f) = pml {
        air.gen_pml(output_f)
    };
}

#[cfg(test)]
mod tests {}
