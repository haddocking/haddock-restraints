use crate::Interactor;
use std::fs::File;
use std::io::BufReader;
use std::process;

/// Reads and processes a JSON file containing Interactor data.
///
/// This function reads a JSON file, deserializes it into a vector of Interactors,
/// and updates the file paths for any structures referenced in the Interactors.
///
/// # Arguments
///
/// * `file_path` - A string slice that holds the path to the JSON file.
///
/// # Returns
///
/// * `Ok(Vec<Interactor>)` - A vector of Interactor objects if successful.
/// * `Err(Box<dyn std::error::Error>)` - An error if file reading or JSON parsing fails.
///
/// # Errors
///
/// This function will return an error if:
/// - The file cannot be opened
/// - The JSON is invalid or cannot be deserialized into `Vec<Interactor>`
/// - Any referenced structure file does not exist
///
/// # Panics
///
/// The function will call `process::exit(1)` if a referenced structure file does not exist.
///
pub fn read_json_file(file_path: &str) -> Result<Vec<Interactor>, Box<dyn std::error::Error>> {
    // Open the file and create a buffered reader
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    // Deserialize the JSON into a vector of Interactors
    let mut data: Vec<Interactor> = serde_json::from_reader(reader)?;

    // Get the directory of the input file
    let wd = std::path::Path::new(file_path).parent().unwrap();

    // Process each Interactor
    data.iter_mut().for_each(|interactor| {
        if !interactor.structure().is_empty() {
            // Construct the full path to the structure file
            let pdb_path = wd.join(interactor.structure());

            if pdb_path.exists() {
                // Update the structure path to the full path
                interactor.set_structure(pdb_path.to_str().unwrap());
            } else {
                // Print an error message and exit if the structure file doesn't exist
                eprintln!("\n### ERROR LOADING STRUCTURE ###");
                eprintln!("# The file `{}` does not exist", interactor.structure());
                eprintln!(
                    "# If you are using relative paths, they should be relative to `{}`\n",
                    file_path
                );
                process::exit(1);
            }

            // Set the structure path (this line seems redundant with the earlier set_structure call)
            interactor.set_structure(pdb_path.to_str().unwrap());
        }
    });

    Ok(data)
}
