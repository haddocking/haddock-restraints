use crate::Interactor;
use std::fs::File;
use std::io::BufReader;

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
        // If the `structure` field was found, load the structure into the
        //  interactor as a PDB object
        if !interactor.structure().is_empty() {
            // Construct the full path to the structure file
            let pdb_path = wd.join(interactor.structure());

            // Update the structure path to the full path
            match interactor.load_structure(pdb_path.to_str().unwrap()) {
                Ok(_) => {}
                Err(e) => {
                    eprintln!("Error opening {:?}", pdb_path);
                    eprintln!("{:?}", e);
                    std::process::exit(1)
                }
            }
        }
    });

    Ok(data)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    #[test]
    fn test_read_json_file_success() {
        // Create a temporary JSON file
        let json_content = r#"
        [
            {
                "id": 1,
                "target": [1, 2, 3],
                "structure": "test.pdb",
                "chain": "A",
                "active": [],
                "passive": []
            }
        ]
        "#;

        let mut json_file = NamedTempFile::new().unwrap();
        write!(json_file, "{}", json_content).unwrap();

        // Create a mock PDB file in the same directory
        let pdb_path = json_file.path().parent().unwrap().join("test.pdb");
        let mut pdb_file = File::create(&pdb_path).unwrap();
        write!(
            pdb_file,
            "ATOM      1  N   MET A   1      11.111  22.222  33.333  1.00  0.00           N"
        )
        .unwrap();

        // Test the function
        let result = read_json_file(json_file.path().to_str().unwrap());
        println!("{:?}", result);
        assert!(result.is_ok());

        let interactors = result.unwrap();
        assert_eq!(interactors.len(), 1);
        assert_eq!(interactors[0].structure(), pdb_path.to_string_lossy());
        assert!(interactors[0].pdb().is_some());
    }

    #[test]
    fn test_read_json_file_invalid_json() {
        // Create invalid JSON file
        let json_content = "invalid json";
        let mut json_file = NamedTempFile::new().unwrap();
        write!(json_file, "{}", json_content).unwrap();

        // Test the function
        let result = read_json_file(json_file.path().to_str().unwrap());
        assert!(result.is_err());
    }

    #[test]
    fn test_read_json_file_no_structure() {
        // Create JSON with no structure field
        let json_content = r#"
        [
            {
                "id": 1,
                "target": [1, 2, 3],
                "chain": "A",
                "active": [],
                "passive": []
            }
        ]
        "#;

        let mut json_file = NamedTempFile::new().unwrap();
        write!(json_file, "{}", json_content).unwrap();

        // Test the function
        let result = read_json_file(json_file.path().to_str().unwrap());
        assert!(result.is_ok());

        let interactors = result.unwrap();
        assert_eq!(interactors.len(), 1);
        assert!(interactors[0].structure().is_empty());
        assert!(interactors[0].pdb().is_none());
    }
}
