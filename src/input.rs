use crate::Interactor;
use std::fs::File;
use std::io::BufReader;
use std::process;

pub fn read_json_file(file_path: &str) -> Result<Vec<Interactor>, Box<dyn std::error::Error>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let mut data: Vec<Interactor> = serde_json::from_reader(reader)?;

    let wd = std::path::Path::new(file_path).parent().unwrap();

    data.iter_mut().for_each(|interactor| {
        if !interactor.structure().is_empty() {
            let pdb_path = wd.join(interactor.structure());

            if pdb_path.exists() {
                interactor.set_structure(pdb_path.to_str().unwrap());
            } else {
                eprintln!("\n### ERROR LOADING STRUCTURE ###");
                eprintln!("# The file `{}` does not exist", interactor.structure());
                eprintln!(
                    "# If you are using relative paths, they should be relative to `{}`\n",
                    file_path
                );
                process::exit(1);
            }

            interactor.set_structure(pdb_path.to_str().unwrap());
        }
    });

    Ok(data)
}
