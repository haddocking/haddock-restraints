use crate::Interactor;
use std::fs::File;
use std::io::BufReader;

pub fn read_json_file(file_path: &str) -> Result<Vec<Interactor>, Box<dyn std::error::Error>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let data: Vec<Interactor> = serde_json::from_reader(reader)?;

    Ok(data)
}

pub fn read_args() -> Result<String, Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        return Err("Please provide a file path".into());
    }
    Ok(args[1].clone())
}
