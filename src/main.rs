mod air;
mod interactor;
mod io;

use air::Air;
use interactor::Interactor;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let restraints_file = io::read_args()?;

    let interactors = io::read_json_file(&restraints_file).unwrap();

    let air = Air::new(interactors);

    let tbl = air.gen_tbl()?;

    println!("{}", tbl);

    Ok(())
}
