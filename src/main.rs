mod air;
mod interactor;
mod io;
mod sasa;
mod structure;

use air::Air;
use interactor::Interactor;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let restraints_file = io::read_args()?;

    let mut interactors = io::read_json_file(&restraints_file).unwrap();

    let wd = std::path::Path::new(restraints_file.as_str())
        .parent()
        .unwrap()
        .to_str()
        .unwrap();

    interactors.iter_mut().for_each(|interactor| {
        if !interactor.structure().is_empty() {
            interactor.set_structure(format!("{}/{}", wd, interactor.structure()).as_str());

            if interactor.passive_from_active() {
                interactor.set_passive_from_active();
            }

            if interactor.surface_as_passive() {
                interactor.set_surface_as_passive();
            }
        }
    });

    let air = Air::new(interactors);

    let _tbl = air.gen_tbl()?;
    println!("{}", _tbl);

    Ok(())
}
