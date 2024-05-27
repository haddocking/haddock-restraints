mod air;
mod interactor;

use air::Air;
use interactor::Interactor;

fn main() {
    let a = Interactor::new(1, "A", vec![1, 2], vec![1, 3], vec![3, 2]);
    let b = Interactor::new(2, "B", vec![10], vec![20], vec![1]);
    let c = Interactor::new(3, "C", vec![100], vec![200], vec![1]);

    let air = Air::new(vec![a, b, c]);

    match air.gen_tbl() {
        Ok(tbl) => println!("{}", tbl),
        Err(e) => println!("{}", e),
    }
}
