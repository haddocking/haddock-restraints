//! Pml-string-building helpers shared between `Air::gen_pml`
//! (config.json-driven visualization) and `tbl2pml` (.tbl-driven
//! visualization), so the two entry points can't drift apart.

use std::fmt::Write;

/// The PyMOL atom selector used to anchor a residue for `distance`/`color`
/// commands. `C1'` is the nucleic-acid equivalent of `CA` (see commit
/// b2763a1, "Fix `--pml` visualization for nucleic-acid structures").
pub fn atom_selector(resnum: i16, chain: &str) -> String {
    format!("resi {} and (name CA or name C1') and chain {}", resnum, chain)
}

/// The PyMOL display settings shared by every generated `.pml` script.
pub fn header() -> String {
    "set label_size, 0\nset dash_gap, 0\nset dash_color, yellow\n".to_string()
}

/// Colors every passive residue green and every active residue red.
pub fn color_footer<'a>(
    passive: impl Iterator<Item = (i16, &'a str)>,
    active: impl Iterator<Item = (i16, &'a str)>,
) -> String {
    let mut footer = String::from("color white\n");
    for (resnum, chain) in passive {
        let _ = writeln!(footer, "color green, (resi {} and chain {})", resnum, chain);
    }
    for (resnum, chain) in active {
        let _ = writeln!(footer, "color red, (resi {} and chain {})", resnum, chain);
    }
    footer
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atom_selector() {
        assert_eq!(
            atom_selector(1, "A"),
            "resi 1 and (name CA or name C1') and chain A"
        );
    }

    #[test]
    fn test_color_footer() {
        let footer = color_footer(vec![(2, "B")].into_iter(), vec![(1, "A")].into_iter());
        assert_eq!(
            footer,
            "color white\ncolor green, (resi 2 and chain B)\ncolor red, (resi 1 and chain A)\n"
        );
    }
}
