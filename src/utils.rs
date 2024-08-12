use std::io::BufReader;
use std::io::Cursor;

pub fn tidy_pdb(pdb_f: &str) -> BufReader<Cursor<Vec<u8>>> {
    let no_remarks_pdb_buf = pdb_handler::remove_remark(pdb_f);

    let padded_lines_pdb_buf = pdb_handler::pad_lines(&no_remarks_pdb_buf);

    padded_lines_pdb_buf
}
