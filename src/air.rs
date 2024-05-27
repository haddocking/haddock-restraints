use crate::Interactor;

pub struct Air(Vec<Interactor>);

impl Air {
    pub fn new(interactors: Vec<Interactor>) -> Self {
        Air(interactors)
    }

    pub fn gen_tbl(&self) -> Result<String, &str> {
        let mut tbl = String::new();

        for (i, interactor) in self.0.iter().enumerate() {
            match interactor.is_valid() {
                Ok(_) => (),
                Err(err) => {
                    println!("Interactor {} is not valid", interactor.id());
                    return Err(err);
                }
            }

            // Find the partners, these are other interactors in which the target point to this id
            let partners: Vec<&Interactor> = self
                .0
                .iter()
                .filter(|x| x.target().contains(&interactor.id()))
                .collect();

            // If there are no partners, skip this interactor
            if partners.is_empty() {
                continue;
            }

            // Find the active and passive of the partners
            let mut target_res = Vec::<(&str, &i16)>::new();

            // Add the active residues of the partners to the target_res
            for partners in &partners {
                target_res.extend(partners.active().iter().map(|x| (partners.chain(), x)));
                target_res.extend(partners.passive().iter().map(|x| (partners.chain(), x)));
            }

            // Make the header
            tbl.push_str(&format!(
                "!========================================!\n! HADDOCK AIR restraints for selection {} !\n!========================================!\n",
                i + 1
            ));

            // Create the assign line
            for resnum in interactor.active() {
                let assign_line = format!(
                    "assign ( resid {} and segid {} )\n       (\n",
                    resnum,
                    interactor.chain()
                );

                let res_lines: Vec<String> = target_res
                    .iter()
                    .enumerate()
                    .map(|(index, res)| {
                        let mut res_line =
                            format!("        ( resid {}  and segid {})\n", res.1, res.0);
                        if index != target_res.len() - 1 {
                            res_line.push_str("     or\n");
                        }
                        res_line
                    })
                    .collect();

                let block = format!(
                    "{}{}       )  2.0 2.0 0.0\n\n",
                    assign_line,
                    res_lines.join("")
                );
                tbl.push_str(&block);
            }
        }
        Ok(tbl)
    }
}
