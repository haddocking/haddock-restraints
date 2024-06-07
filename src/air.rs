use crate::interactor;

pub struct Air(Vec<interactor::Interactor>);

impl Air {
    pub fn new(interactors: Vec<interactor::Interactor>) -> Self {
        Air(interactors)
    }

    pub fn find_partners(
        &self,
        interactor: &interactor::Interactor,
    ) -> Vec<&interactor::Interactor> {
        self.0
            .iter()
            .filter(|x| x.target().contains(&interactor.id()) && x.id() != interactor.id())
            .collect()
    }

    pub fn gen_tbl(&self) -> Result<String, &str> {
        let mut tbl = String::new();

        for interactor in self.0.iter() {
            // println!("# {:?} ", interactor.id());
            // println!("Active: {:?} ", interactor.active());
            // println!("Passive: {:?} ", interactor.passive());
            match interactor.is_valid() {
                Ok(_) => (),
                Err(err) => {
                    eprintln!("## Interactor {} is not valid ##", interactor.id());
                    return Err(err);
                }
            }

            let partners = self.find_partners(interactor);

            if partners.is_empty() {
                continue;
            }

            // let header = append_header(i);
            // tbl.push_str(&header);

            let target_res = interactor::collect_resnums(partners);
            println!("------ {:?}", interactor);
            println!("------ {:?}", target_res);
            let block = interactor.create_block(target_res);
            println!("#### {:?}", block);
            tbl.push_str(&block);
        }
        Ok(tbl)
    }
}

fn _append_header(index: usize) -> String {
    format!(
        "!========================================!\n\
         ! HADDOCK AIR restraints for selection {} !\n\
         !========================================!\n",
        index + 1
    )
}
