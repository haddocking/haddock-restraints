//! Parses HADDOCK `.tbl` AIR restraint files (the *inverse* of
//! `Interactor::create_block`, see `src/core/interactor.rs`).
//!
//! Only `resid <N> and segid <chain>` selections are extracted from each
//! `assign` block — that's all pml rendering (`tbl2pml`) needs. Extra
//! clauses like `and name CA` or `and attr z gt 42.00` are ignored.

/// A single residue selection extracted from a `.tbl` restraint block.
#[derive(Debug, Clone, PartialEq)]
pub struct ResSelection {
    pub resid: i16,
    pub chain: String,
}

/// One `assign` block: one active residue restrained against a list of
/// (possibly OR'd) partner residues.
///
/// ```text
/// assign ( resid 1 and segid A )
///        (
///         ( resid 2 and segid B )
///      or
///         ( resid 3 and segid B )
///        ) 2.0 2.0 0.0
/// ```
/// parses into `active = {1, A}`, `partners = [{2, B}, {3, B}]`.
#[derive(Debug, Clone, PartialEq)]
pub struct ParsedRestraint {
    pub active: ResSelection,
    pub partners: Vec<ResSelection>,
}

/// Parses a HADDOCK `.tbl` restraints file into `ParsedRestraint`s.
///
/// # Errors
/// Returns `Err` if an `assign` block doesn't contain at least one
/// `resid ... and segid ...` selection to use as the active residue.
pub fn parse_tbl(content: &str) -> Result<Vec<ParsedRestraint>, String> {
    let cleaned = strip_comments(content);

    // Splitting on "assign" turns the file into one chunk per restraint
    // block (the first chunk, before the first "assign", is discarded).
    let mut restraints = Vec::new();

    for chunk in cleaned.split("assign").skip(1) {
        let mut selections = extract_selections(chunk).into_iter();

        let active = selections
            .next()
            .ok_or_else(|| format!("assign block has no selections: \"assign{}\"", chunk))?;

        restraints.push(ParsedRestraint {
            active,
            partners: selections.collect(),
        });
    }

    Ok(restraints)
}

/// Strips `!`-prefixed CNS comments (rest-of-line) from tbl content.
fn strip_comments(content: &str) -> String {
    content
        .lines()
        .map(|line| line.split('!').next().unwrap_or(""))
        .collect::<Vec<_>>()
        .join("\n")
}

/// Extracts every `resid <N> ... segid <chain>` selection found in `chunk`,
/// in source order. `resid` and `segid` need not be adjacent — any clauses
/// in between (`and name CA`, `and attr z gt 42.00`, ...) are skipped over —
/// and parens don't need surrounding whitespace (`(resid` / `A)` both work).
fn extract_selections(chunk: &str) -> Vec<ResSelection> {
    // Give every paren its own token, regardless of whether the source had
    // whitespace around it, so `(resid` doesn't hide the `resid` keyword.
    let normalized = chunk.replace('(', " ( ").replace(')', " ) ");
    let tokens: Vec<&str> = normalized.split_whitespace().collect();
    let mut selections = Vec::new();

    let mut i = 0;
    while i < tokens.len() {
        if tokens[i] == "resid"
            && let Some(resid) = tokens.get(i + 1).and_then(|t| t.parse::<i16>().ok())
        {
            // Scan ahead for this selection's `segid`, skipping over any
            // intervening clauses. Stop at the next `resid`/`)` so we
            // never borrow the chain from a different selection.
            let mut j = i + 2;
            while j < tokens.len() && tokens[j] != "segid" && tokens[j] != "resid" && tokens[j] != ")" {
                j += 1;
            }
            if tokens.get(j) == Some(&"segid")
                && let Some(chain) = tokens.get(j + 1)
            {
                selections.push(ResSelection {
                    resid,
                    chain: chain.to_string(),
                });
            }
        }
        i += 1;
    }

    selections
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_oneline() {
        let tbl = "assign ( resid 1 and segid A ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n";
        let restraints = parse_tbl(tbl).unwrap();

        assert_eq!(
            restraints,
            vec![ParsedRestraint {
                active: ResSelection {
                    resid: 1,
                    chain: "A".to_string()
                },
                partners: vec![ResSelection {
                    resid: 2,
                    chain: "B".to_string()
                }],
            }]
        );
    }

    #[test]
    fn test_parse_multiline_or() {
        let tbl = "assign ( resid 1 and segid A )\n       (\n        ( resid 2 and segid B )\n     or\n        ( resid 3 and segid B )\n       ) 2.0 2.0 0.0\n\n";
        let restraints = parse_tbl(tbl).unwrap();

        assert_eq!(
            restraints,
            vec![ParsedRestraint {
                active: ResSelection {
                    resid: 1,
                    chain: "A".to_string()
                },
                partners: vec![
                    ResSelection {
                        resid: 2,
                        chain: "B".to_string()
                    },
                    ResSelection {
                        resid: 3,
                        chain: "B".to_string()
                    },
                ],
            }]
        );
    }

    #[test]
    fn test_parse_ignores_atom_clauses() {
        let tbl = "assign ( resid 1 and segid A and name CA ) ( resid 2 and segid B and name CB ) 2.0 2.0 0.0\n\n";
        let restraints = parse_tbl(tbl).unwrap();

        assert_eq!(restraints[0].active.resid, 1);
        assert_eq!(restraints[0].active.chain, "A");
        assert_eq!(restraints[0].partners[0].resid, 2);
        assert_eq!(restraints[0].partners[0].chain, "B");
    }

    #[test]
    fn test_parse_multiple_blocks() {
        let tbl = "assign ( resid 1 and segid A ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n\
                   assign ( resid 5 and segid A ) ( resid 9 and segid B ) 2.0 2.0 0.0\n\n";
        let restraints = parse_tbl(tbl).unwrap();

        assert_eq!(restraints.len(), 2);
        assert_eq!(restraints[1].active.resid, 5);
    }

    #[test]
    fn test_parse_ignores_comments() {
        let tbl = "! this is a comment\nassign ( resid 1 and segid A ) ( resid 2 and segid B ) 2.0 2.0 0.0 ! inline comment\n\n";
        let restraints = parse_tbl(tbl).unwrap();

        assert_eq!(restraints[0].active.resid, 1);
    }

    #[test]
    fn test_parse_atom_clause_before_segid() {
        // `and name CA` sits between `resid` and `segid` here, unlike the
        // generator's own output which always puts `segid` right after
        // `resid`. Real hand-written .tbl files use both orderings.
        let tbl = "assign ( resid 1 and name CA and segid A ) ( resid 2 and segid B ) 2.0 2.0 0.0\n\n";
        let restraints = parse_tbl(tbl).unwrap();

        assert_eq!(restraints[0].active.resid, 1);
        assert_eq!(restraints[0].active.chain, "A");
    }

    #[test]
    fn test_extract_selections_no_surrounding_whitespace() {
        let selections = extract_selections("(resid 1 and segid A)(resid 2 and segid B)");
        assert_eq!(
            selections,
            vec![
                ResSelection {
                    resid: 1,
                    chain: "A".to_string()
                },
                ResSelection {
                    resid: 2,
                    chain: "B".to_string()
                },
            ]
        );
    }

    #[test]
    fn test_extract_selections() {
        let selections = extract_selections("( resid 1 and segid A ) ( resid 2 and segid B )");
        assert_eq!(
            selections,
            vec![
                ResSelection {
                    resid: 1,
                    chain: "A".to_string()
                },
                ResSelection {
                    resid: 2,
                    chain: "B".to_string()
                },
            ]
        );
    }
}
