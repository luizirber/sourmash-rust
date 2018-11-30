use std::fs::File;
use std::io;

use clap::{load_yaml, App};
use exitfailure::ExitFailure;
use failure::{Error, ResultExt};
use human_panic::setup_panic;
use log::{debug, info, LevelFilter};

use sourmash::index::nodegraph::Nodegraph;
use sourmash::index::sbt::{scaffold, Node, SBT};
use sourmash::index::{Index, Leaf};
use sourmash::Signature;

struct Query<T> {
    data: T,
}

impl Query<Signature> {
    fn ksize(&self) -> u64 {
        // TODO: this might panic
        self.data.signatures[0].ksize as u64
    }

    fn moltype(&self) -> String {
        // TODO: this might panic
        if self.data.signatures[0].is_protein {
            "protein".into()
        } else {
            "DNA".into()
        }
    }

    fn name(&self) -> String {
        self.data.name.clone().unwrap()
    }
}

fn load_query_signature(
    query: &str,
    ksize: usize,
    moltype: Option<&str>,
    scaled: Option<u64>,
) -> Result<Query<Signature>, Error> {
    let mut reader = io::BufReader::new(File::open(query)?);
    let sigs = Signature::load_signatures(&mut reader, ksize, moltype, scaled)?;

    debug!("{:?}", sigs);
    // TODO: what if we have more than one left?
    let data = sigs[0].clone();

    Ok(Query { data })
}

struct Database {
    data: String,
}

fn load_sbts_and_sigs(
    databases: Vec<&str>,
    query: &Query<Signature>,
    containment: bool,
    traverse_directory: bool,
) -> Result<Vec<Database>, Error> {
    let dbs = Vec::default();

    let ksize = query.ksize();
    let moltype = query.moltype();

    Ok(dbs)
}

fn main() -> Result<(), ExitFailure> {
    //setup_panic!();

    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let yml = load_yaml!("sourmash.yml");
    let m = App::from_yaml(yml).get_matches();

    match m.subcommand_name() {
        Some("scaffold") => {
            let cmd = m.subcommand_matches("scaffold").unwrap();
            let sbt_file = cmd.value_of("current_sbt").unwrap();

            let sbt: SBT<Node<Nodegraph>, Leaf<Signature>> = SBT::from_path(sbt_file)?;
            let new_sbt: SBT<Node<Nodegraph>, Leaf<Signature>> = scaffold(sbt.leaves());

            assert_eq!(new_sbt.leaves().len(), 100);
            Ok(())
        }
        Some("search") => {
            let cmd = m.subcommand_matches("search").unwrap();

            if cmd.is_present("quiet") {
                log::set_max_level(LevelFilter::Warn);
            }

            let query = load_query_signature(
                cmd.value_of("query").unwrap(),
                if cmd.is_present("ksize") {
                    cmd.value_of("ksize").unwrap().parse().unwrap()
                } else {
                    0
                },
                Some("dna"), // TODO: select moltype,
                if cmd.is_present("scaled") {
                    Some(cmd.value_of("scaled").unwrap().parse().unwrap())
                } else {
                    None
                },
            )?;

            let databases = load_sbts_and_sigs(
                cmd.values_of("databases")
                    .map(|vals| vals.collect::<Vec<_>>())
                    .unwrap(),
                &query,
                !cmd.is_present("containment"),
                cmd.is_present("traverse-directory"),
            );

            info!(
                "loaded query: {}... (k={}, {})",
                query.name(),
                query.ksize(),
                query.moltype()
            );

            /*
            info!("Shouldn't show up with -q!");
            debug!("{:?}", cmd);
            */

            Ok(())
        }
        _ => {
            println!("{:?}", m);
            Ok(())
        }
    }
}
