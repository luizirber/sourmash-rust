#[macro_use]
extern crate clap;
extern crate env_logger;
extern crate exitfailure;
extern crate failure;
extern crate human_panic;
#[macro_use]
extern crate log;

use exitfailure::ExitFailure;
use failure::{Error, ResultExt};
use human_panic::setup_panic;

extern crate sourmash;

use clap::App;

use sourmash::index::nodegraph::Nodegraph;
use sourmash::index::sbt::{scaffold, Node, SBT};
use sourmash::index::{Index, Leaf};
use sourmash::Signature;

fn main() -> Result<(), ExitFailure> {
    setup_panic!();
    env_logger::init();

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
        _ => Ok(()),
    }
}
