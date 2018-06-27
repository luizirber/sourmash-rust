extern crate serde_json;
extern crate sourmash;

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use sourmash::{Signature, KmerMinHash};

#[test]
fn load_signature() {
    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("tests/data/genome-s10+s11.sig");

    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    let sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

    assert_eq!(sigs.len(), 1);

    let sig = sigs.get(0).unwrap();
    assert_eq!(sig.class, "sourmash_signature");
    assert_eq!(sig.email, "");
    assert_eq!(sig.filename, "-");
    assert_eq!(sig.hash_function, "0.murmur64");
    assert_eq!(sig.name, "s10+s11");
    assert_eq!(sig.signatures.len(), 4);
}
