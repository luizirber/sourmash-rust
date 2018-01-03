extern crate sourmash;

use sourmash::KmerMinHash;

#[test]
fn throws_error() {
    let mut mh = KmerMinHash::new(1, 4, false, 42, 0xffffffffffffffff, false);
    match mh.add_sequence(b"ATGR", false) {
        Ok(_) => assert!(false, "R is not a valid DNA character"),
        Err(_) => assert!(true),
    }
}
