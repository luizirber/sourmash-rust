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

#[test]
fn merge() {
    let mut a = KmerMinHash::new(20, 10, false, 42, 0xffffffffffffffff, false);
    let mut b = KmerMinHash::new(20, 10, false, 42, 0xffffffffffffffff, false);

    a.add_sequence(b"TGCCGCCCAGCA", false);
    b.add_sequence(b"TGCCGCCCAGCA", false);

    a.add_sequence(b"GTCCGCCCAGTGA", false);
    b.add_sequence(b"GTCCGCCCAGTGG", false);

    if let Ok((mins, abunds)) = a.merge(&b) {
    assert_eq!(mins, vec![2996412506971915891, 4448613756639084635, 8373222269469409550, 9390240264282449587, 11085758717695534616, 11668188995231815419, 11760449009842383350, 14682565545778736889]);
    };
}
