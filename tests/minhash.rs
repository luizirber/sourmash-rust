extern crate sourmash;

use sourmash::KmerMinHash;

#[test]
fn throws_error() {
    let mut mh = KmerMinHash {
        num: 1,
        ksize: 4,
        ..Default::default()
    };

    match mh.add_sequence(b"ATGR", false) {
        Ok(_) => assert!(false, "R is not a valid DNA character"),
        Err(_) => assert!(true),
    }
}

#[test]
fn merge() {
    let mut a = KmerMinHash {
        num: 20,
        ksize: 10,
        ..Default::default()
    };
    let mut b = KmerMinHash {
        num: 20,
        ksize: 10,
        ..Default::default()
    };

    a.add_sequence(b"TGCCGCCCAGCA", false);
    b.add_sequence(b"TGCCGCCCAGCA", false);

    a.add_sequence(b"GTCCGCCCAGTGA", false);
    b.add_sequence(b"GTCCGCCCAGTGG", false);

    a.merge(&b);
    assert_eq!(
        a.mins,
        vec![
            2996412506971915891,
            4448613756639084635,
            8373222269469409550,
            9390240264282449587,
            11085758717695534616,
            11668188995231815419,
            11760449009842383350,
            14682565545778736889,
        ]
    );
}

#[test]
fn compare() {
    let mut a = KmerMinHash {
        num: 20,
        ksize: 10,
        ..Default::default()
    };
    let mut b = KmerMinHash {
        num: 20,
        ksize: 10,
        ..Default::default()
    };

    a.add_sequence(b"TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA", false);
    b.add_sequence(b"TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA", false);
    assert_eq!(a.compare(&b).unwrap(), 1.0);
    //    assert_eq!(b.compare(&b).unwrap(), 1.0);
    assert_eq!(b.compare(&a).unwrap(), 1.0);
    //    assert_eq!(a.compare(&a).unwrap(), 1.0);

    b.add_sequence(b"TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA", false);
    assert_eq!(a.compare(&b).unwrap(), 1.0);
    //    assert_eq!(b.compare(&b).unwrap(), 1.0);
    assert_eq!(b.compare(&a).unwrap(), 1.0);
    //    assert_eq!(a.compare(&a).unwrap(), 1.0);

    b.add_sequence(b"GATTGGTGCACACTTAACTGGGTGCCGCGCTGGTGCTGATCCATGAAGTT", false);
    assert!(a.compare(&b).unwrap() >= 0.3);
    assert!(b.compare(&a).unwrap() >= 0.3);
}
