extern crate sourmash;

#[cfg(feature = "from-finch")]
extern crate finch;

use sourmash::KmerMinHash;

#[cfg(feature = "from-finch")]
use finch::minhashes::MinHashKmers;


#[cfg(feature = "from-finch")]
#[test]
fn from_finch() {
    let mut a = KmerMinHash { num: 20, ksize: 10, .. Default::default() };
    let mut b = MinHashKmers::new(20, 42);

    let seq = b"TGCCGCCCAGCACCGGGTGACTAGGTTGAGCCATGATTAACCTGCAATGA";

    a.add_sequence(seq, false);

    for kmer in seq.windows(10) {
        b.push(kmer, 0);
    }

    let b_hashes = b.into_vec();

    println!("{:?}", b_hashes.len());
    println!("{:?}", a.mins.len());

    for item in b_hashes.iter() {
        assert!(a.mins.contains(&(item.hash as u64)));
    }

}
