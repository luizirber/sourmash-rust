#[macro_use]
extern crate criterion;
extern crate rand;
extern crate sourmash;

use std::path::PathBuf;

use criterion::{Bencher, Criterion, Fun};
use sourmash::index::linear::LinearIndexBuilder;
use sourmash::index::sbt::{Node, SBT};
use sourmash::index::search::search_minhashes;
use sourmash::index::{Index, Leaf};

fn find_bench(c: &mut Criterion) {
    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("tests/data/v5.sbt.json");

    let sbt: SBT<Node, Leaf> = SBT::from_path(filename).expect("Loading error");

    let leaf: Leaf = (*sbt.leaves().first().unwrap()).clone();

    let mut linear = LinearIndexBuilder::default()
        .storage(sbt.storage())
        .build()
        .unwrap();
    for l in &sbt.leaves() {
        linear.insert(*l);
    }

    let sbt_find = Fun::new("sbt_find", move |b: &mut Bencher, leaf: &Leaf| {
        b.iter(|| sbt.find(search_minhashes, leaf, 0.1))
    });

    let linear_find = Fun::new("linear_find", move |b: &mut Bencher, leaf: &Leaf| {
        b.iter(|| linear.find(search_minhashes, leaf, 0.1))
    });

    let functions = vec![sbt_find, linear_find];
    c.bench_functions("find", functions, leaf);
}

criterion_group!(benches, find_bench);
criterion_main!(benches);
