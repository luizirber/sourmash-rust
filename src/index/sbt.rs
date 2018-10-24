use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};
use std::rc::Rc;

use failure::Error;
use serde::de::{Deserialize, Deserializer};
use serde::ser::{Serialize, SerializeStruct, Serializer};

use index::nodegraph::Nodegraph;
use index::storage::{FSStorage, ReadData, Storage, StorageInfo};
use index::{Comparable, Index, Leaf, LeafInfo};
use Signature;

#[derive(Builder)]
pub struct SBT<N, L> {
    #[builder(default = "2")]
    d: u32,

    storage: Rc<Storage>,
    factory: Factory,

    #[builder(setter(skip))]
    nodes: HashMap<u64, N>,

    #[builder(setter(skip))]
    leaves: HashMap<u64, L>,
}

impl<N, L> SBT<N, L> {
    #[inline(always)]
    fn parent(&self, pos: u64) -> Option<u64> {
        if pos == 0 {
            None
        } else {
            Some((pos - 1) / (u64::from(self.d)))
        }
    }

    #[inline(always)]
    fn child(&self, parent: u64, pos: u64) -> u64 {
        u64::from(self.d) * parent + pos + 1
    }

    #[inline(always)]
    fn children(&self, pos: u64) -> Vec<u64> {
        (0..u64::from(self.d)).map(|c| self.child(pos, c)).collect()
    }

    pub fn leaves(&self) -> Vec<&L> {
        self.leaves.values().collect()
    }

    pub fn storage(&self) -> Rc<Storage> {
        Rc::clone(&self.storage)
    }

    // combine
}

impl SBT<Node, Leaf> {
    pub fn from_reader<R, P>(rdr: &mut R, path: P) -> Result<SBT<Node, Leaf>, Error>
    where
        R: Read,
        P: AsRef<Path>,
    {
        let sbt: SBTInfo<NodeInfo, LeafInfo> = serde_json::from_reader(rdr)?;

        // TODO: match with available Storage while we don't
        // add a function to build a Storage from a StorageInfo
        let mut basepath = PathBuf::new();
        basepath.push(path);
        basepath.push(&sbt.storage.args["path"]);

        let storage: Rc<Storage> = Rc::new(FSStorage { basepath });

        Ok(SBT {
            d: sbt.d,
            factory: sbt.factory,
            storage: Rc::clone(&storage),
            nodes: sbt
                .nodes
                .into_iter()
                .map(|(n, l)| {
                    let new_node = Node {
                        filename: l.filename,
                        name: l.name,
                        metadata: l.metadata,
                        storage: Some(Rc::clone(&storage)),
                    };
                    (n, new_node)
                }).collect(),
            leaves: sbt
                .leaves
                .into_iter()
                .map(|(n, l)| {
                    let new_node = Leaf {
                        filename: l.filename,
                        name: l.name,
                        metadata: l.metadata,
                        storage: Some(Rc::clone(&storage)),
                    };
                    (n, new_node)
                }).collect(),
        })
    }

    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<SBT<Node, Leaf>, Error> {
        let file = File::open(&path)?;
        let mut reader = BufReader::new(file);

        // TODO: match with available Storage while we don't
        // add a function to build a Storage from a StorageInfo
        let mut basepath = PathBuf::new();
        basepath.push(path);
        basepath.canonicalize()?;

        let sbt = SBT::<Node, Leaf>::from_reader(&mut reader, &basepath.parent().unwrap())?;
        Ok(sbt)
    }
}

impl<N, L> Index for SBT<N, L>
where
    N: Comparable<N> + Comparable<L>,
    L: Comparable<L>,
{
    type Item = L;

    fn find<F>(&self, search_fn: F, sig: &L, threshold: f64) -> Result<Vec<&L>, Error>
    where
        F: Fn(&Comparable<Self::Item>, &Self::Item, f64) -> bool,
    {
        let mut matches = Vec::new();
        let mut visited = HashSet::new();
        let mut queue = vec![0u64];

        while !queue.is_empty() {
            let pos = queue.pop().unwrap();
            if !visited.contains(&pos) {
                visited.insert(pos);

                if let Some(node) = self.nodes.get(&pos) {
                    if search_fn(&node, sig, threshold) {
                        for c in self.children(pos) {
                            queue.push(c);
                        }
                    }
                } else if let Some(leaf) = self.leaves.get(&pos) {
                    if search_fn(leaf, sig, threshold) {
                        matches.push(leaf);
                    }
                }
            }
        }

        Ok(matches)
    }

    fn insert(&mut self, node: &L) {}

    fn save<P: AsRef<Path>>(&self, path: P) -> Result<(), Error> {
        Ok(())
    }

    fn load<P: AsRef<Path>>(path: P) -> Result<(), Error> {
        Ok(())
    }
}

#[derive(Builder, Clone, Default, Deserialize)]
pub struct Factory {
    class: String,
    args: Vec<u64>,
}

pub struct Node {
    filename: String,
    name: String,
    metadata: HashMap<String, u64>,
    storage: Option<Rc<Storage>>,
}

impl Comparable<Node> for Node {
    fn similarity(&self, other: &Node) -> f64 {
        if let Some(storage) = &self.storage {
            let ng: Nodegraph = self.data(&**storage).unwrap();
            let ong: Nodegraph = other.data(&**storage).unwrap();
            ng.similarity(&ong)
        } else {
            // TODO: in this case storage is not set up,
            // so we should throw an error?
            0.0
        }
    }

    fn containment(&self, other: &Node) -> f64 {
        if let Some(storage) = &self.storage {
            let ng: Nodegraph = self.data(&**storage).unwrap();
            let ong: Nodegraph = other.data(&**storage).unwrap();
            ng.containment(&ong)
        } else {
            // TODO: in this case storage is not set up,
            // so we should throw an error?
            0.0
        }
    }
}

impl Comparable<Leaf> for Node {
    fn similarity(&self, other: &Leaf) -> f64 {
        if let Some(storage) = &self.storage {
            let ng: Nodegraph = self.data(&**storage).unwrap();
            let oth: Signature = other.data(&**storage).unwrap();

            // TODO: select the right signatures...
            let sig = &oth.signatures[0];
            if sig.size() == 0 {
                return 0.0;
            }

            let matches: usize = sig.mins.iter().map(|h| ng.get(*h)).sum();

            let min_n_below = self.metadata["min_n_below"] as f64;

            // This overestimates the similarity, but better than truncating too
            // soon and losing matches
            matches as f64 / min_n_below
        } else {
            // TODO: throw error, storage not initialized
            0.0
        }
    }

    fn containment(&self, other: &Leaf) -> f64 {
        if let Some(storage) = &self.storage {
            let ng: Nodegraph = self.data(&**storage).unwrap();
            let oth: Signature = other.data(&**storage).unwrap();

            // TODO: select the right signatures...
            let sig = &oth.signatures[0];
            if sig.size() == 0 {
                return 0.0;
            }

            let matches: usize = sig.mins.iter().map(|h| ng.get(*h)).sum();

            matches as f64 / sig.size() as f64
        } else {
            // TODO: throw error, storage not initialized
            0.0
        }
    }
}

impl<S: Storage + ?Sized> ReadData<Nodegraph, S> for Node {
    fn data(&self, storage: &S) -> Result<Nodegraph, Error> {
        // TODO: cache this call!
        // probably use lazy-init with a field in the struct?
        // https://docs.rs/lazy-init/0.3.0/lazy_init/
        // or lazy_cell
        let raw = storage.load(&self.filename)?;
        Nodegraph::from_reader(&mut &raw[..])
    }
}

#[derive(Deserialize)]
struct NodeInfo {
    filename: String,
    name: String,
    metadata: HashMap<String, u64>,
}

#[derive(Deserialize)]
struct SBTInfo<N, L> {
    d: u32,
    version: u32,
    storage: StorageInfo,
    factory: Factory,
    nodes: HashMap<u64, N>,
    leaves: HashMap<u64, L>,
}

#[cfg(test)]
mod test {
    use super::*;
    use index::linear::{LinearIndex, LinearIndexBuilder};
    use index::search::{search_minhashes, search_minhashes_containment};

    #[test]
    fn load_sbt() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("tests/data/v5.sbt.json");

        let sbt = SBT::<Node, Leaf>::from_path(filename).expect("Loading error");

        assert_eq!(sbt.d, 2);
        //assert_eq!(sbt.storage.backend, "FSStorage");
        //assert_eq!(sbt.storage.args["path"], ".sbt.v5");
        assert_eq!(sbt.factory.class, "GraphFactory");
        assert_eq!(sbt.factory.args, [1, 100000, 4]);

        println!("sbt leaves {:?} {:?}", sbt.leaves.len(), sbt.leaves);

        let leaf = &sbt.leaves[&7];

        let results = sbt.find(search_minhashes, &leaf, 0.5).unwrap();
        assert_eq!(results.len(), 1);
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);

        let results = sbt.find(search_minhashes, &leaf, 0.1).unwrap();
        assert_eq!(results.len(), 2);
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);

        let mut linear = LinearIndexBuilder::default()
            .storage(Rc::clone(&sbt.storage) as Rc<Storage>)
            .build()
            .unwrap();
        for l in &sbt.leaves {
            linear.insert(l.1);
        }

        println!(
            "linear leaves {:?} {:?}",
            linear.leaves.len(),
            linear.leaves
        );

        let results = linear.find(search_minhashes, &leaf, 0.5).unwrap();
        assert_eq!(results.len(), 1);
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);

        let results = linear.find(search_minhashes, &leaf, 0.1).unwrap();
        assert_eq!(results.len(), 2);
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);

        let results = linear
            .find(search_minhashes_containment, &leaf, 0.5)
            .unwrap();
        assert_eq!(results.len(), 2);
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);

        let results = linear
            .find(search_minhashes_containment, &leaf, 0.1)
            .unwrap();
        assert_eq!(results.len(), 4);
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);
    }
}
