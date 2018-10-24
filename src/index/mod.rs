#[macro_use]
pub mod sbt;

#[macro_use]
pub mod storage;

#[macro_use]
pub mod nodegraph;

#[macro_use]
pub mod linear;

#[macro_use]
pub mod search;

use std::path::Path;
use std::rc::Rc;

use failure::Error;

use index::storage::{ReadData, Storage};
use Signature;

pub trait Index {
    type Item;

    fn find<F>(
        &self,
        search_fn: F,
        sig: &Self::Item,
        threshold: f64,
    ) -> Result<Vec<&Self::Item>, Error>
    where
        F: Fn(&Comparable<Self::Item>, &Self::Item, f64) -> bool;

    fn insert(&mut self, node: &Self::Item);

    fn save<P: AsRef<Path>>(&self, path: P) -> Result<(), Error>;

    fn load<P: AsRef<Path>>(path: P) -> Result<(), Error>;
}

// TODO: split into two traits, Similarity and Containment?
pub trait Comparable<O> {
    fn similarity(&self, other: &O) -> f64;
    fn containment(&self, other: &O) -> f64;
}

impl<'a, N, L> Comparable<L> for &'a N
where
    N: Comparable<L>,
{
    fn similarity(&self, other: &L) -> f64 {
        (*self).similarity(&other)
    }

    fn containment(&self, other: &L) -> f64 {
        (*self).containment(&other)
    }
}

#[derive(Deserialize)]
pub struct LeafInfo {
    pub filename: String,
    pub name: String,
    pub metadata: String,
}

#[derive(Builder, Default, Clone)]
pub struct Leaf {
    pub(crate) filename: String,
    pub(crate) name: String,
    pub(crate) metadata: String,

    #[builder(setter(skip))]
    pub(crate) storage: Option<Rc<Storage>>,
}

impl std::fmt::Debug for Leaf {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "Leaf [filename: {}, name: {}, metadata: {}]",
            self.filename, self.name, self.metadata
        )
    }
}

impl<S: Storage + ?Sized> ReadData<Signature, S> for Leaf {
    fn data(&self, storage: &S) -> Result<Signature, Error> {
        // TODO: cache this call!
        let raw = storage.load(&self.filename)?;
        let sigs: Vec<Signature> = serde_json::from_reader(&mut &raw[..])?;
        // TODO: select the right sig?
        Ok(sigs[0].clone())
    }
}

impl Comparable<Leaf> for Leaf {
    fn similarity(&self, other: &Leaf) -> f64 {
        if let Some(storage) = &self.storage {
            let ng: Signature = self.data(&**storage).unwrap();
            let ong: Signature = other.data(&**storage).unwrap();

            // TODO: select the right signatures...
            ng.signatures[0].compare(&ong.signatures[0]).unwrap()
        } else {
            // TODO: in this case storage is not set up,
            // so we should throw an error?
            0.0
        }
    }

    fn containment(&self, other: &Leaf) -> f64 {
        if let Some(storage) = &self.storage {
            let mut ng: Signature = self.data(&**storage).unwrap();
            let ong: Signature = other.data(&**storage).unwrap();

            // TODO: select the right signatures...
            let common = ng.signatures[0].count_common(&ong.signatures[0]).unwrap();
            let size = ng.signatures[0].mins.len();
            common as f64 / size as f64
        } else {
            // TODO: in this case storage is not set up,
            // so we should throw an error?
            0.0
        }
    }
}
