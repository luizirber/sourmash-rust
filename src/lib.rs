extern crate backtrace;
extern crate byteorder;
extern crate fixedbitset;
extern crate md5;
extern crate murmurhash3;

#[macro_use]
extern crate serde_derive;
extern crate serde;
extern crate serde_json;

#[cfg(feature = "from-finch")]
extern crate finch;

#[macro_use]
extern crate lazy_static;

#[cfg(test)]
#[macro_use]
extern crate proptest;

#[macro_use]
pub mod errors;

#[macro_use]
pub mod utils;

#[macro_use]
pub mod ffi;

#[macro_use]
pub mod nodegraph;

#[cfg(feature = "from-finch")]
pub mod from;

use serde::de::{Deserialize, Deserializer};
use serde::ser::{Serialize, SerializeStruct, Serializer};

use std::collections::HashMap;
use std::collections::HashSet;
use std::hash::{BuildHasherDefault, Hasher};
use std::iter::FromIterator;
use std::str;

use murmurhash3::murmurhash3_x64_128;

use errors::{ErrorKind, Result};

pub fn _hash_murmur(kmer: &[u8], seed: u64) -> u64 {
    murmurhash3_x64_128(kmer, seed).0
}

// This comes from finch
pub struct NoHashHasher(u64);

impl Default for NoHashHasher {
    #[inline]
    fn default() -> NoHashHasher {
        NoHashHasher(0x0)
    }
}

impl Hasher for NoHashHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        *self = NoHashHasher(
            (u64::from(bytes[0]) << 24)
                + (u64::from(bytes[1]) << 16)
                + (u64::from(bytes[2]) << 8)
                + u64::from(bytes[3]),
        );
    }
    fn finish(&self) -> u64 {
        self.0
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct KmerMinHash {
    pub num: u32,
    pub ksize: u32,
    pub is_protein: bool,
    pub seed: u64,
    pub max_hash: u64,
    pub mins: Vec<u64>,
    pub abunds: Option<Vec<u64>>,
}

impl Default for KmerMinHash {
    fn default() -> KmerMinHash {
        KmerMinHash {
            num: 1000,
            ksize: 21,
            is_protein: false,
            seed: 42,
            max_hash: 0,
            mins: Vec::with_capacity(1000),
            abunds: None,
        }
    }
}

impl Serialize for KmerMinHash {
    fn serialize<S>(&self, serializer: S) -> ::std::result::Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let n_fields = match &self.abunds {
            Some(_) => 8,
            _ => 7,
        };

        let mut md5_ctx = md5::Context::new();
        md5_ctx.consume(&self.ksize.to_string());
        self.mins
            .iter()
            .map(|x| md5_ctx.consume(x.to_string()))
            .count();

        let mut partial = serializer.serialize_struct("KmerMinHash", n_fields)?;
        partial.serialize_field("num", &self.num)?;
        partial.serialize_field("ksize", &self.ksize)?;
        partial.serialize_field("seed", &self.seed)?;
        partial.serialize_field("max_hash", &self.max_hash)?;
        partial.serialize_field("mins", &self.mins)?;

        partial.serialize_field("md5sum", &format!("{:x}", md5_ctx.compute()))?;

        if let Some(abunds) = &self.abunds {
            partial.serialize_field("abundances", abunds)?;
        }

        partial.serialize_field(
            "molecule",
            match &self.is_protein {
                true => "protein",
                false => "DNA",
            },
        )?;

        partial.end()
    }
}

impl<'de> Deserialize<'de> for KmerMinHash {
    fn deserialize<D>(deserializer: D) -> ::std::result::Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct TempSig {
            num: u32,
            ksize: u32,
            seed: u64,
            max_hash: u64,
            md5sum: String,
            mins: Vec<u64>,
            abundances: Option<Vec<u64>>,
            molecule: String,
        }

        let tmpsig = TempSig::deserialize(deserializer)?;

        let num = if tmpsig.max_hash != 0 { 0 } else { tmpsig.num };

        Ok(KmerMinHash {
            num,
            ksize: tmpsig.ksize,
            seed: tmpsig.seed,
            max_hash: tmpsig.max_hash,
            mins: tmpsig.mins,
            abunds: tmpsig.abundances,
            is_protein: match tmpsig.molecule.as_ref() {
                "protein" => true,
                "DNA" => false,
                _ => false, // TODO: throw error
            },
        })
    }
}

impl KmerMinHash {
    pub fn new(
        num: u32,
        ksize: u32,
        is_protein: bool,
        seed: u64,
        max_hash: u64,
        track_abundance: bool,
    ) -> KmerMinHash {
        let mins: Vec<u64>;
        let abunds: Option<Vec<u64>>;

        if num > 0 {
            mins = Vec::with_capacity(num as usize);
        } else {
            mins = Vec::with_capacity(1000);
        }

        if track_abundance {
            abunds = Some(Vec::with_capacity(mins.capacity()));
        } else {
            abunds = None
        }

        KmerMinHash {
            num,
            ksize,
            is_protein,
            seed,
            max_hash,
            mins,
            abunds,
        }
    }

    pub fn check_compatible(&mut self, other: &KmerMinHash) -> Result<bool> {
        if self.ksize != other.ksize {
            return Err(ErrorKind::MismatchKSizes.into());
        }
        if self.is_protein != other.is_protein {
            return Err(ErrorKind::MismatchDNAProt.into());
        }
        if self.max_hash != other.max_hash {
            return Err(ErrorKind::MismatchMaxHash.into());
        }
        if self.seed != other.seed {
            return Err(ErrorKind::MismatchSeed.into());
        }
        Ok(true)
    }

    pub fn add_hash(&mut self, hash: u64) {
        let current_max = match self.mins.last() {
            Some(&x) => x,
            None => u64::max_value(),
        };

        if hash <= self.max_hash || self.max_hash == 0 {
            // empty? add it, if within range / no range specified.
            if self.mins.is_empty() {
                self.mins.push(hash);
                if let Some(ref mut abunds) = self.abunds {
                    abunds.push(1);
                }
                return;
            } else if hash <= self.max_hash
                || current_max > hash
                || (self.mins.len() as u32) < self.num
            {
                // "good" hash - within range, smaller than current entry, or
                // still have space available
                let pos = match self.mins.binary_search(&hash) {
                    Ok(p) => p,
                    Err(p) => p,
                };

                if pos == self.mins.len() {
                    // at end - must still be growing, we know the list won't
                    // get too long
                    self.mins.push(hash);
                    if let Some(ref mut abunds) = self.abunds {
                        abunds.push(1);
                    }
                } else if self.mins[pos] != hash {
                    // didn't find hash in mins, so inserting somewhere
                    // in the middle; shrink list if needed.
                    self.mins.insert(pos, hash);
                    if let Some(ref mut abunds) = self.abunds {
                        abunds.insert(pos, 1);
                    }

                    // is it too big now?
                    if self.num != 0 && self.mins.len() > (self.num as usize) {
                        self.mins.pop();
                        if let Some(ref mut abunds) = self.abunds {
                            abunds.pop();
                        }
                    }
                } else if let Some(ref mut abunds) = self.abunds {
                    // pos == hash: hash value already in mins, inc count
                    abunds[pos] += 1;
                }
            }
        }
    }

    pub fn add_word(&mut self, word: &[u8]) {
        let hash = _hash_murmur(word, self.seed);
        self.add_hash(hash);
    }

    pub fn add_sequence(&mut self, seq: &[u8], force: bool) -> Result<()> {
        let sequence: Vec<u8> = seq
            .iter()
            .map(|&x| (x as char).to_ascii_uppercase() as u8)
            .collect();
        if sequence.len() >= (self.ksize as usize) {
            if !self.is_protein {
                // dna
                for kmer in sequence.windows(self.ksize as usize) {
                    if _checkdna(kmer) {
                        let rc = revcomp(kmer);
                        if kmer < &rc {
                            self.add_word(kmer);
                        } else {
                            self.add_word(&rc);
                        }
                    } else if !force {
                        return Err(ErrorKind::InvalidDNA(
                            String::from_utf8(kmer.to_vec()).unwrap(),
                        ).into());
                    }
                }
            } else {
                // protein
                let rc = revcomp(&sequence);
                let aa_ksize = self.ksize / 3;

                for i in 0..3 {
                    let substr: Vec<u8> = sequence
                        .iter()
                        .cloned()
                        .skip(i)
                        .take(sequence.len() - i)
                        .collect();
                    let aa = to_aa(&substr);

                    aa.windows(aa_ksize as usize)
                        .map(|n| self.add_word(n))
                        .count();

                    let rc_substr: Vec<u8> =
                        rc.iter().cloned().skip(i).take(rc.len() - i).collect();
                    let aa_rc = to_aa(&rc_substr);

                    aa_rc
                        .windows(aa_ksize as usize)
                        .map(|n| self.add_word(n))
                        .count();
                }
            }
        }
        Ok(())
    }

    pub fn merge(&mut self, other: &KmerMinHash) -> Result<()> {
        self.check_compatible(other)?;
        let max_size = self.mins.len() + other.mins.len();
        let mut merged: Vec<u64> = Vec::with_capacity(max_size);
        let mut merged_abunds: Vec<u64> = Vec::with_capacity(max_size);

        {
            let mut self_iter = self.mins.iter();
            let mut other_iter = other.mins.iter();

            let mut self_abunds_iter: Option<std::slice::Iter<u64>>;
            if let Some(ref mut abunds) = self.abunds {
                self_abunds_iter = Some(abunds.iter());
            } else {
                self_abunds_iter = None;
            }

            let mut other_abunds_iter: Option<std::slice::Iter<u64>>;
            if let Some(ref abunds) = other.abunds {
                other_abunds_iter = Some(abunds.iter());
            } else {
                other_abunds_iter = None;
            }

            let mut self_value = self_iter.next();
            let mut other_value = other_iter.next();
            while self_value.is_some() {
                let value = self_value.unwrap();
                match other_value {
                    None => {
                        merged.push(*value);
                        merged.extend(self_iter);
                        if let Some(sai) = self_abunds_iter {
                            merged_abunds.extend(sai);
                        }
                        break;
                    }
                    Some(x) if x < value => {
                        merged.push(*x);
                        other_value = other_iter.next();

                        if let Some(ref mut oai) = other_abunds_iter {
                            if let Some(v) = oai.next() {
                                merged_abunds.push(*v)
                            }
                        }
                    }
                    Some(x) if x == value => {
                        merged.push(*x);
                        other_value = other_iter.next();
                        self_value = self_iter.next();

                        if let Some(ref mut oai) = other_abunds_iter {
                            if let Some(v) = oai.next() {
                                if let Some(ref mut sai) = self_abunds_iter {
                                    if let Some(s) = sai.next() {
                                        merged_abunds.push(*v + *s)
                                    }
                                }
                            }
                        }
                    }
                    Some(x) if x > value => {
                        merged.push(*value);
                        self_value = self_iter.next();

                        if let Some(ref mut sai) = self_abunds_iter {
                            if let Some(v) = sai.next() {
                                merged_abunds.push(*v)
                            }
                        }
                    }
                    Some(_) => {}
                }
            }
            if let Some(value) = other_value {
                merged.push(*value);
            }
            merged.extend(other_iter);
            if let Some(oai) = other_abunds_iter {
                merged_abunds.extend(oai);
            }
        }

        if merged.len() < (self.num as usize) || (self.num as usize) == 0 {
            self.mins = merged;
            self.abunds = Some(merged_abunds);
        } else {
            self.mins = merged
                .iter()
                .map(|&x| x as u64)
                .take(self.num as usize)
                .collect();
            self.abunds = Some(merged_abunds) // TODO: reduce this one too
        }
        Ok(())
    }

    pub fn add_from(&mut self, other: &KmerMinHash) -> Result<()> {
        for min in &other.mins {
            self.add_hash(*min);
        }
        Ok(())
    }

    pub fn add_many(&mut self, hashes: &[u64]) -> Result<()> {
        for min in hashes {
            self.add_hash(*min);
        }
        Ok(())
    }

    pub fn add_many_with_abund(&mut self, hashes: &[(u64, u64)]) -> Result<()> {
        for item in hashes {
            for _i in 0..item.1 {
                self.add_hash(item.0);
            }
        }
        Ok(())
    }

    pub fn count_common(&mut self, other: &KmerMinHash) -> Result<u64> {
        self.check_compatible(other)?;
        let s1: HashSet<&u64, BuildHasherDefault<NoHashHasher>> =
            HashSet::from_iter(self.mins.iter());
        let s2: HashSet<&u64, BuildHasherDefault<NoHashHasher>> =
            HashSet::from_iter(other.mins.iter());
        Ok(s1.intersection(&s2).count() as u64)
    }

    pub fn intersection(&mut self, other: &KmerMinHash) -> Result<(Vec<u64>, u64)> {
        self.check_compatible(other)?;

        let mut combined_mh = KmerMinHash::new(
            self.num,
            self.ksize,
            self.is_protein,
            self.seed,
            self.max_hash,
            self.abunds.is_some(),
        );

        combined_mh.merge(&self)?;
        combined_mh.merge(&other)?;

        let s1: HashSet<_, BuildHasherDefault<NoHashHasher>> = HashSet::from_iter(self.mins.iter());
        let s2: HashSet<_, BuildHasherDefault<NoHashHasher>> =
            HashSet::from_iter(other.mins.iter());
        let s3: HashSet<_, BuildHasherDefault<NoHashHasher>> =
            HashSet::from_iter(combined_mh.mins.iter());

        let i1 = &s1 & &s2;
        let i2 = &i1 & &s3;

        let common: Vec<u64> = i2.into_iter().cloned().collect();
        Ok((common, combined_mh.mins.len() as u64))
    }

    pub fn intersection_size(&mut self, other: &KmerMinHash) -> Result<(u64, u64)> {
        self.check_compatible(other)?;

        let mut combined_mh = KmerMinHash::new(
            self.num,
            self.ksize,
            self.is_protein,
            self.seed,
            self.max_hash,
            self.abunds.is_some(),
        );

        combined_mh.merge(&self)?;
        combined_mh.merge(&other)?;

        let s1: HashSet<_, BuildHasherDefault<NoHashHasher>> = HashSet::from_iter(self.mins.iter());
        let s2: HashSet<_, BuildHasherDefault<NoHashHasher>> =
            HashSet::from_iter(other.mins.iter());
        let s3: HashSet<_, BuildHasherDefault<NoHashHasher>> =
            HashSet::from_iter(combined_mh.mins.iter());

        let i1 = &s1 & &s2;
        let i2 = &i1 & &s3;

        Ok((i2.into_iter().count() as u64, combined_mh.mins.len() as u64))
    }

    pub fn compare(&mut self, other: &KmerMinHash) -> Result<f64> {
        self.check_compatible(other)?;
        if let Ok((common, size)) = self.intersection_size(other) {
            return Ok(common as f64 / u64::max(1, size) as f64);
        } else {
            return Ok(0.0);
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Signature {
    #[serde(default = "default_class")]
    pub class: String,

    #[serde(default)]
    pub email: String,
    pub hash_function: String,

    pub filename: Option<String>,
    pub name: Option<String>,

    #[serde(default = "default_license")]
    pub license: String,

    pub signatures: Vec<KmerMinHash>,

    #[serde(default = "default_version")]
    pub version: f64,
}

fn default_license() -> String {
    "CC0".to_string()
}

fn default_class() -> String {
    "sourmash_signature".to_string()
}

fn default_version() -> f64 {
    0.4
}

impl Default for Signature {
    fn default() -> Signature {
        Signature {
            class: default_class(),
            email: "".to_string(),
            hash_function: "0.murmur64".to_string(),
            license: default_license(),
            filename: None,
            name: None,
            signatures: Vec::<KmerMinHash>::new(),
            version: default_version(),
        }
    }
}

impl PartialEq for Signature {
    fn eq(&self, other: &Signature) -> bool {
        let metadata = self.class == other.class
            && self.email == other.email
            && self.hash_function == other.hash_function
            && self.filename == other.filename
            && self.name == other.name;

        let mh = &self.signatures[0];
        let other_mh = &other.signatures[0];
        metadata && (mh == other_mh)
    }
}

#[inline]
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|n| match *n as char {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            x => x,
        } as u8) // TODO: error?
        .collect()
}

lazy_static! {
    static ref CODONTABLE: HashMap<&'static str, u8> = {
        let mut m = HashMap::new();

        m.insert("TTT", b'F');
        m.insert("TTC", b'F');
        m.insert("TTA", b'L');
        m.insert("TTG", b'L');

        m.insert("TCT", b'S');
        m.insert("TCC", b'S');
        m.insert("TCA", b'S');
        m.insert("TCG", b'S');

        m.insert("TAT", b'Y');
        m.insert("TAC", b'Y');
        m.insert("TAA", b'*');
        m.insert("TAG", b'*');

        m.insert("TGT", b'C');
        m.insert("TGC", b'C');
        m.insert("TGA", b'*');
        m.insert("TGG", b'W');

        m.insert("CTT", b'L');
        m.insert("CTC", b'L');
        m.insert("CTA", b'L');
        m.insert("CTG", b'L');

        m.insert("CCT", b'P');
        m.insert("CCC", b'P');
        m.insert("CCA", b'P');
        m.insert("CCG", b'P');

        m.insert("CAT", b'H');
        m.insert("CAC", b'H');
        m.insert("CAA", b'Q');
        m.insert("CAG", b'Q');

        m.insert("CGT", b'R');
        m.insert("CGC", b'R');
        m.insert("CGA", b'R');
        m.insert("CGG", b'R');

        m.insert("ATT", b'I');
        m.insert("ATC", b'I');
        m.insert("ATA", b'I');
        m.insert("ATG", b'M');

        m.insert("ACT", b'T');
        m.insert("ACC", b'T');
        m.insert("ACA", b'T');
        m.insert("ACG", b'T');

        m.insert("AAT", b'N');
        m.insert("AAC", b'N');
        m.insert("AAA", b'K');
        m.insert("AAG", b'K');

        m.insert("AGT", b'S');
        m.insert("AGC", b'S');
        m.insert("AGA", b'R');
        m.insert("AGG", b'R');

        m.insert("GTT", b'V');
        m.insert("GTC", b'V');
        m.insert("GTA", b'V');
        m.insert("GTG", b'V');

        m.insert("GCT", b'A');
        m.insert("GCC", b'A');
        m.insert("GCA", b'A');
        m.insert("GCG", b'A');

        m.insert("GAT", b'D');
        m.insert("GAC", b'D');
        m.insert("GAA", b'E');
        m.insert("GAG", b'E');

        m.insert("GGT", b'G');
        m.insert("GGC", b'G');
        m.insert("GGA", b'G');
        m.insert("GGG", b'G');

        m
    };
}

#[inline]
fn to_aa(seq: &[u8]) -> Vec<u8> {
    let mut converted: Vec<u8> = Vec::with_capacity(seq.len() / 3);

    for chunk in seq.chunks(3) {
        if chunk.len() != 3 {
            break;
        }
        if let Some(codon) = CODONTABLE.get(str::from_utf8(chunk).unwrap()) {
            converted.push(*codon);
        }
    }

    converted
}

#[inline]
fn _checkdna(seq: &[u8]) -> bool {
    for n in seq.iter() {
        match *n as char {
            'A' | 'a' | 'C' | 'c' | 'G' | 'g' | 'T' | 't' => (),
            _ => return false,
        }
    }
    true
}
