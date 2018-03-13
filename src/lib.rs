extern crate backtrace;
extern crate murmurhash3;
extern crate ordslice;

#[macro_use]
extern crate error_chain;

#[macro_use]
extern crate lazy_static;

#[macro_use]
pub mod errors;

#[macro_use]
pub mod utils;

#[macro_use]
pub mod ffi;

use std::collections::HashMap;
use std::collections::HashSet;
use std::iter::FromIterator;
use std::str;

use murmurhash3::murmurhash3_x64_128;
use ordslice::Ext;

use errors::{ErrorKind, Result};

pub fn _hash_murmur(kmer: &[u8], seed: u64) -> u64 {
    murmurhash3_x64_128(kmer, seed).0
}

#[derive(Debug)]
pub struct KmerMinHash {
    pub num: u32,
    pub ksize: u32,
    pub is_protein: bool,
    pub seed: u64,
    pub max_hash: u64,
    pub mins: Vec<u64>,
    pub abunds: Option<Vec<u64>>,
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

        if (self.max_hash != 0 && hash <= self.max_hash) || self.max_hash == 0 {
            // empty? add it, if within range / no range specified.
            if self.mins.is_empty() {
                self.mins.push(hash);
                if let Some(ref mut abunds) = self.abunds {
                    abunds.push(1);
                }
                return;
            } else if hash <= self.max_hash || current_max > hash
                || (self.mins.len() as u32) < self.num
            {
                // "good" hash - within range, smaller than current entry, or
                // still have space available
                let pos = self.mins.lower_bound(&hash);

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
                } else { // pos == hash: hash value already in mins, inc count
                    if let Some(ref mut abunds) = self.abunds {
                        abunds[pos] += 1;
                    }
                }
            }
        }
    }

    pub fn add_word(&mut self, word: &[u8]) {
        let hash = _hash_murmur(word, self.seed);
        self.add_hash(hash);
    }

    pub fn add_sequence(&mut self, sequence: &[u8], force: bool) -> Result<()> {
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
                    } else {
                        if !force {
                            let last = vec![kmer.last().cloned().unwrap()];
                            return Err(ErrorKind::InvalidDNA(String::from_utf8(last).unwrap())
                                .into());
                        }
                    }
                }
            } else {
                // protein
                let rc = revcomp(sequence);
                let aa_ksize = self.ksize / 3;

                for i in 0..3 {
                    let substr: Vec<u8> = sequence.iter()
                                         .cloned()
                                         .skip(i)
                                         .take(sequence.len() - i)
                                         .collect();
                    let aa = to_aa(&substr);

                    aa.windows(aa_ksize as usize)
                      .map(|n| self.add_word(n))
                      .count();
                }

                for i in 0..3 {
                    let substr: Vec<u8> = rc.iter()
                                         .cloned()
                                         .skip(i)
                                         .take(sequence.len() - i)
                                         .collect();
                    let aa = to_aa(&substr);

                    aa.windows(aa_ksize as usize)
                      .map(|n| self.add_word(n))
                      .count();
                }






                for (kmer, rc) in sequence
                    .windows((self.ksize * 3) as usize)
                    .zip(rc.windows((self.ksize * 3) as usize))
                {
                    let aa_kmer = to_aa(kmer);
                    let aa_rc = to_aa(rc);

                    //println!("{:?} {:?}, {:?} {:?}", kmer, aa_kmer, rc, aa_rc);

                    if aa_kmer < aa_rc {
                        self.add_word(&aa_kmer);
                    } else {
                        self.add_word(&aa_rc);
                    }
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
            while !self_value.is_none() {
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
            self.mins = merged.iter()
                           .map(|&x| x as u64)
                           .take(self.num as usize)
                           .collect();
            self.abunds = Some(merged_abunds)  // TODO: reduce this one too
        }
        Ok(())
    }

    pub fn count_common(&mut self, other: &KmerMinHash) -> Result<u64> {
        self.check_compatible(other)?;
        let s1: HashSet<&u64> = HashSet::from_iter(self.mins.iter());
        let s2: HashSet<&u64> = HashSet::from_iter(other.mins.iter());
        Ok(s1.intersection(&s2).count() as u64)
    }

    pub fn intersection(&mut self, other: &KmerMinHash) -> Result<(Vec<u64>, u64)> {
        self.check_compatible(other)?;

        let mut combined_mh = KmerMinHash::new(self.num, self.ksize,
            self.is_protein, self.seed, self.max_hash,
            !self.abunds.is_none());

        combined_mh.merge(&self);
        combined_mh.merge(&other);

        let s1: HashSet<_> = HashSet::from_iter(self.mins.iter());
        let s2: HashSet<_> = HashSet::from_iter(other.mins.iter());
        let s3: HashSet<_> = HashSet::from_iter(combined_mh.mins.iter());

        let i1 = &s1 & &s2;
        let i2 = &i1 & &s3;

        let common: Vec<u64> = i2.into_iter().cloned().collect();
        Ok((common, combined_mh.mins.len() as u64))
    }

    pub fn compare(&mut self, other: &KmerMinHash) -> Result<f64> {
        self.check_compatible(other)?;
        if let Ok((common, size)) = self.intersection(other) {
            return Ok(common.len() as f64 / u64::max(1, size) as f64)
        } else {
            return Ok(0.0)
        }
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
             x => x
            } as u8) // TODO: error?
       .collect()
}

lazy_static! {
    static ref CODONTABLE: HashMap<&'static str, u8> = {
        let mut m = HashMap::new();

        m.insert("TTT", 'F' as u8);
        m.insert("TTC", 'F' as u8);
        m.insert("TTA", 'L' as u8);
        m.insert("TTG", 'L' as u8);

        m.insert("TCT", 'S' as u8);
        m.insert("TCC", 'S' as u8);
        m.insert("TCA", 'S' as u8);
        m.insert("TCG", 'S' as u8);

        m.insert("TAT", 'Y' as u8);
        m.insert("TAC", 'Y' as u8);
        m.insert("TAA", '*' as u8);
        m.insert("TAG", '*' as u8);

        m.insert("TGT", 'C' as u8);
        m.insert("TGC", 'C' as u8);
        m.insert("TGA", '*' as u8);
        m.insert("TGG", 'W' as u8);

        m.insert("CTT", 'L' as u8);
        m.insert("CTC", 'L' as u8);
        m.insert("CTA", 'L' as u8);
        m.insert("CTG", 'L' as u8);

        m.insert("CCT", 'P' as u8);
        m.insert("CCC", 'P' as u8);
        m.insert("CCA", 'P' as u8);
        m.insert("CCG", 'P' as u8);

        m.insert("CAT", 'H' as u8);
        m.insert("CAC", 'H' as u8);
        m.insert("CAA", 'Q' as u8);
        m.insert("CAG", 'Q' as u8);

        m.insert("CGT", 'R' as u8);
        m.insert("CGC", 'R' as u8);
        m.insert("CGA", 'R' as u8);
        m.insert("CGG", 'R' as u8);

        m.insert("ATT", 'I' as u8);
        m.insert("ATC", 'I' as u8);
        m.insert("ATA", 'I' as u8);
        m.insert("ATG", 'M' as u8);

        m.insert("ACT", 'T' as u8);
        m.insert("ACC", 'T' as u8);
        m.insert("ACA", 'T' as u8);
        m.insert("ACG", 'T' as u8);

        m.insert("AAT", 'N' as u8);
        m.insert("AAC", 'N' as u8);
        m.insert("AAA", 'K' as u8);
        m.insert("AAG", 'K' as u8);

        m.insert("AGT", 'S' as u8);
        m.insert("AGC", 'S' as u8);
        m.insert("AGA", 'R' as u8);
        m.insert("AGG", 'R' as u8);

        m.insert("GTT", 'V' as u8);
        m.insert("GTC", 'V' as u8);
        m.insert("GTA", 'V' as u8);
        m.insert("GTG", 'V' as u8);

        m.insert("GCT", 'A' as u8);
        m.insert("GCC", 'A' as u8);
        m.insert("GCA", 'A' as u8);
        m.insert("GCG", 'A' as u8);

        m.insert("GAT", 'D' as u8);
        m.insert("GAC", 'D' as u8);
        m.insert("GAA", 'E' as u8);
        m.insert("GAG", 'E' as u8);

        m.insert("GGT", 'D' as u8);
        m.insert("GGC", 'D' as u8);
        m.insert("GGA", 'E' as u8);
        m.insert("GGG", 'E' as u8);

        m
    };
}

#[inline]
fn to_aa(seq: &[u8]) -> Vec<u8> {
    let dna_size = seq.len() / 3;
    seq.chunks(3)
       .map(|n| { //TODO: better error handling
           CODONTABLE.get(str::from_utf8(n).unwrap()).unwrap()
       })
       .cloned()
       .collect()
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
