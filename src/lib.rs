extern crate backtrace;
extern crate murmurhash3;
extern crate ordslice;

#[macro_use]
extern crate error_chain;

#[macro_use]
pub mod errors;

#[macro_use]
pub mod utils;

#[macro_use]
pub mod ffi;

use std::collections::HashSet;
use std::iter::FromIterator;

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
        n: u32,
        k: u32,
        prot: bool,
        seed: u64,
        mx: u64,
        track_abundance: bool,
    ) -> KmerMinHash {
        let mins: Vec<u64>;
        let abunds: Option<Vec<u64>>;

        if n > 0 {
            mins = Vec::with_capacity(n as usize);
        } else {
            mins = Vec::with_capacity(1000);
        }

        if track_abundance {
            abunds = Some(Vec::with_capacity(mins.capacity()));
        } else {
            abunds = None
        }

        KmerMinHash {
            num: n,
            ksize: k,
            is_protein: prot,
            seed: seed,
            max_hash: mx,
            mins: mins,
            abunds: abunds,
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
            if self.mins.len() == 0 {
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
                            self.add_word(&kmer);
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
                for (kmer, rc) in sequence
                    .windows((self.ksize * 3) as usize)
                    .zip(revcomp(sequence).windows((self.ksize * 3) as usize))
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

    pub fn merge(&mut self, other: &KmerMinHash) -> Result<(Vec<u64>, Option<Vec<u64>>)> {
        self.check_compatible(other)?;
        let max_size = self.mins.len() + other.mins.len();
        let mut merged: Vec<u64> = Vec::with_capacity(max_size);
        let mut merged_abunds: Vec<u64> = Vec::with_capacity(max_size);

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

        if merged.len() < (self.num as usize) || (self.num as usize) == 0 {
            return Ok((merged, Some(merged_abunds)));
        } else {
            return Ok((merged
                .iter()
                .map(|&x| x as u64)
                .take(self.num as usize)
                .collect(),
                Some(merged_abunds)));
        }
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

        if let Ok((mins, abunds)) = self.merge(&other) {
            combined_mh.mins = mins;
            combined_mh.abunds = abunds
        }

        let s1: HashSet<_> = HashSet::from_iter(self.mins.iter());
        let s2: HashSet<_> = HashSet::from_iter(other.mins.iter());
        let s3: HashSet<_> = HashSet::from_iter(combined_mh.mins.iter());

        let i1 = &s1 & &s2;
        let i2 = &i1 & &s3;

        let common: Vec<u64> = i2.into_iter().map(|n| *n).collect();
        Ok((common, combined_mh.mins.len() as u64))
    }

    pub fn compare(&mut self, other: &KmerMinHash) -> Result<f64> {
        self.check_compatible(other)?;
        if let Ok((common, size)) = self.intersection(other) {
            return Ok(common.len() as f64 / size as f64)
        } else {
            return Ok(0.0)
        }
    }
}

#[inline]
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
       .map(|n| match *n as char {
             'A' | 'a' => 'T',
             'T' | 't' => 'A',
             'C' | 'c' => 'G',
             'G' | 'g' => 'C',
             x => x
            } as u8) // TODO: error?
       .rev()
       .collect()
}

#[inline]
fn to_aa(seq: &[u8]) -> Vec<u8> {
    seq.iter()
       .map(|n| match *n as char {
             'A' | 'a' => 'T',
             'T' | 't' => 'A',
             'C' | 'c' => 'G',
             'G' | 'g' => 'C',
             x => x
            } as u8) // TODO: error?
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
