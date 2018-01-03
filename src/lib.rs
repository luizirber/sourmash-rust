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
    pub abunds: Vec<u64>,
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
        let abunds: Vec<u64>;

        if n > 0 {
            mins = Vec::with_capacity(n as usize);
        } else {
            mins = Vec::with_capacity(1000);
        }

        if track_abundance {
            abunds = Vec::with_capacity(mins.capacity());
        } else {
            abunds = Vec::new();
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
            if self.mins.len() == 0 {
                self.mins.push(hash);
                return;
            } else if hash <= self.max_hash || current_max > hash
                || (self.mins.len() as u32) < self.num
            {
                let pos = self.mins.lower_bound(&hash);

                if pos == self.mins.len() {
                    self.mins.push(hash);
                } else if self.mins[pos] != hash {
                    self.mins.insert(pos, hash);
                    if self.num != 0 && self.mins.len() > (self.num as usize) {
                        self.mins.pop();
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
                      if ! force {
                          let last = vec![kmer.last().cloned().unwrap()];
                          return Err(ErrorKind::InvalidDNA(String::from_utf8(last).unwrap()).into());
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

    pub fn merge(&mut self, other: &KmerMinHash) -> Result<Vec<u64>> {
        self.check_compatible(other)?;
        let mut merged: Vec<u64> = Vec::with_capacity(self.mins.len() + other.mins.len());
        let mut self_iter = self.mins.iter();
        let mut other_iter = other.mins.iter();

        let mut self_value: u64 = match self_iter.next() {
            Some(x) => *x,
            None => {
                merged.extend(other_iter);
                return Ok(merged);
            }
        };

        loop {
            match other_iter.next() {
                None => {
                    merged.push(self_value);
                    merged.extend(self_iter);
                    // TODO: copy abunds too
                    break;
                }
                Some(x) if *x < self_value => merged.push(*x),
                Some(x) if *x == self_value => {
                    merged.push(*x);
                    // TODO: sum abunds
                    self_value = match self_iter.next() {
                        None => break,
                        Some(x) => *x,
                    }
                }
                Some(x) if *x > self_value => {
                    merged.push(self_value);
                    self_value = match self_iter.next() {
                        None => break,
                        Some(x) => *x,
                    }
                }
                Some(_) => {}
            }
        }
        merged.extend(other_iter);

        if merged.len() < (self.num as usize) || (self.num as usize) == 0 {
            return Ok(merged);
        } else {
            return Ok(merged
                .iter()
                .map(|&x| x as u64)
                .take(self.num as usize)
                .collect());
        }
    }

    pub fn count_common(&mut self, other: &KmerMinHash) -> Result<u64> {
        self.check_compatible(other)?;
        let s1: HashSet<&u64> = HashSet::from_iter(self.mins.iter());
        let s2: HashSet<&u64> = HashSet::from_iter(other.mins.iter());
        Ok(s1.intersection(&s2).count() as u64)
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
