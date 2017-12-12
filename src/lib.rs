extern crate murmurhash3;
extern crate ordslice;

use murmurhash3::murmurhash3_x64_128;
use ordslice::Ext;

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
    pub abunds: Vec<u64>
}

impl KmerMinHash {
    pub fn new(n: u32, k: u32, prot: bool, seed: u64,
               mx: u64, track_abundance: bool) -> KmerMinHash {
        let mins: Vec<u64>;
        let abunds: Vec<u64>;

        if n > 0 {
            mins = Vec::with_capacity(n as usize);
        } else {
            mins = Vec::with_capacity(1000);
        }

        if track_abundance { abunds = Vec::with_capacity(mins.capacity());
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
            abunds: abunds
        }
    }

    pub fn add_hash(&mut self, hash: u64) {
       let current_max = match self.mins.last() {
         Some(&x) => x,
         None => u64::max_value(),
       };
       println!("{}", current_max);

       if (self.max_hash != 0 && hash <= self.max_hash) || self.max_hash == 0 {
           println!("adding {}", current_max);
           if self.mins.len() == 0 {
             self.mins.push(hash);
             return
           } else if hash <= self.max_hash ||
                     current_max > hash ||
                     (self.mins.len() as u32) < self.num {
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

    pub fn add_sequence(&mut self, sequence: &[u8], force: bool) {
        if sequence.len() < (self.ksize as usize) {
            return
        }

        if !self.is_protein { // dna
            for kmer in sequence.windows(self.ksize as usize) {
                let rc = revcomp(kmer);
                if kmer < &rc {
                    self.add_word(&kmer);
                } else {
                    self.add_word(&rc);
                }
            }
        } else {  // protein
            for (kmer, rc) in sequence.windows((self.ksize * 3) as usize).zip(
                              revcomp(sequence).windows((self.ksize * 3) as usize)) {
                let aa_kmer = to_aa(kmer);
                let aa_rc = to_aa(rc);

                println!("{:?} {:?}, {:?} {:?}", kmer, aa_kmer, rc, aa_rc);

                if aa_kmer < aa_rc {
                    self.add_word(&aa_kmer);
                } else {
                    self.add_word(&aa_rc);
                }
            }
        }
    }

    pub fn merge(&mut self, other: KmerMinHash) -> Vec<u64> {
        // self.check_compatible(other);
        let mut merged: Vec<u64> = Vec::with_capacity(self.mins.len() + other.mins.len());
        let mut self_iter = self.mins.iter();
        let mut other_iter = other.mins.iter();

        let mut self_value = self_iter.next().unwrap();
        loop {
            if Some(self_value) == None {
                break;
            }
            match other_iter.next() {
                None => {
                    merged.extend(self_iter);
                    // TODO: copy abunds too
                    break },
                Some(x) if x < self_value => merged.push(*x),
                Some(x) if x == self_value => {
                    merged.push(*x);
                    // TODO: sum abunds
                    self_value = self_iter.next().unwrap()},
                Some(x) if x > self_value => {
                    merged.push(*self_value);
                    self_value = self_iter.next().unwrap()},
                Some(_) => {}
            }
        }
        merged.extend(other_iter);

        merged
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
