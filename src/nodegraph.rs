use std::fs::File;
use std::io;
use std::path::Path;

use byteorder::{BigEndian, LittleEndian, ReadBytesExt, WriteBytesExt};
use fixedbitset::FixedBitSet;

use errors::Result;

type HashIntoType = u64;

#[derive(Debug)]
pub struct Nodegraph {
    bs: Vec<FixedBitSet>,
    ksize: usize,
    occupied_bins: usize,
    unique_kmers: usize,
}

impl Nodegraph {
    pub fn new(tablesizes: &[usize], ksize: usize) -> Nodegraph {
        let mut bs = Vec::with_capacity(tablesizes.len());
        for size in tablesizes.iter() {
            bs.push(FixedBitSet::with_capacity(*size));
        }

        Nodegraph {
            bs,
            ksize,
            occupied_bins: 0,
            unique_kmers: 0,
        }
    }

    pub fn count(&mut self, hash: HashIntoType) -> bool {
        let mut is_new_kmer = false;

        for bitset in &mut self.bs {
            let bin = hash % bitset.len() as u64;
            if !bitset.put(bin as usize) {
                self.occupied_bins += 1;
                is_new_kmer = true;
            }
        }

        if is_new_kmer {
            self.unique_kmers += 1
        }
        return is_new_kmer;
    }

    pub fn get(&self, hash: HashIntoType) -> usize {
        for bitset in &self.bs {
            let bin = hash % bitset.len() as u64;
            if !bitset.contains(bin as usize) {
                return 0;
            }
        }
        return 1;
    }

    // update
    pub fn update(&mut self, other: &Nodegraph) {}

    // save
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        self.save_to_writer(&mut File::open(path)?)?;
        Ok(())
    }

    pub fn save_to_writer<W>(&self, wtr: &mut W) -> Result<()>
    where
        W: io::Write,
    {
        wtr.write(b"OXLI")?;
        wtr.write_u8(4)?; // version
        wtr.write_u8(2)?; // ht_type
        wtr.write_u32::<LittleEndian>(self.ksize as u32)?; // ksize
        wtr.write_u8(self.bs.len() as u8)?; // n_tables
        wtr.write_u64::<LittleEndian>(self.occupied_bins as u64)?; // n_occupied
        for count in &self.bs {
            wtr.write_u64::<LittleEndian>(count.len() as u64)?;
            for (i, chunk) in count.as_slice().iter().enumerate() {
                let next = (i + 1) * 32;
                if next <= count.len() {
                    wtr.write_u32::<LittleEndian>(*chunk).unwrap()
                } else {
                    let rem = count.len() - (i * 32);
                    let remainder = if rem % 8 != 0 { rem / 8 + 1 } else { rem / 8 };

                    if remainder == 0 {
                        wtr.write_u8(0).unwrap()
                    } else {
                        for pos in 0..remainder {
                            let byte: u8 = (chunk.wrapping_shr(pos as u32 * 8) & 0xff) as u8;
                            wtr.write_u8(byte).unwrap()
                        }
                    }
                }
            }
        }
        Ok(())
    }

    pub fn from_reader<R>(rdr: &mut R) -> Result<Nodegraph>
    where
        R: io::Read,
    {
        let signature = rdr.read_u32::<BigEndian>()?;
        assert_eq!(signature, 0x4f584c49);

        let version = rdr.read_u8()?;
        assert_eq!(version, 0x04);

        let ht_type = rdr.read_u8()?;
        assert_eq!(ht_type, 0x02);

        let ksize = rdr.read_u32::<LittleEndian>()?;
        let n_tables = rdr.read_u8()?;
        let occupied_bins = rdr.read_u64::<LittleEndian>()? as usize;

        let mut _n_occupied = 0;
        let mut bs = Vec::with_capacity(n_tables as usize);
        for _i in 0..n_tables {
            let tablesize: usize = rdr.read_u64::<LittleEndian>()? as usize;
            let byte_size = tablesize / 8 + 1;

            let mut counts = FixedBitSet::with_capacity(tablesize);
            for pos in 0..byte_size {
                let mut byte = rdr.read_u8()?;
                if byte == 0 {
                    continue;
                }

                _n_occupied += 1;
                for i in 0..8u32 {
                    if byte & (1 << i) != 0 {
                        counts.insert(pos * 8 + i as usize);
                    }
                }
            }

            bs.push(counts);
        }

        //assert_eq!(occupied_bins, _n_occupied);
        Ok(Nodegraph {
            bs,
            ksize: ksize as usize,
            occupied_bins,
            unique_kmers: 0, // This is a khmer issue, it doesn't save unique_kmers
        })
    }

    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Nodegraph> {
        Ok(Nodegraph::from_reader(&mut File::open(path)?)?)
    }

    pub fn tablesizes(&self) -> Vec<usize> {
        self.bs.iter().map(|x| x.len()).collect()
    }

    pub fn n_occupied_bins(&self) -> usize {
        //self.bs.iter().map(|x| x.count_ones(..)).sum::<usize>() / 8
        self.occupied_bins
    }

    pub fn unique_kmers(&self) -> usize {
        self.unique_kmers
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use proptest::num::u64;

    proptest!{
      #[test]
      fn count_and_get(hash in u64::ANY) {
          let mut ng: Nodegraph = Nodegraph::new(&[10], 3);
          ng.count(hash);
          assert_eq!(ng.get(hash), 1);
      }
    }
}
