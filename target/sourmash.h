#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>

enum SourmashErrorCode {
  NoError = 0,
  Panic = 1,
  Internal = 2,
  Msg = 3,
  Unknown = 4,
  MismatchKSizes = 101,
  MismatchDNAProt = 102,
  MismatchMaxHash = 103,
  MismatchSeed = 104,
  InvalidDNA = 1101,
  InvalidProt = 1102,
  Io = 100001,
  Utf8Error = 100002,
  ParseInt = 100003,
};
typedef uint32_t SourmashErrorCode;

typedef struct KmerMinHash KmerMinHash;

uint64_t hash_murmur(const char *kmer, uint64_t seed);

void kmerminhash_add_hash(KmerMinHash *ptr, uint64_t h);

void kmerminhash_add_sequence(KmerMinHash *ptr, const char *sequence, bool force);

void kmerminhash_add_word(KmerMinHash *ptr, const char *word);

uint64_t kmerminhash_count_common(KmerMinHash *ptr, const KmerMinHash *other);

void kmerminhash_free(KmerMinHash *ptr);

uint64_t kmerminhash_get_min_idx(KmerMinHash *ptr, uint64_t idx);

size_t kmerminhash_get_mins_size(KmerMinHash *ptr);

bool kmerminhash_is_protein(KmerMinHash *ptr);

uint32_t kmerminhash_ksize(KmerMinHash *ptr);

uint64_t kmerminhash_max_hash(KmerMinHash *ptr);

void kmerminhash_merge(KmerMinHash *ptr, const KmerMinHash *other);

KmerMinHash *kmerminhash_new(uint32_t n,
                             uint32_t k,
                             bool prot,
                             uint64_t seed,
                             uint64_t mx,
                             bool track_abundance);

uint32_t kmerminhash_num(KmerMinHash *ptr);

uint64_t kmerminhash_seed(KmerMinHash *ptr);

/*
 * Clears the last error.
 */
void sourmash_err_clear();

/*
 * Returns the last error code.
 *
 * If there is no error, 0 is returned.
 */
SourmashErrorCode sourmash_err_get_last_code();

/*
 * Initializes the library
 */
void sourmash_init();
