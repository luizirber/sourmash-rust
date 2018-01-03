/* c bindings to the sourmash library */

#ifndef SOURMASH_H_INCLUDED
#define SOURMASH_H_INCLUDED

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>

enum SourmashErrorCode {
  SOURMASH_ERROR_CODE_NO_ERROR = 0,
  SOURMASH_ERROR_CODE_PANIC = 1,
  SOURMASH_ERROR_CODE_INTERNAL = 2,
  SOURMASH_ERROR_CODE_MSG = 3,
  SOURMASH_ERROR_CODE_UNKNOWN = 4,
  SOURMASH_ERROR_CODE_MISMATCH_K_SIZES = 101,
  SOURMASH_ERROR_CODE_MISMATCH_D_N_A_PROT = 102,
  SOURMASH_ERROR_CODE_MISMATCH_MAX_HASH = 103,
  SOURMASH_ERROR_CODE_MISMATCH_SEED = 104,
  SOURMASH_ERROR_CODE_INVALID_D_N_A = 1101,
  SOURMASH_ERROR_CODE_INVALID_PROT = 1102,
  SOURMASH_ERROR_CODE_IO = 100001,
  SOURMASH_ERROR_CODE_UTF8_ERROR = 100002,
  SOURMASH_ERROR_CODE_PARSE_INT = 100003,
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

#endif /* SOURMASH_H_INCLUDED */
