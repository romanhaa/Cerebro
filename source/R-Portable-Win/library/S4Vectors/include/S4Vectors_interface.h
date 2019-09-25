/*****************************************************************************
 S4Vectors C interface: prototypes
 ---------------------------------

   The S4Vectors C interface is split in 2 files:
     1. S4Vectors_defines.h (in this directory): contains the typedefs and
        defines of the interface.
     2. S4Vectors_interface.h (this file): contains the prototypes of the
        S4Vectors C routines that are part of the interface.

 *****************************************************************************/
#include "S4Vectors_defines.h"


/*
 * Safe signed integer arithmetic.
 * (see safe_arithm.c)
 */

void reset_ovflow_flag();

int get_ovflow_flag();

int safe_int_add(
	int x,
	int y
);

int safe_int_mult(
	int x,
	int y
);

int as_int(
	const char *val,
	int val_len
);

/*
 * Low-level sorting utilities.
 * (see sort_utils.c)
 */

void sort_int_array(
	int *x,
	size_t nelt,
	int desc
);

void get_order_of_int_array(
	const int *x,
	int nelt,
	int desc,
	int *out,
	int out_shift
);

int sort_ints(
	int *base,
	int base_len,
	const int *x,
	int desc,
	int use_radix,
	unsigned short int *rxbuf1,
	int *rxbuf2
);

void get_order_of_int_pairs(
	const int *a,
	const int *b,
	int nelt,
	int a_desc,
	int b_desc,
	int *out,
	int out_shift
);

int sort_int_pairs(
	int *base,
	int base_len,
	const int *a,
	const int *b,
	int a_desc,
	int b_desc,
	int use_radix,
	unsigned short int *rxbuf1,
	int *rxbuf2
);

void get_matches_of_ordered_int_pairs(
	const int *a1,
	const int *b1,
	const int *o1,
	int nelt1,
	const int *a2,
	const int *b2,
	const int *o2,
	int nelt2,
	int nomatch,
	int *out,
	int out_shift
);

void get_order_of_int_quads(
	const int *a,
	const int *b,
	const int *c,
	const int *d,
	int nelt,
	int a_desc,
	int b_desc,
	int c_desc,
	int d_desc,
	int *out,
	int out_shift
);

int sort_int_quads(
	int *base,
	int base_len,
	const int *a,
	const int *b,
	const int *c,
	const int *d,
	int a_desc,
	int b_desc,
	int c_desc,
	int d_desc,
	int use_radix,
	unsigned short int *rxbuf1,
	int *rxbuf2
);

void get_matches_of_ordered_int_quads(
	const int *a1,
	const int *b1,
	const int *c1,
	const int *d1,
	const int *o1,
	int nelt1,
	const int *a2,
	const int *b2,
	const int *c2,
	const int *d2,
	const int *o2,
	int nelt2,
	int nomatch,
	int *out,
	int out_shift
);

/*
 * Hash table management.
 * (see hash_utils.c)
 */

struct htab new_htab(int n);

int get_hbucket_val(
	const struct htab *htab,
	int bucket_idx
);

void set_hbucket_val(
	struct htab *htab,
	int bucket_idx,
	int val
);

/*
 * Low-level manipulation of the Auto-Extending buffers.
 * (see AEbufs.c)
 */

size_t increase_buflength(size_t buflength);

size_t IntAE_get_nelt(const IntAE *ae);

size_t IntAE_set_nelt(
	IntAE *ae,
	size_t nelt
);

void IntAE_set_val(
	const IntAE *ae,
	int val
);

void IntAE_extend(
	IntAE *ae,
	size_t new_buflength
);

void IntAE_insert_at(
	IntAE *ae,
	size_t at,
	int val
);

IntAE *new_IntAE(
	size_t buflength,
	size_t nelt,
	int val
);

void IntAE_append(
	IntAE *ae,
	const int *newvals,
	size_t nnewval
);

void IntAE_delete_at(
	IntAE *ae,
	size_t at,
	size_t nelt
);

void IntAE_shift(
	const IntAE *ae,
	size_t offset,
	int shift
);

void IntAE_sum_and_shift(
	const IntAE *ae1,
	const IntAE *ae2,
	int shift
);

void IntAE_qsort(
	const IntAE *ae,
	size_t offset,
	int desc
);

void IntAE_uniq(
	IntAE *ae,
	size_t offset
);

SEXP new_INTEGER_from_IntAE(const IntAE *ae);

IntAE *new_IntAE_from_INTEGER(SEXP x);

IntAE *new_IntAE_from_CHARACTER(
	SEXP x,
	int keyshift
);

size_t IntAEAE_get_nelt(const IntAEAE *aeae);

size_t IntAEAE_set_nelt(
	IntAEAE *aeae,
	size_t nelt
);

void IntAEAE_extend(
	IntAEAE *aeae,
	size_t new_buflength
);

void IntAEAE_insert_at(
	IntAEAE *aeae,
	size_t at,
	IntAE *ae
);

IntAEAE *new_IntAEAE(
	size_t buflength,
	size_t nelt
);

void IntAEAE_pappend(
	const IntAEAE *aeae1,
	const IntAEAE *aeae2
);

void IntAEAE_shift(
	const IntAEAE *aeae,
	int shift
);

void IntAEAE_sum_and_shift(
	const IntAEAE *aeae1,
	const IntAEAE *aeae2,
	int shift
);

SEXP new_LIST_from_IntAEAE(
	const IntAEAE *aeae,
	int mode
);

IntAEAE *new_IntAEAE_from_LIST(SEXP x);

SEXP IntAEAE_toEnvir(
	const IntAEAE *aeae,
	SEXP envir,
	int keyshift
);

size_t IntPairAE_get_nelt(const IntPairAE *ae);

size_t IntPairAE_set_nelt(
	IntPairAE *ae,
	size_t nelt
);

void IntPairAE_extend(
	IntPairAE *ae,
	size_t new_buflength
);

void IntPairAE_insert_at(
	IntPairAE *ae,
	size_t at,
	int a,
	int b
);

IntPairAE *new_IntPairAE(
	size_t buflength,
	size_t nelt
);

size_t IntPairAEAE_get_nelt(const IntPairAEAE *aeae);

size_t IntPairAEAE_set_nelt(
	IntPairAEAE *aeae,
	size_t nelt
);

void IntPairAEAE_extend(
	IntPairAEAE *aeae,
	size_t new_buflength
);

void IntPairAEAE_insert_at(
	IntPairAEAE *aeae,
	size_t at,
	IntPairAE *ae
);

IntPairAEAE *new_IntPairAEAE(
	size_t buflength,
	size_t nelt
);

size_t LLongAE_get_nelt(const LLongAE *ae);

size_t LLongAE_set_nelt(
	LLongAE *ae,
	size_t nelt
);

void LLongAE_set_val(
	const LLongAE *ae,
	long long val
);

void LLongAE_extend(
	LLongAE *ae,
	size_t new_buflength
);

void LLongAE_insert_at(
	LLongAE *ae,
	size_t at,
	long long val
);

LLongAE *new_LLongAE(
	size_t buflength,
	size_t nelt,
	long long val
);

size_t CharAE_get_nelt(const CharAE *ae);

size_t CharAE_set_nelt(
	CharAE *ae,
	size_t nelt
);

void CharAE_extend(
	CharAE *ae,
	size_t new_buflength
);

void CharAE_insert_at(
	CharAE *ae,
	size_t at,
	char c
);

CharAE *new_CharAE(size_t buflength);

CharAE *new_CharAE_from_string(const char *string);

void CharAE_append_string(
	CharAE *ae,
	const char *string
);

void CharAE_delete_at(
	CharAE *ae,
	size_t at,
	size_t nelt
);

SEXP new_CHARSXP_from_CharAE(const CharAE *ae);

SEXP new_RAW_from_CharAE(const CharAE *ae);

SEXP new_LOGICAL_from_CharAE(const CharAE *ae);

size_t CharAEAE_get_nelt(const CharAEAE *aeae);

size_t CharAEAE_set_nelt(
	CharAEAE *aeae,
	size_t nelt
);

void CharAEAE_extend(
	CharAEAE *aeae,
	size_t new_buflength
);

void CharAEAE_insert_at(
	CharAEAE *aeae,
	size_t at,
	CharAE *ae
);

CharAEAE *new_CharAEAE(
	size_t buflength,
	size_t nelt
);

void CharAEAE_append_string(
	CharAEAE *aeae,
	const char *string
);

SEXP new_CHARACTER_from_CharAEAE(const CharAEAE *aeae);

/*
 * SEXP_utils.c
 */

const char *get_classname(SEXP x);

/*
 * LLint_class.c
 */

int is_LLint(SEXP x);

R_xlen_t get_LLint_length(SEXP x);

long long int *get_LLint_dataptr(SEXP x);

SEXP alloc_LLint(const char *classname, R_xlen_t length);

/*
 * subsetting_utils.c
 */

long long int copy_vector_block(
	SEXP dest,
	long long int dest_offset,
	SEXP src,
	long long int src_offset,
	long long int block_nelt
);

int copy_vector_positions(
	SEXP dest,
	int dest_offset,
	SEXP src,
	const int *pos,
	int npos
);

int copy_vector_ranges(
	SEXP dest,
	int dest_offset,
	SEXP src,
	const int *start,
	const int *width,
	int nranges
);

/*
 * vector_utils.c
 */

int vector_memcmp(
	SEXP x1,
	int x1_offset,
	SEXP x2,
	int x2_offset,
	int nelt
);

SEXP list_as_data_frame(
	SEXP x,
	int nrow
);

/*
 * integer_utils.c
 */

int check_integer_pairs(
	SEXP a,
	SEXP b,
	const int **a_p,
	const int **b_p,
	const char *a_argname,
	const char *b_argname
);

SEXP find_interv_and_start_from_width(
	const int *x,
	int x_len,
	const int *width,
	int width_len
);

/*
 * Low-level manipulation of Hits objects.
 * (see Hits_class.c)
 */

SEXP new_Hits(
	const char *Class,
	int *from,
	const int *to,
	int nhit,
	int nLnode,
	int nRnode,
	int already_sorted
);

int get_select_mode(SEXP select);

/*
 * Low-level manipulation of Rle objects.
 * (see Rle_class.c)
 */

SEXP construct_logical_Rle(
	R_xlen_t nrun_in,
	const int *values_in,
	const void *lengths_in,
	int lengths_in_is_L
);

SEXP construct_integer_Rle(
	R_xlen_t nrun_in,
	const int *values_in,
	const void *lengths_in,
	int lengths_in_is_L
);

SEXP construct_numeric_Rle(
	R_xlen_t nrun_in,
	const double *values_in,
	const void *lengths_in,
	int lengths_in_is_L
);

SEXP construct_complex_Rle(
	R_xlen_t nrun_in,
	const Rcomplex *values_in,
	const void *lengths_in,
	int lengths_in_is_L
);

SEXP construct_character_Rle(
	SEXP values_in,
	const void *lengths_in,
	int lengths_in_is_L
);

SEXP construct_raw_Rle(
	R_xlen_t nrun_in,
	const Rbyte *values_in,
	const void *lengths_in,
	int lengths_in_is_L
);

SEXP construct_Rle(
	SEXP values_in,
	const void *lengths_in,
	int lengths_in_is_L
);

/*
 * Low-level manipulation of Vector objects.
 * (see List_class.c)
 */

const char *get_List_elementType(SEXP x);

void set_List_elementType(SEXP x, const char *type);

/*
 * Low-level manipulation of SimpleList objects.
 * (see SimpleList_class.c)
 */

SEXP new_SimpleList(const char *classname, SEXP listData);

/*
 * Low-level manipulation of DataFrame objects.
 * (see DataFrame_class.c)
 */

SEXP new_DataFrame(const char *classname, SEXP vars, SEXP rownames, SEXP nrows);

