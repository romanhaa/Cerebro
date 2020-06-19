#include "S4Vectors_interface.h"

#define DEFINE_CCALLABLE_STUB(retT, stubname, Targs, args) \
typedef retT(*__ ## stubname ## _funtype__)Targs; \
retT stubname Targs \
{ \
	static __ ## stubname ## _funtype__ fun = NULL; \
	if (fun == NULL) \
		fun = (__ ## stubname ## _funtype__) R_GetCCallable("S4Vectors", "_" #stubname); \
	return fun args; \
}

/*
 * Using the above macro when retT (the returned type) is void will make Sun
 * Studio 12 C compiler unhappy. So we need to use the following macro to
 * handle that case.
 */
#define DEFINE_NOVALUE_CCALLABLE_STUB(stubname, Targs, args) \
typedef void(*__ ## stubname ## _funtype__)Targs; \
void stubname Targs \
{ \
	static __ ## stubname ## _funtype__ fun = NULL; \
	if (fun == NULL) \
		fun = (__ ## stubname ## _funtype__) R_GetCCallable("S4Vectors", "_" #stubname); \
	fun args; \
	return; \
}


/*
 * Stubs for callables defined in safe_arithm.c
 */

DEFINE_NOVALUE_CCALLABLE_STUB(reset_ovflow_flag,
	(),
	()
)

DEFINE_CCALLABLE_STUB(int, get_ovflow_flag,
	(),
	()
)

DEFINE_CCALLABLE_STUB(int, safe_int_add,
	(int x, int y),
	(    x,     y)
)

DEFINE_CCALLABLE_STUB(int, safe_int_mult,
	(int x, int y),
	(    x,     y)
)

DEFINE_CCALLABLE_STUB(int, as_int,
	(const char *val, int val_len),
	(            val,     val_len)
)

/*
 * Stubs for callables defined in sort_utils.c
 */

DEFINE_NOVALUE_CCALLABLE_STUB(sort_int_array,
	(int *x, size_t nelt, int desc),
	(     x,        nelt,     desc)
)

DEFINE_NOVALUE_CCALLABLE_STUB(get_order_of_int_array,
	(const int *x, int nelt, int desc, int *out, int out_shift),
	(           x,     nelt,     desc,      out,     out_shift)
)

DEFINE_CCALLABLE_STUB(int, sort_ints,
	(int *base, int base_len, const int *x, int desc, int use_radix, unsigned short int *rxbuf1, int *rxbuf2),
	(     base,     base_len,            x,     desc,     use_radix,                     rxbuf1,      rxbuf2)
)

DEFINE_NOVALUE_CCALLABLE_STUB(get_order_of_int_pairs,
	(const int *a, const int *b, int nelt, int a_desc, int b_desc, int *out, int out_shift),
	(           a,            b,     nelt,     a_desc,     b_desc,      out,     out_shift)
)

DEFINE_CCALLABLE_STUB(int, sort_int_pairs,
	(int *base, int base_len, const int *a, const int *b, int a_desc, int b_desc, int use_radix, unsigned short int *rxbuf1, int *rxbuf2),
	(     base,     base_len,            a,            b,     a_desc,     b_desc,     use_radix,                     rxbuf1,      rxbuf2)
)

DEFINE_NOVALUE_CCALLABLE_STUB(get_matches_of_ordered_int_pairs,
	(const int *a1, const int *b1, const int *o1, int nelt1, const int *a2, const int *b2, const int *o2, int nelt2, int nomatch, int *out, int out_shift),
	(           a1,            b1,            o1,     nelt1,            a2,            b2,            o2,     nelt2,     nomatch,      out,     out_shift)
)

DEFINE_NOVALUE_CCALLABLE_STUB(get_order_of_int_quads,
	(const int *a, const int *b, const int *c, const int *d, int nelt, int a_desc, int b_desc, int c_desc, int d_desc, int *out, int out_shift),
	(           a,            b,            c,            d,     nelt,     a_desc,     b_desc,     c_desc,     d_desc,      out,     out_shift)
)

DEFINE_CCALLABLE_STUB(int, sort_int_quads,
	(int *base, int base_len, const int *a, const int *b, const int *c, const int *d, int a_desc, int b_desc, int c_desc, int d_desc, int use_radix, unsigned short int *rxbuf1, int *rxbuf2),
	(     base,     base_len,            a,            b,            c,            d,     a_desc,     b_desc,     c_desc,     d_desc,     use_radix,                     rxbuf1,      rxbuf2)
)

DEFINE_NOVALUE_CCALLABLE_STUB(get_matches_of_ordered_int_quads,
	(const int *a1, const int *b1, const int *c1, const int *d1, const int *o1, int nelt1, const int *a2, const int *b2, const int *c2, const int *d2, const int *o2, int nelt2, int nomatch, int *out, int out_shift),
	(           a1,            b1,            c1,            d1,            o1,     nelt1,            a2,            b2,            c2,            d2,            o2,     nelt2,     nomatch,      out,     out_shift)
)

/*
 * Stubs for callables defined in hash_utils.c
 */

DEFINE_CCALLABLE_STUB(struct htab, new_htab,
	(int n),
	(    n)
)

DEFINE_CCALLABLE_STUB(int, get_hbucket_val,
	(const struct htab *htab, int bucket_idx),
	(                   htab,     bucket_idx)
)

DEFINE_NOVALUE_CCALLABLE_STUB(set_hbucket_val,
	(struct htab *htab, int bucket_idx, int val),
	(             htab,     bucket_idx,     val)
)

/*
 * Stubs for callables defined in AEbufs.c
 */

DEFINE_CCALLABLE_STUB(size_t, increase_buflength,
	(size_t buflength),
	(       buflength)
)

DEFINE_CCALLABLE_STUB(size_t, IntAE_get_nelt,
	(const IntAE *ae),
	(             ae)
)

DEFINE_CCALLABLE_STUB(size_t, IntAE_set_nelt,
	(IntAE *ae, size_t nelt),
	(       ae,        nelt)
)

DEFINE_NOVALUE_CCALLABLE_STUB(IntAE_set_val,
	(const IntAE *ae, int val),
	(             ae,     val)
)

DEFINE_NOVALUE_CCALLABLE_STUB(IntAE_extend,
	(IntAE *ae, size_t new_buflength),
	(       ae,        new_buflength)
)

DEFINE_NOVALUE_CCALLABLE_STUB(IntAE_insert_at,
	(IntAE *ae, size_t at, int val),
	(       ae,        at,     val)
)

DEFINE_CCALLABLE_STUB(IntAE *, new_IntAE,
	(size_t buflength, size_t nelt, int val),
	(       buflength,        nelt,     val)
)

DEFINE_NOVALUE_CCALLABLE_STUB(IntAE_append,
	(IntAE *ae, const int *newvals, size_t nnewval),
	(       ae,            newvals,        nnewval)
)

DEFINE_NOVALUE_CCALLABLE_STUB(IntAE_delete_at,
	(IntAE *ae, size_t at, size_t nelt),
	(       ae,        at,        nelt)
)

DEFINE_NOVALUE_CCALLABLE_STUB(IntAE_shift,
	(const IntAE *ae, size_t offset, int shift),
	(             ae,        offset,     shift)
)

DEFINE_NOVALUE_CCALLABLE_STUB(IntAE_sum_and_shift,
	(const IntAE *ae1, const IntAE *ae2, int shift),
	(             ae1,              ae2,     shift)
)

DEFINE_NOVALUE_CCALLABLE_STUB(IntAE_qsort,
	(const IntAE *ae, size_t offset, int desc),
	(             ae,        offset,     desc)
)

DEFINE_NOVALUE_CCALLABLE_STUB(IntAE_uniq,
	(IntAE *ae, size_t offset),
	(       ae,        offset)
)

DEFINE_CCALLABLE_STUB(SEXP, new_INTEGER_from_IntAE,
	(const IntAE *ae),
	(             ae)
)

DEFINE_CCALLABLE_STUB(IntAE *, new_IntAE_from_INTEGER,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(IntAE *, new_IntAE_from_CHARACTER,
	(SEXP x, int keyshift),
	(     x,     keyshift)
)

DEFINE_CCALLABLE_STUB(size_t, IntAEAE_get_nelt,
	(const IntAEAE *aeae),
	(               aeae)
)

DEFINE_CCALLABLE_STUB(size_t, IntAEAE_set_nelt,
	(IntAEAE *aeae, size_t nelt),
	(         aeae,        nelt)
)

DEFINE_NOVALUE_CCALLABLE_STUB(IntAEAE_extend,
	(IntAEAE *aeae, size_t new_buflength),
	(         aeae,        new_buflength)
)

DEFINE_NOVALUE_CCALLABLE_STUB(IntAEAE_insert_at,
	(IntAEAE *aeae, size_t at, IntAE *ae),
	(         aeae,        at,        ae)
)

DEFINE_CCALLABLE_STUB(IntAEAE *, new_IntAEAE,
	(size_t buflength, size_t nelt),
	(       buflength,        nelt)
)

DEFINE_NOVALUE_CCALLABLE_STUB(IntAEAE_pappend,
	(const IntAEAE *aeae1, const IntAEAE *aeae2),
	(               aeae1,                aeae2)
)

DEFINE_NOVALUE_CCALLABLE_STUB(IntAEAE_shift,
	(const IntAEAE *aeae, int shift),
	(               aeae,     shift)
)

DEFINE_NOVALUE_CCALLABLE_STUB(IntAEAE_sum_and_shift,
	(const IntAEAE *aeae1, const IntAEAE *aeae2, int shift),
	(               aeae1,                aeae2,     shift)
)

DEFINE_CCALLABLE_STUB(SEXP, new_LIST_from_IntAEAE,
	(const IntAEAE *aeae, int mode),
	(               aeae,     mode)
)

DEFINE_CCALLABLE_STUB(IntAEAE *, new_IntAEAE_from_LIST,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(SEXP, IntAEAE_toEnvir,
	(const IntAEAE *aeae, SEXP envir, int keyshift),
	(               aeae,      envir,     keyshift)
)

DEFINE_CCALLABLE_STUB(size_t, IntPairAE_get_nelt,
	(const IntPairAE *ae),
	(                 ae)
)

DEFINE_CCALLABLE_STUB(size_t, IntPairAE_set_nelt,
	(IntPairAE *ae, size_t nelt),
	(           ae,     nelt)
)

DEFINE_NOVALUE_CCALLABLE_STUB(IntPairAE_extend,
	(IntPairAE *ae, size_t new_buflength),
	(           ae,        new_buflength)
)

DEFINE_NOVALUE_CCALLABLE_STUB(IntPairAE_insert_at,
	(IntPairAE *ae, size_t at, int a, int b),
	(           ae,        at,     a,     b)
)

DEFINE_CCALLABLE_STUB(IntPairAE *, new_IntPairAE,
	(size_t buflength, size_t nelt),
	(       buflength,        nelt)
)

DEFINE_CCALLABLE_STUB(size_t, IntPairAEAE_get_nelt,
	(const IntPairAEAE *aeae),
	(                   aeae)
)

DEFINE_CCALLABLE_STUB(size_t, IntPairAEAE_set_nelt,
	(IntPairAEAE *aeae, size_t nelt),
	(             aeae,        nelt)
)

DEFINE_NOVALUE_CCALLABLE_STUB(IntPairAEAE_extend,
	(IntPairAEAE *aeae, size_t new_buflength),
	(             aeae,        new_buflength)
)

DEFINE_NOVALUE_CCALLABLE_STUB(IntPairAEAE_insert_at,
	(IntPairAEAE *aeae, size_t at, IntPairAE *ae),
	(             aeae,        at,            ae)
)

DEFINE_CCALLABLE_STUB(IntPairAEAE *, new_IntPairAEAE,
	(size_t buflength, size_t nelt),
	(       buflength,        nelt)
)

DEFINE_CCALLABLE_STUB(size_t, LLongAE_get_nelt,
	(const LLongAE *ae),
	(               ae)
)

DEFINE_CCALLABLE_STUB(size_t, LLongAE_set_nelt,
	(LLongAE *ae, size_t nelt),
	(         ae,        nelt)
)

DEFINE_NOVALUE_CCALLABLE_STUB(LLongAE_set_val,
	(const LLongAE *ae, long long val),
	(               ae,           val)
)

DEFINE_NOVALUE_CCALLABLE_STUB(LLongAE_extend,
	(LLongAE *ae, size_t new_buflength),
	(         ae,        new_buflength)
)

DEFINE_NOVALUE_CCALLABLE_STUB(LLongAE_insert_at,
	(LLongAE *ae, size_t at, long long val),
	(         ae,     at,           val)
)

DEFINE_CCALLABLE_STUB(LLongAE *, new_LLongAE,
	(size_t buflength, size_t nelt, long long val),
	(       buflength,        nelt,           val)
)

DEFINE_CCALLABLE_STUB(size_t, CharAE_get_nelt,
	(const CharAE *ae),
	(              ae)
)

DEFINE_CCALLABLE_STUB(size_t, CharAE_set_nelt,
	(CharAE *ae, size_t nelt),
	(        ae,        nelt)
)

DEFINE_NOVALUE_CCALLABLE_STUB(CharAE_extend,
	(CharAE *ae, size_t new_buflength),
	(        ae,        new_buflength)
)

DEFINE_NOVALUE_CCALLABLE_STUB(CharAE_insert_at,
	(CharAE *ae, size_t at, char c),
	(        ae,        at,      c)
)

DEFINE_CCALLABLE_STUB(CharAE *, new_CharAE,
	(size_t buflength),
	(       buflength)
)

DEFINE_CCALLABLE_STUB(CharAE *, new_CharAE_from_string,
	(const char *string),
	(            string)
)

DEFINE_NOVALUE_CCALLABLE_STUB(CharAE_append_string,
	(CharAE *ae, const char *string),
	(        ae,             string)
)

DEFINE_NOVALUE_CCALLABLE_STUB(CharAE_delete_at,
	(CharAE *ae, size_t at, size_t nelt),
	(        ae,        at,        nelt)
)

DEFINE_CCALLABLE_STUB(SEXP, new_CHARSXP_from_CharAE,
	(const CharAE *ae),
	(              ae)
)

DEFINE_CCALLABLE_STUB(SEXP, new_RAW_from_CharAE,
	(const CharAE *ae),
	(              ae)
)

DEFINE_CCALLABLE_STUB(SEXP, new_LOGICAL_from_CharAE,
	(const CharAE *ae),
	(              ae)
)

DEFINE_CCALLABLE_STUB(size_t, CharAEAE_get_nelt,
	(const CharAEAE *aeae),
	(                aeae)
)

DEFINE_CCALLABLE_STUB(size_t, CharAEAE_set_nelt,
	(CharAEAE *aeae, size_t nelt),
	(          aeae,        nelt)
)

DEFINE_NOVALUE_CCALLABLE_STUB(CharAEAE_extend,
	(CharAEAE *aeae, size_t new_buflength),
	(          aeae,        new_buflength)
)

DEFINE_NOVALUE_CCALLABLE_STUB(CharAEAE_insert_at,
	(CharAEAE *aeae, size_t at, CharAE *ae),
	(          aeae,        at,         ae)
)

DEFINE_CCALLABLE_STUB(CharAEAE *, new_CharAEAE,
	(size_t buflength, size_t nelt),
	(       buflength,        nelt)
)

DEFINE_NOVALUE_CCALLABLE_STUB(CharAEAE_append_string,
	(CharAEAE *aeae, const char *string),
	(          aeae,             string)
)

DEFINE_CCALLABLE_STUB(SEXP, new_CHARACTER_from_CharAEAE,
	(const CharAEAE *aeae),
	(                aeae)
)

/*
 * Stubs for callables defined in SEXP_utils.c
 */

DEFINE_CCALLABLE_STUB(const char *, get_classname,
	(SEXP x),
	(     x)
)

/*
 * Stubs for callables defined in LLint_class.c
 */

DEFINE_CCALLABLE_STUB(int, is_LLint,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(R_xlen_t, get_LLint_length,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(long long int *, get_LLint_dataptr,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(SEXP, alloc_LLint,
	(const char *classname, R_xlen_t length),
	(            classname,          length)
)

/*
 * Stubs for callables defined in subsetting_utils.c
 */

DEFINE_CCALLABLE_STUB(long long int, copy_vector_block,
	(SEXP dest, long long int dest_offset, SEXP src, long long int src_offset, long long int block_nelt),
	(     dest,               dest_offset,      src,               src_offset,               block_nelt)
)

DEFINE_CCALLABLE_STUB(int, copy_vector_positions,
	(SEXP dest, int dest_offset, SEXP src, const int *pos, int npos),
	(     dest,     dest_offset,      src,            pos,     npos)
)

DEFINE_CCALLABLE_STUB(int, copy_vector_ranges,
	(SEXP dest, int dest_offset, SEXP src, const int *start, const int *width, int nranges),
	(     dest,     dest_offset,      src,            start,            width,     nranges)
)

/*
 * Stubs for callables defined in vector_utils.c
 */

DEFINE_CCALLABLE_STUB(int, vector_memcmp,
	(SEXP x1, int x1_offset, SEXP x2, int x2_offset, int nelt),
	(     x1,     x1_offset,      x2,     x2_offset,     nelt)
)

DEFINE_CCALLABLE_STUB(SEXP, list_as_data_frame,
	(SEXP x, int nrow),
	(     x,     nrow)
)

/*
 * Stubs for callables defined in integer_utils.c
 */

DEFINE_CCALLABLE_STUB(int, check_integer_pairs,
	(SEXP a, SEXP b, const int **a_p, const int **b_p, const char *a_argname, const char *b_argname),
	(     a,      b,             a_p,             b_p,             a_argname,             b_argname)
)

DEFINE_CCALLABLE_STUB(SEXP, find_interv_and_start_from_width,
	(const int *x, int x_len, const int *width, int width_len),
	(           x,     x_len,            width,     width_len)
)

/*
 * Stubs for callables defined in Hits_class.c
 */

DEFINE_CCALLABLE_STUB(SEXP, new_Hits,
	(const char *Class, int *from, const int *to, int nhit, int nLnode, int nRnode, int already_sorted),
	(            Class,      from,            to,     nhit,     nLnode,     nRnode,     already_sorted)
)

DEFINE_CCALLABLE_STUB(int, get_select_mode,
	(SEXP select),
	(     select)
)

/*
 * Stubs for callables defined in Rle_class.c
 */

DEFINE_CCALLABLE_STUB(SEXP, construct_logical_Rle,
	(R_xlen_t nrun_in, const int *values_in, const void *lengths_in, int lengths_in_is_L),
	(         nrun_in,            values_in,             lengths_in,     lengths_in_is_L)
)

DEFINE_CCALLABLE_STUB(SEXP, construct_integer_Rle,
	(R_xlen_t nrun_in, const int *values_in, const void *lengths_in, int lengths_in_is_L),
	(         nrun_in,            values_in,             lengths_in,     lengths_in_is_L)
)

DEFINE_CCALLABLE_STUB(SEXP, construct_numeric_Rle,
	(R_xlen_t nrun_in, const double *values_in, const void *lengths_in, int lengths_in_is_L),
	(         nrun_in,               values_in,             lengths_in,     lengths_in_is_L)
)

DEFINE_CCALLABLE_STUB(SEXP, construct_complex_Rle,
	(R_xlen_t nrun_in, const Rcomplex *values_in, const void *lengths_in, int lengths_in_is_L),
	(         nrun_in,                 values_in,             lengths_in,     lengths_in_is_L)
)

DEFINE_CCALLABLE_STUB(SEXP, construct_character_Rle,
	(SEXP values_in, const void *lengths_in, int lengths_in_is_L),
	(     values_in,             lengths_in,     lengths_in_is_L)
)

DEFINE_CCALLABLE_STUB(SEXP, construct_raw_Rle,
	(R_xlen_t nrun_in, const Rbyte *values_in, const void *lengths_in, int lengths_in_is_L),
	(         nrun_in,              values_in,             lengths_in,     lengths_in_is_L)
)

DEFINE_CCALLABLE_STUB(SEXP, construct_Rle,
	(SEXP values_in, const void *lengths_in, int lengths_in_is_L),
	(     values_in,             lengths_in,     lengths_in_is_L)
)

/*
 * Stubs for callables defined in List_class.c
 */

DEFINE_CCALLABLE_STUB(const char *, get_List_elementType,
	(SEXP x),
	(     x)
)

DEFINE_NOVALUE_CCALLABLE_STUB(set_List_elementType,
	(SEXP x, const char *type),
	(     x,             type)
)

/*
 * Stubs for callables defined in SimpleList_class.c
 */

DEFINE_CCALLABLE_STUB(SEXP, new_SimpleList,
	(const char *classname, SEXP listData),
	(            classname,      listData)
)

/*
 * Stubs for callables defined in DataFrame_class.c
 */

DEFINE_CCALLABLE_STUB(SEXP, new_DataFrame,
	(const char *classname, SEXP vars, SEXP rownames, SEXP nrows),
	(            classname,      vars,      rownames,      nrows)
)

