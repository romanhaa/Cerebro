/*****************************************************************************
 S4Vectors C interface: typedefs and defines
 -------------------------------------------

   The S4Vectors C interface is split in 2 files:
     1. S4Vectors_defines.h (this file): contains the typedefs and defines
        of the interface.
     2. S4Vectors_interface.h (in this directory): contains the prototypes
        of the S4Vectors C routines that are part of the interface.

   Please consult S4Vectors_interface.h for how to use this interface in your
   package.

 *****************************************************************************/
#ifndef S4VECTORS_DEFINES_H
#define S4VECTORS_DEFINES_H

#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#define	NA_LLINT	LLONG_MIN

/* Get or set i-th element from int or long long int array 'x'.
   GET_INT_OR_LLINT() always returns a long long int.
   SET_INT_OR_LLINT() can take a value 'v' that is int or long long int but
   is not safe if 'is_llint' is 0 and 'v' is a long long int. */
#define	GET_INT_OR_LLINT(x, is_llint, i)			\
	((is_llint) ? ((const long long int *)(x))[i]		\
		    : (long long int) ((const int *)(x))[i])
#define	SET_INT_OR_LLINT(x, is_llint, i, v)			\
{								\
        if (is_llint)						\
                ((long long int *)(x))[i] = (v);		\
        else							\
                ((int *)(x))[i] = (v);				\
}

/* Hash table -- modified from R_HOME/src/main/unique.c */
struct htab {
	int K, M;
	unsigned int Mminus1;
	int *buckets;
};


/*
 * Auto-Extending buffers used for temporary storage of incoming data whose
 * size is not known in advance:
 *
 *   o IntAE:       Auto-Extending buffer of ints;
 *   o IntAEAE:     Auto-Extending buffer of Auto-Extending buffers of ints;
 *   o IntPairAE:   Auto-Extending buffer of pairs of ints;
 *   o IntPairAEAE: Auto-Extending buffer of Auto-Extending buffers of pairs
 *                  of ints;
 *   o LLongAE:     Auto-Extending buffer of long long ints;
 *   o CharAE:      Auto-Extending buffer of chars;
 *   o CharAEAE:    Auto-Extending buffer of Auto-Extending buffers of chars.
 *
 * Some differences between AE buffers and SEXP: (a) AE buffers auto-extend
 * i.e. they automatically reallocate when more room is needed to add a new
 * element, (b) they are faster, and (c) they don't require any
 * PROTECT/UNPROTECT mechanism.
 */

typedef struct int_ae {
	size_t _buflength;
	size_t _nelt;
	int *elts;
} IntAE;

typedef struct int_aeae {
	size_t _buflength;
	size_t _nelt;
	IntAE **elts;
} IntAEAE; 

typedef struct intpair_ae {
	IntAE *a;
	IntAE *b;
} IntPairAE;

typedef struct intpair_aeae {
	size_t _buflength;
	size_t _nelt;
	IntPairAE **elts;
} IntPairAEAE;

typedef struct longlong_ae {
	size_t _buflength;
	size_t _nelt;
	long long int *elts;
} LLongAE;

typedef struct char_ae {
	size_t _buflength;
	size_t _nelt;
	char *elts;
} CharAE; 

typedef struct char_aeae {
	size_t _buflength;
	size_t _nelt;
	CharAE **elts;
} CharAEAE; 


/*
 * Holder structs.
 */

typedef struct chars_holder {
	const char *ptr;
	int length;
} Chars_holder;

typedef struct ints_holder {
	const int *ptr;
	int length;
} Ints_holder;

typedef struct doubles_holder {
	const double *ptr;
	int length;
} Doubles_holder;


/*
 * Hit selection modes.
 */

#define ALL_HITS	1
#define FIRST_HIT	2
#define LAST_HIT	3
#define ARBITRARY_HIT	4
#define COUNT_HITS	5

#endif
