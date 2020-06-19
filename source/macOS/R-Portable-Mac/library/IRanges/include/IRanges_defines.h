/*****************************************************************************
 IRanges C interface: typedefs and defines
 -----------------------------------------

   The IRanges C interface is split in 2 files:
     1. IRanges_defines.h (this file): contains the typedefs and defines
        of the interface.
     2. IRanges_interface.h (in this directory): contains the prototypes
        of the IRanges C routines that are part of the interface.

   Please consult IRanges_interface.h for how to use this interface in your
   package.

 *****************************************************************************/
#ifndef IRANGES_DEFINES_H
#define IRANGES_DEFINES_H

#include "S4Vectors_defines.h"

#include <Rdefines.h>
#include <R_ext/Rdynload.h>


/*
 * *_holder structs.
 */

typedef struct compressed_chars_list_holder {
	int length;
	const char *unlisted;
	const int *breakpoints;
} CompressedCharsList_holder;

typedef struct compressed_ints_list_holder {
	int length;
	const int *unlisted;
	const int *breakpoints;
} CompressedIntsList_holder;

typedef struct compressed_doubles_list_holder {
	int length;
	const double *unlisted;
	const int *breakpoints;
} CompressedDoublesList_holder;

typedef struct iranges_holder {
	const char *classname;
	int is_constant_width;
	int length;
	const int *width;
	const int *start;
	const int *end;
	int SEXP_offset;  /* offset in 'names' member below */
	SEXP names;
} IRanges_holder;

typedef struct compressed_iranges_list_holder {
	const char *classname;
	int length;
	const int *end;
	IRanges_holder unlistData_holder;
} CompressedIRangesList_holder;

#endif
