#include "IRanges_interface.h"

#define DEFINE_CCALLABLE_STUB(retT, stubname, Targs, args) \
typedef retT(*__ ## stubname ## _funtype__)Targs; \
retT stubname Targs \
{ \
	static __ ## stubname ## _funtype__ fun = NULL; \
	if (fun == NULL) \
		fun = (__ ## stubname ## _funtype__) R_GetCCallable("IRanges", "_" #stubname); \
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
		fun = (__ ## stubname ## _funtype__) R_GetCCallable("IRanges", "_" #stubname); \
	fun args; \
	return; \
}


/*
 * Stubs for callables defined in IPosRanges_comparison.c
 */

DEFINE_CCALLABLE_STUB(int, overlap_code,
	(int x_start, int x_width, int y_start, int y_width),
	(    x_start,     x_width,     y_start,     y_width)
)

DEFINE_CCALLABLE_STUB(int, invert_overlap_code,
	(int code),
	(    code)
)

/*
 * Stubs for callables defined in IRanges_class.c
 */

DEFINE_CCALLABLE_STUB(SEXP, get_IRanges_start,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(SEXP, get_IRanges_width,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(SEXP, get_IRanges_names,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(int, get_IRanges_length,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(IRanges_holder, hold_IRanges,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(int, get_length_from_IRanges_holder,
	(const IRanges_holder *x_holder),
	(                      x_holder)
)

DEFINE_CCALLABLE_STUB(int, get_width_elt_from_IRanges_holder,
	(const IRanges_holder *x_holder, int i),
	(                      x_holder,     i)
)

DEFINE_CCALLABLE_STUB(int, get_start_elt_from_IRanges_holder,
	(const IRanges_holder *x_holder, int i),
	(                      x_holder,     i)
)

DEFINE_CCALLABLE_STUB(int, get_end_elt_from_IRanges_holder,
	(const IRanges_holder *x_holder, int i),
	(                      x_holder,     i)
)

DEFINE_CCALLABLE_STUB(SEXP, get_names_elt_from_IRanges_holder,
	(const IRanges_holder *x_holder, int i),
	(                      x_holder,     i)
)

DEFINE_CCALLABLE_STUB(IRanges_holder, get_linear_subset_from_IRanges_holder,
	(const IRanges_holder *x_holder, int offset, int length),
	(                      x_holder,     offset,     length)
)

DEFINE_NOVALUE_CCALLABLE_STUB(set_IRanges_names,
	(SEXP x, SEXP names),
	(     x,      names)
)

DEFINE_NOVALUE_CCALLABLE_STUB(copy_IRanges_slots,
	(SEXP x, SEXP x0),
	(     x,      x0)
)

DEFINE_CCALLABLE_STUB(SEXP, new_IRanges,
	(const char *classname, SEXP start, SEXP width, SEXP names),
	(            classname,      start,      width,      names)
)

DEFINE_CCALLABLE_STUB(SEXP, new_IRanges_from_IntPairAE,
	(const char *classname, const IntPairAE *intpair_ae),
	(            classname,                  intpair_ae)
)

DEFINE_CCALLABLE_STUB(SEXP, new_list_of_IRanges_from_IntPairAEAE,
	(const char *element_type, const IntPairAEAE *intpair_aeae),
	(            element_type,                    intpair_aeae)
)

DEFINE_CCALLABLE_STUB(SEXP, alloc_IRanges,
	(const char *classname, int length),
	(            classname,     length)
)

/*
 * Stubs for callables defined in Grouping_class.c
 */

DEFINE_CCALLABLE_STUB(SEXP, get_H2LGrouping_high2low,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(SEXP, get_H2LGrouping_low2high,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(SEXP, get_Partitioning_names,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(SEXP, get_PartitioningByEnd_end,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(SEXP, new_PartitioningByEnd,
	(const char *classname, SEXP end, SEXP names),
	(            classname,      end,      names)
)

/*
 * Stubs for callables defined in CompressedList_class.c
 */

DEFINE_CCALLABLE_STUB(SEXP, get_CompressedList_unlistData,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(SEXP, get_CompressedList_partitioning,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(int, get_CompressedList_length,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(SEXP, get_CompressedList_names,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(SEXP, new_CompressedList,
	(const char *classname, SEXP unlistData, SEXP partitioning),
	(            classname,      unlistData,      partitioning)
)

DEFINE_CCALLABLE_STUB(CompressedIntsList_holder, hold_CompressedIntegerList,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(int, get_length_from_CompressedIntsList_holder,
	(const CompressedIntsList_holder *x_holder),
	(                                 x_holder)
)

DEFINE_CCALLABLE_STUB(Ints_holder, get_elt_from_CompressedIntsList_holder,
	(const CompressedIntsList_holder *x_holder, int i),
	(                                 x_holder,     i)
)

/*
 * Stubs for callables defined in CompressedIRangesList_class.c
 */

DEFINE_CCALLABLE_STUB(CompressedIRangesList_holder, hold_CompressedIRangesList,
	(SEXP x),
	(     x)
)

DEFINE_CCALLABLE_STUB(int, get_length_from_CompressedIRangesList_holder,
	(const CompressedIRangesList_holder *x_holder),
	(                                    x_holder)
)

DEFINE_CCALLABLE_STUB(IRanges_holder, get_elt_from_CompressedIRangesList_holder,
	(const CompressedIRangesList_holder *x_holder, int i),
	(                                    x_holder,     i)
)

DEFINE_CCALLABLE_STUB(int, get_eltNROWS_from_CompressedIRangesList_holder,
	(const CompressedIRangesList_holder *x_holder, int i),
	(                                    x_holder,     i)
)

