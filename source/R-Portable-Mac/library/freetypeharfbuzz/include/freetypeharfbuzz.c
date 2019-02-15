#define R_NO_REMAP
#include <Rinternals.h>

#include <freetypeharfbuzz.h>


int (*fthb_get_font_info)(const char*, double, struct fthb_font_info*) = NULL;
int (*fthb_calc_string_info)(const char*, const char*, double, struct fthb_string_info*) = NULL;
int (*fthb_calc_string_width)(const char*, const char*, double, double*) = NULL;

void fthb_init() {
  fthb_get_font_info =
    (int (*)(const char*, double, struct fthb_font_info*))
    R_GetCCallable("freetypeharfbuzz", "fthb_get_font_info");

  fthb_calc_string_info =
    (int (*)(const char*, const char*, double, struct fthb_string_info*))
    R_GetCCallable("freetypeharfbuzz", "fthb_calc_string_info");

  fthb_calc_string_width =
    (int (*)(const char*, const char*, double, double*))
    R_GetCCallable("freetypeharfbuzz", "fthb_calc_string_width");
}
