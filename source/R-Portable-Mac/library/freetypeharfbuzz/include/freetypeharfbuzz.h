#ifndef FREETYPEHARFBUZZ_H
#define FREETYPEHARFBUZZ_H


struct fthb_font_info {
  double ascent;
  double descent;
  double linegap;
};

struct fthb_string_info {
  double width;
  double height;
  double ascent;
  double descent;
};

extern
int (*fthb_get_font_info)(const char* font_path,
                          double font_size,
                          struct fthb_font_info* metrics_out);

extern
int (*fthb_calc_string_info)(const char* string,
                             const char* font_path,
                             double font_size,
                             struct fthb_string_info* metrics_out);

extern
int (*fthb_calc_string_width)(const char* string,
                              const char* font_path,
                              double font_size,
                              double* width_out);


#endif
