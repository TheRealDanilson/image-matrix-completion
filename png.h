#ifndef PNG_CUSTOM_H
#define PNG_CUSTOM_H
#include <png.h>
typedef struct {
    png_uint_32 height;
    png_uint_32 width;
    png_bytepp rows;
} image;

image read_png(char* filename);

int write_png(image img, char* filename);
#endif