#include <png.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "png.h"

image read_png(char* filename) {
    png_structp png_ptr;
    png_infop info_ptr;
    png_uint_32 width, height;
    int bit_depth, color_type, interlace_type;

    FILE *fp;
    if ((fp = fopen(filename, "rb")) == NULL) {
        fprintf(stderr, "something went opening file %s\n", filename);
        exit(1);
    }

    png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

    if (png_ptr == NULL) {
        fclose(fp);
        fprintf(stderr, "something went creating read struct\n");
        exit(1); 
    }

    info_ptr = png_create_info_struct(png_ptr);
    if (info_ptr == NULL) {
        fclose(fp);
        fprintf(stderr, "something went creating info struct\n");
        png_destroy_read_struct(&png_ptr, NULL, NULL);
        exit(1);
    }

    if (setjmp(png_jmpbuf(png_ptr))) {
      png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
      fclose(fp);
      fprintf(stderr, "something went wrong reading the png file\n");
      exit(1);
   }

   png_init_io(png_ptr, fp);

   png_read_info(png_ptr, info_ptr);
   png_get_IHDR(png_ptr, info_ptr, &width, &height, &bit_depth, &color_type, &interlace_type, NULL, NULL);

    if (color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
        png_set_gray_to_rgb(png_ptr);
    }
    png_set_strip_alpha(png_ptr);
    png_set_scale_16(png_ptr);
    //png_set_alpha_mode(png_ptr, PNG_ALPHA_STANDARD, PNG_DEFAULT_sRGB);

    png_read_update_info(png_ptr, info_ptr);
    png_bytepp row_pointers = calloc(height, sizeof(png_bytep));
    for (png_uint_32 row = 0; row < height; row++) {
        row_pointers[row] = png_malloc(png_ptr, png_get_rowbytes(png_ptr, info_ptr));
    }
    png_read_image(png_ptr, row_pointers);
    png_read_end(png_ptr, info_ptr);

    fclose(fp);

    png_destroy_read_struct(&png_ptr, &info_ptr, NULL);
    image img = {.height = height, .width = width, .rows = row_pointers};
    return img;
}

int write_png(image img, char* filename) {
    png_structp w_png_ptr;
    png_infop w_info_ptr;
    FILE *w_fp;
    if ((w_fp = fopen(filename, "wb")) == NULL) {
        fprintf(stderr, "something went opening file %s\n", filename);
        exit(1);
    }
    w_png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (w_png_ptr == NULL) {
        fclose(w_fp);
        fprintf(stderr, "something went creating write struct\n");
        exit(1); 
    }

    w_info_ptr = png_create_info_struct(w_png_ptr);
    if (w_info_ptr == NULL) {
        fclose(w_fp);
         png_destroy_write_struct(&w_png_ptr, NULL);
        fprintf(stderr, "something went wrong creating the info struct\n");
        exit(1);
    }

if (setjmp(png_jmpbuf(w_png_ptr)))
   {
      fclose(w_fp);
      png_destroy_write_struct(&w_png_ptr, &w_info_ptr);
      fprintf(stderr, "something went wrong writing the png file\n");
      exit(1);
   }

   png_init_io(w_png_ptr, w_fp);
   png_set_IHDR(w_png_ptr, w_info_ptr, img.width, img.height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

   png_write_info(w_png_ptr, w_info_ptr);

   png_write_image(w_png_ptr, img.rows);
   png_write_end(w_png_ptr, w_info_ptr);

   fclose(w_fp);

   png_destroy_write_struct(&w_png_ptr, &w_info_ptr);

   return 0;
}