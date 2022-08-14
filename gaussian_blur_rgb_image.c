#include "header.h"
#include "proto.h"

void gaussian_blur_rgb_image(
 int *inp_image_arr,
 int xdim,
 int ydim,
 double sigma,
 int precision,
 int *out_image_arr
)

{

 int *I;
 int i;
 int j;
 int ind;
 int *I_out;
 int cind;

 /*
 Allocate memory to store image intensity for each channel
 */

 I= (int *)calloc(xdim*ydim,sizeof(int));

 /*
 Allocate memory to store blurred image intensity for each channel
 */

 I_out= (int *)calloc(xdim*ydim,sizeof(int));

 /*
 Process one channel at a time
 */

 for ( cind= 0 ; cind< 3 ; cind++ ) {

    for ( i= 0 ; i< ydim ; i++ ) {
       for ( j= 0 ; j< xdim ; j++ ) {
          ind= i*xdim+j;
          I[ind]= inp_image_arr[3*ind+cind];
       }
    }

    gaussian_blur_image(
     I,
     xdim,
     ydim,
     sigma,
     precision,
     I_out
    );

    for ( i= 0 ; i< ydim ; i++ ) {
       for ( j= 0 ; j< xdim ; j++ ) {
          ind= i*xdim+j;
          out_image_arr[3*ind+cind]= I_out[ind];
       }
    }
 }

 /*
 Free memory to store image intensity for each channel
 */

 free(I);

 /*
 Free memory to store blurred image intensity for each channel
 */

 free(I_out);

}
