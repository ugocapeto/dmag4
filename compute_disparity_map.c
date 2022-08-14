#include "header.h"
#include "proto.h"

void compute_disparity_map(
 int *I,
 int *disp_rgb_arr,
 int *disp_alph_arr,
 int *edge_alph_arr,
 int width,
 int height,
 int beta,
 int maxiter,
 int scale_nbr,
 int con_level,
 int con_level2,
 int **pdisp_arr
)

{

 int *disp_arr;
 int i;
 int j;
 int ind;
 int r;
 int g;
 int b;
 double *v_arr;
 double dval;
 int intensity;
 int pixel;
 int disp;
 int ngb_disp;
 int max_ngb_disp;
 int i2;
 int j2;
 int ngb_i;
 int ngb_j;
 int ngb_pixel;

 disp_arr= (int *)calloc(width*height,sizeof(int));

 for ( i= 0 ; i< height ; i++ ) {
    for ( j= 0 ; j< width ; j++ ) {
       ind= i*width+j;
       r= disp_rgb_arr[3*ind+0];
       g= disp_rgb_arr[3*ind+1];
       b= disp_rgb_arr[3*ind+2];
       intensity= .2989*(double)r+.5870*(double)g+.1140*(double)b;
       if ( intensity < 0 )
        intensity= 0;
       if ( intensity > 255 )
        intensity= 255;
       disp_arr[ind]= intensity;
    }
 }

 /*
 Transform the disparities into voltage
 */

 v_arr= (double *)calloc(width*height,sizeof(double));

 for ( i= 0 ; i< height ; i++ ) {
    for ( j= 0 ; j< width ; j++ ) {
       ind= i*width+j;
       v_arr[ind]= (double)disp_arr[ind]/255.;
    }
 }

 /*
 Ready to use the Random Walk algorithm on the disparities
 */

 random_walks_scale_space(
  I,
  v_arr,
  disp_alph_arr,
  edge_alph_arr,
  width,
  height,
  beta,
  maxiter,
  scale_nbr,
  con_level,
  con_level2
 );

 /*
 Transform the voltage back into disparities
 */

 for ( i= 0 ; i< height ; i++ ) {
    for ( j= 0 ; j< width ; j++ ) {
       ind= i*width+j;
       dval= v_arr[ind]*255.;
       disp_arr[ind]= (int)(dval+.5);
       if ( disp_arr[ind] < 0 )
        disp_arr[ind]= 0;
       if ( disp_arr[ind] > 255 )
        disp_arr[ind]= 255;
    }
 }

 free(v_arr);

 /*
 Any pixel in the edge image
 should take the disparity of its neighbor
 if the neighbor is not in the edge image and is more in the foreground
 */

 for ( i= 0 ; i< height ; i++ ) {
    for ( j= 0 ; j< width ; j++ ) {
       pixel= i*width+j;
       if ( edge_alph_arr[pixel] == 0 )
        continue;
       if ( edge_alph_arr[pixel] != 255 ) {
          error_handler("compute_disparity_map"); 
       }
       disp= disp_arr[pixel];
       for ( i2= -1 ; i2<= +1 ; i2++ ) {
          for ( j2= -1 ; j2<= +1 ; j2++ ) {
             ngb_i= i+i2;
             ngb_j= j+j2;
             if ( ngb_i < 0 )
              continue;
             if ( ngb_i > height-1 )
              continue;
             if ( ngb_j < 0 )
              continue;
             if ( ngb_j > width-1 )
              continue;
             ngb_pixel= ngb_i*width+ngb_j; 
             if ( ngb_pixel == pixel )
              continue;
             if ( edge_alph_arr[ngb_pixel] == 255 )
              continue;
             if ( edge_alph_arr[ngb_pixel] != 0 ) {
                error_handler("compute_disparity_map");
             }
             ngb_disp= disp_arr[ngb_pixel];
             if ( !(ngb_disp > disp) )
              continue;
             disp_arr[pixel]= ngb_disp;
             disp= disp_arr[pixel];
          }
       }
    }
 }

 (*pdisp_arr)= disp_arr;

}
