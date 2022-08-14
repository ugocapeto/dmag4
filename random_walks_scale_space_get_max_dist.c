#include "header.h"
#include "proto.h"

void random_walks_scale_space_get_max_dist(
 int **ref_image,
 double *v_arr,
 int *alph_arr,
 int width,
 int height,
 int scale_nbr,
 int color_space,
 int radius,
 int con_level,
 int con_level2,
 double *pmax_dist
)

/*
Copy of random_walks_scale_space
to get maximum distance for all edges in graph
*/

{

 int scale_ind;
 int i;
 int j;
 int ind;
 double v;
 int alph;
 int r;
 int g;
 int b;
 double x;
 double y;
 double z;
 double CIEL;
 double CIEa;
 double CIEb;
 int radius_i;
 int radius_j;
 int i2;
 int j2;
 int ind2;
 double v2;
 int alph2;
 int r2;
 int g2;
 int b2;
 double x2;
 double y2;
 double z2;
 double CIEL2;
 double CIEa2;
 double CIEb2;
 int dist2_int;
 double dist2;
 double dist;
 double max_dist;
 int scale_ind2;

 /*
 Compute the max distance in color space
 so that the distances can be normalized between 0 and 1
 considering all edges in the graph
 */

 max_dist= 0;

 for ( scale_ind= 0 ; scale_ind< scale_nbr ; scale_ind++ ) {

    for ( i= 0 ; i< height ; i++ ) {
       for ( j= 0 ; j< width ; j++ ) {
          ind= i*width+j;
          v= v_arr[ind];
          alph= alph_arr[ind];
          r= ref_image[scale_ind][3*ind+0];
          g= ref_image[scale_ind][3*ind+1];
          b= ref_image[scale_ind][3*ind+2];
          if ( color_space == 1 ) {
             rgb2xyz(r,g,b,&x,&y,&z);
             xyz2Lab(x,y,z,&CIEL,&CIEa,&CIEb);
          }

          if ( alph == 255 ) {

             /*
             This is a seed
             Solution is obviously already known
             */

             /*
             No edges to consider for that pixel
             */

             continue;
          }

          /*
          If here, this is not a seed
          */

          if ( alph != 0 ) {
             error_handler("random_walks");
          }

          /*
          Consider previous, current and next scale
          */

          for ( scale_ind2= scale_ind-1 ;
                scale_ind2<= scale_ind+1 ;
                scale_ind2++ ) {
             if ( !(scale_ind2 >= 0) )
              continue;
             if ( !(scale_ind2 < scale_nbr) )
              continue;

             /*
             Consider neighboring pixels
             */

             for ( radius_i= -radius ; radius_i<= radius ; radius_i++ ) {
                for ( radius_j= -radius ; radius_j<= radius ; radius_j++ ) {
                   if ( scale_ind2 == scale_ind ) {
                      if ( con_level == 1 ) {
                         if ( abs(radius_j) == abs(radius_i)  ) {
                            continue;
                         }
                      }
                      else if ( con_level == 2 ) {
                         if ( radius_i == 0 && radius_j == 0 ) {
                            continue;
                         }
                      }
                      else {
                         error_handler("random_walks_scale_space");
                      }
                   }
                   else {
                      if ( con_level2 == 1 ) {
                         if ( !(radius_i == 0 && radius_j == 0) ) {
                            continue;
                         }
                      }
                      else if ( con_level2 == 2 ) {
                         if ( abs(radius_j) == abs(radius_i) &&
                              !(radius_i == 0 && radius_j == 0) ) {
                            continue;
                         }
                      }
                      else if ( con_level2 == 3 ) {
                      }
                      else {
                         error_handler("random_walks_scale_space");
                      }
                   }

                   i2= i+radius_i;
                   j2= j+radius_j;
                   if ( i2 < 0 ) {
                      continue;
                   }
                   if ( i2 >= height ) {
                      continue;
                   }
                   if ( j2 < 0 ) {
                      continue;
                   }
                   if ( j2 >= width ) {
                      continue;
                   }
                   ind2= i2*width+j2;
                   v2= v_arr[ind2];
                   alph2= alph_arr[ind2];
                   r2= ref_image[scale_ind2][3*ind2+0];
                   g2= ref_image[scale_ind2][3*ind2+1];
                   b2= ref_image[scale_ind2][3*ind2+2];
                   if ( color_space == 1 ) {
                      rgb2xyz(r2,g2,b2,&x2,&y2,&z2);
                      xyz2Lab(x2,y2,z2,&CIEL2,&CIEa2,&CIEb2);
                   }
                   dist2_int= 0;
                   dist2_int+= (r2-r)*(r2-r);
                   dist2_int+= (g2-g)*(g2-g);
                   dist2_int+= (b2-b)*(b2-b);
                   dist2= (double)dist2_int;
                   dist= sqrt(dist2);
                   if ( color_space == 1 ) {
                      dist2= (CIEL2-CIEL)*(CIEL2-CIEL)+
                             (CIEa2-CIEa)*(CIEa2-CIEa)+
                             (CIEb2-CIEb)*(CIEb2-CIEb);
                      dist= sqrt(dist2);
                   }

                   if ( dist > max_dist ) {
                      max_dist= dist;
                   }
                }
             }
          } /* loop on previous, current, and next scale */
       } /* loop on j */
    } /* loop on i */
 } /* loop on scale_ind */

 (*pmax_dist)= max_dist;

}
