#include "header.h"
#include "proto.h"

void random_walks_sor(
 int *ref_image,
 double *v_arr,
 int *alph_arr,
 int width,
 int height,
 int beta
)

/*
See chapter 3 of
"Image segmentation through the scale-space random walker"
by Richard Rzeszutek, Ryerson University
and
"Semi-automatic 2d to 3d image conversion
using a hybrid random walks and graph cuts based approad"
by Raymond Phan, Richard Rzeszutek, and Dimitrios Androutsos
and
"Image segmentation uisng scale-space random walks"
by Richard Rzeszutek, Thomas El-Maraghi, and Dimitrios Androutsos
*/

{

 double *ref_image_Lab;
 int mask;
 int window_size;
 double sigma;
 int i;
 int j;
 int ind;
 double CIEL;
 double CIEa;
 double CIEb;
 double CIEL2;
 double CIEa2;
 double CIEb2;
 double sigmav2;
 int mask_i;
 int mask_j;
 int i2;
 int j2;
 int ind2;
 double v2;
 double dist2;
 double dist;
 double G;
 double alph2;
 double v;
 double alph;
 double max_dist;
 char **argv_duh;
 int r;
 int g;
 int b;
 double x;
 double y;
 double z;
 int r2;
 int g2;
 int b2;
 double x2;
 double y2;
 double z2;
 int iter;
 double *new_v_arr;
 double GSIGMA;
 double GU;
 double GS;
 double new_v;
 double max_diff;
 double diff;
 double tol;
 double sum_diff;
 double sor_w;

 mask= 1;
 window_size= (2*mask+1)*(2*mask+1);

 /*
 Figure out max dist
 so that all dist can be normalized from 0 to 1
 */

 max_dist= 0.0;

 for ( i= 0 ; i< height ; i++ ) {
    for ( j= 0 ; j< width ; j++ ) {
       ind= i*width+j;
       r= ref_image[3*ind+0];
       g= ref_image[3*ind+1];
       b= ref_image[3*ind+2];
       rgb2xyz(r,g,b,&x,&y,&z);
       xyz2Lab(x,y,z,&CIEL,&CIEa,&CIEb);

       for ( mask_i= -mask ; mask_i<= mask ; mask_i++ ) {
          for ( mask_j= -mask ; mask_j<= mask ; mask_j++ ) {
             if ( mask_i == 0 && mask_j == 0 ) {
                continue;
             }
             i2= i+mask_i;
             j2= j+mask_j;
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
             r2= ref_image[3*ind2+0];
             g2= ref_image[3*ind2+1];
             b2= ref_image[3*ind2+2];
             rgb2xyz(r2,g2,b2,&x2,&y2,&z2);
             xyz2Lab(x2,y2,z2,&CIEL2,&CIEa2,&CIEb2);
             dist2= (CIEL2-CIEL)*(CIEL2-CIEL)+
                    (CIEa2-CIEa)*(CIEa2-CIEa)+
                    (CIEb2-CIEb)*(CIEb2-CIEb);
             dist= sqrt(dist2);
             if ( dist > max_dist ) {
                max_dist= dist;
             }
          }
       }
    }
 }

 /*
 We're gonna use an SOR scheme
 */

 iter= 0;

 SOR_ITER:

 iter++;

 fprintf(stdout,"SOR iteration = %d\n",iter);

 /*
 Allocate memory for the new voltages
 */

 new_v_arr= (double *)calloc(width*height,sizeof(double));

 /*
 Initialize the new voltage array with the current voltage array
 */

 for ( i= 0 ; i< height ; i++ ) {
    for ( j= 0 ; j< width ; j++ ) {
       ind= i*width+j;
       new_v_arr[ind]= v_arr[ind];
    }
 }

 /*
 Let's do the SOR iteration
 */

 for ( i= 0 ; i< height ; i++ ) {
    for ( j= 0 ; j< width ; j++ ) {
       ind= i*width+j;
       v= v_arr[ind];
       alph= alph_arr[ind];
       r= ref_image[3*ind+0];
       g= ref_image[3*ind+1];
       b= ref_image[3*ind+2];
       rgb2xyz(r,g,b,&x,&y,&z);
       xyz2Lab(x,y,z,&CIEL,&CIEa,&CIEb);

       if ( alph == 255 ) {

          /*
          This is a seed
          */

          continue;
       }

       /*
       If here, this is not a seed
       */

       if ( alph != 0 ) {
          error_handler("random_walks");
       }

       GSIGMA= 0.0;
       GU= 0.0;
       GS= 0.0;

       for ( mask_i= -mask ; mask_i<= mask ; mask_i++ ) {
          for ( mask_j= -mask ; mask_j<= mask ; mask_j++ ) {
             if ( mask_i == 0 && mask_j == 0 ) {
                continue;
             } 
             i2= i+mask_i;
             j2= j+mask_j;
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
             ind2= i2*width+j2;
             r2= ref_image[3*ind2+0];
             g2= ref_image[3*ind2+1];
             b2= ref_image[3*ind2+2];
             rgb2xyz(r2,g2,b2,&x2,&y2,&z2);
             xyz2Lab(x2,y2,z2,&CIEL2,&CIEa2,&CIEb2);
             dist2= (CIEL2-CIEL)*(CIEL2-CIEL)+
                    (CIEa2-CIEa)*(CIEa2-CIEa)+
                    (CIEb2-CIEb)*(CIEb2-CIEb);
             dist= sqrt(dist2);
             /*
             Normalize the distance in the color space
             */
             dist/= max_dist;
             /*
             Compute the conductance
             */
             /*
             Instead of using an exponential like Grady,
             we're gonna use a sigmoid
             */
             /*
             G= exp(-(double)beta*dist);
             */
             G= 2/(1+exp((double)beta*dist));

             GSIGMA+= G;

             if ( alph2 == 0 ) {

                /*
                This is not a seed
                */

                GU+= G*v2;
             }
             else if ( alph2 == 255 ) {

                /*
                This is a seed
                */

                GS+= G*v2;
             }
             else {
                error_handler("random_walks");
             }
          }
       }

       /*
       Compute the new voltage
       using the SOR relaxation parameter
       */

       /*
       sor_w= 1.9;
       */
       sor_w= 1.0;
       new_v= (1-sor_w)*v+sor_w*(GS+GU)/GSIGMA;

       /*
       Store the new voltage in the new voltage array
       */

       new_v_arr[ind]= new_v;
    }
 }

 /*
 Compute the maximum variation (from current voltage to new voltage)
 so that we know when to stop iterating
 */

 /*
 Recall that voltages vary from 0 to 1
 */

 max_diff= 0.0;
 sum_diff= 0.0;
 for ( i= 0 ; i< height ; i++ ) {
    for ( j= 0 ; j< width ; j++ ) {
       ind= i*width+j;
       v= v_arr[ind];
       new_v= new_v_arr[ind];
       diff= fabs(new_v-v);
       if ( diff > max_diff )
        max_diff= diff;
       sum_diff+= diff;
    }
 }
 sum_diff/= (double)(width*height);

 fprintf(stdout,"Average voltage variation = %f\n",sum_diff);
 fprintf(stdout,"Max voltage variation = %f\n",max_diff);

 /*
 Update v_arr
 */

 for ( i= 0 ; i< height ; i++ ) {
    for ( j= 0 ; j< width ; j++ ) {
       ind= i*width+j;
       v_arr[ind]= new_v_arr[ind];
    }
 }

 /*
 Free the new voltage array
 */

 free(new_v_arr);

 /*
 Check for convergence
 before doing another SOR iteration
 */

 tol= 0.01;
 if ( max_diff > tol )
  goto SOR_ITER;

}
