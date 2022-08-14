#include "header.h"
#include "proto.h"

void random_walks(
 int *ref_image,
 double *v_arr,
 int *alph_arr,
 int width,
 int height,
 int beta,
 int maxiter
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
 LIS_MATRIX matA;
 LIS_VECTOR vecb;
 LIS_VECTOR vecx;
 LIS_VECTOR vecy;
 LIS_SOLVER solver;
 LIS_INT nnz_row;
 LIS_INT ierr;
 LIS_SCALAR value;
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
 LIS_INT argc_duh;
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
 char option_string[80];

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

 argc_duh= 1;
 argv_duh= (char **)calloc(1,sizeof(char *));
 argv_duh[0]= (char *)calloc(80,sizeof(char));
 argv_duh[0]="duh";
 ierr= lis_initialize(&argc_duh,&argv_duh);

 ierr= lis_solver_create(&solver);
 sprintf(option_string,"-i sor -p none -maxiter %d -print out",maxiter);
 ierr= lis_solver_set_option(option_string,solver);

 ierr= lis_matrix_create(0,&matA);
 ierr= lis_matrix_set_size(matA,0,width*height);

 nnz_row= window_size;

 /*
 Although this is the exact number of non zero terms per row,
 lis still wants to call lis_realloc when setting values
 Apparently, it's a very good idea to over-estimate nnz
 */

 nnz_row*= 2;

 ierr= lis_matrix_malloc(matA,nnz_row,0);

 ierr= lis_vector_create(0,&vecx);
 ierr= lis_vector_set_size(vecx,0,width*height);

 ierr= lis_vector_create(0,&vecb);
 ierr= lis_vector_set_size(vecb,0,width*height);

 ierr= lis_vector_create(0,&vecy);
 ierr= lis_vector_set_size(vecy,0,width*height);

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
          Solution is obviously already known
          */

          lis_matrix_set_value(LIS_INS_VALUE,ind,ind,1,matA);
          lis_vector_set_value(LIS_INS_VALUE,ind,v,vecb);

          continue;
       }

       /*
       If here, this is not a seed
       */

       if ( alph != 0 ) {
          error_handler("random_walks");
       } 

       /*
       Initialize
       sigma = sum of the weights at all edges
       sigmav2 = sum of the weights (times voltage) at only the edges
       that connect to a seed pixel (input voltage)
       */

       sigma= 0;
       sigmav2= 0;

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
             if ( max_dist > 0 )
              dist/= max_dist;
             /*
             Instead of using an exponential like Grady,
             we're gonna use a sigmoid
             */
             /*
             G= exp(-(double)beta*dist);
             */
             G= 2/(1+exp((double)beta*dist));

             sigma+= G;

             if ( alph2 == 0 ) {

                /*
                This is not a seed
                */

                lis_matrix_set_value(LIS_INS_VALUE,ind,ind2,-G,matA);
             }
             else if ( alph2 == 255 ) {

                /*
                This is a seed
                */

                sigmav2+= G*v2; 
             }
             else {
                error_handler("random_walks");
             }
          }
       }

       lis_matrix_set_value(LIS_INS_VALUE,ind,ind,sigma,matA);
       lis_vector_set_value(LIS_INS_VALUE,ind,sigmav2,vecb);
    }
 }

 ierr= lis_matrix_set_type(matA,LIS_MATRIX_CSR);
 ierr= lis_matrix_assemble(matA);

 /*
 ierr= lis_output_matrix(matA,LIS_FMT_MM,"matA.txt");
 ierr= lis_output_vector(vb,LIS_FMT_MM,"vecb.txt");
 */

 ierr= lis_solve(matA,vecb,vecx,solver);

 /*
 ierr= lis_output_vector(vecx,LIS_FMT_MM,"vecx.txt");
 */

 /*
 Let's compute y = Ax-b
 to see if convergence was actually obtained
 */

 /* uncomment if need be
 ierr= lis_matvec(matA,vecx,vecy);
 ierr= lis_vector_axpy(-1,vecb,vecy);

 ierr= lis_output_vector(vecy,LIS_FMT_MM,"vecy.txt");
 */ 

 /*
 Update v_arr
 */

 for ( i= 0 ; i< height ; i++ ) {
    for ( j= 0 ; j< width ; j++ ) {
       ind= i*width+j;
       v= v_arr[ind];
       alph= alph_arr[ind];
       ierr= lis_vector_get_value(vecx,ind,&value);

       if ( alph != 0 ) { 

          /*
          This is a seed
          */

          continue;
       } 

       v_arr[ind]= value;
    }
 }

 ierr= lis_matrix_destroy(matA);
 ierr= lis_vector_destroy(vecx);
 ierr= lis_vector_destroy(vecb);
 ierr= lis_vector_destroy(vecy);
 ierr= lis_solver_destroy(solver);

 ierr= lis_finalize();

}
