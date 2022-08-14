#include "header.h"
#include "proto.h"

void random_walks_scale_space(
 int *inp_ref_image,
 double *v_arr,
 int *disp_alph_arr,
 int *edge_alph_arr,
 int width,
 int height,
 int beta,
 int maxiter,
 int scale_nbr,
 int con_level,
 int con_level2
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
 int radius;
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
 int radius_i;
 int radius_j;
 int i2;
 int j2;
 int ind2;
 double v2;
 double dist2;
 int dist2_int;
 double dist;
 double G;
 double disp_alph2;
 double v;
 double disp_alph;
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
 int **ref_image;
 int scale_ind;
 int precision;
 int matA_ind;
 int nbr_pixels;
 int vecb_ind;
 int sigma_n;
 int matA_ind2;
 int scale_ind2;
 int vecx_ind;
 char filename[80];
 int cind;
 int color_space;
 int edge_alph;
 int edge_alph2;

 /*
 Default color space is rgb (0)
 */

 color_space= 0;

 /*
 Switch color space to CIELab
 Comment out if don't want that
 */

 color_space= 1;

 /*
 Set the radius
 */

 radius= 1;

 /*
 Set the number of pixels
 */

 nbr_pixels= width*height;

 /*
 Allocate memory to hold the reference images at each scale
 */

 ref_image= (int **)calloc(scale_nbr,sizeof(int *));

 /*
 Generate the scale space
 */

 /*
 The input image is always at scale 0
 */

 scale_ind= 0;

 /*
 Allocate memory to hold the reference image at that scale
 */

 ref_image[scale_ind]= (int *)calloc(3*width*height,sizeof(int));

 for ( i= 0 ; i< height ; i++ ) {
    for ( j= 0 ; j< width ; j++ ) {
       ind= i*width+j;
       for ( cind= 0 ; cind< 3 ; cind++ ) {
          ref_image[scale_ind][3*ind+cind]= inp_ref_image[3*ind+cind];
       }
    }
 }

 for ( scale_ind= 1 ; scale_ind< scale_nbr ; scale_ind++ ) {

    /*
    Compute sigma for that scale
    */

    sigma_n= pow(2,(double)(scale_ind-1));

    /*
    Allocate memory to hold the reference image at that scale
    */

    ref_image[scale_ind]= (int *)calloc(3*width*height,sizeof(int));

    /*
    Apply a Gaussian blur to the input reference image
    */

    precision= 5;
    gaussian_blur_rgb_image(
     inp_ref_image,
     width,
     height,
     sigma_n,
     precision,
     ref_image[scale_ind]
    );

    /*
    sprintf(filename,"duh%d.tiff",scale_ind);
    write_rgb_image_arr(
     filename,
     width,
     height,
     ref_image[scale_ind]
    );
    */
 }

 /*
 Compute the max distance in color space
 so that the distances can be normalized between 0 and 1
 considering all edges in the graph
 */

 random_walks_scale_space_get_max_dist(
  ref_image,
  v_arr,
  disp_alph_arr,
  width,
  height,
  scale_nbr,
  color_space,
  radius,
  con_level,
  con_level2,
  &max_dist
 );

 argc_duh= 1;
 argv_duh= (char **)calloc(1,sizeof(char *));
 argv_duh[0]= (char *)calloc(80,sizeof(char));
 argv_duh[0]="duh";
 ierr= lis_initialize(&argc_duh,&argv_duh);

 ierr= lis_solver_create(&solver);
 sprintf(option_string,"-i sor -p none -maxiter %d -print out",maxiter);
 ierr= lis_solver_set_option(option_string,solver);

 ierr= lis_matrix_create(0,&matA);
 ierr= lis_matrix_set_size(matA,0,width*height*scale_nbr);

 nnz_row= 1;
 if ( con_level == 1 ) {
    nnz_row+= 4;
 }
 else if ( con_level == 2 ) {
    nnz_row+= 8;
 }
 else {
    error_handler("random_walks_scale_space");
 }
 if ( con_level2 == 1 ) {
    nnz_row+= 2*1;
 }
 else if ( con_level2 == 2 ) {
    nnz_row+= 2*5;
 }
 else if ( con_level2 == 3 ) {
    nnz_row+= 2*9;
 }
 else {
    error_handler("random_walks_scale_space");
 }

 /*
 Although this is the exact number of non zero terms per row,
 lis still wants to call lis_realloc when setting values
 Apparently, it's a very good idea to over-estimate nnz
 */

 nnz_row*= 2;

 ierr= lis_matrix_malloc(matA,nnz_row,0);

 ierr= lis_vector_create(0,&vecx);
 ierr= lis_vector_set_size(vecx,0,width*height*scale_nbr);

 ierr= lis_vector_create(0,&vecb);
 ierr= lis_vector_set_size(vecb,0,width*height*scale_nbr);

 ierr= lis_vector_create(0,&vecy);
 ierr= lis_vector_set_size(vecy,0,width*height*scale_nbr);

 for ( scale_ind= 0 ; scale_ind< scale_nbr ; scale_ind++ ) {

    /*
    Fill matrix and right-handside vector for that scale
    */

    for ( i= 0 ; i< height ; i++ ) {
       for ( j= 0 ; j< width ; j++ ) {
          ind= i*width+j;
          v= v_arr[ind];
          disp_alph= disp_alph_arr[ind];
          edge_alph= edge_alph_arr[ind];
          r= ref_image[scale_ind][3*ind+0];
          g= ref_image[scale_ind][3*ind+1];
          b= ref_image[scale_ind][3*ind+2];
          if ( color_space == 1 ) {
             rgb2xyz(r,g,b,&x,&y,&z);
             xyz2Lab(x,y,z,&CIEL,&CIEa,&CIEb);
          }

          if ( disp_alph == 255 ) {

             /*
             This is a seed
             Solution is obviously already known
             */

             matA_ind= scale_ind*nbr_pixels+ind;
             lis_matrix_set_value(
              LIS_INS_VALUE,
              matA_ind,
              matA_ind,
              1,
              matA
             );
             vecb_ind= scale_ind*nbr_pixels+ind;
             lis_vector_set_value(
              LIS_INS_VALUE,
              vecb_ind,
              v,
              vecb
             );

             continue;
          }

          /*
          If here, this is not a seed
          */

          if ( disp_alph != 0 ) {
             error_handler("random_walks");
          }

          if ( edge_alph == 255 ) {

             /*
             Pixel is opaque in the segmented reference image
             Do as if solution was already known even though it's not
             Disparity will have to be computed later
             */

             matA_ind= scale_ind*nbr_pixels+ind;
             lis_matrix_set_value(
              LIS_INS_VALUE,
              matA_ind,
              matA_ind,
              1,
              matA
             );
             vecb_ind= scale_ind*nbr_pixels+ind;
             lis_vector_set_value(
              LIS_INS_VALUE,
              vecb_ind,
              v,
              vecb
             );

             continue;
          }

          /*
          If here, pixel is not opaque in the segmented reference image

          if ( edge_alph != 0 ) {
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
                   disp_alph2= disp_alph_arr[ind2];
                   edge_alph2= edge_alph_arr[ind2];
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

                   /*
                   Normalize the Euclidian distance
                   so that it's between 0 and 1
                   */

                   dist/= max_dist;
                   if ( dist > 1. )
                    error_handler("random_walks_scale_space");

                   if ( edge_alph2 == 255 ) {

                      /*
                      Neighboring pixel is opaque in the segmented reference image
                      */

                      /*
                      Do not consider this neighboring pixel in the equation that involves
                      the pixel under consideration and the neighboring pixels
                      */

                      continue;
                   }

                   /*
                   Instead of using an exponential like Grady,
                   we're gonna use a sigmoid
                   */
                   /*
                   G= exp(-(double)beta*dist);
                   */
                   G= 2/(1+exp((double)beta*dist));

                   sigma+= G;

                   if ( disp_alph2 == 0 ) {

                      /*
                      This is not a seed
                      */

                      matA_ind= scale_ind*nbr_pixels+ind;
                      matA_ind2= scale_ind2*nbr_pixels+ind2;
                      lis_matrix_set_value(
                       LIS_INS_VALUE,
                       matA_ind,
                       matA_ind2,
                       -G,
                       matA
                      );
                   }
                   else if ( disp_alph2 == 255 ) {

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
          } /* loop on previous, current, and next scale */

          if ( sigma == 0.0 ) {

             /*
             error_handler("random_walks_scale_space");
             */

             /*
             All neighboring pixels are opaque in the segmented reference image
             Do as if solution was already known even though it's not
             */

             matA_ind= scale_ind*nbr_pixels+ind;
             lis_matrix_set_value(
              LIS_INS_VALUE,
              matA_ind,
              matA_ind,
              1,
              matA
             );
             vecb_ind= scale_ind*nbr_pixels+ind;
             lis_vector_set_value(
              LIS_INS_VALUE,
              vecb_ind,
              v,
              vecb
             );

             continue; 
          }
         
          matA_ind= scale_ind*nbr_pixels+ind;
          lis_matrix_set_value(
           LIS_INS_VALUE,
           matA_ind,
           matA_ind,
           sigma,
           matA
          );
          vecb_ind= scale_ind*nbr_pixels+ind;
          lis_vector_set_value(
           LIS_INS_VALUE,
           vecb_ind,
           sigmav2,
           vecb
          );
       }
    }
 } /* loop on scale_ind */

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
       disp_alph= disp_alph_arr[ind];

       v= 1.0;
       for ( scale_ind= 0 ; scale_ind< scale_nbr ; scale_ind++ ) {
          vecx_ind= scale_ind*nbr_pixels+ind;
          ierr= lis_vector_get_value(vecx,vecx_ind,&value);
          if ( value < 1.0e-16 )
           value= 0;
          v*= value;
       }
       v= pow(v,1/(double)scale_nbr);

       if ( disp_alph != 0 ) { 

          /*
          This is a seed
          */

          continue;
       } 

       v_arr[ind]= v;
    }
 }

 ierr= lis_matrix_destroy(matA);
 ierr= lis_vector_destroy(vecx);
 ierr= lis_vector_destroy(vecb);
 ierr= lis_vector_destroy(vecy);
 ierr= lis_solver_destroy(solver);

 ierr= lis_finalize();

 /*
 Free memory to hold the reference images at each scale
 */

 for ( scale_ind= 0 ; scale_ind< scale_nbr ; scale_ind++ ) {
    free(ref_image[scale_ind]);
 }
 free(ref_image);

}
