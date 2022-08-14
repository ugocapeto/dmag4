#include "header.h"
#include "proto.h"

void gaussian_blur_image(
 int *I,
 int xdim,
 int ydim,
 double sigma,
 int precision,
 int *I_out
)

{

 int size;
 double *G;
 int i;
 double x;
 double pi=acos(-1);
 double *row;
 int j;
 double norm;
 int k;
 double val;
 double *col;
 double *I_dbl;
 double val_dbl;
 int val_int;

 /*
 Allocate memory for output image in double form
 */

 I_dbl= (double *)calloc(xdim*ydim,sizeof(double));

 /*
 Copy input image into output image
 */

 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {
       I_dbl[i*xdim+j]= (double)I[i*xdim+j];
    }
 }

 /*
 We are gonna apply a 1d kernel
 (defined by sigma and precision)
 in both directions
 */

 size= (int)(precision*sigma)+1;

 G= (double *)calloc(size,sizeof(double));
 for ( i= 0 ; i< size ; i++ ) {
    x= (double)i;
    G[i]= 1/(sqrt(2*pi)*sigma)*exp(-x*x/(2*sigma*sigma));
 }

 /*
 We are gonna copy the image row by row,
 convolute each row, and
 put the results back in I
 */

 row= (double *)calloc(xdim,sizeof(double));
 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {
       norm= 0.0;
       k= 0;
       val= G[k]*I_dbl[i*xdim+j];
       norm+= G[k];
       for ( k= 1 ; k< size ; k++ ) {
          if ( !( (j+k)<xdim ) )
           continue;
          val+= G[k]*I_dbl[i*xdim+(j+k)];
          norm+= G[k];
       }
       for ( k= 1 ; k< size ; k++ ) {
          if ( !( (j-k)>=0 ) )
           continue;
          val+= G[k]*I_dbl[i*xdim+(j-k)];
          norm+= G[k];
       }
       row[j]= val/norm;
    }
    for ( j= 0 ; j< xdim ; j++ ) {
       I_dbl[i*xdim+j]= row[j];
    }
 }
 free(row);

 /*
 We are gonna copy the image col by col,
 convolute each col, and
 put the results back in I
 */

 col= (double *)calloc(ydim,sizeof(double));
 for ( j= 0 ; j< xdim ; j++ ) {
    for ( i= 0 ; i< ydim ; i++ ) {
       norm= 0.0;
       k= 0;
       val= G[k]*I_dbl[i*xdim+j];
       norm+= G[k];
       for ( k= 1 ; k< size ; k++ ) {
          if ( !( (i+k)<ydim ) )
           continue;
          val+= G[k]*I_dbl[(i+k)*xdim+j];
          norm+= G[k];
       }
       for ( k= 1 ; k< size ; k++ ) {
          if ( !( (i-k)>=0 ) )
           continue;
          val+= G[k]*I_dbl[(i-k)*xdim+j];
          norm+= G[k];
       }
       col[i]= val/norm;
    }
    for ( i= 0 ; i< ydim ; i++ ) {
       I_dbl[i*xdim+j]= col[i];
    }
 }
 free(col);

 for ( i= 0 ; i< ydim ; i++ ) {
    for ( j= 0 ; j< xdim ; j++ ) {
       val_dbl= I_dbl[i*xdim+j];
       val_int= (int)(val_dbl+0.5);
       if ( val_int < 0 )
        val_int= 0;
       if ( val_int > 255 )
        val_int= 255;
       I_out[i*xdim+j]= val_int;
    }
 }

 free(G);

 /*
 Free memory for output image in double form
 */

 free(I_dbl);

}
