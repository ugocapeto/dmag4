#include "header.h"
#include "proto.h"

int main(
 int argc,
 char* argv[]
)

{

 FILE *fp;
 char filename[80];
 char filename_I[80];
 char filename_edge[80];
 char filename_disp[80];
 char filename_disp2[80];
 int width;
 int height;
 int *I;
 int *disp_arr;
 int *disp_alph_arr;
 double *v_arr;
 int i;
 int j;
 int ind;
 int beta;
 double dval;
 int maxiter;
 int scale_nbr;
 int con_level;
 int con_level2;
 int err_flag;
 int r;
 int g;
 int b;
 int *disp_rgb_arr;
 int a;
 int *edge_rgb_arr;
 int *edge_alph_arr;

 /*
 Let's read in the input file
 */

 fp= fopen("dmag4_input.txt","r");

 /*
 Get filename for reference image
 */

 fscanf(fp,"%s",filename);

 fprintf(stdout,"reference image = %s\n",filename);

 strcpy(filename_I,filename);

 /*
 Get filename for input scribbled disparity map
 */

 fscanf(fp,"%s",filename);

 fprintf(stdout,"scribbled disparity map = %s\n",filename);

 strcpy(filename_disp,filename);

 /*
 Get filename for edge image
 */

 fscanf(fp,"%s",filename);

 fprintf(stdout,"edge image = %s\n",filename);

 strcpy(filename_edge,filename);

 /*
 Get filename for output disparity map
 */

 fscanf(fp,"%s",filename);

 fprintf(stdout,"output disparity map = %s\n",filename);

 strcpy(filename_disp2,filename);

 /*
 Get beta
 */

 fscanf(fp,"%d",&beta);

 fprintf(stdout,"beta = %d\n",beta);

 /*
 Get maximum number of iterations
 */

 fscanf(fp,"%d",&maxiter);

 fprintf(stdout,"maxiter = %d\n",maxiter);

 /*
 Get number of scales
 */

 fscanf(fp,"%d",&scale_nbr);

 fprintf(stdout,"scale_nbr = %d\n",scale_nbr);

 /*
 Get connection level between graph nodes within a scale
 */

 fscanf(fp,"%d",&con_level);

 fprintf(stdout,"con_level = %d\n",con_level);

 /*
 Get connection level between graph nodes across scales
 */

 fscanf(fp,"%d",&con_level2);

 fprintf(stdout,"con_level2 = %d\n",con_level2);

 /*
 Done reading the input file
 */

 fclose(fp);

 /*
 Check con_level
 */

 err_flag= 1;

 if ( con_level == 1 ) {
    err_flag= 0;
 }
 else if ( con_level == 2 ) {
    err_flag= 0;
 }
 else {
    err_flag= 1;
 }

 if ( err_flag == 1 ) {
    fprintf(stdout,"con_level should be either 1 or 2.\n");
    return 1;
 }

 /*
 Check con_level2
 */

 err_flag= 1;

 if ( con_level2 == 1 ) {
    err_flag= 0;
 }
 else if ( con_level2 == 2 ) {
    err_flag= 0;
 }
 else {
    err_flag= 1;
 }

 if ( err_flag == 1 ) {
    fprintf(stdout,"con_level2 should be either 1 or 2.\n");
    return 1;
 }

 /*
 Load reference image
 */

 err_flag= load_rgb_image(
  filename_I,
  &I,
  &width,
  &height
 );

 if ( err_flag == 1 )
  return 1;

 /*
 Load the scribbled disparities
 (includes the alpha channel)
 alpha=  0 if transparent (unknown disparity)
 alpha=255 if 100% opaque (known disparity=seed)
 */

 err_flag= load_rgba_image(
  filename_disp,
  &disp_rgb_arr,
  &disp_alph_arr,
  &width,
  &height
 );

 if ( err_flag == 1 )
  return 1;

 /*
 check for the presence of semi-transparent pixels
 */
 
 err_flag= 0;
 for ( i= 0 ; i< height ; i++ ) {
    for ( j= 0 ; j< width ; j++ ) {
       ind= i*width+j;
       a= disp_alph_arr[ind];
       if ( a > 0 && a < 255 )
        err_flag= 1;
    }
 }
 if ( err_flag == 1 ) {
    fprintf(stdout,"Sparse depth map contains semi-transparent pixels. Use the threshold alpha command in Gimp or Photoshop to make them either transparent or opaque.\n");
    return 1;
 }

 /*
 Check for the presence of an edge image
 */

 fp= fopen(filename_edge,"r");

 if ( fp == 0 ) {

   /*
   There is no edge image provided
   */

   /*
   Make all pixels of the edge image transparent
   */ 

   edge_rgb_arr= (int *)calloc(3*width*height,sizeof(int));
   edge_alph_arr= (int *)calloc(width*height,sizeof(int));
 }
 else {

    fclose(fp);

    /*
    There is an edge image provided
    */

    if ( con_level == 2 ) {
       fprintf(stdout,"Can't have con_level equal to 2 if an edge image is provided.\n");
       return 1;
    }

    /*
    Load edge image
    */

    err_flag= load_rgba_image(
     filename_edge,
     &edge_rgb_arr,
     &edge_alph_arr,
     &width,
     &height
    );

    if ( err_flag == 1 )
     return 1;
 }

 /*
 check for the presence of semi-transparent pixels
 */

 err_flag= 0;
 for ( i= 0 ; i< height ; i++ ) {
    for ( j= 0 ; j< width ; j++ ) {
       ind= i*width+j;
       a= edge_alph_arr[ind];
       if ( a > 0 && a < 255 )
        err_flag= 1;
    }
 }
 if ( err_flag == 1 ) {
    fprintf(stdout,"Edge image contains semi-transparent pixels. Use the threshold alpha command in Gimp or Photoshop to make them either transparent or opaque.\n");
    return 1;
 }

 fprintf(stdout,"Computing the disparity map ...\n");

 compute_disparity_map(
  I,
  disp_rgb_arr,
  disp_alph_arr,
  edge_alph_arr,
  width,
  height,
  beta,
  maxiter,
  scale_nbr,
  con_level,
  con_level2,
  &disp_arr
 );

 fprintf(stdout,"Computing the disparity map ... done.\n");

 /*
 Let's dump the output disparity map into an image
 */

 err_flag= write_image(
  filename_disp2,
  disp_arr,
  width,
  height
 );

 if ( err_flag == 1 ) {
    return 1;
 }

 /*
 Free memory
 */

 free(I);

 free(disp_rgb_arr);
 free(disp_alph_arr);

 free(edge_rgb_arr);
 free(edge_alph_arr);

 free(disp_arr);

}
