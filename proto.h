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
);

void gaussian_blur_image(
 int *I,
 int xdim,
 int ydim,
 double sigma,
 int precision,
 int *I_out
);

void gaussian_blur_rgb_image(
 int *inp_image_arr,
 int xdim,
 int ydim,
 double sigma,
 int precision,
 int *out_image_arr
);

void random_walks(
 int *ref_image,
 double *v_arr,
 int *alph_arr,
 int width,
 int height,
 int beta,
 int maxiter
);

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
);

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
);

void random_walks_sor(
 int *ref_image,
 double *v_arr,
 int *alph_arr,
 int width,
 int height,
 int beta
);
