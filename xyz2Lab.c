#include "header.h"
#include "proto.h"

void xyz2Lab(
 double x,
 double y,
 double z,
 double *L,
 double *a,
 double *b
)

{

 /*
 double xo=244.66128;
 double yo=255.0;
 double zo=277.63227;

 (*L)= 116*F(y/yo)-16;
 (*a)= 500*(F(x/xo)-F(y/yo));
 (*b)= 200*(F(y/yo)-F(z/zo));
 */

 /*
 double xn;
 double yn;
 double zn;

 rgb2xyz(255,255,255,&xn,&yn,&zn);
 */

 double xn=242.36628;
 double yn=255;
 double zn=277.63227;

 (*L)= 116*F(y/yn)-16;
 (*a)= 500*(F(x/xn)-F(y/yn));
 (*b)= 200*(F(y/yn)-F(z/zn));

}
