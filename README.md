# dmag4

Written and tested on Ubuntu 22.04.

dmag4 is a semi-automatic 2d to 3d converter. It creates a depth map given depth clues and a so-called edge map which segments the 2d image.

The executable needs lis.

To create the lis stuff, mkdir lis, tar xvf lis-1.4.11.tar, cd lis-1.4.11, ./configure --prefix=LIS_PATH (where LIS_PATH is where your lis directory is located), make, and make install.

To create the executable, compile the code in directory "dmag4" using "make -f Makefile_g/Makefile_O" and then go into the "main" directory and create the exec using "make".

Test cases are given in the "test" directory under "main".

Info about dmag4 (theory behind it and how to use it) can be found here:

[Depth Map Automatic Generator 4 (DMAG4) ](https://3dstereophoto.blogspot.com/2014/02/depth-map-automatic-generator-4-dmag4.html)

[Semi-Automatic 2D to 3D Image Conversion using Random Walks ](https://3dstereophoto.blogspot.com/2014/02/semi-automatic-2d-to-3d-conversion.html)

Dependencies (check the Makefile in main):

common repo
