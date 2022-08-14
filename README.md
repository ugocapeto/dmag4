# sfm10

Written and tested on Ubuntu 22.04.

sfm10 (structure from motion 10) builds a sparse 3d scene reconstruction given a set of still frames and a focal length. Of course, it also computes the camera poses (rotation and translation) for each camera/view.

The Structure from Motion pipeline in sfm10 is as follows:
- For each camera/view, extract the features using SIFT.
- For each camera pair, compute the matches.
- For each camera pair, compute the good matches by removing outliers and keeping only the inliers. This is done using A Contrario RANSAC aka AC-RANSAC.
- Initialize the 3d scene with an adequate initial pair. The initial pair should be chosen so that the corresponding two views are not related by a homography, in other words, the baseline should be large enough (while having a good number of matches).
- Remove the spurious (low-confidence) 3d points. If the separation angle between the rays emanating from a given 3d point to the camera centers is below some threshold, the 3d point is considered unreliable and is therefore removed from the 3d reconstruction.
- Perform Bundle Adjustment on camera poses and the 3d points. The non-linear cost function is minimized with LBFGS instead of the more widely used Levenberg-Marquardt (just to be different).
- For each remaining camera:
-- Compute the camera pose using EPnP and add measurements to existing 3d points. This is called "resectioning" in Computer Vision lingo.
-- Add the 3d points seen by the camera.
-- Remove the spurious (low confidence) 3d points.
-- Perform Bundle Adjustment.

The executable needs lbfgs and pba.

To create the lbfgs stuff, mkdir liblbfgs, tar xvf liblbfgs-1.10.tar, cd liblbfgs-1.10, ./configure --prefix=LIBLBFGS_PATH (where LIBLBFGS_PATH is where your liblbfgs is located), make, and make install.

To create the pba stuff, unzip pba-master.zip, rename pba-master to pba, cd pba, and make -f makefile_no_gpu.

To create the executable, compile the code in directory "sfm10" using "make -f Makefile_g/Makefile_O" and then go into the "main" directory and create the exec using "make".

Test cases are given in the "test" directory under "main".

Info about sfm10 (theory behind it and how to use it) can be found here:

[Structure from Motion 10 (SfM10)](http://3dstereophoto.blogspot.com/2016/04/structure-from-motion-10-sfm10.html)

Dependencies (check the Makefile in main):

common repo

er9b repo
