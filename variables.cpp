#include "opencv/highgui.h"
#include "opencv/cv.h"
#include "variables.h"
#include <stdio.h>




 IplImage *normalized=0;
 IplImage *greynormalized=0;
 IplImage *eqnormalized=0;

 IplImage *real_gabor=cvCreateImage   (cvSize(360,60), IPL_DEPTH_32F, 1);
 IplImage *imag_gabor=cvCreateImage   (cvSize(360,60), IPL_DEPTH_32F, 1);
 IplImage *mag_gabor =cvCreateImage   (cvSize(360,60), IPL_DEPTH_32F, 1);
 //IplImage *phase_gabor=cvCreateImage   (cvSize(360,60), IPL_DEPTH_32F, 1);

IplImage *quantized_real = cvCreateImage   (cvSize(360,60), IPL_DEPTH_32F, 1);
IplImage *quantized_imag = cvCreateImage   (cvSize(360,60), IPL_DEPTH_32F, 1);
IplImage *code            =cvCreateImage   (cvSize(720,60), IPL_DEPTH_8U, 1);
IplImage * scaled_code    = cvCreateImage(cvSize(720,60), IPL_DEPTH_8U, 1);


 CvMat* kernel_real =0;
 CvMat* kernel_imag =0;
 CvMat* kernel_mag  =0;
 cv::Mat Lkernel;




IplImage * captured_iris =cvCreateImage(cvSize(640,480),IPL_DEPTH_8U, 3 );

IplImage *src    = cvLoadImage("d:/test.jpg");

IplImage *dst    = cvCreateImage( cvSize( src->width, src->height ), IPL_DEPTH_8U, 3 );


IplImage *grayscale    = cvCreateImage( cvSize( src->width, src->height ), IPL_DEPTH_8U, 1 );

IplImage *smooth = cvCreateImage( cvSize( src->width, src->height ), IPL_DEPTH_8U, 1 );
IplImage *pupil  = cvCreateImage( cvSize( src->width, src->height ), IPL_DEPTH_8U, 1 );
IplImage *iris   = cvCreateImage( cvSize( src->width, src->height ), IPL_DEPTH_8U, 1 );
IplImage *pedge  = cvCreateImage( cvSize( src->width, src->height ), IPL_DEPTH_8U, 1 );
IplImage *iedge  = cvCreateImage( cvSize( src->width, src->height ), IPL_DEPTH_8U, 1 );

//IplImage *contour = cvCreateImage( cvSize( src->width, src->height ), IPL_DEPTH_8U, 1 );
//IplImage *con     = cvCreateImage( cvSize( src->width, src->height ), IPL_DEPTH_8U, 3 );

   
IplImage *eyelid_mask  = cvCreateImage( cvGetSize(src), IPL_DEPTH_8U, 1 );
IplImage *iris_mask    = cvCreateImage(cvGetSize(src), 8, 1);
IplImage *mask         = cvCreateImage(cvGetSize(src), IPL_DEPTH_8U, 3 );
IplImage *res          = cvCreateImage(cvGetSize(src), 8, 3);

IplImage *compare_result=cvCreateImage   (cvSize(720,60), IPL_DEPTH_8U, 1);






//cvZero(mask);
//cvZero(eyelid_mask);
//cvZero(iris_mask);
//cvZero(res);


int xroi,yroi;

float xp=0;
float yp=0;
float rp=0;
float xi=0;
float yi=0;
float ri=0;

float shift = (float)(src->width/5);

