#include "opencv/highgui.h"
#include "opencv/cv.h"
#include <math.h>
#include <stdio.h>

#include "variables.h"
#include "display_image_values.h"
#include "process_image.h"
#include "run_camera.h"

#include "normalize.h"
#include "create_gabor.h"
#include "quantize_gabor.h"
#include "generate_code.h"
#include "compare_code.h"
#include "save_image_code.h"


void Smooth(IplImage *src , IplImage *dst);




void main(){



Gabor gabor;

run_camera();


cvCvtColor( src, grayscale, CV_RGB2GRAY );
cvCopy(src,dst);
Smooth(src,smooth);
//cvShowImage( "smooth",smooth);


process_image();

//normalize(float xp,float xi, float yp,float yi, float rp,float ri ,IplImage *src)
normalize(xp, yp, rp, xi, yi,ri ,src);

cvCvtColor(normalized,greynormalized,CV_RGB2GRAY);
//cvShowImage("greynormalized",greynormalized);

cvEqualizeHist(greynormalized,eqnormalized);
//cvShowImage("eqnormalized",eqnormalized);



//1-void create_gabor_kernel_v1( float sig, float thet, float lm, float gamma , float ps , float bw)
//gabor.create_gabor_kernel_v1(0,45,8,1,0,1);   //sigma = 1.6 , gamma = 1 for kernel of size 9



//2-void create_gabor_kernel_v2(int kernel_size, float sig, float thet, float lm, float gamma , float ps);
//gabor.create_gabor_kernel_v2(21,5,45,10,1,90); 



//3-void create_gabor_kernel_v3(int kernel_size, float sig, float thet, float freq, float gamma , float ps)
//gabor.create_gabor_kernel_v3(9,2,0,0.2,1,0); 



//4-void create_gabor_kernel_v4(int kernel_size, double sigma, double theta, double lambd, double gamma, double psi, int ktype);
//gabor.create_gabor_kernel_v4(9,2,0,35,1,0,CV_32F); 



//5-void create_gabor_kernel_v5 (int kernel_size, double sigma,double theta ,  double freq ,double gamma)
//gabor.create_gabor_kernel_v5(9,2,0,0.25,1);


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////


cvFilter2D(eqnormalized,real_gabor, kernel_real);
cvShowImage("real wavelet",real_gabor);


cvFilter2D(eqnormalized,imag_gabor, kernel_imag);
cvShowImage("imag wavelet",imag_gabor);

//cvFilter2D(eqnormalized,mag_gabor, kernel_mag);
//cvShowImage("mag wavelet",mag_gabor);





//quantize_gabor_image_real();
//quantize_gabor_image_imag();



generate_code();

//save_image_code("code_4.bmp",code);       //scales code before saving

//display_image_values(0,kernel_real);     // FOR DEBUGING PURPOSES
//display_image_values(code,0);            // FOR DEBUGING PURPOSES



float HD =compare_code(code,"code_4.bmp");
printf("Hamming Distance = %1.2f",HD);


cvScale(code,scaled_code,255);
cvShowImage("code",scaled_code);

cvNamedWindow( "Source");
cvShowImage( "Source",dst);

//printf ( "xp:%f ,xi:%f, yp:%f, yi:%f, rp:%f, ri:%f",xp,xi,yp,yi,rp,ri);



cvWaitKey();
cvReleaseImage(&dst);
cvReleaseImage(&smooth);
cvReleaseImage(&pupil);
cvReleaseImage(&iris);
cvReleaseImage(&pedge);
cvReleaseImage(&iedge);
cvReleaseImage(&src);

cvReleaseImage(&eqnormalized);
cvReleaseImage(&greynormalized);
cvReleaseImage(&normalized);
cvReleaseImage(&real_gabor);
cvReleaseImage(&imag_gabor);
cvReleaseImage(&mag_gabor);
//cvReleaseImage(&phase_gabor);
cvReleaseImage(&quantized_real);
cvReleaseImage(&quantized_imag);
cvReleaseImage(&code);




}








void Smooth(IplImage *src , IplImage *dst){
cvCvtColor( src, dst, CV_RGB2GRAY );
cvSmooth( dst, dst, CV_MEDIAN,9,9);

}





