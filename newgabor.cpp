 
#include <cv.h>
#include <highgui.h>
#include <stdlib.h>
#include <stdio.h>


IplImage* src = 0;
IplImage* image = 0;
IplImage* dest = 0;
CvMat * kernel=0;

int kernel_size =9;

float var = 32;
float w = 10;
float phase = 45*CV_PI/180;

int pos_var = 0;
int pos_w = 0;
int pos_phase = 0;

void Process(int pos)   
{
	//CvFileStorage *fs = cvOpenFileStorage("kernel.xml",NULL,CV_STORAGE_WRITE);

	int x,y;
	float kernel_val;
	var = (float)pos_var;
	w = (float)pos_w;
	phase = (float) pos_phase*CV_PI/180;

	cvZero(kernel);
	for (x = -kernel_size/2+1;x<=kernel_size/2; x++) {
		for (y = (-kernel_size/2+1);y<=(kernel_size/2); y++) {

			kernel_val = 1/(2*CV_PI*var)*exp( -((x*x)+(y*y))/(2*var))*cos( w*x*cos(phase)+w*y*sin(phase));

			cvSet2D(kernel,y+kernel_size/2-1,x+kernel_size/2-1,cvScalar(kernel_val));
			}
		}
	//cvWrite( fs, "kernel", kernel, cvAttrList(0,0) );
	//cvReleaseFileStorage(&fs);

	cvFilter2D(src, dest,kernel);
    cvShowImage("Process window",dest);
}   





int main( int argc, char** argv )
{
    char* filename = argv[1];
    if( (image = cvLoadImage("5.bmp",1)) == 0 )
        return -1;
    
		
	kernel = cvCreateMat(kernel_size,kernel_size,CV_32FC1);

	src = cvCreateImage(cvSize(image->width,image->height),IPL_DEPTH_8U,1);

	cvCvtColor(image,src,CV_BGR2GRAY);
    dest = cvCloneImage(src);
	
    cvNamedWindow("Process window",1);

    cvShowImage("Process window",src);

    cvCreateTrackbar("Variance","Process window", &pos_var,50,Process);
    cvCreateTrackbar("Pulsation","Process window",&pos_w ,50,Process);
	cvCreateTrackbar("Phase","Process window",&pos_phase ,180,Process);

    cvWaitKey(0);
     
    cvReleaseImage(&src);
    cvReleaseImage(&image);
    cvReleaseImage(&dest);
    
    cvDestroyWindow("Process window"); 
    return 0;
}

  

