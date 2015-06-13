#include "opencv/highgui.h"
#include "opencv/cv.h"
#include <math.h>
#include <stdio.h>
#include "run_camera.h"
#include "variables.h"




void run_camera(){


CvMemStorage* eyestorage = cvCreateMemStorage(0);
CvHaarClassifierCascade* cascade = (CvHaarClassifierCascade*)cvLoad("D:/PROGRAMS\OPENCV 2.4/opencv/sources/data/haarcascades/haarcascade_mcs_righteye.xml",0,0,0 );
CvRect *eye;

cvNamedWindow( "Live", CV_WINDOW_AUTOSIZE );
CvCapture* capture = cvCreateCameraCapture( 0);


IplImage* frame;
while(1) {
frame = cvQueryFrame( capture );
if( !frame ) break;


cvClearMemStorage( eyestorage );
CvSeq* eyes = cvHaarDetectObjects(frame,cascade,eyestorage,1.1,1,0 ,cvSize(30, 30));

/* draw a rectangle for each detected eye */
    for( int i = 0; i < (eyes ? eyes->total : 0); i++ ) 
	{
        /* get one eye */
       eye = (CvRect*)cvGetSeqElem(eyes, i);
       
        /* draw a red rectangle */
		cvRectangle(frame,cvPoint((eye->x)-eye->width, (eye->y)-eye->height),cvPoint(eye->x + eye->width*2, eye->y + eye->height*2),CV_RGB(255, 0, 0),1,8,0 );
    }


cvShowImage( "Live", frame );

char c=cvWaitKey(10);

if(c==99){

CvRect frame_roi = cvRect(   eye->x-eye->width, eye->y-eye->height,  eye->width*3,   eye->height*3   ) ;
cvSetImageROI(frame, frame_roi  );
IplImage* captured = frame;

cvResize(captured,captured_iris);
cvShowImage("Captured IRIS",captured_iris);
break;
}


char e = cvWaitKey(33);
if( e == 27 )
break;




}



cvReleaseCapture( &capture );
cvDestroyWindow( "Live" );
}