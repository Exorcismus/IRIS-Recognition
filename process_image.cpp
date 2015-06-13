
#include "opencv/highgui.h"
#include "opencv/cv.h"
#include <math.h>
#include <stdio.h>
#include "variables.h"
#include "draw_histogram.h"
#include "process_image.h"

void process_image(){







///////////////////////////////////////////////////////
//////////////////// PUPIL/////////////////////////////
///////////////////////////////////////////////////////

int numBins = 256;
float range[] = {0, 255};
float *ranges[] = { range };
 
 CvHistogram *hist = cvCreateHist(1, &numBins, CV_HIST_ARRAY, ranges, 1);
 cvClearHist(hist);


	cvCalcHist(&smooth, hist, 0, 0);
    IplImage* imgHist = DrawHistogram(hist,1,1);
    cvClearHist(hist);
	

//cvShowImage("hist", imgHist);



cvThreshold(smooth,pupil,50,255,CV_THRESH_BINARY);
//cvShowImage( "pupi_binary",pupil);

cvCanny(pupil,pedge,40,50);
//cvShowImage( "pupil_edge",pedge);


//////////////////////////////////////////////////////////
//////////////////////IRIS////////////////////////////////
//////////////////////////////////////////////////////////

//cvEqualizeHist(smooth,smooth);
//cvShowImage("Equalized",smooth);

cvThreshold(smooth,iris,100,255,CV_THRESH_BINARY); //115
//cvShowImage( "iris_binary",iris);


//cvSobel(iris,iedge,1,0,3);
cvCanny(iris,iedge,1,255);
//cvShowImage( "iris_edge",iedge);


/////////////////////////////////////////////////////////
///////////////////////Eyelids///////////////////////////
/////////////////////////////////////////////////////////



cvThreshold(smooth,eyelid_mask,150,255,CV_THRESH_OTSU);
cvNot(eyelid_mask,eyelid_mask);
//cvShowImage("eyelid",eyelid_mask);



//cvAdaptiveThreshold(smooth,contour,255,CV_ADAPTIVE_THRESH_MEAN_C,CV_THRESH_BINARY,9,1);

//cvThreshold(smooth,contour,130,255,CV_THRESH_BINARY);
//cvShowImage( "contour",contour);


//CvSeq* firstContour = NULL;
//CvMemStorage* cstorage = cvCreateMemStorage(0);
//cvFindContours(con, cstorage, &firstContour,sizeof(CvContour), CV_RETR_LIST,CV_CHAIN_APPROX_SIMPLE);
//cvDrawContours(dst,firstContour,CV_RGB(0,255,0),CV_RGB(0,0,255),10,2,8);



CvMemStorage* storage_pupil = cvCreateMemStorage(0);

CvSeq* presults = cvHoughCircles(pedge,storage_pupil,CV_HOUGH_GRADIENT,2,src->width,255,1);
for( int i = 0; i < presults->total; i++ )
{
float* p = (float*) cvGetSeqElem( presults, i );
CvPoint pt = cvPoint( cvRound( p[0] ),cvRound( p[1] ) );

xp=cvRound( p[0] );
yp=cvRound( p[1] );
rp=p[2];

cvCircle(dst,pt,cvRound( p[2] ),CV_RGB(0,255,255),1,400);


xroi= xp-shift;
yroi= yp-shift;

cvRectangle(dst,cvPoint(( p[0] )-shift,p[1]-shift),cvPoint(( p[0] )+shift,p[1]+shift),CV_RGB(255,0,255), 1);

CvRect roi= cvRect(xroi  ,yroi,shift*2,shift*2);



cvSetImageROI( iedge, roi );


//cvShowImage("ROI",iedge);


}
////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

CvMemStorage* storage_iris = cvCreateMemStorage(0);

CvSeq* iresults = cvHoughCircles(iedge,storage_iris,CV_HOUGH_GRADIENT,2,src->width,1,50,50);
for( int i = 0; i < iresults->total; i++ )
{
float* p = (float*) cvGetSeqElem( iresults, i );

CvPoint pt = cvPoint( cvRound( p[0] )+xroi,cvRound( p[1] )+yroi );
cvCircle(dst,pt,cvRound( p[2] ),CV_RGB(255,0,0),1,400);


xi=cvRound( p[0] )+xroi;
yi=cvRound( p[1] )+yroi;
ri=(p[2]);


cvCircle(iris_mask,pt,cvRound( p[2] ),CV_RGB(255, 255, 255),-1, 8, 0);
//cvShowImage("iris_mask",iris_mask);


}
///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
///////////////////////////////////////////


cvResetImageROI(iedge);

cvAnd(dst,dst,res,iris_mask);
//cvShowImage("iris_mask",res);

cvAnd(res,res, mask, eyelid_mask);
//cvShowImage("Mask",mask);



//cvLogPolar(mask,finalres,cvPoint2D32f (xp,yp),100, CV_INTER_LINEAR  );
//cvShowImage("Final Result",finalres);


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/*




*/



}