
#include "opencv/highgui.h"
#include "opencv/cv.h"
#include <math.h>
#include <stdio.h>
#include "variables.h"
#include "display_image_values.h"



float compare_code(IplImage * temp1,IplImage * temp2){




	CvScalar val;
	int   number_of_ones = 0;
	float HD;
    int   N = 43200;



	cvXor(temp1,temp2,compare_result);


	for(int i=0; i<compare_result->height;i++){
		for(int j=0; j<compare_result->width;j++){
		
		
			val=	cvGet2D(compare_result,i,j);
		if(val.val[0]==1)
			number_of_ones+=1;
		
			
		
		}}





HD = number_of_ones/N;









	return HD;


}









////////////////////////////////////////////////////


////////////////////////////////////////////////////
float compare_code(IplImage * temp1,  char * file_name){



	IplImage *gray_temp2 = cvCreateImage   (cvSize(720,60), IPL_DEPTH_8U, 1);
	IplImage *temp2      = cvCreateImage   (cvSize(720,60), IPL_DEPTH_8U, 1);

	IplImage *temp_scaled_code = cvLoadImage(file_name);

	cvCvtColor(temp_scaled_code,gray_temp2,CV_RGB2GRAY);  //loaded image is 3 channels , convert to 1


	cvScale(gray_temp2,temp2,(float)1/255);
	


	CvScalar val;
	float   number_of_ones = 0.0;
	float HD;
    float   N = 43200.0;


	cvXor(temp1,temp2,compare_result);


	for(int i=0; i<compare_result->height;i++){
		for(int j=0; j<compare_result->width;j++){
		
		
			val=	cvGet2D(compare_result,i,j);
		if(val.val[0]==1)
			number_of_ones+=1;
		
		
		}}

	//display_image_values(temp2,0);             // FOR DEBUGING PURPOSES

HD = number_of_ones/N;

	return HD;


}