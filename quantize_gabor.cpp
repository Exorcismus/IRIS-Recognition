#include "opencv/highgui.h"
#include "opencv/cv.h"
#include <math.h>
#include <stdio.h>
#include "variables.h"
#include "quantize_gabor.h"







void quantize_gabor_image_real(){

CvScalar  old_val;
CvScalar   new_val;

	for(int i=0; i<60;i++){
		for(int j=0;j<360;j++){
		
		old_val =	cvGet2D(real_gabor,i,j);
		
		if (old_val.val[0]>=0)
			new_val.val[0]=1;
		else
            new_val.val[0]=0;

		cvSet2D(quantized_real,i,j,new_val);
		
		}}



	
	cvShowImage("quantized_real",quantized_real);


}

void quantize_gabor_image_imag(){



CvScalar   old_val;
CvScalar   new_val;

	for(int i=0; i<60;i++){
		for(int j=0;j<360;j++){
		
		
	old_val =	cvGet2D(imag_gabor,i,j);
		
		if (old_val.val[0]>=0)
			new_val.val[0]=1;
		else
			new_val.val[0]=0;
		cvSet2D(quantized_imag,i,j,new_val);
		
		
		}}

	cvShowImage("quantized_imag",quantized_imag);


}