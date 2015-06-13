
#include "opencv/highgui.h"
#include "opencv/cv.h"
#include <math.h>
#include <stdio.h>
#include "variables.h"
#include "generate_code.h"




void generate_code(){


	CvScalar qreal;
	CvScalar qimag;

	int step = 1;


	for(int i=0;i<60;i++){
		for(int j=0;j<360;j=j++){
		
		
		
			qreal = cvGet2D(quantized_real,i,j);
			qimag = cvGet2D(quantized_imag,i,j);
		
		
		
		
			if (qreal.val[0]==0 && qimag.val[0]==0)
	  {
			
			cvSet2D(code,i,2*j,qreal);
			cvSet2D(code,i,2*j+step,qimag);

			
	   }
	   else if (qreal.val[0]==0 && qimag.val[0]==1)
	   {

            cvSet2D(code,i,2*j,qreal);
			cvSet2D(code,i,2*j+step,qimag);
			

	   }

       else if (qreal.val[0]==1 && qimag.val[0]==0)
	   {
	   
	        cvSet2D(code,i,2*j,qreal);
			cvSet2D(code,i,2*j+step,qimag);
			
	   
	   }
	   else if (qreal.val[0]==1 && qimag.val[0]==1)
	   {
	   
	        cvSet2D(code,i,2*j,qreal);
			cvSet2D(code,i,2*j+step,qimag);
			
	   
	   }
		
		





		
		
		}}

	










}