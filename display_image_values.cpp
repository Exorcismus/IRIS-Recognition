#include "opencv/highgui.h"
#include "opencv/cv.h"
#include <math.h>
#include <stdio.h>
#include "display_image_values.h"




void display_image_values(IplImage *img,CvMat* kernel){

	if(kernel==0){

for(int i=0;i<=img->height/2;i++){
		for(int j=0;j<=img->width/2;j++){
		
			

			CvScalar value;
			value=cvGet2D(img,i,j);
			printf("   <%1.2f>   ",value.val[0]);
	
		}}
	}

	else if (img==0){
	
	
	
	for(int i=0;i<=kernel->height/2;i++){
		for(int j=0;j<=kernel->width/2;j++){
		
			

			CvScalar value;
			value=cvGet2D(kernel,i,j);
			printf("  <%1.2f>  ",value.val[0]);

		
		}}
	
	}

}
