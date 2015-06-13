
#include "opencv/highgui.h"
#include "opencv/cv.h"
#include <math.h>
#include <stdio.h>
#include "variables.h"
#include "create_gabor.h"










	

	
//////////////////////////////////////////
///               VERSION 1            ///
///         WEB MATLAB VERSION         ///
///           USING BANDWIDTH          ///
//////////////////////////////////////////
void Gabor:: create_gabor_kernel_v1( float sig, float thet, float lm, float gamma , float ps , float bw){

   

    float    sigma=sig ;
    float    lmbda=lm ;
    float    slratio;
	float    N;
	int      kernel_size;
   

	slratio = (1/CV_PI)*sqrt(log(2.0)/2)*( (pow(2,bw)+1) /(pow(2,bw)-1));

	if  (sigma == 0)
  sigma  = slratio * lmbda;
   else if (lmbda == 0)
  lmbda = sigma / slratio;


    float    theta = (thet*CV_PI/180);
    float    psi   = (ps  *CV_PI/180);    

      int      x = 0;
      int      y = 0;
      float    X = 0;
      float    Y = 0;
      

	  if (gamma <= 1 && gamma > 0)

    N = cvRound(2.5*sigma/gamma);  //sigma/gamma = 1.6 so N=4
	  else
    N = cvRound(2.5*sigma);
   
		

kernel_size = (int)(2*N+1);


kernel_real =cvCreateMat(kernel_size,kernel_size, CV_32FC1);
 kernel_imag =cvCreateMat(kernel_size,kernel_size, CV_32FC1);
 kernel_mag  =cvCreateMat(kernel_size,kernel_size, CV_32FC1);

     


      float    X2,Y2;
      float    sigma2;
      float    gamma2;
      float    F;

      float    real = 0;
      float    imag = 0;
      float    mag  = 0;
      

      float    K;
      float    temp1;
      float    temp2;
      float    temp3;

      

  int   kernel_mid    = (int) (kernel_size-1)/2;     //kernel anchor


     sigma2=  sigma*sigma;
     gamma2=  gamma*gamma;

	 K    = (    1/(2*(CV_PI)*sigma2) );
     F    = (2*CV_PI )/(float) lmbda ;          // Frequency


	printf("Sigma = %f ; Lmbda = %f , kernel_size = %d ",sigma,lmbda,kernel_size);


 for(int i=0;i<kernel_size;i++){   //row

      for(int j=0;j<kernel_size;j++){  //col


          x = j-(kernel_mid);        //col
          y = i-(kernel_mid);        //row


       
         //temp1 = (pow(K,2.0)/pow(Sigma,2.0))*exp(-(pow((double)x,2.0)+pow((double)y,2.0))*pow(K,2.0)/(2.0*pow(Sigma,2.0)));

         // temp2 = cos(   K*cos(Phi*(CV_PI/180))*x +  K*sin(Phi*(CV_PI/180))*y  );              //Re (s(x, y)) = cos (2PI(u0 x + v0 y) + P)
         // temp3 = sin(   K*cos(Phi*(CV_PI/180))*x +  K*sin(Phi*(CV_PI/180))*y  );              //Im (s(x, y)) = sin (2PI(u0 x + v0 y) + P)
 


          X =    (x*cos(theta))   + ( y*sin(theta))  ;
          Y =    (-x*sin(theta))  + ( y*cos(theta))  ;
 

     Y2    =  Y*Y;
     X2    =  X*X;
     

    temp1=K*exp( -0.5* (  (X2+(Y2*gamma2)) / (sigma2) )  );  //checked

	temp2= cos(  ((X*F)+ psi)  );     //checked

    temp3= sin(  ((X*F)+ psi)  );




          real = (temp1*temp2);                         //checked
          imag = (temp1*temp3);

          mag  = ((real*real)+(imag*imag));


       //  printf ( "x: %d , y: %d  ,real:%2f\n",x,y,real);

              //gan_mat_set_el(pmReal, i, j, dReal);
              //cvmSet( (CvMat*)mReal, i, j, dReal );

                cvSetReal2D((CvMat*)kernel_real, i, j, real );   // NOTE : REMEMBER TO CHANGE
                cvSetReal2D((CvMat*)kernel_imag, i, j, imag );
                cvSetReal2D((CvMat*)kernel_mag , i, j,mag  );


      }}

  

 //cvNormalize(kernel_real,kernel_real,1,0,CV_C,NULL);



//cvScale(kernel_real,kernel_real,0.5);
//cvAddS(kernel_real,cvScalar(0.5),kernel_real,NULL);

//cvShowImage("Kernel_real",kernel_real);
//cvShowImage("Kernel_imag",kernel_imag);
//cvShowImage("Kernenl_magnitude",kernel_mag);



 //Lkernel =cvCreateMat(kernel_size*10, kernel_size*10, CV_32FC1);
//cv::resize(kernel_real, Lkernel, Lkernel.size());
//Lkernel /= 2.;
//Lkernel += 0.5;
 // cv::imshow("kernel",Lkernel);


}





















///////////////////////////////////////////
//               VERSION 2              ///
//-          DIFFERENT SIGMAs           ///
///////////////////////////////////////////
void Gabor::create_gabor_kernel_v2(int kernel_size, float sig, float thet, float lm, float gamma , float ps){

    kernel_real =cvCreateMat(kernel_size,kernel_size, CV_32FC1);
    kernel_imag =cvCreateMat(kernel_size,kernel_size, CV_32FC1);
    kernel_mag  =cvCreateMat(kernel_size,kernel_size, CV_32FC1);

  


  
    float    theta = (thet*CV_PI/180);
    float    psi   = (ps  *CV_PI/180);    

      int      x = 0;
      int      y = 0;
      float    X = 0;
      float    Y = 0;
      
      
   

      float    F;
      float    lmbda = lm ; 
	  float    sigma =sig;
  
      float    X2,Y2;
   
      float    gamma2;
	  float    sigma_x2;
	  float    sigma_y2;
    

      float    real = 0;
      float    imag = 0;
      float    mag  = 0;
	  float    phase= 0;
      

      float    K;
      float    temp1;
      float    temp2;
      float    temp3;

      float    sigma_x = sigma;
	  float    sigma_y = sigma/gamma;

  int   kernel_mid    = (kernel_size-1)/2;     //kernel anchor


    
     gamma2=  gamma*gamma;
	
	 K    = (  1/(2*(CV_PI)*sigma_x*sigma_y) );
     F    = (2*CV_PI )/lmbda ;          // Frequency


	 printf("Sigma = %f ; Lmbda= %f",sigma,lmbda);

 for(int i=0;i<kernel_size;i++){

      for(int j=0;j<kernel_size;j++){



		  y = i-(kernel_mid);        //row
          x = j-(kernel_mid);        //col
                


       

          X =    (x*cos(theta))   + ( y*sin(theta))  ;
          Y =    (-x*sin(theta))  + ( y*cos(theta))  ;
 

     Y2        =  Y*Y;
     X2        =  X*X;
     sigma_x2  =  pow(sigma_x,2);
     sigma_y2  =  pow(sigma_y,2);
     

    temp1=K*exp( -0.5 *(     (X2/sigma_x2)+(Y2/sigma_y2)      ) );  //checked

	temp2= cos(  ((X*F)+ psi)    );     //checked

    temp3= sin(  ((X*F)+ psi)  );




          real = (temp1*temp2);                         //checked
          imag = (temp1*temp3);

          mag  = (real*real)+(imag*imag);
		  phase = atan(imag/real);


        //printf ( "x: %d , y: %d  ,real:%2f\n",x,y,real);
		  

              

                cvSetReal2D((CvMat*)kernel_real, i, j, real );   // NOTE : REMEMBER TO CHANGE
                cvSetReal2D((CvMat*)kernel_imag, i, j, imag );
                cvSetReal2D((CvMat*)kernel_mag , i, j,mag  );
				//cvSetReal2D(phase_gabor, i, j,phase );




      }}



cvShowImage("Real",kernel_real);
cvShowImage("imag",kernel_imag);
cvShowImage("magnitude",kernel_mag);
//cvShowImage("Phase asd",phase_gabor);



Lkernel =cvCreateMat(kernel_size*10, kernel_size*10, CV_32FC1);
cv::resize(kernel_real, Lkernel, Lkernel.size());
  Lkernel /= 2.;
  Lkernel += 0.5;
  cv::imshow("kernel",Lkernel);


}









///////////////////////////////////////////
//               VERSION 3              ///
//-           USING FREQUENCY           ///
///////////////////////////////////////////
void Gabor::create_gabor_kernel_v3(int kernel_size, float sig, float thet, float freq, float gamma , float ps){

    kernel_real =cvCreateMat(kernel_size,kernel_size, CV_32FC1);
    kernel_imag =cvCreateMat(kernel_size,kernel_size, CV_32FC1);
    kernel_mag  =cvCreateMat(kernel_size,kernel_size, CV_32FC1);



  
    float    theta = (thet*CV_PI/180);
    float    psi   = (ps  *CV_PI/180);    

      int      x = 0;
      int      y = 0;
      float    X = 0;
      float    Y = 0;
      
      
   

      
      float    lmbda; 
	  float    F=freq;
	  float    sigma =sig;
  
      float    X2,Y2;
   
     
	 
    

      float    real = 0;
      float    imag = 0;
      float    mag  = 0;
      

      float    K;
      float    temp1;
      float    temp2;
      float    temp3;

     

  int   kernel_mid    = (kernel_size-1)/2;     //kernel anchor


    
   float  gamma2=  gamma*gamma;
   float  sigma2=  sigma*sigma;
	
	 K     = (  1/(sqrt(2*CV_PI*sigma)) );  // (  1/(2*sqrt(gamma*CV_PI)*sigma) )
	 lmbda = (2*CV_PI)/F ;          


	 printf("Sigma = %f ; Lmbda= %f",sigma,lmbda);

 for(int i=0;i<kernel_size;i++){

      for(int j=0;j<kernel_size;j++){


          x = j-(kernel_mid);        //col
          y = i-(kernel_mid);        //row


       
         


          X =    (x*cos(theta))   + ( y*sin(theta))  ;
          Y =    (-x*sin(theta))  + ( y*cos(theta))  ;
 

     Y2        =  pow(Y,2);
     X2        =  pow(X,2);
   
     

    temp1= K* exp( -0.5 *(     (X2/sigma2)+(Y2/sigma2)      ) );  //checked

	temp2= cos(  ((X*F)+ psi)   );     //checked

    temp3= sin(  ((X*F)+ psi)   );




          real = (temp1*temp2);                         //checked
          imag = (temp1*temp3);

          mag  = (real*real)+(imag*imag);


        //printf ( "x: %d , y: %d  ,real:%2f\n",x,y,real);
		  

             

                 cvSetReal2D((CvMat*)kernel_real, i, j, real );   // NOTE : REMEMBER TO CHANGE
                cvSetReal2D((CvMat*)kernel_imag, i, j, imag );
                cvSetReal2D((CvMat*)kernel_mag , i, j,mag  );


			


      }}



cvShowImage("Real",kernel_real);
cvShowImage("imag",kernel_imag);
cvShowImage("magnitude",kernel_mag);



Lkernel =cvCreateMat(kernel_size*10, kernel_size*10, CV_32FC1);
cv::resize(kernel_real, Lkernel, Lkernel.size());
  Lkernel /= 2.;
  Lkernel += 0.5;
  cv::imshow("kernel",Lkernel);



}





////////////////////////////////////////////
///               VERSION 4              ///
///          OPENCV ONLINE VERSION       ///
///            DIFFERENT SIGMAS          ///
////////////////////////////////////////////
void Gabor::create_gabor_kernel_v4 (int kernel_size, double sigma, double theta, double lambd, double gamma, double psi, int ktype){


    double sigma_x = sigma;
    double sigma_y = sigma/gamma;

    int nstds = 3;

    int xmin, xmax, ymin, ymax;

	double c = cos(theta*(CV_PI/180));
	double s = sin(theta*(CV_PI/180));

    if( kernel_size > 0 )

        xmax = kernel_size-1/2;

    else

        xmax = std::max(fabs(nstds*sigma_x*c), fabs(nstds*sigma_y*s));

    if( kernel_size > 0 )

        ymax = kernel_size-1/2;

    else

        ymax = std::max(fabs(nstds*sigma_x*s), fabs(nstds*sigma_y*c));

        xmin = -xmax;
        ymin = -ymax;

    CV_Assert( ktype == CV_32FC1 || ktype == CV_64FC1 );

	//kernel_real =cvCreateMat(kernel_size,kernel_size, CV_32F);

	kernel_real = cvCreateMat(ymax - ymin + 1, xmax - xmin + 1, ktype);
	kernel_imag = cvCreateMat(ymax - ymin + 1, xmax - xmin + 1, ktype);
	
	
    double scale = 1/(2*CV_PI*sigma_x*sigma_y);

    double ex = -0.5/(sigma_x*sigma_x);

    double ey = -0.5/(sigma_y*sigma_y);

    double cscale = CV_PI*2/lambd;

    for( int y = ymin; y <= ymax; y++ )  //row
	{
    for( int x = xmin; x <= xmax; x++ )  //col
 {



            double xr = x*c + y*s;

            double yr = -x*s + y*c;

			double real = scale*exp(ex*xr*xr + ey*yr*yr)*cos((cscale*xr + psi));
			double imag = scale*exp(ex*xr*xr + ey*yr*yr)*sin((cscale*xr + psi));

			if( ktype == CV_32F ){

          cvSetReal2D((CvMat*)kernel_real,ymax - y,xmax - x,  (float)real );
		  cvSetReal2D((CvMat*)kernel_imag,ymax - y,xmax - x, (float)imag);
            //kernel_real.at<float>(ymax - y, xmax - x) = (float)v;
			}
			else{
           cvSetReal2D((CvMat*)kernel_real,ymax - y,xmax - x,real );
		   cvSetReal2D((CvMat*)kernel_imag,ymax - y,xmax - x, imag);
            // kernel_real.at<double>(ymax - y, xmax - x) = v;
			}
     	   

 }

	}


	cvShowImage("Real",kernel_real);
	cvShowImage("Imag",kernel_imag);

 Lkernel =cvCreateMat(kernel_size*10, kernel_size*10, CV_32FC1);
cv::resize(kernel_real, Lkernel, Lkernel.size());
  Lkernel /= 2.;
  Lkernel += 0.5;
  cv::imshow("kernel",Lkernel);



}






/////////////////////////////////////////////
///               VERSION 5               ///
///          ONLINE BOOK VERSION          ///
///            DIFFERENT SIGMAS           ///
///                                       ///
/////////////////////////////////////////////
void Gabor::create_gabor_kernel_v5 (int kernel_size, double sigma, double theta,  double freq, double gamma ){

    kernel_real =cvCreateMat(kernel_size,kernel_size, CV_32FC1);
    kernel_imag =cvCreateMat(kernel_size,kernel_size, CV_32FC1);
    kernel_mag  =cvCreateMat(kernel_size,kernel_size, CV_32FC1);




      int      x = 0;
      int      y = 0;
      float    X = 0;
      float    Y = 0;
      
      
   

      float    F = freq;
      
	 
  
      float    X2,Y2;
   
      float    gamma2;
	  float    sigma_x2;
	  float    sigma_y2;
    

      float    real = 0;
      float    imag = 0;
      float    mag  = 0;
      

      float    K;
      float    temp1;
      float    temp2;
      float    temp3;

      float    sigma_x = sigma;
	  float    sigma_y = sigma/gamma;

  int   kernel_mid    = (kernel_size-1)/2;     //kernel anchor


    
   
	
	 K    = (  1/(2*(CV_PI)*sigma_x*sigma_y) );
     


	

 for(int i=0;i<kernel_size;i++){

      for(int j=0;j<kernel_size;j++){



		  y = i-(kernel_mid);        //row
          x = j-(kernel_mid);        //col
                


       

          X =    (x*cos(theta))   + ( y*sin(theta))  ;
          Y =    (-x*sin(theta))  + ( y*cos(theta))  ;
 

     Y2        =  Y*Y;
     X2        =  X*X;
     sigma_x2  =  pow(sigma_x,2);
     sigma_y2  =  pow(sigma_y,2);
     

    temp1=K*exp( -0.5 *(     (X2/sigma_x2)+(Y2/sigma_y2)      ) );  

	temp2= cos( 2*CV_PI*F*sqrt(X2+Y2) );   

    temp3= sin( 2*CV_PI*F*sqrt(X2+Y2) );




          real = (temp1*temp2);                         //checked
          imag = (temp1*temp3);

          mag  = (real*real)+(imag*imag);


        //printf ( "x: %d , y: %d  ,real:%2f\n",x,y,real);
		  

              

               cvSetReal2D((CvMat*)kernel_real, i, j, real );   // NOTE : REMEMBER TO CHANGE
                cvSetReal2D((CvMat*)kernel_imag, i, j, imag );
                cvSetReal2D((CvMat*)kernel_mag , i, j,mag  );



      }}



cvShowImage("Real",kernel_real);
cvShowImage("imag",kernel_imag);
cvShowImage("magnitude",kernel_mag);



Lkernel =cvCreateMat(kernel_size*10, kernel_size*10, CV_32FC1);
cv::resize(kernel_real, Lkernel, Lkernel.size());
  Lkernel /= 2.;
  Lkernel += 0.5;
  cv::imshow("kernel",Lkernel);





}
