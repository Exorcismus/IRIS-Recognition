
#include "opencv/highgui.h"
#include "opencv/cv.h"
#include <math.h>
#include <stdio.h>
#include "normalize.h"
#include "variables.h"




void normalize(float &xp,float &yp, float &rp,float &xi, float &yi,float &ri ,IplImage *src){



        float xpupil;   // point on pupil x
        float ypupil;   // point on pupil y
        float xiris;    // point on iris x
        float yiris;    // point on iris y
        float xpa=0;
        float ypa=0;

        float theta = 0 ;  //col
        float K = 0.0;     
	    
        
      //  int N = cvRound(ri-rp);       // height
		int N = 60; 
        //int r = N-1;

		normalized    = cvCreateImage   (cvSize(360,N), IPL_DEPTH_8U, 3);
       greynormalized= cvCreateImage   (cvSize(360,N), IPL_DEPTH_8U, 1);
       eqnormalized  = cvCreateImage   (cvSize(360,N), IPL_DEPTH_8U, 1);

       
       
        // x = p1.x + blend * (p2.x - p1.x);
        // y = p1.y + blend * (p2.y - p1.y);
            
        
        
            for ( theta=0;theta<360;theta++) // col
        {
                
                  // point on the pupil
                    xpupil = (rp * cos(theta*(CV_PI/180)  )) + xp;  
                    ypupil = (rp * sin(theta*(CV_PI/180)  )) + yp;
             
                 // point on the iris
                    
                    
                    
                   xiris =  (ri *cos(theta*(CV_PI/180))   ) + xi;  
                   yiris =  (ri * sin(theta*(CV_PI/180))  ) + yi;
        
                  //  xiris =  ((N+rp) *cos(theta*(CV_PI/180))   ) + (N+rp);  
                   // yiris =  ((N+rp)* sin(theta*(CV_PI/180))  ) + (N+rp);
             
            

             // WHERE N = HEIGHT OF THE IMAGE
             
         for(K=0;K<N;K++)                  //row
			 // for(K=N-1;K>=0;K--)
				 { 
 
					 {
            // xpa= xpupil* (   1-  (((float)K+1)/(float)N)   )+ (xiris *(   ((float)K+1)/(float)N));
            // ypa= ypupil* (   1-  (((float)K+1)/(float)N)   )+ (yiris *(   ((float)K+1)/(float)N));
					 }
					 {
            // xpa = xpupil +    ((float)K/(float)N)*(     (xiris - xpupil)   );
           // ypa = (ypupil) +    ((float)K/(float)N)*(     (yiris  - ypupil)  );
					 }
					 {
            //xps=(    xp+ (rp*cos(theta*(3.14159265/180))  )           )
        //yps=(    yp+ (rp*sin(theta*(3.14159265/180))  )           )

        //xis=(xi+ (ri*cos(theta*(3.14159265/180))))
        //yis=(yi+ (ri*sin(theta*(3.14159265/180))))

        // r =    (((float)K+1)/(float)N)
            
        // xpa=(1-r)*xps + r * xis
        // ypa=(1-r)*yps + r * yis
         			 }

					 xpa=  ( (1-(((float)K+1)/(float)N) ) * (xp+ (rp*cos(theta*(CV_PI/180))))   ) + (  (((float)K+1)/(float)N) *  (xi+ (ri*cos(theta*(CV_PI/180))))  );
					 ypa=  ( (1-(((float)K+1)/(float)N) ) * (yp+ (rp*sin(theta*(CV_PI/180))))  ) + (  (((float)K+1)/(float)N) *  (yi+ (ri*sin(theta*(CV_PI/180))))   );


   

            // printf ( "xpa:%f , ypa:%f  \n",xpa,ypa);
            cvSet2D( normalized, K, theta ,cvGet2D(src,(ypa),(xpa)));    //Nomralized(theta,k)= image(xpa,ypa);
         

 }

             

}

 cvNamedWindow( "normalized");
cvShowImage( "normalized",normalized);

}






