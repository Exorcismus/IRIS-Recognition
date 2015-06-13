class gabor{


void gabor(){}
/*!
    \fn CvGabor::Init(int iMu, int iNu, double dSigma, double dF)
Initilize the.gabor

Parameters:
        iMu     The orientations which is iMu*PI.8
        iNu     The scale can be from -5 to infinit
        dSigma  The Sigma value of gabor, Normally set to 2*PI
        dF      The spatial frequence , normally is sqrt(2)

Returns:

Initilize the.gabor with the orientation iMu, the scale iNu, the sigma dSigma, the frequency dF, it will call the function creat_kernel(); So a gabor is created.
 */
void CvGabor::Init(int iMu, int iNu, double dSigma, double dF)
{
  //Initilise the parameters
    bInitialised = false;
    bKernel = false;

    Sigma = dSigma;
    F = dF;
   
    Kmax = PI/2;
   
    // Absolute value of K
    K = Kmax / pow(F, (double)iNu);
    Phi = PI*iMu/8;
    bInitialised = true;
    Width = mask_width();
    Real = cvCreateMat( Width, Width, CV_32FC1);
    Imag = cvCreateMat( Width, Width, CV_32FC1);
    creat_kernel();  
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*!
    \fn CvGabor::Init(double dPhi, int iNu, double dSigma, double dF)
Initilize the.gabor

Parameters:
        dPhi    The orientations
        iNu     The scale can be from -5 to infinit
        dSigma  The Sigma value of gabor, Normally set to 2*PI
        dF      The spatial frequence , normally is sqrt(2)

Returns:
        None

Initilize the.gabor with the orientation dPhi, the scale iNu, the sigma dSigma, the frequency dF, it will call the function creat_kernel(); So a gabor is created.filename  

    The name of the image file
        file_format     The format of the file, e.g. GAN_PNG_FORMAT
        image   The image structure to be written to the file
        octrlstr        Format-dependent control structure

 */
void CvGabor::Init(double dPhi, int iNu, double dSigma, double dF)
{

    bInitialised = false;
    bKernel = false;
    Sigma = dSigma;
    F = dF;
   
    Kmax = PI/2;
   
    // Absolute value of K
    K = Kmax / pow(F, (double)iNu);
    Phi = dPhi;
    bInitialised = true;
    Width = mask_width();
    Real = cvCreateMat( Width, Width, CV_32FC1);
    Imag = cvCreateMat( Width, Width, CV_32FC1);
    creat_kernel();  
}






/*!
    \fn CvGabor::mask_width()
Give out the width of the mask

Parameters:
        None

Returns:
        The long type show the width.

Return the width of mask (should be NxN) by the value of Sigma and iNu.
 */
long CvGabor::mask_width()
{

    long lWidth;
    if (IsInit() == false)  {
       perror ("Error: The Object has not been initilised in mask_width()!\n");
       return 0;
    }
    else {
       //determine the width of Mask

      double dModSigma = Sigma/K;
      double dWidth = round(dModSigma*6 + 1);

      //test whether dWidth is an odd.
      if (fmod(dWidth, 2.0)==0.0)
		  dWidth++;
      lWidth = (long)dWidth;

      return lWidth;
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*!
    \fn CvGabor::creat_kernel()
Create gabor kernel

Parameters:
        None

Returns:
        None

Create 2 gabor kernels - REAL and IMAG, with an orientation and a scale
 */
void CvGabor::creat_kernel()
{
   
    if (IsInit() == false) {perror("Error: The Object has not been initilised in creat_kernel()!\n");}
    else {

      CvMat *mReal, *mImag;

      mReal = cvCreateMat( Width, Width, CV_32FC1);
      mImag = cvCreateMat( Width, Width, CV_32FC1);
     
      /**************************** Gabor Function ****************************/
      int x, y;
      double dReal;
      double dImag;
      double dTemp1, dTemp2, dTemp3;

      for (int i = 0; i < Width; i++)
      {
          for (int j = 0; j < Width; j++)
          {
              x = i-(Width-1)/2;   
              y = j-(Width-1)/2;

              dTemp1 = (pow(K,2)/pow(Sigma,2))*exp(-(pow((double)x,2)+pow((double)y,2))*pow(K,2)/(2*pow(Sigma,2)));

              dTemp2 = cos(   K*cos(Phi)*x +  K*sin(Phi)*y  )- exp(-(pow(Sigma,2)/2));              //Re (s(x, y)) = cos (2PI(u0 x + v0 y) + P)
              dTemp3 = sin(   K*cos(Phi)*x +  K*sin(Phi)*y  );                                       //Im (s(x, y)) = sin (2PI(u0 x + v0 y) + P)
            
			  dReal = dTemp1*dTemp2;
              dImag = dTemp1*dTemp3;

              //gan_mat_set_el(pmReal, i, j, dReal);
              //cvmSet( (CvMat*)mReal, i, j, dReal );
                cvSetReal2D((CvMat*)mReal, i, j, dReal );

              //gan_mat_set_el(pmImag, i, j, dImag);
              //cvmSet( (CvMat*)mImag, i, j, dImag );
                cvSetReal2D((CvMat*)mImag, i, j, dImag );

          }
       }
       /**************************** Gabor Function ****************************/
       bKernel = true;
       cvCopy(mReal, Real, NULL);
       cvCopy(mImag, Imag, NULL);
      //printf("A %d x %d Gabor kernel with %f PI in arc is created.\n", Width, Width, Phi/PI);
       cvReleaseMat( &mReal );
       cvReleaseMat( &mImag );
     }
}
















/*!
    \fn CvGabor::conv_img_a(IplImage *src, IplImage *dst, int Type)
 */
///////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////// ////////IMPORTANT///////////////////
////////////////////////////////////////////////
///////////////////////////////////////////////
void CvGabor::conv_img_a(IplImage *src, IplImage *dst, int Type)
{
    double ve, re,im;
 
    int width = src->width;
    int height = src->height;
    CvMat *mat = cvCreateMat(src->width, src->height, CV_32FC1);
   
    for (int i = 0; i < width; i++)
    {
       for (int j = 0; j < height; j++)
       {
              ve = cvGetReal2D((IplImage*)src, j, i);
              cvSetReal2D( (CvMat*)mat, i, j, ve );
       }
    }

    CvMat *rmat = cvCreateMat(width, height, CV_32FC1);
    CvMat *imat = cvCreateMat(width, height, CV_32FC1);

    CvMat *kernel = cvCreateMat( Width, Width, CV_32FC1 );

    switch (Type)
    {
      case CV_GABOR_REAL:
        cvCopy( (CvMat*)Real, (CvMat*)kernel, NULL );
        cvFilter2D( (CvMat*)mat, (CvMat*)mat, (CvMat*)kernel, cvPoint( (Width-1)/2, (Width-1)/2));
        break;
      case CV_GABOR_IMAG:
        cvCopy( (CvMat*)Imag, (CvMat*)kernel, NULL );
        cvFilter2D( (CvMat*)mat, (CvMat*)mat, (CvMat*)kernel, cvPoint( (Width-1)/2, (Width-1)/2));
        break;





      case CV_GABOR_MAG:

        /* Real Response */
        cvCopy( (CvMat*)Real, (CvMat*)kernel, NULL );
        cvFilter2D( (CvMat*)mat, (CvMat*)rmat, (CvMat*)kernel, cvPoint( (Width-1)/2, (Width-1)/2));
       
		
		/* Imag Response */
        cvCopy( (CvMat*)Imag, (CvMat*)kernel, NULL );


        cvFilter2D( (CvMat*)mat, (CvMat*)imat, (CvMat*)kernel, cvPoint( (Width-1)/2, (Width-1)/2));
        /* Magnitude response is the square root of the sum of the square of real response and imaginary response */
        for (int i = 0; i < width; i++)
        {
           for (int j = 0; j < height; j++)
           {
               re = cvGetReal2D((CvMat*)rmat, i, j);
               im = cvGetReal2D((CvMat*)imat, i, j);
               ve = sqrt(re*re + im*im);
               cvSetReal2D( (CvMat*)mat, i, j, ve );
           }
        }      
        break;



      case CV_GABOR_PHASE:
        break;
    }
   
    if (dst->depth == IPL_DEPTH_8U)
    {
        cvNormalize((CvMat*)mat, (CvMat*)mat, 0, 255, CV_MINMAX, NULL);
        for (int i = 0; i < width; i++)
        {
            for (int j = 0; j < height; j++)
            {
                ve = cvGetReal2D((CvMat*)mat, i, j);
                ve = cvRound(ve);
                cvSetReal2D( (IplImage*)dst, j, i, ve );
            }
        }
     }

     if (dst->depth == IPL_DEPTH_32F)
     {
         for (int i = 0; i < width; i++)
         {
            for (int j = 0; j < height; j++)
            {
                ve = cvGetReal2D((CvMat*)mat, i, j);
                cvSetReal2D( (IplImage*)dst, j, i, ve );
            }
         }
     }

    cvReleaseMat(&kernel);
    cvReleaseMat(&imat);
    cvReleaseMat(&rmat);
    cvReleaseMat(&mat);
}




};