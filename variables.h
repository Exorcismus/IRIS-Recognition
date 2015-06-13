#ifndef variables_H
#define variables_H





extern IplImage* captured_iris;
extern IplImage *normalized;
extern IplImage *greynormalized;
extern IplImage *eqnormalized;

extern IplImage *real_gabor;
extern IplImage *imag_gabor;
extern IplImage *mag_gabor ;
//extern IplImage *phase_gabor;
extern IplImage *quantized_real;
extern IplImage *quantized_imag;

extern IplImage *code      ;
extern IplImage * scaled_code;


extern CvMat* kernel_real ;
extern CvMat* kernel_imag ;
extern CvMat* kernel_mag  ;

extern cv::Mat Lkernel;


extern IplImage *src    ;
extern IplImage *dst   ;


extern IplImage *grayscale   ;

extern IplImage *smooth ;
extern IplImage *pupil  ;
extern IplImage *iris   ;
extern IplImage *pedge  ;
extern IplImage *iedge  ;

//IplImage *contour ;
//IplImage *con     ;

   
extern IplImage *eyelid_mask  ;
extern IplImage *iris_mask    ;
extern IplImage *mask         ;
extern IplImage *res         ;
extern IplImage *compare_result     ;
extern IplImage *normalized ;








extern int xroi,yroi;

extern float xp;
extern float yp;
extern float rp;
extern float xi;
extern float yi;
extern float ri;

extern float shift ;


#endif 