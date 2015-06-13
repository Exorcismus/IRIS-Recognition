#include "opencv/highgui.h"
#include "opencv/cv.h"
#include "variables.h"




void save_image_code(char *file_name ,IplImage* code){



cvScale(code,scaled_code,255);
cvSaveImage(file_name,scaled_code);





}