//////////// ///////// ///////// ///////// ///////// ///////// /////////
//////
> > //////////// ////////Linkage des libraries/// ///////// ///////// /////////
/
> > //////////// ///////// ///////// ///////// ///////// ///////// /////////
/////
> >
> > #pragma comment(lib, "cv.lib")
> > #pragma comment(lib, "cvaux.lib")
> > #pragma comment(lib, "cvcam.lib")
> > #pragma comment(lib, "cxcore.lib" )
> > #pragma comment(lib, "cxts.lib")
> > #pragma comment(lib, "highgui.lib" )
> >
> > //////////// ///////// //////INCLUDE/ ///////// ///////// ///////// ////////
> > //////////// ///////// ///////// ///////// ///////// ///////// ///////// ///
> > #include <stdio.h>
> > #include <cv.h>
> > #include <highgui.h>
> > //#include "cvgabor.h"
> >
> >
> >
> > //////////// ///////// ///////// ///////// ///////// ///////// /////////
/////
> > //////////// ////Fonction qui charge et affiche une image/////// ////////
> > //////////// ///////// ///////// ///////// ///////// ///////// /////////
////
> >
> > IplImage* funcImage(const char* filename)
> > {
> > // Ouvre un fichier image
> > IplImage* image = cvLoadImage( filename, -1);
> >
> > // En cas de probl�me, fin de la fonction
> > if(!image)
> > {
> > printf("L'image est introuvable ou incorrect\n" );
> > return NULL;
> > }
> >
> > // Renvoi l'image
> > return image;
> > }
> >
> >
> >
> //////////// ///////// ///////// ///////// ///////// ///////// /////////
///////
> > //////////// //// Fonction qui transforme en gris
> ///"manuellement" ///////// /
> > //////////// ///////// ///////// ///////// ///////// ///////// /////////
//////
> >
> > IplImage* funcGray(IplImage* img_color)
> > {
> > // Nouvelle image avec 1 canal et 8 bits par canaux
> > IplImage* img_gray = cvCreateImage(
> > cvSize(img_color- >width, img_color->height) ,
> > 8, 1);
> >
> > // Pour chaque pixel
> > int i,color,size = img_gray->width * img_gray->height;
> > unsigned char* data = (unsigned char*)img_color- >imageData;
> > for(i=0; i<size; i++)
> > {
> > // Applique la formule : R*0.299+G*0. 587+B*0.114
> > color = data[i*3 + 0] * 299;
> > color += data[i*3 + 1] * 587;
> > color += data[i*3 + 2] * 114;
> > color /= 1000;
> >
> > // Place la couleur
> > img_gray->imageData [i] = color;
> >
> > }
> > // Applique un egaliseur d'histogramme
> > cvEqualizeHist( img_gray, img_gray);
> >
> > //Applique un filtre gaussian
> > cvSmooth( img_gray, img_gray,
> > CV_GAUSSIAN,
> > 1, 1, 0, 0 );
> >
> > // Renvoi l'image grise
> > return img_gray;
> > }
> >
> > //////////// ///////// ///////// ///////// ///////// ///////// /////////
////
> > //////////// ///Fonction qui applique un filtre de Gabor/////// ////////
> > //////////// ///////// ///////// ///////// ///////// ///////// /////////
////
> > /*
> > IplImage* funcGabor(IplImage* img_gray)
> > {
> >
> >
> > double Sigma = 2*PI;
> > double F = sqrt(2.0);
> > CvGabor* gabor1;
> > gabor1->Init( PI/4, 3, Sigma, F);
> >
> > IplImage *img_gabor =
> cvCreateImage( cvSize(img_ gray->width, img_gray- >height), IPL_DEPTH_8U, 1);
> >
> > gabor1->conv_ img(img_gray, img_gabor, 1);
> >
> > // Renvoi l'image gabor
> > return img_gray;
> > }
> >
> > */
> >
> //////////// ///////// ///////// ///////// ///////// ///////// /////////
/////////
> > //////////// // Fonction qui applique un seuil pour la binarisation
> //////
> >
> //////////// ///////// ///////// ///////// ///////// ///////// /////////
////////
> >
> > IplImage* funcThreshold( IplImage* img_gray, unsigned char thres)
> > {
> > // Cr�e une image avec 1 canal et 8 bits par canaux
> > IplImage* img_thres = cvCreateImage(
> > cvSize(img_gray- >width, img_gray->height) ,
> > 8, 1);
> >
> > // Pour chaque pixel
> > int i, size = img_thres->width * img_thres->height, color;
> > for(i=0; i<size; i++)
> > {
> >
> > color = ((unsigned char)img_gray- >imageData[ i]);
> >
> > // Selon le seuil
> > if(color < thres)
> > color = 0;
> > else
> > color = 255;
> >
> > // Place le pixel
> > img_thres->imageDat a[i] = color;
> > }
> > //Applique un seuil d'erosion
> >
> > cvErode(img_ thres,img_ thres,NULL, 1);
> > //cvDilate(img_ thres,img_ thres,NULL, 1);
> >
> > // Renvoi l'image grise
> > return img_thres;
> >
> > }
> >
> >
> //////////// ///////// ///////// ///////// ///////// ///////// /////////
///////// ////
> > //////////// ///////// ////fonction
> Erode/////// ///////// ///////// ///////// //////
> >
> //////////// ///////// ///////// ///////// ///////// ///////// /////////
///////// ////
> >
> >
> >
> >
> //////////// ///////// ///////// ///////// ///////// ///////// /////////
///////// //////
>
> > //////////// ///////// /////Fonction
> MAIN//////// ///////// ///////// ///////// ///////
> >
> //////////// ///////// ///////// ///////// ///////// ///////// /////////
///////// /////
> >
> > int main()//int argc, char* argv[])
> > {
> >
> > // R�cup�re l'image
> > IplImage* image = funcImage("fingerpr int_grey. jpg");
> > if(!image)
> > {
> > printf("Impossible de trouver l'image");
> > return 0;
> > }
> >
> >
> > // Transforme en gris
> > IplImage* gray = funcGray(image) ;
> >
> > //Applique le filtre de Gabor
> > // IplImage* gabor = funcGabor(gray) ;
> >
> > // Applique un seuil
> > IplImage* thres = funcThreshold( gray, 95);
> >
> >
> >
> > //Affichage de Gabor
> > /*cvNamedWindow( "Real Response", 1);
> > cvShowImage( "Real Response",gabor) ;
> > cvWaitKey(0) ;
> > cvDestroyWindow( "Real Response");
> > */
> >
> >
> >
> >
> > // Affiche les r�sultats
> > cvNamedWindow( "Fingerprint_ Image");
> > cvShowImage( "Fingerprint_ Image", image);
> > cvWaitKey(1000) ;
> >
> > cvNamedWindow( "Fingerprint_ Grey");
> > cvShowImage( "Fingerprint_ Grey", gray);
> > cvWaitKey(1000) ;
> >
> > cvNamedWindow( "Fingerprint_ Binarized" );
> > cvShowImage( "Fingerprint_ Binarized" , thres);
> > cvSaveImage( "Fingerprint_ Binarized. jpg", thres);
> > cvWaitKey(000) ;
> >
> >
> > // Supprime les images de la m�moire
> > cvReleaseImage( &image);
> > cvReleaseImage( &gray);
> > cvReleaseImage( &thres);
> >
> > return 0;
> > }
> >