#ifndef create_gabor_H
#define create_gabor_H


class Gabor
{
public:


	

void create_gabor_kernel_v1( float sig, float thet, float lm, float gamma , float ps , float bw);
void create_gabor_kernel_v2(int kernel_size, float sig, float thet, float lm, float gamma , float ps);
void create_gabor_kernel_v3(int kernel_size, float sig, float thet, float freq, float gamma , float ps);
void create_gabor_kernel_v4 (int kernel_size, double sigma, double theta, double lambd, double gamma, double psi, int ktype);
void create_gabor_kernel_v5 (int kernel_size, double sigma, double theta,  double freq, double gamma );
};
#endif 