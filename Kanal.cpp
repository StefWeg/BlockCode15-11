#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#define PI 3.141592654
 
float gauss(float mean, float sigma);

//************************************************************************

// Function kanal changes binary values into bipolar ones (-1/+1) and adds noise
// *wej - Input vector of binary values (0/1)
// *wyj - Output vector of real numbers
// es_n0 - Es/N0
// dl_kan - the number of input bits
void kanal(float es_n0, long dl_kan, int *wej, float *wyj)
{
  float mean=0;
  float es=1;
  float sygnal;
  float sigma;
  float s_n;
  long y;

  s_n=(float) pow(10, (es_n0/10));
  sigma=(float) sqrt (es/(2*s_n));

  for (y=0; y<dl_kan; y++)
    {
      sygnal = 2 * *(wej+y)-1; // change the binary value (0/1) into symbol (-1/+1) 
      *(wyj+y)=sygnal+gauss(mean,sigma);  // noise addition
    }
}
//**********************************************************************
float gauss(float mean, float sigma)
{
  double x;
  double z;
  
  z=(double)rand()/RAND_MAX;
  if (z==1.0) z=0.9999999;
  x=sigma*sqrt(2.0*log(1.0/(1.0-z)));

  z=(double)rand()/RAND_MAX;
  if (z==1.0) z=0.9999999;
  return((float)(mean+x*cos(2*PI*z)));
}
//*******************************************************************

