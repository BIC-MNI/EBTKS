#include <iostream>
#include <math.h>

#include <EBTKS/Matrix.h>
#include <complex.h>

  dcomplex fft_linear[64]={2016+0*I,
-32+651.37496399959*I,
-32+324.901452403484*I,
-32+215.72647697328*I,
-32+160.874863748027*I,
-32+127.751161080643*I,
-32+105.489862686026*I,
-32+89.4340087196953*I,
-32+77.254833995939*I,
-32+67.6583154415565*I,
-32+59.8677891772605*I,
-32+53.3887745786723*I,
-32+47.8913844052956*I,
-32+43.147005231575*I,
-32+38.9921128188152*I,
-32+35.3065592234713*I,
-32+32*I,
-32+29.0031094086127*I,
-32+26.2617213065171*I,
-32+23.7328174807052*I,
-32+21.3817164134176*I,
-32+19.1800618778216*I,
-32+17.1043563504253*I,
-32+15.1348728285222*I,
-32+13.254833995939*I,
-32+11.4497830820648*I,
-32+9.7070938754349*I,
-32+8.0155827261217*I,
-32+6.365195756149*I,
-32+4.7467516012271*I,
-32+3.1517249074292*I,
-32+1.5720591926229*I,
-32+0*I,
-32-1.5720591926229*I,
-32-3.1517249074293*I,
-32-4.7467516012271*I,
-32-6.3651957561491*I,
-32-8.0155827261218*I,
-32-9.7070938754349*I,
-32-11.4497830820648*I,
-32-13.254833995939*I,
-32-15.1348728285222*I,
-32-17.1043563504253*I,
-32-19.1800618778216*I,
-32-21.3817164134176*I,
-32-23.7328174807052*I,
-32-26.2617213065171*I,
-32-29.0031094086127*I,
-32-32*I,
-32-35.3065592234712*I,
-32-38.9921128188152*I,
-32-43.147005231575*I,
-32-47.8913844052956*I,
-32-53.3887745786722*I,
-32-59.8677891772605*I,
-32-67.6583154415565*I,
-32-77.254833995939*I,
-32-89.4340087196952*I,
-32-105.489862686026*I,
-32-127.751161080643*I,
-32-160.874863748027*I,
-32-215.72647697328*I,
-32-324.901452403484*I,
-32-651.37496399959*I
};


#ifdef USE_COMPMAT
int test_fft(void)
{
  Mat<double>   real_matrix(64,1,0.0);
  Mat<dcomplex> complex_matrix(64,1,0.0);
  double sigma=8;
  
  //make a distribution
  for(int i=0;i<64;i++)
  {
    real_matrix(i,0) = i;
  }
  
  //perform FFT
  complex_matrix=fft(real_matrix,0,1);
  
  double diff=0.0;
  for(int i=0;i<64;i++)
  {
    diff+=abs(fft_linear[i]-complex_matrix(i,0))*abs(fft_linear[i]-complex_matrix(i,0));
  }
  diff/=64;
  if(diff>0.001) //threshold
  {
    std::cerr<<"Discrepancy in forward FFT:"<<diff<<std::endl;
    std::cout<<"Result\tExpected\t"<<std::endl;
    for(int i=0;i<64;i++)
    {
      std::cout<<complex_matrix(i,0)<<"\t"<<fft_linear[i]<<std::endl;
    }
    
    return 1;
  }
  
  Mat<dcomplex> real_matrix_2(64,1,0.0);
  real_matrix_2=ifft(complex_matrix,0,1);
  diff=0.0;
  for(int i=0;i<64;i++)
  {
    diff+=abs(real_matrix_2(i,0) - real_matrix(i,0))*abs(real_matrix_2(i,0)-real_matrix(i,0));
  }
  diff/=64;
  if(diff>0.001) //threshold
  {
    std::cerr<<"Discrepancy in inverse FFT:"<<diff<<std::endl;
    return 2;
  }
 
  return 0;
}
#endif //USE_COMPMAT

#ifdef USE_FCOMPMAT

int test_ffft(void)
{
  Mat<float>    real_matrix(64,1,0.0);
  Mat<fcomplex> complex_matrix(64,1,0.0);
  float sigma=8;
  
  //make a distribution
  for(int i=0;i<64;i++)
  {
    real_matrix(i,0) = i;
  }
  //perform FFT
  complex_matrix=ffft(real_matrix,0,1);
  
  double diff=0.0;
  for(int i=0;i<64;i++)
  {
    dcomplex d=(fcomplex)(fft_linear[i])-complex_matrix(i,0);
    diff+=abs(d)*abs(d);
  }
  diff/=64;
  if(diff>0.001) //threshold
  {
    std::cerr<<"Discrepancy in forward FFFT:"<<diff<<std::endl;
    std::cout<<"Result\tExpected\t"<<std::endl;
    for(int i=0;i<64;i++)
    {
      std::cout<<complex_matrix(i,0)<<"\t"<<fft_linear[i]<<std::endl;
    }
    return 1;
  }
  
  Mat<fcomplex> real_matrix_2(64,1,0.0);
  real_matrix_2=iffft(complex_matrix,0,1);
  diff=0.0;
  for(int i=0;i<64;i++)
  {
    diff+=abs(real_matrix_2(i,0)-real_matrix(i,0))*abs(real_matrix_2(i,0)-real_matrix(i,0));
  }
  diff/=64;
  if(diff>0.001) //threshold
  {
    std::cerr<<"Discrepancy in inverse FFFT:"<<diff<<std::endl;
    return 2;
  }
 
  return 0;
}
#endif 

int  main(int argc, char *argv[])
{
  int r=0;
  #ifdef USE_COMPMAT
  test_fft();
  if(r>0) return r;
  #endif 
  
  #ifdef USE_FCOMPMAT
  r=test_ffft();
  if(r>0) return r;
  #endif 
  
  return 0;
}
