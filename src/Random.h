#ifndef RANDOM_H
#define RANDOM_H

#include <math.h>
#include <time.h>
#include <algorithm>
#include <vector>

#include "DenseVector.h"

template <typename TPrecision>
class Random {

  public:
    Random(){
      srand( time(NULL) );
    };



    //Marsaglia-polar algorithm for sampling from Normal(0, 1)
    void Normal(TPrecision &s1, TPrecision &s2){
        double x1, x2, w;
         do {
           x1 = 2.0 * rand()/(double) RAND_MAX - 1.0;
           x2 = 2.0 * rand()/(double) RAND_MAX - 1.0;
           w = x1 * x1 + x2 * x2;
         } while ( w >= 1.0 );

         w = sqrt( (-2.0 * log( w ) ) / w );
         s1 =(TPrecision)( x1 * w );
         s2 =(TPrecision)( x2 * w );
    };
    

    //Marsaglia- Normal(0, sigma)
    void Normal(TPrecision &s1, TPrecision &s2, TPrecision sigma){
      Normal(s1, s2);
      s1*=sigma;
      s2*=sigma;
    };

    //Marsaglia
    TPrecision Normal(){
      TPrecision s1, s2;
      Normal(s1, s2);
      return s1;
    };



    TPrecision Uniform(){
      return rand()/(TPrecision) RAND_MAX; 
    };



    std::vector<unsigned int> Permutation(unsigned int N){
      std::vector<unsigned int> a;

      for( unsigned int i=0; i<N; i++){
        a.push_back(i); 
      }

      std::random_shuffle(a.begin(), a.end());

      return a;
    };

  
};

#endif
