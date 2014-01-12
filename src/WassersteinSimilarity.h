#ifndef WASSERSTEINSIMILARITY_H
#define WASSERSTEINSIMILARITY_H

#include "DenseVector.h"
#include "DenseMatrix.h"

#include "EuclideanMetric.h"

#include "SVD.h"
#include "RandomSVD.h"

#include <algorithm>
#include <map>

#include "GWT.h"
#include "Wasserstein.h"

template <typename TPrecision>
class WassersteinSimilarity : public SimilarityVisitor<TPrecision>{
private:

  std::map<int, GWTInfo<TPrecision> * > nodes;
  int maxID;
  EuclideanMetric<TPrecision> l2;


public:
  
   WassersteinSimilarity<TPrecision>(){
     maxID = 0;
   };

   void visit(GMRANode<TPrecision> *node){
    GWTInfo<TPrecision> *gwt = (GWTInfo<TPrecision>*)node->getNodeInfo();
    
    int ID = gwt->ID;
    if(ID > maxID){
      maxID = ID;
    }
    nodes[ID] = gwt;  
  };

  std::map<int, GWTInfo<TPrecision> * > &getNodes(){
    return nodes;
  };

  FortranLinalg::DenseMatrix<TPrecision> distances(){
    using namespace FortranLinalg;

    DenseMatrix<TPrecision> D(maxID+1, maxID+1);
    for(int i=0; i< D.M(); i++){
       //D(i, i) = 0;
       GWTInfo<TPrecision> *v1 = nodes[i];
       for(int j=0; j<D.N(); j++){
         GWTInfo<TPrecision> *v2 = nodes[j];
         TPrecision tmp = distance(v1, v2);
         if(tmp < 0){
           tmp = 0;
         }
         D(i, j) = tmp;
         D(j, i) = D(i, j);
      }
    }
    return D;
  };


   TPrecision distance(GWTInfo<TPrecision> *n1, GWTInfo<TPrecision> *n2){
    using namespace FortranLinalg;
    /*DenseVector<TPrecision> m1 = computeMean(n1);
    DenseVector<TPrecision> m2 = computeMean(n2);

    DenseVector<bool> flip( std::max(m1.N(), m2.N()) );
    Linalg<bool>::Set(flip, false);
    
    DenseMatrix<TPrecision> cov1 = computeCov(n1, flip);

    for(int i=0; i< std::min(m1.N(), m2.N()); i++){
      if( fabs( m1(i) - m2(i) ) > fabs( m1(i) + m2(i) ) ){
        flip(i) = true;
        m1(i) = - m1(i);
      }
    }
    DenseMatrix<TPrecision> cov2 = computeCov(n2, flip);


    TPrecision dMean;
    if(m1.N() < m2.N() ){
      Linalg<TPrecision>::Subtract(m2, m1, m2);
      dMean = Linalg<TPrecision>::Length(m2);
    }
    else{
      Linalg<TPrecision>::Subtract(m1, m2, m1);
      dMean = Linalg<TPrecision>::Length(m1);
    }
 
    TPrecision dCov = Wasserstein<TPrecision>::distance2(cov1, cov2);
    
    m1.deallocate();
    m2.deallocate();
    cov1.deallocate();
    cov2.deallocate();
*/
    TPrecision dCov = Wasserstein<TPrecision>::distance2(n1->phi, n1->sigma,
        n2->phi, n2->sigma);
    return dCov;
  };




private:


  static FortranLinalg::DenseVector<TPrecision> computeMean(GWTInfo<TPrecision> *gwt){
    using namespace FortranLinalg;
    if(gwt->parent != NULL){
      DenseVector<TPrecision> mean =
        Linalg<TPrecision>::Multiply(gwt->parent->phi, gwt->center, true);
      for(int i=0; i< mean.N(); i++){
        mean(i) /= gwt->parent->sigma(i);
      }
      return mean;
    }
    else{
      DenseVector<TPrecision> mean(gwt->phi.N());
      Linalg<TPrecision>::Zero(mean);
      return mean;
    }
  };



  static FortranLinalg::DenseMatrix<TPrecision> computeCov(GWTInfo<TPrecision> *gwt, FortranLinalg::DenseVector<bool> flip){
    using namespace FortranLinalg;
   

    if(gwt->parent != NULL){
      DenseMatrix<TPrecision> cov = Linalg<TPrecision>::Multiply(gwt->parent->phi, gwt->phi, true);
      for(int i=0; i<gwt->parent->sigma.N(); i++){
        TPrecision s1 = gwt->parent->sigma(i);
        
        if(flip(i)){
          s1 *= -s1;
        }
        else{
          s1 *= s1;
        }
        Linalg<TPrecision>::ScaleRow( cov, i, 1.0 / s1 );
      }
      for(int i=0; i<gwt->sigma.N(); i++){
        TPrecision s2 = gwt->sigma(i);
        Linalg<TPrecision>::ScaleColumn( cov, i, (s2*s2) );
      }
      return cov;
    }
    else{
      DenseMatrix<TPrecision> cov( gwt->sigma.N(), gwt->sigma.N() );
      Linalg<TPrecision>::Zero(cov);
      for(int i=0; i<cov.N(); i++){
        cov(i, i) = gwt->sigma(i) * gwt->sigma(i);
      }
      return cov;
   }
  
    
  };


};




#endif
