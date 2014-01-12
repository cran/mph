#ifndef RELATIVESUBSPACESIMILARITY_H
#define RELATIVESUBSPACESIMILARITY_H

#include "DenseVector.h"
#include "DenseMatrix.h"

#include "EuclideanMetric.h"

#include "SVD.h"
#include "RandomSVD.h"

#include <algorithm>
#include <map>

#include "GWT.h"


template <typename TPrecision>
class RelativeSubspaceSimilarity : public SimilarityVisitor<TPrecision>{
private:

  std::map<int, DenseVector<TPrecision> > s_psi;
  int maxID;
  EuclideanMetric<TPrecision> l2;

public:
  
   RelativeSubspaceSimilarity<TPrecision>(){
     maxID = 0;
   };

   void visit(GMRANode *node){
    using namespace FortranLinalg;
    GWTNode<TPrecision> *gwt = (GWTNode<TPrecision>*)node->getGWTInfo();
    
    int ID = gwt->ID;
    if(ID > maxID){
      maxID = ID;
    }
    
    DenseVector<TPrecision> sigma;
    if(gwt->parent != NULL){
      
      DenseMatrix<TPrecision> phi_p = Linalg<TPrecision>::Multiply(gwt->parent->phi,
		    gwt->phi, true);
      DenseMatrix<TPrecision> psi = Linalg<TPrecision>::Multiply(gwt->parent->phi, phi_p);
      Linalg<TPrecision>::Subtract(gwt->phi, psi, psi);
      SVD<TPrecision> svd(psi);
      sigma = Linalg<TPrecision>::Copy(svd.S);
      svd.deallocate();

      phi_p.deallocate();
      psi.deallocate();
    }
    else{
      sigma = DenseVector<TPrecision>( gwt->phi.N() );
      Linalg<TPrecision>::Zero( sigma );
    }

    s_psi[ID] = sigma;
  };


  std::map<int, FortranLinalg::DenseVector<TPrecision> > &getSigmaPsi(){
    return s_psi;
  };

  FortranLinalg::DenseMatrix<TPrecision> distances(){
    using namespace FortranLinalg;
    DenseMatrix<TPrecision> D(maxID+1, maxID+1);
    for(int i=0; i< D.M(); i++){
       DenseVector<TPrecision> v1 = s_psi[i];
       for(int j=0; j<D.N(); j++){
         DenseVector<TPrecision> v2 = s_psi[j]; 
         D(i, j) = l2.distance(v1, v2);
         D(j, i) = D(i, j);
      }
    }
    return D;
  };


};




#endif
