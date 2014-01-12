#ifndef GRASSMANNIANNODESIMILARITY_H
#define GRASSMANNIANNODESIMILARITY_H

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
class GrassmannianNodeSimilarity : public SimilarityVisitor<TPrecision>{
private:

  std::map<int, GWTInfo<TPrecision> * > nodes;
  int maxID;
  EuclideanMetric<TPrecision> l2;


public:
  
   GrassmannianNodeSimilarity<TPrecision>(){
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
      D(i, i) = 0;
      
      GWTInfo<TPrecision> *v1 = nodes[i];
      DenseMatrix<TPrecision> A1 = Linalg<TPrecision>::Multiply(v1->phi, v1->phi,
          false, true);
      
      for(int j=i+1; j<D.N(); j++){
        GWTInfo<TPrecision> *v2 = nodes[j];
        DenseMatrix<TPrecision> A2 = Linalg<TPrecision>::Multiply(v2->phi, v2->phi,
            false, true);
        Linalg<TPrecision>::Subtract(A1, A2, A2);
        TPrecision d = Linalg<TPrecision>::SquaredFrobeniusNorm(A2);
        A2.deallocate();
        if(d < 0){
          d = 0;
        }
        D(i, j) = d;
        D(j, i) = D(i, j);
      }
      A1.deallocate();

    }
    return D;
  };







};




#endif
