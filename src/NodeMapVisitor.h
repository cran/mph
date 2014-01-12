#ifndef NODEMAPVISITOR_H
#define NODEMAPVISITOR_H

#include "DenseVector.h"
#include "DenseMatrix.h"

#include "EuclideanMetric.h"

#include "SVD.h"
#include "RandomSVD.h"

#include <algorithm>
#include <map>

#include "GMRATree.h"
#include "IPCATree.h"
#include "Wasserstein.h"

template <typename TPrecision>
class NodeMapVisitor : public Visitor{
private:

  std::map<int, GMRANode * > nodes;
  int maxID;
  EuclideanMetric<TPrecision> l2;


public:
  
   NodeMapVisitor<TPrecision>(){
     maxID = 0;
   };

   void visit(GMRANode *node){
    nodes[maxID] = node; 
    maxID++; 
  };

  std::map<int, GMRANode * > &getNodes(){
    return nodes;
  };




};




#endif
