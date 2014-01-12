#ifndef SIMILARITYVISITOR_H
#define SIMILARITYVISITOR_H


#include "GMRATree.h"
#include "Linalg.h"

template <typename TPrecision>
class SimilarityVisitor : public Visitor<TPrecision>{


public:


  virtual FortranLinalg::DenseMatrix<TPrecision> distances() = 0;

  //virtual TPrecision distance(GWTNode<TPrecision> *n1, GWTNode<TPrecision> *n2) = 0;



};




#endif
