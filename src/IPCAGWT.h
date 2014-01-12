#ifndef IPCAGWT_H
#define IPCAGWT_H

#include "DenseVector.h"
#include "DenseMatrix.h"

#include "SVD.h"
#include "RandomSVD.h"

#include "GWT.h"
#include "IPCATree.h"



template <typename TPrecision> 
class IPCAGWTNode : public GWTNode<TPrecision>{
  private:
    IPCANode<TPrecision> *node;

  public:

    IPCAGWTNode<TPrecision>(int id, GWTNode<TPrecision> *p, IPCANode<TPrecision>
        *n ) : GWTNode<TPrecision>(n->phi, n->center, id, p, n), node(n) {
    };

    virtual ~IPCAGWTNode(){};


    virtual FortranLinalg::DenseMatrix<TPrecision> &getPhi(){
      return node->phi;
    };

    virtual FortranLinalg::DenseVector<TPrecision> &getSigma(){
      return node->sigma;
    };

    virtual FortranLinalg::DenseVector<TPrecision> &getCenter(){
      return node->center;
    };

};





template <typename TPrecision> 
class IPCAGWT : public GWT<TPrecision>{

  private:

    int nodeID; 

    GWTNode<TPrecision> *decorate(GMRANode<TPrecision> *gnode,
        GWTNode<TPrecision> *parent){
      using namespace FortranLinalg;
      
      IPCANode<TPrecision> *node = dynamic_cast< IPCANode<TPrecision> *>( gnode );
      GWTNode<TPrecision> *gwt = new IPCAGWTNode<TPrecision>(nodeID, parent, node);
      ++nodeID;
      gwt->setParent(parent);

      std::vector< GMRANode<TPrecision>* > &children = node->getChildren();
      for( int i=0; i < children.size(); i++ ){ 
        children[i] =  decorate( children[i], gwt );
      }

      return gwt;

    };



  protected:

    void setupGWT(GMRATree<TPrecision> *tree){
      nodeID= 0;
      tree->setRoot( decorate(tree->getRoot(), NULL) );
    };  


};


#endif
