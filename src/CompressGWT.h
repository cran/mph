#ifndef COMPRESSGWT_H
#define COMPRESSGWT_H

#include "DenseVector.h"
#include "DenseMatrix.h"
#include "PlaneAngles.h"

#include "SVD.h"

#include <algorithm>

//Interface for GWT
template <typename TPrecision> 
class CompressGWT : public Visitor{

public:
    std::map<int, DenseMatrix<TPrecision> > dict;
    typename std::map<int, DenseMatrix<TPrecision> >::iterator dictIt;
    
    std::map<int, std::list<GWTNode<TPrecision> *> > dict2node;   
    
    DenseMatrix<TPrecision> &X; 
    TPrecision threshold;



public:
  
  CompressGWT(DenseMatrix<TPrecision> &pts, TPrecision t) : X(pts), threshold(t){
  };

  void visit(GMRANode *gmraNode){
 
      GWTNode<TPrecision> *node = (GWTNode<TPrecision> *) gmraNode->getGWTInfo();
      
      DenseMatrix<TPrecision> phi1 = node->phi;
      bool newElement = true;
      for(dictIt = dict.begin(); dictIt != dict.end(); ++dictIt){
        DenseMatrix<TPrecision> phi2 = dictIt->second;
        DenseVector<TPrecision> a = PlaneAngles<TPrecision>::angles(phi1, phi2);
        std::cout << phi1.N() << " : " << a(0) << std::endl;
        if(a(0) < threshold){
          newElement = false;
          std::list< GMRANode *> nodes;
          nodes.push_back(gmraNode);
          DenseMatrix<TPrecision> M = createMatrix(nodes);
          SVD<TPrecision> svd(M);
          dictIt->second = Linalg<TPrecision>::ExtractColumns(svd.U, 0, std::max( phi1.N(), phi2.N() ) ); 
          svd.deallocate();
          M.deallocate();

          for(std::list<GMRANode *>::iterator it = nodes.begin(); it != nodes.end(); ++it){
            GWTNode<TPrecision> *gn = (GWTNode<TPrecision> *) (*it)->getGWTInfo();
            gn->phi.deallocate();
            gn->phi = dictIt->second;
          } 
          break;
        } 
        a.deallocate();
      }  
      if(newElement){
        int dictID = dict.size();
        dict[dictID] = phi1;
        dict2node[dictID].push_back(node);
      }
      

  };



private:

  DenseMatrix<TPrecision> createMatrix(std::list< GMRANode*> nodes){
   
    int nPoints = 0;
    for(std::list< GMRANode *>::iterator it = nodes.begin(); it != nodes.end(); ++it){
      std::vector<int> &pts = (*it)->getPoints();
      nPoints += pts.size();
    }

    DenseMatrix<TPrecision> M(X.M(), nPoints);
    int index = 0;
    for(std::list< GMRANode *>::iterator it = nodes.begin(); it != nodes.end(); ++it){
      std::vector<int> &pts = (*it)->getPoints();
      for(std::vector<int>::iterator pit = pts.begin(); pit != pts.end(); ++pit, ++index){
        Linalg<TPrecision>::SetColumn(M, index, X, *pit);
        GWTNode<TPrecision> *gn = (GWTNode<TPrecision> *) (*it)->getGWTInfo();
        Linalg<TPrecision>::SubtractColumn(M, index, gn->center, M);
      }
    }
    return M;
  };

};

#endif



