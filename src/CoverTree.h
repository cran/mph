#ifndef COVERTREE_H
#define COVERTREE_H

#include "GMRATree.h"

#include <queue>
#include <list>
#include <iostream>
#include <fstream>

#include "Cover.H"
//#include "DisjointLists.H"
#include "EnlargeThreads.H"


typedef TPrecision REAL;

//Dimension for cover tree
int dim;


//Wrapper for cover class
template <typename TPrecision>
class CoverTree : public GMRATree<TPrecision>{

  private:
     Cover *cover; 


     void buildTree(CoverNodeAdapter *n){
        DLPtrList<CoverNode>* tmp = n->coverNode()->getChildren();
        for(DLPtrListNode<CoverNode> *current = tmp->first(); current != NULL; current =
            current->next() ){
          CoverNodeAdapter<TPrecision> *child = new
            CoverNodeAdapter<TPrecision>(current->getPtr() );
          n->addChild(child);
          buildTree(child);
        }

     };

  public:


    CoverTree(Cover *c, int d) : cover(c) {
      dim = d;
 
      //Construct GMRA tree structure based on cover 
      CoverNodeAdapter<TPrecision> *root = new CoverNodeAdapter<TPrecision>( c->getRoot() );
      buildTree( root );
      setRoot( root );
      
    };


    //
    void construct(TPrecision *X, int m, int n){
      throw("Not supported operation");
    };



    void add(TPrecision *X, int n){
      throw("Not supported operation");
    };



    std::vector<GMRANode *> getLeafPath(double *x) {
     
      Vector pt(x);
      std::vector<DescendNodePtr> nodes; 
      Cover::DescendList dlist;
      const Point **ptarr = new const Point*[1];
      double *dists = new double[1];

      cover->findNearest(pt, nodes, 1, 1, 0, dlist, ptarr, dists);
      
    
      delete [] ptarr;
      delete [] dists;

      std::vector<GMRANode *> path;
      for(DescenNode *current = nodes.first(); current != NULL;
          current=current->next()){
        CoverNode *cNode = current->getCoverNode();
        path.push_back( new CoverNodeAdapter(cNode) );
      }
      
      return path;
    }; 



};



template <typename TPrecision>
class CoverNodeAdapter : public GMRANode<TPrecision>{
  
  
  private:
  
    CoverNode *node;
    DenseVector<TPrecision> center;
  
  
  public:
  
    CoverNodeAdapter(CoverNode * cNode) : node(cNode) {
      Vector *v = (Vector*)node->getPoint();
      center = DenseVector<TPrecision>(dim, v->getPoint() );
    };
  

    ~CoverNodeAdapter(){
      center.deallocate();
    };
    
    std::vector<int> &getPoints(){
      std::vector<int> points;
      std::list<CoverNode *> kids;
      kids.insert( node );
      while(!kids.empty() ){
        CoverNode *n = kids.front();
        kids.pop_front();

        //Obsolete
        IndexedPoint *ip = (IndexedPoint*) n->getPoint();
        points.push_back(ip->getIndex());

        DLPtrList<CoverNode>* tmp = n->getChildren();
        for(DLPtrListNode<CoverNode> *current = tmp->first(); current != NULL; current =
            current->next() ){
          kids.push_back(current->getPtr());
        }

      }
      return points;
    };
   


    //Partition representative and radius 
    virtual FortranLinalg::DenseVector<TPrecision> &getCenter(){
      return center;      
    };

    virtual TPrecision getRadius(){
      return node->getRadius();
    };


};



#endif
