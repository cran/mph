#ifndef POTENTIALGMRANEIGHBORHOOD_H
#define POTENTIALGMRANEIGHBORHOOD_H

#include "GMRATree.h"

#include <list>
#include <set>
#include <vector>


#include "EuclideanMetric.h"
#include "NodeDistance.h"


template <typename TPrecision>
class PotentialInfo : public NodeInfo{
  public:
    PotentialInfo(){
      piMin = std::numeric_limits<TPrecision>::max();
      piMax = std::numeric_limits<TPrecision>::min();
    };

    TPrecision piMin;
    TPrecision piMax;
}


//Neighborhood search with potentials minima and maxima for each node and GMRA tree structure
template <typename TPrecision>
class PotentialGMRANeighborhood : public GMRANeighborhood<TPrecision>{
  
  private:
    typedef typename GMRANeighborhood<TPrecision>::Neighbor Neighbor;
    typedef typename GMRANeighborhood<TPrecision>::NeighborList NeighborList;
    typedef typename NeighborList::iterator NeighborListIterator;
    
    typedef typename GMRANode<TPrecision>::NodeVector  NodeVector;
    typedef typename NodeVector::iterator NodeVectorIterator;
    
    GMRATree<TPrecision> &tree;

    //Minimal and maximal node potential in each node
    std::vector<TPrecision> piMax;
    std::vector<TPrecision> piMin;
    
    //Distance measure between nodes
    NodeDistance<TPrecision> *dist;


  public:

    PotentialGMRANeighborhood(GMRATree<TPrecision> &t) : tree(t) {
      dist = new CenterNodeDistance<TPrecision>( new EuclideanMetric<TPrecision>() );
    };
        
    
    PotentialGMRANeighborhood(GMRATree<TPrecision> &t, NodeDistance<TPrecision> *d) : tree(t),
      dist(d){
    };
    
/*
    //find all points within epsilon
    std::list<int> epsilon(GMRANode<TPrecision> &x, TPrecision eps){
      std::list<int> &result;
      NeighborList candidates;
      collectNodes(x, eps, candidates);
      for(NeighborListIterator it = candidates.begin(); it !=  candidates.end();
          ++it){
        std::list<int> &points = it->getPoints();
        for(std:list<int>::iterator pit = points.begin(); pit<points.end();
            ++pit){
          DenseVector<TPrecision> &p = tree.getPoint(*pit);
          dist->distance(p, x);
          if(d < eps){
            result.push_back(*pit);
          }
        }
      }

      return result;
    }; 
  */ 


    //Find all nodes within epsilon at a certain scale
    void neighbors(GMRANode<TPrecision> *x, TPrecision eps,
        const std::set<GMRANode<TPrecision> *> &stopNodes,
        NeighborList &result) const{
      collectNodes(x, eps, result, stopNodes);
    }; 


    GMRATree<TPrecision> &getTree(){
      return tree;
    };


    void updateNodePotentials(DenseVector<TPrecision> &potentials){
      GMRANode<TPrecision> *node = tree.getRoot();
      
    };

  private:

    TPrecision reducedCost(GMRANode<TPrecision> *n1, GMRANode<TPrecision> *n2){
      TPrecision d = dist->distance(n1, n2);
      TPrecision e1 = n1->getRadius();
      TPrecision e2 = n2->getRadius();

      PotentialInfo *p1 = (PotentialInfo*) n1->getNodeInfo();
      PotentialInfo *p2 = (PotentialInfo*) n2->getNodeInfo();

      TPrecision delta1 = p1->piMax - p2->piMin;
      //TPrecision delta2 = p2->piMax - p1->piMin;


      //Minimal possibly distance - maximal possible potential difference
      return d - e1 - e2 - delta1;



    };



    void updatePotentials(GMRANode<TPrecision> *node, DenseVector<TPrecision> &potentials){
      std::vector< GMRANode<TPrecision>* > kids = node->getChildren();
      PotentialInfo *pi = new PotentialInfo();
      for(int i=0; i< kids.size(); i++){
        updatePotentials( kids[i] );
        PotentialInfo *p = (PotentialInfo *) kids[i]->getNodeInfo();
        pi->piMax = std::max(pi->piMax, p->piMax);
        pi->piMin = std::max(pi->piMin, p->piMin);
      }
      if(kids.size() == 0){
        std::vector<int> points = node->getPoints(); 
        pi->piMax = potentials[points[0]];
        pi->piMin = pi->piMax;
      }
      NodeInfo *tmp = node->getNodeInfo();
      if(tmp != NULL){
        delete tmp;
      }
      node->setNodeInfo(pi);
    };
       

    void collectNodes(GMRANode<TPrecision> *x, TPrecision eps, NeighborList
        &collected, const std::set<GMRANode<TPrecision> *> &stopNodes) const{

      
      NeighborList nodes;
      GMRANode<TPrecision> *root = tree.getRoot();
      TPrecision d = reducedCost(x, root);
      nodes.push_back( Neighbor(d, root) );
     
      int nVisited =0; 
      while( !nodes.empty() && nVisited < stopNodes.size() ){
        
        Neighbor &n = nodes.front();

        NodeVector kids = n.second->getChildren();
        
        if( stopNodes.find(n.second) != stopNodes.end() ){
          if(n.first <= eps){
            collected.push_back(n);
          }
          nVisited++;
        }
        else{
          //For each kid check if nearest neighbors within epsilon are possible.
          for(NodeVectorIterator it = kids.begin(); it != kids.end(); ++it){
            GMRANode<TPrecision> *kid = *it;
            TPrecision d = reducedCost(x, k);
            if(d <= kid->getRadius() + eps){
              nodes.push_back( Neighbor(d, kid) );
            }
          }
        }
        
        nodes.pop_front();
      }


    };




};









#endif
