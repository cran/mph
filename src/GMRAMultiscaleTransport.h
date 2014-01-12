#ifndef GMRAMULTISCALETRANSPORT_H
#define GMRAMULTIScALETRANSPORT_H

#include "MultiscaleTransport.h"
#include "GMRATree.h"
#include "NodeDistance.h"
#include "GMRANeighborhood.h"

//#include "boost/functional.hpp"


template <typename TPrecision>
class GMRATransportNodeDecorator : public GMRANodeDecorator<TPrecision>{
  
  public:
    
    GMRATransportNodeDecorator(GMRANode<TPrecision> * node) : GMRANodeDecorator<TPrecision>(node){
      
    };

    virtual ~GMRATransportNodeDecorator(){};


    std::map<int, TransportNode<TPrecision> *> nodemap;

};



template <typename TPrecision>
class GMRATransportNodeBase : public TransportNode<TPrecision>{
  
  public:
    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;

    typedef typename GMRANode<TPrecision>::NodeVector NodeVector;
    typedef typename NodeVector::iterator NodeVectorIterator;
  
  
  protected:

    NodeDistance<TPrecision> *dist;
  

  private:
    GMRATransportNodeDecorator<TPrecision> *node;
    TPrecision localRadius;



  public:
  
    GMRATransportNodeBase(GMRATransportNodeDecorator<TPrecision> *n, NodeDistance<TPrecision> *d, int
        id) :  TransportNode<TPrecision>(id), node(n), dist(d){
      localRadius = -1;
    };


    virtual ~GMRATransportNodeBase(){};


    virtual TPrecision getLocalNodeRadius() {

      if(localRadius == -1){

        localRadius = 0;
        //static EuclideanMetric<TPrecision> metric;
        const TransportNodeVector &kids = this->getChildren();
        for(TransportNodeVectorCIterator it = kids.begin(); it!=kids.end(); it++){
          //TPrecision tmp = metric.distance(node->getCenter(),
          //    (*it)->getCenter() );
          TPrecision tmp = this->getTransportCost(*it, 1);
          if(localRadius < tmp){
            localRadius = tmp;
          } 
        }

      }
      return localRadius;
    
    };

   
    
    virtual TPrecision getNodeRadius() {
      return node->getRadius();
    };


    GMRATransportNodeDecorator<TPrecision> *getGMRANode() const{
      return node;
    };




    virtual std::vector<int> &getPoints(){
      return node->getPoints();
    };



/*
    virtual bool operator == (const TransportNode<TPrecision> &other) const{
       const GMRATransportNode<TPrecision> &o = (GMRATransportNode<TPrecision>&) other;
       return o.node == node;
    };
   

    virtual bool operator <  (const TransportNode<TPrecision> &other) const{
       GMRATransportNode<TPrecision> &o = (GMRATransportNode<TPrecision>&) other;
       return node < o.node;
    
    };
   

    virtual bool operator >  (const TransportNode<TPrecision> &other) const{
       const GMRATransportNode<TPrecision> &o = (GMRATransportNode<TPrecision>&) other;
       return node > o.node;
    };
*/

    virtual void print(){
      std::cout << "GMRATransportNode: " << node << " - " << node->getScale() << std::endl;
    };

     /*
     virtual size_t hashKey() const{
       //static boost::hash< GMRANode<TPrecision>* > hasher;
       return hash_pointer( (uint64_t) node);
     }
     */

  private:

     /*
     uint32_t hash_pointer( uint32_t a) const{
       a = (a ^ 61) ^ (a >> 16);
       a = a + (a << 3);
       a = a ^ (a >> 4);
       a = a * 0x27d4eb2d;
       a = a ^ (a >> 15);
       return a;
     };

     uint64_t hash_pointer(uint64_t key) const
     {
       key += ~(key << 32);
       key ^= (key >> 22);
       key += ~(key << 13);
       key ^= (key >> 8);
       key += (key << 3);
       key ^= (key >> 15);
       key += ~(key << 27);
       key ^= (key >> 31);
       return key;
     };
     */

};







template <typename TPrecision>
class GMRATransportNode : public GMRATransportNodeBase<TPrecision>{
  
  public:
    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;

    typedef typename GMRANode<TPrecision>::NodeVector NodeVector;
    typedef typename NodeVector::iterator NodeVectorIterator;
  
  

  public:
  
    GMRATransportNode(GMRATransportNodeDecorator<TPrecision> *n, NodeDistance<TPrecision> *d, int
        id) :  GMRATransportNodeBase<TPrecision>(n, d, id){
    };



    virtual TPrecision getTransportCost(TransportNode<TPrecision> *other, double p) const{
       GMRATransportNode<TPrecision> *o = dynamic_cast<
         GMRATransportNode<TPrecision>* >(other);

       GMRATransportNodeDecorator<TPrecision> *n1 = this->getGMRANode();
       GMRATransportNodeDecorator<TPrecision> *n2 = o->getGMRANode();
       return pow( this->dist->distance( n1->getDecoratedNode(), n2->getDecoratedNode() ), p );
    };

};

template <typename TPrecision>
class GMRATransportNodeMS : public GMRATransportNodeBase<TPrecision>{
  
  public:
    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename TransportNodeVector::const_iterator TransportNodeVectorCIterator;

    typedef typename GMRANode<TPrecision>::NodeVector NodeVector;
    typedef typename NodeVector::iterator NodeVectorIterator;
  
  

  public:
  
    GMRATransportNodeMS(GMRATransportNodeDecorator<TPrecision> *n, NodeDistance<TPrecision> *d, int
        id) :  GMRATransportNodeBase<TPrecision>(n, d, id){
    };


    virtual TPrecision getTransportCost(TransportNode<TPrecision> *other, double p) const{
      using namespace FortranLinalg;
       GMRATransportNodeBase<TPrecision> *o = dynamic_cast<
         GMRATransportNodeBase<TPrecision>* >(other);

       GMRATransportNodeBase<TPrecision> *p1 = dynamic_cast<
         GMRATransportNodeBase<TPrecision>* >( this->getParent() );
       GMRATransportNodeBase<TPrecision> *p2 = dynamic_cast<
         GMRATransportNodeBase<TPrecision>* >( o->getParent() );
       if(p1 == NULL || p2 == NULL){
         return pow( this->dist->distance(
               this->getGMRANode()->getDecoratedNode(),
               o->getGMRANode()->getDecoratedNode() ), p );
       }

       //TPrecision pCost = 0;
       //if(recursive){
       TPrecision pCost = p1->getTransportCost(p2, p);
       //}

       DenseVector<TPrecision> x1 = p1->getGMRANode()->getCenter();
       DenseVector<TPrecision> x2 = p2->getGMRANode()->getCenter();
       DenseVector<TPrecision> delta = Linalg<TPrecision>::Subtract(x2, x1);

       this->getGMRANode()->translate(delta);
       return pCost + pow(this->dist->distance(
             this->getGMRANode()->getDecoratedNode(),
             o->getGMRANode()->getDecoratedNode() ), p );

       Linalg<TPrecision>::Scale(delta, -1, delta);
       this->getGMRANode()->translate(delta);

       delta.deallocate();

    };  



};


template <typename TPrecision>
class GMRAMultiscaleTransportLevel : public MultiscaleTransportLevel<TPrecision>{


  public:

    typedef typename TransportNode<TPrecision>::TransportNodeVector TransportNodeVector;
    typedef typename TransportNodeVector::iterator TransportNodeVectorIterator;
    typedef typename GMRANeighborhood<TPrecision>::NeighborList NeighborList;
    typedef typename NeighborList::iterator NeighborListIterator;




  private:

    GMRATree<TPrecision> &tree;
    GMRANeighborhood<TPrecision> &neighborhood;
    NodeDistance<TPrecision> *dist;
    int scale;




  public:


    GMRAMultiscaleTransportLevel(GMRANeighborhood<TPrecision> &nh, int s,
        NodeDistance<TPrecision> *d ) : MultiscaleTransportLevel<TPrecision>(s),
    tree(nh.getTree()), neighborhood(nh),
        dist(d){
      scale = s;
    };





    virtual TransportNodeVector getNeighborhood(TransportNode<TPrecision> *node, TPrecision eps) const{
      GMRATransportNode<TPrecision> *n =
        dynamic_cast<GMRATransportNode<TPrecision> *>( node );
      GMRATransportNodeDecorator<TPrecision> *tNode = n->getGMRANode();
      //GMRANode<TPrecision> *gNode = tNode->getDecoratedNode();
    
       
      NeighborList gnodes;
      int nn = neighborhood.neighbors( tNode, eps, gnodes, scale );
     
      TransportNodeVector neighbors;
      neighbors.reserve(nn);
      for(NeighborListIterator it = gnodes.begin(); it != gnodes.end(); ++it){
        GMRATransportNodeDecorator<TPrecision> *tNode =
          dynamic_cast< GMRATransportNodeDecorator<TPrecision> *>( it->second );
        neighbors.push_back( tNode->nodemap[scale] );
      }
      return neighbors;

    };




    static std::vector< MultiscaleTransportLevel<TPrecision> *>
      buildTransportLevels(GMRANeighborhood<TPrecision> &nh,
          NodeDistance<TPrecision> *d, std::vector<TPrecision> &weights, bool
          multiscaleCost ){


        class MaxScale : public Visitor<TPrecision>{
          public:
            int maxScale;
            MaxScale(){ maxScale = 0; };

            void visit(GMRANode<TPrecision> *node){
              if(maxScale < node->getScale() ){
                maxScale = node->getScale();
              }
            };
        };


        //Decorate tree with transport node decorator
        GMRATree<TPrecision> &t = nh.getTree();
        GMRATransportNodeDecorator<TPrecision> *root = decorate( t.getRoot(), NULL );
        t.setRoot( root );

        //Compte maximal scale
        MaxScale ms;
        t.breadthFirstVisitor(&ms);
        
        std::vector< MultiscaleTransportLevel<TPrecision> *>
          levels(ms.maxScale+1);

        std::vector< int > idCounter(ms.maxScale+1, 0);

        for(int i=0; i < ms.maxScale+1; i++){
          levels[i] = new GMRAMultiscaleTransportLevel<TPrecision>(nh, i, d);
        }


        //populate multiscale transport levels 
        std::list< GMRANode<TPrecision> * > queue;
        std::list< TransportNode<TPrecision> * > parents;
        std::list< int > scales;
        scales.push_back(0);
        queue.push_back( root );
        parents.push_back( NULL );


        TPrecision nPoints = root->getPoints().size();
        while( !queue.empty() ){
          GMRANode<TPrecision> *node = queue.front();
          queue.pop_front();
          TransportNode<TPrecision> *parent = parents.front();
          parents.pop_front();
          int scale = scales.front();
          scales.pop_front();

          TPrecision mass = 0;
          if(weights.empty()){
            mass = node->getPoints().size() / nPoints;
          }
          else{
            std::vector<int> pts = node->getPoints();
            for(int i=0; i<pts.size(); i++){
              mass += weights[pts[i]];
            }
          }
          if(mass == 0){
            continue;
          }

          GMRATransportNodeDecorator<TPrecision> *dec =
            dynamic_cast<GMRATransportNodeDecorator<TPrecision> *>(node);
          GMRATransportNodeBase<TPrecision> *tNode;
         
          if(multiscaleCost){  
            tNode = new GMRATransportNodeMS<TPrecision>( dec, d, idCounter[scale] );
          }
          else{
            tNode = new GMRATransportNode<TPrecision>( dec, d, idCounter[scale] );
          }

          idCounter[scale] += 1;

          tNode->setMass(mass);
          levels[scale]->addNode(tNode);

          dec->nodemap[scale] = tNode;

          tNode->setParent(parent);
          if( parent != NULL ) {
            parent->addChild( tNode );
          }

          std::vector< GMRANode<TPrecision> * > &kids = node->getChildren();
          for(int i = 0; i<kids.size(); i++){
            queue.push_back(   kids[i] );
            parents.push_back( tNode   );
            scales.push_back(  scale+1 );
          }
         

          if( kids.empty() ){
            if(scale < ms.maxScale){
              queue.push_back(   node    );
              parents.push_back( tNode   );
              scales.push_back(  scale+1 );
            }
            else{
              tNode->addChild(tNode);
            }
          }
          
        }

        return levels; 
    };



    static GMRATransportNodeDecorator<TPrecision> *decorate(GMRANode<TPrecision>
        *node, GMRATransportNodeDecorator<TPrecision> *parent){
      using namespace FortranLinalg;
      
      GMRATransportNodeDecorator<TPrecision> *tNode = new GMRATransportNodeDecorator<TPrecision>( node );
      tNode->setParent(parent);

      std::vector< GMRANode<TPrecision>* > &children = node->getChildren();
      for(int i=0; i<children.size(); i++){ 
        children[i] =  decorate(children[i], tNode);
      }

      return tNode;

    };

};








#endif
