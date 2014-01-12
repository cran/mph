//author: Samuel Gerber


#ifndef MULTISCALEMAPPINGRIPS_H
#define MULTISCALEMAPPINGRIPS_H


#include "MappingRips2.h"
#include "GMRATree.h"
#include "GMRANeighborhood.h"
#include "EuclideanMetric.h"
#include "LinalgIO.h"

#include <map>
#include <vector>
#include <utility>
#include <limits>
#include <list>
#include <vector>


#include "RelativeGMRANeighborhood.h"
#include "MappingPersistentHomology.h"


template<typename TPrecision>
class MultiscaleMappingRips{


  public:

    typedef typename MappingRips2<TPrecision>::IFiltration IFiltration;
    typedef typename MappingRips2<TPrecision>::IFiltrationIterator IFiltrationIterator;

    typedef typename MappingRips2<TPrecision>::Filtration Filtration;
    typedef typename MappingRips2<TPrecision>::FiltrationIterator FiltrationIterator;

    typedef typename MappingPersistentHomology<TPrecision>::Event Event;
    typedef typename MappingPersistentHomology<TPrecision>::Events Eventss;
    typedef typename MappingPersistentHomology<TPrecision>::EventsIterator EventsIterator;


    typedef typename MappingRips2<TPrecision>::NeighborMap NeighborMap;
    typedef typename MappingRips2<TPrecision>::Neighbors Neighbors;




  private:


    FortranLinalg::DenseMatrix<TPrecision> X;
    GMRATree<TPrecision> &gmra;
    //std::map< int, Filtration > filtrations;
    std::map< int, FortranLinalg::DenseMatrix<TPrecision> > diagrams;




  public:

    MultiscaleMappingRips(FortranLinalg::DenseMatrix<TPrecision> &Xin, GMRATree<TPrecision> &tree) : X(Xin), gmra(tree){ 
    };



    void run( int maxD, bool singleScale = false){

      using namespace FortranLinalg;
      //Find all leave nodes
      SetupRips setup;
      gmra.breadthFirstVisitor( &setup );


      GenericGMRANeighborhood<TPrecision> nh( gmra );

      std::set< GMRANode<TPrecision> * > current = setup.leaves;
      std::map< GMRANode<TPrecision> *, int > currentIndexes;
      DenseMatrix<TPrecision> Xcurrent = getCenters(current, currentIndexes);
      std::vector<int> mapping;
      IFiltration mapped;

      TPrecision prevTime = 0;
      TPrecision prevEpsilon = 0; 
      int scale = 0;
      while( current.size() > 1){



        clock_t t1 = clock();


        //Rips at current scale
        //TPrecision radius = getNthRadius(current, std::min(current.size()/2.0,
        //     10.0) );
        //TPrecision radius = 1.2 * getMinRadius(current);
        TPrecision radius = getNthRadius(current, std::min(1.0 +
              current.size()/2.5, current.size()-1.0) );
        //TPrecision moveRadius = 1.5 * getMinRelativeRadius(current);
        //radius = std::max(moveRadius, radius);

        TPrecision epsilon = radius;

        if( singleScale ){
          radius = std::numeric_limits<TPrecision>::max() / 3;
        }



        for(typename std::set<GMRANode<TPrecision> *>::iterator it = current.begin(); it
            != current.end(); ++it){
          (*it)->setStop( true );
        }

        //Compute neighborhood info
        Neighbors N(current.size());
        int index = 0;
        for(typename std::set<GMRANode<TPrecision> *>::iterator it = current.begin(); it
            != current.end(); ++it, ++index){

          NeighborMap &nl = N[index];

          typename GMRANeighborhood<TPrecision>::NeighborList nList;
          nh.neighbors(*it, 3*radius, nList);

          for(typename GMRANeighborhood<TPrecision>::NeighborListIterator nIt
              =nList.begin(); nIt != nList.end(); ++nIt){
              int nInd = currentIndexes[nIt->second];
              nl[nInd] = nIt->first;
              N[nInd][index] = nIt->first;
          } 
        }

        for(typename std::set<GMRANode<TPrecision> *>::iterator it = current.begin(); it
            != current.end(); ++it){
          (*it)->setStop( false );
        }



        clock_t t2 =  clock();
#if VERBOSE
        std::cout << "Find neighboors time: " << ((double)t2 - t1) / CLOCKS_PER_SEC << std::endl; 
#endif


        //Run mapping rips with mapped filtration from previous scale
        MappingRips2<TPrecision> rips(Xcurrent);
        rips.setMappedFiltration( mapped );
        int maxSize = 100000000;
        rips.run( N, prevTime, 3*radius, maxSize, maxD, scale);
        Filtration filt = rips.getFiltration();

        clock_t t3 =  clock();

#if VERBOSE
        std::cout << "Rips time: " << ((double)t3 - t2) / CLOCKS_PER_SEC << std::endl; 
#endif


        
/* 
        //Debug
        std::stringstream ss;
        ss << "vertices-scale-" << scale << ".data";
        LinalgIO<TPrecision>::writeMatrix(ss.str(), Xcurrent);


        //Debug
        std::stringstream ss1;
        ss1 << "edges-scale-" << scale << "-ghost.data";
        writeEdges(mapped, ss1.str());


        //Debug
        std::stringstream ss2;
        ss2 << "edges-scale-" << scale << ".data";
        writeEdges( filt, ss2.str() );
*/



        //Compute mapping to next scale
        DenseMatrix<TPrecision> Xnext; 
        std::set< GMRANode<TPrecision> * > next;
        std::map< GMRANode<TPrecision>*, int > nextIndexes;
          
        TPrecision mapTime = radius;
        if( rips.isTruncated() ){
          TPrecision ripsTime = rips.getMaxTime();
#if VERBOSE
          std::cout << "Rips max length: " << ripsTime << std::endl;
#endif
          mapTime = std::min( ripsTime/3.0, radius ); 
        }

        if( !singleScale ){
          collectNodes( current, next, mapTime );
          Xnext = getCenters( next, nextIndexes );
          mapVertices( currentIndexes, nextIndexes, mapping, mapTime );
        }
        else{
          mapTime = std::numeric_limits<TPrecision>::max()/2.0;
        }

        clock_t t4 =  clock();
#if VERBOSE
        std::cout << "Compute mapping time: " << ((double)t4 - t3) /CLOCKS_PER_SEC << std::endl; 
#endif



        //Map filtration to next scale scale
        if( !singleScale ){
          mapped.clear();
          mapFiltration(filt, mapping, mapped, scale);
        }

        clock_t t5 =  clock();
#if VERBOSE
        std::cout << "Map filtration  time: " << ((double)t5 - t4) /CLOCKS_PER_SEC << std::endl; 
#endif





        //Do reduction according to mapped filtration        
        MappingPersistentHomology<TPrecision> homology;
        homology.run( filt, 2*mapTime, maxD );
        diagrams[scale] = homology.getDiagram(); 



        clock_t t6 =  clock();


        TPrecision maxTime = homology.getMaximalTime(); 

#if VERBOSE
        std::cout << "Reduction time: " << ((double)t6 - t5) /CLOCKS_PER_SEC << std::endl; 

        std::cout << "#Vertices : " << Xcurrent.N() << std::endl;
        std::cout << "Filtration size: " << filt.size() << std::endl;
        std::cout << "Radius: " << radius << std::endl;
        std::cout << "Max Time: " << maxTime << std::endl << std::endl;
#endif





        //Setup for next scale
        Xcurrent.deallocate();
        Xcurrent = Xnext; 

        current = next;
        currentIndexes = nextIndexes;
        prevTime = maxTime;
        prevEpsilon = epsilon;
        scale++;

        if(singleScale){
          break;
        }



      }

      Xcurrent.deallocate();

    };


    std::map<int, FortranLinalg::DenseMatrix<TPrecision> > &getDiagrams(){
      return diagrams;
    };


    //  std::map<int, IFiltration> &getInverseFiltrations(){
    //    return filtrations;
    //  };






  private:



    class SetupRips : public Visitor<TPrecision>{

      public:

        std::set<GMRANode<TPrecision> *> leaves;

        void visit(GMRANode<TPrecision> *node){
          //node->setNodeInfo( new RipsInfo() );
          if(node->getChildren().size() == 0){
            leaves.insert(node);
          }
        };

    };







    void collectNodes(std::set< GMRANode<TPrecision> * > &nodes,
        std::set<GMRANode<TPrecision> *> &pnodes, TPrecision r){

      for(typename std::set< GMRANode<TPrecision> * >::iterator it = nodes.begin(); it!=
          nodes.end(); it++){
        pnodes.insert( getMappedNode(*it, r) );
      }

    };




    GMRANode<TPrecision> *getMappedNode(GMRANode<TPrecision> *node, TPrecision
        r){

      GMRANode<TPrecision> *current = node;
      while( getRadius(current) <= r ){
        GMRANode<TPrecision> *parent = current->getParent();
        if(parent != NULL){
          current = parent;
        }
        else{
          break;
        }
      }

      return current;
    };


    TPrecision getRadius(GMRANode<TPrecision> *node){
      //return getParentRadius(node);

      TPrecision r = 0;
      GMRANode<TPrecision> *p = node->getParent();
      if(p == NULL){
        r = node->getRadius();
      }
      else{
        std::vector< GMRANode<TPrecision>* > &kids = p->getChildren();

        for(int i=0; i<kids.size(); i++){
          r = std::max(r, kids[i]->getRadius() );
        }
      }

      return std::max( r, getRelativeRadius(node) );
    };




    TPrecision getParentRadius(GMRANode<TPrecision> *node){
      GMRANode<TPrecision> *p = node->getParent();
      if(p == NULL){
        return node->getRadius();
      }
      return p->getRadius();
    };



    TPrecision getRelativeRadius(GMRANode<TPrecision> *from,
        GMRANode<TPrecision> *to){
      using namespace FortranLinalg;

      GMRANode<TPrecision> *p = to->getParent();
      if(p == NULL){
        return std::numeric_limits<TPrecision>::max();
      }

      static EuclideanMetric<TPrecision> l2norm;
      DenseVector<TPrecision> x1 = p->getCenter();
      DenseVector<TPrecision> x2 = from->getCenter();
      TPrecision r = l2norm.distance(x1, x2);

      return r;
    };


    TPrecision getRelativeRadius(GMRANode<TPrecision> *node){
      using namespace FortranLinalg;

      GMRANode<TPrecision> *p = node->getParent();
      if(p == NULL){
        return std::numeric_limits<TPrecision>::max();
      }

      TPrecision r = 0;

      std::vector< GMRANode<TPrecision>* > &kids = p->getChildren();

      static EuclideanMetric<TPrecision> l2norm;
      DenseVector<TPrecision> x1 = p->getCenter();
      for(int i=0; i<kids.size(); i++){
          DenseVector<TPrecision> x2 = kids[i]->getCenter();
          TPrecision tmp = l2norm.distance(x1, x2);
          r = std::max(r, tmp);
      }

      return r;
    };




    TPrecision getMinRadius( std::set< GMRANode<TPrecision> * > &nodes){
      TPrecision r = std::numeric_limits<TPrecision>::max();
      for(typename std::set< GMRANode<TPrecision> * >::iterator it = nodes.begin(); it!=
          nodes.end(); it++){
        GMRANode<TPrecision> *n = *it;
        r = std::min(r, getRadius(n) );
      }
      return r;
    };



    TPrecision getMinRelativeRadius( std::set< GMRANode<TPrecision> * > &nodes){
      TPrecision r = std::numeric_limits<TPrecision>::max();
      for(typename std::set< GMRANode<TPrecision> * >::iterator it = nodes.begin(); it!=
          nodes.end(); it++){
        GMRANode<TPrecision> *n = *it;
        r = std::min(r, getRelativeRadius(n) );
      }
      return r;
    };





    TPrecision getNthRadius( std::set< GMRANode<TPrecision> * > &nodes, int n){
      std::vector<double> radii(nodes.size());
      int index = 0;
      for(typename std::set< GMRANode<TPrecision> * >::iterator it = nodes.begin(); it!=
          nodes.end(); it++, index++){
        GMRANode<TPrecision> *n = *it;
        radii[index] =  getRadius(n);
      }
      std::sort( radii.begin(), radii.end() );
      return radii[n];
    };





    void mapVertices(std::map<GMRANode<TPrecision> *, int> &X,
        std::map<GMRANode<TPrecision> *,int> &Y, std::vector<int>
        &mapping, TPrecision r){

      mapping.clear();
      mapping.resize(X.size());

      for(typename std::map<GMRANode<TPrecision> *, int>::iterator it = X.begin(); it !=
          X.end(); ++it){ 
        mapping[it->second] = Y[ getMappedNode(it->first, r) ];
      }

    };





    FortranLinalg::DenseMatrix<TPrecision> getCenters(std::set<GMRANode<TPrecision> *> current,
        std::map<GMRANode<TPrecision> *, int> &indexes){
      using namespace FortranLinalg;
      //Collect centers for this scale
      DenseMatrix<TPrecision> Xs(X.M(), current.size());
      int index= 0;
      for(typename std::set<GMRANode<TPrecision> *>::iterator it = current.begin(); it !=
          current.end(); ++it, ++index){
        GMRANode<TPrecision> *node = (*it);
        indexes[node] = index;        
        DenseVector<TPrecision> center = node->getCenter();
        Linalg<TPrecision>::SetColumn(Xs, index, center);
      }
      return Xs;
    };



    //Mapp filtration according to ampping
    void mapFiltration(Filtration &filt, std::vector<int> &mapping, IFiltration
        &mapped, int scale){

      for( FiltrationIterator it = filt.begin(); it != filt.end(); it++){
        Simplex sNew(scale);

        //
        Simplex &s = it->simplex;

        for(std::set<int>::iterator vIt = s.vertices.begin(); vIt !=
            s.vertices.end(); ++vIt){
          int v = *vIt;
          sNew.vertices.insert( mapping[v] );
        }
        if(s.vertices.size() == sNew.vertices.size() ){
          s.setMapsTo(sNew.vertices);
        }

        if(sNew.vertices.size() > 0 ){

          IFiltrationIterator it2 = mapped.find(sNew);
          if(it2 == mapped.end()){
            sNew.setScale(s.getScale());
            mapped[sNew] = it->time;
          }
          else{
            const Simplex &s2 = it2->first;
            if(s.getScale() < s2.getScale()){
              mapped.erase(sNew);
              sNew.setScale( s.getScale() );
              mapped[sNew] = it->time;
            }
            else if(s.getScale() == s2.getScale() ){
              it2->second =  std::max(it->time, it2->second);
            }
            //else leave old one
          }
        }
        else{
#if VERBOSE
          std::cout << "Huuuummmmm" << std::endl;
#endif
        }
      }

      clock_t t2 = clock();

      //Make sure the order of simplicies is preserved. It can happen that a
      //lower order Simplex has a larger time stamp
      makeConsistent(mapped);
#if VERBOSE
      std::cout << "Mapped size: " << mapped.size() << std::endl;
#endif

    };



    void makeConsistent(IFiltration &ifilt){

      for(IFiltrationIterator it = ifilt.begin(); it != ifilt.end(); ++it){

        Simplex s = it->first;
        TPrecision tMax = it->second;

        if(s.vertices.size() > 1){
          //check maximal time of face addition
          std::list<Simplex> sList = s.getFaces();
          for(std::list<Simplex>::iterator sIt = sList.begin(); sIt !=
              sList.end(); ++sIt){

            Simplex &f = *sIt;

            IFiltrationIterator fIt = ifilt.find(f);
            if(fIt != ifilt.end() ){
              const Simplex &f2 = fIt->first;
              if( f2.getScale() <= s.getScale() ){
                if(tMax < fIt->second){
                  tMax = fIt->second;
                }
              }
              else{
#if VERBOSE
                std::cout << "woof" << std::endl;
#endif
              }
            }
            else{
#if VERBOSE
              std::cout << "bark" << std::endl;
#endif
            }
          }

          it->second = tMax;
        }
      }

    };





    //----Debug ----//

    struct Triple{
      Triple(int s, int index1, int index2):scale(s), i1(index1), i2(index2){};
      int scale;
      int i1;
      int i2;
    };


    //Debug method
    void writeEdges(IFiltration &ifilt, std::string filename){
      using namespace FortranLinalg;
      std::list< Triple  > edges;
      for(IFiltrationIterator it = ifilt.begin(); it != ifilt.end(); ++it){
        const Simplex &s = it->first;
        if(s.vertices.size() == 2){
          int i1 = *s.vertices.begin();
          int i2 = *s.vertices.rbegin();
          Triple edge(s.getScale(), i1, i2);
          edges.push_back(edge);
        }
      }

      DenseMatrix<TPrecision> E(3, edges.size());
      int index = 0;
      for(typename std::list< Triple >::iterator it = edges.begin(); it!=
          edges.end(); ++it, ++index){
        E(0, index) = it->i1;
        E(1, index) = it->i2;
        E(2, index) = it->scale;
      }
      LinalgIO<TPrecision>::writeMatrix(filename, E);
      E.deallocate();

    };


    //Debug method
    void writeEdges(Filtration &filt, std::string filename){
      using namespace FortranLinalg;
      std::list< Triple  > edges;
      for(FiltrationIterator it = filt.begin(); it != filt.end(); ++it){
        const Simplex &s = it->simplex;
        if(s.vertices.size() == 2){
          int i1 = *s.vertices.begin();
          int i2 = *s.vertices.rbegin();
          Triple edge(s.getScale(), i1, i2);
          edges.push_back(edge);
        }
      }

      DenseMatrix<TPrecision> E(3, edges.size());
      int index = 0;
      for(typename std::list< Triple >::iterator it = edges.begin(); it!=
          edges.end(); ++it, ++index){
        E(0, index) = it->i1;
        E(1, index) = it->i2;
        E(2, index) = it->scale;
      }
      LinalgIO<TPrecision>::writeMatrix(filename, E);
      E.deallocate();
    };


};

#endif 
