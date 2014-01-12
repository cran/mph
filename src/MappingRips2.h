//author: Samuel Gerber

#ifndef MAPPINGRIPS2_H
#define MAPPINGRIPS2_H


#include <map>
#include <set>
#include <vector>
#include <utility>
#include <limits>
#include <list>

#include "Simplex.h"

#include "DenseMatrix.h"
#include "Linalg.h"
#include "EuclideanMetric.h"



//Computes a Rips filtration F2 given a metric space X2, a filtration F1 on
//X1 and a mapping g from X1 to X2. The filtration is computed by first
//inserting the simplicies from F1 that map to F2 according to their time in
//filtration F1 and then adding according to Rips on F2 any remaining simplicies
//S with times max( max( time( F1 ), time(S) ) )

template<typename TPrecision>
class MappingRips2{


  public:
    typedef typename std::map<Simplex, TPrecision> IFiltration;
    typedef typename IFiltration::iterator IFiltrationIterator;
    typedef typename IFiltration::reverse_iterator IFiltrationRIterator;


    typedef typename std::map<int, TPrecision> NeighborMap;
    typedef typename NeighborMap::iterator NeighborMapIterator;
    typedef typename std::vector<NeighborMap> Neighbors; 

    class FiltrationEntry{
      public:
        Simplex simplex;
        TPrecision time;


        FiltrationEntry(const Simplex &s, TPrecision t):simplex(s), time(t){
        }; 

        bool operator == (const FiltrationEntry& other) const{
          return time== other.time && simplex == other.simplex; 
        };

        bool operator < (const FiltrationEntry& other) const{
          if(time < other.time){
            return true;
          }
          else if(time > other.time){
            return false;
          }
          else{
            return simplex < other.simplex;
          }
        };

        bool operator > (const FiltrationEntry& other) const{
          return other < *this;
        };
    };



    typedef typename std::vector< FiltrationEntry  > Filtration;
    typedef typename Filtration::iterator FiltrationIterator;




  private:
    //Simplex and adding time
    Filtration filtration;
    Filtration tmpFiltration;
    IFiltration mappedFiltration;


    //Data
    FortranLinalg::DenseMatrix<TPrecision> X;
    EuclideanMetric<TPrecision> metric;

    bool truncated; 


  public:


    MappingRips2(FortranLinalg::DenseMatrix<TPrecision> &Xin) : X(Xin){ 
      truncated = false;
    };


    ~MappingRips2(){

    };



    //Compute the filtration based on the nearest neighbor structure on X2, and
    //the minimal inserttion tim minTime = max(F1). Compute the Rips up to time
    //maxTime and with simplicial dimension up to maxD. Truncate the Rips at
    //maxSize. Annotate the added the newly added simplicies, i.e. not mapped
    //from F1, as coming from Filtration scale. This information is required to
    //run a multicsale persistent homology.
    //Before this is run setMappedFiltration (setting F1) should be called. If
    //no mapped filtration is set this is simply doing Rips with some overhead
    //computations.
    void run(Neighbors &N, TPrecision minTime, TPrecision maxTime, int maxSize,
        int maxD, int scale){

      clock_t t0 = clock();
      tmpFiltration.clear();
      tmpFiltration.reserve(N.size()*maxD);

      std::vector< std::set<int> > edges;
      //Build rips based on nn graph and mapped filtration
      for(unsigned int i=0; i < N.size(); i++){
        Simplex simplex(scale);
        simplex.vertices.insert(i);
        FiltrationEntry fe(simplex, 0);

        NeighborMap &knn = N[i];
        NeighborMap lower;
        lowerNeighbors(i, knn, lower);
        addCofaces(fe, N, lower, maxD, minTime, maxTime);    
      }

      /*
         std::cout << " rips inc. duplicates size: " << nTotal << std::endl;
         clock_t t1 = clock();
         std::sort(tmpFiltration.begin(), tmpFiltration.end() );
         clock_t t2 = clock();
         FiltrationIterator it = std::unique (tmpFiltration.begin(), tmpFiltration.end());
         tmpFiltration.erase( it, tmpFiltration.end() );
         clock_t t3 = clock();
       */ 

      //Combine the mapped and the new filtration
      filtration.clear();
      filtration.reserve( mappedFiltration.size() + tmpFiltration.size() );

      for(IFiltrationIterator it = mappedFiltration.begin(); it !=
          mappedFiltration.end(); ++it){
        Simplex s = it->first;
        s.setMapped(true);
        filtration.push_back( FiltrationEntry(s, it->second) );
      }

      filtration.insert(filtration.end(), tmpFiltration.begin(),
          tmpFiltration.end() );

      clock_t t1 = clock();
      std::sort(filtration.begin(), filtration.end() );
      clock_t t2 = clock();
      //FiltrationIterator it = std::unique (filtration.begin(), filtration.end());
      //filtration.erase( it, filtration.end() );
      clock_t t3 = clock();


#if VERBOSE
      std::cout << " rips inc. duplicates size: " << filtration.size() << std::endl;
      std::cout << " rips unique size: " << filtration.size() << std::endl;

      std::cout << " build time: " << ((double) t1-t0) / CLOCKS_PER_SEC << std::endl;
      std::cout << " sort time: " << ((double) t2-t1) / CLOCKS_PER_SEC << std::endl;
      std::cout << " unique time: " << ((double) t3-t2) / CLOCKS_PER_SEC << std::endl;
#endif

      if(filtration.size() > maxSize + X.N()){
        FiltrationIterator fIt  = filtration.begin();
        std::advance(fIt,maxSize+X.N());
        filtration.erase( fIt, filtration.end() );
        truncated = true;
      }

    };


    bool isTruncated(){
      return truncated;
    };

    Filtration &getFiltration(){
      return filtration;
    };


    TPrecision getMaxTime(){
      return filtration.back().time;
    };

    void setMappedFiltration(IFiltration &f){
      mappedFiltration=f;
    };



  private:


    void lowerNeighbors(int index, NeighborMap &nn, NeighborMap &lower){
      lower.clear();
      for(NeighborMapIterator it = nn.begin(); it != nn.end(); ++it){
        if(it->first < index){
          lower.insert(*it);
        }
      } 
    };




    void addCofaces(FiltrationEntry &f, Neighbors &N, NeighborMap &nn, int maxD,
        double minTime, double maxTime){

      if(f.time > maxTime){
        return;
      }

      if(mappedFiltration.find(f.simplex) == mappedFiltration.end() ){
        tmpFiltration.push_back(f);
      }
      if(f.simplex.vertices.size() > maxD ){
        return;
      }


      for(NeighborMapIterator it = nn.begin(); it != nn.end(); ++it){

        //Create new simplex a new entrance time
        Simplex s = f.simplex;
        double time = f.time;
        NeighborMap &knn = N[it->first];
        for(std::set<int>::iterator vit = f.simplex.vertices.begin(); vit !=
            f.simplex.vertices.end(); ++vit){
          double tmp = knn[*vit];
          if(tmp > time){
            time = tmp;
          }
        }
        if(time < minTime){
          time = minTime;
        }

        s.vertices.insert(it->first);
        FiltrationEntry fe(s, time);

        //Compute lower neighbors intersections
        NeighborMap lower;
        lowerNeighbors(it->first, knn, lower);
        NeighborMap res;
        std::set_intersection(lower.begin(), lower.end(), nn.begin(), nn.end(),
            std::inserter( res, res.begin()), lower.value_comp() );


        addCofaces(fe, N, res, maxD, minTime, maxTime);

      }

    };



};

#endif 

