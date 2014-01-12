//author: Samuel Gerber

#ifndef MAPPINGRIPS_H
#define MAPPINGRIPS_H


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
class MappingRips{


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


    MappingRips(FortranLinalg::DenseMatrix<TPrecision> &Xin) : X(Xin){ 
      truncated = false;
    };


    ~MappingRips(){

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

      int nAdded = 0;
      int nTotal = 0;
      //Build rips based on nn graph and potential previous filtration
      for(unsigned int i=0; i < N.size(); i++){
        //Add all simplicies formed by the neighbors of this vertex if they are
        //within the nearest neighbor graph and within distance alpha

        int nAdd = add(i, N, minTime, maxTime, maxD, scale);
        nAdded += nAdd;
        nTotal += nAdd;
        if(nAdded > 500000){
          
          std::sort(tmpFiltration.begin(), tmpFiltration.end());
          
          FiltrationIterator it = std::unique (tmpFiltration.begin(), tmpFiltration.end());
          tmpFiltration.erase( it, tmpFiltration.end() );
          
          nAdded = 0;
        }
      
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
      FiltrationIterator it = std::unique (filtration.begin(), filtration.end());
      filtration.erase( it, filtration.end() );
      clock_t t3 = clock();
       

#if VERBOSE
      std::cout << " rips inc. duplicates size: " << nTotal << std::endl;
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


   TPrecision getInsertTime(Simplex &s, Neighbors &NN){

     TPrecision time = 0;

     if(s.vertices.size() == 2){
       std::set<int>::iterator vIt = s.vertices.begin();
       int i1 = *vIt;
       ++vIt;
       int i2 = *vIt;

       NeighborMapIterator nIt = NN[i1].find(i2);
       //simplex is not in nearest neighbor map
       if(nIt ==  NN[i1].end()){
         return std::numeric_limits<TPrecision>::max();
       }
       time = nIt->second;
     }
     else{

       std::list<Simplex> faces = s.getFaces(1);
       for(std::list<Simplex>::iterator it = faces.begin(); it != faces.end();
           it++){

         Simplex &f = *it;
         std::set<int>::iterator vIt = f.vertices.begin();
         int i1 = *vIt;
         ++vIt;
         int i2 = *vIt;
         NeighborMapIterator nIt = NN[i1].find(i2);
         //simplex is not in nearest neighbor map
         if(nIt ==  NN[i1].end()){
           return std::numeric_limits<TPrecision>::max();
         }
         TPrecision d = nIt->second;
         if(time < d){
           time = d;
         }
       
       }
     }

     return time;

   };




   int add(int index, Neighbors &NN, TPrecision minTime, TPrecision maxTime, int maxD, int
       scale){
     int nAdded = 0;

     NeighborMap &nn = NN[index];
     Simplex s(scale);
     for(NeighborMapIterator it = nn.begin(); it != nn.end(); ++it){
       s.vertices.insert(it->first);
     } 
     s.vertices.erase(index);

     Simplex point(scale);
     point.vertices.insert(index);
     IFiltrationIterator fIt = mappedFiltration.find(point);
     if( fIt == mappedFiltration.end() ){
       tmpFiltration.push_back( FiltrationEntry( point, 0 ) );
       nAdded++; 
     }



     for(int i = 0; i < maxD; i++){
       std::list<Simplex> faces = s.getFaces(i);
       for(std::list<Simplex>::iterator it = faces.begin(); it != faces.end();
           ++it){
         Simplex &f = *it;

         f.vertices.insert(index);
         IFiltrationIterator ifIt = mappedFiltration.find(f);
         if( ifIt == mappedFiltration.end() ){

             TPrecision t = getInsertTime(f, NN);
             if(t <= maxTime){
               if(t < minTime){
                 t = minTime;
               }
               tmpFiltration.push_back( FiltrationEntry( f, t ) );
               nAdded++; 
             }
         }
       }
     }
     return nAdded;

   };



};

#endif 

