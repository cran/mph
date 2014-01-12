//author: Samuel Gerber

#ifndef MAPPINGPERSISTENTHOMOLOGY_H
#define MAPPINGPERSISTENTHOMOLOGY_H


#include <map>
#include <set>
#include <vector>
#include <utility>
#include <limits>
#include <list>
#include <stdlib.h>

#include "Simplex.h"
#include "Linalg.h"
#include "MappingRips2.h"


#include "DenseMatrix.h"



template <typename TPrecision>
class MappingPersistentHomology{

  public:

    typedef typename MappingRips2<TPrecision>::FiltrationEntry FiltrationEntry; 
    typedef typename MappingRips2<TPrecision>::Filtration Filtration;

    typedef typename std::map<Simplex, TPrecision> IFiltration;

    typedef typename IFiltration::iterator IFiltrationIterator;
    typedef typename Filtration::iterator FiltrationIterator;

    class Event{
      public:
        Event(const TPrecision &b, const TPrecision &d, int dim) : birth(b),
        death(d), dimension(dim){};

        TPrecision birth;
        TPrecision death;
        int dimension;
    };

    typedef std::vector< Event > Events;
    typedef typename Events::iterator EventsIterator;


    typedef std::map<Simplex, int > IDMap;
    typedef typename IDMap::iterator IDMapIterator;

    class Column{
      public:
        Column(){
          time = 0;
        }
        std::list<int> entries;
        TPrecision time;
        int dim;
        //Simplex simplex;
    };



  private:

    IDMap ids;
    IDMap mids;
    TPrecision maximalEntranceTime;

    Events events;


    void getColumn(const Simplex &s, std::list<int> &col, int filtSize){

      std::list<Simplex> faces = s.getFaces();
      for(std::list<Simplex>::iterator it = faces.begin(); it != faces.end();
          ++it){
        if(ids[*it] == 0){
#ifdef VERBOSE
          std::cout << "moo" << std::endl;
#endif
        }
        col.push_back( getID(*it) );
      }

      //Add mapped simplex id at top of column with negative id
      Simplex mapped;
      mapped.vertices = s.getMapsTo();
      if( mapped.vertices.size() == s.vertices.size() ){
        col.push_back( getMappedID(mapped) - filtSize - 1 );
      }
      
      
      col.sort();
     // std::sort(col.begin(), col.end());
    };




    int getLow(std::list<int> &col){
      if( col.empty() ){
        return 0;
      }
      return col.back();
    };




    int getID(const Simplex &s){
      IDMapIterator it = ids.find(s);
      if(it == ids.end()){
        int id = ids.size() + 1;
        ids[s] = id;
        return id;
      }
      return it->second;
    };




    int getMappedID(const Simplex &s){
      IDMapIterator it = mids.find(s);
      if(it == mids.end()){
        int id = mids.size() + 1;
        mids[s] = id;
        return id;
      }
      return it->second;
    };





  public:

    MappingPersistentHomology(){ 
    };



    //The filtration should not contain any vertices but start at the edges
    //This methods expects and "inverted filtration", i.e. each simplex mapping
    //to a time, the sorting is done as a preprocessing step here
    // (Note: this conversion might not be necessary since the simplicies in the
    // set are in lexicographical order)
    void run(IFiltration &ifilt, TPrecision mapTime, int maxDim){

#ifdef VERBOSE
      std::cout << "IFilt size: " << ifilt.size() << std::endl;
#endif
      Filtration filt;
      for(IFiltrationIterator it = ifilt.begin(); it != ifilt.end();
          ++it){
        filt.push_back( FiltrationEntry(it->first, it->second) );
      }
      filt.sort();

      run(filt, mapTime, maxDim);
    };




    //The filtration should not contain any vertices but start at the edges
    void run(Filtration &filt, TPrecision mapTime, int maxDim){

      //clear any possible leftovers
      ids.clear();
      mids.clear();
      events.clear();

      int filtSize = filt.size();
      //columns, id to a set of ids
      std::vector< Column > columns( filtSize + 1 );
      std::vector< int > low2id( filtSize + 1 );
      std::vector< int > mLow2id( filtSize + 1 );
      
      
      std::set< int > zeroes;

      //Do matrix reduction
      TPrecision currentTime = 0;
      //FortranLinalg::DenseVector<int> delta(filtSize);
      //int index = 0;
      for(FiltrationIterator it = filt.begin(); it != filt.end(); ++it){

        const Simplex &current = it->simplex;
        int id = getID(current);
        Column &c1 = columns[id];
        //c1.simplex = current;
        c1.time = it->time;
        c1.dim = current.vertices.size();
        currentTime = c1.time;


        //Do reduction with column c1
        //clock_t t1 =  clock();
        getColumn(current, c1.entries, filtSize);
        int low = getLow( c1.entries );
        bool collision = low != 0;
        while( collision ){
         
          int lowInd = 0;
          if(low < 0){
            lowInd = mLow2id[-low];
          }
          else{
            lowInd = low2id[low];
          } 
          Column &c2 = columns[ lowInd ];
          

          /*
          std::list<int> diff;
          std::set_symmetric_difference(c1.entries.begin(), c1.entries.end(),
             c2.entries.begin(), c2.entries.end(), std::inserter(diff, diff.end()) );
          c1.entries = diff;
          */
          symmetric_diff_inplace(c1.entries, c2.entries);

          low = getLow(c1.entries);
          collision = ( !c2.entries.empty() && low != 0 );
        }
        //clock_t t2 =  clock();
        //delta(index) = (t2-t1);



        //Check possible cases of deaths
        if(current.vertices.size() > 1){
          if( low > 0 ){

            //remove from death in mapping
            zeroes.erase(low);
            
            low2id[ low ] = id;
            Column &cLow = columns[low];
            if( cLow.time < it->time  && !current.isMapped() ){//!maps[id] ){
              events.push_back( Event( cLow.time, it->time,
                    current.vertices.size() ) );
            }
          }
          else if(low == 0){
            //potentially dies in mapping
            Column &c1 = columns[id];
            if(c1.time < mapTime && current.vertices.size() < maxDim+1){
              zeroes.insert(id);
            }
          }
          else{
            //it will be either transfered or die within this scale
            mLow2id[-low] = id;
          }
        }
      
      };


      maximalEntranceTime = currentTime;
      
      //Detect deaths in mapping
      for(std::set<int>::iterator it = zeroes.begin(); it != zeroes.end();
          ++it){
        int id = *it;
        Column &c1 = columns[id];
        if(c1.time < mapTime ){
          int dim = c1.dim + 1;
          events.push_back( Event( c1.time, mapTime, dim  ) );
        }
      }


      //FortranLinalg::LinalgIO<int>::writeVector("times.data", delta);
      //std::cout << "#zeros: " << zeroes.size() << std::endl;
#if VERBOSE
      std::cout << "#Events: " << events.size() << std::endl;
#endif
    };



    Events &getEvents(){
      return events; 
    };



    FortranLinalg::DenseMatrix<TPrecision> getDiagram(){
      FortranLinalg::DenseMatrix<TPrecision> ph(3, events.size());
      int index = 0;
      for(typename Events::iterator it = events.begin(); it != events.end();
          ++it, ++index){
        Event &e = (*it);
        //ph(0, index) =  e.s.scale;
        ph(0, index) =   e.birth;
        ph(1, index) =   e.death;
        ph(2, index) =   e.dimension;
      }
      return ph;
    };

    TPrecision getMaximalTime(){
      return maximalEntranceTime;
    };


  private:

    void symmetric_diff_inplace(std::list<int> &c1, std::list<int> &c2){
      std::list<int>::iterator it1 = c1.begin();
      std::list<int>::iterator it2 = c2.begin();

      while (true)
      {
        if( it1 == c1.end() ){
          c1.insert(it1, it2, c2.end() );
          return;
        }
        if( it2 == c2.end() ){
          return;
        }

        if( *it1 < *it2 ) { 
          ++it1; 
        }
        else if( *it2 < *it1 ) { 
          c1.insert(it1, *it2);
          ++it2;
        }
        else { 
          it1 = c1.erase(it1); 
          ++it2; 
        }
      }
    };



};

#endif 

