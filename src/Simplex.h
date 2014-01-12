//author: Samuel Gerber

#ifndef SIMPLEX_H
#define SIMPLEX_H


#include <set>
#include <list>




//typedef std::set<vertices> Simplex;


class Simplex{

private:
  //Time stamp for simplex for mutliscale rips
  int scale;
  std::set<int> mapsTo;
  bool mapped;
  //std::vector<Simplex *> faces;
  //int id;

public:

  Simplex(int s = 0, bool m=false) :  scale(s), mapped(m) {

    
  };
  
  std::set<int> vertices;

  /*
  std::set<int> getVertices(){
    std::set<int> vertices;
    if(faces.empty()){
      vertices.insert(id);
    }
    else{
      for(std::vector<Simplex *>::iterator it= faces.begin(); it != faces.end();
          ++it){
      vertices.insert(
    }
    return vertices;
  };
  */


  bool operator == (const Simplex& other) const{
    return this->vertices == other.vertices;
  };
  

  bool operator < (const Simplex& other) const{
    int s1 = this->vertices.size();
    int s2 = other.vertices.size();
    if(s1 !=  s2){
      return s1 < s2;
    }
    
    return this->vertices < other.vertices;
  };



  bool operator > (const Simplex& other) const{

    int s1 = this->vertices.size();
    int s2 = other.vertices.size();
    if(s1 !=  s2){
      return s1 > s2;
    }

    return this->vertices > other.vertices;
  };



  
  std::list<Simplex> getFaces() const{
    std::list<Simplex> faces;
    if(vertices.size() == 1){
      return faces;
    }

    for(std::set<int>::iterator it=vertices.begin(); it != vertices.end();
        ++it){
      Simplex s = *this;
      s.vertices.erase(*it);
      faces.push_back(s);
    }
    return faces;
  };



/*
  std::list<Simplex> getFaces(unsigned int d){
    std::list<Simplex> faces;
    if( d+1 > vertices.size() ){
      return faces;
    }

    if( (d+1) == vertices.size()){
      faces.push_back(*this);
      return faces;
    }


    std::vector< std::set<int>::iterator > its(d+1);
    its[0] = vertices.begin();
    for(unsigned int i=1; i<=d; i++){
      its[i] = its[i-1];
      ++its[i];
    }

    while(its[d] != vertices.end() ){
      Simplex s(scale);
      for(unsigned int i=0; i <= d; i++){
        s.vertices.insert( *its[i] );
      }
      faces.push_back(s);
      increment(its, d);
    }
    
    return faces;

  };
*/

  int getScale() const{
    return scale;
  };


  void setScale(int s){
    if(s > scale){
#ifdef VERBOSE
      std::cout << "increased scale" << std::endl;
#endif 
    }
    scale = s;
  };

  void setMapped(bool m){
    mapped = m;
  };

  bool isMapped() const{
    return mapped;
  };

  void setMapsTo(std::set<int> &vertices){
    mapsTo = vertices;
  };

  std::set<int> getMapsTo() const{
    return mapsTo;    
  };


  void print() const{
    
#ifdef VERBOSE
    for(std::set<int>::iterator it = vertices.begin(); it !=
           vertices.end(); ++it){
      std::cout << *it << ",";
    }
    std::cout << std::endl;
#endif
  };



private:

  bool increment( std::vector< std::set<int>::iterator > &its, int i){
    if(i < 0){
      return false;
    }
    ++its[i];
    
    if(its[i] == vertices.end()){
      increment(its, i-1);
    }
    else{
      for(unsigned int j=i+1; j < its.size(); j++){
        its[j] = its[j-1];
        ++its[j];
        if(its[j] == vertices.end()){
          increment(its, i-1);
        }
      }
    }
    return true;
  };


};

/*
struct SimplexCompare : public std::binary_function<const char*, const char*, bool> {
public:
    bool operator() (const Simplex* s1, const Simplex* s2) const{ 
      return s1->vertices < s2->vertices; 
    };
};
*/


#endif 

