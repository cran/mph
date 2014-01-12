#ifndef LPDISTANCE_H
#define LPDISTANCE_H

#include "DenseVector.h"
#include "DenseMatrix.h"

#include "GWT.h"
#include "GMRATree.h"
#include "Wasserstein.h"

#include <glpk.h>

#include <time.h>

template < typename TPrecision >
class NodeDistance{

  public:

    virtual TPrecision distance(GMRANode *n1, GMRANode *n2) = 0;

};



template < typename TPrecision >
class WassersteinNodeDistance : public NodeDistance<TPrecision> {

  public:

    TPrecision distance(GMRANode *n1, GMRANode *n2){
      GWTNode<TPrecision> *gwt1 = (GWTNode<TPrecision> *) n1->getGWTInfo();
      GWTNode<TPrecision> *gwt2 = (GWTNode<TPrecision> *) n2->getGWTInfo();

      DenseVector<TPrecision> sigma12 = Linalg<TPrecision>::Copy(gwt1->sigma);
      DenseVector<TPrecision> sigma22 = Linalg<TPrecision>::Copy(gwt2->sigma);
      for(int i=0; i< sigma12.N(); i++){
        sigma12(i) *= sigma12(i);
      }
      for(int i=0; i< sigma22.N(); i++){
        sigma22(i) *= sigma22(i);
      }

      TPrecision d = Wasserstein<TPrecision>::distance2( gwt1->phi, sigma12, gwt1->center, 
                                                         gwt2->phi, sigma22, gwt2->center );
      sigma12.deallocate();
      sigma22.deallocate();
      return d;
    };

};


class LPSolution{
  
  public:
    double objective;  //objective function value
    int optimizationStatus;
    DenseMatrix<double> W; //solution weights of costs
    DenseMatrix<double> C; //costs
    DenseVector<int> pts1;		
    DenseVector<int> pts2;
    unsigned long time;
    double errorBound;

    void deallocate(){
      W.deallocate();
      C.deallocate();
      pts1.deallocate();
      pts2.deallocate();
    };  

};

template <typename TPrecision>
class LPDistance{


  public:

    //Compares the distance between the densities described by two gmra trees. The
    //metric is how to measure distance between the base distirbution between two
    //nodes in the gmra tree. Depending on he metric the trees need to have
    //specific GWTInfo attached,e.g. for gaussian base measures the GWTInfo needs
    //to be a GWTNode with the low d gaussian measure info (phi,sigma and
    //centers)
    static std::list<LPSolution> distance(GMRATree<TPrecision> *tree1,
        GMRATree<TPrecision> *tree2, NodeDistance<TPrecision> *metric){

      //compute distances at each tree level
      std::list<LPSolution> solutions; 

      //data strcutures to keep track of current nodes to be matched ...
      std::list<GMRANode *> nodes1;
      std::list<GMRANode *> nodes2;

      //... and their corresponding probability weights
      std::list< double > w1;
      std::list< double > w2;


      GMRANode* root1 = tree1->getRoot();
      GMRANode* root2 = tree2->getRoot();
      
      int nPoints1 = root1->getPoints().size();
      int nPoints2 = root2->getPoints().size();

      nodes1.push_back( root1 );
      nodes2.push_back( root2 );

      bool new1 = true;
      bool new2 = true;
      while( new1 || new2 ){
        //Solve linear program
        clock_t start = clock();
        LPSolution sol = match(nodes1, nodes2, metric, nPoints1, nPoints2);
        clock_t end = clock();
        sol.time = (end - start) * 1000 / CLOCKS_PER_SEC;
        solutions.push_back(sol);

        //Compute weights and nodes involved in next level
        std::list<GMRANode *> nn1;
        new1 = nextLevel(nodes1, nn1);
        nodes1 = nn1;

        std::list<GMRANode *> nn2;
        new2 = nextLevel(nodes2, nn2);
        nodes2 = nn2;

      }

      return solutions;
    };




  private:


    static LPSolution match( std::list<GMRANode *> &n1, std::list<GMRANode *> &n2,  
                  NodeDistance<TPrecision> *metric, int nPoints1, int nPoints2 ){
  
        
       std::cout << "Solving LP: " << n1.size() << " x " << n2.size() << std::endl; 
        DenseMatrix<TPrecision> D(n1.size(), n2.size());
        int i1 = 0;
        for( std::list<GMRANode *>::iterator it1 = n1.begin(); it1 != n1.end(); ++it1, ++i1){
          int i2 = 0;
          for( std::list<GMRANode *>::iterator it2 = n2.begin(); it2 != n2.end(); ++it2, ++i2){
             D(i1, i2) = metric->distance(*it1, *it2);
          }
        }

        //setup linear programm
        glp_prob *lp = glp_create_prob();
        glp_set_obj_dir(lp, GLP_MIN);
        
        //rows - contraints:
        // M of type sum of weights going from node A tree1 to nodes in tree2 =
        // weight of node A in tree1
        // N of type sum of weights from nodes in tree1 going to node A in tree2
        // = weight of node A in tree2
        glp_add_rows(lp, D.M() + D.N() );
        //columns - coefficents for each combination of nodes from tree1 to
        //tree2 = M*N
        glp_add_cols(lp, D.M() * D.N() );

        //set coefficients equal to cost of moving mass from node i in tree1 to
        //node j in tree2
        for(int i=0; i<D.M(); i++){
          for(int j=0; j<D.N(); j++){
            int col = 1+i*D.N() + j;
            glp_set_col_bnds(lp, col, GLP_DB, 0, 1);
            glp_set_obj_coef(lp, col, D(i, j) ); 
          }
        }

        //set constraints bounds
        int row=1;
        //The weights of each nodes in tree1
        for(std::list<GMRANode *>::iterator it = n1.begin(); it !=n1.end(); ++it, ++row){
          double w = (*it)->getPoints().size()/((double)nPoints1);
          glp_set_row_bnds(lp, row, GLP_FX, w, w);
        }
        //the weight of each node in tree2
        for(std::list<GMRANode *>::iterator it = n2.begin(); it !=n2.end(); ++it, ++row){
          double w = (*it)->getPoints().size()/((double)nPoints2);
          glp_set_row_bnds(lp, row, GLP_FX, w, w);
        }

        //set constraint coefficents for the M nodes giving weight from tree1 to
        //tree2:
        //the sum of the N weights leaving from i in tree1 has to equal weight
        //i in tree2 
        for(int i=0; i<D.M(); i++){
          int ind[D.N()+1];
          double val[D.N()+1];
          for(int j=0; j<D.N(); j++){
            ind[j+1] = 1 + i*D.N() + j;
            val[j+1] = 1;
          }
          glp_set_mat_row(lp, 1+i, D.N(), ind, val);
        }        

        //set constraint coefficents for the N nodes of tree2 receving weight
        //from tree1:
        //the sum of the M weights ending up in node i of tree2 has to equal the
        //weight of node i in tree2
        for(int i=0; i<D.N(); i++){
          int ind[D.M()+1];
          double val[D.M()+1];
          for(int j=0; j<D.M(); j++){
            ind[j+1] = 1 + j*D.N() + i;
            val[j+1] = 1;
          }
          glp_set_mat_row(lp, 1+i + D.M(), D.M(), ind, val);
        }
   
        //solve lp with simplex algorithm
        int ret = glp_simplex(lp, NULL);
        //int ret = glp_interior(lp, NULL);

        LPSolution sol;
        sol.optimizationStatus = ret;
        sol.objective = glp_get_obj_val(lp);
        
        sol.W = DenseMatrix<double>( D.M(), D.N() );        

        int rIndex =0;
        sol.errorBound = 0;
        for(std::list<GMRANode *>::iterator it1 = n1.begin(); it1 !=n1.end(); ++it1, ++rIndex){
          int cIndex = 0;
          for(std::list<GMRANode *>::iterator it2 = n2.begin(); it2 !=n2.end(); ++it2, ++cIndex){
             int col = 1+rIndex*D.N() + cIndex;
             sol.W(rIndex, cIndex) = glp_get_col_prim(lp, col);
             //sol.W(rIndex, cIndex) = glp_ipt_col_prim(lp, col);
             sol.errorBound += sol.W(rIndex, cIndex) * ( (*it1)->r + (*it2)->r );
          }
        }
         
        sol.C = D;
        
        int index = 0;
        sol.pts1 = DenseVector<int>(nPoints1);
        for( std::list<GMRANode *>::iterator it = n1.begin(); it != n1.end(); ++it, ++index){
           std::vector<int> pts = (*it)->getPoints();
           for(int i=0; i<pts.size(); i++){
             sol.pts1(pts[i]) = index;
           }
        }

        index = 0;
        sol.pts2 = DenseVector<int>(nPoints2);
        for( std::list<GMRANode *>::iterator it = n2.begin(); it != n2.end(); ++it, ++index){
           std::vector<int> pts = (*it)->getPoints();
           for(int i=0; i<pts.size(); i++){
             sol.pts2(pts[i]) = index;
           }
        }

        glp_delete_prob(lp);

        return sol;
    };



    static bool nextLevel(std::list< GMRANode * > &nodes, std::list<GMRANode *> &nn){
      bool newNodes = false;
      while( !nodes.empty() ){
        GMRANode *n = nodes.front();
        nodes.pop_front();

        std::vector< GMRANode * > children = n->getChildren();
        if(children.empty() ){
          nn.push_back(n);
        }
        else{
          newNodes = true;
          double nPts = n->getPoints().size();
          for(std::vector< GMRANode * >::iterator it = children.begin(); it !=
              children.end(); ++it){
            GMRANode *ch = (*it);
            nn.push_back(ch);
          }
        }
      }
      return newNodes;
    };

};




#endif
