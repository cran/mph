#ifndef HIERARCHICALLPDISTANCE_H
#define HIERARCHICALLPDISTANCE_H

#include "DenseVector.h"
#include "DenseMatrix.h"

#include "GWT.h"
#include "GMRATree.h"
#include "Wasserstein.h"

#include "LPDistance.h"


#include <set>
#include <map>
#include <list>
#include <vector>


template <typename TPrecision>
class HierarchicalLPDistance{


  public:

    //Compares the distance between the densities described by two gmra trees. The
    //metric is how to measure distance between the base distirbution between two
    //nodes in the gmra tree. Depending on he metric the trees need to have
    //specific GWTInfo attached,e.g. for gaussian base measures the GWTInfo needs
    //to be a GWTNode with the low d gaussian measure info (phi,sigma and
    //centers)
    static std::list<LPSolution> distance(GMRATree<TPrecision> *tree1,
        GMRATree<TPrecision> *tree2, NodeDistance<TPrecision> *metric, int maxVar = 0){

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
      
      nodes1.push_back( root1 );
      nodes2.push_back( root2 );
      
      std::vector< int > parent1;
      std::vector< int > parent2;
      
      parent1.push_back( 0 );
      parent2.push_back( 0 );

      std::vector<int> prevParent1 = parent1;
      std::vector<int> prevParent2 = parent2;
      
      int nPoints1 = root1->getPoints().size();
      int nPoints2 = root2->getPoints().size();
      
      DenseMatrix<TPrecision> WInit(1,1);
      WInit(0,0) = 1;

      DenseMatrix<TPrecision> W = WInit;
      DenseMatrix<TPrecision> prevW =W;

      bool new1 = true;
      bool new2 = true;
      LPSolution *prevSol = NULL;
      while( new1 || new2 ){
        //Solve linear program       
        clock_t start = clock();
        LPSolution sol = match(nodes1, nodes2, metric, nPoints1, nPoints2,
                   parent1, parent2, W, prevParent1, prevParent2, prevW, maxVar);
        clock_t end = clock();
        sol.time = (end - start) * 1000 / CLOCKS_PER_SEC;
        solutions.push_back(sol);

        prevParent1 = parent1;
        prevParent2 = parent2;
        prevW = W;
        W = sol.W;

        //Compute weights and nodes involved in next level
        std::list<GMRANode *> nn1;
        new1 = nextLevel(nodes1, nn1, parent1);
        nodes1 = nn1;

        std::list<GMRANode *> nn2;
        new2 = nextLevel(nodes2, nn2, parent2);
        nodes2 = nn2;

      }

      WInit.deallocate();

      return solutions;
    };




  private:


    static LPSolution match( std::list<GMRANode *> &n1,std::list<GMRANode *> &n2,  
                  NodeDistance<TPrecision> *metric, int nPoints1, int nPoints2,
                  std::vector<int> parent1, std::vector<int> parent2,
                  DenseMatrix<TPrecision> pW, std::vector<int> pparent1,
                  std::vector<int> pparent2,
                  DenseMatrix<TPrecision> ppW, int maxVar ){
  
        
        DenseMatrix<TPrecision> D(n1.size(), n2.size());
        int nVar = n1.size()*n2.size();
        std::map< int, std::list<int> > var1;
        std::map< int, std::list<int> > var2;
        int nCombinations = 0;
        int i1 = 0;
        for( std::list<GMRANode *>::iterator it1 = n1.begin(); it1 != n1.end(); ++it1, ++i1){
          int i2 = 0;
          for( std::list<GMRANode *>::iterator it2 = n2.begin(); it2 != n2.end(); ++it2, ++i2){
             //if(prevW(parent1[i1], parent2[i2]) != 0 ){
             if(maxVar > nVar || ppW(pparent1[parent1[i1]], pparent2[parent2[i2]]) != 0 ){
               var1[i1].push_back(i2);
               nCombinations++;
               var2[i2].push_back(i1);
               D(i1, i2) = metric->distance(*it1, *it2);
             }
          }
        }
        std::cout << "Solving LP for tree sizes: " << n1.size() << " x " << n2.size() << std::endl; 
        std::cout << "Reduced to " << nCombinations << " variables" << std::endl; 

        //setup linear programm
        glp_prob *lp = glp_create_prob();
        glp_set_obj_dir(lp, GLP_MIN);
        
        //rows - constraints:
        // M of type sum of weights going from node A tree1 to nodes in tree2 =
        // weight of node A in tree1
        // N of type sum of weights from nodes in tree1 going to node A in tree2
        // = weight of node A in tree2
        glp_add_rows(lp, var1.size() + var2.size() );
        std::cout << "Sanity check: " << var1.size() << " x " << var2.size() << std::endl; 
        //columns - coefficents for each combination of nodes from tree1 to
        //tree2 = M*N
        glp_add_cols(lp, nCombinations );

        //set coefficients equal to cost of moving mass from node i in tree1 to
        //node j in tree2
        int index = 1;
        for(std::map< int, std::list<int> >::iterator it = var1.begin(); it !=
            var1.end(); ++it){
          int from = it->first;
          std::list<int> &to = it->second;
          for(std::list<int>::iterator toIt = to.begin(); toIt!=to.end(); ++toIt){
            glp_set_col_bnds(lp, index, GLP_DB, 0, 1);
            glp_set_obj_coef(lp, index, D(from, *toIt) ); 
            ++index;
          }
        }

        //set constraints bounds
        int row=1;
        //The weights of each nodes in tree1
        for(std::list<GMRANode *>::iterator it = n1.begin(); it !=n1.end(); ++it, ++row){
          double w = (*it)->getPoints().size()/((double) nPoints1);
          glp_set_row_bnds(lp, row, GLP_FX, w, w);
        }
        //the weight of each node in tree2
        for(std::list<GMRANode *>::iterator it = n2.begin(); it !=n2.end(); ++it, ++row){
          double w = (*it)->getPoints().size()/((double) nPoints2);
          glp_set_row_bnds(lp, row, GLP_FX, w, w);
        }


        //set constraint coefficents for the M nodes giving weight from tree1 to
        //tree2:
        //the sum of the N weights leaving from i in tree1 has to equal weight
        //i in tree2 
        index = 1;
        int colIndex = 1;
        std::map< std::pair<int, int>, int> varmap;
        for(std::map< int, std::list<int> >::iterator it = var1.begin(); it !=
            var1.end(); ++it){
          int from = it->first;
          std::list<int> &to = it->second;
          int ind[to.size()+1];
          double val[to.size()+1];
          int j=1;
          for(std::list<int>::iterator toIt = to.begin(); toIt!=to.end(); ++toIt,
              ++j){
            ind[j] = colIndex;
            val[j] = 1;
            varmap[ std::pair<int, int>(from, *toIt) ] = colIndex;
            colIndex++; 
          }
          glp_set_mat_row(lp, index, to.size(), ind, val);
          index++;
        }

        //set constraint coefficents for the N nodes of tree2 receving weight
        //from tree1:
        //the sum of the M weights ending up in node i of tree2 has to equal the
        //weight of node i in tree2
        for(std::map< int, std::list<int> >::iterator it = var2.begin(); it !=
            var2.end(); ++it){
          int to = it->first;
          std::list<int> &from = it->second;
          int ind[from.size()+1];
          double val[from.size()+1];
          int j=1;
          for(std::list<int>::iterator fromIt = from.begin(); fromIt!=from.end();
              ++fromIt, ++j){
            val[j] = 1;
            ind[j]= varmap[ std::pair<int, int>(*fromIt, to) ];
          }
          glp_set_mat_row(lp, index, from.size(), ind, val);
          index++;
        }

        //solve lp with simplex algorithm
        //int ret = glp_simplex(lp, NULL);
        int ret = glp_interior(lp, NULL);

        LPSolution sol;
        sol.optimizationStatus = ret;
        sol.objective = glp_get_obj_val(lp);
        
        sol.W = DenseMatrix<double>( D.M(), D.N() );
        Linalg<TPrecision>::Zero(sol.W);
        index = 1;
        for(std::map< int, std::list<int> >::iterator it = var1.begin(); it !=
            var1.end(); ++it){
          int from = it->first;
          std::list<int> &to = it->second;
          for(std::list<int>::iterator toIt = to.begin(); toIt!=to.end(); ++toIt){
             //sol.W(from, *toIt) = glp_get_col_prim(lp, index);
             sol.W(from, *toIt) = glp_ipt_col_prim(lp, index);
             index++;
          }
        }
         
        sol.C = D;
        
        index = 0;
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



    static bool nextLevel(std::list< GMRANode * >
        &nodes, std::list<GMRANode *> &nn, std::vector<int> &parents){
      bool newNodes = false;
      int pID = 0;
      nn.clear();
      parents.clear();
      while( !nodes.empty() ){

        GMRANode *n = nodes.front();
        nodes.pop_front();

        std::vector< GMRANode * > children = n->getChildren();
        if(children.empty() ){
          nn.push_back(n);
          parents.push_back(pID);
        }
        else{
          newNodes = true;
          double nPts = n->getPoints().size();
          for(std::vector< GMRANode * >::iterator it = children.begin(); it !=
              children.end(); ++it){
            GMRANode *ch = (*it);
            nn.push_back(ch);
            parents.push_back(pID);
          }
        }
        pID++;
      }
      return newNodes;
    };

};




#endif
