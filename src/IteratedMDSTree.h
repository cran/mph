#ifndef ITERATEDMDSTREE_H
#define ITERATEDMDSTREE_H

#include "GMRATree.h"
#include "SVD.h"
#include "RandomSVD.h"

#include <queue>
#include <list>
#include <iostream>



template <typename TPrecision> 
class IMDSNode : public GMRANode{

  public:
    DenseVector<TPrecision> dir;
    TPrecision a; 



    IMDSNode<TPrecision>(){};



    IMDSNode<TPrecision>(DenseMatrix<TPrecision> &Xi, int d, std::vector<int>
        ind):indices(ind){

        DenseMatrix<TPrecision> Xc = createMatrix(Xi);

        //Center matrix m
        //remove row mean
        DenseVector<TPrecision> rowMean = Linalg<TPrecision>::SumRows(m);
        Linalg<TPrecision>::Scale(rowMean, (TPrecision) 1.0/m.M(), rowMean);
        Linalg<TPrecision>::SubtractRowwise(Xc, rowMean, Xc);
        rowMean.deallocate();

        //remove column mean
        DenseVector<TPrecision> colMean = Linalg<TPrecision>::SumColumns(m);
        Linalg<TPrecision>::Scale(colMean, (TPrecision) 1.0/m.N(), colMean);
        Linalg<TPrecision>::SubtractColumnwise(Xc, colMean, Xc);
        colMean.deallocate();
        

        DenseVector<TPrecision> l = Linalg<TPrecision>::ColumnwiseSquaredNorm(Xc);
        r = sqrt( Linalg<TPrecision>::Max(l) );

        DenseMatrix<TPrecision> U;
        DenseVector<TPrecision> S;
        DenseMatrix<TPrecision> Vt;
        /*if(Xc.N() > 1000 && Xc.M() > d+20){
          RandomSVD<TPrecision> svd(Xc, d+10, 1);
          U = svd.U;
          S = svd.S;
          }
          else{*/
        SVD<TPrecision> svd(Xc);
        U = svd.U;
        S = svd.S;
        Vt = svd.Vt;
        //}

        if(U.N() < d){
          phi = DenseMatrix<TPrecision>(U.M(), d);
          Linalg<TPrecision>::Zero(phi);
          Linalg<TPrecision>::SetColumns(phi, 0, U.N(), U, 0);
        }
        else{
          phi = Linalg<TPrecision>::ExtractColumns(U, 0, d);
        }

        sigma = DenseVector<TPrecision>(d);
        Linalg<TPrecision>::Zero(sigma);
        for(int i=0; i < std::min(d, (int)S.N()); i++){
          sigma(i) = S(i);
        }

        sigma2 = DenseVector<TPrecision>( sigma.N() );
        for(int i=0; i<sigma.N(); i++){
          sigma2(i) = sigma(i) * sigma(i);
        } 
        Linalg<TPrecision>::Scale(sigma , 1.0 /  sqrt( Xc.N() ) , sigma );
        Linalg<TPrecision>::Scale(sigma2 , 1.0 /  Xc.N() , sigma2 );


        mse = DenseVector<TPrecision>( d+1 );
        Linalg<Precision>::Zero( mse );
        mse(0) = Linalg<TPrecision>::Sum(l) / Xc.N();

        for(int i = 0; i <  std::min( (int) sigma2.N(), d ); i++){
          mse(i+1) = mse(i) - sigma2(i);
          if(mse(i+1) < 0 ){
            mse(i+1) = 0;
          }
        }

        U.deallocate();
        S.deallocate();
        Vt.deallocate();

        Xc.deallocate();

        dir = Linalg<TPrecision>::ExtractColumn(phi, 0);
        a = Linalg<TPrecision>::Dot(center, dir);
      };




    ~IMDSNode(){
      dir.deallocate();
      phi.deallocate();
      mse.deallocate();
      sigma.deallocate();
      sigma2.deallocate();
    }; 




    std::list< GMRANode * > &getChildren(){
      return children;
    };

   
    std::vector<int> &getPoints(){
      return indices;
    };


    DenseVector<TPrecision> sigma;
    DenseVector<TPrecision> sigma2;
    DenseMatrix<TPrecision> phi;

    DenseVector<TPrecision> mse;
    TPrecision r;

  private:
    std::list< GMRANode* > children;
   
    std::vector<int> indices;

    DenseMatrix<TPrecision> createMatrix(DenseMatrix<TPrecision> &X){
      DenseMatrix<TPrecision> M(X.M(), indices.size());
      
      std::vector<int>::iterator it = indices.begin();
      for(int i=0; i < M.N(); i++, ++it){
        Linalg<TPrecision>::SetColumn(M, i, X, *it);
      }

      return M;
    };

};






template <typename TPrecision>
class IteratedMDSTree : public GMRATree<TPrecision>{

  private:
    unsigned int d;
    TPrecision epsilon;
    int minPoints;


    IMDSNode<TPrecision> *root;


    DenseMatrix<TPrecision> createMatrix(std::vector<int> ind, DenseMatrix<TPrecision> X){
      DenseMatrix<TPrecision> M(X.M(), ind.size());
      int index = 0;
      for(std::vector<int>::iterator it = ind.begin(); it != ind.end(); ++it){
        int i = *it;
        Linalg<TPrecision>::SetColumn(M, index, X, i);
        ++index;
      }
      return M;
    };



    void buildTreeRecursive(IMDSNode<TPrecision> *node, DenseMatrix<TPrecision> X){


      std::vector<int> m1;
      std::vector<int> m2;
      DenseVector<TPrecision> tmp(X.M());
      for(std::vector<int>::iterator it = node->getPoints().begin();
          it!=node->getPoints().end(); ++it){
        int i = *it;
        Linalg<TPrecision>::ExtractColumn(X, i, tmp);
        if(node->split(tmp) > 0){
          m1.push_back( node->getPoints()[i] );
        }
        else{
          m2.push_back( node->getPoints()[i] );
        }
      }
      tmp.deallocate();

      if( node->getPoints().size() < std::max(minPoints, (int)d) || node->mse(d) <= epsilon ){
        return;
      }

      std::cout << m1.size() << " , " << m2.size() << ": " << node->mse(d) << ", "
        << node->r << std::endl;

      if(m1.size() == 0 || m2.size() == 0){
        return;
      }

      IMDSNode<TPrecision> *n1;
      //if(X1.N() > nPoints){
      //  n1 = new IMDSNode<TPrecision>(X1, node);
      //  node->children.push_back(n1);
      //}
      //else{
      if(m1.size() > 0){
        n1 = new IMDSNode<TPrecision>(X, d, m1);
        node->getChildren().push_back(n1);
      }
      // }

      IMDSNode<TPrecision> *n2;
      //if(X2.N() > nPoints){
      //  n2 = new IMDSNode<TPrecision>(X2, node);
      //  node->children.push_back(n2);
      //}
      //else{
      if(m2.size() > 0){
        n2 = new IMDSNode<TPrecision>(X, d, m2);
        node->getChildren().push_back(n2);
      }
      //}

      if(m1.size() >= d){
        buildTreeRecursive(n1, X);
      }

      if(m2.size() >= d){
        buildTreeRecursive(n2, X);
      }

    };





  public:


    IteratedMDSTree(){
    };



    IteratedMDSTree(int dim, TPrecision eps, int minLeafSize=0){
      d = dim;
      epsilon = eps;
      minPoints = minLeafSize;
    };



    //
    void construct(TPrecision *X, int m, int n){

      DenseMatrix<TPrecision> Xd(m, n, X);

      std::vector<int> all;
      for(int i=0; i<n; i++){
        all.push_back(i);
      }
      
      //if(Xc.N() < nPoints){
      root = new IMDSNode<TPrecision>(Xd, d, all);
      //}
      //else{
      //  root = new IMDSNode<TPrecision>(Xc);
      //}
      buildTreeRecursive((IMDSNode<TPrecision>*)root, Xd); 
    };



    void add(TPrecision *x){
      throw "not implemented yet";
    };



    std::list<GMRANode *> getLeafPath(double *x) {
      DenseVector<TPrecision> xv(root->center.N(), x);

      std::list<GMRANode *> path;
      GMRANode *node = root;
      while( !node->getChildren().empty() ){
        path.push_back(node);
        
        IMDSNode<TPrecision> *inode = (IMDSNode<TPrecision> *) node; 
        TPrecision d = inode->split(xv);
        if(d > 0){
          node = node->getChildren().front(); 
        }
        else{
          node = node->getChildren().back(); 
        }  
      }
      return path;
    }; 



    GMRANode *getRoot(){
      return (GMRANode*) root;
    };



    void flatten(std::ofstream &file){
      std::list<IMDSNode<TPrecision> *> nodes;
      nodes.push_back((IMDSNode<TPrecision>*) root);

      file.write( (char*) &epsilon, sizeof(TPrecision) ) ;
      file.write( (char*) &d, sizeof(unsigned int) )  ;   
      unsigned int m = root->phi.M();
      file.write( (char*) &m, sizeof(unsigned int) ) ;   
      file.write( (char*) &minPoints, sizeof(int) );   

      while( !nodes.empty() ){
        IMDSNode<TPrecision> *node =  nodes.front();
        nodes.pop_front();

        file.write((char*)node->phi.data(), node->phi.M()*node->phi.N()*sizeof(TPrecision) );
        file.write((char*)node->sigma.data(), node->sigma.N()*sizeof(TPrecision) );
        file.write((char*)node->mse.data(), node->mse.N()*sizeof(TPrecision) );
        file.write((char*)node->dir.data(), node->dir.N()*sizeof(TPrecision) );
        TPrecision r = node->r;
        file.write((char*) &r, sizeof(TPrecision) );

        for(typename std::list< GMRANode* >::iterator it =
            node->getChildren().begin(); it != node->getChildren().end(); ++it){
          nodes.push_back((IMDSNode<TPrecision>*) *it);
        }

      }
    };





    void unflatten(std::ifstream &file){
      file.read( (char*) &epsilon, sizeof(TPrecision) ) ;
      file.read( (char*) &d, sizeof(unsigned int) )  ;
      unsigned int m;   
      file.read( (char*) &m, sizeof(unsigned int) );   
      file.read( (char*) &minPoints, sizeof(int) );   

      std::list<IMDSNode<TPrecision> *> nodes;

      IMDSNode<TPrecision> *cur = NULL;
      while( !file.eof() ){
        DenseMatrix<TPrecision> phi(m, d);
        file.read((char*)phi.data(), d*m*sizeof(TPrecision));

        DenseVector<TPrecision> sigma(d);
        file.read((char*)sigma.data(), d*sizeof(TPrecision));


        DenseVector<TPrecision> mse(d+1);
        file.read((char*)mse.data(), (d+1)*sizeof(TPrecision));

        DenseVector<TPrecision> dir(m);
        file.read((char*)dir.data(), m*sizeof(TPrecision));

        IMDSNode<TPrecision> *node;
        if(cur == NULL){
          root = new IMDSNode<TPrecision>();
          node = (IMDSNode<TPrecision> *) root;
        }
        else{
          node = new IMDSNode<TPrecision>();
        }

        TPrecision r;
        file.read((char*) &r, sizeof(TPrecision) );
        node->r = r;
        node->phi = phi;
        node->sigma = sigma;
        node->dir = dir;
        node->mse = mse;

        node->sigma2 = DenseVector<TPrecision>(sigma.N());
        for(int i=0; i<sigma.N(); i++){
          node->sigma2(i) = sigma(i) * sigma(i);
        }

        if(cur != NULL){
          nodes.push_back(node);
          if(cur->getChildren().size() == 2){
            cur = nodes.front();
            nodes.pop_front();
          }
          cur->getChildren().push_back(node);
        }
        else{
          cur = node;
        }
      }

    };


};

#endif
