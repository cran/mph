#ifndef IPCATREE_H
#define IPCATREE_H

#include "GMRATree.h"
#include "SVD.h"
#include "RandomSVD.h"

#include <queue>
#include <map>
#include <list>
#include <iostream>
#include <fstream>




//Node factory decides on the dimensionality of each node
template <typename TPrecision>
class NodeFactory{
  protected:
    FortranLinalg::DenseVector<TPrecision> sigma;
    FortranLinalg::DenseVector<TPrecision> sigma2;
    FortranLinalg::DenseVector<TPrecision> mse;
    FortranLinalg::DenseMatrix<TPrecision> phi;
    TPrecision radius;

  public:

    virtual void createNode(FortranLinalg::DenseMatrix<TPrecision> X) = 0;

    FortranLinalg::DenseMatrix<TPrecision> getPhi(){
      return phi;
    };

    FortranLinalg::DenseVector<TPrecision> getSigma(){
      return sigma;
    };

    FortranLinalg::DenseVector<TPrecision> getSigma2(){
      return sigma2;
    };

    FortranLinalg::DenseVector<TPrecision> getMSE(){
      return mse;
    };

    TPrecision getRadius(){
      return radius;
    };

};


template <typename TPrecision>
class FixedNodeFactory : public NodeFactory<TPrecision>{
  private:
    int d;
  public:
    FixedNodeFactory(int dim = 5){
      d = std::max( 1, dim);
    };

    virtual void createNode(FortranLinalg::DenseMatrix<TPrecision> X){
      using namespace FortranLinalg; 
      DenseMatrix<TPrecision> U;
      DenseVector<TPrecision> S;
      DenseMatrix<TPrecision> Vt;
      if(X.N() == 1){
        U = DenseMatrix<TPrecision>(X.M(), 1);
        Linalg<TPrecision>::Zero(U);
        S = DenseVector<TPrecision>(1);
        S(0) = 0;
      }
      else if(X.N() > 1000 && X.M() > d+4 ){
        RandomSVD<TPrecision> svd(X, d+4, 1);
        U = svd.U;
        S = svd.S;
      }
      else{
        SVD<TPrecision> svd(X);
        U = svd.U;
        S = svd.S;
        Vt = svd.Vt;
      }

      if(U.N() < d){
        this->phi = Linalg<TPrecision>::Copy(U);
      }
      else{
        this->phi = Linalg<TPrecision>::ExtractColumns(U, 0, d);
      }

      this->sigma = DenseVector<TPrecision>( this->phi.N() );
      Linalg<TPrecision>::Zero(this->sigma);
      for(unsigned int i=0; i < this->sigma.N(); i++){
        this->sigma(i) = S(i);
      }
      Linalg<TPrecision>::Scale(this->sigma , 1.0 /  sqrt( (TPrecision) X.N() ) , this->sigma );


      this->sigma2 = DenseVector<TPrecision>( this->phi.N() );
      for(unsigned int i=0; i<this->sigma.N(); i++){
        this->sigma2(i) = this->sigma(i) * this->sigma(i);
      } 


      DenseVector<TPrecision> l = Linalg<TPrecision>::ColumnwiseSquaredNorm(X);
      this->radius = sqrt( Linalg<TPrecision>::Max(l) );

      this->mse = DenseVector<TPrecision>( this->phi.N() + 1 );
      Linalg<TPrecision>::Zero( this->mse );
      this->mse(0) = Linalg<TPrecision>::Sum(l) / X.N();

      for(unsigned int i = 0; i <  this->sigma2.N(); i++){
        this->mse(i+1) = this->mse(i) - this->sigma2(i);
        if(this->mse(i+1) < 0 ){
          this->mse(i+1) = 0;
        }
      }

      l.deallocate(); 
      U.deallocate();
      S.deallocate();
      Vt.deallocate();

    };

};





template <typename TPrecision>
class RelativePrecisionNodeFactory : public NodeFactory<TPrecision>{
  private:
    int maxD;
    TPrecision t;
  public:
    RelativePrecisionNodeFactory(int maxDim, TPrecision rel) : t(rel){
      maxD = std::max(1, maxDim);
    };

    virtual void createNode(FortranLinalg::DenseMatrix<TPrecision> X){
      using namespace FortranLinalg; 

      DenseMatrix<TPrecision> U;
      DenseVector<TPrecision> S;
      DenseMatrix<TPrecision> Vt;
      if(X.N() == 1){
        U = DenseMatrix<TPrecision>(X.M(), 1);
        Linalg<TPrecision>::Zero(U);
        S = DenseVector<TPrecision>(1);
        S(0) = 0;
      }
      else if(X.N() > 1000 && X.M() > maxD+4 ){
        RandomSVD<TPrecision> svd(X, maxD+4, 1);
        U = svd.U;
        S = svd.S;
      }
      else{
        SVD<TPrecision> svd(X);
        U = svd.U;
        S = svd.S;
        Vt = svd.Vt;
      }


      DenseVector<TPrecision> l = Linalg<TPrecision>::ColumnwiseSquaredNorm(X);
      this->radius = sqrt( Linalg<TPrecision>::Max(l) );
      
      int d = 1;
      TPrecision mse0 = Linalg<TPrecision>::Sum(l);
      TPrecision mseTmp = 0;
      for(int i=0; i<S.N(); i++){
        mseTmp += S(i)*S(i);
        if(mseTmp/mse0 > t){
          break;
        }
        d++;
        if(d==maxD){
          break;
        }
      }

      if(U.N() < d){
        this->phi = Linalg<TPrecision>::Copy(U);
      }
      else{
        this->phi = Linalg<TPrecision>::ExtractColumns(U, 0, d);
      }

      this->sigma = DenseVector<TPrecision>( this->phi.N() );
      Linalg<TPrecision>::Zero(this->sigma);
      for(unsigned int i=0; i < this->sigma.N(); i++){
        this->sigma(i) = S(i);
      }
      Linalg<TPrecision>::Scale(this->sigma , 1.0 /  sqrt( X.N() ) , this->sigma );


      this->sigma2 = DenseVector<TPrecision>( this->phi.N() );
      for(unsigned int i=0; i<this->sigma.N(); i++){
        this->sigma2(i) = this->sigma(i) * this->sigma(i);
      } 


      this->mse = DenseVector<TPrecision>( this->phi.N()+1 );
      Linalg<TPrecision>::Zero( this->mse );
      this->mse(0) = mse0/X.N();

      for(unsigned int i = 0; i <  this->sigma2.N(); i++){
        this->mse(i+1) = this->mse(i) - this->sigma2(i);
        if(this->mse(i+1) < 0 ){
          this->mse(i+1) = 0;
        }
      }

      l.deallocate(); 
      U.deallocate();
      S.deallocate();
      Vt.deallocate();



    };

};


template <typename TPrecision>
class RelativeRatioNodeFactory : public NodeFactory<TPrecision>{
  private:
    int maxD;
    TPrecision t;
  public:
    RelativeRatioNodeFactory(int maxDim, TPrecision ratio) : t(ratio){
      maxD= std::max(1, maxDim);
    };

    virtual void createNode(FortranLinalg::DenseMatrix<TPrecision> X){
      using namespace FortranLinalg; 

      DenseMatrix<TPrecision> U;
      DenseVector<TPrecision> S;
      DenseMatrix<TPrecision> Vt;
      if(X.N() == 1){
        U = DenseMatrix<TPrecision>(X.M(), 1);
        Linalg<TPrecision>::Zero(U);
        S = DenseVector<TPrecision>(1);
        S(0) = 0;
      }
      else if(X.N() > 1000 && X.M() > maxD+4 ){
        RandomSVD<TPrecision> svd(X, maxD+4, 1);
        U = svd.U;
        S = svd.S;
      }
      else{
        SVD<TPrecision> svd(X);
        U = svd.U;
        S = svd.S;
        Vt = svd.Vt;
      }

      DenseVector<TPrecision> l = Linalg<TPrecision>::ColumnwiseSquaredNorm(X);
      this->radius = sqrt( Linalg<TPrecision>::Max(l) );

      int d = 1;
      for(int i=1; i<S.N(); i++){
        double mse1 = S(i)*S(i);
        double mse0 = S(i-1)*S(i-1);
        if(mse0/mse1 < t){
          break;
        }
        d++;
        if(d==maxD){
          break;
        }
      }

      if(U.N() < d){
        this->phi = Linalg<TPrecision>::Copy(U);
      }
      else{
        this->phi = Linalg<TPrecision>::ExtractColumns(U, 0, d);
      }

      this->sigma = DenseVector<TPrecision>( this->phi.N() );
      Linalg<TPrecision>::Zero(this->sigma);
      for(unsigned int i=0; i < this->sigma.N(); i++){
        this->sigma(i) = S(i);
      }
      Linalg<TPrecision>::Scale(this->sigma , 1.0 /  sqrt( X.N() ) , this->sigma );


      this->sigma2 = DenseVector<TPrecision>( this->phi.N() );
      for(unsigned int i=0; i<this->sigma.N(); i++){
        this->sigma2(i) = this->sigma(i) * this->sigma(i);
      } 


      TPrecision mse0 = Linalg<TPrecision>::Sum(l);
      this->mse = DenseVector<TPrecision>( this->phi.N()+1 );
      Linalg<TPrecision>::Zero( this->mse );
      this->mse(0) = mse0/X.N();

      for(unsigned int i = 0; i <  this->sigma2.N(); i++){
        this->mse(i+1) = this->mse(i) - this->sigma2(i);
        if(this->mse(i+1) < 0 ){
          this->mse(i+1) = 0;
        }
      }

      l.deallocate(); 
      U.deallocate();
      S.deallocate();
      Vt.deallocate();



    };

};










//GMRANode subclass
template <typename TPrecision> 
class IPCANode : public GMRANodeBase<TPrecision>{
  public:
    typedef typename GMRANode<TPrecision>::NodeVector NodeVector;
    typedef typename NodeVector::iterator NodeVectorIterator;


  private:
    NodeVector children;
    std::map<int, int> childmap;

  public:
    enum SplitDirectionStrategy {PC, RANDOM_PC};
    enum SplitStrategy {MEAN, MIDPOINT, MEDIAN, RANDOM_VARIANCE};

    FortranLinalg::DenseMatrix<TPrecision> dir;
    FortranLinalg::DenseVector<TPrecision> a; 



    virtual ~IPCANode(){
      dir.deallocate();    
      sigma.deallocate();
      sigma2.deallocate();
      phi.deallocate();
      center.deallocate();
      mse.deallocate();
    };



    IPCANode<TPrecision>(){
    };



    IPCANode<TPrecision>(FortranLinalg::DenseMatrix<TPrecision> &Xi, std::vector<int>
        &ind, NodeFactory<TPrecision> &nf, SplitDirectionStrategy splitDir = PC,
        SplitStrategy splitPoint = MEAN) : indices(ind){

      using namespace FortranLinalg;


      DenseMatrix<TPrecision> Xc = createMatrix(Xi);
      center = Linalg<TPrecision>::SumColumns(Xc);
      Linalg<TPrecision>::Scale(center, 1.0/Xc.N(), center);
      Linalg<TPrecision>::SubtractColumnwise(Xc, center, Xc);

      nf.createNode(Xc);
      phi = nf.getPhi();
      sigma = nf.getSigma();
      sigma2 = nf.getSigma2();
      mse = nf.getMSE();
      radius = nf.getRadius();

      if(splitDir == PC){
        dir = Linalg<TPrecision>::Copy(phi);
      }
      else if(splitDir == RANDOM_PC){
        static Random<TPrecision> random;
        DenseVector<TPrecision> tmp = DenseVector<TPrecision>(phi.M());
        for(unsigned int i=0; i<phi.N(); i++){
          double w = random.Uniform() * sigma(i);
          Linalg<TPrecision>::AddScale(tmp, w, phi, i, tmp);
        }
        Linalg<TPrecision>::Normalize(tmp);
        dir = DenseMatrix<TPrecision>(phi.M(), 1, tmp.data());
      }

      DenseVector<TPrecision> splitCenter;
      if(splitPoint == MEAN){
        splitCenter = Linalg<TPrecision>::Copy(center);
      }
      else if(splitPoint == MIDPOINT){
        DenseMatrix<TPrecision> Xp = Linalg<TPrecision>::Multiply(phi, Xc, true);
        DenseVector<TPrecision> mid = Linalg<TPrecision>::RowMin(Xp);
        DenseVector<TPrecision> maxP = Linalg<TPrecision>::RowMax(Xp);
        Linalg<TPrecision>::Add(mid, maxP, mid);
        Linalg<TPrecision>::Scale(mid, 0.5, mid);
        splitCenter = Linalg<TPrecision>::Multiply(phi, mid);
        Linalg<TPrecision>::Add(splitCenter, center, splitCenter);

        Xp.deallocate();
        mid.deallocate();
        maxP.deallocate();

      }
      else if(splitPoint == MEDIAN){
        /*
           DenseVector<TPrecision> a = Linalg<TPrecision>::Min(Vt);
           DenseVector<TPrecision> b = Linalg<TPrecision>::Max(Vt);
           Linalg<TPrecision>::Subtract(b, a, b);
           Linalg<TPrecision>::Scale(b, 0.5, b);
           Linalg<TPrecision>::Add(a, b, a);
           splitCenter = Linalg<TPrecision>::Multiply(U, a);
           Linalg<TPrecision>::Add(splitCenter, center, splitCenter);
           a.deallocate();
           b.deallocate();
         */
      }
      a = Linalg<TPrecision>::Multiply(dir, splitCenter, true);
      splitCenter.deallocate();

      Xc.deallocate();



    };


    void addChild(GMRANode<TPrecision> *node, int childIndex){
      children.push_back(node);
      childmap[childIndex] = children.size()-1;
    };




    GMRANode<TPrecision> *getChild(int childIndex){
      std::map<int, int>::iterator it = childmap.find(childIndex);
      if(it == childmap.end()){
        return NULL;
      }
      return children[it->second];
    };




    int getChildIndex(FortranLinalg::DenseVector<TPrecision> x){
      using namespace FortranLinalg;
      DenseVector<TPrecision> s = Linalg<TPrecision>::Multiply(dir, x, true);
      int childIndex = 0;
      int factor=1;
      for(int i=0; i<s.N(); i++){
        if( s(i) > a(i) ){
          childIndex += factor;
        }
        factor *= 2;
      }
      s.deallocate();
      return childIndex;

    };



    virtual GMRANode<TPrecision> *findDescendant( FortranLinalg::DenseVector<TPrecision> &x ){
       int index = getChildIndex(x);
       return getChild(index);
    };



    virtual FortranLinalg::DenseVector<TPrecision> affine(FortranLinalg::DenseVector<TPrecision> x){
      using namespace FortranLinalg;

      DenseVector<TPrecision> xm = Linalg<TPrecision>::Subtract(x, center);
      DenseVector<TPrecision> p = Linalg<TPrecision>::Multiply(phi, xm, true);
      xm.deallocate();
      return p;
    };



    FortranLinalg::DenseVector<TPrecision> reconstruct(FortranLinalg::DenseVector<TPrecision> p){
      return reconstruct(p, phi.N());
    };



    FortranLinalg::DenseVector<TPrecision> reconstruct(FortranLinalg::DenseVector<TPrecision> p, int d){
      using namespace FortranLinalg;
      for(int i=d; i<p.N(); i++){
        p(i) = 0;
      }
      DenseVector<TPrecision> x = Linalg<TPrecision>::Multiply(phi, p);
      Linalg<TPrecision>::Add(x, center, x);
      return x;
    };



    FortranLinalg::DenseVector<TPrecision> project(FortranLinalg::DenseVector<TPrecision> x){
      using namespace FortranLinalg;
      DenseVector<TPrecision> p = affine(x);
      DenseVector<TPrecision> xr = reconstruct(p);
      p.deallocate();
      return xr;
    };


    virtual NodeVector &getChildren(){
       return children;
    };


    std::vector<int> &getPoints(){
      return indices;
    };

    FortranLinalg::DenseVector<TPrecision> &getCenter(){
      return center;
    };

    TPrecision getRadius(){
      return radius;
    }



    void flatten(std::ostream &file){
      int nPhi = phi.N();
      file.write((char*)&nPhi, sizeof(int) );

      int nKids = children.size();
      file.write((char*)&nKids, sizeof(int) );

      for(std::map<int, int>::iterator it = childmap.begin(); it !=
          childmap.end(); it++){
        file.write((char*)&it->first , sizeof(int) );
        file.write((char*)&it->second, sizeof(int) );
      }

      file.write((char*)phi.data(), phi.M()*nPhi * sizeof(TPrecision) );
      file.write((char*)sigma.data(), sigma.N() * sizeof(TPrecision) );
      file.write((char*)center.data(), center.N() * sizeof(TPrecision) );
      file.write((char*)mse.data(), mse.N() * sizeof(TPrecision) );
      file.write((char*)dir.data(), dir.N() * dir.M()*sizeof(TPrecision) );
      file.write((char*)a.data(), a.N() * sizeof(TPrecision) );

      int nPoints = indices.size();
      file.write((char*) &nPoints, sizeof(int) );
      file.write((char*) indices.data(), indices.size()*sizeof(int) );

      file.write((char*) &radius, sizeof(TPrecision) );

    };




    int unflatten(std::istream &file, int m){
      using namespace FortranLinalg;
      
      int nPhi = 0;
      file.read( (char*) &nPhi, sizeof(int) );

      int nKids = 0;
      file.read( (char*) &nKids, sizeof(int) );

      for(int i=0; i<nKids; i++){
        int i1 = 0;
        file.read( (char*) &i1, sizeof(int) );
        int i2 = 0;
        file.read( (char*) &i2, sizeof(int) );
        childmap[i1] = i2;
      }

      
      phi =  DenseMatrix<TPrecision>(m, nPhi);
      file.read((char*)phi.data(), nPhi*m*sizeof(TPrecision));

      sigma = DenseVector<TPrecision> (nPhi);
      file.read((char*)sigma.data(), nPhi*sizeof(TPrecision));

      center = DenseVector<TPrecision>(m);
      file.read((char*)center.data(), m*sizeof(TPrecision));

      mse = DenseVector<TPrecision>(nPhi+1);
      file.read((char*)mse.data(), (nPhi+1)*sizeof(TPrecision));

      dir = DenseMatrix<TPrecision>(m, nPhi);
      file.read((char*)dir.data(), m*nPhi*sizeof(TPrecision));

      a = DenseVector<TPrecision>(nPhi);
      file.read((char*)a.data(), sizeof(TPrecision)*nPhi );

      int nPoints;
      file.read( (char*) &nPoints, sizeof(int) );
      indices.resize(nPoints);
      file.read((char*)indices.data(), nPoints*sizeof(int) );

      file.read((char*) &radius, sizeof(TPrecision) );



      sigma2 = DenseVector<TPrecision>( sigma.N() );
      for(int i=0; i<sigma.N(); i++){
        sigma2(i) = sigma(i) * sigma(i);
      } 
      return nKids;
    };

    virtual void translate(FortranLinalg::DenseVector<TPrecision> &x){
      using namespace FortranLinalg;
      Linalg<TPrecision>::Add(center, x, center);
    };


    virtual void affine(FortranLinalg::DenseMatrix<TPrecision> &A){
      using namespace FortranLinalg;
     
      DenseVector<TPrecision> centerN = Linalg<TPrecision>::Multiply(A, center);
      center.deallocate();
      center = centerN;
      
      DenseMatrix<TPrecision> phiN = Linalg<TPrecision>::Multiply(A, phi);
      phi.deallocate();
      phi = phiN;
    };



    FortranLinalg::DenseVector<TPrecision> sigma;
    FortranLinalg::DenseVector<TPrecision> sigma2;
    FortranLinalg::DenseMatrix<TPrecision> phi;
    FortranLinalg::DenseVector<TPrecision> center;
    TPrecision radius;


    FortranLinalg::DenseVector<TPrecision> mse;
    std::vector<int> indices;



  private:

    FortranLinalg::DenseMatrix<TPrecision> createMatrix(FortranLinalg::DenseMatrix<TPrecision> &X){
      FortranLinalg::DenseMatrix<TPrecision> M(X.M(), indices.size());

      std::vector<int>::iterator it = indices.begin();
      for(unsigned int i=0; i < M.N(); i++, ++it){
        FortranLinalg::Linalg<TPrecision>::SetColumn(M, i, X, *it);
      }

      return M;
    };

};






template <typename TPrecision>
class IPCATree : public GMRATree<TPrecision>{

  public:
    enum StoppingCriterium {TOTAL_R2, NODE_R2, NODE_MSE, NODE_RADIUS,
      RELATIVE_NODE_RADIUS};

  private:

    TPrecision epsilon;
    unsigned int minPoints;
    FortranLinalg::DenseMatrix<TPrecision> Xref;
    StoppingCriterium stop;
    NodeFactory<TPrecision> *nf;
    typedef typename IPCANode<TPrecision>::SplitStrategy SplitStrategy;
    typedef typename IPCANode<TPrecision>::SplitDirectionStrategy SplitDirectionStrategy;
    SplitStrategy splitStrategy;
    SplitDirectionStrategy splitDirectionStrategy;





    void buildTreeRecursive(IPCANode<TPrecision> *node, FortranLinalg::DenseMatrix<TPrecision>
        X, TPrecision rootMSE, TPrecision rootRadius){

      using namespace FortranLinalg;




      TPrecision mse = node->mse( node->mse.N() - 1);
      TPrecision relativeR2 = mse / rootMSE;
      TPrecision relativeRadius = node->getRadius() / rootRadius;
      TPrecision R2 = mse / node->mse(0);



#ifdef VERBOSE
      std::cout << "Node R^2 : " << R2 << std::endl;
      std::cout << "Node total R^2 : " << relativeR2 << std::endl;
      std::cout << "Node MSE : " << mse << std::endl;
      std::cout << "Node radius : " << node->getRadius() << std::endl;
      std::cout << "Node dim : " << node->phi.N() << std::endl << std::endl;
#endif


      if( node->getPoints().size() < std::max(minPoints, (unsigned int) 2 ) ){ 
        return;
      }
      if(stop == TOTAL_R2 && relativeR2 <= epsilon){
        return;
      }
      if(stop == NODE_R2 && R2 <= epsilon){
        return;
      }
      if(stop == NODE_MSE && mse <= epsilon){
        return;
      }
      if(stop == NODE_RADIUS && node->getRadius() <= epsilon){
        return;
      }
      if(stop == RELATIVE_NODE_RADIUS && relativeRadius <= epsilon){
        return;
      }

      int size = pow( (double) 2, (int) node->phi.N());
      std::vector< std::vector<int> > children(size);
      DenseVector<TPrecision> tmp(X.M());
      std::vector<int> &nodePts = node->getPoints();
      for(std::vector<int>::iterator it = nodePts.begin();
          it!=nodePts.end(); ++it){
        int i = *it;
        Linalg<TPrecision>::ExtractColumn(X, i, tmp);
        int childIndex = node->getChildIndex(tmp);
        children[childIndex].push_back(i);
      }
      tmp.deallocate();

      for(int i=0; i< children.size(); i++){
        if(children[i].size() > 0){
          IPCANode<TPrecision> *n = new IPCANode<TPrecision>(X, children[i], *nf,
              splitDirectionStrategy, splitStrategy);
          node->addChild(n, i);
          buildTreeRecursive(n, X, rootMSE, rootRadius);
        }
      }

    };





  public:


    IPCATree(){};



    IPCATree(TPrecision eps, StoppingCriterium stopC, NodeFactory<TPrecision>
        *factory, SplitStrategy ss, int minLeafSize=1): stop(stopC),
        epsilon(eps), splitStrategy(ss), minPoints(minLeafSize) 
    { 
          splitDirectionStrategy = IPCANode<TPrecision>::PC;
          nf = factory;
    };



    ~IPCATree(){
      
    };




    //
    void construct(TPrecision *X, int m, int n){

      FortranLinalg::DenseMatrix<TPrecision> Xd(m, n, X);
      Xref = Xd;

      std::vector<int> all;
      for(int i=0; i<n; i++){
        all.push_back(i);
      }

      //if(Xc.N() < nPoints){
      IPCANode<TPrecision> *root = new IPCANode<TPrecision>(Xd, all, *nf,
          splitDirectionStrategy, splitStrategy);
      //}
      //else{
      //  root = new IPCANode<TPrecision>(Xc);
      //}
      buildTreeRecursive( (IPCANode<TPrecision>*) root, Xd, root->mse(0),
          root->getRadius() ); 

      this->setRoot(root);
      this->setupParents();
    };




    FortranLinalg::DenseVector<TPrecision> getPoint(int index){
      FortranLinalg::DenseVector<TPrecision> x = FortranLinalg::Linalg<TPrecision>::ExtractColumn(Xref, index);
      return x;
    };




    void add(TPrecision *x){
      throw "not implemented yet";
    };




    std::vector<GMRANode<TPrecision> *> getLeafPath( double *x ) {
      using namespace FortranLinalg;
   
      GMRANode<TPrecision> *node = this->getRoot();
      DenseVector<TPrecision> xv( node->getCenter().N(), x);

      std::vector<GMRANode<TPrecision> *> path;
      while(node != NULL ){
        path.push_back( node );
        node = node->findDescendant( xv );
      }

      return path;
    }; 






    void flatten( std::ofstream &file ){
      using namespace FortranLinalg;
  
      std::list<IPCANode<TPrecision> *> nodes;
      IPCANode<TPrecision>* root = (IPCANode<TPrecision>*) this->getRoot();
      nodes.push_back(root);

      file.write( (char*) &epsilon, sizeof(TPrecision) ) ;
      unsigned int m = root->phi.M();
      file.write( (char*) &m, sizeof(unsigned int) ) ;   
      file.write( (char*) &minPoints, sizeof(unsigned int) );   

      while( !nodes.empty() ){
        IPCANode<TPrecision> *node =  nodes.front();
        nodes.pop_front();

        node->flatten(file);

        for(typename std::vector< GMRANode<TPrecision> * >::iterator it =
            node->getChildren().begin(); it != node->getChildren().end(); ++it){
          nodes.push_back((IPCANode<TPrecision>*) *it);
        }
      }

    };






    void unflatten( std::ifstream &file ){
      using namespace FortranLinalg;
    
      file.read( (char*) &epsilon, sizeof(TPrecision) ) ;
      unsigned int m;   
      file.read( (char*) &m, sizeof(unsigned int) );   
      file.read( (char*) &minPoints, sizeof(unsigned int) );   

      std::list<IPCANode<TPrecision> *> nodes;
      std::list<int> nKids;

      IPCANode<TPrecision> *cur = NULL;
      int nK = 0;
      while( !file.eof() ){
        
        IPCANode<TPrecision> *node = new IPCANode<TPrecision>();
        int n = node->unflatten(file, m);
        if(cur == NULL){
          this->setRoot(node);
          cur = node;
          nK = n;
        }
        else{ 
          if( n > 0 ){
            nodes.push_back(node);
            nKids.push_back(n);
          }
          if( cur->getChildren().size() == nK ){
            cur = nodes.front();
            nodes.pop_front();

            nK = nKids.front();
            nKids.pop_front();
          }
          cur->getChildren().push_back(node);

        }
        file.peek();
      }

      this->setupParents();

    };


};


#endif
