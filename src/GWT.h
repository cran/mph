#ifndef GWT_H
#define GWT_H

#include "DenseVector.h"
#include "DenseMatrix.h"

#include "SVD.h"
#include "RandomSVD.h"

#include <algorithm>
#include <map>

#include "IPCATree.h"
#include "GMRATree.h"
//Interface for GWT

template <typename TPrecision> 
class GWTNode : public GMRANodeDecorator<TPrecision>{
  
  private:
    FortranLinalg::DenseVector<TPrecision> w;

    FortranLinalg::DenseMatrix<TPrecision> psi;
    FortranLinalg::DenseMatrix<TPrecision> psi_t_phi;
   
    int ID;
  

  public:


    //Caller is responsible for cleaning up passed in variables properly
    //This class stores pointers to cx, s, T and p
    GWTNode<TPrecision>(FortranLinalg::DenseMatrix<TPrecision> phi,
        FortranLinalg::DenseVector<TPrecision> &center, int id,
        GWTNode<TPrecision> *p, GMRANode<TPrecision> *node ) :
      GMRANodeDecorator<TPrecision>(node) {
        
      using namespace FortranLinalg;
      
      ID = id;
      //phi_t_phi = Linalg<TPrecision>::Multiply(phi, phi, true);

      if(p != NULL){
        DenseMatrix<TPrecision> &p_phi = p->getPhi();

        DenseMatrix<TPrecision> phi_p =
          Linalg<TPrecision>::Multiply(p_phi, phi, true);
        psi = Linalg<TPrecision>::Multiply(p_phi, phi_p);
        Linalg<TPrecision>::Subtract(phi, psi, psi);
        SVD<TPrecision> svd(psi);
        int cut = svd.S.N();
        for(int i=0; i<svd.S.N(); i++){
          if( svd.S(i) < 0.00000000001 ){
            cut = i;
            break;
          }
        }
        psi.deallocate();
        if(cut == 0){
          psi = DenseMatrix<TPrecision>(0, 0);
        }
        else{
          psi = Linalg<TPrecision>::ExtractColumns(svd.U, 0, cut);
        }
        svd.deallocate();

        
        if(cut == 0){
          psi_t_phi = DenseMatrix<TPrecision>(0, 0);
        }
        else{
          psi_t_phi = Linalg<TPrecision>::Multiply(psi, phi, true);
        }

        DenseVector<TPrecision> t = Linalg<TPrecision>::Subtract(center,
            p->getCenter() );
        DenseVector<TPrecision> pt = Linalg<TPrecision>::Multiply(p_phi, t, true);
        w = Linalg<TPrecision>::Multiply(p_phi, pt);
        Linalg<TPrecision>::Subtract(t, w, w);

        t.deallocate();
        pt.deallocate();
        phi_p.deallocate();
      }

    };


    //Onlyd eallocate newly created objects
    virtual  ~GWTNode(){
      psi.deallocate();
      psi_t_phi.deallocate();
      w.deallocate();

    };


    //FortranLinalg::DenseMatrix<TPrecision> phi;
    //FortranLinalg::DenseVector<TPrecision> center;
    //FortranLinalg::DenseVector<TPrecision> sigma;
    
    //The variables below are deallcoated in the deconstructor

    //DenseMatrix<TPrecision> phi_t_phi;

    FortranLinalg::DenseMatrix<TPrecision> &getPsi(){
      return psi;
    };

    FortranLinalg::DenseMatrix<TPrecision> &getPsitPhi(){
      return psi_t_phi;
    };


    FortranLinalg::DenseVector<TPrecision> &getW(){
      return w;
    };

    int getID(){
      return ID;
    };

    virtual FortranLinalg::DenseMatrix<TPrecision> &getPhi() = 0;
    virtual FortranLinalg::DenseVector<TPrecision> &getSigma() = 0;
    virtual FortranLinalg::DenseVector<TPrecision> &getCenter() = 0;

};






template <typename TPrecision> 
class GWTCoefficients{

  public:
    GWTCoefficients(){
      maxD = 0;
    };

    //Matrix of coefficents at each 
    //scale from finest (0) to coarsest (coeff.M()-1)
    std::vector< FortranLinalg::DenseVector<TPrecision> > coeff;
    std::vector< int > ids;
    int maxD;

    //Node ids for each scale
    GMRANode<TPrecision> *leaf;
};





template <typename TPrecision>
class GWT{

  private:
    GMRATree<TPrecision> *tree;


  protected:
    //Subclass needs to decorate ech GMRATree node info to each GMRANode in the tree.
    virtual void setupGWT(GMRATree<TPrecision> *tree) = 0;

    
  public:

    
    //Setup tree for GWT, i.e. add covariance structur to each node -> see
    //setupGWT
    void setTree(GMRATree<TPrecision> *gTree){
      tree = gTree;
      setupGWT(tree);
    };

    

    //Encoding of GMRA coefficents 
    GWTCoefficients<TPrecision> encode(FortranLinalg::DenseVector<TPrecision> x){
        using namespace FortranLinalg;
      GMRANode<TPrecision> *node = tree->getLeafPath(x.data()).back();
      GWTNode<TPrecision> *gwt = dynamic_cast<GWTNode<TPrecision> *>( node );
      GWTNode<TPrecision> *parent = dynamic_cast<GWTNode<TPrecision> *>(
          node->getParent() );
      
      GWTCoefficients<TPrecision> c;
      c.leaf = node;

      std::list< DenseVector<TPrecision> > ctmp;


      DenseVector<TPrecision> delta = Linalg<TPrecision>::Subtract(x,
          gwt->getCenter());
      DenseVector<TPrecision> p = Linalg<TPrecision>::Multiply(gwt->getPhi(), delta, true);
      DenseVector<TPrecision> x_J = Linalg<TPrecision>::Multiply(gwt->getPhi(), p);
      Linalg<TPrecision>::Add(x_J, gwt->getCenter(), x_J);

      int d = 0;
      while(parent != NULL){
        DenseVector<TPrecision> q;
        if(gwt->getPsi().N() != 0){
          q = Linalg<TPrecision>::Multiply( gwt->getPsitPhi(), p );
        }
        c.coeff.push_back( q );
        c.ids.push_back( gwt->getID() );
        if(q.N() > c.maxD){
          c.maxD = q.N();
        }

        Linalg<TPrecision>::Subtract(x_J, parent->getCenter(), delta);
        p.deallocate();
        p = Linalg<TPrecision>::Multiply(parent->getPhi(), delta, true);

        gwt = parent;
        parent = dynamic_cast< GWTNode<TPrecision> *>( gwt->getParent() );
      }

      c.coeff.push_back(p);
      c.ids.push_back(gwt->getID());
      std::reverse(c.coeff.begin(), c.coeff.end() );
      std::reverse(c.ids.begin(), c.ids.end() );
      if(p.N() > c.maxD){
        c.maxD = p.N();
      }


      delta.deallocate();
      x_J.deallocate();
      return c;
    };




    //Decoding from GMRA coefficents at specific scale (0 = coarsest,
    //a.coeff.M() =finest also -1 = finest)
    FortranLinalg::DenseVector<TPrecision> decode(GWTCoefficients<TPrecision> &a, int scale =-1){
        using namespace FortranLinalg;
      GMRANode<TPrecision> *gmraNode = a.leaf; 
      GWTNode<TPrecision> *gwt = (GWTNode<TPrecision> *) gmraNode;
      GWTNode<TPrecision> *parent = (GWTNode<TPrecision> *)
        gmraNode->getParent();
      
      DenseVector<TPrecision> sumq( gwt->getCenter().N() );  

      int nScales = 0;
      GWTNode<TPrecision> *gwtTmp = gwt;
      while(gwtTmp->getParent() != NULL){
        nScales++;
        gwtTmp = dynamic_cast< GWTNode<TPrecision> *>( gwtTmp->getParent() );
      }
       
      if(scale < 0 || scale > nScales){
        scale = nScales;
      }
      int scaleIndex = nScales;
      while(scaleIndex != scale && gwt->getParent() != NULL){
          scaleIndex--;
          gwt = dynamic_cast< GWTNode<TPrecision> *>( gwt->getParent() );
      }
      if(parent != NULL){
        DenseVector<TPrecision> q = a.coeff[scaleIndex];
        if(gwt->getPsi().N() != 0){
          q.shorten(gwt->getPsi().N());
          Linalg<TPrecision>::Multiply(gwt->getPsi(), q, sumq); 
        }
        else{
          Linalg<TPrecision>::Zero(sumq);
        }
        Linalg<TPrecision>::Add(sumq, gwt->getW(), sumq); 

        DenseVector<TPrecision> qtmp1( sumq.N() );
        DenseVector<TPrecision> qtmp2( sumq.N() );
        
        gwt = parent;
        parent = (GWTNode<TPrecision> *) gwt->getParent();

        scaleIndex--;
        while(parent != NULL){
          q = a.coeff[scaleIndex];

          if(gwt->getPsi().N() != 0){
            q.shorten(gwt->getPsi().N());  
            Linalg<TPrecision>::Multiply( gwt->getPsi(), q, qtmp1); 
          }
          else{
            Linalg<TPrecision>::Zero(qtmp1);
          }
          Linalg<TPrecision>::Add(qtmp1, gwt->getW(), qtmp1);
          
          q = Linalg<TPrecision>::Multiply( parent->getPhi(), sumq, true);
          Linalg<TPrecision>::Multiply( parent->getPhi(), q, qtmp2);
          q.deallocate();
        
          Linalg<TPrecision>::Subtract(qtmp1, qtmp2, qtmp1);
          Linalg<TPrecision>::Add(sumq, qtmp1, sumq);
          //Linalg<TPrecision>::Add(qtmp2, qtmp1, sumq);
          
          gwt = parent;
          parent = dynamic_cast< GWTNode<TPrecision> *>( gwt->getParent() );
          scaleIndex--;
        }
        qtmp1.deallocate();
        qtmp2.deallocate();
      }
      else{
        Linalg<TPrecision>::Set(sumq, 0);
      }

      DenseVector<TPrecision> p = a.coeff[scaleIndex];
      DenseVector<TPrecision> x = Linalg<TPrecision>::Multiply(gwt->getPhi(), p);
      
      Linalg<TPrecision>::Add(x, gwt->getCenter(), x);
      Linalg<TPrecision>::Add(x, sumq, x);
      

      sumq.deallocate();

      return x;
    };



};

#endif
