#ifndef NULL
#define NULL 0
#endif

#define R_NO_REMAP


#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdio.h>

#include "DenseMatrix.h"

#include "IPCATree.h"
#include "GMRATree.h"
#include "GMRANeighborhood.h"

#include "MultiscaleMappingRips.h"

#include "NodeDistance.h"

#include "Linalg.h"



extern "C" {


  /**/  
  /* Multiscale Topology Calls */

  //Multiscale Rips from GMRA tree
  SEXP runMultiscaleRips(FortranLinalg::DenseMatrix<double> &X, GMRATree<double>
      &tree,  int maxD, bool single){

    using namespace FortranLinalg;

    MultiscaleMappingRips<double> mrips(X, tree);
    mrips.run(maxD, single);


    std::map<int, DenseMatrix<double> > &dgms = mrips.getDiagrams(); 

    int size = 0;
    for(std::map<int, DenseMatrix<double> >::iterator it = dgms.begin(); it !=
        dgms.end(); ++it){
      size += it->second.N();
    };

    DenseMatrix<double> dgm(size, 4);
    int index = 0;
    for(std::map<int, DenseMatrix<double> >::iterator it = dgms.begin(); it !=
        dgms.end(); ++it){
      DenseMatrix<double> dg = it->second;
      for(unsigned int i=0; i< dg.N(); i++, index++){
        dgm(index, 0) = dg(0, i);
        dgm(index, 1) = dg(1, i);
        dgm(index, 2) = it->first;
        dgm(index, 3) = dg(2, i);
      }
      dg.deallocate();
    }

    SEXP Rdgm;
    PROTECT( Rdgm = Rf_allocMatrix(REALSXP, dgm.M(), dgm.N()) );
    memcpy( REAL(Rdgm), dgm.data(), dgm.N()*dgm.M()*sizeof(double) );
    UNPROTECT(1);

    dgm.deallocate();

    return Rdgm;

  };

/*
  //Multiscale Rips from GMRA tree
  SEXP runMultiscaleRipsHarer(FortranLinalg::DenseMatrix<double> &X, GMRATree<double> &tree,
      int maxD, bool single){
    using namespace FortranLinalg;

    MultiscaleRipsHarer<double> mrips(X, tree);
    mrips.run( maxD, single );

    MultiscalePersistentHomology<double> persistence;
    persistence.reduce( mrips.getFiltrations(), mrips.getMappings(),
        mrips.getAlphas(), maxD ); 

    DenseMatrix<double> dgmT = persistence.getDiagram();
    DenseMatrix<double> dgm = Linalg<double>::Transpose(dgmT);

    SEXP Rdgm;
    PROTECT( Rdgm = Rf_allocMatrix(REALSXP, dgm.M(), dgm.N()) );
    memcpy( REAL(Rdgm), dgm.data(), dgm.N()*dgm.M()*sizeof(double) );
    UNPROTECT(1);

    dgm.deallocate();
    dgmT.deallocate();

    return Rdgm;

  };

*/


  SEXP multiscale_rips_ipca(SEXP Rx, SEXP Rn, SEXP Rm, SEXP RepsGMRA, 
      SEXP Rd, SEXP Rsingle) {
    using namespace FortranLinalg;

    int d = *INTEGER(Rd);
    int m = *INTEGER(Rm);
    int n = *INTEGER(Rn);
    double *x = REAL(Rx);
    double epsGMRA = *REAL(RepsGMRA);
   // int type = *INTEGER(Rtype);
    bool single = *INTEGER(Rsingle) != 0;

    DenseMatrix<double> X(m, n, x);

    FixedNodeFactory<double> factory(1);
    //RelativePrecisionNodeFactory<double> factory(5, 0.9);
    IPCATree<double> ipcaTree(epsGMRA, IPCATree<double>::RELATIVE_NODE_RADIUS,
        &factory, IPCANode<double>::MIDPOINT );
    ipcaTree.construct(X.data(), m, n);


   // Rprintf("Tree constructed \n");



    //SEXP Rdgm;
    //if(type == 0){
    SEXP Rdgm = runMultiscaleRips(X, ipcaTree, d, single);
    //}
    //else{
    //  Rdgm = runMultiscaleRipsHarer(X, ipcaTree, d, single);
    //}

    return Rdgm;  
  };


/*
  //Multiscale Rips using stored GMRA tree
  SEXP multiscale_rips_gmra(SEXP Rx, SEXP Rn, SEXP Rm,
      SEXP Rd, SEXP RtreeID, SEXP Rtype) {
    using namespace FortranLinalg;

    int d = *INTEGER(Rd);
    int m = *INTEGER(Rm);
    int n = *INTEGER(Rn);
    double *x = REAL(Rx);
    int treeID = *INTEGER(RtreeID);
    int type = *INTEGER(Rtype);

    DenseMatrix<double> X(m, n, x);

    GMRATree<double> *tree = trees[treeID];
    if(tree == NULL){
      return R_NilValue;
    }



    SEXP Rdgm;
    if(type == 0){
      Rdgm = runMultiscaleRips(X, *tree, d);
    }
    else{
      Rdgm = runMultiscaleRipsHarer(X, *tree,  d);
    }

    return Rdgm;  
  };

*/




}//end extern C
