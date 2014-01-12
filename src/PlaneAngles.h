#ifndef PLANEANGLES_H
#define PLANEANGLES_H

#include "DenseVector.h"
#include "DenseMatrix.h"

#include "EuclideanMetric.h"

#include "SVD.h"
#include "RandomSVD.h"

#include <algorithm>
#include <map>

#include "GWT.h"


template <typename TPrecision>
class PlaneAngles{


public:
  

   static DenseVector<TPrecision> angles(DenseMatrix<TPrecision> &phi1, DenseMatrix<TPrecision> &phi2){
     DenseMatrix<TPrecision> phi_p = Linalg<TPrecision>::Multiply(phi1, phi2, true);
     DenseMatrix<TPrecision> psi = Linalg<TPrecision>::Multiply(phi1, phi_p);
     Linalg<TPrecision>::Subtract(phi2, psi, psi);
     SVD<TPrecision> svd(psi);
   
     svd.U.deallocate();
     svd.Vt.deallocate();
     phi_p.deallocate();
     psi.deallocate();
   
     return svd.S;
   };

};




#endif
