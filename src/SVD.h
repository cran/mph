#ifndef SVD_H
#define SVD_H

#include "Linalg.h"
#include "LapackDefs.h"

namespace FortranLinalg{

template <typename TPrecision>
class SVD{
#define isDoubleTPrecision() sizeof(TPrecision) == sizeof(double)

  public:
    DenseMatrix<TPrecision> U;
    DenseVector<TPrecision> S;
    DenseMatrix<TPrecision> Vt;
    DenseVector<TPrecision> c;


    SVD(DenseMatrix<TPrecision> Xo, bool center=false){

      DenseMatrix<TPrecision> Xc = Linalg<TPrecision>::Copy(Xo);

      if(center){
        c = Linalg<TPrecision>::SumColumns(Xc);
        Linalg<TPrecision>::Scale(c, 1.0/Xc.N(), c);
        Linalg<TPrecision>::SubtractColumnwise(Xc, c, Xc);
      }

      int n = Xc.N();
      int m = Xc.M();
      int s = std::min(m,n);
      
      S  = DenseVector<TPrecision>( s );
      U  = DenseMatrix<TPrecision>( m, s );
      Vt = DenseMatrix<TPrecision>( s, n );
      
      TPrecision *work = new TPrecision[1];
      int lwork = -1;
      char jobu = 'S';
      char jobv = 'S';
      int info = 0;

      if(isDoubleTPrecision()){
        lapack::dgesvd_(&jobu, &jobv, &m, &n, (double*)Xc.data(), &m,
          (double*)S.data(), (double*)U.data(), &m, (double*)Vt.data(), &s,
          (double*)work, &lwork, &info);

        //std::cout << info << std::endl;
        int lwork = work[0];
        delete[] work;
        work =new TPrecision[lwork];        
        lapack::dgesvd_(&jobu, &jobv, &m, &n, (double*)Xc.data(), &m,
          (double*)S.data(), (double*)U.data(), &m, (double*)Vt.data(), &s,
          (double*)work, &lwork, &info);
      }
      else{
        lapack::sgesvd_(&jobu, &jobv, &m, &n, (float*)Xc.data(), &m,
          (float*)S.data(), (float*)U.data(), &m, (float*)Vt.data(), &s,
          (float*)work, &lwork, &info);

        int lwork = work[0];
        delete[] work;
        work =new TPrecision[lwork];        
        lapack::sgesvd_(&jobu, &jobv, &m, &n, (float*)Xc.data(), &m,
          (float*)S.data(), (float*)U.data(), &m, (float*)Vt.data(), &s,
          (float*)work, &lwork, &info);
      } 
      Xc.deallocate();
      delete[] work;
    };

    void deallocate(){
      U.deallocate();
      S.deallocate();
      Vt.deallocate();
      c.deallocate();
    };

};

}

#endif 
