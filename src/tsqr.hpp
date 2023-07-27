#ifndef __TSQR_HPP__
#define __TSQR_HPP__

#include "ops.hpp"
#include "svd.hpp"
#include "qr_decomposition.hpp"
#include "processor.hpp"



template<typename CT, typename T>
void performTSQR(matrix<CT>& A, matrix<CT>& U_local, matrix<CT>& V_T_r, matrix<CT>& Sigma) {
    int ierr;
    using size_type = typename matrix<CT>::size_type;

    auto std_inner = [](T v1, T v2) {
        return v1 * v2;
    };

    pout << "Doing QR on A1 in TSQR: \n";
//    std::pair<matrix<T>, matrix<T>>  qr_pair = qr_decomposition(A, std_inner);
#ifdef CheckQR
//		matrix<CT>B(4,3);
//		CT a(-1,0);
//		B(0,0)=a;
//		CT b(-1,0);
//		B(0,1)=b;
//		CT c(1,0);
//		B(0,2)=c;
//		CT d(1,0);
//		B(1,0)=d;
//		CT e(3,0);
//		B(1,1)=e;
//		CT f(3,0);
//		B(1,2)=f;
//		CT g(-1,0);
//		B(2,0)=g;
//		CT h(-1,0);
//		B(2,1)=h;
//		CT i(5,0);
//		B(2,2)=i;
//		CT j(1,0);
//		B(3,0)=j;
//		CT k(3,0);
//		B(3,1)=k;
//		CT l(7,0);
//		B(3,2)=l;
#endif		
		
    std::pair<matrix<CT>, matrix<CT>>  qr_pair =computeQRFactorizationUsingLAPACKOfAComplexMatrix(A);

    MPI_Datatype local_RMatType;
    MPI_Datatype global_RMatType;

    MPI_Type_vector( qr_pair.second.colCount(), qr_pair.second.rowCount(), qr_pair.second.rowCount(), MPI_DOUBLE_COMPLEX, &local_RMatType );
    MPI_Type_vector( qr_pair.second.colCount(), qr_pair.second.rowCount(), qr_pair.second.rowCount(), MPI_DOUBLE_COMPLEX, &global_RMatType);

    MPI_Type_commit( &local_RMatType);
    MPI_Type_commit( &global_RMatType);

    matrix<CT> r_global(qr_pair.second.rowCount()*size_type(MPI_PROC_TOTAL_NUM), qr_pair.second.colCount());
    ierr = MPI_Allgather(&qr_pair.second(0,0), 1, local_RMatType, &r_global(0,0), 1, global_RMatType, MPI_COMM_WORLD);
    MPI_Type_free( &local_RMatType);
    MPI_Type_free( &global_RMatType);

    qr_pair.second.destroy();

#ifdef Print
    pout << "r_global after all gather: \n";
    pout << r_global << "\n";
#endif

    pout << "Doing QR on r_global: \n";
//    std::pair<matrix<T>, matrix<T>> qr_pair_fromR = qr_decomposition(r_global, std_inner);
    std::pair<matrix<CT>, matrix<CT>> qr_pair_fromR =	computeQRFactorizationUsingLAPACKOfAComplexMatrix(r_global);

    matrix<CT> q_local = matMatMul_QoldQnew_ForParallelQR(qr_pair.first, qr_pair_fromR.first);
    qr_pair.first.destroy();
    qr_pair_fromR.first.destroy();
#ifdef Print
    pout << "Q_local: \n";
    pout << q_local ;
#endif

    pout << "Feeding R matrix to SVD: \n";
    size_type R_rows = qr_pair_fromR.second.rowCount();
    size_type R_cols = qr_pair_fromR.second.colCount();
    matrix<CT> U_r(R_rows, R_rows);
//	  matrix<CT> qr_pair_fromR_second_complex(qr_pair_fromR.second.rowCount(),qr_pair_fromR.second.colCount());
//  	constructComplexMatrixFromRealAndImagParts(qr_pair_fromR_second_complex,qr_pair_fromR.second,qr_pair_fromR.second, 0);

    int info =  computeComplexSVDUsingLAPACK<CT,T>(qr_pair_fromR.second, U_r, Sigma, V_T_r);
    pout << "info: " << info << "\n";

		qr_pair_fromR.second.destroy();

//		matrix<CT> q_local_complex(q_local.rowCount(), q_local.colCount());
//  	constructComplexMatrixFromRealAndImagParts(q_local_complex,q_local, q_local, 0);

    U_local = q_local*U_r;

		q_local.destroy();
    U_r.destroy();


}
#endif //TSQR
