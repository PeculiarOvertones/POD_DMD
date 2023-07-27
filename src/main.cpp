#include <fstream>
#include "readIO.hpp"
#include "nvector.hpp"
#include "ops.hpp"
#include "processor.hpp"
#include "svd.hpp"
#include "tsqr.hpp"
#include <math.h>
#include <complex>

#include <mpi.h>
int main(int argc, char *argv[]) {
	using std::cout;
	using Type = double;
	using ComplexType = std::complex<Type>;
	using size_type = unsigned long int;
	int ierr;
	defineMPIVarsAndProcFiles(argc, argv);
	MPI_Barrier(MPI_COMM_WORLD);

	container<nvector<ComplexType>> DataVector; 
#ifdef OzgurN2Case
    container<Type> data;
	matrix<Type> A1(40000, numSnapShots);
    std::string iFd  ="../InputFiles_N2/NOTvib_180_670K_500Steps/SplitFiles/Tvib_NO_Split_";
	std::string pID;
	if(EntireData){
		std::cout << "EntireData directive: true \n";
    	pID =static_cast<std::ostringstream*>( &(std::ostringstream() << MPI_PROC_ID) )->str();
	} else {
		std::cout << "EntireData directive: false \n";
    	pID =static_cast<std::ostringstream*>( &(std::ostringstream() << PlaneID) )->str();
	}
    std::string fE=".dat";
    std::string fname = iFd + pID + fE;
    std::ifstream infile(fname.c_str(), std::ios_base::in);
//    std::ifstream infile("../InputFiles_N2/NOTvib_180_670K_500Steps/Tvib_NO.dat", std::ios_base::in);
//    std::cout << "reading: " << fileName << "\n";
    if(infile.fail()){
        pout << "Problem reading file! \n";
        throw std::invalid_argument("Problem reading file in readInput");
    }else{
		pout << "reading file: " << fname.c_str() << "\n"; 
        std::stringstream buffer;
        buffer << infile.rdbuf();
        infile.close();
        data = split<Type>(buffer.str(),'\t');
//		std::cout << data << "\n";
		int fileCounter=0;
		int counter=0;
        pout << "Done reading data, gonna put it in A1! \n";
//		container<Type> tempV;
		for (auto it = data.begin(); it != data.end(); ++it)
   		{
				size_type col = counter%numSnapShots;
				size_type row = size_type(counter/numSnapShots);
				A1(row, col) = *it;
   		    counter++;
//   		    tempV.push_back(*it);
//   		    if(counter%(numSnapShots) ==0){
//				fileCounter++;
//				if(tempV.size() > 0){
//					DataVector.push_back(tempV);
//				}
//				tempV.clear();
//   		    }
   		}
	    data.clear();
    	container<Type>(data).swap(data);
    }
    pout << "Done forming A1! \n";
//	matrix<Type> A1 = construct_from_row_vectors(DataVector);
//	pout << A1 << "\n";
#else
	//Data Snapshots
#ifdef DMD
	size_type TimeSnapsToRead;
	if(MPI_PROC_TOTAL_NUM==1){
		TimeSnapsToRead = numSnapShots-TimeShift;
	} else {
		TimeSnapsToRead = numSnapShots-TimeShift-1;
	}
#else
	size_type TimeSnapsToRead = (numSnapShots-TimeShift);
#endif

    for(size_type i=0; i< TimeSnapsToRead; ++i) { 
		size_type snapShot = snapShotStart + snapShotInterval*i;
		if(EntireData){
			size_type fileIDToRead = MPI_PROC_ID%NFiles;
			if(MPI_PROC_ID >= NFiles){
				snapShot += size_type(MPI_PROC_ID/NFiles);
			} 
			std::cout  << MPI_PROC_ID << " fileIDToRead: " << fileIDToRead << " snapShot: " << snapShot << " \n";
			std::cout << "EntireData directive: true \n";
			DataVector.push_back(nvector<ComplexType>(readInput<ComplexType>(snapShot, fileIDToRead)));
		} else {
			std::cout << "EntireData directive: false \n";
			DataVector.push_back(nvector<ComplexType>(readInput<ComplexType>(snapShot, PlaneID)));
		}
	}

#ifdef SingleProc
		#ifdef DMD
		matrix<ComplexType> A1 = time_shifted_construct_from_column_vectors(DataVector,0, DataVector.size()-TimeShift);
		#else
		matrix<ComplexType> A1 = construct_from_column_vectors(DataVector,0, DataVector.size());
		#endif
#else
		matrix<ComplexType> A1 = construct_from_column_vectors(DataVector,0, DataVector.size());
#endif
	pout << "A1 rows/cols: " << A1.rowCount() << " " << A1.colCount() << "\n";

//	matrix<ComplexType> A1(4,3);
//	ComplexType a(1,0);
//	ComplexType b(2,0);
//	ComplexType c(3,0);
//	ComplexType d(4,0);
//	ComplexType e(5,0);
//	ComplexType f(6,0);
//	ComplexType g(7,0);
//	ComplexType h(8,0);
//	ComplexType i(9,0);
//	ComplexType j(10,0);
//	ComplexType k(11,0);
//	ComplexType l(12,0);
//	A1(0,0) = a;
//	A1(0,1) = b;
//	A1(0,2) = c;
//	A1(1,0) = d;
//	A1(1,1) = e;
//	A1(1,2) = f;
//	A1(2,0) = g;
//	A1(2,1) = h;
//	A1(2,2) = i;
//	A1(3,0) = j;
//	A1(3,1) = k;
//	A1(3,2) = l;
//	pout << "A1 dummy: \n";
//	pout << A1 << "\n";
	//Free container<nvector<T>>
    for (size_type i = 0; i < DataVector.size(); ++i) {
        DataVector[i].destroy();
    }
    DataVector.clear();
    container<nvector<ComplexType>>(DataVector).swap(DataVector);
#endif

#ifdef CovariancePOD 
	size_type NModes = SpatialModesToOP;
//	matrix<Type> AT_A = A1*transpose(A1); /*1st way to feed: Only use for small A1s*/
	matrix<Type> AT_A = transpose(A1)*A1; /*2nd way to feed*/
	matrix<Type> L_r(AT_A.rowCount(), AT_A.rowCount());
	matrix<Type> L_i(AT_A.rowCount(), AT_A.rowCount());
	matrix<Type> V_R(AT_A.rowCount(), AT_A.rowCount());
	matrix<Type> V_I(AT_A.rowCount(), AT_A.rowCount());
	computeEigenValuesUsingLAPACK(AT_A, L_r, L_i, V_R, V_I);	
	AT_A.destroy();
	
#ifdef Print
	pout << "L_r\n";
	pout << L_r << "\n";
#endif
//	matrix<Type>SpatialMode(A1.rowCount(),NModes); /*Use for first way to feed*/
//    for(size_type i=0; i<NModes; ++i){
//    	for(size_type k=0; k<SpatialMode.rowCount(); ++k){
//			SpatialMode(k,i) = V_R(k,i);
//		}
//	}
	matrix<Type>SpatialMode(A1.rowCount(),NModes); /*Use for 2nd way to feed*/
    for(size_type i=0; i<NModes; ++i){
		container<Type> tempV = A1*V_R(":",i);
		Type norm = l2_norm(tempV);
    	for(size_type k=0; k<SpatialMode.rowCount(); ++k){
			SpatialMode(k,(i)) = tempV[k]/norm; /*Check sometimes eigenvalues are arranged in ascending order in L_r*/
		}
	}
//	reverseArrangeDiagonal(L_r);
	
#ifdef Print
	pout << "L_r\n";
	pout << L_r << "\n";
	pout << "SpatialMode: " << "\n";
	pout << SpatialMode << "\n";
#endif
	matrix<Type>TemporalMode_T = transpose(SpatialMode)*A1;
//	matrix<Type>TemporalMode = transpose(transpose(A1)*SpatialMode);
	Type TotalModeEnergy = 0.0;
    for(size_type i=0; i<L_r.rowCount(); ++i){
		TotalModeEnergy += L_r(i,i);
    }
    pout << TotalModeEnergy << "\n";
	matrix<Type> RelativeModeEnergy(NModes,NModes);
    for(size_type i=0; i<L_r.rowCount(); ++i){
		RelativeModeEnergy(i,i) = L_r(i,i)/TotalModeEnergy;
    }
#ifdef Print
	pout << "TemporalMode_T: " << "\n";
	pout << TemporalMode_T << "\n";
#endif
	if(MPI_PROC_ID == 0) {
		PrintModesAndEigValForCovariancePOD(SpatialMode, TemporalMode_T, L_r, RelativeModeEnergy, "CovariancePOD");
	}

#ifdef Print
	pout << "Spatial*Temporal: \n";
	pout << SpatialMode*TemporalMode_T << "\n";
#endif
	L_r.destroy();	
	L_i.destroy();
	V_R.destroy();	
	V_I.destroy();	
	TemporalMode_T.destroy();
	SpatialMode.destroy();	



#endif

//Converting data matrix to complex type 
//	matrix<ComplexType> A1_complex(A1.rowCount(), A1.colCount());
//  constructComplexMatrixFromRealAndImagParts(A1_complex,A1,A1, 0);
//	A1.destroy();


#ifdef DirectSVD
  //Direct SVD
	pout << "Feeding A1 matrix to SVD in DirectSVD: \n";
	matrix<ComplexType> U_svd(A1.rowCount(), A1.rowCount());
	matrix<ComplexType> V_T_svd(A1.colCount(), A1.colCount());

	matrix<ComplexType> Sigma_svd(MAX(A1.rowCount(), A1.colCount()),MAX(A1.rowCount(), A1.colCount()), 0 );
	int info = computeComplexSVDUsingLAPACK<ComplexType,Type>(A1, U_svd, Sigma_svd, V_T_svd);

	PrintModesAndEigVal(U_svd, V_T_svd, Sigma_svd, "DirSVD");
#elif Print
	pout << "A1 data matrix: \n";
	pout << A1 << "\n";
#endif

MPI_Barrier(MPI_COMM_WORLD);
	
#ifdef TSQR
	/*Do QR factorization of A*/
    matrix<ComplexType> U_local_full(A1.rowCount(), A1.colCount());
    matrix<ComplexType>   V_T_r_full(A1.colCount(), A1.colCount());
    matrix<ComplexType>   Sigma_full(A1.colCount(), A1.colCount(), 0 );
		performTSQR<ComplexType, Type>(A1, U_local_full, V_T_r_full, Sigma_full);

    matrix<ComplexType> U_local(A1.rowCount(), size_type(RANK));
		for (size_type c=0; c<size_type(RANK); ++c){
			for (size_type r=0; r<U_local.rowCount(); ++r){
	  		U_local(r,c) = U_local_full(r,c);	
			}
		}
		U_local_full.destroy();

    matrix<ComplexType>   V_T_r(size_type(RANK), A1.colCount());
		for (size_type c=0; c<V_T_r.colCount(); ++c){
			for (size_type r=0; r<size_type(RANK); ++r){
	  	  V_T_r(r,c) = V_T_r_full(r,c);	
			}
		}

		V_T_r_full.destroy();

    matrix<ComplexType> Sigma(size_type(RANK), size_type(RANK), 0 );
		for (size_type c=0; c<size_type(RANK); ++c){
			for (size_type r=0; r<size_type(RANK); ++r){
	  	  Sigma(r,c) = Sigma_full(r,c);	
			}
		}
		Sigma_full.destroy();

		A1.destroy();
#ifdef Print
	pout << "V_T_r: \n";
	pout << V_T_r << "\n";
	pout << "Sigma: \n";
	pout << Sigma << "\n";
	pout << "U_local: \n";
	pout << U_local << "\n";
#endif
    matrix<ComplexType> U_global(U_local.rowCount()*size_type(MPI_PROC_TOTAL_NUM), U_local.colCount());

    MPI_Datatype local_UMatType;
    MPI_Datatype global_UMatType;

    MPI_Type_vector( U_local.colCount(), U_local.rowCount(), U_local.rowCount(), MPI_DOUBLE_COMPLEX, &local_UMatType );
    MPI_Type_vector( U_local.colCount(), U_local.rowCount(), U_local.rowCount(), MPI_DOUBLE_COMPLEX, &global_UMatType);

    MPI_Type_commit( &local_UMatType);
    MPI_Type_commit( &global_UMatType);

    ierr = MPI_Allgather(&U_local(0,0), 1, local_UMatType, &U_global(0,0), 1, global_UMatType, MPI_COMM_WORLD);
    MPI_Type_free( &local_UMatType);
    MPI_Type_free( &global_UMatType);

#ifdef Print
    pout << "U_global: after all gather: \n";
    pout << U_global << "\n";
    pout << "U_global*Sigma*V_T_r: \n";
    pout << U_global*Sigma*V_T_r << "\n";
#endif

	if(MPI_PROC_ID == 0) {
		PrintModesAndEigVal(U_global, V_T_r, Sigma, "POD");
	}
	U_global.destroy();
#endif //TSQR


	
#ifdef DMD

    for(size_type i=1; i< TimeSnapsToRead+1; ++i) { 
		size_type snapShot = snapShotStart + snapShotInterval*i;
		if(EntireData){
			size_type fileIDToRead = MPI_PROC_ID%NFiles;
			if(MPI_PROC_ID >= NFiles){
				snapShot += size_type(MPI_PROC_ID/NFiles);
			} 
			std::cout << "EntireData directive: true \n";
			DataVector.push_back(nvector<ComplexType>(readInput<ComplexType>(snapShot, fileIDToRead)));
		} else {
			std::cout << "EntireData directive: false \n";
			DataVector.push_back(nvector<ComplexType>(readInput<ComplexType>(snapShot, PlaneID)));
		}
	}

#ifdef SingleProc
	matrix<ComplexType>	A2 = time_shifted_construct_from_column_vectors(DataVector,0, DataVector.size()-TimeShift);
#else
	matrix<ComplexType>	A2 = construct_from_column_vectors(DataVector,0, DataVector.size());
#endif
	pout << "A2 rows/cols: " << A2.rowCount() << " " << A2.colCount() << "\n";

//	matrix<ComplexType> A2(4,3);
//	ComplexType a2(4,0);
//	ComplexType b2(5,0);
//	ComplexType c2(6,0);
//	ComplexType d2(7,0);
//	ComplexType e2(8,0);
//	ComplexType f2(9,0);
//	ComplexType g2(10,0);
//	ComplexType h2(11,0);
//	ComplexType i2(12,0);
//	ComplexType j2(13,0);
//	ComplexType k2(14,0);
//	ComplexType l2(15,0);
//	A2(0,0) = a2;
//	A2(0,1) = b2;
//	A2(0,2) = c2;
//	A2(1,0) = d2;
//	A2(1,1) = e2;
//	A2(1,2) = f2;
//	A2(2,0) = g2;
//	A2(2,1) = h2;
//	A2(2,2) = i2;
//	A2(3,0) = j2;
//	A2(3,1) = k2;
//	A2(3,2) = l2;
//	pout << "A2 dummy: \n";
//	pout << A2 << "\n";

    for (size_type i = 0; i < DataVector.size(); ++i) {
        DataVector[i].destroy();
    }
    DataVector.clear();
    container<nvector<ComplexType>>(DataVector).swap(DataVector);

//Converting data matrix to complex type 
//	matrix<ComplexType> A2_complex(A2.rowCount(), A2.colCount());
//  constructComplexMatrixFromRealAndImagParts(A2_complex,A2,A2, 0);
//	A2.destroy();

#ifdef Print
	pout << "A2 data matrix: \n";
	pout << A2 << "\n";
#endif

	//Computing Atilde which is rxr projection of matrix A on POD modes
	matrix<ComplexType> UT_x_A2 = transpose(U_local)*A2;

#ifdef Print
	pout << "UT_x_A2_local: " << "\n";	
	pout << UT_x_A2  << "\n";	
#endif
	MPI_Datatype local_UT_x_A2;
	MPI_Type_vector( UT_x_A2.colCount(), UT_x_A2.rowCount(), UT_x_A2.rowCount(), MPI_DOUBLE_COMPLEX, &local_UT_x_A2 );
	MPI_Type_commit( &local_UT_x_A2);

    MPI_Op vector_add;
    MPI_Op_create(add_double_vector, 1, &vector_add);
    MPI_Allreduce(MPI_IN_PLACE, &UT_x_A2(0,0), 1, local_UT_x_A2, vector_add, MPI_COMM_WORLD);
    MPI_Op_free(&vector_add);
	MPI_Type_free( &local_UT_x_A2);
#ifdef Print
	pout << "UT_x_A2_global: " << "\n";	
	pout << UT_x_A2  << "\n";	
#endif

	matrix<ComplexType> UT_x_A2_x_V = UT_x_A2*transpose(V_T_r);
	UT_x_A2.destroy();
#ifdef Print
	pout << "UT_x_A2_x_V: " << "\n";	
	pout << UT_x_A2_x_V  << "\n";	
	pout << "Inverted Sigma: " << "\n";
	pout << invertDiagonalVector(Sigma) << "\n";
#endif
	matrix<ComplexType> Atilde = UT_x_A2_x_V*invertDiagonalVector(Sigma);
	UT_x_A2_x_V.destroy();

#ifdef Print
//	matrix<Type> Atilde(2,2);
//	Atilde(0,0) = 1;
//	Atilde(0,1) = 1;
//	Atilde(1,0) = -1;
//	Atilde(1,1) = 1;
	pout << "Atilde: " << "\n"; //Symmetric	
	pout << Atilde  << "\n";	
#endif

	matrix<ComplexType> L(Atilde.rowCount(), Atilde.rowCount());
	matrix<ComplexType> eV(Atilde.rowCount(), Atilde.rowCount());
	computeEigenValuesOfAComplexMatUsingLAPACK<ComplexType,Type>(Atilde, L, eV);	

	matrix<ComplexType> Omega(Atilde.rowCount(), Atilde.rowCount());

	pout << "Continuous time eigenvalues: \n";
	ComplexType dt(Type(TMAX)/Type(numSnapShots-1),0.0);
	pout << "dt: " << dt << "\n";
	for(unsigned int i=0; i<L.rowCount(); ++i){
		Omega(i,i) = log(L(i,i))/dt;
	}

//	matrix<Type> phi_local = U_local*eV;
	U_local.destroy();
	matrix<ComplexType> factor = A2*transpose(V_T_r)*invertDiagonalVector(Sigma);
	matrix<ComplexType> phi_local = factor*eV;
#ifdef Print
	pout << "L: \n";
	pout << L << "\n";
	pout << "eV: \n";
	pout << eV << "\n";
	pout << "factor: \n";
	pout << factor << "\n";
	pout << "phi_local: factor*eV: \n";
	pout << phi_local << "\n";
#endif

	factor.destroy();
//  L.destroy();
  eV.destroy();
	A2.destroy();
	V_T_r.destroy();
	Sigma.destroy();

  MPI_Datatype local_phiMatType;
  MPI_Datatype global_phiMatType;

  MPI_Type_vector( phi_local.colCount(), phi_local.rowCount(), phi_local.rowCount(), MPI_DOUBLE_COMPLEX, &local_phiMatType );
  MPI_Type_vector( phi_local.colCount(), phi_local.rowCount(), phi_local.rowCount(), MPI_DOUBLE_COMPLEX, &global_phiMatType);

  MPI_Type_commit( &local_phiMatType);
  MPI_Type_commit( &global_phiMatType);

  matrix<ComplexType> phi_global(phi_local.rowCount()*size_type(MPI_PROC_TOTAL_NUM), phi_local.colCount());
  ierr = MPI_Allgather(&phi_local(0,0), 1, local_phiMatType, &phi_global(0,0), 1, global_phiMatType, MPI_COMM_WORLD);
  MPI_Type_free( &local_phiMatType);
  MPI_Type_free( &global_phiMatType);

#ifdef Print
  pout << "phi_global after all gather: \n";
  pout << phi_global << "\n";
#endif


	matrix<ComplexType> U_phi_local(phi_local.rowCount(), phi_local.colCount());
	matrix<ComplexType> V_T_phi(phi_local.colCount(), phi_local.colCount());
	matrix<ComplexType> Sigma_phi(phi_local.colCount(), phi_local.colCount(), 0);

//	matrix<ComplexType> phi_global_test(3,2);
//	ComplexType a(0,1);
//	ComplexType b(0,1);
//	ComplexType c(0,0);
//	ComplexType d(0,1);
//	ComplexType e(0,1);
//	ComplexType f(-1,0);
//	phi_global_test(0,0) = a;
//	phi_global_test(0,1) = b;
//	phi_global_test(1,0) = c;
//	phi_global_test(1,1) = d;
//	phi_global_test(2,0) = e;
//	phi_global_test(2,1) = f;
//	pout << "phi_global_test: \n";
//	pout << phi_global_test << "\n";
//	matrix<ComplexType> U_phi_global(3,3);
//	matrix<ComplexType> V_T_phi(2,2);
//	matrix<ComplexType> Sigma_phi(3,2,0);

//	int infoA1 = computeComplexSVDUsingLAPACK<ComplexType,Type>(phi_global, U_phi_global, Sigma_phi, V_T_phi);

	performTSQR<ComplexType, Type>(phi_local, U_phi_local,  V_T_phi, Sigma_phi);
	phi_local.destroy();


	MPI_Datatype local_UphiMatType;
	MPI_Datatype global_UphiMatType;
	
	MPI_Type_vector( U_phi_local.colCount(), U_phi_local.rowCount(), U_phi_local.rowCount(), MPI_DOUBLE_COMPLEX, &local_UphiMatType );
	MPI_Type_vector( U_phi_local.colCount(), U_phi_local.rowCount(), U_phi_local.rowCount(), MPI_DOUBLE_COMPLEX, &global_UphiMatType);
	
	MPI_Type_commit( &local_UphiMatType);
	MPI_Type_commit( &global_UphiMatType);
	
	matrix<ComplexType> U_phi_global(U_phi_local.rowCount()*size_type(MPI_PROC_TOTAL_NUM), U_phi_local.colCount());
	ierr = MPI_Allgather(&U_phi_local(0,0), 1, local_UphiMatType, &U_phi_global(0,0), 1, global_UphiMatType, MPI_COMM_WORLD);
	MPI_Type_free( &local_UphiMatType);
	MPI_Type_free( &global_UphiMatType);

	U_phi_local.destroy();

#ifdef Print
	pout << "U_phi_global after all gather: \n";
	pout << U_phi_global << "\n";
#endif

	matrix<ComplexType> U_phi_global_star(U_phi_global.colCount(), U_phi_global.rowCount());
	computeComplexConjugate(U_phi_global_star, U_phi_global);
	U_phi_global.destroy();
	matrix<ComplexType> V_T_phi_star(V_T_phi.colCount(), V_T_phi.rowCount());
	computeComplexConjugate(V_T_phi_star, V_T_phi);
	V_T_phi.destroy();
	matrix<ComplexType> Sigma_phi_plus(Sigma_phi.colCount(), Sigma_phi.rowCount(), 0);
	computeSigmaPlus<ComplexType,Type>(Sigma_phi_plus, Sigma_phi);
	Sigma_phi.destroy();


	matrix<ComplexType> Psuedoinverse_phi(V_T_phi_star.rowCount(), U_phi_global_star.colCount());
	Psuedoinverse_phi = V_T_phi_star*Sigma_phi_plus*U_phi_global_star;
#ifdef Print
	pout<< "Psuedoinverse_phi: \n";
	pout << Psuedoinverse_phi << "\n";
#endif

	U_phi_global_star.destroy();
	V_T_phi_star.destroy();
	Sigma_phi_plus.destroy();


	size_type TimeSnapsToReadForFirstVector;
	if(MPI_PROC_TOTAL_NUM==1){
		TimeSnapsToReadForFirstVector = 1 + TimeShift;
	} else {
		TimeSnapsToReadForFirstVector = 1;
	}

    for(size_type i=0; i < TimeSnapsToReadForFirstVector; ++i) { 
		size_type snapShot = snapShotStart + snapShotInterval*i;
		if(EntireData){
			size_type fileIDToRead = MPI_PROC_ID%NFiles;
			std::cout  << MPI_PROC_ID << " fileIDToRead: " << fileIDToRead << " \n";
			if(MPI_PROC_ID >= NFiles){
				snapShot += size_type(MPI_PROC_ID/NFiles);
			} 
			std::cout << "EntireData directive: true \n";
			DataVector.push_back(nvector<ComplexType>(readInput<ComplexType>(snapShot, fileIDToRead)));
		} else {
			std::cout << "EntireData directive: false \n";
			DataVector.push_back(nvector<ComplexType>(readInput<ComplexType>(snapShot, PlaneID)));
		}
	}

#ifdef SingleProc
	matrix<ComplexType>  X0_local = time_shifted_construct_from_column_vectors(DataVector,0, DataVector.size()-TimeShift);
#else
	matrix<ComplexType>	X0_local = construct_from_column_vectors(DataVector,0, DataVector.size());
#endif
	MPI_Datatype local_X0Type;
	MPI_Datatype global_X0Type;
	
	MPI_Type_vector( X0_local.colCount(), X0_local.rowCount(), X0_local.rowCount(), MPI_DOUBLE_COMPLEX, &local_X0Type );
	MPI_Type_vector( X0_local.colCount(), X0_local.rowCount(), X0_local.rowCount(), MPI_DOUBLE_COMPLEX, &global_X0Type);
	
	MPI_Type_commit( &local_X0Type);
	MPI_Type_commit( &global_X0Type);
	
	matrix<ComplexType> X0_global(X0_local.rowCount()*size_type(MPI_PROC_TOTAL_NUM), X0_local.colCount());
	ierr = MPI_Allgather(&X0_local(0,0), 1, local_X0Type, &X0_global(0,0), 1, global_X0Type, MPI_COMM_WORLD);
	MPI_Type_free( &local_X0Type);
	MPI_Type_free( &global_X0Type);
	X0_local.destroy();


//	matrix<ComplexType> firstSnapShotVector(4,1);
//	firstSnapShotVector(0,0) = a;
//	firstSnapShotVector(1,0) = d;
//	firstSnapShotVector(2,0) = g;
//	firstSnapShotVector(3,0) = j;
//	pout << "firstSnapShotVector dummy: \n";
//	pout << firstSnapShotVector << "\n";

    for (size_type i = 0; i < DataVector.size(); ++i) {
        DataVector[i].destroy();
    }
    DataVector.clear();
    container<nvector<ComplexType>>(DataVector).swap(DataVector);

//	matrix<ComplexType>firstSnapShotVector_complex(firstSnapShotVector.rowCount(),firstSnapShotVector.colCount());
//  constructComplexMatrixFromRealAndImagParts(firstSnapShotVector_complex,firstSnapShotVector,firstSnapShotVector, 0);
//	firstSnapShotVector.destroy();
#ifdef Print
	pout<< "X0_global: \n";
	pout << X0_global << "\n";
#endif

	matrix<ComplexType> bVector = Psuedoinverse_phi*X0_global;
	Psuedoinverse_phi.destroy();
	X0_global.destroy();


#ifdef Print
	pout<< "bVector: \n";
	pout << bVector << "\n";
#endif

	matrix<ComplexType> t(size_type(numSnapShots)-1,1);
    for(size_type i=0; i<size_type(numSnapShots)-1; ++i) { 
			t(i,0) = dt*ComplexType(i);
	}
#ifdef Print
	pout<< "time vector: \n";
	pout << t << "\n";
	pout<< "Omega: \n";
	pout << Omega << "\n";
#endif
//	matrix<ComplexType> time_dynamics_T((size_type(numSnapShots)-2),bVector.rowCount());
//
//    for (size_type r = 0; r < time_dynamics_T.rowCount(); ++r) {
//    	for (size_type  k= 0; k < bVector.rowCount(); ++k) {
//				ComplexType val = exp(Omega(k,k)*t(r,0));
//				time_dynamics_T(r,k) = bVector(k,0)*val;
//				pout << bVector(k,0) << " " << Omega(k,k) << " " << t(r,0) << " " << Omega(k,k)*t(r,0) << " " << exp(Omega(k,k)*t(r,0)) << " " << time_dynamics_T(r,k) << "\n";
//		}
//	}
	matrix<ComplexType> time_dynamics(bVector.rowCount(),(size_type(numSnapShots)-2));

    for (size_type c = 0; c < time_dynamics.colCount(); ++c) {
    	for (size_type  r= 0; r < bVector.rowCount(); ++r) {
				ComplexType val = exp(Omega(r,r)*t(c,0));
				time_dynamics(r,c) = bVector(r,0)*val;
				pout << bVector(r,0) << " " << Omega(r,r) << " " << t(r,0) << " " << Omega(r,r)*t(c,0) << " " << exp(Omega(r,r)*t(c,0)) << " " << time_dynamics(r,c) << "\n";
		}
	}
//	matrix<ComplexType> time_dynamics_T=transpose(time_dynamics);
	if(MPI_PROC_ID == 0) {
		PrintModesAndEigValForDMD<ComplexType,Type>(phi_global, time_dynamics, Omega, L, bVector, dt, "DMD");
	}
//	time_dynamics_T.destroy();

  if(MPI_PROC_ID == 0){
     reconstructFromDMDAndWriteToTec<ComplexType, Type>(phi_global, time_dynamics);
	}

#ifdef Print
	pout<< "time_dynamics: \n";
	pout << time_dynamics << "\n";
	pout<< "phi_global*time_dynamics: \n";
	pout << phi_global*time_dynamics << "\n";
#endif

	t.destroy();
	Omega.destroy();
	phi_global.destroy();
	time_dynamics.destroy();
	bVector.destroy();
#endif //DMD


#ifndef DMD
  U_local.destroy();
	V_T_r.destroy();
	Sigma.destroy();
#endif

	MPI_Barrier(MPI_COMM_WORLD);
	pout.close();
	ierr = MPI_Finalize();
	return 0;
}
