
#ifndef __PROCESSOR_HPP__  
#define __PROCESSOR_HPP__


#include "mpi.h"
// alias template for a container type

    int MPI_PROC_TOTAL_NUM, MPI_PROC_ID;

	void defineMPIVarsAndProcFiles(int argc, char *argv[] )	{
    	int ierr = MPI_Init(&argc, &argv);
    	ierr = MPI_Comm_size(MPI_COMM_WORLD, &MPI_PROC_TOTAL_NUM);
    	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &MPI_PROC_ID);
    	std::cout << "Processor " << MPI_PROC_ID << " is active.\n";
	    pout.open((ProcessorPrintOutPath(MPI_PROC_ID)).c_str());
		pout << "MPI_PROC_TOTAL_NUM: " << MPI_PROC_TOTAL_NUM << "\n";
    	pout << "Processor " << MPI_PROC_ID << " is active.\n";
	}

	void add_double_vector(void *in, void *inout, int *len, MPI_Datatype *dtype)
	{
//		double *invec = (double*) in;
//		double *inoutvec = (double*) inout;
		std::complex<double> *invec = (std::complex<double>*) in;
		std::complex<double> *inoutvec = (std::complex<double>*) inout;
		int nints, naddresses, ntypes;
		int combiner;
		if (*len != 1) {
		    throw "my_add: len>1 not implemented.";
		    return;
		} 
		MPI_Type_get_envelope(*dtype, &nints, &naddresses, &ntypes, &combiner); 
		if (combiner != MPI_COMBINER_VECTOR) {
		    throw "my_add: do not understand composite datatype.";
		    return;
		} 
		int vecargs [nints];
		MPI_Aint vecaddrs[naddresses];
		MPI_Datatype vectypes[ntypes];
		
		MPI_Type_get_contents(*dtype, nints, naddresses, ntypes, 
		        vecargs, vecaddrs, vectypes);
		
//		if (vectypes[0] != MPI_DOUBLE) {
//		    throw "my_add: not a vector of DOUBLEs.";
//		}
		
		int count    = vecargs[0];
		int blocklen = vecargs[1];
		int stride   = vecargs[2];
		
		for ( int i=0; i<count; i++ ) {
		    for ( int j=0; j<blocklen; j++) {
		        inoutvec[i*stride+j] += invec[i*stride+j]; 
		    } 
		}
	}

//	MPI_Datatype RMatType;
//
//	struct MPI_RMat {
//		matrix<T> rmat;
//	};
//
//	void setMPIVectorTypeForRMatrix(size_type r_size) {
//		MPI_RMat v_rmat;
//		MPI_Datatype type[1] = {MPI_DOUBLE};
//		int blocklen[1] = {r_size};
//		MPI_Aint typedisp[1];
//		MPI_Aint startaddress, address;
//
//		typedisp[0] = 0;
//		MPI_Get_address(&v_rmat.rmat, &startaddress);
//		MPI_Type_create_struct(1, blocklen, typedisp, type, &RMatType);
//		MPI_Type_commit(&RMatType);
//
//	}

#endif 



