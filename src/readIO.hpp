#ifndef __READIOHPP__ 
#define __READIOHPP__ 

#include "input.hpp"
#include "common.hpp"
#include "matrix.hpp"

#include <iostream>
#include <string> 
#include <sstream>
#include <algorithm>

template<typename T>
void split(const std::string &s, char delim, std::vector<T> &elems) {
//	std::cout << s << "\n";
    std::stringstream ss(s);
    std::string item;
//	pout << "reading: " << std::getline(ss,item,delim) << "\n";
	int counter=0;
    while (ss >> item) {
//    while (std::getline(ss, item, delim)) {
        if(RealInput) {
					T read(std::stod(item),0.0);
        	elems.push_back(read);
				} 
		++counter;
    }
	std::cout << "elements read per line: " << counter << "\n";
}

template<typename T>
std::vector<T> split(const std::string &s, char delim) {
    std::vector<T> elems;
    split(s, delim, elems);
    return elems;
}

std::string getInputFileName(int snapShot, int planeID){
    using std::string;
#ifdef OzgurN2Case
    string iFd  ="../InputFiles_N2/NOTvib_180_670K_500Steps/SplitData/Tvib_NO_Split_";
    string pID =static_cast<std::ostringstream*>( &(std::ostringstream() << planeID) )->str();
    string fE=".dat";
    string fname = iFd + pID + fE;
#else
    string iFd  ="../InputFiles/";
    string sFd =static_cast<std::ostringstream*>( &(std::ostringstream() << snapShot) )->str();
    string fS  ="/SampledInstaDump_";
	string under = "_";
    string pID =static_cast<std::ostringstream*>( &(std::ostringstream() << planeID) )->str();
    string ReadVarString =static_cast<std::ostringstream*>( &(std::ostringstream() << ReadVar) )->str();
    string fE=".dat";
    string fname = iFd + sFd + fS + ReadVarString + under + pID + fE;
#endif
   return fname;
}
template<typename T>
std::vector<T> readInput(int snapShot, int planeID) {
	std::vector<T> data;
	std::string fileName = getInputFileName(snapShot, planeID);
	std::ifstream infile(fileName.c_str(), std::ios_base::in);
	std::cout << "reading: " << fileName << "\n";
	if(infile.fail()){
	    pout << "Problem reading file! \n";
        throw std::invalid_argument("Problem reading file in readInput");
	}else{
	    std::stringstream buffer;
	    buffer << infile.rdbuf();
	    infile.close();
	   	data = split<T>(buffer.str(),'\t');
	}
//	std::cout << "snapShot: " << snapShot << " size: " << data.size() << "\n";
	std::cout << "snapShot: " << snapShot << " size: " << data.size() << "\n";
//	std::cout << "Data: " << data << "\n";
	return data;
}


std::string ProcessorPrintOutPath(int procID){
    using std::string;
	string fID = "ProcFiles/Processor_";
	string pID =static_cast<std::ostringstream*>( &(std::ostringstream() << procID) )->str();
    std::string OPString =static_cast<std::ostringstream*>( &(std::ostringstream() << OPFolder) )->str();
	string fE = ".out";
	string path = OPString + fID+pID+fE;
	return path;
}

template<typename CT>
void PrintModesAndEigVal(matrix<CT>& U, matrix<CT>& V_T, matrix<CT>& Sigma, std::string type){

    using size_type = typename matrix<CT>::size_type;

    for(size_type k=0; k<std::min(size_type(SpatialModesToOP),U.colCount()); ++k) {
        std::string ModeNo =static_cast<std::ostringstream*>( &(std::ostringstream() << k) )->str();
        std::string OPString =static_cast<std::ostringstream*>( &(std::ostringstream() << OPFolder) )->str();
        std::string filename = OPString + "SpatialMode_" + type + "_"+ ModeNo + ".dat";
        std::ofstream spatialModes(filename.c_str());
        size_type counter=1;
		counter=1;
        for(size_type i=0; i<U.rowCount(); ++i) {
            spatialModes << std::setw(20) << std::real(U(i,k));
            if(counter%3 == 0){
                spatialModes << "\n";
            } 
            counter++;
        } 
        spatialModes.close();
    }

    std::string OPString =static_cast<std::ostringstream*>( &(std::ostringstream() << OPFolder) )->str();
    std::string filename = OPString + "TemporalModes_" + type + ".dat";
    std::ofstream temporalModes(filename.c_str());

	for(size_type i=0; i<V_T.colCount(); ++i) {
    	temporalModes << std::setw(15) << i;
    	for(int k=0; k<std::min(size_type(TemporalModesToOP),V_T.rowCount()); ++k) {
			temporalModes << std::setw(15) << std::real(V_T(k,i));
    	}
		temporalModes << "\n";
	}
    temporalModes.close();
	
	  matrix<CT> s_V_T = Sigma*V_T;

    filename = OPString + "s_TemporalModes_" + type + ".dat";
    std::ofstream s_temporalModes(filename.c_str());

	for(size_type i=0; i<s_V_T.colCount(); ++i) {
    	s_temporalModes << std::setw(15) << i;
    	for(int k=0; k<std::min(size_type(TemporalModesToOP),s_V_T.rowCount()); ++k) {
			s_temporalModes << std::setw(15) << std::real(s_V_T(k,i));
    	}
		s_temporalModes << "\n";
	}
    s_temporalModes.close();
	s_V_T.destroy();


    filename = OPString + "SingularValues_" + type + ".dat";
    std::ofstream eigenValues(filename.c_str());
		CT sum(0.0,0.0);
    for(size_type i=0; i<Sigma.colCount(); ++i) {
       sum +=  pow(Sigma(i,i),2.0);
    }
	std::cout << "Total Sigma sum: " << sum << "\n";
    for(size_type i=0; i<Sigma.colCount(); ++i) {
//        eigenValues << std::setw(20) << Sigma(i,i)/sum<< "\n";
        eigenValues << std::setw(20) << std::real(Sigma(i,i)) << std::setw(20) << std::real(pow(Sigma(i,i),2.0)/sum)<< "\n";
    }
    eigenValues.close();

}

template<typename CT, typename T>
void PrintModesAndEigValForDMD(matrix<CT>& U, matrix<CT>& V_T, matrix<CT>& Omega, matrix<CT>& L, matrix<CT>& bVector, CT dt, std::string type){

    using size_type = typename matrix<CT>::size_type;
    std::string OPString =static_cast<std::ostringstream*>( &(std::ostringstream() << OPFolder) )->str();


    for(size_type k=0; k<std::min(size_type(SpatialModesToOP),U.colCount()); ++k) {
        std::string ModeNo =static_cast<std::ostringstream*>( &(std::ostringstream() << k) )->str();
        std::string filename = OPString  + "SpatialMode_" + type + "_Real_"+ ModeNo + ".dat";
        std::ofstream spatialModes(filename.c_str());
        size_type counter=1;
        for(size_type i=0; i<U.rowCount(); ++i) {
            spatialModes << std::setw(20) << std::real(U(i,k));
            if(counter%3 == 0){
                spatialModes << "\n";
            } 
            counter++;
        } 
        spatialModes.close();
    }
    for(size_type k=0; k<std::min(size_type(SpatialModesToOP),U.colCount()); ++k) {
        std::string ModeNo =static_cast<std::ostringstream*>( &(std::ostringstream() << k) )->str();
        std::string filename = OPString  + "SpatialMode_" + type + "_Imag_"+ ModeNo + ".dat";
        std::ofstream spatialModes(filename.c_str());
        size_type counter=1;
        for(size_type i=0; i<U.rowCount(); ++i) {
            spatialModes << std::setw(20) << std::imag(U(i,k));
            if(counter%3 == 0){
                spatialModes << "\n";
            } 
            counter++;
        } 
        spatialModes.close();
    }

    std::string filename = OPString  + "TemporalModes_Real_" + type + ".dat";
    std::ofstream temporalModesR(filename.c_str());
	for(size_type i=0; i<V_T.colCount(); ++i) {
    	temporalModesR << std::setw(15) << i;
    	for(int k=0; k<std::min(size_type(TemporalModesToOP),V_T.rowCount()); ++k) {
				temporalModesR << std::setw(15) << std::real(V_T(k,i));
    	}
		temporalModesR << "\n";
	}
    temporalModesR.close();

    filename = OPString  + "TemporalModes_Imag_" + type + ".dat";
    std::ofstream temporalModesI(filename.c_str());
	for(size_type i=0; i<V_T.colCount(); ++i) {
    	temporalModesI << std::setw(15) << i;
    	for(int k=0; k<std::min(size_type(TemporalModesToOP),V_T.rowCount()); ++k) {
				temporalModesI << std::setw(15) << std::imag(V_T(k,i));
    	}
		temporalModesI << "\n";
	}
    temporalModesI.close();

    filename = OPString + "EigenValues_" + type + ".dat";
    std::ofstream eigenValues(filename.c_str());
    for(size_type i=0; i< size_type(RANK); ++i) {
        eigenValues << std::setw(20) << std::real(L(i,i))
                    << std::setw(20) << std::imag(L(i,i))
                    << std::setw(20) << std::real(Omega(i,i)) 
										<< std::setw(20) << std::imag(Omega(i,i))
//		                << std::setw(20) << std::real(exp(Sigma(i,i))) 
//							      << std::setw(20) << std::imag(exp(Sigma(i,i))) 
							      << std::setw(20) <<  std::abs(bVector(i,0)) 
									  << std::setw(20) << std::arg(L(i,i))/(2*M_PI*std::real(dt))<< "\n";
    }
    eigenValues.close();

    filename = OPString + "UnitCircle.dat";
    std::ofstream unitCircle(filename.c_str());
		T dTheta=(2*M_PI-0)/200.0;
    for(size_type i=0; i<200; ++i) {
				T theta = i*dTheta;
			  CT arg(0,theta);
        unitCircle << std::setw(20) << std::real(exp(arg)) << std::setw(20) << std::imag(exp(arg)) << "\n";
    }
    unitCircle.close();

}

template<typename T>
void PrintModesAndEigValForCovariancePOD(matrix<T>& Phi, matrix<T>& aT, matrix<T>& Lambda, matrix<T>& RelE, std::string type){

    using size_type = typename matrix<T>::size_type;
    std::string OPString =static_cast<std::ostringstream*>( &(std::ostringstream() << OPFolder) )->str();


    for(size_type k=0; k<std::min(size_type(SpatialModesToOP),Phi.colCount()); ++k) {
        std::string ModeNo =static_cast<std::ostringstream*>( &(std::ostringstream() << k) )->str();
        std::string filename = OPString + "SpatialMode_" + type +"_"+ ModeNo + ".dat";
        std::ofstream spatialModes(filename.c_str());
        size_type counter=1;
        counter=1;
        for(size_type i=0; i<Phi.rowCount(); ++i) {
            spatialModes << std::setw(20) << Phi(i,k);
            if(counter%3 == 0){
                spatialModes << "\n";
            } 
            counter++;
        } 
        spatialModes.close();
    }

    std::string filename = OPString + "TemporalModes_" + type + ".dat";
    std::ofstream temporalModes(filename.c_str());

	for(size_type i=0; i<aT.colCount(); ++i) {
    	temporalModes << std::setw(15) << i;
    	for(int k=0; k<std::min(size_type(TemporalModesToOP),aT.rowCount()); ++k) {
			temporalModes << std::setw(15) << aT(k,i);
    	}
		temporalModes << "\n";
	}
    temporalModes.close();


    filename = OPString  + "EigenValues_" + type + ".dat";
    std::ofstream eigenValues(filename.c_str());
    for(size_type i=0; i<RelE.colCount(); ++i) {
        eigenValues << std::setw(20) << Lambda(i,i) << std::setw(20) << RelE(i,i)<< "\n";
    }
    eigenValues.close();

}

template<typename CT, typename T>
void reconstructFromDMDAndWriteToTec(matrix<CT>& U, matrix<CT>& V_T){
    using size_type = typename matrix<T>::size_type;
    std::string OPString =static_cast<std::ostringstream*>( &(std::ostringstream() << OPFolder) )->str();

	matrix<CT> f_DMD = U*V_T;
#define TMIN 0.0
#define XMIN 0.0
#define XMAX 80.0

	T DataPoints = U.rowCount()/2.0-1;
  T TimeDurationPerSnap = (T(TMAX)-T(TMIN))/numSnapShots;
  T ResolutionPerDataPoint = (T(XMAX)-T(XMIN))/(DataPoints);

  using std::string;
  string ff= OPString + "Realfreq_Reconstructed.tec";
  std::ofstream data(ff.c_str());
  data << "Title = f_reconstructed_dmd_real\n";
  data << "VARIABLES= \"x\", \"f_DMD(x,t)\"\n";
  for(size_type t=0; t<numSnapShots; ++t) {
    data << "ZONE F=BLOCK, T=\"" << t*TimeDurationPerSnap << "\", I=" << DataPoints << "\n";
    int counter=0;
    for(size_type d=0; d<DataPoints; ++d) {
      data << std::setw(15) << T(XMIN)+d*ResolutionPerDataPoint;
      counter++;
      if(counter%3==0) {data << "\n";}
    }
    data << "\n";
    for(size_type d=0; d<DataPoints; ++d) {
      T spat=T(XMIN)+d*ResolutionPerDataPoint;
      T time=T(TMIN)+t*TimeDurationPerSnap;
      data << std::setw(15) << std::real(f_DMD(d,t));
      counter++;
      if(counter%3==0) {data << "\n";}
    }
    data << "\n";
  }
  data.close();
  ff= OPString + "Imagfreq_Reconstructed.tec";
  std::ofstream dataI(ff.c_str());
  dataI << "Title = f_reconstructed_dmd_real\n";
  dataI << "VARIABLES= \"x\", \"f_DMD_"<< f_DMD << "(x,t)\"\n";
  for(size_type t=0; t<numSnapShots; ++t) {
    dataI << "ZONE F=BLOCK, T=\"" << t*TimeDurationPerSnap << "\", I=" << DataPoints << "\n";
    int counter=0;
    for(size_type d=0; d<U.rowCount(); ++d) {
      dataI << std::setw(15) << T(XMIN)+d*ResolutionPerDataPoint;
      counter++;
      if(counter%3==0) {data << "\n";}
    }
    dataI << "\n";
    for(size_type d=0; d<DataPoints; ++d) {
      T spat=T(XMIN)+d*ResolutionPerDataPoint;
      T time=T(TMIN)+t*TimeDurationPerSnap;
      dataI << std::setw(15) << std::real(f_DMD(d,t));
      counter++;
      if(counter%3==0) {data << "\n";}
    }
    dataI << "\n";
  }
  dataI.close();
	f_DMD.destroy();
}


#endif
