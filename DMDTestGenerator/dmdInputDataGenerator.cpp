#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <stdexcept>

#include <math.h>
#include <complex>
#include <iomanip>
#include <sys/stat.h>

//#define SnapShots 200
//#define DataPoints 400
//#define NFiles 1
//#define TMIN 0.0
//#define TMAX 4*M_PI
//#define XMIN -10.0
//#define XMAX 10.0
//#define NFreqs 2
#define SnapShots 10
#define DataPoints 10
#define NFiles 2
#define TMIN 0.0
#define TMAX 20.0
#define XMIN 0.0
#define XMAX 80.0
#define NFreqs 1
#define visc 0.075
using size_type = unsigned long int;


void generateFolders(){
    using std::string;
	string f1= "DMD_Input";
    const int dir_err = mkdir(f1.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

	for(size_type n=0; n<size_type(SnapShots); ++n) {
		string f2=static_cast<std::ostringstream*>( &(std::ostringstream() << n) )->str();
		string ff = f1+"/"+f2;
    	const int dir2_err = mkdir(ff.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	} 
}
template<typename T, typename MyFunc>
void writeToTec(MyFunc func){

	T DataPointsPerFile = T(DataPoints)/NFiles;
	T TimeDurationPerSnap = (T(TMAX)-T(TMIN))/SnapShots;
	T ResolutionPerDataPoint = (T(XMAX)-T(XMIN))/DataPoints;
	std::cout << "ResolutionPerDataPoint: " << ResolutionPerDataPoint << "\n"; 
	using std::string;
	string f1="Realfreq_";
	for(size_type n=0; n<size_type(NFreqs)+1; ++n){

		string f2=static_cast<std::ostringstream*>( &(std::ostringstream() << n) )->str();
		if(n==NFreqs){
			f2="composed";
		}
		string ff=f1+f2+".tec";
		std::ofstream data(ff.c_str());
		data << "Title=" << f1+f2 << "\n";
		data << "VARIABLES= \"x\", \"fr_"<< f2 << "(x,t)\"\n";
		for(size_type t=0; t<SnapShots; ++t) {
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
				data << std::setw(15) << std::real(func(spat,time,n));
				counter++;
				if(counter%3==0) {data << "\n";}
			}
			data << "\n";
		}
		data.close();
	}

	f1="Imagfreq_";
	for(size_type n=0; n<size_type(NFreqs)+1; ++n){

		string f2=static_cast<std::ostringstream*>( &(std::ostringstream() << n) )->str();
		if(n==NFreqs){
			f2="composed";
		}
		string ff=f1+f2+".tec";
		std::ofstream data(ff.c_str());
		data << "Title=" << f1+f2 << "\n";
		data << "VARIABLES= \"x\", \"fi_"<< f2 << "(x,t)\"\n";
		for(size_type t=0; t<SnapShots; ++t) {
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
				data << std::setw(15) << std::imag(func(spat,time,n));
				counter++;
				if(counter%3==0) {data << "\n";}
			}
			data << "\n";
		}
		data.close();
	}

}

template<typename T>
void generateAnalyticalModes() {


		T DataPointsPerFile = T(DataPoints)/NFiles;
		T ResolutionPerDataPoint = (T(XMAX)-T(XMIN))/DataPoints;
		T TimeDurationPerSnap = (T(TMAX)-T(TMIN))/SnapShots;
	  T alpha = 1.0/(4.0*T(visc))*(sqrt(2+2*sqrt(1+16*pow(T(visc),2)))-2);
	  T beta = 1.0/(4.0*T(visc))*sqrt(-2+2*sqrt(1+16*pow(T(visc),2)));

    auto phi_1 = [](T x, T alpha, T beta) {
    		T phi_1; 

    		phi_1 = sqrt(2*alpha)/beta*exp(-alpha*x)*(beta*cos(beta*x) +
          (-alpha + sqrt(pow(alpha,2)+pow(beta,2))*sin(beta*x)));

        return phi_1;
    };
    auto phi_2 = [](T x, T alpha, T beta) {
    		T phi_2; 

    		phi_2 = -sqrt(2*alpha)/beta*exp(-alpha*x)*(beta*cos(beta*x)+
          (-alpha - sqrt(pow(alpha,2)+pow(beta,2))*sin(beta*x)));

        return phi_2;
    };
    auto at_1 = [](T t, T alpha, T beta) {
    		T at_1; 

    		at_1 = sqrt(2)/(4*sqrt(alpha)*sqrt(pow(alpha,2)+pow(beta,2)))*((alpha+
          sqrt(pow(alpha,2)+pow(beta,2)))*cos(t)+beta*sin(t));

        return at_1;
    };

    auto at_2 = [](T t, T alpha, T beta) {
    		T at_2; 

    		at_2 = -sqrt(2)/(4*sqrt(alpha)*sqrt(-1*pow(alpha,2)+pow(beta,2)))*((alpha+
          sqrt(pow(alpha,2)+pow(beta,2)))*cos(t)-beta*sin(t));

        return at_2;
    };

	  std::ofstream analyticalSModes("AnalyticalSpatialModesOzgurTest.dat");
    analyticalSModes << "title = \"SUGAR TECPLOT\"" << '\n';
    analyticalSModes << "variables = \"x\", \"Phi1\", \"Phi2\"" << '\n';
    analyticalSModes << '\n';
		for(size_type d=0; d<DataPoints; ++d) {
				T spat=T(XMIN)+d*ResolutionPerDataPoint;
        analyticalSModes << std::setw(15) << spat  << std::setw(15) << phi_1(spat,alpha, beta) << std::setw(15) << phi_2(spat,alpha, beta)  << "\n";
    }
    analyticalSModes.close();

  	std::ofstream analyticalTModes("AnalyticalTemporalModesOzgurTest.dat");
    analyticalTModes << "title = \"SUGAR TECPLOT\"" << '\n';
    analyticalTModes << "variables = \"y\", \"at1\", \"at2\"" << '\n';
    analyticalTModes << '\n';
		for(size_type t=0; t<SnapShots; ++t) {
				T time=T(TMIN)+t*TimeDurationPerSnap;
        analyticalTModes << std::setw(15) << time  << std::setw(15) << at_1(time,alpha, beta)  << std::setw(15) << at_2(time,alpha, beta) << "\n";
    }
    analyticalTModes.close();
}

template<typename T, typename MyFunc>
void generateDataFiles(MyFunc func){

	T DataPointsPerFile = T(DataPoints)/NFiles;
	T TimeDurationPerSnap = (T(TMAX)-T(TMIN))/SnapShots;
	T ResolutionPerDataPoint = (T(XMAX)-T(XMIN))/DataPoints;
	using std::string;
	string f1= "DMD_Input/";
	for(size_type t=0; t<SnapShots; ++t) {
		string f2=static_cast<std::ostringstream*>( &(std::ostringstream() << t) )->str();
		for(size_type f=0; f<NFiles; ++f) {
			string f3=static_cast<std::ostringstream*>( &(std::ostringstream() << f) )->str();
			string ff=f1+f2+"/SampledInstaDump_X_"+f3+".dat";
			std::ofstream data(ff.c_str());
			for(size_type d=0; d<DataPointsPerFile; ++d) {
                				
				T spat=T(XMIN)+(T(d)+T(f)*T(DataPointsPerFile))*ResolutionPerDataPoint;
				T time=T(TMIN)+(T(t))*TimeDurationPerSnap;
				data << std::setw(15) << std::real(func(spat,time,NFreqs));
				
			}
			data.close();
		}
	}
}


int main() {
using std::cout;
using Type = double;

    auto func = [](Type x, Type t, size_type k) {
    
    std::complex<Type> arg[NFreqs]; 
    std::complex<Type> freq[NFreqs+1]; 
#ifdef KutzTest
//    std::complex<Type> omega1_t(-0.1*t,t*2.3); 
//    std::complex<Type> omega2_t(-0.2*t,t*2.8); 
    std::complex<Type> omega1_t(0.0,t*2.3); 
    std::complex<Type> omega2_t(0.0,t*2.8); 
    arg[0] = omega1_t;
    arg[1] = omega2_t;
		cout << "arg[0] : " << arg[0] << "\n";
		cout << "arg[1] : " << arg[1] << "\n";

		freq[0] = (1.0/std::cosh(x+3))*exp(arg[0]);
		freq[1] = (2*std::tanh(x)/std::cosh(x))*exp(arg[1]);

#elif OzgurTest

	  Type alpha = 1.0/(4.0*Type(visc))*(sqrt(2+2*sqrt(1+16*pow(Type(visc),2)))-2);
	  Type beta = 1.0/(4.0*Type(visc))*sqrt(-2+2*sqrt(1+16*pow(Type(visc),2)));
//	  std::cout << "alpha: " << alpha << "\n";
//	  std::cout << "beta: " << beta << "\n";

    freq[0] = exp(-alpha*x)*cos(beta*x-t);
#endif


		std::complex<Type> dummy(0.0,0.0);
		freq[NFreqs]=dummy;
		for(size_type n=0; n<NFreqs; ++n){
			freq[NFreqs] += freq[n];
		}
        return freq[k];
    };

	generateFolders();
	writeToTec<Type>(func);
	generateDataFiles<Type>(func);
#ifdef OzgurTest
	generateAnalyticalModes<Type>();
#endif
	return 0;
}
