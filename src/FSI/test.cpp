#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sstream>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/multi_array.hpp>
using namespace std;

	extern "C"{
        void pass_(const int*, double*, double*, double*,
		   const int*, int*,
                   const int*, int*, double*,
                   const double*, const double*, const double*, const double*,
                   int*, const double*, const double*, const double*, 
                   const int*, double*, double*, double*,
                   double*, double*, double*, 
                   const double*, int*, 
                   const int*, const int*, const int*,
                   const double*, const double*, const double*, const double*,
                   const double*, const double*, const double*, 
                   double*,double*,double*,double*,double*);
	}

   int main(){

	namespace ublas = boost::numeric::ublas;
 
       	const int LociSTREAMNumNodes = 257 ;
       	const int LociSTREAMNumElems = 444 ;
       	const int LociSTREAMNumBC = 54 ;
       	int LociSTREAMIntScheme = 1 ;
       	const int LociSTREAMAnsSize = LociSTREAMNumNodes * 6 ;
//       int LociSTREAMtimeStep = 5;
       	int LociSTREAMRotType = 2 ;
       	int LociSTREAMExtype  = 1 ;
       	int LociSTREAMPlugType = 2 ;

	ublas::matrix<int,ublas::column_major> LociSTREAMConnect;
        LociSTREAMConnect.resize(LociSTREAMNumElems,3) ;

       	const double LociSTREAME1 = 70e+9 ;
       	const double LociSTREAMNU12 = 0.3 ;
       	const double LociSTREAMStructDens= 2700. ;
       	const double LociSTREAMStruThickness = 2.e-4 ;
       	const double LociSTREAMBeta = 0.25 ;
       	const double LociSTREAMGamma = 0.5 ;
       	const double LociSTREAMGenalpha = 0.4;
       	const double LociSTREAMdeltaT= 5.e-6 ;
       	const double LociSTREAMFreq = 30.;
        const double LociSTREAMPitchAmp = 0.0 ;
        const double LociSTREAMFlapAmp = 17.0;
        const double LociSTREAMLagAmp = 0.0 ;
        const double LociSTREAMXfiAmp = 0.0 ;
	const double LociSTREAMYfiAmp = 0.0 ;
	const double LociSTREAMZfiAmp = 0.0 ;
	
       	ublas::vector<double> LociSTREAMXgl(LociSTREAMNumNodes);
       	ublas::vector<double> LociSTREAMYgl(LociSTREAMNumNodes);
       	ublas::vector<double> LociSTREAMZgl(LociSTREAMNumNodes);
//
       	ublas::matrix<double,ublas::column_major> LociSTREAMdisp(LociSTREAMAnsSize,1);
       	ublas::matrix<double,ublas::column_major> LociSTREAMvelo(LociSTREAMAnsSize,1);
       	ublas::matrix<double,ublas::column_major> LociSTREAMacce(LociSTREAMAnsSize,1);
      
	ublas::vector<int> LociSTREAMBCdof(LociSTREAMNumBC);
        ublas::vector<double> LociSTREAMBCval(LociSTREAMNumBC); 
       	ublas::matrix<double,ublas::column_major> LociSTREAMaeroforce(LociSTREAMAnsSize,1);
       	ublas::matrix<double,ublas::column_major> LociSTREAMaeroforce_pre(LociSTREAMAnsSize,1);
        boost::multi_array<double,4>  LociSTREAMNodal_CSYS(boost::extents[3][3][LociSTREAMNumElems][3]);
//        boost::multi_array<double,4>  LociSTREAMNodal_CSYS(boost::extents[3][LociSTREAMNumElems][3][3],boost::fortran_storage_order());
       

 
       	ublas::matrix<double,ublas::column_major>  NLAMSdispRemesh ;
	NLAMSdispRemesh.resize(LociSTREAMNumNodes,3);

        ublas::matrix<double,ublas::column_major> NLAMSoutdisp(LociSTREAMAnsSize,1);
        ublas::matrix<double,ublas::column_major> NLAMSoutvelo(LociSTREAMAnsSize,1);
        ublas::matrix<double,ublas::column_major> NLAMSoutacce(LociSTREAMAnsSize,1);
        boost::multi_array<double,4>  NLAMSNodal_CSYS(boost::extents[3][3][LociSTREAMNumElems][3]);
//        boost::multi_array<double,4>  NLAMSNodal_CSYS(boost::extents[3][LociSTREAMNumElems][3][3],boost::fortran_storage_order());
// read grid file

	ifstream gridfile;
	gridfile.open("trimesh.dat");
        if(!gridfile){
		cerr<<"Error: file could not be opened"<<endl;
		exit(1);
	}
	cout<<" grid file is successfully opened "<< endl;

        int counter = 0 ;
        int tempInt;
        ublas::vector<double> tempVector(3) ;
          
        for(unsigned i=0; i<LociSTREAMXgl.size();++i){
		gridfile >> tempInt >> tempVector(0) >> tempVector(1) >> tempVector(2);
		LociSTREAMXgl(i) = tempVector(0);
                LociSTREAMYgl(i) = tempVector(1);
		LociSTREAMZgl(i) = tempVector(2);
		
//		cout << LociSTREAMXgl(i)<<" "<<LociSTREAMYgl(i)<<" "<<LociSTREAMZgl(i)<<endl;

	}
       
        gridfile.close();

// real conectivity file
	ifstream connectivityfile;
	connectivityfile.open("connect.dat");
	 if(!connectivityfile){
                cerr<<"Error: file could not be opened"<<endl;
                exit(1);
        }
        cout<<" connectivity file is successfully opened "<< endl;

// initialization
        counter = 0;
        ublas::vector<int> tempVectorInt(3) ;

// Create nodes by reading the conectivity information until EOF found
        for(unsigned i=0; i<LociSTREAMConnect.size1();++i){
                connectivityfile >> tempInt >> tempInt >> tempVectorInt(0) >> tempVectorInt(1) >> tempVectorInt(2);
                LociSTREAMConnect(i,0) = tempVectorInt(0);
                LociSTREAMConnect(i,1) = tempVectorInt(1);
                LociSTREAMConnect(i,2) = tempVectorInt(2);

//                cout << i << " " << LociSTREAMConnect(i,0) << " " << LociSTREAMConnect(i,1) << " " << LociSTREAMConnect(i,2)<< endl;
        }

        connectivityfile.close();

// read boundary condition file:
        ifstream bcData ;
        bcData.open("bc.dat") ;
        if (!bcData) {
                cerr << "[E] CSD connectivity file couldn't be opened" << endl;
                exit(1);
        }
        cout << " boundary condition file is successfully opened " << endl;

// initialization
        counter = 0 ;
        double tempReal ;

// Create nodes by reading the mesh until EOF found
         for(unsigned i=0; i<LociSTREAMNumBC; ++i){
		bcData >> tempInt >> tempReal;
		LociSTREAMBCdof(i) = tempInt;
		LociSTREAMBCval(i) = tempReal;
//		cout << LociSTREAMBCdof(i) << " " << LociSTREAMBCval(i) << endl;
	}
		bcData.close();
      

        for(int i=0; i<50000; ++i){
        int LociSTREAMtimeStep = i + 1;
//        cout << "i: "<< i << ", LociSTREMtimeStep: "<< LociSTREAMtimeStep << endl;

        pass_(&LociSTREAMNumNodes, &LociSTREAMXgl(0), &LociSTREAMYgl(0), &LociSTREAMZgl(0),
	      &LociSTREAMNumElems, &LociSTREAMConnect(0,0),
              &LociSTREAMNumBC, &LociSTREAMBCdof(0), &LociSTREAMBCval(0),
              &LociSTREAME1, &LociSTREAMNU12, &LociSTREAMStructDens, &LociSTREAMStruThickness,
              &LociSTREAMIntScheme, &LociSTREAMBeta, &LociSTREAMGamma, &LociSTREAMGenalpha,
              &LociSTREAMAnsSize, &LociSTREAMdisp(0,0), &LociSTREAMvelo(0,0), &LociSTREAMacce(0,0),
              &LociSTREAMaeroforce(0,0), &LociSTREAMaeroforce_pre(0,0),&LociSTREAMNodal_CSYS[0][0][0][0],
              &LociSTREAMdeltaT, &LociSTREAMtimeStep,
              &LociSTREAMRotType, &LociSTREAMExtype, &LociSTREAMPlugType,
              &LociSTREAMFreq, &LociSTREAMPitchAmp, &LociSTREAMFlapAmp, &LociSTREAMLagAmp,
              &LociSTREAMXfiAmp, &LociSTREAMYfiAmp, &LociSTREAMZfiAmp,
              &NLAMSdispRemesh(0,0),&NLAMSoutdisp(0,0),&NLAMSoutvelo(0,0),&NLAMSoutacce(0,0),&NLAMSNodal_CSYS[0][0][0][0]);
      
	
/*          cout<< "in test.cpp" << endl;
          cout<< NLAMSNodal_CSYS[0][0][0][0] << " , " << NLAMSNodal_CSYS[0][0][0][1]<< " , " << NLAMSNodal_CSYS[0][0][0][2] << endl;
          cout<< NLAMSNodal_CSYS[0][0][1][0] << " , " << NLAMSNodal_CSYS[0][0][1][1]<< " , " << NLAMSNodal_CSYS[0][0][1][2]<< endl;
          cout<< NLAMSNodal_CSYS[0][0][2][0] << " , " << NLAMSNodal_CSYS[0][0][2][1]<< " , " << NLAMSNodal_CSYS[0][0][2][2] << endl;
*/
              LociSTREAMdisp = NLAMSoutdisp;
              LociSTREAMvelo = NLAMSoutvelo;
              LociSTREAMacce = NLAMSoutacce;
              LociSTREAMNodal_CSYS = NLAMSNodal_CSYS;
       }
       return 0;
   }
