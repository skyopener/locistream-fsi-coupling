
#include "flamelet_table.h"
#include "fstream"
#include <iomanip>
#include <iostream>
#include <string>
#include <stdlib.h>
using namespace std;

	
// Loci includes.
#include <Loci.h>

void Flamelet_Table::read(string fileName) {
	string tag;
	int nVar;
 
	if (Loci::MPI_rank==0) cout << "[I] Reading the flamelet table" << endl;
	
	//Open file
	ifstream chemTable;
	chemTable.open(fileName.c_str(), ios::in);
	if (!chemTable.is_open()) {
		if (Loci::MPI_rank==0) cerr << "[E] Flamelet table file " << fileName  << " couldn't be opened" << endl;
		exit(1);
	}

	chemTable >> nZ; if (Loci::MPI_rank==0) cout << "[I] ...Number of Z points: " << nZ << endl;
	chemTable >> nS; if (Loci::MPI_rank==0) cout << "[I] ...Number of S points: " << nS << endl;
	chemTable >> nChi; if (Loci::MPI_rank==0) cout << "[I] ...Number of Chi points: " << nChi << endl;
	chemTable >> nVar; if (Loci::MPI_rank==0) cout << "[I] ...Number of variables: " << nVar << endl;
	
	Z.resize(nZ);
	S.resize(nS);
	Chi.resize(nChi);
	weights.resize(6);

	// A fix if S has only one entry 
	bool Sflag=false;
 	if (S.size()==1) {
		Sflag=true;
		S.resize(2);
		nS=2;
	}

	rho.resize(nZ);
	T.resize(nZ);
	laminarViscosity.resize(nZ);
	diffusivity.resize(nZ);
	Fchi.resize(nZ);
	//Y.resize(nZ);
	for (int i=0;i<nZ;++i) {
		rho[i].resize(nS);
		T[i].resize(nS);
		laminarViscosity[i].resize(nS);
		diffusivity[i].resize(nS);
		Fchi[i].resize(nS);
		//Y[i].resize(nS);
		for (int j=0;j<nS;++j) {
			rho[i][j].resize(nChi);
			T[i][j].resize(nChi);
			laminarViscosity[i][j].resize(nChi);
			diffusivity[i][j].resize(nChi);
			Fchi[i][j].resize(nChi);
			//Y[i][j].resize(nChi);
		}
	}

	if (Sflag) nS=1;	
	
	for (int i=0;i<nZ;++i) chemTable >> Z[i];
	for (int i=0;i<nS;++i) chemTable >> S[i];
	for (int i=0;i<nChi;++i) chemTable >> Chi[i];
	
	double dummy;
	
	for (int i=0;i<nVar;++i){
		chemTable >> tag;
		if (Loci::MPI_rank==0) cout << "[I] ...Found variable: " << tag << endl;
		if (tag=="RHO") {
			for (int j=0;j<nZ;++j) for (int k=0;k<nS;++k) for (int m=0;m<nChi;++m) chemTable >> rho[j][k][m];
		} else if (tag=="T") {
			for (int j=0;j<nZ;++j) for (int k=0;k<nS;++k) for (int m=0;m<nChi;++m) chemTable >> T[j][k][m]; 
		} else if (tag=="MU") {
			for (int j=0;j<nZ;++j) for (int k=0;k<nS;++k) for (int m=0;m<nChi;++m) chemTable >> laminarViscosity[j][k][m]; 
		} else if (tag=="DIFF") {
			for (int j=0;j<nZ;++j) for (int k=0;k<nS;++k) for (int m=0;m<nChi;++m) chemTable >> diffusivity[j][k][m];
		} else if (tag=="F(Z)-ChiSt") {
			for (int j=0;j<nZ;++j) for (int k=0;k<nS;++k) for (int m=0;m<nChi;++m) chemTable >> Fchi[j][k][m]; // This one isn't actually a function of Chi
		} else {
			if (Loci::MPI_rank==0) cout << "[I] ...Skipping" << endl;
			for (int j=0;j<nZ;++j) for (int k=0;k<nS;++k) for (int m=0;m<nChi;++m) chemTable >> dummy;
		}		
      	}
	
	if (Sflag) {
		for (int j=0;j<nZ;++j) for (int m=0;m<nChi;++m) {
			rho[j][1][m]=rho[j][0][m];
			T[j][1][m]=T[j][0][m];
			laminarViscosity[j][1][m]=laminarViscosity[j][0][m];
			diffusivity[j][1][m]=diffusivity[j][0][m];
			Fchi[j][1][m]=Fchi[j][0][m];
		}
		nS=2;
	}





	return; 
} // end Flamelet_Table::read

void Flamelet_Table::get_weights(const double &Z_in, const double &Zvar_in, const double &Chi_in) {
	
	// Convert Zvar_in to S_in
	double S_in;
	
	if ( fabs(Z_in)<1.e-8 || fabs(Z_in-1.)<1.e-8) {
		S_in=0.;
	} else {
		S_in=Zvar_in/(Z_in*(1.-Z_in));
	}
		
	i1=-1;
	if (Z_in<=Z[0]) {
		i1=0;
		weights[0]=1.0;
	} else if (Z_in>=Z[nZ-1]) {
		i1=nZ-2;
		weights[0]=0.0;
	} else {
		for (int i=0;i<nZ-1;++i) {
			if (Z_in<=Z[i+1]) { 
				i1=i; break;
			}
		}
		weights[0]=1.0-(Z_in-Z[i1])/(Z[i1+1]-Z[i1]);	
	}
	weights[1]=1.0-weights[0];
	
	i2=-1;
	if (S_in<=S[0]) {
		i2=0;
		weights[2]=1.0;
	} else if (S_in>=S[nS-1]) {
		i2=nS-2;
		weights[2]=0.0;
	} else {
		for (int i=0;i<nS-1;++i) {
			if (S_in<=S[i+1]) { 
				i2=i; break;
			}
		}
		weights[2]=1.0-(S_in-S[i2])/(S[i2+1]-S[i2]);	
	}

	if (S.size()==2) weights[2]=1.; i2=0;
	weights[3]=1.0-weights[2];

	// Now first convert Chi_in to Chi_st
	// Fchi doesn't depend on Chi, so
	i3=0;
	double F=weights[2]*( weights[0]*Fchi[i1][i2][i3]+weights[1]*Fchi[i1+1][i2][i3] )
			+weights[3]*( weights[0]*Fchi[i1][i2+1][i3]
			+weights[1]*Fchi[i1+1][i2+1][i3] );
	F=max(1.e-8,F); // Clip
	// Note that Chi_in is actually log10(Chi_in).
	// So the conversion goes like:
	double Chi_st=Chi_in-log10(F); 

	//*******************************************************
	// Comment out the following to turn-of non-eq. effects
	// i.e. use the smallest chi (need a dirac pdf table)
 	// Chi_st=Chi[0];
	//*******************************************************

	i3=-1;
	if (Chi_st<=Chi[0]) {
		i3=0;
		weights[4]=1.0;
	} else if (Chi_st>=Chi[nChi-1]) {
		i3=nChi-2;
		weights[4]=0.0;
	} else {
		for (int i=0;i<nChi-1;++i) {
			if (Chi_st<=Chi[i+1]) {
				i3=i; break;
			}
		}
		weights[4]=1.0-(Chi_st-Chi[i3])/(Chi[i3+1]-Chi[i3]);	
	}

	weights[5]=1.0-weights[4];

	// A little check
	if ((i1+1)*(i2+1)*(i3+1)==0){
		cout << "Error in Lookup table" <<endl;
		exit(1);
	}
	
	return;
}

double Flamelet_Table::get_rho(const double &Z_in, const double &Zvar_in, const double &Chi_in,bool refreshWeights) {
	
	if (refreshWeights) get_weights(Z_in,Zvar_in,Chi_in);
	
	return weights[4]*( weights[2]*( weights[0]*rho[i1][i2][i3]+weights[1]*rho[i1+1][i2][i3] )
			+weights[3]*( weights[0]*rho[i1][i2+1][i3]
			+weights[1]*rho[i1+1][i2+1][i3] ) ) 
			+weights[5]*( weights[2]*( weights[0]*rho[i1][i2][i3+1]
			+weights[1]*rho[i1+1][i2][i3+1] ) 
			+weights[3]*( weights[0]*rho[i1][i2+1][i3+1] 
			+weights[1]*rho[i1+1][i2+1][i3+1] ) );

}

double Flamelet_Table::get_temperature(const double &Z_in, const double &Zvar_in, const double &Chi_in,bool refreshWeights) {
	
	if (refreshWeights) get_weights(Z_in,Zvar_in,Chi_in);
	
	return weights[4]*( weights[2]*( weights[0]*T[i1][i2][i3]
			+weights[1]*T[i1+1][i2][i3] )
			+weights[3]*( weights[0]*T[i1][i2+1][i3]
			+weights[1]*T[i1+1][i2+1][i3] ) ) 
			+weights[5]*( weights[2]*( weights[0]*T[i1][i2][i3+1]
			+weights[1]*T[i1+1][i2][i3+1] ) 
			+weights[3]*( weights[0]*T[i1][i2+1][i3+1] 
			+weights[1]*T[i1+1][i2+1][i3+1] ) );
}

double Flamelet_Table::get_laminarViscosity(const double &Z_in, const double &Zvar_in, const double &Chi_in,bool refreshWeights) {
	
	if (refreshWeights) get_weights(Z_in,Zvar_in,Chi_in);
	
	return  weights[4]*( weights[2]*( weights[0]*laminarViscosity[i1][i2][i3]
			+weights[1]*laminarViscosity[i1+1][i2][i3] )
			+weights[3]*( weights[0]*laminarViscosity[i1][i2+1][i3]
			+weights[1]*laminarViscosity[i1+1][i2+1][i3] ) ) 
			+weights[5]*( weights[2]*( weights[0]*laminarViscosity[i1][i2][i3+1]   
			+weights[1]*laminarViscosity[i1+1][i2][i3+1] ) 
			+weights[3]*( weights[0]*laminarViscosity[i1][i2+1][i3+1] 
			+weights[1]*laminarViscosity[i1+1][i2+1][i3+1] ) );
}

double Flamelet_Table::get_diffusivity(const double &Z_in, const double &Zvar_in, const double &Chi_in,bool refreshWeights) {
	
	if (refreshWeights) get_weights(Z_in,Zvar_in,Chi_in);
	
	return  weights[4]*( weights[2]*( weights[0]*diffusivity[i1][i2][i3]
			+weights[1]*diffusivity[i1+1][i2][i3] )
			+weights[3]*( weights[0]*diffusivity[i1][i2+1][i3]
			+weights[1]*diffusivity[i1+1][i2+1][i3] ) ) 
			+weights[5]*( weights[2]*( weights[0]*diffusivity[i1][i2][i3+1]   
			+weights[1]*diffusivity[i1+1][i2][i3+1] ) 
			+weights[3]*( weights[0]*diffusivity[i1][i2+1][i3+1] 
			+weights[1]*diffusivity[i1+1][i2+1][i3+1] ) );
}
