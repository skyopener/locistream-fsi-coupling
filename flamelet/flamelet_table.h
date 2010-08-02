
#ifndef FLAMELET_TABLE_H
#define FLAMELET_TABLE_H

// Standard library includes.
#include <vector>
#include <string>
using std::string ;
using std::vector ;

class Flamelet_Table {
	public:
	int nZ,nS,nChi;
	int i1,i2,i3;  // indices for Z,Zvar and Chi
	vector<double> weights;
	vector<double> Z;
	vector<double> S;
	vector<double> Chi; // This is actually Chi stoichiometric (Chi_st)
	vector<vector<vector<double> > > rho;
	vector<vector<vector<double> > > T;
	vector<vector<vector<double> > > laminarViscosity;
	vector<vector<vector<double> > > diffusivity;
	vector<vector<vector<double> > > Fchi; // The function to convert Chi to Chi_st
	//vector<vector<vector<vector<double> > > > Y; // mass fractions
	void read(string fileName);
	void get_weights(const double &Z_in,const double &Zvar_in,const double &Chi_in);
	double get_rho(const double &Z_in,const double &Zvar_in,const double &Chi_in,bool refreshWeights=true);
	double get_temperature(const double &Z_in,const double &Zvar_in,const double &Chi_in,bool refreshWeights=true);
	double get_laminarViscosity(const double &Z_in,const double &Zvar_in,const double &Chi_in,bool refreshWeights=true);
	double get_diffusivity(const double &Z_in,const double &Zvar_in,const double &Chi_in,bool refreshWeights=true);
};

#endif
