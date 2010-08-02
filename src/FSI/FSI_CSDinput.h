
#ifndef CSDINPUT_H
#define CSDINPUT_H

// Standard library includes.
#include <vector>
#include <string>
using std::string ;
using std::vector ;

class CSDinput {
	public:
	int numberOfNodes ;
	int numberOfElements ;
	int numberOfBCs ;
	vector<vector<double, 3> > csdMesh;
	vector<vector<int, 3> > csdConnectivity;
	vector<int> csdBCdof ;
	vector<real> csdBCZeroConstraint ;
	
	void readMesh(string fileName);
	void readConnectivity(string fileName);
	void readBC(string fileName);
};

#endif
