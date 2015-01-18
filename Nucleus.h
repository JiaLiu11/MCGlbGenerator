#ifndef Nucleus_h
#define Nucleus_h

#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include "stdlib.h"
#include "Nucleon.h"
#include "arsenal.h"

using namespace std;


class Nucleus
{
protected:
	int A;    //atom number
	double ws_r, ws_a;  //Wood-Saxon model parameters
	double ws_prob0;
	double nS;  //nucleon size
	double mS;  //minimum separation
	vector<Nucleon*> nucleons;
	double nucleon_radius;

	//culmulative distribution function look-up table
	double *cdf_table;  
	double cdf_max;   //last elements of cdf table which used to do normalization
	double tbl_min, tbl_max;  //lower and upper limit for nucleon position r in CDF table
	double tbl_step;   //spacing for position r
	long int max_table;  //length of the CDF lookup table

	void wsInitializion(void);  //calculate ws_r, ws_d from a given atom number A
	void prepareCDFtable(void);  //generate CDF look up table
	void getWSCoordinates(int atom_num); //get nucleon coordinates
											//by invert CDF
	double getWoodsSaxonModel(double distance);

public:
	Nucleus(int A_num, double NS=0.4, double MS=0.4);
	~Nucleus();

	void generateConfiguration(void);  //generate nuleus configuration
	void shiftNucleus(double x_ctr, double y_ctr=0.);//shift the nucleus down in the x-y plane
													 //to centered in(x_ctr, y_ctr)
	double getNucleonSize(void) {return nucleon_radius;}	
	void getNucleonCoordinates(int idx, double* x, double* y, double* z) {
		*x= nucleons[idx]->getX();  *y=nucleons[idx]->getY();	*z=nucleons[idx]->getZ();
	}

	void setNucleonBinaryCollision(int idx){
		nucleons[idx]->setBinaryCollision(); 
	} //# of binary collision for a nucleon ++	

	int getNucleonBCNum(int idx) {
		return nucleons[idx]->getBCNum();
	} //get the # of collisions for a nucleon					
	void dumpNucleonsCoordinates(string filename);  //output the positions of nucleons
	 
};

#endif
