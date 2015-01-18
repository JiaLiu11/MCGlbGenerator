#ifndef mc_glauber_h
#define mc_glauber_h

#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <cmath>
#include "Coordinates.h"
#include "Nucleon.h"
#include "Nucleus.h"

using namespace std;

class mc_glauber
{
protected:
	int atom_num;    //assign atomic number
	double impact_parameter;   //impact parameter for collision
	double glauber_entropy_width;  //the width of entropy deposited in the fireball
	double alpha;     //weight for wounded nucleon
	Nucleus* Nuc1;    //declare two nuclei
	Nucleus* Nuc2;
	vector<Coordinates*> wn_coordinates; //coordinates of wounded nucleons
	vector<Coordinates*> bc_coordinates; //coordinates of binary collision positions
	double** entropy_density; //table for entropy density: dS/(tau_0d^2rd\eta_s)|\eta_s=0

	double sd_tbl_lower, sd_tbl_upper, sd_tbl_step;  //parameters for entropy density table
	int max_sd_tbl;

	bool hit(double rp, double x0, double y0, double x1, double y1);   //if the collision happens
	void distEntropy();     //calculate entropy density in the in the transverse plane
							//sd = (1-alpha)*wn + alpha*bc
	void findSdCM(double* xcm, double *ycm);   //find the coordinate of center of entropy density

public:
	mc_glauber(int Atom_num, double Impact_parameter, 
			double Sd_tbl_min, double Sd_tbl_max, double Sd_tbl_step) ;
	~mc_glauber() ;
	void overlap();  //count wounded nucleons and binary collisions
	void dumpSdTable(string filename);  //dump entropy density table
	double getEccentricity(int order);   //calculate encentricity at specific order
};

#endif
