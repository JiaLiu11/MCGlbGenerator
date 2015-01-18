/*
Programmed by: Jia Liu

Contact information: liu.2053@osu.edu

Code Name: Event-by-Event Monte-Carlo Glauber(MCG) Generator

Version: 1.0

Purpose: Generate entropy density and calculate spatial eccentricity.

Input: (1)Atomic numbers for two nuclei, which should be the same;
	   (2)Impact parameter, which controls the centrality. There is a 
	      table in this folder tells the conversion between impact 
	      parameters and centrality for various colliding nuclei and 
	      energy. 
	   (3)Cross-section of Nucleon-Nucleon collision. In LHC, it is
	      around 60mb; This parameter is specified in Nucleus.cpp
       (4)Size, spacing of the final entropy density table;
       (5)Event number for MCG generator.

Result check: 
(1) For impact parameter equals 6.0fm, (arxiv: 1012.1657) gives the number
	of participants is in the range 220~270. This code reproduce it;    
(2) Checked with superMC, get the similar profile. The difference is
    on the overall scale, which should be adjusted with final multiplicity;
(3) Random test: Ran the code for several times, got different output.  

Need to do:
After generating the entropy density, the profile should be scaled 
according to final multiplicity before putting it to hydrodynamics simulation.

Revise history:
Apr.29, 2013 add a loop in the main program, which enables
			 generating multiple profiles.
Apr.26, 2013 hit() routine and count binary collision
		     and wounded nucleon;
Apr.16, 2013 Switch to invert CDF sampling method to 
             have higher efficiency;
Apr.09, 2013 Finish setWSCoordinates();
Apr.04, 2013 Finish generating position according to 
             Woods-Saxon distribution;
Mar.23, 2013 Creat header file and source file;
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "mc_glauber.h"
#include "time.h"
using namespace std;

int main()
{
	//parameters for generating nuclei configurations
	int atom_num = 208;  //atomic number of colliding nucleus
	double impact_parameter = 6.;  //specify impact parameter

    //parameters for entropy density table
	double sd_tbl_min = -13.;  
	double sd_tbl_max = 13.;  // unit: fm. It should be large enough
							  // to counts all the collisions. 
	double sd_tbl_step = 0.1;  // choose according to the precision 
							   // requirment and computer speed

	//parameters for the main program
	int nevents = 10;   //specify the total events of colllision
	int ecc_order = 2;  //specify the order of eccentricity


	//open file for dumping eccentricity
	ostringstream ecc_filename_stream;
	ecc_filename_stream.str("");  //clean before using it
	ecc_filename_stream << "data/Ecc_A_" << atom_num
	                << "_order_" << ecc_order << ".dat";
	ofstream ecc_of;
	ecc_of.open(ecc_filename_stream.str().c_str(),std::ios_base::app);

	//composee file names for entropy density profiles
	ostringstream sd_filename_stream;

	for(int i=0;i<nevents;i++)
	{
		mc_glauber* glauber_sim;   //create a MCG generator
		glauber_sim = new mc_glauber(atom_num, impact_parameter, 
			sd_tbl_min, sd_tbl_max, sd_tbl_step);

		glauber_sim->overlap();  //get binary collision

		//prepare file name of the entropy density profile
		sd_filename_stream.str("");
		sd_filename_stream << "data/Sd_A_"<<atom_num
		 				   << "_event_" << i+1 << ".dat";
		//dump entropy density table 				   
		glauber_sim->dumpSdTable(sd_filename_stream.str().c_str());  

		//dump eccentricity
		ecc_of << setw(8) << setprecision(5) << ecc_order
		       << setw(15)<< setprecision(8) << glauber_sim->getEccentricity(ecc_order)
		       << endl;  

		//clean up before next loop
		delete glauber_sim;
		cout << "Loop " << i+1 << " completed!" << endl << endl << endl;
	}

	ecc_of.close();  //finish eccentricity output file

	return 0;
}