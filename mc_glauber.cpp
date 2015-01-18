/*
Programmed by: Jia Liu

Contact information: liu.2053@osu.edu

Owned by Code: Event-by-Event Monte-Carlo Glauber(MCG) Generator

Version: 1.0

Purpose: Manipulate nucleus and Simulate collision
1. Use Nucleus class to create nucleus;
2. overlap() can simulate colllisions and count the number of wounded
   nucleons and binary collisions; 
3. hit() function controls collision;
4. distEntropy() collects entropy generated by collisions. The radius, 
   glauber_entropy_width, is specify by user in this code. While superMC 
   chooses this parameters in a way to reproduce the nucleon-nucleon collision 
   cross-section;
5. dumpSdTable() dumps entropy profile;
6. findSdCM() finds the center of the profile;
6. getEccentricity() firstly calls findSdCM() to find the center of the profile,
   recenter it, then calculates eccentricity to any given order.
*/


#include <iostream>
#include <fstream>
#include <iomanip>
#include "mc_glauber.h"

using namespace std; 

mc_glauber::mc_glauber(int Atom_num, double Impact_parameter, 
			double Sd_tbl_min, double Sd_tbl_max, double Sd_tbl_step)
{	
	atom_num = Atom_num;    //read in atomic number
	impact_parameter = Impact_parameter;    //assign impact parameters
	//parameter for entropy density profile
	alpha = 0.3;   //weight of wounded nucleon

	//parameters for entropy density table
	sd_tbl_lower = Sd_tbl_min;
	sd_tbl_upper = Sd_tbl_max;
	sd_tbl_step = Sd_tbl_step;
	max_sd_tbl = (int)((sd_tbl_upper-sd_tbl_lower)/sd_tbl_step+0.1)+1;
	entropy_density = 0; //not assigned value
	glauber_entropy_width = 0.7;  //width for collecting entropy
								  //

	cout << "***********************************************" << endl
	     << "Monte-Carlo Glauber Model" << endl;

	//construct new nuclei
	Nuc1 = new Nucleus(atom_num);
	Nuc2 = new Nucleus(atom_num);
}

mc_glauber::~mc_glauber()
{
	for(int i=0;i<(int)wn_coordinates.size();i++)
		delete wn_coordinates[i];
    wn_coordinates.clear();

 	for(int i=0;i<(int)bc_coordinates.size();i++)
		delete bc_coordinates[i];
    bc_coordinates.clear();  

    if(entropy_density)
    {
    	for(int i=0;i<max_sd_tbl;i++)
     		delete [] entropy_density[i];
     	delete [] entropy_density;    	
    }

	delete Nuc1;
	delete Nuc2;
	cout << "***********************************************" << endl;
}

bool mc_glauber::hit(double rp, double x0, double y0, double x1, double y1)
{
	double distance = sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)* (y0-y1));
	if(distance <= 2.*rp)   //rule for hit: distance between two nucleons
							//is not greater than twice of the nucleon radius
		return true;
	return false;
} 


void mc_glauber::overlap()
{
	double nuc_size_1 = Nuc1->getNucleonSize();
	double nuc_size_2 = Nuc2->getNucleonSize();

	if(fabs(nuc_size_1 - nuc_size_2) > 1e-18)
	{
		cout << "Two different kinds of nucleus, cannot proceed" << endl;
		exit(0);
	}

	// cout << "Nucleon size is: " << nuc_size_1 << endl;

	long int binary_collision_num=0;
	//generate nucleus configuration
	Nuc1->generateConfiguration();
	Nuc2->generateConfiguration();

	//shift centers of nuclei in the x-direction
	Nuc1->shiftNucleus(impact_parameter/2.);
	Nuc2->shiftNucleus(-impact_parameter/2.);

	for(int i=0;i<atom_num;i++)
		for(int j=0;j<atom_num;j++)
		{
			//get the position of two 
			double x0, y0, z0;
			double x1, y1, z1;
			Nuc1->getNucleonCoordinates(i, &x0, &y0, &z0);
			Nuc2->getNucleonCoordinates(j, &x1, &y1, &z1);

			bool hit_here;
			hit_here = hit(nuc_size_1, x0, y0, x1, y1);
			if(hit_here == true)
			{
				binary_collision_num++;
				Nuc1->setNucleonBinaryCollision(i);
				Nuc2->setNucleonBinaryCollision(j);

				//position of binary collison
				double bc_x = (x0 + x1)/2.;
				double bc_y = (y0 + y1)/2.;

				//store binary collision positions debug
				Coordinates* ptr;
				ptr = new Coordinates(bc_x, bc_y);
				bc_coordinates.push_back(ptr);				
			}
		}//<-> for i=0:atom_num-1		

	//loop over to find all wounded nucleons
	long int counts1=0; 
	long int counts2=0;
	for(int i=0;i<atom_num;i++)
    {
    	if(Nuc1->getNucleonBCNum(i)>0)
    	{
    		counts1++;   //# of wounded nucleon in nucleus1 +1

    		//store the coordinates of wouned nucleons
    		Coordinates* ptr;
    		double wn_x, wn_y, z;
    		Nuc1->getNucleonCoordinates(i, &wn_x, &wn_y, &z);
			ptr = new Coordinates(wn_x, wn_y);
			wn_coordinates.push_back(ptr);
    	}

    	if(Nuc2->getNucleonBCNum(i)>0)
    	{
    		counts2++;   //# of wounded nucleon in nucleus2 +1

    		//store the coordinates of wouned nucleons
    		Coordinates* ptr;
    		double wn_x, wn_y, z;
    		Nuc2->getNucleonCoordinates(i, &wn_x, &wn_y, &z);
			ptr = new Coordinates(wn_x, wn_y);
			wn_coordinates.push_back(ptr);
    	}
	}
	//no collision at all
	if(binary_collision_num ==0)
	{
		cout << "No binary collsion!" << endl;
		exit(0);
	}		    

	// cout << "Collison process complete!" << endl
	//      << "Number of participants in nucleus 1: "<< counts1 << endl
	//      << "Number of participants in nucleus 2: "<< counts2 << endl;
	cout << "Number of participants: " << counts1+counts2 << endl
		 << "Total binary collision: " << binary_collision_num<<endl;
    distEntropy();
}

void mc_glauber::distEntropy()
{
/*
distribute entropy density across the transverse plane,
assume the weight of wounded nucleons and binary collisions are 
the same. Normalization needs to be fit by multiplicity, so it is 
not done in this function
*/
//	cout << "start to distribute entropy" << endl;
	//initialize entropy density table
	entropy_density = new double* [max_sd_tbl];
	for(int i=0;i<max_sd_tbl;i++)
	{
		entropy_density[i] = new double [max_sd_tbl];
		for(int j=0;j<max_sd_tbl;j++)
			entropy_density[i][j] = 0.;		
	}

	for(int i=0;i<max_sd_tbl;i++)
		for(int j=0;j<max_sd_tbl;j++)
		{
			double x_tbl = sd_tbl_lower + i*sd_tbl_step;
			double y_tbl = sd_tbl_lower + j*sd_tbl_step;

			//find contribution from wounded nucleons
			for(int k=0;k<(int)wn_coordinates.size();k++)
			{
				double wn_x = wn_coordinates[k]->getX();
				double wn_y = wn_coordinates[k]->getY();

				double distance = sqrt((x_tbl - wn_x)*(x_tbl - wn_x)
						   +(y_tbl - wn_y)*(y_tbl - wn_y));
				if(distance <= glauber_entropy_width)
					entropy_density[i][j]+=alpha;
			}

			//find contribution from binary collisions
			for(int k=0;k<(int)bc_coordinates.size();k++)
			{
				double bc_x = bc_coordinates[k]->getX();
				double bc_y = bc_coordinates[k]->getY();

				double distance = sqrt((x_tbl - bc_x)*(x_tbl - bc_x)
						   +(y_tbl - bc_y)*(y_tbl - bc_y));
				if(distance <= glauber_entropy_width)
					entropy_density[i][j]+=(1.-alpha);
			}
		}
	cout << "Entropy profile is generated!" << endl
	     << "Tips: fit to final multiplicity before put it into hydro!"
	     << endl << endl;
}


void mc_glauber::dumpSdTable(string filename)
{
	//safety check
	if(entropy_density == 0)
    {
    	cout << "No entropy density table" << endl;
    	exit(0);
    }

    ofstream of;
    of.open(filename.c_str(), std::ios_base::out);
    of << "% x, y from: " << sd_tbl_lower << " to " << sd_tbl_upper
       << ", with step: " << sd_tbl_step <<endl;

    of << "% # of wounded nucleons: "<< (int)wn_coordinates.size()
       << "; # of binary collisions: "<< (int)bc_coordinates.size()
       << endl;

    for(int i=0;i<max_sd_tbl;i++)
    {
    	for(int j=0;j<max_sd_tbl;j++)
    	{
    		of << setw(16) << setprecision(8) << entropy_density[i][j];
    	}
    	of << endl;
    }
    cout << "entropy density table dumped to file: "
         << filename << endl;
    cout << "Run Matlab script sd_plot.m to see the entropy density profile" << endl;

}

void mc_glauber::findSdCM(double* xcm, double* ycm)
{
/*find the weighted center of the profile, now use entropy density
as weighting function: xcm=(\int dxdy sd*x)/(\int dxdy sd) 
*/
    double x_ave=0., y_ave=0.;
    double weight=0.;    //use entropy density as weight
    double sd_total=0.;
    for(int i=0;i<max_sd_tbl;i++)
		for(int j=0;j<max_sd_tbl;j++)
		{
			double x= sd_tbl_lower + i*sd_tbl_step;
			double y= sd_tbl_lower + j*sd_tbl_step;
			weight = entropy_density[i][j];
			sd_total+=weight*sd_tbl_step*sd_tbl_step;  //total entropy

			x_ave+= weight*x*sd_tbl_step*sd_tbl_step;
			y_ave+= weight*y*sd_tbl_step*sd_tbl_step;
		}
	*xcm = x_ave/(sd_total + 1e-18);
	*ycm = y_ave/(sd_total + 1e-18);
}


double mc_glauber::getEccentricity(int order)
{
/*
Get the eccentricity of the profile at various order
*/
	double x_cm, y_cm;
	double ecc = 0.;
	double ecc_nu_real = 0.;
	double ecc_nu_img = 0.;
	double ecc_dn = 0.;

	findSdCM(&x_cm, &y_cm);
	//debug
	cout << "Current profile centered at: "
	     << "x=" << x_cm << ", "
	     << "y=" << y_cm << endl;
  
    for(int i=0; i<max_sd_tbl; i++)
    {
        for(int j=0; j<max_sd_tbl; j++)
        {
            double x = sd_tbl_lower + sd_tbl_step*i - x_cm;  //recenter the profile
            double y = sd_tbl_lower + sd_tbl_step*j - y_cm;
            double phi = atan2(y,x);

            double sd = entropy_density[i][j];

            ecc_nu_real += pow(x*x + y*y, double(order)/2.) 
                  * cos(order * phi) * sd * sd_tbl_step * sd_tbl_step;
            ecc_nu_img += pow(x*x + y*y, double(order)/2.) 
                  * sin(order * phi) * sd * sd_tbl_step * sd_tbl_step;
            ecc_dn += sd * pow(y*y + x*x, double(order)/2.)
                  * sd_tbl_step * sd_tbl_step;
        }       
    }//<-> for i=0:max_sd_tbl
  
    ecc = sqrt( ecc_nu_real * ecc_nu_real + ecc_nu_img * ecc_nu_img )
       /(ecc_dn + 1e-18);

    cout << "Spatial Eccentricity at " << order << "th order is: "
         << ecc << endl;

	return ecc;
}