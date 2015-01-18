/*
Programmed by: Jia Liu

Contact information: liu.2053@osu.edu

Owned by Code: Event-by-Event Monte-Carlo Glauber(MCG) Generator

Purpose: Generate and Store information of one nucleus 
1. The positions of nucleons satisfy Woods-Saxon distribution
2. The position of each nucleon is sampled by invert CDF method;
3. shiftNucleus() is called by the mc_glauber class to shift the nucleus
   to a position specified by the impact parameter. It should be done 
   before collision
4. dumpNucleonsCoordinates() is designed for testing the positions of 
   sampled nucleons. Use scatter3() in matlab to plot nucleons; 
   hist() in Matlab can plot histogram, use it to see if the distribution
    of generated nucleons satisfies Woods-Saxon distribution.
*/

#include <cmath>
#include <iomanip>
#include "Nucleus.h"
#include "arsenal.h"

using namespace std;

extern unsigned long int random_seed ();   // routine to generate a seed
 											
Nucleus::Nucleus(int A_num, double NS, double MS)
{
	A  = A_num;
	nS = NS;
	mS = MS;

  //parameters for CDF table
  tbl_min = 0.;
  tbl_max = 20.;
  tbl_step = 0.01;
  max_table = (long int)((tbl_max-tbl_min)/tbl_step+0.1)+1; //find the length of the CDF lookup table	
  cdf_max = 0.;  //Initial value, will be assigned a value later
  cdf_table = 0; //CDF table has not been initialized

  double sigma_nn = 60.;  //nucleon-nucleon cross section unit: mb
  nucleon_radius = sqrt(0.1/(2.*M_PI) * sigma_nn)/2.; //effective radius = sqrt(sigma_nn/2/pi)/2
                                                //0.1 for convert from sqrt(barn) to fm
                                                //a factor of 2 since sigma_nn is effective x-section
  wsInitializion();   //get parameters from Wood-Saxon Model
}

Nucleus::~Nucleus()
{
	for(int i=0;i<(int)nucleons.size();i++) //clear up nucleons
	{
    	delete nucleons[i];
  	}
  nucleons.clear();
  
  if(cdf_table)  //clean up CDF lookup table
    delete [] cdf_table;

}


void Nucleus::wsInitializion()  
//calculate ws_r, ws_d from a given A
{
	ws_r = 1.25 * pow(double(A), 1./3.);  //unit: fm
	ws_a = 0.5;		//unit: fm
	ws_prob0 = 0.16;  // rho0

	cout<<"Atom number is: "<<A<<endl;  //debug
}


void Nucleus::generateConfiguration()
{
  prepareCDFtable();  //generate CDF table and shift it
  //begin invert CDF sampling
  getWSCoordinates(A);
}


void Nucleus::getWSCoordinates(int atom_num)
{
  //invert CDF to get the coordinates

  if(cdf_table == 0 || cdf_max ==0.)  //no cdf table to invert
                                      //or cdf table is wrongly found
  {
    cout<< "No CDF table, or CDF is wrong! Exit...." << endl;
    exit(-0);
  }

  // cout << "Start to get nucleon coordinates:" << endl;

  //prepare CDF table for invert, since it cannot be used by external functions
  vector<double>* cdf_table_copy=new vector<double>(max_table,0.0);    //copy CDF table
  for(long int i=0;i<max_table;i++)
    (*cdf_table_copy)[i]=cdf_table[i];

  srand48(random_seed ());  //random seed
  for(int count = 0; count < atom_num; count ++)
  {
    double t_rand=drand(0., cdf_max);  //get a random number between 0 ~ max value of CDF
    double cdf_prob = t_rand;  //position probability for CDF table
    double r_sampled = 0.;  //sampled spherical coordinate r

    int r_idx = binarySearch(cdf_table_copy, cdf_prob, false);  //find the index of CDF element
    r_sampled = tbl_min + tbl_step* r_idx;   //covert to coordinate r

    //generate theta and phi
    double cos_theta = drand(-1., 1.);
    double sin_theta = sqrt(1 - cos_theta * cos_theta);
    double phi = drand(0., 2.*M_PI);
    //transform to Cartisan coordinates
    double x = r_sampled * sin_theta * cos(phi);
    double y = r_sampled * sin_theta * sin(phi);
    double z = 0.;   //z=0 due to lorentz contraction
  
    //construct a new nucleon and put the pointer to the vector "nucleons"
    Nucleon* ptr;
    ptr = new Nucleon(nucleon_radius, x, y, z); 
    nucleons.push_back(ptr);
  }
  cout << "Nucleus Configuration has been generated!" << endl << endl;

  //clear up before leave
  delete cdf_table_copy;
}


void Nucleus::prepareCDFtable(void)
{
  // cout << "Start to generate CDF table " << endl;

  //initalize CDF look-up table
  cdf_table = new double[max_table];

  //generate CDF look-up table
  cdf_table[0] = 0. ;
  for (long int i=1;i<max_table; i++)
  {
    double r_step = tbl_min +  i * tbl_step;  //current position
    cdf_table[i] = cdf_table[i-1] + getWoodsSaxonModel(r_step);
  }
  //Assign value for the "normalization"
  cdf_max = cdf_table[max_table-1];
  // cout << "CDF table is ready!" << endl << endl;
}



double Nucleus::getWoodsSaxonModel(double distance)
{
// Woods-Saxon model
	double rho0 = ws_prob0;
  double weight = distance * distance;   //geometry factor for 3D position sampling
	double result =  weight * rho0 /(1 + exp((distance - ws_r)/ws_a));

  return result;
}



void Nucleus::shiftNucleus(double x_ctr, double y_ctr)
{
  // cout << "start to shift nucleus to a new center: "
  //      << "(" << x_ctr <<", " << y_ctr << ")" <<endl;
  for(int i=0;i<(int)nucleons.size();i++)
  {
    double x = nucleons[i]->getX() - x_ctr;
    double y = nucleons[i]->getY() - y_ctr;

    //set the new center
    nucleons[i]->setX(x);
    nucleons[i]->setY(y);
  }
}



void Nucleus::dumpNucleonsCoordinates(string filename)
{
  ofstream of;
  of.open(filename.c_str(), std::ios_base::out);
  for(int i=0;i<(int)nucleons.size();i++)
  {
    double x = nucleons[i]->getX();
    double y = nucleons[i]->getY();
    double z = nucleons[i]->getZ();

    of << setw(16) << setprecision(10) << x
       << setw(16) << setprecision(10) << y
       << setw(16) << setprecision(10) << z 
       << endl;
  }
  of.close();
  // cout << "Nucleus has been shifted!" << endl << endl;
}
