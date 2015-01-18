/*
Programmed by: Jia Liu

Contact information: liu.2053@osu.edu

Revised time: Apr.26, 2013

Owned by Code: Event-by-Event Monte-Carlo Glauber(MCG) Generator

Purpose: Store the information of one nucleon 
1. Called by Nucleus class
2. Save and give the position of one nucleon
3. Remember the number of collision experienced by this nucleon
*/

#ifndef Nucleon_h
#define Nucleon_h

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;

class Nucleon
{
protected:
	double x, y, z;    //position of a nucleon
	double radius;  // for disk-like nucleon, radius of one nucleon
	int binary_collision_num;  //total number of binary collisions

public:
	Nucleon(double r0, double x0, double y0, double z0) {	
		x=x0; y=y0; z=z0;
	    radius = r0; 
	    binary_collision_num = 0;
	};
	~Nucleon() {};

	void setX(double x0)  { x=x0; }
	void setY(double y0)  { y=y0; }
	void setZ(double z0)  { z=z0; }
	void setRadius(double r0) { radius = r0; }

	double getX()  { return x; }
	double getY()	{ return y; }
	double getZ()  { return z; }
	double getRadius()  {return radius;}
	int getBCNum() {return binary_collision_num;}

	double distanceTo(Nucleon* nucleon1) {
		double x1=nucleon1->getX();
		double y1=nucleon1->getY();
		double z1=nucleon1->getZ();
		return sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1));
	}
	void setBinaryCollision() {binary_collision_num++;}  //# of binary collision ++
};

#endif
