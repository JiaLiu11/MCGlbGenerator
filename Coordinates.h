/*
Coded by      Jia Liu 
              liu.2053@osu.edu
Revised time: Apr.26, 2013

Owned by Code: Event-by-Event Monte-Carlo Glauber(MCG) Generator

Purpose: store coordinates of binary collision and wouneded 
         which are located in the transverse plane
*/


#ifndef Coordinates_h
#define Coordinates_h

class Coordinates
{
protected:
	double x, y;

public:
	Coordinates(double x0, double y0) {	
		x=x0; y=y0; 
	};
	~Coordinates() {};

	void setX(double x0)  { x=x0; }
	void setY(double y0)  { y=y0; }

	double getX()  { return x; }
	double getY()	{ return y; }

};

#endif

