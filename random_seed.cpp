//*********************************************************************
//
//  random_seed returns a starting seed for random number generators
//
//  Based on an contribution to the gsl-discuss@sources.redhat.com
//   discussion list by Robert Brown of Duke (rgb@phy.duke.edu)
//   on 01-Mar-2004.
//
//  Using /dev/urandom to avoid blocking, but the random numbers
//   may not be secure.
//
//  Revision history:
//   05-Mar-2004  original version based on Robert Brown's code snippet
//
//**********************************************************************   

// include files
#include <iostream>		// cout and cin
#include <iomanip>		// manipulators like setprecision
#include <fstream>		// file input and output
using namespace std;		// we need this when .h is omitted

#include <sys/time.h>		// for the random number seed

// function prototypes
unsigned long int random_seed ();	// routine to generate a seed

//********************************************************************
unsigned long int
random_seed ()
{
  unsigned long int seed = 0;
  unsigned long int *seed_ptr = &seed;
  struct timeval tv;

  // open a stream to read from /dev/random as binary
  ifstream dev_urandom ("/dev/urandom", ios::in | ios::binary);

  if (dev_urandom.bad ())
    {				// failed to open /dev/urandom
      cout << "Using time from gettimeofday() for seed . . ." << endl;
      gettimeofday (&tv, 0);
      seed = tv.tv_sec + tv.tv_usec;  // build a pseudo-random seed
      seed_ptr = &seed;
    }
  else
    {
      // get bits from /dev/urandom to fill the seed
      dev_urandom.read ((char *) seed_ptr, sizeof (long int));
    }

  dev_urandom.close ();		//  close the stream
  return (*seed_ptr);		// send back the seed!
}

