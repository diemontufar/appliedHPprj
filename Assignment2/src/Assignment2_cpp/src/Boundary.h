/*
 * Boundary.h
 *
 */

#include <fstream>
#include <iostream>
#include <math.h>
#include <string>

using namespace std;

class	Boundary
{
public:
	Boundary()
	{

	}
    string  name_;
    string  type_;
	int		N_;
    int*    indices_;
	double	value_;
};


