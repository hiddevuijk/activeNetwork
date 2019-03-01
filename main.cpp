
#include "ConfigFile.h"


#include "xyz.h"
#include "system.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <vector>
#include <string>


using namespace std;



int main()
{

	// read input into config
	ConfigFile config("input.txt");

	// read integration parameters
	Integration int_params(config);

	// read system parameters
	System system(config);

	for(unsigned int ti=0; ti<int_params.Nt; ++ti) {
		system.step();
	}

	for(unsigned int i =0;i<system.N;++i)	
		cout << system.r[i] << endl;

	return 0;
}



