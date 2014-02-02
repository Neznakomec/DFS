#include <mpi.h>

#include <stdio.h>
#include "yavnMPI3D.h"

int main()
{
	int count = 9;
	char* callString[] = {"very_cool_software", "10", "10", "10",
					 "99", "99", "99", 
					 "0.1", "100"};
	run(count, callString);

}

