#include <iostream>
#include <mpi.h>
#include <string>
#include "dg_param.h"

//variable you must change--------------------------------------------

/// mesh file
namespace fileinfo{
	const std::string meshfile = "example.txt";
	const std::string fileplace = "../gmsh_files/";
}

//---------------------------------------------------------------------


// variables you do not need to change----------------------------------

/// mpi variables
namespace mpi{
	int rank;    //!> process rank
	int num_proc;    //!> total number of processor
}

//-----------------------------------------------------------------------
