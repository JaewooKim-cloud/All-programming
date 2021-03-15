/*This is a skeleton code file for use with the Finite Element Method for Problems in Physics.
It uses the deal.II FEM library, dealii.org*/

//Include files
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "FEM1.h"
//#include "writeSolutions.h"

using namespace dealii; // it is like using vecor::std something likt that

//The main program, using the FEM class
int main (){
  try{
    deallog.depth_console (0);

		//Specify the basis function order: 1, 2, or 3
		unsigned int order = 1;

		//Specify the subproblem: 1 or 2
		unsigned int problem = 2;

    FEM<1> problemObject(order,problem); // the name of class is FEM. FEM(1) means one dimensional problem
    
    //Define the number of elements as an input to "generate_mesh"
    problemObject.generate_mesh(10); //e.g. a 10 element mesh, all the element has equal length
    problemObject.setup_system(); // it initializes the data objects,  it resizes force vector, solution vector, global stiffness matrix. it also set up basis.
    problemObject.assemble_system();// create k local, f local and assemble those into global stiffness matrix and apply bcs
    problemObject.solve(); // solve d=k^-1f
    std::cout << problemObject.l2norm_of_error() << std::endl;
    
    //write output file in vtk format for visualization
    problemObject.output_results();
    
    //write solutions to h5 file
    char tag[21];
    sprintf(tag, "CA1_Order%d_Problem%d",order,problem);
   // writeSolutionsToFileCA1(problemObject.D, problemObject.l2norm_of_error(), tag);
  }
  catch (std::exception &exc){
    std::cerr << std::endl << std::endl
	      << "----------------------------------------------------"
	      << std::endl;
    std::cerr << "Exception on processing: " << std::endl
	      << exc.what() << std::endl
	      << "Aborting!" << std::endl
	      << "----------------------------------------------------"
	      << std::endl;

    return 1;
  }
  catch (...){
    std::cerr << std::endl << std::endl
	      << "----------------------------------------------------"
	      << std::endl;
    std::cerr << "Unknown exception!" << std::endl
	      << "Aborting!" << std::endl
	      << "----------------------------------------------------"
	      << std::endl;
    return 1;
  }

  return 0;
}
