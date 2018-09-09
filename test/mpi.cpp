#include <chrono>
#include <stdexcept>
#include <iostream>
#include <string>
#include <mpi.h>

#include "../src/SWPL.hpp"

int main(int argc, char* argv[]){
	int xsize = 0;
	int ysize = 0;
	double pitch = 0;
	double wvl = 0;
	double aptD = 0;
	double propz = 0;
	std::string ofpath;

	if(argc != 8){

		MPI_Init(&argc, &argv);

	   int rank;
	   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	   std::cout << "Hello world from " << rank << "-th process" << std::endl;

	   MPI_Barrier(MPI_COMM_WORLD);
	   if (rank == 0) {
	       std::cout << "Press enter to close" << std::endl;
	       std::cin.ignore(10000, '\n');
	   }

	   	MPI_Finalize();
		return 0;
	}
	try{
		xsize = std::stoi(std::string(argv[1]));
		ysize = std::stoi(std::string(argv[2]));
		pitch = std::stod(std::string(argv[3]));
		wvl = std::stod(std::string(argv[4]));
		aptD = std::stod(std::string(argv[5]));
		propz = std::stod(std::string(argv[6]));
		ofpath = std::string(argv[7]);
		std::cout<<ofpath<<std::endl;
	}catch(std::invalid_argument errorcode){
		std::cout<<"error. unexcepted argument."<<std::endl;
		return 0;
	}
	
	Wavefield wf(xsize,ysize,pitch,wvl);
	//wf.setCircAparture(1);
	wf.setCircAparture(aptD);
	//wf.showWFcons();
	//wf.saveWfield_bin("source.bin");
	std::chrono::system_clock::time_point start = std::chrono::system_clock::now(); 
	try{
		//wf.RSprop(propz);
		//offsetObsScr(wf, propz, pitch*(xsize/2), 0);//);
		//differntGridSpace(wf, propz, 4);
		//differntSizeObsScr(wf, propz, 128, 128);
	}catch(std::invalid_argument errorcode){
		std::cout<<"error. unexcepted argument."<<std::endl;
		return 0;
	}
	std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
	std::chrono::system_clock::duration calcdur = end-start;
	wf.saveWfield_bin(ofpath);
	std::cout<< "using eigen :" <<std::chrono::duration_cast<std::chrono::milliseconds>(calcdur).count() <<std::endl;

	return 0;
}