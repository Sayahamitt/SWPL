#include <chrono>
#include <stdexcept>
#include <iostream>
#include <string>

//#define SWPL_UNUSE_OPENMP
#include "../src/SWPL.hpp"

int RStest(int argc, char* argv[]){
	int xsize = 0;
	int ysize = 0;
	double pitch = 0;
	double wvl = 0;
	double aptD = 0;
	double propz = 0;
	std::string ofpath;

	if(argc != 8){
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
	wf.saveWfield_bin("source.bin");
	std::chrono::system_clock::time_point start = std::chrono::system_clock::now(); 
	try{
		wf.RSprop(propz);
	}catch(std::invalid_argument errorcode){
		std::cout<<"error. unexcepted argument."<<std::endl;
		return 0;
	}
	std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
	std::chrono::system_clock::duration calcdur = end-start;
	wf.saveWfield_bin(ofpath);
	std::cout<< "using eigen :" <<std::chrono::duration_cast<std::chrono::milliseconds>(calcdur).count() <<std::endl;

	wf.setCircAparture(aptD);
	start = std::chrono::system_clock::now(); 
	try{
		wf.RSprop_nosimd(propz);
	}catch(std::invalid_argument errorcode){
		std::cout<<"error. unexcepted argument."<<std::endl;
		return 0;
	}
	end = std::chrono::system_clock::now();
	calcdur = end-start;
	//std::cout<<double(1)/std::complex<double>(0,1)<<std::endl;
	std::cout<< "only openmp :" <<std::chrono::duration_cast<std::chrono::milliseconds>(calcdur).count() <<std::endl;

	return 0;
}

int main(int argc, char* argv[]){
	//et();
	RStest(argc, argv);
	return 0;
}
