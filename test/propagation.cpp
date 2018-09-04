#include <chrono>
#include <stdexcept>
#include <iostream>
#include <string>

#include "../src/SWPL.hpp"

void offsetObsScr(Wavefield& srcwf, double propDist, double xoffset, double yoffset){
	Wavefield::FieldAxis xaxis_dst = srcwf.getXaxis();
	Wavefield::FieldAxis yaxis_dst = srcwf.getYaxis();
	xaxis_dst.shift(xoffset);
	yaxis_dst.shift(yoffset);

	srcwf.RSprop(propDist, xaxis_dst, yaxis_dst);
}

void differntGridSpace(Wavefield& srcwf, double propDist, double gridRatio){
	Wavefield::FieldAxis xaxis_dst = srcwf.getXaxis();
	Wavefield::FieldAxis yaxis_dst = srcwf.getYaxis();
/*
	for(std::vector<double>::iterator itr=xaxis_dst.begin();itr<xaxis_dst.end();++itr){
		std::cout<<*itr<<", "<<std::flush;
	}std::cout<<std::endl;
*/
	xaxis_dst.scale(gridRatio);
	yaxis_dst.scale(gridRatio);
	srcwf.savevectord_bin(xaxis_dst, "xaxis_dst.bin");
/*
	for(std::vector<double>::iterator itr=xaxis_dst.begin();itr<xaxis_dst.end();++itr){
		std::cout<<*itr<<", "<<std::flush;
	}std::cout<<std::endl;
*/
	srcwf.RSprop(propDist, xaxis_dst, yaxis_dst);
}

void differntSizeObsScr(Wavefield srcwf, double propDist, int dstsizeEx_x, int dstsizeEx_y){
	Wavefield::FieldAxis xaxis_dst = srcwf.getXaxis();
	Wavefield::FieldAxis yaxis_dst = srcwf.getYaxis();
	xaxis_dst.extend_front(dstsizeEx_x/2);
	xaxis_dst.extend_back(dstsizeEx_x/2);
	yaxis_dst.extend_front(dstsizeEx_y/2);
	yaxis_dst.extend_back(dstsizeEx_y/2);

	srcwf.RSprop(propDist, xaxis_dst, yaxis_dst);
}

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
		//wf.RSprop(propz);
		//offsetObsScr(wf, propz, pitch*(xsize/2), 0);//);
		differntGridSpace(wf, propz, 0.5);
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

int main(int argc, char* argv[]){
	RStest(argc, argv);
	return 0;
}