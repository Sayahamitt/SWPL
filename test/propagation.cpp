#include <chrono>
#include <stdexcept>
#include <iostream>
#include <string>

//#define SWPL_UNUSE_OPENMP
#include "../src/SWPL.hpp"

void offsetObsScr(SWPL::Wavefield& srcwf, double propDist, double xoffset, double yoffset){
    SWPL::FieldAxis xaxis_dst = srcwf.getXaxis();
    SWPL::FieldAxis yaxis_dst = srcwf.getYaxis();
    xaxis_dst.shift(xoffset);
    yaxis_dst.shift(yoffset);

    srcwf.RSprop(propDist, xaxis_dst, yaxis_dst);
}

void differntGridSpace(SWPL::Wavefield& srcwf, double propDist, double gridRatio){
    SWPL::FieldAxis xaxis_dst = srcwf.getXaxis();
    SWPL::FieldAxis yaxis_dst = srcwf.getYaxis();

    xaxis_dst.scale(gridRatio);
    yaxis_dst.scale(gridRatio);

    srcwf.RSprop(propDist, xaxis_dst, yaxis_dst);
}

void differntSizeObsScr(SWPL::Wavefield& srcwf, double propDist, int dstsizeEx_x, int dstsizeEx_y){
    SWPL::FieldAxis xaxis_dst = srcwf.getXaxis();
    SWPL::FieldAxis yaxis_dst = srcwf.getYaxis();

    xaxis_dst.extend_front(dstsizeEx_x/2);
    xaxis_dst.extend_back(dstsizeEx_x/2);
    yaxis_dst.extend_front(dstsizeEx_y/2);
    yaxis_dst.extend_back(dstsizeEx_y/2);
    
    srcwf.RSprop(propDist, xaxis_dst, yaxis_dst);
}

void differntAxis(SWPL::Wavefield& srcwf, double propDist){
    double xf = srcwf.getXaxis().front()/2;
    double xb = srcwf.getXaxis().back()/2;
    double yf = srcwf.getYaxis().front()/1;
    double yb = srcwf.getYaxis().back()/1;

    SWPL::FieldAxis xaxis_dst(xf,xb,srcwf.getXaxis().Pitch());
    //Wavefield::FieldAxis xaxis_dst = srcwf.getXaxis();
    SWPL::FieldAxis yaxis_dst(yf,yb,srcwf.getYaxis().Pitch());
    //Wavefield::FieldAxis yaxis_dst = srcwf.getYaxis();
    try{
        srcwf.RSprop(propDist, xaxis_dst, yaxis_dst);
    }catch(const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
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
    }catch(...){
        std::cout<<"error in argument parsing"<<std::endl;
    }
    
    SWPL::Wavefield wf(xsize,ysize,pitch,wvl);
    //wf.setCircAparture(1);
    wf.setCircAparture(aptD);
    std::chrono::system_clock::time_point start = std::chrono::system_clock::now(); 
    try{
        //wf.RSprop(propz);
        //offsetObsScr(wf, propz, pitch*(xsize/2), 0);//);
        //differntGridSpace(wf, propz, 4);
        //differntSizeObsScr(wf, propz, 128, 128);
        differntAxis(wf,propz);
    }catch(std::invalid_argument errorcode){
        std::cout<<"error. unexcepted argument."<<std::endl;
        return 0;
    }catch(...){
        std::cout<<"something error in calling"<<std::endl;
    }
    std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
    std::chrono::system_clock::duration calcdur = end-start;
    SWPL::binWriteVCD(wf.getField(), ofpath);
    std::cout<< "using eigen :" <<std::chrono::duration_cast<std::chrono::milliseconds>(calcdur).count() <<std::endl;
    return 0;
}

int main(int argc, char* argv[]){
    RStest(argc, argv);
    return 0;
}