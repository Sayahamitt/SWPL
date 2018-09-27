#include <mpi.h>
#include <chrono>
#include <stdexcept>
#include <iostream>
#include <string>

#include <Eigen/Core>
#include "../src/SWPL.hpp"

std::vector<std::complex<double>> mpitest(SWPL::Wavefield& src, double wvl, const SWPL::FieldAxis& dst_xaxis,const SWPL::FieldAxis& dst_zaxis){
    const double xaxisPitch = dst_xaxis.Pitch();
    const double xaxisFirst = dst_xaxis.front();
    const double xaxisEnd = dst_xaxis.back();
    const size_t xaxisNum = dst_xaxis.size();
    const double yaxisPitch = dst_zaxis.Pitch();
    const double yaxisFirst = dst_zaxis.front();
    const double yaxisEnd = dst_zaxis.back();
    const size_t yaxisNum = dst_zaxis.size();
    
    int mpi_procRank;
    int mpi_commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_commSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_procRank);

    //自スレッドの担当範囲を決める
    //二次元の行列で定義される計算範囲をY軸についてのみ分割する
    //分割数で行数を割った余りの分の計算量をさらに各スレッドで1行分ずつ分担する
    int quotient = dst_zaxis.size()/mpi_commSize;//C++でint型の割り算は切り捨て
    int modulo = dst_zaxis.size()%mpi_commSize;
    bool isModZero = (modulo==0 ? true : false);
    int startPos = quotient*mpi_procRank;
    int endPos = quotient*(mpi_procRank+1)-1;
    if(modulo-mpi_procRank > 0){
        startPos += mpi_procRank;
        endPos += mpi_procRank+1;
    }else if(isModZero==false){
        startPos += modulo;
        endPos = startPos+quotient-1;
    }
    Wavefield::FieldAxis thZaxis(dst_zaxis[startPos],dst_zaxis.Pitch(),endPos-startPos+1);
    std::cout<< "quotient: "<< quotient<< ", modulo: "<<modulo<< ", ThID:  " << mpi_procRank << ", StartPos: "<< startPos << ", endPos: "<< endPos << ", Num: "<< endPos-startPos+1<< std::endl;
    std::cout<< "quotient: "<< quotient<< ", modulo: "<<modulo<< ", ThID:  " << mpi_procRank << ", Start: "<< thZaxis.front() << ", End: "<< thZaxis.back() << ", Num: "<< thZaxis.size()<< std::endl;

    
    std::vector<std::complex<double>> bs(dst_xaxis.size()*thZaxis.size(),std::complex<double>(0));
    Wavefield::FieldAxis src_xaxis= src.getXaxis();
    Wavefield::FieldAxis src_yaxis = src.getYaxis();

    const Eigen::ArrayXXcd srcwf = Eigen::Map<Eigen::ArrayXXcd>(src.getField().data(),src_yaxis.size(),src_xaxis.size());
    const Eigen::ArrayXXd srcax_y = Eigen::Map<Eigen::ArrayXd>(src_yaxis.data(),src_yaxis.size()).replicate(1, src_xaxis.size());
    const Eigen::ArrayXXd srcax_x = Eigen::Map<Eigen::ArrayXd>(src_xaxis.data(),src_xaxis.size()).transpose().replicate(src_yaxis.size(), 1);

    /*
    std::complex<double> jK = Wavefield::imunt*(2*SWPL_PI/wvl);
    Eigen::ArrayXXd r01 = Eigen::ArrayXXd::Zero(src_yaxis.size(),src_xaxis.size());
    double dst_y = src_yaxis[src_yaxis.size()/2];
    #pragma omp parallel for private(r01)
    for(int ii=0;ii<thZaxis.size();ii++){
        for(int jj=0;jj<dst_xaxis.size();jj++){
            r01 = (std::pow(thZaxis[ii],2)+(srcax_x-dst_xaxis[jj]).square()+(srcax_y-dst_y).square()).sqrt();
            bs[dst_xaxis.size()*ii+jj] = (srcwf*((jK*r01).exp()/r01.square())).sum()*(thZaxis[ii]/(Wavefield::imunt*wvl))*(dst_xaxis.Pitch()*dst_zaxis.Pitch())*(src_xaxis.Pitch()/dst_xaxis.Pitch())*(src_yaxis.Pitch()/dst_zaxis.Pitch());

        }
    }*/

    
    #pragma omp parallel for
    for(int ii=0;ii<thZaxis.size();ii++){
        for(int jj=0;jj<dst_xaxis.size();jj++){
            bs[dst_xaxis.size()*ii+jj] = src.RSpoint(thZaxis[ii], dst_xaxis[jj], src_yaxis[src_yaxis.size()/2], dst_xaxis.Pitch(), dst_zaxis.Pitch(),
             srcwf, srcax_y, srcax_x);
        }
    }
    return bs;
}

int main(int argc, char* argv[]){
    int xsize = 0;
    double zpitch = 0;
    double xpitch = 0;
    double wvl = 0;
    double aptD = 0;
    double propz = 0;
    std::string ofpath;
    bool unifiedFileOut = false;

    if(argc != 8){
        return 0;
    }
    try{
        xsize = std::stoi(std::string(argv[1]));
        xpitch = std::stod(std::string(argv[2]));
        wvl = std::stod(std::string(argv[3]));
        aptD = std::stod(std::string(argv[4]));
        propz = std::stod(std::string(argv[5]));
        zpitch = std::stod(std::string(argv[6]));
        ofpath = std::string(argv[7]);
    }catch(std::invalid_argument errorcode){
        std::cout<<"error. unexcepted argument."<<std::endl;
        return 0;
    }

    //int mpi_argc=0;char*** mpi_argv;
    MPI_Init(&argc,&argv);
    int mpi_procRank;
    int mpi_commSize;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_commSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_procRank);

    SWPL::Wavefield wf(xsize,xsize,xpitch,wvl);
    wf.setCircAparture(aptD);
    std::chrono::system_clock::time_point start;
    if(mpi_procRank==0){
        start = std::chrono::system_clock::now(); 
    }
    
    std::vector<std::complex<double>> bs;
    try{
        //std::cout<<propz<<","<<zpitch<<","<<int(propz/zpitch)<<std::endl;
        Wavefield::FieldAxis zaxis = Wavefield::FieldAxis(zpitch,zpitch,int(propz/zpitch));
        bs = mpitest(wf,wvl,wf.getXaxis(),zaxis);
        MPI_Barrier(MPI_COMM_WORLD);
    }catch(std::invalid_argument errorcode){
        std::cout<<"error. unexcepted argument."<<std::endl;
        return 0;
    }

    if(mpi_procRank==0){
        std::chrono::system_clock::time_point end = std::chrono::system_clock::now();
        std::chrono::system_clock::duration calcdur = end-start;
        std::cout<< "using eigen :" <<std::chrono::duration_cast<std::chrono::milliseconds>(calcdur).count() <<std::endl;
    }

    //ファイル出力
    if(unifiedFileOut){
        std::string foutpath = ofpath+std::string("th_obs.bin");
        std::vector<std::complex<double>> wfvec = bs;
        std::ofstream fout;
        if(mpi_procRank==0){
            fout.open(foutpath.c_str(), std::ios::out|std::ios::binary|std::ios::trunc);
            fout.close();
        }
        for(int ii=0;ii<mpi_commSize;ii++){
            MPI_Barrier(MPI_COMM_WORLD);
            if (mpi_procRank != ii){
                continue;
            }
            
            fout.open(foutpath.c_str(), std::ios::out|std::ios::binary|std::ios::app);
            std::for_each(wfvec.begin(),wfvec.end(),[&fout](std::complex<double> num){
                    double r = std::real(num);
                    double i = std::imag(num);
                    fout.write((char *) &r,sizeof( double ) );
                    fout.write((char *) &i,sizeof( double ) );
                });
            fout.close();
        }
    }else{
        std::string foutpath = std::to_string(mpi_procRank)+std::string("th_obs.bin");
        SWPL::binWriteVCD(wf.getField(), foutpath);
    }
    MPI_Finalize();
    return 0;
}