#include <chrono>
#include <stdexcept>
#include <iostream>
#include <string>
#include <mpi.h>

#include "../src/SWPL.hpp"

void mpitest(SWPL::Wavefield& src,const SWPL::FieldAxis& dst_xaxis,const SWPL::FieldAxis& dst_yaxis, double propz){
	const double xaxisPitch = dst_xaxis.Pitch();
	const double xaxisFirst = dst_xaxis.front();
	const double xaxisEnd = dst_xaxis.back();
	const size_t xaxisNum = dst_xaxis.size();
	const double yaxisPitch = dst_yaxis.Pitch();
	const double yaxisFirst = dst_yaxis.front();
	const double yaxisEnd = dst_yaxis.back();
	const size_t yaxisNum = dst_yaxis.size();
	
	int mpi_procRank;
	int mpi_commSize;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_commSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_procRank);

	//自スレッドの担当範囲を決める
	//二次元の行列で定義される計算範囲をY軸についてのみ分割する
	//分割数で行数を割った余りの分の計算量をさらに各スレッドで1行分ずつ分担する
	int quotient = dst_yaxis.size()/mpi_commSize;//C++でint型の割り算は切り捨て
	int modulo = dst_yaxis.size()%mpi_commSize;
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

	//Wavefield::FieldAxis thAxis(dst_yaxis[startPos],dst_yaxis[endPos],dst_yaxis.Pitch());
	SWPL::FieldAxis thAxis(dst_yaxis[startPos],dst_yaxis.Pitch(),endPos-startPos+1);
	//std::cout<< "quotient: "<< quotient<< ", modulo: "<<modulo<< ", ThID:  " << mpi_procRank << ", srcStart: "<< dst_yaxis[startPos] << ", srcEnd: "<< dst_yaxis[endPos] << ", Num: "<< thAxis.size()<< std::endl;
	//std::cout<< "quotient: "<< quotient<< ", modulo: "<<modulo<< ", ThID:  " << mpi_procRank << ", Start: "<< thAxis.front() << ", End: "<< thAxis.back() << ", Num: "<< thAxis.size()<< std::endl;
	//std::cout<< "quotient: "<< quotient<< ", modulo: "<<modulo<< ", ThID:  " << mpi_procRank << ", StartPos: "<< startPos << ", endPos: "<< endPos << ", Num: "<< endPos-startPos+1<< std::endl;
	src.RSprop(propz, dst_xaxis, thAxis);

	// std::string foutpath = std::to_string(mpi_procRank)+std::string("th_obs.bin");
	//src.saveWfield_bin(foutpath);
}

int main(int argc, char* argv[]){
	int xsize = 0;
	int ysize = 0;
	double pitch = 0;
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
		ysize = std::stoi(std::string(argv[2]));
		pitch = std::stod(std::string(argv[3]));
		wvl = std::stod(std::string(argv[4]));
		aptD = std::stod(std::string(argv[5]));
		propz = std::stod(std::string(argv[6]));
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

	SWPL::Wavefield wf(xsize,ysize,pitch,wvl);
	wf.setCircAparture(aptD);
	std::chrono::system_clock::time_point start;
	if(mpi_procRank==0){
		start = std::chrono::system_clock::now(); 
	}
	
	try{
		mpitest(wf,wf.getXaxis(),wf.getYaxis(),propz);
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
		std::vector<std::complex<double>> wfvec = wf.getField();
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