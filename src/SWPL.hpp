#define EIGEN_NO_DEBUG
#define SWPL_PI 3.14159265358979323846

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <complex>
#include <vector>
#include <algorithm>
#include <string>
#include <cmath>
#include <Eigen/Core>

class wavefield{
	std::vector<std::complex<double>> field;//1次元配列に2次元データを行優先で充填する
	std::vector<double> xaxis;
	std::vector<double> yaxis;
	int sizex;
	int sizey;
	double pitchx;
	double pitchy;
	double wvl;
	double wvnum;
	std::complex<double> jK;
	constexpr static std::complex<double> imunt = std::complex<double>(0,1);

public:
	wavefield(int xsize, int ysize, double gridspace, double wavelength):
	field(ysize*xsize, std::complex<double>(0)),
	sizex(xsize),
	sizey(ysize),
	pitchy(gridspace),
	pitchx(gridspace),
	wvl(wavelength),
	wvnum(2*SWPL_PI/wavelength),
	jK(imunt*(2*SWPL_PI/wavelength)){
		//X,Y軸の生成
		for(double ii=-(xsize/2);xaxis.size()<xsize; ii++){
			xaxis.push_back(ii*pitchx);
		}
		for(long long ii=-(xsize/2);yaxis.size()<ysize; ii++){
			yaxis.push_back(ii*pitchy);
		}
	}

	void setCircAparture(double D,double centx=0, double centy=0){
		double r = D/2;
		for(int ii=0;ii<sizey;ii++){
			for(int jj=0;jj<sizex;jj++){
				if(std::sqrt(yaxis[ii]*yaxis[ii]+xaxis[jj]*xaxis[jj]) <= r){
					field[ii*sizex+jj] = std::complex<double>(1,0);
				}else{
					field[ii*sizex+jj] = std::complex<double>(0);
				}
			}
		}
	}

	void setField(const std::vector<std::complex<double>>& srcfield){
		std::vector<std::complex<double>>::iterator itr_field = field.begin();
		std::for_each(srcfield.begin(),srcfield.end(),[&itr_field](const std::complex<double>& fieldelement){
				*itr_field = fieldelement;
				++itr_field;
			});
	}

	const std::vector<std::complex<double>>& getField(){
		return field;
	}

	void RSprop(const double z){
		field = RSprop(z, xaxis, yaxis).getField();
	}

	wavefield RSprop(const double z,const std::vector<double> distax_x,const std::vector<double> distax_y){
		const double dst_pitchy = distax_x[1]-distax_y[0];
		const double dst_pitchx = distax_x[1]-distax_x[0];
		const double ZZ = std::pow(z,2);
		std::vector<std::complex<double>> dst_field = std::vector<std::complex<double>>(distax_x.size()*distax_y.size());
		
		const Eigen::ArrayXXcd srcwf = Eigen::Map<Eigen::ArrayXXcd>(field.data(),sizey,sizex);
		const Eigen::ArrayXXd srcax_y = Eigen::Map<Eigen::ArrayXd>(yaxis.data(),yaxis.size()).replicate(1, sizex);
		const Eigen::ArrayXXd srcax_x = Eigen::Map<Eigen::ArrayXd>(xaxis.data(),xaxis.size()).transpose().replicate(sizey, 1);
		Eigen::ArrayXXcd e_dst_field = Eigen::ArrayXXcd::Zero(distax_y.size(),distax_x.size());
		Eigen::ArrayXXd r01 = Eigen::ArrayXXd::Zero(distax_y.size(),distax_x.size());

		#pragma omp parallel for private(r01)
		for(int ii=0;ii<distax_y.size();ii++){
			for(int jj=0;jj<distax_x.size();jj++){
				r01 = (ZZ+(srcax_x-distax_x[jj]).square()+(srcax_y-distax_y[ii]).square()).sqrt();
				e_dst_field(ii,jj) = (srcwf*((jK*r01).exp()/r01.square())).sum();
			}
		}
		e_dst_field *= (z/(imunt*wvl))*(dst_pitchy*dst_pitchx);

		Eigen::Map<Eigen::Array<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>(dst_field.data(),sizey,sizex) = e_dst_field;
		wavefield dst_wf(distax_x.size(),distax_y.size(),dst_pitchx,wvl);
		dst_wf.setField(dst_field);

		return dst_wf;
	}

	void RSprop_nosimd(double z){
		std::vector<std::complex<double>> distwf(sizey*sizex, std::complex<double>(0));
		std::vector<double> distax_x = xaxis;
		std::vector<double> distax_y = yaxis;

		#pragma omp parallel for
		for(int ii=0;ii<sizey;ii++){
			for(int jj=0;jj<sizex;jj++){
				for(int kk=0;kk<sizey;kk++){
					for(int ll=0;ll<sizex;ll++){
						double r01 = std::sqrt(std::pow(z,2) + std::pow(distax_x[jj]-xaxis[ll],2) + std::pow(distax_y[ii]-yaxis[kk],2));
						distwf[ii*sizex+jj] += (z/(imunt*wvl))*field[kk*sizex+ll]*(std::exp(imunt*wvnum*r01)/(r01*r01))*(pitchy*pitchx);
					}
				}
			}
		}
		field = distwf;
	}

	void showWFcons(){
		int ii=0;
		std::for_each(field.begin(),field.end(),[&ii,this](std::complex<double> num){
				std::cout<<num<<", "<<std::flush;
				if((ii+1)%sizex==0){
					std::cout<<std::endl;
				}ii++;
			});
	}

	void saveWfield_bin(std::string path){
		std::ofstream fout;
		fout.open(path.c_str(), std::ios::out|std::ios::binary|std::ios::trunc);
		std::for_each(field.begin(),field.end(),[&fout](std::complex<double> num){
				double r = std::real(num);
				double i = std::imag(num);
				fout.write((char *) &r,sizeof( double ) );
				fout.write((char *) &i,sizeof( double ) );
			});
		fout.close();
	}

	void savevectord_bin(std::vector<double> vec,std::string path){
		std::ofstream fout;
		fout.open(path.c_str(), std::ios::out|std::ios::binary|std::ios::trunc);

		std::for_each(vec.begin(),vec.end(),[&fout](double num){
				double r = num;
				fout.write((char *) &r,sizeof( double ) );
			});
		fout.close();
	}
};