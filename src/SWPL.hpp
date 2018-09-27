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
#include <memory>
#include <Eigen/Core>

namespace SWPL{

	class FieldAxis:public std::vector<double>{
		double pitch;
	public:
		FieldAxis(int size, double arg_pitch):
		std::vector<double>(),
		pitch(arg_pitch){
			//軸の生成
			for(double ii=-(size/2);this->size()<size; ii++){
				this->push_back(ii*arg_pitch);
			}
		}
		FieldAxis(double lower, double upper, double arg_pitch):
		std::vector<double>(),
		pitch(arg_pitch){
			//軸の生成
			for(double value=lower;value<=upper; value+=arg_pitch){
				push_back(value);
				if(value==upper){
				}
			}
		}
		FieldAxis(double lower, double arg_pitch, int size):
		std::vector<double>(),
		pitch(arg_pitch){
			//軸の生成
			for(double ii=0; ii<size; ii++){
				push_back(lower+ii*arg_pitch);
			}
		}

		void swap(FieldAxis& obj){
			double tmp_objPitch = obj.Pitch();
			obj.pitch = pitch;
			pitch = tmp_objPitch;

			std::vector<double>::swap(obj);
		}

		void shift(double deltashift){
			for(std::vector<double>::iterator itr=this->begin();itr<this->end();++itr){
				*itr+=deltashift;
			}
		}
		void scale(double ratio){
			this->pitch*=ratio;
			for(std::vector<double>::iterator itr=this->begin();itr<this->end();++itr){
				*itr*=ratio;
			}
		}
		void extend_back(const size_t& num){
			const int cap = capacity();
			const double endValue = back();
			if(cap<num){
				reserve(cap+num);
			}

			for(size_t ii=1; ii<=num; ii++){
				push_back(endValue+ii*pitch);
			}
		}
		void extend_front(const size_t& num ){
			extend_back(num);
			shift(-int(num)*pitch);
		}

		double Pitch() const{return pitch;}
	};


	class Wavefield{
		std::vector<std::complex<double>> field;//1次元配列に2次元データを行優先で充填する
		double wvl;
		double wvnum;
		std::complex<double> jK;
		FieldAxis xaxis;
		FieldAxis yaxis;
	public:
		static const std::complex<double> imunt;
	
		Wavefield(int xsize, int ysize, double gridspace, double wavelength):
		field(std::vector<std::complex<double>>(ysize*xsize, std::complex<double>(0))),
		wvl(wavelength),
		wvnum(2*SWPL_PI/wavelength),
		jK(imunt*(2*SWPL_PI/wavelength)),
		xaxis(xsize,gridspace),
		yaxis(ysize,gridspace){
		}
		Wavefield(std::vector<std::complex<double>> arg_field, FieldAxis arg_xaxis, FieldAxis arg_yaxis, double wavelength):
		field(arg_field),
		wvl(wavelength),
		wvnum(2*SWPL_PI/wavelength),
		jK(imunt*(2*SWPL_PI/wavelength)),
		xaxis(arg_xaxis),
		yaxis(arg_yaxis){
		}
		//Copy Constructor
		Wavefield(const Wavefield& rhs):
		field(rhs.field),
		wvl(rhs.wvl),
		wvnum(rhs.wvnum),
		jK(rhs.jK),
		xaxis(rhs.xaxis),
		yaxis(rhs.yaxis){}
	
		int sizex() const{return xaxis.size();}
		int sizey() const{return yaxis.size();}
		double pitchx() const{return xaxis.Pitch();}
		double pitchy() const{return yaxis.Pitch();}
		double wavelength() const{return wvl;}
		void setWavelength(double wavelength){wvl = wavelength;}
		std::vector<std::complex<double>> getField() const{return field;}
	
		void setCircAparture(double diamater,double centx=0, double centy=0){
			double r = diamater/2;
			for(int ii=0;ii<sizey();ii++){
				for(int jj=0;jj<sizex();jj++){
					if(std::sqrt(yaxis[ii]*yaxis[ii]+xaxis[jj]*xaxis[jj]) <= r){
						field[ii*sizex()+jj] = std::complex<double>(1,0);
					}else{
						field[ii*sizex()+jj] = std::complex<double>(0);
					}
				}
			}
		}
	
		//インスタンスの複素振幅分布から任意の点における光波の複素振幅を求める関数
		std::complex<double> RSpoint(const double z, const double dst_x, const double dst_y, const double dst_xPitch, const double dst_yPitch){
			const double ZZ = std::pow(z,2);
	
			const Eigen::ArrayXXcd srcwf = Eigen::Map<Eigen::ArrayXXcd>(field.data(),sizey(),sizex());
			const Eigen::ArrayXXd srcax_y = Eigen::Map<Eigen::ArrayXd>(yaxis.data(),yaxis.size()).replicate(1, sizex());
			const Eigen::ArrayXXd srcax_x = Eigen::Map<Eigen::ArrayXd>(xaxis.data(),xaxis.size()).transpose().replicate(sizey(), 1);
	
			Eigen::ArrayXXd r01 = Eigen::ArrayXXd::Zero(yaxis.size(),xaxis.size());
	
			//レイリーゾンマーフェルト積分の計算
			r01 = (ZZ+(srcax_x-dst_x).square()+(srcax_y-dst_y).square()).sqrt();
			std::complex<double> U0 = (srcwf*((jK*r01).exp()/r01.square())).sum();
			U0 *= (z/(imunt*wvl))*(dst_xPitch*dst_yPitch)*(this->xaxis.Pitch()/dst_xPitch)*(this->yaxis.Pitch()/dst_yPitch);
	
			return U0;
		}
	
		std::complex<double> RSpoint(const double& z, const double& dst_x, const double& dst_y, const double& dst_xPitch, const double& dst_yPitch, const Eigen::ArrayXXcd& srcwf, const Eigen::ArrayXXd& srcax_y, const Eigen::ArrayXXd& srcax_x) const{
			Eigen::ArrayXXd r01 = Eigen::ArrayXXd::Zero(yaxis.size(),xaxis.size());
	
			//レイリーゾンマーフェルト積分の計算
			r01 = (std::pow(z,2)+(srcax_x-dst_x).square()+(srcax_y-dst_y).square()).sqrt();
			std::complex<double> U0 = (srcwf*((jK*r01).exp()/r01.square())).sum();
			U0 *= (z/(imunt*wvl))*(dst_xPitch*dst_yPitch)*(this->xaxis.Pitch()/dst_xPitch)*(this->yaxis.Pitch()/dst_yPitch);
	
			return U0;
		}
	
		//伝搬先スクリーンの計算後には計算後にfieldを上書きする。
		void RSprop(const double z){
			RSprop(z, xaxis, yaxis);
		}
	
		void RSprop(const double z, const FieldAxis& distax_x, const FieldAxis& distax_y){
			std::vector<std::complex<double>> dst_field(distax_x.size()*distax_y.size(),std::complex<double>(0));
			const double ZZ = std::pow(z,2);
	
			const Eigen::ArrayXXcd srcwf = Eigen::Map<Eigen::ArrayXXcd>(field.data(),sizey(),sizex());
			const Eigen::ArrayXXd srcax_y = Eigen::Map<Eigen::ArrayXd>(yaxis.data(),yaxis.size()).replicate(1, sizex());
			const Eigen::ArrayXXd srcax_x = Eigen::Map<Eigen::ArrayXd>(xaxis.data(),xaxis.size()).transpose().replicate(sizey(), 1);
	
			Eigen::ArrayXXcd m_dst_field = Eigen::ArrayXXcd::Zero(distax_y.size(),distax_x.size());
			Eigen::ArrayXXd r01 = Eigen::ArrayXXd::Zero(yaxis.size(),xaxis.size());
	
			//レイリーゾンマーフェルト積分の計算
			const int dsize_x = distax_x.size();
			const int dsize_y = distax_y.size();
			#ifndef SWPL_UNUSE_OPENMP
				#pragma omp parallel for private(r01)
			#endif
			for(int ii=0;ii<dsize_y;ii++){
				for(int jj=0;jj<dsize_x;jj++){
					r01 = (ZZ+(srcax_x-distax_x[jj]).square()+(srcax_y-distax_y[ii]).square()).sqrt();
					m_dst_field(ii,jj) = (srcwf*((jK*r01).exp()/r01.square())).sum();
					//dst_field[dsize_x*ii+jj] = (srcwf*((jK*r01).exp()/r01.square())).sum();
				}
			}
			//積分の外側の係数 z/jλ と 積分単位面積dxdy を乗算する. 入力面と観測面のピクセルサイズの比を乗算し単位面積当たりのエネルギー保存則を満たす.
			m_dst_field *= (z/(imunt*wvl))*(distax_x.Pitch()*distax_y.Pitch())*(this->xaxis.Pitch()/distax_x.Pitch())*(this->yaxis.Pitch()/distax_y.Pitch());
			#ifndef SWPL_UNUSE_OPENMP
				#pragma omp parallel for
			#endif
			for(int ii=0; ii<dst_field.size(); ii++){
				dst_field[ii] *= (z/(imunt*wvl))*(distax_x.Pitch()*distax_y.Pitch())*(this->xaxis.Pitch()/distax_x.Pitch())*(this->yaxis.Pitch()/distax_y.Pitch());	
			}
			Eigen::Map<Eigen::Array<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>(dst_field.data(),dsize_y,dsize_x) = m_dst_field;
			field.swap(dst_field);
			xaxis=distax_x;
			yaxis=distax_y;
		}
	
		void RSprop_nosimd(double z){
			std::vector<std::complex<double>> distwf(sizey()*sizex(), std::complex<double>(0));
			std::vector<double> distax_x = xaxis;
			std::vector<double> distax_y = yaxis;
	
			#ifndef SWPL_UNUSE_OPENMP
				#pragma omp parallel for
			#endif
			for(int ii=0;ii<sizey();ii++){
				for(int jj=0;jj<sizex();jj++){
					for(int kk=0;kk<sizey();kk++){
						for(int ll=0;ll<sizex();ll++){
							double r01 = std::sqrt(std::pow(z,2) + std::pow(distax_x[jj]-xaxis[ll],2) + std::pow(distax_y[ii]-yaxis[kk],2));
							distwf[ii*sizex()+jj] += (z/(imunt*wvl))*field[kk*sizex()+ll]*(std::exp(imunt*wvnum*r01)/(r01*r01))*(pitchy()*pitchx());
						}
					}
				}
			}
			field = distwf;
		}
	
		FieldAxis getXaxis() const{return xaxis;}
		FieldAxis getYaxis() const{return yaxis;}
	
		void showWFcons(){
			int ii=0;
			std::for_each(field.begin(),field.end(),[&ii,this](std::complex<double> num){
				std::cout<<num<<", "<<std::flush;
				if((ii+1)%sizex()==0){
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
	
		void saveCompDouble_bin(std::vector<std::complex<double>> vec,std::string path){
			std::ofstream fout;
			fout.open(path.c_str(), std::ios::out|std::ios::binary|std::ios::trunc);
			std::for_each(vec.begin(),vec.end(),[&fout](std::complex<double> num){
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

//static const メンバの初期化
const std::complex<double> Wavefield::imunt = std::complex<double>(0,1);
}