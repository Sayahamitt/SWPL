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

class Wavefield{
	std::shared_ptr<std::vector<std::complex<double>>> field;//1次元配列に2次元データを行優先で充填する
	//std::vector<double> xaxis;
	//std::vector<double> yaxis;
	int sizex;
	int sizey;
	double pitchx;
	double pitchy;
	double wvl;
	double wvnum;
	std::complex<double> jK;
	constexpr static std::complex<double> imunt = std::complex<double>(0,1);
public:
	class FieldParamater{
		unsigned int sizex;
		unsigned int sizey;
		double pitchx;
		double pitchy;
		double wvl;
		double wvnum;
		std::complex<double> jK;
		//std::vector<double> xaxis;
		//std::vector<double> yaxis;
	public:
		FieldParamater(unsigned int arg_sizex, unsigned int arg_sizey, double arg_pitchx, double arg_pitchy, double arg_wvl):
		sizex(arg_sizex),
		sizey(arg_sizey),
		pitchx(arg_pitchx),
		pitchy(arg_pitchy),
		wvl(arg_wvl),
		wvnum(2*SWPL_PI/arg_wvl),
		jK(imunt*(2*SWPL_PI/arg_wvl)){/*
			//X,Y軸の生成
			for(double ii=-(sizex/2);xaxis.size()<sizex; ii++){
				xaxis.push_back(ii*pitchx);
			}
			for(long long ii=-(sizey/2);yaxis.size()<sizey; ii++){
				yaxis.push_back(ii*pitchy);
			}*/
		}

		unsigned int SizeX() const{
			return sizex;
		}
		unsigned int SizeY() const{
			return sizey;
		}
		double PtithX() const{
			return pitchx;
		}
		double PtithY() const{
			return pitchy;
		}
		double Wavelength() const{
			return wvl;
		}
		double WaveNumber() const{
			return wvnum;
		}/*
		const std::vector<double> AxisX(){
			return xaxis;
		}
		const std::vector<double> AxisY(){
			return yaxis;
		}*/
	};

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
			const int cap = this->capacity();
			const double endValue = *(--(this->end()));
			if(cap<num){
				this->reserve(cap+num);
			}

			for(size_t ii=1; ii<=num; ii++){
				this->push_back(endValue+ii*pitch);
			}
		}
		void extend_front(const unsigned int& num ){
			this->extend_back(num);
			this->shift(-num*pitch);
		}

		double Pitch() const{return pitch;}
	};

private:
	FieldAxis xaxis;
	FieldAxis yaxis;

public:
	Wavefield(int xsize, int ysize, double gridspace, double wavelength):
	field(new std::vector<std::complex<double>>(ysize*xsize, std::complex<double>(0))),
	sizex(xsize),
	sizey(ysize),
	pitchy(gridspace),
	pitchx(gridspace),
	wvl(wavelength),
	wvnum(2*SWPL_PI/wavelength),
	jK(imunt*(2*SWPL_PI/wavelength)),
	xaxis(xsize,gridspace),
	yaxis(ysize,gridspace){
	/*
		//X,Y軸の生成
		for(double ii=-(xsize/2);xaxis.size()<xsize; ii++){
			xaxis.push_back(ii*pitchx);
		}
		for(long long ii=-(xsize/2);yaxis.size()<ysize; ii++){
			yaxis.push_back(ii*pitchy);
		}*/
	}

	void setCircAparture(double D,double centx=0, double centy=0){
		double r = D/2;
		for(int ii=0;ii<sizey;ii++){
			for(int jj=0;jj<sizex;jj++){
				if(std::sqrt(yaxis[ii]*yaxis[ii]+xaxis[jj]*xaxis[jj]) <= r){
					(*field)[ii*sizex+jj] = std::complex<double>(1,0);
				}else{
					(*field)[ii*sizex+jj] = std::complex<double>(0);
				}
			}
		}
	}

	void setField(const std::vector<std::complex<double>>& srcfield){
		std::vector<std::complex<double>>::iterator itr_field = (*field).begin();
		std::for_each(srcfield.begin(),srcfield.end(),[&itr_field](const std::complex<double>& fieldelement){
				*itr_field = fieldelement;
				++itr_field;
			});
	}

	const std::vector<std::complex<double>>& getField(){
		return *field;
	}

	//伝搬計算では計算後にfieldを上書きする。
	void RSprop(const double z){
		RSprop(z, xaxis, yaxis);
	}

	void RSprop(const double z, const FieldAxis distax_x, const FieldAxis distax_y){
		std::shared_ptr<std::vector<std::complex<double>>> dst_field(new std::vector<std::complex<double>>(distax_x.size()*distax_y.size(),std::complex<double>(1)));
		const double ZZ = std::pow(z,2);
		

		const Eigen::ArrayXXcd srcwf = Eigen::Map<Eigen::ArrayXXcd>((*field).data(),sizey,sizex);
		const Eigen::ArrayXXd srcax_y = Eigen::Map<Eigen::ArrayXd>(yaxis.data(),yaxis.size()).replicate(1, sizex);
		const Eigen::ArrayXXd srcax_x = Eigen::Map<Eigen::ArrayXd>(xaxis.data(),xaxis.size()).transpose().replicate(sizey, 1);

		Eigen::ArrayXXcd m_dst_field = Eigen::Map<Eigen::ArrayXXcd>((*dst_field).data(),distax_y.size(),distax_x.size());
		Eigen::ArrayXXd r01 = Eigen::ArrayXXd::Zero(distax_y.size(),distax_x.size());

		//レイリーゾンマーフェルト積分の計算
		#pragma omp parallel for private(r01)
		for(int ii=0;ii<distax_y.size();ii++){
			for(int jj=0;jj<distax_x.size();jj++){
				r01 = (ZZ+(srcax_x-distax_x[jj]).square()+(srcax_y-distax_y[ii]).square()).sqrt();
				m_dst_field(ii,jj) = (srcwf*((jK*r01).exp()/r01.square())).sum();
			}
		}
		//積分の外側の係数 z/jλ と 積分単位面積dxdy を乗算する. 入力面と観測面のピクセルサイズの比を乗算し単位面積当たりのエネルギー保存則を満たす.
		m_dst_field *= (z/(imunt*wvl))*(distax_x.Pitch()*distax_y.Pitch())*(this->xaxis.Pitch()/distax_x.Pitch())*(this->yaxis.Pitch()/distax_y.Pitch());

		Eigen::Map<Eigen::Array<std::complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>>((*dst_field).data(),sizey,sizex) = m_dst_field;

		dst_field.swap(field);

		std::cout<<"Eigen src data: "<<srcwf.data()<<std::endl;
		std::cout<<"Vector src data: "<<(*field).data()<<std::endl;
		std::cout<<"Eigen dst data: "<<m_dst_field.data()<<std::endl;
		std::cout<<"Vector dst data: "<<(*dst_field).data()<<std::endl;

		std::string path = "deb.bin";
		std::ofstream fout;
		fout.open(path.c_str(), std::ios::out|std::ios::binary|std::ios::trunc);
		for(int kk=0;kk<sizex;kk++){
			for(int ll=0;ll<sizey;ll++){
				double r = std::real((*dst_field)[kk*sizey+ll]);
				double i = std::imag((*dst_field)[kk*sizey+ll]);
				//double r = std::real(m_dst_field(kk,ll));
				//double i = std::imag(m_dst_field(kk,ll));
				fout.write((char *) &r,sizeof( double ) );
				fout.write((char *) &i,sizeof( double ) );
			}
		}
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
						distwf[ii*sizex+jj] += (z/(imunt*wvl))*(*field)[kk*sizex+ll]*(std::exp(imunt*wvnum*r01)/(r01*r01))*(pitchy*pitchx);
					}
				}
			}
		}
		*field = distwf;
	}

	FieldAxis const getXaxis(){return xaxis;}
	FieldAxis const getYaxis(){return yaxis;}

	void showWFcons(){
		int ii=0;
		std::for_each((*field).begin(),(*field).end(),[&ii,this](std::complex<double> num){
				std::cout<<num<<", "<<std::flush;
				if((ii+1)%sizex==0){
					std::cout<<std::endl;
				}ii++;
			});
	}

	void saveWfield_bin(std::string path){
		std::ofstream fout;
		fout.open(path.c_str(), std::ios::out|std::ios::binary|std::ios::trunc);
		std::for_each((*field).begin(),(*field).end(),[&fout](std::complex<double> num){
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