#include <functional>
#include <cmath>
#include <random>
#include <fstream>
#include <Eigen/Dense>
#include <string>
#include <iostream>
#include <boost/progress.hpp>

using namespace Eigen;
using namespace std;

typedef unsigned int uint;

/*
	LIST OF PARAMETERS:
	nos: #of Spins
	g: parameter that represents the velocity of the decay of the coupling
	hCentral, hBath: a priori stepsize
	givenError: error that sets the precision in the stepsize adaptation
	s0t0: the central-spin at t=0
	x: either the overhauserfield or the central-spin (depending on the used function CS/BS)
	iter: #iterations the simulation calculates
	dim: dimension of timeC/autoC, equals the maximum time value on the x axis in the plot
	counter: global position of coupling constant vector Ji
*/

namespace{
/**************************************************************************************/
double timeOfMeasurement, h, relError;
uint counter = 0;
const uint nos = 1+1e3;
const uint iter = 1e3;
const uint dim = 100*1e1;
const uint tPulse = 10;
const double g = 1e-2;
const double givenError = 1e-6;
Vector3d hCentral(0, 0, 10);
Vector3d hBath(0, 0, 0);
Vector3d s0, si, s01, B, s0t0, delta, delta1, delta2, sOrg;
VectorXd autoCx(dim), autoCy(dim), autoCz(dim), senv(dim), Ji(nos-1);
MatrixXd spinContainer(3, nos-1);
fstream data;
/**************************************************************************************/
}

Vector3d RK4(double h, Vector3d yn, Vector3d x, function<Vector3d(Vector3d, Vector3d)> f){
	
	Vector3d yn1, k1, k2, k3, k4;

	k1 = h * f(yn, x);
	k2 = h * f(yn + 0.5 * k1, x);
	k3 = h * f(yn + 0.5 * k2, x);
	k4 = h * f(yn + k3, x);
	yn1 = yn + (1./6.) * (k1 + 2.*k2 + 2.*k3 + k4);

	return yn1;
}

double stepsize(double &h, Vector3d delta, Vector3d yn1){
	
	h *= pow(givenError * yn1.norm()/delta.norm(), 0.2);

	return h;
}

//initializes one spin
Vector3d init(){
	random_device seed;
	mt19937 gen(seed());
	normal_distribution<double> distr(0., 0.5);
	Vector3d spin;
	spin(0) = distr(gen);
	spin(1) = distr(gen);
	spin(2) = distr(gen);
	return spin;
}

void expCoupling(){

        static const double A0 = sqrt((1.-exp(-2.*g))/(1.-exp(-2.*g*(nos-1))));
        for(uint i = 0; i<nos-1; i++){    
            Ji(i) = A0 * exp(-g*double(i));
        }
}

void linCoupling(){
    
        static const double A0 = sqrt(6.*(nos-1)/(2.*(nos-1)*(nos-1)+3.*(nos-1)+1.));
        for(uint i = 1; i<nos-1; i++){
            Ji(i) = A0 * (nos-double(i))/(nos-1);
        }
}

//actual eom for central- and bathspin
Vector3d CS(Vector3d s0, Vector3d B){
	
	return B.cross(s0) - hCentral.cross(s0);
}

Vector3d BS(Vector3d si, Vector3d s0){

	return Ji(counter)*(s0.cross(si)) - 0.001*hBath.cross(si);
}

void pulse(){
	
	s0(0) = s0.norm();
	s0(1) = 0.;
	s0(2) = 0.;
}

//one RK4 time-step
void timeEvol(){

	sOrg = s0;
	s01 = RK4(h, s0, B, CS);
	delta2 = RK4(2.*h, s0, B, CS);
	delta1 = RK4(h, s01, B, CS);
	relError = delta.norm()/((h*s01).norm());

	if(relError > givenError){
		h = stepsize(h, delta, h*s01);
		s0 = RK4(h, s0, B, CS);
	}
	if(relError < givenError){
		s0 = s01;
		h = (stepsize(h, delta, h*s01)<=1e-1) ? stepsize(h, delta, h*s01) : 1e-1;
	}
	else{
		s0 = s01;
	}
        B << 0, 0, 0;
	for(uint i = 0; i<nos-1; i++){
				
		counter = i;

		si(0) = spinContainer(0, i);
		si(1) = spinContainer(1, i);
		si(2) = spinContainer(2, i);
					
		si = RK4(h, si, sOrg, BS);
				
		spinContainer(0, i) = si(0);
		spinContainer(1, i) = si(1);
		spinContainer(2, i) = si(2);
	}
	B = spinContainer * Ji;
}

int main(){

        expCoupling();
	//linCoupling();

	data.open("./data/cs_iter=1e"+to_string(int(log(iter)/log(10))+1)+"_N=1e"+to_string(int(log(nos-1)/log(10))+1)+"_h0z="+to_string(int(hBath(2)))+".dat", fstream::out);
	data << "#t\t<s0_z(t)*s0_z(0)>" << '\n';

	cout << "*************************\n";
	cout << "Simulation started with:\n";
	cout << iter << "\tEnsembles\n";
	cout << nos-1 << "\tBathspins\n";
	cout << dim << "\tDatapoints\n";
	cout << "gamma = " << g << '\n';
	cout << "*************************\n";
	cout << "\nWriting data in file...\n";

	boost::progress_display show_progress(iter);
	for(uint i = 0; i<iter; i++){
					
		s0 = init();
		s0t0 = s0;
		timeOfMeasurement = 0.;
		h = 1e-2;
		B << 0, 0, 0;
		
                for(uint i = 0; i<nos-1; i++){
                    
                    si = init();
                    spinContainer(0, i) = si(0);
                    spinContainer(1, i) = si(1);
                    spinContainer(2, i) = si(2);
                }
		B = spinContainer * Ji;
		
		for(uint j = 1; j<=autoCz.rows(); j++){
			
			if((j-1)%(tPulse*10)==0){
				pulse();
			}
			
			//calculating mean of autocorrelationfunction for every component
			autoCx(j-1) += s0t0(0)*s0(0)/iter;
			autoCy(j-1) += s0t0(0)*s0(1)/iter;
			//autoCz(j-1) += s0t0(2)*s0(2)/iter;

			while(timeOfMeasurement<j*1e-1){

				timeEvol();
				timeOfMeasurement += h;
			}
		}
		++show_progress;
	}

	//calculating S_env(t)
	for(uint i = 0; i<nos-1; i++){
		autoCx(i) *= autoCx(i);
		autoCy(i) *= autoCy(i);
		senv(i) = sqrt(autoCx(i) + autoCy(i));
	}

	for(uint j = 0; j<autoCz.rows(); j++){
		data << j*0.1 << '\t' << senv(j) << '\n';
	}
	data.close();

	return 0;
}
