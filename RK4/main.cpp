#include <functional>
#include <cmath>
#include <random>
#include <fstream>
#include <Eigen/Dense>
#include <string>
#include <omp.h>

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
*/

namespace{
/**************************************************************************************/
double timeOfMeasurement = 0.;
double hCentral, hBath, relError; 
const uint nos = 1e2+1;
const uint iter = 1e4;
const uint dim = 10*1e2;
const double g = 3e-2;
const double givenError = 1e-6;
Vector3d s0, si, s01, si1, B, s0t0, delta, save;
Vector3d h0(0, 0, 0);
VectorXd timeC(dim), autoC(dim);
MatrixXd spinContainer(3, nos-1);
fstream data;
/**************************************************************************************/
}

Vector3d RK4(double h, Vector3d yn, Vector3d x, function<Vector3d(Vector3d, Vector3d)> f){
	
	Vector3d yn1, k1, k2, k3, k4;
	
	if(hCentral==h){
		save = f(yn, x);
		k1 = h * save;
	}
	else{
		k1 = h * save;
	}
	//k1 = h * f(yn, x);
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

inline double coupling(uint i){

	return sqrt((1.-exp(-2.*g))/(1.-exp(-2.*g*nos))) * exp(-g*(i-1));
}

Vector3d CS(Vector3d s0, Vector3d B){
	
	return B.cross(s0) - h0.cross(s0);
}

Vector3d BS(Vector3d s0, Vector3d si){

	return s0.cross(si) - h0.cross(si);
}

//one RK4 time-step
void timeEvol(double &h, Vector3d &s, Vector3d &s1, Vector3d x, function<Vector3d(Vector3d, Vector3d)> f){

	s1 = RK4(h, s, x, f);
	delta = RK4(2.*h, s, x, f) - RK4(h, s1, x, f);
	relError = delta.norm()/((h*s1).norm());

	if(relError > givenError){
		h = stepsize(h, delta, h*s1);
		s = RK4(h, s, x, f);
	}
	if(relError < givenError){
		s = s1;
		h = (stepsize(h, delta, h*s1)<=1e-2) ? stepsize(h, delta, h*s1) : 1e-2;
	}
	else{
		s = s1;
	}
}

int main(){

	data.open("./data/cs_"+to_string(iter)+".dat", fstream::out);
	data << "#t\t<s0_z(t)*s0_z(0)>" << '\n';

	#pragma omp parallel for
	for(uint i = 0; i<iter; i++){
					
		s0 = init();
		s0t0 = s0;
		timeOfMeasurement = 0.;
		hCentral = 1e-2;

		for(uint i = 0; i<nos-1; i++){
		
			si = init();
			B += coupling(i+1) * si;
		
			spinContainer(0, i) = si(0);
			spinContainer(1, i) = si(1);
			spinContainer(2, i) = si(2);
		}
		
		for(uint j = 1; j<=timeC.rows(); j++){
			
			timeC(j-1) += timeOfMeasurement;
			autoC(j-1) += s0t0(2)*s0(2);

			while(timeOfMeasurement < j*1e-1){
				timeEvol(hCentral, s0, s01, B, CS);
				/*
				B << 0, 0, 0; 
				//one RK4 time-step for all bath-spins
				for(uint i = 0; i<nos-1; i++){
		
					si(0) = spinContainer(0, i);
					si(1) = spinContainer(1, i);
					si(2) = spinContainer(2, i);
		
					timeEvol(hBath, si, si1, s0, BS);
		
					B += coupling(i+0) * si;
					
					spinContainer(0, i) = si(0);
					spinContainer(1, i) = si(1);
					spinContainer(2, i) = si(2);
				}
				*/
				timeOfMeasurement += hCentral;
			}
		}
	}

	#pragma omp parallel for
	for(uint j = 0; j<timeC.rows(); j++){
		data << timeC(j)/iter << '\t' << autoC(j)/iter << '\n';
	}
	data.close();

	return 0;
}
