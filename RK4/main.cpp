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
	h, hBath: a priori stepsize
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
const uint nos = 1+1e2;
const uint iter = 1e1;
const uint dim = 100*1e1;
const double g = 1e-2;
const double givenError = 1e-6;
Vector3d h0(0, 0, 0);
Vector3d s0, si, s1, sOrg, B, s0t0, delta;
VectorXd timeC(dim), autoC0(dim), autoC1(dim), autoC2(dim), Ji(nos-1);
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
    
        static const double A0 = sqrt(6*(nos-1)/(2*(nos-1)*(nos-1)+3*(nos-1)+1));
        for(uint i = 1; i<nos-1; i++){
            Ji(i) = A0 * (nos-double(i))/(nos-1);
        }
}

//actual eom for central- and bathspin
Vector3d CS(Vector3d s0, Vector3d B){
	
	return B.cross(s0) - h0.cross(s0);
}

Vector3d BS(Vector3d si, Vector3d s0){

	return Ji(counter)*(s0.cross(si)) - 0.001*h0.cross(si);
}

void pulse(Vector3d s){
	
	s(0) = s0.norm();
	s(1) = 0;
	s(2) = 0;
}

//one RK4 time-step
void timeEvol(){

	sOrg = s0;
	s1 = RK4(h, s0, B, CS);
	delta = RK4(2.*h, s0, B, CS) - RK4(h, s1, B, CS);
	relError = delta.norm()/((h*s1).norm());

	if(relError > givenError){
		h = stepsize(h, delta, h*s1);
		s0 = RK4(h, s0, B, CS);
	}
	if(relError < givenError){
		s0 = s1;
		h = (stepsize(h, delta, h*s1)<=1e-2) ? stepsize(h, delta, h*s1) : 1e-2;
	}
	else{
		s0 = s1;
	}
        
	B << 0, 0, 0;
	for(uint i = 0; i<nos-1; i++){
			
		counter = i;

		si(0) = spinContainer(0, i);
		si(1) = spinContainer(1, i);
		si(2) = spinContainer(2, i);
		
		si = RK4(h, si, sOrg, BS);
		B += Ji(i) * si;
				
		spinContainer(0, i) = si(0);
		spinContainer(1, i) = si(1);
		spinContainer(2, i) = si(2);
	}

}

int main(){

        expCoupling();

	data.open("./data/cs_iter=1e"+to_string(int(log(iter)/log(10)))+"_N=1e"+to_string(int(log(nos-1)/log(10)))+"_h0z=0"+".dat", fstream::out);
	//data << "#t\t<s0_x(t)*s0_x(0)>\t<s0_y(t)*s0_y(0)>\t<s0_z(t)*s0_z(0)>" << '\n';
	data << "#t\t<s0_z(t)*s0_z(0)>" << '\n';
	
	cout << "************************\n";
	cout << "Simulation started with:\n";
	cout << iter << "\tEnsembles\n";
	cout << nos-1 << "\tBathspins\n";
	cout << dim << "\tDatapoints\n";
	cout << "gamma = " << g << '\n';
	cout << "Magneticfield: (" << h0(0) << ", " << h0(1) << ", " << h0(2) << ")\n";
	cout << "************************\n";
        cout << "\nWriting data in file...\n";

	boost::progress_display show_progress(iter);	
        for(uint i = 0; i<iter; i++){
					
		s0 = init();
		s0t0 = s0;
		timeOfMeasurement = 0.;
		h = 1e-3;
		B << 0, 0, 0;
		
                for(uint i = 0; i<nos-1; i++){
                    
                    si = init();
                    spinContainer(0, i) = si(0);
                    spinContainer(1, i) = si(1);
                    spinContainer(2, i) = si(2);
                    B += Ji(i) * si;
                }
		
		for(uint j = 1; j<=timeC.rows(); j++){
			
			timeC(j-1) += timeOfMeasurement;
			//autoC0(j-1) += s0t0(0)*s0(0);
			//autoC1(j-1) += s0t0(1)*s0(1);
			autoC2(j-1) += s0t0(2)*s0(2);

			while(timeOfMeasurement<j*1e-2){
				/*
				if(timeOfMeasurement==((j-1)*10)){
					pulse(s0);
				}
				*/
				timeEvol();
				timeOfMeasurement += h;
			}
		}
		++show_progress;	
	}

	for(uint j = 0; j<timeC.rows(); j++){
		//data << timeC(j)/iter << '\t' << autoC0(j)/iter << '\t' << autoC1(j)/iter << '\t' << autoC2(j)/iter << '\n';
		data << timeC(j)/iter << '\t' << autoC2(j)/iter << '\n';
	}
	data.close();

	return 0;
}
