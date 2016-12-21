/*
	ES FEHLT DAS SCHREIBEN DER AUSGABEN IN EINE .dat DATEI
	ES FEHLEN DIE TATSÄCHLICHGEN DGLN
	ES FEHLEN KONKRETE WERTE FÜR ZEITSCHRITTE UND FEHLER SOWIE MAXIMALE SCHRITTWEITE
 	DIE INITIALISIERTEN BADSPINS PRO ZEITSCHRITT IN EINE MATRIX ALS SPALTEN EINTRAGEN
*/

#include <omp.h>
#include <Eigen/Eigenvalues>
#include <functional>
#include <cmath>
#include <random>

using namespace Eigen;
using namespace std;

VectorXd RK4(double h, VectorXd yn, function<VectorXd(VectorXd)> f){
	
	VectorXd yn1(3), k1(3), k2(3), k3(3), k4(3);
	
	k1 = h * f(yn);
	k2 = h * f(yn + 0.5 * k1);
	k3 = h * f(yn + 0.5 * k2);
	k4 = h * f(yn + k3);
	yn1 = yn + (1/6) * (k1 + 2*k2 + 2*k3 + k4);

	return yn1;
}

double stepsize(double h, double givenError, VectorXd deltaY, VectorXd yn1){
	
	h *= pow(((givenError * yn1.norm())/deltaY.norm()), 0.2);

	return h;
}

VectorXd init(){
	random_device seed;
	mt19937_64 gen(seed);
	normal_distribution<double> distr(0., 1.);
	VectorXd bathspin(3);
	bathspin << distr(gen), distr(gen), distr(gen);
	return bathspin;
}

//includes the external magnetic field h0
VectorXd BS(VectorXd si){
	return si; //left side of the CS-DE for the bath-spins
}

int main(){
	//N is the number of steps the simulation has to calculate
	double hAdaptive = 0, givenError, hMax, h, n; 
	unsigned int nos; //nos := Number of (Bath-)Spins
	VectorXd yn1(3), yn(3);

	//simulation of all bathspins
	for(int j = 0; j<nos; j++){
		//simulation of one bathspin
		for(int i = 0; i<n; i++){
			//one RK4 time-step for one bathspin
			yn1 = RK4(h, yn, BS);
			hAdaptive = stepsize(h, givenError, RK4(2*h, yn, BS) - yn1, yn1);
	
			if(hAdaptive < h){
				yn = RK4(hAdaptive, yn, BS);
			}
			if(hAdaptive > h){
				yn = yn1;
				h = hAdaptive;
			}
		}
	}

	return 0;
}
