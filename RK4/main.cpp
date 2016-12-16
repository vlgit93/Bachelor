/*
	ES FEHLT DAS SCHREIBEN DER AUSGABEN IN EINE .dat DATEI
	ES FEHLEN DIE TATSÄCHLICHGEN DGLN
	ES FEHLEN KONKRETE WERTE FÜR ZEITSCHRITTE UND FEHLER SOWIE MAXIMALE SCHRITTWEITE
	ES FEHLT DIE IMPLEMENTIERUNG EINER WAHRSCHEINLICHKEITSVERTEILUNG FÜR DIE ANFANGSBEDINGUNGEN (GAUẞ-VERTEILUNG)
*/

#include <omp.h>
#include <vector>
#include <functional>
#include <cmath>

using namespace std;

vector<double> vectorProduct(vector<double> x, vector<double> y){

	vector<double> z;

	z[0] = x[1] * y[2] - x[2] * y[1];
	z[1] = x[0] * y[2] - x[2] * y[0];
	z[2] = x[0] * y[1] - x[1] * y[0];
	
	return z;
}

inline vector<double> scalarVectorProd(vector<double> x, double a){
	
	for(int i = 0; i<3; i++){
		x[i] *= a;
	}

	return x;
}

vector<double> RK4(double h, vector<double> yn, function<vector<double>(vector<double>)> f){
	
	vector<double> yn1, k1, k2, k3, k4;
	
	k1 = scalarVectorProd(f(yn), h);
	for(int i = 0; i<3; i++){
		yn[i] += 0.5 * k1[i];
	}
	k2 = scalarVectorProd(f(yn), h);
	for(int i = 0; i<3; i++){
		yn[i] += 0.5 * k2[i];
	}
	k3 = scalarVectorProd(f(yn), h);
	for(int i = 0; i<3; i++){
		yn[i] += k3[i];
	}
	k4 = scalarVectorProd(f(yn), h);
	
	for(int i = 0; i<3; i++){
		yn1[i] = yn[i] + (1/6) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
	}

	return yn1;
}

double stepsize(double h, double givenError, vector<double> deltaY, vector<double> yn1, function<vector<double>(vector<double>)> f){
	
	double absDeltaY, absY = 0;
	for(int i = 0; i<3; i++){
		absDeltaY += deltaY[i]*deltaY[i];
		absY += f(yn1)[i]*f(yn1)[i];
	}
	absDeltaY = sqrt(absDeltaY);
	absY = sqrt(absY);
	h *= pow(((givenError * absY)/absDeltaY), 0.2);

	return h;
}

/*
//includes the external magnetic field h0
vector<double> BS(vector<double> si){
	return ; //left side of the CS-DE for the bath-spins
}
*/

int main(){
	/*
	double hAdaptive = 0, givenError, hMax, h, N; //N is the number of steps the simulation has to calculate
	vector<double> yn1;

	//actual simulation
	for(int i = 0; i<N; i++){
		yn1 = RK4(h, yn, BS);
		hAdaptive = stepzise(h, givenError, RK4(2*h, yn, BS) - yn1, yn1, BS);

		if(hAdaptive < h){
			RK4(h, yn, BS);
		}
		if(hAdaptive > h){
			RK4(h, yn, BS);
			h = hAdaptive;
		}
	}

	*/
	return 0;
}
