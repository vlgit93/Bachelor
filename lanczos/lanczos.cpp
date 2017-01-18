/*
 * lanczos.cpp
 *
 *  Created on: 01.12.2015
 *      Author: huedepohl
 *
 *  Modified on: 17.11.2016
 *      Author: schering
 */

#include <float.h>

#define DEBUG 0
#include <fstream>
#include <iostream>
#define _USE_MATH_DEFINES
#include <boost/program_options.hpp>
#include <cassert>
#include <chrono>
#include <functional>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <iomanip>
#include <math.h>
#include <random>
#include <sstream>
#include <string>

using namespace std;
namespace po = boost::program_options;


unsigned seed1 = (unsigned)chrono::system_clock::now().time_since_epoch().count();
// default_random_engine generator;
// random_device rd; //Seed fuer Mersenne-Twister
mt19937 gen(10); // Mersenne-Twister
// uniform_real_distribution<double> gleich(0, 1); //Gleichverteilung
normal_distribution<double> normal(0.0, 1. / 2.); // Gaussverteilung
// mu	sigma

struct parameter {
    int level;      // level of approx.
    double B;       // B-Field value
    int B_function; // function of B-field
    double* alpha;
    double* beta;
};

double Bfield(const double t, const double B0, const int function)
{
    assert(function > 0);
    switch(function) {
        case 1:
            return B0; // B=const
    }
}

void out(const double* v, const int v_size)
{
    cout << endl;
    for(int i = 0; i < v_size; i++) {
        cout << v[i] << endl;
    }
}

double norm(const double* v, const int v_size)
{
    double result = 0.;
    for(int i = 0; i < v_size; i++) {
        result += v[i] * v[i];
    }
    return sqrt(result);
}

void couplings_gamma(double* J, const int J_size, const double gamma)
{
    assert(gamma < 1);
    assert(J_size > 1);
    const double A0 = sqrt((1. - exp(-2. * gamma)) / (1. - exp(-2. * gamma * J_size)));
    for(int i = 0; i < J_size; i++) {
        J[i] = A0 * exp(-double(i) * gamma);
    }
}

void couplings_posneg(double* J, const int J_size, const double gamma)
{
    assert(gamma < 1);
    assert(J_size % 2 == 0);
    assert(J_size > 1);
    const double A0 = sqrt((1. - exp(-2. * gamma)) / (1. - exp(-gamma * J_size)));
    for(int i = 0; i < J_size / 2; i++) {
        J[i] = A0 / sqrt(2) * exp(-double(i) * gamma);
    }
    for(int i = 0; i < J_size / 2; i++) {
        J[i + J_size / 2] = -A0 / sqrt(2) * exp(-double(i) * gamma);
    }
}

void couplings_equidist(double* J, const int N)
{
    for(int i = 0; i < N; i++) {
        J[i] = sqrt(6. * N / (2. * N * N + 3. * N + 1.)) * (N - i) / N;
    }
}

// evaluate polynom p at x with horner scheme
double poly_val(const double* p, const int p_size, const double x)
{
    if(p_size == 2) return p[1] * x + p[0];
    double res = x * p[p_size - 1];
    for(int i = p_size - 2; i >= 1; i--) {
        res += p[i];
        res *= x;
    }
    return res + p[0];
}

// vectorproduct for vectors of arbitrary length with given indices
inline void crossprod(const double* v1, const double* v2, double* dest, const int ind1,
                      const int ind2, const int inddest)
{
    dest[inddest]     = v1[ind1 + 1] * v2[ind2 + 2] - v1[ind1 + 2] * v2[ind2 + 1];
    dest[inddest + 1] = v1[ind1 + 2] * v2[ind2] - v1[ind1] * v2[ind2 + 2];
    dest[inddest + 2] = v1[ind1] * v2[ind2 + 1] - v1[ind1 + 1] * v2[ind2];
}


// right side of hierarchy with lancz alg. for RK4
int rightside(double t, const double* y0, double* dest, void* params)
{
    parameter* p             = (parameter*)params;
    static int lvl           = p->level;
    double H                 = p->B;
    double HBath             = 0; // 0.001 * H;
    static double* alpha     = p->alpha;
    static double* beta      = p->beta;
    static double B_function = p->B_function;
    H                        = Bfield(t, H, B_function);
    int y0_size;
    lvl ? y0_size   = 3 * lvl + 3 : 6;
	static int lvl3 = 3 * lvl;
	int i3;
    double* y0_temp = new double[y0_size];
	double* P_temp = new double[3];

    // dynamic of CS around P1
    //crossprod(y0, y0, dest, 3, 0, 0);
    dest[0] = y0[4] * y0[2] - y0[5] * y0[1];
    dest[1] = y0[5] * y0[0] - y0[3] * y0[2];
    dest[2] = y0[3] * y0[1] - y0[4] * y0[0];
    // contribution of b field
    dest[0] += H * y0[1];
    dest[1] -= H * y0[0];
    /*switch(lvl) {
        case 0:
            dest[3] = 0;
            dest[4] = 0;
            dest[5] = 0;
            break;
        case 1:
            for(unsigned i = 0; i < y0_size; ++i) {
                y0_temp[i] = alpha[0] * y0[i];
            }
            crossprod(y0, y0_temp, dest, 0, 3, 3);
            break;
        default:*/
			P_temp[0] = beta[0] * y0[6]  + alpha[0] * y0[3];
			P_temp[1] = beta[0] * y0[7]  + alpha[0] * y0[4];
			P_temp[2] = beta[0] * y0[8]  + alpha[0] * y0[5];
			//crossprod(y0, P_temp, dest, 0, 0, 3);
			dest[3] = y0[1] * P_temp[2] - y0[2] * P_temp[1];
			dest[4] = y0[2] * P_temp[0] - y0[0] * P_temp[2];
			dest[5] = y0[0] * P_temp[1] - y0[1] * P_temp[0];
                        
                        dest[3] += HBath * y0[4];
                        dest[4] -= HBath * y0[3];
			
			for(unsigned i = 2; i < lvl; i++) {
					i3 = 3 * i;
					P_temp[0] = beta[i - 1] * y0[i3 + 3] + alpha[i - 1] * y0[i3] + beta[i - 2] * y0[i3 - 3];
					P_temp[1] = beta[i - 1] * y0[i3 + 4] + alpha[i - 1] * y0[i3 + 1] + beta[i - 2] * y0[i3 - 2];
					P_temp[2] = beta[i - 1] * y0[i3 + 5] + alpha[i - 1] * y0[i3 + 2] + beta[i - 2] * y0[i3 - 1];
					//crossprod(y0, P_temp, dest, 0, 0, 3 * i);
					dest[i3] = y0[1] * P_temp[2] - y0[2] * P_temp[1];
					dest[i3 + 1] = y0[2] * P_temp[0] - y0[0] * P_temp[2];
					dest[i3 + 2] = y0[0] * P_temp[1] - y0[1] * P_temp[0];
                                        
                                        dest[i3] += HBath * y0[i3 + 1];
                                        dest[i3 + 1] -= HBath * y0[i3];
			}			
			
			P_temp[0] = alpha[lvl - 1] * y0[lvl3] + beta[lvl - 2] * y0[lvl3 - 3];
			P_temp[1] = alpha[lvl - 1] * y0[lvl3 + 1] + beta[lvl - 2] * y0[lvl3 - 2];
			P_temp[2] = alpha[lvl - 1] * y0[lvl3 + 2] + beta[lvl - 2] * y0[lvl3 - 1];
			//crossprod(y0, P_temp, dest, 0, 0, 3 * lvl);
			dest[lvl3] = y0[1] * P_temp[2] - y0[2] * P_temp[1];
			dest[lvl3 + 1] = y0[2] * P_temp[0] - y0[0] * P_temp[2];
			dest[lvl3 + 2] = y0[0] * P_temp[1] - y0[1] * P_temp[0];
                        
                        dest[lvl3] += HBath * y0[lvl3 + 1];
                        dest[lvl3 + 1] -= HBath * y0[lvl3];
    //}

    delete[] y0_temp;
	delete[] P_temp;

    return GSL_SUCCESS;
}

// calc lanczos coeff. alpha beta
// void lanczos(lanczosvals& res,const int lvl,const vector<double>&J,const double gamma){
//	double numzero=1e-12; //lower limit for lanczos
//	//base to orthogonalize
//	vector<vector<double>> p(lvl+1,vector<double>());
//	//vectors for lanczos
//	vector<double> p1;
//	vector<double> p2;
//	p[0] = vector<double>(2, 0.0);
//	//coefficients
//	vector<double> alpha(lvl,0.);
//	vector<double> beta(lvl,0.);
//	//P1=a
//	p[0][1] = 1.0;
//
//	//lanczos algorithm
//	for(unsigned n=0;n<lvl;n++){
//		p1=Poly_mult_x(p[n]);
//		alpha[n]=scalarprod(p1,p[n],J,gamma);
//		if(n==0) p2=p1-alpha[n]*p[n];
//		else p2=p1-alpha[n]*p[n]-beta[n-1]*p[n-1];
//		beta[n]=poly_norm(p2,J,gamma);
//		if(fabs(beta[n])>numzero){
//			p[n+1]=1./beta[n]*p2;
//		}
//		else{
//			break;
//		}
//		GramSchmidt(p,J,n,gamma);
//	}
//	res.alpha=alpha;
//	res.beta=beta;
//}

void flip_spin(double* y0, const int dir) //Pulse 2
{
    // flip spin in direction dir  --> pulse
    double spinval = sqrt(y0[0] * y0[0] + y0[1] * y0[1] + y0[2] * y0[2]);
    y0[0]          = 0.;
    y0[1]          = 0.;
    y0[2]          = 0.;
    y0[dir]        = spinval;
}

void pi_pulse(double* y0)
{
    y0[0] = -y0[0];
    y0[2] = -y0[2];
}

void abs_pulse(double* y0, const int dir) //Pulse 1
{
    y0[dir] = abs(y0[dir]);
}

void around_z_pulse(double* y0)
{
    if(y0[0] < 0) {
        y0[0] = -y0[0];
        y0[1] = -y0[1];
    }
}

double hierarchy(const int N, const int max_configs, const int lvl, const int tmax,
                 const int coupling, const float gamma, const string name, const int datapoints,
                 const double H, const int B_function, const double t_pulse)
{
    // lvl represents the level of approximation
    const double delta = double(tmax) / datapoints;
	const double h_step = 1e-2;
	const int number_steps = delta/h_step;
    double t, ti = 0.;
    fstream daten(name, ios::out);
    unsigned dim = lvl ? (3 * (lvl + 1)) : 6;

    // txts for coefficients
    string file_a = "alpha_";
    string file_b = "beta_";

    stringstream sdoub;
    sdoub << setprecision(2) << gamma;

    switch(coupling) {
        case 1:
            file_a = "equi/" + file_a;
            file_b = "equi/" + file_b;
            break;
        case 2:
            file_a = "gamma/" + file_a;
            file_b = "gamma/" + file_b;
            file_a += "gamma" + sdoub.str() + "_";
            file_b += "gamma" + sdoub.str() + "_";
            break;
    }

    // coeff of tridiagonal matrix
    parameter params;
    params.alpha      = new double[lvl];
    params.beta       = new double[lvl];
    params.B          = H;
    params.level      = lvl;
    params.B_function = B_function;
    parameter* p      = &params;

    // import exact lanczos coefficients alpha, beta from txt (maple)
    long double a = 0.;
    long double b = 0.;
    if(N > 1) {
        file_a += "N" + to_string(N) + "_";
        file_b += "N" + to_string(N) + "_";
    } else if(N == 1) {
        file_a += "Ninfty_";
        file_b += "Ninfty_";
    }
    file_a += "lvl=" + to_string(lvl) + ".txt";
    file_b += "lvl=" + to_string(lvl) + ".txt";

    cout << file_a << endl;

    //	file_a="gamma/alpha_gamma0.01_N1000_lvl=32.txt";
    //	file_b="gamma/beta_gamma0.01_N1000_lvl=32.txt";

    fstream coeffa(file_a, ios::in);
    fstream coeffb(file_b, ios::in);
    for(unsigned i = 0; i < lvl; ++i) {
        coeffa >> a;
        coeffb >> b;
        params.alpha[i] = a;
        params.beta[i]  = b;
    }

    /*cout << "lanczos coefficients:" << endl;
    cout << "alpha:" << endl;
    out(params.alpha,lvl);
    cout << endl;
    cout << "beta:" << endl;
    out(params.beta,lvl);
    cout << endl;
    */

    // P-vectors with P[0,1,2]=CS
    double* P = new double[dim];

    // starting value of S0_z
    double S0[3];

    // steps before pulse
    int pulsestep = int(t_pulse / delta);
    // ensemble mean
    double* sx_mean   = new double[datapoints];
    double* sxsx_mean = new double[datapoints];
    double sx_oldmean;
    double* sx_err    = new double[datapoints];
    double* sy_mean   = new double[datapoints];
    double* sysy_mean = new double[datapoints];
    double sy_oldmean;
    double* sy_err    = new double[datapoints];
    double* sz_mean   = new double[datapoints];
    double* szsz_mean = new double[datapoints];
    double sz_oldmean;
    double* sz_err = new double[datapoints];
    double* sysx0  = new double[datapoints];

    double sxsx = 0.;
    double sysy = 0.;
    double szsz = 0.;

    for(unsigned i = 0; i < datapoints; i++) {
        sx_mean[i]   = 0.;
        sxsx_mean[i] = 0.;
        sx_err[i]    = 0.;
        sy_mean[i]   = 0.;
        sysy_mean[i] = 0.;
        sy_err[i]    = 0.;
        sz_mean[i]   = 0.;
        szsz_mean[i] = 0.;
        sz_err[i]    = 0.;
        sysx0[i]     = 0.;
    }

    // binning of component of overhauser
    //	int noofpulses=int(tmax/t_pulse); //number of pulses in simulation
    //	double bin_max=2.;
    //	double bin_width=0.01;
    //	int noofbins=int(2*bin_max/bin_width+1);
    //	int *bin=new int[noofbins*(noofpulses-1)]; //2d-array auf 1d gemapt
    int pulse_no=0; //number of current pulse
    //	fstream binning_daten("binning_"+name,ios::out);
    //	int bin_index =0;

    // initialize odes
    gsl_odeiv2_system ode = {rightside, NULL, dim, p};

    // rk4 = classic 4th order RK, rkf45 = explicit embedded Runge-Kutta-Fehlberg (4, 5) method.
    gsl_odeiv2_driver* d =
        gsl_odeiv2_driver_alloc_y_new(&ode, gsl_odeiv2_step_rkf45, 1e-6, 1e-6, 0);

    int status;

    //	gsl_odeiv2_driver_set_hmin(d,1e-6);
    //	gsl_odeiv2_driver_set_hmax(d,delta);
    //	gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&ode,gsl_odeiv2_step_rk4,0.002,1e-6,0.)

    double *puls0 = new double[max_configs];
    double *puls50 = new double[max_configs];
    double *puls200 = new double[max_configs];
    double *puls1000 = new double[max_configs];

    // ensemble average
    for(int config = 1; config <= max_configs; config++) {
        t  = 0.;
        ti = 0.;
        if(config % (max_configs / 10) == 0) {
            cout << endl
                 << "Fortschritt: " << name << " bei " << int(config * 10 / max_configs) << "/10";
        }


        // starting values for P-fields

        // random
        for(unsigned i = 0; i < dim; i++) {
            P[i] = normal(gen);
        }
        // turn spin to x-dir at t=0
        flip_spin(P, 0);

        // starting values for correlation
        S0[0] = P[0];
        S0[1] = P[1];
        S0[2] = P[2];
        
        // DGLs with RK4
        for(unsigned step = 0; step < datapoints; step++) {
            ti = (step + 1) * delta;

            if(step == 0){
                puls0[config - 1] = P[5];
            }
            
            if(pulsestep) {
                if(step % pulsestep == 0 && step != 0) {

                    pulse_no=int(step/pulsestep)-1; //first pulse at t=t_pulse
                    //					//vor jedem Puls z-komponente von pseudo-Overhauser (=P[5])
                    switch(pulse_no){
                        case 49:
                            puls50[config - 1] = P[5];
                            break;
                        case 199:
                            puls200[config - 1] = P[5];
                            break;
                        case 999:
                            puls1000[config - 1] = P[5];
                    }
                    //binnen

                    //					bin_index
                    //=int((noofbins-1)/2+pulse_no*noofbins+P[5]/bin_width);
                    //					if(bin_index>=0&&bin_index<noofbins*(noofpulses-1)){
                    //						bin[bin_index]++;
                    //					}

                    // pulse S0:
                    flip_spin(P, 0); // Pulse 2
                    //					pi_pulse(P);
                    //					abs_pulse(P,0);	//Pulse 1

                    //					around_z_pulse(P);
                }
            }
            
            sxsx = P[0] * S0[0];
            sysy = P[1] * S0[1];
            szsz = P[2] * S0[2];

            // S_0 mitteln ohne Fehler
            sx_mean[step] += P[0]; // /max_configs;
            sy_mean[step] += P[1]; // /max_configs;
            sz_mean[step] += P[2]; // /max_configs;

            sxsx_mean[step] += sxsx; // /max_configs;
            sysy_mean[step] += sysy; // /max_configs;
            szsz_mean[step] += szsz; // /max_configs;

            sysx0[step] += P[1] * S0[0]; // /max_configs;


            // S_0 mitteln mit Fehler
            //			sx_oldmean=sxsx_mean[step];
            //			sxsx_mean[step]+=(sxsx-sxsx_mean[step])/config;
            //			sx_err[step]+=(sxsx-sxsx_mean[step])*(sxsx-sx_oldmean);
            //			sy_oldmean=sysy_mean[step];
            //			sysy_mean[step]+=(sysy-sysy_mean[step])/config;
            //			sy_err[step]+=(sysy-sysy_mean[step])*(sysy-sy_oldmean);
            //			sz_oldmean=szsz_mean[step];
            //			szsz_mean[step]+=(szsz-szsz_mean[step])/config;
            //			sz_err[step]+=(szsz-szsz_mean[step])*(szsz-sz_oldmean);

            // time evolution:
            status = gsl_odeiv2_driver_apply(d, &t, ti, P);
			//status = gsl_odeiv2_driver_apply_fixed_step(d, &t, h_step, number_steps, P);

			
            if(status != GSL_SUCCESS) {
                cerr << "GSL Error: return value = " << status << endl;
                exit(1);
            }
        }
    }

    fstream PFieldDist;
    PFieldDist.open("./data/overhauserDist.dat", fstream::out);
    PFieldDist << "#p0\tp50\tp200\tp1000\n";
    for(int i = 0; i<max_configs; i++){
        PFieldDist << puls0[i] << '\t' << puls50[i] << '\t' << puls200[i] << '\t' << puls1000[i] << '\n';
    }
    PFieldDist.close();
    
    delete[] puls0;
    delete[] puls50;
    delete[] puls200;
    delete[] puls1000;

    // Mittelwert bilden:
    for(unsigned i = 0; i < datapoints; i++) {
        sx_mean[i] = sx_mean[i] / max_configs;
        sy_mean[i] = sy_mean[i] / max_configs;
        sz_mean[i] = sz_mean[i] / max_configs;

        sxsx_mean[i] = sxsx_mean[i] / max_configs;
        sysy_mean[i] = sysy_mean[i] / max_configs;
        szsz_mean[i] = szsz_mean[i] / max_configs;

        sysx0[i] = sysx0[i] / max_configs;
    }


    gsl_odeiv2_driver_free(d);
    
    /*
    //save data:
    daten << "#t\tSxSx\tSySy\tSzSz\tSx\tSy\tSz\tSySx0\n";
    for(int step = 0; step < datapoints; step++) {
        daten << step * delta << "\t" << sxsx_mean[step] << "\t";
        daten << sysy_mean[step] << "\t";
        daten << szsz_mean[step];
        daten << "\t" << sx_mean[step] << "\t";
        daten << sy_mean[step] << "\t";
        daten << sz_mean[step] << "\t";
        daten << sysx0[step];
        //
        ////			//Fehler:
        ////			daten << "\t" << sqrt(sx_err[step]/(max_configs*(max_configs-1))) << "\t";
        ////			daten << sqrt(sy_err[step]/(max_configs*(max_configs-1))) << "\t";
        ////			daten << sqrt(sz_err[step]/(max_configs*(max_configs-1)));
        ////
        daten << endl;
    }
    */

    // save binning data
    //	for(unsigned i=0;i<noofbins*(noofpulses-1);i++){
    //		if(i%noofbins==0&&i>0) binning_daten << endl;
    //		binning_daten << bin[i] << "\t";
    //	}

    delete[] P;
    delete[] params.alpha;
    delete[] params.beta;
    delete[] sx_mean;
    delete[] sxsx_mean;
    delete[] sx_err;
    delete[] sy_mean;
    delete[] sysy_mean;
    delete[] sy_err;
    delete[] sz_mean;
    delete[] szsz_mean;
    delete[] sz_err;
    delete[] sysx0;
    //	delete[] bin;
}


// int main(int argN, char** args){
int main(int argC, char* argV[])
{

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()("help", "produce help message")("N", po::value<int>()->default_value(1),
                                                       "set amount of bath spins, 1 := infty")(
        "ensembles", po::value<int>()->default_value(1000),
        "set ensemble size")("lvl", po::value<int>()->default_value(8), "set lanczos level")(
        "tmax", po::value<int>()->default_value(100), "set maximal time")(
        "datapoints", po::value<int>()->default_value(1000),
        "set amount of data points")("coupling", po::value<int>()->default_value(2),
                                     "set coupling distribution, 1 := equi, 2 := exp, 3 := posneg")(
        "gamma", po::value<double>()->default_value(0.01), "set gamma value")(
        "H", po::value<double>()->default_value(1.), "set B-field strength in z-direction")(
        "Hfunction", po::value<int>()->default_value(1), "set function of H-field, 1 := H = const")(
        "tpulse", po::value<double>()->default_value(10.), "set pulse period")(
        "add", po::value<string>()->default_value(""), "add data-file text extension");

    po::variables_map vm;
    po::store(po::parse_command_line(argC, argV, desc), vm);
    po::notify(vm);

    if(vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }

    int N          = vm["N"].as<int>();
    int ensembles  = vm["ensembles"].as<int>();
    int lvl        = vm["lvl"].as<int>();
    int tmax       = vm["tmax"].as<int>();
    int datapoints = vm["datapoints"].as<int>();
    int coupling   = vm["coupling"].as<int>();
    double gamma   = vm["gamma"].as<double>();
    double H       = vm["H"].as<double>();
    int B_function = vm["Hfunction"].as<int>();
    double t_pulse = vm["tpulse"].as<double>();
    string zusatz  = vm["add"].as<string>();

    stringstream sdoub;
    sdoub << setprecision(2) << gamma;
    string dateiname = "data/lanczos_";

    switch(coupling) {
        case 1:
            dateiname += "equi_";
            break;
        case 2:
            dateiname += "gamma" + sdoub.str() + "_";
            break;
        case 3:
            dateiname += "posneg";
            break;
    }

    dateiname += "lvl=";
    dateiname += to_string(lvl);
    dateiname += "_";
    if(N == 1)
        dateiname += "INFTY";
    else
        dateiname += to_string(N);
    dateiname += "spins_";
    dateiname += zusatz;
    dateiname += ".txt";
    cout << endl << dateiname << endl;
    cout << "Hierarchy with Lanczos started with: " << endl;
    cout << ensembles << " Ensembles" << endl;
    cout << datapoints << " Datapoints" << endl;
    cout << "Gamma=" << gamma << endl;
    cout << "level=" << lvl << endl;
    cout << "Daten schreiben in: " << dateiname << endl;

    hierarchy(N, ensembles, lvl, tmax, coupling, gamma, dateiname, datapoints, H, B_function,
              t_pulse);
    cout << endl << endl;


    return 0;
}
