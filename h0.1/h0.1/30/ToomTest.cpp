// Toom's Rule
#include <fstream>
#include <iostream>
#include <random>
#include <math.h>
#include <complex>

using namespace std;

typedef complex<double> dcomp;

#define N 30
#define TempPoints 40
#define MAXWIN 20
double pi = 2*asin(1);

//random seed generators
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> dis (0.0, 1.0);
//double kT[TempPoints] = {0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.10, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135, 0.14, 0.145, 0.15, 0.155, 0.16, 0.165, 0.17, 0.175, 0.18, 0.185, 0.19, 0.195, 0.2};
//double kT[TempPoints] = {0.002, 0.004, 0.006, 0.008, 0.01, 0.012, 0.014, 0.016, 0.018, 0.020, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026, 0.027, 0.028, 0.029, 0.03, 0.031, 0.032, 0.033, 0.034, 0.035, 0.036, 0.037, 0.038, 0.039, 0.04};
//double kT[TempPoints] = {0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8};
double kT[TempPoints] = {0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14, 0.16, 0.18, 0.2, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.34, 0.36, 0.38, 0.4, 0.42, 0.44, 0.46, 0.48, 0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8};
short neighbors[2];
double mag[3];	//store moments of magnetization- <mag>, <mag^2>, <mag^4>
double chiCorel;
double WinMag[3];	//store window average of the moments of magnetization
double WinChiCorel;
double energy1;
double WinEnergy1;
double energy2;
double WinEnergy2;

void create (short* arr)
{
	for (int i=0; i<N*N; i++)
	{	
		arr[i]=1;
		// if (dis(gen)<0.5)
		// 	arr[i]=-1;
		// else
		// 	arr[i]=+1;
	}
}

short* findAdjacent (int x, short* arr)
{
	int row = x/N, col = x % N;	
	//neighbors of lattice[x]
	int r= N*row+ col+1;
	int u= N*(row-1)+ col;
	
	// set boundaries
	if (row == 0)
		u = N*(N-1)+col;
	if (col == N-1)
		r = N*(row);

	neighbors[0] = arr[u];
	neighbors[1] = arr[r];
	return neighbors;
} 

void update (short* lattice, short* lattice2, double T)
{	
	double p, q;
	p = 11*T/20;
	q = 9*T/20;	//h=0.1

	for (int x = 0; x < N*N; x++)
	{
		findAdjacent (x, lattice);
		if (neighbors[0]==neighbors[1])
			lattice2[x]=neighbors[0];
		else
			lattice2[x]=lattice[x];

		if (lattice2[x]==1)
		{
			if (dis(gen)<=q)
				lattice2[x]=-lattice2[x];	//If spin up, flip with prob q
		}
		else if (lattice2[x]==-1)
		{
			if (dis(gen)<=p)
				lattice2[x]=-lattice2[x];	//If spin down, flip with prob p
		}
	}
}

double* magnetization (short* arr)
{
	double sum=0;
	for (int i=0; i<N*N; i++)
		sum += arr[i];		
	mag[0] = abs(sum/(N*N));	// magnetization of system
	mag[1] = pow (mag[0], 2); //square of magnetization of system
	mag[2] = pow (mag[0], 4); // fourth power of magnetization of system
	return mag;
}

double Chi (short* arr)
{
	dcomp i, chi = 0, b;
	double exponent;
	i = -1;
	i = sqrt(i);
	for (int iter=0; iter<N*N; iter++)
	{
		exponent = 2*(pi/N)*((iter)%N);
		//exponent = 0;	// k=0 gives X(0) = N*N*<m^2>
		b = exp(exponent*i);
		b *= arr[iter];
		b /= N*N;
		chi +=  b;
	}
	return N*N*pow(abs(chi),2);
	return 0;
}

double* ToomEnergy (short* arr)
{
	double sum1 = 0;
	double sum2 = 0;

	for (int iter=0; iter<N*N; iter++)
	{
		findAdjacent(iter, arr);
		if (neighbors[0]==neighbors[1] and neighbors[0]==arr[iter])
		{
			sum1 += -1;
			sum2 += -1;
		}
		else if (neighbors[0]==neighbors[1] and neighbors[0]==-arr[iter])
		{
			sum2 += 1;
		}
	}

	energy1 = sum1/(N*N);
	energy2 = sum2/(N*N);
}

int main ()
{
	ofstream opfile, opfile2;
  	opfile.open ("Toom30h0_1.ods");
  	opfile2.open ("Toom30h0_1evolve.ods");

  	opfile << "N=30" << endl;
  	opfile2 << "N=30" << endl;
	
	int iter=0, T=0, window=0;
	short* temp;
	short* lattice = new short[N*N];
	short* lattice2 = new short[N*N];

	for (T=0; T<TempPoints; T++) 
	{
		//int stepCount = 0;
	 	create(lattice);
		opfile2 << "T =" << kT[T] << endl;
		for (int i=0; i<N; i++)
		{
			for (int j=0; j<N; j++)
				opfile2 << lattice[N*i + j] << "\t";
			opfile2 << endl;
		}	// print initial lattice

		opfile2 << endl;
		opfile << "T \t window \t <mag> \t <mag^2> \t <mag^4> \t <chi(k)> \t correl length \t Binder \t Energy1 \t Energy2" << endl;
		for (window=0; window < MAXWIN ; window++)
		{
			int windowsize = pow(2,window);
			for (int k=0;k<3;k++)
				WinMag[k]=0;
			WinEnergy1=0;
			WinEnergy2=0;		
			WinChiCorel = 0;	// Initialize window averages to 0
			for (int i=0; i<windowsize; i++)
			{
				update(lattice, lattice2, kT[T]);	//1 time step

				temp = lattice;
				lattice = lattice2;
				lattice2 = temp;	// swap lattice with lattice 2 to work on updated lattice

				//stepCount++;
				// if (stepCount%10000 == 0)
				// {
				// 	for (int i=0; i<N; i++)
				// 	{
				// 		for (int j=0; j<N; j++)
				// 			opfile2 << lattice[N*i + j] << "\t";
				// 		opfile2 << endl;
				// 	}	// print lattice after every 10,000 timesteps
				// 	opfile2 << endl;
				// }
				magnetization(lattice);
				chiCorel = Chi(lattice);
				ToomEnergy(lattice);
				for (int k=0; k<3; k++)
					WinMag[k] += mag[k]/windowsize;	//storing Window average of magnetization and energy
				WinEnergy1 += energy1/windowsize;
				WinEnergy2 += energy2/windowsize;	
				WinChiCorel += chiCorel/windowsize; //storing Window average of chi

			}

			opfile2 << window << endl;

			for (int i=0; i<N; i++)
			{
				for (int j=0; j<N; j++)
					opfile2 << lattice[N*i + j] << "\t";
				opfile2 << endl;
			}	// also print lattice after every window
			opfile2 << endl;

			opfile << kT[T] << "\t" << window << "\t" << WinMag[0] << "\t" << WinMag[1] << "\t" << WinMag[2] << "\t" << WinChiCorel <<  "\t" << 0.5*sqrt(N*N*WinMag[1]/WinChiCorel-1)/sin(pi/N)  << "\t" << 1-WinMag[2]/(3*pow(WinMag[1],2)) << "\t" << WinEnergy1 << "\t" << WinEnergy2 << endl;	// Window avg of chi to be used for correlation length
		}
	}	
	
	delete(lattice);
	delete(lattice2);
	opfile.close();
	opfile2.close();		
	return 0;
}
