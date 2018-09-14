// Toom's Rule
#include <fstream>
#include <iostream>
#include <random>
#include <math.h>
#include <complex>

using namespace std;

typedef complex<double> dcomp;

#define N 50
#define TempPoints 20
#define MAXWIN 20
double pi = 2*asin(1);

//random seed generators
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<> dis (0.0, 1.0);
//double kT [TempPoints] = {0.1, 0.5, 1.0, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.22, 2.24, 2.26, 2.27, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.5, 4.0, 4.5, 4.9};
double kT[TempPoints] = {0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26, 0.28, 0.30, 0.32, 0.34, 0.36, 0.38, 0.40};
//double kT[TempPoints] = {1.02, 1.04, 1.06, 1.08, 1.10, 1.12, 1.14, 1.16, 1.18, 1.20, 1.22, 1.24, 1.26, 1.28, 1.30, 1.32, 1.34, 1.36, 1.38, 1.40, 1.42, 1.44, 1.46, 1.48, 1.50};
short neighbors[2];
double mag[3];	//store moments of magnetization- <mag>, <mag^2>, <mag^4>
double chiCorel;
double WinMag[3];	//store window average of the moments of magnetization
double WinChiCorel;

void create (short* arr)
{
	for (int i=0; i<N*N; i++)
		if (dis(gen)<0.5)
			arr[i]=-1;
		else
			arr[i]=+1;
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
	p = T/2;
	q = T/2;

	for (int x = 0; x < N*N; x++)
	{
		findAdjacent (x, lattice);
		if (neighbors[0]==neighbors[1])
			lattice2[x]=neighbors[0];
		else
			lattice2[x]=lattice[x];

		if (lattice2[x]==1)
			if (dis(gen)<=q)
				lattice2[x]=-lattice2[x];	//If spin up, flip with prob q
		if (lattice2[x]==-1)
			if (dis(gen)<=p)
				lattice2[x]=-lattice2[x];	//If spin down, flip with prob p
	}
}

double* magnetization (short* arr)
{
	double sum=0;
	for (int i=0; i<N*N; i++)
		sum += arr[i];		
	mag[0] = sum/(N*N);	// magnetization of system
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

int main ()
{
	ofstream opfile, opfile2;
  	opfile.open ("Toom50h0.ods");
  	//opfile2.open ("Toom10evolve.ods");

  	opfile << "N=50" << endl;
  	//opfile2 << "N=10" << endl;
	
	int iter=0, T=0, window=0;
	short* temp;
	short* lattice = new short[N*N];
	short* lattice2 = new short[N*N];

	for (T=0; T<TempPoints; T++) 
	{
	 	create(lattice);
		//opfile2 << "T =" << kT[T] << endl;
		// for (int i=0; i<N; i++)
		// {
		// 	for (int j=0; j<N; j++)
		// 		opfile2 << lattice[N*i + j] << "\t";
		// 	opfile2 << endl;
		// }	// print initial lattice

		//opfile2 << endl;
		opfile << "T \t window \t <mag> \t |<mag>| \t <mag^2> \t <mag^4> \t <chi(k)> \t correl length \t Binder" << endl;
		for (window=0; window < MAXWIN ; window++)
		{
			int windowsize = pow(2,window);
			for (int k=0;k<3;k++)
				WinMag[k]=0;	
			WinChiCorel = 0;	// Initialize window averages to 0
			for (int i=0; i<windowsize; i++)
			{
				//for (iter=0; iter < N*N; iter++)
				update(lattice, lattice2, kT[T]);	//1 time step
				magnetization(lattice);
				chiCorel = Chi(lattice);
				for (int k=0; k<3; k++)
					WinMag[k] += mag[k]/windowsize;	//storing Window average of magnetization and energy
				WinChiCorel += chiCorel/windowsize; //storing Window average of chi
				
				temp = lattice;
				lattice = lattice2;
				lattice2 = temp;	// swap lattice with lattice 2 to work on updated lattice

				// for (int i=0; i<N; i++)
				// {
				// 	for (int j=0; j<N; j++)
				// 		opfile2 << lattice[N*i + j] << "\t";
				// 	opfile2 << endl;
				// }	// print lattice after every timestep
				// opfile2 << endl;
			}
			
			opfile << kT[T] << "\t" << window << "\t" << WinMag[0] << "\t"  << abs(WinMag[0]) << "\t" << WinMag[1] << "\t" << WinMag[2] << "\t" << WinChiCorel <<  "\t" << 0.5*sqrt(N*N*WinMag[1]/WinChiCorel-1)/sin(pi/N)  << "\t" << 1-WinMag[2]/(3*pow(WinMag[1],2)) <<endl;	// Window avg of chi to be used for correlation length
		}
	}	
	
	delete(lattice);
	delete(lattice2);
	opfile.close();
	//opfile2.close();		
	return 0;
}