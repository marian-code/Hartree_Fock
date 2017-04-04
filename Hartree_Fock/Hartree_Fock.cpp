// Hartree_Fock.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <fstream>
#include <Eigen\Core>
#include <Eigen\Dense>
#include <Eigen\Eigenvalues>

using namespace Eigen;
using namespace std;

double Integral(VectorXd Wave_function, const double h_step, int i, const int matrix_size, double density_norm)
{
	int j, k;
	double density;
	double integral = 0;

	if (i == 0) i++;

	for (j = i; j < matrix_size; j++)
	{
		density = 0;
		//spocitat hustotu pravdepodobnosti
		for (k = 1; k < j; k++)
			density += h_step*pow(Wave_function(k), 2);
		
		//integral += (h_step*density)*(1.0/pow(double(j*h_step), 2)) / density_norm;
		integral += (density) / (j*j*h_step);
	}
	integral /= density_norm;

	return integral;
}

double Density_norm(VectorXd Wave_function, const double h_step, const int matrix_size)
{
	double density_norm; 
	int k;
	
	//spocitat normovanie hustoty pravdepodobnosti
	density_norm = 0;
	for (k = 1; k < matrix_size; k++)
		density_norm += h_step*pow(Wave_function(k), 2);
	return density_norm;
}

MatrixXd Hamiltonian_compute(VectorXd Wave_function, const double h_step, const int matrix_size, double density_norm)
{
	int j, k;
	MatrixXd Hamiltonian(matrix_size, matrix_size);

	auto time_start = chrono::high_resolution_clock::now();
	cout << "computing new hamiltonian... " << endl;

	for (j = 0; j < matrix_size; j++)
	{
		for (k = 0; k < matrix_size; k++)
		{
			if (j == k) Hamiltonian(j, k) = 1.0 / (h_step*h_step) - 2.0 / (h_step*(j+1)) + Integral(Wave_function, h_step, j, matrix_size, density_norm);
			else if (j == k - 1) Hamiltonian(j, k) = -1.0 / (2*h_step*h_step);
			else if (j == k + 1) Hamiltonian(j, k) = -1.0 / (2*h_step*h_step);
			else Hamiltonian(j, k) = 0;		
		}
	}

	auto time_end = chrono::high_resolution_clock::now();
	cout << "Hamiltoniam computed in: " << chrono::duration_cast<chrono::nanoseconds>(time_end - time_start).count()*10E-10 << "s" << endl << endl;

	return Hamiltonian;
}

double Energy(VectorXd Wave_function, double energy, const int matrix_size, const double h_step, double density_norm)
{
	int j, k;
	double density;
	double total_energy = 0;

	for (j = 1; j < matrix_size; j++)
	{
		density = 0;
		//spocitat hustotu pravdepodobnosti
		for (k = 1; k < j; k++)
			density += h_step*pow(Wave_function(k), 2);

		//total_energy += h_step*(1.0 / double(j*h_step))*pow(Wave_function(j), 2)*density;
		total_energy += pow(Wave_function(j), 2)*density / double(j);
	}

	total_energy = 2 * energy - 2 * total_energy / pow(density_norm, 2);

	return total_energy;
}

int main()
{
	const double pi = 3.1415926535897932384626433;
	const double cut_off = 14; //14 je rozumny cutoff
	const int matrix_size = 1001;
	const double h_step = cut_off / double(301);
	const double delta_epsilon = 10E-8;

	MatrixXd Hamiltonian; //pouzivaju sa hodnoty typu double k�li r�chlosti
	VectorXd Wave_function;
	SelfAdjointEigenSolver<MatrixXd> es; //selfadjoint je pre hermitovske matice
	
	double epsilon_new, epsilon_old;
	double norm;
	double total_energy;
	int i;

#pragma region inicializacia

	//nastavenie presnosti vypisu
	cout.precision(7);
	cout.setf(std::ios::fixed, std::ios::floatfield);

	cout << "delta epsilon: " << delta_epsilon << endl;
	cout << "cut off: " << cut_off << endl;
	cout << "space step: " << h_step << endl << endl;


	//inicializacia casovych premennych
	auto time_start = chrono::high_resolution_clock::now();
	auto time_end = chrono::high_resolution_clock::now();

	auto time_total_start = chrono::high_resolution_clock::now();
	auto time_total_end = chrono::high_resolution_clock::now();

	auto time_cycle_start = chrono::high_resolution_clock::now();

	//wavefunction vector resize
	Wave_function.resize(matrix_size);

	//inicializacia epsilonu z predoslej iteracie 
	epsilon_old = 10E10;
	epsilon_new = 10;

#pragma endregion

	//initial wavefunction guess
	for (i = 0; i < matrix_size; i++)
		//Wave_function(i) = h_step*i*exp((-1)*pow((i*h_step -0.98) / (double(matrix_size) / 1000.0), 2));
		Wave_function(i) = sin(pi*i / (matrix_size - 1));

	//cout << "initial wavefunction is: " << endl << Wave_function << endl << endl;

	//cyklus
	i = 0;
	while (abs(epsilon_new - epsilon_old) > delta_epsilon)
	{
		cout << "start of iteration number: " << i + 1 << "!!!!!!!!!!!!!!!!!!!!!" << endl << endl;

		epsilon_old = epsilon_new;

		time_cycle_start = chrono::high_resolution_clock::now();

		//spocitat normalizaciu vlnovej funkcie pre zrychlenie vypoctu
		norm = Density_norm(Wave_function, h_step, matrix_size);

		//naplnit Hamiltonin hodnotami
		Hamiltonian = Hamiltonian_compute(Wave_function, h_step, matrix_size, norm);
		//cout << "hamiltonian: " << endl << Hamiltonian << endl;
		/*for (j = 0; j < matrix_size; j++)
		{
			if (j - 1 >= 0 && j + 1 <= matrix_size - 1)
				cout << Hamiltonian(j, j - 1) << "\t" << Hamiltonian(j, j) << "\t" << Hamiltonian(j, j + 1) << endl;
			else if (j - 1 < 0)
				cout << "\t" << "\t" << Hamiltonian(j, j) << "\t" << Hamiltonian(j, j + 1) << endl;
			else if (j + 1 > matrix_size - 1)
				cout << Hamiltonian(j, j - 1) << "\t" << Hamiltonian(j, j) << endl;
		}*/

		time_start = chrono::high_resolution_clock::now();
		cout << "diagonalizing hamiltonian... " << endl;

		//diagonalizacia hamiltonianu
		es.compute(Hamiltonian);


		time_end = chrono::high_resolution_clock::now();
		cout << "hamiltonian diagonalized sucesfully in: " << chrono::duration_cast<chrono::nanoseconds>(time_end - time_start).count()*10E-10 << "s" << endl << endl;

		//vlastne hodnoty by mali byt usporiadane od najmensej po najv��iu
		epsilon_new = es.eigenvalues()[0];
		cout << "Lowest eigenvalue is:" << es.eigenvalues()[0] << endl;
		cout << "difference of eigenvalues: " << abs(epsilon_new - epsilon_old) << endl << endl;

		//cout << "eigenvalues are: " << endl << es.eigenvalues() << endl << endl;
		//cout << "Eigenvector with lowest eigenvalule is:" << endl << es.eigenvectors().col(iterator) << endl << endl;

		//priradenie novej vlnovej funkcie
		Wave_function = es.eigenvectors().col(0);

		time_total_end = chrono::high_resolution_clock::now();
		cout << "actuall iteration time: " << chrono::duration_cast<chrono::nanoseconds>(time_total_end - time_cycle_start).count()*10E-10 << "s" << endl;
		cout << "total CPU time till now: " << chrono::duration_cast<chrono::nanoseconds>(time_total_end - time_total_start).count()*10E-10 << "s" << endl;
		cout << "end of iteration number: " << i + 1 << "!!!!!!!!!!!!!!!!!!!" << endl << endl;
		i++;
	}

	//celkova energia
	total_energy = Energy(es.eigenvectors().col(0), es.eigenvalues()[0], matrix_size, h_step, norm);
	cout << "total energy is: " << total_energy << endl << endl;

	ofstream outfile1;
	outfile1.open("function.dat");
	for (i = 1; i < matrix_size; i++)
	{
		outfile1 << i*h_step << "\t" << Wave_function(i) << endl;
	}
	outfile1.close();

	system("PAUSE");
}

