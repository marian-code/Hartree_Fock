// Hartree_Fock.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <Eigen\Core>
#include <Eigen\Dense>
#include <Eigen\Eigenvalues>

using namespace Eigen;
using namespace std;

float Integral(VectorXf Wave_function, const float h_step, int i, const int matrix_size ) //pozor na const premenne  - neoverene spravanie!!!!
{
	int j, k;
	float density, density_norm;
	float integral = 0;

	if (i == 0) i++;

	//spocitat normovanie hustoty pravdepodobnosti
	density_norm = 0;
	for (k = 1; k < matrix_size; k++)
		density_norm += h_step*pow(Wave_function(k), 2);

	for (j = i; j < matrix_size; j++)
	{
		density = 0;
		//spocitat hustotu pravdepodobnosti
		for (k = 1; k < j; k++)
			density += h_step*pow(Wave_function(k), 2);
		
		integral += (h_step*density)*(1.0/pow(float(j*h_step), 2)) / density_norm;
	}
	return integral;
}

MatrixXf Hamiltonian_compute(VectorXf Wave_function, const float h_step, const int matrix_size)
{
	int j, k;
	MatrixXf Hamiltonian(matrix_size, matrix_size);

	auto time_start = chrono::high_resolution_clock::now();
	cout << "computing new hamiltonian... " << endl;

	for (j = 0; j < matrix_size; j++)
	{
		for (k = 0; k < matrix_size; k++)
		{
			if (j == k)
			{
				if (j == 0) Hamiltonian(j, k) = 2.0 / (h_step*h_step) - 4.0 / (h_step*j + h_step) + 2 * Integral(Wave_function, h_step, j, matrix_size); //toto nie je tak ako by malo byt indexy sú poposúvané
				else Hamiltonian(j, k) = 2.0 / (h_step*h_step) - 4.0 / (h_step*j) + 2 * Integral(Wave_function, h_step, j, matrix_size);
			}
			else if (j == k - 1) Hamiltonian(j, k) = -1.0 / (h_step*h_step);
			else if (j == k + 1) Hamiltonian(j, k) = -1.0 / (h_step*h_step);
			else Hamiltonian(j, k) = 0;		
		}
	}

	auto time_end = chrono::high_resolution_clock::now();
	cout << "Hamiltoniam computed in: " << chrono::duration_cast<chrono::nanoseconds>(time_end - time_start).count()*10E-10 << "s" << endl << endl;

	return Hamiltonian;
}

float Energy(VectorXf Wave_function, float energy, const int matrix_size, const float h_step)
{
	int j, k;
	float density, density_norm;
	float total_energy = 0;

	//spocitat normovanie hustoty pravdepodobnosti
	density_norm = 0;
	for (k = 1; k < matrix_size; k++)
		density_norm += h_step*pow(Wave_function(k), 2);

	for (j = 1; j < matrix_size; j++)
	{
		density = 0;
		//spocitat hustotu pravdepodobnosti
		for (k = 1; k < j; k++)
			density += h_step*pow(Wave_function(k), 2);

		total_energy += h_step*(1.0 / double(j*h_step))*pow(Wave_function(j), 2)*density;
	}

	total_energy = 2 * energy - 2 * total_energy / pow(density_norm, 2);

	return total_energy;
}

int main()
{
	const double pi = 3.1415926535897932384626433;
	const int matrix_size = 301;
	MatrixXf Hamiltonian; //pouzivaju sa hodnoty typu float kôli rýchlosti
	VectorXf Wave_function;
	SelfAdjointEigenSolver<MatrixXf> es; //selfadjoint je pre hermitovske matice
	const float h_step = 0.01; //elektron by mal byt niekde v 0.58 borh cca
	const float epsilon = 10E-8;
	float energy_minimum;
	float total_energy;
	int i, j, k;
	int iterator;

#pragma region inicializacia

	//nastavenie presnosti vypisu
	cout.precision(7);
	cout.setf(std::ios::fixed, std::ios::floatfield);

	cout << "epsilon " << epsilon << endl;
	cout << "space step: " << h_step << endl;

	//inicializacia casovych premennych
	auto time_start = chrono::high_resolution_clock::now();
	auto time_end = chrono::high_resolution_clock::now();

	auto time_total_start = chrono::high_resolution_clock::now();
	auto time_total_end = chrono::high_resolution_clock::now();

	auto time_cycle_start = chrono::high_resolution_clock::now();

	//wavefunction vector resize
	Wave_function.resize(matrix_size);

#pragma endregion

	//initial wavefunction guess
	for (i = 0; i < matrix_size; i++)
		//Wave_function(i) = h_step*i*exp((-1)*pow((i*h_step -0.98) / (double(matrix_size) / 1000.0), 2));
		Wave_function(i) = sin(pi*i / (matrix_size - 1));

	cout << "initial wavefunction is: " << endl << Wave_function << endl << endl;

	//cyklus
	i = 0;
	while (i < 13)
	{
		time_cycle_start = chrono::high_resolution_clock::now();

		//naplnit Hamiltonin hodnotami
		Hamiltonian = Hamiltonian_compute(Wave_function, h_step, matrix_size);
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
		cout << "hamiltonian diagonalized sucesfully in: " << chrono::duration_cast<chrono::nanoseconds>(time_end - time_start).count()*10E-10 << "s" << endl;

		
		//vlastne hodnoty by mali byt usporiadane od najmensej po najväèšiu
		iterator = 0;
		/*energy_minimum = 10E10;
		for (j = 0; j < matrix_size; j++)
		{
			if (es.eigenvalues()[j] < energy_minimum)
			{
				iterator = j;
				energy_minimum = es.eigenvalues()[j];
			}
		}*/

		//cout << "Eigenvector with lowest eigenvalule is:" << endl << es.eigenvectors().col(iterator) << endl << endl;
		cout << "Lowest eigenvalue is:" << endl << es.eigenvalues()[iterator] << endl << endl;
		//cout << "eigenvalues are: " << endl << es.eigenvalues() << endl << endl;

		//celkova energia
		total_energy = Energy(es.eigenvectors().col(iterator), es.eigenvalues()[iterator], matrix_size, h_step);

		cout << "total energy is: " << total_energy << endl << endl;

		//priradenie novej vlnovej funkcie
		Wave_function = es.eigenvectors().col(iterator);

		time_total_end = chrono::high_resolution_clock::now();
		cout << "actuall iteration time: " << chrono::duration_cast<chrono::nanoseconds>(time_total_end - time_cycle_start).count()*10E-10 << "s" << endl;
		cout << "total CPU time till now: " << chrono::duration_cast<chrono::nanoseconds>(time_total_end - time_total_start).count()*10E-10 << "s" << endl << endl;
		cout << "end of iteration number: " << i + 1 << "!!!!!!!!!!" << endl;
		i++;
	}

	system("PAUSE");
}

