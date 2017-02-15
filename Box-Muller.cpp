// Includes libraries 
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;

//Méthode de Monte-Carlo pour le calcul de surface et de volume en plusieurs dimensions
//Méthode de Box-Muller pour le calcul d'aire à échantillonage non uniforme avec des gaussiennes

int box(void) {

	double M_PI = 3.14159265358979323846;
	int N = 100;
	double S = 0;
	double sigma = 0.2;

	for (int i = 0; i < N / 2; i++) {
		double u = rand() / (double)RAND_MAX;
		double v = rand() / (double)RAND_MAX;
		double x = sqrt(-2 * log(u))*cos(2 * M_PI*v);
		double y = sqrt(-2 * log(u))*sin(2 * M_PI*v);

		if (x > -M_PI / 2 && x < M_PI / 2) { // Si x est compris dans l'interval ou f(x) est non nul
			S += pow(cos(x), 20) / exp(-x*x / (2 * sigma * sigma)) / (sigma*sqrt(2 * M_PI));
		}
		if (y > -M_PI / 2 && y < M_PI / 2) { // si y est compris dans l'interval ou f(y) est non nul
			S += pow(cos(y), 20) / exp(-y*y / (2 * sigma * sigma)) / (sigma*sqrt(2 * M_PI));
		}
	}

	return S / N;
}