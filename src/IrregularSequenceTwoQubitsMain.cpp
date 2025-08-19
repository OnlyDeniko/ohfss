#include <iostream>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <vector>
#include <fstream>
#include <math.h>
#include <complex>
#include <mkl.h>
#include <omp.h>
#include <map>
#include "TwoQubitsKernel.h"
#include "ConstantsDescriptor.h"

using namespace std;

int main() {
	cout.precision(15);
	double time = omp_get_wtime();
	int N = 20; // кол-во уровней кубита
	int M = 3;
	double val = 2 * PI * 1e9;
	double tstep = 1e-14; // time grid step
	// main qubit frequencies
	double w1 = 5.4 * (2 * PI) * 1e9; // Частота внешнего управляющего поля
	double w2 = 5.0 * (2 * PI) * 1e9; // Частота внешнего управляющего поля
	// anharmonicities
	double mu1 = 0.35 * (2 * PI) * 1e9; // Параметр нелинейности первого кубита
	double mu2 = 0.25 * (2 * PI) * 1e9; // Параметр нелинейности первого кубита
	double g = 0.02 * (2 * PI) * 1e9; // параметр взаимодействия между кубитами
	
	// qubit capacities
	double Cq1 = 1e-12;
	double Cq2 = 1e-12;

	// connection capacities
	double Cc1 = 1e-15;
	double Cc2 = 1e-15;

	// pulse generation frequencies
	double wg1 = w2;
	double wg2 = w2;
	double tau = 4e-12; // Длительность импульса
	double phi = 0; // phase (number of grid steps paused on Q2)

	// wait time after pulse
	int waitq1 = 0;
	int waitq2 = 0;

	string init = "00"; // initial condition
	double Coeffs = 1;
	string operation = "PY"; // required operation (for fidelity calculation)
	int M1 = 0, M2 = 60;
	TwoQubitsConstantsDescriptor desc(N, M, M1, M2, val, tstep, w1, w2, mu1, mu2, g, Cq1, Cq2,
		Cc1, Cc2, wg1, wg2, tau, phi, waitq1, waitq2, init, Coeffs, operation, 3);

	string seqs(M2, '1');
	vector<int> seq;
	for (int i = 0; i < seqs.size(); i++) {
		if (seqs[i] == '-') {
			i++;
			seq.push_back(-1);
		}
		else if (seqs[i] == '0') {
			seq.push_back(0);
		}
		else {
			seq.push_back(1);
		}
	}
	assert(seq.size() == M1 + M2);
	double start = omp_get_wtime();
	TwoQubitsKernel kernel(desc);
	auto res = kernel.Fidelity(seq);
	start = omp_get_wtime() - start;
	cout << "TIME ELAPSED: " << start << '\n';
	cout << "\nF: " << res.fidelity << '\n';
	for (int i = 0; i < res.spec.Energies.size(); i++) {
		cout << res.spec.Energies[i] << '\t'
			<< res.spec.State_q1[i] << '\t'
			<< res.spec.State_q2[i] << '\t'
			<< res.Probabilities[i] << '\n';
	}
	return 0;
}



