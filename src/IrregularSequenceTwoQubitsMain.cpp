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
	double time = omp_get_wtime();
	int N = 3; // кол-во уровней кубита
	double val = 2 * PI * 1e9;
	double tstep = 1e-14; // time grid step
	// main qubit frequencies
	double w1 = 5.0 * (2 * PI) * 1e9; // Частота внешнего управляющего поля
	double w2 = 5.2 * (2 * PI) * 1e9; // Частота внешнего управляющего поля
	// anharmonicities
	double mu1 = 0.25 * (2 * PI) * 1e9; // Параметр нелинейности первого кубита
	double mu2 = 0.4 * (2 * PI) * 1e9; // Параметр нелинейности первого кубита
	double g = 0.02 * (2 * PI) * 1e9; // параметр взаимодействия между кубитами
	
	// qubit capacities
	double Cq1 = 1e-12;
	double Cq2 = 1e-12;

	// connection capacities
	double Cc1 = 4.9e-16;
	double Cc2 = 4e-16;

	// pulse generation frequencies
	double wg1 = w2;
	double wg2 = w2;
	double tau = 4 * 1e-12; // Длительность импульса
	double phi = 0; // phase (number of grid steps paused on Q2)

	// wait time after pulse
	int waitq1 = 0;
	int waitq2 = 0;

	string init = "00"; // initial condition
	string operation = "h0"; // required operation (for fidelity calculation)
	TwoQubitsConstantsDescriptor desc(N, 1000, 0, val, tstep, w1, w2, mu1, mu2, g, Cq1, Cq2,
		Cc1, Cc2, wg1, wg2, tau, phi, waitq1, waitq2, init, operation, 3);

	vector<int> seq(1000, 1);
	
	TwoQubitsKernel kernel(desc);
	kernel.Fidelity(seq);
	double start = omp_get_wtime();
	auto res = kernel.Fidelity(seq);
	start = omp_get_wtime() - start;
	cout << "TIME ELAPSED: " << start << '\n';
	cout.precision(20);
	cout << "F: " << res.fidelity << '\n';
	for (auto& i : res.probs) {
		cout << fixed << setprecision(20) << i.first << ' ' << i.second << '\n';
	}
	return 0;
}



