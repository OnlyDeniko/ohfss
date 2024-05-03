#include <fstream>
#include "ArgsPreprocessor.h"
#include "TwoQubitsGeneticAlgorithm.h"

using namespace std;

void Genetic(
	double CrossoverProbability = 0.9, 
	double MutationProbability = 0.5, 
	int MaxIter = 500
) {
	
	int N = 3; // кол-во уровней кубита
	double val = 2 * PI * 1e9;
	double tstep = 1e-14; // time grid step
	// main qubit frequencies
	double w1 = 5.0 * (2 * PI) * 1e9; // „астота внешнего управл€ющего пол€
	double w2 = 5.2 * (2 * PI) * 1e9; // „астота внешнего управл€ющего пол€
	// anharmonicities
	double mu1 = 0.25 * (2 * PI) * 1e9; // ѕараметр нелинейности первого кубита
	double mu2 = 0.4 * (2 * PI) * 1e9; // ѕараметр нелинейности первого кубита
	double g = 0.02 * (2 * PI) * 1e9; // параметр взаимодействи€ между кубитами

	// qubit capacities
	double Cq1 = 1e-12;
	double Cq2 = 1e-12;

	// connection capacities
	double Cc1 = 4.9e-16;
	double Cc2 = 4e-16;

	// pulse generation frequencies
	double wg1 = w2;
	double wg2 = w2;
	double tau = 4 * 1e-12; // ƒлительность импульса
	double phi = 0; // phase (number of grid steps paused on Q2)

	// wait time after pulse
	int waitq1 = 0;
	int waitq2 = 0;

	string init = "00"; // initial condition
	string operation = "h0"; // required operation (for fidelity calculation)
	int N1 = 100, N2 = 0;
	TwoQubitsConstantsDescriptor config(N, N1, N2, val, tstep, w1, w2, mu1, mu2, g, Cq1, Cq2,
		Cc1, Cc2, wg1, wg2, tau, phi, waitq1, waitq2, init, operation, 3);

	vector<vector<int>> seqs(2 * (N1 + N2));
	uniform_int_distribution<> dist(-1, 1);
	random_device rd;
	mt19937 gen(rd());
	for (auto& seq : seqs) {
		seq.resize(N1 + N2);
		for (auto& j : seq) {
			j = dist(gen);
		}
	}
	GeneticHyperParameters hyperParams(CrossoverProbability, MutationProbability, MaxIter);
	TwoQubitsGeneticAlgorithm algo(seqs, config, hyperParams);
	auto exec_time = algo.run();

	string filename = "tmp.txt";

	ofstream fout;
	fout.open(filename, std::ios::app);

	fout << N1 << '\t' << N2 << '\t';
	auto sequences = algo.getSequences();
	for (auto& i : sequences[0]) fout << i;
	fout << '\t';
	for (auto& i : sequences[1]) fout << i;
	fout << '\t';
	// fout << algo.getNumberOfCycles() << '\t';
	fout << algo.getBestIteration() << '\t';
	fout << algo.getFidelity() << '\t';
	fout << exec_time << '\n';
	fout.close();
}

int main(int argc, char** argv) {
	omp_set_num_threads(4);
	auto mp = ArgsPreprocessor::run(argc, argv);
	Genetic();
	return 0;
}