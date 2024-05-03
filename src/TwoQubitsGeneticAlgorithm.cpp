#include "TwoQubitsGeneticAlgorithm.h"

TwoQubitsGeneticAlgorithm::TwoQubitsGeneticAlgorithm(const std::vector<std::vector<int>>& _sequences, TwoQubitsConstantsDescriptor _config, GeneticHyperParameters _hyperParams) : 
	config(_config), kernel(_config) {
	hyperParams = _hyperParams;
	populationSize = _sequences.size();
	for (size_t i = 0; i < _sequences.size(); ++i) {
		population.push_back(CreateIndividual(_sequences[i]));
	}
}

double TwoQubitsGeneticAlgorithm::getLeak() {
	return population[0].leak;
}

double TwoQubitsGeneticAlgorithm::getFidelity() {
	return population[0].fidelity;
}

array<vector<int>, 2> TwoQubitsGeneticAlgorithm::getSequences() {
	array<vector<int>, 2> ans;
	for (int i = 0; i < population[0].sequence.size(); ++i) {
		ans[i >= config.N1].push_back(population[0].sequence[i]);
	}
	return ans;
}

int TwoQubitsGeneticAlgorithm::getNumberOfCycles() {
	return population[0].numberOfCycles;
}

BaseIndividual TwoQubitsGeneticAlgorithm::CreateIndividual(const std::vector<int>& sequence) {
	auto res = _compute_fidelity(sequence);
	
	return BaseIndividual(sequence, res.fidelity, 0, 1, 0);
}

TwoQubitsKernel::FidelityResult TwoQubitsGeneticAlgorithm::_compute_fidelity(const std::vector<int>& sequence) {
	return kernel.Fidelity(sequence);
}

void TwoQubitsGeneticAlgorithm::CrossoverImpl(std::vector<int>& ls, std::vector<int>& rs) {
	assert(ls.size() == rs.size());
	int index = GenerateInt(2, (int)rs.size() - 1);
	for (size_t i = index; i < ls.size(); ++i) {
		swap(ls[i], rs[i]);
	}
}

void TwoQubitsGeneticAlgorithm::MutationImpl(std::vector<int>& sequence) {
	double p = 1.0 / sequence.size();
	vector<int> tmp;
	for (int i = 0; i < sequence.size(); ++i) {
		if (GenerateProbability() < p) {
			if (config.type == 2) tmp = { 0, 1 };
			else tmp = { -1, 0, 1 };
			do {
				shuffle(tmp.begin(), tmp.end(), randomGenerator);
			} while (tmp[0] == sequence[i]);
			sequence[i] = tmp[0];
		}
	}
}

bool TwoQubitsGeneticAlgorithm::CheckStopCondition() {
	return population[0].fidelity < 0.0001;
}
