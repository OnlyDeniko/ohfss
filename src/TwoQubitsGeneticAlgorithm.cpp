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
	int ones_count = 0;
	for (int i = config.N1; i < config.N2 + config.N1; ++i) {
		ones_count += sequence[i] == 1;
	}
	assert(ones_count <= config.SecondSequenceOnesLimit);
	auto res = _compute_fidelity(sequence);
	return BaseIndividual(sequence, res.fidelity, 0, 1, 0);
}

TwoQubitsKernel::FidelityResult TwoQubitsGeneticAlgorithm::_compute_fidelity(const std::vector<int>& sequence) {
	return kernel.Fidelity(sequence);
}

void TwoQubitsGeneticAlgorithm::CrossoverImpl(std::vector<int>& ls, std::vector<int>& rs) {
	assert(ls.size() == rs.size());
	int len = config.N1 + config.N2;
	array<vector<int>, 2> ones_count = { vector<int>(config.N2, 0), vector<int>(config.N2, 0) };
	ones_count[0][config.N2 - 1] = (ls[len - 1] == 1);
	ones_count[1][config.N2 - 1] = (rs[len - 1] == 1);
	for (int i = len - 2; i >= config.N1; --i) {
		ones_count[0][i - config.N1] = ones_count[0][i + 1 - config.N1] + (ls[i] == 1);
		ones_count[1][i - config.N1] = ones_count[1][i + 1 - config.N1] + (rs[i] == 1);
	}
	assert(ones_count[0][0] <= config.SecondSequenceOnesLimit);
	assert(ones_count[1][0] <= config.SecondSequenceOnesLimit);
		
	while (true) {
		int index = GenerateInt(2, len - 1);
		if (index < config.N1) {
			for (size_t i = index; i < ls.size(); ++i) {
				swap(ls[i], rs[i]);
			}
			return;
		}
		int new_ones_count1 = ones_count[0][0] - ones_count[0][index - config.N1] + ones_count[1][index - config.N1];
		int new_ones_count2 = ones_count[1][0] - ones_count[1][index - config.N1] + ones_count[0][index - config.N1];
		if (new_ones_count1 <= config.SecondSequenceOnesLimit && new_ones_count2 <= config.SecondSequenceOnesLimit) {
			for (size_t i = index; i < len; ++i) {
				swap(ls[i], rs[i]);
			}
			return;
		}
	}
}

void TwoQubitsGeneticAlgorithm::MutationImpl(std::vector<int>& sequence) {
	int ones_count = 0;
	for (int i = config.N1; i < config.N2 + config.N1; ++i) {
		ones_count += sequence[i] == 1;
	}
	double p = 1.0 / sequence.size();
	
	vector<int> tmp;
	if (config.type == 2) tmp = { 0, 1 };
	else tmp = { -1, 0, 1 };

	for (int i = 0; i < sequence.size(); ++i) {
		if (GenerateProbability() < p) {
			if (i >= config.N1 && ones_count > config.SecondSequenceOnesLimit) {
				ones_count -= sequence[i] == 1;
				sequence[i] = 0;
			}
			else {
				do {
					shuffle(tmp.begin(), tmp.end(), randomGenerator);
				} while (tmp[0] == sequence[i]);
				if (i >= config.N1) ones_count -= sequence[i] == 1;
				sequence[i] = tmp[0];
				if (i >= config.N1) ones_count += sequence[i] == 1;
			}
		}
	}
	for (int i = config.N1; i < config.N1 + config.N2 && ones_count > config.SecondSequenceOnesLimit; i++) {
		ones_count -= sequence[i] == 1;
		sequence[i] = 0;
	}
}

bool TwoQubitsGeneticAlgorithm::CheckStopCondition() {
	return 1 - population[0].fidelity < 0.0001;
}
