#pragma once
#include <vector>
#include <random>
#include <chrono>
#include <cassert>
#include <omp.h>
#include <limits.h>
#include "GeneticAlgorithmBase.h"
#include "ConstantsDescriptor.h"
#include "TwoQubitsKernel.h"

class TwoQubitsGeneticAlgorithm : public GeneticAlgorithmBase<BaseIndividual> {
public:
	TwoQubitsGeneticAlgorithm(
		const std::vector<std::vector<int>>& _sequences,
		TwoQubitsConstantsDescriptor _config,
		GeneticHyperParameters _hyperParams
	);
	double getLeak();
	double getFidelity();
	array<vector<int>, 2> getSequences();
	int getNumberOfCycles();
private:
	TwoQubitsConstantsDescriptor config;
	TwoQubitsKernel kernel;

	BaseIndividual CreateIndividual(const std::vector<int>& sequence);
	TwoQubitsKernel::FidelityResult _compute_fidelity(const std::vector<int>& sequence);
	void CrossoverImpl(std::vector<int>& ls, std::vector<int>& rs);
	void MutationImpl(std::vector<int>& sequence);
	bool CheckStopCondition();
};
