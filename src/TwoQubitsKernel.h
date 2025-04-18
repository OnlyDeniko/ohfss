#pragma once
#include <algorithm>
#include <cassert>
#include <complex>
#include <numeric>
#include <map>
#include <set>
#include <vector>
#include "ConstantsDescriptor.h"
#include "linalg.h"

using namespace std;

static const double PI = acos(-1.0);
static const double h = 1.054571817e-34;// planck constant
static const double F0 = 2.0678338e-15;// magnetic flux quantum

class TwoQubitsKernel {
private:
	int L, MM;
	const complex<double> zero = { 0, 0 };
	vector<complex<double>>
		a1, a2, Identity, MMIdentity, V1, V2, HQ1, HQ2, aa;
	vector<complex<double>> tmp1, tmp2, tmp3, tmp4, tmp5;
	vector<int> IndexEigValuesAndVectors;
	vector<complex<double>> E10, E01;
	vector<complex<double>> WF_init, WF00, WF10, WF01, WF11, WF20, WF02;
	vector<complex<double>> H10, H01, H11, Hm10, H0m1, Hm11, H1m1, Hm1m1;
	vector<complex<double>> EigVectorsL, EigVectorsR, EigValues;
	vector<complex<double>> ReducedEigVectorsL, ReducedEigVectorsR, ReducedEigValues;
	vector<complex<double>> mtmp1, mtmp2, mtmp3;
	vector<complex<double>> ltmp1, ltmp2, ltmp3, ltmp4, ltmp5, H00;
	vector<complex<double>> Hint, Hred;
	vector<complex<double>> ReducingEigVectors, ReducingEigVectorsT;
	vector<vector<complex<double>>> operators;
	TwoQubitsConstantsDescriptor config;

	void fillIdentity(vector<complex<double>>& A, int dim);
	void matadd(const vector<complex<double>>& A, const vector<complex<double>>& B, vector<complex<double>>& Res);
	void matsub(const vector<complex<double>>& A, const vector<complex<double>>& B, vector<complex<double>>& Res);
	void vsMul(const vector<complex<double>>& A, complex<double> b, vector<complex<double>>& Res);
	void mvMul(const vector<complex<double>>& A, const vector<complex<double>>& B, vector<complex<double>>& Res);
	void kMul(const vector<complex<double>>& A, const vector<complex<double>>& B, vector<complex<double>>& Res, int dim);
	void ctranspose(const vector<complex<double>>& A, vector<complex<double>>& Res, int dim);
	void transpose(const vector<complex<double>>& A, vector<complex<double>>& Res, int dim);
	complex<double> trace(const vector<complex<double>>& A, int dim);
	void _separate_sequence(const vector<int>& sequence, vector<int>& seq1, vector<int>& seq2);
	vector<complex<double>> _prepare_control_operator(
		vector<complex<double>>& V,
		int Vdim,
		vector<complex<double>>& WFs,
		int WFdim1,
		int WFdim2
	);
	vector<pair<int, int>> _combine_sequences(
		const vector<int>& seq1,
		const vector<int>& seq2
	);
public:
	TwoQubitsKernel(const TwoQubitsConstantsDescriptor& _config);

	struct Spectrum {
		vector<double> Energies;
		vector<int> State_q1, State_q2;

		Spectrum(
			const vector<double>& _energies = {},
			const vector<int>& _state_q1 = {},
			const vector<int>& _state_q2 = {}
		) : Energies(_energies), State_q1(_state_q1), State_q2(_state_q2) {}
	};
	Spectrum spectrum, reducedSpectrum;
	void GetSpectrum(
		vector<complex<double>>& a1,
		vector<complex<double>>& a2
	);
	void ReduceBasis();

	struct FidelityResult {
		double fidelity;
		vector<double> Probabilities;
		Spectrum spec;
		vector<complex<double>> U;

		FidelityResult(double _fidelity, const vector<double>& _probabilities, Spectrum _spec, const vector<complex<double>>& _U) :
			fidelity(_fidelity), Probabilities(_probabilities), spec(_spec), U(_U) {}
	};

	FidelityResult Fidelity(const vector<int>& sequence);
};
