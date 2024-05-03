#pragma once
#include <algorithm>
#include <cassert>
#include <complex>
#include <numeric>
#include <map>
#include <vector>
#include "ConstantsDescriptor.h"
#include "linalg.h"

using namespace std;

static const double PI = acos(-1.0);
static const double h = 1.054e-34;// planck constant
static const double F0 = 2.06e-15;// magnetic flux quantum

class TwoQubitsKernel {
private:
	int L;
	const complex<double> zero = { 0, 0 };
	vector<complex<double>>
		a1, a2, Identity, lIdentity, V1, V2, HQ1, HQ2, aa;
	vector<complex<double>> tmp1, tmp2, tmp3, tmp4, tmp5;
	vector<int> IndexEigValuesAndVectors;
	vector<complex<double>> WF00, WF10, WF01, WF11, WF20, WF02;
	vector<complex<double>> H10, H01, H11, Hm10, H0m1, Hm11, H1m1, Hm1m1;
	vector<complex<double>> EigVectorsL, EigVectorsR, EigValues;
	vector<complex<double>> ltmp1, ltmp2, ltmp3, ltmp4, ltmp5, H00;
	vector<complex<double>> Hint;
	vector<vector<complex<double>>> operators;
	TwoQubitsConstantsDescriptor config;

	void fillIdentity(vector<complex<double>>& A, int dim);
	void matadd(const vector<complex<double>>& A, const vector<complex<double>>& B, vector<complex<double>>& Res);
	void matsub(const vector<complex<double>>& A, const vector<complex<double>>& B, vector<complex<double>>& Res);
	void vsMul(const vector<complex<double>>& A, complex<double> b, vector<complex<double>>& Res);
	void mvMul(const vector<complex<double>>& A, const vector<complex<double>>& B, vector<complex<double>>& Res);
	void kMul(const vector<complex<double>>& A, const vector<complex<double>>& B, vector<complex<double>>& Res, int dim);
	void ctranspose(const vector<complex<double>>& A, vector<complex<double>>& Res, int dim);
	complex<double> trace(const vector<complex<double>>& A, int dim);
	void _separate_sequence(const vector<int>& sequence, vector<int>& seq1, vector<int>& seq2);
	vector<pair<int, int>> _combine_sequences(
		const vector<int>& seq1,
		const vector<int>& seq2
	);
public:
	TwoQubitsKernel(const TwoQubitsConstantsDescriptor& _config);

	struct FidelityResult {
		double fidelity;
		map<string, double> probs;

		FidelityResult(double _fidelity = 0, map<string, double> _probs = {}) :
			fidelity(_fidelity), probs(_probs) {}
	};

	FidelityResult Fidelity(const vector<int>& sequence);
};
