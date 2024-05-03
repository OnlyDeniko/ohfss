#include "TwoQubitsKernel.h"

void TwoQubitsKernel::fillIdentity(vector<complex<double>>& A, int dim) {
	for (size_t i = 0; i < dim; ++i) {
		A[i * dim + i] = 1;
	}
}

void TwoQubitsKernel::matadd(const vector<complex<double>>& A, const vector<complex<double>>& B, vector<complex<double>>& Res) {
	for (size_t i = 0; i < A.size(); ++i) {
		Res[i] = A[i] + B[i];
	}
}

void TwoQubitsKernel::matsub(const vector<complex<double>>& A, const vector<complex<double>>& B, vector<complex<double>>& Res) {
	for (size_t i = 0; i < A.size(); ++i) {
		Res[i] = A[i] - B[i];
	}
}

void TwoQubitsKernel::vsMul(const vector<complex<double>>& A, complex<double> b, vector<complex<double>>& Res) {
	for (size_t i = 0; i < A.size(); ++i) {
		Res[i] = A[i] * b;
	}
}

void TwoQubitsKernel::mvMul(const vector<complex<double>>& A, const vector<complex<double>>& B, vector<complex<double>>& Res) {
	for (size_t i = 0; i < Res.size(); ++i) {
		Res[i] = 0;
		for (size_t j = 0; j < B.size(); ++j) {
			Res[i] += A[i * B.size() + j] * B[j];
		}
	}
}

void TwoQubitsKernel::kMul(const vector<complex<double>>& A, const vector<complex<double>>& B, vector<complex<double>>& Res, int dim) {
	for (size_t i = 0; i < dim * dim; ++i) {
		for (size_t j = 0; j < dim * dim; ++j) {
			// Res[i][j] = A[i / Size][j / Size] * B[i % Size][j % Size];
			Res[i * dim * dim + j] = A[(i / dim) * dim + j / dim] * B[(i % dim) * dim + j % dim];
		}
	}
}

void TwoQubitsKernel::ctranspose(const vector<complex<double>>& A, vector<complex<double>>& Res, int dim) {
	for (int i = 0; i < dim; i++) {
		for (int j = i; j < dim; j++) {
			Res[i * dim + j] = conj(A[j * dim + i]);
			Res[j * dim + i] = conj(A[i * dim + j]);
		}
	}
}

complex<double> TwoQubitsKernel::trace(const vector<complex<double>>& A, int dim) {
	complex<double> trace = 0;
	for (int i = 0; i < dim; ++i) {
		trace += A[i * dim + i];
	}
	return trace;
}

void TwoQubitsKernel::_separate_sequence(const vector<int>& sequence, vector<int>& seq1, vector<int>& seq2) {
	assert(sequence.size() == config.N1 + config.N2);
	for (int i = 0; i < sequence.size(); ++i) {
		if (i < config.N1) {
			seq1.push_back(sequence[i]);
		}
		else {
			seq2.push_back(sequence[i]);
		}
	}
}

vector<pair<int, int>> TwoQubitsKernel::_combine_sequences(const vector<int>& seq1, const vector<int>& seq2) {
	/*
	% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
	% combines two grids and decides which combined operator  %
	% should be applied on each grid step                     %
	% 1 grid step = 1 grs                                     %
	% % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
	*/
	double T1 = 2 * PI / config.wg1, T2 = 2 * PI / config.wg2;

	int Nw = ceil(config.tau / config.tstep); // pulse width(in grs)
	int NT1 = ceil(T1 / config.tstep) - Nw; // distance b / w pulses on Q1(in grs)
	int NT2 = ceil(T2 / config.tstep) - Nw; // distance b / w pulses on Q2(in grs)
	vector<int> sz(Nw, 0); // array for 0 pulse
	vector<int> sp(Nw, 1); // array for + 1 pulse
	vector<int> sm(Nw, -1); // array for - 1 pulse
	vector<int> s0_1(NT1, 0); // array for distance b / w pulses on Q1
	vector<int> s0_2(NT2, 0); // array for distance b / w pulses on Q2
	vector<pair<int, int>> arr1, arr2;	

	// phase between pulses(in grs)
	/*if (phi > 0)
		phi_arr = zeros(1, phi);
	end*/

	// placing pulsesand the distance b / w them
	int arr1_size = Nw * seq1.size() + NT1 * max(0, (int)seq1.size() - 1);
	int arr2_size = Nw * seq2.size() + NT2 * max(0, (int)seq2.size() - 1);
	for (int j = 0; j < seq1.size(); ++j) {
		if (seq1[j] == 0) {
			arr1.push_back({ 0, Nw });
		}
		else if (seq1[j] == 1) {
			arr1.push_back({ 1, Nw });
		}
		else {
			arr1.push_back({ -1, Nw });
		}
		if (j + 1 != seq1.size()) {
			arr1.push_back({ 0, NT1 });
		}
	}
	for (int j = 0; j < seq2.size(); ++j) {
		if (seq2[j] == 0) {
			arr2.push_back({ 0, Nw });
		}
		else if (seq2[j] == 1) {
			arr2.push_back({ 1, Nw });
		}
		else {
			arr2.push_back({ -1, Nw });
		}
		if (j + 1 != seq2.size()) {
			arr2.push_back({ 0, NT2 });
		}
	}

	// wait time after pulses(in grs)
	if (config.waitq1 != 0) arr1.push_back({ 0, config.waitq1 });
	if (config.waitq2 != 0) arr2.push_back({ 0, config.waitq2 });
	arr1_size += config.waitq1;
	arr2_size += config.waitq2;

	// equalizing pulse strings by adding zeros to the lesser one
	if (arr1_size > arr2_size) {
		arr2.push_back({ 0, arr1_size - arr2_size });
	}
	if (arr1_size < arr2_size) {
		arr1.push_back({ 0, arr2_size - arr1_size });
	}
	
	auto _compress = [&](vector<pair<int, int>>& arr) {
		int index = 0;
		for (int i = 1; i < arr.size(); i++) {
			if (arr[i].first == arr[index].first) {
				arr[index].second += arr[i].second;
			}
			else {
				arr[++index] = arr[i];
			}
		}
		arr.resize(index + 1);
	};
	_compress(arr1);
	_compress(arr2);

	vector<int> arr;
	/*
	combining into one string
	legend:
			0 = -1-1
			1 = 0-1
			2 = 1-1
			3 = -10
			4 = 00
			5 = 10
			6 = -11
			7 = 01
			8 = 11
	*/
	vector<pair<int, int>> freq;

	int index1 = 0, index2 = 0;
	while (index1 < arr1.size() || index2 < arr2.size()) {
		int gg = min(arr1[index1].second, arr2[index2].second);
		freq.push_back({(arr1[index1].first + 1) + (arr2[index2].first + 1) * 3, gg});
		arr1[index1].second -= gg;
		arr2[index2].second -= gg;
		if (arr1[index1].second == 0) index1 += 1;
		if (arr2[index2].second == 0) index2 += 1;
	}
	return freq;
}

TwoQubitsKernel::TwoQubitsKernel(const TwoQubitsConstantsDescriptor& _config) :
	config(_config),
	L(_config.N* _config.N),
	a1(_config.N* _config.N), a2(_config.N* _config.N),
	Identity(_config.N* _config.N), lIdentity(L* L),
	V1(_config.N* _config.N), V2(_config.N* _config.N),
	HQ1(_config.N* _config.N), HQ2(_config.N* _config.N),
	tmp1(_config.N* _config.N), tmp2(_config.N* _config.N),
	tmp3(_config.N* _config.N), tmp4(_config.N* _config.N),
	tmp5(_config.N* _config.N),
	IndexEigValuesAndVectors(_config.N* _config.N),
	H10(L* L),
	H01(L* L),
	H11(L* L),
	Hm10(L* L),
	H0m1(L* L),
	Hm11(L* L),
	H1m1(L* L),
	Hm1m1(L* L),
	EigVectorsL(L* L),
	EigVectorsR(L* L),
	EigValues(L* L),
	ltmp1(L* L), ltmp2(L* L),
	ltmp3(L* L), ltmp4(L* L), ltmp5(L* L),
	H00(L* L),
	Hint(L* L) {

	fill(a1.begin(), a1.end(), zero);
	fill(a2.begin(), a2.end(), zero);
	// second order quantizaion
	for (int i = 1; i < config.N; i++) {
		a1[i * config.N + i - 1] = { sqrt(i), 0 };
		a2[(i - 1) * config.N + i] = { sqrt(i), 0 };
	}
	fillIdentity(Identity, config.N);
	fillIdentity(lIdentity, config.N * config.N);
	aa = linalg::matmul(a1, a2, config.N, config.N, config.N, config.N, config.N, config.N);

	// field operator
	double V0 = F0 / config.tau; // Voltage
	double Amp1 = config.Cc1 * V0 * sqrt(h * config.w1 / (2 * config.Cq1)); //first generator amplitude
	double Amp2 = config.Cc2 * V0 * sqrt(h * config.w2 / (2 * config.Cq2)); //second generator amplitude

	matsub(a2, a1, tmp1);
	vsMul(tmp1, complex<double>{0, Amp1}, V1);
	vsMul(tmp1, complex<double>{0, Amp2}, V2);


	// now to the hamiltonians
	int L = config.N * config.N;
	vsMul(aa, h * config.w1, tmp1);
	vsMul(aa, h * config.w2, tmp3);

	matsub(aa, Identity, tmp5);
	tmp4 = linalg::matmul(aa, tmp5, config.N, config.N, config.N, config.N, config.N, config.N);
	vsMul(tmp4, h * config.mu1 / 2, tmp2);
	matsub(tmp1, tmp2, HQ1);
	vsMul(tmp4, h * config.mu2 / 2, tmp2);
	matsub(tmp3, tmp2, HQ2);

	matadd(a1, a2, tmp2);
	kMul(tmp2, tmp2, ltmp1, config.N);
	vsMul(ltmp1, h * config.g, Hint);


	// no field hamiltonian
	kMul(HQ1, Identity, ltmp2, config.N);
	kMul(Identity, HQ2, ltmp3, config.N);
	matadd(ltmp2, ltmp3, ltmp1);
	matadd(ltmp1, Hint, H00);

	// eigens
	linalg::eig(H00, EigVectorsL, EigVectorsR, EigValues, L);

	iota(IndexEigValuesAndVectors.begin(), IndexEigValuesAndVectors.end(), 0);
	sort(IndexEigValuesAndVectors.begin(), IndexEigValuesAndVectors.end(), [&](int el1, int el2) {
		return EigValues[el1].real() < EigValues[el2].real();
	});
	auto getEigVector = [&](int index) {
		vector<complex<double>> vec(L);
		for (int i = 0; i < L; ++i) {
			vec[i] = EigVectorsR[i * L + IndexEigValuesAndVectors[index]];
		}
		return vec;
	};

	WF00 = getEigVector(0);
	WF10 = getEigVector(1);
	WF01 = getEigVector(2);
	if (config.N == 2) {
		WF11 = getEigVector(3);
	}
	else if (config.N == 3) {
		WF20 = getEigVector(3);
		WF02 = getEigVector(4);
		WF11 = getEigVector(5);
	}

	kMul(V1, Identity, ltmp1, config.N);
	kMul(Identity, V2, ltmp2, config.N);
	vsMul(V1, -1, tmp1);
	kMul(tmp1, Identity, ltmp3, config.N);
	vsMul(V2, -1, tmp2);
	kMul(Identity, tmp2, ltmp4, config.N);

	matadd(H00, ltmp1, H10);
	matadd(H00, ltmp2, H01);
	matadd(H10, ltmp2, H11);
	matadd(H00, ltmp3, Hm10);
	matadd(H00, ltmp4, H0m1);
	matadd(Hm10, ltmp2, Hm11);
	matadd(H10, ltmp4, H1m1);
	matadd(Hm10, ltmp4, Hm1m1);

	operators = {
		linalg::getUMatrix(lIdentity, Hm1m1, config.tstep, h, L),
		linalg::getUMatrix(lIdentity, H0m1, config.tstep, h, L),
		linalg::getUMatrix(lIdentity, H1m1, config.tstep, h, L),
		linalg::getUMatrix(lIdentity, Hm10, config.tstep, h, L),
		linalg::getUMatrix(lIdentity, H00, config.tstep, h, L),
		linalg::getUMatrix(lIdentity, H10, config.tstep, h, L),
		linalg::getUMatrix(lIdentity, Hm11, config.tstep, h, L),
		linalg::getUMatrix(lIdentity, H01, config.tstep, h, L),
		linalg::getUMatrix(lIdentity, H11, config.tstep, h, L)
	};
}

TwoQubitsKernel::FidelityResult TwoQubitsKernel::Fidelity(const vector<int>& sequence) {
	vector<int> seq1, seq2;
	_separate_sequence(sequence, seq1, seq2);

	vector<pair<int, int>> ComressedPulseString = _combine_sequences(seq1, seq2); // {type, frequency}

	vector<complex<double>> U(L * L);
	fillIdentity(U, L);
	vector<complex<double>> WF, updatedWF(L);
	auto getProbability = [&](const vector<complex<double>>& eigVector) {
		complex<double> dot_product = 0;
		for (int j = 0; j < L; ++j) {	
			dot_product += conj(eigVector[j]) * WF[j];
		}
		return norm(dot_product);
	};

	map<pair<int, int>, vector<complex<double>>> compress2matrix;
	for (auto& i : ComressedPulseString) {
		if (compress2matrix.find(i) == compress2matrix.end()) {
			compress2matrix[i] = linalg::matpow(operators[i.first], i.second, L);
		}
	}

	for (int i = 0; i < ComressedPulseString.size(); i++) {
		auto UPulse = compress2matrix[ComressedPulseString[i]];
		U = linalg::matmul(UPulse, U, L, L, L, L, L, L);
		if (i + 1 == ComressedPulseString.size()) {
			if (config.init == "00") {
				WF = WF00;
			}
			else if (config.init == "01") {
				WF = WF01;
			}
			else if (config.init == "10") {
				WF = WF10;
			}
			else if (config.init == "11") {
				WF = WF11;
			}
			else if (config.init == "20") {
				WF = WF20;
			}
			else if (config.init == "02") {
				WF = WF02;
			}
			mvMul(U, WF, updatedWF);
			swap(WF, updatedWF);
		}
	}

	map<string, double> probs = {
		{"00", getProbability(WF00)},
		{"10", getProbability(WF10)},
		{"01", getProbability(WF01)},
		{"11", getProbability(WF11)}
	};
	if (config.N == 3) {
		probs["20"] = getProbability(WF20);
		probs["02"] = getProbability(WF02);
	}

	// THIS IS THE PART THAT CALCULATES FIDELITY(MAY BE WRONG)
	// 1Q ideal gate matrices
	double dth = PI / 2;
	vector<complex<double>> Ypi2 = {
		cos(dth / 2), -sin(dth / 2), 0,
		sin(dth / 2), cos(dth / 2), 0,
		0, 0, 1
	};
	vector<complex<double>> Ypi =
		linalg::matmul(Ypi2, Ypi2, config.N, config.N, config.N, config.N, config.N, config.N);
	// 2Q ideal gate matrices
	vector<complex<double>>
		Y00(L * L), Yh0(L * L), Y0h(L * L),
		Y10(L * L), Y01(L * L), Yhh(L * L),
		Y1h(L * L), Yh1(L * L), Y11(L * L);
	kMul(Identity, Identity, Y00, config.N);
	kMul(Ypi2, Identity, Yh0, config.N);
	kMul(Identity, Ypi2, Y0h, config.N);
	kMul(Ypi, Identity, Y10, config.N);
	kMul(Identity, Ypi, Y01, config.N);
	kMul(Ypi2, Ypi2, Yhh, config.N);
	kMul(Ypi, Ypi2, Y1h, config.N);
	kMul(Ypi2, Ypi, Yh1, config.N);
	kMul(Ypi, Ypi, Y11, config.N);

	vector<complex<double>> Uid;
	if (config.operation == "00") Uid = Y00;
	else if (config.operation == "01") Uid = Y01;
	else if (config.operation == "10") Uid = Y10;
	else if (config.operation == "11") Uid = Y11;
	else if (config.operation == "h0") Uid = Yh0;
	else if (config.operation == "0h") Uid = Y0h;
	else if (config.operation == "h1") Uid = Yh1;
	else if (config.operation == "1h") Uid = Y1h;
	else if (config.operation == "hh") Uid = Yhh;

	vector<complex<double>> conj_Uid(Uid.size());
	ctranspose(Uid, conj_Uid, L);
	auto M = linalg::matmul(conj_Uid, U, L, L, L, L, L, L);
	vector<complex<double>> conj_M(Uid.size());
	ctranspose(M, conj_M, L);
	auto MMconj = linalg::matmul(M, conj_M, L, L, L, L, L, L);
	double F = (abs(trace(MMconj, L)) + norm(trace(M, L))) / (config.N * config.N * (config.N * config.N + 1));
	//cout << F << endl;
	return FidelityResult(F, probs);
}

