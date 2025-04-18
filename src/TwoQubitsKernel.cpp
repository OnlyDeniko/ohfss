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

void TwoQubitsKernel::transpose(const vector<complex<double>>& A, vector<complex<double>>& Res, int dim) {
	for (int i = 0; i < dim; i++) {
		for (int j = i; j < dim; j++) {
			Res[i * dim + j] = A[j * dim + i];
			Res[j * dim + i] = A[i * dim + j];
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
	vector<pair<int, int>> arr1, arr2;	
	
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

vector<complex<double>> TwoQubitsKernel::_prepare_control_operator(
	vector<complex<double>>& V,
	int Vdim,
	vector<complex<double>>& WFs,
	int WFdim1,
	int WFdim2
) {
	assert(V.size() == Vdim * Vdim);
	assert(WFs.size() == WFdim1 * WFdim2);
	assert(Vdim == WFdim2);

	vector<complex<double>> WFsT(WFs.size());
	for (int i = 0; i < WFdim1; ++i) {
		for (int j = 0; j < WFdim2; ++j) {
			WFsT[j * WFdim1 + i] = WFs[i * WFdim2 + j];
		}
	}
	auto right_part = linalg::matmul(V, WFsT, Vdim, WFdim1, Vdim, Vdim, WFdim1, WFdim1);

	vector<complex<double>> ConjWFs(WFs.size());
	for (int i = 0; i < WFdim1; ++i) {
		for (int j = 0; j < WFdim2; ++j) {
			ConjWFs[i * WFdim2 + j] = conj(WFs[i * WFdim2 + j]);
		}
	}
	return linalg::matmul(ConjWFs, right_part, WFdim1, WFdim1, WFdim2, WFdim2, WFdim1, WFdim1);
}

void TwoQubitsKernel::GetSpectrum(vector<complex<double>>& a1, vector<complex<double>>& a2){
	double E00_est, E10_est, E01_est, E11_est, E20_est, E02_est, E21_est, E12_est, E22_est;
	vector<double> Energies;
	vector<string> States;
	if (config.N == 2){
		E00_est = 0;
		E10_est = config.w1 / (2 * PI);
		E01_est = config.w2 / (2 * PI);
		E11_est = (config.w1 + config.w2) / (2 * PI);
		Energies = {E00_est, E10_est, E01_est, E11_est};
		States = {"00", "10", "01", "11"};
	} else if (config.N == 3){
		E00_est = 0;
		E10_est = config.w1 / (2 * PI);
		E01_est = config.w2 / (2 * PI);
		E11_est = (config.w1 + config.w2) / (2 * PI);
		E20_est = (2 * config.w1 - config.mu1) / (2 * PI);
		E02_est = (2 * config.w2 - config.mu2) / (2 * PI);
		E21_est = (2 * config.w1 + config.w2 - config.mu1) / (2 * PI);
		E12_est = (2 * config.w2 + config.w1 - config.mu2) / (2 * PI);
		E22_est = (2 * config.w1 + 2 * config.w2 - config.mu1 - config.mu2) / (2 * PI);
		Energies = {E00_est, E10_est, E01_est, E11_est, E20_est, E02_est, E21_est, E12_est, E22_est};
		States = {"00", "10", "01", "11", "20", "02", "21", "12", "22"};
	} else {
		assert(0);
	}

	// now to the hamiltonians
	aa = linalg::matmul(a1, a2, config.N, config.N, config.N, config.N, config.N, config.N);
	// linalg::print_matrix("aa", config.N, config.N, aa, config.N);
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

	linalg::eig(H00, EigVectorsL, EigVectorsR, EigValues, L);

	vector<int> indices(L);
	iota(indices.begin(), indices.end(), 0);
	std::sort(indices.begin(), indices.end(), [&](int a, int b) {
		return EigValues[a].real() > EigValues[b].real();
	});
	
	for (int i = 0; i < L; ++i) {
		for (int j = 0; j < L; ++j) {
			ltmp1[i + j * L] = EigVectorsL[indices[i] + j * L];
		}
	}
	for (int i = 0; i < L; ++i) {
		for (int j = 0; j < L; ++j) {
			ltmp2[i + j * L] = EigVectorsR[indices[i] + j * L];
		}
	}
	for (int i = 0; i < L; ++i) {
		for (int j = 0; j < L; ++j) {
			ltmp3[i + j * L] = EigValues[indices[i] + j * L];
		}
	}
	swap(EigVectorsL, ltmp1);
	swap(EigVectorsR, ltmp2);
	swap(EigValues, ltmp3);

	for (int i = 0; i < Energies.size(); ++i) {
		Energies[i] = EigValues[i].real() / (2 * PI * h);
	}
	transpose(EigVectorsL, ltmp3, L);

	kMul(aa, Identity, ltmp1, config.N);
	kMul(Identity, aa, ltmp2, config.N);

	vector<int> State_q1(L), State_q2(L);
	vector<complex<double>> wf1(L), wf2(L);
	for (int i = 0; i < L; ++i) {
		for (int j = 0; j < L; ++j) {
			wf1[j] = ltmp3[i * L + j];
			wf2[j] = EigVectorsL[j * L + i];
		}
		mvMul(ltmp1, wf1, tmp1);
		complex<double> check1 = 0;
		for (int j = 0; j < L; ++j) {
			check1 += tmp1[j] * wf2[j];
		}
		mvMul(ltmp2, wf1, tmp1);
		complex<double> check2 = 0;
		for (int j = 0; j < L; ++j) {
			check2 += tmp1[j] * wf2[j];
		}
		State_q1[i] = std::round(check1.real());
		State_q2[i] = std::round(check2.real());
	}

	spectrum = Spectrum(Energies, State_q1, State_q2);
	swap(EigVectorsL, ltmp3);
}

void TwoQubitsKernel::ReduceBasis() {
	vector<int> indices, State_q1, State_q2;
	for (int i = 0; i < L; ++i) {
		if (spectrum.State_q1[i] < config.M && spectrum.State_q2[i] < config.M) {
			indices.push_back(i);
			State_q1.push_back(spectrum.State_q1[i]);
			State_q2.push_back(spectrum.State_q2[i]);
		}
	}
	assert(indices.size() == MM);
	for (int i = 0; i < MM;++i) {
		for (int j = 0; j < L; ++j) {
			ReducingEigVectors[i * L + j] = EigVectorsL[indices[i] * L + j];
		}
	}
	Hred = _prepare_control_operator(H00, L, ReducingEigVectors, MM, L);
	linalg::eig(Hred, ReducedEigVectorsL, ReducedEigVectorsR, ReducedEigValues, MM);

	iota(indices.begin(), indices.end(), 0);
	std::sort(indices.begin(), indices.end(), [&](int a, int b) {
		return ReducedEigValues[a].real() > ReducedEigValues[b].real();
	});

	for (int i = 0; i < MM; ++i) {
		for (int j = 0; j < MM; ++j) {
			mtmp1[i * MM + j] = ReducedEigVectorsL[indices[i] + j * MM];
			mtmp2[i * MM + j] = ReducedEigVectorsR[indices[i] + j * MM];
			mtmp3[i + j * MM] = ReducedEigValues[indices[i] + j * MM];
		}
	}
	swap(ReducedEigVectorsL, mtmp1);
	swap(ReducedEigVectorsR, mtmp2);
	swap(ReducedEigValues, mtmp3);

	vector<double> ReducedEnergies(MM);
	for (int i = 0; i < ReducedEnergies.size(); ++i) {
		ReducedEnergies[i] = ReducedEigValues[i].real() / (2 * PI * h);
	}
	reducedSpectrum = Spectrum(ReducedEnergies, State_q1, State_q2);
}

TwoQubitsKernel::TwoQubitsKernel(const TwoQubitsConstantsDescriptor& _config) :
	config(_config),
	L(_config.N* _config.N), MM(_config.M* _config.M),
	a1(_config.N* _config.N), a2(_config.N* _config.N),
	Identity(_config.N* _config.N), MMIdentity(MM * MM),
	V1(_config.N* _config.N), V2(_config.N* _config.N),
	HQ1(_config.N* _config.N), HQ2(_config.N* _config.N),
	tmp1(_config.N* _config.N), tmp2(_config.N* _config.N),
	tmp3(_config.N* _config.N), tmp4(_config.N* _config.N),
	tmp5(_config.N* _config.N),
	IndexEigValuesAndVectors(_config.N* _config.N),
	H10(MM * MM),
	H01(MM* MM),
	H11(MM* MM),
	Hm10(MM* MM),
	H0m1(MM* MM),
	Hm11(MM* MM),
	H1m1(MM* MM),
	Hm1m1(MM* MM),
	EigVectorsL(L* L),
	EigVectorsR(L* L),
	EigValues(L* L),
	ReducedEigVectorsL(MM * MM),
	ReducedEigVectorsR(MM * MM),
	ReducedEigValues(MM * MM),
	ReducingEigVectors(L * MM),
	ReducingEigVectorsT(L* MM),
	Hred(MM * MM),
	ltmp1(L* L), ltmp2(L* L), ltmp3(L* L), ltmp4(L* L), ltmp5(L* L),
	mtmp1(MM* MM), mtmp2(MM* MM), mtmp3(MM* MM),
	H00(L* L),
	Hint(L* L),
	WF_init(MM) {

	fill(a1.begin(), a1.end(), zero);
	fill(a2.begin(), a2.end(), zero);
	// second order quantizaion
	for (int i = 1; i < config.N; i++) {
		a1[i * config.N + i - 1] = { sqrt(i), 0 };
		a2[(i - 1) * config.N + i] = { sqrt(i), 0 };
	}
	fillIdentity(Identity, config.N);

	GetSpectrum(a1, a2);
	ReduceBasis();

	int _match_index = 0;
	for (int i = 0; i < MM; ++i) {
		string gg = to_string(reducedSpectrum.State_q1[i]) + to_string(reducedSpectrum.State_q2[i]);
		if (gg == config.init) {
			_match_index = i;
			break;
		}
	}
	for (int i = 0; i < MM; ++i) {
		WF_init[i] = ReducedEigVectorsL[_match_index * MM + i];
	}
	// field operator
	double V0 = F0 / config.tau; // Voltage
	double V1 = config.Cc1 * V0 * sqrt(h * config.w1 / (2 * config.Cq1)); //first generator amplitude
	double V2 = config.Cc2 * V0 * sqrt(h * config.w2 / (2 * config.Cq2)); //second generator amplitude

	matsub(a2, a1, tmp1);

	kMul(tmp1, Identity, ltmp1, config.N);
	kMul(Identity, tmp1, ltmp2, config.N);

	E10 = _prepare_control_operator(ltmp1, L, ReducingEigVectors, MM, L); // [MM x MM]
	E01 = _prepare_control_operator(ltmp2, L, ReducingEigVectors, MM, L); // [MM x MM]
	
	vsMul(E10, complex<double>{0, V1}, mtmp1);
	vsMul(E01, complex<double>{0, V2}, mtmp2);
	
	matadd(Hred, mtmp1, H10);
	matadd(Hred, mtmp2, H01);
	matadd(H10, mtmp2, H11);
	matsub(Hred, mtmp1, Hm10);
	matadd(Hred, mtmp2, H0m1);
	matadd(Hm10, mtmp2, Hm11);
	matsub(H10, mtmp2, H1m1);
	matsub(Hm10, mtmp2, Hm1m1);

	fillIdentity(MMIdentity, MM);
	operators = {
		linalg::getUMatrix(MMIdentity, Hm1m1, config.tstep, h, MM),
		linalg::getUMatrix(MMIdentity, H0m1, config.tstep, h, MM),
		linalg::getUMatrix(MMIdentity, H1m1, config.tstep, h, MM),
		linalg::getUMatrix(MMIdentity, Hm10, config.tstep, h, MM),
		linalg::getUMatrix(MMIdentity, Hred, config.tstep, h, MM),
		linalg::getUMatrix(MMIdentity, H10, config.tstep, h, MM),
		linalg::getUMatrix(MMIdentity, Hm11, config.tstep, h, MM),
		linalg::getUMatrix(MMIdentity, H01, config.tstep, h, MM),
		linalg::getUMatrix(MMIdentity, H11, config.tstep, h, MM)
	};
}

TwoQubitsKernel::FidelityResult TwoQubitsKernel::Fidelity(const vector<int>& sequence) {
	vector<int> seq1, seq2;
	_separate_sequence(sequence, seq1, seq2);

	vector<pair<int, int>> ComressedPulseString = _combine_sequences(seq1, seq2); // {type, frequency}
	vector<complex<double>> U(MM * MM);
	fillIdentity(U, MM);
	vector<complex<double>> WF(MM);

	map<pair<int, int>, vector<complex<double>>> compress2matrix;
	for (auto& i : ComressedPulseString) {
		if (compress2matrix.find(i) == compress2matrix.end()) {
			compress2matrix[i] = linalg::matpow(operators[i.first], i.second, MM);
		}
	}

	for (int i = 0; i < ComressedPulseString.size(); i++) {
		auto UPulse = compress2matrix[ComressedPulseString[i]];
		U = linalg::matmul(UPulse, U, MM, MM, MM, MM, MM, MM);
	}
	mvMul(U, WF_init, WF);
	
	ctranspose(ReducedEigVectorsL, mtmp1, MM);
	vector<complex<double>> mtmp(MM);
	mvMul(mtmp1, WF, mtmp);

	vector<double> probs(MM);
	for (int i = 0; i < MM; ++i) {
		probs[i] = norm(mtmp[i]);
	}

	vector<int> indices;
	set<string> basic_states = {
		"00",
		"01",
		"10",
		"11"
	};

	for (int i = 0; i < reducedSpectrum.State_q1.size(); ++i) {
		string cur = to_string(reducedSpectrum.State_q1[i]) + to_string(reducedSpectrum.State_q2[i]);
		if (basic_states.find(cur) != basic_states.end()) {
			indices.push_back(i);
		}
	}
	assert(indices.size() == basic_states.size());

	vector<complex<double>> Ured(4 * 4);
	for (int i = 0; i < 4; ++i) {
		for(int j = 0;j < 4;++j){
			Ured[i * 4 + j] = U[indices[i] * MM + indices[j]];
		}
	}
	// fidelity calculations

	auto calculate_fidelity = [&](vector<complex<double>>& U, vector<complex<double>>& Uid) {
		vector<complex<double>> conj_U(U.size());
		ctranspose(U, conj_U, 4);
		auto conj_U_U = linalg::matmul(conj_U, U, 4, 4, 4, 4, 4, 4);
		auto conj_U_Uid = linalg::matmul(conj_U, Uid, 4, 4, 4, 4, 4, 4);
		double F = (abs(trace(conj_U_U, 4)) + norm(trace(conj_U_Uid, 4))) / 20;
		return F;
	};

	vector<complex<double>> Xid, Yid, Zid, Pid;
	Xid = {
		{0, 0}, {1, 0},
		{1, 0}, {0, 0}
	};
	Yid = {
		{0, 0}, {0, -1},
		{0, 1}, {0, 0}
	};
	Zid = {
		{1, 0}, {0, 0},
		{0, 0}, {-1, 0}
	};

	vector<complex<double>> Id2(2 * 2), Uid(4 * 4);
	double F;
	fillIdentity(Id2, 2);
	if (config.operation == "IX") {
		kMul(Id2, Xid, Uid, 2);
		F = calculate_fidelity(Ured, Uid);
	}
	else if (config.operation == "XI") {
		kMul(Xid, Id2, Uid, 2);
		F = calculate_fidelity(Ured, Uid);
	}
	else if (config.operation == "IY") {
		kMul(Id2, Yid, Uid, 2);
		F = calculate_fidelity(Ured, Uid);
	}
	else if (config.operation == "YI") {
		kMul(Yid, Id2, Uid, 2);
		F = calculate_fidelity(Ured, Uid);
	}
	else if (config.operation == "ZX") {
		kMul(Zid, Xid, Uid, 2);
		F = calculate_fidelity(Ured, Uid);
	}
	else if (config.operation == "XZ") {
		kMul(Xid, Zid, Uid, 2);
		F = calculate_fidelity(Ured, Uid);
	}
	else if (config.operation == "ZY") {
		kMul(Zid, Yid, Uid, 2);
		F = calculate_fidelity(Ured, Uid);
	}
	else if (config.operation == "YZ") {
		kMul(Yid, Zid, Uid, 2);
		F = calculate_fidelity(Ured, Uid);
	}
	else if (config.operation == "YP" || config.operation == "PY" || config.operation == "XP" || config.operation == "PX") {
		double best_F = -1;
		int iterations = 7200;
		for (int i = 0; i < iterations; ++i) {
			Pid = {
				{1, 0}, {0, 0},
				{0, 0}, {cos(2 * PI / iterations * i), sin(2 * PI / iterations * i)}
			};
			if (config.operation == "PY") kMul(Pid, Yid, Uid, 2);
			else if (config.operation == "YP") kMul(Yid, Pid, Uid, 2);
			else if (config.operation == "PX") kMul(Pid, Xid, Uid, 2);
			else if (config.operation == "XP") kMul(Xid, Pid, Uid, 2);

			double F = calculate_fidelity(Ured, Uid);
			if (best_F < F) {
				best_F = F;
			}
		}
		F = best_F;
	}
	return FidelityResult(F, probs, reducedSpectrum, Ured);
}

