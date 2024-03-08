#pragma once
#include<string>

using namespace std;

struct ConstantsDescriptor {
public:
	double w01, w12, wt, w, T, Tstep, Theta;
	int numberOfCycles;
	double neededAngle;
	double angleUpperBound;
	int type;
	int regularLen;
	ConstantsDescriptor(
		double _w01, double _w12, double _wt, double _w,
		double _T, double _Tstep, double _Theta,
		double _neededAngle, int _numberOfCycles,
		double _angleUpperBound, int _type, int _regular_len) :
		w01(_w01), w12(_w12), wt(_wt), w(_w),
		T(_T), Tstep(_Tstep), Theta(_Theta),
		neededAngle(_neededAngle), numberOfCycles(_numberOfCycles),
		angleUpperBound(_angleUpperBound), type(_type), regularLen(_regular_len) {}
};

struct TwoQubitsConstantsDescriptor {
public:
	int N = 3; // кол-во уровней кубита
	int N1, N2; // длины двух последовательностей
	double val;
	double tstep; // time grid step
	// main qubit frequencies
	double w1; // „астота внешнего управл€ющего пол€
	double w2; // „астота внешнего управл€ющего пол€
	// anharmonicities
	double mu1; // ѕараметр нелинейности первого кубита
	double mu2; // ѕараметр нелинейности первого кубита
	double g; // параметр взаимодействи€ между кубитами

	// qubit capacities
	double Cq1;
	double Cq2;

	// connection capacities
	double Cc1;
	double Cc2;

	// pulse generation frequencies
	double wg1;
	double wg2;
	double tau; // ƒлительность импульса
	double phi; // phase (number of grid steps paused on Q2)

	// wait time after pulse
	int waitq1;
	int waitq2;

	std::string init; // initial condition
	string operation; // required operation (for fidelity calculation)

	int type; // 2 or 3
	TwoQubitsConstantsDescriptor(
		int _N,
		int _N1,
		int _N2,
		double _val,
		double _tstep,
		double _w1,
		double _w2,
		double _mu1,
		double _mu2,
		double _g,
		double _Cq1,
		double _Cq2,
		double _Cc1,
		double _Cc2,
		double _wg1,
		double _wg2,
		double _tau,
		double _phi,
		int _waitq1,
		int _waitq2,
		string _init,
		string _operation,
		int _type
	) :
	N(_N), N1(_N1), N2(_N2),
	val(_val), tstep(_tstep),
	w1(_w1), w2(_w2), mu1(_mu1), mu2(_mu2),
	g(_g), Cq1(_Cq1), Cq2(_Cq2), Cc1(_Cc1), Cc2(_Cc2),
	wg1(_wg1), wg2(_wg2), tau(_tau), phi(_phi),
	waitq1(_waitq1), waitq2(_waitq2),
	init(_init), operation(_operation), type(_type) {}
};