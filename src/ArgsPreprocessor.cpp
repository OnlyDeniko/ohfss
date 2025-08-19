#include "ArgsPreprocessor.h"

std::map<std::string, double> ArgsPreprocessor::run(int argc, char** argv) {
	std::map<std::string, double> mp = {
		{"N1", 100},
		{"N2", 100},
		{"type", 2},
		{"cp", 0.9},
		{"mp", 0.5},
		{"max_iter", 500},
		{"ones_limit", 10},
	};
	for (int i = 1; i < argc; i += 2) {
		std::string name = argv[i];
		name = name.substr(2, std::string::npos);
		double value;
		if (name == "type") {
			std::string val = argv[i + 1];
			value = val == "bipolar" ? 3 : 2;
		}
		else {
			value = atof(argv[i + 1]);
		}
		mp[name] = value;
	}
	return mp;
}
