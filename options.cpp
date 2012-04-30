#include "options.h"
#include <cstdlib>
#include <iostream>
#include <utility>

using std::pair;
using std::string;
using std::cout;
using std::endl;

Options::Options() {}

void Options::readOptions(int argc, const char * argv[]) {
	for (int i = 1; i < argc; ++i) {
		string option = argv[i];
		if (option.size() < 3) {
			cout << "ERROR: Unknown option: " << option << "!" << endl;
			exit(EXIT_FAILURE);
		}
		if (option[0] != '-' || option[1] != '-') {
			cout << "ERROR: Unknown option: " << option << "!" << endl;
			exit(EXIT_FAILURE);
		}
		option.erase(0, 2);
		size_t pos = option.find('=');
		if (pos != string::npos) {
			string name = option.substr( 0, pos );
			string value = option.substr( pos + 1, option.size() );
			stringOptions.push_back( pair<string, string>(name, value) );
		} else
			boolOptions.push_back( option );
	}
}

bool Options::checkBoolOption(const std::string & name) const {
	for (unsigned i = 0; i < boolOptions.size(); ++i) {
		if (boolOptions[i] == name)
			return true;
	}
	return false;
}

const string & Options::checkStringOption(const std::string & name) const {
	for (unsigned i = 0; i < stringOptions.size(); ++i) {
		if (stringOptions[i].first == name)
			return stringOptions[i].second;
	}
	cout << "ERROR: Unknown option: " << name << "!" << endl;
	exit(EXIT_FAILURE);
}
