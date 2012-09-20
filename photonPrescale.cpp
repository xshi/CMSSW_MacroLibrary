#include "photonPrescale.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>

using std::map;
using std::ifstream;
using std::stringstream;
using std::string;

bool PhotonPrescale::RunLumi::operator<(const RunLumi & rl) const {
	if (rl.run == run)
		return lumi < rl.lumi;
	return run < rl.run;
}

unsigned PhotonPrescale::PhotonTrigger::getPrescale(const RunLumi & rl) {
	auto pos = prescales.find(rl);
	if (pos != prescales.end()) {
		return pos->second;
	} else {
		std::cout << getTriggerName() << " : " << getThreshold() << " : " << rl.getRun() << " " << rl.getLumi() << std::endl;
		std::cin.get();
		return 0;
		throw string("ERROR: Can't find prescale!");
	}
}

void PhotonPrescale::PhotonTrigger::readInPrescales(const std::string & inFN) {
	ifstream in( inFN.c_str() );
	if (!in.is_open())
		throw string("ERROR: Can't open input file: " + inFN + "!");

	while (in.good()) {
		string line;
		getline(in, line);
		if (!line.size() || line[0] == '#')
			continue;
		stringstream ss;
		ss << line;
		unsigned r;
		unsigned l;
		unsigned p;
		ss >> r >> l >> p;
		prescales[ RunLumi(r, l) ] = p;
	}

	std::cout << "Done reading in prescales for " << name << "! # line: " << prescales.size() << std::endl;
}

void PhotonPrescale::addTrigger(const std::string & tN, double th, const std::string & fileName) {
	triggers.push_back( PhotonTrigger(tN, th) );
	triggers[ triggers.size() - 1 ].readInPrescales(fileName);
	sort( triggers.begin(), triggers.end(), [](const PhotonTrigger & trg1, const PhotonTrigger & trg2) {
			return trg1.getThreshold() < trg2.getThreshold();
			});
}

unsigned PhotonPrescale::getPrescale(unsigned run, unsigned lumi, double pt) {
	for (unsigned i = 0; i < triggers.size(); ++i) {
		if ( (pt >= triggers[i].getThreshold()) && ((i == triggers.size() - 1) || pt < triggers[i + 1].getThreshold()) )
			return triggers[i].getPrescale( RunLumi(run, lumi) );
	}
	std::cout << "WARNING: Can't find proper prescale for p_T = " << pt << "!" << std::endl; 
	return 0;
}
