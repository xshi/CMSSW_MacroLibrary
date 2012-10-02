#include "eventPrinter.h"
#include "event.h"
#include "electron.h"
#include "toolbox.h"
#include "muon.h"
#include <fstream>
#include "TLorentzVector.h" 
#include <iomanip>
#include <sstream>

using std::stringstream;
using std::setw;
using std::ofstream;
using std::string;
using std::cout;
using std::endl;

EventPrinter::EventPrinter(const Event & ev, const string & fN) {
	run = ev.getSVA<int>("run");
	lumi = ev.getSVA<int>("lumi");
	event = ev.getSVA<int>("event");
	if (fN != "") {
		outputFileName = fN;
		output = new ofstream(outputFileName.c_str());
	} else
		output = &cout;
	floatVars.push_back( ev.getSVA<float>("rho") );
	floatVars.push_back( ev.getSVA<float>("rho25") );
	intVars.push_back( ev.getSVA<int>("cat") );
	floatVarsNames.push_back( "rho" );
	floatVarsNames.push_back( "rho25" );
	intVarsNames.push_back( "cat" );
	printEle = false;
	electrons = 0;
	printMu = false;
	muons = 0;
	useList = false;
}

EventPrinter::~EventPrinter() {
	if (outputFileName == "")
		delete output;
}

void EventPrinter::printHeader() const {
	(*output) << setw(10) << "RUN" << setw(10) << "LUMI" << setw(10) << "EVENT";
	for (unsigned i = 0; i < floatVarsNames.size(); ++i)
		(*output) << setw(10) << floatVarsNames[i];
	for (unsigned i = 0; i < intVarsNames.size(); ++i)
		(*output) << setw(10) << intVarsNames[i];
	if (printEle) {
		(*output) << setw(10) << "Ele Pt";
		(*output) << setw(10) << "Ele Eta";
		(*output) << setw(10) << "Ele E";
		(*output) << setw(10) << "Ele px";
		(*output) << setw(10) << "Ele py";
		(*output) << setw(10) << "Ele pz";
		(*output) << setw(10) << "sceta";
		(*output) << setw(10) << "isInCrack";
		(*output) << setw(10) << "MediumID";
		(*output) << setw(10) << "pfIso";
		(*output) << setw(10) << "gIso";
		(*output) << setw(10) << "chIso";
		(*output) << setw(10) << "nhIso";
		(*output) << setw(10) << "isEB";
		(*output) << setw(10) << "detain";
		(*output) << setw(10) << "dphiin";
		(*output) << setw(10) << "sihih";
		(*output) << setw(10) << "hoe";
		(*output) << setw(10) << "d0";
		(*output) << setw(10) << "dZ";
		(*output) << setw(10) << "ooemoop";
		(*output) << setw(10) << "VFP";
		(*output) << setw(10) << "trkLostIn";
	}
	if (printMu) {
		(*output) << setw(10) << "Mu PT";
		(*output) << setw(10) << "Mu ETA";
		(*output) << setw(10) << "Mu PHI";
		(*output) << setw(10) << "nMatches";
		(*output) << setw(10) << "innerTrac";
	}
	(*output) << endl;
}

void EventPrinter::print() const {
	if ( useList && selectedEvents.find( EventAddress(*run, *lumi, *event) ) == selectedEvents.end() )
		return;
	output->precision(2);
	(*output) << std::scientific;
	int lineLength = 0;
	(*output) << setw(10) << (*run) << setw(10) << (*lumi) << setw(10) << (*event);
	lineLength += 30;
	for (unsigned i = 0; i < floatVars.size(); ++i) {
		(*output) << setw(10) << (*floatVars[i]);
		lineLength += 10;
	}
	for (unsigned i = 0; i < intVars.size(); ++i) {
		(*output) << setw(10) << (*intVars[i]);
		lineLength += 10;
	}
	unsigned maxObjects = 0;
	int electronsLength = 0;
	int muonsLength = 0;
	if (printEle) {
		if (!electrons)
			throw string("ERROR: Empty electrons collection!");
		maxObjects = max( maxObjects, electrons->size() );
		electronsLength = 50;
	}
	if (printMu) {
		if (!muons)
			throw string("ERROR: Empty muons collection!");
		maxObjects = max( maxObjects, muons->size() );
		muonsLength = 50;
	}

	for (unsigned i = 0; i < maxObjects; ++i) {
		if (i > 0)
			(*output) << setw(lineLength) << "";
		if (printEle) {
			if (i < electrons->size() ) {
				(*output) << setw(10) << (*electrons)[i].lorentzVector().Pt();
				(*output) << setw(10) << (*electrons)[i].lorentzVector().Eta();
				(*output) << setw(10) << (*electrons)[i].en;
				(*output) << setw(10) << (*electrons)[i].px;
				(*output) << setw(10) << (*electrons)[i].py;
				(*output) << setw(10) << (*electrons)[i].pz;
				(*output) << setw(10) << (*electrons)[i].sceta;
				(*output) << setw(10) << (*electrons)[i].isInCrack();
				(*output) << setw(10) << (*electrons)[i].passesMediumID();
				(*output) << setw(10) << (*electrons)[i].pfIsolation(*floatVars[0], false);
				(*output) << setw(10) << (*electrons)[i].gIso;
				(*output) << setw(10) << (*electrons)[i].chIso;
				(*output) << setw(10) << (*electrons)[i].nhIso;
				(*output) << setw(10) << (*electrons)[i].isEB();
				(*output) << setw(10) << (*electrons)[i].detain;
				(*output) << setw(10) << (*electrons)[i].dphiin;
				(*output) << setw(10) << (*electrons)[i].sihih;
				(*output) << setw(10) << (*electrons)[i].hoe;
				(*output) << setw(10) << (*electrons)[i].d0;
				(*output) << setw(10) << (*electrons)[i].dZ;
				(*output) << setw(10) << (*electrons)[i].ooemoop;
				(*output) << setw(10) << (((*electrons)[i].idbits >> 5) & 0x1);
				(*output) << setw(10) << (*electrons)[i].trkLostInnerHits;
			} else
				(*output) << setw(electronsLength) << "";
		}
		if (printMu) {
			if (i < muons->size() ) {
				(*output) << setw(10) << (*muons)[i].lorentzVector().Pt();
				(*output) << setw(10) << (*muons)[i].lorentzVector().Eta();
				(*output) << setw(10) << (*muons)[i].lorentzVector().Phi();
				(*output) << setw(10) << (*muons)[i].nMatches;
				(*output) << setw(10) << (*muons)[i].innerTrackChi2;
			} else
				(*output) << setw(muonsLength) << "";
		}
		if (i < maxObjects - 1)
			(*output) << endl;
	}
	(*output) << endl;
}

void EventPrinter::readInEvents(const string & inFileName) {
	ifstream input(inFileName.c_str());

	if (!input.is_open())
		throw string("ERROR: Can't open input file: " + inFileName + "!");

	while (input.good()) {
		string tmpLine;
		getline(input, tmpLine);
		if (!tmpLine.size() || tmpLine[0] == '#')
			continue;
		stringstream ss;
		ss << tmpLine;
		unsigned r;
		ss >> r;
		unsigned l;
		ss >> l;
		unsigned e;
		ss >> e;
		selectedEvents.insert( EventAddress(r, l, e) );
	}
	useList = true;
}

void EventPrinter::printSelectedEvents() const {
	for (auto iter = selectedEvents.begin(); iter != selectedEvents.end(); ++iter)
		cout << *iter << endl;
}

EventPrinter::EventAddress::EventAddress() : run(0), lumi(0), event(0) {};

EventPrinter::EventAddress::EventAddress(unsigned r, unsigned l, unsigned e)
	: run(r), lumi(l), event(e) {};

bool EventPrinter::EventAddress::operator==(const EventAddress & ev) const {
	return run == ev.run && lumi == ev.lumi && event == ev.event;
}

bool EventPrinter::EventAddress::operator<(const EventAddress & ev) const {
	if (event == ev.event) {
		if (lumi == ev.lumi)
			return run < ev.run;
		else
			return lumi < ev.lumi;
	} else
		return event < ev.event;
	return run < ev.run && lumi == ev.lumi && event == ev.event;
}

ostream & operator<<(ostream & os, const EventPrinter::EventAddress & ev) {
	os << setw(10) << ev.run << setw(10) << ev.lumi << setw(10) << ev.event;
	return os;
}
