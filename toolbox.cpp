#include "toolbox.h"
#include <cmath>
#include <sstream>
#include "event.h"
#include <string>
#include <iomanip>

using std::setw;
using std::string;
using std::vector;

double min( double a, double b ) {
	if (a < b)
		return a;
	else
		return b;
}

double max( double a, double b ) {
	if (a > b)
		return a;
	else
		return b;
}

double deltaPhi( double phi1, double phi2 ) {
	double delPhi = phi1 - phi2;
	while (delPhi > M_PI)
		delPhi -= 2 * M_PI;
	while (delPhi <= -M_PI)
		delPhi += 2 * M_PI;
	return std::fabs(delPhi);
}

double deltaR( double eta1, double phi1, double eta2, double phi2) {
	double delPhi = deltaPhi(phi1, phi2);
	double delEta = std::fabs(eta1 - eta2);
	return std::sqrt(delPhi*delPhi + delEta*delEta);
}

std::string double2string(double num) {
	std::stringstream temp;
	temp << num;
	return temp.str();
}

double string2double(const std::string & num) {
	std::stringstream temp;
	temp << num;
	double n;
	temp >> n;
	return n;
}

double calculateParA(double a1, double b1, double c1, double b2, double c2, double x) {
	double logx = log(x);
	return a1*b1*sin(b1*(logx + c1)) / ( b2 * sin(b2*(logx + c2)) );
}

double calculateParD(double a1, double b1, double c1, double d1, double a2, double b2, double c2, double x) {
	double logx = log(x);
	return a1*cos(b1*(logx + c1)) + d1 - a2*cos(b2*(logx + c2));
}

double myFunc(double a, double b, double c, double d, double x) {
	return exp(a * cos( b * (log(x) + c) ) + d);
}

double ptFunc(double a1p, double b1, double c1, double d1, double b2, double c2,
		double b3, double c3, double a4, double b4, double c4, double d4,
		double zpt, bool merge) {

	double x1 = 90;
	double x2 = 160;
	double x3 = 265;

	double a1 = a1p / b1;
	if (zpt < x1) {
		return myFunc(a1, b1, c1, d1, zpt);
	} else {
		double b2p = b1 + b2;
		double c2p = c1 + c2;
		double a2 = calculateParA(a1, b1, c1, b2p, c2p, x1);
		double d2 = calculateParD(a1, b1, c1, d1, a2, b2p, c2p, x1);
		if (zpt < x2) {
			return myFunc(a2, b2p, c2p, d2, zpt);
		} else {
			double b3p = b2p + b3;
			double c3p = c2p + c3;
			double a3 = calculateParA(a2, b2p, c2p, b3p, c3p, x2);
			double d3 = calculateParD(a2, b2p, c2p, d2, a3, b3p, c3p, x2);
			if (zpt < x3) {
				return myFunc(a3, b3p, c3p, d3, zpt);
			} else {
				if (merge) {
					double b4p = b3p + b4;
					double c4p = c3p + c4;
					double a4p = calculateParA(a3, b3p, c3p, b4p, c4p, x3);
					double d4p = calculateParD(a3, b3p, c3p, d3, a4p, b4p, c4p, x3);
					return myFunc(a4p, b4p, c4p, d4p, zpt);
				} else
					return myFunc(a4, b4, c4, d4, zpt);
			}
		}
	}
}

unsigned evCategory(int nhardjet, int nsoftjet, double delEta, double mjj, bool isPhotonSample) {
	bool vbf = (nhardjet >= 2) && (delEta > 4.0) && (mjj > 500.0); 
	bool cat1 = (nhardjet == 0);
	if (isPhotonSample)
		cat1 = cat1 && (nsoftjet > 0);
	bool cat2 = (nhardjet >= 1) && !vbf;
	if ( cat1 )
		return 1;
	else if ( cat2 )
		return 2;
	else if ( vbf )
		return 3;
	else {
		std::cout << nhardjet << std::endl;
		std::cout << nsoftjet << std::endl;
		std::cout << delEta << std::endl;
		std::cout << mjj << std::endl;
		throw string("ERROR: Unrecognized category!");
	}
}

std::vector<string> tokenize(std::string text, char token) {
	vector<string> fields;
	auto pos = text.find_first_of(token);
	while (pos != string::npos) {
		fields.push_back( text.substr(0, pos) );
		text.erase(0, pos + 1);
		pos = text.find_first_of(token);
	}
	fields.push_back( text );
	return fields;
}

std::string encode(const std::string & text) {
	string str(text);
	for (unsigned i = 0; i < str.size(); ++i) {
		if (str[i] == ' ' || str[i] == '(' || str[i] == ')')
			str[i] = '_';
	}
	auto pos = str.find("_");
	while (pos != string::npos) {
		str.erase(pos, 1);
		pos = str.find("_");
	}
	pos = str.find("&&");
	while (pos != string::npos) {
		str.replace(pos, 2, "_AND_");
		pos = str.find("&&");
	}
	pos = str.find("||");
	while (pos != string::npos) {
		str.replace(pos, 2, "_OR_");
		pos = str.find("||");
	}
	pos = str.find(">");
	while (pos != string::npos) {
		str.replace(pos, 1, "_GT_");
		pos = str.find(">");
	}
	pos = str.find("<");
	while (pos != string::npos) {
		str.replace(pos, 1, "_LT_");
		pos = str.find("<");
	}

	if (str.size() > 128)
		str.erase(128);
	return str;
}

size_t std::hash<EventAdr>::operator()(EventAdr ev) const {
	return hash<unsigned>()(ev.run) ^ hash<unsigned>()(ev.lumi) ^ hash<unsigned>()(ev.event);
}

ostream & operator<<(ostream & os, const EventAdr & ev) {
	os << setw(15) << ev.run << setw(15) << ev.lumi << setw(15) << ev.event;
	return os;
}

double metCutValue(double mh) {
	if (mh < 300)
		return 90;
	else if (mh < 450)
		return 100;
	else if (mh < 500)
		return 120;
	else
		return 130;
}

std::pair<double, double> mtCutValue(double mh) {
	if (mh < 250)
		return std::make_pair(150, 225);
	else if (mh < 300)
		return std::make_pair(175, 275);
	else if (mh < 350)
		return std::make_pair(275, 350);
	else if (mh < 400)
		return std::make_pair(325, 375);
	else if (mh < 450)
		return std::make_pair(350, 450);
	else if (mh < 500)
		return std::make_pair(400, 500);
	else if (mh < 550)
		return std::make_pair(400, 700);
	else if (mh < 600)
		return std::make_pair(450, 700);
	else
		return std::make_pair(500, -1);
}	
