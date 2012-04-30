#include "toolbox.h"
#include <cmath>
#include <sstream>

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

