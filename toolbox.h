#ifndef TOOLBOX_H
#define TOOLBOX_H

#include <string>
#include <vector>

double min( double a, double b );
double max( double a, double b );
double deltaPhi( double phi1, double phi2 );
double deltaR( double eta1, double phi1, double eta2, double phi2);
std::string double2string(double num);
template <typename T>
T min( const std::vector<T> values) {
	T min = 999999999;
	for (unsigned i = 0; i < values.size(); ++i)
		if ( values[i] < min )
			min = values[i];
	return min;
}
template <typename T>
T max(T a, T b) {
	if (a > b)
		return a;
	else
		return b;
}

#endif // EVENT_H
