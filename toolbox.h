#ifndef TOOLBOX_H
#define TOOLBOX_H

#include <string>
#include <vector>
#include <iostream>
#include <functional>

class Event;

double min( double a, double b );
double max( double a, double b );
double deltaPhi( double phi1, double phi2 );
double deltaR( double eta1, double phi1, double eta2, double phi2);
std::string double2string(double num);
double string2double(const std::string & num);
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
double calculateParA(double a1, double b1, double c1, double b2, double c2, double x);
double calculateParD(double a1, double b1, double c1, double d1, double a2, double b2, double c2, double x);
double myFunc(double a, double b, double c, double d, double x);
double ptFunc(double a1, double b1, double c1, double d1, double b2, double c2,
		double b3, double c3, double b4, double c4, double zpt);
unsigned evCategory(int nhardjet, int nsoftjet, double delEta, double mjj, bool isPhotonSample);
std::vector<std::string> tokenize(std::string text, char token);
std::string encode(const std::string & str);

class EventAdr;
namespace std {
	template <>
		class hash<EventAdr>{
			public :
				size_t operator()(EventAdr ev) const;
		};
}

class EventAdr {
	private :
		unsigned run;
		unsigned lumi;
		unsigned event;
	public :
		EventAdr() : run(0), lumi(0), event(0) {};
		EventAdr(unsigned r, unsigned l, unsigned e) : run(r), lumi(l), event(e) {};
		bool operator==(const EventAdr & ev) const {
			return run == ev.run && lumi == ev.lumi && event == ev.event;
		}
		bool operator<(const EventAdr & ev) const {
			if (run == ev.run) {
				if (lumi == ev.lumi)
					return event < ev.event;
				else
					return lumi < ev.lumi;
			} else
				return run < ev.run;
		}
		friend std::ostream & operator<<(std::ostream & os, const EventAdr & ev);
		friend size_t std::hash<EventAdr>::operator()(EventAdr) const;
};

std::ostream & operator<<(std::ostream & os, const EventAdr & ev);

#endif // EVENT_H
