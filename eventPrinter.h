#ifndef EVENTPRINTER_H
#define EVENTPRINTER_H

#include<string>
#include<iostream>
#include<memory>
#include<vector>
#include<set>

class Event;
class Electron;
class Muon;

class EventPrinter {
	private :
		class EventAddress {
			private :
				unsigned run;
				unsigned lumi;
				unsigned event;
			public :
				EventAddress();
				EventAddress(unsigned r, unsigned l, unsigned e);
				bool operator==(const EventAddress & ev) const;
				bool operator<(const EventAddress & ev) const;
				friend std::ostream & operator<<(std::ostream & os, const EventAddress & ev);
		};
		const int * run;
		const int * lumi;
		const int * event;
		std::string outputFileName;
		std::ostream * output;
		std::vector<const float *> floatVars;
		std::vector<std::string> floatVarsNames;
		bool printEle;
		const std::vector<Electron> * electrons;
		bool printMu;
		const std::vector<Muon> * muons;
		bool useList;
		std::set<EventAddress> selectedEvents;
		EventPrinter(const EventPrinter &);
		bool operator=(const EventPrinter &);
		friend std::ostream & operator<<(std::ostream & os, const EventAddress & ev);
	public :
		EventPrinter(const Event & ev, const std::string & fN = "");
		~EventPrinter();
		void print() const;
		void printHeader() const;
		void printElectrons() {
			printEle = true;
		}
		void setElectronCollection(const std::vector<Electron> & ele) {
			electrons = &ele;
		}
		void printMuons() {
			printMu = true;
		}
		void setMuonCollection(const std::vector<Muon> & mu) {
			muons = &mu;
		}
		void readInEvents(const std::string & inFileName);
		void printSelectedEvents() const;
};

#endif // EVENTPRINTER_H
