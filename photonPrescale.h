#ifndef PHOTONPRESCALE
#define PHOTONPRESCALE

#include <string>
#include <map>
#include <vector>

class PhotonPrescale {
	public :
		class RunLumi {
			public :
				RunLumi(unsigned r, unsigned l) {
					run = r;
					lumi = l;
				}
				unsigned getRun() const {
					return run;
				}
				unsigned getLumi() const {
					return lumi;
				}
				bool operator<(const RunLumi & rl) const;
			private :
				unsigned run;
				unsigned lumi;
		};
		class PhotonTrigger {
			public :
				PhotonTrigger(const std::string & trgName, double t) : name(trgName), threshold(t) {};
				std::string getTriggerName() const {
					return name;
				}
				double getThreshold() const {
					return threshold;
				}
				unsigned getPrescale(const RunLumi & rl, std::ostream & os);
				void readInPrescales(const std::string & inputFileName);
			private :
				std::string name;
				double threshold;
				std::map<RunLumi, unsigned> prescales;
		};
		PhotonPrescale() {};
		void addTrigger(const std::string & tN, double th, const std::string & fileName);
		void addTrigger(const std::string & tN, double th);
		unsigned getPrescale(unsigned run, unsigned lumi, double pt, std::ostream & os);
		double nextThreshold(double pt) const;
	private :
		std::vector<PhotonTrigger> triggers;
};

#endif // PHOTONPRESCALE
