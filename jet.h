#ifndef JET_H
#define JET_H

#include "lepton.h"

class TLorentzVector;

class Jet : public Lepton {
	public :
		virtual TLorentzVector lorentzVector() const;
		bool passesPUID() const;
};

#endif // JET_H
