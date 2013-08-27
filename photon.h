#ifndef PHOTON_H
#define PHOTON_H

#include "lepton.h"
#include <vector>


class TLorentzVector;
class Jet;

class Photon : public Lepton {
	public :
		virtual TLorentzVector lorentzVector() const;
		bool isEB() const;
		bool isEE() const;
		bool isInCrack() const;
		double eta() const;
		bool isSelected(double rho);
};

#endif // PHOTON_H
