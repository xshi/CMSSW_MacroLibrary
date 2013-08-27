#include "lepton.h"

#ifndef ELECTRON_H
#define ELECTRON_H

class Electron : public Lepton {
	public :
		bool isEB() const;
		bool isEE() const;
		bool isInCrack() const;
		double effAreaMC() const;
		double effAreaDATA() const;
		double pfIsolation(double rho, bool isData) const;
		bool isPFIsolatedVeto(double rho, bool isData) const;
		bool isPFIsolatedLoose(double rho, bool isData) const;
		bool isPFIsolatedMedium(double rho, bool isData) const;
		bool isPFIsolatedTight(double rho, bool isData) const;
		bool passesVetoID() const;
		bool passesLooseID() const;
		bool passesMediumID() const;
//		bool passesTightID() const;
};

#endif // ELECTRON_H
