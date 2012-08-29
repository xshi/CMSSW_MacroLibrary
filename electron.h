#include "lepton.h"

#ifndef ELECTRON_H
#define ELECTRON_H

class Electron : public Lepton {
	public :
		int idbits;
		float hoe;
		float dphiin;
		float detain;
		float sihih;
		float sipip;
		float r9;
		float sce;
		float sceta;
		float scphi;
		float e2x5max;
		float e1x5;
		float e5x5;
		float h2te;
		float h2tebc;
		float ooemoop;
		float fbrem;
		float eopin;
		float dEtaCalo;
		float kfchi2;
		float kfhits;
		float etawidth;
		float phiwidth;
		float e1x5e5x5;
		float preShowerOverRaw;
		float eopout;

		Electron( float px_, float py_, float pz_, float en_, float ptErr_, float ecalIso_, float hcalIso_,
				float trkIso_, float gIso_, float chIso_, float puchIso_, float nhIso_, int id_, int genid_,
				float ensf_, float ensferr_, float d0_, float dZ_, float ip3d_, float trkpt_, float trketa_, float trkphi_,
				float trkchi2_, float trkValidPixelHits_, float trkValidTrackerHits_, float trkLostInnerHits_,
				int idbits_, float hoe_, float dphiin_, float detain_, float sihih_, float sipip_, float r9_,
				float sce_, float sceta_, float scphi_, float e2x5max_, float e1x5_, float e5x5_, float h2te_,
				float h2tebc_, float ooemoop_, float fbrem_, float eopin_, float dEtaCalo_, float kfchi2_,
				float kfhits_, float etawidth_, float phiwidth_, float e1x5e5x5_, float preShowerOverRaw_,
				float eopout_ );
		bool isEB() const;
		bool isEE() const;
		bool isInCrack() const;
		static double effAreaMC(double eta);
		static double effAreaDATA(double eta);
		double pfIsolation(double rho, bool isData) const;
		bool isPFIsolatedVeto(double rho, bool isData) const;
		bool isPFIsolatedLoose(double rho, bool isData) const;
		bool isPFIsolatedMedium(double rho, bool isData) const;
		bool isPFIsolatedTight(double rho, bool isData) const;
		bool passesVetoID() const;
		bool passesLooseID() const;
		bool passesMediumID() const;
		bool passesTightID() const;
		bool passesTightTriggerID() const;
		bool passes2011ID() const;
		bool passesMvaTriggerPreselection() const;
};

#endif // ELECTRON_H
