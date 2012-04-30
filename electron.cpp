#include "electron.h"
#include <cmath>
#include "TLorentzVector.h"

Electron::Electron( float px_, float py_, float pz_, float en_, float ptErr_, float ecalIso_, float hcalIso_,
		float trkIso_, float gIso_, float chIso_, float puchIso_, float nhIso_, int id_, int genid_,
		float ensf_, float ensferr_, float d0_, float dZ_, float trkpt_, float trketa_, float trkphi_,
		float trkchi2_, float trkValidPixelHits_, float trkValidTrackerHits_, float trkLostInnerHits_,
		int idbits_, float hoe_, float dphiin_, float detain_, float sihih_, float sipip_, float r9_,
		float sce_, float sceta_, float scphi_, float e2x5max_, float e1x5_, float e5x5_, float h2te_,
		float h2tebc_, float ooemoop_, float fbrem_, float eopin_ ) :
		Lepton( px_, py_, pz_, en_, ptErr_, ecalIso_, hcalIso_, trkIso_, gIso_, chIso_, puchIso_, nhIso_,
			id_, genid_, ensf_, ensferr_, d0_, dZ_, trkpt_, trketa_, trkphi_, trkchi2_,
			trkValidPixelHits_, trkValidTrackerHits_, trkLostInnerHits_ ),
		idbits(idbits_),
		hoe(hoe_),
		dphiin(dphiin_),
		detain(detain_),
		sihih(sihih_),
		sipip(sipip_),
		r9(r9_),
		sce(sce_),
		sceta(sceta_),
		scphi(scphi_),
		e2x5max(e2x5max_),
		e1x5(e1x5_),
		e5x5(e5x5_),
		h2te(h2te_),
		h2tebc(h2tebc_),
		ooemoop(ooemoop_),
		fbrem(fbrem_),
		eopin(eopin_) {
}

bool Electron::isEB() const {
	return fabs(sceta) < 1.479;
}

bool Electron::isEE() const {
	return fabs(sceta) > 1.479 && fabs(sceta) < 2.5;
}

bool Electron::isInCrack() const {
	double abseta = fabs(lorentzVector().Eta());
	return (abseta > 1.4442 && abseta < 1.566);
}
double Electron::effAreaMC(double eta) {
	double abseta = fabs(eta);
	if (abseta < 1.0)
		return 0.11;
	else if (abseta < 1.479)
		return 0.13;
	else if (abseta < 2.0)
		return 0.089;
	else if (abseta < 2.2)
		return 0.13;
	else if (abseta < 2.3)
		return 0.15;
	else if (abseta < 2.4)
		return 0.16;
	else
		return 0.19;
}

double Electron::effAreaDATA(double eta) {
	double abseta = fabs(eta);
	if (abseta < 1.0)
		return 0.10;
	else if (abseta < 1.479)
		return 0.12;
	else if (abseta < 2.0)
		return 0.085;
	else if (abseta < 2.2)
		return 0.11;
	else if (abseta < 2.3)
		return 0.12;
	else if (abseta < 2.4)
		return 0.12;
	else
		return 0.13;
}

bool Electron::isPFIsolatedVeto(double rho) const {
	double pfIsol = (chIso + nhIso + gIso) / lorentzVector().Pt() - effAreaMC(lorentzVector().Eta()) * rho;
//	if (isEB())
		return pfIsol < 0.15;
//	else if (isEE())
//		return pfIsol < 0.15;
//	return false;
}

bool Electron::isPFIsolatedLoose(double rho) const {
	double pfIsol = (chIso + nhIso + gIso) / lorentzVector().Pt() - effAreaMC(lorentzVector().Eta()) * rho;
//	if (isEB())
		return pfIsol < 0.15;
//	else if (isEE())
//		return pfIsol < 0.15;
//	return false;
}

bool Electron::isPFIsolatedMedium(double rho) const {
	double pfIsol = (chIso + nhIso + gIso) / lorentzVector().Pt() - effAreaMC(lorentzVector().Eta()) * rho;
//	if (isEB())
		return pfIsol < 0.15;
//	else if (isEE())
//		return pfIsol < 0.15;
//	return false;
}

bool Electron::isPFIsolatedTight(double rho) const {
	double pfIsol = (chIso + nhIso + gIso) / lorentzVector().Pt() - effAreaMC(lorentzVector().Eta()) * rho;
//	if (isEB())
		return pfIsol < 0.10;
//	else if (isEE())
//		return pfIsol < 0.10;
//	return false;
}

bool Electron::passesTightID() const {
	if (isEB()) {
		if (detain < 0.004 && dphiin < 0.03 && sihih < 0.01 && hoe < 0.12 && d0 < 0.02 && dZ < 0.1
				&& fabs( ooemoop ) < 0.05 && ((idbits >> 5) & 0x1) && trkLostInnerHits == 0)

			return true;
	} else if (isEE()) {
		if (detain < 0.005 && dphiin < 0.02 && sihih < 0.03 && hoe < 0.10 && d0 < 0.02 && dZ < 0.1
				&& fabs( ooemoop ) < 0.05 && ((idbits >> 5) & 0x1) && trkLostInnerHits == 0)
			return true;
	}
	return false;
}
