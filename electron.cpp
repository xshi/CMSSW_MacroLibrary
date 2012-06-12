#include "electron.h"
#include <cmath>
#include "TLorentzVector.h"

using std::max;

Electron::Electron( float px_, float py_, float pz_, float en_, float ptErr_, float ecalIso_, float hcalIso_,
		float trkIso_, float gIso_, float chIso_, float puchIso_, float nhIso_, int id_, int genid_,
		float ensf_, float ensferr_, float d0_, float dZ_, float ip3d_, float trkpt_, float trketa_, float trkphi_,
		float trkchi2_, float trkValidPixelHits_, float trkValidTrackerHits_, float trkLostInnerHits_,
		int idbits_, float hoe_, float dphiin_, float detain_, float sihih_, float sipip_, float r9_,
		float sce_, float sceta_, float scphi_, float e2x5max_, float e1x5_, float e5x5_, float h2te_,
		float h2tebc_, float ooemoop_, float fbrem_, float eopin_, float dEtaCalo_, float kfchi2_,
		float kfhits_, float etawidth_, float phiwidth_, float e1x5e5x5_, float preShowerOverRaw_,
		float eopout_ ) :
	Lepton( px_, py_, pz_, en_, ptErr_, ecalIso_, hcalIso_, trkIso_, gIso_, chIso_, puchIso_, nhIso_,
			id_, genid_, ensf_, ensferr_, d0_, dZ_, ip3d_, trkpt_, trketa_, trkphi_, trkchi2_,
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
		eopin(eopin_),
		dEtaCalo(dEtaCalo_),
		kfchi2(kfchi2_),
		kfhits(kfhits_),
		etawidth(etawidth_),
		phiwidth(phiwidth_),
		e1x5e5x5(e1x5e5x5_),
		preShowerOverRaw(preShowerOverRaw_),
		eopout(eopout_) {
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

double Electron::pfIsolation(double rho) const {
	double rhoPr = max(rho, 0.0);
	double eta = fabs(sceta);
	double nIso = max(nhIso + gIso - rhoPr * effAreaMC(eta), 0.0);
	return (nIso + chIso) / lorentzVector().Pt();
}

bool Electron::isPFIsolatedVeto(double rho) const {
	double pfIsol = pfIsolation(rho);
	if (isEB())
		return pfIsol < 0.15;
	else if (isEE())
		return pfIsol < 0.15;
	return false;
}

bool Electron::isPFIsolatedLoose(double rho) const {
	TLorentzVector lv = lorentzVector();
	double pfIsol = pfIsolation(rho);
	if (isEB())
		return pfIsol < 0.15;
	else if (isEE()) {
		if (lv.Pt() > 20)
			return pfIsol < 0.15;
		else
			return pfIsol < 0.10;
	}
	return false;
}

bool Electron::isPFIsolatedMedium(double rho) const {
	TLorentzVector lv = lorentzVector();
	double pfIsol = pfIsolation(rho);
	if (isEB())
		return pfIsol < 0.15;
	else if (isEE()) {
		if (lv.Pt() > 20)
			return pfIsol < 0.15;
		else
			return pfIsol < 0.10;
	}
	return false;
}

bool Electron::isPFIsolatedTight(double rho) const {
	TLorentzVector lv = lorentzVector();
	double pfIsol = pfIsolation(rho);
	if (isEB())
		return pfIsol < 0.10;
	else if (isEE()) {
		if (lv.Pt() > 20)
			return pfIsol < 0.10;
		else
			return pfIsol < 0.07;
	}
	return false;
}

bool Electron::passesVetoID() const {
	if (isEB()) {
		if (fabs(detain) < 0.007 && fabs(dphiin) < 0.8 && sihih < 0.01 && hoe < 0.15 && d0 < 0.04 && dZ < 0.2)
			return true;
	} else if (isEE()) {
		if (fabs(detain) < 0.01 && fabs(dphiin) < 0.7 && sihih < 0.03 && d0 < 0.04 && dZ < 0.4)
			return true;
	}
	return false;
}

bool Electron::passesLooseID() const {
	if (isEB()) {
		if (fabs(detain) < 0.007 && fabs(dphiin) < 0.15 && sihih < 0.01 && hoe < 0.12 && d0 < 0.02 && dZ < 0.2
				&& fabs( ooemoop ) < 0.05 && ((idbits >> 5) & 0x1) && trkLostInnerHits <= 1)
			return true;
	} else if (isEE()) {
		if (fabs(detain) < 0.009 && fabs(dphiin) < 0.10 && sihih < 0.03 && hoe < 0.10 && d0 < 0.02 && dZ < 0.2
				&& fabs( ooemoop ) < 0.05 && ((idbits >> 5) & 0x1) && trkLostInnerHits <= 1)
			return true;
	}
	return false;
}

bool Electron::passesMediumID() const {
	if (isEB()) {
		if (fabs(detain) < 0.004 && fabs(dphiin) < 0.06 && sihih < 0.01 && hoe < 0.12 && d0 < 0.02 && dZ < 0.1
				&& fabs( ooemoop ) < 0.05 && ((idbits >> 5) & 0x1) && trkLostInnerHits <= 1)
			return true;
	} else if (isEE()) {
		if (fabs(detain) < 0.007 && fabs(dphiin) < 0.06 && sihih < 0.03 && hoe < 0.10 && d0 < 0.02 && dZ < 0.1
				&& fabs( ooemoop ) < 0.05 && ((idbits >> 5) & 0x1) && trkLostInnerHits <= 1)
			return true;
	}
	return false;
}

bool Electron::passesTightID() const {
	if (isEB()) {
		if (fabs(detain) < 0.004 && fabs(dphiin) < 0.03 && sihih < 0.01 && hoe < 0.12 && d0 < 0.02 && dZ < 0.1
				&& fabs( ooemoop ) < 0.05 && ((idbits >> 5) & 0x1) && trkLostInnerHits <= 0)
			return true;
	} else if (isEE()) {
		if (fabs(detain) < 0.005 && fabs(dphiin) < 0.02 && sihih < 0.03 && hoe < 0.10 && d0 < 0.02 && dZ < 0.1
				&& fabs( ooemoop ) < 0.05 && ((idbits >> 5) & 0x1) && trkLostInnerHits <= 0)
			return true;
	}
	return false;
}

bool Electron::passesTightTriggerID() const {
	double pt = lorentzVector().Pt();
	if (isEB()) {
		if (fabs(detain) < 0.007 && fabs(dphiin) < 0.15 && sihih < 0.01 && hoe < 0.12
				&& (ecalIso / pt) < 0.2 && (hcalIso / pt) < 0.2 && (trkIso / pt) < 0.2 )
			return true;
	} else if (isEE()) {
		if (fabs(detain) < 0.009 && fabs(dphiin) < 0.10 && sihih < 0.03 && hoe < 0.10
				&& (ecalIso / pt) < 0.2 && (hcalIso / pt) < 0.2 && (trkIso / pt) < 0.2 )
			return true;
	}
	return false;
}

bool Electron::passesMvaTriggerPreselection() const {
	double pt = lorentzVector().Pt();
	if (isEB()) {
		if (sihih < 0.014 && hoe < 0.15 && trkLostInnerHits == 0
				&& (ecalIso / pt) < 0.2 && (hcalIso / pt) < 0.2 && (trkIso / pt) < 0.2 )
			return true;
	} else if (isEE()) {
		if (sihih < 0.035 && hoe < 0.10 && trkLostInnerHits == 0
				&& (ecalIso / pt) < 0.2 && (hcalIso / pt) < 0.2 && (trkIso / pt) < 0.2 )
			return true;
	}
	return false;

}
bool Electron::passesMvaIdWP70(double mvaVal, double rho) {
	double eta = fabs(sceta);
	double pfIsol = pfIsolation(rho);
	if (!passesTightTriggerID())
		return false;
	if (eta < 0.8) {
		if (mvaVal > 0.977 && pfIsol < 0.093)
			return true;
	} else if (eta < 1.479) {
		if (mvaVal > 0.956 && pfIsol < 0.095)
			return true;
	} else {
		if (mvaVal > 0.966 && pfIsol < 0.171)
			return true;
	}
	return false;
}

bool Electron::passesMvaIdWP80(double mvaVal, double rho) {
	double eta = fabs(sceta);
	double pfIsol = pfIsolation(rho);
	if (!passesTightTriggerID())
		return false;
	if (eta < 0.8) {
		if (mvaVal > 0.913 && pfIsol < 0.105)
			return true;
	} else if (eta < 1.479) {
		if (mvaVal > 0.964 && pfIsol < 0.178)
			return true;
	} else {
		if (mvaVal > 0.899 && pfIsol < 0.150)
			return true;
	}
	return false;
}

bool Electron::passesMvaIdWP85(double mvaVal, double rho) {
	double eta = fabs(sceta);
	double pfIsol = pfIsolation(rho);
	if (!passesTightTriggerID())
		return false;
	if (eta < 0.8) {
		if (mvaVal > 0.929 && pfIsol < 0.135)
			return true;
	} else if (eta < 1.479) {
		if (mvaVal > 0.931 && pfIsol < 0.159)
			return true;
	} else {
		if (mvaVal > 0.805 && pfIsol < 0.155)
			return true;
	}
	return false;
}

bool Electron::passesMvaIdWP90(double mvaVal, double rho) {
	double eta = fabs(sceta);
	double pfIsol = pfIsolation(rho);
	if (!passesTightTriggerID())
		return false;
	if (eta < 0.8) {
		if (mvaVal > 0.877 && pfIsol < 0.177)
			return true;
	} else if (eta < 1.479) {
		if (mvaVal > 0.794 && pfIsol < 0.180)
			return true;
	} else {
		if (mvaVal > 0.846 && pfIsol < 0.244)
			return true;
	}
	return false;
}

bool Electron::passesMvaIdWP95(double mvaVal, double rho) {
	double eta = fabs(sceta);
	double pfIsol = pfIsolation(rho);
	if (!passesTightTriggerID())
		return false;
	if (eta < 0.8) {
		if (mvaVal > 0.858 && pfIsol < 0.253)
			return true;
	} else if (eta < 1.479) {
		if (mvaVal > 0.425 && pfIsol < 0.225)
			return true;
	} else {
		if (mvaVal > 0.759 && pfIsol < 0.308)
			return true;
	}
	return false;
}
