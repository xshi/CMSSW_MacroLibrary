#include "muon.h"
#include "TLorentzVector.h"

Muon::Muon( float px_, float py_, float pz_, float en_, float ptErr_, float ecalIso_, float hcalIso_,
		float trkIso_, float gIso_, float chIso_, float puchIso_, float nhIso_, int id_, int genid_,
		float ensf_, float ensferr_, float d0_, float dZ_, float trkpt_, float trketa_, float trkphi_,
		float trkchi2_, float trkValidPixelHits_, float trkValidTrackerHits_, float trkLostInnerHits_,
		int idbits_, int nMatches_, int validMuonHits_ ) :
		Lepton( px_, py_, pz_, en_, ptErr_, ecalIso_, hcalIso_, trkIso_, gIso_, chIso_, puchIso_, nhIso_,
			id_, genid_, ensf_, ensferr_, d0_, dZ_, trkpt_, trketa_, trkphi_, trkchi2_,
			trkValidPixelHits_, trkValidTrackerHits_, trkLostInnerHits_ ),
		idbits(idbits_),
		nMatches(nMatches_),
		validMuonHits(validMuonHits_) {
}

bool Muon::isTrackerMuon() const {
	return (0x1 << 4) & idbits;
}

bool Muon::isGlobalMuon() const {
	return (0x1 << 5) & idbits;
}

bool Muon::isLooseMuon() const {
	return isTrackerMuon() || isGlobalMuon();
}

bool Muon::isTightMuon() const {
	if (isGlobalMuon() && trkchi2 < 10.0 && validMuonHits > 0 && nMatches > 1 && d0 < 0.2 && dZ < 0.5
			&& trkValidPixelHits > 0 && trkValidTrackerHits > 10)
		return true;
	return false;
}

bool Muon::isTrackIsolatedLoose() const {
	return trkIso / lorentzVector().Pt() < 0.1;
}

bool Muon::isTrackIsolatedTight() const {
	return trkIso / lorentzVector().Pt() < 0.05;
}

bool Muon::isPFIsolatedLoose(double effArea, double rho) const {
	double pfIsol = (chIso + nhIso + gIso) / lorentzVector().Pt() - effArea * rho;
	return pfIsol < 0.2;
}

bool Muon::isPFIsolatedTight(double effArea, double rho) const {
	double pfIsol = (chIso + nhIso + gIso) / lorentzVector().Pt() - effArea * rho;
	return pfIsol < 0.12;
}