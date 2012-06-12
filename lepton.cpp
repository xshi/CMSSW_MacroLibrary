#include "lepton.h"
#include "TLorentzVector.h"

Lepton::Lepton( float px_, float py_, float pz_, float en_, float ptErr_, float ecalIso_, float hcalIso_,
				float trkIso_, float gIso_, float chIso_, float puchIso_, float nhIso_, int id_, int genid_,
				float ensf_, float ensferr_, float d0_, float dZ_, float ip3d_, float trkpt_, float trketa_, float trkphi_,
				float trkchi2_, float trkValidPixelHits_, float trkValidTrackerHits_, float trkLostInnerHits_ ) :
	px(px_),
	py(py_),
	pz(pz_),
	en(en_),
	ptErr(ptErr_),
	ecalIso(ecalIso_),
	hcalIso(hcalIso_),
	trkIso(trkIso_),
	gIso(gIso_),
	chIso(chIso_),
	puchIso(puchIso_),
	nhIso(nhIso_),
	id(id_),
	genid(genid_),
	ensf(ensf_),
	ensferr(ensferr_),
	d0(d0_),
	dZ(dZ_),
	ip3d(ip3d_),
	trkpt(trkpt_),
	trketa(trketa_),
	trkphi(trkphi_),
	trkchi2(trkchi2_),
	trkValidPixelHits(trkValidPixelHits_),
	trkValidTrackerHits(trkValidTrackerHits_),
	trkLostInnerHits(trkLostInnerHits_) {
}

TLorentzVector Lepton::lorentzVector() const {
	return TLorentzVector(px, py, pz, en);
}
