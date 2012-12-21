#include "jet.h"
#include "TLorentzVector.h"

Jet::Jet( float px_, float py_, float pz_, float en_, float btag_, float ptgen_, int idbits_ ) :
	px(px_),
	py(py_),
	pz(pz_),
	en(en_),
	btag(btag_),
	genpt(ptgen_),
	idbits(idbits_) {}

TLorentzVector Jet::lorentzVector() const {
	return TLorentzVector(px, py, pz, en);
}

bool Jet::passesPUID() const {
	return (idbits & (0x1 << 14));
}
