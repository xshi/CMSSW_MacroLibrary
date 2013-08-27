#include "photon.h"
#include "toolbox.h"
#include "jet.h"
#include "TLorentzVector.h"
#include <iostream>

using std::vector;
using std::cout;
using std::endl;
using std::cin;

TLorentzVector Photon::lorentzVector() const {
	return TLorentzVector(getVarF("gn_px"), getVarF("gn_py"), getVarF("gn_pz"), getVarF("gn_en"));
}

bool Photon::isEB() const {
	return fabs(getVarF("egn_sceta")) < 1.479;
}

bool Photon::isEE() const {
	return fabs(getVarF("egn_sceta")) > 1.479 && fabs(getVarF("egn_sceta")) < 2.5;
}

bool Photon::isInCrack() const {
	double abseta = fabs(getVarF("egn_sceta"));
	return (abseta > 1.4442 && abseta < 1.566);
}

bool Photon::isSelected(double rho) {
	TLorentzVector vec = lorentzVector();
	bool isMedium = (getVarI("gn_idbits") & (0x1 << 1));
	if ( vec.Pt() > 30 && isEB() && !isInCrack() && isMedium ) {
		return true;
	}
	return false;

//	const vector<float> & pt = *ev.getVectorFloatAdr("Photons_PT");
//	const vector<float> & eta = *ev.getVectorFloatAdr("Photons_superClusterETA");
//	const vector<float> & r9 = *ev.getVectorFloatAdr("Photons_R9");
//	const vector<float> & IsoSumOverEt = *ev.getVectorFloatAdr("Photons_IsoSumOverEt");
//	const vector<float> & IsoSumOverEtWorst = *ev.getVectorFloatAdr("Photons_IsoSumOverEtWorst");
//	const vector<float> & TrkIsoOverEt = *ev.getVectorFloatAdr("Photons_TrkIsoOverEt");
//	const vector<float> & sigIeIe = *ev.getVectorFloatAdr("Photons_sigmaIetaIeta");
//	const vector<float> & hOe = *ev.getVectorFloatAdr("Photons_hadronicOverEm");
//	const vector<float> & delR = *ev.getVectorFloatAdr("Photons_dREleTrack");
//	const vector<float> & nLostHits = *ev.getVectorFloatAdr("Photons_nLostHitsEleTrack");
//
//	photons.clear();
//	const unsigned nPhot = pt.size();
//	for ( unsigned i = 0; i < nPhot; ++i ) {
//		bool cat1 = (fabs(eta[i]) < 1.4442) && (r9[i] > 0.94);
//		bool cat2 = (fabs(eta[i]) < 1.4442) && (r9[i] <= 0.94);
//		bool cat3 = (fabs(eta[i]) > 1.566 && fabs(eta[i]) < 2.5) && (r9[i] > 0.94);
//		bool cat4 = (fabs(eta[i]) > 1.566 && fabs(eta[i]) < 2.5) && (r9[i] <= 0.94);
//		
//		bool passID = false;
//		if ( pt[i] > 55 ) {
//			if ( cat1 ) {
//				if ( IsoSumOverEt[i] < 3.8 &&
//						IsoSumOverEtWorst[i] < 11.7 &&
//						TrkIsoOverEt[i] < 3.5 &&
//						sigIeIe[i] < 0.0106 &&
//						hOe[i] < 0.082 &&
//						r9[i] > 0.94 )
//					passID = true;
//			}
//			if ( cat2 ) {
//				if ( IsoSumOverEt[i] < 2.2 &&
//						IsoSumOverEtWorst[i] < 3.4 &&
//						TrkIsoOverEt[i] < 2.2 &&
//						sigIeIe[i] < 0.0097 &&
//						hOe[i] < 0.062 &&
//						r9[i] > 0.36
//						// && delR[i] < 0.062
//						)
//					passID = true;
//			}
//			if ( cat3 ) {
//				if ( IsoSumOverEt[i] < 1.77 &&
//						IsoSumOverEtWorst[i] < 3.9 &&
//						TrkIsoOverEt[i] < 2.3 &&
//						sigIeIe[i] < 0.028 &&
//						hOe[i] < 0.065 &&
//						r9[i] > 0.94 )
//					passID = true;
//			}
//			if ( cat4 ) {
//				if ( IsoSumOverEt[i] < 1.29 &&
//						IsoSumOverEtWorst[i] < 1.84 &&
//						TrkIsoOverEt[i] < 1.45 &&
//						sigIeIe[i] < 0.027 &&
//						hOe[i] < 0.048 &&
//						r9[i] > 0.32 )
//					passID = true;
//			}
//		}
//		bool passEleVeto = true;
//		if ( delR[i] < 0.5 && nLostHits[i] == 0 )
//			passEleVeto = false;
//		if (passID && passEleVeto)
//			photons.push_back( i );
//	}
}
