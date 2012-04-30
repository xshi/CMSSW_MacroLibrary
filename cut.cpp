// Standard Libraries
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <iomanip>
// ROOT Libraries
#include <TLorentzVector.h>
// Other
#include "cut.h"
#include "event.h"

using std::setw;
using std::cout;
using std::endl;
using std::string;
using std::stringstream;

string DCut::GetName() const {
	stringstream temp;
	temp << varName_;
	if ( cT_ == GT ) {
		temp << "_GT_";
		temp << GetCutValue();
	} else if ( cT_ == LT ) {
		temp << "_LT_";
		temp << GetCutValue();
	} else if ( cT_ == EQ ) {
		temp << "_EQ_";
		temp << GetCutValue();
	}
	return temp.str();
}

bool Cut::operator()() const {
	if ( desc_.GetType() == DCut::GT )
		return ( var_->getValue() > desc_.GetCutValue() );
	else if ( desc_.GetType() == DCut::LT )
		return ( var_->getValue() < desc_.GetCutValue() );
	else if ( desc_.GetType() == DCut::EQ )
		return ( var_->getValue() == desc_.GetCutValue() );
	else
		return true;
}

string Cut::GetName() const {
	return desc_.GetName();
}

/*
double Cut::Optimize(const Samples & samples, OptConf & conf) {
	cout << "--> Starting optimization of " << GetVarName() << endl;
	unsigned nPoints = conf.GetNpoints();
	double inc = conf.GetIncrement();
	double lPoint = GetCutValue() - ( nPoints / 2 ) * inc;
	if ( !nPoints ) {
		cout << "ERROR: Number of points for optimization has to be greater than 0!" << endl;
		exit( EXIT_FAILURE );
	}
	if (!samples.getNsignalSamples() ) {
		cout << "ERROR: Can't optimize anything without signal sample!" << endl;
		exit( EXIT_FAILURE );
	}
	if (!samples.getNbackgroundSamples() ) {
		cout << "ERROR: Can't optimize anything without background samples!" << endl;
		exit( EXIT_FAILURE );
	}

	double sign = samples.significance();
	double iniVal = GetCutValue();
	double iniSign = sign;
	for (unsigned i = 0; i < nPoints; ++i) {
		Samples tempSamples(samples);
		double tempCV = lPoint + i * inc;
		SetCutValue( tempCV );
		tempSamples.selectEvents( GetVarName(), false );
		if (!i) {
			conf.SetBest( lPoint );
			sign = tempSamples.significance();
		}
		if ( tempSamples.significance() > sign ) {
			sign = tempSamples.significance();
			conf.SetBest( tempCV );
			cout << "----> Found better significance: " << tempSamples.significance() << endl;
		}
	}
	double hPoint = lPoint + (nPoints - 1) * inc;
	bool atLowB = conf.GetBest() == lPoint;
	bool atHighB = conf.GetBest() == hPoint;
	if (sign != iniSign && (atLowB || atHighB) ) {
		double oldSign = -1;
		while ( sign > oldSign ) {
			oldSign = sign;
			Samples tempSamples(samples);
			double tempCV;
			if (atLowB)
				tempCV = conf.GetBest() - inc;
			else
				tempCV = conf.GetBest() + inc;
			SetCutValue( tempCV );
			tempSamples.selectEvents( GetVarName(), false );
			if ( tempSamples.significance() > sign ) {
				sign = tempSamples.significance();
				conf.SetBest( tempCV );
				cout << "----> Found better significance: " << tempSamples.significance() << endl;
			}
		}
	}
	if (iniSign > sign) {
		cout << "WARNING: This cut decreases significance!" << endl;
	}
	if (sign != iniSign)
		SetCutValue( conf.GetBest() );
	else
		SetCutValue( iniVal );
	cout.precision(8);
	cout << "Initial significance : " << iniSign << endl;
	cout << "Best significance    : " << sign << endl;
	cout << "Finished optimization of " << GetVarName() << " with best significance for cut value: " << GetCutValue() << endl;
	return conf.GetBest();
}

bool JetSelector::operator() (const Event & ev, unsigned index) const {
	if (index > max_index) {
		cout << "ERROR: Index out of range!" << endl;
		exit( EXIT_FAILURE );
	}
	double jetPt = ev.getDouble25ArrayValue("RECOJET_PT", index);
	double jetEta = ev.getDouble25ArrayValue("RECOJET_ETA", index);
	double delLep1 = ev.getDouble25ArrayValue("RECOJET_DELTAR1LEP", index);
	double delLep2 = ev.getDouble25ArrayValue("RECOJET_DELTAR2LEP", index);
	if ( jetPt > pt_min && std::fabs(jetEta) < eta_max && delLep1 > 0.5 && delLep2 > 0.5)
		return true;
	return false;
}


bool NlepCut::operator()(const Event & ev) const {
	int isEle1 = ev.getInt5ArrayValue("RECO_isElectron", 3);
	int isEle2 = ev.getInt5ArrayValue("RECO_isElectron", 4);
	int isMu1 = ev.getInt5ArrayValue("RECO_isGlobalMu", 3);
	int isMu2 = ev.getInt5ArrayValue("RECO_isGlobalMu", 4);
	int nEle = ev.getIntVariableValue("RECO_NELE");
	int nMu = ev.getIntVariableValue("RECO_NMU");
	if ( isEle1 == 1 && isEle2 == 1 ) {
		if ( nEle > int(GetValue()) || nMu > 0 )
			return false;
	} else if ( isMu1 == 1 && isMu2 == 1 ) {
		if ( nMu > int(GetValue()) || nEle > 0 )
			return false;
	} else {
		std::cout << "NMU = " << nMu << std::endl;
		std::cout << "NELE = " << nEle << std::endl;
		std::cout << "l1 is GLB = " << isMu1 << std::endl;
		std::cout << "l2 is GLB = " << isMu2 << std::endl;
		std::cout << "l1 is Ele = " << isEle1 << std::endl;
		std::cout << "l2 is Ele = " << isEle2 << std::endl;
		std::cin.get();
	}
	return true;
}

bool ZmassMinCut::operator()(const Event & ev) const {
	if (ev.getDouble5ArrayValue("RECO_MASS", 1) < GetValue())
		return false;
	return true;
}

bool ZmassMaxCut::operator()(const Event & ev) const {
	if (ev.getDouble5ArrayValue("RECO_MASS", 1) > GetValue())
		return false;
	return true;
}

bool Lep1PtCut::operator()(const Event & ev) const {
	if (ev.getDouble5ArrayValue("RECO_PT", 3) < GetValue())
		return false;
	return true;
}

bool Lep2PtCut::operator()(const Event & ev) const {
	if (ev.getDouble5ArrayValue("RECO_PT", 4) < GetValue())
		return false;
	return true;
}

bool ProjMetCut::operator()(const Event & ev) const {
	double dPhiLep1Met = ev.getDouble5ArrayValue("RECO_DPHIMET", 3);
	double dPhiLep2Met = ev.getDouble5ArrayValue("RECO_DPHIMET", 4);
	double minDelPhi = ( dPhiLep1Met < dPhiLep2Met ) ? dPhiLep1Met : dPhiLep2Met;
	double pfMet = ev.getDoubleVariableValue("RECO_PFMET");
	double projMET = ( minDelPhi > ( M_PI / 2 ) ) ? ( pfMet ) : ( pfMet * std::sin ( minDelPhi ) );
	if (projMET < GetValue())
		return false;
	return true;
}

bool BtagCut::operator()(const Event & ev) const {
	const double * btags = ev.getDouble25ArrayAdr("RECOJET_BTAG_TrackCountingHighEff", 0);
	for (int i = 0; i < 25; ++i) {
		if (*(btags + i) > GetValue())
			return false;
	}
	return true;
}

bool ZptMinCut::operator()(const Event & ev) const {
	if (ev.getDouble5ArrayValue("RECO_PT", 1) < GetValue())
		return false;
	return true;
}

void ZptMinCut::SetValue(double val) {
	if ( val < 0 )
		Cut::SetValue(0);
	else
		Cut::SetValue(val);
}

bool ZptMaxCut::operator()(const Event & ev) const {
	if (ev.getDouble5ArrayValue("RECO_PT", 1) > GetValue())
		return false;
	return true;
}

bool PtMinPtMaxCut::operator()(const Event & ev) const {
	double ptMin = ev.getDouble5ArrayValue("RECO_PT", 4);
	double ptMax = ev.getDouble5ArrayValue("RECO_PT", 3);
	double ratio = ptMin / ptMax;
	if (ratio < GetValue())
		return false;
	return true;
}

void PtMinPtMaxCut::SetValue(double val) {
	if ( val < 0 )
		Cut::SetValue(0);
	else
		Cut::SetValue(val);
}

bool HmassMinCut::operator()(const Event & ev) const {
	if (ev.getDouble5ArrayValue("RECO_MASS", 0) < GetValue())
		return false;
	return true;
}

bool HmassMaxCut::operator()(const Event & ev) const {
	if (ev.getDouble5ArrayValue("RECO_MASS", 0) > GetValue())
		return false;
	return true;
}

bool MetMinCut::operator()(const Event & ev) const {
	if (ev.getDouble5ArrayValue("RECO_PT", 2) < GetValue())
		return false;
	return true;
}

bool MetMaxCut::operator()(const Event & ev) const {
	if (ev.getDouble5ArrayValue("RECO_PT", 2) > GetValue())
		return false;
	return true;
}

bool JetMetPhiCut::operator()(const Event & ev) const {
	const double * angle = ev.getDouble25ArrayAdr("RECOJET_DELTAPHIMET", 0);
	for (int i = 0; i < 25; ++i) {
		double val = *(angle + i);
		if (val < GetValue() && val > -998)
			return false;
	}
	return true;
}

bool MetSignCut::operator()(const Event & ev) const {
	if (ev.getDoubleVariableValue("RECO_PFMET_SIGN") < GetValue())
		return false;
	return true;
}

bool ThetaStarMinCut::operator()(const Event & ev) const {
	double zPt = ev.getDouble5ArrayValue("RECO_PT", 1);
	double zEta = ev.getDouble5ArrayValue("RECO_ETA", 1);
	double zPhi = ev.getDouble5ArrayValue("RECO_PHI", 1);
	double zE = ev.getDouble5ArrayValue("RECO_E", 1);
	TLorentzVector zV;
	zV.SetPtEtaPhiE( zPt, zEta, zPhi, zE);
	
	int posLepIdx = -1;
	if ( ev.getDouble5ArrayValue("RECO_CHARGE", 3) == 1 )
		posLepIdx = 3;
	else
		posLepIdx = 4;

	double lPt = ev.getDouble5ArrayValue("RECO_PT", posLepIdx);
	double lEta = ev.getDouble5ArrayValue("RECO_ETA", posLepIdx);
	double lPhi = ev.getDouble5ArrayValue("RECO_PHI", posLepIdx);
	double lE = ev.getDouble5ArrayValue("RECO_E", posLepIdx);
	TLorentzVector lV;
	lV.SetPtEtaPhiE( lPt, lEta, lPhi, lE);

	TVector3 boostV = - zV.BoostVector();
	lV.Boost(boostV);
	double thetaS = lV.Angle(zV.Vect());
	if ( thetaS < GetValue())
		return false;
	return true;
}

bool ThetaStarMaxCut::operator()(const Event & ev) const {
	double zPt = ev.getDouble5ArrayValue("RECO_PT", 1);
	double zEta = ev.getDouble5ArrayValue("RECO_ETA", 1);
	double zPhi = ev.getDouble5ArrayValue("RECO_PHI", 1);
	double zE = ev.getDouble5ArrayValue("RECO_E", 1);
	TLorentzVector zV;
	zV.SetPtEtaPhiE( zPt, zEta, zPhi, zE);
	
	int posLepIdx = -1;
	if ( ev.getDouble5ArrayValue("RECO_CHARGE", 3) == 1 )
		posLepIdx = 3;
	else
		posLepIdx = 4;

	double lPt = ev.getDouble5ArrayValue("RECO_PT", posLepIdx);
	double lEta = ev.getDouble5ArrayValue("RECO_ETA", posLepIdx);
	double lPhi = ev.getDouble5ArrayValue("RECO_PHI", posLepIdx);
	double lE = ev.getDouble5ArrayValue("RECO_E", posLepIdx);
	TLorentzVector lV;
	lV.SetPtEtaPhiE( lPt, lEta, lPhi, lE);

	TVector3 boostV = - zV.BoostVector();
	lV.Boost(boostV);
	double thetaS = lV.Angle(zV.Vect());
	if ( thetaS > GetValue())
		return false;
	return true;
}

extern JetSelector sel_jet;

bool Met_max_Cut::operator()(const Event & ev) const {
	if (ev.RECO_PT[2] > GetValue())
		return false;
	return true;
}
bool ZMet_acop_Cut::operator()(const Event & ev) const {
	if (ev.RECO_ZMETACOP > GetValue())
		return false;
	return true;
}

bool dPhiMET2Lep_Cut::operator()(const Event & ev) const {
	if (ev.RECO_DPHIMET[4] < GetValue())
		return false;
	return true;
}

bool Isol_Cut::operator()(const Event & ev) const {
	double temp_isolation = 2.0 * ev.RECO_HCALISO[3] + ev.RECO_ECALISO[3] + 2.0 * ev.RECO_HCALISO[4] + ev.RECO_ECALISO[4];
	if (temp_isolation > GetValue())
		return false;
	return true;
}

bool ZMetPt_Cut::operator()(const Event & ev) const {
	double temp =  ev.RECO_PT[1] + ev.RECO_PT[2];
	if (temp > GetValue())
		return false;
	return true;
}

bool DPhiLep_Cut::operator()(const Event & ev) const {
	if (std::fabs(ev.RECO_PHI[3] - ev.RECO_PHI[4]) < GetValue())
		return false;
	return true;
}


bool TrackIso_Cut::operator()(const Event & ev) const {
	if ( (ev.RECO_TRACKISO[3] / ev.RECO_PT[3]) > GetValue() || (ev.RECO_TRACKISO[4] / ev.RECO_PT[4]) > GetValue())
		return false;
	return true;
}

bool Nvtx_Cut::operator()(const Event & ev) const {
	if ( ev.RECO_NVTX != GetValue() )
		return false;
	return true;
}
*/
