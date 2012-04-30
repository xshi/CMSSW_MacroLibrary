// ROOT Libraries
#include <TFile.h>
#include <TLorentzVector.h>
// Standard Libraries
//#include <cmath>
#include <iostream>
#include <iomanip>
//#include <vector>
//#include <string>
//#include <sstream>
// Other
#include "preselectionCMG.h"
#include "event.h"
#include "options.h"
//#include "triggerinfo.h"
#include "toolbox.h"
#include "muon.h"
#include "electron.h"

using std::cout;
using std::cin;
using std::endl;
using std::string;
using std::vector;
using std::stringstream;
using std::setw;

void LeptonPreselectionCMG( const Options & opt, PreselType type ) {
	if (type == ELE)
		cout << "Entering ElectronPreselection() ..." << endl;
	else if (type == MU)
		cout << "Entering MuonPreselection() ..." << endl;

	TFile * file = new TFile( opt.checkStringOption("inpath").c_str() );
	if (!file->IsOpen())
		throw string("ERROR: Can't open the file: " + opt.checkStringOption("inpath") + "!");
	TDirectory * dir = (TDirectory *) file->Get("evAnalyzer");
	TTree * tree = ( TTree * ) dir->Get( "data" );
	Event ev( tree );

	string outputFile( opt.checkStringOption("outpath") );
	size_t pos = outputFile.find(".root");
	if ( pos != string::npos )
		outputFile.erase(pos);

	if (type == ELE)
		outputFile += "_elePresel.root";
	else if (type == MU)
		outputFile += "_muPresel.root";
	
	cout << outputFile << endl;
	TFile * out = new TFile( outputFile.c_str(), "recreate" );

	unsigned run;
	unsigned lumi;
	unsigned event;
	double pfmet;
	int nele;
	int nmu;
	int nsoftmu;
	double l1pt;
	double l1eta;
	double l1phi;
	double l2pt;
	double l2eta;
	double l2phi;
	double zmass;
	double zpt;
	double zeta;
	double mt;
	int njet;
	int nsoftjet;
	double maxJetBTag;
	double minDeltaPhiJetMet;
	double minDeltaPhiSoftJetMet;
	int nvtx;
	int ni;

	TTree * smallTree = new TTree("HZZ2l2nuAnalysis", "HZZ2l2nu Analysis Tree");
	smallTree->Branch( "Run", &run, "Run/i" );
	smallTree->Branch( "Lumi", &lumi, "Lumi/i" );
	smallTree->Branch( "Event", &event, "Event/i" );
	smallTree->Branch( "PFMET", &pfmet, "PFMET/D" );
	smallTree->Branch( "NELE", &nele, "NELE/I" );
	smallTree->Branch( "NMU", &nmu, "NMU/I" );
	smallTree->Branch( "NSOFTMU", &nsoftmu, "NSOFTMU/I" );
	smallTree->Branch( "L1PT", &l1pt, "L1PT/D" );
	smallTree->Branch( "L1ETA", &l1eta, "L1ETA/D" );
	smallTree->Branch( "L1PHI", &l1phi, "L1PHI/D" );
	smallTree->Branch( "L2PT", &l2pt, "L2PT/D" );
	smallTree->Branch( "L2ETA", &l2eta, "L2ETA/D" );
	smallTree->Branch( "L2PHI", &l2phi, "L2PHI/D" );
	smallTree->Branch( "ZMASS", &zmass, "ZMASS/D" );
	smallTree->Branch( "ZPT", &zpt, "ZPT/D" );
	smallTree->Branch( "ZETA", &zeta, "ZETA/D" );
	smallTree->Branch( "MT", &mt, "MT/D" );
	smallTree->Branch( "NJET", &njet, "NJET/I" );
	smallTree->Branch( "NSOFTJET", &nsoftjet, "NSOFTJET/I" );
	smallTree->Branch( "MAXJETBTAG", &maxJetBTag, "MAXJETBTAG/D" );
	smallTree->Branch( "MINDPJETMET", &minDeltaPhiJetMet, "MINDPJETMET/D" );
	smallTree->Branch( "MINDPSOFTJETMET", &minDeltaPhiSoftJetMet, "MINDPSOFTJETMET/D" );
	smallTree->Branch( "NVTX", &nvtx, "NVTX/I" );
	smallTree->Branch( "nInter" , &ni, "nInter/I" );

	// bool isData = opt.checkBoolOption("isData");

	// unsigned long nentries = tree->GetEntries();
	unsigned long nentries = 100000;
	for ( unsigned long i = 0; i < nentries; i++ ) {
		if ( i % 10000 == 0) {
			cout << string(40, '\b');
			cout << setw(10) << i << " / " << setw(10) << nentries << " done ..." << std::flush;
		}

		tree->GetEntry( i );

		run = -999;
		lumi = -999;
		event = -999;
		pfmet = -999;
		nele = -999;
		nmu = -999;
		nsoftmu = -999;
		l1pt = -999;
		l1eta = -999;
		l1phi = -999;
		l2pt = -999;
		l2eta = -999;
		l2phi = -999;
		zmass = -999;
		zpt = -999;
		zeta = -999;
		mt = -999;
		njet = -999;
		nsoftjet = -999;
		maxJetBTag = -999;
		minDeltaPhiJetMet = -999;
		minDeltaPhiSoftJetMet = -999;
		nvtx = -999;
		ni = -999;

		run = ev.getSVV<int>("run");
		lumi = ev.getSVV<int>("lumi");
		event = ev.getSVV<int>("event");

//		cout << run << ":" << lumi << ":" << event << endl;
//		if (isData) {
//			if (!triggerAccept(ev, type))
//				continue;
//		}

		pfmet = ev.getAVV<float>("met_pt", 0);

		vector<Electron> electrons = buildElectronCollection(ev);
		vector<Muon> muons = buildMuonCollection(ev);

		float rho = ev.getSVV<float>("rho25");

		vector<Electron> selectedElectrons;
		for (unsigned i = 0; i < electrons.size(); ++i) {
			TLorentzVector lv = electrons[i].lorentzVector();
			if ( lv.Pt() > 20 && fabs(lv.Eta()) < 2.5 && !electrons[i].isInCrack() && electrons[i].passesTightID()
					&& electrons[i].isPFIsolatedTight(rho) )
				selectedElectrons.push_back(electrons[i]);
		}

		vector<Electron> looseElectrons;
		for (unsigned i = 0; i < electrons.size(); ++i) {
			TLorentzVector lv = electrons[i].lorentzVector();
			if ( lv.Pt() > 10 && fabs(lv.Eta()) < 2.5 && !electrons[i].isInCrack() && electrons[i].passesTightID()
					&& electrons[i].isPFIsolatedTight(rho) )
				looseElectrons.push_back(electrons[i]);
		}

		vector<Muon> selectedMuons;
		for (unsigned i = 0; i < muons.size(); ++i) {
			TLorentzVector lv = muons[i].lorentzVector();
			if ( lv.Pt() > 20 && fabs(lv.Eta()) < 2.4 && muons[i].isTightMuon()
					&& muons[i].isPFIsolatedTight(Electron::effAreaMC(lv.Eta()), rho) ) {
				selectedMuons.push_back(muons[i]);
			}
		}

		vector<Muon> looseMuons;
		for (unsigned i = 0; i < muons.size(); ++i) {
			TLorentzVector lv = muons[i].lorentzVector();
			if ( lv.Pt() > 10 && fabs(lv.Eta()) < 2.4 && muons[i].isTightMuon()
					&& muons[i].isPFIsolatedTight(Electron::effAreaMC(lv.Eta()), rho) )
				looseMuons.push_back(muons[i]);
		}

		int nLeptons = looseElectrons.size() + looseMuons.size();
		if ( nLeptons > 2 )
			continue;

		string leptonsType;
		Lepton * selectedLeptons[2] = {0};
		if (type == ELE) {
			if (selectedElectrons.size() < 2)
				continue;
		   else {
			   selectedLeptons[0] = &selectedElectrons[0];
			   selectedLeptons[1] = &selectedElectrons[1];
		   }
		} else if (type == MU) {
			if (selectedMuons.size() < 2)
				continue;
			else {
			   selectedLeptons[0] = &selectedMuons[0];
			   selectedLeptons[1] = &selectedMuons[1];
			}
		}
		
//		cout << selectedLeptons[0]->id << endl;
//		cout << selectedLeptons[1]->id << endl;

//		if ( selectedLeptons[0]->id != -selectedLeptons[1]->id )
//			continue;

//		cout << "# loose electrons : " << looseElectrons.size() << endl;
//		cout << "# muons : " << selectedMuons.size() << endl;
//		cout << "# loose muons : " << looseMuons.size() << endl;

		nele = looseElectrons.size();
		nmu = looseMuons.size();
		nsoftmu = 0;

		TLorentzVector lep1 = selectedLeptons[0]->lorentzVector();
		TLorentzVector lep2 = selectedLeptons[1]->lorentzVector();

		if (lep1.Pt() < 20 || lep2.Pt() < 20)
			continue;

		if (lep2.Pt() > lep1.Pt()) {
			TLorentzVector temp = lep1;
			lep1 = lep2;
			lep2 = temp;
		}

		l1pt = lep1.Pt();
		l1eta = lep1.Eta();
		l1phi = lep1.Phi();

		l2pt = lep2.Pt();
		l2eta = lep2.Eta();
		l2phi = lep2.Phi();

		TLorentzVector Zcand = lep1 + lep2;
		zpt = Zcand.Pt();
		zeta = Zcand.Eta();
		zmass = Zcand.M();
		if (zmass < 76.1876 || zmass > 106.1876 || zpt < 55)
			continue;

		double metphi = ev.getAVV<float>("met_phi", 0);
		double metPx = pfmet * cos( metphi );
		double metPy = pfmet * sin( metphi );
		double px = metPx + Zcand.Px();
		double py = metPy + Zcand.Py();
		double pt2 = px * px + py * py;
		double e = sqrt(zpt * zpt + zmass * zmass) + sqrt(pfmet * pfmet + zmass * zmass);
		double mt2 = e * e - pt2;
		mt = (mt2 > 0) ? sqrt(mt2) : 0;

		vector<unsigned> jets;
		selectJetsCMG( ev, jets );
		njet = jets.size();

		vector<unsigned> softjets;
		selectJetsCMG( ev, softjets, 15 );
		nsoftjet = softjets.size();

		const float * jnPx = ev.getAVA<float>("jn_px");
		const float * jnPy = ev.getAVA<float>("jn_py");
		const float * jnPz = ev.getAVA<float>("jn_pz");
		const float * jnEn = ev.getAVA<float>("jn_en");
		const float * btag = ev.getAVA<float>("jn_btag1");

		maxJetBTag = -999;
		for ( int j = 0; j < njet; ++j ) {
			unsigned jetIdx = jets[j];
			TLorentzVector jet( jnPx[jetIdx], jnPy[jetIdx], jnPz[jetIdx], jnEn[jetIdx] );
			if ( btag[jetIdx] > maxJetBTag && fabs(jet.Eta()) < 2.4 )
				maxJetBTag = btag[jetIdx];
		}
		/*
		unsigned temprun = 160956;
		unsigned templumi = 53;
		unsigned tempevent = 27115287;
		if (run == temprun && lumi == templumi  && event == tempevent) {
			cout << "Found!" << endl;
			cout << "Trigger: " << triggerAccept(ev, type) << endl;
			cout << "Nlep20: " << leptons20.size() << endl;
			cout << "Nlep: " << nmu + nele << endl;
			cout << "Nsoftmu: " << nsoftmu << endl;
			cout << "Same charge: " << (lepCharge[ leptons20[0] ] == lepCharge[ leptons20[1] ]) << endl;
			cout << "zmass: " << zmass << endl;
			cout << "zpt: " << zpt << endl;
			cout << "maxJetBTag: " << maxJetBTag << endl;
			cin.get();
		}
		*/
		if (maxJetBTag > 2.0)
			continue;

		minDeltaPhiJetMet = 999;
		for ( int j = 0; j < njet; ++j ) {
			unsigned jetIdx = jets[j];
			TLorentzVector jet( jnPx[jetIdx], jnPy[jetIdx], jnPz[jetIdx], jnEn[jetIdx] );
			double tempDelPhiJetMet = deltaPhi(metphi, jet.Phi());
			if ( tempDelPhiJetMet < minDeltaPhiJetMet )
				minDeltaPhiJetMet = tempDelPhiJetMet;
		}
		minDeltaPhiSoftJetMet = 999;
		for ( int j = 0; j < nsoftjet; ++j ) {
			unsigned jetIdx = softjets[j];
			TLorentzVector jet( jnPx[jetIdx], jnPy[jetIdx], jnPz[jetIdx], jnEn[jetIdx] );
			double tempDelPhiSoftJetMet = deltaPhi(metphi, jet.Phi());
			if ( tempDelPhiSoftJetMet < minDeltaPhiSoftJetMet )
				minDeltaPhiSoftJetMet = tempDelPhiSoftJetMet;
		}

		nvtx = ev.getSingleVariableValue<int>("nvtx");

		ni = ev.getSVV<int>("ngenITpu");
		
		smallTree->Fill();
	}
	cout << endl;
	delete file;
	smallTree->Write("", TObject::kOverwrite);
	delete smallTree;
	delete out;
}

vector<Muon> buildMuonCollection( const Event & ev) {
	vector<Muon> muons;
	const int * m_idbits	= ev.getAVA<int>("mn_idbits");
	const float * m_nMatches	= ev.getAVA<float>("mn_nMatches");
	const float * m_validMuonHits	= ev.getAVA<float>("mn_validMuonHits");
	
	if (fabs(ev.getSVV<int>("l1_id")) == 13) {
		float px = ev.getSVV<float>("l1_px");
		float py = ev.getSVV<float>("l1_py");
		float pz = ev.getSVV<float>("l1_pz");
		float en = ev.getSVV<float>("l1_en");
		float ptErr = ev.getSVV<float>("l1_ptErr");
		float ecalIso = ev.getSVV<float>("l1_ecalIso");
		float hcalIso = ev.getSVV<float>("l1_hcalIso");
		float trkIso = ev.getSVV<float>("l1_trkIso");
		float gIso = ev.getSVV<float>("l1_gIso");
		float chIso = ev.getSVV<float>("l1_chIso");
		float puchIso = ev.getSVV<float>("l1_puchIso");
		float nhIso = ev.getSVV<float>("l1_nhIso");
		int id = ev.getSVV<int>("l1_id");
		int genid = ev.getSVV<int>("l1_genid");
		float ensf = ev.getSVV<float>("l1_ensf");
		float ensferr = ev.getSVV<float>("l1_ensferr");
		float d0 = ev.getSVV<float>("l1_d0");
		float dZ = ev.getSVV<float>("l1_dZ");
		float trkpt = ev.getSVV<float>("l1_trkpt");
		float trketa = ev.getSVV<float>("l1_trketa");
		float trkphi = ev.getSVV<float>("l1_trkphi");
		float trkchi2 = ev.getSVV<float>("l1_trkchi2");
		float trkValidPixelHits = ev.getSVV<float>("l1_trkValidPixelHits");
		float trkValidTrackerHits = ev.getSVV<float>("l1_trkValidTrackerHits");
		float trkLostInnerHits = ev.getSVV<float>("l1_trkLostInnerHits");
		int pid = ev.getSVV<int>("l1_pid");
		int idbits = m_idbits[pid];
		int nMatches = m_nMatches[pid];
		int validMuonHits = m_validMuonHits[pid];

		Muon tmp(px, py, pz, en, ptErr, ecalIso, hcalIso, trkIso, gIso, chIso, puchIso, nhIso, id, genid,
				ensf, ensferr, d0, dZ, trkpt, trketa, trkphi, trkchi2, trkValidPixelHits,
				trkValidTrackerHits, trkLostInnerHits, idbits, nMatches, validMuonHits);
		muons.push_back(tmp);
	}
	if (fabs(ev.getSVV<int>("l2_id")) == 13) {
		float px = ev.getSVV<float>("l2_px");
		float py = ev.getSVV<float>("l2_py");
		float pz = ev.getSVV<float>("l2_pz");
		float en = ev.getSVV<float>("l2_en");
		float ptErr = ev.getSVV<float>("l2_ptErr");
		float ecalIso = ev.getSVV<float>("l2_ecalIso");
		float hcalIso = ev.getSVV<float>("l2_hcalIso");
		float trkIso = ev.getSVV<float>("l2_trkIso");
		float gIso = ev.getSVV<float>("l2_gIso");
		float chIso = ev.getSVV<float>("l2_chIso");
		float puchIso = ev.getSVV<float>("l2_puchIso");
		float nhIso = ev.getSVV<float>("l2_nhIso");
		int id = ev.getSVV<int>("l2_id");
		int genid = ev.getSVV<int>("l2_genid");
		float ensf = ev.getSVV<float>("l2_ensf");
		float ensferr = ev.getSVV<float>("l2_ensferr");
		float d0 = ev.getSVV<float>("l2_d0");
		float dZ = ev.getSVV<float>("l2_dZ");
		float trkpt = ev.getSVV<float>("l2_trkpt");
		float trketa = ev.getSVV<float>("l2_trketa");
		float trkphi = ev.getSVV<float>("l2_trkphi");
		float trkchi2 = ev.getSVV<float>("l2_trkchi2");
		float trkValidPixelHits = ev.getSVV<float>("l2_trkValidPixelHits");
		float trkValidTrackerHits = ev.getSVV<float>("l2_trkValidTrackerHits");
		float trkLostInnerHits = ev.getSVV<float>("l2_trkLostInnerHits");
		int pid = ev.getSVV<int>("l2_pid");
		int idbits = m_idbits[pid];
		int nMatches = m_nMatches[pid];
		int validMuonHits = m_validMuonHits[pid];

		Muon tmp(px, py, pz, en, ptErr, ecalIso, hcalIso, trkIso, gIso, chIso, puchIso, nhIso, id, genid,
				ensf, ensferr, d0, dZ, trkpt, trketa, trkphi, trkchi2, trkValidPixelHits,
				trkValidTrackerHits, trkLostInnerHits, idbits, nMatches, validMuonHits);
		muons.push_back(tmp);
	}

	int ln = ev.getSVV<int>("ln");
	const float * px = ev.getAVA<float>("ln_px");
	const float * py = ev.getAVA<float>("ln_py");
	const float * pz = ev.getAVA<float>("ln_pz");
	const float * en = ev.getAVA<float>("ln_en");
	const float * ptErr = ev.getAVA<float>("ln_ptErr");
	const float * ecalIso = ev.getAVA<float>("ln_ecalIso");
	const float * hcalIso = ev.getAVA<float>("ln_hcalIso");
	const float * trkIso = ev.getAVA<float>("ln_trkIso");
	const float * gIso = ev.getAVA<float>("ln_gIso");
	const float * chIso = ev.getAVA<float>("ln_chIso");
	const float * puchIso = ev.getAVA<float>("ln_puchIso");
	const float * nhIso = ev.getAVA<float>("ln_nhIso");
	const int * id = ev.getAVA<int>("ln_id");
	const int * genid = ev.getAVA<int>("ln_genid");
	const float * ensf = ev.getAVA<float>("ln_ensf");
	const float * ensferr = ev.getAVA<float>("ln_ensferr");
	const float * d0 = ev.getAVA<float>("ln_d0");
	const float * dZ = ev.getAVA<float>("ln_dZ");
	const float * trkpt = ev.getAVA<float>("ln_trkpt");
	const float * trketa = ev.getAVA<float>("ln_trketa");
	const float * trkphi = ev.getAVA<float>("ln_trkphi");
	const float * trkchi2 = ev.getAVA<float>("ln_trkchi2");
	const float * trkValidPixelHits = ev.getAVA<float>("ln_trkValidPixelHits");
	const float * trkValidTrackerHits = ev.getAVA<float>("ln_trkValidTrackerHits");
	const float * trkLostInnerHits = ev.getAVA<float>("ln_trkLostInnerHits");
	const int * pid = ev.getAVA<int>("ln_pid");

	for (int i = 0; i < ln; ++i) {
		if (fabs(id[i]) == 13) {
			Muon tmp(px[i], py[i], pz[i], en[i], ptErr[i], ecalIso[i], hcalIso[i], trkIso[i], gIso[i],
					chIso[i], puchIso[i], nhIso[i], id[i], genid[i], ensf[i], ensferr[i], d0[i], dZ[i],
					trkpt[i], trketa[i], trkphi[i], trkchi2[i], trkValidPixelHits[i],
					trkValidTrackerHits[i], trkLostInnerHits[i], m_idbits[pid[i]], m_nMatches[pid[i]],
					m_validMuonHits[pid[i]]);
			muons.push_back(tmp);
		}
	}
	return muons;
}

vector<Electron> buildElectronCollection(const Event & ev) {
	vector<Electron> electrons;
	const int * e_idbits	= ev.getAVA<int>("en_idbits");
	const float * e_hoe	= ev.getAVA<float>("en_hoe");
	const float * e_dphiin = ev.getAVA<float>("en_dphiin");
	const float * e_detain = ev.getAVA<float>("en_detain");
	const float * e_sihih = ev.getAVA<float>("en_sihih");
	const float * e_sipip = ev.getAVA<float>("en_sipip");
	const float * e_r9 = ev.getAVA<float>("en_r9");
	const float * e_sce = ev.getAVA<float>("en_sce");
	const float * e_sceta = ev.getAVA<float>("en_sceta");
	const float * e_scphi = ev.getAVA<float>("en_scphi");
	const float * e_e2x5max = ev.getAVA<float>("en_e2x5max");
	const float * e_e1x5 = ev.getAVA<float>("en_e1x5");
	const float * e_e5x5 = ev.getAVA<float>("en_e5x5");
	const float * e_h2te = ev.getAVA<float>("en_h2te");
	const float * e_h2tebc = ev.getAVA<float>("en_h2tebc");
	const float * e_ooemoop = ev.getAVA<float>("en_ooemoop");
	const float * e_fbrem = ev.getAVA<float>("en_fbrem");
	const float * e_eopin = ev.getAVA<float>("en_eopin");
	
	if (fabs(ev.getSVV<int>("l1_id")) == 11) {
		float px = ev.getSVV<float>("l1_px");
		float py = ev.getSVV<float>("l1_py");
		float pz = ev.getSVV<float>("l1_pz");
		float en = ev.getSVV<float>("l1_en");
		float ptErr = ev.getSVV<float>("l1_ptErr");
		float ecalIso = ev.getSVV<float>("l1_ecalIso");
		float hcalIso = ev.getSVV<float>("l1_hcalIso");
		float trkIso = ev.getSVV<float>("l1_trkIso");
		float gIso = ev.getSVV<float>("l1_gIso");
		float chIso = ev.getSVV<float>("l1_chIso");
		float puchIso = ev.getSVV<float>("l1_puchIso");
		float nhIso = ev.getSVV<float>("l1_nhIso");
		int id = ev.getSVV<int>("l1_id");
		int genid = ev.getSVV<int>("l1_genid");
		float ensf = ev.getSVV<float>("l1_ensf");
		float ensferr = ev.getSVV<float>("l1_ensferr");
		float d0 = ev.getSVV<float>("l1_d0");
		float dZ = ev.getSVV<float>("l1_dZ");
		float trkpt = ev.getSVV<float>("l1_trkpt");
		float trketa = ev.getSVV<float>("l1_trketa");
		float trkphi = ev.getSVV<float>("l1_trkphi");
		float trkchi2 = ev.getSVV<float>("l1_trkchi2");
		float trkValidPixelHits = ev.getSVV<float>("l1_trkValidPixelHits");
		float trkValidTrackerHits = ev.getSVV<float>("l1_trkValidTrackerHits");
		float trkLostInnerHits = ev.getSVV<float>("l1_trkLostInnerHits");
		int pid = ev.getSVV<int>("l1_pid");

		Electron tmp(px, py, pz, en, ptErr, ecalIso, hcalIso, trkIso, gIso, chIso, puchIso, nhIso, id, genid,
				ensf, ensferr, d0, dZ, trkpt, trketa, trkphi, trkchi2, trkValidPixelHits,
				trkValidTrackerHits, trkLostInnerHits, e_idbits[pid], e_hoe[pid], e_dphiin[pid], e_detain[pid],
				e_sihih[pid], e_sipip[pid], e_r9[pid], e_sce[pid], e_sceta[pid], e_scphi[pid], e_e2x5max[pid],
				e_e1x5[pid], e_e5x5[pid], e_h2te[pid], e_h2tebc[pid], e_ooemoop[pid], e_fbrem[pid], e_eopin[pid]);
		electrons.push_back(tmp);
	}
	if (fabs(ev.getSVV<int>("l2_id")) == 11) {
		float px = ev.getSVV<float>("l2_px");
		float py = ev.getSVV<float>("l2_py");
		float pz = ev.getSVV<float>("l2_pz");
		float en = ev.getSVV<float>("l2_en");
		float ptErr = ev.getSVV<float>("l2_ptErr");
		float ecalIso = ev.getSVV<float>("l2_ecalIso");
		float hcalIso = ev.getSVV<float>("l2_hcalIso");
		float trkIso = ev.getSVV<float>("l2_trkIso");
		float gIso = ev.getSVV<float>("l2_gIso");
		float chIso = ev.getSVV<float>("l2_chIso");
		float puchIso = ev.getSVV<float>("l2_puchIso");
		float nhIso = ev.getSVV<float>("l2_nhIso");
		int id = ev.getSVV<int>("l2_id");
		int genid = ev.getSVV<int>("l2_genid");
		float ensf = ev.getSVV<float>("l2_ensf");
		float ensferr = ev.getSVV<float>("l2_ensferr");
		float d0 = ev.getSVV<float>("l2_d0");
		float dZ = ev.getSVV<float>("l2_dZ");
		float trkpt = ev.getSVV<float>("l2_trkpt");
		float trketa = ev.getSVV<float>("l2_trketa");
		float trkphi = ev.getSVV<float>("l2_trkphi");
		float trkchi2 = ev.getSVV<float>("l2_trkchi2");
		float trkValidPixelHits = ev.getSVV<float>("l2_trkValidPixelHits");
		float trkValidTrackerHits = ev.getSVV<float>("l2_trkValidTrackerHits");
		float trkLostInnerHits = ev.getSVV<float>("l2_trkLostInnerHits");
		int pid = ev.getSVV<int>("l2_pid");

		Electron tmp(px, py, pz, en, ptErr, ecalIso, hcalIso, trkIso, gIso, chIso, puchIso, nhIso, id, genid,
				ensf, ensferr, d0, dZ, trkpt, trketa, trkphi, trkchi2, trkValidPixelHits,
				trkValidTrackerHits, trkLostInnerHits, e_idbits[pid], e_hoe[pid], e_dphiin[pid], e_detain[pid],
				e_sihih[pid], e_sipip[pid], e_r9[pid], e_sce[pid], e_sceta[pid], e_scphi[pid], e_e2x5max[pid],
				e_e1x5[pid], e_e5x5[pid], e_h2te[pid], e_h2tebc[pid], e_ooemoop[pid], e_fbrem[pid], e_eopin[pid]);
		electrons.push_back(tmp);
	}

	int ln = ev.getSVV<int>("ln");
	const float * px = ev.getAVA<float>("ln_px");
	const float * py = ev.getAVA<float>("ln_py");
	const float * pz = ev.getAVA<float>("ln_pz");
	const float * en = ev.getAVA<float>("ln_en");
	const float * ptErr = ev.getAVA<float>("ln_ptErr");
	const float * ecalIso = ev.getAVA<float>("ln_ecalIso");
	const float * hcalIso = ev.getAVA<float>("ln_hcalIso");
	const float * trkIso = ev.getAVA<float>("ln_trkIso");
	const float * gIso = ev.getAVA<float>("ln_gIso");
	const float * chIso = ev.getAVA<float>("ln_chIso");
	const float * puchIso = ev.getAVA<float>("ln_puchIso");
	const float * nhIso = ev.getAVA<float>("ln_nhIso");
	const int * id = ev.getAVA<int>("ln_id");
	const int * genid = ev.getAVA<int>("ln_genid");
	const float * ensf = ev.getAVA<float>("ln_ensf");
	const float * ensferr = ev.getAVA<float>("ln_ensferr");
	const float * d0 = ev.getAVA<float>("ln_d0");
	const float * dZ = ev.getAVA<float>("ln_dZ");
	const float * trkpt = ev.getAVA<float>("ln_trkpt");
	const float * trketa = ev.getAVA<float>("ln_trketa");
	const float * trkphi = ev.getAVA<float>("ln_trkphi");
	const float * trkchi2 = ev.getAVA<float>("ln_trkchi2");
	const float * trkValidPixelHits = ev.getAVA<float>("ln_trkValidPixelHits");
	const float * trkValidTrackerHits = ev.getAVA<float>("ln_trkValidTrackerHits");
	const float * trkLostInnerHits = ev.getAVA<float>("ln_trkLostInnerHits");
	const int * pid = ev.getAVA<int>("ln_pid");

	for (int i = 0; i < ln; ++i) {
		if (fabs(id[i]) == 11) {
			Electron tmp(px[i], py[i], pz[i], en[i], ptErr[i], ecalIso[i], hcalIso[i], trkIso[i], gIso[i],
					chIso[i], puchIso[i], nhIso[i], id[i], genid[i], ensf[i], ensferr[i], d0[i], dZ[i],
					trkpt[i], trketa[i], trkphi[i], trkchi2[i], trkValidPixelHits[i],
					trkValidTrackerHits[i], trkLostInnerHits[i], e_idbits[pid[i]], e_hoe[pid[i]],
					e_dphiin[pid[i]], e_detain[pid[i]],	e_sihih[pid[i]], e_sipip[pid[i]], e_r9[pid[i]],
					e_sce[pid[i]], e_sceta[pid[i]], e_scphi[pid[i]], e_e2x5max[pid[i]], e_e1x5[pid[i]],
					e_e5x5[pid[i]], e_h2te[pid[i]], e_h2tebc[pid[i]], e_ooemoop[pid[i]], e_fbrem[pid[i]],
					e_eopin[pid[i]]);
			electrons.push_back(tmp);
		}
	}
	return electrons;
}

void selectJetsCMG(const Event & ev, vector<unsigned> & jets, double ptMin, double etaMax) {
	const int nJets = ev.getSVV<int>("jn");
	const float * jnPx = ev.getAVA<float>("jn_px");
	const float * jnPy = ev.getAVA<float>("jn_py");
	const float * jnPz = ev.getAVA<float>("jn_pz");
	const float * jnEn = ev.getAVA<float>("jn_en");

	jets.clear();
	for ( int i = 0; i < nJets; ++i ) {
		TLorentzVector jet(jnPx[i], jnPy[i], jnPz[i], jnEn[i]);
		if ( jet.Pt() > ptMin && fabs(jet.Eta()) < etaMax )
			jets.push_back( i );
	}
}

/*
void selectPhotons(const Event & ev, vector<unsigned> & photons) {
	const vector<float> & pt = *ev.getVectorFloatAdr("Photons_PT");
	const vector<float> & eta = *ev.getVectorFloatAdr("Photons_superClusterETA");
	const vector<int> & pixSeed = *ev.getVectorIntAdr("Photons_hasPixelSeed");
	const vector<float> & sigIeIe = *ev.getVectorFloatAdr("Photons_sigmaIetaIeta");
	const vector<float> & hOe = *ev.getVectorFloatAdr("Photons_hadronicOverEm");
	const vector<float> & trkIso = *ev.getVectorFloatAdr("Photons_TRACKISO");
	const vector<float> & ecalIso = *ev.getVectorFloatAdr("Photons_ECALISO");
	const vector<float> & hcalIso = *ev.getVectorFloatAdr("Photons_HCALISO");

	photons.clear();
	const unsigned nPhot = pt.size();
	for ( unsigned i = 0; i < nPhot; ++i ) {
		if ( pt[i] > 55 &&
				fabs(eta[i]) < 1.4442 &&
				pixSeed[i] == 0 &&
				sigIeIe[i] > 0.001 &&
				sigIeIe[i] < 0.013 &&
				hOe[i] < 0.05 &&
				trkIso[i] < (2.0 + 0.001 * pt[i]) &&
				ecalIso[i] < (4.2 + 0.006 * pt[i]) &&
				hcalIso[i] < (2.2 + 0.0025 * pt[i])
		) {
			photons.push_back( i );
		}
	}
	const vector<float> & pt = *ev.getVectorFloatAdr("Photons_PT");
	const vector<float> & eta = *ev.getVectorFloatAdr("Photons_superClusterETA");
	const vector<float> & r9 = *ev.getVectorFloatAdr("Photons_R9");
	const vector<float> & IsoSumOverEt = *ev.getVectorFloatAdr("Photons_IsoSumOverEt");
	const vector<float> & IsoSumOverEtWorst = *ev.getVectorFloatAdr("Photons_IsoSumOverEtWorst");
	const vector<float> & TrkIsoOverEt = *ev.getVectorFloatAdr("Photons_TrkIsoOverEt");
	const vector<float> & sigIeIe = *ev.getVectorFloatAdr("Photons_sigmaIetaIeta");
	const vector<float> & hOe = *ev.getVectorFloatAdr("Photons_hadronicOverEm");
	const vector<float> & delR = *ev.getVectorFloatAdr("Photons_dREleTrack");
	const vector<float> & nLostHits = *ev.getVectorFloatAdr("Photons_nLostHitsEleTrack");

	photons.clear();
	const unsigned nPhot = pt.size();
	for ( unsigned i = 0; i < nPhot; ++i ) {
		bool cat1 = (fabs(eta[i]) < 1.4442) && (r9[i] > 0.94);
		bool cat2 = (fabs(eta[i]) < 1.4442) && (r9[i] <= 0.94);
		bool cat3 = (fabs(eta[i]) > 1.566 && fabs(eta[i]) < 2.5) && (r9[i] > 0.94);
		bool cat4 = (fabs(eta[i]) > 1.566 && fabs(eta[i]) < 2.5) && (r9[i] <= 0.94);
		
		bool passID = false;
		if ( pt[i] > 55 ) {
			if ( cat1 ) {
				if ( IsoSumOverEt[i] < 3.8 &&
						IsoSumOverEtWorst[i] < 11.7 &&
						TrkIsoOverEt[i] < 3.5 &&
						sigIeIe[i] < 0.0106 &&
						hOe[i] < 0.082 &&
						r9[i] > 0.94 )
					passID = true;
			}
			if ( cat2 ) {
				if ( IsoSumOverEt[i] < 2.2 &&
						IsoSumOverEtWorst[i] < 3.4 &&
						TrkIsoOverEt[i] < 2.2 &&
						sigIeIe[i] < 0.0097 &&
						hOe[i] < 0.062 &&
						r9[i] > 0.36
						// && delR[i] < 0.062
						)
					passID = true;
			}
			if ( cat3 ) {
				if ( IsoSumOverEt[i] < 1.77 &&
						IsoSumOverEtWorst[i] < 3.9 &&
						TrkIsoOverEt[i] < 2.3 &&
						sigIeIe[i] < 0.028 &&
						hOe[i] < 0.065 &&
						r9[i] > 0.94 )
					passID = true;
			}
			if ( cat4 ) {
				if ( IsoSumOverEt[i] < 1.29 &&
						IsoSumOverEtWorst[i] < 1.84 &&
						TrkIsoOverEt[i] < 1.45 &&
						sigIeIe[i] < 0.027 &&
						hOe[i] < 0.048 &&
						r9[i] > 0.32 )
					passID = true;
			}
		}
		bool passEleVeto = true;
		if ( delR[i] < 0.5 && nLostHits[i] == 0 )
			passEleVeto = false;
		if (passID && passEleVeto)
			photons.push_back( i );
	}
}

bool triggerAccept( const Event & ev, double pt, double & weight ) {
	bool hlt50 = (ev.getSingleVariableValue<int>("HLT_Photon50_CaloIdVL_accept") == 1);
	bool hlt50i = (ev.getSingleVariableValue<int>("HLT_Photon50_CaloIdVL_IsoL_accept") == 1);
	bool hlt75 = (ev.getSingleVariableValue<int>("HLT_Photon75_CaloIdVL_accept") == 1);
	bool hlt75i = (ev.getSingleVariableValue<int>("HLT_Photon75_CaloIdVL_IsoL_accept") == 1);
	bool hlt90 = (ev.getSingleVariableValue<int>("HLT_Photon90_CaloIdVL_accept") == 1);
	bool hlt90i = (ev.getSingleVariableValue<int>("HLT_Photon90_CaloIdVL_IsoL_accept") == 1);
	int prescale50 = ev.getSingleVariableValue<int>("HLT_Photon50_CaloIdVL_prescale");
	int prescale50i = ev.getSingleVariableValue<int>("HLT_Photon50_CaloIdVL_IsoL_prescale");
	int prescale75 = ev.getSingleVariableValue<int>("HLT_Photon75_CaloIdVL_prescale");
	int prescale75i = ev.getSingleVariableValue<int>("HLT_Photon75_CaloIdVL_IsoL_prescale");
	int prescale90 = ev.getSingleVariableValue<int>("HLT_Photon90_CaloIdVL_prescale");
	int prescale90i = ev.getSingleVariableValue<int>("HLT_Photon90_CaloIdVL_IsoL_prescale");

	if ( pt >= 55 && pt < 80 ) {
		if ( hlt50 || hlt50i ) {
			vector<int> prescales;
			if (hlt50)
				prescales.push_back( prescale50 );
			if (hlt50i)
				prescales.push_back( prescale50i );
			weight = min<int>( prescales );
			return true;
		}
	} else if ( pt >= 80 && pt < 95 ) {
		if ( hlt50 || hlt50i || hlt75 || hlt75i ) {
			vector<int> prescales;
			if (hlt50)
				prescales.push_back( prescale50 );
			if (hlt50i)
				prescales.push_back( prescale50i );
			if (hlt75)
				prescales.push_back( prescale75 );
			if (hlt75i)
				prescales.push_back( prescale75i );
			weight = min<int>( prescales );
			return true;
		}
	} else if ( pt >= 95 ) {
		if ( hlt50 || hlt50i || hlt75 || hlt75i || hlt90 || hlt90i ) {
			vector<int> prescales;
			if (hlt50)
				prescales.push_back( prescale50 );
			if (hlt50i)
				prescales.push_back( prescale50i );
			if (hlt75)
				prescales.push_back( prescale75 );
			if (hlt75i)
				prescales.push_back( prescale75i );
			if (hlt90)
				prescales.push_back( prescale90 );
			if (hlt90i)
				prescales.push_back( prescale90i );
			weight = min<int>( prescales );
			return true;
		}
	} else {
		cout << "ERROR: Unknown preselection type!" << endl;
		exit( EXIT_FAILURE );
	}
	return false;
}

bool triggerAccept( const Event & ev, PreselType type ) {
	int run = ev.getSingleVariableValue<unsigned>("Run");
	if ( type == ELE ) {
		if ( run >= 160329 && run <= 170064 ) {
			if ( ev.getSingleVariableValue<int>("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_prescale") != 1 ) {
				cout << "ERROR: Prescale different from 1!" << endl;
				exit( EXIT_FAILURE );
			}
			bool hltEle17cEle8c = (ev.getSingleVariableValue<int>("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_accept") == 1);
			if ( hltEle17cEle8c )
				return true;
		} else if ( run >= 170065 && run <= 180296 ) {
			if ( ev.getSingleVariableValue<int>("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_prescale") != 1 ) {
				cout << "ERROR: Prescale different from 1!" << endl;
				exit( EXIT_FAILURE );
			}
			bool hltEle17ctEle8ct = (ev.getSingleVariableValue<int>("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_accept") == 1);
			if ( hltEle17ctEle8ct )
				return true;
		} else {
			cout << "ERROR: Unknown run: " << run << "!" << endl;
			exit( EXIT_FAILURE );
		}
	} else if ( type == MU ) {
		if ( run >= 160329 && run <= 163869 ) {
			if ( ev.getSingleVariableValue<int>("HLT_DoubleMu7_prescale") != 1 ) {
				cout << "ERROR: Prescale different from 1!" << endl;
				exit( EXIT_FAILURE );
			}
			bool hltDoubleMu7 = (ev.getSingleVariableValue<int>("HLT_DoubleMu7_accept") == 1);
			if ( hltDoubleMu7 )
				return true;
		} else if ( run >= 165071 && run <= 178410 ) {
			if ( ev.getSingleVariableValue<int>("HLT_Mu13_Mu8_prescale") != 1 ) {
				cout << "ERROR: Prescale different from 1!" << endl;
				exit( EXIT_FAILURE );
			}
			bool hltMu13Mu8 = (ev.getSingleVariableValue<int>("HLT_Mu13_Mu8_accept") == 1);
			if ( hltMu13Mu8 )
				return true;
		} else if ( run >= 178411 && run <= 180296 ) {
			if ( ev.getSingleVariableValue<int>("HLT_Mu17_Mu8_prescale") != 1 ) {
				cout << "ERROR: Prescale different from 1!" << endl;
				exit( EXIT_FAILURE );
			}
			bool hltMu17Mu8 = (ev.getSingleVariableValue<int>("HLT_Mu17_Mu8_accept") == 1);
			if ( hltMu17Mu8 )
				return true;
		} else {
			cout << "ERROR: Unknown run: " << run << "!" << endl;
			exit( EXIT_FAILURE );
		}
	} else {
		cout << "ERROR: Unknown preselection type!" << endl;
		exit( EXIT_FAILURE );
	}
	return false;
}

*/
