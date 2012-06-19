#include "lepton.h"

#ifndef MUON_H
#define MUON_H

class Muon : public Lepton {
	public :
		int idbits;
		float nMatches;
		float nMatchedStations;
		float validMuonHits;
		float innerTrackChi2;
		float trkLayersWithMeasurement;
		float pixelLayersWithMeasurement;

		Muon( float px_, float py_, float pz_, float en_, float ptErr_, float ecalIso_, float hcalIso_,
				float trkIso_, float gIso_, float chIso_, float puchIso_, float nhIso_, int id_, int genid_,
				float ensf_, float ensferr_, float d0_, float dZ_, float ip3d_, float trkpt_, float trketa_, float trkphi_,
				float trkchi2_, float trkValidPixelHits_, float trkValidTrackerHits_, float trkLostInnerHits_,
				int idbits_, float nMatches_, float nMatchedStations_, float validMuonHits_, float innerTrackChi2_,
				float trkLayersWithMeasurement_, float pixelLayersWithMeasurement_ );

		bool isTrackerMuon() const;
		bool isGlobalMuon() const;
		bool isPFMuon() const;
		bool isTMOneStationTight() const;
		bool isLooseMuon() const;
		bool isSoftMuon() const;
		bool isTightMuon() const;
		bool isTrackIsolatedLoose() const;
		bool isTrackIsolatedTight() const;
		bool isPFIsolatedLoose() const;
		bool isPFIsolatedTight() const;
};

#endif // MUON_H
