#ifndef LEPTON_H
#define LEPTON_H

class TLorentzVector;

class Lepton {
	public :
		float px;
		float py;
		float pz;
		float en;
		float ptErr;
		float ecalIso;
		float hcalIso;
		float trkIso;
		float gIso;
		float chIso;
		float puchIso;
		float nhIso;
		int id;
		int genid;
		float ensf;
		float ensferr;
		float d0;
		float dZ;
		float ip3d;
		float trkpt;
		float trketa;
		float trkphi;
		float trkchi2;
		float trkValidPixelHits;
		float trkValidTrackerHits;
		float trkLostInnerHits;

		Lepton( float px_, float py_, float pz_, float en_, float ptErr_, float ecalIso_, float hcalIso_,
				float trkIso_, float gIso_, float chIso_, float puchIso_, float nhIso_, int id_, int genid_,
				float ensf_, float ensferr_, float d0_, float dZ_, float ip3d_, float trkpt_, float trketa_, float trkphi_,
				float trkchi2_, float trkValidPixelHits_, float trkValidTrackerHits_, float trkLostInnerHits_ );		
		TLorentzVector lorentzVector() const;
		double detIsolation(double rho) const;
};

#endif // LEPTON_H
