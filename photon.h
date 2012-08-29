#ifndef PHOTON_H
#define PHOTON_H

class TLorentzVector;

class Photon {
	public :
		float px;
		float py;
		float pz;
		float en;
		float iso1;
		float iso2;
		float iso3;
		float sihih;
		float sipip;
		float r9;
		float hoe;
		float htoe;
		float corren;
		float correnerr;
		int idbits;

		Photon( float px_, float py_, float pz_, float en_, float iso1_, float iso2_, float iso3_, float sihih_, float sipip_,
			 float r9_, float hoe_, float htoe_, float corren_, float correnerr_, int idbits_ );
		TLorentzVector lorentzVector() const;
		bool isSelected();
};

#endif // PHOTON_H
