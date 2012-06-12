#ifndef JET_H
#define JET_H

class TLorentzVector;

class Jet {
	public :
		float px;
		float py;
		float pz;
		float en;
		float btag;
		float genpt;

		Jet( float px_, float py_, float pz_, float en_, float btag_, float ptgen_ );		
		TLorentzVector lorentzVector() const;
};

#endif // JET_H
