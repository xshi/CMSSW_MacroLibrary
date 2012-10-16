#include <Riostream.h> 
#include "RooZPtPdf.h" 
#include "toolbox.h"

ClassImp(RooZPtPdf);

RooZPtPdf::RooZPtPdf(const char *name, const char *title, RooAbsReal & _zpt,
		RooAbsReal & _a1, RooAbsReal & _b1,
		RooAbsReal & _c1, RooAbsReal & _d1,
		RooAbsReal & _b2, RooAbsReal & _c2,
		RooAbsReal & _b3, RooAbsReal & _c3,
		RooAbsReal & _b4, RooAbsReal & _c4) :
	RooAbsPdf(name,title), zpt("zpt", "zpt", this, _zpt), a1("a1", "a1", this, _a1),
	b1("b1", "b1", this, _b1), c1("c1", "c1", this, _c1), d1("d1", "d1", this, _d1),
	b2("b2", "b2", this, _b2), c2("c2", "c2", this, _c2),
	b3("b3", "b3", this, _b3), c3("c3", "c3", this, _c3),
	b4("b4", "b4", this, _b4), c4("c4", "c4", this, _c4) {}

RooZPtPdf::RooZPtPdf(const RooZPtPdf& other, const char* name) : 
		RooAbsPdf(other,name), zpt("zpt",this,other.zpt), a1("a1",this,other.a1),
		b1("b1",this,other.b1), c1("c1",this,other.c1), d1("d1",this,other.d1),
		b2("b2",this,other.b2), c2("c2",this,other.c2),
		b3("b3",this,other.b3), c3("c3",this,other.c3),
		b4("b4",this,other.b4), c4("c4",this,other.c4) {}



Double_t RooZPtPdf::evaluate() const {
	return ptFunc(a1, b1, c1, d1, b2, c2, b3, c3, b4, c4, zpt);
}
