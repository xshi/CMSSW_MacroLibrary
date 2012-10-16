/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ROOZPTPDF
#define ROOZPTPDF

#include <RooAbsPdf.h>
#include <RooRealProxy.h>
#include <RooCategoryProxy.h>
#include <RooAbsReal.h>
#include <RooAbsCategory.h>

class RooZPtPdf : public RooAbsPdf {
	public:
		RooZPtPdf() {} ; 
		RooZPtPdf(const char *name, const char *title,
				RooAbsReal & _zpt,
				RooAbsReal & _a1,
				RooAbsReal & _b1,
				RooAbsReal & _c1,
				RooAbsReal & _d1,
				RooAbsReal & _b2,
				RooAbsReal & _c2,
				RooAbsReal & _b3,
				RooAbsReal & _c3,
				RooAbsReal & _b4,
				RooAbsReal & _c4
				);
		RooZPtPdf(const RooZPtPdf& other, const char* name = 0) ;
		inline virtual TObject* clone(const char* newname) const { return new RooZPtPdf(*this, newname); }
		inline virtual ~RooZPtPdf() { }
		Double_t evaluate() const ;
		ClassDef(RooZPtPdf,2);
	protected:
		RooRealProxy zpt ;
		RooRealProxy a1;
		RooRealProxy b1;
		RooRealProxy c1;
		RooRealProxy d1;
		RooRealProxy b2;
		RooRealProxy c2;
		RooRealProxy b3;
		RooRealProxy c3;
		RooRealProxy b4;
		RooRealProxy c4;
};

#endif
