#ifndef _KPixel
#define _KPixel

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// KPixel                                                               //
//                                                                      //
// Description of the pixel detector                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////


#include "TArray.h"
#include "TArrayI.h"
#include "TArrayF.h"
#include "TGraph.h"
#include <stdio.h>
#include <stdlib.h>
#include "KStruct.h"
#include "TMinuit.h"
#include "KDetector.h"

class KPixel : public KDetector {

private:
public:
  Int_t Pix;
  Float_t CellZ;
  Float_t CellX;
  Float_t CellY;

  Float_t *PSx; //[Pix]
  Float_t *PSy; //[Pix]
  Float_t *PSWx;//[Pix]
  Float_t *PSWy;//[Pix]
  Float_t *PSd; //[Pix]
  Short_t *PSW; //[Pix]
  

  KPixel(Int_t, Float_t=200, Float_t=50, Float_t=125);
  ~KPixel();
  void SetUpVolume(Float_t, Float_t, Float_t);
  void SetUpPixel(Int_t, Float_t, Float_t, Float_t, Float_t, Float_t, Short_t);
  void SetUpElectrodes(Int_t=0);
  

  ClassDef(KPixel,1) 
};


#endif

