#include "KStruct.h"

KStruct::KStruct()
{
Clear();
}

void KStruct::Clear()
{
Int_t i=0;
for(i=0;i<MAXPOINT;i++) 
  {
  Xtrack[i]=0; Ytrack[i]=0; Ztrack[i]=0; 
  Charge[i]=0; Time[i]=0; 
  Efield[i]=0; MulCar[i]=0;
  }
Ylenght=0; Xlenght=0; TTime=0; Steps=0; DStrip=0; TCharge=0;

}

void KStruct::Info()
{
printf("\nParticle Cahrge= %d\n",PCharge);
printf("Number of Steps= %d\n",Steps);
printf("Drift Strip    = %d\n",DStrip);
printf("Pathlenght X   = %f\n",Xlenght);
printf("Pathlenght Y   = %f\n",Ylenght);
printf("Pathlenght Z   = %f\n",Zlenght);
printf("Total Time     = %f\n",TTime*1e9);
printf("Total Charge   = %f\n",TCharge);
}

void KStruct::Draw(Char_t *option)
{
  Char_t *name="Vector Plot";
  Float_t *x,*y;
  if(! strcmp(option,"xy") || ! strcmp(option,"yx")) {x=Xtrack; y=Ytrack;}
  if(! strcmp(option,"xt") || ! strcmp(option,"tx")) {x=Time; y=Xtrack;}
  if(! strcmp(option,"xc") || ! strcmp(option,"cx")) {x=Xtrack; y=Charge;}
  if(! strcmp(option,"yt") || ! strcmp(option,"ty")) {x=Time; y=Ytrack;}
  if(! strcmp(option,"yc") || ! strcmp(option,"cy")) {x=Ytrack; y=Charge;}
  if(! strcmp(option,"tc") || ! strcmp(option,"ct")) {x=Time; y=Charge;}
  
   TGraph *gr = new TGraph(Steps-1,&x[3],&y[3]);
  gr->SetTitle(name);
  gr->Draw("AL"); gr->GetHistogram()->Draw(); gr->Draw("AL*");
  
}

void KStruct::GetCH(TH1F *histo,Int_t Update,Float_t Mult, Float_t tau)
{
  //
  //
  //   Float_t tau;  trapping time [s]

Int_t i;
TH1F *his=new TH1F("ch","Charge Histogram",histo->GetNbinsX(),histo->GetXaxis()->GetXmin(),histo->GetXaxis()->GetXmax());
Double_t *ch=new Double_t[Steps+1];
Axis_t *ti=new Axis_t[Steps+1]; // Changed when migrating from 2.23 to 2.25 or higher
 for(i=1;i<Steps+1;i++) {ch[i]=Charge[i]; ti[i]=Time[i];} // Changed when migrating from 2.23 to 2.25
his->FillN(Steps,&ti[1],&ch[1]); // Changed when migrating from 2.23 to 2.25

// Trappping is included if tau>0   // added 20.1.2010
if(tau>0)
  {
    for(i=1;i<his->GetNbinsX();i++) 
      his->SetBinContent(i,his->GetBinContent(i)*TMath::Exp(-(his->GetBinCenter(i)-Time[0])/tau));
  }
////////////////////////////////////////////////////////


if(Update) {his->Scale(Mult); histo->Add(his); } else his->Copy(*histo);

delete his;
delete [] ti; // Changed when migrating from 2.23 to 2.25
delete [] ch;  // Changed when migrating from 2.23 to 2.25
}

Float_t KStruct::GetCHMult(TH1F *histo,Int_t Update,Float_t Mult, Float_t tau)
{
  // 
  //
  //   Float_t tau;  trapping time [s]
  Double_t dx,dy,dif,dift,sum=1,summ=1.,mulf,traf;
Int_t i;
TH1F *his=new TH1F("ch","Charge Histogram",histo->GetNbinsX(),histo->GetXaxis()->GetXmin(),histo->GetXaxis()->GetXmax());
Double_t *ch=new Double_t[Steps+1];
Double_t *Multi=new Double_t[Steps+1];
Axis_t *ti=new Axis_t[Steps+1]; // Changed when migrating from 2.23 to 2.25 or higher


 for(i=1;i<Steps+1;i++) 
   {
     if(PCharge<0) dif=KAlpha(0.5*(Efield[i+1]+Efield[i]),-1,KMaterial::ImpactIonization);
     else  dif=KAlpha(0.5*(Efield[i+1]+Efield[i]),1,KMaterial::ImpactIonization);
   
   
 dift=Time[i+1]-Time[i];
 dx=TMath::Sqrt(TMath::Power((Xtrack[i+1]-Xtrack[i]),2)+TMath::Power((Ytrack[i+1]-Ytrack[i]),2));

 //---- Calculation of multiplication and trapping factors in given step //
 mulf=(1+dif*dx);
 traf=(1-dift/tau);
 MulCar[i]=(mulf-1)*sum*traf; 
 summ*=mulf; sum*=mulf;  // multiplication 8.3.2011
 if(tau>0) sum*=traf;    // trapping 8.3.2011
 //---- END ---- //

 ch[i]=Charge[i]*sum; ti[i]=Time[i];
// Trappping is included if tau>0   // added 20.1.2010
////////////////////////////////////////////////////////

// printf("%d :: X=%4.1f , Y=%4.1f :: E=%4.2e ::  Time=%4.1e ; Charge=%4.1e ; dif=%5.2e ; MultT=%5.4e Mult=%5.4f hole=%5.3e\n",i,Xtrack[i],Ytrack[i],Efield[i],ti[i],ch[i],dif,sum,summ,MulCar[i]);
   } 


// Changed ti time when migrating from 2.23 to 2.25
his->FillN(Steps,&ti[1],&ch[1]); // Changed when migrating from 2.23 to 2.25

if(Update) {his->Scale(Mult); histo->Add(his); } else his->Copy(*histo);




delete his;
delete [] ti; // Changed when migrating from 2.23 to 2.25
delete [] ch;  // Changed when migrating from 2.23 to 2.25
delete [] Multi; // 8.3.2011
return summ;
}


//  void KStruct::GetGraph(TGraph *gr,Int_t What)
//  {
//  Int_t i;
//  gr=new TGraph(Steps,Xtrack,Ytrack);
//  gr->SetLineColor(1);
//  }


TH1D *KStruct::GetElFieldAlongTheDrift()
{
  Double_t Step=TMath::Sqrt(TMath::Power(Xtrack[1]-Xtrack[0],2)+TMath::Power(Ytrack[1]-Ytrack[0],2)+TMath::Power(Ztrack[1]-Ztrack[0],2));

  TH1D *his=new TH1D("E-Track","E-Track",Steps,0,Steps*Step);
  for(Int_t i=0;i<Steps;i++)
    his->SetBinContent(i,Efield[i]);
  return his;
}
