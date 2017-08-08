#include "KPixel.h"

ClassImp(KPixel)
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// KPixel                                                               //
//                                                                      //
// Description of the pixel detector                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

KPixel::KPixel(Int_t x1, Float_t x, Float_t y, Float_t z)
{
  Int_t i;
  Pix=x1;
  PSx=new Float_t [Pix];
  PSy=new Float_t [Pix];
  PSd=new Float_t [Pix];
  PSWx=new Float_t [Pix];
  PSWy=new Float_t [Pix];
  PSW=new Short_t [Pix];

  for(i=0;i<Pix;i++) 
    { PSx[i]=0; PSy[i]=0; PSd[i]=0; PSWx[i]=0; PSWy[i]=0; PSW[i]=0;}

  CellZ=z;
  CellX=x;
  CellY=y;

}

KPixel::~KPixel()
{
    delete PSx; 
    delete PSy; 
    delete PSd;
    delete PSWx;
    delete PSW;
}


void KPixel::SetUpVolume(Float_t St1, Float_t St2, Float_t St3)
{
  nx=(int)(CellX/St1);
  ny=(int)(CellY/St2);
  nz=(int)(CellZ/St3);
  
 EG=new TH3I("EG","EG",nx,0,CellX,ny,0,CellY,nz,0,CellZ);
  EG->GetXaxis()->SetTitle("x [#mum]");
  EG->GetYaxis()->SetTitle("y [#mum]");
  EG->GetZaxis()->SetTitle("z [#mum]");

 DM=new TH3I("DM","DM",nx,0,CellX,ny,0,CellY,nz,0,CellZ);
  DM->GetXaxis()->SetTitle("x [#mum]");
  DM->GetYaxis()->SetTitle("y [#mum]");
  DM->GetZaxis()->SetTitle("z [#mum]");


}

void KPixel::SetUpElectrodes(Int_t Material)
{
  Int_t i,j,k,q;
  Int_t xpl,ypl,zpl,xpr,ypr,zpr;
  for(k=1;k<=nz;k++)
    for(j=1;j<=ny;j++)
      for(i=1;i<=nx;i++)
	{
	if(k==1) 
	  {
	    EG->SetBinContent(i,j,k,2); 
	    DM->SetBinContent(i,j,k,100); 
	  }
	else 
	  EG->SetBinContent(i,j,k,0);
  	  DM->SetBinContent(i,j,k,Material);
	}

for(Int_t q=0;q<Pix;q++)
  {

           xpl=EG->GetXaxis()->FindBin(PSx[q]-PSWx[q]);
	   ypl=EG->GetYaxis()->FindBin(PSy[q]-PSWy[q]);
	   zpl=nz;
	   if(xpl<1) xpl=1; if(ypl<1) ypl=1;


           xpr=EG->GetXaxis()->FindBin(PSWx[q]+PSx[q]);
	   ypr=EG->GetYaxis()->FindBin(PSWy[q]+PSy[q]);
	   zpr=EG->GetZaxis()->FindBin(CellZ-PSd[q]);
	   if(xpl>nx) xpr=nx; if(ypl>ny) ypr=ny; 

	   // printf("Hole %d:: Bins:  X(%d %d) Y(%d %d) Z(%d %d)\n",q,xpl,xpr,ypl,ypr,zpl,zpr);

  for(k=zpl;k>=zpr;k--)	   
    for(j=ypl;j<=ypr;j++)
      for(i=xpl;i<=xpr;i++)
	{
	  //	  if(x0+i>1 && x0+i<=nx && y0+i>1 && y0+i<=ny)
	  EG->SetBinContent(i,j,k,PSW[q]);
	}
	
 enp[0]=CellX/2;  exp[0]=enp[0];
 enp[1]=1;        exp[1]=CellY;
 enp[2]=CellZ/2;  exp[2]=CellZ/2;

  }



  

}

void KPixel::SetUpPixel(Int_t n, Float_t posX, Float_t posY, Float_t WX, Float_t WY, Float_t Depth, Short_t Weigth)
{
  PSx[n]=posX;
  PSy[n]=posY;
  PSWx[n]=WX;
  PSWy[n]=WY;
  PSd[n]=Depth;
  PSW[n]=Weigth;
}
