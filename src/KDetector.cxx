#include "TPolyLine3D.h"
#include "KDetector.h"
#include "TFile.h"

#define ABS(x) x>0?x:-x
#define PREDZNAK(x) x>0?1:-1

// Declaration of macros used to setup the solver of partial differential equations

// #define C1(x,y,z) y3[n]=-2./(x*x)-2./(y*y)-2./(z*z);

// #define L1(x) y2[n]=1./(x*x);
// #define R1(x) y4[n]=1./(x*x);
// #define U1(x) y5[n]=1./(x*x);
// #define D1(x) y6[n]=1./(x*x);
// #define I1(x) y7[n]=1./(x*x);
// #define O1(x) y8[n]=1./(x*x);

// #define L2(x) y2[n]=2./(x*x);
// #define R2(x) y4[n]=2./(x*x);
// #define U2(x) y5[n]=2./(x*x);
// #define D2(x) y6[n]=2./(x*x);
// #define I2(x) y7[n]=2./(x*x);
// #define O2(x) y8[n]=2./(x*x);


#define C1(x,y,z) y3[n]=x+y+z;

#define L1(x) y2[n]=x;
#define R1(x) y4[n]=x;
#define U1(x) y5[n]=x;
#define D1(x) y6[n]=x;
#define I1(x) y7[n]=x;
#define O1(x) y8[n]=x;

#define L2(x) y2[n]=2*x;
#define R2(x) y4[n]=2*x;
#define U2(x) y5[n]=2*x;
#define D2(x) y6[n]=2*x;
#define I2(x) y7[n]=2*x;
#define O2(x) y8[n]=2*x;

#define C0 y3[n]=1.;
#define U0 y5[n]=0.;
#define D0 y6[n]=0.;
#define R0 y4[n]=0.;
#define L0 y2[n]=0.;
#define I0 y7[n]=0.;
#define O0 y8[n]=0.;

#define  PI  3.1415927
#define EPS 1.0e-14
#define STEP_DET_CUR 25e-9

ClassImp(KDetector)
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// KDetector                                                            //
//                                                                      //
// The base class for all detectors. It incorporates the calculation of //
// electric and weithing field as well as simulation of drift.          //
//                                                                      //
// Calculation of the drift:                                            //
// BEGIN_LATEX                                                          
//  #frac{d^{2}f(x)}{dx^{2}}=#frac{f(x+h1)-2*f(x+0.5(h1-h2))+f(x-h2)}{(0.5 h1+ 0.5 h2)^{2}}
//  h1=h2=h  
//  #frac{d^{2}f(x)}{dx^{2}}=#frac{f(x+h)-2*f(x)+f(x-h)}{h^{2}}
// END_LATEX                                                            
//////////////////////////////////////////////////////////////////////////



double **a,*b,*y6,*y2,*y3,*y4,*y5,*y7,*y8;
double *dvector(long, long);
void free_dvector(double*,long,long);
/**********************************************************
 snrm: Calculates the norm of vector 
 unsigned long n - dimension of the vector
 double sx[]     - components
 itol            - <=3 real norm , otherwise max component.
************************************************************/

double snrm(unsigned long n, double sx[], int itol)
{
	unsigned long i,isamax;
	double ans;

	if (itol <= 3) {
		ans = 0.0;
		for (i=1;i<=n;i++) ans += sx[i]*sx[i];
		return sqrt(ans);
	} else {
		isamax=1;
		for (i=1;i<=n;i++) {
			if (fabs(sx[i]) > fabs(sx[isamax])) isamax=i;
		}
		return fabs(sx[isamax]);
	}
}
/* (C) Copr. 1986-92 Numerical Recipes Software &c&):)+!. */
/**********************************************************
 linbcg: Solves linear eqution sparse system Ax=b
 unsigned long n - number of equations
 int nx - # of points in x direction
 int ny - # of points in y direction
 double b[] - right side vector
 double x[] - solution of the system
  int itol  - a way to calculte norm
 double tol
 int itmax  - max # of iterations
 int *iter -
 double *err -
************************************************************/


/***************************************************** 
atimes : How to multiply the vector with the matrices 
Modified in this form by GK 4.10.2012
******************************************************/

void atimes(unsigned long n,int dim[], double x[],double r[],int itrnsp)
{
  // This is a function used to multuply the vectors with matrices!
  // Used for calculation of the electric and ramo field
  int nx,ny,nz;
  int i,j,k,q=0;
  double C,L,D,O,R,U,I;

  nx=dim[0]; ny=dim[1]; nz=dim[2];

  for(k=1; k<=nz; k++)
   for(j=1; j<=ny; j++)          /*mnozenje po stolpcu*/ 
     for(i=1; i<=nx; i++)           /*mnozenje po vrstici*/ 
       { 
	 q++;
	                 C=y3[q]*x[q];
         if(q-1>1)       L=y2[q]*x[q-1]; else L=0;
         if(q-nx>1)      D=y6[q]*x[q-nx]; else D=0;
	 if((q-nx*ny)>1) I=y7[q]*x[q-ny*nx]; else I=0;

         if(q+1<=n)       R=y4[q]*x[q+1]; else R=0;
         if(q+nx<=n)      U=y5[q]*x[q+nx]; else U=0;
	 if((q+nx*ny)<=n) O=y8[q]*x[q+ny*nx]; else O=0;
	 
	 r[q]=C+L+D+O+R+U+I;
       }
   if(n!=q) printf("\n Error in matrix solving!");
  return;
}


/***************************************************
Solve a very simple system of equations 
with the diagonal elements to be the only ones 
****************************************************/
void asolve(unsigned long n, double b[], double x[], int itrnsp)
{
  //Solve a very simple system of n equations
  //with the diagonal elements to be the only ones!
	unsigned long i;

	for(i=1;i<=n;i++) x[i]=(y3[i] != 0.0 ? b[i]/y3[i] : b[i]);
}

void linbcg(unsigned long n, int dim[], double b[], double x[], int itol, double tol,
 	int itmax, int *iter, double *err)
{
  // The main function for electric field calcualtion
	void asolve(unsigned long n, double b[], double x[], int itrnsp);
	void atimes(unsigned long n, int dim[],double x[], double r[], int itrnsp);
	double snrm(unsigned long n, double sx[], int itol);
	double *dvector(long, long);
	void free_dvector(double *, long, long);
	void nrerror(char error_text[]);
	unsigned long j;
	double ak,akden,bk,bkden,bknum,bnrm,dxnrm,xnrm,zm1nrm,znrm;
	double *p,*pp,*r,*rr,*z,*zz;

	p=dvector(1,n);
	pp=dvector(1,n);
	r=dvector(1,n);
	rr=dvector(1,n);
	z=dvector(1,n);
	zz=dvector(1,n);

	*iter=0;
	atimes(n,dim,x,r,0);
	for (j=1;j<=n;j++) {
		r[j]=b[j]-r[j];
		rr[j]=r[j];
	}
	atimes(n,dim,r,rr,0); // minimal residual invariant
	znrm=1.0;
	if (itol == 1) bnrm=snrm(n,b,itol);
	else if (itol == 2) {
		asolve(n,b,z,0);
		bnrm=snrm(n,z,itol);
	}
	else if (itol == 3 || itol == 4) {
		asolve(n,b,z,0);
		bnrm=snrm(n,z,itol);
		asolve(n,r,z,0);
		znrm=snrm(n,z,itol);
	} else nrerror("illegal itol in linbcg");
	asolve(n,r,z,0);
	while (*iter <= itmax) {
		++(*iter);
		zm1nrm=znrm;
		asolve(n,rr,zz,1);
		for (bknum=0.0,j=1;j<=n;j++) bknum += z[j]*rr[j];
		if (*iter == 1) {
			for (j=1;j<=n;j++) {
				p[j]=z[j];
				pp[j]=zz[j];
			}
		}
		else {
			bk=bknum/bkden;
			for (j=1;j<=n;j++) {
				p[j]=bk*p[j]+z[j];
				pp[j]=bk*pp[j]+zz[j];
			}
		}
		bkden=bknum;
		atimes(n,dim,p,z,0);
		for (akden=0.0,j=1;j<=n;j++) akden += z[j]*pp[j];
		ak=bknum/akden;
		atimes(n,dim,pp,zz,1);
		for (j=1;j<=n;j++) {
			x[j] += ak*p[j];
			r[j] -= ak*z[j];
			rr[j] -= ak*zz[j];
		}
		asolve(n,r,z,0);
		if (itol == 1 || itol == 2) {
			znrm=1.0;
			*err=snrm(n,r,itol)/bnrm;
		} else if (itol == 3 || itol == 4) {
			znrm=snrm(n,z,itol);
			if (fabs(zm1nrm-znrm) > EPS*znrm) {
				dxnrm=fabs(ak)*snrm(n,p,itol);
				*err=znrm/fabs(zm1nrm-znrm)*dxnrm;
			} else {
				*err=znrm/bnrm;
				continue;
			}
			xnrm=snrm(n,x,itol);
			if (*err <= 0.5*xnrm) *err /= xnrm;
			else {
				*err=znrm/bnrm;
				continue;
			}
		}
		if(*iter==itmax) {printf("\n Number of iterations exceeded the max. value \n"); 
                                  printf("iter=%4d err=%12.6f\n",*iter,*err);}

	if (*err <= tol) break;
	}

	free_dvector(p,1,n);
	free_dvector(pp,1,n);
	free_dvector(r,1,n);
    	free_dvector(rr,1,n);
	free_dvector(z,1,n);
	free_dvector(zz,1,n);
}


double KDetector::V(int val, int dowhat)
{
double voltage;
int k=0;
  if(dowhat==0) 
    {
    if(val&1) voltage=0;  
    if(val&2) voltage=Voltage;
    if(val & 32768) 
    //if bit 15 is on - many voltages
    voltage=Voltages[val>>16];
    }
  else
    {
    // numerical calculation converges faster if 1000 is used instead of 1
    // therefore the potential is scaled after calculation to 1
    if(val&16384) voltage=10000; else voltage=0;
    }
  return voltage;
}


Double_t KDetector::kappa(int i,int j, int k,  int dowhat )
{
 //Sets the effective space charge values for given point in the mesh!
  Double_t x,y,z,ret;
 
  //  if(NeffF!=NULL && NeffH!=NULL) printf("Warning:: Histogram values will be taken for Neff!\n");

    //Position in space
    x=EG->GetXaxis()->GetBinCenter(i);
    y=EG->GetYaxis()->GetBinCenter(j);
    z=EG->GetZaxis()->GetBinCenter(k);

    if(DM!=NULL) KMaterial::Mat=DM->GetBinContent(i,j,k); else KMaterial::Mat=0;

 if (dowhat==0) 
   {
     if(NeffF!=NULL)  // Neff=v enotah [um-3]
       //            ret=(NeffF->Eval(x,y,z)*1e6*e_0)/(KMaterial::Perm()*perm0); /*printf("i=%d,j=%d,y=%e\n",i,j,y);*/
            ret=(NeffF->Eval(x,y,z)*1e6*e_0)/(perm0); /*printf("i=%d,j=%d,y=%e\n",i,j,y);*/

     if(NeffH!=NULL)
       //   ret=(NeffH->GetBinContent(i,j,k)*1e6*e_0)/(KMaterial::Perm()*perm0);
            ret=(NeffH->GetBinContent(i,j,k)*1e6*e_0)/(perm0);

   }

//if (dowhat==0) if(j>nc) y=(Step*Step*1e-12)/(Si_mue*Ro*perm*perm0); else y=-(Step*Step*1e-12)/(Si_mue*Ro*perm*perm0);
 else 
    ret=0.;
return ret;
}


void KDetector::Declaration(Int_t dowhat)
{
  // New declaration for a general detector class 
  Int_t i,j,k,val;
  Double_t Rd,Ld,Dd,Ud,Od,Id;
  Double_t PRd,PLd,PDd,PUd,POd,PId;

  Double_t Xr=0,Yr=0,Zr=0;
  Double_t Xc=0,Yc=0,Zc=0;
  Double_t Xl=0,Yl=0,Zl=0;
  Int_t ii,jj,kk;
  Double_t fac;

  long n=0;
  Int_t num=nx*ny*nz;

   for (k=1;k<=nz;k++)
        for (j=1;j<=ny;j++)
		for(i=1;i<=nx;i++)
		  {
 		    n=(k-1)*nx*ny+(j-1)*nx+i; //Get index of the matrix element
  		    if(j-1<1) jj=1;  else jj=j-1;
		    if(i-1<1) ii=1;  else ii=i-1; 
		    if(k-1<1) kk=1;  else kk=k-1; 
   
		 /////////// DEFINE STEPS IN X //////////////////////////////////////
		 Rd=fabs(EG->GetXaxis()->GetBinCenter(i+1)-EG->GetXaxis()->GetBinCenter(i));
		 Ld=fabs(EG->GetXaxis()->GetBinCenter(i)-EG->GetXaxis()->GetBinCenter(i-1));
		 if(i+1>nx) Rd=Ld; if(i-1<1) Ld=Rd;

		 ////////// DEFINE PEMITIVITY IN X - normal surface ////////////////////////////
		 PRd=Perm(DM->GetBinContent(i,j,k))+Perm(DM->GetBinContent(i,jj,k));

		 if(nz!=1) 
		   {
		   PRd+=Perm(DM->GetBinContent(i,j,kk))+Perm(DM->GetBinContent(i,jj,kk));
		   PRd/=4;
		   } else PRd/=2;

		 PLd=Perm(DM->GetBinContent(ii,j,k))+Perm(DM->GetBinContent(ii,jj,k));
		 if(nz!=1) 
		   {
		     PLd+=Perm(DM->GetBinContent(ii,j,kk))+Perm(DM->GetBinContent(ii,jj,kk));
		     PLd/=4;
		   } else PLd/=2;

		 /////////// DEFINE STEPS IN Y //////////////////////////////////////

		 Ud=fabs(EG->GetYaxis()->GetBinCenter(j+1)-EG->GetYaxis()->GetBinCenter(j));
		 Dd=fabs(EG->GetYaxis()->GetBinCenter(j)-EG->GetYaxis()->GetBinCenter(j-1));
		 if(j+1>ny) Ud=Dd; if(j-1<1) Dd=Ud;

		 ////////// DEFINE PEMITIVITY IN Y ////////////////////////////
		 PUd=Perm(DM->GetBinContent(i,j,k))  +Perm(DM->GetBinContent(ii,j,k));
		 if(nz!=1) 
		   {
		     PUd+=Perm(DM->GetBinContent(i,j,kk))+Perm(DM->GetBinContent(ii,j,kk));
		     PUd/=4;
		   } else PUd/=2;

		 PDd=Perm(DM->GetBinContent(i,jj,k))+Perm(DM->GetBinContent(ii,jj,k));
		   if(nz!=1) 
		   {
		     PDd+=Perm(DM->GetBinContent(i,jj,kk))+Perm(DM->GetBinContent(ii,jj,kk));
		     PDd/=4;
		   } else PDd/=2;

		 /////////// DEFINE STEPS IN Z //////////////////////////////////////

		 Od=fabs(EG->GetZaxis()->GetBinCenter(k+1)-EG->GetZaxis()->GetBinCenter(k));
		 Id=fabs(EG->GetZaxis()->GetBinCenter(k)-EG->GetZaxis()->GetBinCenter(k-1));
		 if(k+1>nz) Od=Id; if(k-1<1) Id=Od;

		 //////////DEFINE PEMITIVITY IN Z ////////////////////////////
		 if(nz!=1)
		   {
		 POd=Perm(DM->GetBinContent(i,jj,k))+Perm(DM->GetBinContent(i,j,k))+
		     Perm(DM->GetBinContent(ii,j,k))+Perm(DM->GetBinContent(ii,jj,k));

		 PId=Perm(DM->GetBinContent(i,jj,kk))+Perm(DM->GetBinContent(i,j,kk))+
		     Perm(DM->GetBinContent(ii,j,kk))+Perm(DM->GetBinContent(ii,jj,kk));

		 POd/=4;
		 PId/=4;

		   }
		 ///////////
		 
		 if(dowhat==1) {PRd=1; PLd=1; PUd=1; PDd=1; POd=1; PId=1;}
		 
		 Xr=PRd/(0.5*Rd*(Rd+Ld));
		 Xl=PLd/(0.5*Ld*(Rd+Ld));  Xc=-(Xr+Xl);
		 Yr=PUd/(0.5*Ud*(Ud+Dd));
		 Yl=PDd/(0.5*Dd*(Ud+Dd));  Yc=-(Yr+Yl);

		 if(nz!=1)
		   {
		 Zr=POd/(0.5*Od*(Od+Id));
		 Zl=PId/(0.5*Id*(Od+Id));  Zc=-(Zr+Zl);
		   } 
		 		 
		 //////////

       		    b[n]=0.;
		    val=EG->GetBinContent(i,j,k);

		    if(nz==1) { C1(Xc,Yc,0) I0 O0 } 
		    else      { C1(Xc,Yc,Zc) I1(Zl) O1(Zr) }

		    R1(Xr) U1(Yr) L1(Xl) D1(Yl) 
		      
		      
  	            if(val&4)            {D0 b[n]-=V(EG->GetBinContent(i,j-1,k),dowhat)*Yl;}
		    if(val&8)            {U0 b[n]-=V(EG->GetBinContent(i,j+1,k),dowhat)*Yr;}
		    if(val&16)           {L0 b[n]-=V(EG->GetBinContent(i-1,j,k),dowhat)*Xl;}
		    if(val&32)           {R0 b[n]-=V(EG->GetBinContent(i+1,j,k),dowhat)*Xr;}
		    if(val&1024)         {I0 b[n]-=V(EG->GetBinContent(i,j,k-1),dowhat)*Zl;}
		    if(val&2048)         {O0 b[n]-=V(EG->GetBinContent(i,j,k+1),dowhat)*Zr;}


		    if(val&64)           {U2(Yr) D0 if(val&8)    {U0 b[n]-=V(EG->GetBinContent(i,j+1,k),dowhat)*Yr;}}
		    if(val&128)          {D2(Yl) U0 if(val&4)    {D0 b[n]-=V(EG->GetBinContent(i,j-1,k),dowhat)*Yl;}}
		    if(val&256)          {R2(Xr) L0 if(val&32)   {R0 b[n]-=V(EG->GetBinContent(i+1,j,k),dowhat)*Xr;}}
		    if(val&512)          {L2(Xl) R0 if(val&16)   {L0 b[n]-=V(EG->GetBinContent(i-1,j,k),dowhat)*Xl;}}
		    if(val&4096)         {O2(Zr) I0 if(val&2048) {O0 b[n]-=V(EG->GetBinContent(i,j,k+1),dowhat)*Zr;}}
		    if(val&8192)         {I2(Zl) O0 if(val&1024) {I0 b[n]-=V(EG->GetBinContent(i,j,k-1),dowhat)*Zl;}}


		    b[n]-=kappa(i,j,k,dowhat);
		    if(val&1 || val&2 || val>=32768)   {U0 D0 L0 R0 C0 O0 I0 b[n]=V(val,dowhat);}   
		    //if(j<=2 && i<=2) printf("stevilki: i=%d, j=%d, k=%d X=(%f %f ::%f %f), Y(%f %f :: %f %f), Z(%f %f :: %f %f) y[2,3,4,5,6,7,8]=%f %f %f %f %f %f %f :: b[n]=%f :: %d\n",i,j,k,Xr,Xl,Ld,Rd,Yr,Yl,Dd,Ud,Zr,Zl,Id,Od,y2[n],y3[n],y4[n],y5[n],y6[n],y7[n],y8[n],b[n],Mat);
		    //      if(k==nz && (j==2 || j==ny-1)) printf("stevilki: i=%d, j=%d, k=%d X=(%f %f ::%f %f), Y(%f %f :: %f %f), Z(%f %f :: %f %f) y[2,3,4,5,6,7,8]=%f %f %f %f %f %f %f :: b[n]=%f\n",i,j,k,Xr,Xl,Ld,Rd,Yr,Yl,Dd,Ud,Zr,Zl,Id,Od,y2[n],y3[n],y4[n],y5[n],y6[n],y7[n],y8[n],b[n]);
		  }

}


// void KDetector::Declaration(Int_t dowhat)
// {
//   // New declaration for a general detector class 
//   Int_t i,j,k,val;
//   Double_t Rd,Ld,Dd,Ud,Od,Id;
//   Double_t Xr=0,Yr=0,Zr=0;
//   Double_t Xc=0,Yc=0,Zc=0;
//   Double_t Xl=0,Yl=0,Zl=0;
  
//   Double_t p1,p2,p3,p4;
//   long n=0;
//   Int_t num=nx*ny*nz;

//    for (k=1;k<=nz;k++)
//         for (j=1;j<=ny;j++)
// 		for(i=1;i<=nx;i++)
// 		  {
//  		    n=(k-1)*nx*ny+(j-1)*nx+i; //Get index of the matrix element

// 		 Rd=fabs(EG->GetXaxis()->GetBinCenter(i+1)-EG->GetXaxis()->GetBinCenter(i));
// 		 Ld=fabs(EG->GetXaxis()->GetBinCenter(i)-EG->GetXaxis()->GetBinCenter(i-1));
// 		 if(i+1>nx) Rd=Ld; if(i-1<1) Ld=Rd;
// 		 Ud=fabs(EG->GetYaxis()->GetBinCenter(j+1)-EG->GetYaxis()->GetBinCenter(j));
// 		 Dd=fabs(EG->GetYaxis()->GetBinCenter(j)-EG->GetYaxis()->GetBinCenter(j-1));
// 		 if(j+1>ny) Ud=Dd; if(j-1<1) Dd=Ud;
// 		 Od=fabs(EG->GetZaxis()->GetBinCenter(k+1)-EG->GetZaxis()->GetBinCenter(k));
// 		 Id=fabs(EG->GetZaxis()->GetBinCenter(k)-EG->GetZaxis()->GetBinCenter(k-1));
// 		 if(k+1>nz) Od=Id; if(k-1<1) Id=Od;
//                  // Rd*=0.1; Ld*=0.1;
// 		 // Ud*=0.1; Dd*=0.1;
		 		 
// 		 if(Rd>=Ld) { Xr=1/(Rd*Ld); Xl=1/(Ld*Ld); } else
//           		    { Xr=1/(Rd*Rd); Xl=1/(Rd*Ld); };  Xc=-(Xr+Xl);
// 		 if(Ud>=Dd) { Yr=1/(Ud*Dd); Yl=1/(Dd*Dd); } else
// 		            { Yr=1/(Ud*Ud); Yl=1/(Ud*Dd); };  Yc=-(Yr+Yl);
// 		 if(nz!=1)  {
// 		 if(Od>=Id) { Zr=1/(Od*Id); Zl=1/(Id*Id); } else
// 		            { Zr=1/(Od*Od); Zl=1/(Od*Id); }   Zc=-(Zr+Zl);
// 		            }
		 
// 		//  //epsilon
//   		 KMaterial::Mat=DM->GetBinContent(i,j,k); p1=KMaterial::Perm();
//  		 KMaterial::Mat=DM->GetBinContent(i,j-1,k); p2=KMaterial::Perm(); Xr*=0.5*(p1+p2);
   	        
//                  KMaterial::Mat=DM->GetBinContent(i-1,j,k); p1=KMaterial::Perm();
//  		 KMaterial::Mat=DM->GetBinContent(i-1,j-1,k); p2=KMaterial::Perm(); Xl*=0.5*(p1+p2);

//                  KMaterial::Mat=DM->GetBinContent(i,j,k); p1=KMaterial::Perm();
//  		 KMaterial::Mat=DM->GetBinContent(i-1,j,k); p2=KMaterial::Perm(); Yr*=0.5*(p1+p2);

//                  KMaterial::Mat=DM->GetBinContent(i-1,j-1,k); p1=KMaterial::Perm();
//  		 KMaterial::Mat=DM->GetBinContent(i,j-1,k); p2=KMaterial::Perm(); Yl*=0.5*(p1+p2);

//                  KMaterial::Mat=DM->GetBinContent(i-1,j-1,k); p1=KMaterial::Perm();
//  		 KMaterial::Mat=DM->GetBinContent(i,j-1,k); p2=KMaterial::Perm(); 
//  		 KMaterial::Mat=DM->GetBinContent(i-1,j,k); p3=KMaterial::Perm(); 
//  		 KMaterial::Mat=DM->GetBinContent(i,j,k); p4=KMaterial::Perm(); 


	      
//        		    b[n]=0.;
// 		    val=EG->GetBinContent(i,j,k);

// 		    if(nz==1) 
//                           { 
// 			    C1(Xc,Yc,0) I0 O0 
// 			    y3[n]=-1./(Ud*Dd)*(p1+p2+p3+p4);
// 		          } 
// 		    else  
// 		          { C1(Xc,Yc,Zc) I1(Zl) O1(Zr) }

// 		    R1(Xr) U1(Yr) L1(Xl) D1(Yl) 
		      
		      
//   	            if(val&4)            {D0 b[n]-=V(EG->GetBinContent(i,j-1,k),dowhat)*Yl;}
// 		    if(val&8)            {U0 b[n]-=V(EG->GetBinContent(i,j+1,k),dowhat)*Yr;}
// 		    if(val&16)           {L0 b[n]-=V(EG->GetBinContent(i-1,j,k),dowhat)*Xl;}
// 		    if(val&32)           {R0 b[n]-=V(EG->GetBinContent(i+1,j,k),dowhat)*Xr;}
// 		    if(val&1024)         {I0 b[n]-=V(EG->GetBinContent(i,j,k-1),dowhat)*Zl;}
// 		    if(val&2048)         {O0 b[n]-=V(EG->GetBinContent(i,j,k+1),dowhat)*Zr;}


// 		    if(val&64)           {U2(Yr) D0 if(val&8)    {U0 b[n]-=V(EG->GetBinContent(i,j+1,k),dowhat)*Yr;}}
// 		    if(val&128)          {D2(Yl) U0 if(val&4)    {D0 b[n]-=V(EG->GetBinContent(i,j-1,k),dowhat)*Yl;}}
// 		    if(val&256)          {R2(Xr) L0 if(val&32)   {R0 b[n]-=V(EG->GetBinContent(i+1,j,k),dowhat)*Xr;}}
// 		    if(val&512)          {L2(Xl) R0 if(val&16)   {L0 b[n]-=V(EG->GetBinContent(i-1,j,k),dowhat)*Xl;}}
// 		    if(val&4096)         {O2(Zr) I0 if(val&2048) {O0 b[n]-=V(EG->GetBinContent(i,j,k+1),dowhat)*Zr;}}
// 		    if(val&8192)         {I2(Zl) O0 if(val&1024) {I0 b[n]-=V(EG->GetBinContent(i,j,k-1),dowhat)*Zl;}}


// 		    b[n]-=kappa(i,j,k,dowhat);
// 		    if(val&1 || val&2)   {U0 D0 L0 R0 C0 O0 I0 b[n]=V(val,dowhat);}   
// 		    	    if(j>=79 && j<82) printf("stevilki: i=%d, j=%d, k=%d X=(%f %f ::%f %f), Y(%f %f :: %f %f), Z(%f %f :: %f %f) y[2,3,4,5,6,7,8]=%f %f %f %f %f %f %f :: b[n]=%f :: %d\n",i,j,k,Xr,Xl,Ld,Rd,Yr,Yl,Dd,Ud,Zr,Zl,Id,Od,y2[n],y3[n],y4[n],y5[n],y6[n],y7[n],y8[n],b[n],Mat);
// 		    //      if(k==nz && (j==2 || j==ny-1)) printf("stevilki: i=%d, j=%d, k=%d X=(%f %f ::%f %f), Y(%f %f :: %f %f), Z(%f %f :: %f %f) y[2,3,4,5,6,7,8]=%f %f %f %f %f %f %f :: b[n]=%f\n",i,j,k,Xr,Xl,Ld,Rd,Yr,Yl,Dd,Ud,Zr,Zl,Id,Od,y2[n],y3[n],y4[n],y5[n],y6[n],y7[n],y8[n],b[n]);
// 		  }

// }


void KDetector::CalField(Int_t what)
{
double err;
int iteracije,i;
int num=nx*ny*nz,dim[3];
Double_t *x;
//booking memory
 b=dvector(1,num); y6=dvector(1,num); y2=dvector(1,num); 
 y3=dvector(1,num); y4=dvector(1,num); y5=dvector(1,num);
 y7=dvector(1,num); y8=dvector(1,num);
 x=dvector(1,num);    
// Setting up the boundary conditions
 printf("Setting up matrix ... \n");
 Declaration(what);
 printf("Solving matrix ...\n");
// matrix solving
for(i=1;i<=num;i++) x[i]=1.;
dim[0]=nx; dim[1]=ny; dim[2]=nz;
linbcg(num,dim,b,x,1,CalErr,MaxIter,&iteracije,&err);
// Calculating the field
if(!what)
   {
Real.U=MapToGeometry(x);
Real.CalField();
   }
 else
   {
 // Scale back to 1 from 1000 when storing the potential
 Ramo.U=MapToGeometry(x,1e-4);

Ramo.CalField();
   }
// Freeing the memory
 free_dvector(x, 1,num);   free_dvector(b, 1,num);   free_dvector(y2, 1,num); 
 free_dvector(y3, 1,num);  free_dvector(y4, 1,num);  free_dvector(y5, 1,num); 
 free_dvector(y6, 1,num);  free_dvector(y7, 1,num);  free_dvector(y8, 1,num);
}


void KDetector::Drift(Double_t sx, Double_t sy, Double_t sz, Float_t charg, KStruct *seg, Double_t t0)
{
  //Drift simulation for a point charge (Float_t charg;)
  //starting from ( sx,sy, sz)
  //KStruct *seg is the structure where the  the drift paths, drift times and induced cahrges are stored
 
  Double_t Stime=0;                       // Step time
  Double_t difx=0,dify=0,difz=0;          // diffusion steps in all directions
  Double_t sigma;                         // sigma of the diffusion step  
  Float_t *xcv,*ycv,*zcv,*time,*charge;   // drift variables to be stored in KStruct
  Double_t cx=0,cy=0,cz=0;                // current position of the charge bucket
  Double_t vel=0;                         // drift velocity
  Double_t t=0;                           // drift time
  Double_t sumc=0;                        // total induced charge
  Double_t ncx=0,ncy=0,ncz=0;             // next position of the charge bucket   
  Double_t deltacx,deltacy,deltacz;       // drift step due to drift
  
  Int_t st=0;                             // current step
  Int_t ishit=0;                          // local counter of the step
  TVector3 *EE;                           // Set up electric field vector
  TVector3 *EEN;                          // Set up electric field vector
  TVector3 FF;                            // Combined drift field
  Float_t pathlen=0;                      // pathlength
  Float_t WPot;                           // current ramo potential

// Inclusion of Magnetic field 28.8.2001 - revised 15.10.2012
TVector3 BB(B);                           // Create a magnetic field vector
Float_t muhe=1650;                        // parametrization of the magnetic field 
Float_t muhh=310;	                  // is based on simply this mobility parameters

 // Start time in the  absolute domain (used for delayed charge generation in multiplication 
 
  t=t0;

// Intitialize KStruct class and its members 
 
  xcv=seg->Xtrack; 
  ycv=seg->Ytrack;  
  zcv=seg->Ztrack; 
  time=seg->Time; 
  charge=seg->Charge; 
  seg->Clear();
  seg->PCharge=(Int_t) charg;

// start drift

cx=sx; cy=sy; cz=sz;                   // set current coordinates
xcv[st]=cx; ycv[st]=cy; zcv[st]=cz;    // put the first point in the KStruct 
time[st]=t;  charge[st]=0; 
           
EE=Real.CalFieldXYZ(cx,cy,cz);         // Get the electric field vector 
seg->Efield[st]=EE->Mag();             // Store the magnitude of E field

//printf("%d %f : (%f %f %f)\n",st,charg,cx,cy,cz);

do 
    {
      //    printf("Calculate field\n");
      st++;
      if(charg>0)
	FF=(*EE)+muhh*EE->Cross(BB); else 
	FF=(*EE)-muhe*EE->Cross(BB); 
      // "-muhe" stands for the fact that at the same 
      // field the drift direction has changed due to different charge

      //printf("Field : %f %f %f (%f %f %f)  ---- ",FF[0],FF[1],FF[2],(*EE)[0],(*EE)[1],(*EE)[2]);
      if(FF.Mag()!=0)
	{
      deltacy=-SStep*charg*FF.y()/FF.Mag();   // deltay 
      deltacx=-SStep*charg*FF.x()/FF.Mag();   // deltax
      deltacz=-SStep*charg*FF.z()/FF.Mag();   // deltaz
	}
      else { deltacx=0; deltacy=0; deltacz=0; }
	
      //    printf("Calculate velocity \n");

      if(DM!=NULL)
	KMaterial::Mat=DM->GetBinContent(DM->FindBin(cx,cy,cz)); 
      else KMaterial::Mat=0;
      
      EEN=Real.CalFieldXYZ(cx+deltacx,cy+deltacy,cz+deltacz); // get field & velocity at new location
      vel=Real.DriftVelocity( (EEN->Mag()+EE->Mag())/2,charg,Temperature,TMath::Abs(NeffF->Eval(cx,cy,cz)),MobMod());

      //      printf("Calculate vel: %e EEN = %e ::: ",vel, EEN->Mag());
      if(vel==0) {
	          deltacx=0; deltacy=0; deltacz=0; 
		  difx=0; dify=0; difz=0;
	          ishit=10;
                 } 
      else 
       if(diff)                      // is diffusion ON
	{
	  Stime=SStep*1e-4/vel; // calcualte step time  
	  sigma=TMath::Sqrt(2*Kboltz*Real.Mobility(EE->Mag(),Temperature,charg,TMath::Abs(NeffF->Eval(cx,cy,cz)),MobMod())*Temperature*Stime); 
	  dify=ran->Gaus(0,sigma)*1e4; 
	  difx=ran->Gaus(0,sigma)*1e4;
	  if(nz!=1) difz=ran->Gaus(0,sigma)*1e4; else difz=0;
	}

       if((cx+deltacx+difx)>=GetUpEdge(0)) ncx=GetUpEdge(0); else
       if((cx+deltacx+difx)<GetLowEdge(0)) ncx=GetLowEdge(0); else
	 ncx=cx+(deltacx+difx);;
 
       if((cy+deltacy+dify)>=GetUpEdge(1)) ncy=GetUpEdge(1); else
       if((cy+deltacy+dify)<GetLowEdge(1)) ncy=GetLowEdge(1); else
	 ncy=cy+(deltacy+dify);

       if((cz+deltacz+difz)>=GetUpEdge(2)) ncz=GetUpEdge(2); else
       if((cz+deltacz+difz)<GetLowEdge(2)) ncz=GetLowEdge(2); else
	 ncz=cz+(deltacz+difz);

       if(Debug) printf("%d %f : x:%f->%f y:%f->%f z:%f->%f (%f %f %f): Mat=%d :: ",st,charg,cx,ncx,cy,ncy,cz,ncz,deltacx,deltacy,deltacz, KMaterial::Mat);

       charge[st]=charg*(Ramo.CalPotXYZ(ncx,ncy,ncz)-Ramo.CalPotXYZ(cx,cy,cz));
       cx=ncx; cy=ncy; cz=ncz;

    sumc+=charge[st];

    if(vel!=0) 
      {
	t=t+SStep*1e-4/vel; //else t+=Stime;
	pathlen+=SStep;
      }

    xcv[st]=cx;
    ycv[st]=cy;
    zcv[st]=cz;
    time[st]=t;

    EE=Real.CalFieldXYZ(cx,cy,cz);     
    seg->Efield[st]=EE->Mag();

    // Checking for termination of the drift //// 

    WPot=Ramo.CalPotXYZ(cx,cy,cz);
    if(WPot>(1-Deps)) ishit=1;
    if(TMath::Abs(WPot)<Deps) ishit=2;
    if(cx<= GetLowEdge(0)) ishit=3;
    if(cx>=  GetUpEdge(0)) ishit=4;
    if(cy<= GetLowEdge(1)) ishit=5;
    if(cy>= GetUpEdge(1))  ishit=6;
    if(cz<= GetLowEdge(2)) ishit=7;
    if(cz>= GetUpEdge(2))  ishit=8;
    if(st>=MAXPOINT-1) ishit=10;   

    if(Debug) printf("(t=%e, vel=%e) [Ch=%f ChInt=%f] Ishit=%d \n",t,vel,charge[st],sumc,ishit);

    } while (!ishit); // Do until the end of drift

   
   (*seg).Xlenght=pathlen; (*seg).Ylenght=pathlen; 
   (*seg).TTime=t; (*seg).TCharge=sumc; (*seg).Steps=st;
   delete EE; delete EEN;

return;
  }

TH1F *KDetector::Draw1D(Char_t *option, Float_t proj,Int_t axis,Float_t pos)
{
  // Draws the 1D projection of the Field
  // Char_t option;  see ::Draw()
  // Char_t proj;    see ::Draw() 
  // Int_t axis;     0=x, 1=y, 2=z;
  // Float_t pos;    position along the missing coordinate if 2D z=0.5;
  TH2F *his=Draw(option,proj);
  Int_t iter,*x,*y,i,j;
  Float_t low,up;
  if(axis==0)
    {
      iter=his->GetNbinsX(); 
      up=his->GetXaxis()->GetBinUpEdge(iter);
      low=his->GetXaxis()->GetBinLowEdge(1);
    }
  else 
    {
      iter=his->GetNbinsY();
      up=his->GetYaxis()->GetBinUpEdge(iter);
      low=his->GetYaxis()->GetBinLowEdge(1);
    }
  TH1F *his1D=new TH1F("Projection","Projection",iter,low,up);

  
  for(i=1;i<=iter;i++)
    {
      if(axis==0) his1D->SetBinContent(i, his->GetBinContent(i,his->GetYaxis()->FindBin(pos)));
      else
	his1D->SetBinContent(i,  his->GetBinContent( his->GetXaxis()->FindBin(pos),i));
    }
  return his1D;

}



TH2F *KDetector::Draw(Char_t *option, Float_t proj)
{
  //The function draws weighting and electric field and also the event display
  // Float_t proj; position along the axis of projection
  // Char_t *option:
  //           Which field do you want to plot?
  //		  W -weighting  
  // 		  E -electric 
  //              G -geometry
  //              M -material  
  //              N -Neff = space charge
  //           What do do want to plot?
  //		  P - potential 
  // 		  F - |field|
  //              X - x component of the field
  //              Y - y component of the field
  //              Z - z component of the field
  //          In case of 3D which plane?
  //              yz - the cross section is along the x value of proj
  //              xz - the cross section is along the y value of proj  
  //              xy - the cross section is along the z value of proj


KField *cf;
TH3F   *ch;
TH2F *histo=new TH2F();
TString opt(option);

 Int_t i=0;
 Int_t iproj=0;
 Int_t What=1;      // default is E
 Int_t Which=0;     // default is potential
 Int_t Which3D=3;   // default is XY

 if(opt.Contains("N")) What=5;  
 if(opt.Contains("M")) What=4;  
 if(opt.Contains("G")) What=3;  
 if(opt.Contains("W")) What=2;  
 if(opt.Contains("E")) What=1;

 if(nz>1)
   {
     if(opt.Contains("xy")) Which3D=3;  
     if(opt.Contains("xz")) Which3D=2;
     if(opt.Contains("yz")) Which3D=1;
   }

   
  if(opt.Contains("P")) Which=0;
  if(opt.Contains("F")) Which=1;
  if(opt.Contains("X")) Which=2;
  if(opt.Contains("Y")) Which=3;
  if(opt.Contains("Z")) Which=4;
  
  if(What<=2) 
    {
    if(What==2) cf=&Ramo; else cf=&Real;    
       switch(Which)
	 {
	 case 0: ch=cf->U;  break;
	 case 1: ch=cf->E;  break;
	 case 2: ch=cf->Ex; break;
	 case 3: ch=cf->Ey; break;
	 case 4: ch=cf->Ez; break;
	 default: ch=cf->U;  break;
	 }
    } 
  else
    {
      if(What==3) ch=(TH3F *)EG;  
      if(What==4) ch=(TH3F *)DM;  
      if(What==5) ch=(TH3F *)NeffH;  
    }

 if(nz==1) histo=KHisProject(ch,3,1); 
 else
   {
     switch(Which3D)
       {
       case 3: iproj=ch->GetZaxis()->FindBin(proj); break;
       case 2: iproj=ch->GetYaxis()->FindBin(proj); break;
       case 1: iproj=ch->GetXaxis()->FindBin(proj); break;
       }
     histo=KHisProject(ch,Which3D,iproj);
   }

 //histo->Draw("COLZ");
return histo;
}

void KDetector::MipIR(Int_t div, Float_t lambda)
{
  // The simulation of the drift for the laser induced carriers - attenuation of the light. 
  // A track is devided into Int_ div buckets. Each bucket is drifted in the field. The
  // induced currents for each carrier is calculated as the sum  all buckets. 
  //	Int_t MobMod; mobility model
  //	Float_t B; magnetic field


Double_t data[4];
 Float_t sp[3];
 Float_t Io,len=0,lent=0,lenlam,scalef=0,pscalef=0,mule,mulh;
 int i,j,k,e;
 KStruct seg,segmul;

TH1F *histop  = new TH1F("ch-","charge+",pos->GetNbinsX(),pos->GetXaxis()->GetXmin(),pos->GetXaxis()->GetXmax()); //2.23 -> 3.0
TH1F *histon  = new TH1F("ch+","charge-",neg->GetNbinsX(),neg->GetXaxis()->GetXmin(),neg->GetXaxis()->GetXmax()); 
sum->Reset(); pos->Reset(); neg->Reset();

/////////////////////////////////////////////////////////////////
// If there is no attenuation coefficient we consider that as mip
/////////////////////////////////////////////////////////////////

 if(lambda!=0)
   {
     for(k=0;k<3;k++) lent+=TMath::Power(exp[k]-enp[k],2); lent=TMath::Sqrt(lent);
     lenlam=lent/lambda;
     Io=div/(lambda*(1-TMath::Exp(-lenlam)));
     
     //     if(Debug) printf("Calculated length=%f um, lenDIVlamda=%f, I0=%f \n", lent,lenlam,I0 ); 
   }
/////////////////////////////////////////////////////////////////


for(i=0;i<div;i++) 
  { 
    for(j=0;j<3;j++) {sp[j]=((exp[j]-enp[j])/div)*i+enp[j]+(exp[j]-enp[j])/(2*div);}
    //      printf("#i=%d div=%d, pointx=%f, pointy=%f pointz=%f \n",i,div,sp[0],sp[1],sp[2]); 
    for(j=0;j<average;j++)
      { 
      
      Drift(sp[0],sp[1],sp[2],1,&seg);
      seg.GetCH(histop,1,1,tauh); 
      Drift(sp[0],sp[1],sp[2],-1,&seg);

      if(MTresh>1)  
	{
	mule=seg.GetCHMult(histon,1,1,taue);  // performing multiplication
	//      if(Debug)  printf(":: Mstep = %f ::",mule);
	}
      else
	seg.GetCH(histon,1,1,taue); 
	
      if(mule>MTresh && MTresh>1) //if the multiplication is large enough then do the hole tracking
 	{
 	  for(e=1;e<seg.Steps+1;e++)
	    {
	      Drift(seg.Xtrack[e],seg.Ytrack[e],seg.Ztrack[e],1,&segmul,seg.Time[e]); 
	      mulh=segmul.GetCHMult(histop,1,seg.MulCar[e],tauh);
	      if(mulh>BDTresh) {printf("HOLE MULTIPLICATION - BREAKDOWN\n"); BreakDown=1;}
	    }
 	}
      }

    
      histop->Scale(1/((Float_t)average)); histon->Scale(1/((Float_t)average)); 

/////////////////////////////////////////////////////////////////
// If there is no attenuation coefficient we consider that as mip
/////////////////////////////////////////////////////////////////
      if(lambda!=0)
	{
      len=0; for(k=0;k<3;k++) len+=TMath::Power(sp[k]-enp[k],2); len=TMath::Sqrt(len); 
      scalef=Io*TMath::Exp(-len/lambda);
      //      printf("step=%d ::  len=%f um, scalef=%f pscalef=%f\n",i,len,scalef,pscalef); 
      histon->Scale(scalef); histop->Scale(scalef);
      pscalef+=scalef;
      histop->Scale(lent/div); histon->Scale(lent/div); 
	}
/////////////////////////////////////////////////////////////////////////////////

      pos->Add(histop);neg->Add(histon);
      histop->Reset();
      histon->Reset();
  
  }

sum->Add(neg);
sum->Add(pos);

delete histop;
delete histon;
}


void KDetector::ShowMipIR(Int_t div, Int_t color,Int_t how)
{
  // The simulation of the drift for the minimum ionizing particles. 
  // A track is devided into Int_ div buckets. Each bucket is drifted in the field. The
  // induced currents for each carrier is calculated as the sum  all buckets. 

TGraph *gr;
TPolyLine3D *gr3D;
Float_t sp[3];
int i,j;
KStruct seg;
TLine *line;

// Draw histograms 

 if(EG!=NULL) 
   { 
     if(nz==1) KHisProject(EG,3,how)->Draw("COL");
     else
       {
       TH3F *hh=GetGeom();
       hh->SetFillColor(color);
       hh->Draw("iso");
       }
   } 

 // Draw drift paths

for(i=0;i<div;i++) 
  {
    for(j=0;j<3;j++) {sp[j]=((exp[j]-enp[j])/div)*i+enp[j]+(exp[j]-enp[j])/(2*div);}

    // hole drift
    if(Debug)   printf("Entry Point: %f %f %f \n",sp[0],sp[1],sp[2]);
      Drift(sp[0],sp[1],sp[2],1,&seg);

      if(nz==1)
      gr=new TGraph(seg.Steps,&seg.Xtrack[1],&seg.Ytrack[1]); 
      else
      gr3D=new TPolyLine3D(seg.Steps,&seg.Xtrack[1],&seg.Ytrack[1],&seg.Ztrack[1]); 	

      if(nz==1)
	{
	  gr->SetLineColor(2);
	  gr->SetLineStyle(1);
	  gr->Draw("L");
	}
      else
	{
	 gr3D->SetLineStyle(1);  
	 gr3D->SetLineColor(2); 
	 gr3D->Draw("SAME"); 
	}

      // electron drift

      Drift(sp[0],sp[1],sp[2],-1,&seg);

      if(nz==1)
      gr=new TGraph(seg.Steps,&seg.Xtrack[1],&seg.Ytrack[1]); 
      else
      gr3D=new TPolyLine3D(seg.Steps,&seg.Xtrack[1],&seg.Ytrack[1],&seg.Ztrack[1]); 	

      if(nz==1)
	{
        gr->SetLineColor(4);
	gr->SetLineStyle(3); 
	gr->Draw("L");
	}
      else
	{
	 gr3D->SetLineStyle(1);  
	 gr3D->SetLineColor(4); 
	 gr3D->Draw("SAME"); 
	}

  }
}

KDetector::KDetector()
{
//////////////////////////////////////////////////////////////////////
// Author: Gregor Kramberger                                        //
// Default constructor for KDetector class                           //
// Default = no space charge -> default space charge is function     //
//////////////////////////////////////////////////////////////////////

NeffF=new TF3("Profile","x[0]*x[1]*x[2]*0+[0]",0,1000,0,1000,0,1000);
NeffF->SetParameter(0,0);
NeffH=NULL; //Neff histogram se to NULL

//set magnetic field to zero
for(Int_t i=0;i<3;i++) B[i]=0;

// setting up default random generator - diffusion
ran=new TRandom(33);  

// Calculation parameters
 CalErr=1e-6;
 MaxIter=2000;

// histograms for storing the drift
pos=NULL; neg=NULL; sum=NULL;
SetDriftHisto(25e-9);

// setting up general variables
// multiplication 
taue=-1;   //no electron trapping 
tauh=-1;   //no hole trapping 
MTresh=-1; //no multiplication 
BDTresh=-1;
// drift
Deps=1e-5; //precision of tracking
//MobMod=1;  //Mobility parametrization
average=1; //average over waveforms
diff=0;    //diffusion
SStep=1;   //Step size of simulation
Temperature=263; //temperature
BreakDown=0; // no breakdown
Debug=0;     // bo printing of debug information
Voltage2=0;
}


KDetector::~KDetector()
{
  //Destructor of the detector class

  if(NeffF!=NULL) delete NeffF;
  if(NeffH!=NULL) delete NeffH;
  if(ran!=NULL) delete ran;
  if(pos!=NULL) delete pos;
  if(neg!=NULL) delete neg;
  if(sum!=NULL) delete sum;
  
}


void KDetector::CalM(KStruct *seg, Double_t *data, Int_t CarrierType)
    { 

      Int_t i,j,numreg=0;
      Double_t *M=new Double_t [seg->Steps+1];
      Double_t *DIF=new Double_t [seg->Steps+1];
      Double_t dx,dy,sum=1.,sum2=0,sum1=0,dif,dif1,dif2;
      //      printf("Number of Steps=%d\n",seg->Steps);
    
      for(i=1; i<seg->Steps;i++)
	{
	  //	  if(i==seg->Steps) DIF[i]=0;
	  dx=TMath::Sqrt(TMath::Power((seg->Xtrack[i+1]-seg->Xtrack[i]),2)+TMath::Power((seg->Ytrack[i+1]-seg->Ytrack[i]),2));

	  if(CarrierType<0)
	  dif=Real.alpha(0.5*(seg->Efield[i+1]+seg->Efield[i]));
	  else
	  dif=Real.beta(0.5*(seg->Efield[i+1]+seg->Efield[i]));

	  sum*=(1+dif*dx);
	  DIF[i]=sum;
	}

      //      for(i=1;i<seg->Steps;i++) { printf("Step=%d [%f %f], E=%f , a[i]=%f M(i)=%f\n",i,seg->Xtrack[i],seg->Ytrack[i],seg->Efield[i],Real.alpha(seg->Efield[i]),DIF[i]); }

     
      data[0]=DIF[seg->Steps-1]; data[1]=0; data[2]=0; data[3]=0;
      // printf("KKK %f\n",data[0]);
      for(i=1;i<seg->Steps;i++) 
	{
	  if((DIF[i]-1)/(data[0]-1)>0.8 && DIF[i]>1.02) 
	    {
	      data[1]+=seg->Xtrack[i]; data[2]+=seg->Ytrack[i]; 
	      data[3]+=seg->Time[i]; numreg++;
	      //printf("MUL=%f (X=%3.1f Y=%3.1f)::::: %e %e %e %d\n",DIF[i],seg->Xtrack[i],seg->Ytrack[i],data[1],data[2],data[3],numreg); 	      
	    }
	}
      if(numreg!=0) {data[1]/=numreg; data[2]/=numreg; data[3]/=numreg;}
    }






// void KPad::Alpha(Float_t E, Float_t Angle, Float_t *enp)
// {
//   // The simulation of the drift for the minimum ionizing particles. 
//   // A track is devided into Int_ div buckets. Each bucket is drifted in the field. The
//   // induced currents for each carrier is calculated as the sum  all buckets. 
//   //	Int_t MobMod; mobility model
//   //	Float_t B; magnetic field

// TH1F trn[4],trp[4];
// Double_t x[1000],y[1000];
// Float_t sp[3];
// int i,j,div;
// KStruct seg;
// TH1F *histop  = new TH1F("ch-","charge+",pos->GetNbinsX(),pos->GetXaxis()->GetXmin(),pos->GetXaxis()->GetXmax());
// TH1F *histon  = new TH1F("ch+","charge-",neg->GetNbinsX(),neg->GetXaxis()->GetXmin(),neg->GetXaxis()->GetXmax()); 
// sum->Reset(); pos->Reset(); neg->Reset();

// div=(Int_t) KMaterial::dEX(E,x,y,1);

// for(i=0;i<=div;i++) 
//   {
//     for(j=0;j<3;j++) {if(j==0) sp[j]=enp[j]+TMath::Sin(Angle)*x[i]; if(j==1) sp[j]=enp[j]-TMath::Cos(Angle)*x[i];} 
//     //if(+enp[j]+(exp[j]-enp[j])/(2*div);}
//     //printf("#i=%d div=%d, pointx=%f, pointy=%f weight=%f\n",i,div,sp[0],sp[1],y[i]); 
//     for(j=0;j<average;j++){ 
//       Drift(sp[0],sp[1],1,&seg,MobMod);
//       seg.GetCH(histop,1); 
//       Drift(sp[0],sp[1],-1,&seg,MobMod);
//       seg.GetCH(histon,1); 
//                      }
//     histop->Scale(y[i]/((Float_t)average)); histon->Scale(y[i]/((Float_t)average)); 
//     //    histop->Scale(0.05); histon->Scale(0.05); 
//     if(trapping) {
//        ht->Trapping(histop,trp);  
//        et->Trapping(histon,trn);
//        pos->Add(&trp[3]); neg->Add(&trn[3]);     
//     } else {pos->Add(histop);neg->Add(histon);}
//       histop->Reset();
//       histon->Reset();
  
//   }
// sum->Add(neg);
// sum->Add(pos);
// delete histop;
// delete histon;
// }

void KDetector::SetDriftHisto(Float_t x, Int_t numbins)
{
  if(x<0 || x>10000e-9) 
    printf("Selected time range out of scope !\n"); 
  else
    {
    if(pos!=NULL) delete pos;
	pos  = new TH1F("charge+","Positive Charge",numbins,0,x);
    if(neg!=NULL) delete neg;
        neg  = new TH1F("charge-","Negative Charge",numbins,0,x); 	
    if(sum!=NULL) delete sum;
	sum  = new TH1F("charge","Total Charge",numbins,0,x); 
    
    sum->SetXTitle("t [s]");  neg->SetXTitle("t [s]"); pos->SetXTitle("t [s]");
    sum->SetYTitle("I [arb.]");neg->SetYTitle("I [arb.]");pos->SetYTitle("I [arb.]");
    sum->GetYaxis()->SetTitleOffset(1.4); neg->GetYaxis()->SetTitleOffset(1.4); pos->GetYaxis()->SetTitleOffset(1.4);
    pos->SetLineColor(2);
    neg->SetLineColor(4);
    sum->SetLineColor(1);	

    }
}


/////////////////////////////////////////////////////////
void  KDetector::Save(Char_t *name,Char_t *file)
{
  Char_t str[100];
  TFile *fn=new TFile(file,"UPDATE");
  sprintf(str,"E_%s",name);
  if(Real.U!=NULL) Real.U->Write(str); else {printf("Histogram U (Electric) does not exist!\n"); return;}
  sprintf(str,"W_%s",name);
  if(Ramo.U!=NULL) Ramo.U->Write(str); else {printf("Histogram U (Ramo) does not exist!\n"); return;}
  sprintf(str,"G_%s",name);
  if(EG!=NULL) EG->Write(str);         else {printf("Histogram Geometry does not exist!\n"); return;}
  sprintf(str,"M_%s",name);
  if(DM!=NULL) DM->Write(str);         else {printf("Histogram Material does not exist!\n"); return;}
  fn->Close();
}

TFile *KDetector::Read(Char_t *name,Char_t *file)
{
  Char_t str[100];
  TFile *fn=new TFile(file);
  if(Real.U==NULL) Real.U=new TH3F();
  if(Ramo.U==NULL) Ramo.U=new TH3F();
  if(EG==NULL) EG=new TH3I();
  if(DM==NULL) DM=new TH3I();
  sprintf(str,"E_%s",name); Real.U->Read(str);
  sprintf(str,"W_%s",name); Ramo.U->Read(str);
  sprintf(str,"G_%s",name); EG->Read(str);
  sprintf(str,"M_%s",name); DM->Read(str);

   nx=EG->GetNbinsX();
   ny=EG->GetNbinsY();
   nz=EG->GetNbinsZ();

   Real.CalField();
   Ramo.CalField();
  

  return fn;
}



void KDetector::GaussBeam(Int_t div, Float_t Lambda, Float_t w0, Float_t CellX, Float_t B, Float_t lambda)  //user manca
{
  ////////////////////////////////////////////////////////////////////////////////////////
  // The simulation of the drift for Gaussian Beam.                                     //
  // A track is divided into Int_ div buckets. Each bucket is drifted in the field.     //
  // The induced currents for each carrier is calculated as the sum all buckets.        //
  //                                                                                    //
  //	- Int_t MobMod; mobility model                                                  //
  //	- Float_t B; magnetic field                                                     //
  //                                                                                    //
  //    - Lambda [um]                                                                   //
  //    - w0 [um]                                                                       //
  //                                                                                    //
  ////////////////////////////////////////////////////////////////////////////////////////

Double_t data[4];
Float_t sp[3];
Float_t Io,len=0,lent=0,lenlam,scalef=0,pscalef=0,mule,mulh;
Int_t i,j,k,e,l,Kint,KAint;
KStruct seg,segmul;
Float_t pi,w,X,xr,K,KA,PenDepth;

PenDepth=1014; //um

 pi=TMath::Pi();
xr=pi/(Lambda/3.54)*TMath::Power(w0,2);
X=enp[0]-CellX*.5;
w=w0*TMath::Power(1+TMath::Power(X/xr,2),0.5);
K=2*w/(CellX/div);
Kint=int(K);
KAint=Kint;
//printf("K=%f\n",K);

TH1F *histop  = new TH1F("ch-","charge+",pos->GetNbinsX(),pos->GetXaxis()->GetXmin(),pos->GetXaxis()->GetXmax()); //2.23 -> 3.0
TH1F *histon  = new TH1F("ch+","charge-",neg->GetNbinsX(),neg->GetXaxis()->GetXmin(),neg->GetXaxis()->GetXmax()); 
sum->Reset(); pos->Reset(); neg->Reset();

if(lambda!=0)
{
  for(k=0;k<3;k++) lent+=TMath::Power(exp[k]-enp[k],2);
  lent=TMath::Sqrt(lent);
  lenlam=lent/lambda;
  Io=div/(lambda*(1-TMath::Exp(-lenlam)));
  //     printf("Calculated length=%f um, len/lamda=%f, I0=%f\n",lent,lenlam,I0); 
}

////////////////////////////////////////////////////////////////////////////////////

for(i=0;i<div;i++) 
  for(l=0;l<=KAint;l++) 
    { 
    for(j=0;j<3;j++) {
      if ((KAint>0) && ((j==0) || (j==2)) ) {sp[j]=((exp[j]-enp[j])/div)*i+enp[j]+(exp[j]-enp[j])/(2*div);}
      else if ((KAint>0) && (j==1)) 
      {
        KA=K*TMath::Exp(-sp[0]/PenDepth);
        KAint=int(KA);
        X=sp[0]-CellX*.5;
        w=w0*TMath::Power(1+TMath::Power(X/xr,2),0.5);
        if (KAint>0) {sp[j]=((exp[j]-enp[j])/div)*i+enp[j]+(exp[j]-enp[j])/(2*div)-w+w/(.5*KAint)*l;}
        //printf("i=%i,j=%i,sp[j]=%f,w=%f,KA=%f\n",i,j,sp[j],w,KA);
      }
    else { KAint=-1; div=-1; }
    }

    for(j=0;j<average;j++)
    { 
      Drift(sp[0],sp[1],sp[2],1,&seg);
      seg.GetCH(histop,1,1,tauh); 
      Drift(sp[0],sp[1],sp[2],-1,&seg);
  
      if(MTresh>1)  
      mule=seg.GetCHMult(histon,1,1,taue);  // performing multiplication
      else          
      seg.GetCH(histon,1,1,taue); 
	
      if(mule>MTresh && MTresh>1) //if the multiplication is large enough then do the hole tracking
 	    {
        for(e=1;e<seg.Steps+1;e++)
	      {
	        Drift(seg.Xtrack[e],seg.Ytrack[e],seg.Ztrack[e],1,&segmul,seg.Time[e]); 
	        mulh=segmul.GetCHMult(histop,1,seg.MulCar[e],tauh);
	        if(mulh>BDTresh) {printf("HOLE MULTIPLICATION - BREAKDOWN\n"); BreakDown=1;}
        }
      }
    }
    
    histop->Scale(1/((Float_t)average)); histon->Scale(1/((Float_t)average)); 

    if(lambda!=0)
    {
      len=0;
      for(k=0;k<3;k++) len+=TMath::Power(sp[k]-enp[k],2); len=TMath::Sqrt(len); 
      scalef=Io*TMath::Exp(-len/lambda);
      //      printf("step=%d ::  len=%f um, scalef=%f pscalef=%f\n",i,len,scalef,pscalef); 
      histon->Scale(scalef); histop->Scale(scalef);
      pscalef+=scalef;
      histop->Scale(lent/div); histon->Scale(lent/div); 
    }
  
////////////////////////////////////////////////////////////////////////////////////

  pos->Add(histop);neg->Add(histon);
  histop->Reset();
  histon->Reset();
  }

sum->Add(neg);
sum->Add(pos);

delete histop;
delete histon;
}



void KDetector::ShowGaussBeam(Int_t div, Float_t Lambda, Float_t w0, Float_t CellX, Int_t color,Int_t how) //user manca
{
  ///////////////////////////////////////////////////////////////////////////////////////
  // Gaussian Beam for strip detector                                                  //
  //                                                                                   //
  // A track is divided into Int_ div buckets. Each bucket is drifted in the field.    //
  // The induced currents for each carrier is calculated as the sum all buckets.       //
  //                                                                                   //
  //    - Lambda [um]                                                                  //
  //    - w0 [um]                                                                      //
  ///////////////////////////////////////////////////////////////////////////////////////

TGraph *gr;
TPolyLine3D *gr3D;
Float_t sp[3];
int i,j,Kint,KAint;
KStruct seg;
TLine *line;
Float_t pi,w,X,xr,l,K,KA,PenDepth;

PenDepth=1014;

 pi=TMath::Pi();
xr=pi/(Lambda/3.54)*TMath::Power(w0,2);
X=enp[0]-CellX*.5;
w=w0*TMath::Power(1+TMath::Power(X/xr,2),0.5);
K=2*w/(CellX/div);
Kint=int(K);
KAint=Kint;
//printf("K=%f\n",K);


// Draw histograms 

if(EG!=NULL) 
{ 
  if(nz==1) KHisProject(EG,3,how)->Draw("COL");
  else
  {
    TH3F *hh=GetGeom();
    hh->SetFillColor(color);
    hh->Draw("iso");
  }
} 

// Draw drift paths

for(i=0;i<div;i++) 
  for(l=0;l<=KAint;l++) 
    { 
    for(j=0;j<3;j++) {
      if ((KAint>0) && ((j==0) || (j==2)) ) {sp[j]=((exp[j]-enp[j])/div)*i+enp[j]+(exp[j]-enp[j])/(2*div);}
      else if ((KAint>0) && (j==1)) 
      {
        KA=K*TMath::Exp(-sp[0]/PenDepth);
        KAint=int(KA);
        X=sp[0]-CellX*.5;
        w=w0*TMath::Power(1+TMath::Power(X/xr,2),0.5);
        if (KAint>0) {sp[j]=((exp[j]-enp[j])/div)*i+enp[j]+(exp[j]-enp[j])/(2*div)-w+w/(.5*KAint)*l;}
       // printf("i=%i,j=%i,sp[j]=%f,w=%f,KA=%f\n",i,j,sp[j],w,KA);
      }
    else { KAint=-1; div=-1; }
    }

    // electron drift
    if(Debug)   printf("Entry Point: %f %f %f \n",sp[0],sp[1],sp[2]);
    Drift(sp[0],sp[1],sp[2],1,&seg);

    if(nz==1)
    gr=new TGraph(seg.Steps,&seg.Xtrack[1],&seg.Ytrack[1]); 
    else
    gr3D=new TPolyLine3D(seg.Steps,&seg.Xtrack[1],&seg.Ytrack[1],&seg.Ztrack[1]); 	
    if(nz==1)
    {
      gr->SetLineColor(2);
      gr->SetLineStyle(1);
      gr->Draw("L");
    }
    else
    {
      gr3D->SetLineStyle(1);  
      gr3D->SetLineColor(2); 
      gr3D->Draw("SAME"); 
    }

    // hole drift

    Drift(sp[0],sp[1],sp[2],-1,&seg);

    if(nz==1)
    gr=new TGraph(seg.Steps,&seg.Xtrack[1],&seg.Ytrack[1]); 
    else
    gr3D=new TPolyLine3D(seg.Steps,&seg.Xtrack[1],&seg.Ytrack[1],&seg.Ztrack[1]); 	
    if(nz==1)
    {
      gr->SetLineColor(4);
      gr->SetLineStyle(3); 
      gr->Draw("L");
    }
    else
    {
      gr3D->SetLineStyle(1);  
      gr3D->SetLineColor(4); 
      gr3D->Draw("SAME"); 
    }
  }
}



