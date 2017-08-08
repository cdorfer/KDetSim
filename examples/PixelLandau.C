{
  KPixel *det=new KPixel(5,50,120,100);
  det->Voltage=60;

  TF3 *f2=new TF3("f2","x[0]*x[1]*x[2]*0+[0]",0,3000,0,3000,0,3000);
  f2->SetParameter(0,-2);
  det->NeffF=f2;

  det->SetUpVolume(2,2,1);
  det->SetUpPixel(0,25, 60,10,10,2,16385);
  det->SetUpPixel(1,10, 10,10,10,2,1);
  det->SetUpPixel(2,40, 10,10,10,2,1);
  det->SetUpPixel(3,40,110,10,10,2,1);
  det->SetUpPixel(4,10,110,10,10,2,1);

  det->SetUpElectrodes();
  det->SetBoundaryConditions();
     
  det->CalField(0);
  det->CalField(1);

  // track entry and exit point  
  det->enp[0]=25;
  det->enp[1]=60;
  det->enp[2]=1;
//det->exp[0]=25;
  det->exp[0]=25;
  det->exp[1]=60;
  det->exp[2]=100;
  ///////////////////////////// START OF THE LANDAU FLUCT.... //////////////////////
  //// define the number of segments to which you cur your track
  int seg=5;

  // define the length of the track through your sensors 
  float dist=TMath::Sqrt(TMath::Power(det->enp[0]-det->exp[0],2)+
                         TMath::Power(det->enp[1]-det->exp[1],2)+
			 TMath::Power(det->enp[2]-det->exp[2],2));
  // calculate most probable charge per segement in ke
  float lanseg=dist/seg*75/1000; 

  // define the distrubution of of energy loss - can be anything
  TF1 *lan=new TF1("lan","TMath::Landau(x,[0],[1])",0,100);
  lan->SetParameter(0,lanseg);
  lan->SetParameter(1,lanseg/10.);
  
  // book histograms
  TH1F *hisp[10]; 
  TH1F *hisn[10];
  TH1F *hiss[10]; 
  
  // define the track segments after (seg for the whole track)
  Float_t *entx=new Float_t [seg];   Float_t *extx=new Float_t [seg];
  Float_t *enty=new Float_t [seg];   Float_t *exty=new Float_t [seg];
  Float_t *entz=new Float_t [seg];   Float_t *extz=new Float_t [seg];

  // vector of the track line
  Float_t k[3]={((float) det->exp[0]-det->enp[0])/seg,((float)det->exp[1]-det->enp[1])/seg,((float)det->exp[2]-det->enp[2])/seg};

  // book the histogram for the segmets and entry points
  for(int i=0;i<seg;i++)
    {
      entx[i]=k[0]*i+det->enp[0]; extx[i]=k[0]*(i+1)+det->enp[0];
      enty[i]=k[1]*i+det->enp[1]; exty[i]=k[1]*(i+1)+det->enp[1];
      entz[i]=k[2]*i+det->enp[2]; extz[i]=k[2]*(i+1)+det->enp[2];
      hisp[i]=new TH1F(); 
      hiss[i]=new TH1F(); 
      hisn[i]=new TH1F();
      printf("%d %f-%f %f-%f %f-%f\n",i,entx[i],extx[i],enty[i],exty[i], entz[i],extz[i]);
    }

  // calculated the induced currents for the segments
  for(int i=0;i<seg;i++)
    {
     det->sum->Reset();
      det->pos->Reset();
      det->neg->Reset();
 
      
     det->SetEntryPoint(entx[i],enty[i],entz[i]);
      det->SetEntryPoint(extx[i],exty[i],extz[i]);
       det->MipIR(dist/seg);
      det->sum->Copy(*hiss[i]);
      det->pos->Copy(*hisp[i]);
      det->neg->Copy(*hisn[i]);
    }

  // histograms for the individual events
  TH1F *hist[100];
  TH1F *htemp=new TH1F();

  // make a loop over 20 events
 for(int j=0;j<20;j++)
  {
    // book the histogram
    hist[j]=new TH1F("landau","landau",200,0,25e-9); 
    
    // multiply the currents at segments with random number according to landau function 
    for(i=0;i<seg;i++)
      hist[j]->Add(hiss[i],lan->GetRandom()*1000/dist*seg);
      
    // draw everything
   if(j==0) 
       hist[j]->Draw(); 
   else  
     {
     hist[j]->SetLineColor(j+1);
     hist[j]->Draw("SAME");
     }

  printf("%d %f\n",j,hist[j]->Integral());


  }
}


