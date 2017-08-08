{
  KPixel *det=new KPixel(5,50,120,100);
  det->Voltage=200;

  det->SetUpVolume(2,2,1);
  det->SetUpPixel(0,25,60,10,10,2,16385);
  det->SetUpPixel(1,10,10,10,10,2,1);
  det->SetUpPixel(2,40,10,10,10,2,1);
  det->SetUpPixel(3,40,110,10,10,2,1);
  det->SetUpPixel(4,10,110,10,10,2,1);


  det->SetUpElectrodes();
  det->SetBoundaryConditions();

   TF3 *f2=new TF3("f2","x[0]*x[1]*x[2]*0+[0]",0,3000,0,3000,0,3000);
   f2->SetParameter(0,-2);
   det->NeffF=f2;

   
   det->CalField(0);
   det->CalField(1);
   
   det->enp[0]=25;
   det->enp[1]=60;
   det->enp[2]=1;
   det->exp[0]=25;
   det->exp[1]=60;
   det->exp[2]=100;
}
