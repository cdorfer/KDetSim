{
  TF1 *neff=new TF1("neff","[0]+x[0]*0",0,1000);
  neff->SetParameter(0,0.2);
  KPad det(50,300);
  det->Neff=neff;
  det->Voltage=-200;
  det->SetUpVolume(1);
  det->SetUpElectrodes();  
}
