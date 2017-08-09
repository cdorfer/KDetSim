{
KDetector det;
det.Voltage = 500;
det.nx = 225;
det.ny = 50;
det.nz = 225;

Float_t dimX = 450;
Float_t dimY = 50;
Float_t dimZ = 22;


//  Electrode Geometry
det.EG = new TH3I("EG","EG",det.nx,0,dimX, det.ny,0,dimY, det.nz,0,dimZ); 
det.EG->GetXaxis()->SetTitle("x [#mum]"); 
det.EG->GetYaxis()->SetTitle("y [#mum]"); 
det.EG->GetZaxis()->SetTitle("z [#mum]");
//  Detector Material
det.DM = new TH3I("DM","DM",det.nx,0,dimX, det.ny,0,dimY, det.nz,0,dimZ);
det.DM->GetXaxis()->SetTitle("x [#mum]"); 
det.DM->GetYaxis()->SetTitle("y [#mum]"); 
det.DM->GetZaxis()->SetTitle("z [#mum]");
//  Space charge
det.NeffH = new TH3F("Neff","Neff",det.nx,0,dimX, det.ny,0,dimY, det.nz,0,dimZ);
det.NeffH->GetXaxis()->SetTitle("x [#mum]");
det.NeffH->GetYaxis()->SetTitle("y [#mum]");
det.NeffH->GetZaxis()->SetTitle("z [#mum]");

//  Electrodes
Int_t weighting_field = 1;
Int_t material = 10;
//  GND Pad
Float_t GndPos[3]={225, 0, 10}; 
Float_t GndSiz[3]={400,0.1,18}; 
det.ElRectangle(GndPos, GndSiz, weighting_field, material);
//  Top Pad
Float_t PadPos[3]={225, 0, 10};
Float_t PadSiz[3]={400,0.1,18};
det.ElRectangle(PadPos, PadSiz, weighting_field, material);

//  Material Setup
Int_t i, j;
for(int j=0;j<=det.ny;j++)
  for(int i=0;i<=det.nx;i++) 
      det.DM->SetBinContent(i,j,1,0); //change here for diamond?
      det.NeffH->SetBinContent(i,j,1,0.1); //change here for diamond?


det.Voltage = - 400;
det.Voltage2 = 0;

det.SetBoundaryConditions();




det.CalField(0); //calculate weighting field
det.CalField(1); //calculate electric field
det.Draw("EFxz")->Draw("COLZ");

}