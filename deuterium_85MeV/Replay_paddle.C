using namespace std;

void Replay_paddle(){
  gStyle->SetOptStat(0);
  const Double_t angles[4] = {40.0, 55.0, 75.0, 90.0};


  const Double_t kY1 = 0.0; //MeV
  const Double_t kX1 = 0.0; //0.001*Q_i
  const Double_t E_beam = 85.0; //MeV, incident photon energy
  const Double_t M_D = 1876.1; //MeV, deuteron mass
  Double_t kY2[4];
  for(int i=0; i<4; i++){
    kY2[i] = E_beam/(1+(E_beam/M_D)*(1-cos(angles[i]*TMath::Pi()/180.0)));
  }
  const Double_t kX2[4] = {4300, 6700, 4000, 1700};
  Double_t m_cal_[4];
  Double_t b_cal_[4];
  for(int i=0; i<4; i++) {
    m_cal_[i] = (kY2[i] - kY1) / (kX2[i] - kX1);
    b_cal_[i] = kY2[i] - m_cal_[i] * kX2[i];
  }



  TH2F *hpad2D[4];
  for(int i=0; i<4; i++){
    hpad2D[i] = new TH2F(Form("%.0f#circ",angles[i]),
			 Form("%.0f#circ",angles[i]),
			 1000,0,200, 25000,0,10000);
  }





  ifstream database("run_database.dat");
  char number[20];
  char name[200];
  string mode;
  string hash = "#";
  float Qcore[4];
  float Qpaddle[4];
  float energy;

  while(database >> number >> name >> mode){
    if(number[0] == hash) continue; //the runs starting with # should be neglected 
    if(mode != "deuterium") continue;//choose run mode: deuterium/empty
     // if(atoi(number) < 817) continue;//start from this run
     if(atoi(number) > 600) break;//end at this run
    printf("Current run number: %s  mode: %s\n", number, mode.c_str());

    TFile *file = new TFile(Form("./root_files/%s",name));
    TTree *tree = (TTree*)file->Get("flat_tree");

    tree->SetBranchAddress("Q_0", &Qcore[0]);// ALAINA (40 deg.)
    tree->SetBranchAddress("Q_7", &Qcore[1]);// SUSAN (55 deg.)
    tree->SetBranchAddress("Q_3", &Qcore[2]);// JONI (75 deg.)
    tree->SetBranchAddress("Q_6", &Qcore[3]);// ROBERTA (90 deg.)
    tree->SetBranchAddress("Q_16", &Qpaddle[0]);// AMY (40 deg.)
    tree->SetBranchAddress("Q_23", &Qpaddle[1]);// MARY (55 deg.)
    tree->SetBranchAddress("Q_19", &Qpaddle[2]);// SAMANTHA (75 deg.)
    tree->SetBranchAddress("Q_22", &Qpaddle[3]);// EVA (90 deg.)
  
    int events = tree->GetEntries();
    for(int i=0; i<events; i++){
      tree->GetEntry(i);

      for(int j=0; j<4; j++){
	energy = Qcore[j]*m_cal_[j]+b_cal_[j]; //energy calibration
	hpad2D[j]->Fill(energy,Qpaddle[j]);
      }

    }//loop over events(i<events)
    file->Close();

  }//loop over while full/empty runs
  database.close();



  TCanvas *c1 = new TCanvas("c1","paddle");
  c1->Divide(2,2);
  for(int i=0; i<4; i++){
    c1->cd(i+1);
    gPad->SetLogz();
    hpad2D[i]->GetXaxis()->SetTitle("E [MeV]");
    hpad2D[i]->GetXaxis()->SetTitleSize(0.05);
    hpad2D[i]->GetXaxis()->SetLabelSize(0.05);
    hpad2D[i]->GetYaxis()->SetTitle("#DeltaE");
    hpad2D[i]->Draw("color");
  }



  // TFile *temp = new TFile("test_paddle.root","recreate");
  // for(int i=0; i<4; i++){
  //   hpad2D[i]->Write();
  // }
  // c1->Write();


}
