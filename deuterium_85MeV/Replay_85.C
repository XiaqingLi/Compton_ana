using namespace std;

void Replay_85(){

  const int nCores = 8;
  const int nShields = 8;
  const string crNames[nCores] = {"ALAINA", "BROOKE", "CINDY", "JONI",
				  "KRISTA", "LINDA", "ROBERTA", "SUSAN"};
  const string shldNames[nShields] = {"ANDY", "BOBBY", "CARL", "JACK",
				      "KEVIN", "LARRY", "REGGIE", "STEVE"};
  const Double_t angles[nCores] = {40.0, 159.0, 125.0, 75.0, 
				   145.0, 110.0, 90.0, 55.0};



  //============== calibration ==========================
  //             (from M.Sikora)
  const Double_t kY1 = 0.0; //MeV
  const Double_t kX1 = 0.0; //0.001*Q_i
  const Double_t E_beam = 85.0; //MeV, incident photon energy
  const Double_t M_D = 1876.1; //MeV, deuteron mass
  Double_t kY2[nCores];
  for(int i=0; i<nCores; i++){
    kY2[i] = E_beam/(1+(E_beam/M_D)*(1-cos(angles[i]*TMath::Pi()/180.0)));
  }

  const Double_t kX2[nCores] = {4400, 3900, 3400, 3900,
  				2300, 2800, 1700, 6800};// Q_i
  //const Double_t kX2[nCores] = {4400, 4000, 3400, 4000,
  //				2300, 2900, 1700, 6800};// Q_i
  // const Double_t kX2[nCores] = {4.4, 3.9, 3.4, 3.9,
  // 				 2.3, 2.8, 1.7, 6.8};// 0.001*Q_i
  Double_t m_cal_[nCores];
  Double_t b_cal_[nCores];
  for(int i=0; i<nCores; i++) {
    m_cal_[i] = (kY2[i] - kY1) / (kX2[i] - kX1);
    b_cal_[i] = kY2[i] - m_cal_[i] * kX2[i];
  }

  const Double_t cutLow[nCores] = {64.0, 44.0, 74.0, 74.0,
				   71.0, 61.0, 69.0, 65.0}; //ns
  const Double_t cutHigh[nCores] = {69.0, 49.0, 79.0, 79.0,
				    76.0, 65.0, 74.0, 70.0}; //ns
  const Double_t shldCut[nCores] = {4000.0, 3500.0, 3500.0, 3000.0,
				    3500.0, 3000.0, 3000.0, 3500.0}; //Q_i
  const Double_t crCut[nCores] = {9.5, 11.0, 13.0, 9.5,
				  23.0, 14.0, 19.5, 11.0}; //MeV
  const Double_t randLow = 130.; //MeV
  const Double_t randHigh = 170.; //MeV

  gStyle->SetOptStat(0);

  TH1F *hcore_raw[nCores];
  TH1F *hcore[nCores];
  TH1F *htof[nCores];
  TH2F *h2D[nCores];
  TH1F *hcore_tof[nCores];
  TH1F *hcore_shld[nCores];
  TH1F *hcore_tof_shld[nCores];
  TH1F *hcore_tof_shldR[nCores];
  TH1F *hcore_net[nCores];
  TH1F *hshield_raw[nShields];
  TH1F *hrand[nCores];
  TH2F *hpad2D[4];
  TH1F *hpad[4];
  TH2F *hpad2DR[4];
  TH1F *hpadR[4];

  for(int i=0; i<nCores; i++){ 
    hcore_raw[i] = new TH1F(Form("%s_coreRaw",crNames[i].c_str()),
			    Form("%s_%.0f#circ",crNames[i].c_str(),angles[i]),
			    2000,0,10000);
    hcore_raw[i]->Sumw2();
    hcore[i] = new TH1F(Form("%s_core",crNames[i].c_str()),
			    Form("%s_%.0f#circ",crNames[i].c_str(),angles[i]),
			    2000,0,200);
    hcore[i]->Sumw2();
    htof[i] = new TH1F(Form("%s_tof",crNames[i].c_str()),
		       Form("%s_%.0f#circ",crNames[i].c_str(),angles[i]),
		       1000,0,200);
    hcore_tof[i] = new TH1F(Form("%s_TOFCut",crNames[i].c_str()),
			    Form("%s_%.0f#circ",crNames[i].c_str(),angles[i]),
			    2000,0,200);
    hcore_tof[i]->Sumw2();
    hcore_shld[i] = new TH1F(Form("%s_ShieldCut",crNames[i].c_str()),
			     Form("%s_%.0f#circ",crNames[i].c_str(),angles[i]),
			     2000,0,200);
    hcore_shld[i]->Sumw2();
    hcore_tof_shld[i] = new TH1F(Form("%s_TOF+ShldCut",crNames[i].c_str()),
				 Form("%s_%.0f#circ",crNames[i].c_str(),angles[i]),
				 2000,0,200);
    hcore_tof_shld[i]->Sumw2();
    hcore_tof_shldR[i] = new TH1F(Form("%s_TOF+ShldCut_Raw",crNames[i].c_str()),
				  Form("%s_%.0f#circ",crNames[i].c_str(),angles[i]),
				  2000,0,10000);
    hcore_net[i] = new TH1F(Form("%s_Net",crNames[i].c_str()),
			    Form("%s_%.0f#circ",crNames[i].c_str(),angles[i]),
			    2000,0,200);
    hcore_net[i]->Sumw2();
    hrand[i] = new TH1F(Form("%s_Rand",crNames[i].c_str()),
			Form("%s_%.0f#circ",crNames[i].c_str(),angles[i]),
			2000,0,200); 
    hrand[i]->Sumw2();
    h2D[i] = new TH2F(Form("%s_2D",crNames[i].c_str()),
		      Form("%s_%.0f#circ",crNames[i].c_str(),angles[i]),
		      1000,0,200, 2000,0,200);
    hshield_raw[i] = new TH1F(Form("%s_shieldRaw",crNames[i].c_str()),
			      Form("%s_%.0f#circ",shldNames[i].c_str(),angles[i]),
			      2000,0,20000);
  }

  Double_t PadAngles[4] = {40, 55, 75, 90};
  for(int i=0; i<4; i++){
    hpad2D[i] = new TH2F(Form("%.0f#circ",PadAngles[i]),
			 Form("%.0f#circ",PadAngles[i]),
			 1000,0,200, 2500,0,500);
    hpad[i] = new TH1F(Form("%.0f#circ",PadAngles[i]),
		       Form("%.0f#circ",PadAngles[i]),
		       2500,0,500);
    hpad[i]->Sumw2();
    hpad2DR[i] = new TH2F(Form("%.0f#circ_no_cut",PadAngles[i]),
			  Form("%.0f#circ_no_cut",PadAngles[i]),
			  1000,0,200, 2500,0,500);
    hpadR[i] = new TH1F(Form("%.0f#circ_no_cut",PadAngles[i]),
			Form("%.0f#circ_no_cut",PadAngles[i]),
			2500,0,500);
    hpadR[i]->Sumw2();
  }



  ifstream database("run_database.dat");
  char number[20];
  char name[200];
  string mode;
  string hash = "#";
  float Qcore[nCores];
  float Qshield[nCores];
  float Qpaddle[4];
  float tof_raw, tof, energy;
  while(database >> number >> name >> mode){
    if(number[0] == hash) continue; //the runs starting with # should be neglected 
    if(mode != "deuterium") continue;//choose run mode
     // if(atoi(number) < 817) continue;//start from this run
     // if(atoi(number) > 817) break;//end at this run
    printf("Current run number: %s  mode: %s\n", number, mode.c_str());

    TFile *file = new TFile(Form("./root_files/%s",name));
    TTree *tree = (TTree*)file->Get("flat_tree");
    
    tree->SetBranchAddress("Q_0", &Qcore[0]);// ALAINA
    tree->SetBranchAddress("Q_1", &Qcore[1]);// BROOKE
    tree->SetBranchAddress("Q_2", &Qcore[2]);// CINDY
    tree->SetBranchAddress("Q_3", &Qcore[3]);// JONI
    tree->SetBranchAddress("Q_4", &Qcore[4]);// KRISTA
    tree->SetBranchAddress("Q_5", &Qcore[5]);// LINDA
    tree->SetBranchAddress("Q_6", &Qcore[6]);// ROBERTA
    tree->SetBranchAddress("Q_7", &Qcore[7]);// SUSAN
    tree->SetBranchAddress("Q_8", &Qshield[0]);// ANDY
    tree->SetBranchAddress("Q_9", &Qshield[1]);// BOBBY
    tree->SetBranchAddress("Q_10", &Qshield[2]);// CARL
    tree->SetBranchAddress("Q_11", &Qshield[3]);// JACK
    tree->SetBranchAddress("Q_12", &Qshield[4]);// KEVIN
    tree->SetBranchAddress("Q_13", &Qshield[5]);// LARRY
    tree->SetBranchAddress("Q_31", &Qshield[6]);// REGGIE (digitizer B)
    tree->SetBranchAddress("Q_14", &tof_raw);
    tree->SetBranchAddress("Q_15", &Qshield[7]);// STEVE
    tree->SetBranchAddress("Q_16", &Qpaddle[0]);// AMY (40 deg.)
    tree->SetBranchAddress("Q_19", &Qpaddle[1]);// SAMANTHA (75 deg.)
    tree->SetBranchAddress("Q_22", &Qpaddle[2]);// EVA (90 deg.)
    tree->SetBranchAddress("Q_23", &Qpaddle[3]);// MARY (55 deg.)
  
 

    int events = tree->GetEntries();
    for(int i=0; i<events; i++){
      tree->GetEntry(i);
      tof = 0.025*tof_raw*0.5739-9.23582;
     
      for(int j=0; j<nCores; j++){
	energy = Qcore[j]*m_cal_[j]+b_cal_[j]; //energy calibration
	hcore_raw[j]->Fill(Qcore[j]); //core uncalibrated
	hcore[j]->Fill(energy);	//core calibrated (no cut)

	//core thershold cut
      	if(energy>crCut[j]){
      	  htof[j]->Fill(tof); //tof hist
      	  h2D[j]->Fill(tof,energy); //tof vs core energy

      	  //tof cut
      	  if(tof>cutLow[j] && tof<cutHigh[j])
      	    hcore_tof[j]->Fill(energy);
	 
      	  //shield cut
      	  if(Qshield[j]<shldCut[j])
      	    hcore_shld[j]->Fill(energy);

      	  //tof and shield cuts
      	  if(tof>cutLow[j] && tof<cutHigh[j] && Qshield[j]<shldCut[j]){
      	    hcore_tof_shld[j]->Fill(energy);
      	    hcore_net[j]->Fill(energy); //for random subtraction
	    hcore_tof_shldR[j]->Fill(Qcore[j]); //to find kX2
      	  }

      	  //random events with shield cut
      	  if(tof>randLow && tof<randHigh && Qshield[j]<shldCut[j])
      	    hrand[j]->Fill(energy);

	  //veto paddle with cut
	  if(tof>cutLow[j] && tof<cutHigh[j] && Qshield[j]<shldCut[j]){
	    if(angles[j]==40){
	      hpad2D[0]->Fill(energy,Qpaddle[0]);
	      hpad[0]->Fill(Qpaddle[0]);
	    }
	    else if(angles[j]==55){
	      hpad2D[1]->Fill(energy,Qpaddle[1]);
	      hpad[1]->Fill(Qpaddle[1]);
	    }
	    else if(angles[j]==75){
	      hpad2D[2]->Fill(energy,Qpaddle[2]);
	      hpad[2]->Fill(Qpaddle[2]);
	    }
	    else if(angles[j]==90){
	      hpad2D[3]->Fill(energy,Qpaddle[3]);
	      hpad[3]->Fill(Qpaddle[3]);
	    }
	  }

	  //all paddle events(no shield or tof cut)
	  if(angles[j]==40){
	    hpad2DR[0]->Fill(energy,Qpaddle[0]);
	    hpadR[0]->Fill(Qpaddle[0]);
	  }
	  else if(angles[j]==55){
	    hpad2DR[1]->Fill(energy,Qpaddle[1]);
	    hpadR[1]->Fill(Qpaddle[1]);
	  }
	  else if(angles[j]==75){
	    hpad2DR[2]->Fill(energy,Qpaddle[2]);
	    hpadR[2]->Fill(Qpaddle[2]);
	  }
	  else if(angles[j]==90){
	    hpad2DR[3]->Fill(energy,Qpaddle[3]);
	    hpadR[3]->Fill(Qpaddle[3]);
	  }


      	}//if thrsh cut

      }//loop over cores(j<nCores)

      //shield energy uncalibrated
      for(int j=0; j<nShields; j++)
      	hshield_raw[j]->Fill(Qshield[j]);
 
 
    }//loop over events(i<events)
    file->Close();


  }//loop over while full/empty runs
  database.close();


  //random subtraction
  for(int i=0; i<nCores; i++)
    hcore_net[i]->Add(hrand[i],-1.*(cutHigh[i]-cutLow[i])/(randHigh-randLow));
  


   
  //plot hists to canvas
  TCanvas *cCore = new TCanvas("cCore","Core Energy Deposition");
  TCanvas *cCoreRaw = new TCanvas("cCoreRaw","Core Energy Deposition (uncalibrated");
  TCanvas *c2D = new TCanvas("c2D","Energy vs. Tof");
  TCanvas *cShield = new TCanvas("cShield","Shield Energy Deposition");
  //TCanvas *cRand = new TCanvas("cRand","Random events with shield cut");
  cCore->Divide(4,2);
  cCoreRaw->Divide(4,2);
  c2D->Divide(4,2);
  cShield->Divide(4,2);
  //cRand->Divide(4,2);

  for(int i=0; i<nCores; i++){
    cCoreRaw->cd(i+1);
    hcore_raw[i]->Draw();
    cCore->cd(i+1);
    //hcore_net[i]->GetYaxis()->SetRangeUser(-5,70);
    hcore_net[i]->GetXaxis()->SetTitle("Energy [MeV]");
    hcore_net[i]->GetXaxis()->SetTitleSize(0.05);
    hcore_net[i]->GetYaxis()->SetLabelSize(0.045);
    hcore_net[i]->GetXaxis()->SetLabelSize(0.045);
    hcore_net[i]->SetLineColor(kViolet);
    hcore_net[i]->Draw();
    hcore[i]->Draw("same");
    hcore_tof[i]->SetLineColor(kGreen+1);
    hcore_tof[i]->Draw("same");
    hcore_shld[i]->SetLineColor(kAzure+1);
    hcore_shld[i]->Draw("same");
    hcore_tof_shld[i]->SetLineColor(kOrange+1);
    hcore_tof_shld[i]->Draw("same");

    c2D->cd(i+1);
    gPad->SetLogz();
    h2D[i]->GetXaxis()->SetTitle("TOF [ns]");
    h2D[i]->GetXaxis()->SetTitleSize(0.05);
    h2D[i]->GetYaxis()->SetTitle("Energy [MeV]");
    h2D[i]->GetYaxis()->SetTitleSize(0.05);
    h2D[i]->GetYaxis()->SetLabelSize(0.045);
    h2D[i]->GetXaxis()->SetLabelSize(0.045);
    h2D[i]->Draw("col");

    // cRand->cd(i+1);
    // hrand[i]->Draw();

    cShield->cd(i+1);
    //hshield_raw[i]->GetYaxis()->SetRangeUser(0,70);
    hshield_raw[i]->Draw();
  }
  

  //paddle plot
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
  TCanvas *c3 = new TCanvas("c3","paddle (no cut)");
  c3->Divide(2,2);
  for(int i=0; i<4; i++){
    c3->cd(i+1);
    gPad->SetLogz();
    hpad2DR[i]->GetXaxis()->SetTitle("E [MeV]");
    hpad2DR[i]->GetXaxis()->SetTitleSize(0.05);
    hpad2DR[i]->GetXaxis()->SetLabelSize(0.05);
    hpad2DR[i]->GetYaxis()->SetTitle("#DeltaE");
    hpad2DR[i]->Draw("color");
  }
  TCanvas *c2 = new TCanvas("c2","c2");
  c2->Divide(2,2);
  for(int i=0; i<4; i++){
    c2->cd(i+1);
    hpad[i]->GetXaxis()->SetTitle("E [MeV]");
    hpad[i]->GetXaxis()->SetTitleSize(0.05);
    hpad[i]->GetXaxis()->SetLabelSize(0.05);
    hpad[i]->Draw();
    hpadR[i]->SetLineColor(kAzure+1);
    hpadR[i]->Draw("same");
  }
  



  //write and save output files
  TFile *xlfile = new TFile("full_new.root","recreate");
  //TFile *xlfile = new TFile("run_180.root","recreate");
  //TFile *xlfile = new TFile("empty_week1.5.root","recreate");
  cCoreRaw->Write();
  cCore->Write();
  c2D->Write();
  //cRand->Write();
  cShield->Write();  

  c1->Write();
  c2->Write();
  c3->Write();

  for(int i=0; i<nCores; i++){
    hcore_raw[i]->Write();
    hcore[i]->Write();
    hcore_tof[i]->Write();
    hcore_shld[i]->Write();
    hcore_tof_shld[i]->Write();
    hcore_tof_shldR[i]->Write(); 
    hrand[i]->Write();
    hcore_net[i]->Write();
    htof[i]->Write();
    hshield_raw[i]->Write();
    h2D[i]->Write();
  }

  for(int i=0; i<4; i++){
    hpad[i]->Write();
    hpad2D[i]->Write();
    hpadR[i]->Write();
    hpad2DR[i]->Write();
  }

  xlfile->Close();
  printf("file written\n");
  
  

  // TFile *temp = new TFile("test_paddle.root","recreate");
  // for(int i=0; i<4; i++){
  //   hpad[i]->Write();
  //   hpad2D[i]->Write();
  //   hpadR[i]->Write();
  //   hpad2DR[i]->Write();
  // }
  // c1->Write();
  // c2->Write();

}
