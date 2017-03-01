using namespace std;

void Replay_mixcut_hinda_paddle(){
  //  gStyle->SetOptStat(0);

  // const string crNames = "ALAINA";
  // const Double_t angles = 40.0;
  // const Double_t kX2 = 4300;
  const string crNames = "JONI";
  const Double_t angles = 75.0;
  const Double_t kX2 = 4000;

  //const string crNames[nCores] = {"ALAINA", "SUSAN", "JONI", "ROBERTA"};
  //const Double_t kX2[4] = {4300, 6700, 4000, 1700};


  //============== for calibration ==========================
  const Double_t E_beam = 85; //MeV
  const Double_t M_D = 1876.1; //MeV, deuteron mass
  const Double_t kX1 = 0.0; //0.001*Q_0
  const Double_t kY1 = 0.0; //MeV
  Double_t kY2 = E_beam/(1+(E_beam/M_D)*(1-cos(angles*TMath::Pi()/180.0)));
  Double_t a = (kY2 - kY1) / (kX2 - kX1);
  Double_t b = kY2 - a * kX2;




  // ====================== make histograms ================================
  TH2F *hEVpadtof = new TH2F(Form("%s_%.0f_EvPadTof",crNames.c_str(),angles),
			     Form("%s (%.0f#circ) shield cut",crNames.c_str(),angles),
			     1800,0,180, 2000,0,200);
  TH2F *hdEVpadtof1 = new TH2F(Form("%s_%.0f_dEvPadTof1",crNames.c_str(),angles),
			       Form("%s (%.0f#circ) shield+tof cut",crNames.c_str(),angles),
			       1000,0,1000, 4000,0,800);
  TH2F *hdEVpadtof2 = new TH2F(Form("%s_%.0f_dEvPadTof2",crNames.c_str(),angles),
			       Form("%s (%.0f#circ) shield+tof cut, Ecore>100MeV",crNames.c_str(),angles),
			       1000,0,1000, 4000,0,800);
  TH2F *hdEVpadtof3 = new TH2F(Form("%s_%.0f_dEvPadTof3",crNames.c_str(),angles),
			       Form("%s (%.0f#circ) shield+tof cut, Ecore=[80,100]MeV",crNames.c_str(),angles),
			       1000,0,1000, 4000,0,800);
  // TH2F *hdEVpadtof4 = new TH2F(Form("%s_%.0f_dEvPadTof4",crNames.c_str(),angles),
  // 			       Form("%s (%.0f#circ) shield+tof cut, Ecore=[20,60]MeV",crNames.c_str(),angles),
  // 			       1000,0,1000, 4000,0,800);
  TH2F *hdEVpadtof1_rndm = new TH2F(Form("%s_%.0f_dEvPadTof1_rndm",crNames.c_str(),angles),
				    Form("%s (%.0f#circ) shield cut, random events",crNames.c_str(),angles),
				    1000,0,1000, 4000,0,800);
  TH2F *hdEVpadtof2_rndm = new TH2F(Form("%s_%.0f_dEvPadTof2_rndm",crNames.c_str(),angles),
				    Form("%s (%.0f#circ) shield cut, random events, Ecore>100MeV",crNames.c_str(),angles),
				    1000,0,1000, 4000,0,800);
  TH2F *hdEVpadtof3_rndm = new TH2F(Form("%s_%.0f_dEvPadTof3_rndm",crNames.c_str(),angles),
				    Form("%s (%.0f#circ) shield cut, random events, Ecore=[80,100]MeV",crNames.c_str(),angles),
				    1000,0,1000, 4000,0,800);

  TH2F *hdEVpadtof1_net,*hdEVpadtof2_net,*hdEVpadtof3_net;

  TH1F *hdE = new TH1F(Form("%s_%.0f_dE",crNames.c_str(),angles),
		       Form("%s (%.0f#circ) shield cut, tof=[130,170]ns, Ecore>120 MeV, tof_pad=[86,102]",crNames.c_str(),angles),
		       800,0,800);
  TH1F *hdE2 = new TH1F(Form("%s_%.0f_dE2",crNames.c_str(),angles),
		       Form("%s (%.0f#circ) shield cut, tof=[75,83]ns, Ecore>120 MeV, tof_pad=[86,102]",crNames.c_str(),angles),
		       800,0,800);


  // ======================== fill histograms ====================================
  ifstream database("run_database.dat");
  char number[20];
  char name[200];
  string mode;
  string hash = "#";
  float Qcore, Qshield, Qpaddle, Qshield;
  float tof_raw, tof, tof_pad, energy;
  while(database >> number >> name >> mode){
    if(number[0] == hash) continue; //the runs starting with # should be neglected 
    if(mode != "empty") continue;//choose run mode "deuterium/empty"
    //if(atoi(number) < 180) continue;//start from this run
    //if(atoi(number) > 599) break;//end at this run
    printf("Current run number: %s  mode: %s\n", number, mode.c_str());

    TFile *file = new TFile(Form("./root_files/new/%s",name));
    TTree *tree = (TTree*)file->Get("flat_tree");
 
    // tree->SetBranchAddress("Q_0", &Qcore);// ALAINA (40 deg.)
    // tree->SetBranchAddress("Q_16", &Qpaddle);// AMY (40 deg.)
    // tree->SetBranchAddress("Q_8", &Qshield);// ANDY
    // tree->SetBranchAddress("Q_7", &Qcore);// SUSAN (55 deg.)
    // tree->SetBranchAddress("Q_23", &Qpaddle);// MARY (55 deg.)
    tree->SetBranchAddress("Q_3", &Qcore);// JONI (75 deg.)
    tree->SetBranchAddress("Q_11", &Qshield);// JACK 
    tree->SetBranchAddress("Q_19", &Qpaddle);// SAMANTHA (75 deg.)
    tree->SetBranchAddress("Xmax_19", &tof_pad);// SAMANTHA (75 deg.)
    // tree->SetBranchAddress("Q_6", &Qcore);// ROBERTA (90 deg.)
    // tree->SetBranchAddress("Q_22", &Qpaddle);// EVA (90 deg.)
    tree->SetBranchAddress("Q_14", &tof_raw);

    int events = tree->GetEntries();
    for(int i=0; i<events; i++){
      tree->GetEntry(i);

      tof = 0.025*tof_raw*0.5739-9.23582;
      energy = a*Qcore+b;

      // hhind_raw->Fill(tof,energy);
      // hpad_raw->Fill(energy,Qpaddle);

      // if(energy>70 && Qpaddle>150 && Qpaddle<300){
      // 	h2D_1->Fill(tof,energy);
      // 	hpad2D_1->Fill(energy,Qpaddle);
      // }
      
      // if(energy>9){
      // 	if(Qpaddle>0 && Qpaddle<80){
      // 	  h2D_2->Fill(tof,energy);
      // 	  hpad2D_2->Fill(energy,Qpaddle);
      // 	}
      // 	if(Qpaddle>0 && Qpaddle<20){
      // 	  h2D_3->Fill(tof,energy);
      // 	  hpad2D_3->Fill(energy,Qpaddle);
      // 	}
      // 	if(energy<50 && Qpaddle>100){
      // 	  h2D_4->Fill(tof,energy);
      // 	  hpad2D_4->Fill(energy,Qpaddle);
      // 	}
      // }

      /*  ~~~~~~~~~~~~~~~~~~ reference for the TOF and shield cut ~~~~~~~~~~~~~
      const string crNames[nCores] = {"ALAINA", "BROOKE", "CINDY", "JONI",
				      "KRISTA", "LINDA", "ROBERTA", "SUSAN"};
      const Double_t cutLow[nCores] = {64.5, 45.0, 75.0, 75.0,
				       70.5, 61.0, 69.5, 65.5}; //ns
      const Double_t cutHigh[nCores] = {69.0, 50.0, 80.0, 79.0,
					75.5, 65.0, 74.0, 69.5}; //ns
      const Double_t shldCut[nCores] = {4500.0, 3500.0, 3500.0, 3000.0,
				        3500.0, 3000.0, 3000.0, 3300.0};
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

      // if(Qshield<3000.0){
      // 	hhind_shld->Fill(tof,energy);
      // 	hpad_shld->Fill(energy,Qpaddle);
      // }

      // if(tof>75.0 && tof<79.0){
      // 	hhind_tof->Fill(tof,energy);
      // 	hpad_tof->Fill(energy,Qpaddle);

      // 	if(Qshield<3000.0){
      // 	  hhind_tofshld->Fill(tof,energy);
      // 	  hpad_tofshld->Fill(energy,Qpaddle);
      // 	}
      // }
  



      // if(Qshield<3000.0)
      // 	hEVpadtof->Fill(tof, energy);

      // if(Qshield<3000.0 && tof>75.0 && tof<79.0 ){
      // 	hdEVpadtof1->Fill(tof_pad,Qpaddle);
      // 	if(energy>100)
      // 	  hdEVpadtof2->Fill(tof_pad,Qpaddle);
      // 	if(energy>80 && energy<100)
      // 	  hdEVpadtof3->Fill(tof_pad,Qpaddle);
      // 	// if(energy>20 && energy<60)
      // 	//   hdEVpadtof4->Fill(tof_pad,Qpaddle);
      // }

      // if(Qshield<3000.0 && tof>130.0 && tof<170.0){
      // 	hdEVpadtof1_rndm->Fill(tof_pad,Qpaddle);
      // 	if(energy>100)	
      // 	  hdEVpadtof2_rndm->Fill(tof_pad,Qpaddle);
      // 	if(energy>80 && energy<100)
      // 	  hdEVpadtof3_rndm->Fill(tof_pad,Qpaddle);
      // }





      //     if(Qshield<3000.0){
      // 	h2D_1cut->Fill(tof,energy);

      // 	if(tof>75.0 && tof<79.0)
      // 	h1D_2cuts->Fill(energy);
    
      //     // h2D_2cuts->Fill(tof,energy);
      //     if(tof_pad<86 || tof_pad>102)
      // 	h1D_3cuts->Fill(energy);
      // 	// h2D_3cuts->Fill(tof,energy);
      
    

      if(Qshield<3000.0 && tof_pad>86 && tof_pad<102 && energy>120){
	if(tof>130.0 && tof<170.0)
		 hdE->Fill(Qpaddle);
	if(tof>75 && tof<83)
		 hdE2->Fill(Qpaddle);
      }

      // if(tof_pad<83 || tof_pad>102)
      //   h2D_2cut->Fill(tof,energy);


    }//loop over events(i<events)
    file->Close();

  }//loop over while full/empty runs
  database.close();



  // ========= Do random subtraction ========================
  // hdEVpadtof1_net = (TH2F*)hdEVpadtof1->Clone(Form("%s_%.0f_dEvPadTof1_net",crNames.c_str(),angles)); 
  // hdEVpadtof1_net->SetTitle(Form("%s (%.0f#circ) shield+tof cut net",crNames.c_str(),angles));
  // hdEVpadtof1_net->Add(hdEVpadtof1_rndm, -1*(79.0-75.0)/(170.0-130.0));
  // hdEVpadtof2_net = (TH2F*)hdEVpadtof2->Clone(Form("%s_%.0f_dEvPadTof2_net",crNames.c_str(),angles)); 
  // hdEVpadtof2_net->SetTitle(Form("%s (%.0f#circ) shield+tof cut net, Ecore>100MeV",crNames.c_str(),angles));
  // hdEVpadtof2_net->Add(hdEVpadtof2_rndm, -1*(79.0-75.0)/(170.0-130.0));
  // hdEVpadtof3_net = (TH2F*)hdEVpadtof3->Clone(Form("%s_%.0f_dEvPadTof3_net",crNames.c_str(),angles)); 
  // hdEVpadtof3_net->SetTitle(Form("%s (%.0f#circ) shield cut net, Ecore=[80,100]MeV",crNames.c_str(),angles));
  // hdEVpadtof3_net->Add(hdEVpadtof3_rndm, -1*(79.0-75.0)/(170.0-130.0));





  // ============== draw canvas ===================
  TCanvas *c1 = new TCanvas();
  hdE->Draw();
  TCanvas *c2 = new TCanvas();
  hdE2->Draw();

  // TCanvas *c1 = new TCanvas("c1","cut1");
  // c1->Divide(2,1);
  // //gPad->SetLogz();
  // c1->cd(1);
  // hpad2D_1->GetXaxis()->SetTitle("E [MeV]");
  // hpad2D_1->GetXaxis()->SetTitleSize(0.045);
  // hpad2D_1->GetXaxis()->SetLabelSize(0.045);
  // hpad2D_1->GetYaxis()->SetTitle("#DeltaE");
  // hpad2D_1->Draw("color");
  // c1->cd(2);
  // h2D_1->GetXaxis()->SetTitle("TOF [ns]");
  // h2D_1->GetXaxis()->SetTitleSize(0.045);
  // h2D_1->GetXaxis()->SetLabelSize(0.045);
  // h2D_1->GetYaxis()->SetTitle("E [MeV]");
  // h2D_1->GetYaxis()->SetTitleSize(0.045);
  // h2D_1->GetYaxis()->SetLabelSize(0.045);
  // h2D_1->Draw("color");

  // TCanvas *c1 = new TCanvas("c1","1D paddle cut comparison w/ TOF+shld+");
  // //gPad->SetLogz();
  // h1D_2cuts->GetXaxis()->SetTitle("E [MeV]");
  // h1D_2cuts->GetXaxis()->SetTitleSize(0.045);
  // h1D_2cuts->GetXaxis()->SetLabelSize(0.045);
  // h1D_2cuts->GetYaxis()->SetTitle("counts");
  // h1D_2cuts->GetYaxis()->SetTitleSize(0.045);
  // h1D_2cuts->SetLineColor(kBlack);
  // h1D_2cuts->Draw();
  // h1D_3cuts->SetLineColor(kRed);
  // h1D_3cuts->Draw("same");

  // TCanvas *c2 = new TCanvas("c2","2D paddle cut comparison");
  // c2->Divide(1,2);
  // c2->cd(1);
  // gPad->SetLogz();
  // h2D_2cuts->GetXaxis()->SetTitle("TOF [ns]");
  // h2D_2cuts->GetXaxis()->SetTitleSize(0.045);
  // h2D_2cuts->GetXaxis()->SetLabelSize(0.045);
  // h2D_2cuts->GetYaxis()->SetTitle("E [MeV]");
  // h2D_2cuts->GetYaxis()->SetTitleSize(0.045);
  // h2D_2cuts->GetYaxis()->SetLabelSize(0.045);
  // h2D_2cuts->Draw("color,colz");
  // c2->cd(2);
  // gPad->SetLogz();
  // h2D_3cuts->GetXaxis()->SetTitle("TOF [ns]");
  // h2D_3cuts->GetXaxis()->SetTitleSize(0.045);
  // h2D_3cuts->GetXaxis()->SetLabelSize(0.045);
  // h2D_3cuts->GetYaxis()->SetTitle("E [MeV]");
  // h2D_3cuts->GetYaxis()->SetTitleSize(0.045);
  // h2D_3cuts->GetYaxis()->SetLabelSize(0.045);
  // h2D_3cuts->Draw("color,colz");



 

  // ================ save file =============================

  TFile *outfile = new TFile("test_75_multicut_empty.root","recreate");

  hdE->Write();
  hdE2->Write();




  outfile->Close();


}
