using namespace std;

void Replay_hinda_paddle(){
  gStyle->SetOptStat(0);

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




  // TH2F *h2D_high = new TH2F(Form("%s_high_cut",crNames.c_str()),
  // 			    Form("%s_%.0f#circ",crNames.c_str(),angles),
  // 			    1000,0,200, 2000,0,200);
  // TH2F *h2D_low = new TH2F(Form("%s_low_cut",crNames.c_str()),
  // 			   Form("%s_%.0f#circ",crNames.c_str(),angles),
  // 			   1000,0,200, 2000,0,200);
  // TH2F *hpad2D_high = new TH2F(Form("paddle_%.0f_high",angles),
  // 			       Form("paddle_%.0f#circ",angles),
  // 			       1000,0,200, 2500,0,500);
  // TH2F *hpad2D_low = new TH2F(Form("paddle_%.0f_low",angles),
  // 			      Form("paddle_%.0f#circ",angles),
  // 			      1000,0,200, 2500,0,500);


  // TH2F *h2D_1 = new TH2F(Form("%s_cut1",crNames.c_str()),
  // 			 Form("%s_%.0f#circ",crNames.c_str(),angles),
  // 			 1000,0,200, 2000,0,200);
  // TH2F *hpad2D_1 = new TH2F(Form("paddle_%.0f_cut1",angles),
  // 			    Form("paddle_%.0f#circ",angles),
  // 			    1000,0,200, 2500,0,500);
  // TH2F *h2D_2 = new TH2F(Form("%s_cut2",crNames.c_str()),
  // 			 Form("%s_%.0f#circ",crNames.c_str(),angles),
  // 			 1000,0,200, 2000,0,200);
  // TH2F *hpad2D_2 = new TH2F(Form("paddle_%.0f_cut2",angles),
  // 			    Form("paddle_%.0f#circ",angles),
  // 			    1000,0,200, 2500,0,500);
  // TH2F *h2D_3 = new TH2F(Form("%s_cut3",crNames.c_str()),
  // 			 Form("%s_%.0f#circ",crNames.c_str(),angles),
  // 			 1000,0,200, 2000,0,200);
  // TH2F *hpad2D_3 = new TH2F(Form("paddle_%.0f_cut3",angles),
  // 			    Form("paddle_%.0f#circ",angles),
  // 			    1000,0,200, 2500,0,500);
  // TH2F *h2D_4 = new TH2F(Form("%s_cut4",crNames.c_str()),
  // 			 Form("%s_%.0f#circ",crNames.c_str(),angles),
  // 			 1000,0,200, 2000,0,200);
  // TH2F *hpad2D_4 = new TH2F(Form("paddle_%.0f_cut4",angles),
  // 			    Form("paddle_%.0f#circ",angles),
  // 			    1000,0,200, 2500,0,500);
  // TH2F *h2D_5 = new TH2F(Form("hinda_%s_TofCut",crNames.c_str()),
  // 			 Form("hinda_%s_%.0f#circ",crNames.c_str(),angles),
  // 			 1000,0,200, 2000,0,200);
  // TH2F *hpad2D_5 = new TH2F(Form("paddle_%.0f_TofCut",angles),
  // 			    Form("paddle_%.0f#circ",angles),
  // 			    1000,0,200, 2500,0,500);
  // TH2F *h2D_6 = new TH2F(Form("hinda_%s_TofShldCut",crNames.c_str()),
  // 			 Form("hinda_%s_%.0f#circ",crNames.c_str(),angles),
  // 			 1000,0,200, 2000,0,200);
  // TH2F *hpad2D_6 = new TH2F(Form("paddle_%.0f_TofShldCut",angles),
  // 			    Form("paddle_%.0f#circ",angles),
  // 			    1000,0,200, 2500,0,500);
   

  TH2F *hhind_raw = new TH2F(Form("hinda_%s_raw",crNames.c_str()),
			     Form("hinda_%s_%.0f#circ_no_cut",crNames.c_str(),angles),
			     // 1000,0,200, 2000,0,200);
			     1000,0,1000, 2000,0,200);
  TH2F *hpad_raw = new TH2F(Form("paddle_%.0f_raw",angles),
			    Form("paddle_%.0f#circ_no_cut",angles),
			    // 1000,0,200, 2500,0,500);
			    1000,0,1000, 4000,0,800);
  TH2F *hhind_tof = new TH2F(Form("hinda_%s_TofCut",crNames.c_str()),
			     Form("hinda_%s_%.0f#circ_TOF_cut",crNames.c_str(),angles),
			     // 1000,0,200, 2000,0,200);
			     1000,0,1000, 2000,0,200);
  TH2F *hpad_tof = new TH2F(Form("paddle_%.0f_TofCut",angles),
			    Form("paddle_%.0f#circ_TOF_cut",angles),
			    // 1000,0,200, 2500,0,500);
			    1000,0,1000, 4000,0,800);
  TH2F *hhind_shld = new TH2F(Form("hinda_%s_ShldCut",crNames.c_str()),
			      Form("hinda_%s_%.0f#circ_shield_cut",crNames.c_str(),angles),
			      // 1000,0,200, 2000,0,200);
			      1000,0,1000, 2000,0,200);
  TH2F *hpad_shld = new TH2F(Form("paddle_%.0f_ShldCut",angles),
			     Form("paddle_%.0f#circ_shield_cut",angles),
			     // 1000,0,200, 2500,0,500);
			     1000,0,1000, 4000,0,800);
  TH2F *hhind_tofshld = new TH2F(Form("hinda_%s_TofShldCut",crNames.c_str()),
				 Form("hinda_%s_%.0f#circ_TOF+shield_cut",crNames.c_str(),angles),
				 // 1000,0,200, 2000,0,200);
				 1000,0,1000, 2000,0,200);
  TH2F *hpad_tofshld = new TH2F(Form("paddle_%.0f_TofShldCut",angles),
				Form("paddle_%.0f#circ_TOF+shield_cut",angles),
				// 1000,0,200, 2000,0,200);
				1000,0,1000, 4000,0,800);
  TH2F *hhind_rand = new TH2F(Form("hinda_%s_Rndm",crNames.c_str()),
			      Form("hinda_%s_%.0f#circ_Random",crNames.c_str(),angles),
			      // 1000,0,200, 2000,0,200);
			      1000,0,1000, 2000,0,200);
  TH2F *hpad_rand = new TH2F(Form("paddle_%.0f_Rndm",angles),
			     Form("paddle_%.0f#circ_Random",angles),
			     // 1000,0,200, 2000,0,200);
			     1000,0,1000, 4000,0,800);
  TH2F *hhind_net, *hpad_net;      
  double scale = (79.-75.)/(170.-130.);


  ifstream database("run_database.dat");
  char number[20];
  char name[200];
  string mode;
  string hash = "#";
  float Qcore, Qshield, Qpaddle, Qshield;
  float tof_raw, tof, tof_paddle, energy;
  while(database >> number >> name >> mode){
    if(number[0] == hash) continue; //the runs starting with # should be neglected 
    if(mode != "empty") continue;//choose run mode "deuterium/empty"
    //if(atoi(number) < 695) continue;//start from this run
    //if(atoi(number) > 696) break;//end at this run
    printf("Current run number: %s  mode: %s\n", number, mode.c_str());

    TFile *file = new TFile(Form("./root_files/%s",name));
    TTree *tree = (TTree*)file->Get("flat_tree");
 
    // tree->SetBranchAddress("Q_0", &Qcore);// ALAINA (40 deg.)
    // tree->SetBranchAddress("Q_16", &Qpaddle);// AMY (40 deg.)
    // tree->SetBranchAddress("Q_8", &Qshield);// ANDY (40 deg.)
    // tree->SetBranchAddress("Q_7", &Qcore);// SUSAN (55 deg.)
    // tree->SetBranchAddress("Q_23", &Qpaddle);// MARY (55 deg.)
    tree->SetBranchAddress("Q_3", &Qcore);// JONI (75 deg.)
    tree->SetBranchAddress("Q_11", &Qshield);// JACK (75 deg.) 
    tree->SetBranchAddress("Q_19", &Qpaddle);// SAMANTHA (75 deg.)
    tree->SetBranchAddress("Xmax_19", &tof_paddle);//SAMANTHA (75 deg.)
    // tree->SetBranchAddress("Q_6", &Qcore);// ROBERTA (90 deg.)
    // tree->SetBranchAddress("Q_22", &Qpaddle);// EVA (90 deg.)
    tree->SetBranchAddress("Q_14", &tof_raw);

    int events = tree->GetEntries();
    for(int i=0; i<events; i++){
      tree->GetEntry(i);

      tof = 0.025*tof_raw*0.5739-9.23582;
      energy = a*Qcore+b;

      hhind_raw->Fill(tof_paddle,energy);
      hpad_raw->Fill(tof_paddle,Qpaddle);

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
      */

      // if(energy>15 && Qshield<3000.0){
      // 	// hhind_shld->Fill(tof_paddle,energy);
      // 	// hpad_shld->Fill(tof_paddle,Qpaddle);
      // 	if(tof>75.0 && tof<79.0){
      // 	  hhind_tofshld->Fill(tof_paddle,energy);
      // 	  hpad_tofshld->Fill(tof_paddle,Qpaddle);
      // 	}
      // 	if(tof>130.0 && tof<170.0){
      // 	  hhind_rand->Fill(tof_paddle,energy);
      // 	  hpad_rand->Fill(tof_paddle,Qpaddle);
      // 	}
      // }

      if(energy>15 && Qshield<3000.0){
	// hhind_shld->Fill(tof_paddle,energy);
	// hpad_shld->Fill(tof_paddle,Qpaddle);
	if(tof>75.0 && tof<79.0){
	  hhind_tofshld->Fill(tof_paddle,energy);
	  hpad_tofshld->Fill(tof_paddle,Qpaddle);
	}
	if(tof>130.0 && tof<170.0){
	  hhind_rand->Fill(tof_paddle,energy);
	  hpad_rand->Fill(tof_paddle,Qpaddle);
	}
      }
 

    }//loop over events(i<events)
    file->Close();

  }//loop over while full/empty runs
  database.close();

  hhind_net = (TH2F*)hhind_tofshld->Clone(Form("hinda_%s_Net",crNames.c_str()));
  hpad_net = (TH2F*)hpad_tofshld->Clone(Form("paddle_%s_Net",crNames.c_str()));
  hhind_net->Add(hhind_rand,-1.*scale);
  hpad_net->Add(hpad_rand,-1.*scale);







  // ============== draw canvas ===================

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

  // TCanvas *c2 = new TCanvas("c2","cut2");
  // c2->Divide(2,1);
  // //gPad->SetLogz();
  // c2->cd(1);
  // hpad2D_2->GetXaxis()->SetTitle("E [MeV]");
  // hpad2D_2->GetXaxis()->SetTitleSize(0.045);
  // hpad2D_2->GetXaxis()->SetLabelSize(0.045);
  // hpad2D_2->GetYaxis()->SetTitle("#DeltaE");
  // hpad2D_2->Draw("color");
  // c2->cd(2);
  // h2D_2->GetXaxis()->SetTitle("TOF [ns]");
  // h2D_2->GetXaxis()->SetTitleSize(0.045);
  // h2D_2->GetXaxis()->SetLabelSize(0.045);
  // h2D_2->GetYaxis()->SetTitle("E [MeV]");
  // h2D_2->GetYaxis()->SetTitleSize(0.045);
  // h2D_2->GetYaxis()->SetLabelSize(0.045);
  // h2D_2->Draw("color");

  // TCanvas *c3 = new TCanvas("c3","cut3");
  // c3->Divide(2,1);
  // //gPad->SetLogz();
  // c3->cd(1);
  // hpad2D_3->GetXaxis()->SetTitle("E [MeV]");
  // hpad2D_3->GetXaxis()->SetTitleSize(0.045);
  // hpad2D_3->GetXaxis()->SetLabelSize(0.045);
  // hpad2D_3->GetYaxis()->SetTitle("#DeltaE");
  // hpad2D_3->Draw("color");
  // c3->cd(2);
  // h2D_3->GetXaxis()->SetTitle("TOF [ns]");
  // h2D_3->GetXaxis()->SetTitleSize(0.045);
  // h2D_3->GetXaxis()->SetLabelSize(0.045);
  // h2D_3->GetYaxis()->SetTitle("E [MeV]");
  // h2D_3->GetYaxis()->SetTitleSize(0.045);
  // h2D_3->GetYaxis()->SetLabelSize(0.045);
  // h2D_3->Draw("color");

  // TCanvas *c4 = new TCanvas("c4","cut4");
  // c4->Divide(2,1);
  // //gPad->SetLogz();
  // c4->cd(1);
  // hpad2D_4->GetXaxis()->SetTitle("E [MeV]");
  // hpad2D_4->GetXaxis()->SetTitleSize(0.045);
  // hpad2D_4->GetXaxis()->SetLabelSize(0.045);
  // hpad2D_4->GetYaxis()->SetTitle("#DeltaE");
  // hpad2D_4->Draw("color");
  // c4->cd(2);
  // h2D_4->GetXaxis()->SetTitle("TOF [ns]");
  // h2D_4->GetXaxis()->SetTitleSize(0.045);
  // h2D_4->GetXaxis()->SetLabelSize(0.045);
  // h2D_4->GetYaxis()->SetTitle("E [MeV]");
  // h2D_4->GetYaxis()->SetTitleSize(0.045);
  // h2D_4->GetYaxis()->SetLabelSize(0.045);
  // h2D_4->Draw("color");

  // TCanvas *c5 = new TCanvas("TofCut","Tof cut");
  // c5->Divide(2,1);
  // //gPad->SetLogz();
  // c5->cd(1);
  // hpad2D_5->GetXaxis()->SetTitle("E [MeV]");
  // hpad2D_5->GetXaxis()->SetTitleSize(0.045);
  // hpad2D_5->GetXaxis()->SetLabelSize(0.045);
  // hpad2D_5->GetYaxis()->SetTitle("#DeltaE");
  // hpad2D_5->Draw("color");
  // c5->cd(2);
  // h2D_5->GetXaxis()->SetTitle("TOF [ns]");
  // h2D_5->GetXaxis()->SetTitleSize(0.045);
  // h2D_5->GetXaxis()->SetLabelSize(0.045);
  // h2D_5->GetYaxis()->SetTitle("E [MeV]");
  // h2D_5->GetYaxis()->SetTitleSize(0.045);
  // h2D_5->GetYaxis()->SetLabelSize(0.045);
  // h2D_5->Draw("color");

  // TCanvas *c6 = new TCanvas("TofShldCut","Tof+shield cut");
  // c6->Divide(2,1);
  // //gPad->SetLogz();
  // c6->cd(1);
  // hpad2D_6->GetXaxis()->SetTitle("E [MeV]");
  // hpad2D_6->GetXaxis()->SetTitleSize(0.045);
  // hpad2D_6->GetXaxis()->SetLabelSize(0.045);
  // hpad2D_6->GetYaxis()->SetTitle("#DeltaE");
  // hpad2D_6->Draw("color");
  // c6->cd(2);
  // h2D_6->GetXaxis()->SetTitle("TOF [ns]");
  // h2D_6->GetXaxis()->SetTitleSize(0.045);
  // h2D_6->GetXaxis()->SetLabelSize(0.045);
  // h2D_6->GetYaxis()->SetTitle("E [MeV]");
  // h2D_6->GetYaxis()->SetTitleSize(0.045);
  // h2D_6->GetYaxis()->SetLabelSize(0.045);
  // h2D_6->Draw("color");

  TCanvas *c1 = new TCanvas("c1","pad no cut");
  gPad->SetLogz();
  hpad_raw->GetXaxis()->SetTitle("paddle TOF");
  hpad_raw->GetXaxis()->SetTitleSize(0.045);
  hpad_raw->GetXaxis()->SetLabelSize(0.045);
  hpad_raw->GetYaxis()->SetTitle("E [MeV]");
  hpad_raw->GetYaxis()->SetTitleSize(0.045);
  hpad_raw->GetYaxis()->SetLabelSize(0.045);
  hpad_raw->Draw("colz");

  TCanvas *c2 = new TCanvas("c2","hinda no cut");
  gPad->SetLogz();
  hhind_raw->GetXaxis()->SetTitle("paddle TOF");
  hhind_raw->GetXaxis()->SetTitleSize(0.045);
  hhind_raw->GetXaxis()->SetLabelSize(0.045);
  hhind_raw->GetYaxis()->SetTitle("#DeltaE");
  hhind_raw->GetYaxis()->SetTitleSize(0.045);
  hhind_raw->GetYaxis()->SetLabelSize(0.045);
  hhind_raw->Draw("colz");

  TCanvas *c3 = new TCanvas("c3","pad tof+shield cut");
  gPad->SetLogz();
  hpad_tofshld->GetXaxis()->SetTitle("paddle TOF");
  hpad_tofshld->GetXaxis()->SetTitleSize(0.045);
  hpad_tofshld->GetXaxis()->SetLabelSize(0.045);
  hpad_tofshld->GetYaxis()->SetTitle("E [MeV]");
  hpad_tofshld->GetYaxis()->SetTitleSize(0.045);
  hpad_tofshld->GetYaxis()->SetLabelSize(0.045);
  hpad_tofshld->Draw("colz");

  TCanvas *c4 = new TCanvas("c4","hinda tof+shiled cut");
  gPad->SetLogz();
  hhind_tofshld->GetXaxis()->SetTitle("paddle TOF");
  hhind_tofshld->GetXaxis()->SetTitleSize(0.045);
  hhind_tofshld->GetXaxis()->SetLabelSize(0.045);
  hhind_tofshld->GetYaxis()->SetTitle("#DeltaE");
  hhind_tofshld->GetYaxis()->SetTitleSize(0.045);
  hhind_tofshld->GetYaxis()->SetLabelSize(0.045);
  hhind_tofshld->Draw("colz");

  TCanvas *c5 = new TCanvas("c5","pad net");
  gPad->SetLogz();
  hpad_tofshld->GetXaxis()->SetTitle("paddle TOF");
  hpad_tofshld->GetXaxis()->SetTitleSize(0.045);
  hpad_tofshld->GetXaxis()->SetLabelSize(0.045);
  hpad_tofshld->GetYaxis()->SetTitle("E [MeV]");
  hpad_tofshld->GetYaxis()->SetTitleSize(0.045);
  hpad_tofshld->GetYaxis()->SetLabelSize(0.045);
  hpad_tofshld->Draw("colz");

  TCanvas *c6 = new TCanvas("c6","hinda net");
  gPad->SetLogz();
  hhind_tofshld->GetXaxis()->SetTitle("paddle TOF");
  hhind_tofshld->GetXaxis()->SetTitleSize(0.045);
  hhind_tofshld->GetXaxis()->SetLabelSize(0.045);
  hhind_tofshld->GetYaxis()->SetTitle("#DeltaE");
  hhind_tofshld->GetYaxis()->SetTitleSize(0.045);
  hhind_tofshld->GetYaxis()->SetLabelSize(0.045);
  hhind_tofshld->Draw("colz");

  // TCanvas *c2 = new TCanvas("c2","shield cut");
  // gPad->SetLogz();
  // hpad_shld->GetXaxis()->SetTitle("E [MeV]");
  // hpad_shld->GetXaxis()->SetTitleSize(0.045);
  // hpad_shld->GetXaxis()->SetLabelSize(0.045);
  // hpad_shld->GetYaxis()->SetTitle("#DeltaE");
  // hpad_shld->Draw("color,colz");

  // TCanvas *c3 = new TCanvas("c3","tof cut");
  // gPad->SetLogz();
  // hpad_tof->GetXaxis()->SetTitle("E [MeV]");
  // hpad_tof->GetXaxis()->SetTitleSize(0.045);
  // hpad_tof->GetXaxis()->SetLabelSize(0.045);
  // hpad_tof->GetYaxis()->SetTitle("#DeltaE");
  // hpad_tof->Draw("color,colz");

  // TCanvas *c4 = new TCanvas("c4","tof+shield cut");
  // gPad->SetLogz();
  // hpad_tofshld->GetXaxis()->SetTitle("E [MeV]");
  // hpad_tofshld->GetXaxis()->SetTitleSize(0.045);
  // hpad_tofshld->GetXaxis()->SetLabelSize(0.045);
  // hpad_tofshld->GetYaxis()->SetTitle("#DeltaE");
  // hpad_tofshld->Draw("color,colz");
  // ================ save file =============================

  TFile *outfile = new TFile("test_75_pad_tof_empty.root","recreate");
  // hpad2D_1->Write();
  // h2D_1->Write();
  // hpad2D_2->Write();
  // h2D_2->Write();
  // hpad2D_3->Write();
  // h2D_3->Write();
  // hpad2D_4->Write();
  // h2D_4->Write();
  // hpad2D_5->Write();
  // h2D_5->Write();
  // hpad2D_6->Write();
  // h2D_6->Write();
  // c1->Write();
  // c2->Write();
  // c3->Write();
  // c4->Write();
  // c5->Write();
  // c6->Write();

  hhind_raw->Write();
  hpad_raw->Write();
  // hhind_shld->Write();
  // hpad_shld->Write();
  // hhind_tof->Write();
  // hpad_tof->Write();
  hhind_tofshld->Write();
  hpad_tofshld->Write();
  hhind_rand->Write();
  hpad_rand->Write();
  hhind_net->Write();
  hpad_net->Write();
  c1->Write();
  c2->Write();
  c3->Write();
  c4->Write();
  c5->Write();
  c6->Write();


  outfile->Close();


}
