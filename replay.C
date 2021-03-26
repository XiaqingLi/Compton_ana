{
  gStyle->SetOptStat(0);
  const int nHINDA = 8;
  const string crNames[nHINDA] = {"ALAINA", "BROOKE", "CINDY", "JONI",
				  "KRISTA", "LINDA", "ROBERTA", "SUSAN"};
  const string shldNames[nHINDA] = {"ANDY", "BOBBY", "CARL", "JACK",
				    "KEVIN", "LARRY", "REGGIE", "STEVE"};
  const Double_t angles[nHINDA] = {55.0, 90.0, 125.0, 125.0, 
				   125.0, 90.0, 90.0, 55.0};
  const string pos[nHINDA] = {"D", "D", "L", "D", "R", "R", "L", "L"};
  
  const Double_t crCut[nHINDA] = {9, 11.5, 12, 8.5, 21, 12.5, 17, 13}; 
  const Double_t randLow1[nHINDA] = {0, 0, 0, 0, 0, 0, 0, 0};
  const Double_t randHigh1[nHINDA] = {0, 0, 20, 0, 0, 0, 0, 0};
  const Double_t randLow2[nHINDA] = {120, 90, 0, 155, 148, 140, 130, 120};
  const Double_t randHigh2[nHINDA] = {170, 160, 0, 170, 173, 155, 170, 160};


  const Double_t coef[nHINDA] = {22.9167, 10.9375, 12.5, 16.6667, 
    				 29.1667, 4.86111, 27.5, 6.54762}; //coef. for timewalk correction
  //after shield cut optimization 
  const Double_t shldCut[nHINDA] = {3500.0, 14000.0, 3500.0, 3500.0,
  				    3000.0, 1900.0, 1500.0, 2000.0}; //Q_i

  const Double_t tofShift[5] = {8.5, 10.5, 1.9, 1.0, 0};
  const Double_t tofCenter[nHINDA] = {94, 64.5, 118, 111.5, 106, 86.5, 105, 95};
  const Double_t tofWidth[nHINDA] = {14, 11, 14, 12, 16, 13, 11, 9};
  Double_t tofplus_ns = 0;

  Double_t tofCutLow[nHINDA], tofCutHigh[nHINDA];
  float Qcore[nHINDA];
  float Qshield[nHINDA];
  float tof_raw, tofpad_raw, tof, energy, delta;

  //============== set energy calibration parameters ==============
  const Double_t kY1 = 0.0; //MeV
  const Double_t kX1 = 0.0; //0.001*Q_i
  const Double_t E_beam = 84.0; //MeV, incident photon energy
  const Double_t Mass = 3727.4; //MeV, 4He nuclues mass
  Double_t kY2[nHINDA];
  for(int i=0; i<nHINDA; i++){
    kY2[i] = E_beam/(1+(E_beam/Mass)*(1-cos(angles[i]*TMath::Pi()/180.0)));
  }
  const Double_t kX2[nHINDA] = {5000, 4400, 3600, 4100,
  				2900, 3100, 2100, 6800};// Q_i
	
  Double_t m_cal_[nHINDA];
  Double_t b_cal_[nHINDA];
  for(int i=0; i<nHINDA; i++) {
    m_cal_[i] = (kY2[i] - kY1) / (kX2[i] - kX1);
    b_cal_[i] = kY2[i] - m_cal_[i] * kX2[i];
  }


  // ============  make histograms  ===============
  const Double_t Emax = 300.0;
  const Int_t nE = 1500;

  TH1F *hcore_raw[nHINDA];
  TH1F *hcore_2Cuts_raw[nHINDA];
  TH1F *hshield_raw[nHINDA];
  TH1F *hcore[nHINDA];
  TH1F *hcore_shldCut[nHINDA];
  TH1F *hcore_tofCut[nHINDA];
  TH1F *hcore_2Cuts[nHINDA];
  TH1F *hcore_rand[nHINDA];
  TH1F *hcore_net[nHINDA];
  TH2F *hEcr_tof[nHINDA];
  TH2F *hEcr_tof_shldCut[nHINDA];
  TH2F *hEcr_Eshld[nHINDA];
  TH2F *hEcr_Eshld_tofcut[nHINDA];

  for(int i=0; i<nHINDA; i++){ 
    hcore_raw[i] = new TH1F(Form("%s_Raw",crNames[i].c_str()),
			    Form("%s #theta=%0.f %s, Raw",crNames[i].c_str(),angles[i],pos[i].c_str()), 500,0,10000);
    hcore_raw[i]->SetLineColor(kBlue);
    hcore_raw[i]->Sumw2();
    hcore_2Cuts_raw[i] = new TH1F(Form("%s_2Cuts_raw",crNames[i].c_str()),
				  Form("%s #theta=%0.f %s, 2 Cuts (Raw)",crNames[i].c_str(),angles[i],pos[i].c_str()), 500,0,10000);
    hcore[i] = new TH1F(Form("%s",crNames[i].c_str()),
			Form("%s #theta=%0.f %s",crNames[i].c_str(),angles[i],pos[i].c_str()), 1500,0,300);
    hcore[i]->Sumw2();
    hcore_shldCut[i] = new TH1F(Form("%s_shldCut",crNames[i].c_str()),
				Form("%s #theta=%0.f %s, Shield Cut",crNames[i].c_str(),angles[i],pos[i].c_str()), 1500,0,300);
    hcore_shldCut[i]->SetLineColor(kCyan);
    hcore_shldCut[i]->Sumw2();
    hcore_tofCut[i] = new TH1F(Form("%s_tofCut",crNames[i].c_str()),
			       Form("%s #theta=%0.f %s, ToF Cut",crNames[i].c_str(),angles[i],pos[i].c_str()), 1500,0,300);
    hcore_tofCut[i]->SetLineColor(kGreen);
    hcore_tofCut[i]->Sumw2();

    hcore_2Cuts[i] = new TH1F(Form("%s_2Cuts",crNames[i].c_str()),
			      Form("%s #theta=%0.f %s, 2 Cuts",crNames[i].c_str(),angles[i],pos[i].c_str()), 1500,0,300);
    hcore_2Cuts[i]->GetXaxis()->SetTitle("E [MeV]");
    hcore_2Cuts[i]->GetXaxis()->SetTitleSize(0.045);
    hcore_2Cuts[i]->GetXaxis()->SetLabelSize(0.045);
    hcore_2Cuts[i]->GetYaxis()->SetTitle("Counts / 0.2 MeV");
    hcore_2Cuts[i]->GetYaxis()->SetTitleSize(0.045);
    hcore_2Cuts[i]->GetYaxis()->SetLabelSize(0.045);
    hcore_2Cuts[i]->GetXaxis()->SetRangeUser(20,200);
    hcore_2Cuts[i]->SetFillStyle(3003);
    hcore_2Cuts[i]->SetFillColor(kRed+1);
    hcore_2Cuts[i]->SetLineColor(kRed+1);

    hcore_rand[i] = new TH1F(Form("%s_rand",crNames[i].c_str()),
			     Form("%s #theta=%0.f %s, random events",crNames[i].c_str(),angles[i],pos[i].c_str()), 1500,0,300);
    hcore_rand[i]->Sumw2();
    hcore_net[i] = new TH1F(Form("%s_net",crNames[i].c_str()),
			    Form("%s #theta=%0.f %s, random subtracted",crNames[i].c_str(),angles[i],pos[i].c_str()), 1500,0,300);
    hcore_net[i]->Sumw2();
    hEcr_tof[i] = new TH2F(Form("%s_2D",crNames[i].c_str()),
			   Form("%s #theta=%0.f %s",crNames[i].c_str(),angles[i],pos[i].c_str()), 
			   1000,0,200, 1500,0,300);
    hEcr_tof[i]->GetXaxis()->SetTitle("ToF [ns]");
    hEcr_tof[i]->GetXaxis()->SetTitleSize(0.045);
    hEcr_tof[i]->GetXaxis()->SetLabelSize(0.045);
    hEcr_tof[i]->GetYaxis()->SetTitle("E [MeV]");
    hEcr_tof[i]->GetYaxis()->SetTitleSize(0.045);
    hEcr_tof[i]->GetYaxis()->SetLabelSize(0.045);
    hEcr_tof_shldCut[i] = new TH2F(Form("%s_2D_shldCut",crNames[i].c_str()),
				   Form("%s #theta=%0.f %s, w/ Shield Cut",crNames[i].c_str(),angles[i],pos[i].c_str()), 
				   1000,0,200, 1500,0,300);
    hEcr_Eshld[i] = new TH2F(Form("%s_cr_shld",crNames[i].c_str()),
			     Form("%s #theta=%0.f %s",crNames[i].c_str(),angles[i],pos[i].c_str()), 
			     1500,0,300, 10000,0,50000);
    hEcr_Eshld_tofcut[i] = new TH2F(Form("%s_cr_shld_tofcut",crNames[i].c_str()),
				    Form("%s #theta=%0.f %s",crNames[i].c_str(),angles[i],pos[i].c_str()), 
			     1500,0,300, 10000,0,50000);
    hshield_raw[i] = new TH1F(Form("%s_Raw",shldNames[i].c_str()),
			      Form("%s #theta=%0.f %s",crNames[i].c_str(),angles[i],pos[i].c_str()), 
			      500,0,50000);
  }



  // ======================== fill histograms ====================================
  ifstream database("/var/phy/project/mepg/xl79/helium_84MeV/run_database.dat");
  char number[20];
  char filename[200];
  string mode;
  int helicity;
  string hash = "#";
  int runnumber;
  TTree *tree;

  while(database >> number >> filename >> mode >> helicity){
    runnumber = atoi(number);
    if(number[0] == hash) continue; //the runs starting with # should be neglected 
    if(mode != "helium_empty") continue;//choose run mode "helium_full/helium_empty"
    // if(atoi(number) != 942) continue;
    // if(runnumber < 1016) continue;//start from this run
    // if(runnumber > 928) break;//end at this run
    // printf("Current run number: %s  mode: %s\n", number, mode.c_str());
    // cout<<number<<"\t"<<filename<<"\t"<<endl;

    if(runnumber>=928 && runnumber<932)
      delta = tofShift[0];
    else if(runnumber>=932 && runnumber<941)
      delta = tofShift[1];
    else if(runnumber>=941 && runnumber<971)
      delta = tofShift[2];
    else if(runnumber>=971 && runnumber<1016)
      delta = tofShift[3];
    else if(runnumber>=1016 && runnumber<1049)
      delta = tofShift[4];
    else 
      continue;

    TFile *infile = new TFile(Form("/var/phy/project/mepg/xl79/helium_84MeV/root_files/%s",filename));
    tree = (TTree*)infile->Get("flat_tree");

    tree->SetBranchAddress("Q_0", &Qcore[0]);// ALAINA (55 deg. down)
    tree->SetBranchAddress("Q_1", &Qcore[1]);// BROOKE (90 deg. down)
    tree->SetBranchAddress("Q_2", &Qcore[2]);// CINDY (125 deg. left)
    tree->SetBranchAddress("Q_3", &Qcore[3]);// JONI (125 deg. down)
    tree->SetBranchAddress("Q_4", &Qcore[4]);// KRISTA (125 deg. right)
    tree->SetBranchAddress("Q_5", &Qcore[5]);// LINDA (90 deg. right)
    tree->SetBranchAddress("Q_6", &Qcore[6]);// ROBERTA (90 deg. left)
    tree->SetBranchAddress("Q_7", &Qcore[7]);// SUSAN (55 deg. left)
    tree->SetBranchAddress("Q_8", &Qshield[0]);// ANDY
    tree->SetBranchAddress("Q_9", &Qshield[1]);// BOBBY
    tree->SetBranchAddress("Q_10", &Qshield[2]);// CARL
    tree->SetBranchAddress("Q_11", &Qshield[3]);// JACK
    tree->SetBranchAddress("Q_12", &Qshield[4]);// KEVIN
    tree->SetBranchAddress("Q_13", &Qshield[5]);// LARRY
    tree->SetBranchAddress("Q_31", &Qshield[6]);// REGGIE (digitizer B)
    tree->SetBranchAddress("Q_14", &tof_raw);
    tree->SetBranchAddress("Q_15", &Qshield[7]);// STEVE
 

    int events = tree->GetEntries();
    cout<<"current run: "<<runnumber<<"\t total events: "<<events<<endl;
    for(int i=0; i<events; i++){

      // if(i%100000 == 0)
      // 	cout<<i<<endl;

      tree->GetEntry(i);

      tof = 0.025*tof_raw*0.5739-9.23582;

      for(int j=0; j<nHINDA; j++){
	if(runnumber<952 && j==2) continue;

	energy = Qcore[j]*m_cal_[j]+b_cal_[j]; //do energy calibration

	tofCutLow[j] = tofCenter[j]+delta-(tofWidth[j]+tofplus_ns)/2.0;
	tofCutHigh[j] = tofCenter[j]+delta+(tofWidth[j]+tofplus_ns)/2.0;

	hcore_raw[j]->Fill(Qcore[j]); 
	hshield_raw[j]->Fill(Qshield[j]);
	hcore[j]->Fill(energy);

	if(energy>kY2[j])
	  tof = tof-(energy-kY2[j])/coef[j];

	if(energy>crCut[j]){
	  hEcr_tof[j]->Fill(tof,energy);
	  hEcr_Eshld[j]->Fill(energy,Qshield[j]);

	  if(tof>tofCutLow[j] && tof<tofCutHigh[j]){ //tof cut only
	    hcore_tofCut[j]->Fill(energy);
	    hEcr_Eshld_tofcut[j]->Fill(energy,Qshield[j]);
	  }
	  if(Qshield[j]<shldCut[j]*1.0){//shield cut 
	    hcore_shldCut[j]->Fill(energy);
	    hEcr_tof_shldCut[j]->Fill(tof,energy);

	    if(tof>tofCutLow[j] && tof<tofCutHigh[j]){//select prompt events
	      hcore_2Cuts[j]->Fill(energy);
	      hcore_2Cuts_raw[j]->Fill(Qcore[j]);
	      hcore_net[j]->Fill(energy); 
	    }//**************fill prompt events

	    if((tof>randLow1[j] && tof<randHigh1[j])||(tof>randLow2[j] && tof<randHigh2[j])){ 
	      hcore_rand[j]->Fill(energy); 
	    }//**************fill random events
	  }//shield cut
	}//core threshold
      }//loop over HINDA
    }//loop over events
  }//loop over while
  database.close();


  // =============== do random subtraction ==================
  for(int j=0; j<nHINDA; j++){
     hcore_net[j]->Add(hcore_rand[j],-1.*tofWidth[j]/((randHigh1[j]-randLow1[j])+(randHigh2[j]-randLow2[j])));
  }


  // ======================== save file ===================
  TFile *fout = new TFile("Empty.root","recreate");
  for(int i=0; i<nHINDA; i++){
    hcore_raw[i]->Write();
    hcore_2Cuts_raw[i]->Write();
    hshield_raw[i]->Write();
    hcore[i]->Write();
    hcore_shldCut[i]->Write();
    hcore_tofCut[i]->Write();
    hcore_2Cuts[i]->Write();
    hcore_rand[i]->Write();
    hcore_net[i]->Write();
    hEcr_tof[i]->Write();
    hEcr_tof_shldCut[i]->Write();
    hEcr_Eshld[i]->Write();
    hEcr_Eshld_tofcut[i]->Write();
  }

  fout->Close();



}
