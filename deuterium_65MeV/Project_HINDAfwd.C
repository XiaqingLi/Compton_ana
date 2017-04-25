using namespace std;
// const int nCores = 8;
// const int nShields = 8;
// const string crNames[nCores] = {"ALAINA", "BROOKE", "CINDY", "JONI",
// 				  "KRISTA", "LINDA", "ROBERTA", "SUSAN"};
// const string shldNames[nShields] = {"ANDY", "BOBBY", "CARL", "JACK",
// 				      "KEVIN", "LARRY", "REGGIE", "STEVE"};
// const Double_t angles[nCores] = {40.0, 159.0, 125.0, 75.0, 
// 				   145.0, 110.0, 90.0, 55.0};

// const Double_t cutLow[nCores] = {82.0, 64.0, 91.0, 91.5,
//                                  87.5, 78.5, 87.0, 83.0};
// const Double_t cutHigh[nCores] = {86.0, 69.0, 97.0, 95.5,
//                                   93.5, 82.5, 91.0, 86.5};
// const Double_t shldCut[nCores] = {4500.0, 0.0, 3000.0, 3000.0,
// 				    4300.0, 3000.0, 2700.0, 3300.0};
// const Double_t crCut[nCores] = {9.0, 10.0, 12.0, 9.0
// 				  19.0, 14.0, 19.0, 11.0};
// const Double_t randLow = 130.;
// const Double_t randHigh = 170.;

// //============== for calibration ==========================
// const E_beam = 65; //MeV
// const Double_t M_D = 1876.1; //MeV, deuteron mass
// const Double_t kX1 = 0.0; //0.001*Q_0
// const Double_t kY1 = 0.0; //MeV
// const Double_t kX2[nCores] = {3500, 3320, 2815, 3305,
// 				2280, 2420, 1370, 5480};// Q_i
//   tree->SetBranchAddress("Q_0", &Qcore[0]);// ALAINA
//   tree->SetBranchAddress("Q_1", &Qcore[1]);// BROOKE
//   tree->SetBranchAddress("Q_2", &Qcore[2]);// CINDY
//   tree->SetBranchAddress("Q_3", &Qcore[3]);// JONI
//   tree->SetBranchAddress("Q_4", &Qcore[4]);// KRISTA
//   tree->SetBranchAddress("Q_5", &Qcore[5]);// LINDA
//   tree->SetBranchAddress("Q_6", &Qcore[6]);// ROBERTA
//   tree->SetBranchAddress("Q_7", &Qcore[7]);// SUSAN
//   tree->SetBranchAddress("Q_8", &Qshield[0]);// ANDY
//   tree->SetBranchAddress("Q_9", &Qshield[6]);// REGGIE
//   tree->SetBranchAddress("Q_10", &Qshield[2]);// CARL
//   tree->SetBranchAddress("Q_11", &Qshield[3]);// JACK
//   tree->SetBranchAddress("Q_12", &Qshield[4]);// KEVIN
//   tree->SetBranchAddress("Q_13", &Qshield[5]);// LARRY
//   tree->SetBranchAddress("Q_14", &tof_raw);
//   tree->SetBranchAddress("Q_15", &Qshield[7]);// STEVE
 



void Project_HINDAfwd(){

  const string crNames[3] = {"ALAINA","SUSAN","JONI"};
  const Double_t angles[3] = {40.0, 55.0, 75.0};
  const Double_t shldCut[3] = {4500.0, 3300.0, 3000.0};
  const char *crBranch[3] = {"Q_0", "Q_7", "Q_3"};
  const char *shldBranch[3] = {"Q_8", "Q_15", "Q_11"};
  // const Double_t tofLow[3] = {81.0, 83.0, 91.0};
  // const Double_t tofHigh[3] = {91.0, 98.0, 104.0};
  const Double_t tofLow[3] = {0.0, 0.0, 0.0};
  const Double_t tofHigh[3] = {170.0, 170.0, 170.0};

  //============== for calibration ==========================
  const Double_t E_beam = 65.; //MeV
  const Double_t M_D = 1876.1; //MeV, deuteron mass
  const Double_t kX1 = 0.0; //0.001*Q_0
  const Double_t kY1 = 0.0; //MeV
  const Double_t kX2[3] = {3500., 5480., 3305.};

  Double_t kY2[3],a[3],b[3];
  for(int i=0; i<3; i++){
    kY2[i] = E_beam/(1+(E_beam/M_D)*(1-cos(angles[i]*TMath::Pi()/180.0)));
    a[i] = (kY2[i] - kY1) / (kX2[i] - kX1);
    b[i] = kY2[i] - a[i] * kX2[i];
  }

  // const int index = 0;

  // // TH1F *hInTime = new TH1F(Form("%s_ShldTofCut",crNames.c_str()),
  // // 			   Form("%s_%.0f#circ_shield+tof=[64.5,71]",crNames.c_str(),angles),
  // // 			   200,0,200);
  // // TH1F *hRand = new TH1F(Form("%s_Rand",crNames.c_str()),
  // // 			 Form("%s_%.0f#circ_shield+tof=[100,170]",crNames.c_str(),angles),
  // // 			 200,0,200);
  // // TH1F *hNet;


  TH1F *hcore[3], *hcore_rand[3], *hcore_net[3], *hcore_raw[3], *hcore_shldcut[3];
  TH2F *h2D[3];
  for(int l=0; l<3; l++){
    h2D[l] = new TH2F(Form("%s_2D",crNames[l].c_str()),
		      Form("%s_%.0f#circ",crNames[l].c_str(),angles[l]),
		      500,0,500, 20000,0,20000);
    hcore[l] = new TH1F(Form("h%s",crNames[l].c_str()),
    			Form("%s(%.0f#circ),shield>%.0f,tof=[%.0f,%.0f]",crNames[l].c_str(), angles[l],shldCut[l],tofLow[l],tofHigh[l]),
    			1000,0,500);
    hcore_rand[l] = new TH1F(Form("h%s_rand",crNames[l].c_str()),
			     Form("%s(%.0f#circ),shield>%.0f,tof=[100,170]",crNames[l].c_str(), angles[l],shldCut[l]),
			     1000,0,500);
    hcore_raw[l] = new TH1F(Form("h%s_raw",crNames[l].c_str()),
			    Form("%s(%.0f#circ),raw,tof=[0,170]",crNames[l].c_str()),
			    1000,0,500);
    hcore_shldcut[l] = new TH1F(Form("h%s_shldcut",crNames[l].c_str()),
				Form("%s(%.0f#circ),shield>%.0f",crNames[l].c_str(), angles[l],shldCut[l]),
				1000,0,500);
  }



  // TH1F *h_40 = new TH1F("hALAINA","ALAINA(40#circ),shield>4500,tof=[63,75]",1000,0,500);
  // TH1F *h_40_rand = new TH1F("hALAINA_rand","ALAINA(40#circ),shield>4500,tof=[100,170]",1000,0,500);
  // TH1F *h_40_net;// = new TH1F("hALAINA_net","ALAINA(40#circ),net",1000,0,500);

  // TH1F *h_55 = new TH1F("hSUSAN","SUSAN(55#circ),shield>3300,tof=[63,83]",1000,0,500);
  // TH1F *h_55_rand = new TH1F("hSUSAN_rand","SUSAN(55#circ),shield>3300,tof=[100,170]",1000,0,500);
  // TH1F *h_55_net;// = new TH1F("hSUSAN_net","SUSAN(55#circ),net",1000,0,500);

  // TH1F *h_75 = new TH1F("hJONI","JONI(75#circ),shield>3000,tof=[73,90]",1000,0,500);
  // TH1F *h_75_rand = new TH1F("hJONI_rand","JONI(75#circ),shield>3000,tof=[100,170]",1000,0,500);
  // TH1F *h_75_net;// = new TH1F("hJONI_net","JONI(75#circ),net",1000,0,500);

  // TH1F *h_90 = new TH1F("hROBERTA","ROBERTA(90#circ),shield>3000,tof=[]",1000,0,500);
  // TH1F *h_90_rand = new TH1F("hROBERTA_rand","ROBERTA(90#circ),shield>3000,tof=[100,170]",1000,0,500);
  // TH1F *h_90_net = new TH1F("hROBERTA_net","ROBERTA(90#circ),net",1000,0,500);

  TChain *flat_tree = new TChain("flat_tree");

  ifstream database("./run_database.dat");
  char number[20];
  char name[200];
  string mode;
  string hash = "#";
  float Qcore[3], Qshield[3], Qpaddle;
  float tof_raw, tof, energy;
  while(database >> number >> name >> mode){
    if(number[0] == hash) continue; //the runs starting with # should be neglected 
    if(mode != "deuterium") continue;//choose run mode "deuterium/empty"
    // if(atoi(number) < 180) continue;//start from this run
    // if(atoi(number) > 146) break;//end at this run
    // printf("Current run number: %s  mode: %s\n", number, mode.c_str());

    flat_tree->Add(Form("./root_files/%s",name));

  }//loop over while
  database.close();

  
  flat_tree->SetBranchAddress(crBranch[0], &Qcore[0]);
  flat_tree->SetBranchAddress(shldBranch[0], &Qshield[0]);
  flat_tree->SetBranchAddress(crBranch[1], &Qcore[1]);
  flat_tree->SetBranchAddress(shldBranch[1], &Qshield[1]);
  flat_tree->SetBranchAddress(crBranch[2], &Qcore[2]);
  flat_tree->SetBranchAddress(shldBranch[2], &Qshield[2]);
  flat_tree->SetBranchAddress("Q_14", &tof_raw);
  
  int events = flat_tree->GetEntries();
  for(int i=0; i<events; i++){
    flat_tree->GetEntry(i);
    tof = 0.025*tof_raw*0.5739-9.23582;

    for(int index=0; index<3; index++){    
      energy = a[index]*Qcore[index]+b[index];
      if(energy>20. )//&& Qshield[index]<shldCut[index])
  	h2D[index]->Fill(energy, Qshield[index]);
    }

  }

  
  TCanvas *c1 = new TCanvas();
  h2D[0]->Draw("colz");
  TCanvas *c2 = new TCanvas();
  h2D[1]->Draw("colz");  
  TCanvas *c3 = new TCanvas();
  h2D[2]->Draw("colz");  


  /*
  TString shieldCut_s[3], coreCut_s[3], tofCut_s[3], var[3];
  TCut shieldCut[3], tofCut[3], coreCut[3];
  TString randcut_s = "(0.025*Q_14*0.5739-9.23582)>100&&(0.025*Q_14*0.5739-9.23582)<170";
  TCut randcut = randcut_s;
  Double_t scale[3];
  for(int l=0; l<3; l++){
    var[l].Form("%s*%f", crBranch[l],a[l]); 
    shieldCut_s[l].Form("%s<%f", shldBranch[l],shldCut[l]); shieldCut[l]=shieldCut_s[l];
    tofCut_s[l].Form("(0.025*Q_14*0.5739-9.23582)>%f&&(0.025*Q_14*0.5739-9.23582)<%f",tofLow[l], tofHigh[l]); tofCut[l]=tofCut_s[l];
    coreCut_s[l].Form("%s*%f>10",crBranch[l],a[l]);  coreCut[l]=coreCut_s[l];

    // flat_tree->Project(Form("h%s",crNames[l].c_str()), var[l], shieldCut[l]+tofCut[l]+coreCut[l]);
    // flat_tree->Project(Form("h%s_rand",crNames[l].c_str()), var[l], shieldCut[l]+randcut+coreCut[l]);

    //compare with and without shield cut
    flat_tree->Project(Form("h%s_raw",crNames[l].c_str()), var[l], tofCut[l]+coreCut[l]);
    flat_tree->Project(Form("h%s_shldcut",crNames[l].c_str()), var[l], shieldCut[l]+coreCut[l]);
    flat_tree->Project(Form("h%s_rand",crNames[l].c_str()), var[l], randcut+coreCut[l]);

    //random subtraction
    hcore_net[l] = (TH1F*)hcore_raw[l]->Clone(Form("h%s_randsub",crNames[l].c_str()));
    hcore_net[l]->SetTitle(Form("%s(%.0f#circ),no cut, random subtraction only",crNames[l].c_str(),angles[l]));
    scale[l] = 170./70.;
    hcore_net[l]->Add(hcore_rand[l],-1.*scale[l]);
  }




  // flat_tree->Project("hALAINA","Q_0*0.0195601","Q_8<4500&&(0.025*Q_14*0.5739-9.23582)>64&&(0.025*Q_14*0.5739-9.23582)<73&&Q_0*0.0195601>10");
  // flat_tree->Project("hALAINA_rand","Q_0*0.0195601","Q_8<4500&&(0.025*Q_14*0.5739-9.23582)>100&&(0.025*Q_14*0.5739-9.23582)<170&&Q_0*0.0195601>10");

  // flat_tree->Project("hSUSAN","Q_7*0.0124461","Q_15<3300&&(0.025*Q_14*0.5739-9.23582)>65&&(0.025*Q_14*0.5739-9.23582)<80&&Q_7*0.0124461>10");
  // flat_tree->Project("hSUSAN_rand","Q_7*0.0124461","Q_15<3300&&(0.025*Q_14*0.5739-9.23582)>100&&(0.025*Q_14*0.5739-9.23582)<170&&Q_7*0.0124461>10");

  // flat_tree->Project("hJONI","Q_3*0.0205596","Q_11<3300&&(0.025*Q_14*0.5739-9.23582)>75&&(0.025*Q_14*0.5739-9.23582)<95&&Q_3*0.0205596>10");
  // flat_tree->Project("hJONI_rand","Q_3*0.0205596","Q_11<3300&&(0.025*Q_14*0.5739-9.23582)>100&&(0.025*Q_14*0.5739-9.23582)<170&&Q_3*0.0205596>10");

  // //random subtraction
  // h_40_net = (TH1F*)h_40->Clone("hALAINA_net");
  // h_40_net->SetTitle("ALAINA(40#circ),net");
  // Double_t scale_40 = (71.-64.)/(170.-100.);
  // h_40_net->Add(h_40_rand,-1.*scale_40);

  // h_55_net = (TH1F*)h_55->Clone("hSUSAN_net");
  // h_55_net->SetTitle("SUSAN(55#circ),net");
  // Double_t scale_55 = (90.-65.)/(170.-100.);
  // h_55_net->Add(h_55_rand,-1.*scale_55);

  // h_75_net = (TH1F*)h_75->Clone("hJONI_net");
  // h_75_net->SetTitle("JONI(75#circ),net");
  // Double_t scale_75 = (95.-75.)/(170.-100.);
  // h_75_net->Add(h_75_rand,-1.*scale_75);



  // //================ draw and save hists =================================
  TCanvas *c1 = new TCanvas();
  h2D[0]->Draw("colz");
  TCanvas *c2 = new TCanvas();
  h2D[1]->Draw("colz");  
  TCanvas *c3 = new TCanvas();
  h2D[2]->Draw("colz");  
  // TCanvas *c4 = new TCanvas();
  // hcore[0]->Draw();
  // TCanvas *c5 = new TCanvas();
  // hcore_rand[0]->Draw();
  // TCanvas *c6 = new TCanvas();
  // hcore_net[0]->Draw();

  // TCanvas *c1 = new TCanvas("c1","hALAINA_net");
  // // gPad->SetLogz();
  // hNet->GetXaxis()->SetTitle("E [MeV]");
  // hNet->GetXaxis()->SetTitleSize(0.045);
  // hNet->GetXaxis()->SetLabelSize(0.045);
  // hNet->GetYaxis()->SetTitle("counts / MeV");
  // hNet->GetYaxis()->SetTitleSize(0.045);
  // hNet->GetYaxis()->SetLabelSize(0.045);
  // hNet->Draw();

  // // TCanvas *c2 = new TCanvas("c2","HINDA core");
  // // //gPad->SetLogz();
  // // h2D->GetXaxis()->SetTitle("TOF [ns]");
  // // h2D->GetXaxis()->SetTitleSize(0.05);
  // // h2D->GetXaxis()->SetLabelSize(0.05);
  // // h2D->GetYaxis()->SetTitle("E [MeV]");
  // // h2D->GetYaxis()->SetTitleSize(0.05);
  // // h2D->GetYaxis()->SetLabelSize(0.05);
  // // h2D->Draw("color");

  */




  TFile *outfile = new TFile("hindaFWD.root","update");
  for(int l=0; l<3; l++){
    // hcore[l]->Write();
    // hcore_rand[l]->Write();
    // hcore_net[l]->Write();
    h2D[l]->Write();
    // hcore_raw[l]->Write();
    // hcore_shldcut[l]->Write();
    // hcore_net[l]->Write("",TObject::kOverwrite);

  }

  outfile->Close();

}
