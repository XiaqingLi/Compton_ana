using namespace std;

void xlReplay(){

  const int nCores = 8;
  const int nShields = 8;
  const string crNames[nCores] = {"ALAINA", "BROOKE", "CINDY", "JONI",
				  "KRISTA", "LINDA", "ROBERTA", "SUSAN"};
  const string shldNames[nShields] = {"ANDY", "BOBBY", "CARL", "JACK",
				      "KEVIN", "LARRY", "REGGIE", "STEVE"};
  const Double_t angles[nCores] = {40.0, 159.0, 125.0, 75.0, 
				   145.0, 110.0, 90.0, 55.0};

  const Double_t cutLow[nCores] = {82.0, 64.0, 91.0, 91.5,
                                   87.5, 78.5, 87.0, 83.0};
  const Double_t cutHigh[nCores] = {86.0, 69.0, 97.0, 95.5,
                                    93.5, 82.5, 91.0, 86.5};
  const Double_t shldCut[nCores] = {4500.0, 0.0, 3000.0, 3000.0,
				    4300.0, 3000.0, 2700.0, 3300.0};
 

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
    hcore_tof_shld[i] = new TH1F(Form("%s_TOFShldCut",crNames[i].c_str()),
                                 Form("%s_%.0f#circ",crNames[i].c_str(),angles[i]),
                                 2000,0,200);
    hcore_tof_shld[i]->Sumw2();
    hcore_tof_shldR[i] = new TH1F(Form("%s_TOFShldCut_Raw",crNames[i].c_str()),
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


  ifstream database("run_database.dat");
  char number[20];
  char name[200];
  string mode;
  string hash = "#";
  float Qcore[nCores];
  float Qshield[nCores];
  float tof_raw, tof, energy;
  while(database >> number >> name >> mode){
    if(number[0] == hash) continue; //the runs starting with # should be neglected 
    if(mode != "deuterium") continue;//choose run mode
    //if(atoi(number) < 180) continue;//start from this run
    //if(atoi(number) > 180) break;//end at this run
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
    tree->SetBranchAddress("Q_9", &Qshield[6]);// REGGIE
    tree->SetBranchAddress("Q_10", &Qshield[2]);// CARL
    tree->SetBranchAddress("Q_11", &Qshield[3]);// JACK
    tree->SetBranchAddress("Q_12", &Qshield[4]);// KEVIN
    tree->SetBranchAddress("Q_13", &Qshield[5]);// LARRY
    tree->SetBranchAddress("Q_14", &tof_raw);
    tree->SetBranchAddress("Q_15", &Qshield[7]);// STEVE
   
 

    int events = tree->GetEntries();
    for(int i=0; i<events; i++){
      tree->GetEntry(i);
      tof = 0.025*tof_raw*0.5739-9.23582;
      Qshield[1] = -10.0; //for BROOKE

      for(int j=0; j<nCores; j++){
	hcore_raw[j]->Fill(energy);

	if(tof>cutLow[j] && tof<cutHigh[j] && Qshield[j]<shldCut[j]){
	  hcore_tof_shldR[j]->Fill(Qcore[j]); //to find kX2
	}

      }//loop over for(j<nCores)

    }//loop over for(i<events)
    file->Close();

  }//loop over while full/empty runs
  database.close();


  TFile *xlfile = new TFile("xlFull.root","recreate");
  for(int i=0; i<nCores; i++){
    hcore_raw[i]->Write();
    hcore_tof_shldR[i]->Write();
  }





  xlfile->Close();
  printf("file written\n");








}
