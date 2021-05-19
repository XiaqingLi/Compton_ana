{
  const int nCores = 8;
  const int nShields = 8;
  const string crNames[nCores] = {"ALAINA", "BROOKE", "CINDY", "JONI",
				  "KRISTA", "LINDA", "ROBERTA", "SUSAN"};
  const string shldNames[nShields] = {"ANDY", "BOBBY", "CARL", "JACK",
				      "KEVIN", "LARRY", "REGGIE", "STEVE"};
  const double angles[nCores] = {55.0, 90.0, 125.0, 125.0, 
				 125.0, 90.0, 90.0, 55.0};
  const string pos[nCores] = {"D", "D", "L", "D", "R", "R", "L", "L"};

  //**** change input and output file names here! ****
  TString infile_full = "/var/phy/project/mepg/xl79/helium_84MeV/Full.root";
  TString infile_empty = "/var/phy/project/mepg/xl79/helium_84MeV/Empty.root";
  TString outfile = "Net.root";
  //**************************************************


  //calculate scale factor
  char line[200];
  char str[100];
  Int_t run_num, num, tmp;
  string R = "R";
  string o = "0";
  string e = "5";
  string f = "6";
  string g = "7";

  char number[20];
  char name[200];
  string mode;
  string hash = "#";
  Int_t helicity;

  Int_t clock;
  Double_t five_pad,veto_pad,BPM;
  Float_t n,n_err;

  Float_t n_full = 0.0, n_full_CINDY = 0.0;
  Float_t n_empty = 0.0;
  Float_t n_full_err = 0.0, n_full_err_CINDY = 0.0;
  Float_t n_empty_err = 0.0;
  Int_t clock_full = 0, clock_full_CINDY = 0;
  Int_t clock_empty = 0;

  Float_t kFactor;
  Float_t kFactorError;
  Double_t scale_factor, scale_factor_CINDY;
  Float_t flux, flux_err;
  ofstream flux_out("flux.dat"); //output file for flux information


	
  //process full target runs
  FILE *scaler = fopen("/var/phy/project/mepg/xl79/helium_84MeV/compton_scalers","r");
  ifstream database("/var/phy/project/mepg/xl79/helium_84MeV/run_database.dat");

  while(database >> number >> name >> mode >> helicity){
    if(number[0] == hash) continue; //the runs starting with # should be neglected 
    if(mode != "helium_full") continue; //choose run mode
    printf("\n");
    printf("Run#%d\n",atoi(number));
 
    while(!feof(scaler)){
      fgets(line, 200, scaler);
      if(line[0] == R){
	sscanf(line, "%s %i", &str, &run_num);
	if(run_num == atoi(number)){
	  printf(line);
	  while(!feof(scaler)){
	    fgets(line, 200, scaler);
	    if(line[0] == o){
	      printf(line);
	      sscanf(line, "%i %i", &num, &tmp);
	      clock = tmp;
	      printf("clock = %d\n",clock);
	      continue;
	    }//if find clock
	    if(line[0] == e){
	      printf(line);
	      sscanf(line, "%i %i", &num, &tmp);
	      BPM = (Double_t)tmp*10000.0;
	      printf("BPM = %.1f\n",BPM);
	      continue;
	    }//if find BPM
	    else if(line[0] == f){
	      printf(line);
	      sscanf(line, "%i %i", &num, &tmp);
	      five_pad = (Double_t)tmp;
	      printf("5pad = %.1f\n",five_pad);
	      continue;
	    }//if find 5pad
	    else if(line[0] == g){
	      printf(line);
	      sscanf(line, "%i %i", &num, &tmp);
	      veto_pad = (Double_t)tmp;
	      printf("veto_pad = %.1f\n",veto_pad);
	      break;
	    }//if find veto pad
	  }
	  break;
	}//if find good run
      }//if line starts with "Run"
    }//while loop scaler

	  
    //calculate flux for every run and add up n_full
    gROOT->ProcessLine(".L /var/phy/project/mepg/xl79/helium_84MeV/fluxcal.C");
    if(atoi(number)>=928 && atoi(number)<=951) {
      kFactor = 214.66;
      kFactorError = 3.16;
    }else if(atoi(number)>=952 && atoi(number)<1049) {
      kFactor = 126.31;
      kFactorError = 0.66;
    }else if(atoi(number)>=1292) {
      Factor = 101.44;
      FactorError = 2.74;
    }

    fluxcal(kFactor, kFactorError, five_pad, veto_pad, BPM, &n, &n_err);
    if(atoi(number)>951){
      n_full_CINDY += n;
      clock_full_CINDY += clock;
      n_full_err_CINDY += n_err*n_err;
    }
    n_full += n;
    clock_full += clock;
    n_full_err += n_err*n_err;

    
    flux = n / ((Float_t)clock*0.1);
    flux_err = n_err / ((Float_t)clock*0.1);
    flux_out << atoi(number) <<"\t"<< (int)flux <<"\t"<< (int)flux_err << endl;
  }//while loop database

  fclose(scaler);
  database.close();
  n_full_err = sqrt(n_full_err);
  n_full_err_CINDY = sqrt(n_full_err_CINDY);



  //process empty target runs
  FILE *scaler = fopen("/var/phy/project/mepg/xl79/helium_84MeV/compton_scalers","r");
  ifstream database("/var/phy/project/mepg/xl79/helium_84MeV/run_database.dat");
  while(database >> number >> name >> mode >> helicity){
    if(number[0] == hash) continue; //the runs starting with # should be neglected 
    if(mode != "helium_empty") continue;//choose run mode
    printf("\n");
    printf("Run#%d\n",atoi(number));
 
    while(!feof(scaler)){
      fgets(line, 200, scaler);
      if(line[0] == R){
	sscanf(line, "%s %i", &str, &run_num);
	if(run_num == atoi(number)){
	  printf(line);
	  while(!feof(scaler)){
	    fgets(line, 200, scaler);
	    if(line[0] == o){
	      printf(line);
	      sscanf(line, "%i %i", &num, &tmp);
	      clock = tmp;
	      printf("clock = %d\n",clock);
	      continue;
	    }//if find clock
	    if(line[0] == e){
	      printf(line);
	      sscanf(line, "%i %i", &num, &tmp);
	      BPM = (Double_t)tmp*10000.0;
	      printf("BPM = %.1f\n",BPM);
	      continue;
	    }//if find BPM
	    else if(line[0] == f){
	      printf(line);
	      sscanf(line, "%i %i", &num, &tmp);
	      five_pad = (Double_t)tmp;
	      printf("5pad = %.1f\n",five_pad);
	      continue;
	    }//if find 5pad
	    else if(line[0] == g){
	      printf(line);
	      sscanf(line, "%i %i", &num, &tmp);
	      veto_pad = (Double_t)tmp;
	      printf("veto_pad = %.1f\n",veto_pad);
	      break;
	    }//if find BPM
	  }
	  break;
	}//if find good run
      }//if line starts with "Run"
    }//while loop scaler

	  
    //calculate flux for every run and add up n_empty
    gROOT->ProcessLine(".L /var/phy/project/mepg/xl79/helium_84MeV/fluxcal.C");
    if(atoi(number)>=928 && atoi(number)<=951) {
      kFactor = 214.66;
      kFactorError = 3.16;
    }else if(atoi(number)>=952 && atoi(number)<1049) {
      kFactor = 126.31;
      kFactorError = 0.66;
    }else if(atoi(number)>=1292) {
      Factor = 101.44;
      FactorError = 2.74;
    }
    fluxcal(kFactor, kFactorError, five_pad, veto_pad, BPM, &n, &n_err);
    n_empty += n;
    clock_empty += clock;
    n_empty_err += n_err*n_err;

    
    flux = n / ((Float_t)clock*0.1);
    flux_err = n_err / ((Float_t)clock*0.1);
    flux_out << atoi(number) <<"\t"<< (int)flux <<"\t"<< (int)flux_err << endl;
    
  }//while loop database

  fclose(scaler);
  database.close();
  n_empty_err = sqrt(n_empty_err);
  flux_out.close();


	
  // //flux calculation code from Rob
  // gROOT->ProcessLine(".L fluxcal.C");
  // Double_t scale_factor = (Double_t)full_5pad / empty_5pad;
  // cout << "5pad: full / empty = " << scale_factor << endl;
  // scale_factor = (Double_t)full_clock / empty_clock;
  // cout << "clock: full / empty = " << scale_factor << endl;
  // Float_t n_full,n_full_err;
  // Float_t n_empty,n_empty_err;
  // const Float_t kFactor = 55.17;// FINAL value from Rob!
  // const Float_t kFactorError = 0.34;// FINAL value from Rob!
  // fluxcal(kFactor, kFactorError, full_5pad, full_veto_pad, full_BPM, &n_full, &n_full_err);
  // fluxcal(kFactor, kFactorError, empty_5pad, empty_veto_pad, empty_BPM, &n_empty, &n_empty_err);
  scale_factor = n_full / n_empty;// scale by number of photons
  cout << "n_full / n_empty = " << scale_factor << endl;
  // cout << n_full/(0.1*full_clock)<<endl;

  //print the total number of photons of all the full runs
  cout<<"N_gamma = "<<n_full<<endl;
  cout<<"sigma(N_gamma) = "<<n_full_err<<endl;

  cout<<"full clock = "<<clock_full<<endl;
  cout<<"empty clock = "<<clock_empty<<endl;
  cout<<"full clock /empty clock = "<<(double)clock_full/clock_empty<<endl;

  cout<<"\naveraged flux = "<<(n_full+n_empty)/((double)(clock_full+clock_empty)*0.1)<<endl;

  //****************special for CINDY**********************
  scale_factor_CINDY = n_full_CINDY / n_empty;
  cout << "\nn_full_CINDY / n_empty = " << scale_factor_CINDY << endl;
  cout<<"N_gamma_CINDY = "<<n_full_CINDY<<endl;
  cout<<"sigma(N_gamma_CINDY) = "<<n_full_err_CINDY<<endl;
  //*******************************************************


	
  //subtract empty target contributions from full target data
  TFile *ffull = new TFile(infile_full);
  TFile *fempty = new TFile(infile_empty);

  int rebin = 1;

  TH1F *hfull[nCores]; 
  TH1F *hempty[nCores]; 
  TH1F *hsub[nCores];
  TLine *zero[nCores];
  TCanvas *c1 = new TCanvas("c1","c1");
  c1->Divide(4,2);
  TCanvas *c2 = new TCanvas("c2","c2");
  c2->Divide(4,2);
  for(int i=0; i<nCores; i++){
    hfull[i] = (TH1F*)((TH1F*)ffull->Get(Form("%s_net",crNames[i].c_str())))->Clone(Form("%s_full",crNames[i].c_str()));
    hempty[i] = (TH1F*)((TH1F*)fempty->Get(Form("%s_net",crNames[i].c_str())))->Clone(Form("%s_empty",crNames[i].c_str()));
    hsub[i] = (TH1F*)hfull[i]->Clone(Form("%s_Net",crNames[i].c_str()));
    hsub[i]->SetTitle(Form("%s #theta=%0.f %s",crNames[i].c_str(),angles[i],pos[i].c_str()));

    if(i==2)
      hsub[i]->Add(hempty[i], -1.0*scale_factor_CINDY);
    else
      hsub[i]->Add(hempty[i], -1.0*scale_factor);
 
    c1->cd(i+1);
    hfull[i]->SetLineColor(kBlack);
    hfull[i]->GetXaxis()->SetTitle("Energy [MeV]");
    hfull[i]->GetXaxis()->SetTitleSize(0.045);
    hfull[i]->GetYaxis()->SetLabelSize(0.045);
    hfull[i]->GetXaxis()->SetLabelSize(0.045);
    hfull[i]->GetXaxis()->SetRangeUser(20,200);
    // hfull[i]->GetYaxis()->SetRangeUser(-10,200);
    hfull[i]->Rebin(rebin);
    hfull[i]->Draw();
    hsub[i]->Rebin(rebin);
    hsub[i]->SetLineColor(kRed);
    hsub[i]->Draw("same");

    zero[i] = new TLine(20,0,200,0);
    // zero[i]->SetLineWidth(2);
    zero[i]->SetLineStyle(2);
    zero[i]->SetLineColor(kBlack);
    zero[i]->Draw("same");

    c2->cd(i+1);
    hfull[i]->SetLineColor(kBlack);
    hfull[i]->GetXaxis()->SetTitle("Energy [MeV]");
    hfull[i]->GetXaxis()->SetTitleSize(0.045);
    hfull[i]->GetYaxis()->SetLabelSize(0.045);
    hfull[i]->GetXaxis()->SetLabelSize(0.045);
    hfull[i]->GetXaxis()->SetRangeUser(20,200);
    // hfull[i]->GetYaxis()->SetRangeUser(-10,200);
    hfull[i]->Draw();
    hempty[i]->Scale(scale_factor);
    hempty[i]->Rebin(rebin);
    hempty[i]->SetLineColor(kGreen);
    hempty[i]->Draw("same");

    zero[i] = new TLine(20,0,200,0);
    // zero[i]->SetLineWidth(2);
    zero[i]->SetLineStyle(2);
    zero[i]->SetLineColor(kBlack);
    zero[i]->Draw("same");
  }



  TFile *fout = new TFile(outfile,"recreate");
  // TFile *fout = new TFile("timeOfFlight.root","recreate");
  for(int i=0; i<nCores; i++){
    hfull[i]->Write();
    hempty[i]->Write();
    hsub[i]->Write();
  }
  c1->Write();
  c2->Write();
  fout->Close();



}
