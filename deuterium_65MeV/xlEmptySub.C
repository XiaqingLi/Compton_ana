using namespace std;

void xlEmptySub(){
  const int nCores = 8;
  const int nShields = 8;
  const string crNames[nCores] = {"ALAINA", "BROOKE", "CINDY", "JONI",
				  "KRISTA", "LINDA", "ROBERTA", "SUSAN"};
  const string shldNames[nShields] = {"ANDY", "BOBBY", "CARL", "JACK",
				      "KEVIN", "LARRY", "REGGIE", "STEVE"};
  const double angles[nCores] = {40.0, 159.0, 125.0, 75.0, 
				 145.0, 110.0, 90.0, 55.0};

//****change input and output file names here!****
  TString infile_full = "xlFull.root";
  TString infile_empty = "xlEmpty.root";
  TString outfile = "xlNet.root"; 
//************************************************


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

  Int_t clock,five_pad,veto_pad,prescale_BPM;
  Int_t full_clock = 0;
  Int_t empty_clock = 0;
  Double_t full_5pad = 0.0;
  Double_t empty_5pad = 0.0;
  Double_t full_veto_pad = 0.0;
  Double_t empty_veto_pad = 0.0;
  Double_t full_BPM = 0.0;
  Double_t empty_BPM = 0.0;

  FILE *scaler = fopen("compton_scalers","r");
  ifstream database("run_database.dat");
  while(database >> number >> name >> mode){
    if(number[0] == hash) continue; //the runs starting with # should be neglected 
    if(mode != "deuterium") continue;//choose run mode
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
	      full_clock += tmp;
	      printf("clock = %d\n",full_clock);
	      continue;
	    }//if find clock
	    if(line[0] == e){
	      printf(line);
	      sscanf(line, "%i %i", &num, &tmp);
	      full_BPM += (Double_t)tmp;
	      printf("BPM = %.1f\n",full_BPM);
	      continue;
	    }//if find BPM
	    else if(line[0] == f){
	      printf(line);
	      sscanf(line, "%i %i", &num, &tmp);
	      full_5pad += (Double_t)tmp;
	      printf("5pad = %.1f\n",full_5pad);
	      continue;
	    }//if find 5pad
	    else if(line[0] == g){
	      printf(line);
	      sscanf(line, "%i %i", &num, &tmp);
	      full_veto_pad += (Double_t)tmp;
	      printf("veto_pad = %.1f\n",full_veto_pad);
	      break;
	    }//if find BPM
	  }
	  break;
	}//if find good run
      }//if line starts with "Run"
    }
  }//while loop database
  fclose(scaler);
  database.close();
  full_BPM = full_BPM*10000.0;


  FILE *scaler = fopen("compton_scalers","r");
  ifstream database("run_database.dat");
  while(database >> number >> name >> mode){
    if(number[0] == hash) continue; //the runs starting with # should be neglected 
    if(mode != "empty") continue;//choose run mode
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
	      empty_clock += tmp;
	      printf("clock = %d\n",empty_clock);
	      continue;
	    }//if find clock
	    if(line[0] == e){
	      printf(line);
	      sscanf(line, "%i %i", &num, &tmp);
	      empty_BPM += (Double_t)tmp;
	      printf("BPM = %.1f\n",empty_BPM);
	      continue;
	    }//if find BPM
	    else if(line[0] == f){
	      printf(line);
	      sscanf(line, "%i %i", &num, &tmp);
	      empty_5pad += (Double_t)tmp;
	      printf("5pad = %.1f\n",empty_5pad);
	      continue;
	    }//if find 5pad
	    else if(line[0] == g){
	      printf(line);
	      sscanf(line, "%i %i", &num, &tmp);
	      empty_veto_pad += (Double_t)tmp;
	      printf("veto_pad = %.1f\n",empty_veto_pad);
	      break;
	    }//if find BPM
	  }
	  break;
	}//if find good run
      }//if line starts with "Run"
    }
  }//while loop database
  fclose(scaler);
  database.close();
  empty_BPM = empty_BPM*10000.0;



  // flux calculation code from Rob
  gROOT->ProcessLine(".L fluxcal.C");
  Double_t scale_factor = (Double_t)full_5pad / empty_5pad;
  cout << "5pad: full / empty = " << scale_factor << endl;
  scale_factor = (Double_t)full_clock / empty_clock;
  cout << "clock: full / empty = " << scale_factor << endl;
  Float_t n_full,n_full_err;
  Float_t n_empty,n_empty_err;
  const Float_t kFactor = 55.17;// FINAL value from Rob!
  const Float_t kFactorError = 0.34;// FINAL value from Rob!
  fluxcal(kFactor, kFactorError, full_5pad, full_veto_pad, full_BPM, &n_full, &n_full_err);
  fluxcal(kFactor, kFactorError, empty_5pad, empty_veto_pad, empty_BPM, &n_empty, &n_empty_err);
  scale_factor = n_full / n_empty;// scale by number of photons
  cout << "n_full / n_empty = " << scale_factor << endl;
  cout << n_full/(0.1*full_clock)<<endl;

  //print the total number of photons of all the full runs
  cout<<"N_gamma = "<<n_full<<endl;








  //Do empty subtraction
  // TFile *ffull = new TFile("full.root");
  // TFile *fempty = new TFile("empty.root");
  TFile *ffull = new TFile(infile_full);
  TFile *fempty = new TFile(infile_empty);

  TH1F *hfull[nCores]; 
  TH1F *hempty[nCores]; 
  TH1F *hsub[nCores];
  TLine *zore[nCores];
  TCanvas *canvas = new TCanvas();
  canvas->Divide(4,2);
  for(int i=0; i<nCores; i++){
    hfull[i] = (TH1F*)ffull->Get(Form("%s_Net",crNames[i].c_str()));
    hempty[i] = (TH1F*)fempty->Get(Form("%s_Net",crNames[i].c_str()));
    hsub[i] = (TH1F*)hfull[i]->Clone(Form("%s_NetSUB",crNames[i].c_str()));
    hsub[i]->Add(hempty[i], -1.0*scale_factor);
    //hfull[i]->Rebin(2);
 /*   
    if(i == 0){
      hfull[i]->GetXaxis()->SetRangeUser(40,100);
      hfull[i]->GetYaxis()->SetRangeUser(-10,300);
    }
    else if(i == 1){
      hfull[i]->GetXaxis()->SetRangeUser(25,100);
      hfull[i]->GetYaxis()->SetRangeUser(-20,300);
    }
    else if(i == 7){
      hfull[i]->GetXaxis()->SetRangeUser(25,100);
      hfull[i]->GetYaxis()->SetRangeUser(-10,300);
    }
    else{
      hfull[i]->GetXaxis()->SetRangeUser(25,100);
      hfull[i]->GetYaxis()->SetRangeUser(-10,250);
    }
*/

    canvas->cd(i+1);
    hfull[i]->SetLineColor(kBlue);
    hfull[i]->GetXaxis()->SetTitle("Energy [MeV]");
    hfull[i]->GetXaxis()->SetTitleSize(0.045);
    hfull[i]->GetYaxis()->SetLabelSize(0.045);
    hfull[i]->GetXaxis()->SetLabelSize(0.045);
    hfull[i]->Draw();
    hsub[i]->SetLineColor(kRed);
    hsub[i]->Draw("same");

    // zore[i] = new TLine(25.0, 0.0, 100.0, 0.0);
    // zore[i]->SetLineWidth(2);
    // zore[i]->SetLineStyle(4);
    // zore[i]->SetLineColor(kViolet);
    // zore[i]->Draw("same");
  }

  TFile *fout = new TFile(outfile,"recreate");
  // TFile *fout = new TFile("timeOfFlight.root","recreate");
  for(int i=0; i<nCores; i++){
    hfull[i]->Write();
    hsub[i]->Write();
  }
  canvas->Write();
  fout->Close();

}
