using namespace std;

void xlEmptySub_paddle(){
  const int nCores = 4;
  //const int nShields = 8;
  //const string crNames[nCores] = {"ALAINA", "BROOKE", "CINDY", "JONI",
//				  "KRISTA", "LINDA", "ROBERTA", "SUSAN"};
  //const string shldNames[nShields] = {"ANDY", "BOBBY", "CARL", "JACK",
//				      "KEVIN", "LARRY", "REGGIE", "STEVE"};
  const int angles[nCores] = {40, 55, 75, 90};

//****change input and output file names here!****
  TString infile_full = "test_hinda_75_randsub_padTOFcut_LD2_new.root";
  TString infile_empty = "test_hinda_75_randsub_padTOFcut_empty_new.root";
  TString outfile = "test_hinda_75_randsub_padTOFcut_net_new.root"; 
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
  TFile *ffull = new TFile(infile_full);
  TFile *fempty = new TFile(infile_empty);


  //~~~~~~~~~~~~~~~~~~~~ for single paddle use ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  TH1F *hfull, *hempty, *hempty_scaled, *hsub;
  TH1F *hfull2, *hempty2, *hempty_scaled2, *hsub2;
  const Int_t rebin = 10;

  hfull = (TH1F*)ffull->Get("JONI_TofShldCut_net")->Clone("JONI_TofShldCut_net_full");
  hempty = (TH1F*)fempty->Get("JONI_TofShldCut_net");
  hempty_scaled = (TH1F*)hempty->Clone("JONI_TofShldCut_netempty_scaled");
  hempty_scaled->Scale(scale_factor);
  hsub = (TH1F*)hfull->Clone("JONI_TofShldCut_net_sub");
  hsub->Add(hempty, -1.0*scale_factor);
  hsub->Rebin(rebin);

  hfull2 = (TH1F*)ffull->Get("JONI_TofShldPadtofCut_net")->Clone("JONI_TofShldPadtofCut_net_full");
  hempty2 = (TH1F*)fempty->Get("JONI_TofShldPadtofCut_net");
  hempty_scaled2 = (TH1F*)hempty2->Clone("JONI_TofShldPadtofCut_net_empty_scaled");
  hempty_scaled2->Scale(scale_factor);
  hsub2 = (TH1F*)hfull2->Clone("JONI_TofShldPadtofCut_net_sub");
  hsub2->Add(hempty2, -1.0*scale_factor);
  hsub2->Rebin(rebin);

  TCanvas *c1 = new TCanvas("c1","tof+shld cut net (with rand sub)");
  //gPad->SetLogz();
  hsub->GetXaxis()->SetRangeUser(20,180);
  hsub->GetXaxis()->SetTitle("E [MeV]");
  hsub->GetXaxis()->SetLabelSize(0.045);
  hsub->GetXaxis()->SetTitleSize(0.045);
  hsub->GetYaxis()->SetRangeUser(0,550);
  hsub->GetYaxis()->SetTitle("counts");
  hsub->GetYaxis()->SetTitleSize(0.045);
  hsub->GetYaxis()->SetLabelSize(0.045);
  hsub->SetTitle("JONI_75#circ_tof+shld_cut_random_sub");
  hsub->SetLineColor(kBlack);
  hsub->Draw();
  hsub2->SetFillStyle(3003);
  hsub2->SetFillColor(kRed+1);
  hsub2->SetLineColor(kRed+1);
  hsub2->Draw("same");

  TLegend *l1 = new TLegend(0.5,0.6,0.85,0.75);
  l1->AddEntry(hsub,"no paddle cut");
  l1->AddEntry(hsub2,"with paddle TOF cut");
  l1->Draw();



  // TH2F *hfull, *hempty, *hempty_scaled, *hsub;
  // hfull = (TH2F*)ffull->Get("hinda_JONI_Net")->Clone("hinda_75deg_ld2");
  // hempty = (TH2F*)fempty->Get("hinda_JONI_Net");
  // hempty_scaled = (TH2F*)hempty->Clone("hinda_75deg_empty_scaled");
  // hempty_scaled->Scale(scale_factor);
  // hsub = (TH2F*)hfull->Clone("hinda_75deg_net");
  // hsub->Add(hempty, -1.0*scale_factor);
  // hsub->SetTitle("hinda_JONI_75#circ_TOF+Shield_cut_w/_Random+Empty_sub");
  // hsub->GetXaxis()->SetTitle("paddle TOF");
  // hsub->GetXaxis()->SetTitleSize(0.045);
  // hsub->GetYaxis()->SetTitle("E [MeV]");
  // hsub->GetYaxis()->SetTitleSize(0.045);
  // hsub->GetYaxis()->SetLabelSize(0.045);
  // TCanvas *c1 = new TCanvas("c1","tof+shld cut & rndm+empty sub");
  // gPad->SetLogz();
  // hsub->Draw("colz");




  TFile *fout = new TFile(outfile, "recreate");
  hfull->Write();
  hempty->Write();
  hempty_scaled->Write();
  hsub->Write();
  hfull2->Write();
  hempty2->Write();
  hempty_scaled2->Write();
  hsub2->Write();
  c1->Write();
  fout->Close();



  //~~~~~~~~~~~~~~~~~~~~ for all 4 paddles use ~~~~~~~~~~~~~~~~~~~~
  // TH2F *hfull[nCores]; 
  // TH2F *hempty[nCores]; 
  // TH2F *hempty_scaled[nCores]; 
  // TH2F *hsub[nCores];
  // //TLine *zore[nCores];
  // TCanvas *canvas1 = new TCanvas("c1","net");
  // canvas1->Divide(2,2);
  // TCanvas *canvas2 = new TCanvas("c2","full");
  // canvas2->Divide(2,2);  
  // TCanvas *canvas3 = new TCanvas("c3","empty scaled");
  // canvas3->Divide(2,2);  
  // TCanvas *canvas4 = new TCanvas("c4","empty unscaled");
  // canvas4->Divide(2,2);

  // for(int i=0; i<nCores; i++){
  //   hfull[i] = (TH2F*)ffull->Get(Form("%i#circ",angles[i]));
  //   hempty[i] = (TH2F*)fempty->Get(Form("%i#circ",angles[i]));
  //   hsub[i] = (TH2F*)hfull[i]->Clone(Form("%i#circ_SUB",angles[i]));
  //   hempty_scaled[i] = (TH2F*)hempty[i]->Clone(Form("%i#circ_scaled",angles[i]));
  //   hempty_scaled[i]->Scale(scale_factor);
  //   hsub[i]->Add(hempty[i], -1.0*scale_factor);

  //   canvas1->cd(i+1);
  //   gPad->SetLogz();
  //   //hfull[i]->SetLineColor(kBlue);
  //   //hfull[i]->GetXaxis()->SetTitle("Energy [MeV]");
  //   //hfull[i]->GetXaxis()->SetTitleSize(0.045);
  //   //hfull[i]->GetYaxis()->SetLabelSize(0.045);
  //   //hfull[i]->GetXaxis()->SetLabelSize(0.045);
  //   //hfull[i]->Draw("col");
  //   //hsub[i]->SetLineColor(kRed);
  //   hsub[i]->Draw("col");

  //   canvas2->cd(i+1);
  //   gPad->SetLogz();
  //   hfull[i]->Draw("col");

  //   canvas3->cd(i+1);
  //   gPad->SetLogz();
  //   hempty_scaled[i]->Draw("col");

  //   canvas4->cd(i+1);
  //   gPad->SetLogz();
  //   hempty[i]->Draw("col");
  //   // zore[i] = new TLine(25.0, 0.0, 100.0, 0.0);
  //   // zore[i]->SetLineWidth(2);
  //   // zore[i]->SetLineStyle(4);
  //   // zore[i]->SetLineColor(kViolet);
  //   // zore[i]->Draw("same");
  // }




  // TFile *fout = new TFile(outfile,"recreate");
  // // TFile *fout = new TFile("timeOfFlight.root","recreate");
  // for(int i=0; i<nCores; i++){
  //   hfull[i]->Write();
  //   hempty[i]->Write(); 
  //   hsub[i]->Write();
  // }
  // canvas1->Write();
  // canvas2->Write();
  // canvas3->Write();

  // fout->Close();

}
