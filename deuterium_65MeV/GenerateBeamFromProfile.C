#include <iostream>
#include <fstream>
#include <cmath>
#include "TROOT.h"
#include "TLine.h"
#include "TAxis.h"
#include "TApplication.h"
#include "Rtypes.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TMatrixD.h"
#include "TEventList.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TError.h"
#include "TMath.h"
#include "TStyle.h"
#include "TRandom.h"

using namespace std;


int main(){
  Double_t energy,dcount;
  Int_t count;
  TH1F *hist = new TH1F("hist","hist",700,0,7e7);
  TRandom rand;
  //  TRandom1 Rx;
  //  TRandom2 Ry;

  ifstream infile("det_resolution/July2016_spectradata/gspectra_65MeV_espd0.25_8mm.txt");
  while(infile >> energy >> dcount){
    count = (Int_t)dcount;
    hist->Fill(energy, count);
  }

  Double_t bincont, bincont_right, bincont_left;
  Int_t bin_right, bin_left;
  Int_t nbins = hist->GetNbinsX();
  Int_t bin_min;
  for(int i=0; i<nbins; i++){
    bincont = hist->GetBinContent(i);
    if(bincont != 0){
      bin_min = i;
      break;
    }  
  }


  for(int i=bin_min; i<nbins; i++){
    bincont = hist->GetBinContent(i);

    if(bincont == 0){
      bin_right = i+1;
      bin_left = i-1;
      bincont_right = hist->GetBinContent(bin_right);
      bincont_left = hist->GetBinContent(bin_left);
      if(bincont_right == 0){
	bin_right = i+2;
	bincont_right = hist->GetBinContent(bin_right);
      }
      if(bincont_left == 0){
	bin_left = i-2;
	bincont_left = hist->GetBinContent(bin_left);
      }
      bincont = (bincont_right-bincont_left)*(i-bin_left)/(bin_right-bin_left)+bincont_left;

      hist->SetBinContent(i, bincont);
    }
  }


  hist->SetLineColor(kBlue);
  hist->Draw();

  Double_t max = hist->GetBinContent(hist->GetMaximumBin());
  Double_t val, prob, E;
  Bool_t found;

  ofstream outfile("det_resolution/July2016_spectradata/beamData_new.dat");
  for(int i=0; i<50000000; i++){
    found = 0;
    while(!found){
      val = rand.Rndm()*7e7;
      prob = rand.Rndm()*max;
      // val = Rx.Rndm()*7e7;
      // prob = Ry.Rndm()*max;
      bincont = hist->GetBinContent(hist->FindBin(val));

      if(prob <= bincont){
	E = val;
	found = true;
      }

    }
  
    outfile << E << endl;

  }
  outfile.close();
  cout<<"File written"<<endl;

  // TH1F *outhist = new TH1F("outh","Beam profile",700,0,7e7);
  // ifstream readfile("beamData.dat");
  // Double_t value;
  // while(readfile >> value)
  //   outhist->Fill(value);



  TH1F *outhist = new TH1F("outh","Beam profile",700,0,7e7);
  Double_t maxSim = outhist->GetBinContent(outhist->GetMaximumBin());
  Double_t scale = max/maxSim;

  outhist->Scale(scale);
  outhist->SetLineColor(kRed);
  outhist->Draw("same");

  return 0;
}
