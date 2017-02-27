#include <iostream>
#include <fstream>
#include <cmath>

#include "TROOT.h"
#include "TLine.h"
#include "TText.h"
#include "TMathText.h"
#include "TLatex.h"
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

using namespace std;


//****change input and output file names here!****
TString infile = "xlNet.root";
TString outfile = "xlFitCorePrime.root";
//************************************************

const Int_t kNcores = 8;
const TString coreName[kNcores] = {"ALAINA", "SUSAN", "JONI", "ROBERTA",
				   "LINDA", "CINDY", "KRISTA", "BROOKE"};
const Double_t Efit_low[kNcores] = {35.0, 35.0, 35.0, 35.0,
				    45.0, 45.0, 45.0, 45.0};
const Double_t Efit_high[kNcores] = {74.0, 74.0, 74.0, 74.0,
				     71.0, 71.0, 71.0, 71.0};
const Double_t E_beam = 65.0;
const Int_t rebin = 12;
const Double_t drawlow = 30;
const Double_t drawhigh = 90;

TH1F *hgeant;
TH1F *hdata; 
TH1F *hfit;
TH1F *hist_data[kNcores];
TH1F *hist_fit[kNcores];
TH1F *hconv[kNcores];
Int_t nbins;
Double_t xlow,xhigh;
TF1 *exp_func;
TF1 *fexp[kNcores];


void GetData(Int_t index);
Double_t GetFitVal(TH1F *geant_hist, const Double_t *par, Int_t j);
Double_t DoConvolution(TH1F *geant_hist, const Double_t *par, Int_t j);
Double_t GetChi2_fwd(const Double_t *par);
Double_t GetChi2_bwd(const Double_t *par);

int main(){
  for(int k=0; k<kNcores; k++){
    GetData(k);

    //****************setup Minuit*****************
    ROOT::Math::Minimizer *minuit = ROOT::Math::Factory::CreateMinimizer("Minuit","");
    minuit->SetMaxFunctionCalls(10000); 
    minuit->SetTolerance(0.1);
    minuit->SetPrintLevel(1);
      
    if(k<4){ //fit for forward angles
      ROOT::Math::Functor f(&GetChi2_fwd,6);//Fitting function
      minuit->SetFunction(f);    
 
      if(coreName[k] == "ALAINA"){   
	minuit->SetVariable(0,"scale", 3.40504e-02, 0.0001);
	minuit->SetVariable(1,"sigma", 1.43309e-01, 0.0001);
	minuit->SetVariable(2,"mean", 1.38267e+00, 0.0001);
	minuit->SetVariable(3,"exp_scale", 1.68444e+00, 0.0001);  
	minuit->SetVariable(4,"exp_decay", 1.83201e-03, 0.0001);
	minuit->SetFixedVariable(5,"index",(Double_t)k);
      } 
      else if(coreName[k] == "SUSAN"){   
	minuit->SetVariable(0,"scale", 2.87215e-02, 0.0001);
	minuit->SetVariable(1,"sigma", 1.31849e-01, 0.0001);
	minuit->SetVariable(2,"mean", 3.26263e+00, 0.0001);
	minuit->SetVariable(3,"exp_scale", 2.30166e-01, 0.0001);  
	minuit->SetVariable(4,"exp_decay", 1.11768e-03, 0.0001);
	minuit->SetFixedVariable(5,"index",(Double_t)k);
      } 
      else if(coreName[k] == "JONI"){   
	minuit->SetVariable(0,"scale", 2.80960e-02, 0.0001);
	minuit->SetVariable(1,"sigma", 1.64846e-01, 0.0001);
	minuit->SetVariable(2,"mean", 1.62840e+00, 0.0001);
	minuit->SetVariable(3,"exp_scale", 1.05347e-01, 0.0001);  
	minuit->SetVariable(4,"exp_decay", 1.76686e-03, 0.0001);
	minuit->SetFixedVariable(5,"index",(Double_t)k);
      } 
      else if(coreName[k] == "ROBERTA"){   
	minuit->SetVariable(0,"scale", 2.52297e-02, 0.0001);
	minuit->SetVariable(1,"sigma", 1.84266e-01, 0.0001);
	minuit->SetVariable(2,"mean", 2.97019e+00, 0.0001);
	minuit->SetVariable(3,"exp_scale", 4.22329e-02, 0.0001);  
	minuit->SetVariable(4,"exp_decay", 1.43182e-03, 0.0001);
	minuit->SetFixedVariable(5,"index",(Double_t)k);
      }  
      
      minuit->Minimize(); //has to be in the same layer as ROOT::Math::Functor f()

    }//if forward angles
    


    else{  //fit for backward angles
      ROOT::Math::Functor f(&GetChi2_bwd,4);//Fitting function
      minuit->SetFunction(f);

      if(coreName[k] == "LINDA"){   
    	minuit->SetVariable(0,"scale", 2.50040e-02, 0.0001);
    	minuit->SetVariable(1,"sigma", 1.98657e-01, 0.0001);
    	minuit->SetVariable(2,"mean", 0, 0.0001);
    	minuit->SetFixedVariable(3,"index",(Double_t)k);
      } 
      else if(coreName[k] == "CINDY"){   
    	minuit->SetVariable(0,"scale", 2.50040e-02, 0.0001);
    	minuit->SetVariable(1,"sigma", 1.98657e-01, 0.0001);
    	minuit->SetVariable(2,"mean", 0, 0.0001);
    	minuit->SetFixedVariable(3,"index",(Double_t)k);
      } 
      else if(coreName[k] == "KRISTA"){   
    	minuit->SetVariable(0,"scale", 2.50040e-02, 0.0001);
    	minuit->SetVariable(1,"sigma", 1.98657e-01, 0.0001);
    	minuit->SetVariable(2,"mean", 0, 0.0001);
    	minuit->SetFixedVariable(3,"index",(Double_t)k);
      } 
      else if(coreName[k] == "BROOKE"){   
    	minuit->SetVariable(0,"scale", 2.50040e-02, 0.0001);
    	minuit->SetVariable(1,"sigma", 1.98657e-01, 0.0001);
    	minuit->SetVariable(2,"mean", 0, 0.0001);
    	minuit->SetFixedVariable(3,"index",(Double_t)k);
      }

      minuit->Minimize(); //has to be in the same layer as ROOT::Math::Functor f()

    }//if backward angles





    //***************execute Minuit**************
    //minuit->Minimize();
    cout<<"haha"<<coreName[k]<<endl;

    const Double_t *result = minuit->X();
    //cout<<GetChi2(result)<<endl;
    hist_data[k] = hdata;
    hist_fit[k] = hfit;
    hconv[k] = new TH1F(coreName[k]+"_conv","convolution response", nbins, xlow, xhigh);
    for(int n=0; n<nbins; n++){
      hconv[k]->SetBinContent(n,DoConvolution(hgeant,result,n));
    }
    if(k<4)
      fexp[k] = exp_func;

  }//loop over k






  //************plot and save histograms************* 
  gStyle->SetOptStat(0);
  TCanvas *c = new TCanvas("c", "c", 1600, 1000);
  c->Divide(4,2);
  for(int k=0; k<kNcores; k++){
    c->cd(k+1);
    hist_data[k]->SetLineColor(kBlack);
    hist_data[k]->SetTitleSize(1.5);
    hist_data[k]->GetXaxis()->SetRangeUser(drawlow,drawhigh);
    hist_data[k]->GetXaxis()->SetLabelSize(0.055);
    hist_data[k]->GetXaxis()->SetTitleSize(0.05);
    hist_data[k]->GetYaxis()->SetLabelSize(0.05);
    hist_data[k]->GetYaxis()->SetNdivisions(510);
    hist_data[k]->Draw("E0");  
    hconv[k]->SetLineColor(kRed); 
    hconv[k]->SetLineWidth(2);
    hconv[k]->Draw("c,same");
    if(k<4){
      hist_fit[k]->SetLineColor(kBlue); 
      hist_fit[k]->SetLineWidth(2);
      hist_fit[k]->Draw("c,same");
      fexp[k]->SetLineColor(kGreen);
      fexp[k]->SetLineWidth(2);
      fexp[k]->Draw("same");
    }

    TLine *zero = new TLine(drawlow,0,drawhigh,0);
    zero->SetLineStyle(3);
    zero->Draw("same");


    cout<<hconv[k]->Integral(hconv[k]->FindBin(50),hconv[k]->FindBin(75))<<endl;

  }


  c->SaveAs("a.pdf");

  // TFile *fout = new TFile(outfile,"recreate");
  // c->Write();
  // for(int k=0; k<kNcores; k++){
  //   hist_data[k]->Write();
  //   hconv[k]->Write();
  //   if(k<4)
  //     fexp[k]->Write();
  //   hist_fit[k]->Write();
  // }






}







void GetData(Int_t index){

  //load data histogram
  TFile *file = new TFile(infile);
  hdata = (TH1F*)file->Get(coreName[index]+"_NetSUB")->Clone(coreName[index]+"_data");
  hdata->Sumw2();
  hdata->Rebin(rebin);

  //load geant histogram(has same binning as data)
  TFile *gfile = new TFile("geant_LD2_65MeV_50mil.root");
  hgeant = (TH1F*)gfile->Get(coreName[index]);
  hgeant->Rebin(rebin);

  //set fitting hist as the same binning as data hist
  nbins =hdata->GetNbinsX();
  xlow = hdata->GetXaxis()->GetBinLowEdge(1);
  xhigh = hdata->GetXaxis()->GetBinUpEdge(nbins);
  hfit = new TH1F(coreName[index]+"_fit","total fit",nbins,xlow,xhigh);
  exp_func = new TF1(coreName[index]+"_exp_func","[0]*exp(-[1]*(x*x-27*27))",10,100);

}







Double_t GetFitVal(TH1F *geant_hist, const Double_t *par, Int_t j){

  Double_t E, val, sigma_gauss, exp, gauss;
  Double_t con_val = 0.0;
  Double_t exp_val, fit_val;
  Double_t E_conv = geant_hist->GetXaxis()->GetBinCenter(j);

  for(int i=1; i<geant_hist->GetNbinsX(); i++){
    E = geant_hist->GetBinCenter(i);
    val = E_conv-E-par[2];
    val = val*val;
    sigma_gauss = par[1]*E_beam/2.354*sqrt(E_beam/E);
    exp = -1.0*(val/(2.0*sigma_gauss*sigma_gauss));
    gauss = TMath::Exp(exp)/sigma_gauss;

    con_val += par[0]*geant_hist->GetBinContent(i)*gauss;
  }

  //  exp_func = new TF1("exp_func","[0]*exp(-[1]*(x*x-27*27))",10,100); 
  exp_func->SetParameter(0, par[3]*1.0e4);
  exp_func->SetParameter(1, par[4]);
  exp_val = exp_func->Eval(E_conv);

  fit_val = con_val + exp_val /*+ par[5]*/;
  //cout << fit_val << endl;
 
  return fit_val;
  
}


Double_t DoConvolution(TH1F *geant_hist, const Double_t *par, Int_t j){

  Double_t E, val, sigma_gauss, exp, gauss;
  Double_t con_val = 0.0;
  Double_t E_conv = geant_hist->GetXaxis()->GetBinCenter(j);


  for(int i=1; i<geant_hist->GetNbinsX(); i++){
    E = geant_hist->GetBinCenter(i);
    val = E_conv-E-par[2];
    val = val*val;
    sigma_gauss = par[1]*E_beam/2.354*sqrt(E_beam/E);
    exp = -1.0*(val/(2.0*sigma_gauss*sigma_gauss));
    gauss = TMath::Exp(exp)/sigma_gauss;

    con_val += par[0]*geant_hist->GetBinContent(i)*gauss;
  }
 
  return con_val;
  
}




Double_t GetChi2_fwd(const Double_t *par){
  for(int j=0; j<nbins; j++){
    hfit->SetBinContent(j,GetFitVal(hgeant,par,j));
  }
  Double_t sigma_chi, chi;
  Double_t chi2 = 0.0;
  Int_t fit_low = hfit->FindBin(Efit_low[(Int_t)par[5]]);
  Int_t fit_high = hfit->FindBin(Efit_high[(Int_t)par[5]]);
  for(int i=fit_low; i<fit_high; i++){
    //if (i > 40 && i < 55) continue;
    chi = hfit->GetBinContent(i)-hdata->GetBinContent(i);
    //sigma_chi = hdata->GetBinError(i);
    sigma_chi = sqrt(hdata->GetBinContent(i));
    chi = chi/sigma_chi;
    chi2 += chi*chi;
  }

  return chi2;
}


Double_t GetChi2_bwd(const Double_t *par){
  for(int j=0; j<nbins; j++){
    hfit->SetBinContent(j,DoConvolution(hgeant,par,j));
  }
  Double_t sigma_chi, chi;
  Double_t chi2 = 0.0;
  Int_t fit_low = hfit->FindBin(Efit_low[(Int_t)par[3]]);
  Int_t fit_high = hfit->FindBin(Efit_high[(Int_t)par[3]]);
  for(int i=fit_low; i<fit_high; i++){
    //if (i > 40 && i < 55) continue;
    chi = hfit->GetBinContent(i)-hdata->GetBinContent(i);
    //sigma_chi = hdata->GetBinError(i);
    sigma_chi = sqrt(hdata->GetBinContent(i));
    chi = chi/sigma_chi;
    chi2 += chi*chi;
  }

  return chi2;
}





