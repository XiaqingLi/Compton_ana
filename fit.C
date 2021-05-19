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
#include "TGraph.h"
#include "TGraphErrors.h"

using namespace std;

const Int_t nHINDA = 8;
const Double_t E_beam = 84;
Double_t Efit_low;
Double_t Efit_high;
const Double_t Lshift = 0;
const Double_t Hshift = 0;
const Int_t kRebin = 4;
Double_t binWidth = 0.2*kRebin;
TString outfile = "Fit.root";
TString infile = "/home/xl79/mepg/helium_84MeV/Net.root";
const TString coreNames[nHINDA] = {"ALAINA","BROOKE","CINDY","JONI","KRISTA","LINDA","ROBERTA","SUSAN"};
const Double_t angles[nHINDA] = {55.0, 90.0, 125.0, 125.0, 125.0, 90.0, 90.0, 55.0};

const Double_t kEfit_low[nHINDA] = {34, 26, 69, 69, 69, 26, 26, 35};
const Double_t kEfit_high[nHINDA] = {109, 110, 109, 109, 109, 109, 109, 109};

Int_t rebin = 1;
const Double_t Emax = 145.;

TH1F *hgeant;
TH1F *hdata; 
TH1F *hconv;
TH1F *hfit;
TH1F *hbg[nHINDA], *hbg_sub[nHINDA];
TF1 *exp_func = new TF1("exp_func","[0]*exp(-[1]*(x*x-27*27))+[2]",Emax/10,Emax); 
TF1 *exp_bg[nHINDA];
Int_t bin_num;
Int_t l;
Double_t xlow,xhigh;
Double_t xE;
Double_t chiPerDof[8];
Double_t sigma_r[8], scale_r[8], mean_r[8], exp_scale_r[8], exp_decay_r[8], const_r[8],sigma_err_r[8], scale_err_r[8], mean_err_r[8], exp_scale_err_r[8], exp_decay_err_r[8], const_err_r[8];

TFile *file = new TFile(infile);
TFile *gfile = new TFile("/home/xl79/mepg/git_Compton/4He_Oct2017/HindaSim_revisedMode7_D16mm_EM_QGSP_BERT_HP_5e7.root");
TFile *fout = new TFile(outfile, "recreate");


TCanvas *c1 = new TCanvas("c1","c1",1800,1200);
TCanvas *c2 = new TCanvas("c2","c2",1800,1200);
TCanvas *c3 = new TCanvas("c3","c3",1800,1200);
Double_t E, val, sigma_gauss, expo, gauss, con_val, E_conv, fit_val, exp_val;
Double_t sigma_chi, chi, chi2, dof, chi_per_dof;
Int_t fit_low, fit_high;

void GetData(int index);
Double_t GetFitVal(TH1F *geant_hist, const Double_t *par, Int_t j);
Double_t DoConvolution(TH1F *geant_hist, const Double_t *par, Int_t j);
Double_t GetChi2(const Double_t *par);
Double_t GetChiPerDof(const Double_t *par);



int main(){
  // gStyle->SetOptStat(0);
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadRightMargin(0.03);
  c1->Divide(4,2);
  c2->Divide(4,2);
  c3->Divide(4,2);

  for(l=0; l<nHINDA; l++){
    Efit_low = kEfit_low[l]+Lshift*binWidth;
    Efit_high = kEfit_high[l]+Hshift*binWidth;

    GetData(l);

    ROOT::Math::Minimizer *minuit = ROOT::Math::Factory::CreateMinimizer("Minuit","");
    minuit->SetMaxFunctionCalls(10000000); 
    minuit->SetTolerance(0.0000001);
    minuit->SetPrintLevel(0);
    ROOT::Math::Functor f(&GetChi2,6);//Fitting function
    minuit->SetFunction(f);

    if(l>=2 && l<=4){
      minuit->SetVariable(0,"scale", 0.01, 0.0001);
      // minuit->SetVariableLimits(0,0,1);
      minuit->SetVariable(1,"sigma", 0.25, 0.0001);
      // minuit->SetVariableLimits(1,0.01,1);
      minuit->SetVariable(2,"mean", 0.95, 0.001);
      minuit->SetFixedVariable(3,"exp_scale", 0);  
      minuit->SetFixedVariable(4,"exp_decay", 0);
      minuit->SetVariable(5,"const", 20.0, 0.0001);
    }
    else{
      minuit->SetVariable(0,"scale", 0.01, 0.0001);
      // minuit->SetVariableLimits(0,0,1);
      minuit->SetVariable(1,"sigma", 0.03, 0.0001);
      // minuit->SetVariableLimits(1,0.01,1);
      minuit->SetVariable(2,"mean", 0.95, 0.001);
      minuit->SetVariable(3,"exp_scale", 0.05, 0.0001);  
      minuit->SetVariable(4,"exp_decay", 0.01, 0.0001);
      minuit->SetVariable(5,"const", 50.0, 0.0001);
    }


    minuit->Minimize();
    cout<<coreNames[l]<<"\t";
    minuit->PrintResults();


    //*********** get fit result ************
    const Double_t *result = minuit->X();
    const Double_t *error = minuit->Errors();
    scale_r[l] = result[0];
    sigma_r[l] = result[1];
    mean_r[l] = result[2];
    exp_scale_r[l] = result[3];
    exp_decay_r[l] = result[4];
    const_r[l] = result[5];
    scale_err_r[l] = error[0];
    sigma_err_r[l] = error[1];
    mean_err_r[l] = error[2];
    exp_scale_err_r[l] = error[3];
    exp_decay_err_r[l] = error[4];
    const_err_r[l] = error[5];
    chiPerDof[l] = GetChiPerDof(result);
    // cout<<GetChi2(result)<<endl;
    // cout<<result[0]<<"\t"<<result[1]<<"\t\t"<<result[2]<<"\t"<<GetChiPerDof(result)<<endl;

    for(int j=0; j<bin_num; j++){
      hconv->SetBinContent(j,DoConvolution(hgeant,result,j));
      hfit->SetBinContent(j,GetFitVal(hgeant,result,j));
    }

    exp_bg[l] = new TF1("exp_bg_"+coreNames[l],"[0]*exp(-[1]*(x*x-27*27))+[2]",10,Emax); 
    exp_bg[l]->SetParameter(0, result[3]*1.0e4);
    exp_bg[l]->SetParameter(1, result[4]);
    exp_bg[l]->SetParameter(2, result[5]);

    for(int j=0; j<bin_num; j++){
      xE = hbg[l]->GetBinCenter(j);
      hbg[l]->SetBinContent(j,exp_bg[l]->Eval(xE));

      double sigma1 = error[3]*1.0e4;
      double sigma2 = error[4];
      double sigma3 = error[5];
      double dirv1 = TMath::Exp(-result[4]*(xE*xE-27*27));
      double dirv2 = result[3]*1.0e4*TMath::Exp(-result[4]*(xE*xE-27*27))*(27*27-xE*xE);

      double hbgError = sqrt(dirv1*dirv1*sigma1*sigma1+dirv2*dirv2*sigma2*sigma2+sigma3*sigma3);
 
      hbg[l]->SetBinError(j,hbgError);
    }
    hbg_sub[l] ->Add(hbg[l], -1);


    c1->cd(l+1);
    hdata->SetLineColor(kBlack);
    hdata->GetXaxis()->SetRangeUser(25, 110);
    hdata->GetXaxis()->SetLabelSize(0.045);
    hdata->GetXaxis()->SetTitle("E [MeV]");
    hdata->GetXaxis()->SetTitleSize(0.045);
    // hdata->GetXaxis()->CenterTitle();
    // hdata->GetYaxis()->SetRangeUser(0, 60*binWidth/0.4);
    hdata->GetYaxis()->SetTitle(Form("Counts / %1.1f MeV",binWidth));
    hdata->GetYaxis()->SetTitleSize(0.045);
    hdata->GetYaxis()->SetTitleOffset(1.5);
    hdata->GetYaxis()->SetLabelSize(0.045);
    hdata->GetYaxis()->SetLabelOffset(0.003);
    // hdata->GetYaxis()->CenterTitle();
    // hdata->SetTitle("");
    hdata->Draw("E0");  
    hconv->SetLineColor(kBlue);
    hconv->SetLineWidth(2);
    hconv->Draw("c, same");
    hfit->SetLineColor(kRed);
    hfit->SetLineWidth(2);
    hfit->Draw("c, same");
    exp_bg[l]->SetLineColor(kGreen);
    exp_bg[l]->SetLineWidth(2);
    exp_bg[l]->Draw("same");
    hbg_sub[l]->SetLineColor(kMagenta);
    hbg_sub[l]->Draw("same");


    TLine *low = new TLine(Efit_low,0,Efit_low,500);
    low->SetLineWidth(1);
    low->SetLineStyle(2);
    // low->SetLineColor(kBlue);
    low->Draw("same");
    TLine *high = new TLine(Efit_high,0,Efit_high,500);
    high->SetLineWidth(1);
    high->SetLineStyle(2);
    // high->SetLineColor(kBlue);
    high->Draw("same");

    // TLatex *t = new TLatex();
    // t->SetNDC();
    // //  t->SetTextColor(1);
    // t->SetTextSize(0.065);
    // //  t.SetTextAlign(13);
    // //  t.SetTextAngle(0);
    // t->DrawLatex(.55,.75,"#theta = 55#circ");


    c2->cd(l+1);
    hdata->Draw("E0");  
    hconv->Draw("same");
    exp_bg[l]->Draw("same");
    hfit->Draw("same");
    low->Draw("same");
    high->Draw("same");

    hconv->Write();
    hfit->Write();
    exp_bg[l]->SetTitle("exp_"+coreNames[l]);
    exp_bg[l]->Write();
    hdata->Write();
    hbg[l]->Write(); 
    hbg_sub[l]->Write();
  }

  Double_t gr_x[8] = {1,2,3,4,5,6,7,8};
  TGraph *gr = new TGraph(8, gr_x, chiPerDof);
  gr->SetName("ChiPerDOF");
  gr->Write();
  TGraphErrors *gre_scale = new TGraphErrors(8, gr_x, scale_r, 0, scale_err_r);
  gre_scale->SetName("FitScale");
  gre_scale->Write();
  TGraphErrors *gre_sigma = new TGraphErrors(8, gr_x, sigma_r, 0, sigma_err_r);
  gre_sigma->SetName("FitSigma");
  gre_sigma->Write();
  TGraphErrors *gre_mean = new TGraphErrors(8, gr_x, mean_r, 0, mean_err_r);
  gre_mean->SetName("FitMean");
  gre_mean->Write();
  TGraphErrors *gre_exp_scale = new TGraphErrors(8, gr_x, exp_scale_r, 0, exp_scale_err_r);
  gre_exp_scale->SetName("FitExpScale");
  gre_exp_scale->Write();
  TGraphErrors *gre_exp_decay = new TGraphErrors(8, gr_x, exp_decay_r, 0, exp_decay_err_r);
  gre_exp_decay->SetName("FitDecay");
  gre_exp_decay->Write();
  TGraphErrors *gre_const = new TGraphErrors(8, gr_x, const_r, 0, const_err_r);
  gre_const->SetName("FitConst");
  gre_const->Write();


  c2->SaveAs("FitAll.pdf");
  c1->Write();
  c2->Write();
  fout->Close();

  return 0;
}






void GetData(int index){
  hdata = (TH1F*)file->Get(coreNames[index] + "_Net")->Clone("hist_data_"+coreNames[index]);
  hdata->Rebin(kRebin);

  bin_num = hdata->GetNbinsX();
  xlow = hdata->GetXaxis()->GetBinLowEdge(1);
  xhigh = hdata->GetXaxis()->GetBinUpEdge(bin_num);
  hfit = new TH1F("hist_fit_"+coreNames[index],"total fit", bin_num, xlow, xhigh);
  hconv = new TH1F("hist_conv_"+coreNames[index],"convolution response", bin_num, xlow, xhigh);
  hbg[index] = new TH1F("hist_bg_"+coreNames[index],"background",bin_num, xlow, xhigh);
  hbg_sub[index] = (TH1F*)hdata->Clone("hist_bg_sub_"+coreNames[index]);

  hgeant = (TH1F*)gfile->Get(coreNames[index]);
  rebin = (Int_t)(hdata->GetBinWidth(10)/hgeant->GetBinWidth(20));
  if(rebin>1)  
    hgeant->Rebin(rebin);
 
  if(hgeant->GetBinWidth(10) != hdata->GetBinWidth(10)){
    cout<<"\nERROR: BIN WIDTH DOES NOT MATCH!!!"<<endl;
    exit(1);
  }

}

Double_t GetFitVal(TH1F *geant_hist, const Double_t *par, Int_t j){
  fit_val = 0.0;
  E_conv = par[2]*geant_hist->GetXaxis()->GetBinCenter(j);

  for(int i=1; i<geant_hist->GetNbinsX(); i++){
    E = geant_hist->GetBinCenter(i);
    val = E_conv-E;
    sigma_gauss = par[1]*E/2.354*sqrt(E_beam/E);
    expo = -1.0*(val*val/(2.0*sigma_gauss*sigma_gauss));
    gauss = TMath::Exp(expo)/(sigma_gauss*sqrt(2*TMath::Pi()));

    fit_val += par[0]*geant_hist->GetXaxis()->GetBinWidth(i)*geant_hist->GetBinContent(i)*gauss;
  }
  exp_func->SetParameter(0, par[3]*1.0e4);
  exp_func->SetParameter(1, par[4]);
  exp_func->SetParameter(2, par[5]);
  exp_val = exp_func->Eval(E_conv);
  fit_val = fit_val + exp_val;

  return fit_val;
}



Double_t DoConvolution(TH1F *geant_hist, const Double_t *par, Int_t j){
  con_val = 0.0;
  E_conv = par[2]*geant_hist->GetXaxis()->GetBinCenter(j);

  for(int i=1; i<geant_hist->GetNbinsX(); i++){
    E = geant_hist->GetBinCenter(i);
    val = E_conv-E;
    sigma_gauss = par[1]*E/2.354*sqrt(E_beam/E);
    expo = -1.0*(val*val/(2.0*sigma_gauss*sigma_gauss));
    gauss = TMath::Exp(expo)/(sigma_gauss*sqrt(2*TMath::Pi()));

    con_val += par[0]*geant_hist->GetXaxis()->GetBinWidth(i)*geant_hist->GetBinContent(i)*gauss;
  }
 
  return con_val;
}


Double_t GetChi2(const Double_t *par){
  for(int j=0; j<bin_num; j++){
    hfit->SetBinContent(j,GetFitVal(hgeant,par,j));
  }

  chi2 = 0.0;
  fit_low = hfit->FindBin(Efit_low);
  fit_high = hfit->FindBin(Efit_high);
  for(int i=fit_low; i<fit_high; i++){
    if(hdata->GetBinContent(i)>1){
      //if (i > 40 && i < 55) continue;
      chi = hfit->GetBinContent(i)-hdata->GetBinContent(i);
      sigma_chi = hdata->GetBinError(i);
      chi = chi/sigma_chi;
      chi2 += chi*chi;
    }
  }
  return chi2;
}


Double_t GetChiPerDof(const Double_t *par){
  chi2 = 0.0;
  dof = 0.0;
  fit_low = hfit->FindBin(Efit_low);
  fit_high = hfit->FindBin(Efit_high);
  for(int i=fit_low; i<fit_high; i++){
    if(hdata->GetBinContent(i)>1){
      //if (i > 40 && i < 55) continue;
      chi = hfit->GetBinContent(i)-hdata->GetBinContent(i);
      sigma_chi = hdata->GetBinError(i);
      chi = chi/sigma_chi;
      chi2 += chi*chi;
      dof += 1;
    }
  }
  if(par[3]==0 && par[4]==0)
    chi_per_dof = chi2/(dof-4);
  else
    chi_per_dof = chi2/(dof-6);

  // cout<<dof<<endl;

  return chi_per_dof;
}
