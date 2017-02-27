using namespace std;

const Double_t kNgamma = 2.32e13;
const Double_t kTarThick = 9.75e23; //nuclei/cm^2
const Double_t kAttenuation = 0.9789; //photon attenuation factor
// const Int_t kNCores = 4;
// const TString kCoreName[kNCores] = {"LINDA","CINDY","KRISTA","BROOKE"};
// const Double_t kCoreAngle[kNCores] = {110.0, 125.0, 145.0, 159.0};
// const Double_t kEfit_high[kNCores] = {72.0, 73.0, 72.0, 71.0};
// const Double_t kEfit_low = 50.0;
// const Double_t solidAng = 0.043; //[sr]
const Int_t kNCores = 8;
const TString kCoreName[kNCores] = {"ALAINA", "SUSAN", "JONI", "ROBERTA",
				   "LINDA", "CINDY", "KRISTA", "BROOKE"};
const Double_t kCoreAngle[kNCores] = {40.0, 55.0, 75.0, 90.0, 
				      110.0, 125.0, 145.0, 159.0};
// const Double_t kEfit_high[kNCores] = {72.0, 72.0, 74.0, 74.0,
// 				      72.0, 73.0, 72.0, 71.0};
// const Double_t kEfit_low[kNCores] = {57.0, 57.0, 57.0, 53.0,
// 				     50.0, 50.0, 50.0, 50.0};
const Double_t kEfit_low[kNCores] = {53.0, 53.0, 53.0, 53.0,
				     55.0, 55.0, 55.0, 55.0};
const Double_t kEfit_high[kNCores] = {72.0, 72.0, 72.0, 70.0,
				      70.0, 70.0, 70.0, 70.0};
const Double_t solidAng[kNCores] = {0.0438, 0.0415, 0.0426, 0.0396,
				    0.0423, 0.0430, 0.0422, 0.0455};
const Double_t Ebg_low = 80;
const Double_t Ebg_high = 100;

//const Double_t kEfit_low = 57.0;

TH1F *hfit;
TH1F *hdata;
Double_t yield[kNCores];
Double_t yield_err[kNCores];
Double_t dXS[kNCores];
Double_t dXS_err[kNCores];
Double_t yield_sub[kNCores];
Double_t yield_err_sub[kNCores];
Double_t dXS_sub[kNCores];
Double_t dXS_err_sub[kNCores];
Double_t ratio1,ratio2;

//****change input and output file names here!****
TString infile = "xlFitCore.root";
//TString infile = "FitCoreConst.root";
//TString outfile = "dXS09_combine_rebin12.root";
TString outfile = "xlDXS_BremSub.root";
//************************************************



void xlCalDXS(){

  TFile *file = new TFile(infile);

  //calculate yield
  for(int i=0; i<kNCores; i++){
    //    hfit = (TH1F*)file->Get(kCoreName[i]+"_fit");
    hfit = (TH1F*)file->Get(kCoreName[i]+"_conv");
    hdata = (TH1F*)file->Get(kCoreName[i]+"_data");
    Int_t binfit_low = hfit->FindBin(kEfit_low[i]);
    Int_t binfit_high = hfit->FindBin(kEfit_high[i]);
    Int_t bin_low = hfit->FindBin(10);
    Int_t bin_high = hfit->FindBin(75);
    Int_t Nfit_w = hfit->Integral(binfit_low, binfit_high);
    Int_t Nfit_t = hfit->Integral(bin_low, bin_high);
    //Int_t Ndata_w = hdata->Integral(binfit_low, binfit_high);
    Int_t Ndata_w = hdata->IntegralAndError(binfit_low, binfit_high,yield_err[i]);
 
    yield[i] = (Double_t)Nfit_t*Ndata_w/Nfit_w;
    //yield[i] = (Double_t)Nfit_t;
    yield_err[i] = (Double_t)Nfit_t*yield_err[i]/Nfit_w;
 
    dXS[i] = 1e33*yield[i]/(kNgamma*solidAng[i]*kTarThick*kAttenuation);//[nb/sr]
    dXS_err[i] = dXS[i]*yield_err[i]/yield[i];
 
   //for high energy background crection 
    Int_t bg_low = hdata->FindBin(Ebg_low);
    Int_t bg_high = hdata->FindBin(Ebg_high);
    Int_t Nbackground = hdata->Integral(bg_low, bg_high);
    Double_t factor = (75.0-10.0)/(Ebg_high-Ebg_low);
    yield_sub[i] = yield[i]-Nbackground*factor;
    yield_err_sub[i] = sqrt(yield_err[i]*yield_err[i]+factor*factor*Nbackground);

    dXS_sub[i] = 1e33*yield_sub[i]/(kNgamma*solidAng[i]*kTarThick*kAttenuation);//[nb/sr]
    dXS_err_sub[i] = dXS_sub[i]*yield_err_sub[i]/yield_sub[i];


    ratio1 = (Double_t)1.0-dXS_sub[i]/dXS[i];
    ratio2 = (Double_t)hdata->Integral(bg_low,bg_high)*factor/hfit->Integral(bin_low, bin_high);
    cout<<ratio1<<", "<<ratio2<<", "<<ratio1/ratio2<<endl;

  } 


  // TCanvas *c1 = new TCanvas("cYield","Yield");
  // //TGraph *gr1 = new TGraph(kNCores, kCoreAngle, yield);
  // TGraphErrors *gr1 = new TGraphErrors(kNCores, kCoreAngle, yield, 0, yield_err);
  // gr1->SetTitle("Yield");
  // gr1->GetXaxis()->SetTitle("#theta [deg]");
  // gr1->GetXaxis()->CenterTitle();
  // gr1->GetXaxis()->SetTitleSize(0.045);
  // gr1->GetXaxis()->SetLabelSize(0.045);
  // gr1->SetMarkerStyle(kFullCircle);
  // gr1->SetMarkerSize(0.7);
  // gr1->SetMarkerColor(kBlue);
  // gr1->Draw("ap");

  TCanvas *c2 = new TCanvas("cdXS","Differential cross section");
  c2->SetGrid();
  //TGraph *gr2 = new TGraph(kNCores, kCoreAngle, dXS);
  TGraphErrors *gr2 = new TGraphErrors(kNCores, kCoreAngle, dXS, 0, dXS_err);
  gr2->SetTitle("Differential cross section, LD2 @ 65 MeV");
  gr2->GetXaxis()->SetTitle("#theta [deg]");
  gr2->GetXaxis()->SetTitleSize(0.045);
  gr2->GetXaxis()->CenterTitle();
  gr2->GetXaxis()->SetLabelSize(0.045);
  gr2->GetYaxis()->SetTitle("d#sigma/d#Omega [nb/sr]");
  gr2->GetYaxis()->SetTitleSize(0.045);
  gr2->GetYaxis()->CenterTitle();
  gr2->GetYaxis()->SetLabelSize(0.045);
  gr2->GetYaxis()->SetRangeUser(0, 35);
  gr2->SetMarkerStyle(kFullCircle);
  gr2->SetMarkerSize(0.7);
  //gr2->SetMarkerColor(kRed);
  //  gr2->Draw("ap");
  //  TGraph *theory = new TGraph("Harald.dat","%lg %lg");
  // theory->SetLineWidth(2);
  // theory->SetLineColor(kRed);
  // theory->Draw("csame");

  TGraphErrors *gr3 = new TGraphErrors(kNCores, kCoreAngle, dXS_sub, 0, dXS_err_sub);
  gr3->SetTitle("");
  gr3->GetXaxis()->SetTitle("#theta_{lab} [deg]");
  gr3->GetXaxis()->SetTitleSize(0.045);
  gr3->GetXaxis()->CenterTitle();
  gr3->GetXaxis()->SetLabelSize(0.045);
  gr3->GetYaxis()->SetTitle("d#sigma/d#Omega [nb/sr]");
  gr3->GetYaxis()->SetTitleSize(0.045);
  gr3->GetYaxis()->CenterTitle();
  gr3->GetYaxis()->SetLabelSize(0.045);
  gr3->GetYaxis()->SetRangeUser(0, 26);
  gr3->SetMarkerStyle(kFullCircle);
  gr3->SetMarkerSize(0.7);
  // gr3->SetMarkerColor(kGreen+1);
  // gr3->SetLineColor(kGreen+1);
  gr3->Draw("apE");

  TFile *fmark = new TFile("lH2_xs_65mev.root");
  TGraph *theory = (TGraph*)fmark->Get("gr_central");
  theory->Draw("same");

  TLegend *leg = new TLegend(0.57, 0.2, 0.83, 0.35);
  //  leg->AddEntry(gr2,"HI#gammaS data","ep");
  leg->AddEntry(gr3,"HI#gammaS, 2016","ep");
  leg->AddEntry(theory, "#chiEFT prediction","l");
  leg->Draw();

  
  //*************for comparison **************
  // TCanvas *c3 = new TCanvas("cComp","Differential cross section comparison");
  // c3->SetGrid();
  // TGraphErrors *gMark = (TGraphErrors*)fmark->Get("gre_xs");
  // gMark->SetMarkerColor(kAzure+1);
  // gMark->SetLineColor(kAzure+1);
  // gMark->SetMarkerStyle(kFullCircle);
  // gMark->SetMarkerSize(0.7);

  // TLegend *leg2 = new TLegend(0.57, 0.2, 0.83, 0.35);
  // leg2->AddEntry(gr2,"HI#gammaS data (Xiaqing)","ep");
  // leg2->AddEntry(gMark,"HI#gammaS data (Mark)","ep");
  // leg2->AddEntry(theory, "Griesshammer","l");

  // gr2->Draw("ap");
  // gMark->Draw("p,same");
  // theory->Draw("same");
  // leg2->Draw();

  


  // TFile *fout = new TFile(outfile,"recreate"); 
  // //  gr1->Write("gYield");
  // gr2->Write("gDXS");
  // gr3->Write("gDXS_Sub");
  // //  c1->Write();
  // c2->Write();
  // //  c3->Write();



}
