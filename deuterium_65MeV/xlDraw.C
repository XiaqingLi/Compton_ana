{
  gStyle->SetOptStat(0);
  TFile *file1  = new TFile("./old_analysis/full05_sumw2.root");
  TFile *file2 = new TFile("xlFull.root");
  TFile *file3 = new TFile("xlEmpty.root");
  TFile *file4 = new TFile("xlNet.root");

  TH1F *hraw, *htofshld, *hnet, *hempty, *hlineshape;

  hraw = (TH1F*)file1->Get("LINDA_coreRaw");
  htofshld = (TH1F*)file2->Get("LINDA_TOFShldCut");
  hnet = (TH1F*)file2->Get("LINDA_Net");
  hempty = (TH1F*)file3->Get("LINDA_Net");
  hlineshape = (TH1F*)file4->Get("LINDA_NetSUB");

  hraw->Rebin(4);
  htofshld->Rebin(4);
  hnet->Rebin(4);
  hempty->Rebin(4);
  hlineshape->Rebin(4);

  TCanvas *c = new TCanvas();
  hraw->SetLineColor(kAzure+1);
  htofshld->SetLineColor(kMagenta);
  hnet->SetLineColor(kBlack);
  hempty->SetLineColor(kGreen+1);
  hlineshape->SetLineColor(kRed+1);


  //################ draw raw and shield+tof cut ############
  // c.SetLogy();
  // htofshld->SetTitle("Energy spectrum in HINDA core detector");
  // htofshld->GetXaxis()->SetTitle("Energy [MeV]");
  // htofshld->GetXaxis()->CenterTitle();
  // htofshld->GetXaxis()->SetTitleSize(0.045);
  // htofshld->GetXaxis()->SetLabelSize(0.045); 
  // htofshld->GetYaxis()->SetTitle("counts / 0.4 MeV");
  // htofshld->GetYaxis()->CenterTitle();
  // htofshld->GetYaxis()->SetTitleSize(0.045);
  // htofshld->GetYaxis()->SetLabelSize(0.045);
  // htofshld->GetXaxis()->SetRangeUser(35,80);
  // //htofshld->GetYaxis()->SetRangeUser(0,500000);

  // htofshld->Draw();
  // hraw->Draw("same");

  // TLegend *leg = new TLegend(0.3,0.6,0.7,0.8);
  // leg->AddEntry(hraw,"No cut","lep");
  // leg->AddEntry(htofshld,"Shield+TOF cut","lep");
  // leg->Draw();

  //##############  draw random subtraction effect ###########
  // htofshld->SetTitle("Energy spectrum in HINDA core detector");
  // htofshld->GetXaxis()->SetTitle("Energy [MeV]");
  // htofshld->GetXaxis()->CenterTitle();
  // htofshld->GetXaxis()->SetTitleSize(0.045);
  // htofshld->GetXaxis()->SetLabelSize(0.045); 
  // htofshld->GetYaxis()->SetTitle("counts / 0.4 MeV");
  // htofshld->GetYaxis()->CenterTitle();
  // htofshld->GetYaxis()->SetTitleSize(0.045);
  // htofshld->GetYaxis()->SetLabelSize(0.045);
  // htofshld->GetXaxis()->SetRangeUser(35,80);

  // htofshld->Draw();
  // hnet->Draw("same");

  // TLegend *leg = new TLegend(0.3,0.6,0.7,0.8);
  // leg->AddEntry(htofshld,"Before rand sub","lep");
  // leg->AddEntry(hnet,"After rand sub","lep");

  // leg->Draw();

  //############## draw full and empty comparison #########
  hnet->SetTitle("Energy spectrum in HINDA core detector");
  hnet->GetXaxis()->SetTitle("Energy [MeV]");
  hnet->GetXaxis()->CenterTitle();
  hnet->GetXaxis()->SetTitleSize(0.045);
  hnet->GetXaxis()->SetLabelSize(0.045); 
  hnet->GetYaxis()->SetTitle("counts / 0.4 MeV");
  hnet->GetYaxis()->CenterTitle();
  hnet->GetYaxis()->SetTitleSize(0.045);
  hnet->GetYaxis()->SetLabelSize(0.045);
  hnet->GetXaxis()->SetRangeUser(30,80);

  hnet->Draw();
  hempty->Draw("same");

  TLegend *leg = new TLegend(0.3,0.6,0.7,0.8);
  leg->AddEntry(hnet,"Full","lep");
  leg->AddEntry(hempty,"Empty","lep");

  leg->Draw();



  //##############  draw empty subtraction effect ###########
  // hnet->SetTitle("Energy spectrum in HINDA core detector");
  // hnet->GetXaxis()->SetTitle("Energy [MeV]");
  // hnet->GetXaxis()->CenterTitle();
  // hnet->GetXaxis()->SetTitleSize(0.045);
  // hnet->GetXaxis()->SetLabelSize(0.045); 
  // hnet->GetYaxis()->SetTitle("counts / 0.4 MeV");
  // hnet->GetYaxis()->CenterTitle();
  // hnet->GetYaxis()->SetTitleSize(0.045);
  // hnet->GetYaxis()->SetLabelSize(0.045);
  // hnet->GetXaxis()->SetRangeUser(30,80);

  // hnet->Draw();
  // hlineshape->Draw("same");

  // TLegend *leg = new TLegend(0.3,0.6,0.7,0.8);
  // leg->AddEntry(hnet,"Full","lep");
  // leg->AddEntry(hlineshape,"Subtracted","lep");

  // leg->Draw();


}
