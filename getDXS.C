//helium-4 target, 84 MeV circ.
void getDXS(){

  const Double_t kNgamma = 3.12047e+12;
  const Double_t kNgamma_Cindy = 2.47755e+12;
  const Double_t kTarThick = 4.17e23; //nuclei/cm^2
  const Double_t kAttenuation = 0.989084; //photon attenuation factor

  const Int_t kNCores = 8;
  const TString kCoreName[kNCores] = {"ALAINA", "BROOKE", "CINDY", "JONI", 
				      "KRISTA", "LINDA", "ROBERTA", "SUSAN"};
  const Double_t kCoreAngle[kNCores] = {55.0, 90.0, 125.0, 125.0, 
					125.0, 90.0, 90.0, 55.0};
  const Double_t coreCuts[8] = {7.5, 9.5, 11, 8, 17.5, 11, 15.3, 10.5};
  const Double_t E_high[8] = {98, 95, 99, 95, 98, 92, 95, 93}; //for 84 MeV analysis
  const Double_t E_low[8] = {74, 74, 72, 72, 72, 72, 72, 74}; //for 84 MeV analysis


  TH1F *hconv[kNCores], *hbg_sub[kNCores], *htotal[kNCores];
  TF1 *exp_bg[kNCores];
  Double_t fit_N0[kNCores], fit_N1[kNCores], data_N1[kNCores], data_err[kNCores], yield[kNCores], yield_err[kNCores];
  Int_t sum_low[8], sum_high[8], bin_low, bin_high; 
  Double_t N_sim = 5e7;
  Double_t effSolAng[kNCores];
  Double_t dxs[kNCores], dxs_err[kNCores], bg_const[kNCores];

  TFile *flineshape = new TFile("Fit.root");
  TGraph *gr_scale = (TGraph*)flineshape->Get("FitScale");
  TGraph *gr_mean = (TGraph*)flineshape->Get("FitMean");


  for(int i=0; i<kNCores; i++){

    hconv[i] = (TH1F*)flineshape->Get("hist_conv_"+kCoreName[i]);
    hbg_sub[i] = (TH1F*)flineshape->Get("hist_bg_sub_"+kCoreName[i]);
    htotal[i] = (TH1F*)flineshape->Get("hist_data_"+kCoreName[i]);
    exp_bg[i] = (TF1*)flineshape->Get("exp_bg_"+kCoreName[i]);
    bin_low = hconv[i]->FindBin(coreCuts[i]);
    bin_high = hconv[i]->FindBin(100);
    sum_high[i] = hconv[i]->FindBin(E_high[i]);
    sum_low[i] = hconv[i]->FindBin(E_low[i]);


    fit_N0[i] = hconv[i]->Integral(bin_low, bin_high);
    fit_N1[i] = hconv[i]->Integral(sum_low[i], sum_high[i]);
    data_N1[i] = hbg_sub[i]->IntegralAndError(sum_low[i], sum_high[i], data_err[i]);

    yield[i] = data_N1[i]*fit_N0[i]/fit_N1[i];
    yield_err[i] = data_err[i]*fit_N0[i]/fit_N1[i];
    effSolAng[i] = 4*TMath::Pi()*hconv[i]->Integral(bin_low, bin_high)*gr_mean->Eval(i+1.0)/(gr_scale->Eval(i+1.0)*N_sim);

//    cout<<effSolAng[i]*1000<<endl;
    //////Equivalently:
    // yield[i] = hbg_sub[i]->IntegralAndError(sum_low[i], sum_high[i], yield_err[i]);
    // effSolAng[i] = 4*TMath::Pi()*hconv[i]->Integral(sum_low[i], sum_high[i])/(gr_scale->Eval(i+1.0)*N_sim);

    dxs[i] = 1e33*yield[i]/(kNgamma*kTarThick*kAttenuation*effSolAng[i]); 
    if(i==2)
      dxs[i] = 1e33*yield[i]/(kNgamma_Cindy*kTarThick*kAttenuation*effSolAng[i]); 
    dxs_err[i] = yield_err[i]*dxs[i]/yield[i];
  }

  for(int i=0; i<8; i++)
    cout<<dxs[i]<<"\n";
  cout<<endl;
  for(int i=0; i<8; i++)
    cout<<dxs_err[i]<<"\n";
  cout<<endl;

  TGraphErrors *gr = new TGraphErrors(kNCores, kCoreAngle, dxs, 0, dxs_err);

	
  //////////////// separate phi angles //////////////
  Double_t angle[3] = {55.0, 90.0, 125.0};
  Double_t angle2[2] = {90,125};
  Double_t dxs_R[2] = {dxs[5], dxs[4]};
  Double_t dxs_L[3] = {dxs[7], dxs[6], dxs[2]};
  Double_t dxs_D[3] = {dxs[0], dxs[1], dxs[3]};
  Double_t dxs_stat_R[2] = {dxs_err[5], dxs_err[4]};
  Double_t dxs_stat_L[3] = {dxs_err[7], dxs_err[6], dxs_err[2]};
  Double_t dxs_stat_D[3] = {dxs_err[0], dxs_err[1], dxs_err[3]};
  Double_t dxs_syst_R[2] = {1.0, 1.2}; //point-to-point syst
  Double_t dxs_syst_L[3] = {1.1, 1.0, 1.2}; //point-to-point syst
  Double_t dxs_syst_D[3] = {2.2, 0.7, 1.1}; //point-to-point syst
  Double_t dxs_err_R[2];
  Double_t dxs_err_L[3];
  Double_t dxs_err_D[3];
  for(int i=0; i<3; i++){
    if(i<2)
      dxs_err_R[i] = sqrt(dxs_stat_R[i]*dxs_stat_R[i]+dxs_syst_R[i]*dxs_syst_R[i]);
    dxs_err_D[i] = sqrt(dxs_stat_D[i]*dxs_stat_D[i]+dxs_syst_D[i]*dxs_syst_D[i]);
    dxs_err_L[i] = sqrt(dxs_stat_L[i]*dxs_stat_L[i]+dxs_syst_L[i]*dxs_syst_L[i]);
  }
  TGraphErrors *gr_R = new TGraphErrors(2, angle2, dxs_R, 0, dxs_err_R);
  TGraphErrors *gr_L = new TGraphErrors(3, angle, dxs_L, 0, dxs_err_L);
  TGraphErrors *gr_D = new TGraphErrors(3, angle, dxs_D, 0, dxs_err_D);


  ///////////////// calculate weighted average ////////////
  Double_t weighted_dxs[3],weighted_err[3];
  weighted_dxs[0] = calculate_dxs(dxs[0], dxs_err[0], dxs[7], dxs_err[7], 0, 0);
  weighted_dxs[1] = calculate_dxs(dxs[1], dxs_err[1], dxs[5], dxs_err[5], dxs[6], dxs_err[6]);
  weighted_dxs[2] = calculate_dxs(dxs[2], dxs_err[2], dxs[3], dxs_err[3], dxs[4], dxs_err[4]);
  weighted_err[0] = calculate_err(dxs_err[0], dxs_err[7], 0);
  weighted_err[1] = calculate_err(dxs_err[1], dxs_err[5], dxs_err[6]);
  weighted_err[2] = calculate_err(dxs_err[2], dxs_err[3], dxs_err[4]);

  TGraphErrors *gr_weighted = new TGraphErrors(3, angle, weighted_dxs, 0, weighted_err);


  for(int i=0; i<3; i++){
    cout<<weighted_dxs[i]<<endl;
  }
  cout<<endl;
  for(int i=0; i<3; i++){
    cout<<weighted_err[i]<<endl;
  }
  cout<<endl;


  cout<<"\n55deg\t"<<weighted_dxs[0]<<"\t"<< weighted_err[0]<<"\t"<<weighted_err[0]/weighted_dxs[0]<<endl;
  cout<<"90deg\t"<<weighted_dxs[1]<<"\t"<< weighted_err[1]<<"\t"<<weighted_err[1]/weighted_dxs[1]<<endl;
  cout<<"125deg\t"<<weighted_dxs[2]<<"\t"<< weighted_err[2]<<"\t"<<weighted_err[2]/weighted_dxs[2]<<endl;



  ////////////////////////  HIGS 61 MeV result///////////////
  Double_t Higs_angle[7] = {40, 55, 75, 110, 125, 145, 159};
  Double_t Higs_dxs[7] = {119.8, 115.9, 89.8, 85.8, 95.3, 120.0, 150.2};
  Double_t Higs_stat[7] = {6.8, 5.4, 3.2, 3.3, 2.9, 3.2, 3.5};
  Double_t Higs_sys[7] = {4.2, 1.5, 1.5, 1.3, 1.5, 1.1, 1.5};
  Double_t Higs_error[7];
  for(int i=0; i<7; i++){
    Higs_error[i] = Higs_stat[i]*Higs_stat[i] + Higs_sys[i]*Higs_sys[i]+Higs_dxs[i]*0.022*Higs_dxs[i]*0.022;
    Higs_error[i] = sqrt(Higs_error[i]);
  }
  TGraphErrors *gr_Higs = new TGraphErrors(7, Higs_angle, Higs_dxs, 0, Higs_error);


  ////////////////////////   Lund 87 MeV result//////////////////
  Double_t Lund_angle[3] = {60.0, 90.0, 150.0};
  Double_t Lund_dxs[3] = {64.0, 51.0, 132.0};
  Double_t Lund_stat[3] = {10.0, 7.0, 15.0};
  Double_t Lund_error[3];
  for(i=0; i<3; i++)
    Lund_error[i] = sqrt(Lund_stat[i]*Lund_stat[i] + Lund_dxs[i]*0.15*Lund_dxs[i]*0.15);
  TGraphErrors *gr_Lund = new TGraphErrors(3, Lund_angle, Lund_dxs, 0, Lund_error);


	
  ///////////////////////// plot graph ////////////////////////
  TCanvas *c1 = new TCanvas("cDXS","cDXS",900,700);
  gPad->SetBottomMargin(0.13);
  // gPad->SetLeftMargin(0.13);
  // gPad->SetRightMargin(0.05);
  //  c1->SetGrid();
  gr->SetTitle("");
  gr->GetXaxis()->SetTitle("#theta_{lab} (deg)");
  gr->GetXaxis()->SetLimits(15, 165);
  gr->GetYaxis()->SetRangeUser(20, 140);
  gr->GetXaxis()->SetTitleSize(0.05);
  gr->GetXaxis()->SetTitleOffset(1);
  gr->GetXaxis()->CenterTitle();
  gr->GetXaxis()->SetLabelSize(0.05);
  gr->GetYaxis()->SetTitle("d#sigma/d#Omega (nb/sr)");
  gr->GetYaxis()->SetTitleSize(0.05);
  gr->GetYaxis()->SetTitleOffset(1);
  gr->GetYaxis()->CenterTitle();
  gr->GetYaxis()->SetLabelSize(0.05);
  gr->SetMarkerStyle(kFullCircle);
  gr->SetMarkerColor(kWhite);
  gr->SetLineColor(kWhite);
  // gr->SetMarkerSize(1.5);
  gr->Draw("ap");

  gr_Higs->SetMarkerStyle(26);
  // gr_Higs->SetMarkerColor(kRed+1);
  // gr_Higs->SetLineColor(kRed+1);
  // gr_Higs->SetMarkerSize(1.5);
  // gr_Higs->Draw("psame");

  gr_Lund->SetMarkerStyle(25);
  gr_Lund->SetMarkerSize(1);
  gr_Lund->SetMarkerColor(kGreen+2);
  gr_Lund->SetLineColor(kGreen+2);
  // gr_Lund->Draw("psame");

  gr_R->SetMarkerSize(2);
  gr_L->SetMarkerSize(2);
  gr_D->SetMarkerSize(2);
  gr_R->SetLineWidth(2);
  gr_L->SetLineWidth(2);
  gr_D->SetLineWidth(2);
  gr_R->SetMarkerColor(kGreen+1);
  gr_L->SetMarkerColor(kBlue);
  gr_D->SetMarkerColor(kRed+1);
  gr_R->SetLineColor(kGreen+1);
  gr_L->SetLineColor(kBlue);
  gr_D->SetLineColor(kRed+1);
  gr_R->SetMarkerStyle(22);
  gr_L->SetMarkerStyle(23);
  gr_D->SetMarkerStyle(20);
  gr_L->Draw("psame");
  gr_D->Draw("psame");
  gr_R->Draw("psame");

  // gr_weighted->SetMarkerStyle(28);
  // gr_weighted->SetMarkerColor(kGreen);
  // gr_weighted->SetLineColor(kGreen);
  gr_weighted->SetMarkerStyle(kFullCircle);
  gr_weighted->SetMarkerSize(1.3);
  gr_weighted->SetLineWidth(1);
  gr_weighted->SetMarkerColor(kRed+1);
  gr_weighted->SetLineColor(kRed+1);
  // gr_weighted->Draw("psame");

  
  //////////////////// plot systematic errors //////////////////
  TGraphErrors *sys_out = new TGraphErrors();
  TGraphErrors *sys_in_right = new TGraphErrors();
  TGraphErrors *sys_in_left = new TGraphErrors();
  TGraphErrors *sys_avr = new TGraphErrors();
  TGraphErrors *sys_onPoint = new TGraphErrors();
  sys_out->SetPoint(1,55,31);
  sys_out->SetPoint(2,90,31);
  sys_out->SetPoint(3,125,31);
  sys_out->SetPointError(1,0,1.40);
  sys_out->SetPointError(2,0,1.30);
  sys_out->SetPointError(3,0,1.74);
  sys_out->SetLineWidth(1.7);
  sys_out->SetLineColor(kRed+1);
  sys_out->SetMarkerColor(kRed+1);
  // sys_out->Draw("pesame");
  sys_in_left->SetPoint(1,55,26);
  sys_in_left->SetPoint(2,90,26);
  sys_in_left->SetPoint(3,125,26);
  sys_in_left->SetPointError(1,0,1.65);
  sys_in_left->SetPointError(2,0,1.12);
  sys_in_left->SetPointError(3,0,1.15);
  sys_in_left->SetLineWidth(1.7);
  sys_in_left->SetLineColor(kBlue);
  sys_in_left->SetMarkerColor(kBlue);
  // sys_in_left->Draw("pesame");
  sys_in_right->SetPoint(1,90,21);
  sys_in_right->SetPoint(2,125,21);
  sys_in_right->SetPointError(1,0,2.06);
  sys_in_right->SetPointError(2,0,1.17);
  sys_in_right->SetLineWidth(1.7);
  sys_in_right->SetLineColor(kGreen+2);
  sys_in_right->SetMarkerColor(kGreen+2);
  // sys_in_right->Draw("pesame");
  sys_avr->SetPoint(1,55,21);
  sys_avr->SetPoint(2,90,21);
  sys_avr->SetPoint(3,125,21);
  sys_avr->SetPointError(1,0,1.87);
  sys_avr->SetPointError(2,0,1.10);
  sys_avr->SetPointError(3,0,1.11);
  sys_avr->SetLineWidth(1.7);
  sys_avr->SetLineColor(kRed+1);
  sys_avr->SetMarkerColor(kRed+1);
  // sys_avr->Draw("pesame");
  sys_onPoint->SetPoint(1,55,weighted_dxs[0]);
  sys_onPoint->SetPoint(2,90,weighted_dxs[1]);
  sys_onPoint->SetPoint(3,125,weighted_dxs[2]);
  sys_onPoint->SetPointError(1,0,sqrt(weighted_err[0]*weighted_err[0]+1.87*1.87+weighted_dxs[0]*0.022*weighted_dxs[0]*0.022));
  sys_onPoint->SetPointError(2,0,sqrt(weighted_err[1]*weighted_err[1]+1.10*1.10+weighted_dxs[1]*0.022*weighted_dxs[1]*0.022));
  sys_onPoint->SetPointError(3,0,sqrt(weighted_err[2]*weighted_err[2]+1.11*1.11+weighted_dxs[2]*0.022*weighted_dxs[2]*0.022));
  sys_onPoint->SetMarkerColor(kRed+1);
  //sys_onPoint->SetMarkerSize(1.5);
  sys_onPoint->SetLineColor(kRed+1);
  // sys_onPoint->Draw("pesame");


  ///////////////// plot legend ///////////////////////
  TLegend *leg = new TLegend(0.3, 0.7, 0.7, 0.9);
  leg->AddEntry(gr_L,"#phi = 0#circ (in plane)","p");
  leg->AddEntry(gr_R,"#phi = 180#circ (in plane)","p");
  leg->AddEntry(gr_D,"#phi = 270#circ (out of plane)","p");
  // leg->AddEntry(gr_R,"#parallel (Right)","p");
  // leg->AddEntry(gr_D,"#perp","p");
  // leg->AddEntry(gr_L,"#parallel (Left)","p");
  // leg->AddEntry(gr_weighted,"This work, E_{#gamma}=84 MeV","p");
  // leg->AddEntry(gr_Higs,"HI#gammaS, E_{#gamma}=61 MeV","p");
  // leg->AddEntry(gr_Lund, "Lund, E_{#gamma}=87 MeV  ","p");
  leg->SetFillColor(0);
  leg->Draw();

  c1->SaveAs("dxs_new.pdf");

  //////////////////   draw fittings and summing windows  /////////////////
  TCanvas *c0 = new TCanvas("cFit","cFit");
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadRightMargin(0.03);
  TLine *low_sum, *high_sum;
  c0->Divide(4,2);
  TH1F *hfit[8];
  for(int i=0; i<8; i++){
    hfit[i] = (TH1F*)flineshape->Get("hist_fit_"+kCoreName[i]);

    c0->cd(i+1);
    htotal[i]->GetXaxis()->SetRangeUser(50,110);
    htotal[i]->Draw();
    hconv[i]->Draw("csame");
    exp_bg[i]->Draw("same");
    hfit[i]->Draw("csame");
   
    low_sum = new TLine(E_low[i],0,E_low[i],170);
    high_sum = new TLine(E_high[i],0,E_high[i],170);
    low_sum->SetLineColor(kMagenta);
    high_sum->SetLineColor(kMagenta);
    low_sum->SetLineStyle(2);
    high_sum->SetLineStyle(2);
    low_sum->Draw();
    high_sum->Draw();
  }


  // TFile *fout = new TFile("all.root","update");
  //c1->SaveAs("dxs.pdf");
}



Double_t calculate_dxs(Double_t var1, Double_t err1, Double_t var2, Double_t err2, Double_t var3, Double_t err3){
  if(var3 == 0){
    Double_t weightSum = 1/(err1*err1) + 1/(err2*err2);
    Double_t avg = var1/(err1*err1)+var2/(err2*err2);
  } 
  else{
    Double_t weightSum = 1/(err1*err1) + 1/(err2*err2) + 1/(err3*err3);
    Double_t avg = var1/(err1*err1)+var2/(err2*err2)+var3/(err3*err3);
  }
  avg = avg/weightSum;

  return avg;
}

Double_t calculate_err(Double_t err1, Double_t err2, Double_t err3){
  if(err3 == 0){
    Double_t weighted_error = 1/(err1*err1) + 1/(err2*err2);
  } 
  else{
    Double_t weighted_error = 1/(err1*err1) + 1/(err2*err2) + 1/(err3*err3);
  }
  weighted_error = sqrt(1/weighted_error);

  return weighted_error;
}
