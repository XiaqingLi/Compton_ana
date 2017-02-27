// Mohammad Ahmed
// -------------------------------------------------------------------------
// Need these headers if you want to run compiled code, otherwise ignore ---
#include "TApplication.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TF1.h"
#include "TF2.h"
#include "TH2F.h"
#include "TFormula.h"
#include "TFile.h"
#include "TList.h"
#include "TMath.h"
#include "TROOT.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TTree.h"
#include "TStyle.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "Math/Polynomial.h"
#include "Math/Interpolator.h"
#include <stdint.h>
#include <sys/stat.h>
using namespace std;


// -----------------------------------------------------------------------


//int flashQDC(int runno, bool draw_traces = false) {
int flashQDC(int numchannels) {
 
  bool draw_traces=true;

  //  load_run_database();
  // ROOT Global Variables ------------------------------
  //gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  //gStyle->SetTitleFillColor(10);
  //gStyle->SetTitleBorderSize(0);
  // ---------------------------------------------------
 
  //********* digitizer configuration ***********
  const int kSampleRate = 500; // in MHz 
  const double kTimeStep = 1000.0 / (double)kSampleRate; // in ns
  const int kNumBits = 14;
  const int kRecordBuffer = 1030; 
  const double kPostTrigger = 0.8; // 70%
 
  double ZTH = 1.0;
  // set baseline low/high values
  int BaseLow = 0;
  int BaseHigh = 100;
  cout << "BaseLow = " << BaseLow << "   BaseHigh = " << BaseHigh << endl; 
  // set integration gate width for qdc
  int StartSignal = 130;
  int EndSignal = 800;
  //for beam test*********
  BaseHigh=100;
  StartSignal=120;
  EndSignal=1000;
  //for beam test, Feb 9,2016
  cout << "StartSignal = " << StartSignal << "   EndSignal = " << EndSignal << endl;
  //********************************************
  short dataword;
  FILE * infile[numchannels];
  int samples;
  int eventnum = 0;
  long lSize;
  int voltage[kRecordBuffer];
  int BaselineSub[kRecordBuffer];
  int signal[kRecordBuffer];
  Double_t xi[kRecordBuffer];
  Double_t yi[kRecordBuffer];
  //****************************************
  float Q[numchannels];
  float Qinter;
  float pedestal[numchannels];
  float Baseline,BaselineSamples;
  float BaselineEntry[numchannels];
  int BaseLowEntry[numchannels];
  int BaseHighEntry[numchannels];
  int StartSignalEntry[numchannels];
  int EndSignalEntry[numchannels];
  Int_t timestamp[numchannels];
  float Xmax[numchannels];
  //*****************************************
  
  
  // ROOT OBJECTs --------------------------------------
  TFile* Rootfile = new TFile("traces.root", "RECREATE");
  TTree* tree = new TTree("flat_tree", "Digitizer Tree");
  TTree** temp_trees = new TTree*[numchannels];//fill each det. separately then merge later
  for (int i = 0; i < numchannels; ++i) temp_trees[i] = new TTree(Form("temp_tree_ch%i", i), "temp tree");
  //----Loop over dets -----------------------------
  int thedet = 0;
  TString filename;
  while (thedet < numchannels) {
    if (thedet < 16) {
      filename = Form("data_in%i", thedet);
    } else if(thedet > 15 && thedet < 32) {
      filename = Form("data_in_b%i", thedet - 16);
    } else {
      filename = "";
    }
    infile[thedet] = fopen(filename.Data(), "rb");
    if (infile[thedet] == NULL) {
      ++thedet;
      break;
    }
    cout << "Found file " << filename << endl;
    
    temp_trees[thedet]->Branch(Form("Q_%i", thedet), &Q[thedet], Form("Q_%i/F", thedet));
    temp_trees[thedet]->Branch(Form("pedestal_%i", thedet), &pedestal[thedet], Form("pedestal_%i/F", thedet));
    temp_trees[thedet]->Branch(Form("Baseline_%i", thedet), &BaselineEntry[thedet], Form("Baseline_%i/F", thedet));// store baseline average
    temp_trees[thedet]->Branch(Form("BaseStart_%i", thedet), &BaseLowEntry[thedet], Form("BaseStart_%i/I", thedet));// store baseline start record
    temp_trees[thedet]->Branch(Form("BaseEnd_%i", thedet), &BaseHighEntry[thedet], Form("BaseEnd_%i/I", thedet));// store baseline end record
    temp_trees[thedet]->Branch(Form("SignalStart_%i", thedet), &StartSignalEntry[thedet], Form("SignalStart_%i/I", thedet));// store signal start
    temp_trees[thedet]->Branch(Form("SignalEnd_%i", thedet), &EndSignalEntry[thedet], Form("SignalEnd_%i/I", thedet));// store signal end
    temp_trees[thedet]->Branch(Form("TimeStamp_%i", thedet), &timestamp[thedet], Form("TimeStamp_%i/I", thedet));// store timestamp
    temp_trees[thedet]->Branch(Form("Xmax_%i", thedet), &Xmax[thedet], Form("Xmax_%i/F", thedet));// store max bin position
    tree->Branch(Form("Q_%i", thedet), &Q[thedet], Form("Q_%i/F", thedet));
    tree->Branch(Form("pedestal_%i", thedet), &pedestal[thedet], Form("pedestal_%i/F", thedet));
    tree->Branch(Form("Baseline_%i", thedet), &BaselineEntry[thedet], Form("Baseline_%i/F", thedet));// store baseline average
    tree->Branch(Form("BaseStart_%i", thedet), &BaseLowEntry[thedet], Form("BaseStart_%i/I", thedet));// store baseline start record
    tree->Branch(Form("BaseEnd_%i", thedet), &BaseHighEntry[thedet], Form("BaseEnd_%i/I", thedet));// store baseline end record
    tree->Branch(Form("SignalStart_%i", thedet), &StartSignalEntry[thedet], Form("SignalStart_%i/I", thedet));// store signal start
    tree->Branch(Form("SignalEnd_%i", thedet), &EndSignalEntry[thedet], Form("SignalEnd_%i/I", thedet));// store signal end
    tree->Branch(Form("TimeStamp_%i", thedet), &timestamp[thedet], Form("TimeStamp_%i/I", thedet));// store timestamp
    tree->Branch(Form("Xmax_%i", thedet), &Xmax[thedet], Form("Xmax_%i/F", thedet));// store max bin position
    ++thedet;
  }


  //make histograms to draw traces
  const int ntraces = 500;
  TH1F*** trace = new TH1F**[numchannels];
  for (int i = 0; i < numchannels; ++i) trace[i] = new TH1F*[ntraces];
  for (int x = 0; x < numchannels; ++x) {
    for (int y = 0; y < ntraces; ++y) {
      trace[x][y] = new TH1F(Form("ch%i_trace_%i", x, y),";Record [#]",
			     kRecordBuffer, 0, kRecordBuffer);
      trace[x][y]->SetDirectory(0);
    }
  }

  // ***********
  // struct stat st0,st16;
  // stat("data_in0", &st0);
  // stat("data_in16", &st16);
  // cout << st0.st_size << "   " << st16.st_size << endl;
  // if (st0.st_size == st16.st_size) {
  //   cout << st0.st_size << "   " << st16.st_size << endl;
  //   cout << "@@@@@digitizers synched!" << endl;
  // }
  fseek (infile[0] , 0 , SEEK_END);
  lSize = ftell (infile[0]);
  rewind (infile[0]);
  int samples0 = (lSize/2.0)/((float) kRecordBuffer+12);
  fseek (infile[16] , 0 , SEEK_END);
  lSize = ftell (infile[16]);
  rewind (infile[16]);
  int samples16 = (lSize/2.0)/((float) kRecordBuffer+12);
  cout << samples0 << "   " << samples16 << endl;
  if (samples0 != samples16) cout << "@@@@@@@@@@@@@@@@@@ digitizers not synched!" << endl;
  // ***************


  uint32_t header[6];
  // loop over all input files ---------------------------------------------------
  for(int f = 0; f < thedet; ++f) {
    fseek (infile[f] , 0 , SEEK_END);
    lSize = ftell (infile[f]);
    rewind (infile[f]);
    cout << "File size is " << lSize / 1024 / 1024 << " MB" << endl;
    samples = (lSize/2.0)/((float) kRecordBuffer+12);
    cout << "Total Samples are : " << samples << endl;
    eventnum = 0;
    
    // set baseline and signal ranges
    if (f == 16 || f == 19 || f == 22 || f == 23) {// paddles 
      BaseLow = 0;
      BaseHigh = 70;
      StartSignal = 80;
      EndSignal = 1000;
    } else if (f == 14 || f == 30) {// TACs
      BaseLow = 0;
      BaseHigh = 300;
      StartSignal = 400;
      EndSignal = 1000;
    } else {// core + shields
      BaseLow = 0;
      BaseHigh = 80;
      StartSignal = 90;
      EndSignal = 1000;
    }

    if (f == 17 || f == 18 || f == 20 || f == 21 || (f > 23 && f < 30)) continue;
    
    for (int n = 0; n < samples; ++n) {
      if (n % 50000 == 0) cout << "Reading and Processing Event " << n<< endl;
      Baseline = 0.0;
      BaselineSamples = 0.0; 
      Q[f] = 0.0;
      Qinter = 0.0;
      pedestal[f] = 0.0;
	
      BaseLowEntry[f] = BaseLow;
      BaseHighEntry[f] = BaseHigh;
      StartSignalEntry[f] = StartSignal;
      EndSignalEntry[f] = EndSignal;

      fread(&header[0],4,6,infile[f]);
      timestamp[f] = header[5];
      //for (int x=0;x<6;++x)cout<<header[x]<<endl;
      // The Event is processed here ---------------------------------
      for (int i = 0; i < kRecordBuffer; ++i) {
	fread(&dataword,2,1,infile[f]);
	voltage[i] = 0;
	voltage[i] = dataword;
	xi[i] = (double) i;
	yi[i] = (double) voltage[i];
	
	// establish baseline
	if (i > BaseLow && i < BaseHigh) {
	  BaselineSamples = BaselineSamples + 1.0;
	  Baseline = Baseline + voltage[i];
	} 
      } // over i, all records

      Baseline = Baseline / BaselineSamples; // Average Baseline
      BaselineEntry[f] = Baseline;
      
      //fill traces
      for (int i = 0; i < kRecordBuffer; ++i) {
      	// baseline subtraction
	if (f == 14 || f == 30) {// timing
	  BaselineSub[i] = voltage[i] - BaselineEntry[f];
	  //BaselineSub[i] = voltage[i];
	} else {
	  BaselineSub[i] = Baseline - voltage[i];
	}

      	signal[i] = 0.0;
      	if (BaselineSub[i] > ZTH) {
	  signal[i] = BaselineSub[i];
	}

      	if (eventnum < ntraces) {
	  trace[f][eventnum]->Fill(i, signal[i]);  // filling few samples
	}

	// QDC integration
	if (i > StartSignal && i < EndSignal && f!=14 && f!=30) {//timing
	  Q[f] = Q[f] + kTimeStep * signal[i];
	  pedestal[f] = pedestal[f] + pow(2.0,(double)kNumBits);
	} 
	
	
      } // over i, all records
      ++eventnum;

      Xmax[f] = 0;          
      if (f == 14 || f == 30) {//timing
	Q[f] = signal[TMath::LocMax(kRecordBuffer, signal)]; 
      } else {
	Q[f] = 2.0 * Q[f] / ((EndSignal - StartSignal));
	if (f == 0 || f == 3 || f == 6 || f == 7 || f == 16 || f == 19 || f == 22 || f == 23 )Xmax[f] = TMath::LocMax(kRecordBuffer, signal);
      }
     
      temp_trees[f]->Fill();

    } // over n, samples in file
    cout << "finished channel " << f << endl;
  } // over f, all channels
  
  for (int f = 0; f < thedet; ++f) {
    fclose(infile[f]);
    temp_trees[f]->SetBranchAddress(Form("Q_%i", f), &Q[f]);
    temp_trees[f]->SetBranchAddress(Form("pedestal_%i", f), &pedestal[f]);
    temp_trees[f]->SetBranchAddress(Form("Baseline_%i", f), &BaselineEntry[f]);// store baseline average
    temp_trees[f]->SetBranchAddress(Form("BaseStart_%i", f), &BaseLowEntry[f]);// store baseline start record
    temp_trees[f]->SetBranchAddress(Form("BaseEnd_%i", f), &BaseHighEntry[f]);// store baseline end record
    temp_trees[f]->SetBranchAddress(Form("SignalStart_%i", f), &StartSignalEntry[f]);// store signal start
    temp_trees[f]->SetBranchAddress(Form("SignalEnd_%i", f), &EndSignalEntry[f]);// store signal end
    temp_trees[f]->SetBranchAddress(Form("TimeStamp_%i", f), &timestamp[f]);// store timestamp
    temp_trees[f]->SetBranchAddress(Form("Xmax_%i", f), &Xmax[f]);// store timestamp
  }

  // merge trees into one hist
  int n_temp_events = temp_trees[0]->GetEntries();
  bool IsEqual = true;
  for (int i = 1; i < numchannels; ++i) {
    if (i == 17 || i == 18 || i == 20 || i == 21 || (i > 23 && i < 30)) continue;
    if (temp_trees[i]->GetEntries() != n_temp_events) {
      IsEqual = false;
      break;
    }
  }

  if (IsEqual) {
    for (int i = 0; i < n_temp_events; ++i) {
      for (int f = 0; f < numchannels; ++f) {
	if (f == 17 || f == 18 || f == 20 || f == 21 || (f > 23 && f < 30)) continue;
	temp_trees[f]->GetEntry(i);
      }
      tree->Fill();
    }
  } else {
    cout << "error: trees are not of equal length!" << endl;
    exit(1);
  }

  //--------------------------------------------------------------------------
  Rootfile->cd();
  tree->Write();
  if (draw_traces) {
    cout << "drawing traces" << endl;

    TCanvas* tC[numchannels];
    for (int f = 0; f < thedet; ++f) {
      tC[f] = new TCanvas(Form("ch%i_tC", f), Form("Channel %i Traces", f), 600, 600);
      //tC[f]->SetFrameFillColor(6);
      //tC[f]->cd();
      trace[f][1]->SetLineStyle(1);
      trace[f][1]->SetLineWidth(2);
      trace[f][1]->SetLineColor(1);
      trace[f][1]->Draw();
      trace[f][1]->GetYaxis()->SetRangeUser(0, pow(2, kNumBits));
      for (int j = 2; j < ntraces; ++j) {
	
	trace[f][j]->SetLineStyle(1);
	trace[f][j]->SetLineWidth(2);
	trace[f][j]->SetLineColor(1);
	trace[f][j]->Draw("same");
	trace[f][j]->Write(Form("ch%i_trace_%i",f,j));
      } // over j, number of waveforms to save
      tC[f]->Write();
    } // over f, number of detectors
    
  } // draw_traces
  Rootfile->Close();
  
  cout << "saved output to " << Rootfile->GetName() << endl;
  //gApplication->Terminate();
  return 0;
} // end of function
