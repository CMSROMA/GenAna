void drawHisto(string output, std::vector< TH1F* > thisVector, std::vector<float> Vec_M, std::vector<float> Vec_R, std::vector<string> Vec_label,int rebin, string xaxisTitle, string yaxisTitle, float xmin, float xmax, bool logX, bool norm){

  TCanvas c1;

  std::vector< TPaveStats* > paveVec;

  if(logX && thisVector.size()>0)
    {
      c1.SetLogx();
      thisVector[0]->GetXaxis()->SetRangeUser(1,xmax);
    }

  //TLegend *legend=new TLegend(0.55,0.74,0.95,0.98);
  TLegend *legend=new TLegend(0.15,0.74,0.55,0.98);
  legend->SetTextFont(72);
  legend->SetTextSize(0.035);

  float ymax = 0;

  for(int k=0;k<thisVector.size();k++)
    {
      cout << k << endl;      

      float Nevents = thisVector[k]->GetEntries(); 

      char name [100];
      sprintf(name, "M=%.1f #lambda=%.1f %s", Vec_M[k] , Vec_R[k], Vec_label[k].c_str() ); 
      thisVector[k]->SetTitle(0);
      thisVector[k]->Rebin(rebin);
      thisVector[k]->SetLineColor(k+1);
      thisVector[k]->GetXaxis()->SetTitle(xaxisTitle.c_str());
      thisVector[k]->GetYaxis()->SetTitle(yaxisTitle.c_str());
      thisVector[k]->GetXaxis()->SetRangeUser(xmin,xmax);

      if(norm == true)
	thisVector[k]->Scale(float(1/Nevents));

      if(thisVector[k]->GetMaximum()>ymax)
	ymax=thisVector[k]->GetMaximum();

      legend->AddEntry(thisVector[k],name,"l");

      if(k==0)
	{
	  thisVector[k]->Draw("HIST");
	}
      else
	{
	  thisVector[k]->Draw("HISTsames");
	}

      c1.Update();

      TPaveStats *stats = (TPaveStats*)c1.GetPrimitive("stats");
      //TPaveStats *stats =(TPaveStats*)(thisVector[k]->GetListOfFunctions()->FindObject("stats"));
      stats->SetName(Vec_label[k].c_str());
      stats->SetY1NDC(.8-k*0.2);
      stats->SetY2NDC(1.0-k*0.2);
      stats->SetTextColor(k+1);      
      paveVec.push_back(stats);
    }

  for(int k=0;k<paveVec.size();k++)
    {
      paveVec[k]->Print();
      paveVec[k]->Draw();
    }

  thisVector[0]->GetYaxis()->SetRangeUser(0,ymax*1.3);
  thisVector[0]->GetXaxis()->SetTicks("+-");
  thisVector[0]->GetXaxis()->SetTitleOffset(1.5);

  legend->Draw();

  char outputname[100];
  sprintf(outputname,"%s.png",output.c_str());
  c1.SaveAs(outputname);

  c1.Close();

}

/*
void drawEff(string output, std::vector< TH1F* > DenVector, std::vector< TH1F* > NumVector, std::vector<float> Vec_M, std::vector<float> Vec_R, int rebin, string xaxisTitle, string yaxisTitle,string selection){

  TCanvas c1;

  TH2F *h2_efficiency = new TH2F("h2_efficiency",selection.c_str(),50,0,10000,21,-0.025,1.025);

  for(int k=0;k<DenVector.size();k++)
    {
      float den = DenVector[k]->GetEntries();
      float num = NumVector[k]->GetEntries();

      float eff = 0;
      if(den>0)
	eff = num/den;

      cout << "M,R = " << Vec_M[k] << "," << Vec_R[k] << " : Efficiency="<< eff << endl;

      h2_efficiency->Fill(Vec_M[k],Vec_R[k],eff);
    }

  h2_efficiency->GetZaxis()->SetRangeUser(0,1);
  h2_efficiency->GetZaxis()->SetTitleOffset(0.7);
  h2_efficiency->GetZaxis()->SetTitle("Selection efficiency");
  h2_efficiency->GetXaxis()->SetTitle(xaxisTitle.c_str());
  h2_efficiency->GetYaxis()->SetTitle(yaxisTitle.c_str());
  h2_efficiency->SetMarkerSize(1.5);

  h2_efficiency->Draw("colztext");

  char outputname[100];
  sprintf(outputname,"%s.png",output.c_str());
  c1.SaveAs(outputname);

  c1.Close();
}
*/

void makePlots_singleLQ(){

  // general options

  gROOT->Reset();
  gStyle->SetOptStat(1111111);

  //load histograms
  TFile *fileHwOld = TFile::Open("Output_HerwigOld/LeptonInducedLQ_umu_18_11_22_M3000_Lambda1p0_mod__1_ANALYSIS.root");
  TFile *fileHwNew = TFile::Open("Output_HerwigNew/LeptonInducedLQ_umu_18_11_22_M3000_Lambda1p0_mod__1_ANALYSIS.root");
  TFile *fileHwNewLO = TFile::Open("Output_HerwigNewLO/LeptonInducedLQ_mod_umu_M3000_Lambda1p0_numEvent1000_ANALYSIS.root");
  TFile *fileHwNewLuca = TFile::Open("Output_HerwigNew-Luca-NLO/LeptonInducedLQ_umu_18_11_22_M3000_Lambda1p0-Luca-NLO_mod_numEvent1000_ANALYSIS.root");

  string myoutputdir = "/afs/cern.ch/user/s/santanas/www/SingleLQ/2022_11_19_SingleLQ_CompareHerwig"; 

  /*
  std::vector< TFile* > Vec_Files;
  Vec_Files.push_back(fileHwOld);
  Vec_Files.push_back(fileHwNew);

  std::vector< float > Vec_Mvalues;
  std::vector< float > Vec_Rvalues;  
  Vec_Mvalues = {2000,3000};
  Vec_Lvalues = {1.0,1.0}; 
  */  

  //make plots
  std::vector< TH1F* > histVec_Mall_Lambda1p0_LQ_mass;
  string variable = "genAnalyzer/h1_LQ_mass";
  histVec_Mall_Lambda1p0_LQ_mass.push_back( (TH1F*)fileHwOld->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_LQ_mass.push_back( (TH1F*)fileHwNew->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_LQ_mass.push_back( (TH1F*)fileHwNewLO->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_LQ_mass.push_back( (TH1F*)fileHwNewLuca->Get(variable.c_str())->Clone() );
  std::vector< float > Vec_M = {3000, 3000, 3000, 3000};
  std::vector< float > Vec_L = {1.0, 1.0, 1.0, 1.0};
  std::vector< string > Vec_label = {"Hw Old NLO", "Hw New NLO", "Hw New LO", "Hw New NLO (Luca)"};
  drawHisto(myoutputdir+"/histVec_Mall_Lambda1p0_LQ_mass",histVec_Mall_Lambda1p0_LQ_mass,Vec_M,Vec_L,Vec_label,2,"LQ mass [GeV]","Events",0,5000,false,false);
  Vec_M.clear(); Vec_L.clear();

  std::vector< TH1F* > histVec_Mall_Lambda1p0_lep_pt;
  variable = "genAnalyzer/h1_lep_pt";
  histVec_Mall_Lambda1p0_lep_pt.push_back( (TH1F*)fileHwOld->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_lep_pt.push_back( (TH1F*)fileHwNew->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_lep_pt.push_back( (TH1F*)fileHwNewLO->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_lep_pt.push_back( (TH1F*)fileHwNewLuca->Get(variable.c_str())->Clone() );
  Vec_M = {3000, 3000, 3000, 3000};
  Vec_L = {1.0, 1.0, 1.0, 1.0};
  Vec_label = {"Hw Old NLO", "Hw New NLO", "Hw New LO", "Hw New NLO (Luca)"};
  drawHisto(myoutputdir+"/histVec_Mall_Lambda1p0_lep_pt",histVec_Mall_Lambda1p0_lep_pt,Vec_M,Vec_L,Vec_label,2,"Lepton pt [GeV]","Events",0,3000,false,false);
  Vec_M.clear(); Vec_L.clear();

  std::vector< TH1F* > histVec_Mall_Lambda1p0_q_pt;
  variable = "genAnalyzer/h1_q_pt";
  histVec_Mall_Lambda1p0_q_pt.push_back( (TH1F*)fileHwOld->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_q_pt.push_back( (TH1F*)fileHwNew->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_q_pt.push_back( (TH1F*)fileHwNewLO->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_q_pt.push_back( (TH1F*)fileHwNewLuca->Get(variable.c_str())->Clone() );
  Vec_M = {3000, 3000, 3000, 3000};
  Vec_L = {1.0, 1.0, 1.0, 1.0};
  Vec_label = {"Hw Old NLO", "Hw New NLO", "Hw New LO", "Hw New NLO (Luca)"};
  drawHisto(myoutputdir+"/histVec_Mall_Lambda1p0_q_pt",histVec_Mall_Lambda1p0_q_pt,Vec_M,Vec_L,Vec_label,2,"Quark pt [GeV]","Events",0,3000,false,false);
  Vec_M.clear(); Vec_L.clear();

  std::vector< TH1F* > histVec_Mall_Lambda1p0_lep_eta;
  variable = "genAnalyzer/h1_lep_eta";
  histVec_Mall_Lambda1p0_lep_eta.push_back( (TH1F*)fileHwOld->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_lep_eta.push_back( (TH1F*)fileHwNew->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_lep_eta.push_back( (TH1F*)fileHwNewLO->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_lep_eta.push_back( (TH1F*)fileHwNewLuca->Get(variable.c_str())->Clone() );
  Vec_M = {3000, 3000, 3000, 3000};
  Vec_L = {1.0, 1.0, 1.0, 1.0};
  Vec_label = {"Hw Old NLO", "Hw New NLO", "Hw New LO", "Hw New NLO (Luca)"};
  drawHisto(myoutputdir+"/histVec_Mall_Lambda1p0_lep_eta",histVec_Mall_Lambda1p0_lep_eta,Vec_M,Vec_L,Vec_label,2,"Lepton #eta","Events",-5,5,false,false);
  Vec_M.clear(); Vec_L.clear();

  std::vector< TH1F* > histVec_Mall_Lambda1p0_q_eta;
  variable = "genAnalyzer/h1_q_eta";
  histVec_Mall_Lambda1p0_q_eta.push_back( (TH1F*)fileHwOld->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_q_eta.push_back( (TH1F*)fileHwNew->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_q_eta.push_back( (TH1F*)fileHwNewLO->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_q_eta.push_back( (TH1F*)fileHwNewLuca->Get(variable.c_str())->Clone() );
  Vec_M = {3000, 3000, 3000, 3000};
  Vec_L = {1.0, 1.0, 1.0, 1.0};
  Vec_label = {"Hw Old NLO", "Hw New NLO", "Hw New LO", "Hw New NLO (Luca)"};
  drawHisto(myoutputdir+"/histVec_Mall_Lambda1p0_q_eta",histVec_Mall_Lambda1p0_q_eta,Vec_M,Vec_L,Vec_label,2,"Quark #eta","Events",-5,5,false,false);
  Vec_M.clear(); Vec_L.clear();

  std::vector< TH1F* > histVec_Mall_Lambda1p0_xip;
  variable = "genAnalyzer/h1_xip";
  histVec_Mall_Lambda1p0_xip.push_back( (TH1F*)fileHwOld->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_xip.push_back( (TH1F*)fileHwNew->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_xip.push_back( (TH1F*)fileHwNewLO->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_xip.push_back( (TH1F*)fileHwNewLuca->Get(variable.c_str())->Clone() );
  Vec_M = {3000, 3000, 3000, 3000};
  Vec_L = {1.0, 1.0, 1.0, 1.0};
  Vec_label = {"Hw Old NLO", "Hw New NLO", "Hw New LO", "Hw New NLO (Luca)"};
  drawHisto(myoutputdir+"/histVec_Mall_Lambda1p0_xip",histVec_Mall_Lambda1p0_xip,Vec_M,Vec_L,Vec_label,2,"Proton #xi = #Delta p_{z} / p_{z}","Events",0,1,false,false);
  Vec_M.clear(); Vec_L.clear();

  std::vector< TH1F* > histVec_Mall_Lambda1p0_xip_over_xiLQ_minus_one__inPPS;
  variable = "genAnalyzer/h1_xip_over_xiLQ_minus_one__inPPS";
  histVec_Mall_Lambda1p0_xip_over_xiLQ_minus_one__inPPS.push_back( (TH1F*)fileHwOld->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_xip_over_xiLQ_minus_one__inPPS.push_back( (TH1F*)fileHwNew->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_xip_over_xiLQ_minus_one__inPPS.push_back( (TH1F*)fileHwNewLO->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_xip_over_xiLQ_minus_one__inPPS.push_back( (TH1F*)fileHwNewLuca->Get(variable.c_str())->Clone() );
  Vec_M = {3000, 3000, 3000, 3000};
  Vec_L = {1.0, 1.0, 1.0, 1.0};
  Vec_label = {"Hw Old NLO", "Hw New NLO", "Hw New LO", "Hw New NLO (Luca)"};
  drawHisto(myoutputdir+"/histVec_Mall_Lambda1p0_xip_over_xiLQ_minus_one__inPPS",histVec_Mall_Lambda1p0_xip_over_xiLQ_minus_one__inPPS,Vec_M,Vec_L,Vec_label,4,"#xi(proton)/#xi(LQ) -1 (within PPS acceptance, #xi<0.2)","Events",-1.,5.,false,false);
  Vec_M.clear(); Vec_L.clear();

  std::vector< TH1F* > histVec_Mall_Lambda1p0_xip_over_xiLQ_minus_one__inPPS__plusLep;
  variable = "genAnalyzer/h1_xip_over_xiLQ_minus_one__inPPS__plusLep";
  histVec_Mall_Lambda1p0_xip_over_xiLQ_minus_one__inPPS__plusLep.push_back( (TH1F*)fileHwOld->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_xip_over_xiLQ_minus_one__inPPS__plusLep.push_back( (TH1F*)fileHwNew->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_xip_over_xiLQ_minus_one__inPPS__plusLep.push_back( (TH1F*)fileHwNewLO->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_xip_over_xiLQ_minus_one__inPPS__plusLep.push_back( (TH1F*)fileHwNewLuca->Get(variable.c_str())->Clone() );
  Vec_M = {3000, 3000, 3000, 3000};
  Vec_L = {1.0, 1.0, 1.0, 1.0};
  Vec_label = {"Hw Old NLO", "Hw New NLO", "Hw New LO", "Hw New NLO (Luca)"};
  drawHisto(myoutputdir+"/histVec_Mall_Lambda1p0_xip_over_xiLQ_minus_one__inPPS__plusLep",histVec_Mall_Lambda1p0_xip_over_xiLQ_minus_one__inPPS__plusLep,Vec_M,Vec_L,Vec_label,1,"#xi(proton)/#xi(LQ+Lepton) -1 (within PPS acceptance, #xi<0.2)","Events (normalized)",-0.5,0.5,false,true);
  Vec_M.clear(); Vec_L.clear();

  std::vector< TH1F* > histVec_Mall_Lambda1p0_gamma_pt__pos;
  variable = "genAnalyzer/h1_gamma_pt__pos";
  histVec_Mall_Lambda1p0_gamma_pt__pos.push_back( (TH1F*)fileHwOld->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_gamma_pt__pos.push_back( (TH1F*)fileHwNew->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_gamma_pt__pos.push_back( (TH1F*)fileHwNewLO->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_gamma_pt__pos.push_back( (TH1F*)fileHwNewLuca->Get(variable.c_str())->Clone() );
  Vec_M = {3000, 3000, 3000, 3000};
  Vec_L = {1.0, 1.0, 1.0, 1.0};
  Vec_label = {"Hw Old NLO", "Hw New NLO", "Hw New LO", "Hw New NLO (Luca)"};
  drawHisto(myoutputdir+"/histVec_Mall_Lambda1p0_gamma_pt__pos",histVec_Mall_Lambda1p0_gamma_pt__pos,Vec_M,Vec_L,Vec_label,2,"Photon pt [GeV]","Events",0,10,false,false);
  Vec_M.clear(); Vec_L.clear();

  std::vector< TH1F* > histVec_Mall_Lambda1p0_gamma_pz__pos;
  variable = "genAnalyzer/h1_gamma_pz__pos";
  histVec_Mall_Lambda1p0_gamma_pz__pos.push_back( (TH1F*)fileHwOld->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_gamma_pz__pos.push_back( (TH1F*)fileHwNew->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_gamma_pz__pos.push_back( (TH1F*)fileHwNewLO->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_gamma_pz__pos.push_back( (TH1F*)fileHwNewLuca->Get(variable.c_str())->Clone() );
  Vec_M = {3000, 3000, 3000, 3000};
  Vec_L = {1.0, 1.0, 1.0, 1.0};
  Vec_label = {"Hw Old NLO", "Hw New NLO", "Hw New LO", "Hw New NLO (Luca)"};
  drawHisto(myoutputdir+"/histVec_Mall_Lambda1p0_gamma_pz__pos",histVec_Mall_Lambda1p0_gamma_pz__pos,Vec_M,Vec_L,Vec_label,2,"Photon pz [GeV]","Events",-30000,30000,false,false);
  Vec_M.clear(); Vec_L.clear();

  std::vector< TH1F* > histVec_Mall_Lambda1p0_LQreco_mass;
  variable = "genAnalyzer/h1_LQreco_mass";
  histVec_Mall_Lambda1p0_LQreco_mass.push_back( (TH1F*)fileHwOld->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_LQreco_mass.push_back( (TH1F*)fileHwNew->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_LQreco_mass.push_back( (TH1F*)fileHwNewLO->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_LQreco_mass.push_back( (TH1F*)fileHwNewLuca->Get(variable.c_str())->Clone() );
  Vec_M = {3000, 3000, 3000, 3000};
  Vec_L = {1.0, 1.0, 1.0, 1.0};
  Vec_label = {"Hw Old NLO", "Hw New NLO", "Hw New LO", "Hw New NLO (Luca)"};
  drawHisto(myoutputdir+"/histVec_Mall_Lambda1p0_LQreco_mass",histVec_Mall_Lambda1p0_LQreco_mass,Vec_M,Vec_L,Vec_label,2,"LQreco mass [GeV]","Events",0,5000,false,false);
  Vec_M.clear(); Vec_L.clear();

  std::vector< TH1F* > histVec_Mall_Lambda1p0_jet_pt;
  variable = "genAnalyzer/h1_jet_pt";
  histVec_Mall_Lambda1p0_jet_pt.push_back( (TH1F*)fileHwOld->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_jet_pt.push_back( (TH1F*)fileHwNew->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_jet_pt.push_back( (TH1F*)fileHwNewLO->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_jet_pt.push_back( (TH1F*)fileHwNewLuca->Get(variable.c_str())->Clone() );
  Vec_M = {3000, 3000, 3000, 3000};
  Vec_L = {1.0, 1.0, 1.0, 1.0};
  Vec_label = {"Hw Old NLO", "Hw New NLO", "Hw New LO", "Hw New NLO (Luca)"};
  drawHisto(myoutputdir+"/histVec_Mall_Lambda1p0_jet_pt",histVec_Mall_Lambda1p0_jet_pt,Vec_M,Vec_L,Vec_label,2,"Jet pt [GeV]","Events",0,3000,false,false);
  Vec_M.clear(); Vec_L.clear();

  std::vector< TH1F* > histVec_Mall_Lambda1p0_jet_eta;
  variable = "genAnalyzer/h1_jet_eta";
  histVec_Mall_Lambda1p0_jet_eta.push_back( (TH1F*)fileHwOld->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_jet_eta.push_back( (TH1F*)fileHwNew->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_jet_eta.push_back( (TH1F*)fileHwNewLO->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_jet_eta.push_back( (TH1F*)fileHwNewLuca->Get(variable.c_str())->Clone() );
  Vec_M = {3000, 3000, 3000, 3000};
  Vec_L = {1.0, 1.0, 1.0, 1.0};
  Vec_label = {"Hw Old NLO", "Hw New NLO", "Hw New LO", "Hw New NLO (Luca)"};
  drawHisto(myoutputdir+"/histVec_Mall_Lambda1p0_jet_eta",histVec_Mall_Lambda1p0_jet_eta,Vec_M,Vec_L,Vec_label,2,"Jet #eta","Events",-5,5,false,false);
  Vec_M.clear(); Vec_L.clear();

  std::vector< TH1F* > histVec_Mall_Lambda1p0_lepFromGamma_pt;
  variable = "genAnalyzer/h1_lepFromGamma_pt";
  histVec_Mall_Lambda1p0_lepFromGamma_pt.push_back( (TH1F*)fileHwOld->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_lepFromGamma_pt.push_back( (TH1F*)fileHwNew->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_lepFromGamma_pt.push_back( (TH1F*)fileHwNewLO->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_lepFromGamma_pt.push_back( (TH1F*)fileHwNewLuca->Get(variable.c_str())->Clone() );
  Vec_M = {3000, 3000, 3000, 3000};
  Vec_L = {1.0, 1.0, 1.0, 1.0};
  Vec_label = {"Hw Old NLO", "Hw New NLO", "Hw New LO", "Hw New NLO (Luca)"};
  drawHisto(myoutputdir+"/histVec_Mall_Lambda1p0_lepFromGamma_pt",histVec_Mall_Lambda1p0_lepFromGamma_pt,Vec_M,Vec_L,Vec_label,2,"Lepton from Gamma pt [GeV]","Events",0,3000,false,false);
  Vec_M.clear(); Vec_L.clear();

  std::vector< TH1F* > histVec_Mall_Lambda1p0_lepFromGamma_eta;
  variable = "genAnalyzer/h1_lepFromGamma_eta";
  histVec_Mall_Lambda1p0_lepFromGamma_eta.push_back( (TH1F*)fileHwOld->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_lepFromGamma_eta.push_back( (TH1F*)fileHwNew->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_lepFromGamma_eta.push_back( (TH1F*)fileHwNewLO->Get(variable.c_str())->Clone() );
  histVec_Mall_Lambda1p0_lepFromGamma_eta.push_back( (TH1F*)fileHwNewLuca->Get(variable.c_str())->Clone() );
  Vec_M = {3000, 3000, 3000, 3000};
  Vec_L = {1.0, 1.0, 1.0, 1.0};
  Vec_label = {"Hw Old NLO", "Hw New NLO", "Hw New LO", "Hw New NLO (Luca)"};
  drawHisto(myoutputdir+"/histVec_Mall_Lambda1p0_lepFromGamma_eta",histVec_Mall_Lambda1p0_lepFromGamma_eta,Vec_M,Vec_L,Vec_label,2,"Lepton from Gamma #eta","Events",-10,10,false,false);
  Vec_M.clear(); Vec_L.clear();


  /*

  // Efficiencies

  //dijet
  std::vector< TH1F* > histVec_Numerator_dijet;
  std::vector< TH1F* > histVec_Denominator_dijet;
  string variable_num = "genAnalyzer/h1_Mjj_sel";
  string variable_den = "genAnalyzer/h1_Mjj";

  for(int k=0; k<Vec_Files.size(); k++)
    {
      histVec_Denominator_dijet.push_back( (TH1F*)Vec_Files[k]->Get(variable_den.c_str())->Clone() );
      histVec_Numerator_dijet.push_back( (TH1F*)Vec_Files[k]->Get(variable_num.c_str())->Clone() );
    }
  drawEff("/afs/cern.ch/user/s/santanas/www/TrijetResoannces13TeV/GenStudies_ggg_madgraph/eff_dijet",histVec_Denominator_dijet,histVec_Numerator_dijet,Vec_Mvalues,Vec_Rvalues,1,"MRes1 [GeV]","R","Dijet selection");

  //trijet
  std::vector< TH1F* > histVec_Numerator_trijet;
  std::vector< TH1F* > histVec_Denominator_trijet;
  variable_num = "genAnalyzer/h1_Mjjj_sel";
  variable_den = "genAnalyzer/h1_Mjjj";

  for(int k=0; k<Vec_Files.size(); k++)
    {
      histVec_Denominator_trijet.push_back( (TH1F*)Vec_Files[k]->Get(variable_den.c_str())->Clone() );
      histVec_Numerator_trijet.push_back( (TH1F*)Vec_Files[k]->Get(variable_num.c_str())->Clone() );
    }
  drawEff("/afs/cern.ch/user/s/santanas/www/TrijetResoannces13TeV/GenStudies_ggg_madgraph/eff_trijet",histVec_Denominator_trijet,histVec_Numerator_trijet,Vec_Mvalues,Vec_Rvalues,1,"MRes1 [GeV]","R","Trijet selection");

  */
 
}




