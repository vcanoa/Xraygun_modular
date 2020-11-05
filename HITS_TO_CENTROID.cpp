#include <TMath.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include "TCanvas.h"
#include <TProfile.h>
#include <TF1.h>
#include <iostream>
#include <fstream>
#include "ConfigurationCard.C"

using namespace std;

//==========================
// MAP CONTROL VARIABLES 
int map_x[256];
int map_y[256];
const float kPitchX = 0.4;
const float kPitchY = 0.4;
//selecting good cluster
const float x_rms_min=0.5;
const float x_rms_max=3.0;

//==========================
// CLUSTER CONTROL VARIABLES 
TH1F *kCluster_QX = new TH1F("kCluster_QX","kCluster_QX;sumQx;events",1000,0,30000);
TH1F *kCluster_QY = new TH1F("kCluster_QY","kCluster_QY;sumQy;events",1000,0,30000);
TH1F *kCluster_XCENTROID = new TH1F("kCluster_XCENTROID","kCluster_XCENTROID;mm;events",5000,0,100);
TH1F *kCluster_YCENTROID = new TH1F("kCluster_YCENTROID","kCluster_YCENTROID;mm;events",1000,0,100);
TH1F *kCluster_XSTD_DEV = new TH1F("kCluster_XSTD_DEV","kCluster_XSTD_DEV;mm;events",1000,0,10);
TH1F *kCluster_YSTD_DEV = new TH1F("kCluster_YSTD_DEV","kCluster_YSTD_DEV;mm;events",1000,0,10);
TH2F *kCluster_XYCENTROID = new TH2F("kCluster_XYCENTROID","kCluster_XYCENTROID;mm;mm",200,0,100,200,0,100);
TH2F *kCluster_XMOD_STD = new TH2F("kCluster_XMOD_STD","kCluster_XMOD_STD;CENTROID_MOD_PITCH  mm;STD_DEV  mm",200,-kPitchX,+kPitchX,200,0,10);
TH2F *kCluster_YMOD_STD = new TH2F("kCluster_YMOD_STD","kCluster_YMOD_STD;CENTROID_MOD_PITCH  mm;STD_DEV  mm",200,-kPitchY,+kPitchY,200,0,10);

//========================================================================
void load_map() {
  // we should be loading a hit map here instead of hardwire
  // this map, but maybe latter
  for(int i=0; i!=256; ++i) {
    map_x[i] = i*kPitchX - kPitchX/2;
    map_y[i] = i*kPitchY - kPitchY/2;
  }
}
//========================================================================
void cluster_calculation() {
  // making clusters
  double sqx_xx=0;
  double sqx_x=0;
  double sqx=0;
  double sqy_yy=0;
  double sqy_y=0;
  double sqy=0;
  for(int i=0; i!=fEvent.NHits; ++i) {
    int chn = fEvent.Channel[i];
    float ene = fEvent.Energy[i];
    float tim = fEvent.Time[i];
    if( chn < 256 ) { // x-direction
      sqx_xx+= ene * map_x[chn] * map_x[chn]; // chargeX*posX*posX
      sqx_x += ene * map_x[chn]; // chargeX*posX
      sqx   += ene; // chargeX  
    } else { // y-direction
      sqy_yy+= ene * map_y[chn-256] * map_y[chn-256]; // chargeY*posY*posY
      sqy_y += ene * map_y[chn-256]; // chargeY*posY
      sqy   += ene; // chargeY
    }
  }
  double mX = sqx_x / sqx;
  double mY = sqy_y / sqy;
  double mXX = sqx_xx / sqx;
  double mYY = sqy_yy / sqy;
  double sX = sqrt( mXX-mX*mX );
  double sY = sqrt( mYY-mY*mY );

  //selecting clusters
  //if((mX>39 &&mX<49)&& sX>0.5 &&sX<3.0){ //Vero
  //if(sX>0.5 &&sX<3.0){ //Vero
  if(sX>x_rms_min &&sX<x_rms_max){ 
    kCluster_QX->Fill( sqx );
    kCluster_QY->Fill( sqy );
    kCluster_XCENTROID->Fill( mX );
    kCluster_YCENTROID->Fill( mY );
    kCluster_XSTD_DEV->Fill( sX );
    kCluster_YSTD_DEV->Fill( sY );
    kCluster_XYCENTROID->Fill( mX, mY);
    kCluster_XMOD_STD->Fill( fmod(mX,kPitchX), sX );
    kCluster_XMOD_STD->Fill( fmod(mY,kPitchY), sY );
  }
}
//=======================================================================
void resolution_calculation(TString file){
  int peakbin;
  float peakpos;

  peakbin = kCluster_XCENTROID->GetMaximumBin();
  peakpos=kCluster_XCENTROID->GetBinCenter(peakbin);
  TF1 *f1 = new TF1("f1", "gaus");
  f1->SetParLimits(1, peakpos-1.5,peakpos+1.5 ); // search within 1.5 mm
  f1->SetParLimits(2, 0.1 - 0.05, 0.1 + 0.05);
  float minnn = peakpos - 0.5;
  float maxxx = peakpos + 0.5;
  kCluster_XCENTROID->GetXaxis()->SetRangeUser(minnn, maxxx);
  kCluster_XCENTROID->Fit("f1","Q");
  kCluster_XCENTROID->Fit("f1","Q");
  kCluster_XCENTROID->Fit("f1");
  float kResidual_MeanX = f1->GetParameter(1);
  float kResidual_MeanX_er=f1->GetParError(1);
  float kResidual_SigmaX = f1->GetParameter(2);
  float kResidual_SigmaX_er=f1->GetParError(2);

  peakbin = kCluster_YCENTROID->GetMaximumBin();
  peakpos=kCluster_YCENTROID->GetBinCenter(peakbin);
  TF1 *f2 = new TF1("f2", "gaus");
  f1->SetParLimits(1, peakpos-1.5,peakpos+1.5 ); // search within 1.5 mm
  f1->SetParLimits(2, 0.1 - 0.05, 0.1 + 0.05);
  minnn = peakpos - 0.5;
  maxxx = peakpos + 0.5;
  kCluster_YCENTROID->GetXaxis()->SetRangeUser(minnn, maxxx);
  kCluster_YCENTROID->Fit("f1","Q");
  kCluster_YCENTROID->Fit("f1","Q");
  kCluster_YCENTROID->Fit("f1");
  float kResidual_MeanY = f1->GetParameter(1);
  float kResidual_MeanY_er=f1->GetParError(1);
  float kResidual_SigmaY = f1->GetParameter(2);
  float kResidual_SigmaY_er=f1->GetParError(2);

  ofstream file_resul;
  file_resul.open( file.Data() );
  file_resul << file.Data() << " ";
  file_resul << kResidual_MeanX<<" "<<kResidual_MeanX_er<<" "<<kResidual_SigmaX<<" "<<kResidual_SigmaX_er << " ";
  file_resul << kResidual_MeanY<<" "<<kResidual_MeanY_er<<" "<<kResidual_SigmaY<<" "<<kResidual_SigmaY_er<<" \n";
  file_resul.close();
}
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================
int main(int nn, char** arg){
  if(nn!=2) {
    std::cout << "usage:" << std::endl;
    std::cout << "  ./hits2centroid  file.root  [xREF]  [yREF]" << std::endl << std::endl;
    std::cout << "output:" << std::endl;
    std::cout << "  ./file.root_CENTROID.root" << std::endl;
    std::cout << "  ./file.root_CENTROID.txt" << std::endl;
    return 1;
  }
  TString file_in=arg[1];

  //============
  load_map();
  TFile *filein = new TFile( file_in.Data() );
  fOutput = (TTree*) filein->Get("triggered");
  ConnectTree();
  UInt_t nev = fOutput->GetEntries();
  for(UInt_t ev=0; ev!=nev; ++ev) {
    fOutput->GetEntry(ev);
    cluster_calculation();
  }
  resolution_calculation( Form("%s_CENTROID.txt",file_in.Data()) );
  //============
  
  //writing histograms 
  TFile *outhistos = new TFile( Form("%s_CENTROID.root",file_in.Data()), "RECREATE" );
  kCluster_QX->Write();
  kCluster_QY->Write();
  kCluster_XCENTROID->Write();
  kCluster_YCENTROID->Write();
  kCluster_XSTD_DEV->Write();
  kCluster_YSTD_DEV->Write();
  kCluster_XYCENTROID->Write();
  outhistos->Close();

  return 0;

}
