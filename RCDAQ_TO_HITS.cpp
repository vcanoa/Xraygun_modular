#include <TMath.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <pmonitor/pmonitor.h>
#include <TTree.h>
#include "TCanvas.h"
#include <Event/Event.h>
#include <Event/EventTypes.h> 
#include <TProfile.h>
#include "TF1.h"
#include <iostream>
#include <fstream>

#include "ConfigurationCard.C"

using namespace std;
int nEvents = 0;


//==========================
// GENERAL CONTROL VARIABLES
const float kControl_NoiseCut = 80.0; // value used to decide if channel was hit or not (centered RMS = StdDev)
const float kControl_RejectFlag1 = 800.; // hardwire-value for flaging bad events (type1) check code below
const int   kControl_SpyChannel = 113; //113;//95; // spy channel for wave inspection
int   kControl_EventNumber = 0;
bool  kControl_Sgn[kControl_NumberOfChannels]; // Channel shows activity?
TProfile *kControl_QA_RAWMEANpc    = new TProfile( "kControl_QA_RAWMEANpc", "kControl_QA_RAWMEANpc;channel;meanRawADC",
						   kControl_NumberOfChannels, -0.5, kControl_NumberOfChannels-0.5 );
TH2F     *kControl_QA_RAWRMSpc     = new TH2F( "kControl_QA_RAWRMSpc",  "kControl_QA_RAWRMSpc;channel;stddevRawADC",
					       kControl_NumberOfChannels, -0.5, kControl_NumberOfChannels-0.5, 100, 0, 1000 );
TH1F     *kControl_QA_RAWWAVEhits  = new TH1F( "kControl_QA_RAWWAVEhits",  "kControl_QA_RAWWAVEhits;sample",
					       kControl_NumberofSamples, -0.5, kControl_NumberofSamples-0.5 );
TH1F     *kControl_QA_RAWWAVEnoise = new TH1F( "kControl_QA_RAWWAVEnoise", "kControl_QA_RAWWAVEnoise;sample",
					       kControl_NumberofSamples, -0.5, kControl_NumberofSamples-0.5 );  
TH2F     *kControl_QA_RAWWAVEspy   = new TH2F( "kControl_QA_RAWWAVEspy", Form("kControl_QA_RAWWAVEspy chn=%d;event id;sample",
									      kControl_SpyChannel),
					       100, -0.5, 99.5, kControl_NumberofSamples, -0.5, kControl_NumberofSamples-0.5 );

//==========================
// PEDESTAL CONTROL VARIABLES
const int   kPedestal_NumberOfSamplesForFastPedestal = 5;
float kPedestal_MEAN[kControl_NumberOfChannels];
TProfile *kPedestal_QA_MEANpc = new TProfile( "kPedestal_QA_MEANpc",    "kPedestal_QA_MEANpc;channel;meanADC",
					      kControl_NumberOfChannels, -0.5, kControl_NumberOfChannels-0.5 );

//==========================
// HIT CONTROL VARIABLES
float kHit_PI[kControl_NumberOfChannels];
float kHit_PH[kControl_NumberOfChannels];
float kHit_PHT[kControl_NumberOfChannels];
float kHit_PHA[kControl_NumberOfChannels];
float kHit_PHAT[kControl_NumberOfChannels];
TH1F *kHit_QA_PI    = new TH1F( "kHit_QA_PI",   "kHit_QA_PI;pulse integral", 400, -100, 100 );
TH1F *kHit_QA_PH    = new TH1F( "kHit_QA_PH",   "kHit_QA_PH;pulse MAX height", 400, -100, 2000 );
TH1F *kHit_QA_PHT   = new TH1F( "kHit_QA_PHT",  "kHit_QA_PHT;pulse MAX height trigger",
				kControl_NumberofSamples, -0.5, kControl_NumberofSamples-0.5 );
TH2F *kHit_QA_PHPHT = new TH2F( "kHit_QA_PHPHT","kHit_QA_PHPHT;pulse MAX height;pulse MAX height trigger",
				100,-100,2000,kControl_NumberofSamples, -0.5, kControl_NumberofSamples-0.5);
TH1F *kHit_QA_PHA   = new TH1F( "kHit_QA_PHA",  "kHit_QA_PHA;pulse AVG height", 400, -100, 2000 );
TH1F *kHit_QA_PHAT  = new TH1F( "kHit_QA_PHAT", "kHit_QA_PHAT;pulse AVG height trigger",
				kControl_NumberofSamples, -0.5, kControl_NumberofSamples-0.5 );
TProfile *kHit_QA_CHANNELMAP= new TProfile( "kHit_QA_CHANNELMAP", "kHit_QA_CHANNELMAP;channel",
					    kControl_NumberOfChannels, -0.5, kControl_NumberOfChannels-0.5 );

//===================
int pinit() {
  for (int i=0; i<kControl_NumberOfChannels; i++) {
    kPedestal_MEAN[i]=0;
    kHit_PH[i]=0;
    kHit_PHT[i]=0;
    kControl_Sgn[i]=0;
    kHit_PHA[i]=0;
    kHit_PHAT[i]=0;
    kHit_PI[i]=0;
  }
  return 0;
}

//===================
int phsave(char const*) {
  return 0;
}
//===================
int pcontrol() {
  return 0;
}
//===================
int pcontrol(int) {
  return 0;
}
//========================================================================
int check_signals(Packet *p) {
  int flag = 0; // good event
  int nsamp =  p->iValue(0,"NSAMPLES");
  for(int idx=0; idx!=kControl_NumberOfChannels; ++idx) { // loop over all channels
    kControl_Sgn[ idx ] = false; // clear
    int c = idx%128; // unpacking particulars of SRS: channel
    int h = idx/128; // unpacking particulars of SRS: hybrid
    // computing MEAN and STD-DEV
    double sx = 0;
    double sxx = 0;
    for(int i = 0; i < nsamp; i++) { // loop over samples
      float raw = p->iValue ( c, i, 2*h ); // unpack data
      sx += raw;
      sxx += raw*raw;
    }
    double ped_mean = sx/nsamp;
    double ped_rms = sqrt( sxx/nsamp - ped_mean*ped_mean );
    if(ped_rms>kControl_RejectFlag1) flag = 1; // flag event that contains a big spread in any channel
    kControl_QA_RAWMEANpc->Fill( idx, ped_mean );
    kControl_QA_RAWRMSpc->Fill(  idx, ped_rms );
    // if STD-DEV is higher than kControl_NoiseCut then declare activity in this channel
    kControl_Sgn[ idx ] = ped_rms > kControl_NoiseCut;
    for(int i = 0; i < nsamp; i++) { // loop over samples
      float raw = p->iValue ( c, i, 2*h ); // unpack data
      if( kControl_Sgn[ idx ] ) kControl_QA_RAWWAVEhits->Fill( i, raw );
      else kControl_QA_RAWWAVEnoise->Fill( i, raw );
      if(idx==kControl_SpyChannel)
	kControl_QA_RAWWAVEspy->Fill( kControl_EventNumber, i, raw );
    }
  }
  return flag;
}
//========================================================================
void pedestal_calculation(Packet *p) {
  // The following strategy comes from Bob:
  // For each channel we compute the MEAN and RMS of the adc values for all samples
  // based on the StandardDeviation (RMS for TH1F in old code), we decide if the channel was hit or not.
  // the pedestal is the MEAN when there is no hit
  // when there is a hit, the approach is different
  int nsamp =  p->iValue(0,"NSAMPLES");
  for(int idx=0; idx!=kControl_NumberOfChannels; ++idx) { // loop over all channels
    int c = idx%128; // unpacking particulars of SRS: channel
    int h = idx/128; // unpacking particulars of SRS: hybrid
    // computing MEAN and RMS
    double sx = 0;
    // if there is no activity in this channel compute mean over all samples
    // otherwise just use the mean over few samples around the highest fluctuation
    int highest = 0;
    int imax = 0;
    for(int i = 0; i < nsamp; i++) { // loop over samples
      float raw = p->iValue ( c, i, 2*h ); // unpack data
      sx += raw;
      if(highest < raw) {
	highest = raw;
	imax = i; //saves the positions of the highest fluctuation
      }
    }
    double ped_mean = sx/nsamp;
    if(kControl_Sgn[idx]) { // there is activity in this cell so we take the mean over few samples around the highest fluctuation
      sx = 0;
      int ilow = imax;
      if(ilow>nsamp-kPedestal_NumberOfSamplesForFastPedestal) ilow = nsamp-kPedestal_NumberOfSamplesForFastPedestal;
      for(int i = ilow; i < ilow+kPedestal_NumberOfSamplesForFastPedestal; i++) {
	float raw = p->iValue ( c, i, 2*h ); // unpack data
	sx += raw;
      }
      double ped_mean = sx/kPedestal_NumberOfSamplesForFastPedestal;
    }
    kPedestal_MEAN[ h*128 + c ] = ped_mean;
    kPedestal_QA_MEANpc->Fill( idx, ped_mean );
  }
}
//========================================================================
void hit_calculation(Packet *p) {
  int nsamp =  p->iValue(0,"NSAMPLES");

  fEvent.NHits = 0;
  for(int idx=0; idx!=kControl_NumberOfChannels; ++idx) { // loop over all channels
    kHit_PH[idx] = 0;  // pulse amplitude
    kHit_PHT[idx] = 0; // pulse trigger time
    kHit_PHA[idx] = 0; // pulse amplitude wrt three samples
    kHit_PHAT[idx] = 0;// pulse trigger time wrt three samples
    kHit_PI[idx] = 0;  // pulse integral over full window
    if(!kControl_Sgn[idx]) continue; // discard the channels with no activity above noise
    int c = idx%128; // unpacking particulars of SRS: channel
    int h = idx/128; // unpacking particulars of SRS: hybrid
    double ped = kPedestal_MEAN[idx];
    float minRAW = 1E+37; // initialize to a riduculus high value
    int whenHappened = 0;
    for(int i = 0; i < nsamp; i++) { // loop over samples
      float raw = p->iValue ( c, i, 2*h ); // unpack data
      kHit_PI[idx] += (ped - raw);
      if(raw<minRAW) {
	minRAW = raw;
	whenHappened = i;
      }
    }
    kHit_PH[idx] = ped - minRAW;
    kHit_PHT[idx] = whenHappened;
    int iSample = whenHappened;
    if(iSample < 1 || iSample > nsamp - 1 )
      continue; // discard wave where maximum happens to close to edges
    // computing height as average of three samples around peak
    float val_nm1 = ped - p->iValue ( c, iSample-1, 2*h ); // unpack data
    float val_n   = ped - p->iValue ( c, iSample,   2*h ); // unpack data
    float val_np1 = ped - p->iValue ( c, iSample+1, 2*h ); // unpack data
    kHit_PHA[idx] = val_nm1 + val_n + val_np1;
    kHit_PHAT[idx] = ((iSample-1)*val_nm1 + iSample*val_n + (iSample+1)*val_np1) / kHit_PHA[idx];
    kHit_PHA[idx] /= 3;
    kHit_QA_PI->Fill( kHit_PI[idx] );
    kHit_QA_PH->Fill( kHit_PH[idx] );
    kHit_QA_PHA->Fill( kHit_PHA[idx] );
    kHit_QA_PHT->Fill( kHit_PHT[idx] );
    kHit_QA_PHPHT->Fill( kHit_PH[idx], kHit_PHT[idx] );
    kHit_QA_PHAT->Fill( kHit_PHAT[idx] );
    kHit_QA_CHANNELMAP->Fill( idx, kHit_PH[idx] );
    //publish maxAmplitude and maxAmplitude time
    //should decide depending on the system if to change to a diffferent quantification of energy/time
    //for the time being for SRS this is good enough
    fEvent.Channel[ fEvent.NHits ] = idx;
    fEvent.Energy[ fEvent.NHits ] = kHit_PH[idx];
    fEvent.Time[ fEvent.NHits ] = kHit_PHT[idx];
    fEvent.NHits++;
  }
  fOutput->Fill();
}
//=======================================================================
int process_event (Event *e){
  if ( e->getEvtType() == 12) {
    cout<<"N_EvtProc: "<<nEvents<<", Event: "<<e->getEvtSequence()<<endl;
  }
  //-------Reset all N counters at beginning of run--------
  if ( e->getEvtType() == 9) { // This if statement must come before cutting on event type below---  9=begin event, 1=all middle events, 12=last event
    nEvents = 0;
  }
  //************************************************************************************
  if ( e->getEvtType() != 1) return 0;  //9=begin event, 1=all middle events, 12=last event //MOVE DOWN: No other event types will be accepted below this line!!!???
  //************************************************************************************
  if ( e->getEvtSequence() <0) return 0;  //skip first n events: e->getEvtSequence() < n  //Bob
  //************************************************************************************
  kControl_EventNumber = e->getEvtSequence();
  //************************************************************************************
  nEvents++;
  if(nEvents%1000 == 0) {
    cout << "-----Events_processed: " << nEvents << ", Event_Number: " << kControl_EventNumber <<endl;
  }
  //here we fill histograms and do things
  Packet *p = e->getPacket(kControl_Packet);  //this is the package for SRS
  if (p) {
    int ret = check_signals(p);
    if(ret!=0) {
      delete p;
      return 0; // skip event if passes basic cuts!
    }
    pedestal_calculation(p);
    // commom_mode_calculation(p); // CALLS pulse_height(p) already
    hit_calculation(p);
    // DONE
    delete p;
  }
  return 0;
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
    std::cout << "  ./rcdaq2hits /path/to/data/from/rcdaq/(file.evt)" << std::endl << std::endl;
    std::cout << "output:" << std::endl;
    std::cout << "  ./file.evt_HITS.root" << std::endl;
    return 1;
  }
  TString file_in=arg[1];
  TObjArray *objarray = file_in.Tokenize("/");
  TString filename = ((TObjString*)objarray->At( objarray->GetEntries()-1 ))->GetString();
  filename += "_HITS.root";
  std::cout << file_in.Data() << "  =>  " << filename.Data() << std::endl;
  pfileopen( file_in.Data() );

  CreateTree();
  prun();

  //writing histograms 
  TFile *outputfile = new TFile( filename.Data(), "RECREATE" );
  fOutput->Write();
  kControl_QA_RAWMEANpc->Write();
  kControl_QA_RAWRMSpc->Write();
  kControl_QA_RAWWAVEhits->Write();
  kControl_QA_RAWWAVEnoise->Write();
  kControl_QA_RAWWAVEspy->Write();
  kPedestal_QA_MEANpc->Write();
  kHit_QA_PI->Write();
  kHit_QA_PH->Write();
  kHit_QA_PHT->Write();
  kHit_QA_PHPHT->Write();
  kHit_QA_PHA->Write();
  kHit_QA_PHAT->Write();
  kHit_QA_CHANNELMAP->Write();
  outputfile->Close();
 return 0;

}
