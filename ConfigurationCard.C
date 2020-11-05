//============================================
//this should come from configuration file
const int kControl_Packet = 1010;
const int kControl_NumberOfChannels = 512;
const int kControl_NumberofSamples  = 30;
const int kControl_MaxADCvalue      = 4096;
//============================================

//============= OUTPUT ============
TTree *fOutput;
struct event {
  // reference position (external tracker coordinate indexes)
  //Int_t Ref[2];
  Int_t NHits;
  Int_t   Channel[kControl_NumberOfChannels];
  Float_t Energy[kControl_NumberOfChannels];
  Float_t Time[kControl_NumberOfChannels];
} fEvent;
//=================================
void CreateTree() {
  fOutput = new TTree("triggered","Triggered");
  //fOutput->Branch("Ref",     &fEvent.Ref,     "Ref[2]/I");
  fOutput->Branch("NHits",   &fEvent.NHits,   "NHits/I");
  fOutput->Branch("Channel", &fEvent.Channel, "Channel[NHits]/I");
  fOutput->Branch("Energy",  &fEvent.Energy,  "Energy[NHits]/F");
  fOutput->Branch("Time",    &fEvent.Time,    "Time[NHits]/F");
}
//=================================
void ConnectTree() {
  //fOutput->SetBranchAddress("Ref",     &fEvent.Ref);
  fOutput->SetBranchAddress("NHits",   &fEvent.NHits);
  fOutput->SetBranchAddress("Channel", &fEvent.Channel);
  fOutput->SetBranchAddress("Energy",  &fEvent.Energy);
  fOutput->SetBranchAddress("Time",    &fEvent.Time);
}
