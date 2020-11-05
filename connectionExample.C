#include <TFile.h>
#include <TRandom3.h>
#include <TTree.h>
#include <vector>
#include <map>

const int   kControl_NumberOfChannels = 512; // 128*4 (256 + 256)

//============= OUTPUT ============
TTree *fOutput;
struct event {
  // reference position (external tracker coordinate indexes)
  Int_t Ref[2];
  Int_t NHits;
  Int_t   Channel[kControl_NumberOfChannels];
  Float_t Energy[kControl_NumberOfChannels];
  Float_t Time[kControl_NumberOfChannels];
} fEvent;
//=================================
void CreateTree() {
  fOutput = new TTree("out","Triggered");
  fOutput->Branch("Ref",     &fEvent.Ref,     "Ref[2]/I");
  fOutput->Branch("NHits",   &fEvent.NHits,   "NHits/I");
  fOutput->Branch("Channel", &fEvent.Channel, "Channel[NHits]/I");
  fOutput->Branch("Energy",  &fEvent.Energy,  "Energy[NHits]/F");
  fOutput->Branch("Time",    &fEvent.Time,    "Time[NHits]/F");
}
//=================================
int main() {
  CreateTree();
  for(int i=0; i!=100; ++i) {
    fEvent.Ref[0] = 1;
    fEvent.Ref[1] = 2;
    fEvent.NHits = gRandom->Integer(50);
    for(int j=0; j!=fEvent.NHits; ++j) {
      fEvent.Channel[j] = 10;
      fEvent.Energy[j] = 1.0;
      fEvent.Time[j] = 2.0;
    }
    fOutput->Fill();
  }
  TFile *file = new TFile("example.root","RECREATE");
  fOutput->Write();
  file->Close();
  return 0;
}
