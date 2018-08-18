// cd ../../.. && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval `scramv1 runtime -sh` && eval `scramv1 runtime -sh` && scram b -j 10

// system include files
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include<map>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "DataFormats/MuonDetId/interface/GEMDetId.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartition.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TTree.h"

using namespace std;
using namespace edm;

class AodSliceTestAnalysis : public edm::EDAnalyzer {
public:
  explicit AodSliceTestAnalysis(const edm::ParameterSet&);
  ~AodSliceTestAnalysis();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(Run const&, EventSetup const&) override;
  virtual void endRun(Run const&, EventSetup const&) override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;
  edm::EDGetTokenT<reco::VertexCollection> vertexCollection_;
  edm::Service<TFileService> fs;

  MuonServiceProxy* theService_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;
  edm::ESHandle<MagneticField> bField_; 
  
  TH2D* h_firstStrip[36][3];
  TH2D* h_allStrips[36][3];
  TH2D* h_globalPosOnGem;
  TH1D* h_clusterSize, *h_totalStrips, *h_bxtotal;
  TH1D* h_inEta[36][3];
  TH1D* h_hitEta[36][3];
  TH1D* h_trkEta[36][3];

  TH1D* h_res_x, *h_res_y, *h_pull_x, *h_pull_y;

  int b_run, b_lumi, b_event;
  int b_firstStrip, b_nStrips, b_chamber, b_layer, b_etaPartition, b_muonQuality, b_bx;
  vector<int> b_strips;
  float b_x, b_y, b_z;
  float b_mu_eta, b_mu_phi, b_mu_pt;
  float b_pull_x, b_pull_y, b_res_x, b_res_y;

  int nEvents, nMuonTotal, nGEMFiducialMuon, nGEMTrackWithMuon, nGEMTotal;
  int b_nMuons, b_nMuonsWithGEMHit;
  int b_valid;

  int b_nGEMHits;

  // muon branches
  int m_nhits, m_nvalidhits;
  int m_nbounds;
  int m_quality;
  float m_pt, m_eta, m_phi;
  // GEMHits included in Muon
  vector<int> m_roll, m_chamber, m_layer; // hit info
  vector<float> m_resx, m_resy, m_pullx, m_pully;
  // Propagation only information
  vector<int> m_in_roll, m_in_chamber, m_in_layer; // propagation bound info
  vector<float> m_in_globalPhi, m_in_globalEta, m_in_nearGemPhi, m_in_nearGemEta; // global info
  vector<float> m_in_x, m_in_y, m_in_gemx, m_in_gemy, m_in_pullx, m_in_pully, m_in_resx, m_in_resy;
  vector<int> m_in_nearGemNStrips, m_in_nearGemFirstStrip;
  
  TTree *t_hit;
  TTree *t_run;
  TTree *t_muon;
  TTree *t_event;
};

AodSliceTestAnalysis::AodSliceTestAnalysis(const edm::ParameterSet& iConfig) :
  nEvents(0),
  nMuonTotal(0),
  nGEMFiducialMuon(0),
  nGEMTrackWithMuon(0),
  nGEMTotal(0)
{ 
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
  muons_ = consumes<View<reco::Muon> >(iConfig.getParameter<InputTag>("muons"));
  
  t_event = fs->make<TTree>("Event", "Event");
  t_event->Branch("nMuons", &b_nMuons, "nMuons/I");
  t_event->Branch("nMuonsWithGEMHit", &b_nMuonsWithGEMHit, "nMuonsWithGEMHit/I");
  t_event->Branch("nGEMHits", &b_nGEMHits, "nGEMHits/I");
  t_event->Branch("run", &b_run, "run/I");
  t_event->Branch("lumi", &b_lumi, "lumi/I");

  t_hit = fs->make<TTree>("Hit", "Hit");
  t_hit->Branch("run", &b_run, "run/I");
  t_hit->Branch("lumi", &b_lumi, "lumi/I");
  t_hit->Branch("bx", &b_bx, "bx/I");
  t_hit->Branch("firstStrip", &b_firstStrip, "firstStrip/I");
  t_hit->Branch("nStrips", &b_nStrips, "nStrips/I");
  t_hit->Branch("strips", &b_strips);
  t_hit->Branch("chamber", &b_chamber, "chamber/I");
  t_hit->Branch("layer", &b_layer, "layer/I");
  t_hit->Branch("etaPartition", &b_etaPartition, "etaPartition/I");
  t_hit->Branch("muonQuality", &b_muonQuality, "muonQuality/I")->SetTitle("muonQuality -1:none 0:noid 1:looseID 2:tightID");
  t_hit->Branch("x", &b_x, "x/F");
  t_hit->Branch("y", &b_y, "y/F");
  t_hit->Branch("z", &b_z, "z/F");
}

AodSliceTestAnalysis::~AodSliceTestAnalysis()
{
  std::cout << "::: GEM Slice Test Results :::" << std::endl;
  std::cout << ": From " << nEvents << " events" << std::endl;
  std::cout << std::endl;
  std::cout << " # GEMHits " << nGEMTotal << std::endl;
  std::cout << " # Muons   " << nMuonTotal << std::endl;
  std::cout << " # FidMu   " << nGEMFiducialMuon << std::endl;
  std::cout << " # GEMMu   " << nGEMTrackWithMuon << std::endl;
  std::cout << std::endl;
}

void
AodSliceTestAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<GEMRecHitCollection> gemRecHits;  
  iEvent.getByToken(gemRecHits_, gemRecHits);

  Handle<View<reco::Muon> > muons;
  iEvent.getByToken(muons_, muons);

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttrackBuilder_);

  nEvents++;

  b_nMuons = 0;
  b_nMuonsWithGEMHit = 0;
  b_nGEMHits = 0;
  
  b_run = iEvent.run();
  b_lumi = iEvent.luminosityBlock();

  for (auto & mu : *muons) {
    b_nMuons++;
    
    if (mu.passed(reco::Muon::Selector::CutBasedIdTight)) m_quality = 2;
    else if (mu.passed(reco::Muon::Selector::CutBasedIdLoose)) m_quality = 1;
    else m_quality = 0;

    m_nhits = 0;
    m_nbounds = 0;

    const reco::Track* muonTrack = 0;  
    if (mu.globalTrack().isNonnull()) muonTrack = mu.globalTrack().get();
    else if (mu.outerTrack().isNonnull()) muonTrack = mu.outerTrack().get();
    if (muonTrack) {
      reco::TransientTrack ttTrack = ttrackBuilder_->build(muonTrack);
      
    }
  }

  for (auto & gem : *gemRecHits) {
    ++b_nGEMHits;
    auto detId = gem.gemId();

    b_strips.clear();
    b_chamber = detId.chamber();
    b_layer = detId.layer();
    b_etaPartition = detId.roll();
    b_firstStrip = gem.firstClusterStrip();
    b_bx = gem.BunchX();
    b_nStrips = gem.clusterSize();
    for (int nstrip = gem.firstClusterStrip(); nstrip < (gem.firstClusterStrip() + gem.clusterSize()); ++nstrip) {
      b_strips.push_back(nstrip);
    }

    t_hit->Fill();
  }

  nMuonTotal += b_nMuons;
  nGEMTotal += b_nGEMHits;
}

void AodSliceTestAnalysis::beginJob(){}
void AodSliceTestAnalysis::endJob(){}

void AodSliceTestAnalysis::beginRun(Run const& run, EventSetup const&){
}
void AodSliceTestAnalysis::endRun(Run const&, EventSetup const&){}

//define this as a plug-in
DEFINE_FWK_MODULE(AodSliceTestAnalysis);
