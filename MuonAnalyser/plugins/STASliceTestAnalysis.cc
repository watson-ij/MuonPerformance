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

#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/MuonData/interface/MuonDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMAMCStatusDigi.h"
#include "DataFormats/GEMDigi/interface/GEMVfatStatusDigi.h"
#include "DataFormats/GEMDigi/interface/GEMGEBStatusDigi.h"
#include "EventFilter/GEMRawToDigi/interface/AMCdata.h"

#include "DataFormats/MuonData/interface/MuonDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMAMC13EventCollection.h"
#include "DataFormats/GEMDigi/interface/GEMAMCdataCollection.h"
#include "DataFormats/GEMDigi/interface/GEMGEBStatusDigiCollection.h"
#include "DataFormats/GEMDigi/interface/GEMVfatStatusDigiCollection.h"

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

#include "TVector2.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

using namespace std;
using namespace edm;

class STASliceTestAnalysis : public edm::EDAnalyzer {
public:
  explicit STASliceTestAnalysis(const edm::ParameterSet&);
  ~STASliceTestAnalysis();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() override;
  virtual void endJob() override;

  virtual void beginRun(Run const&, EventSetup const&) override;
  virtual void endRun(Run const&, EventSetup const&) override;

  // ----------member data ---------------------------
  edm::EDGetTokenT<GEMRecHitCollection> gemRecHits_;
  //  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;
  edm::EDGetTokenT<reco::VertexCollection> vertexCollection_;
  edm::Service<TFileService> fs;
  edm::EDGetTokenT<reco::TrackCollection> staTracks_;

  MuonServiceProxy* theService_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;
  edm::ESHandle<MagneticField> bField_;
  edm::EDGetTokenT<LumiScalersCollection> lumiScalers_;

  edm::EDGetTokenT<GEMAMC13EventCollection> amc13Event_;
  edm::EDGetTokenT<GEMAMCdataCollection> amcData_;
  edm::EDGetTokenT<GEMGEBStatusDigiCollection> gebStatusCol_;
  edm::EDGetTokenT<GEMVfatStatusDigiCollection> vfatStatusCol_;

  bool checkEtaPartitionGood(const GEMEtaPartition* part);
  
  edm::Handle<GEMAMC13EventCollection> amc13Event;
  edm::Handle<GEMAMCdataCollection> amcData;
  edm::Handle<GEMGEBStatusDigiCollection> gebStatusCol;  
  edm::Handle<GEMVfatStatusDigiCollection> vfatStatusCol;  

  ULong64_t b_event;
  int b_run, b_lumi;
  float b_instLumi;
  int b_firstStrip, b_nStrips, b_chamber, b_layer, b_etaPartition, b_muonQuality, b_bx;
  int b_latency;
  vector<int> b_strips;
  float b_x, b_y, b_z;
  float b_mu_eta, b_mu_phi, b_mu_pt;
  float b_pull_x, b_pull_y, b_res_x, b_res_y;

  int b_amcBx;

  int nEvents, nMuonTotal, nGEMFiducialMuon, nGEMTrackWithMuon, nGEMTotal;
  int b_nMuons, b_nMuonsWithGEMHit, b_nMuonsInMuonTree, b_nHitsInHitTree;
  int b_valid;

  int b_nGEMHits;

  // muon branches
  int m_nGEMhits, m_nCSChits, m_nhits, m_nvalidhits, m_nDThits, m_nRPChits, m_nLosthits, m_nBadhits, m_nGEMTrackRecHits;
  float m_innerHitPos_x, m_innerHitPos_y, m_innerHitPos_z;
  int m_nbounds;
  int m_quality, m_charge;
  float m_pt, m_eta, m_phi, m_chi2;
  float m_dxy, m_dz, m_lambda, m_dsz, m_vx, m_vy, m_vz;
  float m_dxyError, m_dzError;
  int m_hasCSCStation1, m_hasCSCRing1, m_hasCSCRing4;

  union {
    struct {
      vector<float> m_resx, m_resy, m_pullx, m_pully;      // GEMHits included in Muon
      vector<float> in_vfat;      // Propagation only information

      vector<float> in_local_x_inner, in_local_y_inner;
      vector<float> in_globalPhi, in_globalEta, in_nearGemPhi, in_nearGemEta; // global info
      vector<float> in_x, in_y, in_local_x, in_local_y, in_gemx, in_gemy, in_local_gemx, in_local_gemy, in_pullx, in_pully, in_resx, in_resy, in_trkextdx, in_gem_dx;
      vector<float> in_local_x_closetsos, in_local_y_closetsos, in_trkextdx_closetsos, in_trkextdx_inner;

      vector<float> rec_x, rec_y;
      vector<float> rec_gem_dx;
      vector<float> rec_resx;
      vector<float> rec_trkextdx;
      vector<float> rec_resy;
      vector<float> rec_pullx;
      vector<float> rec_pully;
    } m_f;
    vector<float> m_fs[(sizeof(m_f))/(sizeof(vector<float>))];
  };
  
  union {
    struct {
      vector<int> roll, chamber, layer; // hit info
      vector<int> in_roll, in_chamber, in_layer, in_goodEta; // propagation bound info
      vector<int> in_gemNStrips, in_gemFirstStrip, in_strip;
      vector<int> in_VfatQual, in_VfatFlag, in_VfatBc, in_nvfat, in_errorc, in_stuckd, in_VfatDBx, in_bx;
      vector<int> rec_roll, rec_chamber, rec_layer;
      vector<int> in_matchingGem, in_GebFu, in_IsGeb;
      
      vector<int> in_assoc; // associated muon hit reconstructed in track?
    } m_i;
    vector<int> m_is[(sizeof(m_i))/(sizeof(vector<int>))];
  };

  // GEM RecHit additional info about quality
  bool b_gebFu = false, b_isGeb = false;
  int b_vfatQual = 0, b_vfatFlag = 0, b_vfatBc = 0, b_nvfat = 0, b_stuckd = -99, b_errorc = -99;

  // dimuon
  int d_mu1, d_mu2; /// indexes
  float d_mass, d_pt, d_eta, d_phi; /// dimuon mass
  float d_dca; // distance of closest approach of the two tracks
  
  TTree *t_hit;
  TTree *t_run;
  TTree *t_muon;
  TTree *t_event;
  TTree *t_dimuon;
};

STASliceTestAnalysis::STASliceTestAnalysis(const edm::ParameterSet& iConfig) :
  nEvents(0),
  nMuonTotal(0),
  nGEMFiducialMuon(0),
  nGEMTrackWithMuon(0),
  nGEMTotal(0)
{
  for (unsigned i = 0; i < (sizeof(m_f)/sizeof(m_fs[0])); ++i)
    new(&m_fs[i]) std::vector<float>();
  for (unsigned i = 0; i < (sizeof(m_i)/sizeof(m_is[0])); ++i)
    new(&m_is[i]) std::vector<int>();

  
  
  gemRecHits_ = consumes<GEMRecHitCollection>(iConfig.getParameter<edm::InputTag>("gemRecHits"));
  staTracks_ = consumes<reco::TrackCollection>(iConfig.getParameter<InputTag>("muons"));
  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService_ = new MuonServiceProxy(serviceParameters);
  lumiScalers_ = consumes<LumiScalersCollection>(iConfig.getParameter<edm::InputTag>("lumiScalers"));
  amc13Event_ = consumes<GEMAMC13EventCollection>(iConfig.getParameter<edm::InputTag>("amc13Event"));
  amcData_ = consumes<GEMAMCdataCollection>(iConfig.getParameter<edm::InputTag>("amcData"));
  gebStatusCol_ = consumes<GEMGEBStatusDigiCollection>(iConfig.getParameter<edm::InputTag>("gebStatusCol"));
  vfatStatusCol_ = consumes<GEMVfatStatusDigiCollection>(iConfig.getParameter<edm::InputTag>("vfatStatusCol"));
 
  t_event = fs->make<TTree>("Event", "Event");
  t_event->Branch("nMuons", &b_nMuons, "nMuons/I");
  t_event->Branch("nMuonsWithGEMHit", &b_nMuonsWithGEMHit, "nMuonsWithGEMHit/I");
  t_event->Branch("nGEMHits", &b_nGEMHits, "nGEMHits/I");
  t_event->Branch("nMuonsInMuonTree", &b_nMuonsInMuonTree, "nMuonsInMuonTree/I");
  t_event->Branch("nHitsInHitTree", &b_nHitsInHitTree, "nHitsInHitTree/I");
  t_event->Branch("event", &b_event, "event/l");
  t_event->Branch("run", &b_run, "run/I");
  t_event->Branch("lumi", &b_lumi, "lumi/I");
  t_event->Branch("instLumi", &b_instLumi, "instLumi/F");
  t_event->Branch("latency", &b_latency, "latency/I");
  t_event->Branch("amcBx", &b_amcBx, "amcBx/I");

  t_dimuon = fs->make<TTree>("DiMuon", "DiMuon");
  t_dimuon->Branch("event", &b_event, "event/l");
  t_dimuon->Branch("run", &b_run, "run/I");
  t_dimuon->Branch("lumi", &b_lumi, "lumi/I");
  t_dimuon->Branch("instLumi", &b_instLumi, "instLumi/F");
  t_dimuon->Branch("latency", &b_latency, "latency/I");
  t_dimuon->Branch("amcBx", &b_amcBx, "amcBx/I");
  t_dimuon->Branch("mu1", &d_mu1, "mu1/I");
  t_dimuon->Branch("mu2", &d_mu2, "mu2/I");
  t_dimuon->Branch("mass", &d_mass, "mass/F");
  t_dimuon->Branch("pt", &d_pt, "pt/F");
  t_dimuon->Branch("eta", &d_eta, "eta/F");
  t_dimuon->Branch("phi", &d_phi, "phi/F");
  t_dimuon->Branch("dca", &d_dca, "dca/F");

  t_muon = fs->make<TTree>("Muon", "Muon");
  t_muon->Branch("event", &b_event, "event/l");
  t_muon->Branch("run", &b_run, "run/I");
  t_muon->Branch("lumi", &b_lumi, "lumi/I");
  t_muon->Branch("instLumi", &b_instLumi, "instLumi/F");
  t_muon->Branch("latency", &b_latency, "latency/I");
  t_muon->Branch("amcBx", &b_amcBx, "amcBx/I");
  t_muon->Branch("chi2", &m_chi2, "chi2/F");
  t_muon->Branch("innerHitPos_x", &m_innerHitPos_x, "innerHitPos_x/F");
  t_muon->Branch("innerHitPos_y", &m_innerHitPos_y, "innerHitPos_y/F");
  t_muon->Branch("innerHitPos_z", &m_innerHitPos_z, "innerHit_z/F");
  t_muon->Branch("nhits", &m_nhits, "nhits/I");
  t_muon->Branch("nLosthits", &m_nLosthits, "nLosthits/I");
  t_muon->Branch("nBadhits", &m_nBadhits, "nBadhits/I");
  t_muon->Branch("nCSChits", &m_nCSChits, "nCSChits/I");
  t_muon->Branch("nRPChits", &m_nRPChits, "nRPChits/I");
  t_muon->Branch("nDThits", &m_nDThits, "nDThits/I");
  t_muon->Branch("nGEMhits", &m_nGEMhits, "nGEMhits/I")->SetTitle("n 'valid' GEM hits associated to muon");
  t_muon->Branch("m_nGEMTrackRecHits", &m_nGEMTrackRecHits, "nGEMTrackRecHits/I")->SetTitle("n GEM hits associated to muon");
  t_muon->Branch("nvalidhits", &m_nvalidhits, "nvalidhits/I")->SetTitle("n GEM hits associated to muon, and muon can propagate to eta partition of hit");
  t_muon->Branch("nbounds", &m_nbounds, "nbounds/I")->SetTitle("times muon is in GEM eta partition bounds");
  t_muon->Branch("pt", &m_pt, "pt/F");
  t_muon->Branch("eta", &m_eta, "eta/F");
  t_muon->Branch("phi", &m_phi, "phi/F");
  t_muon->Branch("charge", &m_charge, "charge/I");

  t_muon->Branch("dxy", &m_dxy, "dxy/F");
  t_muon->Branch("dz", &m_dz, "dz/F");
  t_muon->Branch("lambda", &m_lambda, "lambda/F");
  t_muon->Branch("dsz", &m_dsz, "dsz/F");
  t_muon->Branch("vx", &m_vx, "vx/F");
  t_muon->Branch("vy", &m_vy, "vy/F");
  t_muon->Branch("vz", &m_vz, "vz/F");
  t_muon->Branch("dxyError", &m_dxyError, "dxyError/F");
  t_muon->Branch("dzError", &m_dzError, "dzError/F");
  
  t_muon->Branch("quality", &m_quality, "quality/I")->SetTitle("muon quality :: 0:noid 1:looseID 2:tightID");
  // t_muon->Branch("hasME11", &m_hasME11, "hasME11/I");
  t_muon->Branch("hasCSCStation1", &m_hasCSCStation1, "hasCSCStation1/I");
  t_muon->Branch("hasCSCRing1", &m_hasCSCRing1, "hasCSCRing1/I");
  t_muon->Branch("hasCSCRing4", &m_hasCSCRing4, "hasCSCRing4/I");
  t_muon->Branch("in_strip", &m_i.in_strip);
  t_muon->Branch("in_vfat", &m_f.in_vfat);
  t_muon->Branch("in_roll", &m_i.in_roll);
  t_muon->Branch("in_chamber", &m_i.in_chamber);
  t_muon->Branch("in_layer", &m_i.in_layer);
  t_muon->Branch("in_goodEta", &m_i.in_goodEta);
  t_muon->Branch("in_resx", &m_f.in_resx);
  t_muon->Branch("in_resy", &m_f.in_resy);
  t_muon->Branch("in_trkextdx", &m_f.in_trkextdx);
  t_muon->Branch("in_pullx", &m_f.in_pullx);
  t_muon->Branch("in_pully", &m_f.in_pully);
  t_muon->Branch("in_phi", &m_f.in_globalPhi);
  t_muon->Branch("in_eta", &m_f.in_globalEta);
  t_muon->Branch("in_nearGemPhi", &m_f.in_nearGemPhi);
  t_muon->Branch("in_nearGemEta", &m_f.in_nearGemEta);
  t_muon->Branch("in_x", &m_f.in_x);
  t_muon->Branch("in_y", &m_f.in_y);

  t_muon->Branch("in_local_x_inner", &m_f.in_local_x_inner);
  t_muon->Branch("in_local_y_inner", &m_f.in_local_y_inner);
  t_muon->Branch("in_local_x_closetsos", &m_f.in_local_x_closetsos);
  t_muon->Branch("in_local_y_closetsos", &m_f.in_local_y_closetsos);
  t_muon->Branch("in_trkextdx_closetsos", &m_f.in_trkextdx_closetsos);
  t_muon->Branch("in_trkextdx_inner", &m_f.in_trkextdx_inner);

  t_muon->Branch("in_local_x", &m_f.in_local_x);
  t_muon->Branch("in_local_y", &m_f.in_local_y);
  t_muon->Branch("in_gemx", &m_f.in_gemx);
  t_muon->Branch("in_gem_dx", &m_f.in_gem_dx);
  t_muon->Branch("in_gemy", &m_f.in_gemy);
  t_muon->Branch("in_local_gemx", &m_f.in_local_gemx);
  t_muon->Branch("in_local_gemy", &m_f.in_local_gemy);
  t_muon->Branch("in_gemFirstStrip", &m_i.in_gemFirstStrip);
  t_muon->Branch("in_gemNStrips", &m_i.in_gemNStrips);
  t_muon->Branch("in_matchingGem", &m_i.in_matchingGem);
  t_muon->Branch("in_assoc", &m_i.in_assoc);
  t_muon->Branch("in_IsGeb", &m_i.in_IsGeb);
  t_muon->Branch("in_GebFu", &m_i.in_GebFu);
  t_muon->Branch("in_stuckd", &m_i.in_stuckd);
  t_muon->Branch("in_errorc", &m_i.in_errorc);
  t_muon->Branch("in_NVfat", &m_i.in_nvfat);
  t_muon->Branch("in_VfatQual", &m_i.in_VfatQual)->SetTitle("Counts number of *good* Vfat in partition by quality flag");
  t_muon->Branch("in_VfatFlag", &m_i.in_VfatFlag);
  t_muon->Branch("in_VfatBc", &m_i.in_VfatBc)->SetTitle("Counts number of *good* Vfat in partition by asking if BC is equal to AMC BC");
  t_muon->Branch("in_VfatDbx", &m_i.in_VfatDBx)->SetTitle("Sum of AMC BX - VFAT BX for partition");
  t_muon->Branch("in_bx", &m_i.in_bx)->SetTitle("BX of GEM RecHit");

  t_muon->Branch("rec_roll", &m_i.rec_roll);
  t_muon->Branch("rec_chamber", &m_i.rec_chamber);
  t_muon->Branch("rec_layer", &m_i.rec_layer);
  t_muon->Branch("rec_resx", &m_f.rec_resx);
  t_muon->Branch("rec_resy", &m_f.rec_resy);
  t_muon->Branch("rec_trkextdx", &m_f.rec_trkextdx);
  t_muon->Branch("rec_pullx", &m_f.rec_pullx);
  t_muon->Branch("rec_pully", &m_f.rec_pully);
  t_muon->Branch("rec_gem_dx", &m_f.rec_gem_dx);
  t_muon->Branch("rec_x", &m_f.rec_x);
  t_muon->Branch("rec_y", &m_f.rec_y);
  
  t_hit = fs->make<TTree>("Hit", "Hit");
  t_hit->Branch("event", &b_event, "event/l");
  t_hit->Branch("run", &b_run, "run/I");
  t_hit->Branch("lumi", &b_lumi, "lumi/I");
  t_hit->Branch("instLumi", &b_instLumi, "instLumi/F");
  t_hit->Branch("latency", &b_latency, "latency/I");
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

  t_hit->Branch("gebFu", &b_gebFu, "gebFu/O");
  t_hit->Branch("isGeb", &b_isGeb, "isGeb/O");
  
  t_hit->Branch("vfatQual", &b_vfatQual, "vfatQual/I");
  t_hit->Branch("vfatFlag", &b_vfatFlag, "vfatFlag/I");
  t_hit->Branch("vfatBc", &b_vfatBc, "vfatBc/I");
  t_hit->Branch("nvfat", &b_nvfat, "nvfat/I");
  t_hit->Branch("stuckd", &b_stuckd, "stuckd/I");
  t_hit->Branch("errorc", &b_errorc, "errorc/I");
}

STASliceTestAnalysis::~STASliceTestAnalysis()
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
STASliceTestAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<LumiScalersCollection> lumiScalers;
  iEvent.getByToken(lumiScalers_, lumiScalers);

  edm::ESHandle<GEMGeometry> GEMGeometry_;
  iSetup.get<MuonGeometryRecord>().get(GEMGeometry_);

  edm::Handle<GEMRecHitCollection> gemRecHits;  
  iEvent.getByToken(gemRecHits_, gemRecHits);

  edm::Handle<reco::TrackCollection> muons;
  iEvent.getByToken(staTracks_, muons);

  iEvent.getByToken(amc13Event_, amc13Event);
  iEvent.getByToken(amcData_, amcData);
  iEvent.getByToken(gebStatusCol_, gebStatusCol);
  iEvent.getByToken(vfatStatusCol_, vfatStatusCol);

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttrackBuilder_);
  theService_->update(iSetup);
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");

  b_latency = -1;

  b_amcBx = -99;
  for (auto g : *amcData) {
    for (auto a = g.second.first; a != g.second.second; ++a) {
      b_amcBx = a->bx();
      if (b_latency != -1 && b_latency != a->param1())
 	std::cout << "CHANGING LATENCY - old: " << b_latency << " new: " << a->param1() << std::endl;
      b_latency = a->param1();
      if (b_latency == -1) std::cout << "-1 LATENCY VALUE - " << iEvent.id().event() << " " << iEvent.id().run() << std::endl;
    }
  }

  nEvents++;

  b_nMuons = 0;
  b_nMuonsWithGEMHit = 0;
  b_nGEMHits = 0;
  b_nMuonsInMuonTree = 0;
  b_nHitsInHitTree = 0;
  
  b_run = iEvent.run();
  b_lumi = iEvent.luminosityBlock();
  if (lumiScalers->size() > 0)
    b_instLumi = (lumiScalers->at(0)).instantLumi();
  else
    b_instLumi = -1;
  b_event = iEvent.id().event(); // = event number

 
  //  cout << "---- MUONS --- " << endl;
  vector<TLorentzVector> mu_p;
  vector<int> mu_charge;
  vector<int> mu_idx;
  for (auto & mu : *muons) {
    b_nMuons++;

    
    TLorentzVector p; p.SetPtEtaPhiM(m_pt, m_eta, m_phi, 0.1056583745);
    mu_p.push_back(p);
    mu_charge.push_back(m_charge);
    mu_idx.push_back(-1);
    // only consider muons going in the right direction (toward the gem slice test)
    if (mu.eta() > 0) continue;
    
    m_nGEMhits = 0; m_nGEMTrackRecHits = 0;
    m_nbounds = 0;

    for (unsigned i = 0; i < (sizeof(m_f)/sizeof(m_fs[0])); ++i)
      m_fs[i].clear();
    for (unsigned i = 0; i < (sizeof(m_i)/sizeof(m_is[0])); ++i)
      m_is[i].clear();

    if (m_nGEMhits != 0) std::cout << "NOOOOO!!!" << std::endl;
    
    m_quality = mu.qualityMask();
    m_charge = mu.charge();
    // m_hasME11 = 0;
    m_hasCSCRing4 = 0; m_hasCSCRing1 = 0; m_hasCSCStation1 = 0;

    const reco::Track* muonTrack = &mu;
    if (muonTrack) {
      bool fidMu = false;
      reco::TransientTrack ttTrack = ttrackBuilder_->build(muonTrack);
      for (auto ch : GEMGeometry_->etaPartitions()) {
  	TrajectoryStateOnSurface tsos = propagator->propagate(ttTrack.innermostMeasurementState(),
  							      ch->surface());
  	if (!tsos.isValid()) continue;

  	GlobalPoint tsosGP = tsos.globalPosition();
	
  	const LocalPoint pos = ch->toLocal(tsosGP);
  	const LocalPoint pos2D(pos.x(), pos.y(), 0);
  	const BoundPlane& bps(ch->surface());

  	if (bps.bounds().inside(pos2D)) {
  	  fidMu = true;
  	  m_nbounds++;
	  
  	  auto gemid = ch->id();

  	  m_i.in_strip.push_back(ch->strip(pos2D));
  	  m_f.in_vfat.push_back((float) ch->strip(pos2D)/128.);
  	  m_i.in_roll.push_back(gemid.roll());
  	  m_i.in_chamber.push_back(gemid.chamber());
  	  m_i.in_layer.push_back(gemid.layer());

	  m_i.in_goodEta.push_back(checkEtaPartitionGood(ch));
	  
  	  float in_x, in_y;
  	  in_x = tsosGP.x(); in_y = tsosGP.y();
  	  m_f.in_x.push_back(in_x); m_f.in_y.push_back(in_y);
  	  m_f.in_globalPhi.push_back(tsosGP.phi());
  	  m_f.in_globalEta.push_back(tsosGP.eta());

  	  // propagate from the "closest" tsos
  	  float in_local_x_closetsos = -9999, in_local_y_closetsos = -9999, in_trkextdx_closetsos = -99;
  	  GlobalPoint gemtagp1 = ch->position();
  	  TrajectoryStateOnSurface tsos1 = ttTrack.stateOnSurface(gemtagp1);
  	  if (!tsos1.isValid()) goto out1;
  	  tsos1 = propagator->propagate(tsos1, ch->surface());
  	  if (!tsos1.isValid()) goto out1;
  	  in_local_x_closetsos = tsos1.localPosition().x(); in_local_y_closetsos = tsos1.localPosition().y();	  
  	  in_trkextdx_closetsos = std::sqrt(tsos1.localError().positionError().xx());
  	out1:
  	  m_f.in_local_x_closetsos.push_back(in_local_x_closetsos); m_f.in_local_y_closetsos.push_back(in_local_y_closetsos);
  	  m_f.in_trkextdx_closetsos.push_back(in_trkextdx_closetsos);

  	  // Propagate from the innermost state
  	  float in_local_x_inner = -9999, in_local_y_inner = -9999, in_trkextdx_inner = -99;
  	  TrajectoryStateOnSurface tsos2 = propagator->propagate(ttTrack.innermostMeasurementState(),
  								ch->surface());
  	  if (!tsos2.isValid()) goto out2;	  
  	  in_local_x_inner = tsos2.localPosition().x(); in_local_y_inner = tsos2.localPosition().y();
  	  in_trkextdx_inner = std::sqrt(tsos2.localError().positionError().xx());
  	out2:
  	  m_f.in_local_x_inner.push_back(in_local_x_inner); m_f.in_local_y_inner.push_back(in_local_y_inner);
  	  m_f.in_trkextdx_inner.push_back(in_trkextdx_inner);
	  
  	  bool matchingGem = false;
  	  int nstrips = -1, firststrip = -1;
  	  float in_resx = -99, in_trkextdx, in_resy = -99, in_pullx = -99, in_pully = -99, in_gem_dx = -99;
  	  float gemEta = +99.0, gemPhi = +99.0, in_gemx=-999, in_gemy=-999, in_local_gemx=-999, in_local_gemy=-999, in_local_x=-999, in_local_y=-999;
  	  auto recHitsRange = gemRecHits->get(gemid);
  	  auto gemRecHit = recHitsRange.first;

  	  in_local_x = pos.x();
  	  in_local_y = pos.y();

  	  LocalPoint loc;

	  // check quality of chamber
	  bool gebFu = false, isGeb = false;
	  int vfatQual = 0, vfatFlag = 0, vfatBc = 0, nvfat = 0, stuckd = -99, errorc = -99, vfatdbx, bx;
	  auto gebs = gebStatusCol->get(gemid.chamberId());
	  for (auto geb = gebs.first; geb != gebs.second; ++geb) {
	    isGeb = true;
	    if (int(geb->getInFu()) != 0) gebFu = false;
	    else gebFu = true;
	    stuckd = geb->getStuckd();
	    errorc = geb->getErrorC();
	  }
	  auto vfats = vfatStatusCol->get(gemid); 
	  for (auto vfat = vfats.first; vfat != vfats.second; ++vfat) {
	    nvfat++;
	    vfatdbx += b_amcBx - vfat->bc();
	    if (vfat->bc() == b_amcBx) vfatBc++;
	    if (int(vfat->quality()) == 0) vfatQual++;
	    if (int(vfat->flag()) == 0) vfatFlag++;
	  }
	  
  	  for (auto hit = gemRecHit; hit != recHitsRange.second; ++hit) {
  	    auto gemGlob = ch->toGlobal(hit->localPosition());
  	    if (fabs(gemGlob.phi() - tsosGP.phi()) < gemPhi) {
	      bx = hit->BunchX();
	      
  	      in_gemx = gemGlob.x();
  	      in_gemy = gemGlob.y();

  	      in_local_gemx = hit->localPosition().x();
  	      in_local_gemy = hit->localPosition().y();
	      
  	      matchingGem = true;
  	      firststrip = hit->firstClusterStrip();
  	      nstrips = hit->clusterSize();
  	      gemPhi = gemGlob.phi(); gemEta = gemGlob.eta();

  	      LocalPoint && tsos_localpos = tsos.localPosition();
  	      loc = hit->localPosition();
  	      LocalError && tsos_localerr = tsos.localError().positionError();
  	      LocalPoint && dethit_localpos = hit->localPosition();     
  	      LocalError && dethit_localerr = hit->localPositionError();
  	      in_gem_dx = dethit_localerr.xx();
  	      in_resx = (dethit_localpos.x() - tsos_localpos.x());
  	      in_trkextdx = std::sqrt(tsos.localError().positionError().xx());
  	      in_resy = (dethit_localpos.y() - tsos_localpos.y()); 
  	      in_pullx = (dethit_localpos.x() - tsos_localpos.x()) / 
  		std::sqrt(dethit_localerr.xx() + tsos_localerr.xx());
  	      in_pully = (dethit_localpos.y() - tsos_localpos.y()) / 
  		std::sqrt(dethit_localerr.yy() + tsos_localerr.yy());
  	    }
  	  }

	  // try to find an associated rechit
	  bool in_assoc = false;
	  for (auto hit = mu.recHitsBegin(); hit != mu.recHitsEnd(); ++hit) {
	    if ((*hit)->geographicalId() == gemid) {
	      m_nGEMTrackRecHits++;
	      
	      in_assoc = true;
	    }
	  }

	  m_i.in_assoc.push_back(in_assoc);
  	  m_f.in_gemx.push_back(in_gemx); m_f.in_gemy.push_back(in_gemy);
  	  m_f.in_gem_dx.push_back(in_gem_dx);
  	  m_f.in_local_x.push_back(in_local_x); m_f.in_local_y.push_back(in_local_y);
  	  m_f.in_local_gemx.push_back(in_local_gemx); m_f.in_local_gemy.push_back(in_local_gemy);
  	  m_f.in_resx.push_back(in_resx); m_f.in_trkextdx.push_back(in_trkextdx); m_f.in_resy.push_back(in_resy);
  	  m_f.in_pullx.push_back(in_pullx); m_f.in_pully.push_back(in_pully);
  	  m_i.in_gemNStrips.push_back(nstrips); m_i.in_gemFirstStrip.push_back(firststrip);
  	  m_i.in_matchingGem.push_back(matchingGem);
	  m_i.in_GebFu.push_back(gebFu); m_i.in_IsGeb.push_back(isGeb); m_i.in_errorc.push_back(errorc); m_i.in_stuckd.push_back(stuckd);
	  m_i.in_nvfat.push_back(nvfat); m_i.in_VfatQual.push_back(vfatQual); m_i.in_VfatFlag.push_back(vfatFlag); m_i.in_VfatBc.push_back(vfatBc); m_i.in_VfatDBx.push_back(vfatdbx);
  	  m_f.in_nearGemPhi.push_back(gemPhi); m_f.in_nearGemEta.push_back(gemEta); m_i.in_bx.push_back(bx);
  	}
      }

      m_nhits = mu.hitPattern().numberOfValidMuonHits();
      m_nGEMhits = mu.hitPattern().numberOfValidMuonGEMHits();
      m_nCSChits = mu.hitPattern().numberOfValidMuonCSCHits();
      m_nRPChits = mu.hitPattern().numberOfValidMuonRPCHits();
      m_nDThits = mu.hitPattern().numberOfValidMuonDTHits();
      m_nLosthits = mu.hitPattern().numberOfLostMuonHits();
      m_nBadhits = mu.hitPattern().numberOfBadMuonHits();
      m_innerHitPos_x = mu.innerPosition().x();
      m_innerHitPos_y = mu.innerPosition().y();
      m_innerHitPos_z = mu.innerPosition().z();

      for (auto hit = mu.recHitsBegin(); hit != mu.recHitsEnd(); ++hit) {
	if ((*hit)->geographicalId().det() == DetId::Muon &&
	    (*hit)->geographicalId().subdetId() == MuonSubdetId::GEM) {
	  GEMDetId id((*hit)->geographicalId());
	  m_i.rec_chamber.push_back(id.chamber());
	  m_i.rec_layer.push_back(id.layer());
	  m_i.rec_roll.push_back(id.roll());

	  auto roll = GEMGeometry_->etaPartition(id);
	  // if ((*hit)->det() == 0) std::cout << "oh no!" << std::endl;
	  // TrajectoryStateOnSurface tsos = propagator->propagate(ttTrack.innermostMeasurementState(),
	  // 							*(*hit)->surface());
	  TrajectoryStateOnSurface tsos = propagator->propagate(ttTrack.innermostMeasurementState(), roll->surface());

	  if (!tsos.isValid()) {
	    m_f.rec_resx.push_back(100); m_f.rec_resy.push_back(100);
	    m_f.rec_pullx.push_back(100); m_f.rec_pully.push_back(100);
	    m_f.rec_gem_dx.push_back(100); m_f.rec_trkextdx.push_back(100);
	    m_f.rec_y.push_back(-999); m_f.rec_x.push_back(-999);
	  }

	  m_f.rec_x.push_back(tsos.globalPosition().x());
	  m_f.rec_y.push_back(tsos.globalPosition().y());
	  
	  LocalPoint && tsos_localpos = tsos.localPosition();
	  LocalError && tsos_localerr = tsos.localError().positionError();
	  LocalPoint && dethit_localpos = (*hit)->localPosition();
	  LocalError && dethit_localerr = (*hit)->localPositionError();
	  m_f.rec_gem_dx.push_back(dethit_localerr.xx());
	  m_f.rec_resx.push_back(dethit_localpos.x() - tsos_localpos.x());
	  m_f.rec_trkextdx.push_back(std::sqrt(tsos.localError().positionError().xx()));
	  m_f.rec_resy.push_back(dethit_localpos.y() - tsos_localpos.y());
	  m_f.rec_pullx.push_back((dethit_localpos.x() - tsos_localpos.x()) / 
				  std::sqrt(dethit_localerr.xx() + tsos_localerr.xx()));
	  m_f.rec_pully.push_back((dethit_localpos.y() - tsos_localpos.y()) / 
				  std::sqrt(dethit_localerr.yy() + tsos_localerr.yy()));
	}
      }

      
      if (fidMu) ++nGEMFiducialMuon;
      if (m_nGEMhits) ++nGEMTrackWithMuon;
    }

    if (m_nGEMhits > 0 or m_nbounds > 0) {
      m_pt = mu.pt();
      m_eta = mu.eta();
      m_phi = mu.phi();
      m_chi2 = mu.normalizedChi2();
      m_dxy = mu.dxy();
      m_dz = mu.dz();
      m_dxyError = mu.dxyError();
      m_dzError = mu.dzError();
      m_lambda = mu.lambda();
      m_phi = mu.phi();
      m_dsz = mu.dsz();
      m_vx = mu.vx();
      m_vy = mu.vy();
      m_vz = mu.vz();

      for (auto hit = mu.recHitsBegin(); hit != mu.recHitsEnd(); ++hit) {
	if ((*hit)->geographicalId().det() == DetId::Muon &&
	    (*hit)->geographicalId().subdetId() == MuonSubdetId::CSC) {
  	  CSCDetId cscid = (*hit)->geographicalId();
	  if (cscid.station() == 1) {
	    m_hasCSCStation1 = 1;
	    if (cscid.ring() == 1)
	      m_hasCSCRing1 = 1;
	    if (cscid.ring() == 4)
	      m_hasCSCRing4 = 1;
	  }
	}
      }

      if ( (m_chi2 < 5) &&
	   (fabs(m_dxy/m_dxyError) < 3))
	mu_idx[mu_idx.size()-1] = t_muon->GetEntries();
      t_muon->Fill();
      b_nMuonsInMuonTree++;
    }
  } // Muon Loop

  for (unsigned i1 = 0; i1 < mu_p.size(); ++i1) {
    for (unsigned i2 = i1+1; i2 < mu_p.size(); ++i2) {
      if ((mu_idx[i1] == -1) && (mu_idx[i2] == -1)) continue;
      if ((mu_charge[i1]*mu_charge[i2]) > 0)
	continue;
      d_mu1 = mu_idx[i1];
      d_mu2 = mu_idx[i2];
      TLorentzVector dip = (mu_p[i1]+mu_p[i2]);
      d_mass = dip.M();
      d_pt = dip.Pt();
      d_eta = dip.Eta();
      d_phi = dip.Phi();

      reco::TransientTrack ttTrack1 = ttrackBuilder_->build(muons->at(i1));
      reco::TransientTrack ttTrack2 = ttrackBuilder_->build(muons->at(i2));
      ClosestApproachInRPhi cApp;
      cApp.calculate(ttTrack1.impactPointState(), ttTrack2.impactPointState());

      if (!cApp.status()) {d_dca=99;}
      else {d_dca = cApp.distance();}
      
      t_dimuon->Fill();
    }
  }

  //  cout << "--- GEM HITS ---" << endl;
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

    auto roll = GEMGeometry_->etaPartition(detId);
    auto globalPosition = roll->toGlobal(gem.localPosition());
    b_x = globalPosition.x();
    b_y = globalPosition.y();
    b_z = globalPosition.z();
    

    b_gebFu = false, b_isGeb = false;
    b_vfatQual = 0, b_vfatFlag = 0, b_vfatBc = 0, b_nvfat = 0, b_stuckd = -99, b_errorc = -99;
    auto gebs = gebStatusCol->get(detId.chamberId()); 
    for (auto geb = gebs.first; geb != gebs.second; ++geb) {
      b_isGeb = true;
      if (int(geb->getInFu()) != 0) b_gebFu = false;
      else b_gebFu = true;
      b_stuckd = geb->getStuckd();
      b_errorc = geb->getErrorC();
    }
    auto vfats = vfatStatusCol->get(detId); 
    for (auto vfat = vfats.first; vfat != vfats.second; ++vfat) {
      b_nvfat++;
      if (vfat->bc() == b_amcBx) b_vfatBc++;
      if (int(vfat->quality()) == 0) b_vfatQual++;
      if (int(vfat->flag()) == 0) b_vfatFlag++;
    }

    t_hit->Fill();
    b_nHitsInHitTree++;
    //    cout << "   --- eta: " << globalPosition.eta() << " phi: " << globalPosition.phi()  << " xyz: " << b_x << " " << b_y << " " << b_z << " ch:" << detId.chamber() << " l:" << detId.layer() << " p:" << detId.roll() << endl;
  } // GEM Loop
  //  cout << endl << endl;

  nMuonTotal += b_nMuons;
  nGEMTotal += b_nGEMHits;

  t_event->Fill();
}

bool STASliceTestAnalysis::checkEtaPartitionGood(const GEMEtaPartition* part)
{
  GEMDetId rId = part->id();
  int amcBx = -1;

  for (GEMAMCdataCollection::DigiRangeIterator amcsIt = amcData->begin(); amcsIt != amcData->end(); ++amcsIt){
    auto amcs = (*amcsIt).second;
    for (auto amc = amcs.first; amc != amcs.second; ++amc) {
      amcBx = amc->bx();
    }
  }
  
  auto gebs = gebStatusCol->get(rId.chamberId()); 
  for (auto geb = gebs.first; geb != gebs.second; ++geb) {
    if (int(geb->getInFu()) != 0 ) return false;
  }
  
  auto vfats = vfatStatusCol->get(rId); 
  for (auto vfat = vfats.first; vfat != vfats.second; ++vfat) {
    if (vfat->bc() != amcBx ) return false;
    if (int(vfat->quality()) != 0 ) return false;
    if (int(vfat->flag()) != 0 ) return false;
  }
  return true;
}

void STASliceTestAnalysis::beginJob(){}
void STASliceTestAnalysis::endJob(){}

void STASliceTestAnalysis::beginRun(Run const& run, EventSetup const&){
}
void STASliceTestAnalysis::endRun(Run const&, EventSetup const&){}

//define this as a plug-in
DEFINE_FWK_MODULE(STASliceTestAnalysis);
