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
  edm::EDGetTokenT<LumiScalersCollection> lumiScalers_;
  
  int b_run, b_lumi, b_event;
  float b_instLumi;
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
  vector<int> m_in_gemNStrips, m_in_gemFirstStrip;
  vector<bool> m_in_matchingGem;
  
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
  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService_ = new MuonServiceProxy(serviceParameters);
  lumiScalers_ = consumes<LumiScalersCollection>(iConfig.getParameter<edm::InputTag>("lumiScalers"));
  
  t_event = fs->make<TTree>("Event", "Event");
  t_event->Branch("nMuons", &b_nMuons, "nMuons/I");
  t_event->Branch("nMuonsWithGEMHit", &b_nMuonsWithGEMHit, "nMuonsWithGEMHit/I");
  t_event->Branch("nGEMHits", &b_nGEMHits, "nGEMHits/I");
  t_event->Branch("run", &b_run, "run/I");
  t_event->Branch("lumi", &b_lumi, "lumi/I");
  t_event->Branch("instLumi", &b_instLumi, "instLumi/F");

  t_muon = fs->make<TTree>("Muon", "Muon");
  t_muon->Branch("run", &b_run, "run/I");
  t_muon->Branch("lumi", &b_lumi, "lumi/I");
  t_muon->Branch("instLumi", &b_instLumi, "instLumi/F");
  t_muon->Branch("nhits", &m_nhits, "nhits/I")->SetTitle("n GEM hits associated to muon");
  t_muon->Branch("nvalidhits", &m_nvalidhits, "nvalidhits/I")->SetTitle("n GEM hits associated to muon, and muon can propagate to eta partition of hit");
  t_muon->Branch("nbounds", &m_nbounds, "nbounds/I")->SetTitle("times muon is in GEM eta partition bounds");
  t_muon->Branch("pt", &m_pt, "pt/F");
  t_muon->Branch("eta", &m_eta, "eta/F");
  t_muon->Branch("phi", &m_phi, "phi/F");
  t_muon->Branch("quality", &m_quality, "quality/I")->SetTitle("muon quality :: 0:noid 1:looseID 2:tightID");
  t_muon->Branch("in_roll", &m_in_roll);
  t_muon->Branch("in_chamber", &m_in_chamber);
  t_muon->Branch("in_layer", &m_in_layer);
  t_muon->Branch("in_resx", &m_in_resx);
  t_muon->Branch("in_resy", &m_in_resy);
  t_muon->Branch("in_pullx", &m_in_pullx);
  t_muon->Branch("in_pully", &m_in_pully);
  t_muon->Branch("in_phi", &m_in_globalPhi);
  t_muon->Branch("in_eta", &m_in_globalEta);
  t_muon->Branch("in_nearGemPhi", &m_in_nearGemPhi);
  t_muon->Branch("in_nearGemEta", &m_in_nearGemEta);
  t_muon->Branch("in_x", &m_in_x);
  t_muon->Branch("in_y", &m_in_y);
  t_muon->Branch("in_gemx", &m_in_gemx);
  t_muon->Branch("in_gemy", &m_in_gemy);
  t_muon->Branch("in_gemFirstStrip", &m_in_gemFirstStrip);
  t_muon->Branch("in_gemNStrips", &m_in_gemNStrips);
  t_muon->Branch("in_matchingGem", &m_in_matchingGem);

  t_hit = fs->make<TTree>("Hit", "Hit");
  t_hit->Branch("run", &b_run, "run/I");
  t_hit->Branch("lumi", &b_lumi, "lumi/I");
  t_hit->Branch("instLumi", &b_instLumi, "instLumi/F");
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
  edm::Handle<LumiScalersCollection> lumiScalers;
  iEvent.getByToken(lumiScalers_, lumiScalers);

  edm::ESHandle<GEMGeometry> GEMGeometry_;
  iSetup.get<MuonGeometryRecord>().get(GEMGeometry_);

  edm::Handle<GEMRecHitCollection> gemRecHits;  
  iEvent.getByToken(gemRecHits_, gemRecHits);

  Handle<View<reco::Muon> > muons;
  iEvent.getByToken(muons_, muons);

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", ttrackBuilder_);
  theService_->update(iSetup);
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");

  nEvents++;

  b_nMuons = 0;
  b_nMuonsWithGEMHit = 0;
  b_nGEMHits = 0;
  
  b_run = iEvent.run();
  b_lumi = iEvent.luminosityBlock();
  b_instLumi = (lumiScalers->at(0)).instantLumi();

  //  cout << "---- MUONS --- " << endl;
  for (auto & mu : *muons) {
    b_nMuons++;
    // only consider muons going in the right direction (toward the gem slice test)
    if (mu.eta() > 0) continue;

    
    m_nhits = 0;
    m_nbounds = 0;
    m_in_roll.clear(); m_in_chamber.clear(); m_in_layer.clear();
    m_in_x.clear(); m_in_y.clear(); m_in_gemx.clear(); m_in_gemy.clear();
    m_in_resx.clear(); m_in_resy.clear(); m_in_pullx.clear(); m_in_pully.clear();
    m_in_matchingGem.clear();
    m_in_gemFirstStrip.clear(); m_in_gemNStrips.clear();
    m_in_nearGemPhi.clear(); m_in_nearGemEta.clear();
    m_in_globalPhi.clear(); m_in_globalEta.clear();
    
    if (mu.passed(reco::Muon::Selector::CutBasedIdTight)) m_quality = 2;
    else if (mu.passed(reco::Muon::Selector::CutBasedIdLoose)) m_quality = 1;
    else m_quality = 0;

    const reco::Track* muonTrack = 0;  
    if (mu.globalTrack().isNonnull()) muonTrack = mu.globalTrack().get();
    else if (mu.outerTrack().isNonnull()) muonTrack = mu.outerTrack().get();
    if (muonTrack) {
      bool fidMu = false;
      reco::TransientTrack ttTrack = ttrackBuilder_->build(muonTrack);
      for (auto ch : GEMGeometry_->etaPartitions()) {
	TrajectoryStateOnSurface tsos = propagator->propagate(ttTrack.outermostMeasurementState(),
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
	  m_in_roll.push_back(gemid.roll());
	  m_in_chamber.push_back(gemid.chamber());
	  m_in_layer.push_back(gemid.layer());
	  float in_x, in_y;
	  in_x = tsosGP.x(); in_y = tsosGP.y();
	  m_in_x.push_back(in_x); m_in_y.push_back(in_y);
	  m_in_globalPhi.push_back(tsosGP.phi());
	  m_in_globalEta.push_back(tsosGP.eta());

	  bool matchingGem = false;
	  int nstrips = -1, firststrip = -1;
	  float in_resx = -99, in_resy = -99, in_pullx = -99, in_pully = -99;
	  float gemEta = +99.0, gemPhi = +99.0, in_gemx=-999, in_gemy=-999;
	  auto recHitsRange = gemRecHits->get(gemid);
	  auto gemRecHit = recHitsRange.first;
	  for (auto hit = gemRecHit; hit != recHitsRange.second; ++hit) {
	    auto gemGlob = ch->toGlobal(hit->localPosition());
	    if (fabs(gemGlob.phi() - tsosGP.phi()) < gemPhi) {
	      in_gemx = gemGlob.x();
	      in_gemy = gemGlob.y();
	      matchingGem = true;
	      firststrip = hit->firstClusterStrip();
	      nstrips = hit->clusterSize();
	      gemPhi = gemGlob.phi(); gemEta = gemGlob.eta();

	      LocalPoint && tsos_localpos = tsos.localPosition();
	      LocalError && tsos_localerr = tsos.localError().positionError();
	      LocalPoint && dethit_localpos = hit->localPosition();     
	      LocalError && dethit_localerr = hit->localPositionError();
	      in_resx = (dethit_localpos.x() - tsos_localpos.x());
	      in_resy = (dethit_localpos.y() - tsos_localpos.y()); 
	      in_pullx = (dethit_localpos.x() - tsos_localpos.x()) / 
		std::sqrt(dethit_localerr.xx() + tsos_localerr.xx());
	      in_pully = (dethit_localpos.y() - tsos_localpos.y()) / 
		std::sqrt(dethit_localerr.yy() + tsos_localerr.yy());
	    }
	  }
	  m_in_gemx.push_back(in_gemx); m_in_gemy.push_back(in_gemy);
	  m_in_resx.push_back(in_resx); m_in_resy.push_back(in_resy);
	  m_in_pullx.push_back(in_pullx); m_in_pully.push_back(in_pully);
	  m_in_gemNStrips.push_back(nstrips); m_in_gemFirstStrip.push_back(firststrip);
	  m_in_matchingGem.push_back(matchingGem);
	  m_in_nearGemPhi.push_back(gemPhi); m_in_nearGemEta.push_back(gemEta);
	  // cout << "   --- eta: " << tsosGP.eta() << " phi: " << tsosGP.phi()  << " xyz: " << in_x << " " << in_y << " GEMxyz: " << in_gemx << " " << in_gemy << " qual:" << m_quality << " eta:" << mu.eta() << " phi:" << mu.phi() << " pt:" << mu.pt() << " ch:" << gemid.chamber() << " l:" << gemid.layer() << " p:" << gemid.roll() << endl;
	}
      }
      if (fidMu) ++nGEMFiducialMuon;
    }

    if (m_nhits > 0 or m_nbounds > 0) {
      m_pt = mu.pt();
      m_eta = mu.eta();
      m_phi = mu.phi();
      t_muon->Fill();
    }
  } // Muon Loop


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
    
    t_hit->Fill();

    //    cout << "   --- eta: " << globalPosition.eta() << " phi: " << globalPosition.phi()  << " xyz: " << b_x << " " << b_y << " " << b_z << " ch:" << detId.chamber() << " l:" << detId.layer() << " p:" << detId.roll() << endl;
  } // GEM Loop
  //  cout << endl << endl;

  nMuonTotal += b_nMuons;
  nGEMTotal += b_nGEMHits;

  t_event->Fill();
}

void AodSliceTestAnalysis::beginJob(){}
void AodSliceTestAnalysis::endJob(){}

void AodSliceTestAnalysis::beginRun(Run const& run, EventSetup const&){
}
void AodSliceTestAnalysis::endRun(Run const&, EventSetup const&){}

//define this as a plug-in
DEFINE_FWK_MODULE(AodSliceTestAnalysis);
