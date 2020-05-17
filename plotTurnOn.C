#include "plotTurnOn.h"


void plotTurnOn()
{
  const std::string recFilePath = "/eos/cms/store/group/phys_heavyions/ikucher/MB2_run327237/HIMinimumBias2/Forest_HIMinimumBias2_run327237_merged.root";
  const std::vector<std::pair<std::string, std::string> > hltFilePath = {
      {"Online", "/eos/cms/store/group/phys_heavyions/anstahll/HLT2021/TREE/OpenHLTTree_HIMinBias_103X_dataRun2_HLT_Run327237.root"},
      {"Miscalibrated", "/eos/cms/store/group/phys_heavyions/anstahll/HLT2021/TREE/OpenHLTTree_HIMinBias_103X_dataRun2_HLT_ForHITestsV3_Run327237.root"}
  };
  const std::map<std::string, int> COLOR = { {"Online", kBlack}, {"Miscalibrated", kRed} };

  // get online information
  std::map<std::string, TriggerReader> triggerInfo;
  for (const auto& c : hltFilePath) {
    triggerInfo.emplace(std::piecewise_construct, std::forward_as_tuple(c.first), std::forward_as_tuple(c.second, true));
  }

  // get offline
  RecoReader recoInfo(recFilePath);
  
  // add trigger paths to trigger reader
  std::set<std::string> doParticle;
  for (const auto& t : TRIGLIST) {
    for (const auto& path : t.second) {
      for (auto& c : triggerInfo) {
        c.second.addTrigger(path);
      }
    }
    if      (t.first.rfind("Electron")!=std::string::npos) { doParticle.insert("electron"); }
    else if (t.first.rfind("Photon"  )!=std::string::npos) { doParticle.insert("photon");   }
    else if (t.first.rfind("Muon"    )!=std::string::npos) { doParticle.insert("muon");     }
  }

  // initialize offline information
  for (const auto& p : doParticle) { recoInfo.initBranches(p); }

  // initialize efficiency objects
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, TEfficiency> > > > effMap;
  for (const auto& t : TRIGLIST) {
    for (const auto& path : t.second) {
      for (const auto& b : BINMAP.at(t.first)) {
        const std::string name = "eff"+t.first+"_"+path+"_"+b.first;
        for (const auto& c : hltFilePath) {
          auto& eff = effMap[t.first][path][b.first][c.first];
          eff = TEfficiency("", "", b.second.size()-1, b.second.data());
          eff.SetName((name+"_"+c.first).c_str());
          eff.SetTitle(name.c_str());
        }
      }
    }
  }

  // fill efficiencies
  const auto nEntries = recoInfo.getEntries();
  for (Long64_t iEntry=0; iEntry<nEntries; iEntry++) {
    if ((iEntry%10000)==0) { std::cout << "[INFO] Processing event " << iEntry << " / " << nEntries << std::endl; }
    recoInfo.setEntry(iEntry, true, true);
    bool ignoreEvent = false;
    for (auto& c : triggerInfo) {
      if (!c.second.setEntry(recoInfo.getEventNumber(), true, true)) { ignoreEvent = true; break; }
    }
    if (ignoreEvent) continue;
    // check that event pass event selection
    if (!recoInfo.passEventSelection()) continue;
    const auto cent = recoInfo.getCentrality();
    // loop over particle type
    for (const auto& p : doParticle) {
      std::string parName = p; parName[0] = toupper(parName[0]);
      const auto particles = recoInfo.getParticles(p);
      // loop over single particles
      for (size_t iPar1=0; iPar1<particles.size(); iPar1++) {
        const auto& particle1 = particles[iPar1];
        // extract variables
        const auto var = std::map<std::string, double>({ {"Pt",  particle1.first.Pt()}, {"Eta", particle1.first.Eta()} });
        // loop over single particle triggers
        auto parTag = "Single"+parName;
        if (TRIGLIST.find(parTag)!=TRIGLIST.end()) {
          for (const auto& path : TRIGLIST.at(parTag)) {
            for (auto& c : triggerInfo) {
              // check if trigger matched
              const auto isMatched = c.second.isTriggerMatched(particle1.first, path);
              // fill efficiency
              for (const auto& v : BINMAP.at(parTag)) {
                effMap.at(parTag).at(path).at(v.first).at(c.first).Fill(isMatched, var.at(v.first));
              }
            }
          }
        }
        // loop over double particles
        parTag = "Double"+parName;
        if (TRIGLIST.find(parTag)!=TRIGLIST.end()) {
          for (size_t iPar2=iPar1+1; iPar2<particles.size(); iPar2++) {
            const auto& particle2 = particles[iPar2];
            const auto p4 = particle1.first + particle2.first;
            // extract variables
            const auto varDP = std::map<std::string, double>({ {"Pt",  p4.Pt()}, {"Rapidity", p4.Rapidity()} });
            // loop over ddouble particle triggers
            for (const auto& path : TRIGLIST.at(parTag)) {
              for (auto& c : triggerInfo) {
                // check if trigger matched
                const bool isMatched1 = c.second.isTriggerMatched(particle1.first, path);
                const bool isMatched2 = c.second.isTriggerMatched(particle2.first, path);
                const bool isMatched = isMatched1 && isMatched2;
                // fill efficiency
                for (const auto& v : BINMAP.at(parTag)) {
                  effMap.at(parTag).at(path).at(v.first).at(c.first).Fill(isMatched, varDP.at(v.first));
                }
              }
            }
          }
        }
      }
    }
  }
  

  // set plot style
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // plot efficiencies
  for (const auto& p : effMap) {
    for (const auto& t : p.second) {
      for (const auto& v : t.second) {
        const auto& par = p.first;
        const auto& path = t.first;
        const auto& var = v.first;
        auto& effM = effMap.at(par).at(path).at(var);
        // create Canvas
        TCanvas c("c", "c", 1000, 1000); c.cd();
        // create the text info
        TLatex tex; tex.SetNDC(); tex.SetTextSize(0.028); float dy = 0;
        std::vector< std::string > textToPrint;
        textToPrint.push_back(path);
        if (par.rfind("Electron")!=std::string::npos) {
          textToPrint.push_back("p^{e}_{T} > 20 GeV/c");
          textToPrint.push_back("|#eta^{e}| < 2.1");
        }
        else if (par.rfind("Photon")!=std::string::npos) {
          textToPrint.push_back("p^{#gamma}_{T} > 40 GeV/c");
          textToPrint.push_back("|#eta^{e}| < 2.4");
        }
        else if (par.rfind("Muon")!=std::string::npos) {
          textToPrint.push_back("p^{#mu}_{T} > 1.5 GeV/c");
          textToPrint.push_back("|#eta^{#mu}| < 2.4");
        }
        // format efficiency
        std::vector<std::pair<std::string, TGraphAsymmErrors> > graphM;
        for (const auto& c : hltFilePath) {
          auto& eff = effM.at(c.first);
          eff.Draw(); gPad->Update();
          graphM.push_back({c.first, *eff.GetPaintedGraph()});
          auto& graph = graphM.back().second;
          graph.SetName(eff.GetName());
          formatEff(graph, par, var);
          graph.SetMarkerColor(COLOR.at(c.first));
          graph.SetLineColor(COLOR.at(c.first));
          graph.GetXaxis()->SetLimits(eff.GetTotalHistogram()->GetXaxis()->GetXmin(), eff.GetTotalHistogram()->GetXaxis()->GetXmax());
          if (c.first!=hltFilePath.begin()->first) {
            graph.SetMarkerSize(0);
            graph.SetLineWidth(2);
          }
        }
        auto& graph = graphM[0].second;
        // add to legend
        TLegend leg(0.43, 0.74, 0.6, 0.84);
        for (auto& g : graphM) {
          leg.AddEntry(&g.second, g.first.c_str(), "pel")->SetTextSize(0.032);
        }
        // draw efficiency
        graph.Draw("ap");
        for (auto& g : graphM) {
          if (g.first==graphM[0].first) continue;
          g.second.Draw("samep");
        }
        // draw legend
        leg.Draw("same");
        // draw line
        TLine line(graph.GetXaxis()->GetXmin(), 1.0,  graph.GetXaxis()->GetXmax(), 1.0);
        line.SetLineStyle(2);
        line.Draw("same");
        // Update
        c.Modified(); c.Update();
        // Draw the text
        tex.SetTextSize(0.058); tex.SetTextFont(61); tex.DrawLatex(0.78, 0.84, "CMS"); tex.SetTextFont(62);
        tex.SetTextSize(0.044); tex.SetTextFont(52); tex.DrawLatex(0.69, 0.79, "Preliminary"); tex.SetTextFont(62);
        for (size_t i=0; i<textToPrint.size(); i++) {
          tex.SetTextSize(0.030); tex.DrawLatex(0.20, 0.86-i*0.05, textToPrint[i].c_str());
        }
        c.Modified(); c.Update();
        // set the CMS style
        CMS_lumi(&c, 33, "PbPb run 327237", "#sqrt{s_{NN}} = 5.02 TeV", false, 0.60, false);
        c.Modified(); c.Update();
        // Create Output Directory
        const std::string plotDir = "Plot/" + par;
        makeDir(plotDir + "/png/");
        makeDir(plotDir + "/pdf/");
        // Save Canvas
        const std::string name = effM.begin()->second.GetTitle();
        c.SaveAs((plotDir + "/png/" + name + ".png" ).c_str());
        c.SaveAs((plotDir + "/pdf/" + name + ".pdf" ).c_str());
        // Clean up memory
        c.Clear(); c.Close();
      }
    }
  }
};
