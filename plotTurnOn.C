#include "plotTurnOn.h"


void plotTurnOn()
{
  const std::string recFilePath = "/eos/cms/store/group/phys_heavyions/ikucher/MB2_run327237/HIMinimumBias2/Forest_HIMinimumBias2_run327237_merged.root";
  const std::vector<std::pair<std::string, std::string> > hltFilePath = {
      {"Online", "/eos/cms/store/group/phys_heavyions/anstahll/HLT2021/TREE/OpenHLTTree_HIMinBias_103X_dataRun2_HLT_NoZS_20200521.root"},
      {"Miscalibrated", "/eos/cms/store/group/phys_heavyions/anstahll/HLT2021/TREE/OpenHLTTree_HIMinBias_103X_dataRun2_HLT_ForHITestsV3_NoZS_20200521.root"}
  };
  const std::map<std::string, int> COLOR = { {"Online", kBlack}, {"Miscalibrated", kRed} };

  // initialize efficiency objects
  TH1::AddDirectory(kFALSE);
  std::vector<TEfficiency> effVec;
  std::map<std::string, std::map<std::string, std::map<std::string, std::string> > > effMapTitle;
  std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, size_t> > > > effMapIdx;
  for (const auto& p : TRIGLIST) {
    for (const auto& t : p.second) {
      for (const auto& b : BINMAP.at(p.first)) {
        const std::string name = "eff"+p.first+"_"+t+"_"+b.first;
        effMapTitle[p.first][t][b.first] = name;
        for (const auto& c : hltFilePath) {
          effVec.emplace_back(TEfficiency("", "", b.second.size()-1, b.second.data()));
          effVec.back().SetName((name+"_"+c.first).c_str());
          effMapIdx[p.first][t][b.first][c.first] = effVec.size()-1;
        }
      }
    }
  }

  std::set<std::string> doParticle;
  for (const auto& t : TRIGLIST) {
    if      (t.first.rfind("Electron")!=std::string::npos) { doParticle.insert("electron"); }
    else if (t.first.rfind("Photon"  )!=std::string::npos) { doParticle.insert("photon");   }
    else if (t.first.rfind("Muon"    )!=std::string::npos) { doParticle.insert("muon");     }
  }

  // prepare multi-threading
  const int nCores = doParticle.size()*hltFilePath.size();
  ROOT::EnableImplicitMT();
  ROOT::TProcessExecutor mpe(nCores);

  TH1::AddDirectory(kFALSE);
  auto fillEfficiency = [=](int idx)
  {
    auto c_start = std::clock();
    const auto& n = doParticle.size();
    const auto& hltFileP = hltFilePath[idx/n];
    const auto& pTag = *std::next(doParticle.begin(), idx%n);
    auto parName = pTag; parName[0] = toupper(parName[0]);

    std::set<int> idxS;
    std::vector<TEfficiency> effV(effVec);

    // get offline
    RecoReader recoInfo(recFilePath, false);

    // get online information
    TriggerReader triggerInfo(hltFileP.second, true);

    // add trigger paths to trigger reader
    for (const auto& t : TRIGLIST) {
      for (const auto& path : t.second) {
        triggerInfo.addTrigger(path);
      }
    }

    // initialize offline information
    recoInfo.initBranches(pTag);

    // fill efficiencies
    const auto nEntries = recoInfo.getEntries();
    for (const auto& iEntry : ROOT::TSeqUL(nEntries)) {
      if ((iEntry%100000)==0) { std::cout << "[INFO] Core " << idx << ":  Processing event " << iEntry << " / " << nEntries << std::endl; }
      recoInfo.setEntry(iEntry, false, true);
      if (!triggerInfo.setEntry(recoInfo.getEventNumber(), false, true)) continue;
      // check that event pass event selection
      if (!recoInfo.passEventSelection()) continue;
      //const auto cent = recoInfo.getCentrality();
      // loop over particle type
      const auto particles = recoInfo.getParticles(pTag);
      // loop over single particles
      for (size_t iPar1=0; iPar1<particles.size(); iPar1++) {
        const auto& particle1 = particles[iPar1];
        // extract variables
        const auto var = std::map<std::string, double>({ {"Pt",  particle1.first.Pt()}, {"Eta", particle1.first.Eta()} });
        // loop over single particle triggers
        auto parTag = "Single"+parName;
        if (TRIGLIST.find(parTag)!=TRIGLIST.end()) {
          for (const auto& path : TRIGLIST.at(parTag)) {
            // check if trigger matched
            const auto isMatched = triggerInfo.isTriggerMatched(particle1.first, path);
            // fill efficiency
            for (const auto& v : BINMAP.at(parTag)) {
              const auto& index = effMapIdx.at(parTag).at(path).at(v.first).at(hltFileP.first);
              effV[index].Fill(isMatched, var.at(v.first));
              idxS.insert(index);
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
              // check if trigger matched
              const bool isMatched1 = triggerInfo.isTriggerMatched(particle1.first, path);
              const bool isMatched2 = triggerInfo.isTriggerMatched(particle2.first, path);
              const bool isMatched = isMatched1 && isMatched2;
              // fill efficiency
              for (const auto& v : BINMAP.at(parTag)) {
                const auto& index = effMapIdx.at(parTag).at(path).at(v.first).at(hltFileP.first);
                effV[index].Fill(isMatched, varDP.at(v.first));
                idxS.insert(index);
              }
            }
          }
        }
      }
    }
    return std::make_pair(idxS, effV);
  };

  const auto& res = mpe.Map(fillEfficiency, ROOT::TSeqI(nCores));

  for (const auto& r : res) {
    for (const auto& idx : r.first) {
      effVec[idx] = r.second[idx];
    }
  }

  // set plot style
  setTDRStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  // plot efficiencies
  for (const auto& p : effMapIdx) {
    for (const auto& t : p.second) {
      for (const auto& v : t.second) {
        const auto& par = p.first;
        const auto& path = t.first;
        const auto& var = v.first;
        const auto& idxM = v.second;
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
          auto& eff = effVec[idxM.at(c.first)];
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
        if (graph.GetN()==0) continue;
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
        const auto& name = effMapTitle[p.first][t.first][v.first];
        c.SaveAs((plotDir + "/png/" + name + ".png" ).c_str());
        c.SaveAs((plotDir + "/pdf/" + name + ".pdf" ).c_str());
        // Clean up memory
        c.Clear(); c.Close();
      }
    }
  }
};
