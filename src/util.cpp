#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <streambuf>
#include <algorithm>
#include <assert.h>
#include <jsoncpp/json/value.h>
#include <jsoncpp/json/reader.h>
#include <TFile.h>
#include <TRandom.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2F.h>

#include "util.h"

#ifndef VERBOSE
#define VERBOSE false
#endif

std::vector<std::string> read_file_list(std::vector<std::string> l) {
  std::vector<std::string> filelist;

  for (size_t i=0; i<l.size(); i++) {
    std::ifstream infile(l[i].c_str());
    std::string name;
    while (infile >> name)
    {
      filelist.push_back(name);
    }
  }

  return filelist;
}


TChain* make_tchain(std::vector<std::string> filenames, std::string name) {
  TChain* t = new TChain(name.c_str());

  for (size_t i=0; i<filenames.size(); i++) {
    std::string filename = filenames[i];
    std::cout << filename << std::endl;
    t->Add(filename.c_str());
  }
  std::cout << "created tchain" << std::endl;
  return t;
}


double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    for (size_t i=0;i<v.size(); i++) {
      std::cout << v[i] << " ";
    }
    std::cout << std::endl;
    size_t n = v.size() / 2;
    std::nth_element(v.begin(), v.begin()+n, v.end());
    return v[n];
}


TNtuple* FakeDataGenerator::make_dataset(bool poisson, std::map<std::string, double>* _norms) {
  if (VERBOSE) {
    std::cout << "FakeDataGenerator::make_dataset: Generating dataset..." << std::endl;
  }
  TNtuple* nt = new TNtuple("tev", "events", "r:e");
  double r = 0;
  double e = 0;
  for (auto it=this->pdfs.begin(); it!=this->pdfs.end(); it++) {
    int nexpected = ((_norms && _norms->find(it->first) != _norms->end()) ? (*_norms)[it->first] : this->default_norms[it->first]);

    if (it->second->IsA() == TH2F::Class()) {
      TH2F* ht = dynamic_cast<TH2F*>(it->second);
      //int x1 = it->second->GetXaxis()->FindBin(r_range.min);
      //int x2 = it->second->GetXaxis()->FindBin(r_range.max);
      //int y1 = it->second->GetYaxis()->FindBin(e_range.min);
      //int y2 = it->second->GetYaxis()->FindBin(e_range.max);
      //double integral = ht->Integral(x1, x2, y1, y2);
      //nexpected *= integral;
      int nobserved = nexpected;
      if (poisson && nexpected > 0) {
        nobserved = gRandom->Poisson(nexpected);
      }
      for (int i=0; i<nobserved; i++) {
        do {
          ht->GetRandom2(r, e);
        } while(e > e_range.max || e < e_range.min || r > r_range.max || r < r_range.min);
        nt->Fill(r, e);
      }
      if (VERBOSE) {
        std::cout << "FakeDataGenerator::make_dataset: " << it->first << ": " << nobserved << " events" << std::endl;
      }
    }
    else if (it->second->IsA() == TH1D::Class()) {
      TH1D* ht = dynamic_cast<TH1D*>(it->second);
      //int x1 = it->second->GetXaxis()->FindBin(e_range.min);
      //int x2 = it->second->GetXaxis()->FindBin(e_range.max);
      //double integral = ht->Integral(x1, x2);
      //nexpected *= integral;
      int nobserved = nexpected;
      if (poisson && nexpected > 0) {
        nobserved = gRandom->Poisson(nexpected);
      }
      for (int i=0; i<nobserved; i++) {
        do {
          e = ht->GetRandom();
        } while(e > e_range.max || e < e_range.min);
        nt->Fill(r, e);
      }
      if (VERBOSE) {
        std::cout << "FakeDataGenerator::make_dataset: " << it->first << ": " << nobserved << " events" << std::endl;
      }
    }
    else {
      std::cerr << "FakeDataGenerator::make_dataset: Unknown histogram class " << it->second->ClassName() << std::endl;
      assert(false);
    }
  }

  return nt;
}


FitConfig::FitConfig(std::string const filename) {
  Json::Value root;
  Json::Reader reader;

  std::ifstream t(filename);
  std::string data((std::istreambuf_iterator<char>(t)),
                   std::istreambuf_iterator<char>());

  bool parse_ok = reader.parse(data, root);
  if (!parse_ok) {
    std::cout  << "FitConfig: JSON parse error:" << std::endl
               << reader.getFormattedErrorMessages();
    assert(false);
  }

  // experiment parameters
  const Json::Value experiment_params = root["experiment"];
  assert(experiment_params.isMember("live_time"));
  this->live_time = experiment_params["live_time"].asFloat();
  if (experiment_params.isMember("efficiency")) {
    this->efficiency = experiment_params["efficiency"].asFloat();
  }
  else {
    this->efficiency = 1;
  }
  assert(experiment_params.isMember("confidence"));
  this->confidence = experiment_params["confidence"].asFloat();

  // fit parameters
  const Json::Value fit_params = root["fit"];
  this->mode = static_cast<FitMode>(fit_params["mode"].asInt());
  this->e_range.min = fit_params["energy_range"][0].asFloat();
  this->e_range.max = fit_params["energy_range"][1].asFloat();
  this->r_range.min = fit_params["radius_range"][0].asFloat();
  this->r_range.max = fit_params["radius_range"][1].asFloat();
  assert(fit_params.isMember("mc_trials"));
  this->mc_trials = fit_params["mc_trials"].asInt();
  if (fit_params.isMember("fake_experiments")) {
    this->fake_experiments = fit_params["fake_experiments"].asInt();
  }
  else {
    this->fake_experiments = this->mc_trials;
  }
  if (fit_params.isMember("output_file")) {
    this->output_file = fit_params["output_file"].asString();
  }
  else {
    this->output_file = "fit_spectrum";
  }
  assert(fit_params.isMember("signal_name"));
  this->signal_name = fit_params["signal_name"].asString();

  std::vector<std::string> fit_signal_names;
  for (Json::Value::iterator it=fit_params["signals"].begin(); it!=fit_params["signals"].end(); it++) {
    fit_signal_names.push_back((*it).asString());
  }

  // signal parameters
  const Json::Value signal_names = root["signals"];
  for (auto it=signal_names.begin(); it!=signal_names.end(); it++) {
    if (std::find(fit_signal_names.begin(), fit_signal_names.end(), it.key().asString()) == fit_signal_names.end()) {
      continue;
    }
    
    const Json::Value signal_params = root["signals"][it.key().asString()];

    Signal s;
    s.name = it.key().asString();
    s.title = signal_params.get("title", s.name).asString();
    s.fixed = signal_params.get("fixed", false).asBool();
    s.constraint = signal_params.get("constraint", 0).asFloat();

    assert(signal_params.isMember("rate"));
    s.rate = signal_params["rate"].asFloat() * this->live_time;  // fixme is event expectation value not rate

    assert(signal_params.isMember("filename"));
    std::string filename = signal_params["filename"].asString();

    TH2F* h2d = load_histogram(filename, "pdf");
    h2d->Sumw2();
    h2d->Scale(1.0/h2d->Integral());

    if (this->mode == FitMode::ENERGY) {
      s.histogram = dynamic_cast<TH1*>(project1d(h2d, &r_range));
      int x1 = s.histogram->GetXaxis()->FindBin(this->e_range.min);
      int x2 = s.histogram->GetXaxis()->FindBin(this->e_range.max);
      double integral = s.histogram->Integral(x1, x2);
      std::cout << s.name << ": "  << integral << " " << s.histogram->Integral() << std::endl;
      s.rate *= integral / h2d->Integral();
      s.histogram->Rebin(4);
    }
    else if (this->mode == FitMode::ENERGY_RADIUS) {
      s.histogram = dynamic_cast<TH1*>(h2d);
      int x1 = h2d->GetXaxis()->FindBin(r_range.min);
      int x2 = h2d->GetXaxis()->FindBin(r_range.max);
      int y1 = h2d->GetYaxis()->FindBin(e_range.min);
      int y2 = h2d->GetYaxis()->FindBin(e_range.max);
      double integral = h2d->Integral(x1, x2, y1, y2);
      s.rate *= integral / h2d->Integral();
      dynamic_cast<TH2F*>(s.histogram)->RebinY(4);
    }
    else {
      std::cerr << "Unknown fit mode " << static_cast<int>(this->mode) << std::endl;
      assert(false);
    }

    s.histogram->Scale(1.0/s.histogram->Integral());
    this->signals.push_back(s);
  }
}


TH2F* FitConfig::load_histogram(std::string const filename, std::string const objname) {
  TFile f(filename.c_str());
  assert(!f.IsZombie());
  TH2F* h = dynamic_cast<TH2F*>(f.Get(objname.c_str()));
  h->SetDirectory(0);
  h->Sumw2();
  return h;
}


TH1D* FitConfig::project1d(TH2F* const hist2d, Range<float>* const r_range) {
  int first_bin = 0;
  int last_bin = -1;
  if (r_range != nullptr) {
    first_bin = hist2d->GetXaxis()->FindBin(r_range->min);
    last_bin = hist2d->GetXaxis()->FindBin(r_range->max);
  }
  TH1D* h = hist2d->ProjectionY("pdf_proj", first_bin, last_bin);
  h->SetDirectory(0);
  return h;  
}


void FitConfig::print() const {
  std::cout << "Fit:" << std::endl
            << "  Mode: " << static_cast<int>(this->mode) << std::endl
            << "  MC trials: " << this->mc_trials << std::endl
            << "  Fake experiments: " << this->fake_experiments << std::endl
            << "  Signal name: " << this->signal_name << std::endl
            << "  Energy: (" << this->e_range.min << ", " << this->e_range.max << ") MeV" << std::endl
            << "  Radius: (" << this->r_range.min << ", " << this->r_range.max << ") mm" << std::endl
            << "  Output plot:" << this->output_file << std::endl
            << "Experiment:" << std::endl
            << "  Live time: " << this->live_time << " y" << std::endl
            << "  Confidence level: " << this->confidence << std::endl
            << "Signals:" << std::endl;

  for (auto it=this->signals.begin(); it!=this->signals.end(); it++) {
    std::cout << "  " << it->name << std::endl;
    std::cout << "    Title: \"" << it->title << "\"" << std::endl;
    std::cout << "    Rate: " << it->rate << std::endl;
    std::cout << "    Constraint: ";
    if (it->constraint != 0) {
      std::cout << it->constraint << std::endl;
    }
    else {
      std::cout << "none" << std::endl;
    }
    std::cout << "    Fixed: " << (it->fixed ? "yes" : "no") << std::endl;
  }
}

