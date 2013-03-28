#include <assert.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <streambuf>
#include <algorithm>
#include <jsoncpp/json/value.h>
#include <jsoncpp/json/reader.h>
#include <TFile.h>
#include <TChain.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2F.h>
#include "util.h"

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

double median(std::vector<double> &v) {
    size_t n = v.size() / 2;
    std::nth_element(v.begin(), v.begin()+n, v.end());
    return v[n];
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
               <<reader.getFormattedErrorMessages();
    assert(false);
  }

  // experiment parameters
  const Json::Value experiment_params = root["experiment"];
  assert(experiment_params.isMember("live_time"));
  this->live_time = experiment_params["live_time"].asFloat();
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

    if (this->mode == FitMode::ENERGY) {
      s.histogram = dynamic_cast<TH1*>(project1d(h2d, &r_range));
      s.rate *= (s.histogram->Integral() / h2d->Integral());
    }
    else if (this->mode == FitMode::ENERGY_RADIUS) {
      s.histogram = dynamic_cast<TH1*>(h2d);
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
            << "  Signal name: " << this->signal_name << std::endl
            << "  Energy: (" << this->e_range.min << ", " << this->e_range.max << ") MeV" << std::endl
            << "  Radius: (" << this->r_range.min << ", " << this->r_range.max << ") mm" << std::endl
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

