// loop over file lists to make histograms
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TVector3.h>
#include <TFile.h>
#include <RAT/DS/Root.hh>
#include <RAT/DS/FitVertex.hh>
#include <RAT/DS/PMTCal.hh>
#include "util.h"

// make a TH1F from a selector
TH1F* make_pdf(const std::vector<std::string> filenames,
               const std::string selector, const std::string cut,
               const std::vector<float> bin_edges,
               const std::string name="pdf") {
  size_t ndim = 1;
  if (ndim != 1) {
    std::cerr << "make_pdf: Inconsistent dimensions for 1D PDF" << std::endl;
    return NULL;
  }

  double* edges = new double[bin_edges.size()];
  for (size_t i=0; i<bin_edges.size(); i++) {
    edges[i] = bin_edges[i];
  }

  TH1F* pdf = new TH1F(name.c_str(), name.c_str(), bin_edges.size()-1, edges);

  TChain* t = make_tchain(filenames);
  std::cout << "Events: " << t->GetEntries() << std::endl;

  std::string full_selector = selector + ">>" + name;

  int nevents = t->Draw(full_selector.c_str(), cut.c_str());

  if (nevents == -1) {
    std::cerr << "make_pdf: Error in TTree::Draw" << std::endl;
    return NULL;
  }

  return pdf;
}


// make a TH2F from a selector
TH2F* make_pdf(const std::vector<std::string> filenames,
               const std::string selector, const std::string cut,
               const std::vector<float> x_bin_edges,
               const std::vector<float> y_bin_edges,
               const std::string name="pdf") {
  size_t ndim = 2;
  if (ndim != 2) {
    std::cerr << "make_pdf: Inconsistent dimensions for 2D PDF" << std::endl;
    return NULL;
  }

  double* xedges = new double[x_bin_edges.size()];
  for (size_t i=0; i<x_bin_edges.size(); i++) {
    xedges[i] = x_bin_edges[i];
  }

  double* yedges = new double[y_bin_edges.size()];
  for (size_t i=0; i<y_bin_edges.size(); i++) {
    yedges[i] = y_bin_edges[i];
  }

  TH2F* pdf = new TH2F(name.c_str(), name.c_str(),
                       x_bin_edges.size()-1, xedges,
                       y_bin_edges.size()-1, yedges);

  TChain* t = make_tchain(filenames);
  std::cout << "Events: " << t->GetEntries() << std::endl;

  std::string full_selector = selector + ">>" + name;

  int nevents = t->Draw(full_selector.c_str(), cut.c_str());

  if (nevents == -1) {
    std::cerr << "make_pdf: Error in TTree::Draw" << std::endl;
    return NULL;
  }

  return pdf;
}


// make a TH3F from a selector
TH3F* make_pdf(const std::vector<std::string> filenames,
               const std::string selector, const std::string cut,
               const std::vector<float> x_bin_edges,
               const std::vector<float> y_bin_edges,
               const std::vector<float> z_bin_edges,
               const std::string name="pdf") {
  size_t ndim = 3;
  if (ndim != 3) {
    std::cerr << "make_pdf: Inconsistent dimensions for 3D PDF" << std::endl;
    return NULL;
  }

  double* xedges = new double[x_bin_edges.size()];
  for (size_t i=0; i<x_bin_edges.size(); i++) {
    xedges[i] = x_bin_edges[i];
  }

  double* yedges = new double[y_bin_edges.size()];
  for (size_t i=0; i<y_bin_edges.size(); i++) {
    yedges[i] = y_bin_edges[i];
  }

  double* zedges = new double[z_bin_edges.size()];
  for (size_t i=0; i<z_bin_edges.size(); i++) {
    zedges[i] = z_bin_edges[i];
  }

  TH3F* pdf = new TH3F(name.c_str(), name.c_str(),
                       x_bin_edges.size()-1, xedges,
                       y_bin_edges.size()-1, yedges,
                       z_bin_edges.size()-1, zedges);

  TChain* t = make_tchain(filenames);
  std::cout << "Events: " << t->GetEntries() << std::endl;

  std::string full_selector = selector + ">>" + name;

  int nevents = t->Draw(full_selector.c_str(), cut.c_str());

  if (nevents == -1) {
    std::cerr << "make_pdf: Error in TTree::Draw" << std::endl;
    return NULL;
  }

  return pdf;
}


// make a pdf in fit radius and energy
// can't use selectors since fit results are buried too deep
TH2F* make_fit_e_r_pdf(const std::vector<std::string> filenames,
                       const std::vector<float> x_bin_edges,
                       const std::vector<float> y_bin_edges,
                       const std::string name="pdf") {
  double* xedges = new double[x_bin_edges.size()];
  for (size_t i=0; i<x_bin_edges.size(); i++) {
    xedges[i] = x_bin_edges[i];
  }

  double* yedges = new double[y_bin_edges.size()];
  for (size_t i=0; i<y_bin_edges.size(); i++) {
    yedges[i] = y_bin_edges[i];
  }

  TH2F* pdf = new TH2F(name.c_str(), name.c_str(),
                     x_bin_edges.size()-1, xedges,
                     y_bin_edges.size()-1, yedges);

  TChain* t = make_tchain(filenames);
  std::cout << "Events: " << t->GetEntries() << std::endl;

  RAT::DS::Root* ds = new RAT::DS::Root;

  t->SetBranchStatus("*", 0);
  t->SetBranchStatus("ev", 1);
  t->SetBranchAddress("ds", &ds);

  for (int i=0; i<t->GetEntries(); i++) {
    if (i % 1000 == 0) {
        std::cout << i << "/" << t->GetEntries() << " events" << std::endl;
    }

    t->GetEntry(i);

    if (ds->GetEVCount() == 0) {
      continue;
    }

    bool fit_valid = false;
    try {
      fit_valid = ds->GetEV(0)->GetFitResult("scintFitter").GetValid();
    }
    catch(...) {
      continue;
    }

    if (!fit_valid) {
      continue;
    }

    RAT::DS::FitVertex vertex = ds->GetEV(0)->GetFitResult("scintFitter").GetVertex(0);

    pdf->Fill(vertex.GetPosition().Mag(), vertex.GetEnergy());
  }

  return pdf;
}


void print_usage() {
  std::cout << "PDFify v0.1 -- "
            << "Make PDFs from RAT ROOT files" << std::endl
            << "Usage: pdfify <mode (ratroot|ratlist)> "
            << "(file1.root file2.root ...|filelist.txt) "
            << "outfile.root" << std::endl;
}


int main(int argc, char* argv[]) {
  if (argc < 4) {
    print_usage();
    return 1;
  }

  std::vector<std::string> filenames;
  for (int i=2; i<argc-1; i++) {
    filenames.push_back(argv[i]);
  }

  std::string mode = std::string(argv[1]);
  std::string outfile = std::string(argv[argc-1]);

  if (mode == "ratlist") {
    filenames = read_file_list(filenames);
  }
  else if (mode != "ratroot") {
    print_usage();
    return 2;
  }  

  std::cout << "Processing " << filenames.size() << " files..." << std::endl;
  std::cout << "Mode: " << mode << std::endl;
  std::cout << "Output: " << outfile << std::endl;

  std::vector<float> rbins = {0,    1000, 2000, 3000, 4000, 5000, 5100, 5200,
                              5300, 5400, 5500, 5600, 5700, 5800, 5900, 6005};

  std::vector<float> ebins(501);
  for (unsigned i=0; i<501; i++) {
    ebins[i] = 0.01 * i;  // 10 keV
  }
    
  TFile f(outfile.c_str(), "recreate");
  f.cd();
  TH2F* pdf = make_fit_e_r_pdf(filenames, rbins, ebins);
  std::cout << "Total hits: " << pdf->Integral() << std::endl;
  pdf->Write();
  f.Close();

  return 0;
}

