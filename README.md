Signal Extraction Fit
=====================
A GPU-accelerated fit program which calculates fully frequentist confidence
intervals or limits via the profile likelihood construction.

Building
--------
`sigex` requires the following libraries:

* CUDA (and the nvcc compiler)
* JsonCpp
* ROOT
* RAT (to make PDFs from MC files)

It also requires a CUDA-enabled Nvidia GPU.

To build, run `make` and specify a `CUDA_ROOT` environment variable. E.g.,

    $ make CUDA_ROOT=/opt/cuda-5.0

Usage
-----
1. Create PDFs: PDFs are ROOT TH2Fs with event energy and radius dimensions.
   `bin/pdf` converts a set of RAT ROOT files to an appropriate PDF.

2. Configure fit: Set up the fit parameters are signal PDFs using a JSON-format
   configuration file.

3. Run fit: `$ ./bin/fit config.json`

Configuration Files
-------------------
The fit is controlled via a JSON-format configuration file.

Example:

    {
      "fit": {
        "mode": 0,
        "mc_trials": 100,
        "signal_name": "0vbb",
        "energy_range": [0.0, 5.0],
        "radius_range": [0.0, 5000.0],
        "signals": ["0vbb", "2vbb", "b8"]
      },
      "experiment": {
        "live_time": 2.0,
        "confidence": 0.9
       },
      "signals": {
        "0vbb": {
          "title": "0#nu#beta#beta",
          "filename": "pdfs/zeronu.root",
          "rate": 405,
          "fixed": false
        },
        "2vbb": {
          "title": "2#nu#beta#beta",
          "filename": "pdfs/twonu.root",
          "rate": 225730,
          "constraint": 0.1,
          "fixed": false
        },
        "b8": {
          "title": "^{8}B #nu ES",
          "filename": "pdfs/b8.root",
          "rate": 1045,
          "constraint": 0.05,
          "fixed": false
        },
      }
    }

* `fixed` controls whether a normalization floats in the fit, and defaults to
  `false`
* `constraint` is a fractional uncertainty on the normalization, and
  defaults to 0 (unconstrained).
* `rate` is a number of events per year

Not all signals specified in `signals` need be used in the fit -- only those
listed in `fit.signals` are considered.

