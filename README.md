# Load-dependent propagation of geometry-freezing errors in monopile foundations

MATLAB source code and numerical data for the mechanical, benchmark, field-comparison and auxiliary prediction analyses.

Start MATLAB in a writable copy of this directory and execute `run_supplement('checks')`. See [RUN_GUIDE.md](RUN_GUIDE.md) for the analysis commands and output locations. [DATA_SOURCES.csv](DATA_SOURCES.csv) distinguishes digitized published curves, interpreted soil inputs, numerical reference results and generated auxiliary data.

The primary calculations use MATLAB R2024b. The supplied full-precision results are reference outputs of the recorded numerical analyses. Published field curves are digitized or author-processed observations. Raw instrument records are not included. The PISA design benchmarks are numerical results from the related model family. The auxiliary prediction problem is a controlled numerical study.

The external OpenPile package is used only for optional implementation comparisons. Its source distribution, binaries and third-party full publications are not bundled. Software and source attribution are provided in [DEPENDENCIES.md](DEPENDENCIES.md).
