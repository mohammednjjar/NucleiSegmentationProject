# MERFISH_analysis Folder 06: root deprecated/ Python translation

This package translates the root `deprecated/` folder from the ZhuangLab MERFISH_analysis repository.

Coverage: 17 MATLAB files mapped to 17 Python files. The translated modules contain executable Python logic for the old utilities, probe-construction helpers, and report calculations. External tools such as BLAST/OligoArray are called through Python subprocess wrappers and can be run with execution disabled for command inspection.

| MATLAB file | Python file | Purpose/use |
|---|---|---|
| `misc/PageBreak.m` | `deprecated/misc/PageBreak.py` | Prints/returns a horizontal separator line for console reports. |
| `misc/PlotCorr2.m` | `deprecated/misc/PlotCorr2.py` | Computes Pearson correlation on linear and log10 values and creates a correlation plot. |
| `misc/StringFind.m` | `deprecated/misc/StringFind.py` | Searches cell-string lists for exact or substring matches. |
| `probe_construction/AddCntrlSeqs.m` | `deprecated/probe_construction/AddCntrlSeqs.py` | Adds random no-target controls and blank controls into picked probe-gene structures. |
| `probe_construction/AssembleProbes.m` | `deprecated/probe_construction/AssembleProbes.py` | Assembles a MERFISH probe library from ProbeData, primers, readouts, controls, and SECDED codebook. |
| `probe_construction/BatchLaunchOligoArray.m` | `deprecated/probe_construction/BatchLaunchOligoArray.py` | Splits gene FASTA entries and launches/records OligoArray commands per gene. |
| `probe_construction/BuildBLASTlib.m` | `deprecated/probe_construction/BuildBLASTlib.py` | Builds a nucleotide BLAST database from a FASTA file. |
| `probe_construction/CompileOligoArrayOutput.m` | `deprecated/probe_construction/CompileOligoArrayOutput.py` | Parses OligoArray output files into a ProbeData structure. |
| `probe_construction/GenProbe.m` | `deprecated/probe_construction/GenProbe.py` | Builds oligo sequences and FASTA output from codebook, probes, readouts, and primers. |
| `probe_construction/OligoArrayCmd.m` | `deprecated/probe_construction/OligoArrayCmd.py` | Constructs the Java OligoArray command string and save path. |
| `probe_construction/WriteFasta.m` | `deprecated/probe_construction/WriteFasta.py` | Writes FASTA records from headers/sequences or MATLAB-style structures. |
| `reports/GenerateBitFlipReport.m` | `deprecated/reports/GenerateBitFlipReport.py` | Estimates per-hybridization bit-flip probabilities from observed codeword counts. |
| `reports/GenerateCompositeImage.m` | `deprecated/reports/GenerateCompositeImage.py` | Generates cell word-overlay composite figures from decoded words and image metadata. |
| `reports/GenerateFPKMReport.m` | `deprecated/reports/GenerateFPKMReport.py` | Compares decoded gene counts against FPKM values and calculates correlation. |
| `reports/GenerateHammingSphereReport.m` | `deprecated/reports/GenerateHammingSphereReport.py` | Counts words in Hamming spheres around exact codewords and computes confidence ratios. |
| `reports/GenerateMoleculeStatsReport.m` | `deprecated/reports/GenerateMoleculeStatsReport.py` | Summarizes per-hybridization molecule intensity/background/height and spatial offsets. |
| `reports/GenerateOnBitHistograms.m` | `deprecated/reports/GenerateOnBitHistograms.py` | Creates histograms of number of on-bits per word, globally and by cell. |

Checks included: Python compile check, unit tests, and source scan.
