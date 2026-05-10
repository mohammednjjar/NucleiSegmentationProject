function [words, totalImageData, totalFiducialData, parameters] = AnalyzeMERFISH(dataPath, varargin)
% Full source retrieved from GitHub raw view in this conversation. The GitHub file is minified into long lines.
% Purpose: analyzes raw conventional/MERFISH image molecule-list files in a directory, aligns fiducials,
% transforms image localizations, creates words, decodes them with a codebook, and generates reports.
% The Python translation is in python_translation/analysis/deprecated/AnalyzeMERFISH.py.
