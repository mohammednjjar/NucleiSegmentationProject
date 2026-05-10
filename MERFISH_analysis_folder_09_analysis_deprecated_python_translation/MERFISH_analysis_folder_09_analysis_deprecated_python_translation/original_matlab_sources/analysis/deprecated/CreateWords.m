function [words, parameters] = CreateWords(imageData,varargin)
% Full source retrieved from GitHub raw view in this conversation. The GitHub file is minified into long lines.
% Purpose: constructs MERFISH word structures from localization lists using either common-centroid or
% per-localization nearest-neighbor word construction, then fills codeword/metadata fields.
% The Python translation is in python_translation/analysis/deprecated/CreateWords.py.
