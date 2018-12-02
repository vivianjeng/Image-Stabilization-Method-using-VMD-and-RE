% Steerable pyramid filters.  Transform described  in:
%
% @INPROCEEDINGS{Simoncelli95b,
%	TITLE = "The Steerable Pyramid: A Flexible Architecture for
%		 Multi-Scale Derivative Computation",
%	AUTHOR = "E P Simoncelli and W T Freeman",
%	BOOKTITLE = "Second Int'l Conf on Image Processing",
%	ADDRESS = "Washington, DC", MONTH = "October", YEAR = 1995 }
%
% Filter kernel design described in:
%
%@INPROCEEDINGS{Karasaridis96,
%	TITLE = "A Filter Design Technique for 
%		Steerable Pyramid Image Transforms",
%	AUTHOR = "A Karasaridis and E P Simoncelli",
%	BOOKTITLE = "ICASSP",	ADDRESS = "Atlanta, GA",
%	MONTH = "May",	YEAR = 1996 }

% Eero Simoncelli, 6/96.

function [lo0filt,hi0filt,lofilt,bfilts,steermtx,harmonics] = G3();

% JEBYRNE: load rather than copy into m file
if (exist('G3.mat','file') == 2)
  load('G3.mat'); 
else
  error('''G3.mat'' not found.  Run create_numerical_hilbert.m');
end

