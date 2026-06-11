% startup.m -- project initialisation (run from project root each session)
%
% Adds project code and the read-only external libraries to the path.
% Nothing under external/ may be modified.

proj = fileparts(mfilename('fullpath'));

addpath(fullfile(proj, 'src'));
addpath(fullfile(proj, 'experiments'));
addpath(fullfile(proj, 'tests'));
addpath(fullfile(proj, 'postprocess'));

% Read-only references, on path for interactive comparison runs only
addpath(fullfile(proj, 'external', 'Moll-codes'));
addpath(genpath(fullfile(proj, 'external', 'SparseEcon', 'lib')));

format compact
