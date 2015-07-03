% Examples of simulations ran with the use of CSHS.m and the therein called
% functions.
%
% Required functions:
%   cget_rho        Obtains the connectivity matrix for the KCs according
%                   to the rules set in the main text.
%   cget_clusters   Randomly generates clusters.
%   cget_obs        Randomly generates the observation matrix.
%   UKF             Unscented Kalman Filter algorithm. Written by Sebastian
%                   Bitzer; can be found at: 
%                   https://github.com/sbitzer/BayesFilters
%                   Additionally, the function isonemptyfield is needed,
%                   which can be found at:
%                   https://github.com/sbitzer/Matlab-helpers

%% Example 1: Simple simulation
%==========================================================================
par.N = 50;
par.M = 30;
par.sqlh = 10;
par.tend = 150;
par.num_seqs = 1;
par.Full = 0;
par.Plot = 1;
par.PlotPN = 1;

[data, par, Kal] = CSHS(par);


%% Example 2: Full response
%==========================================================================
% The simulation will be ran twice, changing the precision values to see
% their effects. In the second set of precision values, the observation
% (inverse) precision is lower, thus the switch between phases is not done
% properly.
% The results can be seen side by side for comparison.

clear par;
par.Full = 1;
par.Qn = 0.1;
par.Rn = 0.001;

[data, par, Kal] = CSHS(par);

% Duplicate plot window
%--------------------------------------------------------------------------
temp = findobj('name','Plots of Simulation Results');
h1=temp(1);
h2=figure;
mystr = sprintf('Plots with Qn = %f, Rn = %f',par.Qn,par.Rn);
set(h2,'name',mystr);
objects=allchild(h1);
copyobj(get(h1,'children'),h2);
set(h2,'position',get(h1,'position'));

% Run again
%--------------------------------------------------------------------------
par.Qn = 0.1;
par.Rn = 0.1;

[data, par, Kal] = CSHS(par);


%% Example 3: Overnumerousness
%==========================================================================
% With a high-enough number of embedded sequences in the connectivity
% matrix, depending on the population sizes, the shape of the observed
% sequences might change. Even in this cases, the inversion can still work.

clear par;
par.N = 50;
par.num_seqs = 100;

[data, par, Kal] = CSHS(par);
