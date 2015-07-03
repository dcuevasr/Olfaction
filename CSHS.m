% Main simulation of the Cluster Stable Heteroclinic Sequence through the
% Lotka-Volterra equations for the KC space, a linear projection to the PN
% space and the Unscented Kalman Filter (each section in its corresponding 
% cell).
%
% See the file Examples.m for simulation examples.
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
%
% Inputs
%  par          Struct with the following fields (all are set to default
%               values when they are not provided):
%   N           Size of KC population
%   M           Size of PN population
%   kcpn        KC-PN connection numbers. 1 KC connects to kcpn PNs
%   sqlh        Number of clusters per sequence
%   tini        Initial time. Meaningless for autonomous systems
%   tend        Final time.
%   dt          Time-step size
%   num_seqs    Number of sequences to be embedded in the connectivity
%               matrix.
%   Full        1 or 0. 1 for Full response, with the three phases (first
%               sequential, steady state, second sequential). 0 for just
%               the first sequential phase.
%   Plot        1 or 0. 1 for plotting data and inversion.
%   PlotPN      1 or 0. 1 for plotting the response of every PN in a single
%               window. Not recommended for large population sizes.
%   Qn          Precision of the hidden states
%   Rn          Precision of the observed states
%   ObsNoise    Variance of the white noise in the observation equation
%   HidNoise    Variance of the white noise in the hidden states equation
%   ClusSiz     Number of neurons per Cluster
%   Cutoff      Only for Full=1. Selects when, in time, the first
%               sequential phase ends and the steady state begins.
%   Repeat      Only for Full=1. Determines how long the steady state phase
%               is, in units of dt.
%
% Outputs:
%  data         Structure with the KC (data.x) and PN (data.y) activity
%  par          Structure with all the parameters used in the simulation.
%               Additionally, it contains the following fields:
%   clusters    Clusters used during the simulation. Each row is a cluster
%               and the numbers represent the KCs. Sequences are not
%               separated, however, the first sqlh rows are the first
%               sequence, and so on.
%   rho         Connectivity matrix.
%   obs_matrix  Observation matrix.
%  Kal          Structure with the outputs of the Kalman Filter. For
%               details, see the UKF.m file. The main field of the
%               structure is Kal.mX, which is the inverted KC activity.

function [data, par, Kal] = CSHS(par)
%% Simulation and population parameters
%==========================================================================
try,   N            = par.N;        catch,  N           = 100;      end
try,   M            = par.M;        catch,  M           = 30;       end
try,   kcpn         = par.kcpn;     catch,  kcpn        = 20;       end
try,   sqlh         = par.sqlh;     catch,  sqlh        = 5;        end
try,   tini         = par.tini;     catch,  tini        = 0;        end
try,   tend         = par.tend;     catch,  tend        = 90;       end
try,   dt           = par.dt;       catch,  dt          = 0.2;      end
try,   num_seqs     = par.num_seqs; catch,  num_seqs    = 2;        end
try,   Full         = par.Full;     catch,  Full        = 0;        end
try,   Plot         = par.Plot;     catch,  Plot        = 1;        end
try,   PlotPN       = par.PlotPN;   catch,  PlotPN      = 0;        end
try,   Qn           = par.Qn;       catch,  Qn          = 0.1;      end
try,   Rn           = par.Rn;       catch,  Rn          = 0.001;    end
try,   ObsNoise     = par.ObsNoise; catch,  ObsNoise    = 0.2;      end
try,   HidNoise     = par.HidNoise; catch,  HidNoise    = 0.0;      end
try,   ClusSiz      = par.ClusSiz;  catch,  ClusSiz     = 3;        end
try,   Cutoff       = par.Cutoff;   catch,  Cutoff      = 300;      end
try,   Repeat       = par.Repeat;   catch,  Repeat      = 1000;     end

%% Generate clusters, connectivity matrix and observation matrix
%==========================================================================

try
    [clusters,~]    = cgenerate_clusters(N,num_seqs,sqlh,ClusSiz);
catch err
    if isequal(err.message,'Excessive repetition')
        [clusters,~]    = cgenerate_clusters(N,num_seqs,sqlh,ClusSiz);
    else
        rethrow(err)
    end
end
 
rho             = cget_rho(N, sqlh, clusters);
[~, obs_matrix] = cget_neurons(N,M,clusters,kcpn,1);
    
%% Begin data generation
%==========================================================================
    
seqs = cell(2,1);
if Full
    i2vec = 0:1;
else
    i2vec = 0;
end

% Solve
%--------------------------------------------------------------------------
for i2=i2vec
    
    % Set initial conditions between first and second clusters in sequence
    %----------------------------------------------------------------------
    S = ones(N,1)*1e-7;
    S(clusters(sqlh*i2+1,clusters(sqlh*i2+1,:)~=0)) = 0.5;
    S(clusters(sqlh*i2+2,clusters(sqlh*i2+2,:)~=0)) = 0.5;

    % Differential equation that gives you dS
    %----------------------------------------------------------------------
    dS = @(S,V,dummy) S.*(1 -rho*S) +V.^2;

    xout = zeros(N,(tend-tini)/dt);
    xout(:,1) = S;
    t = 1;
    for i=tini:dt:tend
        V= HidNoise*randn(N,1)+0.0001; %Noise in the hidden states
        S = S + dS(S,V,0)*dt;
        t = t+1;
        xout(:,t) = S;
    end
    seqs{i2+1} = xout; %Save for later
end

if Full
    xout = seqs{1}(:,1:Cutoff);
    xout = [xout, repmat(seqs{1}(:,Cutoff),1,Repeat)];
    xout = [xout, seqs{2}];
else
    xout = seqs{1};
end

t   = numel(xout(1,:));
vt  = (1:t)/(tend-tini)*dt*1500; % Escale time to match experimental results


% Set-up Figure and plot data
%--------------------------------------------------------------------------
if Plot
    fig_old = findobj('name','Plots of Simulation Results');
    if isempty(fig_old)
        fig1 = figure(1);
        set(fig1,'name','Plots of Simulation Results');
        Ss = get(0,'screensize');
        Ss = [Ss(1) + 100,Ss(2) + 100,Ss(3)*0.3, Ss(4)*0.8];
        set(fig1,'position',Ss);
    else
        fig1 = figure(fig_old(1));
        clf(fig1);
    end
    subplot(4,1,1)
    plot(vt,xout');
    title('Generated KC data');
end



%% Projection to the PNs using equation 3 from the main paper
%==========================================================================

% Observation equation
%--------------------------------------------------------------------------
out = @(X,W,dummy) 4*obs_matrix*X + W;

W = ObsNoise*randn(M,t);
Y = out(xout,W,0);

if Plot
    figure(fig1);
    subplot(4,1,2);
    plot(vt,Y');
    title('Generated PN data');
end

%Generate subplot with all neurons
%--------------------------------------------------------------------------
% Note: not recommended for large population sizes
if PlotPN
    fig_old = findobj('name','Firing rates of all PNs');
    if isempty(fig_old)
        fig2 = figure(2);
        set(fig2,'name','Firing rates of all PNs');
        Sss = [Ss(1)+Ss(3)+20, Ss(2), Ss(3)/0.3*0.6, Ss(4)];
        set(fig2,'position',Sss);
    else
        fig2 = figure(fig_old(1));
        clf;
    end
    for i=1:M
        maxY = max(max(Y));
        subplot(ceil(M/3),3,i)
        plot(vt,Y(i,:)');
        line([0 vt(end)],[1 1], 'color', 'r');
        axis([0 vt(end) 0 maxY]);
        
    end
end


%% Kalman Filter
%==========================================================================

Q       =   Qn*eye(N); %best 0.1
R       =   Rn*eye(M); %best 0.001
P0      =   eye(N);
x0      =   zeros(N,1);
% Other possibilities for initial conditions for the inference.
% x0   	=  0.5*ones(N,1);
% x0      =   rand(N,1);

T = (1:t)*dt;

opt.fopt    = [];
opt.gopt    = [];
opt.alpha   = 0.01;
opt.kappa   = 3-N;
opt.beta    = 2;
opt.dt      = dt;

ffun = dS;
gfun = out;

[Kal.mX, Kal.P, Kal.peY,  Kal.peX, Kal.mYpred, Kal.mXpred, Kal.Ppred,...
    Kal.Pypred, Kal.nposdeferr] =...
    UKF(Y, T, x0, P0, ffun, gfun, Q, R, opt);

% Plotting
%--------------------------------------------------------------------------
eucdist = sqrt(sum( (xout(1:N,:)' - Kal.mX(1:N,:)').^2,2));

if Plot
    figure(fig1)
    subplot(4,1,3)
    title('Inferred activity');
    plot(vt,Kal.mX(1:N,:)','--');
    subplot(4,1,4)
    plot(vt,eucdist);
    axis([0,vt(end),0, max(eucdist)]);
    title('Error in prediction')
end

% Writing outputs
%--------------------------------------------------------------------------
data.x = xout;
data.y = Y;

par.N           = N;
par.M           = M;
par.kcpn        = kcpn;
par.sqlh        = sqlh;
par.tini        = tini;
par.tend        = tend;
par.dt          = dt;
par.num_seqs    = num_seqs;
par.Full        = Full;
par.Plot        = Plot;
par.PlotPN      = PlotPN;
par.Qn          = Qn;
par.Rn          = Rn;
par.ObsNoise    = ObsNoise;
par.HidNoise    = HidNoise;
par.ClusSiz     = ClusSiz;
par.Cutoff      = Cutoff;
par.Repeat      = Repeat;

par.clusters    = clusters;
par.rho         = rho;
par.obs_matrix  = obs_matrix;


end
