% [neurons,obs_matrix] = cget_neurons(N,M,clusters,kcpn, rep)
% Function that generates the observation matrix (equation 3 of main text)
% randomly, connecting each KC randomly to a set number of PNs.
%
% Inputs
%   N               Size of KC population
%   M               Size of PN population
%   clusters        Matrix containing all the clusters in the connectivity
%                   matrix.
%   kcpn            Number of PNs that receive input from each KC.
%   rep             1 or 0. Switch to allow or disallow the repetition of
%                   PNs in subsequent clusters. That is, whether a PN can
%                   be connected to KCs who appear in two consecutive
%                   clusters in a sequences. Beware, if kcpn is too big,
%                   repetition might be unavoidable.
%
% Outputs:
%   neurons         [N,kcpn] matrix containing a list of all the PNs that
%                   connect to each KC. For legacy applications.
%   obs_matrix      Observation matrix for the observation equation.


function [neurons,obs_matrix] = cget_neurons(N,M,clusters,kcpn, rep)

neurons = zeros(N,kcpn);
obs_matrix = zeros(M,N);
appear = zeros(N,1);

[clnum, clele] = size(clusters);


for i=1:clnum
    pool = randperm(M);
    negpool = []; 
    for j=1:numel(negpool)
        pool(pool==negpool(j)) = [];
    end

    
    %Remove from pool the neurons from the elements in current cluster
    %----------------------------------------------------------------------
    if rep==0
        for j=1:clele
            for k=1:kcpn
                pool(pool==neurons(clusters(i,j),k)) = [];
            end
        end
    end
    % Stores which PNs have been together in a KC cluster.
    %----------------------------------------------------------------------
    for j=1:numel(clusters(i,clusters(i,:)~=0))
        if appear(clusters(i,j))==0
            neurons(clusters(i,j),:) = pool(1:kcpn);
            if rep==0, pool(1:kcpn) = [];
            else pool = randperm(M);
            end
            appear(clusters(i,j)) = 1;
        end
    end
end

%Fill the rest of the entries with random numbers.
%--------------------------------------------------------------------------
%Note: For legacy applications. Does not affect obs_matrix.
for i=1:N
    if neurons(i,1)==0
        neurons(i,:) = ceil(M*rand(1,kcpn));
    end
end

%Create the obs_matrix for the out function
%--------------------------------------------------------------------------
for i=1:N
    for j=1:kcpn
        obs_matrix( neurons(i,j),i ) = 1;
    end
end
end
