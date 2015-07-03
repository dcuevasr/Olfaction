%Function to obtain the rho connectivity matrix for a CSHS, using the
%clusters as a guide to set the connections. The values are those from
%Figure 4 of the main text.
%
% Inputs:
%  N            Size of the KC population
%  sqlh         Number of clusters per sequence
%  clusters     Matrix containing all the clusters in all the sequences.
%               Each row is a cluster.
%
% Outputs:
%  rho          [N,N] Connectivity matrix.

function rho = cget_rho(N,sqlh,clusters)

cltotal = size(clusters,1);

rho = 2*ones(N); % Base level of high all-to-all inhibition.

% Set the values for all the neurons, cluster by cluster
%--------------------------------------------------------------------------
for i=1:cltotal
    frac = 1/(numel(clusters(i,clusters(i,:)~=0)));

    vals = [frac/10, frac];


    for j=1:numel(clusters(i,clusters(i,:)~=0))
        %How much a cluster inhibits the previous cluster
        %------------------------------------------------------------------
        if mod(i-1,sqlh)>0
            rho(clusters(i-1,clusters(i-1,:)~=0),clusters(i,j) )= 2*vals(2);
        end
        %How much the members of a cluster inhibit each other
        %------------------------------------------------------------------
        rho(clusters(i,clusters(i,:)~=0) ,clusters(i,j)) = vals(2);
        
        %How much the cluster inhibits the next one
        %------------------------------------------------------------------
        if mod(i,sqlh)>0
            rho(clusters(i+1,clusters(i+1,:)~=0),clusters(i,j) )= vals(1);
        end
    end
end

%NOTE: The WITHIN-CLUSTER inhibition must be bigger than  
%TO-PREVIOUS-CLUSTER inhibition

end

