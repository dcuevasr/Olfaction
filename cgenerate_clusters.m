% [mclusters,mconnect] = cgenerate_clusters(N,D,Z,M)
% Generates the requested number of clusters for a CSHS. Each cluster
% contains M elements and they're set so that each cluster will excite
% another cluster, thus creating a sequence. A proto-connectivity matrix is
% also created, which, for each element, contains a list of the elements it
% excites.
%
% This variant creates clusters trying to minimize repetition of neurons.
%   Inputs:
%       N               Dimension of the space
%       D               Number of desired sequences
%       Z               Number of clusters per sequence
%       M               Number of elements per cluster
%
%   Outputs:
%       mclusters       [N,M] Contains the N clusters with M elements each
%       mconnect        Tells you which elements are excited by a given
%                       element. Row is the element, columns are the
%                       elements that it excites. NxM


function [mclusters,mconnect] = cgenerate_clusters(N,D,Z,M)


C = D*Z;                %Total number of clusters to compute
Mpool = M*ones(Z,1);
Mmax = max(Mpool);
Mpool = repmat(Mpool,D,1);


clusters = cell(C,1);
connect  = cell(N,1);

cmates = cell(N,1);
AlreadySeen = zeros(N,1);
partition(1) = 1;


%% Generate clusters minimizing repetition of neurons across sequences
%==========================================================================

pool = randperm(N);
clusters{1} = pool(1:M);
AlreadySeen(clusters{1}) = AlreadySeen(clusters{1})+1;
for ii=2:C
    [SortedAS, SortedIndex] = sort(AlreadySeen);
    zz=1;
    for kk=1:N-1
        if SortedAS(kk)~=SortedAS(kk+1)
            partition(zz) = kk+1;
            zz = zz+1;
        end
    end
    if partition(1)~=0
        partition = [1,partition];
    end
    for kk=1:(numel(partition)-1)
        auxShuffle = randperm(partition(kk+1) - partition(kk));
        auxSortedIndex = SortedIndex(partition(kk):(partition(kk+1)-1));
        SortedIndex(partition(kk):(partition(kk+1)-1)) =...
                                            auxSortedIndex(auxShuffle);
    end
    clusters{ii} = SortedIndex(1:M)';
    flag=0;
    flag2=0;
    AddK=1;
    changeN = 1;
    while flag==0
        for jj=1:ii-1
            if isequal(sort(clusters{ii}),sort(clusters{jj}))
                flag2=1;
                break
            end
        end
        if flag2==1
            clusters{ii}(changeN) = SortedIndex(M+AddK);
            AddK = AddK+1;
            if M+AddK==numel(SortedIndex)
                changeN = changeN+1;
                AddK = 1;
            end
            if M+AddK==numel(SortedIndex) && changeN==M
                error('Excessive repetition',...
                    'Too much repetition led to bad results. Try again');
            end
            flag2 = 0;
        else
            flag=1;
        end
    end
    AlreadySeen(clusters{ii}) = AlreadySeen(clusters{ii})+1;
end


%% Turn output into matrices for legacy reasons.
%==========================================================================
max_clusters_size = 0;
for i=1:C
    max_clusters_size = max(max_clusters_size,numel(clusters{i}));
end
        
max_connect_size = 0;
for i=1:N
    max_connect_size = max(max_connect_size,numel(connect{i}));
end

mconnect = zeros(N,max_connect_size);
mclusters = zeros(Z,max_clusters_size);

for i=1:C
    for j=1:numel(clusters{i})
        mclusters(i,j) = clusters{i}(j);
    end
end

for i=1:N
    for j=1:numel(connect{i})
        mconnect(i,j) = connect{i}(j);
    end
end
end
        
        
        
        
        
