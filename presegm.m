function [adjmat] = presegm( neighb, affinity, k, adj_thr)
%PRESEGM performs a screening of relevant EEG predictors by grouping the candidate set into clusters to improve
%   efficiency and to reduce computational load. 
%  adjmat=PRESEGM(neighb, affinity, k, adj_thr) where the input adjacency matrix 'neighb'
%  exclude non-neighbours and the affinities between all pairs of EEG channels form the 'affinity'
%  matrix neighb. The matrix 'adjmat' returns the channels belong to the same cluster.
%  Parameter  Value:
%  - adj_thr is the threshold of minimum correlation being considered
%  - k is a scale parameter, in that a larger k causes a preference for larger clusters.
%   Example:  Spatial clustering of EEG sensors.
%       x = randn(neighbM, affinityM, 0.2, 0.4);       
%   References:
%      Cross Multivariate Correlation Coefficients as Screening Tool for Analysis of Concurrent EEG-fMRI Recordings
 
if nargin<2
   error('Please provide adjacency matrix and affinity matrix');
end
if size(neighb)~=size(affinity)
  error(' matrix and affinity matrix have to be equal in size');
end
if  ~k 
     k = 0.2;
end
if ~adj_thr
    adj_thr = 0.4;
end
affinitymat = neighb.*affinity;
[W_descend,IND] = sort(affinitymat(:),'descend');
nch = size(neighb,1);
[I,J] = ind2sub([nch,nch],IND);
adjmat = zeros(nch,nch);

for i=1:length(W_descend)
    if W_descend(i)>0
        Ci = [I(i),find(adjmat(I(i),:)~=0)];
        Cj = [J(i),find(adjmat(J(i),:)~=0)];
        if ~isempty(setxor(Ci,Cj))
            if length(Ci)>1 || length(Cj)>1 
                Omega_Ci = sqrt(det(affinity(Ci,Ci)));   %compute the internal MUC of Ci
                Omega_Cj = sqrt(det(affinity(Cj,Cj)));   %compute the internal MUC of Cj
                MInt = min(Omega_Ci+k/length(Ci), Omega_Cj+k/length(Cj));
                Omega = sqrt(1-W_descend(i)^2);
                if Omega<=MInt  % decision-making for merging
                    adjmat(I(i),J(i)) = 1;
                    adjmat(J(i),I(i)) = 1;                    
                end
            else
                if W_descend(i)>adj_thr
                    adjmat(I(i),J(i)) = 1;
                    adjmat(J(i),I(i)) = 1;
                end
            end
        end
    end
end