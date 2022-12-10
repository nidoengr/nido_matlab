function [Centroids, CentroidsUnc] = calc_centroids_from_clusters(ydetc_pxl,xdetc_pxl,IdxClust,IB_hat)

clusters = unique(IdxClust);
nClusters = max(clusters);
Centroids = nan(nClusters,2);
CentroidsUnc = nan(nClusters,2);

sig=10; %read noise TODO! MAKE AN INPUT PARM
% Loop through each cluster
for iCluster = 1:nClusters

   % Get the indecies of the cluster
   idxClustMod = IdxClust; 
   % Find the current cluster index
   idxClustMod(IdxClust~=iCluster) = zeros(size(idxClustMod(IdxClust~=iCluster)));
   idxClustMod(IdxClust==iCluster) = ones(size(idxClustMod(IdxClust==iCluster)));
   idxClustMod = logical(idxClustMod);
   % pull out all the x y coordinates of the current cluster
   xi = xdetc_pxl(idxClustMod);
   yi = ydetc_pxl(idxClustMod);

   % Determine the signal at each x y
   si = IB_hat(xi,yi); % s = sum( s_i ) from i to m (m_eta x m_nu)
   s  = sum(si,"all");
   %stdxi = std(xi);
   %stdyi = std(yi);
   m = numel(xi)*numel(yi);
   % Centroid of the signal density
   xhat_c = sum( xi .* si ) / sum ( si ); % sum(si,'all');
   yhat_c = sum( yi .* si ) / sum ( si ); % sum(si,'all');
   
   Centroids(iCluster,:) = [xhat_c, yhat_c];
   
   CentroidsUnc(iCluster,:)=sqrt([sum(xi.^2,"all"), sum(yi.^2,"all")].*...
      sig^2/(s^2+m*sig^2) + (1/(24*m*sqrt(pi))));


%    [sum(xi), sum(yi)] ./(si.^2 + numel(xi)*numel(yi)*sig_readnoise)
% 
%    sig3 = 1/2 * 1/sqrt(pi) * 1/12 * 1/
   
end

