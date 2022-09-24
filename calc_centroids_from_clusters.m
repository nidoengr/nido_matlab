function Centroids = calc_centroids_from_clusters(ydetc_pxl,xdetc_pxl,IdxClust,IB_hat)

clusters = unique(IdxClust);
nClusters = max(clusters);
Centroids = nan(nClusters,2);
for iCluster = 1:nClusters

   idxClustMod = IdxClust; 
   idxClustMod(IdxClust~=iCluster) = zeros(size(idxClustMod(IdxClust~=iCluster)));
   idxClustMod(IdxClust==iCluster) = ones(size(idxClustMod(IdxClust==iCluster)));
   idxClustMod = logical(idxClustMod);
   xi = xdetc_pxl(idxClustMod);
   yi = ydetc_pxl(idxClustMod);
   si = IB_hat(xi,yi);

   xhat_c = sum( xi .* si ) / sum ( si );
   yhat_c = sum( yi .* si ) / sum ( si );
   Centroids(iCluster,:) = [xhat_c, yhat_c];
   
end

