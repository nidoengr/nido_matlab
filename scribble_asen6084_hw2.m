% [IdxClust C] = kmeans(IB_hat_flat,5);
% [IdxClust C] = imsegkmeans(IB_hat_flat,1);
% IdxClust = reshape(IdxClust,size(IB_hat,1),size(IB_hat,2));
% B = labeloverlay(IB_hat,IdxClust,'Transparency',1);
% B = im2gray(B);
% B(B<0) = zeros(size(B(B<0)));
% B= IB_hat;
% figure("Name","Label Overlay")
% imshow(B)
% imshow(B - uint8(IB_hat));
% imshow(IB_hat - IM_true); % All black
% imshow(B-uint8(IM_true)); % All black
% imshowpair(mat2gray(RGB),IB_hat,'diff')
% IB_hat_clust_flat = IB_hat(:);

%%
% % Figure for dbscan? Doesnt work?
% figure()
% hold all;
% for i=1:size(B,1)
%    if IdxClust(i) == -1
%       plot(B(i,1),B(i,2),'r.','DisplayName','Outliers'); hold on;
%    else
%       if C(i)==1
%          plot(B(i,1),B(i,2),'b.','MarkerSize',5,'DisplayName','Obj'); hold on;
%       else
%          plot(B(i,1),B(i,2),'g.','MarkerSize',5,'DisplayName','Obj'); hold on;
%       end
%    end
% end
%% Plot
% figure()
% % imshow(IB_hat - IM_true)
% % imagesc(IM_true)
% imagesc(IB_hat)

%% Lets assume we got the cluster detections working

% MTWs = [[150 250 320 380]
%         [70 105 255 295]
%         [70 120 180 230]
%         [115 155 110 150]];
%         %[? ? ? ?]]
% 
% for iDetect = 1 : numel(MTWs(:,1))
%    i0 = MTW(iDetect,1);
%    iF = MTW(iDetect,2);
%    j0 = MTW(iDetect,3);
%    jF = MTW(iDetect,4);
%    TW = B(i0:iF,j0:jF);
%    xSum = 0;
% %    for i=i0:iF
% %       for j=j0:jF
% %          xSum = xSum + i
% %       end
% %    end
% end
