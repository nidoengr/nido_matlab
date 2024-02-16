function visualize_eul_rot_seq( s, eulVec_deg, seqChar,  f , startDCM, colSeq)
arguments
   s matlab.graphics.chart.primitive.Surface;
   eulVec_deg(1,3) double = [10 20 30];
   seqChar(1,3) char = '121';
   f matlab.ui.Figure = figure();
   startDCM (3,3)=eye(3);
   colSeq (3,3) double = [[1 0 0];[0 0 1];[0 1 0]];
end
% close all;
R = nan(3,3,3);
axMult = 1.3;
txMult = 1.45;

switch seqChar
   case '121'
      R(1,:,:) = rotx(eulVec_deg(1));
      R(2,:,:) = roty(eulVec_deg(2));
      R(3,:,:) = rotx(eulVec_deg(3));
      C = rotx(eulVec_deg(1)) * ...
          roty(eulVec_deg(2)) * ...
          rotx(eulVec_deg(3)) ;
      tstep=0.5;
      t=0:tstep:eulVec_deg(1);
      xR1C = zeros(size(t));
      yR1C = sind(t);
      zR1C = cosd(t);
   otherwise
      warning('Case not supported. Please develop')
end

%% Draw Rot 1
R1N = startDCM'*squeeze(R(1,:,:)); %* ICFM;
% [isGood, prAng_rd, prUvec] = determine_principal_rot(R1N);
% prAng_deg=prAng_rd*180/pi;
% if isGood
%    rotate(s,prUvec,prAng_deg,[0 0 0])
% end

% This is the axes plotting
AXR1M = R1N*axMult;
plot3([0 AXR1M(1,1)],[0 AXR1M(2,1)],[0 AXR1M(3,1)],'-','LineWidth',4,'Color',colSeq(1,:))
plot3([0 AXR1M(1,2)],[0 AXR1M(2,2)],[0 AXR1M(3,2)],'-','LineWidth',4,'Color',colSeq(1,:))
plot3([0 AXR1M(1,3)],[0 AXR1M(2,3)],[0 AXR1M(3,3)],'-','LineWidth',4,'Color',colSeq(1,:))
ATR1M = R1N*txMult;
text(ATR1M(1,1),AXR1M(2,1),AXR1M(3,1),'^{R1F}X''','FontWeight','bold')
text(AXR1M(1,2),ATR1M(2,2),AXR1M(3,2),'^{R1F}Y''','FontWeight','bold')
text(AXR1M(1,3),AXR1M(2,3),AXR1M(3,3),'^{R1F}Z''','FontWeight','bold')
%% Draw Rot 1 Angles - TODO! FIX ME! not doing this the smartest way
plot_rot_arc(f,startDCM',eulVec_deg(1),'X',1,colSeq(1,:));

%% Draw Rot 2
R2N = startDCM'*[squeeze(R(1,:,:))*squeeze(R(2,:,:))]; %* ICFM;
% R12FM = squeeze(R(2,:,:));
% [isGood, prAng_rd, prUvec] = determine_principal_rot(R2N);
% % [isGood, prAng_rd, prUvec] = determine_principal_rot(squeeze(R(2,:,:)));
% prAng_deg=prAng_rd*180/pi;
% if isGood
%    rotate(s,prUvec,prAng_deg,[0 0 0])
% end

% This is the axes plotting
AXR12M = R2N*axMult;
plot3([0 AXR12M(1,1)],[0 AXR12M(2,1)],[0 AXR12M(3,1)],'b-','LineWidth',4,'Color',colSeq(2,:))
plot3([0 AXR12M(1,2)],[0 AXR12M(2,2)],[0 AXR12M(3,2)],'b-','LineWidth',4,'Color',colSeq(2,:))
plot3([0 AXR12M(1,3)],[0 AXR12M(2,3)],[0 AXR12M(3,3)],'b-','LineWidth',4,'Color',colSeq(2,:))
ATR12M = R2N*txMult;
text(ATR12M(1,1),AXR12M(2,1),AXR12M(3,1),'^{R12F}X''''','FontWeight','bold')
text(AXR12M(1,2),ATR12M(2,2),AXR12M(3,2),'^{R12F}Y''''','FontWeight','bold')
text(AXR12M(1,3),AXR12M(2,3),AXR12M(3,3),'^{R12F}Z''''','FontWeight','bold')

%% Draw Rot 2 Angles - TODO! FIX ME! not doing this the smartest way
plot_rot_arc(f,R1N,eulVec_deg(2),'Y',2,colSeq(2,:));

%% Draw Rot 3
R3N = startDCM'*[squeeze(R(1,:,:))*squeeze(R(2,:,:))*squeeze(R(3,:,:))]; %* ICFM;
[isGood, prAng_rd, prUvec] = determine_principal_rot(startDCM*R3N');%*squeeze(R(1,:,:)));
prAng_deg=prAng_rd*180/pi;
if isGood
   rotate(s,prUvec,prAng_deg,[0 0 0])
end

% This is the axes plotting
AXR123M = R3N*axMult;
plot3([0 AXR123M(1,1)],[0 AXR123M(2,1)],[0 AXR123M(3,1)],'g-','LineWidth',4,'Color',colSeq(3,:))
plot3([0 AXR123M(1,2)],[0 AXR123M(2,2)],[0 AXR123M(3,2)],'g-','LineWidth',4,'Color',colSeq(3,:))
plot3([0 AXR123M(1,3)],[0 AXR123M(2,3)],[0 AXR123M(3,3)],'g-','LineWidth',4,'Color',colSeq(3,:))
ATR123M = R3N*txMult;
text(ATR123M(1,1),ATR123M(2,1),ATR123M(3,1),'^{R123F}X''''''','FontWeight','bold')
text(ATR123M(1,2),ATR123M(2,2),ATR123M(3,2),'^{R123F}Y''''''','FontWeight','bold')
text(ATR123M(1,3),ATR123M(2,3),ATR123M(3,3),'^{R123F}Z''''''','FontWeight','bold')
%% Draw angles - TODO! FIX ME! not doing this the smartest way
plot_rot_arc(f,R2N,eulVec_deg(3),'X',3,colSeq(3,:));

%%
fh = gcf;
fh.WindowState='maximized';
