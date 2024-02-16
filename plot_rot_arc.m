function plot_rot_arc(f,DCM,rotAng_deg,charCase,rotId, col)
arguments
   f matlab.ui.Figure = figure();
   DCM (3,3) double = eye(3);
   rotAng_deg (1,1) double = 20;
   charCase (1,:) char = 'Y'
   rotId (1,1) double = 1;
   col (1,3) = [1 0 0];
end
xyzArc = eul_rot_ang_pnts(rotAng_deg,charCase);
figure(f);
for i=1:size(xyzArc,2)
   XYZArc = DCM* xyzArc(:,i);
   % XYZArc = squeeze(R(2,:,:))'* xyzArc(:,i);
   
   plot3(gca,XYZArc(1),XYZArc(2),XYZArc(3),'.','Color',col);
   if i==floor(size(xyzArc,2)/2)
      mltp=1.2;
      txtx=sprintf('\\theta_%d = %0.2f',rotId,rotAng_deg);
      text(XYZArc(1)*mltp,XYZArc(2)*mltp,XYZArc(3)*mltp,txtx);
   end
end

xyzArc = eul_rot_ang_pnts(rotAng_deg,charCase,90);
for i=1:size(xyzArc,2)
   XYZArc = DCM* xyzArc(:,i);
   % XYZArc = squeeze(R(2,:,:))'* xyzArc(:,i);
   plot3(gca,XYZArc(1),XYZArc(2),XYZArc(3),'.','Color',col);
   if i==floor(size(xyzArc,2)/2)
      mltp=1.2;
      txtx=sprintf('\\theta_%d = %0.2f',rotId,rotAng_deg);
      text(XYZArc(1)*mltp,XYZArc(2)*mltp,XYZArc(3)*mltp,txtx);
   end
end