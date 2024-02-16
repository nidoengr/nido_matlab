function [xyzArc] = eul_rot_ang_pnts(ang_deg,rotCase,shiftOpt_deg)
arguments
   ang_deg(1,1) double = 20;
   rotCase (1,1) char = 'X'; %'Y', 'Z'
   shiftOpt_deg(1,1) double  = 0;
end
tstep=0.5;


switch rotCase
   case 'X'
      t=(0:tstep:ang_deg)+shiftOpt_deg;
      xR1C = zeros(size(t));
      yR1C = cosd(t);
      zR1C = sind(t);
      xyzArc = [xR1C; yR1C; zR1C];

   case 'Y'
      t=(0:tstep:ang_deg)+shiftOpt_deg;
      xR12C = sind(t);
      yR12C = zeros(size(t));
      zR12C = cosd(t);
      xyzArc = [xR12C; yR12C; zR12C];
   case 'Z'
      t=(0:tstep:ang_deg)+shiftOpt_deg;
      xR12C = sind(t);
      yR12C = cosd(t);
      zR12C = zeros(size(t));
      xyzArc = [xR12C; yR12C; zR12C];
   otherwise
      error()
end
