function DCM = eul_vec_to_dcm(eulVec_rd, seqChar)
arguments
   eulVec_rd (1,:) double = [20 -10 120]*pi/180;
   seqChar (1,:) char = '321';
   % gaken from hw2-w3-concepCheck3-principalRotationAddition
end

% Use syms to derive the matrix multiplication equations
% (t1 t2 t3) = (theta1, theta2, theta3)
syms t1 t2 t3

% Define the rotation matrices as anonymous functions
M1_func = @(theta) [1 0 0;0 cos(theta) sin(theta);0 -sin(theta) cos(theta)];
M2_func = @(theta) [cos(theta) 0 -sin(theta);0 1 0;sin(theta) 0 cos(theta)];
M3_func = @(theta) [cos(theta) sin(theta) 0;-sin(theta) cos(theta) 0;0 0 1];

%% Rotation Sequence
caseInput = seqChar;
isUnhandledCase = false;
switch caseInput
   case '321'
      C = M1_func(t3) * M2_func(t2) * M3_func(t1);
   case '313'
      C = M3_func(t3) * M1_func(t2) * M3_func(t1);
   case '323'
      C = M3_func(t3) * M2_func(t2) * M3_func(t1);
   case '232'
      C = M2_func(t3) * M3_func(t2) * M2_func(t1);
   case '231'
      C = M1_func(t3) * M3_func(t2) * M2_func(t1);
   case '123'
      C = M3_func(t3) * M2_func(t2) * M1_func(t1);
   case '121'
      C = M1_func(t3) * M2_func(t2) * M1_func(t1);
   otherwise
      warning('Euler sequence not recognized or handled in switch case. Update if needed.')
      C = nan(3,3);
      isUnhandledCase = true;
end
%%
t1 = eulVec_rd(1);
t2 = eulVec_rd(2);
t3 = eulVec_rd(3);

DCMtemp = vpa(subs(C));

% This is a workaround to get the dcm as double
DCM = nan(3,3);
if ~isUnhandledCase
   for i=1:3, for j=1:3, DCM(i,j)=double(DCMtemp(i,j)); end, end
end