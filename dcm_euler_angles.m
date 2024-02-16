function C_321 = dcm_euler_angles(eulerAngles_rd, seqInt )
arguments
   eulerAngles_rd (1,3) = [0 0 0];
   seqInt (1,:) {mustBeMember(seqInt,{321,313,323})} = 321;
end
t1 = eulerAngles_rd(1);
t2 = eulerAngles_rd(2);
t3 = eulerAngles_rd(3);

ct2 = cos(eulerAngles_rd(2));
c21 = cos(eulerAngles_rd(2));