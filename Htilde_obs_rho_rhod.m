function Htilde_obs_sc = Htilde_obs_rho_rhod(x_sc, x_obs)

% Relative vecs
x_rel = x_sc - x_obs;

r_rel = x_rel(1:3);
v_rel = x_rel(4:6);

x = r_rel(1);
y = r_rel(2);
z = r_rel(3);

xdot = v_rel(1);
ydot = v_rel(2);
zdot = v_rel(3);

rdotv = dot(r_rel,v_rel);

rho = sqrt( dot(r_rel,r_rel) )

% rhodot = dot( (r_rel), (v_rel) ) / rho;

drhodx = x/rho;
drhody = y/rho;
drhodz = z/rho;
% drhodxdot = 0;
% drhodydot = 0;
% drhodzdot = 0;

drhodotdx = ( xdot*rho - rdotv*(x/rho) )/rho^2;
drhodotdy = ( ydot*rho - rdotv*(y/rho) )/rho^2;
drhodotdz = ( zdot*rho - rdotv*(z/rho) )/rho^2;
% drhodotdxdot = x/rho;
% drhodotdydot = y/rho;
% drhodotdzdot = z/rho;

Htilde_obs_sc = [-drhodx -drhody -drhodz ;...
    -drhodotdx -drhodotdy -drhodotdz ];