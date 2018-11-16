% Calculate vector of inverse dynamics base and joint forces/torques for a
% S5RRRRR1
% Use numerical implementation of recursive Newton-Euler Algorithm
%
% Input:
% q [5x1]
%   Joint Angles [rad]
% qD [5x1]
%   Joint Velocities [rad/s]
% qDD [5x1]
%   Joint Accelerations [rad/s]
% phi_base [3x1]
%   Base orientation in world frame. Expressed with XYZ-Euler angles
% xD_base [6x1]
%   time derivative of r_base and phi_base
% xDD_base [6x1]
%   second time derivative of r_base and phi_base
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% tau [(6+5)x1]
%   Base forces and joint torques required to compensate inverse dynamics.
%   Base moments in generalized coordinates (Euler-XYZ)
%
% See also:
% IMES-Robotik-Toolbox: dynamics/robot_tree_invdyn_floatb_eulxyz_nnew_vp1.m

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 14:53
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [tau, v_i_i_ges, w_i_i_ges] = S5RRRRR1_invdyn_floatb_eulxyz_nnew_vp1(q, qD, qDD, phi_base, xD_base, xDD_base, ...
  pkin, m_num, rSges_num_mdh, Icges_num_mdh)

%%Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1)
%$cgargs  zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(q) && all(size(q) == [5 1]), ...
  'S5RRRRR1_invdyn_floatb_eulxyz_nnew_vp1: q has to be [5x1] (double)');
assert(isreal(qD) && all(size(qD) == [5 1]), ...
  'S5RRRRR1_invdyn_floatb_eulxyz_nnew_vp1: qD has to be [5x1] (double)');
assert(isreal(qDD) && all(size(qDD) == [5 1]), ...
  'S5RRRRR1_invdyn_floatb_eulxyz_nnew_vp1: qDD has to be [5x1] (double)');
assert(isreal(phi_base) && all(size(phi_base) == [3 1]), ...
  'S5RRRRR1_invdyn_floatb_eulxyz_nnew_vp1: phi_base has to be [3x1] (double)');
assert(isreal(xD_base) && all(size(xD_base) == [6 1]), ...
  'S5RRRRR1_invdyn_floatb_eulxyz_nnew_vp1: xD_base has to be [6x1] (double)');
assert(isreal(xDD_base) && all(size(xDD_base) == [6 1]), ...
  'S5RRRRR1_invdyn_floatb_eulxyz_nnew_vp1: xDD_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_invdyn_floatb_eulxyz_nnew_vp1: Kinematic parameters pkin have to be [6x1] (double)');
assert(isreal(m_num) && all(size(m_num) == [6 1]), ...
  'S5RRRRR1_invdyn_floatb_eulxyz_nnew_vp1: m_num has to be [6x1] (double)');
assert(isreal(rSges_num_mdh) && all(size(rSges_num_mdh) == [6,3]), ...
  'S5RRRRR1_invdyn_floatb_eulxyz_nnew_vp1: rSges_num_mdh has to be [6x3] (double)');
assert(isreal(Icges_num_mdh) && all(size(Icges_num_mdh) == [6 6]), ...
  'S5RRRRR1_invdyn_floatb_eulxyz_nnew_vp1: Icges_num_mdh has to be [6x6] (double)');

%% Init
tau_J = NaN(5,1);
R_W_0 = eulxyz2r(phi_base);
T_basevel = eulxyzjac(phi_base);
[v_mdh, sigma] = S5RRRRR1_structural_kinematic_parameters();


%% Vorwärts-Iteration
% Positionen. Berechne symbolisch (effizienter).
T_mdh = S5RRRRR1_joint_trafo_rotmat_mdh_sym_varpar(q, pkin);

% Geschwindigkeiten
v_i_i_ges = NaN(3,6);
w_i_i_ges = NaN(3,6);

% Beschleunigungen
vD_i_i_ges = NaN(3,6);
wD_i_i_ges = NaN(3,6);

% Anfangswerte: Geschwindigkeit und Beschleunigung der Basis
v_i_i_ges(:,1) = R_W_0'*xD_base(1:3);
w_i_i_ges(:,1) = R_W_0'*eulxyzD2omega(phi_base, xD_base(4:6));

vD_i_i_ges(:,1) = R_W_0'*xDD_base(1:3);
wD_i_i_ges(:,1) = R_W_0'*eulxyzDD2omegaD(phi_base, xD_base(4:6), xDD_base(4:6));

for i = 2:6
  % Nummer des Vorgänger-Segments
  j = v_mdh(i-1)+1; % Gelenk 1 führt zu Körper 2 (Matlab-Indizes) usw.

  % Temporäre Ausdrücke belegen
  v_j_j = v_i_i_ges(:,j);
  w_j_j = w_i_i_ges(:,j);
  vD_j_j = vD_i_i_ges(:,j);
  wD_j_j = wD_i_i_ges(:,j);
  R_j_i = T_mdh(1:3,1:3,i-1);
  r_j_j_i = T_mdh(1:3,4,i-1);

  % Berechnung
  % [3], Gl. 9.17
  w_i_i = R_j_i'*w_j_j;
  if sigma(i-1) == 0
    w_i_i = w_i_i + [0;0;1]*qD(i-1);
  end
  % [3], Gl. 9.18
  v_i_i = R_j_i'*( v_j_j + cross(w_j_j, r_j_j_i) );
  if sigma(i-1) == 1
    v_i_i = v_i_i + [0;0;1]*qD(i-1);
  end
  wD_i_i = R_j_i'*wD_j_j;
  if sigma(i-1) == 0
    wD_i_i = wD_i_i + [0;0;1]*qDD(i-1) + cross(R_j_i'*w_j_j, [0;0;1]*qD(i-1));
  end
  vD_i_i = R_j_i'*( vD_j_j + cross(wD_j_j, r_j_j_i) +cross(w_j_j, cross(w_j_j, r_j_j_i)) );
  if sigma(i-1) == 1
    vD_i_i = vD_i_i + [0;0;1]*qDD(i-1) + 2*cross(R_j_i'*w_j_j, [0;0;1]*qD(i-1));;
  end

  % Ausgabeausdrücke belegen
  v_i_i_ges(:,i) = v_i_i;
  w_i_i_ges(:,i) = w_i_i;
  vD_i_i_ges(:,i) = vD_i_i;
  wD_i_i_ges(:,i) = wD_i_i;
end


%% Rückwärts-Rekursion
f_i_i_ges = NaN(3,6);
n_i_i_ges = NaN(3,6);

for i = 6:-1:1
  % Temporäre Ausdrücke belegen
  vD_i_i = vD_i_i_ges(:,i);
  w_i_i = w_i_i_ges(:,i);
  wD_i_i = wD_i_i_ges(:,i);
  c_i = rSges_num_mdh(i,:)';
  I_i = inertiavector2matrix(Icges_num_mdh(i,:));

  % Dynamik-Terme
  F_i = m_num(i)*(vD_i_i + cross(wD_i_i, c_i) + cross( w_i_i,cross(w_i_i, c_i)) );
  N_i = I_i*wD_i_i + cross(w_i_i, I_i*w_i_i);

  f_i_i = F_i;
  n_i_i = N_i + cross(c_i, F_i);

  % Suche alle Nachfolger und addiere das Schnittmoment
  I_nf = find( (v_mdh == (i-1)) ) + 1;
   % Wähle diese Konstruktion um Schleifen mit variabler Länge zu vermeiden (Kompilierbarkeit)
  if ~isempty(I_nf)
    for tmp = 1:length(v_mdh) % Index des Nachfolgers
      j = I_nf(tmp);
      R_i_j = T_mdh(1:3,1:3,j-1);
      f_j_j = f_i_i_ges(:,j);
      n_j_j = n_i_i_ges(:,j);
      r_i_i_j = T_mdh(1:3,4,j-1);

      f_i_i = f_i_i + R_i_j*f_j_j;
      n_i_i = n_i_i + R_i_j*n_j_j + cross(r_i_i_j, R_i_j*f_j_j);
      if tmp == length(I_nf)
        break; % Abbruch. Alle Nachfolger untersucht.
      end
    end
  end

  % Ausgabeausdrücke belegen
  f_i_i_ges(:,i) = f_i_i;
  n_i_i_ges(:,i) = n_i_i;
end

%% Projektion auf die Gelenke
for i = 2:6
  if sigma(i-1) == 0 % Drehgelenk
    tau_J(i-1) = [0 0 1] * n_i_i_ges(:,i);
  else % Schubgelenk
    tau_J(i-1) = [0 0 1] * f_i_i_ges(:,i);
  end
end

%% Basis-Kraft
tau_B = [R_W_0*f_i_i_ges(:,1); T_basevel' * R_W_0*n_i_i_ges(:,1)];

%% Ausgabe
tau = [tau_B; tau_J];
