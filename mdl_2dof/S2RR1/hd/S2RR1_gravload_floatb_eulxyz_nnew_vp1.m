% Calculate vector of base and joint gravitational load for
% S2RR1
% Use numerical implementation of recursive Newton-Euler Algorithm
%
% Input:
% q [2x1]
%   Joint Angles [rad]
% phi_base [3x1]
%   Base orientation in world frame. Expressed with XYZ-Euler angles
% g_world [3x1]
%   gravitation vector in world frame [m/s^2]
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
%   mass of all robot links (including the base)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
%
% Output:
% tau_g [(6+2)x1]
%   Base forces and joint torques required to compensate gravitational load.
%   Base moments in generalized coordinates (Euler-XYZ)
%
% See also:
% IMES-Robotik-Toolbox: dynamics/robot_tree_invdyn_floatb_eulxyz_nnew_vp1.m

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 16:45
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tau_g = S2RR1_gravload_floatb_eulxyz_nnew_vp1(q, phi_base, g_world, pkin, m_num, rSges_num_mdh)

%%Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(3,1),zeros(3,1),zeros(1,1),zeros(3,1),zeros(3,3)}
assert(isreal(q) && all(size(q) == [2 1]), ...
  'S2RR1_gravload_floatb_eulxyz_nnew_vp1: q has to be [2x1] (double)');
assert(isreal(phi_base) && all(size(phi_base) == [3 1]), ...
  'S2RR1_gravload_floatb_eulxyz_nnew_vp1: phi_base has to be [3x1] (double)');
assert(isreal(g_world) && all(size(g_world) == [3 1]), ...
  'S2RR1_gravload_floatb_eulxyz_nnew_vp1: g_world has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_gravload_floatb_eulxyz_nnew_vp1: Kinematic parameters pkin have to be [1x1] (double)');
assert(isreal(m_num) && all(size(m_num) == [3 1]), ...
  'S2RR1_gravload_floatb_eulxyz_nnew_vp1: m_num has to be [3x1] (double)');
assert(isreal(rSges_num_mdh) && all(size(rSges_num_mdh) == [3,3]), ...
  'S2RR1_gravload_floatb_eulxyz_nnew_vp1: rSges_num_mdh has to be [3x3] (double)');

%% Init
tau_J = NaN(2,1);
R_W_0 = eulxyz2r(phi_base);
T_basevel = eulxyzjac(phi_base);
vD_i_i_ges = NaN(3,3);
[v_mdh, sigma] = S2RR1_structural_kinematic_parameters();

%% Vorwärts-Iteration
% Positionen. Berechne symbolisch (effizienter).
T_mdh = S2RR1_joint_trafo_rotmat_mdh_sym_varpar(q, pkin);

% Anfangswerte: Geschwindigkeit und Beschleunigung der Basis
vD_i_i_ges(:,1) = -R_W_0'*g_world;

for i = 2:3
  % Nummer des Vorgänger-Segments
  j = v_mdh(i-1)+1; % Gelenk 1 führt zu Körper 2 (Matlab-Indizes) usw.

  % Temporäre Ausdrücke belegen
  vD_j_j = vD_i_i_ges(:,j);
  R_j_i = T_mdh(1:3,1:3,i-1);

  % Berechnung
  vD_i_i = R_j_i'*vD_j_j;

  % Ausgabeausdrücke belegen
  vD_i_i_ges(:,i) = vD_i_i;
end


%% Rückwärts-Rekursion
f_i_i_ges = NaN(3,3);
n_i_i_ges = NaN(3,3);

for i = 3:-1:1
  % Temporäre Ausdrücke belegen
  vD_i_i = vD_i_i_ges(:,i);
  c_i = rSges_num_mdh(i,:)';

  % Dynamik-Terme
  F_i = m_num(i)*vD_i_i;

  f_i_i = F_i;
  n_i_i = cross(c_i, F_i);

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
for i = 2:3
  if sigma(i-1) == 0 % Drehgelenk
    tau_J(i-1) = [0 0 1] * n_i_i_ges(:,i);
  else % Schubgelenk
    tau_J(i-1) = [0 0 1] * f_i_i_ges(:,i);
  end
end

%% Basis-Kraft
tau_B = [R_W_0*f_i_i_ges(:,1); T_basevel' * R_W_0*n_i_i_ges(:,1)];

%% Ausgabe
tau_g = [tau_B; tau_J];
