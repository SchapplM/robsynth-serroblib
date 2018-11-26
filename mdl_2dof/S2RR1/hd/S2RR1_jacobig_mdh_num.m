% Gelenk-Jacobimatrix zu beliebigen Punkten eines Körpers für
% S2RR1
%
% Input:
% qJ [2x1]
%   Joint Angles [rad]
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: S2RR1_fkine_fixb_rotmat_mdh_sym_varpar.m
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
%
% Output:
% Jg_C [6x2]
%   geometric body jacobian for the defined point
%
% Quellen:
% [1] Ortmaier: Robotik I Skript

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 16:45
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover
% Moritz Schappler, schappler@irt.uni-hannover.de, 2016-03
% (C) Institut für Regelungstechnik, Leibniz Universität Hannover

function Jg_C = S2RR1_jacobig_mdh_num(qJ, link_index, r_i_i_C, pkin)
%% Init
%#codegen
%$cgargs {zeros(2,1),uint8(zeros(1,1)),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_jacobig_mdh_num: Joint angles qJ have to be [2x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
  'S2RR1_jacobig_mdh_num: link_index has to be [1x1] uint8');
assert(isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
  'S2RR1_jacobig_mdh_num: Position vector r_i_i_C has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_jacobig_mdh_num: Kinematic parameters pkin have to be [1x1] (double)');

if link_index > 3-1
  error('Index exceeds number of bodies');
end

% Initialisierung. Alle Spalten die nicht gesetzt werden haben keinen
% Einfluss. Fallunterscheidung für symbolische Eingabe.
if isa([qJ;pkin;r_i_i_C], 'double'), Jg_C = zeros(6,2);           % numerisch
else,                                Jg_C = sym('xx', [6,2]); end % symbolisch

if link_index == 0
  % Die Gelenkwinkel haben keinen Einfluss auf die Basis
  return;
end

%% Kinematik berechnen
T_c_mdh = S2RR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin);
[v_mdh, sigma_mdh] = S2RR1_structural_kinematic_parameters();
T_0_i = T_c_mdh(:,:,link_index+1);
R_0_i = T_0_i(1:3,1:3);
r_0_i_C = R_0_i * r_i_i_C;

j = link_index;
for tmp = 1:2
  % Vorgänger-Index
  k = v_mdh(j);

  % Drehachse des Gelenks, das diesen Körper bewegt ist die z-Achse dieses
  % Körpers (bei DH-Notation ist es der vorherige, hier MDH-Notation).
  ax = T_c_mdh(1:3,3,j+1);
  
  % Vektor vom Gelenk zum Punkt
  r_0_j_i = -T_c_mdh(1:3,4,j+1) + T_0_i(1:3,4);
  r_0_j_C = r_0_j_i + r_0_i_C;
  
  if sigma_mdh(j) == 0
    % Drehgelenk
    % [1], Gl. (4.19)
    % Hebelarm vom Gelenk zum Punkt
    jt = cross(ax, r_0_j_C);
    jr = ax;
  else
    % Schubgelenk, [1], Gl. (4.18)
    jt = ax;
    jr = zeros(3,1);
  end
  % Spalte der Jacobi-Matrix eintragen
  Jg_C(:,j) = [jt; jr];
  
  % Indizes tauschen: Kinematische Kette weiter Richtung Basis entlanggehen
  j = k;
  if j == 0
    % An Basis angekommen
    return;
  end
end

