% Zeitableitung der Gelenk-Jacobimatrix zu beliebigen Punkten eines Körpers für
% S4PRRR1
%
% Input:
% qJ [1x4]
%   Joint Angles [rad]
% qJD [1x4]
%   Joint Velocities [rad/s]
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: S4PRRR1_fkine_fixb_rotmat_mdh_sym_varpar.m
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
%
% Output:
% JgD_C [6x4]
%   time derivative of geometric body jacobian for the defined point
%
% Siehe auch: S4PRRR1_jacobig_mdh_num.m
%
% Quelle:
% Berechnungen Moritz Schappler, 21.06.2016

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:45
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover
% Moritz Schappler, schappler@irt.uni-hannover.de, 2016-06
% (C) Institut für Regelungstechnik, Leibniz Universität Hannover

function JgD_C = S4PRRR1_jacobigD_mdh_num(qJ, qJD, link_index, r_i_i_C, pkin)
%% Init
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(zeros(1,1)),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_jacobigD_mdh_num: Joint angles qJ have to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_jacobigD_mdh_num: Joint velocities qJD have to be [4x1] (double)');
assert(isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
  'S4PRRR1_jacobigD_mdh_num: Position vector r_i_i_C has to be [3x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
  'S4PRRR1_jacobigD_mdh_num: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_jacobigD_mdh_num: Kinematic parameters pkin have to be [7x1] (double)');

% Initialisierung. Alle Spalten die nicht gesetzt werden haben keinen
% Einfluss.
JgD_C = zeros(6,4);

if link_index == 0
  % Die Gelenkwinkel haben keinen Einfluss auf die Basis
  return;
end

%% Kinematik berechnen
% direkte Kinematik
T_c_mdh = S4PRRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin);
v_mdh = S4PRRR1_structural_kinematic_parameters();
T_0_i = T_c_mdh(:,:,link_index+1);
R_0_i = T_0_i(1:3,1:3);
r_0_i_C = R_0_i * r_i_i_C;

% Geschwindigkeit des Punktes C
Jg_C = S4PRRR1_jacobig_mdh_num(qJ, link_index, r_i_i_C, pkin);
V_C = Jg_C*qJD;
rD_0_0_C = V_C(1:3);

j = link_index;
for tmp = 1:4 % Schleife mit Dummy-Variable begrenzt maximale Anzahl an Durchgängen
  % Vorgänger-Index
  k = v_mdh(j);

  % Geschwindigkeit des Körpers j
  Jg_j = S4PRRR1_jacobig_mdh_num(qJ, j, zeros(3,1), pkin);
  V_j = Jg_j*qJD;
  rD_0_0_j = V_j(1:3);
  omega_j = V_j(4:6);

  %  Drehachse des Gelenks, das diesen Körper bewegt ist die z-Achse dieses
  %  Körpers (bei DH-Notation ist es der vorherige, hier MDH-Notation).
  ax = T_c_mdh(1:3,3,j+1);
  jrD = cross(omega_j, ax); % Zeitableitung von ax

  % Vektor vom Gelenk zum Punkt
  r_0_j_i = -T_c_mdh(1:3,4,j+1) + T_0_i(1:3,4);
  r_0_j_C = r_0_j_i + r_0_i_C;

  % Hebelarm vom Gelenk zum Punkt
  rD_0_j_C = -rD_0_0_j + rD_0_0_C;
  jtD = cross(cross(omega_j, ax), r_0_j_C) + ... % Zeitableitung von ax
    cross(ax, rD_0_j_C); % Zeitableitung von r_0_j_C

  % Spalte der Jacobi-Matrix eintragen
  JgD_C(:,j) = [jtD; jrD];

  % Indizes tauschen: Kinematische Kette weiter Richtung Basis entlanggehen
  j = k;

  if j == 0
    % An Basis angekommen
    return;
  end
end