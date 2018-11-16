% Jacobimatrix-Zeitableitung zu beliebigen Punkten eines Körpers für
% S5RRRRR1
% Floating-Base-Notation: Jacobi-Matrix bezogen auf Basis und Gelenke
%
% Input:
% phi_base [3x1]
%   Base orientation in world frame. Expressed with XYZ-Euler angles
% phiD_base [3x1]
%   Time Derivative of Base Orientation in world frame.
%   Expressed with XYZ-Euler angles ("rpy")
% qJ [5x1]
%   Joint Angles [rad]
% qJD [5x1]
%   Joint Velocities [rad/s]
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: S5RRRRR1_fkine_fixb_rotmat_mdh_sym_varpar.m
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
%
% Output:
% JgD_C [6x(6+5)]
%   time derivative of geometric body jacobian for the defined point with respect to
%   generalized coordinates:
%
% Quellen:
% [1] Aufzeichnungen Schappler (vom 22.06.2016)
% [2] BouyarmaneKhe2012: On the dynamics modeling of free-floating-base
%     articulated mechanisms and applications to humanoid whole-body dynamics
%     and control

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 14:53
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover
% Moritz Schappler, schappler@irt.uni-hannover.de, 2016-06
% (C) Institut für Regelungstechnik, Leibniz Universität Hannover

function JgD_C = S5RRRRR1_jacobigD_mdh_eulxyz_num(phi_base, phiD_base, qJ, qJD, link_index, r_i_i_C, ...
  pkin)
%% Init
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(5,1),zeros(5,1),
%$cgargs  uint8(zeros(1,1)),zeros(3,1),zeros(6,1)}
assert(isreal(phi_base) && all(size(phi_base) == [3 1]), ...
  'S5RRRRR1_jacobigD_mdh_eulxyz_num: Base Euler angles have to be [3x1] (double)');
assert(isreal(phiD_base) && all(size(phiD_base) == [3 1]), ...
  'S5RRRRR1_jacobigD_mdh_eulxyz_num: Base Euler angles time derivatives have to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_jacobigD_mdh_eulxyz_num: Joint angles qJ have to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_jacobigD_mdh_eulxyz_num: Joint velocities qJD have to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
  'S5RRRRR1_jacobigD_mdh_eulxyz_num: link_index has to be [1x1] uint8');
assert(isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
  'S5RRRRR1_jacobigD_mdh_eulxyz_num: Position vector r_i_i_C has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_jacobigD_mdh_eulxyz_num: Kinematic parameters pkin have to be [6x1] (double)');

if link_index > 6-1
  error('Index exceeds number of bodies');
end

%% Kinematik berechnen
T_c_mdh = S5RRRRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin);
v_mdh = S5RRRRR1_structural_kinematic_parameters();
T_B_i = T_c_mdh(:,:,link_index+1);
% Vektor von der Basis (B) zum Punkt C (auf Körper i)
r_B_B_C = T_B_i(1:3,4) + T_B_i(1:3,1:3) * r_i_i_C; % r_B_B_i + R_B_i_C

% Rotationsmatrix vom Welt- ins Basis-KS
R_W_B = eulxyz2r(phi_base);

% Transformation zwischen Euler-Winkel-Zeitableitung und
% Winkelgeschwindigkeit der Basis im Welt-KS
T_angvel = eulxyzjac(phi_base);
TD_angvel = eulxyzjacD(phi_base, phiD_base);

omega_W_B = T_angvel*phiD_base;

%% Gelenk-Jacobi-Matrix
Jg_JC_B = S5RRRRR1_jacobig_mdh_num(qJ, link_index, r_i_i_C, pkin);
JgD_JC_B = S5RRRRR1_jacobigD_mdh_num(qJ, qJD, link_index, r_i_i_C, pkin);
% Rotation ins Welt-KS und Zeitableitung der Rotation
JgD_JTC = skew(omega_W_B)*R_W_B*Jg_JC_B(1:3,:) + R_W_B*JgD_JC_B(1:3,:);
JgD_JRC = skew(omega_W_B)*R_W_B*Jg_JC_B(4:6,:) + R_W_B*JgD_JC_B(4:6,:);

%% Basis-Jacobi-Matrix
% Siehe [1]
% [2], equ. (14). Zeitableitung durch Produktregel und Eulersche
% Differentiationsregel
Term5_Teil1 = skew(Jg_JC_B(1:3,:)*qJD) * (R_W_B'); % erster Term aus [1] Gl. (5)
Term5_Teil2 = skew(r_B_B_C)*((skew(omega_W_B)*R_W_B)'); % zweiter Term aus [1] Gl. (5)
Term5 = Term5_Teil1 + Term5_Teil2;
Term4_Teil1 = skew(omega_W_B) * R_W_B * skew(r_B_B_C) * (R_W_B'); % erster Term aus [1] Gl. (4)
Term4 = Term4_Teil1 + R_W_B*Term5;
Term3_Teil2 = R_W_B*skew(r_B_B_C)*R_W_B'*TD_angvel; % zweiter Term aus [1] Gl. (3)
Term3 = Term4*T_angvel + Term3_Teil2;

JgD_BTC = [zeros(3,3),     -Term3];
% [2], equ. (17). Zeitableitung
JgD_BRC = [zeros(3,3), TD_angvel];

%% Gesamt-Jacobi-Matrix zusammensetzen
JgD_C = [JgD_BTC, JgD_JTC; ... % [2], equ. (14)
         JgD_BRC, JgD_JRC]; % [2], equ. (17)
