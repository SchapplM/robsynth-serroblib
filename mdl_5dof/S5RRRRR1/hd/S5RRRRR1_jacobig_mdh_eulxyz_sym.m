% Jacobimatrix zu beliebigen Punkten eines Körpers für
% S5RRRRR1
%
% Input:
% phi_base [3x1]
%   Base orientation in world frame. Expressed with XYZ-Euler angles
% qJ [5x1]
%   Joint Angles [rad]
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
% Jg_C [6x(6+5)]
%   geometric body jacobian for the defined point with respect to
%   generalized coordinates:
%
% Quellen:
% [1] Aufzeichnungen Schappler (vom 10.06.2016)
% [2] BouyarmaneKhe2012: On the dynamics modeling of free-floating-base
%     articulated mechanisms and applications to humanoid whole-body dynamics
%     and control
%
% Siehe auch: S5RRRRR1_jacobig_mdh_eulxyz_num.m

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 14:53
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover
% Moritz Schappler, schappler@irt.uni-hannover.de, 2016-06
% (C) Institut für Regelungstechnik, Leibniz Universität Hannover

function Jg_C = S5RRRRR1_jacobig_mdh_eulxyz_sym(phi_base, qJ, link_index, r_i_i_C, pkin)
%% Init
%#codegen
%$cgargs {zeros(3,1),zeros(5,1),uint8(zeros(1,1)),zeros(3,1),zeros(6,1)}
assert(isreal(phi_base) && all(size(phi_base) == [3 1]), ...
  'S5RRRRR1_jacobig_mdh_eulxyz_sym: Base Euler-XYZ angles have to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_jacobig_mdh_eulxyz_sym: Joint angles qJ have to be [5x1] (double)');
assert(isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
  'S5RRRRR1_jacobig_mdh_eulxyz_sym: Position vector r_i_i_C has to be [3x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
  'S5RRRRR1_jacobig_mdh_eulxyz_sym: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_jacobig_mdh_eulxyz_sym: Kinematic parameters pkin have to be [6x1] (double)');
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

% Transformation zwischen Euler-XYZ-Winkel-Zeitableitung und
% Winkelgeschwindigkeit der Basis im Welt-KS
T_angvel = eulxyzjac(phi_base);
%% Basis-Jacobi-Matrix
% Siehe [1]
% [2], equ. (14)
Jg_BTC = [eye(3),     -skew(R_W_B*r_B_B_C)*T_angvel];
% [2], equ. (17)
Jg_BRC = [zeros(3,3), T_angvel];

%% Gelenk-Jacobi-Matrix
Jg_JC_B = S5RRRRR1_jacobig_floatb_twist_sym_varpar(qJ, link_index, r_i_i_C, pkin);
Jg_JTC = R_W_B * Jg_JC_B(1:3,:);
Jg_JRC = R_W_B * Jg_JC_B(4:6,:);

%% Gesamt-Jacobi-Matrix zusammensetzen
Jg_C = [Jg_BTC, Jg_JTC; ... % [2], equ. (14)
        Jg_BRC, Jg_JRC]; % [2], equ. (17)
