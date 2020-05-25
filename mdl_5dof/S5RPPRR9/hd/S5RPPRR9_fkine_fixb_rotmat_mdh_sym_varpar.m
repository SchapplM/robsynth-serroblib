% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% T_c_mdh [4x4x(5+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   6:  mdh base (link 0) -> mdh frame (6-1), link (6-1)
%   ...
%   5+1:  mdh base (link 0) -> mdh frame (5)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPPRR9_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR9_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR9_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:02:13
% EndTime: 2019-12-31 18:02:13
% DurationCPUTime: 0.09s
% Computational Cost: add. (87->36), mult. (124->36), div. (0->0), fcn. (189->8), ass. (0->24)
t16 = sin(qJ(4));
t26 = sin(pkin(8));
t27 = cos(pkin(8));
t30 = sin(qJ(1));
t31 = cos(qJ(1));
t3 = -t30 * t26 - t31 * t27;
t33 = t3 * t16;
t4 = t31 * t26 - t30 * t27;
t32 = t4 * t16;
t15 = sin(qJ(5));
t18 = cos(qJ(4));
t29 = t15 * t18;
t17 = cos(qJ(5));
t28 = t17 * t18;
t14 = pkin(5) + 0;
t25 = t31 * pkin(1) + t30 * qJ(2) + 0;
t8 = -qJ(3) + t14;
t24 = t31 * pkin(2) + t25;
t23 = -pkin(4) * t18 - pkin(7) * t16;
t22 = t30 * pkin(1) - t31 * qJ(2) + 0;
t21 = -t3 * pkin(3) + t4 * pkin(6) + t24;
t20 = t30 * pkin(2) + t22;
t19 = -t4 * pkin(3) - t3 * pkin(6) + t20;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t31, -t30, 0, 0; t30, t31, 0, 0; 0, 0, 1, t14; 0, 0, 0, 1; t31, 0, t30, t25; t30, 0, -t31, t22; 0, 1, 0, t14; 0, 0, 0, 1; -t3, -t4, 0, t24; -t4, t3, 0, t20; 0, 0, -1, t8; 0, 0, 0, 1; -t3 * t18, t33, t4, t21; -t4 * t18, t32, -t3, t19; -t16, -t18, 0, t8; 0, 0, 0, 1; t4 * t15 - t3 * t28, t4 * t17 + t3 * t29, -t33, t23 * t3 + t21; -t3 * t15 - t4 * t28, -t3 * t17 + t4 * t29, -t32, t23 * t4 + t19; -t16 * t17, t16 * t15, t18, -t16 * pkin(4) + t18 * pkin(7) + t8; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,5+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,5+1]); end % symbolisch
for i = 1:5+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
