% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPPRP1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:00
% EndTime: 2020-01-03 11:25:00
% DurationCPUTime: 0.11s
% Computational Cost: add. (121->42), mult. (83->38), div. (0->0), fcn. (125->8), ass. (0->32)
t13 = qJ(1) + pkin(7);
t10 = cos(t13);
t14 = sin(pkin(8));
t33 = t10 * t14;
t17 = sin(qJ(4));
t32 = t14 * t17;
t15 = cos(pkin(8));
t31 = t15 * t17;
t19 = cos(qJ(4));
t30 = t15 * t19;
t29 = pkin(5) + 0;
t18 = sin(qJ(1));
t28 = t18 * pkin(1) + 0;
t9 = sin(t13);
t27 = t9 * pkin(2) + t28;
t26 = -pkin(4) * t17 - qJ(3);
t11 = qJ(2) + t29;
t20 = cos(qJ(1));
t25 = -t20 * pkin(1) + 0;
t24 = pkin(3) * t15 + pkin(6) * t14;
t16 = -qJ(5) - pkin(6);
t8 = t19 * pkin(4) + pkin(3);
t23 = -t14 * t16 + t15 * t8;
t22 = -t10 * qJ(3) + t27;
t21 = -t9 * qJ(3) + t25;
t7 = t14 * t19;
t5 = t9 * t14;
t4 = -t10 * t30 - t9 * t17;
t3 = t10 * t31 - t9 * t19;
t2 = -t10 * t17 + t9 * t30;
t1 = -t10 * t19 - t9 * t31;
t6 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t29; t18, t20, 0, 0; -t20, t18, 0, 0; 0, 0, 0, 1; 0, 0, 1, t11; t9, t10, 0, t28; -t10, t9, 0, t25; 0, 0, 0, 1; t14, t15, 0, t11; t9 * t15, -t5, -t10, t22; -t10 * t15, t33, -t9, -t10 * pkin(2) + t21; 0, 0, 0, 1; t7, -t32, -t15, t14 * pkin(3) - t15 * pkin(6) + t11; t2, t1, t5, t24 * t9 + t22; t4, t3, -t33, (-pkin(2) - t24) * t10 + t21; 0, 0, 0, 1; t7, -t32, -t15, t14 * t8 + t15 * t16 + t11; t2, t1, t5, t26 * t10 + t23 * t9 + t27; t4, t3, -t33, t26 * t9 + (-pkin(2) - t23) * t10 + t25; 0, 0, 0, 1;];
T_ges = t6;
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
