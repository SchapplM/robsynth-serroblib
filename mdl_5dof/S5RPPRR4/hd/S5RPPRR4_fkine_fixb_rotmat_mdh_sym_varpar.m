% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPPRR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:30:16
% EndTime: 2020-01-03 11:30:16
% DurationCPUTime: 0.15s
% Computational Cost: add. (124->61), mult. (117->64), div. (0->0), fcn. (168->10), ass. (0->30)
t15 = sin(pkin(9));
t27 = -t15 * pkin(3) - qJ(2);
t17 = cos(pkin(9));
t6 = t17 * pkin(3) + pkin(2);
t18 = cos(pkin(8));
t20 = sin(qJ(1));
t32 = t20 * t18;
t16 = sin(pkin(8));
t21 = cos(qJ(1));
t31 = t21 * t16;
t30 = t21 * t18;
t19 = -pkin(6) - qJ(3);
t13 = pkin(9) + qJ(4);
t7 = sin(t13);
t29 = -pkin(4) * t7 + t27;
t14 = pkin(5) + 0;
t28 = t20 * pkin(1) + 0;
t26 = -t20 * qJ(2) + 0;
t8 = cos(t13);
t1 = pkin(4) * t8 + t6;
t12 = -pkin(7) + t19;
t25 = t1 * t18 - t12 * t16;
t24 = -t16 * t19 + t18 * t6;
t23 = pkin(2) * t18 + qJ(3) * t16;
t22 = -t21 * qJ(2) + t28;
t9 = qJ(5) + t13;
t5 = cos(t9);
t4 = sin(t9);
t3 = t20 * t16;
t2 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t14; t20, t21, 0, 0; -t21, t20, 0, 0; 0, 0, 0, 1; t16, t18, 0, t14; t32, -t3, -t21, t22; -t30, t31, -t20, -t21 * pkin(1) + t26; 0, 0, 0, 1; t16 * t17, -t16 * t15, -t18, t16 * pkin(2) - t18 * qJ(3) + t14; -t21 * t15 + t17 * t32, -t15 * t32 - t21 * t17, t3, t23 * t20 + t22; -t20 * t15 - t17 * t30, t15 * t30 - t20 * t17, -t31, (-pkin(1) - t23) * t21 + t26; 0, 0, 0, 1; t16 * t8, -t16 * t7, -t18, t16 * t6 + t18 * t19 + t14; -t21 * t7 + t8 * t32, -t21 * t8 - t7 * t32, t3, t24 * t20 + t27 * t21 + t28; -t20 * t7 - t8 * t30, -t20 * t8 + t7 * t30, -t31, 0 + t27 * t20 + (-pkin(1) - t24) * t21; 0, 0, 0, 1; t16 * t5, -t16 * t4, -t18, t16 * t1 + t18 * t12 + t14; -t21 * t4 + t5 * t32, -t21 * t5 - t4 * t32, t3, t25 * t20 + t29 * t21 + t28; -t20 * t4 - t5 * t30, -t20 * t5 + t4 * t30, -t31, 0 + t29 * t20 + (-pkin(1) - t25) * t21; 0, 0, 0, 1;];
T_ges = t2;
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
