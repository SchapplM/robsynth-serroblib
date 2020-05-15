% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5RPPPR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:19:57
% EndTime: 2020-01-03 11:19:57
% DurationCPUTime: 0.11s
% Computational Cost: add. (131->47), mult. (83->48), div. (0->0), fcn. (125->10), ass. (0->31)
t15 = cos(pkin(8));
t11 = qJ(1) + pkin(7);
t5 = sin(t11);
t32 = t5 * t15;
t13 = sin(pkin(8));
t7 = cos(t11);
t31 = t7 * t13;
t30 = t7 * t15;
t12 = sin(pkin(9));
t29 = t12 * t15;
t14 = cos(pkin(9));
t28 = t14 * t15;
t27 = pkin(5) + 0;
t17 = sin(qJ(1));
t26 = t17 * pkin(1) + 0;
t25 = t5 * pkin(2) + t26;
t24 = -pkin(4) * t12 - qJ(3);
t8 = qJ(2) + t27;
t18 = cos(qJ(1));
t23 = -t18 * pkin(1) + 0;
t16 = -pkin(6) - qJ(4);
t3 = t14 * pkin(4) + pkin(3);
t22 = -t13 * t16 + t15 * t3;
t21 = pkin(3) * t15 + qJ(4) * t13;
t20 = -t7 * qJ(3) + t25;
t19 = -t5 * qJ(3) + t23;
t10 = pkin(9) + qJ(5);
t6 = cos(t10);
t4 = sin(t10);
t1 = t5 * t13;
t2 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; 0, 0, 1, t27; t17, t18, 0, 0; -t18, t17, 0, 0; 0, 0, 0, 1; 0, 0, 1, t8; t5, t7, 0, t26; -t7, t5, 0, t23; 0, 0, 0, 1; t13, t15, 0, t8; t32, -t1, -t7, t20; -t30, t31, -t5, -t7 * pkin(2) + t19; 0, 0, 0, 1; t13 * t14, -t13 * t12, -t15, t13 * pkin(3) - t15 * qJ(4) + t8; -t7 * t12 + t5 * t28, -t7 * t14 - t5 * t29, t1, t21 * t5 + t20; -t5 * t12 - t7 * t28, -t5 * t14 + t7 * t29, -t31, (-pkin(2) - t21) * t7 + t19; 0, 0, 0, 1; t13 * t6, -t13 * t4, -t15, t13 * t3 + t15 * t16 + t8; t6 * t32 - t7 * t4, -t4 * t32 - t7 * t6, t1, t22 * t5 + t24 * t7 + t25; -t6 * t30 - t5 * t4, t4 * t30 - t5 * t6, -t31, t24 * t5 + (-pkin(2) - t22) * t7 + t23; 0, 0, 0, 1;];
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
