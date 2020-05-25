% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5PRRRR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:09:26
% EndTime: 2019-12-05 17:09:26
% DurationCPUTime: 0.11s
% Computational Cost: add. (121->48), mult. (92->54), div. (0->0), fcn. (138->10), ass. (0->32)
t15 = qJ(2) + qJ(3);
t11 = cos(t15);
t16 = sin(pkin(9));
t35 = t16 * t11;
t18 = sin(qJ(4));
t34 = t16 * t18;
t20 = cos(qJ(4));
t33 = t16 * t20;
t17 = cos(pkin(9));
t32 = t17 * t11;
t31 = t17 * t18;
t30 = t17 * t20;
t21 = cos(qJ(2));
t7 = t21 * pkin(2) + pkin(1);
t29 = t17 * t7 + 0;
t13 = qJ(1) + 0;
t23 = -pkin(6) - pkin(5);
t28 = t16 * t7 + t17 * t23 + 0;
t19 = sin(qJ(2));
t27 = t19 * pkin(2) + t13;
t9 = sin(t15);
t26 = pkin(3) * t11 + pkin(7) * t9;
t22 = -pkin(8) - pkin(7);
t6 = t20 * pkin(4) + pkin(3);
t25 = t11 * t6 - t22 * t9;
t24 = -t16 * t23 + t29;
t14 = qJ(4) + qJ(5);
t10 = cos(t14);
t8 = sin(t14);
t4 = t17 * t9;
t3 = t16 * t9;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t17, -t16, 0, 0; t16, t17, 0, 0; 0, 0, 1, t13; 0, 0, 0, 1; t17 * t21, -t17 * t19, t16, t17 * pkin(1) + t16 * pkin(5) + 0; t16 * t21, -t16 * t19, -t17, t16 * pkin(1) - t17 * pkin(5) + 0; t19, t21, 0, t13; 0, 0, 0, 1; t32, -t4, t16, t24; t35, -t3, -t17, t28; t9, t11, 0, t27; 0, 0, 0, 1; t11 * t30 + t34, -t11 * t31 + t33, t4, t26 * t17 + t24; t11 * t33 - t31, -t11 * t34 - t30, t3, t26 * t16 + t28; t9 * t20, -t9 * t18, -t11, t9 * pkin(3) - t11 * pkin(7) + t27; 0, 0, 0, 1; t10 * t32 + t16 * t8, t16 * t10 - t8 * t32, t4, t25 * t17 + (pkin(4) * t18 - t23) * t16 + t29; t10 * t35 - t17 * t8, -t17 * t10 - t8 * t35, t3, -pkin(4) * t31 + t25 * t16 + t28; t9 * t10, -t9 * t8, -t11, t11 * t22 + t9 * t6 + t27; 0, 0, 0, 1;];
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
