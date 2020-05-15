% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5PRPRR5_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:17
% EndTime: 2019-12-05 15:53:17
% DurationCPUTime: 0.14s
% Computational Cost: add. (120->55), mult. (117->64), div. (0->0), fcn. (168->10), ass. (0->31)
t18 = sin(pkin(9));
t34 = t18 * pkin(3);
t20 = cos(pkin(9));
t7 = t20 * pkin(3) + pkin(2);
t19 = sin(pkin(8));
t33 = t19 * t18;
t24 = cos(qJ(2));
t32 = t19 * t24;
t21 = cos(pkin(8));
t31 = t21 * t24;
t22 = -pkin(6) - qJ(3);
t17 = pkin(9) + qJ(4);
t30 = t19 * pkin(1) + 0;
t15 = qJ(1) + 0;
t29 = t21 * pkin(1) + t19 * pkin(5) + 0;
t9 = cos(t17);
t1 = pkin(4) * t9 + t7;
t16 = -pkin(7) + t22;
t23 = sin(qJ(2));
t28 = t1 * t24 - t16 * t23;
t27 = -t22 * t23 + t24 * t7;
t26 = pkin(2) * t24 + qJ(3) * t23;
t25 = -t21 * pkin(5) + t30;
t10 = qJ(5) + t17;
t8 = sin(t17);
t6 = cos(t10);
t5 = sin(t10);
t4 = t21 * t23;
t3 = t19 * t23;
t2 = pkin(4) * t8 + t34;
t11 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t21, -t19, 0, 0; t19, t21, 0, 0; 0, 0, 1, t15; 0, 0, 0, 1; t31, -t4, t19, t29; t32, -t3, -t21, t25; t23, t24, 0, t15; 0, 0, 0, 1; t20 * t31 + t33, -t18 * t31 + t19 * t20, t4, t26 * t21 + t29; -t21 * t18 + t20 * t32, -t18 * t32 - t21 * t20, t3, t26 * t19 + t25; t23 * t20, -t23 * t18, -t24, t23 * pkin(2) - t24 * qJ(3) + t15; 0, 0, 0, 1; t19 * t8 + t9 * t31, t19 * t9 - t8 * t31, t4, pkin(3) * t33 + t27 * t21 + t29; -t21 * t8 + t9 * t32, -t21 * t9 - t8 * t32, t3, (-pkin(5) - t34) * t21 + t27 * t19 + t30; t23 * t9, -t23 * t8, -t24, t24 * t22 + t23 * t7 + t15; 0, 0, 0, 1; t19 * t5 + t6 * t31, t19 * t6 - t5 * t31, t4, t19 * t2 + t28 * t21 + t29; -t21 * t5 + t6 * t32, -t21 * t6 - t5 * t32, t3, (-pkin(5) - t2) * t21 + t28 * t19 + t30; t23 * t6, -t23 * t5, -t24, t23 * t1 + t24 * t16 + t15; 0, 0, 0, 1;];
T_ges = t11;
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
