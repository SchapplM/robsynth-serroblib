% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S5PRRRR7
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
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S5PRRRR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:11:38
% EndTime: 2019-12-05 17:11:38
% DurationCPUTime: 0.14s
% Computational Cost: add. (120->55), mult. (117->66), div. (0->0), fcn. (168->10), ass. (0->33)
t24 = -pkin(7) - pkin(6);
t20 = sin(qJ(3));
t36 = t20 * pkin(3);
t22 = cos(qJ(3));
t7 = t22 * pkin(3) + pkin(2);
t18 = sin(pkin(9));
t35 = t18 * t20;
t23 = cos(qJ(2));
t34 = t18 * t23;
t19 = cos(pkin(9));
t33 = t19 * t23;
t32 = t20 * t23;
t31 = t22 * t23;
t17 = qJ(3) + qJ(4);
t30 = t18 * pkin(1) + 0;
t15 = qJ(1) + 0;
t29 = t19 * pkin(1) + t18 * pkin(5) + 0;
t21 = sin(qJ(2));
t28 = pkin(2) * t23 + pkin(6) * t21;
t9 = cos(t17);
t1 = pkin(4) * t9 + t7;
t16 = -pkin(8) + t24;
t27 = t1 * t23 - t16 * t21;
t26 = -t21 * t24 + t23 * t7;
t25 = -t19 * pkin(5) + t30;
t13 = qJ(5) + t17;
t8 = sin(t17);
t6 = cos(t13);
t5 = sin(t13);
t4 = t19 * t21;
t3 = t18 * t21;
t2 = pkin(4) * t8 + t36;
t10 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t19, -t18, 0, 0; t18, t19, 0, 0; 0, 0, 1, t15; 0, 0, 0, 1; t33, -t4, t18, t29; t34, -t3, -t19, t25; t21, t23, 0, t15; 0, 0, 0, 1; t19 * t31 + t35, t18 * t22 - t19 * t32, t4, t28 * t19 + t29; t18 * t31 - t19 * t20, -t18 * t32 - t19 * t22, t3, t28 * t18 + t25; t21 * t22, -t21 * t20, -t23, t21 * pkin(2) - t23 * pkin(6) + t15; 0, 0, 0, 1; t18 * t8 + t9 * t33, t18 * t9 - t8 * t33, t4, pkin(3) * t35 + t26 * t19 + t29; -t19 * t8 + t9 * t34, -t19 * t9 - t8 * t34, t3, (-pkin(5) - t36) * t19 + t26 * t18 + t30; t21 * t9, -t21 * t8, -t23, t21 * t7 + t23 * t24 + t15; 0, 0, 0, 1; t18 * t5 + t6 * t33, t18 * t6 - t5 * t33, t4, t18 * t2 + t27 * t19 + t29; -t19 * t5 + t6 * t34, -t19 * t6 - t5 * t34, t3, (-pkin(5) - t2) * t19 + t27 * t18 + t30; t21 * t6, -t21 * t5, -t23, t21 * t1 + t23 * t16 + t15; 0, 0, 0, 1;];
T_ges = t10;
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
