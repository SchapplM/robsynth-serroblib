% Calculate inertial parameters regressor of joint inertia matrix for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPPRR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t21 = cos(pkin(8));
t23 = sin(qJ(4));
t25 = cos(qJ(4));
t19 = sin(pkin(8));
t20 = cos(pkin(9));
t36 = t19 * t20;
t3 = t21 * t25 + t23 * t36;
t40 = t3 ^ 2;
t24 = cos(qJ(5));
t39 = 0.2e1 * t24;
t38 = t3 * t25;
t18 = sin(pkin(9));
t37 = t18 * t19;
t35 = t23 * t18;
t34 = t24 * t23;
t33 = t25 * t18;
t22 = sin(qJ(5));
t32 = t25 * t22;
t31 = t25 * t24;
t14 = t22 ^ 2;
t16 = t24 ^ 2;
t30 = t14 + t16;
t29 = t30 * t23;
t5 = -t21 * t23 + t25 * t36;
t1 = -t5 * t22 + t24 * t37;
t2 = t22 * t37 + t5 * t24;
t28 = -t1 * t22 + t2 * t24;
t6 = -t18 * t32 - t24 * t20;
t7 = t18 * t31 - t22 * t20;
t27 = -t6 * t22 + t7 * t24;
t17 = t25 ^ 2;
t15 = t23 ^ 2;
t13 = t21 ^ 2;
t12 = t20 ^ 2;
t11 = t19 ^ 2;
t10 = t18 ^ 2;
t9 = t15 * t10;
t8 = t10 * t11;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 + t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t11 + t13 + t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 ^ 2 + t40 + t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2 + t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t23 * t3 + t25 * t5 - t36) * t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 * t6 + t2 * t7 + t3 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 + t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t10 + t12 + t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 ^ 2 + t7 ^ 2 + t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t23 - t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t23 - t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t27 - t33) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 + t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 * t15 + t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t5, 0, 0, 0, 0, 0, 0, 0, 0, -t3 * t24, t3 * t22, t28, -t3 * pkin(4) + t28 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t33, 0, 0, 0, 0, 0, 0, 0, 0, -t18 * t34, t22 * t35, t27, -pkin(4) * t35 + t27 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t23, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t32, t29, t25 * pkin(4) + pkin(6) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t14, t22 * t39, 0, t16, 0, 0, pkin(4) * t39, -0.2e1 * pkin(4) * t22, 0.2e1 * t30 * pkin(6), t30 * pkin(6) ^ 2 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22 * t23, -t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, 0, t24, 0, -t22 * pkin(6), -t24 * pkin(6), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t4;
