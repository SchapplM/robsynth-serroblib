% Calculate inertial parameters regressor of joint inertia matrix for
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
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPPR1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t32 = sin(pkin(8));
t56 = -0.2e1 * t32;
t31 = sin(pkin(9));
t34 = cos(pkin(9));
t37 = sin(qJ(5));
t38 = cos(qJ(5));
t16 = -t37 * t31 + t38 * t34;
t55 = 0.2e1 * t32;
t35 = cos(pkin(8));
t54 = -0.2e1 * t35;
t33 = sin(pkin(7));
t53 = t33 * pkin(1);
t36 = cos(pkin(7));
t52 = t36 * pkin(1);
t25 = qJ(3) + t53;
t51 = t25 * t31;
t28 = t32 ^ 2;
t50 = t28 * t34;
t49 = t31 * t32;
t48 = t31 * t35;
t47 = t34 * t32;
t46 = t34 * t35;
t26 = -pkin(2) - t52;
t15 = -t35 * pkin(3) - t32 * qJ(4) + t26;
t6 = t31 * t15 + t25 * t46;
t27 = t31 ^ 2;
t29 = t34 ^ 2;
t43 = t27 + t29;
t30 = t35 ^ 2;
t42 = t28 + t30;
t41 = t35 * t55;
t13 = t34 * t15;
t5 = -t25 * t48 + t13;
t40 = t6 * t31 + t5 * t34;
t17 = t38 * t31 + t37 * t34;
t24 = t29 * t28;
t23 = t27 * t28;
t22 = t25 ^ 2;
t20 = t32 * t25;
t19 = t28 * t22;
t14 = pkin(4) * t49 + t20;
t11 = t16 * t32;
t9 = t17 * t32;
t8 = t11 ^ 2;
t7 = t9 ^ 2;
t4 = -pkin(6) * t49 + t6;
t3 = -pkin(6) * t47 + t13 + (-pkin(4) - t51) * t35;
t2 = t37 * t3 + t38 * t4;
t1 = t38 * t3 - t37 * t4;
t10 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t52, -0.2e1 * t53, 0, (t33 ^ 2 + t36 ^ 2) * pkin(1) ^ 2, t28, t41, 0, t30, 0, 0, t26 * t54, t26 * t55, 0.2e1 * t42 * t25, t30 * t22 + t26 ^ 2 + t19, t24, -0.2e1 * t31 * t50, t46 * t56, t23, t31 * t41, t30, 0.2e1 * t28 * t51 - 0.2e1 * t5 * t35, 0.2e1 * t25 * t50 + 0.2e1 * t6 * t35, t40 * t56, t5 ^ 2 + t6 ^ 2 + t19, t8, -0.2e1 * t11 * t9, t11 * t54, t7, -t9 * t54, t30, -0.2e1 * t1 * t35 + 0.2e1 * t14 * t9, 0.2e1 * t14 * t11 + 0.2e1 * t2 * t35, -0.2e1 * t1 * t11 - 0.2e1 * t2 * t9, t1 ^ 2 + t14 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t25 * t35 - t31 * t5 + t34 * t6) * t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t9 + t2 * t11 - t14 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24 + t23 + t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 + t7 + t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, t32, 0, t26, 0, 0, 0, 0, 0, 0, -t46, t48, -t43 * t32, t40, 0, 0, 0, 0, 0, 0, -t16 * t35, t17 * t35, -t16 * t11 - t17 * t9, t1 * t16 + t2 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t17 - t9 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 ^ 2 + t17 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t47, 0, t20, 0, 0, 0, 0, 0, 0, t9, t11, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, -t9, -t35, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t17, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t10;
