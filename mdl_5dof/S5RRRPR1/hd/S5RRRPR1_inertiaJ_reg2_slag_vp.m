% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t49 = sin(qJ(3));
t50 = sin(qJ(2));
t52 = cos(qJ(3));
t53 = cos(qJ(2));
t30 = t49 * t50 - t52 * t53;
t32 = t49 * t53 + t50 * t52;
t46 = sin(pkin(9));
t47 = cos(pkin(9));
t17 = t47 * t30 + t32 * t46;
t41 = -pkin(2) * t53 - pkin(1);
t23 = pkin(3) * t30 + t41;
t10 = pkin(4) * t17 + t23;
t66 = 0.2e1 * t10;
t65 = 0.2e1 * t23;
t64 = 0.2e1 * t41;
t63 = 0.2e1 * t53;
t62 = -pkin(7) - pkin(6);
t61 = t46 * pkin(3);
t60 = t49 * pkin(2);
t44 = t50 ^ 2;
t45 = t53 ^ 2;
t59 = t44 + t45;
t58 = t47 * t60;
t43 = t52 * pkin(2);
t40 = t43 + pkin(3);
t27 = t40 * t46 + t58;
t57 = -t27 - t61;
t36 = t62 * t50;
t37 = t62 * t53;
t20 = t52 * t36 + t37 * t49;
t14 = -qJ(4) * t32 + t20;
t21 = t36 * t49 - t37 * t52;
t15 = -qJ(4) * t30 + t21;
t5 = t47 * t14 - t15 * t46;
t25 = t47 * t40 - t46 * t60;
t6 = t14 * t46 + t15 * t47;
t51 = cos(qJ(5));
t48 = sin(qJ(5));
t42 = t47 * pkin(3);
t38 = t42 + pkin(4);
t35 = t51 * t38;
t28 = t48 * t38 + t51 * t61;
t26 = -t48 * t61 + t35;
t24 = pkin(4) + t25;
t22 = t51 * t24;
t19 = -t30 * t46 + t32 * t47;
t13 = t48 * t24 + t51 * t27;
t12 = -t27 * t48 + t22;
t9 = -t17 * t48 + t19 * t51;
t7 = t51 * t17 + t19 * t48;
t4 = -pkin(8) * t17 + t6;
t3 = -pkin(8) * t19 + t5;
t2 = t3 * t48 + t4 * t51;
t1 = t3 * t51 - t4 * t48;
t8 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t44, t50 * t63, 0, t45, 0, 0, pkin(1) * t63, -0.2e1 * pkin(1) * t50, 0.2e1 * t59 * pkin(6), pkin(6) ^ 2 * t59 + pkin(1) ^ 2, t32 ^ 2, -0.2e1 * t32 * t30, 0, t30 ^ 2, 0, 0, t30 * t64, t32 * t64, -0.2e1 * t20 * t32 - 0.2e1 * t21 * t30, t20 ^ 2 + t21 ^ 2 + t41 ^ 2, t19 ^ 2, -0.2e1 * t19 * t17, 0, t17 ^ 2, 0, 0, t17 * t65, t19 * t65, -0.2e1 * t17 * t6 - 0.2e1 * t19 * t5, t23 ^ 2 + t5 ^ 2 + t6 ^ 2, t9 ^ 2, -0.2e1 * t9 * t7, 0, t7 ^ 2, 0, 0, t7 * t66, t9 * t66, -0.2e1 * t1 * t9 - 0.2e1 * t2 * t7, t1 ^ 2 + t10 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, t53, 0, -t50 * pkin(6), -t53 * pkin(6), 0, 0, 0, 0, t32, 0, -t30, 0, t20, -t21, (-t30 * t49 - t32 * t52) * pkin(2), (t20 * t52 + t21 * t49) * pkin(2), 0, 0, t19, 0, -t17, 0, t5, -t6, -t17 * t27 - t19 * t25, t25 * t5 + t27 * t6, 0, 0, t9, 0, -t7, 0, t1, -t2, -t12 * t9 - t13 * t7, t1 * t12 + t13 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t43, -0.2e1 * t60, 0, (t49 ^ 2 + t52 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t25, -0.2e1 * t27, 0, t25 ^ 2 + t27 ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t12, -0.2e1 * t13, 0, t12 ^ 2 + t13 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, -t30, 0, t20, -t21, 0, 0, 0, 0, t19, 0, -t17, 0, t5, -t6, (-t17 * t46 - t19 * t47) * pkin(3), (t46 * t6 + t47 * t5) * pkin(3), 0, 0, t9, 0, -t7, 0, t1, -t2, -t26 * t9 - t28 * t7, t1 * t26 + t2 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t43, -t60, 0, 0, 0, 0, 0, 0, 0, 1, t25 + t42, -t58 + (-pkin(3) - t40) * t46, 0, (t25 * t47 + t27 * t46) * pkin(3), 0, 0, 0, 0, 0, 1, t48 * t57 + t22 + t35, t57 * t51 + (-t24 - t38) * t48, 0, t12 * t26 + t13 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t42, -0.2e1 * t61, 0, (t46 ^ 2 + t47 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t26, -0.2e1 * t28, 0, t26 ^ 2 + t28 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t19, 0, t23, 0, 0, 0, 0, 0, 0, t7, t9, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, -t7, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t12, -t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t26, -t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t8;
