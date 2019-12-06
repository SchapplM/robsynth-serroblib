% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPR6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t42 = sin(pkin(10));
t44 = cos(pkin(10));
t46 = sin(qJ(5));
t49 = cos(qJ(5));
t76 = -t46 * t42 + t49 * t44;
t45 = cos(pkin(5));
t47 = sin(qJ(3));
t50 = cos(qJ(3));
t43 = sin(pkin(5));
t48 = sin(qJ(2));
t68 = t43 * t48;
t17 = -t45 * t50 + t47 * t68;
t16 = t17 ^ 2;
t34 = -t44 * pkin(4) - pkin(3);
t75 = 0.2e1 * t34;
t74 = -0.2e1 * t50;
t40 = t47 ^ 2;
t73 = t40 * pkin(7);
t36 = t47 * pkin(7);
t72 = t17 * t47;
t71 = t42 * t44;
t70 = t42 * t47;
t69 = t42 * t50;
t51 = cos(qJ(2));
t67 = t43 * t51;
t66 = t44 * t47;
t65 = t44 * t50;
t63 = t47 * t50;
t61 = pkin(8) + qJ(4);
t26 = -t50 * pkin(3) - t47 * qJ(4) - pkin(2);
t12 = pkin(7) * t65 + t42 * t26;
t37 = t42 ^ 2;
t39 = t44 ^ 2;
t60 = t37 + t39;
t59 = 0.2e1 * t63;
t58 = t42 * t66;
t19 = t45 * t47 + t50 * t68;
t6 = -t19 * t42 - t44 * t67;
t7 = t19 * t44 - t42 * t67;
t57 = -t6 * t42 + t7 * t44;
t56 = -pkin(3) * t47 + qJ(4) * t50;
t21 = t44 * t26;
t11 = -pkin(7) * t69 + t21;
t55 = -t11 * t42 + t12 * t44;
t54 = t19 * t50 + t72;
t24 = t49 * t42 + t46 * t44;
t53 = pkin(7) ^ 2;
t41 = t50 ^ 2;
t38 = t43 ^ 2;
t35 = t40 * t53;
t31 = t38 * t51 ^ 2;
t28 = t61 * t44;
t27 = t61 * t42;
t25 = pkin(4) * t70 + t36;
t15 = t76 * t47;
t13 = t24 * t47;
t10 = -t46 * t27 + t49 * t28;
t9 = -t49 * t27 - t46 * t28;
t8 = -pkin(8) * t70 + t12;
t5 = -pkin(8) * t66 + t21 + (-pkin(7) * t42 - pkin(4)) * t50;
t4 = t46 * t6 + t49 * t7;
t3 = -t46 * t7 + t49 * t6;
t2 = t46 * t5 + t49 * t8;
t1 = -t46 * t8 + t49 * t5;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38 * t48 ^ 2 + t45 ^ 2 + t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 ^ 2 + t16 + t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6 ^ 2 + t7 ^ 2 + t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 + t4 ^ 2 + t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, -t68, 0, 0, 0, 0, 0, 0, 0, 0, t50 * t67, -t47 * t67, t54, pkin(2) * t67 + t54 * pkin(7), 0, 0, 0, 0, 0, 0, t17 * t70 - t6 * t50, t17 * t66 + t7 * t50, (-t42 * t7 - t44 * t6) * t47, pkin(7) * t72 + t6 * t11 + t7 * t12, 0, 0, 0, 0, 0, 0, t17 * t13 - t3 * t50, t17 * t15 + t4 * t50, -t4 * t13 - t3 * t15, t3 * t1 + t17 * t25 + t4 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t40, t59, 0, t41, 0, 0, 0.2e1 * pkin(2) * t50, -0.2e1 * pkin(2) * t47, 0.2e1 * (t40 + t41) * pkin(7), pkin(2) ^ 2 + t41 * t53 + t35, t39 * t40, -0.2e1 * t40 * t71, -0.2e1 * t44 * t63, t37 * t40, t42 * t59, t41, -0.2e1 * t11 * t50 + 0.2e1 * t42 * t73, 0.2e1 * t12 * t50 + 0.2e1 * t44 * t73, 0.2e1 * (-t11 * t44 - t12 * t42) * t47, t11 ^ 2 + t12 ^ 2 + t35, t15 ^ 2, -0.2e1 * t15 * t13, t15 * t74, t13 ^ 2, -t13 * t74, t41, -0.2e1 * t1 * t50 + 0.2e1 * t25 * t13, 0.2e1 * t25 * t15 + 0.2e1 * t2 * t50, -0.2e1 * t1 * t15 - 0.2e1 * t2 * t13, t1 ^ 2 + t2 ^ 2 + t25 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t19, 0, 0, 0, 0, 0, 0, 0, 0, -t17 * t44, t17 * t42, t57, -t17 * pkin(3) + t57 * qJ(4), 0, 0, 0, 0, 0, 0, -t17 * t76, t17 * t24, -t3 * t24 + t4 * t76, t4 * t10 + t17 * t34 + t3 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, t50, 0, -t36, -t50 * pkin(7), 0, 0, t58, (-t37 + t39) * t47, -t69, -t58, -t65, 0, -pkin(7) * t66 + t56 * t42, pkin(7) * t70 + t56 * t44, t55, -pkin(3) * t36 + t55 * qJ(4), t15 * t24, -t24 * t13 + t15 * t76, -t24 * t50, -t13 * t76, -t76 * t50, 0, t34 * t13 - t25 * t76 - t9 * t50, t10 * t50 + t34 * t15 + t25 * t24, -t1 * t24 - t10 * t13 - t9 * t15 + t2 * t76, t1 * t9 + t2 * t10 + t25 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t37, 0.2e1 * t71, 0, t39, 0, 0, 0.2e1 * pkin(3) * t44, -0.2e1 * pkin(3) * t42, 0.2e1 * t60 * qJ(4), t60 * qJ(4) ^ 2 + pkin(3) ^ 2, t24 ^ 2, 0.2e1 * t24 * t76, 0, t76 ^ 2, 0, 0, -t76 * t75, t24 * t75, 0.2e1 * t10 * t76 - 0.2e1 * t9 * t24, t10 ^ 2 + t34 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t66, 0, t36, 0, 0, 0, 0, 0, 0, t13, t15, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t42, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -t76, t24, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, 0, -t13, -t50, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, t76, 0, t9, -t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t14;
