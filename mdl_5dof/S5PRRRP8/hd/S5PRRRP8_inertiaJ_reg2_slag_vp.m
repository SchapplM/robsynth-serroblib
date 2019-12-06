% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRP8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t43 = sin(qJ(4));
t46 = cos(qJ(4));
t42 = cos(pkin(5));
t44 = sin(qJ(3));
t47 = cos(qJ(3));
t41 = sin(pkin(5));
t45 = sin(qJ(2));
t71 = t41 * t45;
t17 = t42 * t44 + t47 * t71;
t48 = cos(qJ(2));
t70 = t41 * t48;
t5 = t17 * t43 + t46 * t70;
t7 = t17 * t46 - t43 * t70;
t60 = t5 * t43 + t7 * t46;
t37 = t43 ^ 2;
t39 = t46 ^ 2;
t86 = t37 + t39;
t15 = -t42 * t47 + t44 * t71;
t14 = t15 ^ 2;
t55 = -t46 * pkin(4) - t43 * qJ(5);
t22 = -pkin(3) + t55;
t85 = -0.2e1 * t22;
t84 = -0.2e1 * t44;
t83 = 0.2e1 * t44;
t82 = pkin(3) * t43;
t81 = pkin(3) * t46;
t80 = pkin(7) * t43;
t38 = t44 ^ 2;
t79 = t38 * pkin(7);
t78 = t43 * pkin(8);
t77 = t44 * pkin(7);
t76 = t46 * pkin(8);
t74 = t15 * t43;
t73 = t15 * t44;
t72 = t15 * t46;
t69 = t43 * t44;
t68 = t43 * t46;
t67 = t43 * t47;
t66 = t44 * t47;
t32 = t46 * t44;
t65 = t46 * t47;
t23 = -t47 * pkin(3) - t44 * pkin(8) - pkin(2);
t12 = pkin(7) * t65 + t43 * t23;
t64 = t86 * pkin(8) ^ 2;
t63 = t47 * qJ(5);
t62 = t43 * t66;
t61 = t38 * t68;
t59 = t5 ^ 2 + t7 ^ 2 + t14;
t58 = t15 * t69 + t5 * t47;
t57 = t60 * pkin(8);
t8 = -t63 + t12;
t19 = t46 * t23;
t9 = -t19 + (pkin(4) + t80) * t47;
t56 = t9 * t43 + t8 * t46;
t54 = -pkin(4) * t43 + t46 * qJ(5);
t11 = -pkin(7) * t67 + t19;
t53 = -t11 * t43 + t12 * t46;
t52 = t17 * t47 + t73;
t51 = (-t43 * t7 + t46 * t5) * t44;
t50 = pkin(7) ^ 2;
t40 = t47 ^ 2;
t36 = t41 ^ 2;
t34 = t38 * t50;
t31 = t39 * t38;
t30 = t37 * t38;
t28 = t36 * t48 ^ 2;
t26 = pkin(8) * t67;
t25 = t43 * t32;
t24 = t65 * t84;
t21 = 0.2e1 * t86 * pkin(8);
t20 = (t37 - t39) * t44;
t13 = (pkin(7) - t54) * t44;
t1 = t15 * t32 + t7 * t47;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36 * t45 ^ 2 + t42 ^ 2 + t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17 ^ 2 + t14 + t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t71, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t70, -t44 * t70, t52, pkin(2) * t70 + t52 * pkin(7), 0, 0, 0, 0, 0, 0, t58, t1, t51, pkin(7) * t73 - t5 * t11 + t7 * t12, 0, 0, 0, 0, 0, 0, t58, t51, -t1, t15 * t13 + t5 * t9 + t7 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t38, 0.2e1 * t66, 0, t40, 0, 0, 0.2e1 * pkin(2) * t47, pkin(2) * t84, 0.2e1 * (t38 + t40) * pkin(7), pkin(2) ^ 2 + t40 * t50 + t34, t31, -0.2e1 * t61, t24, t30, 0.2e1 * t62, t40, -0.2e1 * t11 * t47 + 0.2e1 * t43 * t79, 0.2e1 * t12 * t47 + 0.2e1 * t46 * t79, (-t11 * t46 - t12 * t43) * t83, t11 ^ 2 + t12 ^ 2 + t34, t31, t24, 0.2e1 * t61, t40, -0.2e1 * t62, t30, 0.2e1 * t13 * t69 + 0.2e1 * t9 * t47, (-t43 * t8 + t46 * t9) * t83, -0.2e1 * t13 * t32 - 0.2e1 * t8 * t47, t13 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t17, 0, 0, 0, 0, 0, 0, 0, 0, -t72, t74, t60, -t15 * pkin(3) + t57, 0, 0, 0, 0, 0, 0, -t72, t60, -t74, t15 * t22 + t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, t47, 0, -t77, -t47 * pkin(7), 0, 0, t25, -t20, -t67, -t25, -t65, 0, t26 + (-pkin(7) * t46 - t82) * t44, pkin(8) * t65 + (t80 - t81) * t44, t53, -pkin(3) * t77 + t53 * pkin(8), t25, -t67, t20, 0, t65, -t25, -t13 * t46 + t22 * t69 + t26, t56, -t13 * t43 + (-pkin(8) * t47 - t22 * t44) * t46, t56 * pkin(8) + t13 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t37, 0.2e1 * t68, 0, t39, 0, 0, 0.2e1 * t81, -0.2e1 * t82, t21, pkin(3) ^ 2 + t64, t37, 0, -0.2e1 * t68, 0, 0, t39, t46 * t85, t21, t43 * t85, t22 ^ 2 + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t7, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, t7, -t5 * pkin(4) + t7 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, -t69, -t47, t11, -t12, 0, 0, 0, t32, 0, -t47, t69, 0, t19 + (-0.2e1 * pkin(4) - t80) * t47, t55 * t44, -0.2e1 * t63 + t12, -t9 * pkin(4) + t8 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, t46, 0, -t78, -t76, 0, 0, 0, t43, 0, 0, -t46, 0, -t78, t54, t76, t54 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, t32, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t2;
