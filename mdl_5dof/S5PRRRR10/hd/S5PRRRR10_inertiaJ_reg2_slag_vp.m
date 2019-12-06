% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR10_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t46 = sin(pkin(5));
t48 = cos(pkin(5));
t51 = sin(qJ(3));
t52 = sin(qJ(2));
t55 = cos(qJ(3));
t47 = cos(pkin(6));
t56 = cos(qJ(2));
t79 = t47 * t56;
t45 = sin(pkin(6));
t84 = t45 * t51;
t15 = t48 * t84 + (t51 * t79 + t52 * t55) * t46;
t81 = t46 * t56;
t25 = -t45 * t81 + t47 * t48;
t50 = sin(qJ(4));
t54 = cos(qJ(4));
t5 = t15 * t50 - t25 * t54;
t104 = t5 ^ 2;
t82 = t46 * t52;
t83 = t45 * t55;
t13 = -t46 * t55 * t79 - t48 * t83 + t51 * t82;
t103 = t13 ^ 2;
t26 = -t54 * t47 + t50 * t84;
t102 = t26 ^ 2;
t28 = t47 * t50 + t54 * t84;
t49 = sin(qJ(5));
t53 = cos(qJ(5));
t16 = t28 * t49 + t53 * t83;
t101 = -0.2e1 * t16;
t100 = -0.2e1 * t28;
t99 = 0.2e1 * t45;
t98 = -0.2e1 * t50;
t97 = 0.2e1 * t54;
t96 = pkin(2) * t51;
t95 = pkin(2) * t55;
t94 = pkin(4) * t53;
t42 = t50 ^ 2;
t93 = t42 * pkin(9);
t92 = t50 * pkin(9);
t67 = pkin(8) * t83;
t23 = t67 + (pkin(9) + t96) * t47;
t24 = (-pkin(3) * t55 - pkin(9) * t51 - pkin(2)) * t45;
t11 = -t23 * t50 + t24 * t54;
t9 = pkin(4) * t83 - t11;
t91 = t9 * t49;
t90 = t9 * t53;
t89 = t16 * t53;
t18 = t28 * t53 - t49 * t83;
t88 = t18 * t49;
t87 = t26 * t54;
t86 = t28 * t50;
t39 = t45 ^ 2;
t85 = t39 * t55;
t80 = t47 * t51;
t78 = t49 * t26;
t77 = t49 * t50;
t76 = t49 * t53;
t75 = t49 * t54;
t74 = t50 * t26;
t73 = t53 * t26;
t72 = t53 * t50;
t71 = t53 * t54;
t41 = t49 ^ 2;
t43 = t53 ^ 2;
t70 = t41 + t43;
t69 = 0.2e1 * t83;
t68 = t50 * t97;
t66 = t50 * t83;
t65 = t49 * t72;
t64 = t54 * t83;
t12 = t23 * t54 + t24 * t50;
t10 = -pkin(10) * t83 + t12;
t35 = pkin(8) * t84;
t22 = t35 + (-pkin(3) - t95) * t47;
t8 = t26 * pkin(4) - t28 * pkin(10) + t22;
t1 = -t10 * t49 + t53 * t8;
t2 = t10 * t53 + t49 * t8;
t63 = -t1 * t49 + t2 * t53;
t7 = t15 * t54 + t25 * t50;
t3 = t13 * t53 - t49 * t7;
t4 = t13 * t49 + t53 * t7;
t62 = -t3 * t49 + t4 * t53;
t61 = t5 * t50 + t54 * t7;
t60 = -t11 * t50 + t12 * t54;
t32 = -pkin(4) * t54 - pkin(10) * t50 - pkin(3);
t19 = -pkin(9) * t75 + t32 * t53;
t20 = pkin(9) * t71 + t32 * t49;
t59 = -t19 * t49 + t20 * t53;
t58 = pkin(9) ^ 2;
t44 = t54 ^ 2;
t38 = t42 * t58;
t36 = t39 * t55 ^ 2;
t30 = pkin(2) * t80 + t67;
t29 = t47 * t95 - t35;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 ^ 2 + (t52 ^ 2 + t56 ^ 2) * t46 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 ^ 2 + t25 ^ 2 + t103, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t103 + t104, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 + t4 ^ 2 + t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, -t82, 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t47 - t25 * t83, -t15 * t47 + t25 * t84, (t13 * t51 + t15 * t55) * t45, -pkin(2) * t25 * t45 - t13 * t29 + t15 * t30, 0, 0, 0, 0, 0, 0, t13 * t26 + t5 * t83, t13 * t28 + t7 * t83, -t26 * t7 + t28 * t5, -t11 * t5 + t12 * t7 + t13 * t22, 0, 0, 0, 0, 0, 0, t16 * t5 + t26 * t3, t18 * t5 - t26 * t4, -t16 * t4 - t18 * t3, t1 * t3 + t2 * t4 + t5 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t39 * t51 ^ 2, 0.2e1 * t51 * t85, t80 * t99, t36, t47 * t69, t47 ^ 2, 0.2e1 * pkin(2) * t85 + 0.2e1 * t29 * t47, -0.2e1 * t30 * t47 - 0.2e1 * t39 * t96, (-t29 * t51 + t30 * t55) * t99, pkin(2) ^ 2 * t39 + t29 ^ 2 + t30 ^ 2, t28 ^ 2, t26 * t100, t83 * t100, t102, t26 * t69, t36, -0.2e1 * t11 * t83 + 0.2e1 * t22 * t26, 0.2e1 * t12 * t83 + 0.2e1 * t22 * t28, -0.2e1 * t11 * t28 - 0.2e1 * t12 * t26, t11 ^ 2 + t12 ^ 2 + t22 ^ 2, t18 ^ 2, t18 * t101, 0.2e1 * t18 * t26, t16 ^ 2, t26 * t101, t102, 0.2e1 * t1 * t26 + 0.2e1 * t16 * t9, 0.2e1 * t18 * t9 - 0.2e1 * t2 * t26, -0.2e1 * t1 * t18 - 0.2e1 * t16 * t2, t1 ^ 2 + t2 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15, 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t54, t13 * t50, t61, -t13 * pkin(3) + pkin(9) * t61, 0, 0, 0, 0, 0, 0, -t3 * t54 + t5 * t77, t4 * t54 + t5 * t72, (-t3 * t53 - t4 * t49) * t50, t19 * t3 + t20 * t4 + t5 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, t83, t47, t29, -t30, 0, 0, t86, t28 * t54 - t74, -t66, -t87, -t64, 0, -pkin(3) * t26 + pkin(9) * t66 - t22 * t54, -pkin(3) * t28 + pkin(9) * t64 + t22 * t50, (t86 - t87) * pkin(9) + t60, -t22 * pkin(3) + pkin(9) * t60, t18 * t72, (-t88 - t89) * t50, -t18 * t54 + t26 * t72, t16 * t77, t16 * t54 - t49 * t74, -t87, -t1 * t54 + t19 * t26 + (pkin(9) * t16 + t91) * t50, t2 * t54 - t20 * t26 + (pkin(9) * t18 + t90) * t50, -t20 * t16 - t19 * t18 + (-t1 * t53 - t2 * t49) * t50, t1 * t19 + t2 * t20 + t9 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t42, t68, 0, t44, 0, 0, pkin(3) * t97, pkin(3) * t98, 0.2e1 * (t42 + t44) * pkin(9), pkin(3) ^ 2 + t44 * t58 + t38, t43 * t42, -0.2e1 * t42 * t76, t71 * t98, t41 * t42, t49 * t68, t44, -0.2e1 * t19 * t54 + 0.2e1 * t49 * t93, 0.2e1 * t20 * t54 + 0.2e1 * t53 * t93, 0.2e1 * (-t19 * t53 - t20 * t49) * t50, t19 ^ 2 + t20 ^ 2 + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t7, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * t53, t5 * t49, t62, -t5 * pkin(4) + pkin(10) * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t26, -t83, t11, -t12, 0, 0, t88, -t16 * t49 + t18 * t53, t78, -t89, t73, 0, -pkin(4) * t16 - pkin(10) * t78 - t90, -pkin(4) * t18 - pkin(10) * t73 + t91, (t88 - t89) * pkin(10) + t63, -t9 * pkin(4) + pkin(10) * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, t54, 0, -t92, -t54 * pkin(9), 0, 0, t65, (-t41 + t43) * t50, -t75, -t65, -t71, 0, -pkin(9) * t72 + (-pkin(4) * t50 + pkin(10) * t54) * t49, pkin(10) * t71 + (pkin(9) * t49 - t94) * t50, t59, -pkin(4) * t92 + pkin(10) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t41, 0.2e1 * t76, 0, t43, 0, 0, 0.2e1 * t94, -0.2e1 * pkin(4) * t49, 0.2e1 * t70 * pkin(10), pkin(10) ^ 2 * t70 + pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t16, t26, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, -t77, -t54, t19, -t20, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, t53, 0, -t49 * pkin(10), -t53 * pkin(10), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t6;
