% Calculate inertial parameters regressor of joint inertia matrix for
% S6PPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PPRRRP1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_inertiaJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t60 = sin(pkin(12));
t62 = sin(pkin(6));
t65 = cos(pkin(6));
t68 = sin(qJ(3));
t71 = cos(qJ(3));
t63 = cos(pkin(12));
t64 = cos(pkin(7));
t87 = t63 * t64;
t61 = sin(pkin(7));
t89 = t61 * t68;
t21 = t65 * t89 + (t60 * t71 + t68 * t87) * t62;
t29 = -t62 * t63 * t61 + t65 * t64;
t67 = sin(qJ(4));
t70 = cos(qJ(4));
t13 = t21 * t67 - t29 * t70;
t105 = t13 ^ 2;
t88 = t61 * t71;
t19 = -t65 * t88 + (t60 * t68 - t71 * t87) * t62;
t104 = t19 ^ 2;
t30 = -t70 * t64 + t67 * t89;
t103 = t30 ^ 2;
t102 = -0.2e1 * t67;
t101 = 0.2e1 * t67;
t100 = 0.2e1 * t70;
t69 = cos(qJ(5));
t99 = pkin(4) * t69;
t66 = sin(qJ(5));
t98 = pkin(9) * t66;
t57 = t67 ^ 2;
t97 = t57 * pkin(9);
t96 = t66 * pkin(5);
t95 = t67 * pkin(9);
t94 = t13 * t30;
t93 = t13 * t67;
t92 = t13 * t69;
t91 = t30 * t67;
t90 = t30 * t69;
t86 = t66 * t67;
t85 = t66 * t69;
t84 = t66 * t70;
t51 = t69 * t67;
t83 = t69 * t70;
t82 = -qJ(6) - pkin(10);
t56 = t66 ^ 2;
t58 = t69 ^ 2;
t81 = t56 + t58;
t80 = qJ(6) * t67;
t79 = t67 * t100;
t78 = pkin(9) * t83;
t37 = -t70 * pkin(4) - t67 * pkin(10) - pkin(3);
t33 = t69 * t37;
t77 = -t69 * t80 + t33;
t15 = t21 * t70 + t29 * t67;
t7 = -t15 * t66 + t19 * t69;
t8 = t15 * t69 + t19 * t66;
t4 = -t7 * t66 + t8 * t69;
t76 = t15 * t70 + t93;
t32 = t67 * t64 + t70 * t89;
t22 = -t66 * t32 - t69 * t88;
t23 = t69 * t32 - t66 * t88;
t11 = -t22 * t66 + t23 * t69;
t25 = -pkin(9) * t84 + t33;
t26 = t66 * t37 + t78;
t75 = -t25 * t66 + t26 * t69;
t74 = t32 * t70 + t91;
t73 = pkin(9) ^ 2;
t59 = t70 ^ 2;
t54 = t61 ^ 2;
t53 = t57 * t73;
t52 = -t69 * pkin(5) - pkin(4);
t50 = t58 * t57;
t49 = t56 * t57;
t47 = t54 * t71 ^ 2;
t46 = 0.2e1 * t85;
t45 = t66 * t51;
t42 = t83 * t102;
t41 = -0.2e1 * t57 * t85;
t40 = t66 * t79;
t39 = t82 * t69;
t38 = t82 * t66;
t36 = (pkin(9) + t96) * t67;
t34 = (-t56 + t58) * t67;
t28 = t30 * t66;
t24 = t78 + (t37 - t80) * t66;
t18 = (-pkin(5) - t98) * t70 + t77;
t17 = t23 * t70 + t30 * t51;
t16 = -t22 * t70 + t30 * t86;
t12 = t13 * t66;
t10 = (-t22 * t69 - t23 * t66) * t67;
t9 = t22 ^ 2 + t23 ^ 2 + t103;
t6 = t13 * t51 + t8 * t70;
t5 = t13 * t86 - t7 * t70;
t3 = (-t66 * t8 - t69 * t7) * t67;
t2 = t7 ^ 2 + t8 ^ 2 + t105;
t1 = t7 * t22 + t8 * t23 + t94;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65 ^ 2 + (t60 ^ 2 + t63 ^ 2) * t62 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 ^ 2 + t29 ^ 2 + t104, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 ^ 2 + t104 + t105, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29 * t64 + (-t19 * t71 + t21 * t68) * t61, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 * t32 - t19 * t88 + t94, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t68 ^ 2 + t64 ^ 2 + t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 ^ 2 + t103 + t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t21, 0, 0, 0, 0, 0, 0, 0, 0, -t19 * t70, t19 * t67, t76, -t19 * pkin(3) + t76 * pkin(9), 0, 0, 0, 0, 0, 0, t5, t6, t3, pkin(9) * t93 + t7 * t25 + t8 * t26, 0, 0, 0, 0, 0, 0, t5, t6, t3, t13 * t36 + t7 * t18 + t8 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, -t89, 0, 0, 0, 0, 0, 0, 0, 0, t70 * t88, -t67 * t88, t74, pkin(3) * t88 + t74 * pkin(9), 0, 0, 0, 0, 0, 0, t16, t17, t10, pkin(9) * t91 + t22 * t25 + t23 * t26, 0, 0, 0, 0, 0, 0, t16, t17, t10, t22 * t18 + t23 * t24 + t30 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t57, t79, 0, t59, 0, 0, pkin(3) * t100, pkin(3) * t102, 0.2e1 * (t57 + t59) * pkin(9), pkin(3) ^ 2 + t59 * t73 + t53, t50, t41, t42, t49, t40, t59, -0.2e1 * t25 * t70 + 0.2e1 * t66 * t97, 0.2e1 * t26 * t70 + 0.2e1 * t69 * t97 (-t25 * t69 - t26 * t66) * t101, t25 ^ 2 + t26 ^ 2 + t53, t50, t41, t42, t49, t40, t59, -0.2e1 * t18 * t70 + 0.2e1 * t36 * t86, 0.2e1 * t24 * t70 + 0.2e1 * t36 * t51 (-t18 * t69 - t24 * t66) * t101, t18 ^ 2 + t24 ^ 2 + t36 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15, 0, 0, 0, 0, 0, 0, 0, 0, -t92, t12, t4, -t13 * pkin(4) + t4 * pkin(10), 0, 0, 0, 0, 0, 0, -t92, t12, t4, t13 * t52 + t7 * t38 - t8 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t32, 0, 0, 0, 0, 0, 0, 0, 0, -t90, t28, t11, -t30 * pkin(4) + t11 * pkin(10), 0, 0, 0, 0, 0, 0, -t90, t28, t11, t22 * t38 - t23 * t39 + t30 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, t70, 0, -t95, -t70 * pkin(9), 0, 0, t45, t34, -t84, -t45, -t83, 0, -pkin(9) * t51 + (-pkin(4) * t67 + pkin(10) * t70) * t66, pkin(10) * t83 + (t98 - t99) * t67, t75, -pkin(4) * t95 + t75 * pkin(10), t45, t34, -t84, -t45, -t83, 0, -t36 * t69 - t38 * t70 + t52 * t86, t36 * t66 - t39 * t70 + t52 * t51 (-t38 * t67 + t24) * t69 + (t39 * t67 - t18) * t66, t18 * t38 - t24 * t39 + t36 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t56, t46, 0, t58, 0, 0, 0.2e1 * t99, -0.2e1 * pkin(4) * t66, 0.2e1 * t81 * pkin(10), t81 * pkin(10) ^ 2 + pkin(4) ^ 2, t56, t46, 0, t58, 0, 0, -0.2e1 * t52 * t69, 0.2e1 * t52 * t66, -0.2e1 * t38 * t66 - 0.2e1 * t39 * t69, t38 ^ 2 + t39 ^ 2 + t52 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t8, 0, t7 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t23, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t23, 0, t22 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, -t86, -t70, t25, -t26, 0, 0, 0, 0, t51, 0, -t86, -t70 (-0.2e1 * pkin(5) - t98) * t70 + t77, -t24, -pkin(5) * t51, t18 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, t69, 0, -t66 * pkin(10), -t69 * pkin(10), 0, 0, 0, 0, t66, 0, t69, 0, t38, t39, -t96, t38 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * pkin(5), 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t51, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, t66, 0, t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t14;
