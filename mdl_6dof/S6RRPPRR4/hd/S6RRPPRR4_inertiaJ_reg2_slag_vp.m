% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPPRR4_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t68 = sin(pkin(6));
t76 = cos(qJ(2));
t112 = t68 * t76;
t73 = sin(qJ(2));
t113 = t68 * t73;
t67 = sin(pkin(11));
t69 = cos(pkin(11));
t38 = -t69 * t112 + t67 * t113;
t70 = cos(pkin(6));
t72 = sin(qJ(5));
t75 = cos(qJ(5));
t26 = -t38 * t75 + t70 * t72;
t143 = t26 ^ 2;
t142 = t38 ^ 2;
t40 = (-t67 * t76 - t69 * t73) * t68;
t37 = t40 ^ 2;
t129 = t67 * pkin(2);
t52 = qJ(4) + t129;
t141 = t52 ^ 2;
t28 = t38 * t72 + t70 * t75;
t71 = sin(qJ(6));
t74 = cos(qJ(6));
t17 = t28 * t71 + t40 * t74;
t140 = -0.2e1 * t17;
t139 = -0.2e1 * t28;
t138 = 0.2e1 * t40;
t137 = 0.2e1 * t52;
t136 = 0.2e1 * t68;
t135 = 0.2e1 * t74;
t134 = 0.2e1 * t75;
t133 = pkin(3) + pkin(9);
t132 = pkin(1) * t73;
t45 = (-pkin(2) * t76 - pkin(1)) * t68;
t79 = t40 * qJ(4) + t45;
t10 = t133 * t38 + t79;
t127 = t70 * pkin(2);
t50 = t70 * t76 * pkin(1);
t99 = pkin(8) + qJ(3);
t29 = -t99 * t113 + t127 + t50;
t93 = t70 * t132;
t32 = t99 * t112 + t93;
t15 = t69 * t29 - t67 * t32;
t8 = -t40 * pkin(4) - t133 * t70 - t15;
t5 = -t72 * t10 + t75 * t8;
t3 = t40 * pkin(5) - t5;
t131 = t3 * t71;
t130 = t3 * t74;
t128 = t69 * pkin(2);
t126 = t72 * pkin(5);
t125 = t75 * pkin(5);
t124 = t17 * t74;
t19 = t28 * t74 - t40 * t71;
t123 = t19 * t71;
t122 = t26 * t72;
t121 = t28 * t72;
t120 = t28 * t75;
t119 = t40 * t70;
t118 = t40 * t72;
t117 = t40 * t75;
t61 = t68 ^ 2;
t116 = t61 * t76;
t63 = t71 ^ 2;
t115 = t63 * t75;
t54 = -pkin(3) - t128;
t51 = -pkin(9) + t54;
t66 = t75 ^ 2;
t114 = t66 * t51;
t111 = t70 * t38;
t110 = t71 * t26;
t56 = t71 * t72;
t109 = t71 * t74;
t108 = t71 * t75;
t107 = t72 * t17;
t106 = t72 * t51;
t105 = t74 * t26;
t104 = t74 * t72;
t59 = t74 * t75;
t103 = t75 * t17;
t102 = t75 * t19;
t21 = t75 * t26;
t101 = t75 * t51;
t100 = t75 * t72;
t16 = t67 * t29 + t69 * t32;
t65 = t74 ^ 2;
t98 = t63 + t65;
t64 = t72 ^ 2;
t97 = t64 + t66;
t96 = t38 * t138;
t95 = t70 * t136;
t94 = -0.2e1 * t100;
t92 = t71 * t102;
t91 = t71 * t21;
t90 = t74 * t21;
t89 = t71 * t59;
t88 = t70 * qJ(4) + t16;
t87 = t98 * pkin(10);
t86 = t98 * t72;
t85 = -pkin(10) * t72 - t125;
t6 = t75 * t10 + t72 * t8;
t4 = -t40 * pkin(10) + t6;
t9 = -t38 * pkin(4) + t88;
t7 = t26 * pkin(5) - t28 * pkin(10) + t9;
t1 = -t71 * t4 + t74 * t7;
t2 = t74 * t4 + t71 * t7;
t84 = -t1 * t71 + t2 * t74;
t83 = t5 * t75 + t6 * t72;
t82 = t123 - t124;
t44 = -t75 * pkin(10) + t126 + t52;
t22 = -t71 * t106 + t74 * t44;
t23 = t51 * t104 + t71 * t44;
t81 = -t22 * t71 + t23 * t74;
t80 = -t120 - t122;
t62 = t70 ^ 2;
t58 = t65 * t75;
t57 = t65 * t66;
t55 = t63 * t66;
t48 = t51 ^ 2;
t46 = t66 * t48;
t43 = pkin(8) * t112 + t93;
t42 = -pkin(8) * t113 + t50;
t36 = t97 * t51;
t20 = t38 * pkin(3) + t79;
t14 = t19 * t72;
t13 = t74 * t103;
t12 = -t70 * pkin(3) - t15;
t11 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t61 * t73 ^ 2, 0.2e1 * t73 * t116, t73 * t95, t61 * t76 ^ 2, t76 * t95, t62, 0.2e1 * pkin(1) * t116 + 0.2e1 * t42 * t70, -0.2e1 * t61 * t132 - 0.2e1 * t43 * t70 (-t42 * t73 + t43 * t76) * t136, t61 * pkin(1) ^ 2 + t42 ^ 2 + t43 ^ 2, t37, t96, -0.2e1 * t119, t142, -0.2e1 * t111, t62, 0.2e1 * t15 * t70 + 0.2e1 * t45 * t38, -0.2e1 * t16 * t70 - 0.2e1 * t45 * t40, 0.2e1 * t15 * t40 - 0.2e1 * t16 * t38, t15 ^ 2 + t16 ^ 2 + t45 ^ 2, t62, 0.2e1 * t119, 0.2e1 * t111, t37, t96, t142, -0.2e1 * t12 * t40 - 0.2e1 * t38 * t88, 0.2e1 * t12 * t70 - 0.2e1 * t20 * t38, 0.2e1 * t20 * t40 + 0.2e1 * t70 * t88, t12 ^ 2 + t20 ^ 2 + t88 ^ 2, t28 ^ 2, t26 * t139, t40 * t139, t143, t26 * t138, t37, 0.2e1 * t9 * t26 - 0.2e1 * t5 * t40, 0.2e1 * t9 * t28 + 0.2e1 * t6 * t40, -0.2e1 * t6 * t26 - 0.2e1 * t5 * t28, t5 ^ 2 + t6 ^ 2 + t9 ^ 2, t19 ^ 2, t19 * t140, 0.2e1 * t19 * t26, t17 ^ 2, t26 * t140, t143, 0.2e1 * t1 * t26 + 0.2e1 * t3 * t17, 0.2e1 * t3 * t19 - 0.2e1 * t2 * t26, -0.2e1 * t1 * t19 - 0.2e1 * t2 * t17, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113, 0, t112, t70, t42, -t43, 0, 0, 0, 0, -t40, 0, -t38, t70, t69 * t127 + t15, -t67 * t127 - t16 (-t38 * t67 + t40 * t69) * pkin(2) (t15 * t69 + t16 * t67) * pkin(2), t70, t40, t38, 0, 0, 0, -t52 * t38 - t54 * t40 (-pkin(3) + t54) * t70 - t15, t52 * t70 + t88, t12 * t54 + t52 * t88, t120, -t21 - t121, -t117, t122, t118, 0, -t101 * t40 + t52 * t26 + t9 * t72, t40 * t106 + t52 * t28 + t9 * t75, t80 * t51 - t83, t83 * t51 + t9 * t52, t74 * t102, -t13 - t92, t14 + t90, t71 * t103, -t91 - t107, t122, t1 * t72 + t22 * t26 + (-t17 * t51 + t131) * t75, -t2 * t72 - t23 * t26 + (-t19 * t51 + t130) * t75, -t23 * t17 - t22 * t19 + (-t1 * t74 - t2 * t71) * t75, t1 * t22 - t101 * t3 + t2 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t128, -0.2e1 * t129, 0 (t67 ^ 2 + t69 ^ 2) * pkin(2) ^ 2, 1, 0, 0, 0, 0, 0, 0, 0.2e1 * t54, t137, t54 ^ 2 + t141, t66, t94, 0, t64, 0, 0, t72 * t137, t52 * t134, -0.2e1 * t36, t64 * t48 + t141 + t46, t57, -0.2e1 * t66 * t109, t100 * t135, t55, t71 * t94, t64, -0.2e1 * t114 * t71 + 0.2e1 * t22 * t72, -0.2e1 * t114 * t74 - 0.2e1 * t23 * t72 (-t22 * t74 - t23 * t71) * t134, t22 ^ 2 + t23 ^ 2 + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t40, 0, t45, 0, 0, 0, 0, 0, 0, 0, -t38, t40, t20, 0, 0, 0, 0, 0, 0, t118, t117, -t21 + t121, -t5 * t72 + t6 * t75, 0, 0, 0, 0, 0, 0, -t91 + t107, t14 - t90, -t13 + t92, t3 * t72 + t75 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t81 - t106) * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57 + t55 + t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t70, 0, t12, 0, 0, 0, 0, 0, 0, -t117, t118, t80, t83, 0, 0, 0, 0, 0, 0, -t26 * t56 - t103, -t104 * t26 - t102, t82 * t72, -t3 * t75 + t72 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t54, 0, 0, 0, 0, 0, 0, 0, 0, -t97, t36, 0, 0, 0, 0, 0, 0, -t97 * t71, -t97 * t74, 0, t72 * t81 + t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-0.1e1 + t98) * t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64 * t98 + t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t26, -t40, t5, -t6, 0, 0, t123, -t71 * t17 + t19 * t74, t110, -t124, t105, 0, -pkin(5) * t17 - pkin(10) * t110 - t130, -pkin(5) * t19 - pkin(10) * t105 + t131, pkin(10) * t82 + t84, -t3 * pkin(5) + pkin(10) * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, -t72, 0, t101, -t106, 0, 0, t89, t58 - t115, t56, -t89, t104, 0, t101 * t74 + t71 * t85, -t101 * t71 + t74 * t85, t81, pkin(5) * t101 + pkin(10) * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t75, 0, 0, 0, 0, 0, 0, 0, 0, -t104, t56, t58 + t115, t75 * t87 - t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t72, 0, 0, 0, 0, 0, 0, 0, 0, t59, -t108, t86, pkin(10) * t86 + t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t63, 0.2e1 * t109, 0, t65, 0, 0, pkin(5) * t135, -0.2e1 * pkin(5) * t71, 0.2e1 * t87, pkin(10) ^ 2 * t98 + pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, -t17, t26, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, -t108, t72, t22, -t23, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, -t59, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t104, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, t74, 0, -t71 * pkin(10), -t74 * pkin(10), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t11;