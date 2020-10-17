% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x33]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 12:24:43
% EndTime: 2019-05-07 12:24:49
% DurationCPUTime: 1.26s
% Computational Cost: add. (2299->181), mult. (5299->348), div. (0->0), fcn. (6352->12), ass. (0->118)
t84 = sin(pkin(6));
t90 = sin(qJ(2));
t114 = t84 * t90;
t86 = cos(pkin(6));
t89 = sin(qJ(3));
t93 = cos(qJ(3));
t58 = t93 * t114 + t86 * t89;
t83 = sin(pkin(12));
t85 = cos(pkin(12));
t72 = t89 * t114;
t97 = t86 * t93 - t72;
t39 = t58 * t83 - t85 * t97;
t133 = -0.2e1 * t39;
t132 = 0.2e1 * t39;
t65 = t83 * t89 - t85 * t93;
t131 = -0.2e1 * t65;
t130 = 0.2e1 * t65;
t77 = -pkin(3) * t85 - pkin(4);
t92 = cos(qJ(5));
t70 = -pkin(5) * t92 + t77;
t129 = 0.2e1 * t70;
t128 = 0.2e1 * t93;
t127 = pkin(1) * t90;
t94 = cos(qJ(2));
t126 = pkin(1) * t94;
t125 = pkin(5) * t39;
t124 = pkin(5) * t65;
t113 = t84 * t94;
t40 = t85 * t58 + t83 * t97;
t88 = sin(qJ(5));
t32 = t92 * t113 + t40 * t88;
t101 = pkin(8) * t113;
t54 = t101 + (pkin(9) + t127) * t86;
t55 = (-pkin(2) * t94 - pkin(9) * t90 - pkin(1)) * t84;
t34 = -t54 * t89 + t93 * t55;
t23 = -pkin(3) * t113 - qJ(4) * t58 + t34;
t35 = t93 * t54 + t89 * t55;
t28 = t97 * qJ(4) + t35;
t14 = t83 * t23 + t85 * t28;
t12 = -pkin(10) * t113 + t14;
t73 = pkin(8) * t114;
t79 = -pkin(3) * t93 - pkin(2);
t44 = t72 * pkin(3) + t73 + (t79 - t126) * t86;
t16 = t39 * pkin(4) - t40 * pkin(10) + t44;
t7 = t12 * t92 + t16 * t88;
t5 = -pkin(11) * t32 + t7;
t91 = cos(qJ(6));
t123 = t5 * t91;
t87 = sin(qJ(6));
t122 = t87 * pkin(5);
t121 = t91 * pkin(5);
t76 = pkin(3) * t83 + pkin(10);
t120 = pkin(11) + t76;
t33 = -t88 * t113 + t40 * t92;
t119 = t33 * t88;
t104 = -qJ(4) - pkin(9);
t71 = t104 * t93;
t98 = t104 * t89;
t50 = -t85 * t71 + t83 * t98;
t118 = t50 * t92;
t69 = t87 * t92 + t88 * t91;
t117 = t69 * t39;
t116 = t69 * t65;
t80 = t84 ^ 2;
t115 = t80 * t94;
t112 = t88 * t39;
t111 = t88 * t65;
t66 = t83 * t93 + t85 * t89;
t110 = t88 * t66;
t109 = t88 * t76;
t108 = t88 * t92;
t47 = pkin(4) * t65 - pkin(10) * t66 + t79;
t21 = t118 + (-pkin(11) * t66 + t47) * t88;
t107 = t91 * t21;
t106 = t92 * t66;
t105 = t92 * t76;
t103 = 0.2e1 * t84 * t86;
t102 = -0.2e1 * t113;
t100 = t89 * t113;
t99 = t93 * t113;
t6 = -t12 * t88 + t92 * t16;
t4 = -pkin(11) * t33 + t125 + t6;
t1 = t91 * t4 - t5 * t87;
t24 = t92 * t47 - t50 * t88;
t20 = -pkin(11) * t106 + t124 + t24;
t9 = t91 * t20 - t21 * t87;
t13 = t23 * t85 - t83 * t28;
t48 = -t83 * t71 - t85 * t98;
t11 = pkin(4) * t113 - t13;
t96 = -t65 * t76 + t66 * t77;
t68 = t87 * t88 - t91 * t92;
t82 = t92 ^ 2;
t81 = t88 ^ 2;
t64 = t66 ^ 2;
t63 = t65 ^ 2;
t62 = t120 * t92;
t61 = t120 * t88;
t60 = t86 * t127 + t101;
t59 = t86 * t126 - t73;
t57 = t92 * t65;
t53 = t73 + (-pkin(2) - t126) * t86;
t51 = t68 * t65;
t46 = -t61 * t87 + t62 * t91;
t45 = -t61 * t91 - t62 * t87;
t42 = t68 * t66;
t41 = t69 * t66;
t38 = t39 ^ 2;
t37 = t92 * t39;
t36 = pkin(5) * t110 + t48;
t31 = t68 * t39;
t30 = t39 * t65;
t25 = t47 * t88 + t118;
t18 = -t32 * t87 + t33 * t91;
t17 = t91 * t32 + t33 * t87;
t10 = t20 * t87 + t107;
t8 = pkin(5) * t32 + t11;
t2 = t4 * t87 + t123;
t3 = [1, 0, 0, t80 * t90 ^ 2, 0.2e1 * t90 * t115, t90 * t103, t94 * t103, t86 ^ 2, 0.2e1 * pkin(1) * t115 + 0.2e1 * t59 * t86, -0.2e1 * t80 * t127 - 0.2e1 * t60 * t86, t58 ^ 2, 0.2e1 * t58 * t97, t58 * t102, t97 * t102, t80 * t94 ^ 2, -0.2e1 * t113 * t34 - 0.2e1 * t53 * t97, 0.2e1 * t113 * t35 + 0.2e1 * t53 * t58, -0.2e1 * t13 * t40 - 0.2e1 * t14 * t39, t13 ^ 2 + t14 ^ 2 + t44 ^ 2, t33 ^ 2, -0.2e1 * t33 * t32, t33 * t132, t32 * t133, t38, 0.2e1 * t11 * t32 + 0.2e1 * t39 * t6, 0.2e1 * t11 * t33 - 0.2e1 * t39 * t7, t18 ^ 2, -0.2e1 * t18 * t17, t18 * t132, t17 * t133, t38, 0.2e1 * t1 * t39 + 0.2e1 * t17 * t8, 0.2e1 * t18 * t8 - 0.2e1 * t2 * t39; 0, 0, 0, 0, 0, t114, t113, t86, t59, -t60, t58 * t89, t58 * t93 + t89 * t97, -t100, -t99, 0, pkin(2) * t97 + pkin(9) * t100 - t53 * t93, -pkin(2) * t58 + pkin(9) * t99 + t53 * t89, -t13 * t66 - t14 * t65 - t39 * t50 + t40 * t48, -t13 * t48 + t14 * t50 + t44 * t79, t33 * t106 (-t32 * t92 - t119) * t66, t106 * t39 + t33 * t65, -t110 * t39 - t32 * t65, t30, t11 * t110 + t24 * t39 + t32 * t48 + t6 * t65, t106 * t11 - t25 * t39 + t33 * t48 - t65 * t7, -t18 * t42, t17 * t42 - t18 * t41, t18 * t65 - t39 * t42, -t17 * t65 - t39 * t41, t30, t1 * t65 + t17 * t36 + t39 * t9 + t41 * t8, -t10 * t39 + t18 * t36 - t2 * t65 - t42 * t8; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t89 ^ 2, t89 * t128, 0, 0, 0, pkin(2) * t128, -0.2e1 * pkin(2) * t89, 0.2e1 * t48 * t66 - 0.2e1 * t50 * t65, t48 ^ 2 + t50 ^ 2 + t79 ^ 2, t82 * t64, -0.2e1 * t64 * t108, t106 * t130, t110 * t131, t63, 0.2e1 * t110 * t48 + 0.2e1 * t24 * t65, 0.2e1 * t106 * t48 - 0.2e1 * t25 * t65, t42 ^ 2, 0.2e1 * t42 * t41, -t42 * t130, t41 * t131, t63, 0.2e1 * t36 * t41 + 0.2e1 * t65 * t9, -0.2e1 * t10 * t65 - 0.2e1 * t36 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, t97, -t113, t34, -t35 (-t39 * t83 - t40 * t85) * pkin(3) (t13 * t85 + t14 * t83) * pkin(3), t119, -t32 * t88 + t33 * t92, t112, t37, 0, -t109 * t39 - t11 * t92 + t32 * t77, -t105 * t39 + t11 * t88 + t33 * t77, t18 * t69, -t17 * t69 - t18 * t68, t117, -t31, 0, t17 * t70 + t39 * t45 + t68 * t8, t18 * t70 - t39 * t46 + t69 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t93, 0, -t89 * pkin(9), -t93 * pkin(9) (-t65 * t83 - t66 * t85) * pkin(3) (-t48 * t85 + t50 * t83) * pkin(3), t88 * t106 (-t81 + t82) * t66, t111, t57, 0, -t48 * t92 + t88 * t96, t48 * t88 + t92 * t96, -t42 * t69, -t41 * t69 + t42 * t68, t116, -t51, 0, t36 * t68 + t41 * t70 + t45 * t65, t36 * t69 - t42 * t70 - t46 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t83 ^ 2 + t85 ^ 2) * pkin(3) ^ 2, t81, 0.2e1 * t108, 0, 0, 0, -0.2e1 * t77 * t92, 0.2e1 * t77 * t88, t69 ^ 2, -0.2e1 * t69 * t68, 0, 0, 0, t68 * t129, t69 * t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, 0, 0, 0, 0, t37, -t112, 0, 0, 0, 0, 0, -t31, -t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, 0, 0, 0, 0, t57, -t111, 0, 0, 0, 0, 0, -t51, -t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32, t39, t6, -t7, 0, 0, t18, -t17, t39, t39 * t121 + t1, -t123 + (-t4 - t125) * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, -t110, t65, t24, -t25, 0, 0, -t42, -t41, t65, t65 * t121 + t9, -t107 + (-t20 - t124) * t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, t92, 0, -t109, -t105, 0, 0, t69, -t68, 0, t45, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, -t88, 0, 0, 0, 0, 0, -t68, -t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t121, -0.2e1 * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t17, t39, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42, -t41, t65, t9, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t68, 0, t45, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t121, -t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t3;
