% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPR1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t91 = sin(pkin(11));
t87 = t91 ^ 2;
t93 = cos(pkin(11));
t88 = t93 ^ 2;
t109 = t87 + t88;
t110 = t109 * qJ(5);
t123 = cos(qJ(4));
t92 = sin(pkin(10));
t94 = cos(pkin(10));
t97 = sin(qJ(2));
t98 = cos(qJ(2));
t68 = -t92 * t97 + t94 * t98;
t71 = t92 * t98 + t94 * t97;
t96 = sin(qJ(4));
t48 = t123 * t71 + t96 * t68;
t137 = -0.2e1 * t48;
t122 = cos(qJ(6));
t95 = sin(qJ(6));
t136 = t122 * t93 - t95 * t91;
t115 = -qJ(3) - pkin(7);
t75 = t115 * t97;
t76 = t115 * t98;
t54 = t94 * t75 + t92 * t76;
t102 = -t71 * pkin(8) + t54;
t55 = t92 * t75 - t94 * t76;
t37 = t68 * pkin(8) + t55;
t22 = -t123 * t102 + t96 * t37;
t135 = t22 ^ 2;
t46 = -t123 * t68 + t96 * t71;
t44 = t46 ^ 2;
t134 = 0.2e1 * t46;
t126 = t93 * pkin(5);
t127 = t92 * pkin(2);
t125 = t94 * pkin(2);
t81 = pkin(3) + t125;
t63 = t123 * t81 - t96 * t127;
t62 = -pkin(4) - t63;
t56 = t62 - t126;
t133 = 0.2e1 * t56;
t83 = -t98 * pkin(2) - pkin(1);
t59 = -t68 * pkin(3) + t83;
t132 = 0.2e1 * t59;
t131 = 0.2e1 * t71;
t82 = -pkin(4) - t126;
t130 = 0.2e1 * t82;
t129 = 0.2e1 * t98;
t117 = t93 * t48;
t21 = t46 * pkin(4) - t48 * qJ(5) + t59;
t24 = t102 * t96 + t123 * t37;
t8 = t93 * t21 - t91 * t24;
t6 = t46 * pkin(5) - pkin(9) * t117 + t8;
t119 = t91 * t48;
t9 = t91 * t21 + t93 * t24;
t7 = -pkin(9) * t119 + t9;
t3 = t122 * t6 - t95 * t7;
t4 = t122 * t7 + t95 * t6;
t72 = t122 * t91 + t95 * t93;
t128 = t136 * t4 - t3 * t72;
t124 = pkin(4) - t62;
t121 = t22 * t93;
t27 = t136 * t48;
t120 = t27 * t136;
t35 = t72 * t46;
t42 = t91 * t46;
t118 = t91 * t93;
t64 = t123 * t127 + t96 * t81;
t61 = qJ(5) + t64;
t52 = (-pkin(9) - t61) * t91;
t86 = t93 * pkin(9);
t53 = t93 * t61 + t86;
t31 = t122 * t52 - t95 * t53;
t32 = t122 * t53 + t95 * t52;
t114 = t136 * t32 - t31 * t72;
t73 = (-pkin(9) - qJ(5)) * t91;
t74 = t93 * qJ(5) + t86;
t50 = t122 * t73 - t95 * t74;
t51 = t122 * t74 + t95 * t73;
t113 = t136 * t51 - t50 * t72;
t112 = t56 + t82;
t111 = t109 * t61;
t89 = t97 ^ 2;
t90 = t98 ^ 2;
t108 = t89 + t90;
t107 = t46 * t137;
t105 = t8 * t93 + t9 * t91;
t5 = -t8 * t91 + t9 * t93;
t104 = -pkin(4) * t48 - qJ(5) * t46;
t103 = -t46 * t61 + t48 * t62;
t77 = 0.2e1 * t118;
t67 = t72 ^ 2;
t66 = t136 ^ 2;
t49 = 0.2e1 * t72 * t136;
t45 = t48 ^ 2;
t43 = t93 * t46;
t40 = t91 * t117;
t34 = t136 * t46;
t30 = (-t87 + t88) * t48;
t25 = t72 * t48;
t20 = t22 * t91;
t18 = t72 * t25;
t17 = t27 * t72;
t16 = t25 * t136;
t13 = pkin(5) * t119 + t22;
t12 = t13 * t72;
t11 = t13 * t136;
t10 = -t18 + t120;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t89, t97 * t129, 0, t90, 0, 0, pkin(1) * t129, -0.2e1 * pkin(1) * t97, 0.2e1 * t108 * pkin(7), pkin(7) ^ 2 * t108 + pkin(1) ^ 2, t71 ^ 2, t68 * t131, 0, t68 ^ 2, 0, 0, -0.2e1 * t83 * t68, t83 * t131, -0.2e1 * t54 * t71 + 0.2e1 * t55 * t68, t54 ^ 2 + t55 ^ 2 + t83 ^ 2, t45, t107, 0, t44, 0, 0, t46 * t132, t48 * t132, 0.2e1 * t22 * t48 - 0.2e1 * t24 * t46, t24 ^ 2 + t59 ^ 2 + t135, t88 * t45, -0.2e1 * t45 * t118, t117 * t134, t87 * t45, t91 * t107, t44, 0.2e1 * t119 * t22 + 0.2e1 * t8 * t46, 0.2e1 * t117 * t22 - 0.2e1 * t9 * t46, t105 * t137, t8 ^ 2 + t9 ^ 2 + t135, t27 ^ 2, -0.2e1 * t27 * t25, t27 * t134, t25 ^ 2, -t25 * t134, t44, 0.2e1 * t13 * t25 + 0.2e1 * t3 * t46, 0.2e1 * t13 * t27 - 0.2e1 * t4 * t46, -0.2e1 * t4 * t25 - 0.2e1 * t3 * t27, t13 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, t98, 0, -t97 * pkin(7), -t98 * pkin(7), 0, 0, 0, 0, t71, 0, t68, 0, t54, -t55 (t68 * t92 - t71 * t94) * pkin(2) (t54 * t94 + t55 * t92) * pkin(2), 0, 0, t48, 0, -t46, 0, -t22, -t24, -t64 * t46 - t63 * t48, -t22 * t63 + t24 * t64, t40, t30, t42, -t40, t43, 0, t103 * t91 - t121, t103 * t93 + t20, t5, t22 * t62 + t5 * t61, t17, t10, t35, -t16, t34, 0, t56 * t25 + t31 * t46 - t11, t56 * t27 - t32 * t46 + t12, -t32 * t25 - t31 * t27 + t128, t13 * t56 + t3 * t31 + t4 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t125, -0.2e1 * t127, 0 (t92 ^ 2 + t94 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t63, -0.2e1 * t64, 0, t63 ^ 2 + t64 ^ 2, t87, t77, 0, t88, 0, 0, -0.2e1 * t62 * t93, 0.2e1 * t62 * t91, 0.2e1 * t111, t109 * t61 ^ 2 + t62 ^ 2, t67, t49, 0, t66, 0, 0, -t136 * t133, t72 * t133, 0.2e1 * t114, t31 ^ 2 + t32 ^ 2 + t56 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, t71, 0, t83, 0, 0, 0, 0, 0, 0, t46, t48, 0, t59, 0, 0, 0, 0, 0, 0, t43, -t42, -t109 * t48, t105, 0, 0, 0, 0, 0, 0, t34, -t35, -t18 - t120, t136 * t3 + t4 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136 * t31 + t32 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67 + t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, -t46, 0, -t22, -t24, 0, 0, t40, t30, t42, -t40, t43, 0, t104 * t91 - t121, t104 * t93 + t20, t5, -t22 * pkin(4) + qJ(5) * t5, t17, t10, t35, -t16, t34, 0, t82 * t25 + t50 * t46 - t11, t82 * t27 - t51 * t46 + t12, -t51 * t25 - t50 * t27 + t128, t13 * t82 + t3 * t50 + t4 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t63, -t64, 0, 0, t87, t77, 0, t88, 0, 0, t124 * t93, -t124 * t91, t110 + t111, -t62 * pkin(4) + t61 * t110, t67, t49, 0, t66, 0, 0, -t112 * t136, t112 * t72, t113 + t114, t31 * t50 + t32 * t51 + t56 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136 * t50 + t72 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t87, t77, 0, t88, 0, 0, 0.2e1 * pkin(4) * t93, -0.2e1 * pkin(4) * t91, 0.2e1 * t110, qJ(5) ^ 2 * t109 + pkin(4) ^ 2, t67, t49, 0, t66, 0, 0, -t136 * t130, t72 * t130, 0.2e1 * t113, t50 ^ 2 + t51 ^ 2 + t82 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, t117, 0, t22, 0, 0, 0, 0, 0, 0, t25, t27, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, t91, 0, t62, 0, 0, 0, 0, 0, 0, -t136, t72, 0, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, t91, 0, -pkin(4), 0, 0, 0, 0, 0, 0, -t136, t72, 0, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, -t25, t46, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, t136, 0, t31, -t32, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, -t72, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, t136, 0, t50, -t51, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t1;