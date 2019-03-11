% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRPP7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t100 = cos(pkin(11));
t102 = sin(qJ(4));
t105 = cos(qJ(4));
t107 = cos(qJ(2));
t99 = sin(pkin(6));
t138 = t99 * t107;
t101 = cos(pkin(6));
t103 = sin(qJ(3));
t106 = cos(qJ(3));
t104 = sin(qJ(2));
t139 = t99 * t104;
t67 = t101 * t103 + t106 * t139;
t44 = t67 * t102 + t105 * t138;
t46 = -t102 * t138 + t67 * t105;
t98 = sin(pkin(11));
t24 = t100 * t44 + t98 * t46;
t183 = t24 ^ 2;
t73 = t100 * t102 + t98 * t105;
t60 = t73 * t103;
t182 = t60 ^ 2;
t65 = -t101 * t106 + t103 * t139;
t64 = t65 ^ 2;
t71 = -t100 * t105 + t98 * t102;
t181 = t71 ^ 2;
t180 = -0.2e1 * t44;
t179 = 0.2e1 * t65;
t178 = -0.2e1 * t67;
t90 = -t105 * pkin(4) - pkin(3);
t177 = 0.2e1 * t90;
t176 = 0.2e1 * t99;
t175 = -0.2e1 * t103;
t174 = 0.2e1 * t106;
t170 = pkin(1) * t107;
t81 = pkin(8) * t139;
t56 = t81 + (-pkin(2) - t170) * t101;
t28 = t65 * pkin(3) - t67 * pkin(10) + t56;
t129 = pkin(8) * t138;
t171 = pkin(1) * t104;
t57 = t129 + (pkin(9) + t171) * t101;
t58 = (-pkin(2) * t107 - pkin(9) * t104 - pkin(1)) * t99;
t33 = t103 * t58 + t106 * t57;
t30 = -pkin(10) * t138 + t33;
t12 = t102 * t28 + t105 * t30;
t10 = -t44 * qJ(5) + t12;
t11 = -t102 * t30 + t105 * t28;
t8 = t65 * pkin(4) - t46 * qJ(5) + t11;
t4 = t100 * t10 + t98 * t8;
t173 = t98 * pkin(4);
t167 = t100 * pkin(4);
t88 = pkin(5) + t167;
t172 = pkin(5) + t88;
t169 = pkin(3) * t105;
t168 = pkin(9) * t102;
t92 = t103 * pkin(9);
t166 = t106 * pkin(4);
t165 = t24 * t60;
t164 = t24 * t71;
t26 = t100 * t46 - t98 * t44;
t163 = t26 * t24;
t155 = -qJ(5) - pkin(10);
t119 = t155 * t102;
t78 = t155 * t105;
t49 = -t100 * t119 - t98 * t78;
t162 = t49 * t65;
t51 = -t100 * t78 + t98 * t119;
t161 = t51 * t65;
t160 = t60 * t71;
t133 = t105 * t103;
t136 = t102 * t103;
t62 = t100 * t133 - t98 * t136;
t159 = t62 * t60;
t158 = t65 * t24;
t157 = t65 * t71;
t156 = t73 * t71;
t77 = -t106 * pkin(3) - t103 * pkin(10) - pkin(2);
t74 = t105 * t77;
t41 = -qJ(5) * t133 + t74 + (-pkin(4) - t168) * t106;
t132 = t105 * t106;
t127 = pkin(9) * t132;
t47 = t127 + (-qJ(5) * t103 + t77) * t102;
t22 = t100 * t47 + t98 * t41;
t76 = pkin(4) * t136 + t92;
t94 = t102 ^ 2;
t96 = t105 ^ 2;
t154 = t94 + t96;
t153 = t102 * t65;
t152 = t105 * t65;
t151 = t106 * t60;
t150 = t106 * t71;
t93 = t99 ^ 2;
t149 = t107 * t93;
t32 = -t103 * t57 + t106 * t58;
t29 = pkin(3) * t138 - t32;
t148 = t29 * t102;
t147 = t29 * t105;
t146 = t44 * t105;
t145 = t46 * t102;
t144 = t49 * t106;
t143 = t51 * t106;
t142 = t65 * t106;
t141 = t67 * t103;
t140 = t73 * t106;
t137 = t101 * t104;
t135 = t102 * t105;
t134 = t102 * t106;
t1 = t65 * qJ(6) + t4;
t131 = t49 ^ 2 + t51 ^ 2;
t130 = 0.2e1 * t138;
t128 = t103 * t174;
t126 = t103 * t138;
t125 = t106 * t138;
t124 = t102 * t133;
t123 = t98 * t10 - t100 * t8;
t122 = -t51 * t24 + t49 * t26;
t121 = t49 * t62 - t51 * t60;
t120 = -t100 * t41 + t98 * t47;
t118 = t62 * t24 + t26 * t60;
t117 = t73 * t24 + t26 * t71;
t116 = -t73 * t60 - t62 * t71;
t115 = t106 * t24 - t65 * t60;
t114 = -t11 * t102 + t12 * t105;
t54 = -pkin(9) * t134 + t74;
t55 = t102 * t77 + t127;
t113 = -t54 * t102 + t55 * t105;
t112 = -t32 * t103 + t33 * t106;
t111 = 0.2e1 * t49 * t73 - 0.2e1 * t51 * t71;
t16 = t44 * pkin(4) + t29;
t109 = pkin(9) ^ 2;
t97 = t106 ^ 2;
t95 = t103 ^ 2;
t91 = t95 * t109;
t85 = t93 * t107 ^ 2;
t84 = qJ(6) + t173;
t70 = t73 ^ 2;
t69 = pkin(1) * t137 + t129;
t68 = t101 * t170 - t81;
t59 = t62 ^ 2;
t53 = -0.2e1 * t62 * t106;
t43 = t73 * t65;
t40 = t62 * t73;
t36 = t71 * pkin(5) - t73 * qJ(6) + t90;
t31 = t60 * pkin(5) - t62 * qJ(6) + t76;
t23 = t26 ^ 2;
t20 = t106 * pkin(5) + t120;
t19 = -t106 * qJ(6) + t22;
t18 = t26 * t73;
t17 = t26 * t179;
t15 = t26 * t62;
t13 = -t26 * t106 + t62 * t65;
t5 = t24 * pkin(5) - t26 * qJ(6) + t16;
t2 = -t65 * pkin(5) + t123;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t93 * t104 ^ 2, 0.2e1 * t104 * t149, t137 * t176, t85, t101 * t130, t101 ^ 2, 0.2e1 * pkin(1) * t149 + 0.2e1 * t68 * t101, -0.2e1 * t69 * t101 - 0.2e1 * t93 * t171 (-t104 * t68 + t107 * t69) * t176, t93 * pkin(1) ^ 2 + t68 ^ 2 + t69 ^ 2, t67 ^ 2, t65 * t178, t138 * t178, t64, t65 * t130, t85, -0.2e1 * t138 * t32 + 0.2e1 * t56 * t65, 0.2e1 * t138 * t33 + 0.2e1 * t56 * t67, -0.2e1 * t32 * t67 - 0.2e1 * t33 * t65, t32 ^ 2 + t33 ^ 2 + t56 ^ 2, t46 ^ 2, t46 * t180, t46 * t179, t44 ^ 2, t65 * t180, t64, 0.2e1 * t11 * t65 + 0.2e1 * t29 * t44, -0.2e1 * t12 * t65 + 0.2e1 * t29 * t46, -0.2e1 * t11 * t46 - 0.2e1 * t12 * t44, t11 ^ 2 + t12 ^ 2 + t29 ^ 2, t23, -0.2e1 * t163, t17, t183, -0.2e1 * t158, t64, -0.2e1 * t123 * t65 + 0.2e1 * t16 * t24, 0.2e1 * t16 * t26 - 0.2e1 * t4 * t65, 0.2e1 * t123 * t26 - 0.2e1 * t4 * t24, t123 ^ 2 + t16 ^ 2 + t4 ^ 2, t23, t17, 0.2e1 * t163, t64, 0.2e1 * t158, t183, -0.2e1 * t2 * t65 + 0.2e1 * t5 * t24, -0.2e1 * t1 * t24 + 0.2e1 * t2 * t26, 0.2e1 * t1 * t65 - 0.2e1 * t5 * t26, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, 0, t138, t101, t68, -t69, 0, 0, t141, -t103 * t65 + t67 * t106, -t126, -t142, -t125, 0, -pkin(2) * t65 + pkin(9) * t126 - t56 * t106, -pkin(2) * t67 + pkin(9) * t125 + t56 * t103 (t141 - t142) * pkin(9) + t112, -t56 * pkin(2) + pkin(9) * t112, t46 * t133 (-t145 - t146) * t103, -t46 * t106 + t133 * t65, t44 * t136, t44 * t106 - t136 * t65, -t142, -t11 * t106 + t54 * t65 + (pkin(9) * t44 + t148) * t103, t12 * t106 - t55 * t65 + (pkin(9) * t46 + t147) * t103, -t55 * t44 - t54 * t46 + (-t102 * t12 - t105 * t11) * t103, t11 * t54 + t12 * t55 + t29 * t92, t15, -t118, t13, t165, t115, -t142, t106 * t123 - t120 * t65 + t16 * t60 + t76 * t24, t4 * t106 + t16 * t62 - t22 * t65 + t76 * t26, t120 * t26 + t123 * t62 - t22 * t24 - t4 * t60, t120 * t123 + t16 * t76 + t4 * t22, t15, t13, t118, -t142, -t115, t165, t2 * t106 - t20 * t65 + t31 * t24 + t5 * t60, -t1 * t60 - t19 * t24 + t2 * t62 + t20 * t26, -t1 * t106 + t19 * t65 - t31 * t26 - t5 * t62, t1 * t19 + t2 * t20 + t5 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t95, t128, 0, t97, 0, 0, pkin(2) * t174, pkin(2) * t175, 0.2e1 * (t95 + t97) * pkin(9), pkin(2) ^ 2 + t97 * t109 + t91, t96 * t95, -0.2e1 * t95 * t135, t132 * t175, t94 * t95, t102 * t128, t97, -0.2e1 * t54 * t106 + 0.2e1 * t168 * t95, 0.2e1 * t95 * pkin(9) * t105 + 0.2e1 * t55 * t106, 0.2e1 * (-t102 * t55 - t105 * t54) * t103, t54 ^ 2 + t55 ^ 2 + t91, t59, -0.2e1 * t159, t53, t182, 0.2e1 * t151, t97, 0.2e1 * t106 * t120 + 0.2e1 * t76 * t60, 0.2e1 * t22 * t106 + 0.2e1 * t76 * t62, 0.2e1 * t120 * t62 - 0.2e1 * t22 * t60, t120 ^ 2 + t22 ^ 2 + t76 ^ 2, t59, t53, 0.2e1 * t159, t97, -0.2e1 * t151, t182, 0.2e1 * t20 * t106 + 0.2e1 * t31 * t60, -0.2e1 * t19 * t60 + 0.2e1 * t20 * t62, -0.2e1 * t19 * t106 - 0.2e1 * t31 * t62, t19 ^ 2 + t20 ^ 2 + t31 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, -t65, -t138, t32, -t33, 0, 0, t145, -t102 * t44 + t46 * t105, t153, -t146, t152, 0, -pkin(3) * t44 - pkin(10) * t153 - t147, -pkin(3) * t46 - pkin(10) * t152 + t148 (t145 - t146) * pkin(10) + t114, -t29 * pkin(3) + pkin(10) * t114, t18, -t117, t43, t164, -t157, 0, t16 * t71 + t90 * t24 - t162, t16 * t73 + t90 * t26 - t161, t123 * t73 - t4 * t71 + t122, t123 * t49 + t16 * t90 + t4 * t51, t18, t43, t117, 0, t157, t164, t36 * t24 + t5 * t71 - t162, -t1 * t71 + t2 * t73 + t122, -t36 * t26 - t5 * t73 + t161, t1 * t51 + t2 * t49 + t5 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, 0, t106, 0, -t92, -t106 * pkin(9), 0, 0, t124 (-t94 + t96) * t103, -t134, -t124, -t132, 0, -pkin(9) * t133 + (-pkin(3) * t103 + pkin(10) * t106) * t102, pkin(10) * t132 + (t168 - t169) * t103, t113, -pkin(3) * t92 + pkin(10) * t113, t40, t116, -t140, t160, t150, 0, t90 * t60 + t76 * t71 + t144, t90 * t62 + t76 * t73 + t143, t120 * t73 - t22 * t71 + t121, t120 * t49 + t22 * t51 + t76 * t90, t40, -t140, -t116, 0, -t150, t160, t31 * t71 + t36 * t60 + t144, -t19 * t71 + t20 * t73 + t121, -t31 * t73 - t36 * t62 - t143, t19 * t51 + t20 * t49 + t31 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t94, 0.2e1 * t135, 0, t96, 0, 0, 0.2e1 * t169, -0.2e1 * pkin(3) * t102, 0.2e1 * t154 * pkin(10), pkin(10) ^ 2 * t154 + pkin(3) ^ 2, t70, -0.2e1 * t156, 0, t181, 0, 0, t71 * t177, t73 * t177, t111, t90 ^ 2 + t131, t70, 0, 0.2e1 * t156, 0, 0, t181, 0.2e1 * t36 * t71, t111, -0.2e1 * t36 * t73, t36 ^ 2 + t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, -t44, t65, t11, -t12, 0, 0, 0, 0, t26, 0, -t24, t65, t167 * t65 - t123, -t173 * t65 - t4 (-t100 * t26 - t24 * t98) * pkin(4) (-t100 * t123 + t4 * t98) * pkin(4), 0, t26, 0, t65, t24, 0, t172 * t65 - t123, -t84 * t24 - t88 * t26, t84 * t65 + t1, t1 * t84 - t2 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, 0, -t136, -t106, t54, -t55, 0, 0, 0, 0, t62, 0, -t60, -t106, -t100 * t166 - t120, t166 * t98 - t22 (-t100 * t62 - t60 * t98) * pkin(4) (-t100 * t120 + t22 * t98) * pkin(4), 0, t62, 0, -t106, t60, 0, -t106 * t172 - t120, -t84 * t60 - t88 * t62 (-qJ(6) - t84) * t106 + t22, t19 * t84 - t20 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, 0, t105, 0, -t102 * pkin(10), -t105 * pkin(10), 0, 0, 0, 0, t73, 0, -t71, 0, -t49, -t51 (-t100 * t73 - t71 * t98) * pkin(4) (-t100 * t49 + t51 * t98) * pkin(4), 0, t73, 0, 0, t71, 0, -t49, -t84 * t71 - t88 * t73, t51, -t49 * t88 + t51 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t167, -0.2e1 * t173, 0 (t100 ^ 2 + t98 ^ 2) * pkin(4) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t88, 0, 0.2e1 * t84, t84 ^ 2 + t88 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t26, 0, t16, 0, 0, 0, 0, 0, 0, t24, 0, -t26, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t62, 0, t76, 0, 0, 0, 0, 0, 0, t60, 0, -t62, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t73, 0, t90, 0, 0, 0, 0, 0, 0, t71, 0, -t73, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t26, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t62, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
