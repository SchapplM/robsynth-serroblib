% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRP10_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t102 = sin(qJ(5));
t171 = cos(qJ(5));
t100 = cos(pkin(11));
t106 = cos(qJ(2));
t99 = sin(pkin(6));
t135 = t99 * t106;
t101 = cos(pkin(6));
t103 = sin(qJ(3));
t105 = cos(qJ(3));
t104 = sin(qJ(2));
t136 = t99 * t104;
t68 = t101 * t103 + t105 * t136;
t98 = sin(pkin(11));
t46 = t100 * t135 + t68 * t98;
t48 = t100 * t68 - t135 * t98;
t24 = t102 * t48 + t171 * t46;
t182 = t24 ^ 2;
t77 = t100 * t102 + t171 * t98;
t62 = t77 * t103;
t181 = t62 ^ 2;
t66 = -t101 * t105 + t103 * t136;
t65 = t66 ^ 2;
t120 = t171 * t100;
t75 = t102 * t98 - t120;
t180 = t75 ^ 2;
t179 = -0.2e1 * t46;
t178 = 0.2e1 * t66;
t177 = -0.2e1 * t68;
t90 = -pkin(4) * t100 - pkin(3);
t176 = 0.2e1 * t90;
t175 = 0.2e1 * t99;
t169 = pkin(1) * t106;
t84 = pkin(8) * t136;
t58 = t84 + (-pkin(2) - t169) * t101;
t30 = t66 * pkin(3) - t68 * qJ(4) + t58;
t126 = pkin(8) * t135;
t170 = pkin(1) * t104;
t59 = t126 + (pkin(9) + t170) * t101;
t60 = (-pkin(2) * t106 - pkin(9) * t104 - pkin(1)) * t99;
t36 = t103 * t60 + t105 * t59;
t31 = -qJ(4) * t135 + t36;
t12 = t100 * t31 + t30 * t98;
t10 = -pkin(10) * t46 + t12;
t11 = t100 * t30 - t31 * t98;
t8 = pkin(4) * t66 - pkin(10) * t48 + t11;
t4 = t10 * t171 + t102 * t8;
t174 = pkin(3) * t98;
t173 = pkin(9) * t98;
t172 = t66 * pkin(5);
t168 = pkin(9) * t100;
t92 = t103 * pkin(9);
t167 = t105 * pkin(5);
t166 = t24 * t62;
t165 = t24 * t75;
t26 = -t102 * t46 + t171 * t48;
t164 = t26 * t24;
t35 = -t103 * t59 + t105 * t60;
t32 = pkin(3) * t135 - t35;
t163 = t32 * t98;
t162 = t48 * t98;
t153 = pkin(10) + qJ(4);
t121 = t153 * t98;
t80 = t153 * t100;
t51 = t102 * t80 + t121 * t171;
t161 = t51 * t66;
t53 = -t102 * t121 + t171 * t80;
t160 = t53 * t66;
t159 = t62 * t75;
t138 = t98 * t103;
t64 = -t102 * t138 + t103 * t120;
t158 = t64 * t62;
t157 = t66 * t24;
t156 = t66 * t75;
t155 = t77 * t75;
t154 = t98 * t66;
t133 = t100 * t103;
t79 = -pkin(3) * t105 - qJ(4) * t103 - pkin(2);
t72 = t100 * t79;
t41 = -pkin(10) * t133 + t72 + (-pkin(4) - t173) * t105;
t132 = t100 * t105;
t57 = pkin(9) * t132 + t79 * t98;
t49 = -pkin(10) * t138 + t57;
t22 = t102 * t41 + t171 * t49;
t78 = pkin(4) * t138 + t92;
t93 = t98 ^ 2;
t95 = t100 ^ 2;
t152 = t93 + t95;
t151 = t100 * t66;
t150 = t105 * t62;
t149 = t105 * t75;
t94 = t99 ^ 2;
t148 = t106 * t94;
t147 = t32 * t100;
t146 = t46 * t100;
t145 = t51 * t105;
t144 = t53 * t105;
t143 = t66 * qJ(6);
t142 = t66 * t105;
t141 = t68 * t103;
t140 = t77 * t105;
t139 = t98 * t100;
t137 = t98 * t105;
t134 = qJ(4) * t105;
t131 = t101 * t104;
t130 = t103 * t105;
t129 = t105 * qJ(6);
t128 = t51 ^ 2 + t53 ^ 2;
t127 = 0.2e1 * t135;
t125 = 0.2e1 * t130;
t124 = t98 * t133;
t123 = t103 * t135;
t122 = t105 * t135;
t119 = t10 * t102 - t171 * t8;
t118 = -t24 * t53 + t26 * t51;
t117 = t51 * t64 - t53 * t62;
t21 = -t102 * t49 + t171 * t41;
t116 = t24 * t64 + t26 * t62;
t115 = t24 * t77 + t26 * t75;
t114 = -t62 * t77 - t64 * t75;
t113 = t12 * t100 - t11 * t98;
t56 = -pkin(9) * t137 + t72;
t112 = t100 * t57 - t56 * t98;
t111 = t105 * t24 - t62 * t66;
t110 = -t35 * t103 + t36 * t105;
t109 = 0.2e1 * t51 * t77 - 0.2e1 * t53 * t75;
t16 = pkin(4) * t46 + t32;
t108 = pkin(9) ^ 2;
t97 = t105 ^ 2;
t96 = t103 ^ 2;
t91 = t96 * t108;
t87 = t94 * t106 ^ 2;
t73 = t77 ^ 2;
t70 = pkin(1) * t131 + t126;
t69 = t101 * t169 - t84;
t61 = t64 ^ 2;
t55 = -0.2e1 * t64 * t105;
t45 = t77 * t66;
t43 = t64 * t77;
t40 = pkin(5) * t75 - qJ(6) * t77 + t90;
t33 = pkin(5) * t62 - qJ(6) * t64 + t78;
t23 = t26 ^ 2;
t20 = -t21 + t167;
t19 = -t129 + t22;
t18 = t26 * t77;
t17 = t26 * t178;
t15 = t26 * t64;
t13 = -t105 * t26 + t64 * t66;
t5 = pkin(5) * t24 - qJ(6) * t26 + t16;
t2 = t119 - t172;
t1 = t143 + t4;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t94 * t104 ^ 2, 0.2e1 * t104 * t148, t131 * t175, t87, t101 * t127, t101 ^ 2, 0.2e1 * pkin(1) * t148 + 0.2e1 * t101 * t69, -0.2e1 * t101 * t70 - 0.2e1 * t170 * t94 (-t104 * t69 + t106 * t70) * t175, pkin(1) ^ 2 * t94 + t69 ^ 2 + t70 ^ 2, t68 ^ 2, t66 * t177, t135 * t177, t65, t66 * t127, t87, -0.2e1 * t135 * t35 + 0.2e1 * t58 * t66, 0.2e1 * t135 * t36 + 0.2e1 * t58 * t68, -0.2e1 * t35 * t68 - 0.2e1 * t36 * t66, t35 ^ 2 + t36 ^ 2 + t58 ^ 2, t48 ^ 2, t48 * t179, t48 * t178, t46 ^ 2, t66 * t179, t65, 0.2e1 * t11 * t66 + 0.2e1 * t32 * t46, -0.2e1 * t12 * t66 + 0.2e1 * t32 * t48, -0.2e1 * t11 * t48 - 0.2e1 * t12 * t46, t11 ^ 2 + t12 ^ 2 + t32 ^ 2, t23, -0.2e1 * t164, t17, t182, -0.2e1 * t157, t65, -0.2e1 * t119 * t66 + 0.2e1 * t16 * t24, 0.2e1 * t16 * t26 - 0.2e1 * t4 * t66, 0.2e1 * t119 * t26 - 0.2e1 * t24 * t4, t119 ^ 2 + t16 ^ 2 + t4 ^ 2, t23, t17, 0.2e1 * t164, t65, 0.2e1 * t157, t182, -0.2e1 * t2 * t66 + 0.2e1 * t24 * t5, -0.2e1 * t1 * t24 + 0.2e1 * t2 * t26, 0.2e1 * t1 * t66 - 0.2e1 * t26 * t5, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, 0, t135, t101, t69, -t70, 0, 0, t141, -t103 * t66 + t105 * t68, -t123, -t142, -t122, 0, -pkin(2) * t66 + pkin(9) * t123 - t105 * t58, -pkin(2) * t68 + pkin(9) * t122 + t103 * t58 (t141 - t142) * pkin(9) + t110, -t58 * pkin(2) + pkin(9) * t110, t48 * t133 (-t146 - t162) * t103, -t105 * t48 + t133 * t66, t46 * t138, t105 * t46 - t138 * t66, -t142, -t11 * t105 + t56 * t66 + (pkin(9) * t46 + t163) * t103, t12 * t105 - t57 * t66 + (pkin(9) * t48 + t147) * t103, -t57 * t46 - t56 * t48 + (-t100 * t11 - t12 * t98) * t103, t11 * t56 + t12 * t57 + t32 * t92, t15, -t116, t13, t166, t111, -t142, t105 * t119 + t16 * t62 + t21 * t66 + t24 * t78, t105 * t4 + t16 * t64 - t22 * t66 + t26 * t78, t119 * t64 - t21 * t26 - t22 * t24 - t4 * t62, -t119 * t21 + t16 * t78 + t22 * t4, t15, t13, t116, -t142, -t111, t166, t105 * t2 - t20 * t66 + t24 * t33 + t5 * t62, -t1 * t62 - t19 * t24 + t2 * t64 + t20 * t26, -t1 * t105 + t19 * t66 - t26 * t33 - t5 * t64, t1 * t19 + t2 * t20 + t33 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t96, t125, 0, t97, 0, 0, 0.2e1 * pkin(2) * t105, -0.2e1 * pkin(2) * t103, 0.2e1 * (t96 + t97) * pkin(9), pkin(2) ^ 2 + t108 * t97 + t91, t95 * t96, -0.2e1 * t96 * t139, -0.2e1 * t100 * t130, t93 * t96, t98 * t125, t97, -0.2e1 * t105 * t56 + 0.2e1 * t173 * t96, 0.2e1 * t105 * t57 + 0.2e1 * t168 * t96, 0.2e1 * (-t100 * t56 - t57 * t98) * t103, t56 ^ 2 + t57 ^ 2 + t91, t61, -0.2e1 * t158, t55, t181, 0.2e1 * t150, t97, -0.2e1 * t105 * t21 + 0.2e1 * t62 * t78, 0.2e1 * t105 * t22 + 0.2e1 * t64 * t78, -0.2e1 * t21 * t64 - 0.2e1 * t22 * t62, t21 ^ 2 + t22 ^ 2 + t78 ^ 2, t61, t55, 0.2e1 * t158, t97, -0.2e1 * t150, t181, 0.2e1 * t105 * t20 + 0.2e1 * t33 * t62, -0.2e1 * t19 * t62 + 0.2e1 * t20 * t64, -0.2e1 * t105 * t19 - 0.2e1 * t33 * t64, t19 ^ 2 + t20 ^ 2 + t33 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, -t66, -t135, t35, -t36, 0, 0, t162, t100 * t48 - t46 * t98, t154, -t146, t151, 0, -pkin(3) * t46 - qJ(4) * t154 - t147, -pkin(3) * t48 - qJ(4) * t151 + t163 (-t146 + t162) * qJ(4) + t113, -t32 * pkin(3) + qJ(4) * t113, t18, -t115, t45, t165, -t156, 0, t16 * t75 + t24 * t90 - t161, t16 * t77 + t26 * t90 - t160, t119 * t77 - t4 * t75 + t118, t119 * t51 + t16 * t90 + t4 * t53, t18, t45, t115, 0, t156, t165, t24 * t40 + t5 * t75 - t161, -t1 * t75 + t2 * t77 + t118, -t26 * t40 - t5 * t77 + t160, t1 * t53 + t2 * t51 + t40 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, 0, t105, 0, -t92, -t105 * pkin(9), 0, 0, t124 (-t93 + t95) * t103, -t137, -t124, -t132, 0, t98 * t134 + (-t168 - t174) * t103, pkin(9) * t138 + (-pkin(3) * t103 + t134) * t100, t112, -pkin(3) * t92 + qJ(4) * t112, t43, t114, -t140, t159, t149, 0, t62 * t90 + t75 * t78 + t145, t64 * t90 + t77 * t78 + t144, -t21 * t77 - t22 * t75 + t117, -t21 * t51 + t22 * t53 + t78 * t90, t43, -t140, -t114, 0, -t149, t159, t33 * t75 + t40 * t62 + t145, -t19 * t75 + t20 * t77 + t117, -t33 * t77 - t40 * t64 - t144, t19 * t53 + t20 * t51 + t33 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t93, 0.2e1 * t139, 0, t95, 0, 0, 0.2e1 * pkin(3) * t100, -0.2e1 * t174, 0.2e1 * t152 * qJ(4), qJ(4) ^ 2 * t152 + pkin(3) ^ 2, t73, -0.2e1 * t155, 0, t180, 0, 0, t75 * t176, t77 * t176, t109, t90 ^ 2 + t128, t73, 0, 0.2e1 * t155, 0, 0, t180, 0.2e1 * t40 * t75, t109, -0.2e1 * t40 * t77, t40 ^ 2 + t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t48, 0, t32, 0, 0, 0, 0, 0, 0, t24, t26, 0, t16, 0, 0, 0, 0, 0, 0, t24, 0, -t26, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, t133, 0, t92, 0, 0, 0, 0, 0, 0, t62, t64, 0, t78, 0, 0, 0, 0, 0, 0, t62, 0, -t64, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, t98, 0, -pkin(3), 0, 0, 0, 0, 0, 0, t75, t77, 0, t90, 0, 0, 0, 0, 0, 0, t75, 0, -t77, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, -t24, t66, -t119, -t4, 0, 0, 0, t26, 0, t66, t24, 0, -t119 + 0.2e1 * t172, -pkin(5) * t26 - qJ(6) * t24, 0.2e1 * t143 + t4, -pkin(5) * t2 + qJ(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, -t62, -t105, t21, -t22, 0, 0, 0, t64, 0, -t105, t62, 0, t21 - 0.2e1 * t167, -pkin(5) * t64 - qJ(6) * t62, -0.2e1 * t129 + t22, -pkin(5) * t20 + qJ(6) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, 0, -t75, 0, -t51, -t53, 0, 0, 0, t77, 0, 0, t75, 0, -t51, -pkin(5) * t77 - qJ(6) * t75, t53, -pkin(5) * t51 + qJ(6) * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, t26, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, t64, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;