% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 19:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRR11_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t99 = sin(qJ(3));
t91 = t99 ^ 2;
t103 = cos(qJ(3));
t94 = t103 ^ 2;
t174 = t91 + t94;
t97 = sin(qJ(6));
t89 = t97 ^ 2;
t101 = cos(qJ(6));
t92 = t101 ^ 2;
t144 = t89 + t92;
t173 = t144 * pkin(11);
t102 = cos(qJ(5));
t162 = -pkin(3) - pkin(4);
t98 = sin(qJ(5));
t67 = t102 * qJ(4) + t98 * t162;
t63 = -pkin(11) + t67;
t115 = t144 * t63;
t100 = sin(qJ(2));
t95 = sin(pkin(6));
t129 = t95 * t100;
t96 = cos(pkin(6));
t49 = -t96 * t103 + t99 * t129;
t51 = t103 * t129 + t96 * t99;
t27 = -t49 * t102 + t51 * t98;
t172 = t27 ^ 2;
t118 = (pkin(9) - pkin(10)) * t99;
t85 = t103 * pkin(9);
t70 = -t103 * pkin(10) + t85;
t31 = -t102 * t118 + t98 * t70;
t171 = t31 ^ 2;
t170 = t49 ^ 2;
t56 = t102 * t103 + t98 * t99;
t169 = t56 ^ 2;
t168 = -0.2e1 * t27;
t167 = 0.2e1 * t56;
t59 = t102 * t99 - t98 * t103;
t166 = 0.2e1 * t59;
t165 = -0.2e1 * t97;
t164 = -0.2e1 * t99;
t163 = 0.2e1 * t101;
t161 = pkin(9) * t51;
t104 = cos(qJ(2));
t77 = t95 * t104;
t121 = pkin(8) * t77;
t157 = pkin(1) * t100;
t41 = t121 + (pkin(9) + t157) * t96;
t42 = (-pkin(2) * t104 - pkin(9) * t100 - pkin(1)) * t95;
t146 = -t103 * t42 + t99 * t41;
t73 = pkin(3) * t77;
t20 = t73 + t146;
t11 = pkin(4) * t77 - t51 * pkin(10) + t20;
t117 = qJ(4) * t77;
t25 = t103 * t41 + t99 * t42;
t19 = -t117 + t25;
t17 = t49 * pkin(10) + t19;
t6 = t102 * t11 - t98 * t17;
t3 = -pkin(5) * t77 - t6;
t160 = t3 * t97;
t159 = t99 * pkin(9);
t65 = t98 * qJ(4) - t102 * t162;
t62 = pkin(5) + t65;
t158 = pkin(5) + t62;
t29 = t51 * t102 + t49 * t98;
t23 = t29 * t101 + t97 * t77;
t156 = t23 * t97;
t155 = t27 * t56;
t154 = t3 * t101;
t153 = t31 * t97;
t152 = t51 * t49;
t151 = t97 * t27;
t150 = t97 * t56;
t149 = t97 * t59;
t148 = t97 * t63;
t147 = t98 * t27;
t7 = t102 * t17 + t98 * t11;
t52 = t96 * t104 * pkin(1) - pkin(8) * t129;
t145 = t174 * pkin(9) ^ 2;
t143 = t101 * t27;
t142 = t101 * t56;
t141 = t101 * t59;
t140 = t101 * t63;
t139 = t101 * t98;
t138 = t102 * t97;
t88 = t95 ^ 2;
t137 = t104 * t88;
t136 = t19 * t103;
t21 = -t101 * t77 + t29 * t97;
t135 = t21 * t101;
t134 = t23 * t101;
t133 = t25 * t103;
t132 = t31 * t101;
t131 = t31 * t102;
t130 = t49 * t103;
t128 = t97 * t101;
t127 = t99 * t103;
t126 = t102 * t101;
t125 = -0.2e1 * t59 * t56;
t124 = -0.2e1 * t77;
t123 = 0.2e1 * t77;
t122 = -0.2e1 * t128;
t68 = -t103 * pkin(3) - t99 * qJ(4) - pkin(2);
t40 = -t96 * pkin(2) - t52;
t120 = t49 * t77;
t119 = t99 * t77;
t116 = t103 * t77;
t58 = t144 * t98;
t54 = t103 * pkin(4) - t68;
t114 = pkin(9) * t116;
t18 = t49 * pkin(3) - t51 * qJ(4) + t40;
t113 = -pkin(5) * t59 - pkin(11) * t56;
t4 = pkin(11) * t77 + t7;
t16 = -t49 * pkin(4) - t18;
t8 = t27 * pkin(5) - t29 * pkin(11) + t16;
t1 = t101 * t8 - t97 * t4;
t2 = t101 * t4 + t97 * t8;
t112 = -t1 * t97 + t2 * t101;
t111 = -t56 * t63 + t59 * t62;
t110 = -pkin(3) * t99 + t103 * qJ(4);
t26 = t56 * pkin(5) - t59 * pkin(11) + t54;
t33 = t102 * t70 + t98 * t118;
t12 = t101 * t26 - t97 * t33;
t13 = t101 * t33 + t97 * t26;
t5 = t13 * t101 - t12 * t97;
t109 = -t135 + t156;
t108 = -t102 * t59 - t98 * t56;
t107 = t51 * t103 - t99 * t49;
t93 = t102 ^ 2;
t90 = t98 ^ 2;
t76 = t88 * t104 ^ 2;
t75 = 0.2e1 * t128;
t69 = pkin(9) * t119;
t64 = 0.2e1 * t174 * pkin(9);
t55 = t59 ^ 2;
t53 = t96 * t157 + t121;
t48 = t51 ^ 2;
t45 = t59 * t128;
t44 = t51 * t99;
t39 = pkin(9) * t130;
t34 = t51 * t124;
t30 = (t89 - t92) * t59;
t9 = -t97 * t21 + t134;
t10 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t88 * t100 ^ 2, 0.2e1 * t100 * t137, 0.2e1 * t96 * t129, t76, t96 * t123, t96 ^ 2, 0.2e1 * pkin(1) * t137 + 0.2e1 * t52 * t96, -0.2e1 * t88 * t157 - 0.2e1 * t53 * t96, 0.2e1 * (-t100 * t52 + t104 * t53) * t95, t88 * pkin(1) ^ 2 + t52 ^ 2 + t53 ^ 2, t48, -0.2e1 * t152, t34, t170, 0.2e1 * t120, t76, 0.2e1 * t146 * t77 + 0.2e1 * t40 * t49, 0.2e1 * t25 * t77 + 0.2e1 * t40 * t51, 0.2e1 * t146 * t51 - 0.2e1 * t25 * t49, t146 ^ 2 + t25 ^ 2 + t40 ^ 2, t48, t34, 0.2e1 * t152, t76, -0.2e1 * t120, t170, 0.2e1 * t18 * t49 + 0.2e1 * t20 * t77, -0.2e1 * t19 * t49 + 0.2e1 * t20 * t51, -0.2e1 * t18 * t51 - 0.2e1 * t19 * t77, t18 ^ 2 + t19 ^ 2 + t20 ^ 2, t29 ^ 2, t29 * t168, t29 * t123, t172, t27 * t124, t76, 0.2e1 * t16 * t27 + 0.2e1 * t6 * t77, 0.2e1 * t16 * t29 - 0.2e1 * t7 * t77, -0.2e1 * t7 * t27 - 0.2e1 * t6 * t29, t16 ^ 2 + t6 ^ 2 + t7 ^ 2, t23 ^ 2, -0.2e1 * t23 * t21, 0.2e1 * t23 * t27, t21 ^ 2, t21 * t168, t172, 0.2e1 * t1 * t27 + 0.2e1 * t3 * t21, -0.2e1 * t2 * t27 + 0.2e1 * t3 * t23, -0.2e1 * t1 * t23 - 0.2e1 * t2 * t21, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, 0, t77, t96, t52, -t53, 0, 0, t44, t107, -t119, -t130, -t116, 0, -pkin(2) * t49 - t40 * t103 + t69, -pkin(2) * t51 + t40 * t99 + t114, t133 - t39 + (t146 + t161) * t99, -t40 * pkin(2) + (t146 * t99 + t133) * pkin(9), t44, -t119, -t107, 0, t116, -t130, -t18 * t103 + t68 * t49 + t69, t136 - t39 + (t20 + t161) * t99, -t18 * t99 - t68 * t51 - t114, t18 * t68 + (t20 * t99 + t136) * pkin(9), t29 * t59, -t59 * t27 - t29 * t56, t59 * t77, t155, -t56 * t77, 0, t16 * t56 + t54 * t27 - t31 * t77, t16 * t59 + t54 * t29 - t33 * t77, -t33 * t27 + t31 * t29 - t7 * t56 - t6 * t59, t16 * t54 - t6 * t31 + t7 * t33, t59 * t134 (-t135 - t156) * t59, t141 * t27 + t23 * t56, t21 * t149, -t149 * t27 - t21 * t56, t155, t1 * t56 + t12 * t27 + t149 * t3 + t31 * t21, -t13 * t27 + t141 * t3 - t2 * t56 + t31 * t23, -t12 * t23 - t13 * t21 + (-t1 * t101 - t2 * t97) * t59, t1 * t12 + t2 * t13 + t3 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t91, 0.2e1 * t127, 0, t94, 0, 0, 0.2e1 * pkin(2) * t103, pkin(2) * t164, t64, pkin(2) ^ 2 + t145, t91, 0, -0.2e1 * t127, 0, 0, t94, -0.2e1 * t68 * t103, t64, t68 * t164, t68 ^ 2 + t145, t55, t125, 0, t169, 0, 0, t54 * t167, t54 * t166, 0.2e1 * t31 * t59 - 0.2e1 * t33 * t56, t33 ^ 2 + t54 ^ 2 + t171, t92 * t55, t55 * t122, t141 * t167, t89 * t55, t97 * t125, t169, 0.2e1 * t12 * t56 + 0.2e1 * t149 * t31, -0.2e1 * t13 * t56 + 0.2e1 * t132 * t59 (-t101 * t12 - t13 * t97) * t166, t12 ^ 2 + t13 ^ 2 + t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, -t49, -t77, -t146, -t25, 0, 0, 0, t51, 0, -t77, t49, 0, -0.2e1 * t73 - t146, -t51 * pkin(3) - t49 * qJ(4), -0.2e1 * t117 + t25, -t20 * pkin(3) + t19 * qJ(4), 0, 0, -t29, 0, t27, -t77, -t65 * t77 - t6, -t67 * t77 + t7, -t67 * t27 + t65 * t29, -t6 * t65 + t7 * t67, -t156, -t9, -t151, t135, -t143, 0, -t148 * t27 + t62 * t21 + t154, -t140 * t27 + t62 * t23 - t160 (t23 * t63 + t1) * t97 + (-t21 * t63 - t2) * t101, t112 * t63 + t3 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, 0, t103, 0, -t159, -t85, 0, 0, 0, t99, 0, 0, -t103, 0, -t159, t110, t85, t110 * pkin(9), 0, 0, -t59, 0, t56, 0, t31, t33, -t67 * t56 + t65 * t59, t31 * t65 + t33 * t67, -t45, t30, -t150, t45, -t142, 0, t111 * t97 + t132, t101 * t111 - t153, -t5, t31 * t62 + t5 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3), 0, 0.2e1 * qJ(4), pkin(3) ^ 2 + qJ(4) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t65, 0.2e1 * t67, 0, t65 ^ 2 + t67 ^ 2, t89, t75, 0, t92, 0, 0, t62 * t163, t62 * t165, -0.2e1 * t115, t144 * t63 ^ 2 + t62 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t51, 0, t20, 0, 0, 0, 0, 0, 0, t102 * t77, -t98 * t77, -t102 * t29 - t147, t6 * t102 + t7 * t98, 0, 0, 0, 0, 0, 0, -t102 * t21 - t147 * t97, -t102 * t23 - t139 * t27, t109 * t98, -t3 * t102 + t112 * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, 0, t159, 0, 0, 0, 0, 0, 0, 0, 0, t108, t33 * t98 - t131, 0, 0, 0, 0, 0, 0, t108 * t97, t108 * t101, 0, t5 * t98 - t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -t102, t98, 0, -t65 * t102 + t67 * t98, 0, 0, 0, 0, 0, 0, -t126, t138, -t58, -t62 * t102 + t58 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90 + t93, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144 * t90 + t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, -t27, t77, t6, -t7, 0, 0, t156, t9, t151, -t135, t143, 0, -pkin(5) * t21 - pkin(11) * t151 - t154, -pkin(5) * t23 - pkin(11) * t143 + t160, pkin(11) * t109 + t112, -t3 * pkin(5) + pkin(11) * t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, -t56, 0, -t31, -t33, 0, 0, t45, -t30, t150, -t45, t142, 0, t113 * t97 - t132, t101 * t113 + t153, t5, -t31 * pkin(5) + pkin(11) * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t65, -t67, 0, 0, -t89, t122, 0, -t92, 0, 0, -t158 * t101, t158 * t97, t115 - t173, -t62 * pkin(5) + pkin(11) * t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, -t98, 0, 0, 0, 0, 0, 0, 0, 0, t126, -t138, t58, t102 * pkin(5) + pkin(11) * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t89, t75, 0, t92, 0, 0, pkin(5) * t163, pkin(5) * t165, 0.2e1 * t173, pkin(11) ^ 2 * t144 + pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, -t21, t27, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, 0, -t149, t56, t12, -t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, 0, -t101, 0, -t148, -t140, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97 * t98, -t139, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, t101, 0, -t97 * pkin(11), -t101 * pkin(11), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t10;
