% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRRP9_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t116 = cos(pkin(6));
t119 = sin(qJ(3));
t123 = cos(qJ(3));
t115 = sin(pkin(6));
t120 = sin(qJ(2));
t149 = t115 * t120;
t77 = -t116 * t123 + t119 * t149;
t76 = t77 ^ 2;
t117 = sin(qJ(5));
t121 = cos(qJ(5));
t118 = sin(qJ(4));
t122 = cos(qJ(4));
t124 = cos(qJ(2));
t148 = t115 * t124;
t79 = t116 * t119 + t123 * t149;
t51 = -t79 * t118 - t122 * t148;
t52 = -t118 * t148 + t79 * t122;
t30 = t117 * t52 - t121 * t51;
t177 = -0.2e1 * t30;
t104 = -t122 * pkin(4) - pkin(3);
t85 = t117 * t118 - t121 * t122;
t64 = t85 * pkin(5) + t104;
t176 = 0.2e1 * t64;
t175 = 0.2e1 * t77;
t174 = -0.2e1 * t79;
t173 = 0.2e1 * t104;
t172 = 0.2e1 * t115;
t171 = -0.2e1 * t119;
t170 = -0.2e1 * t123;
t169 = 0.2e1 * t123;
t168 = -pkin(11) - pkin(10);
t167 = pkin(1) * t120;
t166 = pkin(1) * t124;
t165 = pkin(3) * t122;
t164 = pkin(9) * t118;
t112 = t119 ^ 2;
t163 = t112 * pkin(9);
t162 = t117 * pkin(4);
t108 = t119 * pkin(9);
t109 = t121 * pkin(4);
t161 = t123 * pkin(5);
t160 = t118 * t77;
t159 = t122 * t77;
t137 = pkin(8) * t148;
t66 = t137 + (pkin(9) + t167) * t116;
t67 = (-pkin(2) * t124 - pkin(9) * t120 - pkin(1)) * t115;
t38 = -t119 * t66 + t123 * t67;
t36 = pkin(3) * t148 - t38;
t158 = t36 * t118;
t157 = t36 * t122;
t156 = t51 * t122;
t155 = t52 * t118;
t154 = t77 * t123;
t153 = t79 * t119;
t152 = t85 * t123;
t87 = t117 * t122 + t121 * t118;
t151 = t87 * t123;
t146 = t118 * t119;
t89 = pkin(4) * t146 + t108;
t110 = t115 ^ 2;
t150 = t110 * t124;
t147 = t116 * t120;
t145 = t118 * t122;
t144 = t118 * t123;
t143 = t122 * t119;
t142 = t122 * t123;
t111 = t118 ^ 2;
t113 = t122 ^ 2;
t141 = t111 + t113;
t140 = 0.2e1 * t148;
t139 = t119 * t169;
t138 = pkin(9) * t142;
t136 = t119 * t148;
t135 = t123 * t148;
t134 = t118 * t143;
t96 = pkin(8) * t149;
t65 = t96 + (-pkin(2) - t166) * t116;
t35 = t77 * pkin(3) - t79 * pkin(10) + t65;
t39 = t119 * t67 + t123 * t66;
t37 = -pkin(10) * t148 + t39;
t13 = t118 * t35 + t122 * t37;
t10 = t51 * pkin(11) + t13;
t12 = -t118 * t37 + t122 * t35;
t7 = t77 * pkin(4) - t52 * pkin(11) + t12;
t3 = -t117 * t10 + t121 * t7;
t90 = -t123 * pkin(3) - t119 * pkin(10) - pkin(2);
t82 = t122 * t90;
t45 = -pkin(11) * t143 + t82 + (-pkin(4) - t164) * t123;
t53 = t138 + (-pkin(11) * t119 + t90) * t118;
t26 = -t117 * t53 + t121 * t45;
t91 = t168 * t118;
t92 = t168 * t122;
t55 = t117 * t92 + t121 * t91;
t4 = t121 * t10 + t117 * t7;
t27 = t117 * t45 + t121 * t53;
t56 = t117 * t91 - t121 * t92;
t133 = -t12 * t118 + t13 * t122;
t61 = -pkin(9) * t144 + t82;
t62 = t118 * t90 + t138;
t132 = -t61 * t118 + t62 * t122;
t131 = -t38 * t119 + t39 * t123;
t32 = t117 * t51 + t121 * t52;
t130 = -t32 * qJ(6) + t3;
t73 = -t117 * t146 + t121 * t143;
t129 = -t73 * qJ(6) + t26;
t18 = -t51 * pkin(4) + t36;
t74 = t77 * pkin(5);
t1 = t130 + t74;
t2 = -t30 * qJ(6) + t4;
t71 = t87 * t119;
t17 = -t71 * qJ(6) + t27;
t128 = pkin(4) ^ 2;
t127 = pkin(9) ^ 2;
t125 = 0.2e1 * pkin(5);
t114 = t123 ^ 2;
t107 = t112 * t127;
t106 = t117 ^ 2 * t128;
t105 = -0.2e1 * t162;
t103 = t109 + pkin(5);
t100 = t110 * t124 ^ 2;
t98 = t123 * t162;
t84 = t87 ^ 2;
t83 = t85 ^ 2;
t81 = pkin(1) * t147 + t137;
t80 = t116 * t166 - t96;
t75 = t85 * t162;
t70 = t73 ^ 2;
t69 = t71 ^ 2;
t63 = t77 * t162;
t60 = t73 * t170;
t59 = t71 * t170;
t58 = t71 * t162;
t54 = -0.2e1 * t87 * t85;
t50 = t87 * t77;
t49 = t85 * t77;
t48 = t73 * t87;
t47 = t71 * t85;
t46 = t71 * pkin(5) + t89;
t42 = -0.2e1 * t73 * t71;
t41 = -t85 * qJ(6) + t56;
t40 = -t87 * qJ(6) + t55;
t33 = -t87 * t71 - t73 * t85;
t29 = t32 ^ 2;
t28 = t30 ^ 2;
t25 = t30 * t162;
t24 = t32 * t87;
t23 = t30 * t85;
t22 = t32 * t175;
t21 = t77 * t177;
t20 = t32 * t73;
t19 = t30 * t71;
t16 = t129 - t161;
t15 = -t32 * t123 + t73 * t77;
t14 = t30 * t123 - t71 * t77;
t11 = t32 * t177;
t9 = t30 * pkin(5) + t18;
t8 = -t87 * t30 - t32 * t85;
t5 = -t73 * t30 - t32 * t71;
t6 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t110 * t120 ^ 2, 0.2e1 * t120 * t150, t147 * t172, t100, t116 * t140, t116 ^ 2, 0.2e1 * pkin(1) * t150 + 0.2e1 * t80 * t116, -0.2e1 * t110 * t167 - 0.2e1 * t81 * t116 (-t120 * t80 + t124 * t81) * t172, t110 * pkin(1) ^ 2 + t80 ^ 2 + t81 ^ 2, t79 ^ 2, t77 * t174, t148 * t174, t76, t77 * t140, t100, -0.2e1 * t38 * t148 + 0.2e1 * t65 * t77, 0.2e1 * t39 * t148 + 0.2e1 * t65 * t79, -0.2e1 * t38 * t79 - 0.2e1 * t39 * t77, t38 ^ 2 + t39 ^ 2 + t65 ^ 2, t52 ^ 2, 0.2e1 * t52 * t51, t52 * t175, t51 ^ 2, t51 * t175, t76, 0.2e1 * t12 * t77 - 0.2e1 * t36 * t51, -0.2e1 * t13 * t77 + 0.2e1 * t36 * t52, -0.2e1 * t12 * t52 + 0.2e1 * t13 * t51, t12 ^ 2 + t13 ^ 2 + t36 ^ 2, t29, t11, t22, t28, t21, t76, 0.2e1 * t18 * t30 + 0.2e1 * t3 * t77, 0.2e1 * t18 * t32 - 0.2e1 * t4 * t77, -0.2e1 * t3 * t32 - 0.2e1 * t4 * t30, t18 ^ 2 + t3 ^ 2 + t4 ^ 2, t29, t11, t22, t28, t21, t76, 0.2e1 * t1 * t77 + 0.2e1 * t9 * t30, -0.2e1 * t2 * t77 + 0.2e1 * t9 * t32, -0.2e1 * t1 * t32 - 0.2e1 * t2 * t30, t1 ^ 2 + t2 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149, 0, t148, t116, t80, -t81, 0, 0, t153, -t119 * t77 + t79 * t123, -t136, -t154, -t135, 0, -pkin(2) * t77 + pkin(9) * t136 - t65 * t123, -pkin(2) * t79 + pkin(9) * t135 + t65 * t119 (t153 - t154) * pkin(9) + t131, -t65 * pkin(2) + pkin(9) * t131, t52 * t143 (-t155 + t156) * t119, -t52 * t123 + t77 * t143, -t51 * t146, -t51 * t123 - t77 * t146, -t154, -t12 * t123 + t61 * t77 + (-pkin(9) * t51 + t158) * t119, t13 * t123 - t62 * t77 + (pkin(9) * t52 + t157) * t119, t62 * t51 - t61 * t52 + (-t118 * t13 - t12 * t122) * t119, t36 * t108 + t12 * t61 + t13 * t62, t20, t5, t15, t19, t14, -t154, -t3 * t123 + t18 * t71 + t26 * t77 + t89 * t30, t4 * t123 + t18 * t73 - t27 * t77 + t89 * t32, -t26 * t32 - t27 * t30 - t3 * t73 - t4 * t71, t18 * t89 + t3 * t26 + t4 * t27, t20, t5, t15, t19, t14, -t154, -t1 * t123 + t16 * t77 + t46 * t30 + t9 * t71, t2 * t123 - t17 * t77 + t46 * t32 + t9 * t73, -t1 * t73 - t16 * t32 - t17 * t30 - t2 * t71, t1 * t16 + t2 * t17 + t9 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t112, t139, 0, t114, 0, 0, pkin(2) * t169, pkin(2) * t171, 0.2e1 * (t112 + t114) * pkin(9), pkin(2) ^ 2 + t114 * t127 + t107, t113 * t112, -0.2e1 * t112 * t145, t142 * t171, t111 * t112, t118 * t139, t114, 0.2e1 * t118 * t163 - 0.2e1 * t61 * t123, 0.2e1 * t122 * t163 + 0.2e1 * t62 * t123, 0.2e1 * (-t118 * t62 - t122 * t61) * t119, t61 ^ 2 + t62 ^ 2 + t107, t70, t42, t60, t69, -t59, t114, -0.2e1 * t26 * t123 + 0.2e1 * t89 * t71, 0.2e1 * t27 * t123 + 0.2e1 * t89 * t73, -0.2e1 * t26 * t73 - 0.2e1 * t27 * t71, t26 ^ 2 + t27 ^ 2 + t89 ^ 2, t70, t42, t60, t69, -t59, t114, -0.2e1 * t16 * t123 + 0.2e1 * t46 * t71, 0.2e1 * t17 * t123 + 0.2e1 * t46 * t73, -0.2e1 * t16 * t73 - 0.2e1 * t17 * t71, t16 ^ 2 + t17 ^ 2 + t46 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, -t77, -t148, t38, -t39, 0, 0, t155, t118 * t51 + t52 * t122, t160, t156, t159, 0, pkin(3) * t51 - pkin(10) * t160 - t157, -pkin(3) * t52 - pkin(10) * t159 + t158 (t155 + t156) * pkin(10) + t133, -t36 * pkin(3) + pkin(10) * t133, t24, t8, t50, t23, -t49, 0, t104 * t30 + t18 * t85 + t55 * t77, t104 * t32 + t18 * t87 - t56 * t77, -t3 * t87 - t56 * t30 - t55 * t32 - t4 * t85, t18 * t104 + t3 * t55 + t4 * t56, t24, t8, t50, t23, -t49, 0, t64 * t30 + t40 * t77 + t9 * t85, t64 * t32 - t41 * t77 + t9 * t87, -t1 * t87 - t2 * t85 - t41 * t30 - t40 * t32, t1 * t40 + t2 * t41 + t9 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, 0, t123, 0, -t108, -t123 * pkin(9), 0, 0, t134 (-t111 + t113) * t119, -t144, -t134, -t142, 0, -pkin(9) * t143 + (-pkin(3) * t119 + pkin(10) * t123) * t118, pkin(10) * t142 + (t164 - t165) * t119, t132, -pkin(3) * t108 + pkin(10) * t132, t48, t33, -t151, t47, t152, 0, t104 * t71 - t55 * t123 + t89 * t85, t104 * t73 + t56 * t123 + t89 * t87, -t26 * t87 - t27 * t85 - t55 * t73 - t56 * t71, t89 * t104 + t26 * t55 + t27 * t56, t48, t33, -t151, t47, t152, 0, -t40 * t123 + t46 * t85 + t64 * t71, t41 * t123 + t46 * t87 + t64 * t73, -t16 * t87 - t17 * t85 - t40 * t73 - t41 * t71, t16 * t40 + t17 * t41 + t46 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t111, 0.2e1 * t145, 0, t113, 0, 0, 0.2e1 * t165, -0.2e1 * pkin(3) * t118, 0.2e1 * t141 * pkin(10), pkin(10) ^ 2 * t141 + pkin(3) ^ 2, t84, t54, 0, t83, 0, 0, t85 * t173, t87 * t173, -0.2e1 * t55 * t87 - 0.2e1 * t56 * t85, t104 ^ 2 + t55 ^ 2 + t56 ^ 2, t84, t54, 0, t83, 0, 0, t85 * t176, t87 * t176, -0.2e1 * t40 * t87 - 0.2e1 * t41 * t85, t40 ^ 2 + t41 ^ 2 + t64 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, t51, t77, t12, -t13, 0, 0, 0, 0, t32, 0, -t30, t77, t77 * t109 + t3, -t4 - t63, -t32 * t109 - t25 (t117 * t4 + t121 * t3) * pkin(4), 0, 0, t32, 0, -t30, t77, t103 * t77 + t1, -t2 - t63, -t103 * t32 - t25, t1 * t103 + t2 * t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, 0, -t146, -t123, t61, -t62, 0, 0, 0, 0, t73, 0, -t71, -t123, -t123 * t109 + t26, -t27 + t98, -t73 * t109 - t58 (t117 * t27 + t121 * t26) * pkin(4), 0, 0, t73, 0, -t71, -t123 (-pkin(5) - t103) * t123 + t129, -t17 + t98, -t103 * t73 - t58, t16 * t103 + t17 * t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, 0, t122, 0, -t118 * pkin(10), -t122 * pkin(10), 0, 0, 0, 0, t87, 0, -t85, 0, t55, -t56, -t87 * t109 - t75 (t117 * t56 + t121 * t55) * pkin(4), 0, 0, t87, 0, -t85, 0, t40, -t41, -t103 * t87 - t75, t40 * t103 + t41 * t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t109, t105, 0, t121 ^ 2 * t128 + t106, 0, 0, 0, 0, 0, 1, 0.2e1 * t103, t105, 0, t103 ^ 2 + t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, -t30, t77, t3, -t4, 0, 0, 0, 0, t32, 0, -t30, t77, t130 + 0.2e1 * t74, -t2, -t32 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, -t71, -t123, t26, -t27, 0, 0, 0, 0, t73, 0, -t71, -t123, t129 - 0.2e1 * t161, -t17, -t73 * pkin(5), t16 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, 0, -t85, 0, t55, -t56, 0, 0, 0, 0, t87, 0, -t85, 0, t40, -t41, -t87 * pkin(5), t40 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t109, -t162, 0, 0, 0, 0, 0, 0, 0, 1, t125 + t109, -t162, 0, t103 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t125, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t32, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, t73, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t87, 0, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
