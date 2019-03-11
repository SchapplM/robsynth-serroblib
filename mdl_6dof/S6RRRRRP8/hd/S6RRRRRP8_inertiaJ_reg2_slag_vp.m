% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRRRRP8
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
% Datum: 2019-03-10 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRRP8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_inertiaJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t120 = sin(qJ(5));
t114 = t120 ^ 2;
t124 = cos(qJ(5));
t116 = t124 ^ 2;
t209 = t114 + t116;
t121 = sin(qJ(4));
t192 = t121 * pkin(3);
t106 = pkin(11) + t192;
t180 = t209 * t106;
t118 = sin(pkin(6));
t127 = cos(qJ(2));
t165 = t118 * t127;
t125 = cos(qJ(4));
t119 = cos(pkin(6));
t122 = sin(qJ(3));
t126 = cos(qJ(3));
t123 = sin(qJ(2));
t166 = t118 * t123;
t70 = t119 * t126 - t122 * t166;
t71 = t119 * t122 + t126 * t166;
t50 = t121 * t70 + t125 * t71;
t40 = t120 * t50 + t124 * t165;
t36 = t40 * t124;
t42 = -t120 * t165 + t124 * t50;
t37 = t42 * t120;
t85 = t121 * t126 + t125 * t122;
t208 = (t37 + t36) * t85;
t207 = t40 ^ 2;
t48 = t121 * t71 - t125 * t70;
t47 = t48 ^ 2;
t199 = -pkin(10) - pkin(9);
t146 = t199 * t122;
t92 = t199 * t126;
t60 = -t121 * t92 - t125 * t146;
t206 = t60 ^ 2;
t83 = t121 * t122 - t125 * t126;
t81 = t83 ^ 2;
t205 = 0.2e1 * t85;
t108 = -t126 * pkin(3) - pkin(2);
t204 = 0.2e1 * t108;
t203 = 0.2e1 * t118;
t202 = -0.2e1 * t120;
t201 = -0.2e1 * t124;
t200 = 0.2e1 * t126;
t172 = t48 * qJ(6);
t151 = pkin(8) * t165;
t195 = pkin(1) * t123;
t66 = t151 + (pkin(9) + t195) * t119;
t67 = (-pkin(2) * t127 - pkin(9) * t123 - pkin(1)) * t118;
t44 = t122 * t67 + t126 * t66;
t32 = t70 * pkin(10) + t44;
t176 = t125 * t32;
t153 = pkin(3) * t165;
t43 = -t122 * t66 + t126 * t67;
t27 = -t71 * pkin(10) - t153 + t43;
t15 = t121 * t27 + t176;
t13 = -pkin(11) * t165 + t15;
t194 = pkin(1) * t127;
t99 = pkin(8) * t166;
t65 = t99 + (-pkin(2) - t194) * t119;
t53 = -t70 * pkin(3) + t65;
t20 = t48 * pkin(4) - t50 * pkin(11) + t53;
t7 = t120 * t20 + t124 * t13;
t3 = t172 + t7;
t143 = t120 * t13 - t124 * t20;
t197 = t48 * pkin(5);
t4 = t143 - t197;
t139 = t4 * t120 + t3 * t124;
t198 = pkin(11) * t83;
t196 = t83 * pkin(5);
t193 = t120 * pkin(11);
t191 = t124 * pkin(11);
t190 = t125 * pkin(3);
t189 = t42 * t40;
t188 = t48 * t40;
t39 = t48 * t83;
t5 = t7 * t124;
t184 = t121 * t32 - t125 * t27;
t12 = pkin(4) * t165 + t184;
t8 = t40 * pkin(5) - t42 * qJ(6) + t12;
t187 = t8 * t120;
t186 = t8 * t124;
t107 = -pkin(4) - t190;
t185 = pkin(4) - t107;
t54 = t83 * pkin(4) - t85 * pkin(11) + t108;
t62 = t121 * t146 - t125 * t92;
t31 = t120 * t54 + t124 * t62;
t136 = -t124 * pkin(5) - t120 * qJ(6);
t88 = -pkin(4) + t136;
t80 = t88 - t190;
t183 = -t80 - t88;
t182 = t180 * pkin(11);
t181 = t209 * t106 ^ 2;
t179 = t106 * t83;
t178 = t12 * t124;
t45 = t120 * t48;
t177 = t120 * t85;
t46 = t124 * t48;
t77 = t124 * t85;
t135 = pkin(5) * t120 - t124 * qJ(6);
t34 = t135 * t85 + t60;
t175 = t34 * t120;
t174 = t34 * t124;
t173 = t42 * t124;
t58 = t60 * t120;
t171 = t60 * t124;
t170 = t70 * t126;
t169 = t71 * t122;
t168 = t83 * qJ(6);
t113 = t118 ^ 2;
t167 = t113 * t127;
t164 = t119 * t123;
t163 = t120 * t106;
t162 = t120 * t124;
t161 = t124 * t106;
t160 = t209 * pkin(11) ^ 2;
t159 = t209 * pkin(11);
t115 = t122 ^ 2;
t117 = t126 ^ 2;
t158 = t115 + t117;
t157 = pkin(11) * t45;
t156 = pkin(11) * t46;
t155 = -0.2e1 * t165;
t154 = 0.2e1 * t165;
t152 = t83 * t177;
t150 = t40 * t177;
t149 = t48 * t163;
t148 = t48 * t161;
t82 = t85 ^ 2;
t147 = t82 * t162;
t145 = t122 * t165;
t144 = t126 * t165;
t142 = t120 * t62 - t124 * t54;
t141 = -pkin(4) * t85 - t198;
t140 = -t85 * t88 + t198;
t138 = t120 * t143 + t5;
t137 = -t80 * t85 + t179;
t134 = t107 * t85 - t179;
t24 = t168 + t31;
t25 = t142 - t196;
t9 = t25 * t120 + t24 * t124;
t16 = t120 * t142 + t31 * t124;
t21 = t120 * t40 - t173;
t132 = -t43 * t122 + t44 * t126;
t131 = t48 * t177 + t83 * t40;
t104 = t113 * t127 ^ 2;
t102 = -0.2e1 * t162;
t101 = 0.2e1 * t162;
t87 = 0.2e1 * t159;
t79 = pkin(1) * t164 + t151;
t78 = t119 * t194 - t99;
t76 = t124 * t83;
t74 = t116 * t82;
t73 = t120 * t83;
t72 = t114 * t82;
t69 = t85 * t162;
t64 = 0.2e1 * t180;
t59 = t159 + t180;
t56 = 0.2e1 * t83 * t77;
t55 = (t114 - t116) * t85;
t38 = t42 ^ 2;
t35 = pkin(11) * t36;
t33 = t40 * t161;
t28 = t85 * t173;
t23 = 0.2e1 * t42 * t48;
t19 = t42 * t83 + t48 * t77;
t11 = t12 * t120;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t113 * t123 ^ 2, 0.2e1 * t123 * t167, t164 * t203, t104, t119 * t154, t119 ^ 2, 0.2e1 * pkin(1) * t167 + 0.2e1 * t78 * t119, -0.2e1 * t113 * t195 - 0.2e1 * t79 * t119 (-t123 * t78 + t127 * t79) * t203, t113 * pkin(1) ^ 2 + t78 ^ 2 + t79 ^ 2, t71 ^ 2, 0.2e1 * t71 * t70, t71 * t155, t70 ^ 2, t70 * t155, t104, -0.2e1 * t43 * t165 - 0.2e1 * t65 * t70, 0.2e1 * t44 * t165 + 0.2e1 * t65 * t71, -0.2e1 * t43 * t71 + 0.2e1 * t44 * t70, t43 ^ 2 + t44 ^ 2 + t65 ^ 2, t50 ^ 2, -0.2e1 * t50 * t48, t50 * t155, t47, t48 * t154, t104, 0.2e1 * t165 * t184 + 0.2e1 * t53 * t48, 0.2e1 * t15 * t165 + 0.2e1 * t53 * t50, -0.2e1 * t15 * t48 + 0.2e1 * t184 * t50, t15 ^ 2 + t184 ^ 2 + t53 ^ 2, t38, -0.2e1 * t189, t23, t207, -0.2e1 * t188, t47, 0.2e1 * t12 * t40 - 0.2e1 * t143 * t48, 0.2e1 * t12 * t42 - 0.2e1 * t7 * t48, 0.2e1 * t143 * t42 - 0.2e1 * t7 * t40, t12 ^ 2 + t143 ^ 2 + t7 ^ 2, t38, t23, 0.2e1 * t189, t47, 0.2e1 * t188, t207, -0.2e1 * t4 * t48 + 0.2e1 * t8 * t40, -0.2e1 * t3 * t40 + 0.2e1 * t4 * t42, 0.2e1 * t3 * t48 - 0.2e1 * t8 * t42, t3 ^ 2 + t4 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, 0, t165, t119, t78, -t79, 0, 0, t169, t122 * t70 + t71 * t126, -t145, t170, -t144, 0, pkin(2) * t70 + pkin(9) * t145 - t65 * t126, -pkin(2) * t71 + pkin(9) * t144 + t65 * t122 (t169 + t170) * pkin(9) + t132, -t65 * pkin(2) + t132 * pkin(9), t50 * t85, -t85 * t48 - t50 * t83, -t85 * t165, t39, t83 * t165, 0, t108 * t48 + t60 * t165 + t53 * t83, t108 * t50 + t165 * t62 + t53 * t85, -t15 * t83 + t184 * t85 - t62 * t48 + t60 * t50, t53 * t108 + t15 * t62 + t184 * t60, t28, -t208, t19, t150, -t131, t39, t12 * t177 - t142 * t48 - t143 * t83 + t60 * t40, t12 * t77 - t31 * t48 + t60 * t42 - t7 * t83, t142 * t42 - t31 * t40 + (-t120 * t7 + t124 * t143) * t85, t12 * t60 + t142 * t143 + t7 * t31, t28, t19, t208, t39, t131, t150, t8 * t177 - t25 * t48 + t34 * t40 - t4 * t83, -t24 * t40 + t25 * t42 + (-t120 * t3 + t124 * t4) * t85, t24 * t48 + t3 * t83 - t34 * t42 - t77 * t8, t3 * t24 + t4 * t25 + t8 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t115, t122 * t200, 0, t117, 0, 0, pkin(2) * t200, -0.2e1 * pkin(2) * t122, 0.2e1 * t158 * pkin(9), t158 * pkin(9) ^ 2 + pkin(2) ^ 2, t82, -0.2e1 * t85 * t83, 0, t81, 0, 0, t83 * t204, t85 * t204, 0.2e1 * t60 * t85 - 0.2e1 * t62 * t83, t108 ^ 2 + t62 ^ 2 + t206, t74, -0.2e1 * t147, t56, t72, -0.2e1 * t152, t81, -0.2e1 * t142 * t83 + 0.2e1 * t85 * t58, 0.2e1 * t171 * t85 - 0.2e1 * t31 * t83 (-t120 * t31 + t124 * t142) * t205, t142 ^ 2 + t31 ^ 2 + t206, t74, t56, 0.2e1 * t147, t81, 0.2e1 * t152, t72, 0.2e1 * t175 * t85 - 0.2e1 * t25 * t83 (-t120 * t24 + t124 * t25) * t205, -0.2e1 * t174 * t85 + 0.2e1 * t24 * t83, t24 ^ 2 + t25 ^ 2 + t34 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, t70, -t165, t43, -t44, 0, 0, 0, 0, t50, 0, -t48, -t165, -t125 * t153 - t184, -t176 + (-t27 + t153) * t121 (-t121 * t48 - t125 * t50) * pkin(3) (t121 * t15 - t125 * t184) * pkin(3), t37, -t21, t45, -t36, t46, 0, t107 * t40 - t149 - t178, t107 * t42 + t11 - t148, -t33 + t5 + (t106 * t42 + t143) * t120, t106 * t138 + t12 * t107, t37, t45, t21, 0, -t46, -t36, t80 * t40 - t149 - t186, t163 * t42 + t139 - t33, -t80 * t42 + t148 - t187, t106 * t139 + t8 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t122, 0, t126, 0, -t122 * pkin(9), -t126 * pkin(9), 0, 0, 0, 0, t85, 0, -t83, 0, -t60, -t62 (-t121 * t83 - t125 * t85) * pkin(3) (t121 * t62 - t125 * t60) * pkin(3), t69, -t55, t73, -t69, t76, 0, t120 * t134 - t171, t124 * t134 + t58, t16, t106 * t16 + t60 * t107, t69, t73, t55, 0, -t76, -t69, -t120 * t137 - t174, t9, t124 * t137 - t175, t106 * t9 + t34 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t190, -0.2e1 * t192, 0 (t121 ^ 2 + t125 ^ 2) * pkin(3) ^ 2, t114, t101, 0, t116, 0, 0, t107 * t201, 0.2e1 * t107 * t120, t64, t107 ^ 2 + t181, t114, 0, t102, 0, 0, t116, t80 * t201, t64, t80 * t202, t80 ^ 2 + t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, 0, -t48, -t165, -t184, -t15, 0, 0, t37, -t21, t45, -t36, t46, 0, -pkin(4) * t40 - t157 - t178, -pkin(4) * t42 + t11 - t156, -t35 + t5 + (pkin(11) * t42 + t143) * t120, -t12 * pkin(4) + pkin(11) * t138, t37, t45, t21, 0, -t46, -t36, t88 * t40 - t157 - t186, pkin(11) * t37 + t139 - t35, -t88 * t42 + t156 - t187, pkin(11) * t139 + t8 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, 0, -t83, 0, -t60, -t62, 0, 0, t69, -t55, t73, -t69, t76, 0, t120 * t141 - t171, t124 * t141 + t58, t16, -t60 * pkin(4) + pkin(11) * t16, t69, t73, t55, 0, -t76, -t69, -t120 * t140 - t174, t9, t124 * t140 - t175, pkin(11) * t9 + t34 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t190, -t192, 0, 0, t114, t101, 0, t116, 0, 0, t185 * t124, -t185 * t120, t59, -t107 * pkin(4) + t182, t114, 0, t102, 0, 0, t116, t183 * t124, t59, t183 * t120, t80 * t88 + t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t114, t101, 0, t116, 0, 0, 0.2e1 * pkin(4) * t124, pkin(4) * t202, t87, pkin(4) ^ 2 + t160, t114, 0, t102, 0, 0, t116, t88 * t201, t87, t88 * t202, t88 ^ 2 + t160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, -t40, t48, -t143, -t7, 0, 0, 0, t42, 0, t48, t40, 0, -t143 + 0.2e1 * t197, -pkin(5) * t42 - t40 * qJ(6), 0.2e1 * t172 + t7, -t4 * pkin(5) + t3 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, 0, -t177, t83, -t142, -t31, 0, 0, 0, t77, 0, t83, t177, 0, -t142 + 0.2e1 * t196, t136 * t85, 0.2e1 * t168 + t31, -t25 * pkin(5) + t24 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, 0, t124, 0, -t163, -t161, 0, 0, 0, t120, 0, 0, -t124, 0, -t163, -t135, t161, -t135 * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, 0, t124, 0, -t193, -t191, 0, 0, 0, t120, 0, 0, -t124, 0, -t193, -t135, t191, -t135 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t42, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, t77, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, 0, t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, 0, t193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t1;
