% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRR13_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_inertiaJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t110 = sin(pkin(13));
t113 = cos(pkin(13));
t117 = sin(qJ(5));
t185 = cos(qJ(5));
t92 = t185 * t110 + t117 * t113;
t199 = -0.2e1 * t92;
t118 = sin(qJ(3));
t121 = cos(qJ(3));
t111 = sin(pkin(7));
t114 = cos(pkin(7));
t115 = cos(pkin(6));
t112 = sin(pkin(6));
t119 = sin(qJ(2));
t152 = t112 * t119;
t122 = cos(qJ(2));
t99 = t115 * t122 * pkin(1);
t65 = t115 * pkin(2) + t99 + (-pkin(10) * t114 - pkin(9)) * t152;
t73 = (-pkin(10) * t111 * t119 - pkin(2) * t122 - pkin(1)) * t112;
t127 = t111 * t73 + t114 * t65;
t149 = t114 * t122;
t137 = t112 * t149;
t151 = t112 * t122;
t184 = pkin(1) * t119;
t87 = pkin(9) * t151 + t115 * t184;
t60 = (t111 * t115 + t137) * pkin(10) + t87;
t34 = -t118 * t60 + t127 * t121;
t154 = t111 * t118;
t64 = t115 * t154 + (t118 * t149 + t119 * t121) * t112;
t80 = t111 * t151 - t115 * t114;
t43 = t64 * t110 + t80 * t113;
t45 = -t80 * t110 + t64 * t113;
t31 = t117 * t45 + t185 * t43;
t198 = t31 ^ 2;
t78 = t110 * t154 - t113 * t114;
t82 = t110 * t114 + t113 * t154;
t52 = t117 * t82 + t185 * t78;
t197 = t52 ^ 2;
t153 = t111 * t121;
t62 = -t115 * t153 + t118 * t152 - t121 * t137;
t61 = t62 ^ 2;
t178 = pkin(11) + qJ(4);
t135 = t178 * t110;
t94 = t178 * t113;
t66 = t117 * t94 + t185 * t135;
t196 = t66 ^ 2;
t90 = t117 * t110 - t185 * t113;
t195 = t90 ^ 2;
t194 = -0.2e1 * t31;
t193 = -0.2e1 * t52;
t192 = -0.2e1 * t62;
t191 = 0.2e1 * t62;
t190 = -0.2e1 * t64;
t103 = -t113 * pkin(4) - pkin(3);
t189 = 0.2e1 * t103;
t188 = 0.2e1 * t111;
t187 = 0.2e1 * t112;
t186 = 0.2e1 * t113;
t183 = pkin(2) * t118;
t182 = pkin(2) * t121;
t181 = t31 * t52;
t180 = t31 * t90;
t179 = t52 * t90;
t40 = -t111 * t65 + t114 * t73;
t25 = t62 * pkin(3) - t64 * qJ(4) + t40;
t35 = t127 * t118 + t121 * t60;
t29 = -t80 * qJ(4) + t35;
t13 = t110 * t25 + t113 * t29;
t141 = pkin(10) * t153;
t74 = t141 + (qJ(4) + t183) * t114;
t75 = (-pkin(3) * t121 - qJ(4) * t118 - pkin(2)) * t111;
t50 = t110 * t75 + t113 * t74;
t177 = t110 * t62;
t176 = t113 * t62;
t116 = sin(qJ(6));
t175 = t116 * t31;
t174 = t116 * t52;
t173 = t116 * t90;
t172 = t116 * t92;
t120 = cos(qJ(6));
t28 = t120 * t31;
t51 = t120 * t52;
t171 = t120 * t92;
t33 = -t117 * t43 + t185 * t45;
t16 = t116 * t33 - t62 * t120;
t170 = t16 * t120;
t18 = t62 * t116 + t120 * t33;
t169 = t18 * t116;
t168 = t18 * t120;
t49 = -t110 * t74 + t113 * t75;
t36 = -pkin(4) * t153 - t82 * pkin(11) + t49;
t39 = -t78 * pkin(11) + t50;
t21 = -t117 * t39 + t185 * t36;
t19 = pkin(5) * t153 - t21;
t167 = t19 * t116;
t166 = t19 * t120;
t165 = t43 * t113;
t164 = t45 * t110;
t54 = -t117 * t78 + t185 * t82;
t46 = t116 * t54 + t120 * t153;
t163 = t46 * t120;
t48 = -t116 * t153 + t120 * t54;
t162 = t48 * t116;
t161 = t48 * t120;
t160 = t66 * t116;
t159 = t66 * t120;
t158 = t78 * t113;
t157 = t82 * t110;
t105 = t111 ^ 2;
t156 = t105 * t121;
t106 = t112 ^ 2;
t155 = t106 * t122;
t150 = t114 * t118;
t148 = t116 * t120;
t104 = t110 ^ 2;
t107 = t113 ^ 2;
t147 = t104 + t107;
t108 = t116 ^ 2;
t109 = t120 ^ 2;
t146 = t108 + t109;
t145 = t90 * t199;
t144 = -0.2e1 * t153;
t143 = 0.2e1 * t153;
t142 = t115 * t187;
t140 = t92 * t148;
t139 = t62 * t153;
t138 = t110 * t153;
t136 = t113 * t153;
t12 = -t110 * t29 + t113 * t25;
t134 = -pkin(5) * t92 - pkin(12) * t90;
t8 = t62 * pkin(4) - t45 * pkin(11) + t12;
t9 = -t43 * pkin(11) + t13;
t6 = t117 * t8 + t185 * t9;
t4 = t62 * pkin(12) + t6;
t30 = t80 * pkin(3) - t34;
t14 = t43 * pkin(4) + t30;
t7 = t31 * pkin(5) - t33 * pkin(12) + t14;
t1 = -t116 * t4 + t120 * t7;
t2 = t116 * t7 + t120 * t4;
t133 = t1 * t120 + t2 * t116;
t132 = -t1 * t116 + t2 * t120;
t22 = t117 * t36 + t185 * t39;
t20 = -pkin(12) * t153 + t22;
t98 = pkin(10) * t154;
t77 = t98 + (-pkin(3) - t182) * t114;
t55 = t78 * pkin(4) + t77;
t27 = t52 * pkin(5) - t54 * pkin(12) + t55;
t10 = -t116 * t20 + t120 * t27;
t11 = t116 * t27 + t120 * t20;
t131 = t10 * t120 + t11 * t116;
t130 = -t10 * t116 + t11 * t120;
t129 = -t12 * t110 + t13 * t113;
t128 = -t49 * t110 + t50 * t113;
t57 = t90 * pkin(5) - t92 * pkin(12) + t103;
t68 = -t117 * t135 + t185 * t94;
t37 = -t116 * t68 + t120 * t57;
t38 = t116 * t57 + t120 * t68;
t126 = t38 * t116 + t37 * t120;
t125 = -t37 * t116 + t38 * t120;
t5 = -t117 * t9 + t185 * t8;
t101 = t105 * t121 ^ 2;
t88 = t92 ^ 2;
t86 = pkin(2) * t150 + t141;
t85 = -pkin(9) * t152 + t99;
t84 = t114 * t182 - t98;
t83 = t120 * t90;
t42 = t116 * t46;
t15 = t116 * t16;
t3 = -t62 * pkin(5) - t5;
t17 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t106 * t119 ^ 2, 0.2e1 * t119 * t155, t119 * t142, t106 * t122 ^ 2, t122 * t142, t115 ^ 2, 0.2e1 * pkin(1) * t155 + 0.2e1 * t85 * t115, -0.2e1 * t106 * t184 - 0.2e1 * t87 * t115 (-t119 * t85 + t122 * t87) * t187, t106 * pkin(1) ^ 2 + t85 ^ 2 + t87 ^ 2, t64 ^ 2, t62 * t190, t80 * t190, t61, t80 * t191, t80 ^ 2, -0.2e1 * t34 * t80 + 0.2e1 * t40 * t62, 0.2e1 * t35 * t80 + 0.2e1 * t40 * t64, -0.2e1 * t34 * t64 - 0.2e1 * t35 * t62, t34 ^ 2 + t35 ^ 2 + t40 ^ 2, t45 ^ 2, -0.2e1 * t45 * t43, t45 * t191, t43 ^ 2, t43 * t192, t61, 0.2e1 * t12 * t62 + 0.2e1 * t30 * t43, -0.2e1 * t13 * t62 + 0.2e1 * t30 * t45, -0.2e1 * t12 * t45 - 0.2e1 * t13 * t43, t12 ^ 2 + t13 ^ 2 + t30 ^ 2, t33 ^ 2, t33 * t194, t33 * t191, t198, t31 * t192, t61, 0.2e1 * t14 * t31 + 0.2e1 * t5 * t62, 0.2e1 * t14 * t33 - 0.2e1 * t6 * t62, -0.2e1 * t6 * t31 - 0.2e1 * t5 * t33, t14 ^ 2 + t5 ^ 2 + t6 ^ 2, t18 ^ 2, -0.2e1 * t18 * t16, 0.2e1 * t18 * t31, t16 ^ 2, t16 * t194, t198, 0.2e1 * t1 * t31 + 0.2e1 * t3 * t16, 0.2e1 * t3 * t18 - 0.2e1 * t2 * t31, -0.2e1 * t1 * t18 - 0.2e1 * t2 * t16, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, 0, t151, t115, t85, -t87, 0, 0, t64 * t154 (-t118 * t62 + t121 * t64) * t111, t64 * t114 - t80 * t154, -t139, -t62 * t114 - t153 * t80, -t80 * t114, t34 * t114 - t84 * t80 + (-pkin(2) * t62 - t121 * t40) * t111, -t35 * t114 + t86 * t80 + (-pkin(2) * t64 + t118 * t40) * t111, -t86 * t62 - t84 * t64 + (-t118 * t34 + t121 * t35) * t111, -t40 * t111 * pkin(2) + t34 * t84 + t35 * t86, t45 * t82, -t82 * t43 - t45 * t78, -t153 * t45 + t82 * t62, t43 * t78, t153 * t43 - t78 * t62, -t139, -t12 * t153 + t30 * t78 + t77 * t43 + t49 * t62, t13 * t153 + t30 * t82 + t77 * t45 - t50 * t62, -t12 * t82 - t13 * t78 - t50 * t43 - t49 * t45, t12 * t49 + t13 * t50 + t30 * t77, t33 * t54, -t54 * t31 - t33 * t52, -t153 * t33 + t54 * t62, t181, t153 * t31 - t52 * t62, -t139, t14 * t52 - t153 * t5 + t21 * t62 + t55 * t31, t14 * t54 + t153 * t6 - t22 * t62 + t55 * t33, -t21 * t33 - t22 * t31 - t5 * t54 - t6 * t52, t14 * t55 + t5 * t21 + t6 * t22, t18 * t48, -t48 * t16 - t18 * t46, t18 * t52 + t48 * t31, t16 * t46, -t16 * t52 - t46 * t31, t181, t1 * t52 + t10 * t31 + t19 * t16 + t3 * t46, -t11 * t31 + t19 * t18 - t2 * t52 + t3 * t48, -t1 * t48 - t10 * t18 - t11 * t16 - t2 * t46, t1 * t10 + t2 * t11 + t3 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t105 * t118 ^ 2, 0.2e1 * t118 * t156, t150 * t188, t101, t114 * t143, t114 ^ 2, 0.2e1 * pkin(2) * t156 + 0.2e1 * t84 * t114, -0.2e1 * t105 * t183 - 0.2e1 * t86 * t114 (-t118 * t84 + t121 * t86) * t188, t105 * pkin(2) ^ 2 + t84 ^ 2 + t86 ^ 2, t82 ^ 2, -0.2e1 * t82 * t78, t82 * t144, t78 ^ 2, t78 * t143, t101, -0.2e1 * t153 * t49 + 0.2e1 * t77 * t78, 0.2e1 * t153 * t50 + 0.2e1 * t77 * t82, -0.2e1 * t49 * t82 - 0.2e1 * t50 * t78, t49 ^ 2 + t50 ^ 2 + t77 ^ 2, t54 ^ 2, t54 * t193, t54 * t144, t197, t52 * t143, t101, -0.2e1 * t153 * t21 + 0.2e1 * t55 * t52, 0.2e1 * t153 * t22 + 0.2e1 * t55 * t54, -0.2e1 * t21 * t54 - 0.2e1 * t22 * t52, t21 ^ 2 + t22 ^ 2 + t55 ^ 2, t48 ^ 2, -0.2e1 * t48 * t46, 0.2e1 * t48 * t52, t46 ^ 2, t46 * t193, t197, 0.2e1 * t10 * t52 + 0.2e1 * t19 * t46, -0.2e1 * t11 * t52 + 0.2e1 * t19 * t48, -0.2e1 * t10 * t48 - 0.2e1 * t11 * t46, t10 ^ 2 + t11 ^ 2 + t19 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, -t62, -t80, t34, -t35, 0, 0, t164, -t110 * t43 + t45 * t113, t177, -t165, t176, 0, -pkin(3) * t43 - qJ(4) * t177 - t30 * t113, -pkin(3) * t45 - qJ(4) * t176 + t30 * t110 (t164 - t165) * qJ(4) + t129, -t30 * pkin(3) + qJ(4) * t129, t33 * t92, -t92 * t31 - t33 * t90, t92 * t62, t180, -t90 * t62, 0, t103 * t31 + t14 * t90 - t66 * t62, t103 * t33 + t14 * t92 - t68 * t62, -t68 * t31 + t66 * t33 - t5 * t92 - t6 * t90, t14 * t103 - t5 * t66 + t6 * t68, t92 * t168 (-t169 - t170) * t92, t171 * t31 + t18 * t90, t16 * t172, -t16 * t90 - t172 * t31, t180, t1 * t90 + t66 * t16 + t172 * t3 + t37 * t31, t171 * t3 + t66 * t18 - t2 * t90 - t38 * t31, -t133 * t92 - t38 * t16 - t37 * t18, t1 * t37 + t2 * t38 + t3 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t154, 0, t153, t114, t84, -t86, 0, 0, t157, -t110 * t78 + t82 * t113, -t138, -t158, -t136, 0, -pkin(3) * t78 + qJ(4) * t138 - t77 * t113, -pkin(3) * t82 + qJ(4) * t136 + t77 * t110 (t157 - t158) * qJ(4) + t128, -t77 * pkin(3) + qJ(4) * t128, t54 * t92, -t92 * t52 - t54 * t90, -t92 * t153, t179, t90 * t153, 0, t103 * t52 + t153 * t66 + t55 * t90, t103 * t54 + t153 * t68 + t55 * t92, -t21 * t92 - t22 * t90 - t68 * t52 + t66 * t54, t55 * t103 - t21 * t66 + t22 * t68, t92 * t161 (-t162 - t163) * t92, t171 * t52 + t48 * t90, t46 * t172, -t172 * t52 - t46 * t90, t179, t10 * t90 + t167 * t92 + t37 * t52 + t66 * t46, -t11 * t90 + t166 * t92 - t38 * t52 + t66 * t48, -t131 * t92 - t37 * t48 - t38 * t46, t10 * t37 + t11 * t38 + t19 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t104, t110 * t186, 0, t107, 0, 0, pkin(3) * t186, -0.2e1 * pkin(3) * t110, 0.2e1 * t147 * qJ(4), qJ(4) ^ 2 * t147 + pkin(3) ^ 2, t88, t145, 0, t195, 0, 0, t90 * t189, t92 * t189, 0.2e1 * t66 * t92 - 0.2e1 * t68 * t90, t103 ^ 2 + t68 ^ 2 + t196, t109 * t88, -0.2e1 * t88 * t148, 0.2e1 * t90 * t171, t108 * t88, t116 * t145, t195, 0.2e1 * t160 * t92 + 0.2e1 * t37 * t90, 0.2e1 * t159 * t92 - 0.2e1 * t38 * t90, t126 * t199, t37 ^ 2 + t38 ^ 2 + t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t45, 0, t30, 0, 0, 0, 0, 0, 0, t31, t33, 0, t14, 0, 0, 0, 0, 0, 0, t28, -t175, -t15 - t168, t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t82, 0, t77, 0, 0, 0, 0, 0, 0, t52, t54, 0, t55, 0, 0, 0, 0, 0, 0, t51, -t174, -t42 - t161, t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, t110, 0, -pkin(3), 0, 0, 0, 0, 0, 0, t90, t92, 0, t103, 0, 0, 0, 0, 0, 0, t83, -t173, -t146 * t92, t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, -t31, t62, t5, -t6, 0, 0, t169, -t15 + t168, t175, -t170, t28, 0, -pkin(5) * t16 - pkin(12) * t175 - t3 * t120, -pkin(5) * t18 - pkin(12) * t28 + t3 * t116 (t169 - t170) * pkin(12) + t132, -t3 * pkin(5) + pkin(12) * t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, -t52, -t153, t21, -t22, 0, 0, t162, -t42 + t161, t174, -t163, t51, 0, -pkin(5) * t46 - pkin(12) * t174 - t166, -pkin(5) * t48 - pkin(12) * t51 + t167 (t162 - t163) * pkin(12) + t130, -t19 * pkin(5) + pkin(12) * t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, 0, -t90, 0, -t66, -t68, 0, 0, t140 (-t108 + t109) * t92, t173, -t140, t83, 0, t116 * t134 - t159, t120 * t134 + t160, t125, -t66 * pkin(5) + pkin(12) * t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t108, 0.2e1 * t148, 0, t109, 0, 0, 0.2e1 * pkin(5) * t120, -0.2e1 * pkin(5) * t116, 0.2e1 * t146 * pkin(12), pkin(12) ^ 2 * t146 + pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t16, t31, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, -t46, t52, t10, -t11, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, 0, -t172, t90, t37, -t38, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, -t116, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, 0, t120, 0, -t116 * pkin(12), -t120 * pkin(12), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t17;