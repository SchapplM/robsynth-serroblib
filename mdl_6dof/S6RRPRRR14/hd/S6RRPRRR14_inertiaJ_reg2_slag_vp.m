% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2019-01-03 10:26
% Revision: 5fdbc45bcf2cc60deefd7ac2d71d743ed41bf7e4 (2018-12-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRR14_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_inertiaJ_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t120 = sin(qJ(4));
t124 = cos(qJ(4));
t111 = sin(pkin(8));
t115 = cos(pkin(8));
t116 = cos(pkin(7));
t110 = sin(pkin(14));
t112 = sin(pkin(7));
t155 = t112 * t110;
t188 = pkin(11) * t115;
t114 = cos(pkin(14));
t150 = t114 * t116;
t98 = pkin(2) * t150;
t67 = t116 * pkin(3) + t98 + (-qJ(3) - t188) * t155;
t189 = pkin(11) * t111;
t75 = (-pkin(3) * t114 - t110 * t189 - pkin(2)) * t112;
t133 = t111 * t75 + t115 * t67;
t151 = t114 * t115;
t140 = t112 * t151;
t154 = t112 * t114;
t158 = t110 * t116;
t84 = pkin(2) * t158 + qJ(3) * t154;
t60 = (t111 * t116 + t140) * pkin(11) + t84;
t35 = -t120 * t60 + t124 * t133;
t117 = cos(pkin(6));
t113 = sin(pkin(6));
t125 = cos(qJ(2));
t152 = t113 * t125;
t128 = t112 * t117 + t116 * t152;
t121 = sin(qJ(2));
t191 = pkin(1) * t121;
t89 = pkin(10) * t152 + t117 * t191;
t61 = qJ(3) * t128 + t89;
t153 = t113 * t121;
t99 = t117 * t125 * pkin(1);
t68 = t117 * pkin(2) + t99 + (-qJ(3) * t116 - pkin(10)) * t153;
t76 = (-qJ(3) * t112 * t121 - pkin(2) * t125 - pkin(1)) * t113;
t37 = -t110 * t61 + t68 * t150 + t76 * t154;
t63 = t110 * t128 + t114 * t153;
t82 = -t112 * t152 + t117 * t116;
t29 = t82 * pkin(3) - t63 * t188 + t37;
t47 = -t112 * t68 + t116 * t76;
t62 = -t110 * t153 + t114 * t128;
t31 = -t62 * pkin(3) - t63 * t189 + t47;
t135 = t111 * t31 + t115 * t29;
t180 = t115 * t62;
t134 = t111 * t82 + t180;
t38 = t114 * t61 + t76 * t155 + t68 * t158;
t25 = pkin(11) * t134 + t38;
t13 = -t120 * t25 + t124 * t135;
t119 = sin(qJ(5));
t123 = cos(qJ(5));
t42 = t120 * t134 + t63 * t124;
t49 = -t62 * t111 + t82 * t115;
t26 = t42 * t119 - t49 * t123;
t205 = t26 ^ 2;
t156 = t111 * t124;
t40 = t63 * t120 - t124 * t180 - t82 * t156;
t204 = t40 ^ 2;
t157 = t111 * t120;
t66 = t116 * t157 + (t110 * t124 + t120 * t151) * t112;
t81 = -t111 * t154 + t115 * t116;
t50 = t119 * t66 - t123 * t81;
t203 = t50 ^ 2;
t64 = -t116 * t156 + t120 * t155 - t124 * t140;
t202 = t64 ^ 2;
t85 = -t123 * t115 + t119 * t157;
t201 = t85 ^ 2;
t200 = -0.2e1 * t26;
t199 = -0.2e1 * t40;
t198 = -0.2e1 * t50;
t197 = -0.2e1 * t64;
t196 = 0.2e1 * t82;
t195 = 0.2e1 * t112;
t194 = 0.2e1 * t113;
t193 = -0.2e1 * t119;
t192 = 0.2e1 * t123;
t122 = cos(qJ(6));
t190 = pkin(5) * t122;
t107 = t119 ^ 2;
t187 = t107 * pkin(12);
t186 = t112 * pkin(2);
t185 = t119 * pkin(12);
t184 = t26 * t50;
t118 = sin(qJ(6));
t14 = t120 * t135 + t124 * t25;
t11 = t49 * pkin(12) + t14;
t15 = -t111 * t29 + t115 * t31;
t12 = t40 * pkin(4) - t42 * pkin(12) + t15;
t5 = -t119 * t11 + t123 * t12;
t3 = -t40 * pkin(5) - t5;
t183 = t3 * t118;
t182 = t3 * t122;
t181 = t40 * t64;
t179 = t118 * t26;
t178 = t118 * t50;
t177 = t119 * t40;
t176 = t119 * t64;
t175 = t122 * t26;
t174 = t122 * t50;
t173 = t123 * t40;
t172 = t123 * t64;
t46 = -t111 * t67 + t115 * t75;
t32 = t64 * pkin(4) - t66 * pkin(12) + t46;
t36 = t120 * t133 + t124 * t60;
t34 = t81 * pkin(12) + t36;
t21 = -t119 * t34 + t123 * t32;
t16 = -t64 * pkin(5) - t21;
t171 = t16 * t118;
t170 = t16 * t122;
t28 = t49 * t119 + t42 * t123;
t18 = t28 * t118 - t40 * t122;
t169 = t18 * t122;
t20 = t40 * t118 + t28 * t122;
t168 = t20 * t118;
t167 = t26 * t123;
t166 = t28 * t119;
t52 = t119 * t81 + t123 * t66;
t43 = t118 * t52 - t122 * t64;
t165 = t43 * t122;
t45 = t118 * t64 + t122 * t52;
t164 = t45 * t118;
t163 = t50 * t123;
t162 = t52 * t119;
t161 = t85 * t119;
t104 = t112 ^ 2;
t160 = t104 * t114;
t105 = t113 ^ 2;
t159 = t105 * t125;
t149 = t118 * t119;
t148 = t118 * t122;
t147 = t118 * t123;
t146 = t122 * t119;
t145 = t122 * t123;
t106 = t118 ^ 2;
t108 = t122 ^ 2;
t144 = t106 + t108;
t143 = t116 * t195;
t142 = t117 * t194;
t141 = t119 * t192;
t139 = t118 * t146;
t6 = t123 * t11 + t119 * t12;
t4 = t40 * pkin(13) + t6;
t10 = -t49 * pkin(4) - t13;
t7 = t26 * pkin(5) - t28 * pkin(13) + t10;
t1 = -t118 * t4 + t122 * t7;
t2 = t118 * t7 + t122 * t4;
t138 = -t1 * t118 + t2 * t122;
t22 = t119 * t32 + t123 * t34;
t17 = t64 * pkin(13) + t22;
t33 = -t81 * pkin(4) - t35;
t23 = t50 * pkin(5) - t52 * pkin(13) + t33;
t8 = -t118 * t17 + t122 * t23;
t9 = t118 * t23 + t122 * t17;
t137 = -t8 * t118 + t9 * t122;
t136 = -t5 * t119 + t6 * t123;
t87 = t119 * t115 + t123 * t157;
t71 = -t118 * t87 - t122 * t156;
t72 = -t118 * t156 + t122 * t87;
t132 = -t71 * t118 + t72 * t122;
t91 = -t123 * pkin(5) - t119 * pkin(13) - pkin(4);
t78 = -pkin(12) * t147 + t122 * t91;
t79 = pkin(12) * t145 + t118 * t91;
t131 = -t78 * t118 + t79 * t122;
t130 = -t21 * t119 + t22 * t123;
t129 = t87 * t123 + t161;
t127 = pkin(12) ^ 2;
t109 = t123 ^ 2;
t103 = t111 ^ 2;
t102 = t107 * t127;
t100 = t103 * t124 ^ 2;
t88 = -pkin(10) * t153 + t99;
t83 = -qJ(3) * t155 + t98;
t19 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t105 * t121 ^ 2, 0.2e1 * t121 * t159, t121 * t142, t105 * t125 ^ 2, t125 * t142, t117 ^ 2, 0.2e1 * pkin(1) * t159 + 0.2e1 * t88 * t117, -0.2e1 * t105 * t191 - 0.2e1 * t89 * t117 (-t121 * t88 + t125 * t89) * t194, t105 * pkin(1) ^ 2 + t88 ^ 2 + t89 ^ 2, t63 ^ 2, 0.2e1 * t62 * t63, t63 * t196, t62 ^ 2, t62 * t196, t82 ^ 2, 0.2e1 * t37 * t82 - 0.2e1 * t47 * t62, -0.2e1 * t38 * t82 + 0.2e1 * t47 * t63, -0.2e1 * t37 * t63 + 0.2e1 * t38 * t62, t37 ^ 2 + t38 ^ 2 + t47 ^ 2, t42 ^ 2, t42 * t199, 0.2e1 * t49 * t42, t204, t49 * t199, t49 ^ 2, 0.2e1 * t13 * t49 + 0.2e1 * t15 * t40, -0.2e1 * t14 * t49 + 0.2e1 * t15 * t42, -0.2e1 * t13 * t42 - 0.2e1 * t14 * t40, t13 ^ 2 + t14 ^ 2 + t15 ^ 2, t28 ^ 2, t28 * t200, 0.2e1 * t28 * t40, t205, t26 * t199, t204, 0.2e1 * t10 * t26 + 0.2e1 * t40 * t5, 0.2e1 * t10 * t28 - 0.2e1 * t40 * t6, -0.2e1 * t26 * t6 - 0.2e1 * t28 * t5, t10 ^ 2 + t5 ^ 2 + t6 ^ 2, t20 ^ 2, -0.2e1 * t20 * t18, 0.2e1 * t20 * t26, t18 ^ 2, t18 * t200, t205, 0.2e1 * t1 * t26 + 0.2e1 * t18 * t3, -0.2e1 * t2 * t26 + 0.2e1 * t20 * t3, -0.2e1 * t1 * t20 - 0.2e1 * t18 * t2, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, 0, t152, t117, t88, -t89, 0, 0, t63 * t155 (t110 * t62 + t114 * t63) * t112, t116 * t63 + t82 * t155, t62 * t154, t62 * t116 + t82 * t154, t82 * t116, t37 * t116 + t83 * t82 + (pkin(2) * t62 - t114 * t47) * t112, -t38 * t116 - t84 * t82 + (-pkin(2) * t63 + t110 * t47) * t112, t84 * t62 - t83 * t63 + (-t110 * t37 + t114 * t38) * t112, -t47 * t186 + t37 * t83 + t38 * t84, t42 * t66, -t40 * t66 - t64 * t42, t81 * t42 + t49 * t66, t181, -t81 * t40 - t49 * t64, t49 * t81, t13 * t81 + t15 * t64 + t35 * t49 + t46 * t40, -t14 * t81 + t15 * t66 - t36 * t49 + t46 * t42, -t13 * t66 - t14 * t64 - t35 * t42 - t36 * t40, t13 * t35 + t14 * t36 + t15 * t46, t28 * t52, -t52 * t26 - t28 * t50, t28 * t64 + t52 * t40, t184, -t26 * t64 - t50 * t40, t181, t10 * t50 + t21 * t40 + t33 * t26 + t5 * t64, t10 * t52 - t22 * t40 + t33 * t28 - t6 * t64, -t21 * t28 - t22 * t26 - t5 * t52 - t6 * t50, t10 * t33 + t5 * t21 + t6 * t22, t20 * t45, -t18 * t45 - t20 * t43, t20 * t50 + t45 * t26, t18 * t43, -t18 * t50 - t43 * t26, t184, t1 * t50 + t16 * t18 + t8 * t26 + t3 * t43, t16 * t20 - t2 * t50 - t9 * t26 + t3 * t45, -t1 * t45 - t9 * t18 - t2 * t43 - t8 * t20, t1 * t8 + t16 * t3 + t2 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t104 * t110 ^ 2, 0.2e1 * t110 * t160, t110 * t143, t104 * t114 ^ 2, t114 * t143, t116 ^ 2, 0.2e1 * pkin(2) * t160 + 0.2e1 * t83 * t116, -0.2e1 * t104 * pkin(2) * t110 - 0.2e1 * t84 * t116 (-t110 * t83 + t114 * t84) * t195, t104 * pkin(2) ^ 2 + t83 ^ 2 + t84 ^ 2, t66 ^ 2, t66 * t197, 0.2e1 * t81 * t66, t202, t81 * t197, t81 ^ 2, 0.2e1 * t35 * t81 + 0.2e1 * t46 * t64, -0.2e1 * t36 * t81 + 0.2e1 * t46 * t66, -0.2e1 * t35 * t66 - 0.2e1 * t36 * t64, t35 ^ 2 + t36 ^ 2 + t46 ^ 2, t52 ^ 2, t52 * t198, 0.2e1 * t52 * t64, t203, t50 * t197, t202, 0.2e1 * t21 * t64 + 0.2e1 * t33 * t50, -0.2e1 * t22 * t64 + 0.2e1 * t33 * t52, -0.2e1 * t21 * t52 - 0.2e1 * t22 * t50, t21 ^ 2 + t22 ^ 2 + t33 ^ 2, t45 ^ 2, -0.2e1 * t45 * t43, 0.2e1 * t45 * t50, t43 ^ 2, t43 * t198, t203, 0.2e1 * t16 * t43 + 0.2e1 * t8 * t50, 0.2e1 * t16 * t45 - 0.2e1 * t9 * t50, -0.2e1 * t9 * t43 - 0.2e1 * t8 * t45, t16 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, t63, 0, t47, 0, 0, 0, 0, 0, 0, t115 * t40 + t49 * t156, t115 * t42 - t157 * t49 (-t120 * t40 - t124 * t42) * t111, t15 * t115 + (t120 * t14 + t124 * t13) * t111, 0, 0, 0, 0, 0, 0, -t156 * t26 - t85 * t40, -t156 * t28 - t87 * t40, -t87 * t26 + t85 * t28, -t10 * t156 - t5 * t85 + t6 * t87, 0, 0, 0, 0, 0, 0, t85 * t18 + t71 * t26, t85 * t20 - t72 * t26, -t72 * t18 - t71 * t20, t1 * t71 + t2 * t72 + t3 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t154, t155, 0, -t186, 0, 0, 0, 0, 0, 0, t115 * t64 + t81 * t156, t115 * t66 - t157 * t81 (-t120 * t64 - t124 * t66) * t111, t46 * t115 + (t120 * t36 + t124 * t35) * t111, 0, 0, 0, 0, 0, 0, -t156 * t50 - t85 * t64, -t156 * t52 - t87 * t64, -t87 * t50 + t85 * t52, -t156 * t33 - t21 * t85 + t22 * t87, 0, 0, 0, 0, 0, 0, t85 * t43 + t71 * t50, t85 * t45 - t72 * t50, -t72 * t43 - t71 * t45, t16 * t85 + t8 * t71 + t9 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103 * t120 ^ 2 + t115 ^ 2 + t100, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87 ^ 2 + t100 + t201, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 ^ 2 + t72 ^ 2 + t201; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, -t40, t49, t13, -t14, 0, 0, t166, -t119 * t26 + t28 * t123, t177, -t167, t173, 0, -pkin(4) * t26 - pkin(12) * t177 - t10 * t123, -pkin(4) * t28 - pkin(12) * t173 + t10 * t119 (t166 - t167) * pkin(12) + t136, -t10 * pkin(4) + pkin(12) * t136, t20 * t146 (-t168 - t169) * t119, -t20 * t123 + t146 * t26, t18 * t149, t18 * t123 - t149 * t26, -t167, -t1 * t123 + t78 * t26 + (pkin(12) * t18 + t183) * t119, t2 * t123 - t79 * t26 + (pkin(12) * t20 + t182) * t119, -t79 * t18 - t78 * t20 + (-t1 * t122 - t118 * t2) * t119, t1 * t78 + t3 * t185 + t2 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, -t64, t81, t35, -t36, 0, 0, t162, -t119 * t50 + t52 * t123, t176, -t163, t172, 0, -pkin(4) * t50 - pkin(12) * t176 - t33 * t123, -pkin(4) * t52 - pkin(12) * t172 + t33 * t119 (t162 - t163) * pkin(12) + t130, -t33 * pkin(4) + pkin(12) * t130, t45 * t146 (-t164 - t165) * t119, -t45 * t123 + t146 * t50, t43 * t149, t43 * t123 - t149 * t50, -t163, -t8 * t123 + t78 * t50 + (pkin(12) * t43 + t171) * t119, t9 * t123 - t79 * t50 + (pkin(12) * t45 + t170) * t119, -t79 * t43 - t78 * t45 + (-t118 * t9 - t122 * t8) * t119, t16 * t185 + t8 * t78 + t9 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, -t157, 0, 0, 0, 0, 0, 0, 0, 0, t123 * t156, -t119 * t156, t129, pkin(4) * t156 + pkin(12) * t129, 0, 0, 0, 0, 0, 0, -t71 * t123 + t149 * t85, t72 * t123 + t146 * t85 (-t118 * t72 - t122 * t71) * t119, pkin(12) * t161 + t71 * t78 + t72 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t107, t141, 0, t109, 0, 0, pkin(4) * t192, pkin(4) * t193, 0.2e1 * (t107 + t109) * pkin(12), pkin(4) ^ 2 + t109 * t127 + t102, t108 * t107, -0.2e1 * t107 * t148, t145 * t193, t106 * t107, t118 * t141, t109, 0.2e1 * t118 * t187 - 0.2e1 * t78 * t123, 0.2e1 * t122 * t187 + 0.2e1 * t79 * t123, 0.2e1 * (-t118 * t79 - t122 * t78) * t119, t78 ^ 2 + t79 ^ 2 + t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t26, t40, t5, -t6, 0, 0, t168, -t118 * t18 + t20 * t122, t179, -t169, t175, 0, -pkin(5) * t18 - pkin(13) * t179 - t182, -pkin(5) * t20 - pkin(13) * t175 + t183 (t168 - t169) * pkin(13) + t138, -t3 * pkin(5) + pkin(13) * t138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52, 0, -t50, t64, t21, -t22, 0, 0, t164, -t118 * t43 + t45 * t122, t178, -t165, t174, 0, -pkin(5) * t43 - pkin(13) * t178 - t170, -pkin(5) * t45 - pkin(13) * t174 + t171 (t164 - t165) * pkin(13) + t137, -t16 * pkin(5) + pkin(13) * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t87, 0, 0, 0, 0, 0, 0, 0, 0, -t85 * t122, t85 * t118, t132, -t85 * pkin(5) + pkin(13) * t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, 0, t123, 0, -t185, -t123 * pkin(12), 0, 0, t139 (-t106 + t108) * t119, -t147, -t139, -t145, 0, -pkin(12) * t146 + (-pkin(5) * t119 + pkin(13) * t123) * t118, pkin(13) * t145 + (pkin(12) * t118 - t190) * t119, t131, -pkin(5) * t185 + pkin(13) * t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t106, 0.2e1 * t148, 0, t108, 0, 0, 0.2e1 * t190, -0.2e1 * pkin(5) * t118, 0.2e1 * t144 * pkin(13), pkin(13) ^ 2 * t144 + pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, -t18, t26, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, -t43, t50, t8, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t72, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, 0, -t149, -t123, t78, -t79, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, 0, t122, 0, -t118 * pkin(13), -t122 * pkin(13), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t19;
