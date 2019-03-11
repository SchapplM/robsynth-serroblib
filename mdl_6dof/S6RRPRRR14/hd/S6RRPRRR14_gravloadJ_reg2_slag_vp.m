% Calculate inertial parameters regressor of gravitation load for
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR14_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR14_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_gravloadJ_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t102 = sin(qJ(6));
t103 = sin(qJ(5));
t106 = cos(qJ(5));
t104 = sin(qJ(4));
t198 = cos(qJ(4));
t100 = sin(pkin(8));
t101 = sin(pkin(7));
t107 = cos(qJ(1));
t187 = sin(pkin(6));
t190 = cos(pkin(7));
t157 = t190 * t187;
t148 = t107 * t157;
t191 = cos(pkin(6));
t199 = cos(qJ(2));
t165 = t191 * t199;
t196 = sin(qJ(2));
t197 = sin(qJ(1));
t86 = -t107 * t165 + t197 * t196;
t133 = t86 * t101 - t148;
t189 = cos(pkin(8));
t188 = cos(pkin(14));
t155 = t188 * t187;
t150 = t101 * t155;
t158 = t190 * t188;
t186 = sin(pkin(14));
t164 = t191 * t196;
t87 = t107 * t164 + t197 * t199;
t56 = t107 * t150 + t86 * t158 + t87 * t186;
t203 = -t133 * t100 + t189 * t56;
t154 = t187 * t186;
t149 = t101 * t154;
t156 = t190 * t186;
t55 = -t107 * t149 - t86 * t156 + t87 * t188;
t23 = t203 * t104 - t55 * t198;
t39 = t56 * t100 + t133 * t189;
t5 = t39 * t103 - t23 * t106;
t213 = t5 * t102;
t105 = cos(qJ(6));
t212 = t5 * t105;
t211 = t23 * t103 + t39 * t106;
t210 = t55 * t104;
t129 = t107 * t196 + t197 * t165;
t125 = t129 * t188;
t88 = t107 * t199 - t197 * t164;
t113 = t190 * t125 - t197 * t150 + t88 * t186;
t206 = t129 * t101 + t197 * t157;
t207 = t113 * t100 + t189 * t206;
t205 = -t100 * t206 + t113 * t189;
t138 = t190 * t155;
t175 = t101 * t191;
t115 = -t199 * t138 + t196 * t154 - t188 * t175;
t162 = t199 * t187;
t128 = t101 * t162 - t191 * t190;
t204 = t128 * t100 + t115 * t189;
t174 = t101 * t189;
t63 = -t87 * t158 + t86 * t186;
t47 = -t63 * t100 + t87 * t174;
t124 = t129 * t186;
t65 = -t88 * t158 + t124;
t48 = -t65 * t100 + t88 * t174;
t202 = -g(1) * t48 - g(2) * t47;
t160 = t187 * t196;
t151 = t101 * t160;
t195 = pkin(2) * t162 + qJ(3) * t151;
t161 = t187 * t197;
t194 = t107 * pkin(1) + pkin(10) * t161;
t80 = -t196 * t138 - t199 * t154;
t192 = t80 * t100;
t185 = qJ(3) * t101;
t184 = t101 * t100;
t183 = t102 * t106;
t182 = t105 * t106;
t20 = t203 * t198 + t210;
t181 = -t20 * pkin(4) - pkin(12) * t23;
t58 = -t190 * t124 + t197 * t149 + t88 * t188;
t24 = t58 * t104 + t205 * t198;
t25 = -t205 * t104 + t58 * t198;
t180 = -t24 * pkin(4) + t25 * pkin(12);
t137 = t190 * t154;
t75 = t199 * t137 + t196 * t155 + t186 * t175;
t35 = t75 * t104 + t204 * t198;
t36 = -t204 * t104 + t75 * t198;
t179 = -t35 * pkin(4) + t36 * pkin(12);
t178 = t100 * t198;
t173 = t107 * t187;
t136 = t189 * t151;
t81 = -t196 * t137 + t199 * t155;
t172 = t81 * pkin(3) + pkin(11) * t136 + t195;
t8 = t25 * t103 - t106 * t207;
t171 = g(1) * t211 + g(2) * t8;
t170 = t101 * t178;
t169 = -t197 * pkin(1) + pkin(10) * t173;
t163 = t189 * t198;
t22 = t133 * t178 - t56 * t163 - t210;
t168 = g(1) * t22 + g(2) * t24;
t167 = -t86 * pkin(2) + t87 * t185;
t166 = -t129 * pkin(2) + t88 * t185;
t159 = -pkin(5) * t106 - pkin(13) * t103;
t64 = -t87 * t156 - t86 * t188;
t153 = t64 * pkin(3) + t167;
t66 = -t88 * t156 - t125;
t152 = t66 * pkin(3) + t166;
t52 = t115 * t100 - t128 * t189;
t14 = -t36 * t103 + t52 * t106;
t147 = g(1) * t8 - g(2) * t211 - g(3) * t14;
t15 = t52 * t103 + t36 * t106;
t9 = t103 * t207 + t25 * t106;
t146 = g(1) * t9 + g(2) * t5 + g(3) * t15;
t29 = t64 * t198 + (t87 * t184 + t189 * t63) * t104;
t10 = t29 * t103 - t47 * t106;
t31 = t66 * t198 + (t88 * t184 + t189 * t65) * t104;
t12 = t31 * t103 - t48 * t106;
t140 = t100 * t151;
t46 = t81 * t198 + (t189 * t80 + t140) * t104;
t68 = t136 - t192;
t32 = t46 * t103 - t68 * t106;
t145 = g(1) * t12 + g(2) * t10 + g(3) * t32;
t144 = g(1) * t24 + g(2) * t20 + g(3) * t35;
t143 = g(1) * t25 - g(2) * t23 + g(3) * t36;
t28 = t64 * t104 - t63 * t163 - t87 * t170;
t30 = t66 * t104 - t65 * t163 - t88 * t170;
t45 = t81 * t104 - t198 * t140 - t80 * t163;
t142 = g(1) * t30 + g(2) * t28 + g(3) * t45;
t141 = t46 * pkin(4) + t45 * pkin(12) + t172;
t135 = t29 * pkin(4) + t28 * pkin(12) + t153;
t134 = t31 * pkin(4) + t30 * pkin(12) + t152;
t132 = -t87 * pkin(2) + qJ(3) * t148 - t86 * t185 + t169;
t131 = -g(1) * t88 - g(2) * t87 - g(3) * t160;
t123 = -t55 * pkin(3) - pkin(11) * t39 + t132;
t122 = t23 * pkin(4) + t22 * pkin(12) + t123;
t121 = t88 * pkin(2) + t206 * qJ(3) + t194;
t117 = (g(3) * t192 + t202) * pkin(11);
t110 = t58 * pkin(3) + t207 * pkin(11) + t121;
t109 = t25 * pkin(4) + t24 * pkin(12) + t110;
t41 = -g(1) * t206 - g(2) * t133 + g(3) * t128;
t33 = t68 * t103 + t46 * t106;
t13 = t48 * t103 + t31 * t106;
t11 = t47 * t103 + t29 * t106;
t3 = t24 * t102 + t9 * t105;
t2 = -t9 * t102 + t24 * t105;
t1 = t144 * t103;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t197 - g(2) * t107, g(1) * t107 + g(2) * t197, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t87 - g(2) * t88, -g(1) * t86 + g(2) * t129, -g(1) * t173 - g(2) * t161, -g(1) * t169 - g(2) * t194, 0, 0, 0, 0, 0, 0, g(1) * t55 - g(2) * t58, -g(1) * t56 + g(2) * t113, g(1) * t133 - g(2) * t206, -g(1) * t132 - g(2) * t121, 0, 0, 0, 0, 0, 0, -g(1) * t23 - g(2) * t25, t168, g(1) * t39 - g(2) * t207, -g(1) * t123 - g(2) * t110, 0, 0, 0, 0, 0, 0, g(1) * t5 - g(2) * t9, t171, -t168, -g(1) * t122 - g(2) * t109, 0, 0, 0, 0, 0, 0, -g(1) * (t22 * t102 - t212) - g(2) * t3, -g(1) * (t22 * t105 + t213) - g(2) * t2, -t171, -g(1) * (-pkin(5) * t5 + pkin(13) * t211 + t122) - g(2) * (t9 * pkin(5) + t8 * pkin(13) + t109); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t129 + g(2) * t86 - g(3) * t162, -t131, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t66 - g(2) * t64 - g(3) * t81, -g(1) * t65 - g(2) * t63 - g(3) * t80, t131 * t101, -g(1) * t166 - g(2) * t167 - g(3) * t195, 0, 0, 0, 0, 0, 0, -g(1) * t31 - g(2) * t29 - g(3) * t46, t142, -g(3) * t68 + t202, -g(1) * t152 - g(2) * t153 - g(3) * t172 + t117, 0, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t11 - g(3) * t33, t145, -t142, -g(1) * t134 - g(2) * t135 - g(3) * t141 + t117, 0, 0, 0, 0, 0, 0, -g(1) * (t30 * t102 + t13 * t105) - g(2) * (t28 * t102 + t11 * t105) - g(3) * (t45 * t102 + t33 * t105) -g(1) * (-t13 * t102 + t30 * t105) - g(2) * (-t11 * t102 + t28 * t105) - g(3) * (-t33 * t102 + t45 * t105) -t145, -g(1) * (t13 * pkin(5) + t12 * pkin(13) + t134) - g(2) * (t11 * pkin(5) + t10 * pkin(13) + t135) - g(3) * (t33 * pkin(5) + t32 * pkin(13) + t141) + t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, t143, 0, 0, 0, 0, 0, 0, 0, 0, t144 * t106, -t1, -t143, -g(1) * t180 - g(2) * t181 - g(3) * t179, 0, 0, 0, 0, 0, 0, -g(1) * (t25 * t102 - t24 * t182) - g(2) * (-t102 * t23 - t20 * t182) - g(3) * (t36 * t102 - t35 * t182) -g(1) * (t25 * t105 + t24 * t183) - g(2) * (-t105 * t23 + t20 * t183) - g(3) * (t36 * t105 + t35 * t183) t1, -g(1) * (t159 * t24 + t180) - g(2) * (t159 * t20 + t181) - g(3) * (t159 * t35 + t179); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, t146, 0, 0, 0, 0, 0, 0, 0, 0, t147 * t105, -t147 * t102, -t146, -g(1) * (-t8 * pkin(5) + t9 * pkin(13)) - g(2) * (pkin(5) * t211 + t5 * pkin(13)) - g(3) * (t14 * pkin(5) + t15 * pkin(13)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * (t20 * t105 - t213) - g(3) * (-t15 * t102 + t35 * t105) g(1) * t3 - g(2) * (-t20 * t102 - t212) - g(3) * (-t35 * t102 - t15 * t105) 0, 0;];
taug_reg  = t4;
