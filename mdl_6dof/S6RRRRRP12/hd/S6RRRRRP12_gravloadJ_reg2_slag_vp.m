% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRP12_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP12_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_gravloadJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t121 = sin(qJ(5));
t124 = cos(qJ(5));
t122 = sin(qJ(4));
t125 = cos(qJ(4));
t126 = cos(qJ(1));
t195 = cos(pkin(6));
t204 = cos(qJ(2));
t170 = t195 * t204;
t201 = sin(qJ(2));
t202 = sin(qJ(1));
t104 = -t126 * t170 + t202 * t201;
t120 = sin(pkin(7));
t193 = sin(pkin(6));
t194 = cos(pkin(7));
t160 = t194 * t193;
t209 = -t104 * t120 + t126 * t160;
t169 = t195 * t201;
t105 = t126 * t169 + t202 * t204;
t123 = sin(qJ(3));
t178 = t126 * t193;
t179 = t123 * t194;
t203 = cos(qJ(3));
t64 = -t120 * t123 * t178 - t104 * t179 + t105 * t203;
t32 = -t122 * t209 + t64 * t125;
t166 = t203 * t193;
t168 = t194 * t203;
t63 = t126 * t120 * t166 + t104 * t168 + t105 * t123;
t11 = t32 * t121 - t63 * t124;
t12 = t63 * t121 + t32 * t124;
t135 = t126 * t201 + t202 * t170;
t210 = t135 * t120 + t202 * t160;
t208 = -pkin(4) * t125 - pkin(12) * t122;
t207 = -t64 * t122 - t125 * t209;
t165 = t193 * t202;
t206 = -t120 * t165 + t135 * t194;
t205 = t195 * t120 + t204 * t160;
t199 = pkin(10) * t120;
t191 = t120 * t122;
t190 = t120 * t125;
t189 = t121 * t125;
t188 = t124 * t125;
t164 = t193 * t201;
t157 = t120 * t164;
t167 = t204 * t193;
t187 = pkin(2) * t167 + pkin(10) * t157;
t186 = t126 * pkin(1) + pkin(9) * t165;
t185 = -t63 * pkin(3) + t64 * pkin(11);
t106 = t126 * t204 - t202 * t169;
t67 = t106 * t123 + t206 * t203;
t68 = t106 * t203 - t206 * t123;
t184 = -t67 * pkin(3) + t68 * pkin(11);
t87 = t123 * t164 - t205 * t203;
t88 = t205 * t123 + t203 * t164;
t183 = -t87 * pkin(3) + t88 * pkin(11);
t182 = pkin(4) * t207 + t32 * pkin(12);
t35 = t68 * t122 - t125 * t210;
t36 = t122 * t210 + t68 * t125;
t181 = -t35 * pkin(4) + t36 * pkin(12);
t103 = -t120 * t167 + t195 * t194;
t61 = t103 * t125 - t88 * t122;
t62 = t103 * t122 + t88 * t125;
t180 = t61 * pkin(4) + t62 * pkin(12);
t176 = -t104 * pkin(2) + t105 * t199;
t175 = -t135 * pkin(2) + t106 * t199;
t15 = t36 * t121 - t67 * t124;
t174 = -g(1) * t11 + g(2) * t15;
t173 = g(1) * t207 + g(2) * t35;
t172 = -g(1) * t63 + g(2) * t67;
t171 = -t202 * pkin(1) + pkin(9) * t178;
t163 = t208 * t63 + t185;
t162 = t208 * t67 + t184;
t161 = t208 * t87 + t183;
t159 = pkin(5) * t124 + qJ(6) * t121;
t145 = t201 * t160;
t97 = t123 * t167 + t203 * t145;
t98 = -t123 * t145 + t204 * t166;
t158 = t98 * pkin(3) + t97 * pkin(11) + t187;
t25 = t62 * t121 - t87 * t124;
t1 = g(1) * t15 + g(2) * t11 + g(3) * t25;
t16 = t67 * t121 + t36 * t124;
t26 = t87 * t121 + t62 * t124;
t154 = g(1) * t16 + g(2) * t12 + g(3) * t26;
t17 = -t64 * t124 - t63 * t189;
t19 = -t68 * t124 - t67 * t189;
t37 = -t88 * t124 - t87 * t189;
t153 = g(1) * t19 + g(2) * t17 + g(3) * t37;
t74 = -t104 * t203 - t105 * t179;
t42 = t105 * t191 + t74 * t125;
t73 = -t104 * t123 + t105 * t168;
t21 = t42 * t121 - t73 * t124;
t76 = -t106 * t179 - t135 * t203;
t44 = t106 * t191 + t76 * t125;
t75 = t106 * t168 - t135 * t123;
t23 = t44 * t121 - t75 * t124;
t79 = t122 * t157 + t98 * t125;
t45 = t79 * t121 - t97 * t124;
t152 = g(1) * t23 + g(2) * t21 + g(3) * t45;
t151 = g(1) * t35 - g(2) * t207 - g(3) * t61;
t150 = g(1) * t36 + g(2) * t32 + g(3) * t62;
t41 = -t105 * t190 + t74 * t122;
t43 = -t106 * t190 + t76 * t122;
t78 = t98 * t122 - t125 * t157;
t149 = g(1) * t43 + g(2) * t41 + g(3) * t78;
t148 = g(1) * t67 + g(2) * t63 + g(3) * t87;
t147 = g(1) * t68 + g(2) * t64 + g(3) * t88;
t146 = g(1) * t75 + g(2) * t73 + g(3) * t97;
t142 = t74 * pkin(3) + t73 * pkin(11) + t176;
t141 = t76 * pkin(3) + t75 * pkin(11) + t175;
t140 = t79 * pkin(4) + t78 * pkin(12) + t158;
t139 = -t105 * pkin(2) + t209 * pkin(10) + t171;
t138 = t42 * pkin(4) + t41 * pkin(12) + t142;
t137 = t44 * pkin(4) + t43 * pkin(12) + t141;
t136 = -g(1) * t106 - g(2) * t105 - g(3) * t164;
t133 = -pkin(3) * t64 - pkin(11) * t63 + t139;
t131 = -pkin(4) * t32 + pkin(12) * t207 + t133;
t130 = t106 * pkin(2) + t210 * pkin(10) + t186;
t128 = t68 * pkin(3) + t67 * pkin(11) + t130;
t127 = t36 * pkin(4) + t35 * pkin(12) + t128;
t46 = t97 * t121 + t79 * t124;
t38 = t88 * t121 - t87 * t188;
t24 = t75 * t121 + t44 * t124;
t22 = t73 * t121 + t42 * t124;
t20 = t68 * t121 - t67 * t188;
t18 = t64 * t121 - t63 * t188;
t10 = t148 * t122;
t6 = t151 * t124;
t5 = t151 * t121;
t4 = -g(1) * t24 - g(2) * t22 - g(3) * t46;
t3 = g(1) * t12 - g(2) * t16;
t2 = -g(1) * t20 - g(2) * t18 - g(3) * t38;
t7 = [0, 0, 0, 0, 0, 0, g(1) * t202 - g(2) * t126, g(1) * t126 + g(2) * t202, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t105 - g(2) * t106, -g(1) * t104 + g(2) * t135, -g(1) * t178 - g(2) * t165, -g(1) * t171 - g(2) * t186, 0, 0, 0, 0, 0, 0, g(1) * t64 - g(2) * t68, t172, -g(1) * t209 - g(2) * t210, -g(1) * t139 - g(2) * t130, 0, 0, 0, 0, 0, 0, g(1) * t32 - g(2) * t36, t173, -t172, -g(1) * t133 - g(2) * t128, 0, 0, 0, 0, 0, 0, t3, t174, -t173, -g(1) * t131 - g(2) * t127, 0, 0, 0, 0, 0, 0, t3, -t173, -t174, -g(1) * (-pkin(5) * t12 - qJ(6) * t11 + t131) - g(2) * (t16 * pkin(5) + t15 * qJ(6) + t127); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t135 + g(2) * t104 - g(3) * t167, -t136, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t76 - g(2) * t74 - g(3) * t98, t146, t136 * t120, -g(1) * t175 - g(2) * t176 - g(3) * t187, 0, 0, 0, 0, 0, 0, -g(1) * t44 - g(2) * t42 - g(3) * t79, t149, -t146, -g(1) * t141 - g(2) * t142 - g(3) * t158, 0, 0, 0, 0, 0, 0, t4, t152, -t149, -g(1) * t137 - g(2) * t138 - g(3) * t140, 0, 0, 0, 0, 0, 0, t4, -t149, -t152, -g(1) * (t24 * pkin(5) + t23 * qJ(6) + t137) - g(2) * (t22 * pkin(5) + t21 * qJ(6) + t138) - g(3) * (t46 * pkin(5) + t45 * qJ(6) + t140); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, t147, 0, 0, 0, 0, 0, 0, 0, 0, t148 * t125, -t10, -t147, -g(1) * t184 - g(2) * t185 - g(3) * t183, 0, 0, 0, 0, 0, 0, t2, t153, t10, -g(1) * t162 - g(2) * t163 - g(3) * t161, 0, 0, 0, 0, 0, 0, t2, t10, -t153, -g(1) * (t20 * pkin(5) + t19 * qJ(6) + t162) - g(2) * (t18 * pkin(5) + t17 * qJ(6) + t163) - g(3) * (t38 * pkin(5) + t37 * qJ(6) + t161); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, t150, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t150, -g(1) * t181 - g(2) * t182 - g(3) * t180, 0, 0, 0, 0, 0, 0, t6, -t150, t5, -g(1) * (-t159 * t35 + t181) - g(2) * (t159 * t207 + t182) - g(3) * (t159 * t61 + t180); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t154, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t154, -g(1) * (-t15 * pkin(5) + t16 * qJ(6)) - g(2) * (-pkin(5) * t11 + qJ(6) * t12) - g(3) * (-t25 * pkin(5) + t26 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t7;
