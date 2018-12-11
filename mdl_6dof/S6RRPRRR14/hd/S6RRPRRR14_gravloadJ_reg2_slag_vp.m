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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

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
t118 = sin(qJ(6));
t121 = cos(qJ(6));
t120 = sin(qJ(4));
t115 = sin(pkin(7));
t216 = sin(pkin(6));
t217 = cos(pkin(7));
t176 = t217 * t216;
t224 = cos(qJ(1));
t162 = t224 * t176;
t207 = pkin(6) + qJ(2);
t182 = cos(t207) / 0.2e1;
t208 = pkin(6) - qJ(2);
t199 = cos(t208);
t156 = t199 / 0.2e1 + t182;
t221 = sin(qJ(2));
t222 = sin(qJ(1));
t95 = -t156 * t224 + t221 * t222;
t148 = t115 * t95 - t162;
t205 = pkin(8) + qJ(4);
t178 = sin(t205) / 0.2e1;
t206 = pkin(8) - qJ(4);
t195 = sin(t206);
t152 = t178 + t195 / 0.2e1;
t181 = cos(t205) / 0.2e1;
t198 = cos(t206);
t155 = t198 / 0.2e1 + t181;
t196 = sin(t207);
t179 = t196 / 0.2e1;
t197 = sin(t208);
t154 = t179 - t197 / 0.2e1;
t223 = cos(qJ(2));
t190 = t222 * t223;
t141 = t154 * t224 + t190;
t203 = pkin(7) + pkin(14);
t174 = sin(t203) / 0.2e1;
t204 = pkin(7) - pkin(14);
t193 = sin(t204);
t150 = t174 + t193 / 0.2e1;
t146 = t150 * t216;
t175 = cos(t204) / 0.2e1;
t194 = cos(t203);
t151 = t175 + t194 / 0.2e1;
t215 = sin(pkin(14));
t228 = t141 * t215 + t146 * t224 + t151 * t95;
t100 = t175 - t194 / 0.2e1;
t116 = cos(pkin(14));
t184 = t224 * t216;
t99 = t174 - t193 / 0.2e1;
t60 = -t100 * t184 + t116 * t141 - t95 * t99;
t20 = t120 * t60 - t148 * t152 + t155 * t228;
t119 = sin(qJ(5));
t122 = cos(qJ(5));
t101 = t178 - t195 / 0.2e1;
t102 = t181 - t198 / 0.2e1;
t123 = cos(qJ(4));
t130 = -t101 * t228 - t102 * t148 + t60 * t123;
t114 = sin(pkin(8));
t117 = cos(pkin(8));
t45 = t114 * t228 + t117 * t148;
t5 = t119 * t45 + t122 * t130;
t238 = t118 * t5 - t121 * t20;
t237 = t118 * t20 + t121 * t5;
t234 = -t119 * t130 + t122 * t45;
t142 = t156 * t222 + t221 * t224;
t191 = t224 * t223;
t143 = -t154 * t222 + t191;
t129 = t142 * t151 + t143 * t215 - t146 * t222;
t232 = t115 * t142 + t176 * t222;
t233 = t114 * t129 + t117 * t232;
t212 = t115 * t117;
t180 = t197 / 0.2e1;
t170 = t180 - t196 / 0.2e1;
t96 = t170 * t224 - t190;
t72 = t151 * t96 + t215 * t95;
t50 = -t114 * t72 - t212 * t96;
t97 = -t170 * t222 - t191;
t74 = t142 * t215 + t151 * t97;
t51 = -t114 * t74 - t212 * t97;
t103 = t182 - t199 / 0.2e1;
t153 = t179 + t180;
t84 = t103 * t151 - t153 * t215;
t69 = -t103 * t212 - t114 * t84;
t229 = -g(1) * t51 - g(2) * t50 - g(3) * t69;
t218 = cos(pkin(6));
t214 = qJ(3) * t115;
t213 = t115 * t102;
t211 = t118 * t122;
t210 = t121 * t122;
t183 = t216 * t222;
t209 = pkin(1) * t224 + pkin(10) * t183;
t202 = -pkin(4) * t20 + pkin(12) * t130;
t63 = t100 * t183 + t116 * t143 - t142 * t99;
t124 = -t101 * t129 - t102 * t232 + t123 * t63;
t25 = t120 * t63 + t129 * t155 - t152 * t232;
t201 = -pkin(4) * t25 + pkin(12) * t124;
t133 = t103 * t215 + t150 * t218 + t151 * t153;
t144 = t115 * t153 - t217 * t218;
t76 = t100 * t218 - t103 * t116 + t153 * t99;
t131 = t101 * t133 + t102 * t144 + t123 * t76;
t37 = t120 * t76 - t133 * t155 + t144 * t152;
t200 = -pkin(4) * t37 + pkin(12) * t131;
t8 = t119 * t124 - t122 * t233;
t192 = g(1) * t234 + g(2) * t8;
t189 = -g(1) * t20 + g(2) * t25;
t188 = -pkin(2) * t95 - t214 * t96;
t187 = -pkin(2) * t142 - t214 * t97;
t186 = -pkin(1) * t222 + pkin(10) * t184;
t185 = pkin(2) * t153 - t103 * t214;
t177 = -pkin(5) * t122 - pkin(13) * t119;
t73 = -t116 * t95 + t96 * t99;
t173 = pkin(3) * t73 + t188;
t75 = -t116 * t142 + t97 * t99;
t172 = pkin(3) * t75 + t187;
t85 = t103 * t99 + t116 * t153;
t171 = pkin(3) * t85 + t185;
t57 = -t114 * t133 - t117 * t144;
t14 = -t119 * t131 + t122 * t57;
t169 = g(1) * t8 - g(2) * t234 - g(3) * t14;
t15 = t119 * t57 + t122 * t131;
t9 = t119 * t233 + t122 * t124;
t168 = g(1) * t9 + g(2) * t5 + g(3) * t15;
t33 = t101 * t72 + t123 * t73 + t213 * t96;
t10 = t119 * t33 - t122 * t50;
t35 = t101 * t74 + t123 * t75 + t213 * t97;
t12 = t119 * t35 - t122 * t51;
t44 = t101 * t84 + t103 * t213 + t123 * t85;
t28 = t119 * t44 - t122 * t69;
t167 = g(1) * t12 + g(2) * t10 + g(3) * t28;
t166 = g(1) * t25 + g(2) * t20 + g(3) * t37;
t165 = g(1) * t124 + g(2) * t130 + g(3) * t131;
t147 = t115 * t152;
t32 = t120 * t73 + t147 * t96 - t155 * t72;
t34 = t120 * t75 + t147 * t97 - t155 * t74;
t43 = t103 * t147 + t120 * t85 - t155 * t84;
t164 = g(1) * t34 + g(2) * t32 + g(3) * t43;
t163 = g(1) * t97 + g(2) * t96 + g(3) * t103;
t160 = pkin(4) * t33 + t32 * pkin(12) + t173;
t159 = pkin(4) * t35 + t34 * pkin(12) + t172;
t158 = pkin(4) * t44 + t43 * pkin(12) + t171;
t149 = -pkin(2) * t141 + qJ(3) * t162 - t214 * t95 + t186;
t145 = -pkin(3) * t60 - pkin(11) * t45 + t149;
t140 = -pkin(4) * t130 - pkin(12) * t20 + t145;
t138 = t143 * pkin(2) + qJ(3) * t232 + t209;
t137 = t229 * pkin(11);
t127 = t63 * pkin(3) + pkin(11) * t233 + t138;
t125 = pkin(4) * t124 + t25 * pkin(12) + t127;
t47 = -g(1) * t232 - g(2) * t148 + g(3) * t144;
t29 = t119 * t69 + t122 * t44;
t13 = t119 * t51 + t122 * t35;
t11 = t119 * t50 + t122 * t33;
t3 = t118 * t25 + t121 * t9;
t2 = -t118 * t9 + t121 * t25;
t1 = t166 * t119;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t222 - g(2) * t224, g(1) * t224 + g(2) * t222, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t141 - g(2) * t143, -g(1) * t95 + g(2) * t142, -g(1) * t184 - g(2) * t183, -g(1) * t186 - g(2) * t209, 0, 0, 0, 0, 0, 0, g(1) * t60 - g(2) * t63, -g(1) * t228 + g(2) * t129, g(1) * t148 - g(2) * t232, -g(1) * t149 - g(2) * t138, 0, 0, 0, 0, 0, 0, g(1) * t130 - g(2) * t124, t189, g(1) * t45 - g(2) * t233, -g(1) * t145 - g(2) * t127, 0, 0, 0, 0, 0, 0, g(1) * t5 - g(2) * t9, t192, -t189, -g(1) * t140 - g(2) * t125, 0, 0, 0, 0, 0, 0, g(1) * t237 - g(2) * t3, -g(1) * t238 - g(2) * t2, -t192, -g(1) * (-pkin(5) * t5 + pkin(13) * t234 + t140) - g(2) * (t9 * pkin(5) + t8 * pkin(13) + t125); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t142 + g(2) * t95 - g(3) * t153, -t163, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t75 - g(2) * t73 - g(3) * t85, -g(1) * t74 - g(2) * t72 - g(3) * t84, t163 * t115, -g(1) * t187 - g(2) * t188 - g(3) * t185, 0, 0, 0, 0, 0, 0, -g(1) * t35 - g(2) * t33 - g(3) * t44, t164, t229, -g(1) * t172 - g(2) * t173 - g(3) * t171 + t137, 0, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t11 - g(3) * t29, t167, -t164, -g(1) * t159 - g(2) * t160 - g(3) * t158 + t137, 0, 0, 0, 0, 0, 0, -g(1) * (t118 * t34 + t121 * t13) - g(2) * (t11 * t121 + t118 * t32) - g(3) * (t118 * t43 + t121 * t29) -g(1) * (-t118 * t13 + t121 * t34) - g(2) * (-t11 * t118 + t121 * t32) - g(3) * (-t118 * t29 + t121 * t43) -t167, -g(1) * (t13 * pkin(5) + t12 * pkin(13) + t159) - g(2) * (t11 * pkin(5) + t10 * pkin(13) + t160) - g(3) * (t29 * pkin(5) + t28 * pkin(13) + t158) + t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, t165, 0, 0, 0, 0, 0, 0, 0, 0, t166 * t122, -t1, -t165, -g(1) * t201 - g(2) * t202 - g(3) * t200, 0, 0, 0, 0, 0, 0, -g(1) * (t118 * t124 - t210 * t25) - g(2) * (t118 * t130 - t20 * t210) - g(3) * (t118 * t131 - t210 * t37) -g(1) * (t121 * t124 + t211 * t25) - g(2) * (t121 * t130 + t20 * t211) - g(3) * (t121 * t131 + t211 * t37) t1, -g(1) * (t177 * t25 + t201) - g(2) * (t177 * t20 + t202) - g(3) * (t177 * t37 + t200); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169, t168, 0, 0, 0, 0, 0, 0, 0, 0, t169 * t121, -t169 * t118, -t168, -g(1) * (-pkin(5) * t8 + pkin(13) * t9) - g(2) * (pkin(5) * t234 + pkin(13) * t5) - g(3) * (pkin(5) * t14 + pkin(13) * t15); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 + g(2) * t238 - g(3) * (-t118 * t15 + t121 * t37) g(1) * t3 + g(2) * t237 - g(3) * (-t118 * t37 - t121 * t15) 0, 0;];
taug_reg  = t4;
