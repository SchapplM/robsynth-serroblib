% Calculate inertial parameters regressor of gravitation load for
% S6RRRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 06:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_gravloadJ_reg2_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 18:43:28
% EndTime: 2019-05-08 18:44:07
% DurationCPUTime: 4.84s
% Computational Cost: add. (3051->274), mult. (8823->456), div. (0->0), fcn. (11549->18), ass. (0->160)
t122 = sin(qJ(6));
t123 = sin(qJ(5));
t127 = cos(qJ(5));
t124 = sin(qJ(4));
t231 = cos(qJ(4));
t120 = sin(pkin(8));
t223 = cos(pkin(8));
t129 = cos(qJ(2));
t130 = cos(qJ(1));
t225 = cos(pkin(6));
t229 = sin(qJ(2));
t187 = t225 * t229;
t230 = sin(qJ(1));
t107 = t230 * t129 + t130 * t187;
t125 = sin(qJ(3));
t128 = cos(qJ(3));
t202 = t129 * t225;
t106 = -t130 * t202 + t230 * t229;
t121 = sin(pkin(7));
t222 = sin(pkin(6));
t200 = t130 * t222;
t224 = cos(pkin(7));
t157 = t224 * t106 + t121 * t200;
t234 = t107 * t125 + t157 * t128;
t182 = t224 * t222;
t239 = -t106 * t121 + t130 * t182;
t236 = t120 * t239 + t223 * t234;
t73 = t107 * t128 - t125 * t157;
t25 = -t236 * t124 + t73 * t231;
t52 = -t234 * t120 + t223 * t239;
t7 = t123 * t52 - t127 * t25;
t247 = t122 * t7;
t126 = cos(qJ(6));
t246 = t126 * t7;
t245 = -t123 * t25 - t52 * t127;
t242 = t124 * t73;
t108 = t129 * t130 - t230 * t187;
t155 = t130 * t229 + t230 * t202;
t185 = t222 * t230;
t141 = t121 * t185 - t155 * t224;
t136 = t108 * t125 - t141 * t128;
t240 = t155 * t121 + t230 * t182;
t241 = t136 * t120 + t223 * t240;
t238 = -t120 * t240 + t136 * t223;
t154 = t225 * t121 + t129 * t182;
t184 = t222 * t229;
t143 = t125 * t184 - t154 * t128;
t201 = t129 * t222;
t153 = -t121 * t201 + t225 * t224;
t237 = -t153 * t120 + t143 * t223;
t206 = t121 * t223;
t203 = t128 * t224;
t81 = t106 * t125 - t107 * t203;
t61 = t107 * t206 - t81 * t120;
t83 = -t108 * t203 + t155 * t125;
t62 = t108 * t206 - t83 * t120;
t235 = -g(1) * t62 - g(2) * t61;
t228 = pkin(11) * t121;
t227 = pkin(12) * t120;
t166 = t229 * t182;
t100 = -t125 * t201 - t128 * t166;
t221 = t100 * t120;
t218 = t120 * t121;
t217 = t120 * t123;
t216 = t120 * t127;
t215 = t122 * t127;
t214 = t126 * t127;
t179 = t121 * t184;
t213 = pkin(2) * t201 + pkin(11) * t179;
t212 = t130 * pkin(1) + pkin(10) * t185;
t24 = t236 * t231 + t242;
t211 = -t24 * pkin(4) + pkin(13) * t25;
t76 = t108 * t128 + t141 * t125;
t28 = t124 * t76 + t238 * t231;
t29 = -t238 * t124 + t76 * t231;
t210 = -t28 * pkin(4) + pkin(13) * t29;
t95 = t154 * t125 + t128 * t184;
t47 = t124 * t95 + t237 * t231;
t48 = -t237 * t124 + t95 * t231;
t209 = -t47 * pkin(4) + pkin(13) * t48;
t208 = t120 * t231;
t205 = t124 * t223;
t204 = t125 * t224;
t101 = -t125 * t166 + t128 * t201;
t161 = t223 * t179;
t198 = t101 * pkin(3) + pkin(12) * t161 + t213;
t197 = -pkin(3) * t234 + t73 * t227;
t196 = -t136 * pkin(3) + t76 * t227;
t195 = -t143 * pkin(3) + t95 * t227;
t8 = t123 * t29 - t127 * t241;
t194 = g(1) * t245 + g(2) * t8;
t193 = t121 * t208;
t192 = -t106 * pkin(2) + t107 * t228;
t191 = -t155 * pkin(2) + t108 * t228;
t186 = t223 * t231;
t26 = -t186 * t234 - t208 * t239 - t242;
t190 = g(1) * t26 + g(2) * t28;
t188 = -t230 * pkin(1) + pkin(10) * t200;
t183 = -pkin(5) * t127 - pkin(14) * t123;
t82 = -t106 * t128 - t107 * t204;
t181 = t82 * pkin(3) + t192;
t84 = -t108 * t204 - t155 * t128;
t180 = t84 * pkin(3) + t191;
t68 = t143 * t120 + t153 * t223;
t18 = -t123 * t48 + t127 * t68;
t177 = g(1) * t8 - g(2) * t245 - g(3) * t18;
t19 = t123 * t68 + t127 * t48;
t9 = t123 * t241 + t29 * t127;
t176 = g(1) * t9 - g(2) * t7 + g(3) * t19;
t35 = t82 * t231 + (t107 * t218 + t223 * t81) * t124;
t10 = t123 * t35 - t61 * t127;
t37 = t84 * t231 + (t108 * t218 + t223 * t83) * t124;
t12 = t123 * t37 - t62 * t127;
t168 = t120 * t179;
t60 = t101 * t231 + (t223 * t100 + t168) * t124;
t86 = t161 - t221;
t42 = t123 * t60 - t86 * t127;
t175 = g(1) * t12 + g(2) * t10 + g(3) * t42;
t39 = -t73 * t205 - t231 * t234;
t14 = t123 * t39 - t73 * t216;
t41 = -t136 * t231 - t76 * t205;
t16 = t123 * t41 - t76 * t216;
t55 = -t143 * t231 - t95 * t205;
t44 = t123 * t55 - t95 * t216;
t174 = g(1) * t16 + g(2) * t14 + g(3) * t44;
t173 = g(1) * t28 + g(2) * t24 + g(3) * t47;
t172 = g(1) * t29 + g(2) * t25 + g(3) * t48;
t34 = -t107 * t193 + t82 * t124 - t81 * t186;
t36 = -t108 * t193 + t84 * t124 - t83 * t186;
t59 = -t100 * t186 + t101 * t124 - t231 * t168;
t171 = g(1) * t36 + g(2) * t34 + g(3) * t59;
t38 = -t124 * t234 + t73 * t186;
t40 = -t136 * t124 + t76 * t186;
t54 = -t143 * t124 + t95 * t186;
t170 = g(1) * t40 + g(2) * t38 + g(3) * t54;
t169 = g(1) * t76 + g(2) * t73 + g(3) * t95;
t167 = t60 * pkin(4) + pkin(13) * t59 + t198;
t164 = t39 * pkin(4) + pkin(13) * t38 + t197;
t163 = t41 * pkin(4) + pkin(13) * t40 + t196;
t162 = t55 * pkin(4) + pkin(13) * t54 + t195;
t160 = t35 * pkin(4) + t34 * pkin(13) + t181;
t159 = t37 * pkin(4) + t36 * pkin(13) + t180;
t156 = -t107 * pkin(2) + pkin(11) * t239 + t188;
t152 = -g(1) * t108 - g(2) * t107 - g(3) * t184;
t147 = -t73 * pkin(3) + pkin(12) * t52 + t156;
t145 = t108 * pkin(2) + pkin(11) * t240 + t212;
t144 = -pkin(4) * t25 + t26 * pkin(13) + t147;
t137 = (g(3) * t221 + t235) * pkin(12);
t133 = t76 * pkin(3) + pkin(12) * t241 + t145;
t132 = t29 * pkin(4) + t28 * pkin(13) + t133;
t45 = t127 * t55 + t95 * t217;
t43 = t123 * t86 + t127 * t60;
t17 = t127 * t41 + t76 * t217;
t15 = t127 * t39 + t73 * t217;
t13 = t123 * t62 + t127 * t37;
t11 = t123 * t61 + t127 * t35;
t3 = t122 * t28 + t126 * t9;
t2 = -t122 * t9 + t126 * t28;
t1 = t173 * t123;
t4 = [0, 0, 0, 0, 0, 0, g(1) * t230 - g(2) * t130, g(1) * t130 + g(2) * t230, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t107 - g(2) * t108, -g(1) * t106 + g(2) * t155, -g(1) * t200 - g(2) * t185, -g(1) * t188 - g(2) * t212, 0, 0, 0, 0, 0, 0, g(1) * t73 - g(2) * t76, -g(1) * t234 + g(2) * t136, -g(1) * t239 - g(2) * t240, -g(1) * t156 - g(2) * t145, 0, 0, 0, 0, 0, 0, g(1) * t25 - g(2) * t29, t190, -g(1) * t52 - g(2) * t241, -g(1) * t147 - g(2) * t133, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9, t194, -t190, -g(1) * t144 - g(2) * t132, 0, 0, 0, 0, 0, 0, -g(1) * (t122 * t26 + t246) - g(2) * t3, -g(1) * (t126 * t26 - t247) - g(2) * t2, -t194, -g(1) * (t7 * pkin(5) + pkin(14) * t245 + t144) - g(2) * (t9 * pkin(5) + t8 * pkin(14) + t132); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t155 + g(2) * t106 - g(3) * t201, -t152, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t84 - g(2) * t82 - g(3) * t101, -g(1) * t83 - g(2) * t81 - g(3) * t100, t152 * t121, -g(1) * t191 - g(2) * t192 - g(3) * t213, 0, 0, 0, 0, 0, 0, -g(1) * t37 - g(2) * t35 - g(3) * t60, t171, -g(3) * t86 + t235, -g(1) * t180 - g(2) * t181 - g(3) * t198 + t137, 0, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t11 - g(3) * t43, t175, -t171, -g(1) * t159 - g(2) * t160 - g(3) * t167 + t137, 0, 0, 0, 0, 0, 0, -g(1) * (t122 * t36 + t126 * t13) - g(2) * (t11 * t126 + t122 * t34) - g(3) * (t122 * t59 + t126 * t43) -g(1) * (-t122 * t13 + t126 * t36) - g(2) * (-t11 * t122 + t126 * t34) - g(3) * (-t122 * t43 + t126 * t59) -t175, -g(1) * (t13 * pkin(5) + t12 * pkin(14) + t159) - g(2) * (t11 * pkin(5) + t10 * pkin(14) + t160) - g(3) * (pkin(5) * t43 + pkin(14) * t42 + t167) + t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152 * t125 + (-g(1) * t141 + g(2) * t157 - g(3) * t154) * t128, t169, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t41 - g(2) * t39 - g(3) * t55, t170, -t169 * t120, -g(1) * t196 - g(2) * t197 - g(3) * t195, 0, 0, 0, 0, 0, 0, -g(1) * t17 - g(2) * t15 - g(3) * t45, t174, -t170, -g(1) * t163 - g(2) * t164 - g(3) * t162, 0, 0, 0, 0, 0, 0, -g(1) * (t122 * t40 + t126 * t17) - g(2) * (t122 * t38 + t126 * t15) - g(3) * (t122 * t54 + t126 * t45) -g(1) * (-t122 * t17 + t126 * t40) - g(2) * (-t122 * t15 + t126 * t38) - g(3) * (-t122 * t45 + t126 * t54) -t174, -g(1) * (pkin(5) * t17 + pkin(14) * t16 + t163) - g(2) * (pkin(5) * t15 + pkin(14) * t14 + t164) - g(3) * (pkin(5) * t45 + pkin(14) * t44 + t162); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t172, 0, 0, 0, 0, 0, 0, 0, 0, t173 * t127, -t1, -t172, -g(1) * t210 - g(2) * t211 - g(3) * t209, 0, 0, 0, 0, 0, 0, -g(1) * (t122 * t29 - t214 * t28) - g(2) * (t122 * t25 - t214 * t24) - g(3) * (t122 * t48 - t214 * t47) -g(1) * (t126 * t29 + t215 * t28) - g(2) * (t126 * t25 + t215 * t24) - g(3) * (t126 * t48 + t215 * t47) t1, -g(1) * (t183 * t28 + t210) - g(2) * (t183 * t24 + t211) - g(3) * (t183 * t47 + t209); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, t176, 0, 0, 0, 0, 0, 0, 0, 0, t177 * t126, -t177 * t122, -t176, -g(1) * (-pkin(5) * t8 + pkin(14) * t9) - g(2) * (pkin(5) * t245 - pkin(14) * t7) - g(3) * (pkin(5) * t18 + pkin(14) * t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * (t126 * t24 + t247) - g(3) * (-t122 * t19 + t126 * t47) g(1) * t3 - g(2) * (-t122 * t24 + t246) - g(3) * (-t122 * t47 - t126 * t19) 0, 0;];
taug_reg  = t4;
