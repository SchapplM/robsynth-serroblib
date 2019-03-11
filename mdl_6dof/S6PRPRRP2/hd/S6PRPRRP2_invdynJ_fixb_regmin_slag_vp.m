% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tau_reg [6x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:03:09
% EndTime: 2019-03-08 20:03:18
% DurationCPUTime: 3.48s
% Computational Cost: add. (3724->413), mult. (8592->547), div. (0->0), fcn. (7161->12), ass. (0->214)
t163 = sin(qJ(4));
t166 = cos(qJ(4));
t204 = pkin(4) * t163 - pkin(9) * t166;
t122 = t204 * qJD(4);
t159 = cos(pkin(11));
t164 = sin(qJ(2));
t158 = sin(pkin(6));
t241 = qJD(1) * t158;
t222 = t164 * t241;
t128 = t159 * t222;
t156 = sin(pkin(11));
t167 = cos(qJ(2));
t221 = t167 * t241;
t90 = t156 * t221 + t128;
t301 = t122 - t90;
t240 = qJD(2) * t158;
t215 = qJD(1) * t240;
t227 = qJDD(1) * t158;
t300 = t164 * t227 + t167 * t215;
t238 = qJD(2) * t166;
t299 = qJD(5) - t238;
t231 = qJD(5) * t163;
t298 = qJD(2) * t231 - qJDD(4);
t162 = sin(qJ(5));
t165 = cos(qJ(5));
t226 = t163 * qJDD(2);
t73 = t162 * ((qJD(5) + t238) * qJD(4) + t226) + t298 * t165;
t252 = t162 * t166;
t127 = t156 * t222;
t93 = t159 * t221 - t127;
t297 = -t165 * t301 - t93 * t252;
t196 = pkin(4) * t166 + pkin(9) * t163 + pkin(3);
t284 = pkin(2) * t159;
t109 = -t196 - t284;
t230 = qJD(5) * t165;
t249 = t165 * t166;
t296 = t109 * t230 + t162 * t301 - t93 * t249;
t161 = cos(pkin(6));
t141 = qJD(1) * t161 + qJD(3);
t124 = qJD(2) * pkin(2) + t221;
t85 = t124 * t156 + t128;
t83 = qJD(2) * pkin(8) + t85;
t295 = t141 * t166 - t163 * t83;
t146 = pkin(2) * t156 + pkin(8);
t243 = t109 * t162 + t146 * t249;
t157 = sin(pkin(10));
t160 = cos(pkin(10));
t253 = t161 * t167;
t294 = -t157 * t253 - t160 * t164;
t197 = t156 * t167 + t159 * t164;
t248 = t167 * t159;
t112 = t156 * t164 - t248;
t98 = t112 * t161;
t67 = -t157 * t197 - t160 * t98;
t70 = t157 * t98 - t160 * t197;
t258 = t158 * t164;
t96 = t156 * t258 - t158 * t248;
t185 = -g(1) * t70 - g(2) * t67 + g(3) * t96;
t97 = t197 * t158;
t198 = t161 * t166 - t163 * t97;
t257 = t158 * t166;
t99 = t197 * t161;
t66 = t112 * t157 - t160 * t99;
t69 = t112 * t160 + t157 * t99;
t187 = g(3) * t198 + g(2) * (-t160 * t257 + t163 * t66) + g(1) * (t157 * t257 + t163 * t69);
t229 = t165 * qJD(4);
t239 = qJD(2) * t163;
t114 = t162 * t239 - t229;
t236 = qJD(4) * t162;
t116 = t165 * t239 + t236;
t54 = -qJD(4) * pkin(4) - t295;
t25 = pkin(5) * t114 - qJ(6) * t116 + t54;
t152 = t166 * qJDD(2);
t228 = qJD(2) * qJD(4);
t111 = t163 * t228 + qJDD(5) - t152;
t282 = pkin(9) * t111;
t293 = -t25 * t299 + t282;
t291 = t116 ^ 2;
t283 = pkin(5) * t111;
t232 = qJD(5) * t162;
t235 = qJD(4) * t163;
t261 = t146 * t165;
t281 = (-t146 * t232 - qJD(6)) * t166 + (qJ(6) - t261) * t235 + t296;
t262 = t146 * t162;
t210 = pkin(5) + t262;
t280 = qJD(5) * t243 - t210 * t235 + t297;
t220 = t166 * t229;
t251 = t163 * t165;
t279 = -t114 * t220 - t251 * t73;
t278 = pkin(9) * qJD(5);
t59 = t141 * t163 + t166 * t83;
t55 = qJD(4) * pkin(9) + t59;
t84 = t124 * t159 - t127;
t65 = -qJD(2) * t196 - t84;
t17 = t162 * t65 + t165 * t55;
t11 = qJ(6) * t299 + t17;
t277 = t11 * t299;
t276 = t299 * t17;
t214 = t166 * t228;
t72 = -qJD(5) * t229 + (-t214 - t226) * t165 + t298 * t162;
t275 = t166 * t72;
t274 = t93 * t114;
t273 = t93 * t116;
t202 = pkin(5) * t162 - qJ(6) * t165;
t272 = -t162 * qJD(6) + t202 * t299 - t59;
t121 = t204 * qJD(2);
t271 = t121 * t162 + t165 * t295;
t270 = qJ(6) * t111;
t269 = t111 * t162;
t268 = t114 * t299;
t267 = t116 * t114;
t266 = t116 * t299;
t138 = qJDD(1) * t161 + qJDD(3);
t264 = t138 * t163;
t260 = t157 * t164;
t259 = t158 * t163;
t256 = t158 * t167;
t254 = t161 * t164;
t250 = t165 * t109;
t16 = -t162 * t55 + t165 * t65;
t247 = qJD(6) - t16;
t246 = qJDD(1) - g(3);
t100 = t111 * t251;
t245 = t220 * t299 + t100;
t154 = t163 ^ 2;
t242 = -t166 ^ 2 + t154;
t237 = qJD(4) * t114;
t234 = qJD(4) * t166;
t233 = qJD(5) * t114;
t137 = t167 * t227;
t94 = qJDD(2) * pkin(2) - t164 * t215 + t137;
t51 = -t156 * t300 + t159 * t94;
t24 = qJD(2) * t122 - qJDD(2) * t196 - t51;
t52 = t156 * t94 + t159 * t300;
t48 = qJDD(2) * pkin(8) + t52;
t8 = qJDD(4) * pkin(9) + qJD(4) * t295 + t166 * t48 + t264;
t225 = t162 * t24 + t165 * t8 + t230 * t65;
t223 = t160 * t253;
t219 = t116 * t234;
t218 = t299 * t232;
t217 = t162 * t231;
t216 = t163 * t230;
t211 = t162 * t8 - t165 * t24 + t230 * t55 + t232 * t65;
t209 = t116 * t235 + t275;
t208 = -t72 + t233;
t206 = t116 * t216;
t203 = pkin(5) * t165 + qJ(6) * t162;
t10 = -pkin(5) * t299 + t247;
t201 = t10 * t165 - t11 * t162;
t200 = t10 * t162 + t11 * t165;
t81 = t161 * t163 + t166 * t97;
t37 = t162 * t96 + t165 * t81;
t36 = t162 * t81 - t165 * t96;
t195 = -t138 * t166 + t141 * t235 + t163 * t48 + t234 * t83;
t194 = pkin(4) + t203;
t192 = t146 + t202;
t191 = -t55 * t232 + t225;
t190 = -t230 * t299 - t269;
t18 = t165 * t66 + t252 * t67;
t20 = t165 * t69 + t252 * t70;
t49 = -t165 * t97 - t252 * t96;
t189 = g(1) * t20 + g(2) * t18 + g(3) * t49;
t19 = -t162 * t66 + t249 * t67;
t21 = -t162 * t69 + t249 * t70;
t50 = t162 * t97 - t249 * t96;
t188 = -g(1) * t21 - g(2) * t19 - g(3) * t50;
t42 = -t160 * t259 - t166 * t66;
t44 = t157 * t259 - t166 * t69;
t186 = g(1) * t44 + g(2) * t42 + g(3) * t81;
t92 = t112 * t240;
t33 = qJD(4) * t81 - t92 * t163;
t34 = qJD(4) * t198 - t92 * t166;
t91 = qJD(2) * t97;
t4 = qJD(5) * t37 + t34 * t162 - t91 * t165;
t184 = -t111 * t36 + t114 * t33 - t198 * t73 - t299 * t4;
t183 = -g(3) * t161 + (-g(1) * t157 + g(2) * t160) * t158;
t182 = -t217 + t220;
t9 = -qJDD(4) * pkin(4) + t195;
t181 = t299 * t54 - t282;
t5 = -qJD(5) * t36 + t91 * t162 + t34 * t165;
t180 = t111 * t37 - t116 * t33 - t198 * t72 + t299 * t5;
t147 = -pkin(3) - t284;
t82 = -qJD(2) * pkin(3) - t84;
t179 = -qJDD(4) * t146 + (qJD(2) * t147 + t82 + t93) * qJD(4);
t12 = t162 * t42 + t165 * t67;
t14 = t162 * t44 + t165 * t70;
t178 = g(1) * t14 + g(2) * t12 + g(3) * t36 - t211;
t177 = -t278 * t299 - t187;
t3 = pkin(5) * t73 + qJ(6) * t72 - qJD(6) * t116 + t9;
t176 = t177 - t3;
t175 = -g(1) * t294 - g(3) * t256;
t1 = qJD(6) * t299 + t191 + t270;
t2 = qJDD(6) + t211 - t283;
t174 = qJD(5) * t201 + t1 * t165 + t2 * t162;
t13 = -t162 * t67 + t165 * t42;
t15 = -t162 * t70 + t165 * t44;
t173 = -g(1) * t15 - g(2) * t13 - g(3) * t37 + t191;
t172 = t116 * t25 + qJDD(6) - t178;
t168 = qJD(4) ^ 2;
t171 = -qJD(2) * t90 + t146 * t168 - t185 - t51 + (-pkin(3) + t147) * qJDD(2);
t170 = -t163 * t269 - t166 * t73 + t114 * t235 - (t162 * t234 + t216) * t299;
t169 = qJD(2) ^ 2;
t133 = pkin(2) * t223;
t132 = qJDD(4) * t166 - t163 * t168;
t131 = qJDD(4) * t163 + t166 * t168;
t87 = t192 * t163;
t79 = pkin(5) * t116 + qJ(6) * t114;
t75 = t166 * t210 - t250;
t74 = -qJ(6) * t166 + t243;
t45 = (qJD(5) * t203 - qJD(6) * t165) * t163 + t192 * t234;
t35 = -t72 + t268;
t29 = -pkin(5) * t239 - t121 * t165 + t162 * t295;
t28 = qJ(6) * t239 + t271;
t6 = [t246, 0 (qJDD(2) * t167 - t164 * t169) * t158 (-qJDD(2) * t164 - t167 * t169) * t158, t138 * t161 - t51 * t96 + t52 * t97 - t84 * t91 - t85 * t92 - g(3), 0, 0, 0, 0, 0, -t96 * t152 - qJD(4) * t33 + qJDD(4) * t198 + (-t166 * t91 + t235 * t96) * qJD(2), t96 * t226 - qJD(4) * t34 - qJDD(4) * t81 + (t163 * t91 + t234 * t96) * qJD(2), 0, 0, 0, 0, 0, t184, -t180, t184, -t114 * t5 + t116 * t4 - t36 * t72 - t37 * t73, t180, t1 * t37 + t10 * t4 + t11 * t5 - t198 * t3 + t2 * t36 + t25 * t33 - g(3); 0, qJDD(2), t137 - g(2) * (t223 - t260) + t175, -g(1) * (t157 * t254 - t160 * t167) - g(2) * (-t157 * t167 - t160 * t254) - t246 * t258, -g(2) * t133 + t84 * t90 - t85 * t93 + (g(2) * t260 + t52 * t156 + t51 * t159 + t175) * pkin(2), qJDD(2) * t154 + 0.2e1 * t163 * t214, 0.2e1 * t152 * t163 - 0.2e1 * t228 * t242, t131, t132, 0, t163 * t179 - t166 * t171, t163 * t171 + t166 * t179, t116 * t182 - t251 * t72, -t206 + (-t219 + (t72 + t233) * t163) * t162 + t279, -t217 * t299 + t209 + t245 (-t236 * t299 + t73) * t166 + (t190 - t237) * t163, -t111 * t166 + t235 * t299, t111 * t250 - (t109 * t232 + t297) * t299 + (t54 * t236 + (t190 + t237) * t146 + t211) * t166 + (t54 * t230 - t274 + t146 * t73 + t9 * t162 + (t262 * t299 + t16) * qJD(4)) * t163 + t188, -t243 * t111 - t296 * t299 + ((t146 * t299 - t55) * t232 + (t116 * t146 + t165 * t54) * qJD(4) + t225) * t166 + (-t54 * t232 - t273 - t146 * t72 + t9 * t165 + (t261 * t299 - t17) * qJD(4)) * t163 + t189, -t111 * t75 + t114 * t45 + t73 * t87 + (t236 * t25 + t2) * t166 - t280 * t299 + (-qJD(4) * t10 + t162 * t3 + t230 * t25 - t274) * t163 + t188, -t72 * t75 - t73 * t74 + t280 * t116 - t281 * t114 + t201 * t234 + (-qJD(5) * t200 - t1 * t162 + t165 * t2 + t185) * t163, t111 * t74 - t116 * t45 + t72 * t87 + (-t229 * t25 - t1) * t166 + t281 * t299 + (qJD(4) * t11 - t165 * t3 + t232 * t25 + t273) * t163 - t189, t1 * t74 + t3 * t87 + t2 * t75 - g(1) * (pkin(2) * t294 + pkin(5) * t21 - pkin(8) * t69 + qJ(6) * t20) - g(2) * (-pkin(2) * t260 + pkin(5) * t19 - pkin(8) * t66 + qJ(6) * t18 + t133) - g(3) * (pkin(2) * t256 + pkin(5) * t50 + pkin(8) * t97 + qJ(6) * t49) + (-t163 * t93 + t45) * t25 + t281 * t11 + t280 * t10 + t185 * t196; 0, 0, 0, 0, t183 + t138, 0, 0, 0, 0, 0, t132, -t131, 0, 0, 0, 0, 0, t170, -t182 * t299 - t100 + t209, t170, t206 + (t163 * t208 + t219) * t162 + t279, -t275 + (-qJD(4) * t116 - t218) * t163 + t245 (qJD(4) * t200 - t3) * t166 + (qJD(4) * t25 + t174) * t163 + t183; 0, 0, 0, 0, 0, -t163 * t169 * t166, t242 * t169, t226, t152, qJDD(4), qJD(4) * t59 - t239 * t82 - t187 - t195, -t264 + (-qJD(2) * t82 - t48) * t166 + t186, -t72 * t162 + t165 * t266 (-t72 - t268) * t165 + (-t266 - t73) * t162 (-t116 * t163 - t249 * t299) * qJD(2) - t190, -t218 + t111 * t165 + (t114 * t163 + t252 * t299) * qJD(2), -t299 * t239, -t16 * t239 - pkin(4) * t73 - t59 * t114 + (t295 * t299 + t181) * t162 + (-t9 - (t121 + t278) * t299 - t187) * t165, pkin(4) * t72 + t271 * t299 + t17 * t239 - t59 * t116 + t181 * t165 + (-t177 + t9) * t162, t10 * t239 + t272 * t114 - t162 * t293 + t176 * t165 - t194 * t73 + t29 * t299, t114 * t28 - t116 * t29 + (t1 + t299 * t10 + (qJD(5) * t116 - t73) * pkin(9)) * t165 + (pkin(9) * t208 + t2 - t277) * t162 - t186, -t11 * t239 - t272 * t116 + t176 * t162 + t165 * t293 - t194 * t72 - t28 * t299, -t10 * t29 - t11 * t28 + t272 * t25 + (t174 - t186) * pkin(9) + (-t3 - t187) * t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t267, -t114 ^ 2 + t291, t35, -t73 + t266, t111, -t116 * t54 + t178 + t276, t114 * t54 + t16 * t299 - t173, -t114 * t79 - t172 + t276 + 0.2e1 * t283, pkin(5) * t72 - qJ(6) * t73 + (t11 - t17) * t116 + (t10 - t247) * t114, 0.2e1 * t270 - t114 * t25 + t116 * t79 - (-0.2e1 * qJD(6) + t16) * t299 + t173, t1 * qJ(6) - t2 * pkin(5) - t25 * t79 - t10 * t17 - g(1) * (-pkin(5) * t14 + qJ(6) * t15) - g(2) * (-pkin(5) * t12 + qJ(6) * t13) - g(3) * (-pkin(5) * t36 + qJ(6) * t37) + t247 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111 + t267, t35, -t299 ^ 2 - t291, t172 - t277 - t283;];
tau_reg  = t6;
