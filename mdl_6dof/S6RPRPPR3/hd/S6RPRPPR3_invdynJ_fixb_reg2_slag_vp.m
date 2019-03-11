% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRPPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:45:31
% EndTime: 2019-03-09 02:45:37
% DurationCPUTime: 3.56s
% Computational Cost: add. (3292->446), mult. (6141->513), div. (0->0), fcn. (3366->10), ass. (0->237)
t158 = sin(qJ(3));
t161 = cos(qJ(3));
t153 = sin(pkin(9));
t104 = pkin(1) * t153 + pkin(7);
t79 = t104 * qJD(1);
t223 = -t161 * qJD(2) + t158 * t79;
t311 = qJD(4) + t223;
t38 = -qJD(3) * pkin(3) + t311;
t154 = cos(pkin(9));
t105 = -t154 * pkin(1) - pkin(2);
t80 = qJD(1) * t105;
t260 = qJD(1) * t161;
t316 = qJD(2) * qJD(3) + qJDD(1) * t104;
t128 = t161 * qJDD(1);
t157 = sin(qJ(6));
t160 = cos(qJ(6));
t252 = t160 * qJD(3);
t254 = qJD(6) * t161;
t234 = t157 * t254;
t312 = t158 * t252 + t234;
t25 = qJD(1) * t312 - qJD(6) * t252 - t157 * qJDD(3) - t160 * t128;
t69 = t157 * t260 - t252;
t261 = qJD(1) * t158;
t97 = qJD(6) + t261;
t295 = t69 * t97;
t309 = t25 - t295;
t249 = qJD(1) * qJD(3);
t230 = t158 * t249;
t259 = qJD(3) * t157;
t70 = t160 * t260 + t259;
t277 = qJD(6) * t70;
t26 = -t160 * qJDD(3) + t277 + (t128 - t230) * t157;
t294 = t70 * t97;
t184 = t26 - t294;
t302 = pkin(3) + pkin(4);
t315 = t158 * t302;
t314 = qJDD(3) * qJ(4) + qJD(3) * qJD(4);
t229 = t161 * t249;
t243 = t158 * qJDD(1);
t313 = t229 + t243;
t236 = t302 * qJD(3);
t231 = t302 * qJDD(3);
t149 = qJD(3) * qJ(4);
t253 = t158 * qJD(2);
t36 = -(qJ(5) * qJD(1) - t79) * t161 + t253;
t310 = t149 + t36;
t54 = t161 * t79 + t253;
t42 = t149 + t54;
t281 = pkin(1) * qJDD(1);
t268 = -qJ(5) + t104;
t145 = qJ(1) + pkin(9);
t122 = sin(t145);
t101 = g(2) * t122;
t123 = cos(t145);
t308 = g(1) * t123 + t101;
t307 = 0.2e1 * t314;
t139 = t158 * pkin(5);
t212 = t161 * pkin(8) + t139;
t258 = qJD(3) * t158;
t66 = qJDD(6) + t313;
t285 = t161 * t66;
t306 = -t160 * (t97 * t258 - t285) - t97 * t234;
t220 = -t158 * qJDD(2) - t161 * t316 + t79 * t258;
t18 = -t220 + t314;
t248 = qJD(1) * qJD(5);
t250 = qJ(5) * qJDD(1);
t94 = qJ(5) * t230;
t12 = t161 * (t248 + t250) - t18 - t94;
t10 = qJDD(3) * pkin(5) - t12;
t251 = pkin(8) + t302;
t297 = g(3) * t158;
t305 = -qJD(6) * t251 * t97 + t161 * t308 - t10 + t297;
t244 = qJDD(3) * t104;
t43 = -pkin(3) * t260 - qJ(4) * t261 + t80;
t134 = t158 * qJ(4);
t142 = t161 * pkin(3);
t262 = t142 + t134;
t58 = -t262 + t105;
t304 = qJD(3) * (qJD(1) * t58 + t43) - t244;
t186 = pkin(5) * t161 - t251 * t158;
t178 = t186 * qJD(3);
t130 = t158 * qJD(4);
t108 = t154 * t281;
t78 = -qJDD(1) * pkin(2) - t108;
t194 = pkin(3) * t128 + qJ(4) * t313 + qJD(1) * t130 - t78;
t179 = pkin(4) * t128 + qJDD(5) + t194;
t6 = qJD(1) * t178 + t212 * qJDD(1) + t179;
t32 = pkin(4) * t260 + qJD(5) - t43;
t22 = t212 * qJD(1) + t32;
t112 = qJ(5) * t261;
t196 = -t112 + t311;
t27 = -t251 * qJD(3) + t196;
t7 = -t157 * t27 + t160 * t22;
t257 = qJD(3) * t161;
t219 = -t161 * qJDD(2) + t158 * t316 + t79 * t257;
t204 = -qJDD(4) - t219;
t167 = -qJ(5) * t313 - t158 * t248 - t204;
t9 = -t251 * qJDD(3) + t167;
t1 = qJD(6) * t7 + t157 * t6 + t160 * t9;
t5 = t160 * t6;
t8 = t157 * t22 + t160 * t27;
t2 = -qJD(6) * t8 - t157 * t9 + t5;
t208 = t157 * t8 + t160 * t7;
t303 = -qJD(6) * t208 + t1 * t160 - t2 * t157;
t301 = t7 * t97;
t300 = t8 * t97;
t102 = g(1) * t122;
t299 = g(2) * qJ(5);
t298 = g(2) * t123;
t144 = g(3) * t161;
t141 = t161 * pkin(4);
t296 = t69 * t70;
t293 = t25 * t158 - t70 * t257;
t292 = t312 * t69;
t233 = t160 * t254;
t291 = t157 * t285 + t97 * t233;
t289 = t157 * t66;
t288 = t158 * t26;
t287 = t160 * t66;
t286 = t161 * t310;
t284 = t161 * t69;
t283 = t161 * t70;
t282 = t25 * t157;
t280 = qJ(4) * t161;
t30 = qJD(3) * pkin(5) + t310;
t279 = qJD(3) * t30;
t278 = qJD(3) * t310;
t276 = qJDD(3) * pkin(3);
t275 = t122 * t158;
t274 = t122 * t161;
t273 = t123 * t158;
t272 = t123 * t161;
t151 = t161 ^ 2;
t165 = qJD(1) ^ 2;
t271 = t151 * t165;
t270 = t157 * t158;
t269 = t158 * t160;
t35 = -t112 + t223;
t267 = qJD(4) + t35;
t265 = -qJD(5) - t32;
t263 = qJ(4) * t257 + t130;
t256 = qJD(6) * t157;
t255 = qJD(6) * t160;
t150 = t158 ^ 2;
t246 = qJDD(1) * t150;
t245 = qJDD(1) * t151;
t242 = g(1) * t273 + g(2) * t275 - t144;
t241 = t97 * t270;
t240 = t97 * t269;
t162 = cos(qJ(1));
t239 = t162 * pkin(1) + t123 * pkin(2) + t122 * pkin(7);
t238 = t70 * t258;
t237 = t141 + t262;
t159 = sin(qJ(1));
t232 = -pkin(1) * t159 + t123 * pkin(7);
t228 = -pkin(2) - t134;
t227 = t102 - t298;
t226 = g(1) * t251;
t52 = t141 - t58;
t225 = qJD(1) * t52 + t32;
t222 = t265 * t158;
t218 = -t308 + (t245 + t246) * t104;
t217 = t157 * t238;
t215 = t158 * t236;
t214 = t158 * t229;
t213 = pkin(3) * t272 + t123 * t134 + t239;
t210 = g(1) * t159 - g(2) * t162;
t207 = -t157 * t7 + t160 * t8;
t164 = qJD(3) ^ 2;
t205 = t104 * t164 + t298;
t203 = t237 + t212;
t34 = t203 - t105;
t59 = t268 * t158;
t16 = -t157 * t59 + t160 * t34;
t17 = t157 * t34 + t160 * t59;
t200 = t157 * t97;
t199 = pkin(4) * t272 + t213;
t198 = t228 - t142;
t197 = -qJD(6) * t22 - t144 - t9;
t193 = -t241 - t284;
t192 = t241 - t284;
t190 = -t97 * t255 - t289;
t189 = -t97 * t256 + t287;
t188 = -qJD(3) * t223 + t220;
t187 = -t219 + t242;
t181 = g(1) * (-qJ(5) * t123 + t232);
t180 = -qJDD(4) + t187;
t177 = t251 * t66 - t97 * t30;
t15 = -qJD(1) * t215 + t179;
t39 = -t215 + t263;
t176 = -qJD(1) * t39 - qJDD(1) * t52 - t15 + t298;
t175 = qJDD(1) * t105 + t205 + t78;
t174 = 0.2e1 * t80 * qJD(3) - t244;
t173 = qJD(3) * t54 + t187;
t20 = pkin(3) * t230 - t194;
t57 = pkin(3) * t258 - t263;
t172 = -qJD(1) * t57 - qJDD(1) * t58 - t20 - t205;
t171 = -qJ(5) * t243 - t180;
t19 = -t204 - t276;
t170 = t19 * t158 + t18 * t161 + (-t158 * t42 + t161 * t38) * qJD(3);
t169 = t219 * t158 - t220 * t161 + (-t158 * t54 + t161 * t223) * qJD(3);
t168 = qJD(6) * t207 + t1 * t157 + t160 * t2 - t298;
t156 = qJ(4) + pkin(5);
t129 = t150 * t165;
t115 = qJ(4) * t260;
t96 = t158 * t165 * t161;
t92 = -t129 - t164;
t89 = g(1) * t274;
t88 = g(1) * t275;
t84 = qJ(4) * t272;
t82 = qJ(4) * t274;
t81 = qJDD(3) + t96;
t76 = -t129 + t271;
t75 = qJDD(3) * t161 - t164 * t158;
t74 = qJDD(3) * t158 + t161 * t164;
t71 = pkin(3) * t261 - t115;
t62 = -0.2e1 * t214 + t245;
t61 = 0.2e1 * t214 + t246;
t60 = t268 * t161;
t55 = -t302 * t261 + t115;
t50 = -t122 * t157 + t123 * t269;
t49 = -t122 * t160 - t123 * t270;
t48 = -t122 * t269 - t123 * t157;
t47 = t122 * t270 - t123 * t160;
t44 = t158 * t128 + (-t150 + t151) * t249;
t41 = -t158 * qJD(5) + t268 * t257;
t40 = qJD(5) * t161 + t258 * t268;
t37 = 0.2e1 * t44;
t31 = qJD(1) * t186 + t115;
t29 = -t236 + t196;
t28 = t178 + t263;
t14 = t157 * t31 + t160 * t36;
t13 = -t157 * t36 + t160 * t31;
t11 = -t231 + t167;
t4 = -qJD(6) * t17 - t157 * t41 + t160 * t28;
t3 = qJD(6) * t16 + t157 * t28 + t160 * t41;
t21 = [0, 0, 0, 0, 0, qJDD(1), t210, g(1) * t162 + g(2) * t159, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t108 + t227, -0.2e1 * t153 * t281 + t308, 0 (t210 + (t153 ^ 2 + t154 ^ 2) * t281) * pkin(1), t61, t37, t74, t62, t75, 0, t158 * t174 - t161 * t175 + t89, t158 * t175 + t161 * t174 - t88, t169 + t218, t78 * t105 - g(1) * (-pkin(2) * t122 + t232) - g(2) * t239 + t169 * t104, t61, t74, -0.2e1 * t44, 0, -t75, t62, t158 * t304 + t172 * t161 + t89, t170 + t218, t172 * t158 - t161 * t304 + t88, -g(1) * t232 - g(2) * t213 - t198 * t102 + t104 * t170 + t20 * t58 + t43 * t57, t62, t37, t75, t61, t74, 0, qJDD(3) * t60 + t88 + (t161 * t225 - t40) * qJD(3) - t176 * t158, qJDD(3) * t59 - t89 + (t158 * t225 + t41) * qJD(3) + t176 * t161 (-qJD(3) * t29 - qJDD(1) * t60 + t12 + (-qJD(3) * t59 + t40) * qJD(1)) * t161 + (t278 - qJDD(1) * t59 - t11 + (qJD(3) * t60 - t41) * qJD(1)) * t158 + t308, t11 * t59 + t29 * t41 - t12 * t60 - t310 * t40 + t15 * t52 + t32 * t39 - t181 - g(2) * t199 + (-g(1) * (t198 - t141) + t299) * t122, -t70 * t234 + (-t161 * t25 - t238) * t160, t217 + (t282 + (-t26 - t277) * t160) * t161 + t292, t293 - t306, t69 * t233 + (t161 * t26 - t69 * t258) * t157, -qJD(3) * t192 + t288 + t291, t158 * t66 + t97 * t257, -g(1) * t48 - g(2) * t50 + t16 * t66 - t26 * t60 + t4 * t97 + t40 * t69 + (t30 * t259 + t2) * t158 + (qJD(3) * t7 - t10 * t157 - t30 * t255) * t161, -g(1) * t47 - g(2) * t49 - t17 * t66 + t25 * t60 - t3 * t97 + t40 * t70 + (t252 * t30 - t1) * t158 + (-qJD(3) * t8 - t10 * t160 + t256 * t30) * t161, -t16 * t25 + t161 * t168 + t17 * t26 - t208 * t258 + t3 * t69 + t4 * t70 + t89, t1 * t17 + t8 * t3 + t2 * t16 + t7 * t4 + t10 * t60 - t30 * t40 - t181 - g(2) * (pkin(5) * t273 + pkin(8) * t272 + t199) + (-g(1) * (t228 - t139) + t299 + t161 * t226) * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, t75, -t74, 0, -t158 * t220 - t161 * t219 - g(3) + (t158 * t223 + t161 * t54) * qJD(3), 0, 0, 0, 0, 0, 0, t75, 0, t74, t158 * t18 - t161 * t19 - g(3) + (t158 * t38 + t161 * t42) * qJD(3), 0, 0, 0, 0, 0, 0, t74, -t75, 0, -t11 * t161 - t12 * t158 - g(3) + (t158 * t29 + t286) * qJD(3), 0, 0, 0, 0, 0, 0, qJD(3) * t193 - t288 + t291, t293 + t306, -t217 + (-t282 + (-t26 + t277) * t160) * t161 + t292, -g(3) + (qJD(3) * t207 + t10) * t158 + (t279 - t303) * t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t76, t243, t96, t128, qJDD(3), -t80 * t261 + t173, t297 + (-qJD(1) * t80 + t308) * t161 + t188, 0, 0, -t96, t243, t76, qJDD(3), -t128, t96, 0.2e1 * t276 - qJDD(4) + (-t158 * t43 + t161 * t71) * qJD(1) + t173 (-pkin(3) * t158 + t280) * qJDD(1) (qJD(1) * t71 - g(3)) * t158 + (qJD(1) * t43 - t308) * t161 - t188 + t307, t18 * qJ(4) - t19 * pkin(3) - t43 * t71 - t38 * t54 - g(1) * (-pkin(3) * t273 + t84) - g(2) * (-pkin(3) * t275 + t82) - g(3) * t262 + t311 * t42, t96, -t76, t128, -t96, t243, qJDD(3), qJD(3) * t35 + t94 + (-qJD(1) * t55 - g(3)) * t158 + (t265 * qJD(1) - t250 - t308) * t161 - t220 + t307, -qJD(3) * t36 - 0.2e1 * t231 + ((-qJ(5) * qJD(3) + t55) * t161 + t222) * qJD(1) + t171 (-t280 + t315) * qJDD(1) + (-t267 + t29 + t236) * t260, -g(1) * t84 - g(2) * t82 - g(3) * t237 - t12 * qJ(4) - t11 * t302 + t267 * t310 - t29 * t36 + t308 * t315 - t32 * t55, t160 * t294 - t282 (-t25 - t295) * t160 + (-t26 - t294) * t157 (-t240 + t283) * qJD(1) + t190, -t26 * t160 + t200 * t69, qJD(1) * t192 - t189, -t97 * t260, -t13 * t97 - t156 * t26 + t177 * t157 - t160 * t305 - t7 * t260 - t267 * t69, t14 * t97 + t156 * t25 + t157 * t305 + t177 * t160 + t8 * t260 - t267 * t70, -t13 * t70 - t14 * t69 + (t7 * t261 - t251 * t26 - t1 + (t251 * t70 + t7) * qJD(6)) * t160 + (t8 * t261 - t251 * t25 + t2 + (t251 * t69 + t8) * qJD(6)) * t157 + t242, t10 * t156 - t8 * t14 - t7 * t13 - g(1) * (pkin(5) * t272 + t84) - g(2) * (pkin(5) * t274 + t82) - g(3) * t203 + t267 * t30 + (t251 * t101 + t123 * t226) * t158 - t303 * t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, t243, t92, -qJD(3) * t42 + t43 * t261 - t180 - t276, 0, 0, 0, 0, 0, 0, t92, t81, -t243, -t278 - t231 + (-qJ(5) * t257 + t222) * qJD(1) + t171, 0, 0, 0, 0, 0, 0, -t160 * t97 ^ 2 + qJD(3) * t69 - t289, qJD(3) * t70 + t200 * t97 - t287, t157 * t309 + t184 * t160, -t279 + (t1 - t301) * t160 + (-t2 - t300) * t157 - t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t229 + t243, -t128 + 0.2e1 * t230, -t129 - t271 (t286 + (t29 - t236) * t158) * qJD(1) + t179 + t227, 0, 0, 0, 0, 0, 0, qJD(1) * t193 + t189 (-t240 - t283) * qJD(1) + t190, t184 * t157 - t160 * t309, t102 + (t158 * t207 + t161 * t30) * qJD(1) + t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t296, -t69 ^ 2 + t70 ^ 2, t309, -t296, t184, t66, -g(1) * t49 + g(2) * t47 + t157 * t197 - t255 * t27 + t30 * t70 + t300 + t5, g(1) * t50 - g(2) * t48 - t30 * t69 + t301 + (qJD(6) * t27 - t6) * t157 + t197 * t160, 0, 0;];
tau_reg  = t21;
