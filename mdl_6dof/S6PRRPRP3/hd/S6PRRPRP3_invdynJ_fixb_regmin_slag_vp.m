% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:38:24
% EndTime: 2019-03-08 21:38:37
% DurationCPUTime: 5.73s
% Computational Cost: add. (5714->519), mult. (13188->702), div. (0->0), fcn. (10410->14), ass. (0->247)
t210 = sin(qJ(3));
t212 = cos(qJ(3));
t259 = pkin(3) * t210 - qJ(4) * t212;
t139 = qJD(3) * t259 - qJD(4) * t210;
t205 = sin(pkin(11));
t207 = cos(pkin(11));
t211 = sin(qJ(2));
t300 = qJD(3) * t210;
t293 = pkin(8) * t300;
t206 = sin(pkin(6));
t306 = qJD(1) * t206;
t213 = cos(qJ(2));
t314 = t212 * t213;
t333 = t207 * t139 + t205 * t293 - (-t205 * t314 + t207 * t211) * t306;
t320 = t205 * t211;
t366 = t205 * t139 - (t207 * t314 + t320) * t306;
t315 = t207 * t212;
t254 = pkin(4) * t210 - pkin(9) * t315;
t365 = -qJD(3) * t254 - t333;
t316 = t207 * t210;
t319 = t205 * t212;
t364 = (-pkin(8) * t316 - pkin(9) * t319) * qJD(3) + t366;
t301 = qJD(3) * t207;
t303 = qJD(2) * t210;
t157 = -t205 * t303 + t301;
t158 = qJD(3) * t205 + t207 * t303;
t209 = sin(qJ(5));
t345 = cos(qJ(5));
t82 = -t345 * t157 + t158 * t209;
t363 = t82 ^ 2;
t245 = -t209 * t157 - t345 * t158;
t350 = t245 ^ 2;
t302 = qJD(2) * t212;
t192 = -qJD(5) + t302;
t362 = t192 * t82;
t361 = t192 * t245;
t290 = t345 * t207;
t243 = -t209 * t205 + t290;
t282 = qJD(5) * t345;
t298 = qJD(5) * t209;
t356 = -t205 * t298 + t207 * t282;
t311 = -t243 * t302 + t356;
t162 = t345 * t205 + t209 * t207;
t146 = t162 * qJD(5);
t232 = t212 * t162;
t310 = -qJD(2) * t232 + t146;
t296 = qJD(2) * qJD(3);
t280 = t212 * t296;
t294 = t210 * qJDD(2);
t360 = t280 + t294;
t260 = pkin(3) * t212 + qJ(4) * t210;
t169 = -pkin(2) - t260;
t123 = pkin(8) * t315 + t205 * t169;
t321 = t205 * t210;
t103 = -pkin(9) * t321 + t123;
t156 = t207 * t169;
t89 = -pkin(9) * t316 + t156 + (-pkin(8) * t205 - pkin(4)) * t212;
t359 = -t103 * t298 - t365 * t209 + t89 * t282 + t364 * t345;
t340 = pkin(9) + qJ(4);
t172 = t340 * t205;
t173 = t340 * t207;
t244 = -t345 * t172 - t209 * t173;
t163 = t259 * qJD(2);
t165 = qJD(2) * pkin(8) + t211 * t306;
t331 = cos(pkin(6));
t275 = qJD(1) * t331;
t233 = t210 * t165 - t212 * t275;
t71 = t207 * t163 + t205 * t233;
t49 = qJD(2) * t254 + t71;
t287 = t205 * t302;
t72 = t205 * t163 - t207 * t233;
t59 = -pkin(9) * t287 + t72;
t358 = qJD(4) * t243 + qJD(5) * t244 - t209 * t49 - t345 * t59;
t105 = -t209 * t172 + t345 * t173;
t357 = qJD(4) * t162 + qJD(5) * t105 - t209 * t59 + t345 * t49;
t200 = t212 * qJDD(2);
t279 = t210 * t296;
t355 = t279 - t200;
t325 = qJDD(3) * pkin(3);
t354 = qJDD(4) - t325;
t330 = cos(pkin(10));
t262 = t331 * t330;
t329 = sin(pkin(10));
t142 = t211 * t262 + t329 * t213;
t261 = t331 * t329;
t144 = -t211 * t261 + t330 * t213;
t318 = t206 * t211;
t148 = t210 * t318 - t331 * t212;
t276 = t206 * t329;
t277 = t206 * t330;
t239 = g(3) * t148 + g(2) * (t142 * t210 + t212 * t277) + g(1) * (t144 * t210 - t212 * t276);
t160 = qJDD(5) + t355;
t202 = pkin(11) + qJ(5);
t198 = sin(t202);
t353 = -t105 * t160 - t239 * t198;
t214 = qJD(3) ^ 2;
t297 = qJD(1) * qJD(2);
t281 = t211 * t297;
t317 = t206 * t213;
t258 = -qJDD(1) * t317 + t206 * t281;
t141 = t329 * t211 - t213 * t262;
t143 = t330 * t211 + t213 * t261;
t265 = g(1) * t143 + g(2) * t141;
t352 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t214 + t206 * (-g(3) * t213 + t281) - t258 + t265;
t295 = t205 * qJDD(3);
t222 = t207 * t360 + t295;
t309 = t360 * t205;
t257 = t207 * qJDD(3) - t309;
t27 = -qJD(5) * t245 + t209 * t222 - t345 * t257;
t247 = t345 * t103 + t209 * t89;
t351 = qJD(5) * t247 + t364 * t209 + t365 * t345;
t347 = qJ(6) * t300 - qJD(6) * t212 + t359;
t346 = -pkin(5) * t300 + t351;
t343 = pkin(5) * t160;
t341 = t245 * t82;
t131 = qJDD(2) * pkin(8) + (qJDD(1) * t211 + t213 * t297) * t206;
t272 = qJDD(1) * t331;
t263 = t210 * t272;
t41 = t263 + qJDD(3) * qJ(4) + t212 * t131 + (qJD(4) - t233) * qJD(3);
t57 = qJD(2) * t139 + qJDD(2) * t169 + t258;
t20 = t205 * t57 + t207 * t41;
t186 = t210 * t275;
t120 = t212 * t165 + t186;
t91 = pkin(4) * t287 + t120;
t338 = t310 * pkin(5) - t311 * qJ(6) - qJD(6) * t162 - t91;
t337 = -qJ(6) * t303 + t358;
t336 = pkin(5) * t303 + t357;
t335 = qJD(2) * pkin(2);
t110 = qJD(3) * qJ(4) + t120;
t288 = t213 * t306;
t121 = qJD(2) * t169 - t288;
t50 = -t110 * t205 + t207 * t121;
t31 = -pkin(4) * t302 - pkin(9) * t158 + t50;
t51 = t207 * t110 + t205 * t121;
t37 = pkin(9) * t157 + t51;
t14 = t209 * t31 + t345 * t37;
t334 = t14 * t192;
t270 = t207 * t293;
t332 = -t270 + t366;
t327 = qJ(6) * t160;
t323 = t198 * t212;
t199 = cos(t202);
t322 = t199 * t212;
t13 = -t209 * t37 + t345 * t31;
t313 = qJD(6) - t13;
t312 = qJDD(1) - g(3);
t308 = pkin(2) * t317 + pkin(8) * t318;
t299 = qJD(3) * t212;
t284 = t205 * t299;
t152 = pkin(4) * t284 + pkin(8) * t299;
t164 = pkin(4) * t321 + t210 * pkin(8);
t203 = t210 ^ 2;
t307 = -t212 ^ 2 + t203;
t305 = qJD(1) * t210;
t304 = qJD(2) * t206;
t292 = t198 * t317;
t196 = t207 * pkin(4) + pkin(3);
t291 = pkin(4) * t205 + pkin(8);
t289 = qJ(4) * t300;
t286 = t211 * t304;
t285 = t213 * t304;
t19 = -t205 * t41 + t207 * t57;
t12 = pkin(4) * t355 - t222 * pkin(9) + t19;
t17 = pkin(9) * t257 + t20;
t278 = t345 * t12 - t209 * t17 - t37 * t282 - t31 * t298;
t274 = -qJD(3) * pkin(3) + qJD(4);
t122 = -pkin(8) * t319 + t156;
t273 = qJD(2) * t122 + t50;
t271 = -qJDD(2) * t122 - t19;
t269 = t82 * t288;
t268 = t245 * t288;
t264 = g(1) * t144 + g(2) * t142;
t256 = t196 * t212 + t210 * t340;
t215 = qJD(2) ^ 2;
t255 = qJDD(2) * t213 - t211 * t215;
t149 = t331 * t210 + t212 * t318;
t95 = -t149 * t205 - t207 * t317;
t96 = t149 * t207 - t205 * t317;
t250 = -t209 * t96 + t345 * t95;
t40 = t209 * t95 + t345 * t96;
t249 = t209 * t12 + t345 * t17 + t31 * t282 - t37 * t298;
t248 = -t209 * t103 + t345 * t89;
t114 = -t199 * t318 + t212 * t292;
t62 = -t141 * t323 - t142 * t199;
t64 = -t143 * t323 - t144 * t199;
t241 = -g(1) * t64 - g(2) * t62 - g(3) * t114;
t115 = (t198 * t211 + t199 * t314) * t206;
t63 = -t141 * t322 + t142 * t198;
t65 = -t143 * t322 + t144 * t198;
t240 = -g(1) * t65 - g(2) * t63 - g(3) * t115;
t102 = qJD(3) * t149 + t210 * t285;
t101 = -qJD(3) * t148 + t212 * t285;
t69 = -t101 * t205 + t207 * t286;
t70 = t101 * t207 + t205 * t286;
t11 = t40 * qJD(5) + t209 * t70 - t345 * t69;
t238 = t102 * t82 + t11 * t192 + t148 * t27 + t160 * t250;
t236 = t207 * t294 + t295;
t100 = t144 * t212 + t210 * t276;
t98 = t142 * t212 - t210 * t277;
t235 = g(1) * t100 + g(2) * t98 + g(3) * t149;
t26 = -t157 * t282 + t158 * t298 - t209 * t257 - t345 * t222;
t229 = -qJD(3) * t186 - t210 * t131 - t165 * t299 + t212 * t272;
t43 = -t229 + t354;
t234 = t239 - t43;
t231 = g(3) * t317 - t265;
t230 = -g(3) * t318 - t264;
t228 = t26 - t362;
t106 = t233 + t274;
t10 = t250 * qJD(5) + t209 * t69 + t345 * t70;
t227 = t10 * t192 - t102 * t245 - t148 * t26 - t160 * t40;
t45 = -t141 * t199 + t198 * t98;
t47 = t100 * t198 - t143 * t199;
t87 = t149 * t198 + t199 * t317;
t226 = g(1) * t47 + g(2) * t45 + g(3) * t87 + t278;
t225 = t234 + t325;
t223 = t160 * t244 + t199 * t239;
t166 = -t288 - t335;
t221 = -pkin(8) * qJDD(3) + (t166 + t288 - t335) * qJD(3);
t73 = -t157 * pkin(4) + t106;
t46 = t141 * t198 + t199 * t98;
t48 = t100 * t199 + t143 * t198;
t88 = t149 * t199 - t292;
t220 = -g(1) * t48 - g(2) * t46 - g(3) * t88 + t249;
t24 = t82 * pkin(5) + qJ(6) * t245 + t73;
t219 = -t24 * t245 + qJDD(6) - t226;
t217 = t229 + t239;
t28 = -pkin(4) * t257 + t43;
t216 = t27 + t361;
t3 = t27 * pkin(5) + t26 * qJ(6) + qJD(6) * t245 + t28;
t138 = t143 * pkin(2);
t137 = t141 * pkin(2);
t135 = t243 * t210;
t134 = t162 * t210;
t80 = -pkin(5) * t243 - qJ(6) * t162 - t196;
t75 = qJD(3) * t232 + t210 * t356;
t74 = t210 * t146 + t209 * t284 - t290 * t299;
t58 = pkin(5) * t134 - qJ(6) * t135 + t164;
t36 = -pkin(5) * t245 + qJ(6) * t82;
t35 = t212 * pkin(5) - t248;
t34 = -qJ(6) * t212 + t247;
t21 = pkin(5) * t75 + qJ(6) * t74 - qJD(6) * t135 + t152;
t18 = -t26 - t362;
t9 = -t192 * qJ(6) + t14;
t5 = t192 * pkin(5) + t313;
t2 = qJDD(6) - t278 - t343;
t1 = -qJD(6) * t192 + t249 + t327;
t4 = [t312, 0, t255 * t206 (-qJDD(2) * t211 - t213 * t215) * t206, 0, 0, 0, 0, 0, -qJD(3) * t102 - qJDD(3) * t148 + (t212 * t255 - t213 * t279) * t206, -qJD(3) * t101 - qJDD(3) * t149 + (-t210 * t255 - t213 * t280) * t206, -t95 * t200 - t102 * t157 - t148 * t257 + (-t212 * t69 + t300 * t95) * qJD(2), t96 * t200 + t102 * t158 + t148 * t236 + (t70 * t212 + (t148 * t315 - t210 * t96) * qJD(3)) * qJD(2), t70 * t157 - t69 * t158 - t222 * t95 + t96 * t257, t102 * t106 + t148 * t43 + t19 * t95 + t20 * t96 + t50 * t69 + t51 * t70 - g(3), 0, 0, 0, 0, 0, t238, t227, t238, -t10 * t82 - t11 * t245 + t250 * t26 - t27 * t40, -t227, t1 * t40 + t10 * t9 + t102 * t24 + t11 * t5 + t148 * t3 - t2 * t250 - g(3); 0, qJDD(2), t312 * t317 + t265, -t312 * t318 + t264, qJDD(2) * t203 + 0.2e1 * t212 * t279, 0.2e1 * t210 * t200 - 0.2e1 * t307 * t296, qJDD(3) * t210 + t212 * t214, qJDD(3) * t212 - t210 * t214, 0, t221 * t210 + t212 * t352, -t210 * t352 + t221 * t212, t230 * t205 + (-pkin(8) * t257 + qJD(3) * t273 + t157 * t288 + t43 * t205) * t210 + ((-pkin(8) * t157 + t106 * t205) * qJD(3) - t333 * qJD(2) - t231 * t207 + t271) * t212, t230 * t207 + (-t158 * t288 + t43 * t207 + (-qJD(2) * t123 - t51) * qJD(3) + t236 * pkin(8)) * t210 + (t123 * qJDD(2) + t20 + (pkin(8) * t158 + t106 * t207) * qJD(3) + t231 * t205 + (t270 + t332) * qJD(2)) * t212, t123 * t257 - t122 * t295 - t333 * t158 + t332 * t157 + (-t205 * t51 - t207 * t273) * t299 + (-t20 * t205 + t207 * t271 - t231) * t210, t20 * t123 + t19 * t122 - g(1) * (-t143 * t260 - t138) - g(2) * (-t141 * t260 - t137) - g(3) * t308 + t332 * t51 + t333 * t50 + (-g(3) * t260 - t106 * t305) * t317 + (t106 * t299 + t43 * t210 - t264) * pkin(8), -t135 * t26 + t245 * t74, t134 * t26 - t135 * t27 + t245 * t75 + t74 * t82, t135 * t160 + t192 * t74 + t212 * t26 - t245 * t300, -t134 * t160 + t192 * t75 + t212 * t27 - t300 * t82, -t160 * t212 - t192 * t300, t13 * t300 + t28 * t134 + t152 * t82 + t248 * t160 + t164 * t27 + t192 * t351 - t210 * t269 - t278 * t212 + t73 * t75 + t240, -t247 * t160 + t249 * t212 - t152 * t245 - t164 * t26 + t28 * t135 - t73 * t74 + (-t14 * qJD(3) + t268) * t210 + t359 * t192 - t241, t134 * t3 - t160 * t35 + t2 * t212 + t21 * t82 + t24 * t75 + t27 * t58 + (-qJD(3) * t5 - t269) * t210 + t346 * t192 + t240, -t1 * t134 + t135 * t2 - t231 * t210 - t245 * t346 - t26 * t35 - t27 * t34 - t347 * t82 - t5 * t74 - t75 * t9, -t1 * t212 - t135 * t3 + t160 * t34 + t21 * t245 + t24 * t74 + t26 * t58 + (qJD(3) * t9 - t268) * t210 - t347 * t192 + t241, t1 * t34 + t3 * t58 + t24 * t21 + t2 * t35 - g(1) * (pkin(5) * t65 + qJ(6) * t64 - t143 * t256 + t144 * t291 - t138) - g(2) * (pkin(5) * t63 + qJ(6) * t62 - t141 * t256 + t142 * t291 - t137) - g(3) * (pkin(5) * t115 + qJ(6) * t114 + t308) + t347 * t9 + t346 * t5 + (-g(3) * pkin(4) * t320 + (-g(3) * t256 - t24 * t305) * t213) * t206; 0, 0, 0, 0, -t210 * t215 * t212, t307 * t215, t294, t200, qJDD(3), t120 * qJD(3) - t166 * t303 + t217, -t263 + (-qJD(2) * t166 - t131) * t212 + t235, t205 * qJ(4) * t200 - pkin(3) * t309 + t120 * t157 + t225 * t207 + (-t50 * t210 + t71 * t212 + (-t289 + (qJD(4) - t106) * t212) * t205) * qJD(2), -t120 * t158 - t259 * t207 * qJDD(2) - t225 * t205 + (t51 * t210 - t72 * t212 + (-t289 + (-t106 + t274) * t212) * t207) * qJD(2), -t72 * t157 + t71 * t158 + (qJ(4) * t257 + qJD(4) * t157 + t302 * t50 + t20) * t207 + (qJ(4) * t222 + qJD(4) * t158 + t302 * t51 - t19) * t205 - t235, -t106 * t120 - t50 * t71 - t51 * t72 + (-t205 * t50 + t207 * t51) * qJD(4) + t234 * pkin(3) + (-t19 * t205 + t20 * t207 - t235) * qJ(4), -t162 * t26 - t245 * t311, -t162 * t27 - t243 * t26 + t245 * t310 - t311 * t82, t160 * t162 - t192 * t311 + t245 * t303, t160 * t243 + t192 * t310 + t303 * t82, t192 * t303, -t13 * t303 + t192 * t357 - t196 * t27 - t243 * t28 + t310 * t73 - t91 * t82 + t223, t14 * t303 + t28 * t162 + t192 * t358 + t196 * t26 + t245 * t91 + t311 * t73 + t353, t192 * t336 + t24 * t310 - t243 * t3 + t27 * t80 + t303 * t5 + t338 * t82 + t223, t1 * t243 - t105 * t27 + t162 * t2 + t244 * t26 - t245 * t336 - t310 * t9 + t311 * t5 - t337 * t82 - t235, -t162 * t3 - t192 * t337 - t24 * t311 + t245 * t338 + t26 * t80 - t303 * t9 - t353, t1 * t105 - t2 * t244 - t340 * t235 + t24 * t338 + t3 * t80 + t336 * t5 + t337 * t9 + t239 * (pkin(5) * t199 + qJ(6) * t198 + t196); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t158 * t302 - t257 (-t157 + t301) * t302 + t236, -t157 ^ 2 - t158 ^ 2, -t51 * t157 + t50 * t158 - t217 + t354, 0, 0, 0, 0, 0, t216, -t228, t216, -t350 - t363, t228, t245 * t5 + t82 * t9 - t239 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t341, t350 - t363, t18, -t27 + t361, t160, t245 * t73 + t226 - t334, -t13 * t192 + t73 * t82 - t220, -t36 * t82 - t219 - t334 + 0.2e1 * t343, pkin(5) * t26 - qJ(6) * t27 - (-t14 + t9) * t245 + (t5 - t313) * t82, 0.2e1 * t327 - t24 * t82 - t36 * t245 + (-0.2e1 * qJD(6) + t13) * t192 + t220, t1 * qJ(6) - t2 * pkin(5) - t24 * t36 - t5 * t14 - g(1) * (-pkin(5) * t47 + qJ(6) * t48) - g(2) * (-pkin(5) * t45 + qJ(6) * t46) - g(3) * (-pkin(5) * t87 + qJ(6) * t88) + t313 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160 - t341, t18, -t192 ^ 2 - t350, t192 * t9 + t219 - t343;];
tau_reg  = t4;
