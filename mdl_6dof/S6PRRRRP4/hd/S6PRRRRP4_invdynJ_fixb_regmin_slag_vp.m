% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRRRP4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:17:50
% EndTime: 2019-03-09 00:18:07
% DurationCPUTime: 7.17s
% Computational Cost: add. (7404->564), mult. (16552->758), div. (0->0), fcn. (12816->14), ass. (0->279)
t220 = sin(qJ(3));
t223 = cos(qJ(3));
t282 = pkin(3) * t220 - pkin(9) * t223;
t166 = t282 * qJD(3);
t219 = sin(qJ(4));
t221 = sin(qJ(2));
t222 = cos(qJ(4));
t333 = qJD(3) * t220;
t216 = sin(pkin(6));
t340 = qJD(1) * t216;
t224 = cos(qJ(2));
t349 = t223 * t224;
t385 = pkin(8) * t219;
t417 = (-t219 * t349 + t221 * t222) * t340 - t222 * t166 - t333 * t385;
t173 = -pkin(3) * t223 - pkin(9) * t220 - pkin(2);
t327 = qJD(4) * t222;
t416 = -(t219 * t221 + t222 * t349) * t340 + t219 * t166 + t173 * t327;
t350 = t222 * t223;
t200 = pkin(8) * t350;
t276 = pkin(4) * t220 - pkin(10) * t350;
t415 = -t276 * qJD(3) - (-t200 + (pkin(10) * t220 - t173) * t219) * qJD(4) + t417;
t328 = qJD(4) * t219;
t332 = qJD(3) * t222;
t331 = qJD(3) * t223;
t304 = t219 * t331;
t409 = t220 * t327 + t304;
t414 = -t409 * pkin(10) + (-t220 * t332 - t223 * t328) * pkin(8) + t416;
t336 = qJD(2) * t223;
t307 = t219 * t336;
t225 = -pkin(10) - pkin(9);
t312 = qJD(4) * t225;
t163 = t282 * qJD(2);
t168 = qJD(2) * pkin(8) + t221 * t340;
t217 = cos(pkin(6));
t354 = t217 * t223;
t401 = qJD(1) * t354 - t220 * t168;
t345 = t219 * t163 + t222 * t401;
t413 = pkin(10) * t307 + t219 * t312 - t345;
t150 = t222 * t163;
t412 = -qJD(2) * t276 + t219 * t401 + t222 * t312 - t150;
t337 = qJD(2) * t220;
t159 = -t219 * t337 + t332;
t334 = qJD(3) * t219;
t160 = t222 * t337 + t334;
t218 = sin(qJ(5));
t387 = cos(qJ(5));
t265 = t218 * t159 + t160 * t387;
t112 = -qJD(3) * pkin(3) - t401;
t80 = -pkin(4) * t159 + t112;
t89 = -t387 * t159 + t160 * t218;
t27 = pkin(5) * t89 - qJ(6) * t265 + t80;
t411 = t27 * t89;
t410 = t80 * t89;
t377 = t265 * t89;
t353 = t218 * t219;
t263 = t222 * t387 - t353;
t398 = qJD(4) + qJD(5);
t302 = t387 * qJD(5);
t399 = t387 * qJD(4) + t302;
t367 = t399 * t222 - t263 * t336 - t398 * t353;
t162 = t218 * t222 + t219 * t387;
t100 = t398 * t162;
t346 = -t162 * t336 + t100;
t289 = qJD(4) + t336;
t317 = t220 * qJDD(2);
t245 = qJD(3) * t289 + t317;
t321 = qJD(2) * qJD(4);
t298 = t220 * t321;
t278 = -qJDD(3) + t298;
t254 = t278 * t219;
t230 = t245 * t222 - t254;
t275 = qJD(3) * qJD(4) + t317;
t322 = qJD(2) * qJD(3);
t300 = t223 * t322;
t244 = t275 + t300;
t389 = t265 ^ 2;
t408 = -t89 ^ 2 + t389;
t198 = -qJD(4) + t336;
t186 = -qJD(5) + t198;
t286 = t244 * t219 + t222 * t298;
t249 = t222 * qJDD(3) - t286;
t326 = qJD(5) * t218;
t25 = -t159 * t302 + t160 * t326 - t218 * t249 - t387 * t230;
t17 = -t186 * t89 - t25;
t44 = pkin(5) * t265 + qJ(6) * t89;
t406 = pkin(4) * t219 + pkin(8);
t342 = t219 * t173 + t200;
t352 = t219 * t220;
t109 = -pkin(10) * t352 + t342;
t158 = t222 * t173;
t351 = t220 * t222;
t96 = -pkin(10) * t351 + t158 + (-pkin(4) - t385) * t223;
t405 = -t109 * t326 - t218 * t415 + t96 * t302 + t414 * t387;
t180 = t225 * t219;
t181 = t225 * t222;
t264 = t180 * t387 + t218 * t181;
t404 = qJD(5) * t264 + t218 * t412 + t387 * t413;
t116 = t218 * t180 - t181 * t387;
t403 = qJD(5) * t116 + t218 * t413 - t387 * t412;
t339 = qJD(1) * t220;
t193 = t217 * t339;
t123 = t223 * t168 + t193;
t113 = qJD(3) * pkin(9) + t123;
t311 = t224 * t340;
t125 = qJD(2) * t173 - t311;
t67 = t222 * t113 + t219 * t125;
t43 = pkin(10) * t159 + t67;
t284 = -t123 + (-t307 + t328) * pkin(4);
t208 = t223 * qJDD(2);
t156 = t220 * t322 + qJDD(4) - t208;
t152 = qJDD(5) + t156;
t137 = t152 * qJ(6);
t172 = t186 * qJD(6);
t402 = t137 - t172;
t366 = cos(pkin(11));
t291 = t366 * t221;
t215 = sin(pkin(11));
t358 = t215 * t224;
t140 = t217 * t291 + t358;
t292 = t216 * t366;
t102 = t140 * t223 - t220 * t292;
t290 = t366 * t224;
t359 = t215 * t221;
t142 = -t217 * t359 + t290;
t104 = t215 * t216 * t220 + t142 * t223;
t356 = t216 * t223;
t147 = t217 * t220 + t221 * t356;
t355 = t216 * t224;
t107 = -t147 * t219 - t222 * t355;
t139 = -t217 * t290 + t359;
t141 = t217 * t358 + t291;
t400 = -g(1) * (-t104 * t219 + t141 * t222) - g(2) * (-t102 * t219 + t139 * t222) - g(3) * t107;
t138 = t152 * pkin(5);
t397 = t138 - qJDD(6);
t357 = t216 * t221;
t146 = t220 * t357 - t354;
t396 = g(3) * t146 - g(2) * (-t140 * t220 - t223 * t292) - g(1) * (-t142 * t220 + t215 * t356);
t323 = qJD(1) * qJD(2);
t131 = qJDD(2) * pkin(8) + (qJDD(1) * t221 + t224 * t323) * t216;
t319 = qJDD(1) * t217;
t297 = t220 * t319;
t52 = qJDD(3) * pkin(9) + qJD(3) * t401 + t131 * t223 + t297;
t301 = t221 * t323;
t279 = -qJDD(1) * t355 + t216 * t301;
t75 = qJD(2) * t166 + qJDD(2) * t173 + t279;
t316 = t125 * t327 + t219 * t75 + t222 * t52;
t258 = -t113 * t328 + t316;
t13 = pkin(10) * t249 + t258;
t66 = -t113 * t219 + t222 * t125;
t42 = -pkin(10) * t160 + t66;
t33 = -pkin(4) * t198 + t42;
t74 = t222 * t75;
t293 = -t219 * t52 + t74;
t318 = t219 * qJDD(3);
t6 = -(t318 + (t300 + t317) * t222) * pkin(10) + t156 * pkin(4) - t43 * qJD(4) + t293;
t295 = t218 * t13 + t43 * t302 + t33 * t326 - t387 * t6;
t214 = qJ(4) + qJ(5);
t209 = sin(t214);
t210 = cos(t214);
t61 = t102 * t209 - t139 * t210;
t63 = t104 * t209 - t141 * t210;
t94 = t147 * t209 + t210 * t355;
t242 = g(1) * t63 + g(2) * t61 + g(3) * t94 - t295;
t232 = t265 * t27 - t242 - t397;
t395 = -t80 * t265 + t242;
t26 = qJD(5) * t265 + t230 * t218 - t387 * t249;
t394 = -t186 * t265 - t26;
t393 = -t116 * t152 - t209 * t396;
t226 = qJD(3) ^ 2;
t380 = g(2) * t139;
t382 = g(1) * t141;
t281 = t380 + t382;
t392 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t226 + t216 * (-g(3) * t224 + t301) - t279 + t281;
t248 = g(3) * t355 - t281;
t391 = qJD(4) * (pkin(8) * t198 + t113) - t248;
t269 = t387 * t109 + t218 * t96;
t390 = qJD(5) * t269 + t414 * t218 + t387 * t415;
t388 = qJ(6) * t333 - qJD(6) * t223 + t405;
t376 = -pkin(5) * t333 + t390;
t375 = t346 * pkin(5) - t367 * qJ(6) - qJD(6) * t162 + t284;
t373 = -qJ(6) * t337 + t404;
t372 = pkin(5) * t337 + t403;
t371 = qJD(2) * pkin(2);
t314 = t387 * t43;
t19 = t218 * t33 + t314;
t370 = t186 * t19;
t369 = t218 * t43;
t21 = t387 * t42 - t369;
t368 = pkin(4) * t302 + qJD(6) - t21;
t364 = qJDD(3) * pkin(3);
t362 = t160 * t198;
t361 = t209 * t223;
t360 = t210 * t223;
t18 = t33 * t387 - t369;
t348 = qJD(6) - t18;
t347 = qJDD(1) - g(3);
t167 = pkin(4) * t352 + t220 * pkin(8);
t212 = t220 ^ 2;
t341 = -t223 ^ 2 + t212;
t338 = qJD(2) * t216;
t335 = qJD(3) * t159;
t330 = qJD(4) * t159;
t329 = qJD(4) * t198;
t325 = t112 * qJD(4);
t324 = t160 * qJD(3);
t313 = t209 * t355;
t124 = pkin(4) * t409 + pkin(8) * t331;
t206 = pkin(4) * t222 + pkin(3);
t309 = t221 * t338;
t308 = t224 * t338;
t306 = t220 * t328;
t299 = t224 * t322;
t296 = t387 * t13 + t218 * t6 + t33 * t302 - t43 * t326;
t288 = t89 * t311;
t287 = t265 * t311;
t285 = t387 * t331;
t20 = t218 * t42 + t314;
t283 = pkin(4) * t326 - t20;
t280 = g(1) * t142 + g(2) * t140;
t227 = qJD(2) ^ 2;
t277 = qJDD(2) * t224 - t221 * t227;
t274 = t206 * t223 - t220 * t225 + pkin(2);
t270 = -t218 * t109 + t387 * t96;
t267 = -t147 * t222 + t219 * t355;
t266 = t107 * t387 + t218 * t267;
t48 = t218 * t107 - t267 * t387;
t261 = t156 * t219 - t198 * t327;
t260 = t222 * t156 + t198 * t328;
t106 = qJD(3) * t147 + t220 * t308;
t105 = -qJD(3) * t146 + t223 * t308;
t37 = qJD(4) * t267 - t105 * t219 + t222 * t309;
t38 = qJD(4) * t107 + t105 * t222 + t219 * t309;
t8 = qJD(5) * t48 + t218 * t38 - t37 * t387;
t257 = t106 * t89 + t146 * t26 + t152 * t266 + t186 * t8;
t117 = -t210 * t357 + t223 * t313;
t76 = -t139 * t361 - t140 * t210;
t78 = -t141 * t361 - t142 * t210;
t256 = -g(1) * t78 - g(2) * t76 - g(3) * t117;
t118 = (t209 * t221 + t210 * t349) * t216;
t77 = -t139 * t360 + t140 * t209;
t79 = -t141 * t360 + t142 * t209;
t255 = -g(1) * t79 - g(2) * t77 - g(3) * t118;
t253 = qJDD(2) * t222 - t219 * t321;
t251 = g(1) * t104 + g(2) * t102 + g(3) * t147;
t250 = qJD(3) * t193 + t220 * t131 + t168 * t331 - t223 * t319;
t247 = -g(3) * t357 - t280;
t246 = -pkin(9) * t156 - t112 * t198;
t53 = t250 - t364;
t7 = qJD(5) * t266 + t218 * t37 + t38 * t387;
t243 = t106 * t265 - t146 * t25 - t152 * t48 + t186 * t7;
t62 = t102 * t210 + t139 * t209;
t64 = t104 * t210 + t141 * t209;
t95 = t147 * t210 - t313;
t241 = g(1) * t64 + g(2) * t62 + g(3) * t95 - t296;
t238 = -pkin(9) * t329 - t396 + t53;
t237 = t152 * t264 + t210 * t396;
t169 = -t311 - t371;
t236 = -pkin(8) * qJDD(3) + (t169 + t311 - t371) * qJD(3);
t234 = -t18 * t186 + t241;
t231 = -g(1) * (-t63 * pkin(5) + qJ(6) * t64) - g(2) * (-t61 * pkin(5) + qJ(6) * t62) - g(3) * (-t94 * pkin(5) + qJ(6) * t95);
t28 = -pkin(4) * t249 + t53;
t205 = -pkin(4) * t387 - pkin(5);
t202 = pkin(4) * t218 + qJ(6);
t135 = t263 * t220;
t134 = t162 * t220;
t86 = -pkin(5) * t263 - qJ(6) * t162 - t206;
t71 = pkin(5) * t134 - qJ(6) * t135 + t167;
t55 = t219 * t285 - t218 * t306 - t326 * t352 + (t218 * t331 + t399 * t220) * t222;
t54 = t100 * t220 + t218 * t304 - t222 * t285;
t41 = t223 * pkin(5) - t270;
t40 = -qJ(6) * t223 + t269;
t30 = pkin(4) * t160 + t44;
t16 = -t186 * qJ(6) + t19;
t15 = pkin(5) * t55 + qJ(6) * t54 - qJD(6) * t135 + t124;
t14 = t186 * pkin(5) + t348;
t3 = t26 * pkin(5) + t25 * qJ(6) - qJD(6) * t265 + t28;
t2 = t295 - t397;
t1 = t296 + t402;
t4 = [t347, 0, t277 * t216 (-qJDD(2) * t221 - t224 * t227) * t216, 0, 0, 0, 0, 0, -qJD(3) * t106 - qJDD(3) * t146 + (-t220 * t299 + t223 * t277) * t216, -qJD(3) * t105 - qJDD(3) * t147 + (-t220 * t277 - t223 * t299) * t216, 0, 0, 0, 0, 0, -t106 * t159 + t107 * t156 - t146 * t249 - t37 * t198, t106 * t160 + t230 * t146 + t156 * t267 + t38 * t198, 0, 0, 0, 0, 0, t257, t243, t257, t25 * t266 - t26 * t48 + t265 * t8 - t7 * t89, -t243, t1 * t48 + t106 * t27 + t14 * t8 + t146 * t3 + t16 * t7 - t2 * t266 - g(3); 0, qJDD(2), t347 * t355 + t281, -t347 * t357 + t280, qJDD(2) * t212 + 0.2e1 * t220 * t300, 0.2e1 * t208 * t220 - 0.2e1 * t322 * t341, qJDD(3) * t220 + t223 * t226, qJDD(3) * t223 - t220 * t226, 0, t236 * t220 + t392 * t223, -t392 * t220 + t236 * t223, -t160 * t306 + (-t220 * t254 + t223 * t324 + t244 * t351) * t222 (t222 * t159 - t160 * t219) * t331 + ((-t160 * qJD(4) + t249) * t222 + (-t330 - t230) * t219) * t220 (-t318 + (-t198 - t289) * t332) * t223 + (-t223 * t253 + t260 + t324) * t220 (t198 * t334 - t249) * t223 + (-t261 + t335) * t220, -t156 * t223 - t198 * t333, t158 * t156 + t417 * t198 + (t173 * t329 + t247) * t219 + (-pkin(8) * t335 - t74 + (-pkin(8) * t156 + qJD(3) * t112 + qJD(4) * t125 + t52) * t219 + t391 * t222) * t223 + (-pkin(8) * t249 + t66 * qJD(3) + t159 * t311 + t53 * t219 + t222 * t325) * t220, -t342 * t156 + t416 * t198 + t247 * t222 + ((pkin(8) * t160 + t112 * t222) * qJD(3) - t391 * t219 + t316) * t223 + (-t160 * t311 - t219 * t325 - t67 * qJD(3) + t53 * t222 + (t318 + t253 * t220 + (-t198 + t289) * t332) * pkin(8)) * t220, -t135 * t25 - t265 * t54, t134 * t25 - t135 * t26 - t265 * t55 + t54 * t89, t135 * t152 + t186 * t54 + t223 * t25 + t265 * t333, -t134 * t152 + t186 * t55 + t223 * t26 - t333 * t89, -t152 * t223 - t186 * t333, t124 * t89 + t28 * t134 + t270 * t152 + t167 * t26 + t18 * t333 + t390 * t186 - t220 * t288 + t295 * t223 + t80 * t55 + t255, -t269 * t152 + t296 * t223 + t124 * t265 - t167 * t25 + t28 * t135 - t80 * t54 + (-t19 * qJD(3) - t287) * t220 + t405 * t186 - t256, t134 * t3 + t15 * t89 - t152 * t41 + t2 * t223 + t26 * t71 + t27 * t55 + (-qJD(3) * t14 - t288) * t220 + t376 * t186 + t255, -t1 * t134 + t135 * t2 - t14 * t54 - t16 * t55 - t220 * t248 - t25 * t41 - t26 * t40 + t265 * t376 - t388 * t89, -t1 * t223 - t135 * t3 - t15 * t265 + t152 * t40 + t25 * t71 + t27 * t54 + (qJD(3) * t16 + t287) * t220 - t388 * t186 + t256, t1 * t40 + t3 * t71 + t27 * t15 + t2 * t41 - g(1) * (pkin(5) * t79 + qJ(6) * t78 + t406 * t142) - g(2) * (pkin(5) * t77 + qJ(6) * t76 + t406 * t140) - g(3) * (pkin(5) * t118 + qJ(6) * t117) + t388 * t16 + t274 * t382 + t376 * t14 + t274 * t380 + (-g(3) * t406 * t221 + (-g(3) * t274 - t27 * t339) * t224) * t216; 0, 0, 0, 0, -t220 * t227 * t223, t341 * t227, t317, t208, qJDD(3), qJD(3) * t123 - t169 * t337 - t250 + t396, -t297 + (-qJD(2) * t169 - t131) * t223 + t251, -t278 * t219 ^ 2 + (t219 * t245 - t362) * t222 (-t286 + t362) * t219 + (t330 + 0.2e1 * t318 + t275 * t222 + (-t306 + (-t159 + t332) * t223) * qJD(2)) * t222 (-t160 * t220 + t198 * t350) * qJD(2) + t261 (-t198 * t219 * t223 - t159 * t220) * qJD(2) + t260, t198 * t337, -pkin(3) * t286 + t150 * t198 - t66 * t337 + t123 * t159 + (-t198 * t401 + t246) * t219 + (-t238 + t364) * t222, -t345 * t198 + t67 * t337 - t123 * t160 + (-pkin(3) * t244 + t246) * t222 + (pkin(3) * t278 + t238) * t219, -t162 * t25 + t265 * t367, -t162 * t26 - t25 * t263 - t265 * t346 - t367 * t89, t152 * t162 - t186 * t367 - t265 * t337, t152 * t263 + t186 * t346 + t337 * t89, t186 * t337, -t18 * t337 + t403 * t186 - t206 * t26 - t263 * t28 + t284 * t89 + t346 * t80 + t237, t28 * t162 + t404 * t186 + t19 * t337 + t206 * t25 + t284 * t265 + t367 * t80 + t393, t14 * t337 + t186 * t372 + t26 * t86 - t263 * t3 + t27 * t346 + t375 * t89 + t237, t1 * t263 - t116 * t26 + t14 * t367 - t16 * t346 + t162 * t2 + t25 * t264 + t265 * t372 - t373 * t89 - t251, -t16 * t337 - t162 * t3 - t186 * t373 + t25 * t86 - t265 * t375 - t27 * t367 - t393, t1 * t116 - t2 * t264 + t14 * t372 + t16 * t373 + t225 * t251 + t27 * t375 + t3 * t86 + t396 * (pkin(5) * t210 + qJ(6) * t209 + t206); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160 * t159, -t159 ^ 2 + t160 ^ 2, t159 * t198 + t230, t249 - t362, t156, -t112 * t160 + t293 + (-qJD(4) - t198) * t67 + t400, -t66 * t198 - t112 * t159 - g(1) * (-t104 * t222 - t141 * t219) - g(2) * (-t102 * t222 - t139 * t219) - g(3) * t267 - t258, t377, t408, t17, t394, t152, -t20 * t186 + (t152 * t387 - t160 * t89 + t186 * t326) * pkin(4) + t395, -t21 * t186 + t410 + (-t152 * t218 - t160 * t265 + t186 * t302) * pkin(4) + t241, -t152 * t205 + t186 * t283 - t30 * t89 - t232, -t202 * t26 - t205 * t25 + (t16 + t283) * t265 + (t14 - t368) * t89, t152 * t202 - t186 * t368 + t265 * t30 - t241 + t402 - t411, t1 * t202 + t2 * t205 - t27 * t30 - t14 * t20 + t368 * t16 + (t14 * t326 + t400) * pkin(4) + t231; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t377, t408, t17, t394, t152, -t370 + t395, t234 + t410, -t44 * t89 + t138 - t232 - t370, pkin(5) * t25 - qJ(6) * t26 + (t16 - t19) * t265 + (t14 - t348) * t89, t265 * t44 + 0.2e1 * t137 - 0.2e1 * t172 - t234 - t411, -t2 * pkin(5) + t1 * qJ(6) - t14 * t19 + t16 * t348 - t27 * t44 + t231; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152 + t377, t17, -t186 ^ 2 - t389, t16 * t186 + t232;];
tau_reg  = t4;
