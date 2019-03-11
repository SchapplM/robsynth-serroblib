% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRRPR4
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRRPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:20:32
% EndTime: 2019-03-08 23:20:46
% DurationCPUTime: 6.51s
% Computational Cost: add. (6339->511), mult. (14753->733), div. (0->0), fcn. (11645->16), ass. (0->265)
t245 = cos(qJ(3));
t337 = qJD(2) * t245;
t216 = -qJD(4) + t337;
t209 = -qJD(6) + t216;
t240 = sin(qJ(4));
t244 = cos(qJ(4));
t326 = t244 * qJD(3);
t241 = sin(qJ(3));
t338 = qJD(2) * t241;
t191 = t240 * t338 - t326;
t335 = qJD(3) * t240;
t193 = t244 * t338 + t335;
t234 = sin(pkin(12));
t236 = cos(pkin(12));
t120 = t191 * t234 - t193 * t236;
t239 = sin(qJ(6));
t243 = cos(qJ(6));
t278 = -t191 * t236 - t193 * t234;
t352 = t243 * t278;
t58 = t120 * t239 + t352;
t373 = t209 * t58;
t327 = qJD(6) * t239;
t324 = qJD(2) * qJD(3);
t304 = t245 * t324;
t322 = t241 * qJDD(2);
t330 = qJD(4) * t241;
t403 = -qJD(2) * t330 + qJDD(3);
t112 = qJD(4) * t326 + (t304 + t322) * t244 + t403 * t240;
t113 = (qJD(3) * (qJD(4) + t337) + t322) * t240 - t403 * t244;
t45 = -t112 * t234 - t113 * t236;
t46 = t112 * t236 - t113 * t234;
t8 = qJD(6) * t352 + t120 * t327 + t239 * t45 + t243 * t46;
t421 = t8 + t373;
t402 = -t243 * t120 + t239 * t278;
t420 = t402 * t58;
t287 = pkin(3) * t241 - pkin(9) * t245;
t195 = t287 * qJD(3);
t242 = sin(qJ(2));
t334 = qJD(3) * t241;
t235 = sin(pkin(6));
t342 = qJD(1) * t235;
t246 = cos(qJ(2));
t350 = t245 * t246;
t387 = pkin(8) * t240;
t419 = (-t240 * t350 + t242 * t244) * t342 - t244 * t195 - t334 * t387;
t199 = -pkin(3) * t245 - pkin(9) * t241 - pkin(2);
t329 = qJD(4) * t244;
t418 = -(t240 * t242 + t244 * t350) * t342 + t240 * t195 + t199 * t329;
t417 = t402 ^ 2 - t58 ^ 2;
t237 = cos(pkin(6));
t366 = sin(pkin(11));
t291 = t366 * t246;
t367 = cos(pkin(11));
t294 = t367 * t242;
t167 = t237 * t294 + t291;
t296 = t235 * t367;
t128 = t167 * t245 - t241 * t296;
t292 = t366 * t242;
t293 = t367 * t246;
t169 = -t237 * t292 + t293;
t295 = t235 * t366;
t130 = t169 * t245 + t241 * t295;
t166 = -t237 * t293 + t292;
t168 = t237 * t291 + t294;
t357 = t235 * t242;
t174 = t237 * t241 + t245 * t357;
t229 = qJ(4) + pkin(12) + qJ(6);
t220 = sin(t229);
t221 = cos(t229);
t196 = qJD(2) * pkin(8) + t242 * t342;
t341 = qJD(1) * t241;
t214 = t237 * t341;
t149 = t245 * t196 + t214;
t141 = qJD(3) * pkin(9) + t149;
t340 = qJD(1) * t246;
t315 = t235 * t340;
t150 = qJD(2) * t199 - t315;
t89 = t141 * t244 + t150 * t240;
t68 = -qJ(5) * t191 + t89;
t372 = t236 * t68;
t87 = -t141 * t240 + t244 * t150;
t67 = -qJ(5) * t193 + t87;
t49 = -pkin(4) * t216 + t67;
t24 = t234 * t49 + t372;
t401 = pkin(10) * t278;
t19 = t24 + t401;
t228 = t245 * qJDD(2);
t183 = t241 * t324 + qJDD(4) - t228;
t325 = qJD(1) * qJD(2);
t155 = qJDD(2) * pkin(8) + (qJDD(1) * t242 + t246 * t325) * t235;
t323 = qJDD(1) * t237;
t301 = t241 * t323;
t355 = t237 * t245;
t397 = qJD(1) * t355 - t241 * t196;
t75 = qJDD(3) * pkin(9) + qJD(3) * t397 + t155 * t245 + t301;
t305 = t242 * t325;
t356 = t235 * t246;
t284 = -qJDD(1) * t356 + t235 * t305;
t98 = qJD(2) * t195 + qJDD(2) * t199 + t284;
t97 = t244 * t98;
t253 = -qJD(4) * t89 - t240 * t75 + t97;
t12 = pkin(4) * t183 - qJ(5) * t112 - qJD(5) * t193 + t253;
t321 = t150 * t329 + t240 * t98 + t244 * t75;
t331 = qJD(4) * t240;
t265 = t141 * t331 - t321;
t14 = -qJ(5) * t113 - qJD(5) * t191 - t265;
t4 = t236 * t12 - t14 * t234;
t2 = pkin(5) * t183 - pkin(10) * t46 + t4;
t306 = -t19 * t327 + t239 * t2;
t140 = -qJD(3) * pkin(3) - t397;
t106 = pkin(4) * t191 + qJD(5) + t140;
t51 = -pkin(5) * t278 + t106;
t416 = -t51 * t58 - g(1) * (-t130 * t221 - t168 * t220) - g(2) * (-t128 * t221 - t166 * t220) - g(3) * (-t174 * t221 + t220 * t356) - t306;
t374 = t209 * t402;
t9 = qJD(6) * t402 + t239 * t46 - t243 * t45;
t414 = -t9 - t374;
t351 = t244 * t245;
t218 = pkin(8) * t351;
t388 = pkin(4) * t241;
t275 = -qJ(5) * t351 + t388;
t328 = qJD(5) * t244;
t413 = -t241 * t328 + t275 * qJD(3) + (-t218 + (qJ(5) * t241 - t199) * t240) * qJD(4) - t419;
t353 = t241 * t244;
t412 = -(-pkin(8) * qJD(3) - qJ(5) * qJD(4)) * t353 - (-qJD(5) * t241 + (-pkin(8) * qJD(4) - qJ(5) * qJD(3)) * t245) * t240 - t418;
t238 = -qJ(5) - pkin(9);
t297 = qJD(4) * t238;
t310 = t240 * t337;
t194 = t287 * qJD(2);
t348 = t240 * t194 + t244 * t397;
t411 = qJ(5) * t310 + t240 * t297 + t328 - t348;
t177 = t244 * t194;
t410 = -qJD(2) * t275 + t244 * t297 - t177 + (-qJD(5) + t397) * t240;
t286 = g(1) * t168 + g(2) * t166;
t409 = -g(3) * t356 + t286;
t5 = t234 * t12 + t236 * t14;
t3 = pkin(10) * t45 + t5;
t316 = t243 * t2 - t239 * t3;
t408 = -t51 * t402 - g(1) * (-t130 * t220 + t168 * t221) - g(2) * (-t128 * t220 + t166 * t221) - g(3) * (-t174 * t220 - t221 * t356) + t316;
t407 = pkin(10) * t120;
t185 = t234 * t244 + t236 * t240;
t264 = t185 * t245;
t405 = qJD(2) * t264 - t185 * qJD(4);
t277 = t234 * t240 - t236 * t244;
t404 = t216 * t277;
t285 = g(1) * t169 + g(2) * t167;
t259 = -g(3) * t357 - t285;
t379 = t412 * t234 + t413 * t236;
t378 = t413 * t234 - t412 * t236;
t370 = -t234 * t411 + t236 * t410;
t369 = t234 * t410 + t236 * t411;
t62 = t234 * t68;
t23 = t236 * t49 - t62;
t16 = -pkin(5) * t216 + t23 + t407;
t398 = (qJD(6) * t16 + t3) * t243;
t133 = -t174 * t240 - t244 * t356;
t396 = -g(1) * (-t130 * t240 + t168 * t244) - g(2) * (-t128 * t240 + t166 * t244) - g(3) * t133;
t395 = -t149 + (-t310 + t331) * pkin(4);
t394 = pkin(4) * t113 + qJDD(5);
t247 = qJD(3) ^ 2;
t393 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t247 + t235 * (-g(3) * t246 + t305) - t284 + t286;
t391 = (pkin(8) * t216 + t141) * qJD(4) + t409;
t389 = pkin(4) * t234;
t307 = t245 * t326;
t333 = qJD(3) * t245;
t308 = t240 * t333;
t103 = t185 * t330 + t234 * t308 - t236 * t307;
t381 = -pkin(5) * t334 - pkin(10) * t103 - t379;
t102 = -qJD(3) * t264 + t277 * t330;
t380 = pkin(10) * t102 + t378;
t279 = -t185 * t239 - t243 * t277;
t377 = qJD(6) * t279 + t239 * t405 + t243 * t404;
t119 = t185 * t243 - t239 * t277;
t376 = qJD(6) * t119 + t239 * t404 - t243 * t405;
t30 = t236 * t67 - t62;
t375 = qJD(2) * pkin(2);
t371 = t243 * t16;
t368 = -pkin(5) * t405 + t395;
t364 = qJDD(3) * pkin(3);
t363 = t112 * t240;
t362 = t191 * t216;
t361 = t193 * t216;
t360 = t193 * t244;
t359 = t220 * t245;
t358 = t221 * t245;
t354 = t240 * t241;
t349 = qJDD(1) - g(3);
t187 = t244 * t199;
t124 = -qJ(5) * t353 + t187 + (-pkin(4) - t387) * t245;
t345 = t240 * t199 + t218;
t135 = -qJ(5) * t354 + t345;
t70 = t234 * t124 + t236 * t135;
t200 = t238 * t240;
t201 = t238 * t244;
t138 = t234 * t200 - t236 * t201;
t344 = pkin(4) * t354 + t241 * pkin(8);
t232 = t241 ^ 2;
t343 = -t245 ^ 2 + t232;
t339 = qJD(2) * t235;
t336 = qJD(3) * t191;
t332 = qJD(4) * t216;
t318 = pkin(4) * t308 + pkin(8) * t333 + t329 * t388;
t224 = pkin(4) * t244 + pkin(3);
t313 = t241 * t340;
t312 = t242 * t339;
t311 = t246 * t339;
t309 = t216 * t326;
t303 = t246 * t324;
t29 = -t234 * t67 - t372;
t69 = t236 * t124 - t135 * t234;
t137 = t236 * t200 + t201 * t234;
t107 = -pkin(10) * t185 + t137;
t289 = pkin(10) * t405 + qJD(6) * t107 + t369;
t108 = -pkin(10) * t277 + t138;
t288 = pkin(5) * t338 + pkin(10) * t404 + qJD(6) * t108 - t370;
t7 = t239 * t16 + t243 * t19;
t162 = t277 * t241;
t38 = -pkin(5) * t245 + pkin(10) * t162 + t69;
t161 = t185 * t241;
t39 = -pkin(10) * t161 + t70;
t283 = t239 * t38 + t243 * t39;
t269 = -t174 * t244 + t240 * t356;
t72 = t133 * t236 + t234 * t269;
t73 = t133 * t234 - t236 * t269;
t282 = -t239 * t73 + t243 * t72;
t281 = t239 * t72 + t243 * t73;
t280 = -t243 * t161 + t162 * t239;
t101 = -t161 * t239 - t162 * t243;
t248 = qJD(2) ^ 2;
t276 = qJDD(2) * t246 - t242 * t248;
t222 = pkin(4) * t236 + pkin(5);
t274 = t222 * t239 + t243 * t389;
t273 = t222 * t243 - t239 * t389;
t173 = t241 * t357 - t355;
t267 = t240 * t183 - t216 * t329;
t266 = t244 * t183 + t216 * t331;
t127 = t167 * t241 + t245 * t296;
t129 = t169 * t241 - t245 * t295;
t263 = g(1) * t129 + g(2) * t127 + g(3) * t173;
t262 = g(1) * t130 + g(2) * t128 + g(3) * t174;
t261 = qJD(3) * t214 + t241 * t155 + t196 * t333 - t245 * t323;
t258 = -pkin(9) * t183 - t140 * t216;
t76 = t261 - t364;
t252 = pkin(9) * t332 + t263 - t76;
t197 = -t315 - t375;
t251 = -pkin(8) * qJDD(3) + (t197 + t315 - t375) * qJD(3);
t37 = t76 + t394;
t250 = -t261 + t263;
t179 = qJDD(6) + t183;
t156 = pkin(5) * t277 - t224;
t132 = qJD(3) * t174 + t241 * t311;
t131 = -qJD(3) * t173 + t245 * t311;
t125 = pkin(5) * t161 + t344;
t93 = pkin(4) * t193 - pkin(5) * t120;
t71 = -pkin(5) * t102 + t318;
t61 = qJD(4) * t133 + t131 * t244 + t240 * t312;
t60 = qJD(4) * t269 - t131 * t240 + t244 * t312;
t34 = qJD(6) * t101 - t243 * t102 - t103 * t239;
t33 = qJD(6) * t280 + t102 * t239 - t103 * t243;
t28 = t234 * t60 + t236 * t61;
t26 = -t234 * t61 + t236 * t60;
t22 = -pkin(5) * t45 + t37;
t21 = t30 + t407;
t20 = t29 - t401;
t6 = -t19 * t239 + t371;
t1 = [t349, 0, t276 * t235 (-qJDD(2) * t242 - t246 * t248) * t235, 0, 0, 0, 0, 0, -qJD(3) * t132 - qJDD(3) * t173 + (-t241 * t303 + t245 * t276) * t235, -qJD(3) * t131 - qJDD(3) * t174 + (-t241 * t276 - t245 * t303) * t235, 0, 0, 0, 0, 0, t113 * t173 + t132 * t191 + t133 * t183 - t216 * t60, t112 * t173 + t132 * t193 + t183 * t269 + t216 * t61, t120 * t26 + t278 * t28 + t45 * t73 - t46 * t72, t106 * t132 + t173 * t37 + t23 * t26 + t24 * t28 + t4 * t72 + t5 * t73 - g(3), 0, 0, 0, 0, 0 -(-qJD(6) * t281 - t239 * t28 + t243 * t26) * t209 + t282 * t179 - t132 * t58 + t173 * t9 (qJD(6) * t282 + t239 * t26 + t243 * t28) * t209 - t281 * t179 + t132 * t402 + t173 * t8; 0, qJDD(2), t349 * t356 + t286, -t349 * t357 + t285, qJDD(2) * t232 + 0.2e1 * t241 * t304, 0.2e1 * t228 * t241 - 0.2e1 * t324 * t343, qJDD(3) * t241 + t245 * t247, qJDD(3) * t245 - t241 * t247, 0, t251 * t241 + t393 * t245, -t393 * t241 + t251 * t245, t112 * t353 + (-t240 * t330 + t307) * t193 (-t191 * t244 - t193 * t240) * t333 + (-t363 - t113 * t244 + (t191 * t240 - t360) * qJD(4)) * t241 (-t112 - t309) * t245 + (qJD(3) * t193 + t266) * t241 (t216 * t335 + t113) * t245 + (-t267 - t336) * t241, -t183 * t245 - t216 * t334, t187 * t183 + t419 * t216 + (t199 * t332 + t259) * t240 + (pkin(8) * t336 - t97 + (-pkin(8) * t183 + qJD(3) * t140 + qJD(4) * t150 + t75) * t240 + t391 * t244) * t245 + (pkin(8) * t113 + qJD(3) * t87 + t140 * t329 - t191 * t315 + t76 * t240) * t241, -t345 * t183 + t418 * t216 + t259 * t244 + ((pkin(8) * t193 + t140 * t244) * qJD(3) - t391 * t240 + t321) * t245 + (-t193 * t315 - t140 * t331 - t89 * qJD(3) + t76 * t244 + (t112 - t309) * pkin(8)) * t241, t102 * t24 + t103 * t23 + t120 * t379 - t161 * t5 + t162 * t4 + t241 * t409 + t278 * t378 + t45 * t70 - t46 * t69, t5 * t70 + t4 * t69 + t37 * t344 + t378 * t24 + t379 * t23 + t259 * (pkin(4) * t240 + pkin(8)) + (-t341 * t356 + t318) * t106 + t409 * (t224 * t245 - t238 * t241 + pkin(2)) t101 * t8 + t33 * t402, -t101 * t9 + t280 * t8 + t33 * t58 - t34 * t402, t101 * t179 - t209 * t33 - t245 * t8 + t334 * t402, t179 * t280 + t209 * t34 + t245 * t9 + t334 * t58, -t179 * t245 - t209 * t334 (-t239 * t39 + t243 * t38) * t179 - t316 * t245 + t6 * t334 - t71 * t58 + t125 * t9 - t22 * t280 + t51 * t34 - g(1) * (-t168 * t358 + t169 * t220) - g(2) * (-t166 * t358 + t167 * t220) + (t239 * t380 + t243 * t381) * t209 + (t209 * t283 + t245 * t7) * qJD(6) + (t58 * t313 - g(3) * (t220 * t242 + t221 * t350)) * t235, -t283 * t179 + (t306 + t398) * t245 - t7 * t334 + t71 * t402 + t125 * t8 + t22 * t101 + t51 * t33 - g(1) * (t168 * t359 + t169 * t221) - g(2) * (t166 * t359 + t167 * t221) + ((qJD(6) * t38 + t380) * t243 + (-qJD(6) * t39 - t381) * t239) * t209 + (-t402 * t313 - g(3) * (-t220 * t350 + t221 * t242)) * t235; 0, 0, 0, 0, -t241 * t248 * t245, t343 * t248, t322, t228, qJDD(3), qJD(3) * t149 - t197 * t338 + t250, -t301 + (-qJD(2) * t197 - t155) * t245 + t262, -t216 * t360 + t363 (t112 + t362) * t244 + (-t113 + t361) * t240 (-t193 * t241 + t216 * t351) * qJD(2) + t267 (-t216 * t240 * t245 + t191 * t241) * qJD(2) + t266, t216 * t338, -t87 * t338 - pkin(3) * t113 - t149 * t191 + t177 * t216 + (-t216 * t397 + t258) * t240 + t252 * t244, -pkin(3) * t112 - t149 * t193 - t216 * t348 - t240 * t252 + t244 * t258 + t338 * t89, t370 * t120 - t137 * t46 + t138 * t45 - t185 * t4 - t23 * t404 + t24 * t405 - t277 * t5 + t369 * t278 - t262, t5 * t138 + t4 * t137 - t37 * t224 - g(1) * (-t129 * t224 - t130 * t238) - g(2) * (-t127 * t224 - t128 * t238) - g(3) * (-t173 * t224 - t174 * t238) + t369 * t24 + t370 * t23 + t395 * t106, t119 * t8 + t377 * t402, -t119 * t9 + t279 * t8 - t376 * t402 + t377 * t58, t119 * t179 - t209 * t377 - t338 * t402, t179 * t279 + t209 * t376 - t338 * t58, t209 * t338 (t107 * t243 - t108 * t239) * t179 + t156 * t9 - t22 * t279 - t6 * t338 - t368 * t58 + t376 * t51 + (t239 * t289 + t243 * t288) * t209 + t263 * t221 -(t107 * t239 + t108 * t243) * t179 + t156 * t8 + t22 * t119 + t7 * t338 + t368 * t402 + t377 * t51 + (-t239 * t288 + t243 * t289) * t209 - t263 * t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193 * t191, -t191 ^ 2 + t193 ^ 2, t112 - t362, -t113 - t361, t183, -t140 * t193 - t89 * t216 + t253 + t396, -t87 * t216 + t140 * t191 - g(1) * (-t130 * t244 - t168 * t240) - g(2) * (-t128 * t244 - t166 * t240) - g(3) * t269 + t265 (t234 * t45 - t236 * t46) * pkin(4) + (-t30 + t23) * t278 + (-t24 - t29) * t120, -t23 * t29 - t24 * t30 + (-t106 * t193 + t5 * t234 + t4 * t236 + t396) * pkin(4), -t420, t417, t421, t414, t179, t273 * t179 + (t20 * t243 - t21 * t239) * t209 + t93 * t58 + (t209 * t274 - t7) * qJD(6) + t408, -t274 * t179 - t243 * t3 - (t20 * t239 + t21 * t243) * t209 - t93 * t402 + (t209 * t273 - t371) * qJD(6) + t416; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120 ^ 2 - t278 ^ 2, -t120 * t23 - t24 * t278 - t250 - t364 + t394, 0, 0, 0, 0, 0, t9 - t374, t8 - t373; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t420, t417, t421, t414, t179 (-qJD(6) - t209) * t7 + t408, -t209 * t6 - t398 + t416;];
tau_reg  = t1;
