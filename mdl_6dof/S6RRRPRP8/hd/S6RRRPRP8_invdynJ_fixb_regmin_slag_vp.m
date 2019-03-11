% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRP8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:19:43
% EndTime: 2019-03-09 17:20:01
% DurationCPUTime: 7.99s
% Computational Cost: add. (5892->591), mult. (12660->735), div. (0->0), fcn. (8578->8), ass. (0->268)
t242 = sin(qJ(3));
t246 = cos(qJ(3));
t248 = cos(qJ(1));
t244 = sin(qJ(1));
t247 = cos(qJ(2));
t358 = t244 * t247;
t125 = t242 * t358 + t246 * t248;
t354 = t248 * t242;
t356 = t246 * t247;
t126 = t244 * t356 - t354;
t241 = sin(qJ(5));
t245 = cos(qJ(5));
t416 = t125 * t245 - t126 * t241;
t127 = -t244 * t246 + t247 * t354;
t355 = t247 * t248;
t128 = t242 * t244 + t246 * t355;
t72 = t127 * t245 - t128 * t241;
t439 = -g(1) * t72 - g(2) * t416;
t243 = sin(qJ(2));
t343 = -t247 * pkin(2) - t243 * pkin(8);
t427 = -pkin(1) + t343;
t136 = t427 * qJD(1);
t340 = qJD(1) * t247;
t220 = pkin(7) * t340;
t175 = qJD(2) * pkin(8) + t220;
t88 = t246 * t136 - t242 * t175;
t353 = qJD(4) - t88;
t333 = qJD(3) * t246;
t438 = t246 * t340 - t333;
t322 = t242 * t340;
t335 = qJD(3) * t242;
t437 = -t322 + t335;
t361 = t243 * t246;
t364 = t242 * t245;
t114 = t241 * t361 - t243 * t364;
t195 = -qJD(3) + t340;
t405 = pkin(3) + pkin(4);
t341 = qJD(1) * t243;
t321 = t246 * t341;
t338 = qJD(2) * t242;
t151 = t321 + t338;
t433 = t151 * pkin(9) - t353;
t38 = t195 * t405 - t433;
t180 = t195 * qJ(4);
t330 = t246 * qJD(2);
t149 = t242 * t341 - t330;
t89 = t242 * t136 + t246 * t175;
t54 = pkin(9) * t149 + t89;
t44 = -t180 + t54;
t18 = t241 * t38 + t245 * t44;
t223 = t247 * qJDD(1);
t329 = qJD(1) * qJD(2);
t418 = -t243 * t329 + t223;
t147 = qJDD(3) - t418;
t122 = pkin(7) * t418 + qJDD(2) * pkin(8);
t292 = pkin(2) * t243 - pkin(8) * t247;
t164 = t292 * qJD(2);
t96 = qJD(1) * t164 + qJDD(1) * t427;
t297 = t242 * t122 + t136 * t335 + t175 * t333 - t246 * t96;
t283 = qJDD(4) + t297;
t310 = t247 * t329;
t328 = t243 * qJDD(1);
t334 = qJD(3) * t243;
t429 = qJD(1) * t334 - qJDD(2);
t74 = -qJD(3) * t330 + (-t310 - t328) * t246 + t429 * t242;
t11 = pkin(9) * t74 - t147 * t405 + t283;
t134 = t147 * qJ(4);
t178 = t195 * qJD(4);
t266 = t246 * t122 + t136 * t333 - t175 * t335 + t242 * t96;
t23 = t134 - t178 + t266;
t75 = t242 * (qJD(2) * (qJD(3) + t340) + t328) + t429 * t246;
t13 = pkin(9) * t75 + t23;
t307 = t245 * t11 - t241 * t13;
t254 = -t18 * qJD(5) + t307;
t174 = -qJD(2) * pkin(2) + pkin(7) * t341;
t70 = t149 * pkin(3) - t151 * qJ(4) + t174;
t46 = -pkin(4) * t149 - t70;
t83 = t149 * t241 + t151 * t245;
t436 = g(3) * t114 - t46 * t83 + t254 + t439;
t161 = t292 * qJD(1);
t130 = t242 * t161;
t351 = qJ(4) * t341 + t130;
t363 = t242 * t247;
t404 = pkin(8) - pkin(9);
t431 = t404 * t335 + (-pkin(7) * t361 + pkin(9) * t363) * qJD(1) + t351;
t177 = t404 * t246;
t325 = -pkin(7) * t242 - pkin(3);
t257 = -pkin(9) * t356 + (-pkin(4) + t325) * t243;
t357 = t246 * t161;
t435 = -qJD(1) * t257 + qJD(3) * t177 + t357;
t407 = t83 ^ 2;
t281 = -t245 * t149 + t151 * t241;
t78 = t281 ^ 2;
t434 = t78 - t407;
t432 = qJ(6) * t83;
t331 = qJD(5) * t245;
t332 = qJD(5) * t241;
t378 = t241 * t438 - t242 * t331 + t245 * t437 + t246 * t332;
t153 = t241 * t242 + t245 * t246;
t265 = t153 * t247;
t91 = t153 * qJD(5) - t241 * t335 - t245 * t333;
t375 = -qJD(1) * t265 - t91;
t360 = t243 * t248;
t362 = t243 * t244;
t430 = g(1) * t360 + g(2) * t362;
t184 = qJD(5) + t195;
t21 = qJD(5) * t83 - t241 * t74 - t245 * t75;
t423 = t184 * t83 - t21;
t426 = t83 * t281;
t425 = t435 * t245;
t20 = -t149 * t331 + t151 * t332 - t241 * t75 + t245 * t74;
t424 = -t184 * t281 + t20;
t422 = qJ(6) * t281;
t289 = g(1) * t248 + g(2) * t244;
t411 = t243 * t289;
t302 = -qJ(4) * t241 - t245 * t405;
t421 = -qJD(5) * t302 + t241 * t54 + t245 * t433;
t165 = qJ(4) * t245 - t241 * t405;
t420 = -qJD(5) * t165 + t241 * t433 - t245 * t54;
t176 = t404 * t242;
t349 = t241 * t176 + t245 * t177;
t419 = -t176 * t331 + t177 * t332 - t241 * t435 + t245 * t431;
t417 = -qJ(4) * t438 + t242 * qJD(4) + t220;
t226 = t242 * qJ(4);
t414 = t246 * pkin(3) + pkin(2) + t226;
t392 = g(3) * t247;
t413 = -t392 + t430;
t398 = pkin(8) * t147;
t412 = t195 * t70 + t398;
t294 = -t242 * t405 - pkin(7);
t336 = qJD(2) * t247;
t316 = t247 * t330;
t347 = qJ(4) * t316 + qJD(4) * t361;
t39 = (-t246 * t405 - t226) * t334 + t294 * t336 + t347;
t115 = t153 * t243;
t274 = -t241 * t11 - t245 * t13 - t38 * t331 + t332 * t44;
t373 = t125 * t241;
t282 = t126 * t245 + t373;
t371 = t127 * t241;
t73 = t128 * t245 + t371;
t409 = -g(1) * t73 - g(2) * t282 - g(3) * t115 - t46 * t281 - t274;
t408 = -0.2e1 * pkin(1);
t406 = t151 ^ 2;
t17 = -t241 * t44 + t245 * t38;
t7 = t17 - t432;
t6 = pkin(5) * t184 + t7;
t403 = -t7 + t6;
t400 = pkin(3) * t147;
t399 = pkin(5) * t241;
t397 = g(1) * t244;
t394 = g(2) * t248;
t235 = g(3) * t243;
t216 = pkin(5) * t245 + pkin(4);
t391 = -pkin(3) - t216;
t390 = qJ(6) * t378 - qJD(6) * t153 - t419;
t366 = t241 * t246;
t280 = -t364 + t366;
t389 = pkin(5) * t341 - qJ(6) * t375 - qJD(5) * t349 + t280 * qJD(6) + t241 * t431 + t425;
t197 = pkin(7) * t363;
t231 = t247 * pkin(3);
t76 = t247 * pkin(4) + t197 + t231 + (-pkin(9) * t243 - t427) * t246;
t365 = t242 * t243;
t199 = pkin(7) * t356;
t346 = t242 * t427 + t199;
t98 = -qJ(4) * t247 + t346;
t87 = pkin(9) * t365 + t98;
t386 = t241 * t76 + t245 * t87;
t385 = pkin(8) * qJD(3);
t60 = -t180 + t89;
t382 = t195 * t60;
t381 = t195 * t89;
t380 = t74 * t242;
t323 = qJD(3) * t405;
t379 = -t242 * t323 + t322 * t405 + t417;
t377 = -t421 - t432;
t376 = t420 + t422;
t370 = t149 * t195;
t369 = t151 * t149;
t368 = t151 * t195;
t367 = t151 * t246;
t251 = qJD(1) ^ 2;
t359 = t243 * t251;
t352 = pkin(3) * t437 - t417;
t350 = t242 * t164 + t333 * t427;
t95 = t151 * pkin(3) + t149 * qJ(4);
t348 = t430 * t246;
t236 = t243 ^ 2;
t342 = -t247 ^ 2 + t236;
t339 = qJD(2) * t151;
t337 = qJD(2) * t243;
t327 = qJ(4) * t337 + t350;
t217 = pkin(7) * t328;
t123 = -qJDD(2) * pkin(2) + pkin(7) * t310 + t217;
t326 = g(1) * t355 + g(2) * t358 + t235;
t324 = pkin(3) * t242 + pkin(7);
t319 = t195 * t338;
t318 = t195 * t330;
t317 = t242 * t336;
t315 = t195 * t335;
t314 = t242 * t334;
t308 = qJ(4) + t399;
t287 = -qJD(3) * t199 + t246 * t164 - t335 * t427;
t34 = pkin(9) * t314 + qJD(2) * t257 - t287;
t35 = (-pkin(7) * qJD(2) + pkin(9) * qJD(3)) * t361 + (-qJD(4) + (-pkin(7) * qJD(3) + pkin(9) * qJD(2)) * t242) * t247 + t327;
t306 = -t241 * t35 + t245 * t34;
t303 = -t241 * t87 + t245 * t76;
t300 = t245 * t176 - t177 * t241;
t299 = t246 * t427 - t197;
t298 = t184 ^ 2;
t142 = t246 * pkin(4) + t414;
t295 = pkin(3) * t356 + qJ(4) * t363 - t343;
t55 = -pkin(4) * t151 - t95;
t293 = t325 * t243;
t291 = -g(1) * t125 + g(2) * t127;
t290 = g(1) * t126 - g(2) * t128;
t288 = -t126 * pkin(3) + t248 * pkin(7) - t125 * qJ(4);
t58 = pkin(3) * t195 + t353;
t285 = -t242 * t60 + t246 * t58;
t284 = qJD(3) * t174 - t398;
t277 = -g(1) * t127 - g(2) * t125 - g(3) * t365;
t24 = t75 * pkin(3) + t74 * qJ(4) - t151 * qJD(4) + t123;
t276 = -t195 * t385 + t392;
t273 = -pkin(7) * qJDD(2) + t329 * t408;
t272 = t241 * t34 + t245 * t35 + t76 * t331 - t332 * t87;
t271 = t242 * t147 - t195 * t333;
t270 = t246 * t147 + t315;
t268 = g(3) * t280;
t267 = -t24 - t276;
t250 = qJD(2) ^ 2;
t261 = pkin(7) * t250 + qJDD(1) * t408 - t397;
t260 = t248 * pkin(1) + pkin(2) * t355 + t128 * pkin(3) + t244 * pkin(7) + pkin(8) * t360 + qJ(4) * t127;
t14 = -pkin(4) * t75 - t24;
t256 = -t277 - t297;
t5 = pkin(5) * t21 + qJDD(6) + t14;
t253 = t151 * t70 + qJDD(4) - t256;
t252 = g(1) * t128 + g(2) * t126 + g(3) * t361 - t195 * t88 - t266;
t240 = -qJ(6) - pkin(9);
t207 = g(2) * t360;
t202 = pkin(8) * t355;
t198 = pkin(8) * t358;
t192 = qJ(4) * t361;
t160 = -pkin(5) + t302;
t135 = -qJDD(5) + t147;
t120 = t127 * pkin(3);
t118 = t125 * pkin(3);
t108 = t243 * t324 - t192;
t99 = t231 - t299;
t97 = t243 * t294 + t192;
t94 = qJD(1) * t293 - t357;
t93 = -pkin(7) * t321 + t351;
t62 = -qJ(6) * t153 + t349;
t61 = qJ(6) * t280 + t300;
t50 = qJ(4) * t314 + pkin(7) * t336 + (t243 * t333 + t317) * pkin(3) - t347;
t45 = qJD(2) * t293 - t287;
t43 = -t74 - t370;
t42 = -t247 * qJD(4) + (-t243 * t330 - t247 * t335) * pkin(7) + t327;
t41 = qJD(2) * t265 + (qJD(3) - qJD(5)) * t243 * t280;
t40 = t241 * t316 + t243 * t91 - t245 * t317;
t30 = pkin(5) * t281 + qJD(6) + t46;
t29 = -qJ(6) * t114 + t386;
t28 = pkin(5) * t247 - qJ(6) * t115 + t303;
t27 = t283 - t400;
t8 = t18 - t422;
t4 = -qJ(6) * t40 - qJD(6) * t114 + t272;
t3 = -pkin(5) * t337 - qJ(6) * t41 - qJD(5) * t386 - qJD(6) * t115 + t306;
t2 = -qJ(6) * t21 - qJD(6) * t281 - t274;
t1 = -pkin(5) * t135 + qJ(6) * t20 - qJD(6) * t83 + t254;
t9 = [qJDD(1), -t394 + t397, t289, qJDD(1) * t236 + 0.2e1 * t243 * t310, 0.2e1 * t223 * t243 - 0.2e1 * t329 * t342, qJDD(2) * t243 + t247 * t250, qJDD(2) * t247 - t243 * t250, 0, t273 * t243 + (-t261 - t394) * t247, t243 * t261 + t247 * t273 + t207, -t74 * t361 + (-t314 + t316) * t151 (-t149 * t246 - t151 * t242) * t336 + (t380 - t246 * t75 + (t149 * t242 - t367) * qJD(3)) * t243 (t74 - t318) * t247 + (t270 + t339) * t243 (t75 + t319) * t247 + (-qJD(2) * t149 - t271) * t243, -t147 * t247 - t195 * t337, -t287 * t195 + t299 * t147 + ((pkin(7) * t149 + t174 * t242) * qJD(2) + t297) * t247 + (t174 * t333 + t88 * qJD(2) + t123 * t242 + (t75 - t319) * pkin(7)) * t243 + t290, t350 * t195 - t346 * t147 + (t174 * t330 + (-t315 + t339) * pkin(7) + t266) * t247 + (-t174 * t335 - t89 * qJD(2) + t123 * t246 + (-t74 - t318) * pkin(7)) * t243 + t291, t108 * t75 - t147 * t99 + t149 * t50 + t195 * t45 + (t338 * t70 + t27) * t247 + (-qJD(2) * t58 + t24 * t242 + t333 * t70) * t243 + t290, -t149 * t42 + t151 * t45 - t74 * t99 - t75 * t98 - t207 + t285 * t336 + (t397 - t23 * t242 + t246 * t27 + (-t242 * t58 - t246 * t60) * qJD(3)) * t243, t108 * t74 + t147 * t98 - t151 * t50 - t195 * t42 + (-t330 * t70 - t23) * t247 + (qJD(2) * t60 - t24 * t246 + t335 * t70) * t243 - t291, -g(1) * t288 - g(2) * t260 + t24 * t108 + t23 * t98 + t27 * t99 - t397 * t427 + t60 * t42 + t58 * t45 + t70 * t50, -t115 * t20 + t41 * t83, t114 * t20 - t115 * t21 - t281 * t41 - t40 * t83, -t115 * t135 + t184 * t41 - t20 * t247 - t337 * t83, t114 * t135 - t184 * t40 - t21 * t247 + t281 * t337, -t135 * t247 - t184 * t337, t306 * t184 - t303 * t135 + t307 * t247 - t17 * t337 + t39 * t281 + t97 * t21 + t14 * t114 + t46 * t40 + g(1) * t282 - g(2) * t73 + (-t18 * t247 - t184 * t386) * qJD(5), g(1) * t416 - g(2) * t72 + t14 * t115 + t135 * t386 + t18 * t337 - t184 * t272 - t97 * t20 + t247 * t274 + t39 * t83 + t46 * t41, -g(1) * t362 - t1 * t115 - t114 * t2 + t20 * t28 - t21 * t29 - t281 * t4 - t3 * t83 - t40 * t8 - t41 * t6 + t207, t2 * t29 + t8 * t4 + t1 * t28 + t6 * t3 + t5 * (t114 * pkin(5) + t192) - g(1) * (-pkin(1) * t244 - pkin(2) * t358 - pkin(5) * t373 - t126 * t216 + t288) - g(2) * (pkin(5) * t371 + t128 * t216 + t260) + (t5 * (-pkin(4) * t242 - t324) - t240 * t394 - (-pkin(8) - t240) * t397) * t243 + (pkin(5) * t40 + t39) * t30; 0, 0, 0, -t247 * t359, t342 * t251, t328, t223, qJDD(2), pkin(1) * t359 - t217 + t413 (pkin(1) * t251 - pkin(7) * qJDD(1)) * t247 + t326, -t195 * t367 - t380 (-t74 + t370) * t246 + (-t75 + t368) * t242 (-t151 * t243 + t195 * t356) * qJD(1) + t271 (t149 * t243 - t195 * t363) * qJD(1) + t270, t195 * t341, -pkin(2) * t75 + t284 * t242 + (-t392 - t123 + (t161 + t385) * t195) * t246 + (-t174 * t363 - t88 * t243 + (-t149 * t247 + t195 * t365) * pkin(7)) * qJD(1) + t348, pkin(2) * t74 - t130 * t195 + t284 * t246 + (-t174 * t356 + t89 * t243 + (-t151 * t247 + t195 * t361) * pkin(7)) * qJD(1) + (t123 + t276 - t411) * t242, t352 * t149 - t195 * t94 - t242 * t412 + t267 * t246 + t58 * t341 - t414 * t75 + t348, t149 * t93 - t151 * t94 + (t23 - t195 * t58 + (qJD(3) * t151 - t75) * pkin(8)) * t246 + (t27 + t382 + (qJD(3) * t149 - t74) * pkin(8)) * t242 - t326, -t60 * t341 - t414 * t74 + t195 * t93 - t352 * t151 + t412 * t246 + (t267 + t411) * t242, -t60 * t93 - t58 * t94 - g(1) * t202 - g(2) * t198 - g(3) * t295 + t352 * t70 + (qJD(3) * t285 + t23 * t246 + t27 * t242) * pkin(8) + (-t24 + t411) * t414, t20 * t280 + t375 * t83, t20 * t153 + t21 * t280 - t281 * t375 + t378 * t83, t135 * t280 + t184 * t375 + t341 * t83, t135 * t153 + t184 * t378 - t281 * t341, t184 * t341, -t300 * t135 + t142 * t21 + t379 * t281 - t378 * t46 - g(3) * t265 + (-t177 * t331 + (-qJD(5) * t176 + t431) * t241 + t425) * t184 + t17 * t341 + (t14 + t411) * t153, t349 * t135 - t142 * t20 - t14 * t280 + t379 * t83 + t375 * t46 + t247 * t268 + t419 * t184 + (-t18 * qJD(1) - t280 * t289) * t243, t1 * t280 - t153 * t2 + t20 * t61 - t21 * t62 - t281 * t390 - t375 * t6 + t378 * t8 - t389 * t83 + t326, t2 * t62 + t1 * t61 + t5 * (pkin(5) * t153 + t142) - g(1) * (t240 * t355 + t202) - g(2) * (t240 * t358 + t198) - g(3) * (t216 * t356 + t295) + t390 * t8 + t389 * t6 + (-pkin(5) * t378 + t417) * t30 + (-t30 * t323 + (qJD(1) * t30 * t405 - g(3) * t399) * t247) * t242 + (-g(3) * t240 + t289 * (t242 * t308 - t246 * t391 + pkin(2))) * t243; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t369, -t149 ^ 2 + t406, t43, -t368 - t75, t147, -t151 * t174 + t256 - t381, t149 * t174 + t252, -t149 * t95 - t253 - t381 + 0.2e1 * t400, pkin(3) * t74 - t75 * qJ(4) + (t60 - t89) * t151 + (t58 - t353) * t149, -t149 * t70 + t151 * t95 + 0.2e1 * t134 - 0.2e1 * t178 - t252, t23 * qJ(4) - t27 * pkin(3) - t70 * t95 - t58 * t89 - g(1) * (qJ(4) * t128 - t120) - g(2) * (qJ(4) * t126 - t118) - g(3) * (-pkin(3) * t365 + t192) + t353 * t60, -t426, t434, t424, -t423, t135, -t302 * t135 + t184 * t420 - t55 * t281 - t436, t165 * t135 + t184 * t421 - t55 * t83 + t409, t160 * t20 - t165 * t21 + (-t377 + t6) * t281 + (-t376 - t8) * t83, t2 * t165 + t1 * t160 - t30 * (-pkin(5) * t83 + t55) - g(1) * (-t127 * t216 + t128 * t308 - t120) - g(2) * (-t125 * t216 + t126 * t308 - t118) - g(3) * t192 + t377 * t8 + t376 * t6 - (pkin(5) * t366 + t242 * t391) * t235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147 + t369, t43, -t195 ^ 2 - t406, t253 + t382 - t400, 0, 0, 0, 0, 0, -t245 * t135 - t151 * t281 - t241 * t298, t241 * t135 - t151 * t83 - t245 * t298, t241 * t423 + t245 * t424, -t151 * t30 + (t184 * t8 + t1) * t245 + (-t184 * t6 + t2) * t241 + t277; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t426, -t434, -t424, t423, -t135, t18 * t184 + t436, t17 * t184 - t409, pkin(5) * t20 - t281 * t403, t403 * t8 + (t243 * t268 - t30 * t83 + t1 + t439) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78 - t407, t281 * t8 + t6 * t83 + t413 + t5;];
tau_reg  = t9;
