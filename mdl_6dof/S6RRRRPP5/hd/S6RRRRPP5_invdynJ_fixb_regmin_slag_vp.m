% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRPP5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:09:12
% EndTime: 2019-03-09 21:09:29
% DurationCPUTime: 7.95s
% Computational Cost: add. (9662->661), mult. (20966->773), div. (0->0), fcn. (14561->10), ass. (0->319)
t266 = cos(qJ(3));
t264 = sin(qJ(2));
t267 = cos(qJ(2));
t399 = t266 * t267;
t315 = pkin(3) * t264 - pkin(9) * t399;
t269 = -pkin(9) - pkin(8);
t362 = qJD(3) * t269;
t331 = pkin(2) * t264 - pkin(8) * t267;
t178 = t331 * qJD(1);
t263 = sin(qJ(3));
t386 = qJD(1) * t264;
t353 = t263 * t386;
t389 = pkin(7) * t353 + t266 * t178;
t477 = -qJD(1) * t315 + t266 * t362 - t389;
t158 = t263 * t178;
t405 = t264 * t266;
t408 = t263 * t267;
t476 = -t158 - (-pkin(7) * t405 - pkin(9) * t408) * qJD(1) + t263 * t362;
t381 = qJD(2) * t266;
t172 = -t353 + t381;
t383 = qJD(2) * t263;
t173 = t266 * t386 + t383;
t262 = sin(qJ(4));
t445 = cos(qJ(4));
t308 = t262 * t172 + t173 * t445;
t101 = t308 ^ 2;
t103 = -t445 * t172 + t173 * t262;
t447 = t103 ^ 2;
t475 = t101 - t447;
t385 = qJD(1) * t267;
t360 = t263 * t385;
t361 = t445 * t266;
t411 = t262 * t263;
t452 = qJD(3) + qJD(4);
t351 = t445 * qJD(4);
t453 = t445 * qJD(3) + t351;
t393 = -t262 * t360 - t266 * t453 + t361 * t385 + t411 * t452;
t428 = qJD(2) * pkin(2);
t196 = pkin(7) * t386 - t428;
t312 = pkin(3) * t172 - t196;
t292 = qJ(5) * t308 + t312;
t47 = pkin(4) * t103 - t292;
t474 = t103 * t47;
t339 = qJD(3) + t385;
t369 = t264 * qJDD(1);
t290 = qJD(2) * t339 + t369;
t372 = qJD(1) * qJD(3);
t348 = t264 * t372;
t319 = -qJDD(2) + t348;
t302 = t319 * t263;
t275 = t290 * t266 - t302;
t313 = qJD(2) * qJD(3) + t369;
t373 = qJD(1) * qJD(2);
t349 = t267 * t373;
t468 = t313 + t349;
t338 = t263 * t468 + t266 * t348;
t299 = t266 * qJDD(2) - t338;
t376 = qJD(4) * t262;
t35 = -t172 * t351 + t173 * t376 - t262 * t299 - t445 * t275;
t226 = -qJD(3) + t385;
t211 = -qJD(4) + t226;
t419 = t103 * t211;
t289 = t35 + t419;
t198 = t269 * t263;
t199 = t269 * t266;
t473 = t198 * t351 + t199 * t376 + t262 * t477 + t476 * t445;
t420 = qJ(5) * t103;
t472 = qJ(6) * t103;
t471 = t103 * t312;
t418 = t308 * t103;
t175 = t262 * t266 + t263 * t445;
t118 = t452 * t175;
t392 = -t175 * t385 + t118;
t249 = t267 * qJDD(1);
t454 = -t264 * t373 + t249;
t169 = qJDD(3) - t454;
t163 = qJDD(4) + t169;
t457 = t163 * qJ(5) - t211 * qJD(5);
t247 = pkin(7) * t385;
t378 = qJD(3) * t263;
t332 = -t247 + (-t360 + t378) * pkin(3);
t268 = cos(qJ(1));
t404 = t264 * t268;
t265 = sin(qJ(1));
t406 = t264 * t265;
t470 = g(1) * t404 + g(2) * t406;
t36 = qJD(4) * t308 + t262 * t275 - t445 * t299;
t417 = t308 * t211;
t469 = -t36 - t417;
t380 = qJD(2) * t267;
t359 = t263 * t380;
t377 = qJD(3) * t266;
t467 = t264 * t377 + t359;
t446 = pkin(4) + pkin(5);
t23 = -t103 * t446 + qJD(6) + t292;
t464 = t36 * qJ(6) + t103 * qJD(6);
t466 = t103 * t23 + t464;
t465 = pkin(4) * t308;
t432 = -qJ(5) * t386 + t473;
t463 = t308 * t47;
t122 = t262 * t198 - t445 * t199;
t462 = qJD(4) * t122 + t476 * t262 - t445 * t477;
t461 = qJ(6) * t308;
t460 = t312 * t308;
t459 = t308 * t446;
t191 = -pkin(2) * t267 - pkin(8) * t264 - pkin(1);
t164 = t191 * qJD(1);
t197 = qJD(2) * pkin(8) + t247;
t114 = t263 * t164 + t266 * t197;
t83 = pkin(9) * t172 + t114;
t425 = t262 * t83;
t113 = t266 * t164 - t197 * t263;
t82 = -pkin(9) * t173 + t113;
t73 = -pkin(3) * t226 + t82;
t38 = t445 * t73 - t425;
t394 = qJD(5) - t38;
t458 = 0.2e1 * t457;
t456 = -t211 ^ 2 - t101;
t455 = qJ(5) * t393 - qJD(5) * t175 + t332;
t439 = g(2) * t265;
t328 = g(1) * t268 + t439;
t150 = t163 * pkin(4);
t451 = t150 - qJDD(5);
t450 = t163 * pkin(5) - t35 * qJ(6) + t308 * qJD(6);
t437 = g(3) * t267;
t449 = -t437 + t470;
t260 = qJ(3) + qJ(4);
t250 = sin(t260);
t251 = cos(t260);
t400 = t265 * t267;
t136 = t250 * t400 + t251 * t268;
t397 = t268 * t250;
t138 = -t265 * t251 + t267 * t397;
t181 = t331 * qJD(2);
t119 = qJD(1) * t181 + qJDD(1) * t191;
t109 = t266 * t119;
t146 = pkin(7) * t454 + qJDD(2) * pkin(8);
t370 = t263 * qJDD(2);
t20 = -t263 * t146 + t109 - (t370 + (t349 + t369) * t266) * pkin(9) + t169 * pkin(3) - t83 * qJD(3);
t300 = t263 * t119 + t266 * t146 + t164 * t377 - t197 * t378;
t26 = pkin(9) * t299 + t300;
t345 = -t445 * t20 + t262 * t26 + t83 * t351 + t73 * t376;
t415 = t250 * t264;
t288 = g(1) * t138 + g(2) * t136 + g(3) * t415 - t345;
t283 = t288 + t451;
t276 = -t23 * t308 - t283 - t450;
t448 = -0.2e1 * pkin(1);
t444 = pkin(7) * t263;
t442 = g(1) * t265;
t153 = t263 * t400 + t266 * t268;
t440 = g(2) * t153;
t438 = g(2) * t268;
t257 = g(3) * t264;
t252 = t264 * pkin(7);
t436 = -t392 * t446 - t455;
t174 = -t361 + t411;
t435 = qJ(6) * t392 + qJD(6) * t174 + t432;
t366 = t264 * t446;
t434 = qJ(6) * t393 + qJD(1) * t366 - t175 * qJD(6) + t462;
t433 = pkin(4) * t392 + t455;
t431 = pkin(4) * t386 + t462;
t45 = t445 * t82 - t425;
t80 = t445 * t83;
t39 = t262 * t73 + t80;
t429 = qJ(5) * t36;
t427 = t211 * t39;
t238 = pkin(3) * t262 + qJ(5);
t426 = t238 * t36;
t424 = t45 * t211;
t171 = t266 * t191;
t112 = -pkin(9) * t405 + t171 + (-pkin(3) - t444) * t267;
t230 = pkin(7) * t399;
t388 = t263 * t191 + t230;
t410 = t263 * t264;
t120 = -pkin(9) * t410 + t388;
t423 = t262 * t112 + t445 * t120;
t22 = t45 + t461;
t220 = pkin(3) * t351 + qJD(5);
t422 = t220 - t22;
t421 = t220 - t45;
t261 = qJDD(2) * pkin(2);
t416 = t173 * t226;
t414 = t250 * t267;
t413 = t251 * t264;
t412 = t251 * t267;
t409 = t263 * t265;
t407 = t263 * t268;
t403 = t264 * t269;
t272 = qJD(1) ^ 2;
t402 = t264 * t272;
t401 = t265 * t266;
t244 = pkin(3) * t266 + pkin(2);
t202 = t267 * t244;
t398 = t267 * t268;
t396 = qJ(6) + t269;
t18 = t38 + t461;
t395 = qJD(5) - t18;
t391 = t263 * t181 + t191 * t377;
t382 = qJD(2) * t264;
t390 = t266 * t181 + t382 * t444;
t227 = pkin(3) * t410;
t182 = t252 + t227;
t258 = t264 ^ 2;
t387 = -t267 ^ 2 + t258;
t384 = qJD(2) * t172;
t379 = qJD(3) * t172;
t375 = t173 * qJD(2);
t374 = t196 * qJD(3);
t368 = pkin(3) * t376;
t365 = t262 * t410;
t364 = t263 * t398;
t245 = pkin(7) * t369;
t147 = pkin(7) * t349 + t245 - t261;
t363 = g(1) * t398 + g(2) * t400 + t257;
t124 = pkin(3) * t467 + pkin(7) * t380;
t358 = t226 * t378;
t357 = t264 * t378;
t355 = t226 * t377;
t354 = t211 * t376;
t347 = -pkin(1) - t202;
t346 = t262 * t20 + t445 * t26 + t73 * t351 - t83 * t376;
t344 = t396 * t268;
t137 = t251 * t400 - t397;
t343 = -t136 * pkin(4) + qJ(5) * t137;
t139 = t250 * t265 + t251 * t398;
t342 = -t138 * pkin(4) + qJ(5) * t139;
t341 = -qJ(5) * t250 - t244;
t340 = -qJD(3) * t164 - t146;
t72 = -t299 * pkin(3) + t147;
t243 = -pkin(3) * t445 - pkin(4);
t337 = g(3) * (pkin(4) * t412 + qJ(5) * t414 + t202);
t336 = t445 * t380;
t200 = qJ(5) * t413;
t335 = -pkin(4) * t415 + t200;
t234 = g(2) * t404;
t334 = -g(1) * t406 + t234;
t19 = t39 + t472;
t55 = -qJ(5) * t267 + t423;
t330 = g(1) * t136 - g(2) * t138;
t329 = g(1) * t137 - g(2) * t139;
t326 = pkin(3) * t401 + t342;
t142 = t264 * t361 - t365;
t325 = qJ(5) * t142 - t182;
t324 = t112 * t445 - t262 * t120;
t323 = t197 * t377 - t109;
t322 = g(2) * t343;
t321 = -pkin(3) * t173 - t420;
t320 = -pkin(8) * t169 + t374;
t318 = -g(3) * t414 + t250 * t470;
t317 = -g(3) * t412 + t251 * t470;
t4 = t346 + t457;
t316 = -t250 * t366 + t200;
t314 = qJ(5) * t175 + t244;
t56 = t267 * pkin(4) - t324;
t5 = t345 - t451;
t311 = g(2) * (-t136 * pkin(5) + t343);
t44 = t262 * t82 + t80;
t309 = -pkin(7) * qJDD(2) + t373 * t448;
t121 = -t198 * t445 - t262 * t199;
t59 = t315 * qJD(2) + (-t230 + (pkin(9) * t264 - t191) * t263) * qJD(3) + t390;
t62 = -t467 * pkin(9) + (-t264 * t381 - t267 * t378) * pkin(7) + t391;
t307 = t112 * t376 + t120 * t351 + t262 * t62 - t445 * t59;
t306 = t263 * t169 - t355;
t305 = t266 * t169 + t358;
t304 = t112 * t351 - t120 * t376 + t262 * t59 + t445 * t62;
t6 = t36 * pkin(4) + t35 * qJ(5) - qJD(5) * t308 + t72;
t303 = pkin(3) * t407 - t137 * pkin(4) + t268 * pkin(7) - qJ(5) * t136 + t265 * t403;
t301 = qJDD(1) * t266 - t263 * t372;
t271 = qJD(2) ^ 2;
t297 = pkin(7) * t271 + qJDD(1) * t448 - t442;
t296 = t268 * pkin(1) + pkin(3) * t409 + t139 * pkin(4) + t265 * pkin(7) + qJ(5) * t138 + t244 * t398;
t295 = -t121 * t163 + t317;
t294 = t122 * t163 + t318;
t63 = t118 * t264 + t262 * t359 - t266 * t336;
t293 = -qJ(5) * t63 + qJD(5) * t142 - t124;
t287 = g(1) * t139 + g(2) * t137 + g(3) * t413 - t346;
t286 = -qJD(3) * pkin(8) * t226 + t147 - t261 + t437;
t3 = -pkin(5) * t36 + qJDD(6) - t6;
t285 = -t163 + t418;
t10 = qJ(5) * t382 - qJD(5) * t267 + t304;
t282 = -t44 * t211 + t288;
t281 = -t211 * t38 + t287;
t280 = g(1) * t364 + t440;
t279 = t238 * t163 - t220 * t211 - t287 + t457;
t278 = -t283 + t463;
t237 = -pkin(5) + t243;
t190 = t211 * qJ(5);
t167 = pkin(3) * t354;
t156 = t266 * t398 + t409;
t155 = -t364 + t401;
t154 = -t265 * t399 + t407;
t141 = t175 * t264;
t131 = t138 * pkin(5);
t100 = pkin(4) * t174 - t314;
t86 = qJ(6) * t174 + t122;
t85 = -t175 * qJ(6) + t121;
t81 = -t174 * t446 + t314;
t74 = pkin(4) * t141 - t325;
t64 = t263 * t336 - t262 * t357 - qJD(4) * t365 + (t262 * t380 + t264 * t453) * t266;
t58 = t420 + t465;
t57 = -t141 * t446 + t325;
t50 = -t321 + t465;
t48 = qJ(6) * t141 + t55;
t46 = t267 * pkin(5) - t142 * qJ(6) + t56;
t43 = -t420 - t459;
t32 = -t190 + t39;
t28 = pkin(4) * t211 + t394;
t27 = t321 - t459;
t21 = t44 + t472;
t15 = t19 - t190;
t13 = t211 * t446 + t395;
t12 = pkin(4) * t64 - t293;
t11 = -pkin(4) * t382 + t307;
t9 = -t446 * t64 + t293;
t8 = qJ(6) * t64 + qJD(6) * t141 + t10;
t7 = t63 * qJ(6) - qJD(2) * t366 - t142 * qJD(6) + t307;
t2 = t4 + t464;
t1 = t5 - t450;
t14 = [qJDD(1), -t438 + t442, t328, qJDD(1) * t258 + 0.2e1 * t264 * t349, 0.2e1 * t249 * t264 - 0.2e1 * t373 * t387, qJDD(2) * t264 + t267 * t271, qJDD(2) * t267 - t264 * t271, 0, t309 * t264 + (-t297 - t438) * t267, t264 * t297 + t267 * t309 + t234, -t173 * t357 + (-t264 * t302 + t267 * t375 + t405 * t468) * t266 (t266 * t172 - t173 * t263) * t380 + ((-t173 * qJD(3) + t299) * t266 + (-t379 - t275) * t263) * t264 (-t370 + (-t226 - t339) * t381) * t267 + (-t267 * t301 + t305 + t375) * t264 (t226 * t383 - t299) * t267 + (-t306 + t384) * t264, -t169 * t267 - t226 * t382 -(-t191 * t378 + t390) * t226 + t171 * t169 - g(1) * t154 - g(2) * t156 + ((t355 - t384) * pkin(7) + (-pkin(7) * t169 + qJD(2) * t196 - t340) * t263 + t323) * t267 + (-pkin(7) * t299 + t113 * qJD(2) + t147 * t263 + t266 * t374) * t264, t391 * t226 - t388 * t169 - g(1) * t153 - g(2) * t155 + (t196 * t381 + (-t358 + t375) * pkin(7) + t300) * t267 + (-t263 * t374 - t114 * qJD(2) + t147 * t266 + (t370 + t301 * t264 + (-t226 + t339) * t381) * pkin(7)) * t264, -t142 * t35 - t308 * t63, t103 * t63 + t141 * t35 - t142 * t36 - t308 * t64, t142 * t163 + t211 * t63 + t267 * t35 + t308 * t382, -t103 * t382 - t141 * t163 + t211 * t64 + t267 * t36, -t163 * t267 - t211 * t382, t124 * t103 + t72 * t141 + t163 * t324 + t182 * t36 + t211 * t307 + t267 * t345 - t312 * t64 + t38 * t382 + t329, t124 * t308 + t72 * t142 - t163 * t423 - t182 * t35 + t211 * t304 + t267 * t346 + t312 * t63 - t382 * t39 - t330, t103 * t12 + t11 * t211 + t141 * t6 - t163 * t56 + t267 * t5 - t28 * t382 + t36 * t74 + t47 * t64 + t329, -t10 * t103 + t11 * t308 - t141 * t4 + t142 * t5 - t28 * t63 - t32 * t64 - t35 * t56 - t36 * t55 - t334, -t10 * t211 - t12 * t308 - t142 * t6 + t163 * t55 - t267 * t4 + t32 * t382 + t35 * t74 + t47 * t63 + t330, t4 * t55 + t32 * t10 + t6 * t74 + t47 * t12 + t5 * t56 + t28 * t11 - g(1) * (t265 * t347 + t303) - g(2) * (-t268 * t403 + t296) t1 * t267 - t103 * t9 - t13 * t382 - t141 * t3 - t163 * t46 + t211 * t7 - t23 * t64 - t36 * t57 + t329, t142 * t3 + t15 * t382 + t163 * t48 - t2 * t267 - t211 * t8 - t23 * t63 + t308 * t9 - t35 * t57 + t330, -t1 * t142 + t103 * t8 + t13 * t63 + t141 * t2 + t15 * t64 - t308 * t7 + t35 * t46 + t36 * t48 + t334, t2 * t48 + t15 * t8 + t1 * t46 + t13 * t7 + t3 * t57 + t23 * t9 - g(1) * (-pkin(5) * t137 + t303) - g(2) * (pkin(5) * t139 - t264 * t344 + t296) - (qJ(6) * t264 + t347) * t442; 0, 0, 0, -t267 * t402, t387 * t272, t369, t249, qJDD(2), pkin(1) * t402 - t245 + t449 (pkin(1) * t272 - pkin(7) * qJDD(1)) * t267 + t363, -t319 * t263 ^ 2 + (t263 * t290 - t416) * t266 (-t338 + t416) * t263 + (t379 + 0.2e1 * t370 + t313 * t266 + (-t357 + (-t172 + t381) * t267) * qJD(1)) * t266 (-t173 * t264 + t226 * t399) * qJD(1) + t306 (-t172 * t264 - t226 * t408) * qJD(1) + t305, t226 * t386, -pkin(2) * t338 + t389 * t226 + t320 * t263 + (-t113 * t264 + (pkin(7) * t172 - t196 * t263) * t267) * qJD(1) + (t264 * t328 - t286) * t266, -t158 * t226 + (-t267 * pkin(7) * t173 + t114 * t264) * qJD(1) + (-pkin(2) * t313 + (t226 * t252 + (-t196 - t428) * t267) * qJD(1) + t320) * t266 + ((pkin(2) * t372 - t328) * t264 + t286) * t263, -t175 * t35 - t308 * t393, t103 * t393 + t174 * t35 - t175 * t36 - t308 * t392, t163 * t175 + t211 * t393 - t308 * t386, t103 * t386 - t163 * t174 + t211 * t392, t211 * t386, t332 * t103 + t72 * t174 + t211 * t462 - t244 * t36 - t312 * t392 - t38 * t386 + t295, t72 * t175 + t211 * t473 + t244 * t35 + t332 * t308 + t393 * t312 + t39 * t386 - t294, t100 * t36 + t103 * t433 + t174 * t6 + t211 * t431 + t28 * t386 + t392 * t47 + t295, -t103 * t432 - t121 * t35 - t122 * t36 - t174 * t4 + t175 * t5 - t28 * t393 + t308 * t431 - t32 * t392 - t363, t100 * t35 - t175 * t6 - t211 * t432 - t308 * t433 - t32 * t386 + t393 * t47 + t294, t4 * t122 + t6 * t100 + t5 * t121 - t337 + t433 * t47 + t432 * t32 + t431 * t28 + t328 * t269 * t267 + (g(3) * t269 + t328 * (pkin(4) * t251 - t341)) * t264, -t103 * t436 + t13 * t386 - t163 * t85 - t174 * t3 + t211 * t434 - t23 * t392 - t36 * t81 + t317, -t15 * t386 + t163 * t86 + t175 * t3 - t211 * t435 - t23 * t393 + t308 * t436 - t35 * t81 + t318, -t1 * t175 + t103 * t435 + t13 * t393 + t15 * t392 + t174 * t2 - t308 * t434 + t35 * t85 + t36 * t86 + t363, t2 * t86 + t1 * t85 + t3 * t81 - t337 + t436 * t23 + t435 * t15 + t434 * t13 + (-g(3) * pkin(5) * t251 + g(1) * t344 + t396 * t439) * t267 + (g(3) * t396 + t328 * (t251 * t446 - t341)) * t264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t173 * t172, -t172 ^ 2 + t173 ^ 2, t172 * t226 + t275, t299 - t416, t169, -g(1) * t155 + t440 - t114 * t226 - t173 * t196 + (t340 + t257) * t263 - t323, g(1) * t156 - g(2) * t154 + g(3) * t405 - t113 * t226 - t172 * t196 - t300, t418, t475, -t289, t469, t163, t460 + (-t103 * t173 + t163 * t445 + t354) * pkin(3) + t282, -t471 - t424 + (-t163 * t262 - t173 * t308 + t211 * t351) * pkin(3) + t287, -t103 * t50 - t163 * t243 + t167 + t282 + t451 - t463, -t243 * t35 - t426 + (t32 - t44 + t368) * t308 + (-t421 + t28) * t103, t308 * t50 + t279 + t424 - t474, t4 * t238 + t5 * t243 - t47 * t50 - t28 * t44 - g(1) * t326 - t322 - g(3) * (-t227 + t335) + t421 * t32 + (t28 * t376 + t280) * pkin(3), t103 * t27 - t163 * t237 - t21 * t211 + t167 - t276, t211 * t22 - t27 * t308 + t279 + t466, t237 * t35 + t426 + (-t15 + t21 - t368) * t308 + (t422 - t13) * t103, t2 * t238 + t1 * t237 - t13 * t21 - t23 * t27 - g(1) * (-t131 + t326) - t311 - g(3) * (-t227 + t316) + t422 * t15 + (t13 * t376 + t280) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t418, t475, -t289, t469, t163, t288 - t427 + t460, t281 - t471, -t103 * t58 + t150 - t278 - t427, pkin(4) * t35 - t429 + (t32 - t39) * t308 + (t28 - t394) * t103, t308 * t58 - t281 + t458 - t474, -t5 * pkin(4) - g(1) * t342 - g(3) * t335 + t4 * qJ(5) - t28 * t39 + t32 * t394 - t47 * t58 - t322, t103 * t43 + t163 * t446 - t19 * t211 - t276, t18 * t211 - t308 * t43 - t287 + t458 + t466, t429 - t446 * t35 + (-t15 + t19) * t308 + (-t13 + t395) * t103, t2 * qJ(5) - t1 * t446 - t13 * t19 - t23 * t43 - g(1) * (-t131 + t342) - t311 - g(3) * t316 + t395 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t285, -t289, t456, t211 * t32 + t278, t285, t456, t289, t15 * t211 + t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36 + t417, -t35 + t419, -t101 - t447, -t103 * t15 + t13 * t308 + t3 + t449;];
tau_reg  = t14;
