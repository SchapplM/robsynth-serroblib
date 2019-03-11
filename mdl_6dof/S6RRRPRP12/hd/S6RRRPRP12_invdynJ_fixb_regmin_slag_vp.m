% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPRP12
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRP12_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP12_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP12_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:59:21
% EndTime: 2019-03-09 17:59:49
% DurationCPUTime: 13.43s
% Computational Cost: add. (11157->708), mult. (26793->897), div. (0->0), fcn. (20979->10), ass. (0->313)
t262 = sin(qJ(3));
t266 = cos(qJ(2));
t260 = sin(pkin(6));
t405 = qJD(1) * t260;
t376 = t266 * t405;
t348 = t262 * t376;
t402 = qJD(3) * t262;
t463 = t348 - t402;
t224 = -qJD(3) + t376;
t239 = pkin(8) * t376;
t263 = sin(qJ(2));
t433 = cos(pkin(6));
t361 = t433 * qJD(1);
t346 = pkin(1) * t361;
t183 = t263 * t346 + t239;
t476 = -pkin(3) * t463 - t262 * qJD(4) - t183;
t453 = cos(qJ(3));
t475 = -t476 + t224 * (pkin(10) * t262 - qJ(4) * t453);
t254 = t453 * pkin(9);
t227 = t453 * pkin(4) + t254;
t377 = t263 * t405;
t180 = -pkin(8) * t377 + t266 * t346;
t334 = pkin(2) * t263 - pkin(9) * t266;
t181 = t334 * t405;
t327 = t262 * t180 - t181 * t453;
t380 = t266 * t453;
t457 = pkin(3) + pkin(10);
t474 = t227 * qJD(3) - (pkin(4) * t380 - t263 * t457) * t405 - t327;
t383 = t260 * t453;
t344 = qJD(1) * t383;
t320 = t266 * t344;
t372 = qJD(3) * t453;
t286 = t320 - t372;
t473 = -t262 * qJ(4) - pkin(2);
t318 = t361 + qJD(2);
t302 = t262 * t318;
t157 = t263 * t344 + t302;
t149 = qJD(5) + t157;
t265 = cos(qJ(5));
t261 = sin(qJ(5));
t359 = t149 * t261;
t292 = t453 * t318;
t358 = t433 * qJDD(1);
t314 = t358 + qJDD(2);
t423 = t260 * t263;
t392 = t262 * t423;
t347 = qJD(3) * t392;
t351 = t263 * t383;
t78 = qJD(1) * t347 - qJD(2) * t320 - qJD(3) * t292 - qJDD(1) * t351 - t262 * t314;
t75 = -qJDD(5) + t78;
t472 = -t149 * t359 - t265 * t75;
t138 = pkin(9) * t318 + t183;
t316 = -pkin(2) * t266 - pkin(9) * t263 - pkin(1);
t147 = t316 * t405;
t76 = t138 * t262 - t453 * t147;
t415 = -qJD(4) - t76;
t217 = -pkin(3) * t453 + t473;
t200 = -pkin(10) * t453 + t217;
t456 = pkin(4) + pkin(9);
t226 = t456 * t262;
t400 = qJD(5) * t265;
t401 = qJD(5) * t261;
t471 = t200 * t401 - t226 * t400 - t474 * t261 + t475 * t265;
t264 = sin(qJ(1));
t454 = cos(qJ(1));
t337 = t433 * t454;
t196 = t263 * t337 + t264 * t266;
t127 = t196 * t262 + t383 * t454;
t195 = t263 * t264 - t266 * t337;
t88 = t127 * t261 + t195 * t265;
t89 = t127 * t265 - t195 * t261;
t289 = qJD(3) * t302 - t453 * t314;
t342 = t263 * t372;
t397 = qJDD(1) * t263;
t366 = t262 * t397;
t403 = qJD(2) * t266;
t373 = t262 * t403;
t269 = (qJD(1) * (t342 + t373) + t366) * t260 + t289;
t257 = t260 ^ 2;
t470 = 0.2e1 * t257;
t410 = t265 * t200 + t261 * t226;
t447 = -t286 * pkin(5) - qJD(5) * t410 + t475 * t261 + t474 * t265;
t446 = -t286 * qJ(6) + t262 * qJD(6) - t471;
t145 = t261 * t377 - t265 * t348;
t416 = t262 * t266;
t162 = (t261 * t416 + t263 * t265) * t260;
t146 = qJD(1) * t162;
t325 = pkin(5) * t265 + qJ(6) * t261;
t313 = -pkin(4) - t325;
t379 = t453 * qJ(6);
t382 = t261 * t453;
t411 = t453 * t180 + t262 * t181;
t84 = -qJ(4) * t377 - t411;
t71 = -pkin(4) * t348 - t84;
t443 = -pkin(5) * t145 + qJ(6) * t146 - t71 + qJD(6) * t382 + (-pkin(5) * t382 + t265 * t379) * qJD(5) + (-pkin(9) + t313) * t402;
t155 = t262 * t377 - t292;
t396 = qJDD(1) * t266;
t243 = t260 * t396;
t398 = qJD(1) * qJD(2);
t368 = t263 * t398;
t341 = t260 * t368;
t177 = qJDD(3) - t243 + t341;
t38 = -t155 * t400 - t265 * t177 - t224 * t401 - t261 * t269;
t106 = -t265 * t155 - t224 * t261;
t431 = t106 * t149;
t469 = t38 - t431;
t108 = t155 * t261 - t224 * t265;
t39 = qJD(5) * t108 + t177 * t261 - t265 * t269;
t429 = t108 * t149;
t468 = -t39 + t429;
t467 = t106 * t224 + t472;
t409 = t286 * qJ(4) + t476;
t364 = t266 * t433;
t312 = pkin(1) * t364 - pkin(8) * t423;
t459 = t149 ^ 2;
t296 = t261 * t75 - t265 * t459;
t466 = t108 * t224 + t296;
t465 = (qJDD(2) + 0.2e1 * t358) * t260;
t464 = -t262 * t177 + t224 * t372;
t414 = pkin(4) * t157 - t415;
t384 = t260 * t454;
t128 = t196 * t453 - t262 * t384;
t365 = t263 * t433;
t198 = -t264 * t365 + t266 * t454;
t422 = t260 * t264;
t132 = t198 * t453 + t262 * t422;
t194 = t262 * t433 + t351;
t298 = g(1) * t132 + g(2) * t128 + g(3) * t194;
t336 = t433 * t453;
t350 = t260 * t380;
t124 = qJD(2) * t350 + qJD(3) * t336 - t347;
t421 = t260 * t266;
t407 = pkin(1) * t365 + pkin(8) * t421;
t175 = pkin(9) * t433 + t407;
t408 = pkin(2) * t421 + pkin(9) * t423;
t176 = -pkin(1) * t260 - t408;
t308 = t334 * qJD(2);
t182 = t260 * t308;
t184 = t312 * qJD(2);
t305 = t175 * t372 + t176 * t402 - t182 * t453 + t262 * t184;
t404 = qJD(2) * t263;
t375 = t260 * t404;
t32 = t124 * pkin(4) - t375 * t457 + t305;
t328 = -t262 * t175 + t176 * t453;
t82 = pkin(3) * t421 - t328;
t58 = t194 * pkin(4) + pkin(10) * t421 + t82;
t193 = -t336 + t392;
t448 = t193 * pkin(10);
t174 = -pkin(2) * t433 - t312;
t187 = t193 * pkin(3);
t360 = t194 * qJ(4) - t187;
t80 = t174 - t360;
t65 = t80 + t448;
t322 = t261 * t58 + t265 * t65;
t123 = qJD(3) * t194 + t260 * t373;
t345 = pkin(1) * qJD(2) * t433;
t374 = t260 * t403;
t185 = pkin(8) * t374 + t263 * t345;
t297 = -qJ(4) * t124 - qJD(4) * t194 + t185;
t37 = t123 * t457 + t297;
t461 = -qJD(5) * t322 - t261 * t37 + t265 * t32;
t460 = t108 ^ 2;
t458 = t157 ^ 2;
t268 = qJD(1) ^ 2;
t455 = pkin(5) * t75;
t452 = pkin(3) * t177;
t77 = t453 * t138 + t262 * t147;
t56 = -pkin(4) * t155 + t77;
t432 = qJ(4) * t155;
t69 = t157 * t457 + t432;
t445 = t261 * t56 + t265 * t69;
t442 = qJ(6) * t75;
t44 = t224 * t457 + t414;
t137 = -pkin(2) * t318 - t180;
t272 = -t157 * qJ(4) + t137;
t48 = t155 * t457 + t272;
t19 = t261 * t44 + t265 * t48;
t17 = qJ(6) * t149 + t19;
t441 = t149 * t17;
t440 = t149 * t19;
t439 = t224 * t77;
t437 = t457 * t75;
t436 = t78 * t262;
t435 = qJD(5) * t325 - qJD(6) * t265 - t157 * t313 - t415;
t434 = -t456 * t402 - t71;
t430 = t108 * t106;
t427 = t157 * t155;
t169 = t177 * qJ(4);
t424 = t257 * t268;
t420 = t261 * t262;
t417 = t262 * t265;
t18 = -t261 * t48 + t265 * t44;
t413 = qJD(6) - t18;
t412 = t453 * t175 + t262 * t176;
t258 = t263 ^ 2;
t406 = -t266 ^ 2 + t258;
t399 = qJD(5) * t457;
t394 = qJ(4) * t421;
t393 = t266 * t424;
t391 = t265 * t421;
t386 = t195 * t453;
t390 = -pkin(3) * t386 + t473 * t195;
t197 = t263 * t454 + t264 * t364;
t385 = t197 * t453;
t389 = -pkin(3) * t385 + t473 * t197;
t315 = qJD(1) * t345;
t338 = pkin(1) * t358;
t387 = pkin(8) * t243 + t263 * t338 + t266 * t315;
t381 = t265 * t453;
t378 = t453 * t177;
t371 = t453 * qJD(5);
t370 = pkin(1) * t470;
t369 = t266 * t398;
t367 = t260 * t397;
t100 = (qJD(1) * t308 + qJDD(1) * t316) * t260;
t290 = -pkin(8) * t341 + t387;
t97 = pkin(9) * t314 + t290;
t356 = -t453 * t100 + t138 * t372 + t147 * t402 + t262 * t97;
t321 = qJDD(4) + t356;
t11 = -pkin(4) * t78 - t177 * t457 + t321;
t353 = pkin(8) * t367 + qJD(2) * t239 + t263 * t315 - t266 * t338;
t98 = -pkin(2) * t314 + t353;
t21 = pkin(3) * t269 + t78 * qJ(4) - t157 * qJD(4) + t98;
t15 = pkin(10) * t269 + t21;
t363 = -t265 * t11 + t261 * t15 + t48 * t400 + t44 * t401;
t357 = -t262 * t100 + t138 * t402 - t147 * t372 - t453 * t97;
t354 = pkin(9) * t378;
t352 = pkin(3) * t350 + t262 * t394 + t408;
t131 = t198 * t262 - t264 * t383;
t91 = -t131 * t265 + t197 * t261;
t340 = g(1) * t89 + g(2) * t91;
t92 = t131 * t261 + t197 * t265;
t339 = g(1) * t88 - g(2) * t92;
t335 = t260 * t268 * t433;
t333 = g(1) * t127 - g(2) * t131;
t332 = -g(1) * t128 + g(2) * t132;
t331 = -g(1) * t195 + g(2) * t197;
t330 = g(1) * t198 + g(2) * t196;
t125 = t193 * t265 + t261 * t421;
t126 = -t193 * t261 + t391;
t326 = pkin(5) * t125 - qJ(6) * t126;
t323 = -t261 * t65 + t265 * t58;
t317 = 0.2e1 * t361 + qJD(2);
t209 = t224 * qJ(4);
t49 = -t209 + t56;
t81 = t394 - t412;
t207 = qJD(4) * t224;
t20 = -t169 + t207 + t357;
t214 = pkin(5) * t261 - qJ(6) * t265 + qJ(4);
t3 = t261 * t11 + t265 * t15 + t44 * t400 - t401 * t48;
t309 = t261 * t32 + t265 * t37 + t58 * t400 - t401 * t65;
t28 = pkin(5) * t106 - qJ(6) * t108 + t49;
t307 = t149 * t28 + t437;
t306 = t454 * pkin(1) + t198 * pkin(2) + t132 * pkin(3) + pkin(8) * t422 + qJ(4) * t131;
t304 = -t175 * t402 + t176 * t372 + t262 * t182 + t453 * t184;
t101 = t195 * t417 + t196 * t261;
t103 = t197 * t417 + t198 * t261;
t161 = t261 * t423 - t262 * t391;
t301 = g(1) * t103 + g(2) * t101 + g(3) * t161;
t102 = -t195 * t420 + t196 * t265;
t104 = -t197 * t420 + t198 * t265;
t300 = -g(1) * t104 - g(2) * t102 - g(3) * t162;
t299 = g(1) * t131 + g(2) * t127 + g(3) * t193;
t294 = g(1) * t197 + g(2) * t195 - g(3) * t421;
t66 = -pkin(4) * t193 - t81;
t293 = -pkin(1) * t264 - t196 * pkin(2) - pkin(3) * t128 + pkin(8) * t384 - qJ(4) * t127;
t288 = t261 * t371 + t265 * t402 + t145;
t287 = -t261 * t402 + t265 * t371 + t146;
t285 = g(1) * t91 - g(2) * t89 - g(3) * t125 - t363;
t25 = t321 - t452;
t284 = -t20 * t453 + t25 * t262 - t330;
t283 = t149 * t399 - t298;
t282 = -t298 - t357;
t281 = t299 - t356;
t280 = pkin(9) * qJD(3) * t224 + t294;
t14 = -pkin(4) * t269 - t20;
t5 = t39 * pkin(5) + t38 * qJ(6) - t108 * qJD(6) + t14;
t279 = -t283 - t5;
t43 = -qJ(4) * t375 + qJD(4) * t421 - t304;
t278 = -t155 * t224 - t78;
t277 = -g(1) * t92 - g(2) * t88 + g(3) * t126 + t3;
t33 = -pkin(4) * t123 - t43;
t276 = t108 * t28 + qJDD(6) - t285;
t275 = g(1) * t385 + g(2) * t386 - g(3) * t350;
t64 = t155 * pkin(3) + t272;
t274 = t157 * t64 + qJDD(4) - t281;
t271 = pkin(9) * t464 + t275;
t148 = pkin(5) * t381 + t261 * t379 + t227;
t120 = t131 * pkin(3);
t118 = t127 * pkin(3);
t111 = -pkin(5) * t262 + t200 * t261 - t226 * t265;
t110 = qJ(6) * t262 + t410;
t99 = pkin(3) * t157 + t432;
t86 = -pkin(3) * t377 + t327;
t70 = t209 - t77;
t68 = pkin(3) * t224 - t415;
t62 = qJD(5) * t125 + t123 * t261 + t265 * t375;
t61 = -t123 * t265 - qJD(5) * t391 + (qJD(5) * t193 + t375) * t261;
t54 = pkin(5) * t108 + qJ(6) * t106;
t47 = pkin(3) * t123 + t297;
t46 = -pkin(3) * t375 + t305;
t34 = -t326 + t66;
t27 = -pkin(5) * t194 - t323;
t26 = qJ(6) * t194 + t322;
t24 = pkin(5) * t155 + t261 * t69 - t265 * t56;
t23 = -qJ(6) * t155 + t445;
t16 = -pkin(5) * t149 + t413;
t8 = pkin(5) * t61 - qJ(6) * t62 + qJD(6) * t126 + t33;
t7 = -pkin(5) * t124 - t461;
t6 = qJ(6) * t124 + qJD(6) * t194 + t309;
t2 = qJDD(6) + t363 + t455;
t1 = qJD(6) * t149 + t3 - t442;
t4 = [qJDD(1), g(1) * t264 - g(2) * t454, g(1) * t454 + g(2) * t264 (qJDD(1) * t258 + 0.2e1 * t266 * t368) * t257 (t263 * t396 - t398 * t406) * t470, t263 * t465 + t317 * t374, t266 * t465 - t317 * t375, t314 * t433, -t185 * t318 + t312 * t314 - t353 * t433 + g(1) * t196 - g(2) * t198 + (-t368 + t396) * t370, -t184 * t318 - t407 * t314 - t290 * t433 + (-t369 - t397) * t370 + t331, t124 * t157 - t194 * t78, -t157 * t123 - t124 * t155 + t78 * t193 - t194 * t269, -t124 * t224 + t177 * t194 + (t157 * t404 + t266 * t78) * t260, t123 * t224 - t193 * t177 + (t289 * t266 - t155 * t404 - (-qJD(1) * t342 - t262 * t369 - t366) * t421) * t260 (-t177 * t266 - t224 * t404) * t260, t137 * t123 + t185 * t155 + t174 * t269 + t177 * t328 + t98 * t193 + t224 * t305 + t356 * t421 - t375 * t76 - t332, t304 * t224 - t412 * t177 + t185 * t157 - t174 * t78 + t98 * t194 + t137 * t124 + (-t266 * t357 - t404 * t77) * t260 - t333, t70 * t123 + t68 * t124 + t43 * t155 + t46 * t157 + t20 * t193 + t25 * t194 + t269 * t81 - t82 * t78 - t331, -t64 * t123 - t47 * t155 + t82 * t177 - t21 * t193 - t46 * t224 - t25 * t421 - t269 * t80 + t375 * t68 + t332, -t124 * t64 - t157 * t47 - t177 * t81 - t194 * t21 + t224 * t43 + t78 * t80 + (t20 * t266 - t404 * t70) * t260 + t333, t21 * t80 + t64 * t47 + t20 * t81 + t70 * t43 + t25 * t82 + t68 * t46 - g(1) * (-pkin(9) * t195 + t293) - g(2) * (pkin(9) * t197 + t306) t108 * t62 + t126 * t38, -t106 * t62 - t108 * t61 - t125 * t38 + t126 * t39, t108 * t124 + t126 * t75 + t149 * t62 - t194 * t38, -t106 * t124 - t125 * t75 - t149 * t61 - t194 * t39, t124 * t149 - t194 * t75, t33 * t106 + t18 * t124 - t14 * t125 + t149 * t461 - t194 * t363 - t323 * t75 + t66 * t39 + t49 * t61 + t339, t33 * t108 - t19 * t124 - t14 * t126 - t149 * t309 - t3 * t194 + t322 * t75 - t66 * t38 + t49 * t62 + t340, t106 * t8 - t124 * t16 - t125 * t5 - t149 * t7 - t194 * t2 + t27 * t75 + t28 * t61 + t34 * t39 + t339, t1 * t125 - t106 * t6 + t108 * t7 - t126 * t2 + t16 * t62 - t17 * t61 - t26 * t39 - t27 * t38 - t332, t1 * t194 - t108 * t8 + t124 * t17 + t126 * t5 + t149 * t6 - t26 * t75 - t28 * t62 + t34 * t38 - t340, t1 * t26 + t17 * t6 + t5 * t34 + t28 * t8 + t2 * t27 + t16 * t7 - g(1) * (-pkin(5) * t88 - pkin(10) * t128 + qJ(6) * t89 - t195 * t456 + t293) - g(2) * (pkin(5) * t92 + pkin(10) * t132 + qJ(6) * t91 + t197 * t456 + t306); 0, 0, 0, -t263 * t393, t406 * t424, -t266 * t335 + t367, t263 * t335 + t243, t314, pkin(1) * t263 * t424 + t183 * t318 + t294 - t353, pkin(1) * t393 + t180 * t318 + (pkin(8) * t398 + g(3)) * t423 + t330 - t387, -t157 * t286 - t436, t155 * t286 + t157 * t463 - t262 * t269 - t78 * t453 (-t157 * t263 + t224 * t380) * t405 - t464, t224 * t402 + t378 + (t155 * t263 - t224 * t416) * t405, t224 * t377, -pkin(2) * t269 - t137 * t463 - t183 * t155 - t327 * t224 + t76 * t377 - t98 * t453 + t271, -t354 + pkin(2) * t78 + t137 * t372 - t411 * t224 - t183 * t157 + (-t137 * t380 + t263 * t77) * t405 + (-t280 + t98) * t262, -g(3) * t423 - t84 * t155 - t86 * t157 - t254 * t269 + t284 - t463 * t70 - t286 * t68 + (t155 * t402 + t157 * t372 - t436) * pkin(9), -t155 * t409 + t21 * t453 - t217 * t269 + t86 * t224 - t68 * t377 + t463 * t64 - t271, t354 - t64 * t372 + t217 * t78 - t84 * t224 - t409 * t157 + (t263 * t70 + t380 * t64) * t405 + (-t21 + t280) * t262, t21 * t217 - t70 * t84 - t68 * t86 - g(1) * t389 - g(2) * t390 - g(3) * t352 + t409 * t64 + ((t262 * t70 + t453 * t68) * qJD(3) + t284) * pkin(9), -t108 * t287 + t38 * t382, t106 * t287 + t108 * t288 + t38 * t381 + t382 * t39, -t108 * t286 - t149 * t287 - t38 * t262 + t382 * t75, t106 * t286 + t149 * t288 - t39 * t262 + t381 * t75, -t149 * t286 - t75 * t262, -t49 * t145 + t227 * t39 - t363 * t262 - t286 * t18 + t434 * t106 + (-t49 * t371 + t200 * t75 + (-qJD(5) * t226 + t475) * t149) * t261 + (-t49 * t402 + t14 * t453 - t226 * t75 + (-qJD(5) * t200 + t474) * t149) * t265 + t300, t434 * t108 - t14 * t382 + t471 * t149 + t286 * t19 - t227 * t38 - t3 * t262 - t287 * t49 + t410 * t75 + t301, t106 * t443 + t111 * t75 + t148 * t39 + t149 * t447 + t16 * t286 - t2 * t262 - t28 * t288 + t381 * t5 + t300, -t1 * t381 - t106 * t446 - t108 * t447 - t110 * t39 - t111 * t38 - t16 * t287 + t17 * t288 - t2 * t382 + t275, t1 * t262 - t108 * t443 - t110 * t75 + t148 * t38 + t149 * t446 - t17 * t286 + t28 * t287 + t382 * t5 - t301, t1 * t110 + t5 * t148 + t2 * t111 - g(1) * (t104 * pkin(5) - pkin(10) * t385 + t103 * qJ(6) + t198 * t456 + t389) - g(2) * (t102 * pkin(5) - pkin(10) * t386 + t101 * qJ(6) + t196 * t456 + t390) - g(3) * (t162 * pkin(5) + t161 * qJ(6) + (pkin(4) * t263 + pkin(10) * t380) * t260 + t352) + t443 * t28 + t446 * t17 - t447 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t427, -t155 ^ 2 + t458, t278, -t157 * t224 - t269, t177, -t137 * t157 + t281 - t439, t137 * t155 + t224 * t76 - t282, pkin(3) * t78 - qJ(4) * t269 + (-t70 - t77) * t157 + (t68 + t415) * t155, t155 * t99 + t274 + t439 - 0.2e1 * t452, -t155 * t64 + t157 * t99 + t224 * t415 + 0.2e1 * t169 - t207 + t282, -t20 * qJ(4) - t25 * pkin(3) - t64 * t99 - t68 * t77 - g(1) * (qJ(4) * t132 - t120) - g(2) * (qJ(4) * t128 - t118) - g(3) * t360 + t415 * t70, -t108 * t359 - t265 * t38 (-t39 - t429) * t265 + (t38 + t431) * t261, t108 * t155 + t472, -t106 * t155 + t296, t149 * t155, qJ(4) * t39 + t18 * t155 + t414 * t106 + (t437 + (t49 - t56) * t149) * t265 + (t14 + (t69 + t399) * t149 - t298) * t261, -qJ(4) * t38 + t445 * t149 - t19 * t155 + t414 * t108 + (-t149 * t49 - t437) * t261 + (t14 + t283) * t265, t106 * t435 + t149 * t24 - t155 * t16 + t214 * t39 - t261 * t279 + t265 * t307, t106 * t23 - t108 * t24 + (-t157 * t17 - t457 * t38 + t2 + (t106 * t457 - t17) * qJD(5)) * t265 + (-t157 * t16 + t457 * t39 - t1 + (-t108 * t457 - t16) * qJD(5)) * t261 + t299, -t108 * t435 - t149 * t23 + t155 * t17 + t214 * t38 + t261 * t307 + t265 * t279, -t17 * t23 - t16 * t24 - g(1) * (-pkin(10) * t131 - t120) - g(2) * (-pkin(10) * t127 - t118) - g(3) * (-t187 - t448) + t435 * t28 - (t1 * t261 - t2 * t265 + (t16 * t261 + t17 * t265) * qJD(5)) * t457 + (t5 - t298) * t214; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t278, t177 - t427, -t224 ^ 2 - t458, -t224 * t70 + t274 - t452, 0, 0, 0, 0, 0, t467, t466, t467, t261 * t468 + t265 * t469, -t466, t224 * t28 + (-t2 + t441) * t265 + (t149 * t16 + t1) * t261 - t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t430, -t106 ^ 2 + t460, -t469, t468, -t75, -t108 * t49 + t285 + t440, t106 * t49 + t149 * t18 - t277, -t106 * t54 - t276 + t440 - 0.2e1 * t455, pkin(5) * t38 - qJ(6) * t39 + (t17 - t19) * t108 + (t16 - t413) * t106, -0.2e1 * t442 - t106 * t28 + t108 * t54 + (0.2e1 * qJD(6) - t18) * t149 + t277, t1 * qJ(6) - t2 * pkin(5) - t28 * t54 - t16 * t19 - g(1) * (-pkin(5) * t91 + qJ(6) * t92) - g(2) * (pkin(5) * t89 + qJ(6) * t88) - g(3) * t326 + t413 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75 + t430, -t469, -t459 - t460, t276 - t441 + t455;];
tau_reg  = t4;
