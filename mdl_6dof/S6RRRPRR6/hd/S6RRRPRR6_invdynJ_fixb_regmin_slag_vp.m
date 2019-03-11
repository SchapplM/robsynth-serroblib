% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% tau_reg [6x33]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRPRR6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:29:35
% EndTime: 2019-03-09 18:29:58
% DurationCPUTime: 10.90s
% Computational Cost: add. (13935->588), mult. (32069->796), div. (0->0), fcn. (24690->16), ass. (0->292)
t304 = sin(qJ(6));
t309 = cos(qJ(6));
t306 = sin(qJ(3));
t311 = cos(qJ(3));
t312 = cos(qJ(2));
t381 = qJD(1) * qJD(2);
t367 = t312 * t381;
t307 = sin(qJ(2));
t380 = t307 * qJDD(1);
t382 = t311 * qJD(2);
t389 = qJD(3) * t307;
t462 = -qJD(1) * t389 + qJDD(2);
t173 = qJD(3) * t382 + (t367 + t380) * t311 + t462 * t306;
t396 = qJD(1) * t312;
t174 = t306 * (qJD(2) * (qJD(3) + t396) + t380) - t462 * t311;
t301 = sin(pkin(11));
t302 = cos(pkin(11));
t103 = -t173 * t301 - t174 * t302;
t104 = t173 * t302 - t174 * t301;
t397 = qJD(1) * t307;
t374 = t306 * t397;
t251 = t374 - t382;
t393 = qJD(2) * t306;
t253 = t311 * t397 + t393;
t180 = t251 * t301 - t253 * t302;
t305 = sin(qJ(5));
t310 = cos(qJ(5));
t341 = -t251 * t302 - t253 * t301;
t385 = qJD(5) * t310;
t386 = qJD(5) * t305;
t35 = t305 * t103 + t310 * t104 + t180 * t386 + t341 * t385;
t457 = -t310 * t180 + t305 * t341;
t36 = qJD(5) * t457 - t310 * t103 + t104 * t305;
t124 = t180 * t305 + t310 * t341;
t53 = t124 * t304 + t309 * t457;
t319 = -qJD(6) * t53 - t304 * t35 - t309 * t36;
t279 = -qJD(3) + t396;
t269 = -qJD(5) + t279;
t262 = -qJD(6) + t269;
t428 = t262 * t53;
t456 = t319 - t428;
t57 = t124 * t309 - t304 * t457;
t465 = t53 * t57;
t427 = t262 * t57;
t383 = qJD(6) * t309;
t384 = qJD(6) * t304;
t8 = t124 * t383 - t304 * t36 + t309 * t35 - t384 * t457;
t460 = t8 + t427;
t412 = t311 * t312;
t442 = pkin(3) * t307;
t337 = -qJ(4) * t412 + t442;
t303 = -qJ(4) - pkin(8);
t359 = qJD(3) * t303;
t350 = pkin(2) * t307 - pkin(8) * t312;
t254 = t350 * qJD(1);
t402 = pkin(7) * t374 + t311 * t254;
t484 = -qJD(1) * t337 - qJD(4) * t306 + t311 * t359 - t402;
t236 = t306 * t254;
t387 = qJD(4) * t311;
t414 = t307 * t311;
t415 = t306 * t312;
t483 = t236 + (-pkin(7) * t414 - qJ(4) * t415) * qJD(1) - t306 * t359 - t387;
t461 = t53 ^ 2 - t57 ^ 2;
t296 = qJ(3) + pkin(11) + qJ(5);
t286 = qJ(6) + t296;
t276 = sin(t286);
t277 = cos(t286);
t313 = cos(qJ(1));
t308 = sin(qJ(1));
t413 = t308 * t312;
t197 = t276 * t413 + t277 * t313;
t411 = t312 * t313;
t199 = -t276 * t411 + t277 * t308;
t295 = t312 * qJDD(1);
t448 = -t307 * t381 + t295;
t243 = qJDD(3) - t448;
t240 = qJDD(5) + t243;
t466 = pkin(9) * t180;
t261 = -pkin(2) * t312 - pkin(8) * t307 - pkin(1);
t241 = t261 * qJD(1);
t290 = pkin(7) * t396;
t267 = qJD(2) * pkin(8) + t290;
t190 = t311 * t241 - t267 * t306;
t151 = -qJ(4) * t253 + t190;
t140 = -pkin(3) * t279 + t151;
t191 = t241 * t306 + t267 * t311;
t152 = -qJ(4) * t251 + t191;
t145 = t301 * t152;
t81 = t302 * t140 - t145;
t61 = -pkin(4) * t279 + t466 + t81;
t455 = pkin(9) * t341;
t417 = t302 * t152;
t82 = t301 * t140 + t417;
t64 = t82 + t455;
t27 = t305 * t61 + t310 * t64;
t255 = t350 * qJD(2);
t192 = qJD(1) * t255 + qJDD(1) * t261;
t186 = t311 * t192;
t226 = pkin(7) * t448 + qJDD(2) * pkin(8);
t60 = pkin(3) * t243 - qJ(4) * t173 - qJD(3) * t191 - qJD(4) * t253 - t306 * t226 + t186;
t388 = qJD(3) * t311;
t390 = qJD(3) * t306;
t328 = t306 * t192 + t311 * t226 + t241 * t388 - t267 * t390;
t66 = -qJ(4) * t174 - qJD(4) * t251 + t328;
t24 = -t301 * t66 + t302 * t60;
t12 = pkin(4) * t243 - pkin(9) * t104 + t24;
t25 = t301 * t60 + t302 * t66;
t16 = pkin(9) * t103 + t25;
t364 = t310 * t12 - t305 * t16;
t320 = -qJD(5) * t27 + t364;
t2 = pkin(5) * t240 - pkin(10) * t35 + t320;
t358 = -t305 * t12 - t310 * t16 - t61 * t385 + t64 * t386;
t3 = -pkin(10) * t36 - t358;
t375 = t309 * t2 - t304 * t3;
t436 = g(3) * t307;
t266 = -qJD(2) * pkin(2) + pkin(7) * t397;
t196 = pkin(3) * t251 + qJD(4) + t266;
t133 = -pkin(4) * t341 + t196;
t73 = -pkin(5) * t124 + t133;
t482 = -g(1) * t199 + g(2) * t197 + t276 * t436 - t73 * t53 + t375;
t481 = pkin(10) * t124;
t21 = t27 + t481;
t17 = t21 * t384;
t198 = t276 * t313 - t277 * t413;
t200 = t276 * t308 + t277 * t411;
t459 = g(1) * t200 - g(2) * t198 + t277 * t436 - t73 * t57 + t17;
t245 = t301 * t311 + t302 * t306;
t329 = t245 * t312;
t469 = qJD(1) * t329 - t245 * qJD(3);
t338 = t301 * t306 - t302 * t311;
t478 = t279 * t338;
t422 = t124 * t269;
t480 = t35 + t422;
t479 = t457 * t124;
t407 = t483 * t301 + t484 * t302;
t406 = t484 * t301 - t483 * t302;
t477 = -t124 ^ 2 + t457 ^ 2;
t282 = sin(t296);
t283 = cos(t296);
t204 = t282 * t313 - t283 * t413;
t206 = t282 * t308 + t283 * t411;
t476 = g(1) * t206 - g(2) * t204 - t133 * t124 + t283 * t436 + t358;
t474 = pkin(10) * t457;
t473 = -pkin(4) * t397 - pkin(9) * t478 + t407;
t464 = pkin(9) * t469 + t406;
t423 = t457 * t269;
t472 = -t36 - t423;
t288 = pkin(7) * t380;
t227 = -qJDD(2) * pkin(2) + pkin(7) * t367 + t288;
t349 = g(1) * t313 + g(2) * t308;
t435 = g(3) * t312;
t324 = t307 * t349 - t435;
t468 = qJD(3) * pkin(8) * t279 - t227 + t324;
t203 = t282 * t413 + t283 * t313;
t205 = -t282 * t411 + t283 * t308;
t467 = -g(1) * t205 + g(2) * t203 - t133 * t457 + t282 * t436 + t320;
t342 = -t245 * t305 - t310 * t338;
t410 = qJD(5) * t342 + t469 * t305 + t478 * t310;
t179 = t245 * t310 - t305 * t338;
t409 = qJD(5) * t179 + t478 * t305 - t469 * t310;
t443 = pkin(3) * t306;
t463 = pkin(3) * t390 - t396 * t443 - t290;
t26 = -t305 * t64 + t310 * t61;
t20 = t26 - t474;
t14 = -pkin(5) * t269 + t20;
t426 = t309 * t21;
t5 = t304 * t14 + t426;
t458 = -qJD(6) * t5 + t482;
t454 = t473 * t310;
t263 = t303 * t306;
t264 = t303 * t311;
t194 = t302 * t263 + t264 * t301;
t165 = -pkin(9) * t245 + t194;
t195 = t301 * t263 - t302 * t264;
t166 = -pkin(9) * t338 + t195;
t408 = t305 * t165 + t310 * t166;
t405 = -t469 * pkin(4) + t463;
t284 = pkin(3) * t302 + pkin(4);
t444 = pkin(3) * t301;
t353 = t310 * t284 - t305 * t444;
t89 = -t151 * t301 - t417;
t74 = t89 - t455;
t90 = t302 * t151 - t145;
t75 = t90 + t466;
t451 = t353 * qJD(5) - t305 * t74 - t310 * t75;
t223 = t284 * t305 + t310 * t444;
t450 = t223 * qJD(5) - t305 * t75 + t310 * t74;
t449 = t165 * t385 - t166 * t386 + t473 * t305 + t464 * t310;
t232 = t306 * t413 + t311 * t313;
t234 = -t306 * t411 + t308 * t311;
t446 = -g(1) * t234 + g(2) * t232;
t441 = pkin(7) * t306;
t118 = t179 * t304 - t309 * t342;
t434 = -qJD(6) * t118 - t409 * t304 + t410 * t309;
t119 = t179 * t309 + t304 * t342;
t433 = qJD(6) * t119 + t410 * t304 + t409 * t309;
t430 = t409 * pkin(5) + t405;
t247 = t311 * t261;
t187 = -qJ(4) * t414 + t247 + (-pkin(3) - t441) * t312;
t281 = pkin(7) * t412;
t401 = t306 * t261 + t281;
t416 = t306 * t307;
t193 = -qJ(4) * t416 + t401;
t127 = t302 * t187 - t193 * t301;
t220 = t338 * t307;
t96 = -pkin(4) * t312 + pkin(9) * t220 + t127;
t128 = t301 * t187 + t302 * t193;
t219 = t245 * t307;
t99 = -pkin(9) * t219 + t128;
t429 = t305 * t96 + t310 * t99;
t425 = t451 + t474;
t424 = t450 - t481;
t421 = t173 * t306;
t420 = t251 * t279;
t419 = t253 * t279;
t418 = t253 * t311;
t392 = qJD(2) * t307;
t403 = t311 * t255 + t392 * t441;
t115 = -t307 * t387 + t337 * qJD(2) + (-t281 + (qJ(4) * t307 - t261) * t306) * qJD(3) + t403;
t404 = t306 * t255 + t261 * t388;
t126 = (-pkin(7) * qJD(2) - qJ(4) * qJD(3)) * t414 + (-qJD(4) * t307 + (-pkin(7) * qJD(3) - qJ(4) * qJD(2)) * t312) * t306 + t404;
t51 = t301 * t115 + t302 * t126;
t399 = pkin(3) * t416 + t307 * pkin(7);
t299 = t307 ^ 2;
t398 = -t312 ^ 2 + t299;
t395 = qJD(2) * t251;
t394 = qJD(2) * t253;
t391 = qJD(2) * t312;
t372 = t306 * t391;
t377 = pkin(3) * t372 + pkin(7) * t391 + t388 * t442;
t287 = pkin(3) * t311 + pkin(2);
t376 = pkin(7) + t443;
t373 = t279 * t382;
t371 = t312 * t382;
t370 = t279 * t390;
t369 = t279 * t388;
t365 = qJD(6) * t14 + t3;
t160 = t245 * t389 + t301 * t372 - t302 * t371;
t50 = t302 * t115 - t126 * t301;
t43 = pkin(4) * t392 + pkin(9) * t160 + t50;
t159 = -qJD(2) * t329 + t338 * t389;
t46 = pkin(9) * t159 + t51;
t362 = -t305 * t46 + t310 * t43;
t360 = -t305 * t99 + t310 * t96;
t356 = t310 * t165 - t166 * t305;
t355 = -qJD(3) * t241 - t226;
t188 = pkin(4) * t219 + t399;
t76 = -pkin(10) * t179 + t356;
t352 = -t409 * pkin(10) + qJD(6) * t76 + t449;
t77 = pkin(10) * t342 + t408;
t351 = pkin(5) * t397 + t410 * pkin(10) + t408 * qJD(5) + qJD(6) * t77 + t464 * t305 - t454;
t142 = pkin(3) * t253 - pkin(4) * t180;
t348 = g(1) * t308 - g(2) * t313;
t347 = t267 * t388 - t186;
t345 = -pkin(8) * t243 + qJD(3) * t266;
t158 = -t219 * t305 - t220 * t310;
t343 = -t310 * t219 + t220 * t305;
t92 = t158 * t304 - t309 * t343;
t93 = t158 * t309 + t304 * t343;
t340 = t287 * t312 - t303 * t307;
t210 = pkin(4) * t338 - t287;
t129 = -pkin(4) * t159 + t377;
t335 = pkin(1) + t340;
t334 = -0.2e1 * pkin(1) * t381 - pkin(7) * qJDD(2);
t333 = t305 * t43 + t310 * t46 + t96 * t385 - t386 * t99;
t332 = t306 * t243 - t369;
t331 = t311 * t243 + t370;
t315 = qJD(1) ^ 2;
t326 = pkin(1) * t315 + t349;
t138 = pkin(3) * t174 + qJDD(4) + t227;
t314 = qJD(2) ^ 2;
t321 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t314 + t348;
t67 = -pkin(4) * t103 + t138;
t235 = t306 * t308 + t311 * t411;
t233 = t306 * t313 - t308 * t412;
t229 = qJDD(6) + t240;
t221 = pkin(5) + t353;
t139 = -pkin(5) * t342 + t210;
t111 = -pkin(5) * t343 + t188;
t78 = pkin(5) * t457 + t142;
t72 = qJD(5) * t158 - t310 * t159 - t160 * t305;
t71 = qJD(5) * t343 + t159 * t305 - t160 * t310;
t44 = pkin(5) * t72 + t129;
t40 = pkin(10) * t343 + t429;
t37 = -pkin(5) * t312 - pkin(10) * t158 + t360;
t19 = qJD(6) * t93 + t304 * t71 + t309 * t72;
t18 = -qJD(6) * t92 - t304 * t72 + t309 * t71;
t15 = pkin(5) * t36 + t67;
t7 = -pkin(10) * t72 + t333;
t6 = pkin(5) * t392 - pkin(10) * t71 - qJD(5) * t429 + t362;
t4 = t14 * t309 - t21 * t304;
t1 = [qJDD(1), t348, t349, qJDD(1) * t299 + 0.2e1 * t307 * t367, 0.2e1 * t295 * t307 - 0.2e1 * t381 * t398, qJDD(2) * t307 + t312 * t314, qJDD(2) * t312 - t307 * t314, 0, t307 * t334 + t312 * t321, -t307 * t321 + t312 * t334, t173 * t414 + (-t306 * t389 + t371) * t253 (-t251 * t311 - t253 * t306) * t391 + (-t421 - t174 * t311 + (t251 * t306 - t418) * qJD(3)) * t307 (-t173 - t373) * t312 + (t331 + t394) * t307 (t279 * t393 + t174) * t312 + (-t332 - t395) * t307, -t243 * t312 - t279 * t392 -(-t261 * t390 + t403) * t279 + t247 * t243 - g(1) * t233 - g(2) * t235 + ((t369 + t395) * pkin(7) + (-pkin(7) * t243 + qJD(2) * t266 - t355) * t306 + t347) * t312 + (pkin(7) * t174 + qJD(2) * t190 + t227 * t306 + t266 * t388) * t307, t404 * t279 - t401 * t243 - g(1) * t232 - g(2) * t234 + (t266 * t382 + (-t370 + t394) * pkin(7) + t328) * t312 + (-t266 * t390 - t191 * qJD(2) + t227 * t311 + (t173 - t373) * pkin(7)) * t307, t103 * t128 - t104 * t127 + t159 * t82 + t160 * t81 + t180 * t50 - t219 * t25 + t220 * t24 + t307 * t348 + t341 * t51, t25 * t128 + t82 * t51 + t24 * t127 + t81 * t50 + t138 * t399 + t196 * t377 + (-g(1) * t376 - g(2) * t335) * t313 + (g(1) * t335 - g(2) * t376) * t308, t158 * t35 + t457 * t71, t124 * t71 - t158 * t36 + t343 * t35 - t457 * t72, t158 * t240 - t269 * t71 - t312 * t35 + t392 * t457, t124 * t392 + t240 * t343 + t269 * t72 + t312 * t36, -t240 * t312 - t269 * t392, -t362 * t269 + t360 * t240 - t364 * t312 + t26 * t392 - t129 * t124 + t188 * t36 - t67 * t343 + t133 * t72 - g(1) * t204 - g(2) * t206 + (t269 * t429 + t27 * t312) * qJD(5), -g(1) * t203 - g(2) * t205 + t129 * t457 + t133 * t71 + t67 * t158 + t188 * t35 - t240 * t429 + t269 * t333 - t27 * t392 - t312 * t358, t18 * t53 + t8 * t93, t18 * t57 - t19 * t53 + t319 * t93 - t8 * t92, -t18 * t262 + t229 * t93 - t312 * t8 + t392 * t53, t19 * t262 - t229 * t92 - t312 * t319 + t392 * t57, -t229 * t312 - t262 * t392 -(-t304 * t7 + t309 * t6) * t262 + (-t304 * t40 + t309 * t37) * t229 - t375 * t312 + t4 * t392 - t44 * t57 - t111 * t319 + t15 * t92 + t73 * t19 - g(1) * t198 - g(2) * t200 + (-(-t304 * t37 - t309 * t40) * t262 + t5 * t312) * qJD(6), -t5 * t392 - g(1) * t197 - g(2) * t199 + t111 * t8 + t15 * t93 - t17 * t312 + t73 * t18 + t44 * t53 + ((-qJD(6) * t40 + t6) * t262 - t37 * t229 + t2 * t312) * t304 + ((qJD(6) * t37 + t7) * t262 - t40 * t229 + t365 * t312) * t309; 0, 0, 0, -t307 * t315 * t312, t398 * t315, t380, t295, qJDD(2), t307 * t326 - t288 - t435, t436 + (-pkin(7) * qJDD(1) + t326) * t312, -t279 * t418 + t421 (t173 + t420) * t311 + (-t174 + t419) * t306 (-t253 * t307 + t279 * t412) * qJD(1) + t332 (t251 * t307 - t279 * t415) * qJD(1) + t331, t279 * t397, -pkin(2) * t174 + t402 * t279 + t345 * t306 + (-t190 * t307 + (-pkin(7) * t251 - t266 * t306) * t312) * qJD(1) + t468 * t311, -pkin(2) * t173 - t236 * t279 + t345 * t311 + (-t266 * t412 + t191 * t307 + (-t253 * t312 + t279 * t414) * pkin(7)) * qJD(1) - t468 * t306, t103 * t195 - t104 * t194 + t407 * t180 - t24 * t245 - t338 * t25 - t349 * t312 + t406 * t341 + t469 * t82 - t478 * t81 - t436, t25 * t195 + t24 * t194 - t138 * t287 - g(3) * t340 + t406 * t82 + t407 * t81 + t463 * t196 + t349 * (t287 * t307 + t303 * t312) t179 * t35 + t410 * t457, t124 * t410 - t179 * t36 + t342 * t35 - t409 * t457, t179 * t240 - t269 * t410 - t397 * t457, -t124 * t397 + t240 * t342 + t269 * t409, t269 * t397, t356 * t240 + t210 * t36 - t67 * t342 - t26 * t397 + (t166 * t385 + (qJD(5) * t165 + t464) * t305 - t454) * t269 + t409 * t133 - t405 * t124 + t324 * t283, t410 * t133 + t67 * t179 + t210 * t35 - t408 * t240 + t269 * t449 + t27 * t397 - t282 * t324 + t405 * t457, t119 * t8 + t434 * t53, -t118 * t8 + t119 * t319 - t433 * t53 + t434 * t57, t119 * t229 - t262 * t434 - t397 * t53, -t118 * t229 + t262 * t433 - t397 * t57, t262 * t397 (-t304 * t77 + t309 * t76) * t229 - t139 * t319 + t15 * t118 - t4 * t397 + t433 * t73 - t430 * t57 + (t304 * t352 + t309 * t351) * t262 + t324 * t277 -(t304 * t76 + t309 * t77) * t229 + t139 * t8 + t15 * t119 + t5 * t397 + t434 * t73 + t430 * t53 + (-t304 * t351 + t309 * t352) * t262 - t324 * t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253 * t251, -t251 ^ 2 + t253 ^ 2, t173 - t420, -t174 - t419, t243, -t191 * t279 - t253 * t266 + (t355 + t436) * t306 - t347 + t446, g(1) * t235 - g(2) * t233 + g(3) * t414 - t190 * t279 + t251 * t266 - t328 (t103 * t301 - t104 * t302) * pkin(3) + (-t90 + t81) * t341 + (-t82 - t89) * t180, -t81 * t89 - t82 * t90 + (g(3) * t416 - t196 * t253 + t24 * t302 + t25 * t301 + t446) * pkin(3), -t479, t477, t480, t472, t240, t124 * t142 + t353 * t240 + t269 * t450 + t467, -t142 * t457 - t223 * t240 + t269 * t451 + t476, -t465, t461, t460, t456, t229 (t221 * t309 - t223 * t304) * t229 + t78 * t57 + (t425 * t304 + t424 * t309) * t262 + (-(-t221 * t304 - t223 * t309) * t262 - t5) * qJD(6) + t482, -t78 * t53 + (-t221 * t229 - t2 + (-qJD(6) * t223 - t424) * t262) * t304 + (-t223 * t229 + (qJD(6) * t221 + t425) * t262 - t365) * t309 + t459; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t180 ^ 2 - t341 ^ 2, -t180 * t81 - t341 * t82 + t138 - t324, 0, 0, 0, 0, 0, t36 - t423, t35 - t422, 0, 0, 0, 0, 0, -t319 - t428, t8 - t427; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t479, t477, t480, t472, t240, -t269 * t27 + t467, -t26 * t269 + t476, -t465, t461, t460, t456, t229 (-t20 * t304 - t426) * t262 + (t229 * t309 + t262 * t384 + t457 * t57) * pkin(5) + t458 (t21 * t262 - t2) * t304 + (-t20 * t262 - t365) * t309 + (-t229 * t304 + t262 * t383 - t457 * t53) * pkin(5) + t459; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t465, t461, t460, t456, t229, -t262 * t5 + t458, -t304 * t2 - t262 * t4 - t309 * t365 + t459;];
tau_reg  = t1;
