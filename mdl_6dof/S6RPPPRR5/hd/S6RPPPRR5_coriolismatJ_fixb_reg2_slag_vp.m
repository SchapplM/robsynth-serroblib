% Calculate inertial parameters regressor of coriolis matrix for
% S6RPPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPPPRR5_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:48
% EndTime: 2019-03-09 01:37:54
% DurationCPUTime: 5.08s
% Computational Cost: add. (4251->390), mult. (7394->545), div. (0->0), fcn. (7024->6), ass. (0->298)
t292 = sin(qJ(6));
t476 = 0.2e1 * t292;
t294 = cos(qJ(6));
t289 = sin(pkin(9));
t295 = cos(qJ(5));
t421 = t295 * t289;
t242 = t294 * t421;
t290 = cos(pkin(9));
t428 = t292 * t290;
t349 = -t242 + t428;
t475 = t349 / 0.2e1;
t474 = t349 * t295;
t422 = t294 * t290;
t202 = t292 * t421 + t422;
t473 = t202 / 0.2e1;
t420 = t295 * t290;
t429 = t292 * t289;
t204 = t294 * t420 + t429;
t472 = -t204 / 0.2e1;
t471 = -t242 / 0.2e1;
t273 = t294 * t289;
t357 = -t273 / 0.2e1;
t285 = t292 ^ 2;
t470 = t285 / 0.2e1;
t293 = sin(qJ(5));
t286 = t293 ^ 2;
t469 = t286 / 0.2e1;
t288 = t295 ^ 2;
t468 = -t288 / 0.2e1;
t467 = t288 / 0.2e1;
t466 = t292 / 0.2e1;
t465 = t294 / 0.2e1;
t464 = t293 * pkin(5);
t463 = t295 * pkin(5);
t462 = t295 * pkin(8);
t461 = pkin(3) + qJ(2);
t203 = t292 * t420 - t273;
t362 = t203 * t466;
t318 = t204 * t465 + t362;
t336 = t349 * t294;
t319 = -t336 / 0.2e1;
t359 = t428 / 0.2e1;
t374 = t293 * t420;
t21 = -t289 * t374 + (t202 * t359 + t289 * t318 + t290 * t319) * t293;
t460 = t21 * qJD(5);
t303 = t202 * t466 + t319;
t287 = t294 ^ 2;
t306 = t467 + (t287 / 0.2e1 + t470 - 0.1e1 / 0.2e1) * t286;
t36 = t289 * t306 - t295 * t303;
t459 = t36 * qJD(5);
t37 = t290 * t306 - t295 * t318;
t458 = t37 * qJD(5);
t291 = pkin(1) + qJ(3);
t316 = t289 * t291 + t290 * t461;
t205 = -pkin(4) - t316;
t339 = -t293 * pkin(8) - t463;
t297 = t205 + t339;
t210 = t289 * t461 - t290 * t291;
t206 = pkin(7) + t210;
t437 = t206 * t295;
t377 = t292 * t437;
t79 = -t294 * t297 + t377;
t457 = t79 * t295;
t425 = t294 * t206;
t375 = t295 * t425;
t80 = t292 * t297 + t375;
t456 = t80 * t295;
t313 = -t202 ^ 2 / 0.2e1 - t349 ^ 2 / 0.2e1;
t282 = t290 ^ 2;
t256 = -t286 * t282 / 0.2e1;
t331 = t256 - t203 ^ 2 / 0.2e1 - t204 ^ 2 / 0.2e1;
t281 = t289 ^ 2;
t432 = t286 * t281;
t20 = -t432 / 0.2e1 + t313 + t331;
t455 = qJD(1) * t20;
t433 = t286 * t206;
t38 = -t292 * t433 - t457;
t454 = qJD(1) * t38;
t423 = t294 * t286;
t39 = -t206 * t423 - t456;
t453 = qJD(1) * t39;
t412 = t286 + t288;
t212 = t412 * t289;
t50 = -t205 * t290 + t206 * t212;
t452 = qJD(1) * t50;
t213 = t412 * t290;
t51 = t205 * t289 + t206 * t213;
t451 = qJD(1) * t51;
t353 = -t421 / 0.2e1;
t58 = (t353 + t303) * t293;
t450 = qJD(1) * t58;
t351 = t420 / 0.2e1;
t59 = (t351 - t318) * t293;
t449 = qJD(1) * t59;
t337 = t349 * t292;
t443 = t202 * t294;
t93 = (-t337 - t443) * t293;
t448 = qJD(1) * t93;
t439 = t204 * t292;
t441 = t203 * t294;
t94 = (-t439 + t441) * t293;
t447 = qJD(1) * t94;
t376 = t293 * t425;
t240 = -t462 + t464;
t430 = t292 * t240;
t124 = -t376 + t430;
t444 = t124 * t293;
t427 = t292 * t293;
t190 = t206 * t427;
t424 = t294 * t240;
t123 = t190 + t424;
t445 = t123 * t293;
t12 = (t445 - t457) * t294 + (t444 + t456) * t292;
t446 = t12 * qJD(1);
t442 = t202 * t295;
t440 = t203 * t295;
t438 = t204 * t295;
t22 = t79 * t293 + (t123 - 0.2e1 * t190) * t295;
t436 = t22 * qJD(1);
t23 = t124 * t295 + (-t80 + 0.2e1 * t375) * t293;
t435 = t23 * qJD(1);
t434 = t285 * t295;
t431 = t289 * t290;
t426 = t293 * t295;
t320 = t337 / 0.2e1;
t31 = t295 * t320 + t362 + (t442 / 0.2e1 + t204 / 0.2e1) * t294;
t419 = t31 * qJD(1);
t307 = (-t439 / 0.2e1 + t441 / 0.2e1) * t295;
t33 = t307 - t303;
t418 = t33 * qJD(1);
t361 = -t429 / 0.2e1;
t382 = 0.1e1 / 0.2e1 + t469;
t81 = (t361 + t472) * t295 - t382 * t422;
t417 = t81 * qJD(1);
t340 = t428 + t471;
t82 = t286 * t357 + t340 * t295 + t357;
t416 = t82 * qJD(1);
t83 = t471 + t440 / 0.2e1 + t382 * t428;
t415 = t83 * qJD(1);
t241 = t286 * t429;
t354 = t422 / 0.2e1;
t360 = t429 / 0.2e1;
t84 = t360 + t241 / 0.2e1 + (t354 + t473) * t295;
t414 = t84 * qJD(1);
t413 = t285 + t287;
t265 = t287 - t285;
t266 = t288 - t286;
t127 = -t241 - t442;
t411 = qJD(1) * t127;
t378 = t286 * t428;
t128 = t378 + t440;
t410 = qJD(1) * t128;
t129 = t273 * t286 - t474;
t409 = qJD(1) * t129;
t130 = -t286 * t422 - t438;
t408 = qJD(1) * t130;
t407 = qJD(1) * t295;
t406 = qJD(5) * t292;
t405 = qJD(5) * t294;
t404 = qJD(6) * t292;
t403 = qJD(6) * t294;
t402 = qJD(6) * t295;
t352 = t421 / 0.2e1;
t326 = t292 * t352 - t202 / 0.2e1;
t355 = -t422 / 0.2e1;
t100 = (t355 + t326) * t293;
t401 = t100 * qJD(1);
t342 = t294 * t352;
t101 = (t342 + t340) * t293;
t400 = t101 * qJD(1);
t345 = t292 * t351;
t325 = t345 - t203 / 0.2e1;
t356 = t273 / 0.2e1;
t103 = (t356 + t325) * t293;
t399 = t103 * qJD(1);
t341 = t294 * t351;
t324 = t341 + t472;
t106 = (t361 + t324) * t293;
t398 = t106 * qJD(1);
t254 = t432 / 0.2e1;
t381 = t467 + 0.1e1 / 0.2e1;
t113 = t254 + t381 * t281 + (t469 + t381) * t282;
t397 = t113 * qJD(1);
t115 = t210 * t289 + t290 * t316;
t396 = t115 * qJD(1);
t116 = t210 * t290 - t289 * t316;
t395 = t116 * qJD(1);
t207 = (t470 - t287 / 0.2e1) * t293;
t394 = t207 * qJD(6);
t393 = t212 * qJD(1);
t392 = t213 * qJD(1);
t219 = t266 * t292;
t391 = t219 * qJD(1);
t220 = t288 * t294 - t423;
t390 = t220 * qJD(1);
t251 = t282 + t281;
t389 = t251 * qJD(1);
t388 = t266 * qJD(1);
t279 = t289 * qJD(1);
t280 = t290 * qJD(1);
t387 = t291 * qJD(1);
t386 = t293 * qJD(1);
t385 = t293 * qJD(5);
t384 = t293 * qJD(6);
t383 = t295 * qJD(5);
t380 = t290 * pkin(5) / 0.2e1;
t379 = -t464 / 0.2e1;
t373 = t294 * t386;
t372 = t292 * t405;
t371 = t292 * t402;
t370 = t294 * t402;
t269 = t292 * t385;
t369 = t292 * t403;
t368 = t293 * t383;
t367 = t289 * t386;
t366 = t290 * t386;
t365 = t290 * t385;
t270 = t295 * t386;
t364 = t294 * t385;
t363 = t295 * t279;
t358 = -t427 / 0.2e1;
t350 = t413 * t293;
t348 = -qJD(6) + t407;
t347 = t292 * t364;
t346 = t286 * t369;
t344 = t293 * t357;
t343 = t290 * t358;
t338 = qJD(5) * t350;
t10 = (t444 / 0.2e1 + t456 / 0.2e1) * t294 + (-t445 / 0.2e1 + t457 / 0.2e1) * t292 + (t469 + t468) * t206;
t13 = t206 ^ 2 * t426 - t79 * t123 + t80 * t124;
t335 = t13 * qJD(1) + t10 * qJD(4);
t14 = t79 * t202 + t289 * t433 - t349 * t80;
t334 = t14 * qJD(1) + t58 * qJD(4);
t15 = -t203 * t79 - t204 * t80 - t290 * t433;
t333 = t15 * qJD(1) + t59 * qJD(4);
t332 = -t123 * t292 + t124 * t294;
t330 = t348 * t293;
t329 = qJD(2) * t290 + qJD(3) * t289;
t328 = -qJD(2) * t289 + qJD(3) * t290;
t327 = t462 / 0.2e1 + t379;
t323 = t465 * t80 + t466 * t79;
t314 = -t240 / 0.2e1 + t327;
t125 = t314 * t292;
t322 = pkin(5) * t405 + qJD(1) * t125;
t126 = t314 * t294;
t321 = pkin(5) * t406 - qJD(1) * t126;
t317 = t294 * t330;
t142 = -qJD(1) * t207 + t372;
t315 = -t206 * t374 + t123 * t203 / 0.2e1 + t124 * t472;
t131 = qJD(1) * t292 * t423 + qJD(5) * t207;
t218 = t265 * t286;
t312 = qJD(1) * t218 + 0.2e1 * t347;
t311 = -qJD(5) * t265 + t373 * t476;
t310 = t292 * t384 - t294 * t383;
t309 = t292 * t383 + t294 * t384;
t308 = t318 * pkin(8);
t305 = -qJD(1) * t205 + t328;
t304 = qJD(5) * (-pkin(8) * t350 - t463);
t298 = t303 * pkin(8);
t2 = t298 + (-t289 * pkin(5) / 0.2e1 + t323 * t290) * t293 + t315;
t53 = t203 * t290 * t427 + t204 * t293 * t422 - t282 * t426;
t302 = -t2 * qJD(1) - t53 * qJD(2) - t21 * qJD(3) - t37 * qJD(4);
t300 = -t206 * t293 * t421 + t123 * t473 + t124 * t475;
t3 = -t308 + (t289 * t323 + t380) * t293 + t300;
t52 = -t281 * t426 + (t202 * t427 - t293 * t336) * t289;
t301 = -t3 * qJD(1) - t21 * qJD(2) - t52 * qJD(3) - t36 * qJD(4);
t143 = (-0.1e1 + t413) * t426;
t299 = t10 * qJD(1) - t37 * qJD(2) - t36 * qJD(3) + t143 * qJD(4);
t296 = t203 * t202 - t204 * t349 + t286 * t431;
t284 = qJ(2) * qJD(1);
t283 = qJ(2) * qJD(2);
t275 = t287 * t295;
t274 = t385 / 0.2e1;
t264 = t295 * t280;
t211 = t270 - t384 / 0.2e1;
t199 = t289 * t385 - t264;
t198 = t363 + t365;
t197 = t289 * t383 + t366;
t196 = -t290 * t383 + t367;
t138 = (-0.1e1 + t412) * t431;
t114 = t254 + t256 + (t468 + 0.1e1 / 0.2e1) * t282 + (-0.1e1 / 0.2e1 + t467) * t281;
t105 = (t324 + t360) * t293;
t104 = t293 * t325 + t344;
t102 = t343 + (t342 + t475) * t293;
t99 = (t326 + t354) * t293;
t88 = t438 / 0.2e1 + t286 * t354 + t292 * t353 + t355;
t87 = -t474 / 0.2e1 + t286 * t356 + t345 + t357;
t86 = -t440 / 0.2e1 - t378 / 0.2e1 + t471 + t359;
t85 = -t442 / 0.2e1 - t241 / 0.2e1 + t341 + t360;
t74 = t190 + t424 / 0.2e1 + t327 * t294;
t73 = t376 - t430 / 0.2e1 - t327 * t292;
t32 = t307 + t303;
t30 = (t320 + t443 / 0.2e1) * t295 - t318;
t19 = t254 - t313 + t331;
t5 = t293 * t355 * t80 + t289 * t379 + t343 * t79 + t298 - t315;
t4 = t289 * t358 * t79 + t293 * t380 + t344 * t80 - t300 - t308;
t1 = qJD(2) * t58 + qJD(3) * t59 + qJD(5) * t10;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t283, 0, 0, 0, 0, 0, 0, qJD(2), 0, qJD(3), qJD(3) * t291 + t283, 0, 0, 0, 0, 0, 0, t329, t328, 0, qJD(2) * t115 - qJD(3) * t116, t368, t266 * qJD(5), 0, -t368, 0, 0, t205 * t385 + t295 * t329, t205 * t383 - t293 * t329, qJD(2) * t212 - qJD(3) * t213, qJD(2) * t50 - qJD(3) * t51, t287 * t368 - t346, -t218 * qJD(6) - 0.2e1 * t295 * t347, -t220 * qJD(5) + t293 * t371, t285 * t368 + t346, t219 * qJD(5) + t293 * t370, -t368, -qJD(2) * t127 - qJD(3) * t128 - qJD(5) * t22 - qJD(6) * t39, qJD(2) * t129 + qJD(3) * t130 + qJD(5) * t23 + qJD(6) * t38, -qJD(2) * t93 - qJD(3) * t94 - qJD(5) * t12, qJD(2) * t14 + qJD(3) * t15 + qJD(5) * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), t284, 0, 0, 0, 0, 0, 0, qJD(1), 0, 0, t284, 0, 0, 0, 0, 0, 0, t280, -t279, 0, t396, 0, 0, 0, 0, 0, 0, t264, -t366, t393, qJD(2) * t138 + qJD(3) * t114 + t452, 0, 0, 0, 0, 0, 0, qJD(5) * t104 + qJD(6) * t88 - t411, qJD(5) * t105 + qJD(6) * t86 + t409, qJD(5) * t32 - t448, qJD(2) * t296 + t19 * qJD(3) + t5 * qJD(5) + t334; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), t387, 0, 0, 0, 0, 0, 0, t279, t280, 0, -t395, 0, 0, 0, 0, 0, 0, t363, -t367, -t392, qJD(2) * t114 - qJD(3) * t138 - t451, 0, 0, 0, 0, 0, 0, qJD(5) * t99 + qJD(6) * t87 - t410, qJD(5) * t102 + qJD(6) * t85 + t408, qJD(5) * t30 - t447, t19 * qJD(2) - qJD(3) * t296 + t4 * qJD(5) + t333; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t270, t388, t383, -t270, -t385, 0, t205 * t386 - t206 * t383, t205 * t407 + t206 * t385, 0, 0, -t394 + (t287 * t386 + t372) * t295 (t275 - t434) * qJD(5) + 0.2e1 * (-qJD(6) - t407) * t294 * t427, t269 - t390, t394 + (t285 * t386 - t372) * t295, t364 + t391, -t211, -t436 + t104 * qJD(2) + t99 * qJD(3) + (t292 * t339 - t375) * qJD(5) + t74 * qJD(6), t435 + t105 * qJD(2) + t102 * qJD(3) + (t294 * t339 + t377) * qJD(5) + t73 * qJD(6), t32 * qJD(2) + t30 * qJD(3) + qJD(5) * t332 - t446, t5 * qJD(2) + t4 * qJD(3) + (-pkin(5) * t437 + pkin(8) * t332) * qJD(5) + t335; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, -t312, t292 * t330, t131, t317, t274, qJD(2) * t88 + qJD(3) * t87 + qJD(5) * t74 - qJD(6) * t80 - t453, qJD(2) * t86 + qJD(3) * t85 + qJD(5) * t73 + qJD(6) * t79 + t454, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), -t284, 0, 0, 0, 0, 0, 0, -qJD(1), 0, 0, -qJD(3) - t284, 0, 0, 0, 0, 0, 0, -t280, t279, 0, -qJD(3) * t251 - t396, 0, 0, 0, 0, 0, 0, t199, t197, -t393, -qJD(3) * t113 - t452, 0, 0, 0, 0, 0, 0, qJD(5) * t103 - qJD(6) * t81 + t411, qJD(5) * t106 - qJD(6) * t83 - t409, qJD(5) * t33 + t448, qJD(3) * t20 - qJD(5) * t2 - t334; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t389, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t397, 0, 0, 0, 0, 0, 0, 0, 0, 0, t455 - t460; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t450 - t458; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, t198, 0, 0, 0, 0, 0, 0, 0, 0, t290 * t310 + t399, t290 * t309 + t398, -t290 * t338 + t418, t290 * t304 + t302; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(6) * t204 + t292 * t365 - t417, qJD(6) * t203 + t290 * t364 - t415, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), qJD(2) - t387, 0, 0, 0, 0, 0, 0, -t279, -t280, 0, qJD(2) * t251 + t395, 0, 0, 0, 0, 0, 0, -t198, t196, t392, qJD(2) * t113 + t451, 0, 0, 0, 0, 0, 0, qJD(5) * t100 - qJD(6) * t82 + t410, qJD(5) * t101 - qJD(6) * t84 - t408, qJD(5) * t31 + t447, -qJD(2) * t20 - qJD(5) * t3 - t333; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t389, 0, 0, 0, 0, 0, 0, 0, 0, 0, t397, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t455 - t460; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t449 - t459; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t197, t199, 0, 0, 0, 0, 0, 0, 0, 0, t289 * t310 + t401, t289 * t309 + t400, -t289 * t338 + t419, t289 * t304 + t301; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(6) * t349 + t269 * t289 - t416, qJD(6) * t202 + t289 * t364 - t414, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t450 - t458; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t449 - t459; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t385, -t383, 0, 0, 0, 0, 0, 0, 0, 0, -t364 - t371, t269 - t370 (t275 + t434) * qJD(5) (t413 * t462 - t464) * qJD(5) + t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t309, t310, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t270, -t388, 0, t270, 0, 0, t305 * t293, t305 * t295, 0, 0, -t270 * t287 - t394, t317 * t476, -t370 + t390, -t270 * t285 + t394, t371 - t391, t211, -qJD(2) * t103 - qJD(3) * t100 + qJD(6) * t126 + t436, -qJD(2) * t106 - qJD(3) * t101 - qJD(6) * t125 - t435, -qJD(2) * t33 - qJD(3) * t31 + t446, qJD(2) * t2 + qJD(3) * t3 - t335; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t367, -t363, 0, 0, 0, 0, 0, 0, 0, 0, -t399, -t398, -t418, -t302; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t366, t264, 0, 0, 0, 0, 0, 0, 0, 0, -t401, -t400, -t419, -t301; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t369, t265 * qJD(6), 0, -t369, 0, 0, -pkin(5) * t404, -pkin(5) * t403, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, -t311, -t348 * t294, -t142, t348 * t292, -t386 / 0.2e1, -pkin(8) * t403 - t321, pkin(8) * t404 - t322, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, t312 (-t292 * t386 + t405) * t295, -t131 (-t373 - t406) * t295, t274, qJD(2) * t81 + qJD(3) * t82 - qJD(5) * t126 + t453, qJD(2) * t83 + qJD(3) * t84 + qJD(5) * t125 - t454, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t417, t415, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t416, t414, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t142, t311, t294 * t407, t142, -t292 * t407, t386 / 0.2e1, t321, t322, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t6;
