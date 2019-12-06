% Calculate vector of inverse dynamics joint torques for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRPRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RRPRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:36:20
% EndTime: 2019-12-05 18:36:25
% DurationCPUTime: 3.42s
% Computational Cost: add. (2533->341), mult. (3859->467), div. (0->0), fcn. (2535->14), ass. (0->192)
t385 = qJDD(1) + qJDD(2);
t388 = qJD(1) + qJD(2);
t393 = sin(pkin(9));
t399 = cos(qJ(5));
t400 = cos(qJ(4));
t464 = qJD(4) + qJD(5);
t437 = t464 * t400;
t395 = sin(qJ(5));
t491 = t395 * t400;
t396 = sin(qJ(4));
t503 = t385 * t396;
t470 = qJD(4) * t396;
t448 = t393 * t470;
t492 = t395 * t396;
t454 = t393 * t492;
t520 = -qJD(5) * t454 - t395 * t448;
t271 = t520 * t388 + (t385 * t491 + (t388 * t437 + t503) * t399) * t393;
t397 = sin(qJ(2));
t509 = pkin(1) * qJD(1);
t461 = t397 * t509;
t346 = qJ(3) * t388 + t461;
t386 = t393 ^ 2;
t394 = cos(pkin(9));
t499 = t388 * t394;
t356 = -qJD(4) + t499;
t466 = qJD(4) + t356;
t521 = t346 * (t386 * t388 + t394 * t466);
t401 = cos(qJ(2));
t474 = qJD(1) * t401;
t431 = -pkin(1) * t474 + qJD(3);
t500 = t388 * t393;
t473 = qJD(2) * t397;
t459 = pkin(1) * t473;
t512 = pkin(1) * t401;
t481 = -qJD(1) * t459 + qJDD(1) * t512;
t444 = qJDD(3) - t481;
t511 = pkin(2) * t385;
t318 = t444 - t511;
t392 = qJ(1) + qJ(2);
t381 = sin(t392);
t383 = cos(t392);
t430 = g(2) * t383 + g(3) * t381;
t517 = -t318 + t430;
t349 = -pkin(3) * t394 - pkin(7) * t393 - pkin(2);
t469 = qJD(4) * t400;
t471 = qJD(3) * t394;
t493 = t394 * t401;
t516 = -(t396 * t397 + t400 * t493) * t509 + t349 * t469 + t400 * t471;
t333 = t349 - t512;
t372 = pkin(1) * t397 + qJ(3);
t494 = t394 * t400;
t485 = t396 * t333 + t372 * t494;
t515 = qJ(3) * t494 + t396 * t349;
t513 = qJD(4) * t515 + (-t396 * t493 + t397 * t400) * t509 + t396 * t471;
t510 = g(1) * t393;
t375 = g(2) * t381;
t303 = t349 * t388 + t431;
t457 = t346 * t494;
t418 = -t303 * t396 - t457;
t497 = t393 * t396;
t463 = pkin(8) * t497;
t284 = -t388 * t463 - t418;
t508 = t284 * t399;
t348 = -qJD(5) + t356;
t507 = t348 * t394;
t506 = t381 * t394;
t505 = t383 * t394;
t504 = t385 * t394;
t502 = t385 * t400;
t465 = qJDD(1) * t397;
t472 = qJD(2) * t401;
t311 = qJ(3) * t385 + qJD(3) * t388 + (qJD(1) * t472 + t465) * pkin(1);
t304 = t386 * t311;
t498 = t388 * t396;
t496 = t393 * t400;
t495 = t394 * t396;
t353 = -qJDD(4) + t504;
t490 = t396 * t353;
t489 = t396 * t400;
t488 = t399 * t400;
t484 = g(2) * t505 + g(3) * t506;
t482 = g(3) * t383 - t375;
t387 = t394 ^ 2;
t480 = t386 + t387;
t390 = t400 ^ 2;
t479 = t396 ^ 2 - t390;
t478 = MDP(11) * t386;
t477 = MDP(12) * t386;
t476 = MDP(13) * t393;
t475 = MDP(14) * t393;
t468 = qJD(5) * t395;
t347 = -qJDD(5) + t353;
t340 = t347 * MDP(22);
t467 = t353 * MDP(15);
t462 = pkin(8) * t496;
t458 = qJ(3) * t495;
t455 = t388 * t496;
t453 = t393 * t488;
t295 = t349 * t385 + t444;
t452 = t396 * t295 + t303 * t469 + t311 * t494;
t366 = pkin(1) * t472 + qJD(3);
t451 = t333 * t469 + t366 * t494 + t396 * t459;
t450 = t388 * t473;
t449 = t388 * t469;
t447 = t394 * t470;
t446 = t346 * t470;
t443 = t366 * t480;
t442 = t480 * t311;
t441 = t480 * t385;
t440 = t353 + t504;
t439 = t388 * t466;
t411 = -t394 * t446 + t452;
t417 = t449 + t503;
t265 = -pkin(8) * t393 * t417 + t411;
t298 = t400 * t303;
t283 = -pkin(8) * t455 - t346 * t495 + t298;
t273 = -pkin(4) * t356 + t283;
t438 = qJD(5) * t273 + t265;
t436 = t388 * t461;
t435 = t387 * t311 + t304 - t482;
t434 = -t481 - t430;
t292 = t400 * t295;
t429 = -t311 * t495 + t292;
t337 = t400 * t349;
t296 = -t462 + t337 + (-qJ(3) * t396 - pkin(4)) * t394;
t428 = qJD(5) * t296 + (-t458 - t462) * qJD(4) + t516;
t302 = -t463 + t515;
t361 = pkin(8) * t448;
t427 = qJD(5) * t302 - t361 + t513;
t426 = -t273 * t395 - t508;
t330 = t400 * t333;
t286 = -t462 + t330 + (-t372 * t396 - pkin(4)) * t394;
t294 = -t463 + t485;
t425 = t286 * t399 - t294 * t395;
t424 = t286 * t395 + t294 * t399;
t423 = t396 * t399 + t491;
t422 = -t488 + t492;
t421 = qJD(4) * (t356 + t499);
t362 = t393 * pkin(4) * t469;
t420 = -t431 * t393 - t362;
t322 = t381 * t495 + t383 * t400;
t324 = -t381 * t400 + t383 * t495;
t416 = -g(2) * t324 - g(3) * t322 + t400 * t304 + t411 * t394;
t323 = t381 * t494 - t383 * t396;
t325 = -t381 * t396 - t383 * t494;
t415 = t386 * t346 * t469 - g(2) * t325 + g(3) * t323 + t396 * t304;
t414 = t436 + t511;
t377 = -pkin(2) - t512;
t413 = -pkin(1) * t450 - t377 * t385;
t264 = -t385 * t462 - pkin(4) * t353 + (-t457 + (pkin(8) * t500 - t303) * t396) * qJD(4) + t429;
t278 = t284 * t468;
t285 = (pkin(4) * t417 + t311) * t393;
t406 = t464 * t423;
t287 = t406 * t393;
t306 = (pkin(4) * t498 + t346) * t393;
t391 = qJ(4) + qJ(5);
t380 = sin(t391);
t382 = cos(t391);
t313 = t380 * t506 + t382 * t383;
t315 = t380 * t505 - t381 * t382;
t328 = t422 * t393;
t412 = -g(2) * t315 - g(3) * t313 + (t395 * t264 + t438 * t399 - t278) * t394 - t285 * t328 - t306 * t287;
t270 = t385 * t453 + (-t385 * t492 - t388 * t406) * t393;
t308 = (-t453 + t454) * t388;
t309 = t423 * t500;
t410 = -t308 * t309 * MDP(18) + (-t309 * t348 + t270) * MDP(20) + (t308 * t348 - t271) * MDP(21) + (t308 ^ 2 - t309 ^ 2) * MDP(19) - t340;
t384 = t388 ^ 2;
t409 = -t356 ^ 2 - t384 * t386;
t258 = qJD(5) * t426 + t399 * t264 - t395 * t265;
t288 = t393 * t399 * t437 + t520;
t314 = -t383 * t380 + t382 * t506;
t316 = -t380 * t381 - t382 * t505;
t327 = t423 * t393;
t408 = -g(2) * t316 + g(3) * t314 - t258 * t394 + t285 * t327 + t306 * t288;
t407 = t431 * t480;
t405 = t278 + t382 * t510 + (t284 * t348 - t264) * t395 - g(2) * t314 - g(3) * t316 + t306 * t309;
t404 = (-t270 * t327 + t271 * t328 + t287 * t309 + t288 * t308) * MDP(19) + (-t270 * t394 + t287 * t348 + t328 * t347) * MDP(20) + (t271 * t394 + t288 * t348 + t327 * t347) * MDP(21) + (-t270 * t328 + t287 * t308) * MDP(18) + (t396 * t421 - t400 * t440) * t476 + (t396 * t440 + t400 * t421) * t475 + 0.2e1 * (qJD(4) * t388 * t479 - t385 * t489) * t477 + (t385 * t390 - 0.2e1 * t396 * t449) * t478 + t385 * MDP(4) + (t340 + t467) * t394;
t403 = -g(2) * t313 + g(3) * t315 + t306 * t308 + t380 * t510 + t258;
t402 = cos(qJ(1));
t398 = sin(qJ(1));
t371 = t383 * qJ(3);
t369 = pkin(4) * t497;
t364 = t400 * t459;
t342 = qJ(3) * t393 + t369;
t341 = -pkin(2) * t388 + t431;
t331 = t372 * t393 + t369;
t317 = t366 * t393 + t362;
t312 = t318 * t393;
t280 = -t485 * qJD(4) - t366 * t495 + t361 + t364;
t279 = (-t372 * t495 - t462) * qJD(4) + t451;
t267 = qJD(4) * t418 + t429;
t1 = [(t318 * t377 + t341 * t459 - g(2) * (-pkin(1) * t402 - t383 * pkin(2) - t381 * qJ(3)) - g(3) * (-pkin(1) * t398 - pkin(2) * t381 + t371) + t346 * t443 + t372 * t442) * MDP(10) + ((-t318 + t413) * t394 + t484) * MDP(7) + (t312 + (-t413 - t430) * t393) * MDP(8) + (((-qJDD(1) - t385) * t397 + (-qJD(1) - t388) * t472) * pkin(1) + t482) * MDP(6) + qJDD(1) * MDP(1) + (t372 * t441 + t388 * t443 + t435) * MDP(9) + (-(-t333 * t470 + t364) * t356 - t330 * t353 + (-(-t366 * t396 - t372 * t469) * t356 + t372 * t490 - t267) * t394 + (t366 * t498 + t372 * t417) * t386 + t415) * MDP(16) + (-(-qJD(5) * t424 - t279 * t395 + t280 * t399) * t348 - t425 * t347 + t317 * t309 + t331 * t271 + t408) * MDP(23) + ((-t372 * t447 + t451) * t356 + t485 * t353 + ((t366 * t388 + t372 * t385) * t400 + (-t372 * t388 - t346) * t470) * t386 + t416) * MDP(17) + ((t385 * t401 - t450) * pkin(1) - t434) * MDP(5) + t404 + ((qJD(5) * t425 + t279 * t399 + t280 * t395) * t348 + t424 * t347 - t317 * t308 + t331 * t270 + t412) * MDP(24) + (g(2) * t402 + g(3) * t398) * MDP(2) + (-g(2) * t398 + g(3) * t402) * MDP(3); (-t341 * t461 - g(3) * t371 + t517 * pkin(2) + (t442 + t375) * qJ(3) + t407 * t346) * MDP(10) + (-(t337 - t458) * t353 - t267 * t394 + t513 * t356 + (qJ(3) * t503 + (qJ(3) * t469 + t396 * t431) * t388) * t386 + t415) * MDP(16) + (t515 * t353 + (-qJ(3) * t447 + t516) * t356 + (qJ(3) * t502 - t446 + (-qJ(3) * t470 + t400 * t431) * t388) * t386 + t416) * MDP(17) + ((-t465 + (-qJD(2) + t388) * t474) * pkin(1) + t482) * MDP(6) + ((-t318 + t414) * t394 + t484) * MDP(7) + (t312 + (-t414 - t430) * t393) * MDP(8) + ((t296 * t395 + t302 * t399) * t347 + t342 * t270 + (-t395 * t427 + t399 * t428) * t348 + t420 * t308 + t412) * MDP(24) + (qJ(3) * t441 + t388 * t407 + t435) * MDP(9) + (-(t296 * t399 - t302 * t395) * t347 + t342 * t271 + (t395 * t428 + t399 * t427) * t348 - t420 * t309 + t408) * MDP(23) + t404 + (-t434 + t436) * MDP(5); -MDP(7) * t504 + t393 * t385 * MDP(8) - t480 * t384 * MDP(9) + (-t346 * t388 * t480 - t517) * MDP(10) + (-t353 * t400 + t396 * t409) * MDP(16) + (t400 * t409 + t490) * MDP(17) + (t406 * t348 + t422 * t347 + (-t393 * t309 - t423 * t507) * t388) * MDP(23) + (t308 * t500 + t423 * t347 + (-t464 * t348 + t507 * t388) * t422) * MDP(24); (-t396 * t439 + t502) * t476 + (-t400 * t439 - t503) * t475 - t467 + (-g(2) * t322 + g(3) * t324 + t292 - t400 * t521 + (-t303 * t466 - t394 * t311 + t510) * t396) * MDP(16) + (g(1) * t496 - g(2) * t323 - g(3) * t325 - t298 * t356 + t396 * t521 - t452) * MDP(17) + ((-t283 * t395 - t508) * t348 + (-t309 * t455 - t399 * t347 + t348 * t468) * pkin(4) + t403) * MDP(23) + ((-t283 * t348 - t438) * t399 + (qJD(5) * t399 * t348 + t308 * t455 + t395 * t347) * pkin(4) + t405) * MDP(24) + t410 + (-t477 * t479 + t478 * t489) * t384; (t348 * t426 + t403) * MDP(23) + ((-t265 + (-qJD(5) - t348) * t273) * t399 + t405) * MDP(24) + t410;];
tau = t1;
