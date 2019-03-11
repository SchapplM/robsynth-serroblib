% Calculate vector of inverse dynamics joint torques for
% S6RPRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6RPRPPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:39:45
% EndTime: 2019-03-09 02:39:54
% DurationCPUTime: 6.86s
% Computational Cost: add. (4498->447), mult. (9819->596), div. (0->0), fcn. (7112->18), ass. (0->196)
t451 = cos(qJ(3));
t526 = cos(pkin(10));
t491 = t526 * t451;
t414 = qJD(1) * t491;
t442 = sin(pkin(10));
t448 = sin(qJ(3));
t504 = qJD(1) * t448;
t388 = t442 * t504 - t414;
t385 = qJD(6) + t388;
t402 = t442 * t451 + t448 * t526;
t391 = t402 * qJD(1);
t441 = sin(pkin(11));
t444 = cos(pkin(11));
t378 = t444 * qJD(3) - t391 * t441;
t450 = cos(qJ(6));
t377 = qJD(3) * t441 + t391 * t444;
t447 = sin(qJ(6));
t519 = t377 * t447;
t542 = t378 * t450 - t519;
t543 = t385 * t542;
t432 = t451 * qJDD(2);
t443 = sin(pkin(9));
t420 = pkin(1) * t443 + pkin(7);
t408 = t420 * qJDD(1);
t459 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + qJD(2) * qJD(3) + t408;
t512 = qJ(4) + t420;
t487 = t512 * qJD(1);
t469 = t487 * qJD(3);
t320 = qJDD(3) * pkin(3) - t448 * t459 - t451 * t469 + t432;
t326 = (qJDD(2) - t469) * t448 + t459 * t451;
t290 = t320 * t526 - t442 * t326;
t289 = -qJDD(3) * pkin(4) + qJDD(5) - t290;
t437 = qJ(3) + pkin(10);
t427 = sin(t437);
t430 = cos(t437);
t438 = qJ(1) + pkin(9);
t428 = sin(t438);
t431 = cos(t438);
t485 = g(1) * t431 + g(2) * t428;
t461 = -g(3) * t430 + t427 * t485;
t538 = t289 - t461;
t471 = -t377 * t450 - t378 * t447;
t540 = t385 * t471;
t403 = t441 * t450 + t444 * t447;
t395 = t403 * qJD(6);
t507 = t403 * t388 + t395;
t493 = -g(1) * t428 + g(2) * t431;
t539 = t430 * t493;
t445 = cos(pkin(9));
t422 = -pkin(1) * t445 - pkin(2);
t434 = t451 * pkin(3);
t537 = t422 - t434;
t500 = qJD(1) * qJD(3);
t494 = t448 * t500;
t457 = pkin(3) * t494 + qJDD(1) * t537 + qJDD(4);
t390 = t402 * qJD(3);
t499 = qJDD(1) * t448;
t360 = qJD(1) * t390 - qJDD(1) * t491 + t442 * t499;
t359 = qJDD(6) + t360;
t401 = t441 * t447 - t450 * t444;
t508 = t385 * t401;
t536 = -t359 * t403 + t385 * t508;
t535 = pkin(8) * t444;
t531 = g(3) * t427;
t529 = g(3) * t451;
t416 = pkin(3) * t442 + qJ(5);
t528 = pkin(8) + t416;
t527 = qJD(3) * pkin(3);
t525 = qJ(5) * t427;
t524 = t542 * t391;
t523 = t471 * t391;
t521 = t360 * t441;
t520 = t360 * t444;
t518 = t388 * t441;
t466 = -t442 * t448 + t491;
t393 = t466 * qJD(3);
t517 = t393 * t441;
t516 = t402 * t441;
t515 = t402 * t444;
t514 = t428 * t430;
t513 = t430 * t431;
t380 = qJD(2) * t448 + t451 * t487;
t366 = t442 * t380;
t511 = qJDD(2) - g(3);
t361 = qJD(3) * t414 + qJDD(1) * t402 - t442 * t494;
t348 = -t444 * qJDD(3) + t361 * t441;
t349 = qJDD(3) * t441 + t361 * t444;
t502 = qJD(6) * t450;
t495 = -t447 * t348 + t450 * t349 + t378 * t502;
t503 = qJD(6) * t447;
t281 = -t377 * t503 + t495;
t510 = -t281 * t466 - t390 * t471;
t291 = t442 * t320 + t526 * t326;
t288 = qJDD(3) * qJ(5) + qJD(3) * qJD(5) + t291;
t300 = pkin(4) * t360 - qJ(5) * t361 - qJD(5) * t391 + t457;
t277 = t444 * t288 + t441 * t300;
t307 = t393 * t403 + t502 * t515 - t503 * t516;
t355 = t403 * t402;
t509 = -t307 * t385 - t355 * t359;
t379 = t451 * qJD(2) - t448 * t487;
t369 = t379 + t527;
t492 = t526 * t380;
t325 = t442 * t369 + t492;
t317 = qJD(3) * qJ(5) + t325;
t387 = qJD(1) * t537 + qJD(4);
t339 = pkin(4) * t388 - qJ(5) * t391 + t387;
t294 = t444 * t317 + t441 * t339;
t332 = t379 * t526 - t366;
t352 = pkin(3) * t504 + pkin(4) * t391 + qJ(5) * t388;
t302 = t444 * t332 + t441 * t352;
t496 = t448 * t527;
t330 = pkin(4) * t390 - qJ(5) * t393 - qJD(5) * t402 + t496;
t488 = qJD(3) * t512;
t381 = qJD(4) * t451 - t448 * t488;
t382 = -qJD(4) * t448 - t451 * t488;
t338 = t381 * t526 + t442 * t382;
t299 = t441 * t330 + t444 * t338;
t354 = -pkin(4) * t466 - qJ(5) * t402 + t537;
t398 = t512 * t448;
t399 = t512 * t451;
t358 = -t442 * t398 + t399 * t526;
t309 = t441 * t354 + t444 * t358;
t423 = t434 + pkin(2);
t452 = cos(qJ(1));
t506 = t452 * pkin(1) + t431 * t423;
t439 = t448 ^ 2;
t505 = -t451 ^ 2 + t439;
t411 = qJD(1) * t422;
t324 = t369 * t526 - t366;
t316 = -qJD(3) * pkin(4) + qJD(5) - t324;
t501 = -qJD(5) + t316;
t498 = qJDD(1) * t451;
t276 = -t288 * t441 + t444 * t300;
t272 = pkin(5) * t360 - pkin(8) * t349 + t276;
t275 = -pkin(8) * t348 + t277;
t490 = t450 * t272 - t275 * t447;
t293 = -t317 * t441 + t444 * t339;
t298 = t444 * t330 - t338 * t441;
t301 = -t332 * t441 + t444 * t352;
t489 = t450 * t348 + t447 * t349;
t308 = t444 * t354 - t358 * t441;
t331 = t379 * t442 + t492;
t337 = t381 * t442 - t526 * t382;
t357 = t526 * t398 + t399 * t442;
t421 = -pkin(3) * t526 - pkin(4);
t449 = sin(qJ(1));
t483 = g(1) * t449 - g(2) * t452;
t482 = -t401 * t359 - t385 * t507;
t446 = -qJ(4) - pkin(7);
t481 = -pkin(1) * t449 - t431 * t446;
t480 = -pkin(4) * t430 - t525;
t479 = t272 * t447 + t275 * t450;
t478 = -t276 * t444 - t277 * t441;
t280 = pkin(5) * t388 - pkin(8) * t377 + t293;
t284 = pkin(8) * t378 + t294;
t273 = t280 * t450 - t284 * t447;
t274 = t280 * t447 + t284 * t450;
t282 = -qJD(6) * t471 + t489;
t477 = t282 * t466 + t390 * t542;
t476 = -t293 * t441 + t294 * t444;
t297 = -pkin(5) * t466 - pkin(8) * t515 + t308;
t304 = -pkin(8) * t516 + t309;
t475 = t297 * t450 - t304 * t447;
t474 = t297 * t447 + t304 * t450;
t306 = -t393 * t401 - t395 * t402;
t356 = t401 * t402;
t473 = -t306 * t385 + t356 * t359;
t472 = -t360 * t402 - t388 * t393;
t397 = t528 * t444;
t468 = pkin(5) * t391 + qJD(5) * t441 + qJD(6) * t397 + t388 * t535 + t301;
t396 = t528 * t441;
t467 = pkin(8) * t518 - qJD(5) * t444 + qJD(6) * t396 + t302;
t463 = -qJD(1) * t411 - t408 + t485;
t462 = 0.2e1 * qJD(3) * t411 - qJDD(3) * t420;
t458 = t289 * t402 + t316 * t393 - t485;
t453 = qJD(3) ^ 2;
t456 = -0.2e1 * qJDD(1) * t422 - t420 * t453 - t493;
t436 = pkin(11) + qJ(6);
t429 = cos(t436);
t426 = sin(t436);
t407 = qJDD(3) * t451 - t448 * t453;
t406 = qJDD(3) * t448 + t451 * t453;
t405 = -t444 * pkin(5) + t421;
t386 = t388 ^ 2;
t373 = t426 * t428 + t429 * t513;
t372 = -t426 * t513 + t428 * t429;
t371 = t426 * t431 - t429 * t514;
t370 = t426 * t514 + t429 * t431;
t336 = pkin(5) * t516 + t357;
t313 = pkin(5) * t517 + t337;
t310 = -pkin(5) * t518 + t331;
t305 = -pkin(5) * t378 + t316;
t287 = -pkin(8) * t517 + t299;
t283 = pkin(5) * t390 - t393 * t535 + t298;
t278 = t348 * pkin(5) + t289;
t1 = [qJDD(1) * MDP(1) + t483 * MDP(2) + (g(1) * t452 + g(2) * t449) * MDP(3) + (t483 + (t443 ^ 2 + t445 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (qJDD(1) * t439 + 0.2e1 * t451 * t494) * MDP(5) + 0.2e1 * (t448 * t498 - t500 * t505) * MDP(6) + t406 * MDP(7) + t407 * MDP(8) + (t448 * t462 + t451 * t456) * MDP(10) + (-t448 * t456 + t451 * t462) * MDP(11) + (-t290 * t402 + t291 * t466 - t324 * t393 - t325 * t390 + t337 * t391 - t338 * t388 + t357 * t361 - t358 * t360 - t485) * MDP(12) + (t291 * t358 + t325 * t338 - t290 * t357 - t324 * t337 + t457 * t537 + t387 * t496 - g(1) * (-t423 * t428 + t481) - g(2) * (-t428 * t446 + t506)) * MDP(13) + (-t276 * t466 + t293 * t390 + t298 * t388 + t308 * t360 - t337 * t378 + t357 * t348 + t441 * t458 - t444 * t539) * MDP(14) + (t277 * t466 - t294 * t390 - t299 * t388 - t309 * t360 + t337 * t377 + t357 * t349 + t441 * t539 + t444 * t458) * MDP(15) + (-t298 * t377 + t299 * t378 - t308 * t349 - t309 * t348 - t493 * t427 + t478 * t402 + (-t293 * t444 - t294 * t441) * t393) * MDP(16) + (t277 * t309 + t294 * t299 + t276 * t308 + t293 * t298 + t289 * t357 + t316 * t337 - g(1) * t481 - g(2) * (pkin(4) * t513 + t431 * t525 + t506) + (-g(1) * (-t423 + t480) + g(2) * t446) * t428) * MDP(17) + (-t281 * t356 - t306 * t471) * MDP(18) + (-t281 * t355 + t282 * t356 + t306 * t542 + t307 * t471) * MDP(19) + (-t473 + t510) * MDP(20) + (t477 + t509) * MDP(21) + (-t359 * t466 + t385 * t390) * MDP(22) + ((t283 * t450 - t287 * t447) * t385 + t475 * t359 - t490 * t466 + t273 * t390 - t313 * t542 + t336 * t282 + t278 * t355 + t305 * t307 - g(1) * t371 - g(2) * t373 + (t274 * t466 - t385 * t474) * qJD(6)) * MDP(23) + (-(t283 * t447 + t287 * t450) * t385 - t474 * t359 + t479 * t466 - t274 * t390 - t313 * t471 + t336 * t281 - t278 * t356 + t305 * t306 - g(1) * t370 - g(2) * t372 + (t273 * t466 - t385 * t475) * qJD(6)) * MDP(24); t511 * MDP(4) + t407 * MDP(10) - t406 * MDP(11) + (-t361 * t466 + t390 * t391 + t472) * MDP(12) + (t290 * t466 + t291 * t402 - t324 * t390 + t325 * t393 - g(3)) * MDP(13) + (-t348 * t466 - t378 * t390) * MDP(14) + (-t349 * t466 + t377 * t390) * MDP(15) + (-t289 * t466 + t316 * t390 - g(3)) * MDP(17) + (-t477 + t509) * MDP(23) + (t473 + t510) * MDP(24) + (t472 * MDP(15) + (-t348 * t402 + t378 * t393) * MDP(16) + (t277 * t402 + t294 * t393) * MDP(17)) * t444 + (t472 * MDP(14) + (t349 * t402 + t377 * t393) * MDP(16) + (-t276 * t402 - t293 * t393) * MDP(17)) * t441; MDP(7) * t499 + MDP(8) * t498 + qJDD(3) * MDP(9) + (t448 * t463 + t432 - t529) * MDP(10) + (-t448 * t511 + t451 * t463) * MDP(11) + ((t325 - t331) * t391 + (-t324 + t332) * t388 + (-t360 * t442 - t361 * t526) * pkin(3)) * MDP(12) + (t324 * t331 - t325 * t332 + (t526 * t290 - t529 + t291 * t442 + (-qJD(1) * t387 + t485) * t448) * pkin(3)) * MDP(13) + (-t416 * t521 - t293 * t391 + t331 * t378 + t348 * t421 + (t441 * t501 - t301) * t388 - t538 * t444) * MDP(14) + (-t416 * t520 + t294 * t391 - t331 * t377 + t349 * t421 + (t444 * t501 + t302) * t388 + t538 * t441) * MDP(15) + (-t531 + t301 * t377 - t302 * t378 - t485 * t430 + (qJD(5) * t378 - t293 * t388 - t348 * t416 + t277) * t444 + (qJD(5) * t377 - t294 * t388 + t349 * t416 - t276) * t441) * MDP(16) + (t289 * t421 - t294 * t302 - t293 * t301 - t316 * t331 - g(3) * (t434 - t480) + (-t276 * t441 + t277 * t444) * t416 + t476 * qJD(5) + t485 * (pkin(3) * t448 + pkin(4) * t427 - qJ(5) * t430)) * MDP(17) + (t281 * t403 + t471 * t508) * MDP(18) + (-t281 * t401 - t282 * t403 + t471 * t507 - t508 * t542) * MDP(19) + (t523 - t536) * MDP(20) + (t482 - t524) * MDP(21) - t385 * t391 * MDP(22) + ((-t396 * t450 - t397 * t447) * t359 + t405 * t282 + t278 * t401 - t273 * t391 + t310 * t542 + (t447 * t467 - t450 * t468) * t385 + t507 * t305 + t461 * t429) * MDP(23) + (-(-t396 * t447 + t397 * t450) * t359 + t405 * t281 + t278 * t403 + t274 * t391 + t310 * t471 + (t447 * t468 + t450 * t467) * t385 - t508 * t305 - t461 * t426) * MDP(24) + (-t448 * t451 * MDP(5) + t505 * MDP(6)) * qJD(1) ^ 2; (-t391 ^ 2 - t386) * MDP(12) + (t378 * t391 + t520) * MDP(14) + (-t377 * t391 - t386 * t444 - t521) * MDP(15) + (-t348 * t441 - t349 * t444) * MDP(16) + (-t316 * t391 - t478 + t493) * MDP(17) + (t482 + t524) * MDP(23) + (t523 + t536) * MDP(24) + ((t377 * t441 + t378 * t444) * MDP(16) + t476 * MDP(17) - MDP(14) * t518) * t388 + (t324 * t391 + t325 * t388 + t457 + t493) * MDP(13); (t377 * t388 + t348) * MDP(14) + (t378 * t388 + t349) * MDP(15) + (-t377 ^ 2 - t378 ^ 2) * MDP(16) + (t293 * t377 - t294 * t378 + t538) * MDP(17) + (t282 - t540) * MDP(23) + (t281 + t543) * MDP(24); t471 * t542 * MDP(18) + (t471 ^ 2 - t542 ^ 2) * MDP(19) + (t495 - t543) * MDP(20) + (-t489 - t540) * MDP(21) + t359 * MDP(22) + (-g(1) * t372 + g(2) * t370 + t274 * t385 + t305 * t471 + t426 * t531 + t490) * MDP(23) + (g(1) * t373 - g(2) * t371 + t273 * t385 - t305 * t542 + t429 * t531 - t479) * MDP(24) + (-MDP(20) * t519 + MDP(21) * t471 - MDP(23) * t274 - MDP(24) * t273) * qJD(6);];
tau  = t1;
