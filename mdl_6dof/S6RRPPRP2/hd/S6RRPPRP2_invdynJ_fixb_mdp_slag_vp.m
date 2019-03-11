% Calculate vector of inverse dynamics joint torques for
% S6RRPPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPPRP2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRP2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RRPPRP2_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:32:08
% EndTime: 2019-03-09 08:32:17
% DurationCPUTime: 6.90s
% Computational Cost: add. (4279->490), mult. (9633->592), div. (0->0), fcn. (6778->10), ass. (0->214)
t433 = sin(pkin(9));
t438 = sin(qJ(2));
t509 = qJD(1) * t438;
t434 = cos(pkin(9));
t441 = cos(qJ(2));
t523 = t434 * t441;
t381 = -qJD(1) * t523 + t433 * t509;
t437 = sin(qJ(5));
t440 = cos(qJ(5));
t394 = t433 * t441 + t434 * t438;
t383 = t394 * qJD(2);
t496 = qJDD(1) * t441;
t407 = t434 * t496;
t497 = qJDD(1) * t438;
t472 = -t433 * t497 + t407;
t446 = qJD(1) * t383 - t472;
t504 = qJD(5) * t440;
t505 = qJD(5) * t437;
t313 = qJD(2) * t505 - t440 * qJDD(2) - t381 * t504 - t437 * t446;
t362 = qJD(2) * t437 - t440 * t381;
t384 = t394 * qJD(1);
t558 = qJD(5) + t384;
t554 = t362 * t558;
t562 = -t313 + t554;
t364 = qJD(2) * t440 + t381 * t437;
t561 = t364 * t558;
t498 = qJD(1) * qJD(2);
t489 = t438 * t498;
t404 = t433 * t489;
t488 = t441 * t498;
t355 = qJDD(1) * t394 + t434 * t488 - t404;
t351 = qJDD(5) + t355;
t425 = t441 * pkin(2);
t419 = t425 + pkin(1);
t398 = -qJD(1) * t419 + qJD(3);
t453 = -qJ(4) * t384 + t398;
t544 = pkin(3) + pkin(8);
t316 = t381 * t544 + t453;
t436 = -qJ(3) - pkin(7);
t400 = t436 * t441;
t397 = qJD(1) * t400;
t387 = t433 * t397;
t399 = t436 * t438;
t396 = qJD(1) * t399;
t392 = qJD(2) * pkin(2) + t396;
t353 = t392 * t434 + t387;
t473 = qJD(4) - t353;
t541 = pkin(4) * t384;
t321 = -qJD(2) * t544 + t473 + t541;
t298 = t316 * t440 + t321 * t437;
t495 = pkin(2) * t489 + qJDD(3);
t529 = t355 * qJ(4);
t456 = -t384 * qJD(4) + t495 - t529;
t459 = t394 * pkin(3);
t525 = t433 * t438;
t294 = -t544 * t407 + (t525 * t544 - t419) * qJDD(1) + (pkin(8) * t394 + t459) * t498 + t456;
t486 = qJD(2) * t436;
t379 = -qJD(3) * t438 + t441 * t486;
t350 = qJDD(2) * pkin(2) + qJD(1) * t379 + qJDD(1) * t399;
t378 = qJD(3) * t441 + t438 * t486;
t356 = qJD(1) * t378 - qJDD(1) * t400;
t309 = t350 * t434 - t433 * t356;
t471 = qJDD(4) - t309;
t301 = pkin(4) * t355 - qJDD(2) * t544 + t471;
t483 = -t294 * t437 + t440 * t301;
t448 = -t298 * qJD(5) + t483;
t284 = pkin(5) * t351 + qJ(6) * t313 - qJD(6) * t364 + t448;
t291 = -qJ(6) * t362 + t298;
t560 = t291 * t558 + t284;
t507 = qJD(5) * t364;
t314 = qJDD(2) * t437 - t440 * t446 + t507;
t493 = -t440 * t294 - t437 * t301 - t321 * t504;
t458 = -t316 * t505 - t493;
t285 = -qJ(6) * t314 - qJD(6) * t362 + t458;
t297 = -t316 * t437 + t440 * t321;
t290 = -qJ(6) * t364 + t297;
t289 = pkin(5) * t558 + t290;
t559 = -t289 * t558 + t285;
t482 = t558 ^ 2;
t557 = -MDP(23) * t482 - MDP(24) * t562 + MDP(25) * t560;
t427 = qJ(2) + pkin(9);
t423 = cos(t427);
t439 = sin(qJ(1));
t442 = cos(qJ(1));
t550 = g(1) * t439 - g(2) * t442;
t553 = t550 * t423;
t422 = sin(t427);
t549 = -g(1) * t442 - g(2) * t439;
t552 = t549 * t422;
t393 = -t523 + t525;
t477 = -qJ(4) * t394 - t419;
t330 = t393 * t544 + t477;
t359 = -t434 * t399 - t400 * t433;
t340 = pkin(4) * t394 + t359;
t514 = t440 * t330 + t437 * t340;
t412 = t422 * qJ(4);
t551 = -t423 * pkin(3) - t412;
t358 = t396 * t434 + t387;
t501 = -qJD(4) + t358;
t518 = t440 * t442;
t521 = t437 * t439;
t374 = t422 * t518 - t521;
t519 = t439 * t440;
t520 = t437 * t442;
t376 = t422 * t519 + t520;
t415 = g(3) * t423;
t548 = -g(1) * t374 - g(2) * t376 + t440 * t415;
t524 = t434 * t397;
t354 = t433 * t392 - t524;
t345 = -qJD(2) * qJ(4) - t354;
t542 = pkin(4) * t381;
t326 = -t345 - t542;
t417 = -pkin(2) * t434 - pkin(3);
t410 = -pkin(8) + t417;
t547 = t326 * t558 + t410 * t351;
t545 = t364 ^ 2;
t543 = pkin(2) * t438;
t540 = pkin(5) * t437;
t539 = pkin(5) * t440;
t533 = g(3) * t441;
t532 = qJ(6) * t440;
t531 = qJDD(2) * pkin(3);
t530 = t313 * t440;
t528 = t393 * t437;
t526 = t423 * t442;
t522 = t436 * t442;
t344 = t440 * t351;
t517 = qJ(6) - t410;
t516 = -t290 + t289;
t485 = pkin(2) * t509 + qJ(4) * t381;
t322 = t384 * t544 + t485;
t357 = t396 * t433 - t524;
t331 = t357 - t542;
t515 = t440 * t322 + t437 * t331;
t310 = t433 * t350 + t434 * t356;
t329 = t440 * t331;
t503 = qJD(6) * t440;
t513 = t505 * t517 - t503 + pkin(5) * t381 - t329 - (-qJ(6) * t384 - t322) * t437;
t391 = t517 * t440;
t512 = -qJD(5) * t391 - qJD(6) * t437 - t384 * t532 - t515;
t418 = pkin(4) + t539;
t511 = t418 - t436;
t431 = t438 ^ 2;
t510 = -t441 ^ 2 + t431;
t508 = qJD(2) * t438;
t506 = qJD(5) * t410;
t502 = t364 * MDP(23);
t500 = t541 - t501;
t421 = pkin(2) * t508;
t492 = qJDD(2) * qJ(4) + t310;
t491 = t425 - t551;
t413 = pkin(2) * t433 + qJ(4);
t484 = -qJ(6) * t393 - t330;
t338 = t378 * t433 - t434 * t379;
t481 = t437 * t558;
t480 = t440 * t558;
t429 = qJD(2) * qJD(4);
t306 = -t429 - t492;
t406 = t442 * t419;
t478 = g(2) * (pkin(3) * t526 + t442 * t412 + t406);
t476 = -pkin(3) * t422 - t543;
t339 = t378 * t434 + t379 * t433;
t360 = t399 * t433 - t400 * t434;
t470 = -t550 + t495;
t435 = -qJ(6) - pkin(8);
t469 = t422 * t540 - t423 * t435;
t468 = -t434 * t489 + t407;
t467 = -t419 + t551;
t466 = -0.2e1 * pkin(1) * t498 - pkin(7) * qJDD(2);
t465 = t383 * t437 + t393 * t504;
t464 = -t383 * t440 + t393 * t505;
t386 = qJD(2) * t523 - t433 * t508;
t461 = -qJ(4) * t386 - qJD(4) * t394 + t421;
t460 = pkin(3) * t525 - t419;
t311 = t383 * t544 + t461;
t323 = pkin(4) * t386 + t338;
t457 = t440 * t311 + t437 * t323 - t330 * t505 + t340 * t504;
t324 = -pkin(4) * t383 + t339;
t334 = pkin(3) * t381 + t453;
t455 = t334 * t384 + t415 + t471;
t454 = -g(3) * t422 + t423 * t549;
t443 = qJD(2) ^ 2;
t451 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t443 + t550;
t444 = qJD(1) ^ 2;
t450 = pkin(1) * t444 - pkin(7) * qJDD(1) - t549;
t302 = -pkin(4) * t446 - t306;
t449 = t302 + t454;
t447 = t338 * t384 - t339 * t381 + t359 * t355 + t549;
t445 = -MDP(22) * t482 - t351 * MDP(23) + MDP(25) * t559;
t288 = t314 * pkin(5) + qJDD(6) + t302;
t403 = qJ(4) * t526;
t401 = t439 * t423 * qJ(4);
t390 = t517 * t437;
t380 = t384 ^ 2;
t377 = -t422 * t521 + t518;
t375 = t422 * t520 + t519;
t372 = qJD(2) * t381;
t361 = t362 ^ 2;
t352 = pkin(3) * t393 + t477;
t342 = -qJD(2) * pkin(3) + t473;
t341 = -pkin(4) * t393 + t360;
t337 = pkin(3) * t384 + t485;
t336 = t440 * t340;
t325 = pkin(3) * t383 + t461;
t320 = t440 * t323;
t312 = t440 * t314;
t307 = t471 - t531;
t305 = pkin(5) * t362 + qJD(6) + t326;
t304 = -(-t488 * t433 + t468) * pkin(3) + t460 * qJDD(1) + t456;
t303 = t393 * t532 + t514;
t296 = pkin(5) * t394 + t437 * t484 + t336;
t287 = -qJ(6) * t464 + t393 * t503 + t457;
t286 = pkin(5) * t386 + t320 + t484 * t504 + (-qJ(6) * t383 - qJD(5) * t340 - qJD(6) * t393 - t311) * t437;
t1 = [(qJDD(1) * t431 + 0.2e1 * t438 * t488) * MDP(4) + 0.2e1 * (t438 * t496 - t498 * t510) * MDP(5) + (qJDD(2) * t438 + t441 * t443) * MDP(6) + (qJDD(2) * t441 - t438 * t443) * MDP(7) + (t438 * t466 + t441 * t451) * MDP(9) + (-t438 * t451 + t441 * t466) * MDP(10) + (t360 * (-qJD(2) * t384 + t472) - t310 * t393 - t354 * t383 - t309 * t394 - t353 * t386 + t447) * MDP(11) + (t310 * t360 + t354 * t339 - t309 * t359 - t353 * t338 - (-qJDD(1) * t419 + t495) * t419 + t398 * t421 - g(1) * (-t419 * t439 - t522) - g(2) * (-t436 * t439 + t406)) * MDP(12) + (t306 * t393 + t307 * t394 + t342 * t386 + t345 * t383 - t360 * t446 + t447) * MDP(13) + (-t325 * t381 + t352 * t472 - t304 * t393 - t334 * t383 + t359 * qJDD(2) - t553 + (-t352 * t384 + t338) * qJD(2)) * MDP(14) + (qJD(2) * t339 + qJDD(2) * t360 - t304 * t394 - t325 * t384 - t334 * t386 - t352 * t355 + t422 * t550) * MDP(15) + (t304 * t352 + t334 * t325 - t306 * t360 - t345 * t339 + t307 * t359 + t342 * t338 + g(1) * t522 - t478 + (-g(1) * t467 + g(2) * t436) * t439) * MDP(16) + (-t313 * t528 + t364 * t465) * MDP(17) + ((-t362 * t437 + t364 * t440) * t383 + (-t530 - t314 * t437 + (-t362 * t440 - t364 * t437) * qJD(5)) * t393) * MDP(18) + (-t313 * t394 + t351 * t528 + t364 * t386 + t465 * t558) * MDP(19) + (-t314 * t394 + t344 * t393 - t362 * t386 - t464 * t558) * MDP(20) + (t351 * t394 + t386 * t558) * MDP(21) + ((-t311 * t437 + t320) * t558 + (-t330 * t437 + t336) * t351 + t483 * t394 + t297 * t386 + t324 * t362 + t341 * t314 - g(1) * t377 - g(2) * t375 + (-t302 * t393 - t326 * t383) * t440 + (-t298 * t394 + t326 * t528 - t514 * t558) * qJD(5)) * MDP(22) + (g(1) * t376 - g(2) * t374 - t298 * t386 + t302 * t528 - t341 * t313 + t324 * t364 + t326 * t465 - t351 * t514 - t394 * t458 - t457 * t558) * MDP(23) + (-t286 * t364 - t287 * t362 + t296 * t313 - t303 * t314 + t553 + (-t289 * t437 + t291 * t440) * t383 + (-t284 * t437 + t285 * t440 + (-t289 * t440 - t291 * t437) * qJD(5)) * t393) * MDP(24) + (t285 * t303 + t291 * t287 + t284 * t296 + t289 * t286 + t288 * (-t393 * t418 + t360) + t305 * (pkin(5) * t464 + t324) - t478 + (-g(1) * t511 - g(2) * t469) * t442 + (-g(1) * (t467 - t469) - g(2) * t511) * t439) * MDP(25) + qJDD(1) * MDP(1) + t550 * MDP(2) - t549 * MDP(3); MDP(6) * t497 + MDP(7) * t496 + qJDD(2) * MDP(8) + (t438 * t450 - t533) * MDP(9) + (g(3) * t438 + t441 * t450) * MDP(10) + ((t354 - t357) * t384 + (t358 - t353) * t381 + (-t434 * t355 + ((-t488 - t497) * t433 + t468) * t433) * pkin(2)) * MDP(11) + (t353 * t357 - t354 * t358 + (-t533 + t309 * t434 + t310 * t433 + (-qJD(1) * t398 - t549) * t438) * pkin(2)) * MDP(12) + (-t413 * t446 + t417 * t355 + (-t345 - t357) * t384 + (t342 + t501) * t381) * MDP(13) + (-qJD(2) * t357 + t337 * t381 + t552 + (-pkin(3) + t417) * qJDD(2) + t455) * MDP(14) + (-qJD(2) * t358 + qJDD(2) * t413 - t334 * t381 + t337 * t384 + 0.2e1 * t429 + t454 + t492) * MDP(15) + (-t306 * t413 + t307 * t417 - t334 * t337 - t342 * t357 - g(1) * (t442 * t476 + t403) - g(2) * (t439 * t476 + t401) - g(3) * t491 + t501 * t345) * MDP(16) + (-t364 * t481 - t530) * MDP(17) + (-t312 - t364 * t480 + (t313 + t554) * t437) * MDP(18) + (t364 * t381 - t481 * t558 + t344) * MDP(19) + (-t351 * t437 - t362 * t381 - t480 * t558) * MDP(20) + t558 * t381 * MDP(21) + (t297 * t381 + t413 * t314 - t329 * t558 + t500 * t362 + t547 * t440 + ((t322 - t506) * t558 + t449) * t437) * MDP(22) + (-t413 * t313 + t515 * t558 - t298 * t381 + t500 * t364 - t547 * t437 + (-t506 * t558 + t449) * t440) * MDP(23) + (-t313 * t391 + t314 * t390 - t512 * t362 - t513 * t364 - t559 * t437 - t440 * t560 - t415 - t552) * MDP(24) + (-t285 * t390 - t284 * t391 + t288 * (t413 + t540) - g(1) * t403 - g(2) * t401 - g(3) * (t469 + t491) + (t539 * t558 + t500) * t305 + t512 * t291 + t513 * t289 + t549 * (t423 * t540 - t543 + (-pkin(3) + t435) * t422)) * MDP(25) + (-MDP(4) * t438 * t441 + MDP(5) * t510) * t444; -t380 * MDP(11) + t470 * MDP(12) + t407 * MDP(14) + (t372 + t404) * MDP(15) + (-pkin(3) * t407 + t470 - t529) * MDP(16) - t312 * MDP(24) - t550 * MDP(25) + (t353 * MDP(12) + (-qJD(4) - t342) * MDP(16) - MDP(13) * t384) * t384 + (-t384 * MDP(14) + (-MDP(14) * t394 - MDP(15) * t523 + MDP(16) * t459) * qJD(1)) * qJD(2) + (t354 * MDP(12) - t345 * MDP(16) + t362 * MDP(22) + t502 + t305 * MDP(25) + (-MDP(11) - MDP(13)) * t381) * t381 + (-t419 * MDP(12) - MDP(14) * t525 - t394 * MDP(15) + MDP(16) * t460) * qJDD(1) + (MDP(24) * t561 + t445) * t440 + (-t351 * MDP(22) - t557) * t437; (t372 + t355) * MDP(13) + (-t381 * t384 + qJDD(2)) * MDP(14) + (-t380 - t443) * MDP(15) + (qJD(2) * t345 + t455 - t531) * MDP(16) + (-qJD(2) * t362 + t344) * MDP(22) - qJD(2) * t502 + (-qJD(2) * t305 + t415) * MDP(25) + t557 * t440 + ((t364 * t384 - t314 + t507) * MDP(24) + t445) * t437 - (-MDP(16) - MDP(25)) * t552; t364 * t362 * MDP(17) + (-t361 + t545) * MDP(18) + t562 * MDP(19) + (-t314 + t561) * MDP(20) + t351 * MDP(21) + (t298 * t558 - t326 * t364 + t448 + t548) * MDP(22) + (g(1) * t375 - g(2) * t377 + t297 * t558 + t326 * t362 + (qJD(5) * t316 - t415) * t437 + t493) * MDP(23) + (pkin(5) * t313 - t362 * t516) * MDP(24) + (t516 * t291 + (-t305 * t364 + t284 + t548) * pkin(5)) * MDP(25); (-t361 - t545) * MDP(24) + (t289 * t364 + t291 * t362 + t288 + t454) * MDP(25);];
tau  = t1;
