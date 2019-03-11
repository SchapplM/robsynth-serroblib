% Calculate vector of inverse dynamics joint torques for
% S6PRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRPR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6PRPRPR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:49:45
% EndTime: 2019-03-08 19:49:54
% DurationCPUTime: 6.99s
% Computational Cost: add. (2715->496), mult. (5741->680), div. (0->0), fcn. (4402->14), ass. (0->211)
t432 = sin(pkin(11));
t435 = cos(pkin(11));
t437 = sin(qJ(6));
t440 = cos(qJ(6));
t572 = -t432 * t437 + t440 * t435;
t570 = t572 * qJD(6);
t438 = sin(qJ(4));
t525 = qJD(2) * t438;
t421 = qJD(6) + t525;
t441 = cos(qJ(4));
t524 = qJD(2) * t441;
t503 = t432 * t524;
t519 = t435 * qJD(4);
t390 = t503 - t519;
t500 = t435 * t524;
t523 = qJD(4) * t432;
t392 = t500 + t523;
t467 = t390 * t437 - t392 * t440;
t575 = t421 * t467;
t478 = pkin(4) * t441 + qJ(5) * t438;
t368 = qJD(4) * t478 - qJD(5) * t441 + qJD(3);
t442 = cos(qJ(2));
t443 = -pkin(2) - pkin(8);
t521 = qJD(4) * t441;
t498 = t443 * t521;
t434 = sin(pkin(6));
t530 = qJD(1) * t434;
t439 = sin(qJ(2));
t543 = t438 * t439;
t536 = t432 * t368 + t435 * t498 - (t432 * t442 + t435 * t543) * t530;
t574 = t435 * t368 - (-t432 * t543 + t435 * t442) * t530;
t497 = t442 * t530;
t476 = qJD(3) - t497;
t386 = qJD(2) * t443 + t476;
t436 = cos(pkin(6));
t544 = t436 * t438;
t573 = -qJD(1) * t544 + t386 * t441;
t562 = cos(pkin(10));
t485 = t562 * t442;
t433 = sin(pkin(10));
t548 = t433 * t439;
t376 = -t436 * t485 + t548;
t486 = t562 * t439;
t547 = t433 * t442;
t378 = t436 * t547 + t486;
t545 = t434 * t442;
t453 = g(1) * t378 + g(2) * t376 - g(3) * t545;
t377 = t436 * t486 + t547;
t379 = -t436 * t548 + t485;
t571 = -g(1) * t379 - g(2) * t377;
t396 = t432 * t440 + t435 * t437;
t458 = t396 * qJD(6);
t526 = qJD(2) * t434;
t502 = t439 * t526;
t409 = qJD(1) * t502;
t517 = qJDD(1) * t434;
t491 = t442 * t517;
t465 = qJDD(3) + t409 - t491;
t348 = qJDD(2) * t443 + t465;
t569 = -qJDD(4) * pkin(4) - t348 * t441 + qJDD(5);
t540 = qJDD(1) - g(3);
t546 = t434 * t439;
t568 = -t540 * t546 - t571;
t506 = t439 * t530;
t528 = qJD(2) * qJ(3);
t399 = t506 + t528;
t567 = qJD(4) * (-t399 + t506 - t528) - qJDD(4) * t443;
t566 = pkin(9) * t441;
t563 = pkin(9) + qJ(5);
t561 = qJDD(2) * pkin(2);
t529 = qJD(1) * t441;
t417 = t436 * t529;
t516 = qJDD(1) * t436;
t522 = qJD(4) * t438;
t507 = -qJD(4) * t417 - t386 * t522 - t438 * t516;
t304 = -t507 + t569;
t559 = t304 * t441;
t372 = t440 * t390;
t555 = t392 * t437;
t324 = t372 + t555;
t558 = t324 * t421;
t518 = qJD(2) * qJD(4);
t493 = t441 * t518;
t513 = qJDD(2) * t438;
t459 = t493 + t513;
t394 = qJDD(6) + t459;
t554 = t394 * t572;
t553 = t394 * t396;
t429 = pkin(11) + qJ(6);
t426 = sin(t429);
t552 = t426 * t438;
t427 = cos(t429);
t551 = t427 * t438;
t549 = t433 * t434;
t542 = t438 * t443;
t490 = t441 * t516;
t303 = t490 + qJDD(4) * qJ(5) + t348 * t438 + (qJD(5) + t573) * qJD(4);
t477 = pkin(4) * t438 - qJ(5) * t441;
t400 = qJ(3) + t477;
t492 = t439 * t517;
t310 = t492 + t400 * qJDD(2) + (t368 + t497) * qJD(2);
t296 = t435 * t303 + t432 * t310;
t489 = -t432 * t443 + pkin(5);
t509 = pkin(9) * t435 * t438;
t539 = (t441 * t489 + t509) * qJD(4) + t574;
t499 = t432 * t522;
t538 = -pkin(9) * t499 - t536;
t342 = t438 * t386 + t417;
t332 = qJD(4) * qJ(5) + t342;
t354 = qJD(2) * t400 + t506;
t306 = t435 * t332 + t432 * t354;
t537 = -t432 * t498 + t574;
t398 = t478 * qJD(2);
t313 = t432 * t398 + t435 * t573;
t460 = t572 * t438;
t535 = qJD(2) * t460 + t570;
t457 = t396 * qJD(2);
t534 = t438 * t457 + t458;
t353 = t432 * t400 + t435 * t542;
t512 = qJDD(2) * t441;
t533 = t432 * qJDD(4) + t435 * t512;
t430 = t438 ^ 2;
t431 = t441 ^ 2;
t532 = t430 - t431;
t444 = qJD(4) ^ 2;
t445 = qJD(2) ^ 2;
t531 = -t444 - t445;
t527 = qJD(2) * t399;
t515 = qJDD(2) * qJ(3);
t514 = qJDD(2) * t430;
t511 = qJDD(4) * t438;
t481 = -qJDD(4) * t435 + t432 * t512;
t494 = t438 * t518;
t355 = t432 * t494 - t481;
t356 = t435 * t494 - t533;
t508 = -qJD(6) * t372 + t437 * t355 - t440 * t356;
t505 = t439 * t529;
t504 = t432 * t525;
t501 = t442 * t526;
t496 = g(3) * (pkin(2) * t545 + qJ(3) * t546);
t495 = qJ(5) * t513;
t488 = pkin(5) * t432 - t443;
t487 = t434 * t562;
t295 = -t303 * t432 + t435 * t310;
t291 = pkin(5) * t459 + pkin(9) * t356 + t295;
t294 = pkin(9) * t355 + t296;
t484 = t440 * t291 - t437 * t294;
t305 = -t332 * t432 + t435 * t354;
t312 = t435 * t398 - t432 * t573;
t483 = -t440 * t355 - t437 * t356;
t482 = -t348 + t527;
t475 = t437 * t291 + t440 * t294;
t474 = -t295 * t432 + t296 * t435;
t300 = pkin(5) * t525 - pkin(9) * t392 + t305;
t302 = -pkin(9) * t390 + t306;
t292 = t300 * t440 - t302 * t437;
t293 = t300 * t437 + t302 * t440;
t473 = t305 * t435 + t306 * t432;
t472 = -t305 * t432 + t306 * t435;
t388 = t435 * t400;
t322 = -t435 * t566 + t438 * t489 + t388;
t331 = -t432 * t566 + t353;
t471 = t322 * t440 - t331 * t437;
t470 = t322 * t437 + t331 * t440;
t383 = t436 * t441 - t438 * t545;
t337 = -t383 * t432 + t435 * t546;
t338 = t383 * t435 + t432 * t546;
t469 = t337 * t440 - t338 * t437;
t468 = t337 * t437 + t338 * t440;
t382 = t441 * t545 + t544;
t405 = t563 * t435;
t463 = qJD(5) * t432 + qJD(6) * t405 + (pkin(5) * t441 + t509) * qJD(2) + t312;
t404 = t563 * t432;
t462 = pkin(9) * t504 - qJD(5) * t435 + qJD(6) * t404 + t313;
t461 = t421 * t572;
t298 = -qJD(6) * t555 + t508;
t456 = g(1) * (-t378 * t441 + t438 * t549) - g(2) * (t376 * t441 + t438 * t487) + g(3) * t382;
t334 = t378 * t438 + t441 * t549;
t336 = -t376 * t438 + t441 * t487;
t455 = g(1) * t334 - g(2) * t336 + g(3) * t383;
t329 = -qJD(4) * pkin(4) + qJD(5) - t573;
t452 = g(3) * t546 - t571;
t451 = -t304 + t456;
t450 = -qJ(5) * t521 + (-qJD(5) + t329) * t438;
t449 = t453 + t491;
t448 = t456 + t507;
t447 = qJDD(3) - t449;
t299 = -qJD(6) * t467 + t483;
t349 = t492 + t515 + (qJD(3) + t497) * qJD(2);
t446 = qJD(2) * t476 - t443 * t444 + t349 - t452 + t515;
t428 = qJDD(4) * t441;
t423 = -pkin(5) * t435 - pkin(4);
t397 = -qJD(2) * pkin(2) + t476;
t389 = t488 * t441;
t373 = t488 * t522;
t371 = t378 * pkin(2);
t370 = t376 * pkin(2);
t367 = t572 * t441;
t366 = t396 * t441;
t357 = t465 - t561;
t352 = -t432 * t542 + t388;
t340 = qJD(4) * t383 - t441 * t502;
t339 = -qJD(4) * t382 + t438 * t502;
t321 = -pkin(5) * t504 + t342;
t319 = -t437 * t438 * t519 - t440 * t499 + t441 * t570;
t318 = -qJD(4) * t460 - t441 * t458;
t317 = t339 * t435 + t432 * t501;
t316 = -t339 * t432 + t435 * t501;
t315 = pkin(5) * t390 + t329;
t297 = -pkin(5) * t355 + t304;
t1 = [t540 * MDP(1) + (qJDD(1) * t436 ^ 2 - g(3)) * MDP(7) + (-qJD(4) * t340 - qJDD(4) * t382) * MDP(13) + (-qJD(4) * t339 - qJDD(4) * t383) * MDP(14) + (t337 * t513 + t340 * t390 - t355 * t382) * MDP(15) + (-t338 * t513 + t340 * t392 - t356 * t382) * MDP(16) + (-t316 * t392 - t317 * t390 + t337 * t356 + t338 * t355) * MDP(17) + (t295 * t337 + t296 * t338 + t304 * t382 + t305 * t316 + t306 * t317 + t329 * t340 - g(3)) * MDP(18) + ((-qJD(6) * t468 + t316 * t440 - t317 * t437) * t421 + t469 * t394 + t340 * t324 + t382 * t299) * MDP(24) + (-(qJD(6) * t469 + t316 * t437 + t317 * t440) * t421 - t468 * t394 - t340 * t467 + t382 * t298) * MDP(25) + ((t316 * t438 + t337 * t521) * MDP(15) + (-t317 * t438 - t338 * t521) * MDP(16)) * qJD(2) + ((-MDP(4) + MDP(6)) * (qJDD(2) * t439 + t442 * t445) + (-MDP(3) + MDP(5)) * (-qJDD(2) * t442 + t439 * t445) + ((-t357 + t527) * MDP(7) + (MDP(13) * t438 + MDP(14) * t441) * t445) * t442 + ((qJD(2) * t397 + t349) * MDP(7) + t459 * MDP(13) + (-t494 + t512) * MDP(14)) * t439) * t434; qJDD(2) * MDP(2) + t449 * MDP(3) + t568 * MDP(4) + (t447 - 0.2e1 * t561) * MDP(5) + (0.2e1 * qJD(2) * qJD(3) + 0.2e1 * t515 - t568) * MDP(6) + (t349 * qJ(3) + t399 * qJD(3) - t357 * pkin(2) - g(1) * (qJ(3) * t379 - t371) - g(2) * (qJ(3) * t377 - t370) - t496 + (-t397 * t439 - t399 * t442) * t530) * MDP(7) + (qJDD(2) * t431 - 0.2e1 * t438 * t493) * MDP(8) + 0.2e1 * (-t438 * t512 + t518 * t532) * MDP(9) + (-t438 * t444 + t428) * MDP(10) + (-t441 * t444 - t511) * MDP(11) + (t446 * t438 - t441 * t567) * MDP(13) + (t438 * t567 + t446 * t441) * MDP(14) + (t453 * t432 + (t390 * t506 + t304 * t432 + t443 * t355 + (qJD(2) * t352 + t305) * qJD(4)) * t441 + (t352 * qJDD(2) + t295 + (-t329 * t432 + t390 * t443) * qJD(4) + t537 * qJD(2) - t452 * t435) * t438) * MDP(15) + (t453 * t435 + (t392 * t506 + t304 * t435 + t443 * t356 + (-qJD(2) * t353 - t306) * qJD(4)) * t441 + (-t353 * qJDD(2) - t296 + (-t329 * t435 + t392 * t443) * qJD(4) - t536 * qJD(2) + t452 * t432) * t438) * MDP(16) + (t352 * t356 + t353 * t355 - t537 * t392 - t536 * t390 + t473 * t522 + (-t295 * t435 - t296 * t432 + t452) * t441) * MDP(17) + (t296 * t353 + t295 * t352 - g(1) * (-pkin(8) * t378 - t371) - g(2) * (-pkin(8) * t376 - t370) - t496 + (t329 * t522 - t559) * t443 + t536 * t306 + t537 * t305 + (-g(3) * pkin(8) * t442 + (-g(3) * t477 + t329 * t529) * t439) * t434 + t571 * t400) * MDP(18) + (t298 * t367 - t318 * t467) * MDP(19) + (-t298 * t366 - t299 * t367 - t318 * t324 + t319 * t467) * MDP(20) + (t298 * t438 + t318 * t421 + t367 * t394 - t467 * t521) * MDP(21) + (-t299 * t438 - t319 * t421 - t324 * t521 - t366 * t394) * MDP(22) + (t394 * t438 + t421 * t521) * MDP(23) + (t471 * t394 + t484 * t438 + t292 * t521 - t373 * t324 + t389 * t299 + t297 * t366 + t315 * t319 - g(1) * (-t378 * t426 + t379 * t551) - g(2) * (-t376 * t426 + t377 * t551) + (t538 * t437 + t539 * t440) * t421 + (-t293 * t438 - t421 * t470) * qJD(6) + (t324 * t505 - g(3) * (t426 * t442 + t427 * t543)) * t434) * MDP(24) + (-t470 * t394 - t475 * t438 - t293 * t521 + t373 * t467 + t389 * t298 + t297 * t367 + t315 * t318 - g(1) * (-t378 * t427 - t379 * t552) - g(2) * (-t376 * t427 - t377 * t552) + (-t437 * t539 + t440 * t538) * t421 + (-t292 * t438 - t421 * t471) * qJD(6) + (-t467 * t505 - g(3) * (-t426 * t543 + t427 * t442)) * t434) * MDP(25); qJDD(2) * MDP(5) - t445 * MDP(6) + (t409 + t447 - t527 - t561) * MDP(7) + (t438 * t531 + t428) * MDP(13) + (t441 * t531 - t511) * MDP(14) + (-t432 * t514 + t355 * t441 + (-t435 * t445 + (t390 - 0.2e1 * t503) * qJD(4)) * t438) * MDP(15) + (-t435 * t514 + t356 * t441 + (t432 * t445 + (t392 - 0.2e1 * t500) * qJD(4)) * t438) * MDP(16) + ((qJD(2) * t392 + t355 * t438 - t390 * t521) * t435 + (qJD(2) * t390 - t356 * t438 + t392 * t521) * t432) * MDP(17) + (-t559 + t474 * t438 - t473 * qJD(2) + (t329 * t438 + t441 * t472) * qJD(4) - t453) * MDP(18) + (-qJD(2) * t461 + (-qJD(4) * t396 * t421 - t299) * t441 + (qJD(4) * t324 - t421 * t570 - t553) * t438) * MDP(24) + (t421 * t457 + (-qJD(4) * t461 - t298) * t441 + (-qJD(4) * t467 + t421 * t458 - t554) * t438) * MDP(25); MDP(10) * t512 - MDP(11) * t513 + qJDD(4) * MDP(12) + (qJD(4) * t342 - t441 * t482 + t448) * MDP(13) + (t482 * t438 + t455 - t490) * MDP(14) + (-t432 * t495 + pkin(4) * t355 - t342 * t390 + t451 * t435 + (-t305 * t441 - t312 * t438 + t432 * t450) * qJD(2)) * MDP(15) + (-t435 * t495 + pkin(4) * t356 - t342 * t392 - t451 * t432 + (t306 * t441 + t313 * t438 + t435 * t450) * qJD(2)) * MDP(16) + (t312 * t392 + t313 * t390 + (qJ(5) * t355 - qJD(5) * t390 - t305 * t525 + t296) * t435 + (-qJ(5) * t356 + qJD(5) * t392 - t306 * t525 - t295) * t432 - t455) * MDP(17) + (-t305 * t312 - t306 * t313 - t329 * t342 + t472 * qJD(5) + t451 * pkin(4) + (-t455 + t474) * qJ(5)) * MDP(18) + (t298 * t396 - t467 * t535) * MDP(19) + (t298 * t572 - t299 * t396 - t324 * t535 + t467 * t534) * MDP(20) + (t421 * t535 + t467 * t524 + t553) * MDP(21) + (t324 * t524 - t421 * t534 + t554) * MDP(22) - t421 * MDP(23) * t524 + ((-t404 * t440 - t405 * t437) * t394 + t423 * t299 - t297 * t572 - t292 * t524 - t321 * t324 + (t437 * t462 - t440 * t463) * t421 + t534 * t315 + t456 * t427) * MDP(24) + (-(-t404 * t437 + t405 * t440) * t394 + t423 * t298 + t297 * t396 + t293 * t524 + t321 * t467 + (t437 * t463 + t440 * t462) * t421 + t535 * t315 - t456 * t426) * MDP(25) + (MDP(8) * t438 * t441 - MDP(9) * t532) * t445; ((t392 - t523) * t525 + t481) * MDP(15) + ((-t390 - t519) * t525 + t533) * MDP(16) + (-t390 ^ 2 - t392 ^ 2) * MDP(17) + (t305 * t392 + t306 * t390 - t448 + t569) * MDP(18) + (t299 - t575) * MDP(24) + (t298 - t558) * MDP(25); -t467 * t324 * MDP(19) + (-t324 ^ 2 + t467 ^ 2) * MDP(20) + (t508 + t558) * MDP(21) + (-t483 - t575) * MDP(22) + t394 * MDP(23) + (t293 * t421 + t315 * t467 - g(1) * (-t334 * t426 + t379 * t427) - g(2) * (t336 * t426 + t377 * t427) - g(3) * (-t383 * t426 + t427 * t546) + t484) * MDP(24) + (t292 * t421 + t315 * t324 - g(1) * (-t334 * t427 - t379 * t426) - g(2) * (t336 * t427 - t377 * t426) - g(3) * (-t383 * t427 - t426 * t546) - t475) * MDP(25) + (-MDP(21) * t555 + MDP(22) * t467 - MDP(24) * t293 - MDP(25) * t292) * qJD(6);];
tau  = t1;
