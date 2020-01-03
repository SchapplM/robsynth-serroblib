% Calculate vector of inverse dynamics joint torques for
% S5RRRPP7
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RRRPP7_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP7_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP7_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_invdynJ_fixb_mdp_slag_vp: pkin has to be [7x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S5RRRPP7_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:06:07
% EndTime: 2019-12-31 21:06:16
% DurationCPUTime: 6.34s
% Computational Cost: add. (2829->518), mult. (6105->602), div. (0->0), fcn. (3782->6), ass. (0->205)
t441 = cos(qJ(2));
t425 = t441 * qJDD(1);
t438 = sin(qJ(2));
t504 = qJD(1) * qJD(2);
t574 = -t438 * t504 + t425;
t370 = qJDD(3) - t574;
t519 = qJD(1) * t441;
t581 = qJD(3) - t519;
t582 = t370 * qJ(4) + qJD(4) * t581;
t442 = cos(qJ(1));
t538 = t438 * t442;
t439 = sin(qJ(1));
t560 = g(2) * t439;
t576 = g(1) * t538 + t438 * t560;
t490 = t441 * t504;
t502 = qJDD(1) * t438;
t580 = qJD(2) * qJD(3) + t490 + t502;
t437 = sin(qJ(3));
t440 = cos(qJ(3));
t513 = qJD(3) * t438;
t489 = qJD(1) * t513;
t465 = (-qJDD(2) + t489) * t440;
t320 = t437 * (qJD(2) * (qJD(3) + t519) + t502) + t465;
t542 = t437 * qJ(4);
t566 = pkin(3) + pkin(4);
t462 = -t440 * t566 - t542;
t367 = pkin(2) - t462;
t508 = t440 * qJD(2);
t520 = qJD(1) * t438;
t372 = t437 * t520 - t508;
t579 = t320 * qJ(5) + t372 * qJD(5);
t578 = 0.2e1 * t582;
t422 = pkin(6) * t519;
t575 = qJD(4) * t437 + t422;
t573 = g(1) * t442 + t560;
t427 = t438 * pkin(7);
t430 = t441 * pkin(2);
t499 = -pkin(1) - t430;
t469 = t499 - t427;
t361 = t469 * qJD(1);
t390 = qJD(2) * pkin(7) + t422;
t325 = t440 * t361 - t437 * t390;
t506 = qJD(4) - t325;
t362 = t370 * pkin(3);
t572 = t362 - qJDD(4);
t497 = t440 * t520;
t517 = qJD(2) * t437;
t374 = t497 + t517;
t389 = -qJD(2) * pkin(2) + pkin(6) * t520;
t464 = qJ(4) * t374 - t389;
t317 = pkin(3) * t372 - t464;
t563 = pkin(7) * t370;
t571 = -t317 * t581 + t563;
t305 = -t372 * t566 + qJD(5) + t464;
t536 = t439 * t441;
t354 = t437 * t536 + t440 * t442;
t533 = t442 * t437;
t356 = -t439 * t440 + t441 * t533;
t478 = pkin(2) * t438 - pkin(7) * t441;
t379 = t478 * qJD(2);
t331 = qJD(1) * t379 + qJDD(1) * t469;
t351 = t574 * pkin(6) + qJDD(2) * pkin(7);
t512 = qJD(3) * t440;
t514 = qJD(3) * t437;
t482 = -t440 * t331 + t437 * t351 + t361 * t514 + t390 * t512;
t541 = t437 * t438;
t452 = g(1) * t356 + g(2) * t354 + g(3) * t541 - t482;
t450 = t452 + t572;
t481 = t437 * qJDD(2) + t440 * t580;
t319 = t437 * t489 - t481;
t552 = qJ(5) * t319;
t570 = (qJD(5) + t305) * t374 + t450 - t552;
t568 = -0.2e1 * pkin(1);
t567 = t372 ^ 2;
t369 = t374 ^ 2;
t565 = pkin(3) * t437;
t564 = pkin(4) * t370;
t562 = g(1) * t439;
t559 = g(2) * t442;
t558 = g(3) * t441;
t557 = pkin(7) - qJ(5);
t556 = pkin(7) * qJD(3);
t555 = qJ(4) * t320;
t554 = qJ(4) * t372;
t553 = qJ(4) * t440;
t326 = t437 * t361 + t440 * t390;
t309 = qJ(5) * t372 + t326;
t393 = t581 * qJ(4);
t304 = t309 + t393;
t551 = t304 * t581;
t313 = t393 + t326;
t550 = t313 * t581;
t549 = t319 * t437;
t548 = t326 * t581;
t547 = t372 * t374;
t546 = t372 * t581;
t545 = t374 * t581;
t544 = t374 * t440;
t378 = t478 * qJD(1);
t543 = t378 * t440;
t540 = t437 * t441;
t539 = t438 * t440;
t445 = qJD(1) ^ 2;
t537 = t438 * t445;
t535 = t440 * t441;
t534 = t441 * t442;
t501 = t566 * t437;
t463 = -t501 + t553;
t532 = t463 * t581 + t575;
t509 = qJD(5) * t440;
t358 = t437 * t378;
t528 = qJ(4) * t520 + t358;
t531 = -t514 * t557 - t509 - (-t539 * pkin(6) + qJ(5) * t540) * qJD(1) - t528;
t472 = -t553 + t565;
t530 = t472 * t581 - t575;
t384 = t557 * t440;
t498 = -pkin(6) * t437 - pkin(3);
t484 = -pkin(4) + t498;
t529 = qJD(3) * t384 - qJD(5) * t437 + t543 - (-qJ(5) * t535 + t438 * t484) * qJD(1);
t522 = -t427 - t430;
t381 = -pkin(1) + t522;
t527 = t437 * t379 + t381 * t512;
t404 = pkin(6) * t535;
t526 = qJD(3) * t404 + t381 * t514;
t525 = t576 * t437;
t524 = t576 * t440;
t523 = t437 * t381 + t404;
t434 = t438 ^ 2;
t521 = -t441 ^ 2 + t434;
t518 = qJD(2) * t374;
t516 = qJD(2) * t438;
t515 = qJD(2) * t441;
t510 = qJD(4) * t440;
t308 = qJ(5) * t374 + t325;
t507 = qJD(4) - t308;
t420 = pkin(6) * t502;
t352 = -qJDD(2) * pkin(2) + pkin(6) * t490 + t420;
t500 = g(1) * t534 + g(2) * t536 + g(3) * t438;
t496 = t581 * t517;
t495 = t581 * t508;
t494 = t305 * t514;
t493 = t305 * t512;
t492 = t581 * t514;
t296 = t320 * pkin(3) + t319 * qJ(4) - t374 * qJD(4) + t352;
t294 = -pkin(4) * t320 + qJDD(5) - t296;
t488 = t294 - t558;
t355 = t439 * t535 - t533;
t487 = -t354 * pkin(3) + qJ(4) * t355;
t357 = t437 * t439 + t440 * t534;
t486 = -t356 * pkin(3) + qJ(4) * t357;
t402 = pkin(6) * t540;
t485 = t381 * t440 - t402;
t480 = pkin(3) * t535 + qJ(4) * t540 - t522;
t479 = t438 * t498;
t477 = g(1) * t354 - g(2) * t356;
t476 = g(1) * t355 - g(2) * t357;
t475 = -t355 * pkin(3) + t442 * pkin(6) - qJ(4) * t354;
t333 = -qJ(4) * t441 + t523;
t474 = t379 * t440 - t526;
t473 = pkin(3) * t440 + t542;
t471 = qJD(3) * t389 - t563;
t312 = -pkin(3) * t581 + t506;
t470 = t312 * t440 - t313 * t437;
t468 = pkin(2) + t473;
t467 = -t556 * t581 - t558;
t466 = qJ(4) * t516 - qJD(4) * t441 + t527;
t461 = -pkin(6) * qJDD(2) + t504 * t568;
t297 = t482 - t572;
t460 = t370 * t437 + t512 * t581;
t459 = t370 * t440 - t492;
t458 = -t296 + t467;
t457 = t437 * t331 + t440 * t351 + t361 * t512 - t390 * t514;
t444 = qJD(2) ^ 2;
t455 = pkin(6) * t444 + qJDD(1) * t568 - t562;
t454 = t442 * pkin(1) + pkin(2) * t534 + t357 * pkin(3) + t439 * pkin(6) + pkin(7) * t538 + qJ(4) * t356;
t295 = t457 + t582;
t451 = t319 - t546;
t448 = t317 * t374 - t450;
t447 = g(1) * t357 + g(2) * t355 + g(3) * t539 - t457;
t446 = t325 * t581 + t447;
t429 = t441 * pkin(3);
t412 = g(2) * t538;
t407 = pkin(7) * t534;
t403 = pkin(7) * t536;
t397 = qJ(4) * t539;
t383 = t557 * t437;
t341 = -t397 + (pkin(6) + t565) * t438;
t334 = t429 - t485;
t332 = t397 + (-pkin(6) - t501) * t438;
t330 = pkin(3) * t374 + t554;
t329 = qJD(1) * t479 - t543;
t328 = -pkin(6) * t497 + t528;
t324 = qJ(5) * t541 + t333;
t321 = pkin(4) * t441 + t402 + t429 + (-qJ(5) * t438 - t381) * t440;
t310 = -t374 * t566 - t554;
t307 = (qJD(3) * t473 - t510) * t438 + (pkin(6) + t472) * t515;
t306 = qJD(2) * t479 - t474;
t302 = (-t438 * t508 - t441 * t514) * pkin(6) + t466;
t301 = (qJD(3) * t462 + t510) * t438 + (-pkin(6) + t463) * t515;
t300 = -t566 * t581 + t507;
t299 = (-pkin(6) * qJD(2) + qJ(5) * qJD(3)) * t539 + (qJD(5) * t438 + (-pkin(6) * qJD(3) + qJ(5) * qJD(2)) * t441) * t437 + t466;
t298 = (-qJ(5) * t515 - t379) * t440 + (qJ(5) * t514 + qJD(2) * t484 - t509) * t438 + t526;
t293 = t295 + t579;
t292 = -qJD(5) * t374 + t297 + t552 - t564;
t1 = [qJDD(1) * MDP(1) + (qJDD(1) * t434 + 0.2e1 * t438 * t490) * MDP(4) + 0.2e1 * (t425 * t438 - t504 * t521) * MDP(5) + (qJDD(2) * t438 + t441 * t444) * MDP(6) + (qJDD(2) * t441 - t438 * t444) * MDP(7) + (t461 * t438 + (-t455 - t559) * t441) * MDP(9) + (t438 * t455 + t441 * t461 + t412) * MDP(10) + (-t319 * t539 + (-t437 * t513 + t441 * t508) * t374) * MDP(11) + ((-t372 * t440 - t374 * t437) * t515 + (t549 - t320 * t440 + (t372 * t437 - t544) * qJD(3)) * t438) * MDP(12) + ((t319 + t495) * t441 + (t459 + t518) * t438) * MDP(13) + ((t320 - t496) * t441 + (-qJD(2) * t372 - t460) * t438) * MDP(14) + (-t370 * t441 + t516 * t581) * MDP(15) + (t474 * t581 + t485 * t370 + ((pkin(6) * t372 + t389 * t437) * qJD(2) + t482) * t441 + (t389 * t512 + t325 * qJD(2) + t352 * t437 + (t320 + t496) * pkin(6)) * t438 + t476) * MDP(16) + (-t527 * t581 - t523 * t370 + (t389 * t508 + (t492 + t518) * pkin(6) + t457) * t441 + (-t389 * t514 - t326 * qJD(2) + t352 * t440 + (-t319 + t495) * pkin(6)) * t438 - t477) * MDP(17) + (-t306 * t581 + t307 * t372 + t320 * t341 - t334 * t370 + (t317 * t517 + t297) * t441 + (-qJD(2) * t312 + t296 * t437 + t317 * t512) * t438 + t476) * MDP(18) + (-t302 * t372 + t306 * t374 - t319 * t334 - t320 * t333 - t412 + t470 * t515 + (t562 - t295 * t437 + t297 * t440 + (-t312 * t437 - t313 * t440) * qJD(3)) * t438) * MDP(19) + (t302 * t581 - t307 * t374 + t319 * t341 + t333 * t370 + (-t317 * t508 - t295) * t441 + (qJD(2) * t313 - t296 * t440 + t317 * t514) * t438 + t477) * MDP(20) + (-g(1) * t475 - g(2) * t454 + t295 * t333 + t296 * t341 + t297 * t334 + t313 * t302 + t312 * t306 + t317 * t307 - t469 * t562) * MDP(21) + (-t298 * t581 - t301 * t372 - t320 * t332 - t321 * t370 + (-t305 * t517 + t292) * t441 + (-qJD(2) * t300 - t294 * t437 - t493) * t438 + t476) * MDP(22) + (t299 * t581 + t301 * t374 - t319 * t332 + t324 * t370 + (t305 * t508 - t293) * t441 + (qJD(2) * t304 + t294 * t440 - t494) * t438 + t477) * MDP(23) + (-t298 * t374 + t299 * t372 + t319 * t321 + t320 * t324 + t412 + (-t300 * t440 + t304 * t437) * t515 + (-t562 - t292 * t440 + t293 * t437 + (t300 * t437 + t304 * t440) * qJD(3)) * t438) * MDP(24) + (t293 * t324 + t304 * t299 + t292 * t321 + t300 * t298 + t294 * t332 + t305 * t301 - g(1) * (-pkin(4) * t355 + t475) - g(2) * (pkin(4) * t357 - qJ(5) * t538 + t454) - (-t438 * t557 + t499) * t562) * MDP(25) + (-t559 + t562) * MDP(2) + t573 * MDP(3); -t441 * MDP(4) * t537 + t521 * t445 * MDP(5) + MDP(6) * t502 + MDP(7) * t425 + qJDD(2) * MDP(8) + (pkin(1) * t537 - t420 - t558 + t576) * MDP(9) + ((pkin(1) * t445 - pkin(6) * qJDD(1)) * t441 + t500) * MDP(10) + (t544 * t581 - t549) * MDP(11) + ((-t319 - t546) * t440 + (-t320 - t545) * t437) * MDP(12) + ((-t374 * t438 - t535 * t581) * qJD(1) + t460) * MDP(13) + ((t372 * t438 + t540 * t581) * qJD(1) + t459) * MDP(14) - t581 * MDP(15) * t520 + (-pkin(2) * t320 + t471 * t437 + (-t558 - t352 - (t378 + t556) * t581) * t440 + (-t389 * t540 - t325 * t438 + (-t372 * t441 - t541 * t581) * pkin(6)) * qJD(1) + t524) * MDP(16) + (pkin(2) * t319 + t358 * t581 + t471 * t440 + (t352 - t467) * t437 + (-t389 * t535 + t326 * t438 + (-t374 * t441 - t539 * t581) * pkin(6)) * qJD(1) - t525) * MDP(17) + (t312 * t520 - t320 * t468 + t329 * t581 + t530 * t372 - t437 * t571 + t458 * t440 + t524) * MDP(18) + (t328 * t372 - t329 * t374 + (t295 + t581 * t312 + (qJD(3) * t374 - t320) * pkin(7)) * t440 + (t297 - t550 + (qJD(3) * t372 - t319) * pkin(7)) * t437 - t500) * MDP(19) + (-t313 * t520 - t319 * t468 - t328 * t581 - t530 * t374 + t458 * t437 + t440 * t571 + t525) * MDP(20) + (-t313 * t328 - t312 * t329 - g(1) * t407 - g(2) * t403 - g(3) * t480 + t530 * t317 + (qJD(3) * t470 + t295 * t440 + t297 * t437) * pkin(7) + (t573 * t438 - t296) * t468) * MDP(21) + (-t494 - t320 * t367 - t370 * t383 + t488 * t440 - t529 * t581 - t532 * t372 + (t300 * t438 + t305 * t540) * qJD(1) + t524) * MDP(22) + (t493 - t319 * t367 + t370 * t384 + t488 * t437 + t531 * t581 + t532 * t374 + (-t304 * t438 - t305 * t535) * qJD(1) + t525) * MDP(23) + (t319 * t383 + t320 * t384 - t529 * t374 + t531 * t372 + (-t300 * t581 - t293) * t440 + (-t292 + t551) * t437 + t500) * MDP(24) + (t293 * t384 + t292 * t383 + t294 * t367 - g(1) * (-qJ(5) * t534 + t407) - g(2) * (-qJ(5) * t536 + t403) - g(3) * (pkin(4) * t535 + t480) + t532 * t305 + t531 * t304 + t529 * t300 + (g(3) * qJ(5) + t367 * t573) * t438) * MDP(25); MDP(11) * t547 + (t369 - t567) * MDP(12) - t451 * MDP(13) + (-t320 + t545) * MDP(14) + t370 * MDP(15) + (-t374 * t389 + t452 + t548) * MDP(16) + (t372 * t389 + t446) * MDP(17) + (-t330 * t372 + t362 - t448 + t548) * MDP(18) + (pkin(3) * t319 - t555 + (t313 - t326) * t374 + (t312 - t506) * t372) * MDP(19) + (-t317 * t372 + t330 * t374 - t446 + t578) * MDP(20) + (t295 * qJ(4) - t297 * pkin(3) - t317 * t330 - t312 * t326 - g(1) * t486 - g(2) * t487 - g(3) * (-pkin(3) * t541 + t397) + t506 * t313) * MDP(21) + (t309 * t581 + t310 * t372 + (pkin(4) + t566) * t370 + t570) * MDP(22) + (t305 * t372 - t308 * t581 - t310 * t374 - t447 + t578 + t579) * MDP(23) + (t555 - t319 * t566 + (-t304 + t309) * t374 + (-t300 + t507) * t372) * MDP(24) + (t293 * qJ(4) - t292 * t566 - t300 * t309 - t305 * t310 - g(1) * (-pkin(4) * t356 + t486) - g(2) * (-pkin(4) * t354 + t487) - g(3) * (-t438 * t501 + t397) + t507 * t304) * MDP(25); (t448 - t550) * MDP(21) + (-t551 - t564 - t570) * MDP(25) + (MDP(18) + MDP(22)) * (-t370 + t547) + (MDP(20) + MDP(23)) * (-t581 ^ 2 - t369) + (-MDP(19) + MDP(24)) * t451; (-t465 - t545) * MDP(22) + (t481 - t546) * MDP(23) + (-t369 - t567) * MDP(24) + (t300 * t374 - t304 * t372 + t488 + t576) * MDP(25) + (-MDP(22) * t580 - MDP(23) * t489) * t437;];
tau = t1;
