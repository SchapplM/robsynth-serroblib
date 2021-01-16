% Calculate vector of inverse dynamics joint torques for
% S6PRPRPR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRPRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:08
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PRPRPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:06:34
% EndTime: 2021-01-16 01:06:46
% DurationCPUTime: 7.24s
% Computational Cost: add. (3317->478), mult. (7286->663), div. (0->0), fcn. (5939->22), ass. (0->227)
t608 = qJ(5) + pkin(8);
t473 = sin(qJ(2));
t476 = cos(qJ(2));
t466 = sin(pkin(6));
t553 = qJD(2) * t466;
t528 = qJD(1) * t553;
t546 = qJDD(1) * t466;
t607 = t473 * t546 + t476 * t528;
t460 = qJ(2) + pkin(11);
t455 = cos(t460);
t465 = sin(pkin(10));
t430 = t465 * t455;
t468 = cos(pkin(10));
t431 = t468 * t455;
t453 = sin(t460);
t469 = cos(pkin(6));
t567 = t468 * t469;
t575 = t465 * t469;
t597 = -g(2) * (t453 * t567 + t430) + g(1) * (t453 * t575 - t431);
t463 = sin(pkin(12));
t472 = sin(qJ(4));
t475 = cos(qJ(4));
t590 = cos(pkin(12));
t415 = t463 * t475 + t472 * t590;
t449 = pkin(4) * t475 + pkin(3);
t467 = cos(pkin(11));
t595 = pkin(2) * t467;
t425 = -t449 - t595;
t524 = t590 * t475;
t492 = -t463 * t472 + t524;
t353 = -pkin(5) * t492 - pkin(9) * t415 + t425;
t405 = t415 * qJD(4);
t545 = qJDD(2) * t472;
t508 = -qJDD(2) * t524 + t463 * t545;
t361 = qJD(2) * t405 + t508;
t360 = qJDD(6) + t361;
t459 = qJ(4) + pkin(12);
t454 = cos(t459);
t593 = g(3) * t466;
t543 = t455 * t593;
t606 = t353 * t360 + (t453 * (g(1) * t468 + g(2) * t465) - t543) * t454;
t438 = t469 * qJDD(1) + qJDD(3);
t421 = t475 * t438;
t435 = t476 * t546;
t393 = qJDD(2) * pkin(2) - t473 * t528 + t435;
t464 = sin(pkin(11));
t349 = t464 * t393 + t467 * t607;
t345 = qJDD(2) * pkin(8) + t349;
t440 = qJD(1) * t469 + qJD(3);
t486 = qJ(5) * qJDD(2) + qJD(2) * qJD(5) + qJD(4) * t440 + t345;
t554 = qJD(1) * t466;
t529 = t476 * t554;
t418 = qJD(2) * pkin(2) + t529;
t530 = t473 * t554;
t423 = t467 * t530;
t374 = t464 * t418 + t423;
t515 = qJD(2) * t608 + t374;
t502 = t515 * qJD(4);
t317 = qJDD(4) * pkin(4) - t472 * t486 - t475 * t502 + t421;
t318 = (t438 - t502) * t472 + t486 * t475;
t308 = t590 * t317 - t463 * t318;
t306 = -qJDD(4) * pkin(5) - t308;
t437 = qJD(2) * t524;
t552 = qJD(2) * t472;
t403 = t463 * t552 - t437;
t402 = qJD(6) + t403;
t406 = t415 * qJD(2);
t432 = t469 * t454;
t444 = pkin(4) * t463 + pkin(9);
t452 = sin(t459);
t579 = t453 * t466;
t510 = g(1) * t465 - g(2) * t468;
t602 = t510 * t466;
t605 = t452 * t597 + t454 * t602 + (t552 * pkin(4) + pkin(5) * t406 + pkin(9) * t403 + qJD(6) * t444) * t402 + g(3) * (-t452 * t579 + t432) + t306;
t385 = t464 * t529 + t423;
t591 = qJD(4) * pkin(4);
t601 = t472 * t591 - t385;
t399 = (t464 * t476 + t467 * t473) * t466;
t443 = t469 * t475;
t520 = -t399 * t472 + t443;
t309 = t463 * t317 + t590 * t318;
t307 = qJDD(4) * pkin(9) + t309;
t346 = t475 * t440 - t472 * t515;
t341 = t346 + t591;
t347 = t440 * t472 + t475 * t515;
t577 = t463 * t347;
t323 = t341 * t590 - t577;
t321 = -qJD(4) * pkin(5) - t323;
t422 = t464 * t530;
t373 = t418 * t467 - t422;
t359 = -qJD(2) * t449 + qJD(5) - t373;
t328 = pkin(5) * t403 - pkin(9) * t406 + t359;
t445 = pkin(2) * t464 + pkin(8);
t561 = qJ(5) + t445;
t413 = t561 * t475;
t523 = t561 * t472;
t357 = t413 * t590 - t463 * t523;
t408 = t492 * qJD(4);
t517 = qJD(4) * t561;
t384 = qJD(5) * t475 - t472 * t517;
t387 = t467 * t529 - t422;
t489 = -qJD(5) * t472 - t475 * t517;
t557 = t384 * t590 - t492 * t387 + t463 * t489;
t599 = (qJD(6) * t328 + t307) * t492 + t306 * t415 + t321 * t408 + (t469 * t510 - t593) * t453 + (-qJD(6) * t353 - t557) * t402 - t357 * t360;
t451 = pkin(6) - t460;
t596 = sin(t451);
t594 = pkin(4) * t472;
t592 = g(3) * t469;
t547 = qJD(2) * qJD(4);
t527 = t472 * t547;
t483 = qJDD(2) * t415 - t463 * t527;
t362 = qJD(4) * t437 + t483;
t471 = sin(qJ(6));
t474 = cos(qJ(6));
t548 = t474 * qJD(4);
t531 = qJD(6) * t548 + t471 * qJDD(4) + t474 * t362;
t549 = qJD(6) * t471;
t329 = -t406 * t549 + t531;
t589 = t329 * t471;
t582 = t406 * t471;
t376 = -t548 + t582;
t587 = t376 * t402;
t586 = t376 * t406;
t378 = qJD(4) * t471 + t406 * t474;
t585 = t378 * t402;
t584 = t378 * t406;
t581 = t452 * t469;
t580 = t453 * t465;
t578 = t454 * t455;
t576 = t465 * t466;
t574 = t465 * t473;
t573 = t466 * t468;
t572 = t466 * t471;
t571 = t466 * t474;
t570 = t466 * t475;
t569 = t466 * t476;
t568 = t468 * t453;
t566 = t469 * t471;
t565 = t469 * t473;
t564 = t469 * t474;
t563 = t469 * t476;
t562 = t471 * t360;
t354 = t474 * t360;
t560 = qJDD(1) - g(3);
t559 = -t329 * t492 + t378 * t405;
t339 = t590 * t347;
t324 = t463 * t341 + t339;
t558 = t384 * t463 - t415 * t387 - t489 * t590;
t556 = pkin(5) * t405 - pkin(9) * t408 + t601;
t461 = t472 ^ 2;
t555 = -t475 ^ 2 + t461;
t551 = qJD(4) * t475;
t550 = qJD(6) * t415;
t544 = qJDD(2) * t475;
t542 = cos(t451) / 0.2e1;
t540 = t415 * t562;
t539 = t415 * t354;
t536 = t471 * t578;
t535 = t474 * t578;
t534 = t454 * t566;
t533 = t454 * t564;
t532 = t468 * t563;
t522 = pkin(6) + t460;
t521 = -t474 * qJDD(4) + t362 * t471;
t500 = sin(t522) / 0.2e1;
t409 = t500 - t596 / 0.2e1;
t519 = t409 * t468 + t430;
t518 = -t409 * t465 + t431;
t516 = t402 * t474;
t348 = t393 * t467 - t464 * t607;
t509 = cos(t522);
t322 = qJD(4) * pkin(9) + t324;
t311 = t322 * t474 + t328 * t471;
t507 = t322 * t471 - t328 * t474;
t330 = t378 * qJD(6) + t521;
t506 = t330 * t492 - t376 * t405;
t364 = t399 * t475 + t469 * t472;
t332 = t364 * t590 + t463 * t520;
t503 = t464 * t473 - t467 * t476;
t398 = t503 * t466;
t505 = t332 * t474 + t398 * t471;
t504 = -t332 * t471 + t398 * t474;
t501 = t354 + (-t403 * t471 - t549) * t402;
t495 = -t465 * t563 - t468 * t473;
t494 = -t408 * t471 - t474 * t550;
t493 = -t408 * t474 + t415 * t549;
t491 = -g(1) * t576 + g(2) * t573 - t592;
t488 = t509 / 0.2e1 + t542;
t367 = -t468 * t488 + t580;
t370 = t465 * t488 + t568;
t410 = t596 / 0.2e1 + t500;
t490 = g(1) * t370 + g(2) * t367 - g(3) * t410;
t326 = t346 * t590 - t577;
t487 = -t444 * t360 + (t321 + t326) * t402;
t365 = -qJD(2) * pkin(3) - t373;
t447 = -pkin(3) - t595;
t485 = -qJDD(4) * t445 + (qJD(2) * t447 + t365 + t387) * qJD(4);
t388 = t503 * t553;
t484 = -qJD(4) * t364 + t472 * t388;
t482 = -g(1) * t495 - g(3) * t569;
t481 = g(3) * t579 - qJD(2) * t365 - t345 - t597;
t333 = pkin(4) * t527 - qJDD(2) * t449 + qJDD(5) - t348;
t477 = qJD(4) ^ 2;
t480 = -g(1) * (t455 * t575 + t568) + g(2) * (t455 * t567 - t580) - qJD(2) * t385 + t445 * t477 - t348 + t543 + (-pkin(3) + t447) * qJDD(2);
t478 = qJD(2) ^ 2;
t446 = -pkin(4) * t590 - pkin(5);
t428 = pkin(2) * t532;
t427 = qJDD(4) * t475 - t472 * t477;
t426 = qJDD(4) * t472 + t475 * t477;
t411 = t542 - t509 / 0.2e1;
t392 = t465 * t533 - t468 * t471;
t391 = t465 * t471 + t468 * t533;
t390 = t465 * t534 + t468 * t474;
t389 = t465 * t474 - t468 * t534;
t386 = qJD(2) * t399;
t383 = t454 * t579 + t581;
t381 = t452 * t571 + t455 * t566;
t380 = -t452 * t572 + t455 * t564;
t356 = t413 * t463 + t523 * t590;
t335 = qJD(4) * t520 - t388 * t475;
t331 = t364 * t463 - t520 * t590;
t325 = t346 * t463 + t339;
t320 = t335 * t590 + t463 * t484;
t319 = t335 * t463 - t484 * t590;
t315 = pkin(5) * t361 - pkin(9) * t362 + t333;
t312 = t474 * t315;
t1 = [t560 * MDP(1) + (-t348 * t398 + t349 * t399 - t373 * t386 - t374 * t388 + t438 * t469 - g(3)) * MDP(5) + (t520 * qJDD(4) + (-qJD(2) * t386 - qJDD(2) * t398) * t475 + (-t399 * t551 + (qJD(2) * t398 - qJD(4) * t469 + t388) * t472) * qJD(4)) * MDP(11) + (t398 * t545 - qJD(4) * t335 - qJDD(4) * t364 + (t386 * t472 + t398 * t551) * qJD(2)) * MDP(12) + (-qJD(4) * t319 - qJDD(4) * t331 + t361 * t398 + t386 * t403) * MDP(13) + (-qJD(4) * t320 - qJDD(4) * t332 + t362 * t398 + t386 * t406) * MDP(14) + (t319 * t406 - t320 * t403 + t331 * t362 - t332 * t361) * MDP(15) + (-t308 * t331 + t309 * t332 - t319 * t323 + t320 * t324 + t333 * t398 + t359 * t386 - g(3)) * MDP(16) + ((-qJD(6) * t505 - t320 * t471 + t386 * t474) * t402 + t504 * t360 + t319 * t376 + t331 * t330) * MDP(22) + (-(qJD(6) * t504 + t320 * t474 + t386 * t471) * t402 - t505 * t360 + t319 * t378 + t331 * t329) * MDP(23) + ((qJDD(2) * t476 - t473 * t478) * MDP(3) + (-qJDD(2) * t473 - t476 * t478) * MDP(4)) * t466; qJDD(2) * MDP(2) + (t435 - g(2) * (t532 - t574) + t482) * MDP(3) + (-g(1) * (t465 * t565 - t468 * t476) - g(2) * (-t465 * t476 - t468 * t565) - t560 * t473 * t466) * MDP(4) + (-g(2) * t428 + t373 * t385 - t374 * t387 + (g(2) * t574 + t348 * t467 + t349 * t464 + t482) * pkin(2)) * MDP(5) + (qJDD(2) * t461 + 0.2e1 * t475 * t527) * MDP(6) + 0.2e1 * (t472 * t544 - t555 * t547) * MDP(7) + t426 * MDP(8) + t427 * MDP(9) + (t472 * t485 - t475 * t480) * MDP(11) + (t472 * t480 + t475 * t485) * MDP(12) + (-qJDD(4) * t356 - t333 * t492 + t359 * t405 + t361 * t425 - t385 * t403 + t490 * t454 + (t403 * t594 - t558) * qJD(4)) * MDP(13) + (-qJDD(4) * t357 + t333 * t415 + t359 * t408 + t362 * t425 - t385 * t406 - t490 * t452 + (t406 * t594 - t557) * qJD(4)) * MDP(14) + (-g(1) * t518 - g(2) * t519 - g(3) * t411 - t308 * t415 + t309 * t492 - t323 * t408 - t324 * t405 + t356 * t362 - t357 * t361 - t403 * t557 + t406 * t558) * MDP(15) + (t309 * t357 - t308 * t356 + t333 * t425 - g(1) * (pkin(2) * t495 - t370 * t449 + t518 * t608) - g(2) * (-pkin(2) * t574 - t367 * t449 + t519 * t608 + t428) - g(3) * (pkin(2) * t569 + t410 * t449 + t411 * t608) + t601 * t359 + t557 * t324 - t558 * t323) * MDP(16) + (t329 * t415 * t474 - t378 * t493) * MDP(17) + ((-t376 * t474 - t378 * t471) * t408 + (-t589 - t330 * t474 + (t376 * t471 - t378 * t474) * qJD(6)) * t415) * MDP(18) + (-t402 * t493 + t539 + t559) * MDP(19) + (t402 * t494 + t506 - t540) * MDP(20) + (-t360 * t492 + t402 * t405) * MDP(21) + (-t507 * t405 - t312 * t492 + t356 * t330 + (g(1) * t392 - g(2) * t391) * t455 + t558 * t376 + (t556 * t402 + (t321 * t415 + t322 * t492 - t357 * t402) * qJD(6) + t606) * t474 + t599 * t471) * MDP(22) + (-t311 * t405 + t356 * t329 + (-g(1) * t390 - g(2) * t389) * t455 + t558 * t378 + ((-qJD(6) * t322 + t315) * t492 - t321 * t550 + (qJD(6) * t357 - t556) * t402 - t606) * t471 + t599 * t474) * MDP(23); (t491 + t438) * MDP(5) + t427 * MDP(11) - t426 * MDP(12) + (-qJD(4) * t405 + qJDD(4) * t492) * MDP(13) + (-qJD(4) * t408 - qJDD(4) * t415) * MDP(14) + (-t361 * t415 - t362 * t492 - t403 * t408 + t405 * t406) * MDP(15) + (t308 * t492 + t309 * t415 - t323 * t405 + t324 * t408 + t491) * MDP(16) + (-t506 - t540) * MDP(22) + (-t539 + t559) * MDP(23) + (MDP(22) * t494 + MDP(23) * t493) * t402; MDP(8) * t545 + MDP(9) * t544 + qJDD(4) * MDP(10) + (-g(3) * t443 + t472 * t481 - t475 * t602 + t421) * MDP(11) + ((-t438 + t602 + t592) * t472 + t481 * t475) * MDP(12) + (t325 * qJD(4) - t359 * t406 - g(1) * (-t452 * t518 + t454 * t576) - g(2) * (-t452 * t519 - t454 * t573) - g(3) * (-t411 * t452 + t432) + (qJDD(4) * t590 - t403 * t552) * pkin(4) + t308) * MDP(13) + (t326 * qJD(4) + t359 * t403 - g(1) * (-t452 * t576 - t454 * t518) - g(2) * (t452 * t573 - t454 * t519) - g(3) * (-t411 * t454 - t581) + (-qJDD(4) * t463 - t406 * t552) * pkin(4) - t309) * MDP(14) + ((t324 - t325) * t406 + (-t323 + t326) * t403 + (-t361 * t463 - t362 * t590) * pkin(4)) * MDP(15) + (t323 * t325 - t324 * t326 + (t309 * t463 + t308 * t590 - t359 * t552 - g(1) * (t465 * t570 - t472 * t518) - g(2) * (-t468 * t570 - t472 * t519) - g(3) * (-t411 * t472 + t443)) * pkin(4)) * MDP(16) + (t378 * t516 + t589) * MDP(17) + ((t329 - t587) * t474 + (-t330 - t585) * t471) * MDP(18) + (t402 * t516 + t562 - t584) * MDP(19) + (t501 + t586) * MDP(20) - t402 * t406 * MDP(21) + (-t325 * t376 + t446 * t330 + t507 * t406 + t487 * t471 - t474 * t605) * MDP(22) + (t311 * t406 - t325 * t378 + t446 * t329 + t471 * t605 + t487 * t474) * MDP(23) + (-t472 * t475 * MDP(6) + MDP(7) * t555) * t478; (0.2e1 * qJD(4) * t406 + t508) * MDP(13) + ((t437 - t403) * qJD(4) + t483) * MDP(14) + (-t403 ^ 2 - t406 ^ 2) * MDP(15) + (t323 * t406 + t324 * t403 + t333 - t490) * MDP(16) + (t501 - t586) * MDP(22) + (-t402 ^ 2 * t474 - t562 - t584) * MDP(23); t378 * t376 * MDP(17) + (-t376 ^ 2 + t378 ^ 2) * MDP(18) + (t531 + t587) * MDP(19) + (-t521 + t585) * MDP(20) + t360 * MDP(21) + (-t471 * t307 + t312 + t311 * t402 - t321 * t378 - g(1) * (t380 * t465 + t390 * t453 - t468 * t536) - g(2) * (-t380 * t468 + t389 * t453 - t465 * t536) - g(3) * (-t383 * t471 - t455 * t571)) * MDP(22) + (-t474 * t307 - t471 * t315 - t507 * t402 + t321 * t376 - g(1) * (-t381 * t465 + t392 * t453 - t468 * t535) - g(2) * (t381 * t468 - t391 * t453 - t465 * t535) - g(3) * (-t383 * t474 + t455 * t572)) * MDP(23) + (-MDP(19) * t582 - MDP(20) * t378 - MDP(22) * t311 + MDP(23) * t507) * qJD(6);];
tau = t1;
