% Calculate vector of inverse dynamics joint torques for
% S6PPRRPR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
% MDP [23x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPRRPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(23,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [23 1]), ...
  'S6PPRRPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [23x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:47:44
% EndTime: 2019-03-08 18:47:54
% DurationCPUTime: 9.61s
% Computational Cost: add. (4116->498), mult. (10537->719), div. (0->0), fcn. (9799->18), ass. (0->211)
t464 = sin(pkin(7));
t587 = cos(pkin(6));
t535 = t464 * t587;
t593 = cos(qJ(3));
t465 = sin(pkin(6));
t463 = sin(pkin(12));
t470 = sin(qJ(3));
t467 = cos(pkin(12));
t586 = cos(pkin(7));
t532 = t467 * t586;
t490 = -t463 * t470 + t593 * t532;
t601 = t465 * t490;
t605 = t593 * t535 + t601;
t447 = qJD(1) * t587 + qJD(2);
t543 = t464 * t593;
t479 = -qJD(1) * t601 - t447 * t543;
t491 = t593 * t463 + t470 * t532;
t485 = t491 * t465;
t571 = t464 * t470;
t370 = qJD(1) * t485 + t447 * t571;
t469 = sin(qJ(4));
t472 = cos(qJ(4));
t520 = pkin(4) * t469 - qJ(5) * t472;
t409 = qJD(4) * t520 - qJD(5) * t469;
t604 = -t370 + t409;
t560 = qJD(3) * t472;
t449 = -qJD(6) + t560;
t462 = sin(pkin(13));
t466 = cos(pkin(13));
t554 = t466 * qJD(4);
t561 = qJD(3) * t469;
t424 = t462 * t561 - t554;
t559 = qJD(4) * t462;
t426 = t466 * t561 + t559;
t468 = sin(qJ(6));
t471 = cos(qJ(6));
t511 = t424 * t468 - t426 * t471;
t603 = t449 * t511;
t367 = qJD(3) * pkin(9) + t370;
t546 = t465 * t467 * t464;
t398 = -qJD(1) * t546 + t447 * t586;
t600 = -t469 * t367 + t398 * t472;
t558 = qJD(4) * t469;
t548 = pkin(9) * t558;
t573 = t462 * t472;
t566 = t462 * t548 + t604 * t466 - t479 * t573;
t569 = t466 * t472;
t599 = t604 * t462 + t479 * t569;
t445 = t587 * qJDD(1) + qJDD(2);
t540 = qJD(3) * t571;
t553 = qJD(1) * qJD(3);
t477 = -qJDD(1) * t601 - t445 * t543 + t447 * t540 + t485 * t553;
t585 = cos(pkin(11));
t522 = t587 * t585;
t584 = sin(pkin(11));
t410 = t463 * t522 + t467 * t584;
t483 = t463 * t584 - t467 * t522;
t534 = t465 * t585;
t597 = t464 * t534 + t483 * t586;
t351 = t410 * t470 + t593 * t597;
t521 = t587 * t584;
t411 = -t463 * t521 + t467 * t585;
t484 = t463 * t585 + t467 * t521;
t533 = t465 * t584;
t596 = -t464 * t533 + t484 * t586;
t353 = t411 * t470 + t593 * t596;
t495 = g(1) * t353 + g(2) * t351 - g(3) * t605;
t598 = t370 * qJD(3) - t477 + t495;
t595 = t491 * qJDD(1);
t594 = -qJDD(4) * pkin(4) + qJDD(5);
t589 = pkin(10) + qJ(5);
t588 = qJD(3) * pkin(3);
t581 = t479 * t469;
t577 = t426 * t468;
t373 = t471 * t424 + t577;
t580 = t373 * t449;
t395 = -qJDD(1) * t546 + t445 * t586;
t579 = t395 * t469;
t459 = pkin(13) + qJ(6);
t455 = sin(t459);
t576 = t455 * t472;
t456 = cos(t459);
t575 = t456 * t472;
t574 = t462 * t468;
t570 = t466 * t469;
t539 = qJD(3) * t593;
t340 = qJDD(3) * pkin(9) + (t445 * t470 + t447 * t539) * t464 + (t490 * t553 + t595) * t465;
t311 = qJDD(4) * qJ(5) + t340 * t472 + t579 + (qJD(5) + t600) * qJD(4);
t507 = pkin(4) * t472 + qJ(5) * t469 + pkin(3);
t321 = qJD(3) * t409 - qJDD(3) * t507 + t477;
t306 = t466 * t311 + t462 * t321;
t508 = pkin(5) * t469 - pkin(10) * t569;
t568 = -qJD(4) * t508 - t566;
t567 = (-pkin(9) * t570 - pkin(10) * t573) * qJD(4) + t599;
t348 = t472 * t367 + t469 * t398;
t346 = qJD(4) * qJ(5) + t348;
t357 = -qJD(3) * t507 + t479;
t318 = t466 * t346 + t462 * t357;
t431 = t520 * qJD(3);
t333 = t462 * t431 + t466 * t600;
t565 = -t466 * t548 + t599;
t429 = -t471 * t466 + t574;
t500 = t429 * t472;
t564 = qJD(3) * t500 - t429 * qJD(6);
t430 = t462 * t471 + t466 * t468;
t501 = t430 * t472;
t563 = -qJD(3) * t501 + t430 * qJD(6);
t400 = pkin(9) * t569 - t462 * t507;
t460 = t469 ^ 2;
t562 = -t472 ^ 2 + t460;
t557 = qJD(4) * t472;
t556 = qJD(6) * t469;
t555 = qJD(6) * t471;
t552 = qJD(3) * qJD(4);
t551 = qJDD(3) * t469;
t550 = qJDD(4) * t462;
t549 = t472 * qJDD(3);
t452 = t466 * qJDD(4);
t536 = t472 * t552;
t499 = t536 + t551;
t396 = t462 * t499 - t452;
t397 = t466 * t499 + t550;
t545 = -t468 * t396 + t471 * t397 - t424 * t555;
t544 = pkin(5) * t462 + pkin(9);
t541 = t462 * t560;
t538 = qJ(5) * t549;
t305 = -t311 * t462 + t466 * t321;
t498 = t469 * t552 - t549;
t303 = pkin(5) * t498 - pkin(10) * t397 + t305;
t304 = -pkin(10) * t396 + t306;
t531 = t471 * t303 - t468 * t304;
t317 = -t346 * t462 + t466 * t357;
t332 = t466 * t431 - t462 * t600;
t530 = t544 * t557 + t581;
t529 = t471 * t396 + t468 * t397;
t527 = t464 * t539;
t526 = qJD(4) * t539;
t519 = t468 * t303 + t471 * t304;
t313 = -pkin(5) * t560 - pkin(10) * t426 + t317;
t314 = -pkin(10) * t424 + t318;
t307 = t313 * t471 - t314 * t468;
t308 = t313 * t468 + t314 * t471;
t380 = t470 * t535 + t485;
t489 = t586 * t587 - t546;
t356 = t380 * t472 + t469 * t489;
t326 = -t356 * t462 - t466 * t605;
t327 = t356 * t466 - t462 * t605;
t518 = t326 * t471 - t327 * t468;
t517 = t326 * t468 + t327 * t471;
t423 = t466 * t507;
t378 = -pkin(10) * t570 - t423 + (-pkin(9) * t462 - pkin(5)) * t472;
t387 = -pkin(10) * t462 * t469 + t400;
t515 = t378 * t471 - t387 * t468;
t514 = t378 * t468 + t387 * t471;
t416 = t469 * t586 + t472 * t571;
t383 = -t462 * t416 - t466 * t543;
t384 = t466 * t416 - t462 * t543;
t513 = t383 * t471 - t384 * t468;
t512 = t383 * t468 + t384 * t471;
t506 = t469 * t340 + t367 * t557 - t395 * t472 + t398 * t558;
t442 = t589 * t466;
t505 = qJD(3) * t508 + qJD(5) * t462 + qJD(6) * t442 + t332;
t441 = t589 * t462;
t504 = pkin(10) * t541 + qJD(5) * t466 - qJD(6) * t441 - t333;
t474 = qJD(3) ^ 2;
t502 = qJDD(3) * t593 - t470 * t474;
t335 = -qJD(6) * t577 + t545;
t352 = t410 * t593 - t470 * t597;
t354 = t411 * t593 - t470 * t596;
t355 = t380 * t469 - t472 * t489;
t475 = t464 * t483 - t534 * t586;
t476 = t464 * t484 + t533 * t586;
t497 = g(1) * (t354 * t469 - t472 * t476) + g(2) * (t352 * t469 - t472 * t475) + g(3) * t355;
t329 = t352 * t472 + t469 * t475;
t331 = t354 * t472 + t469 * t476;
t496 = g(1) * t331 + g(2) * t329 + g(3) * t356;
t494 = g(1) * t354 + g(2) * t352 + g(3) * t380;
t344 = -qJD(4) * pkin(4) + qJD(5) - t600;
t415 = t469 * t571 - t472 * t586;
t312 = t506 + t594;
t493 = -t312 + t497;
t492 = -qJ(5) * t558 + (qJD(5) - t344) * t472;
t366 = t479 - t588;
t488 = -pkin(9) * qJDD(4) + (t366 - t479 - t588) * qJD(4);
t487 = -g(1) * t533 + g(2) * t534 - g(3) * t587;
t336 = -qJD(6) * t511 + t529;
t482 = t497 - t506;
t473 = qJD(4) ^ 2;
t478 = 0.2e1 * qJDD(3) * pkin(3) - pkin(9) * t473 + t598;
t451 = -pkin(5) * t466 - pkin(4);
t435 = t544 * t469;
t428 = qJDD(6) + t498;
t407 = t429 * t469;
t406 = t430 * t469;
t399 = -pkin(9) * t573 - t423;
t386 = qJD(4) * t416 + t469 * t527;
t385 = -qJD(4) * t415 + t472 * t527;
t372 = t380 * qJD(3);
t371 = t605 * qJD(3);
t365 = qJD(4) * t501 + t555 * t570 - t556 * t574;
t364 = -qJD(4) * t500 - t430 * t556;
t363 = t385 * t466 + t462 * t540;
t362 = -t385 * t462 + t466 * t540;
t342 = pkin(5) * t541 + t348;
t334 = pkin(5) * t424 + t344;
t325 = -qJD(4) * t355 + t371 * t472;
t324 = qJD(4) * t356 + t371 * t469;
t316 = t325 * t466 + t372 * t462;
t315 = -t325 * t462 + t372 * t466;
t309 = pkin(5) * t396 + t312;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t445 * t587 - g(3) + (t463 ^ 2 + t467 ^ 2) * t465 ^ 2 * qJDD(1)) * MDP(2) + (-qJD(3) * t372 + qJDD(3) * t605) * MDP(4) + (-qJD(3) * t371 - qJDD(3) * t380) * MDP(5) + (t605 * t549 - qJD(4) * t324 - qJDD(4) * t355 + (-t372 * t472 - t558 * t605) * qJD(3)) * MDP(11) + (-t605 * t551 - qJD(4) * t325 - qJDD(4) * t356 + (t372 * t469 - t557 * t605) * qJD(3)) * MDP(12) + (-t326 * t549 + t324 * t424 + t355 * t396 + (-t315 * t472 + t326 * t558) * qJD(3)) * MDP(13) + (t327 * t549 + t324 * t426 + t355 * t397 + (t316 * t472 - t327 * t558) * qJD(3)) * MDP(14) + (-t315 * t426 - t316 * t424 - t326 * t397 - t327 * t396) * MDP(15) + (t305 * t326 + t306 * t327 + t312 * t355 + t315 * t317 + t316 * t318 + t324 * t344 - g(3)) * MDP(16) + (-(-qJD(6) * t517 + t315 * t471 - t316 * t468) * t449 + t518 * t428 + t324 * t373 + t355 * t336) * MDP(22) + ((qJD(6) * t518 + t315 * t468 + t316 * t471) * t449 - t517 * t428 - t324 * t511 + t355 * t335) * MDP(23); (t487 + t445) * MDP(2) + (-t386 * qJD(4) - t415 * qJDD(4)) * MDP(11) + (-t385 * qJD(4) - t416 * qJDD(4)) * MDP(12) + (-t383 * t549 + t386 * t424 + t396 * t415 + (-t362 * t472 + t383 * t558) * qJD(3)) * MDP(13) + (t384 * t549 + t386 * t426 + t397 * t415 + (t363 * t472 - t384 * t558) * qJD(3)) * MDP(14) + (-t362 * t426 - t363 * t424 - t383 * t397 - t384 * t396) * MDP(15) + (t305 * t383 + t306 * t384 + t312 * t415 + t317 * t362 + t318 * t363 + t344 * t386 + t487) * MDP(16) + (-(-qJD(6) * t512 + t362 * t471 - t363 * t468) * t449 + t513 * t428 + t386 * t373 + t415 * t336) * MDP(22) + ((qJD(6) * t513 + t362 * t468 + t363 * t471) * t449 - t512 * t428 - t386 * t511 + t415 * t335) * MDP(23) + (t502 * MDP(4) + (-qJDD(3) * t470 - t474 * t593) * MDP(5) + (-t469 * t526 + t472 * t502) * MDP(11) + (-t469 * t502 - t472 * t526) * MDP(12)) * t464; qJDD(3) * MDP(3) + t598 * MDP(4) + (-t445 * t571 - t465 * t595 + t494) * MDP(5) + (qJDD(3) * t460 + 0.2e1 * t469 * t536) * MDP(6) + 0.2e1 * (t469 * t549 - t552 * t562) * MDP(7) + (qJDD(4) * t469 + t472 * t473) * MDP(8) + (qJDD(4) * t472 - t469 * t473) * MDP(9) + (t469 * t488 + t472 * t478) * MDP(11) + (-t469 * t478 + t472 * t488) * MDP(12) + (-t494 * t462 + (pkin(9) * t396 + t312 * t462 + t479 * t424 + (qJD(3) * t399 + t317) * qJD(4)) * t469 + (-t399 * qJDD(3) - t305 + (pkin(9) * t424 + t344 * t462) * qJD(4) - t566 * qJD(3) + t495 * t466) * t472) * MDP(13) + (-t494 * t466 + (pkin(9) * t397 + t312 * t466 + t479 * t426 + (-qJD(3) * t400 - t318) * qJD(4)) * t469 + (t400 * qJDD(3) + t306 + (pkin(9) * t426 + t344 * t466) * qJD(4) + t565 * qJD(3) - t495 * t462) * t472) * MDP(14) + (-t396 * t400 - t397 * t399 - t566 * t426 - t565 * t424 + (-t317 * t466 - t318 * t462) * t557 + (-t305 * t466 - t306 * t462 + t495) * t469) * MDP(15) + (t344 * t581 + t305 * t399 + t306 * t400 + t565 * t318 + t566 * t317 + (t312 * t469 + t344 * t557 - t494) * pkin(9) + t495 * t507) * MDP(16) + (-t335 * t407 - t364 * t511) * MDP(17) + (-t335 * t406 + t336 * t407 - t364 * t373 + t365 * t511) * MDP(18) + (-t335 * t472 - t364 * t449 - t407 * t428 - t511 * t558) * MDP(19) + (t336 * t472 + t365 * t449 - t373 * t558 - t406 * t428) * MDP(20) + (-t428 * t472 - t449 * t558) * MDP(21) + (t515 * t428 - t531 * t472 + t307 * t558 + t435 * t336 + t309 * t406 + t334 * t365 - g(1) * (-t353 * t575 + t354 * t455) - g(2) * (-t351 * t575 + t352 * t455) - g(3) * (t380 * t455 + t575 * t605) + (t468 * t567 + t471 * t568) * t449 + t530 * t373 + (t308 * t472 + t449 * t514) * qJD(6)) * MDP(22) + (-t514 * t428 + t519 * t472 - t308 * t558 + t435 * t335 - t309 * t407 + t334 * t364 - g(1) * (t353 * t576 + t354 * t456) - g(2) * (t351 * t576 + t352 * t456) - g(3) * (t380 * t456 - t576 * t605) + (-t468 * t568 + t471 * t567) * t449 - t530 * t511 + (t307 * t472 + t449 * t515) * qJD(6)) * MDP(23); MDP(8) * t551 + MDP(9) * t549 + qJDD(4) * MDP(10) + (qJD(4) * t348 - t366 * t561 + t482) * MDP(11) + (-t579 + (-qJD(3) * t366 - t340) * t472 + t496) * MDP(12) + (t462 * t538 - pkin(4) * t396 - t348 * t424 + t493 * t466 + (-t317 * t469 + t332 * t472 + t462 * t492) * qJD(3)) * MDP(13) + (t466 * t538 - pkin(4) * t397 - t348 * t426 - t493 * t462 + (t318 * t469 - t333 * t472 + t466 * t492) * qJD(3)) * MDP(14) + (t332 * t426 + t333 * t424 + (-qJ(5) * t396 - qJD(5) * t424 + t317 * t560 + t306) * t466 + (qJ(5) * t397 + qJD(5) * t426 + t318 * t560 - t305) * t462 - t496) * MDP(15) + (-t317 * t332 - t318 * t333 - t344 * t348 + (-t317 * t462 + t318 * t466) * qJD(5) + t493 * pkin(4) + (-t305 * t462 + t306 * t466 - t496) * qJ(5)) * MDP(16) + (t335 * t430 - t511 * t564) * MDP(17) + (-t335 * t429 - t336 * t430 - t373 * t564 + t511 * t563) * MDP(18) + (t428 * t430 - t449 * t564 + t511 * t561) * MDP(19) + (t373 * t561 - t428 * t429 + t449 * t563) * MDP(20) + t449 * MDP(21) * t561 + ((-t441 * t471 - t442 * t468) * t428 + t451 * t336 + t309 * t429 - t307 * t561 - t342 * t373 + (t468 * t504 + t471 * t505) * t449 + t563 * t334 + t497 * t456) * MDP(22) + (-(-t441 * t468 + t442 * t471) * t428 + t451 * t335 + t309 * t430 + t308 * t561 + t342 * t511 + (-t468 * t505 + t471 * t504) * t449 + t564 * t334 - t497 * t455) * MDP(23) + (-MDP(6) * t469 * t472 + MDP(7) * t562) * t474; (t462 * t551 - t452 + (-t426 + t559) * t560) * MDP(13) + (t466 * t551 + t550 + (t424 + t554) * t560) * MDP(14) + (-t424 ^ 2 - t426 ^ 2) * MDP(15) + (t317 * t426 + t318 * t424 - t482 + t594) * MDP(16) + (t336 + t603) * MDP(22) + (t335 + t580) * MDP(23); -t511 * t373 * MDP(17) + (-t373 ^ 2 + t511 ^ 2) * MDP(18) + (t545 - t580) * MDP(19) + (-t529 + t603) * MDP(20) + t428 * MDP(21) + (-t308 * t449 + t334 * t511 - g(1) * (-t331 * t455 + t353 * t456) - g(2) * (-t329 * t455 + t351 * t456) - g(3) * (-t356 * t455 - t456 * t605) + t531) * MDP(22) + (-t307 * t449 + t334 * t373 - g(1) * (-t331 * t456 - t353 * t455) - g(2) * (-t329 * t456 - t351 * t455) - g(3) * (-t356 * t456 + t455 * t605) - t519) * MDP(23) + (-MDP(19) * t577 + MDP(20) * t511 - MDP(22) * t308 - MDP(23) * t307) * qJD(6);];
tau  = t1;
