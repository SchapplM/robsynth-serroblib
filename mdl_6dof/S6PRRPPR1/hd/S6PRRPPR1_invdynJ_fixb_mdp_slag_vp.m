% Calculate vector of inverse dynamics joint torques for
% S6PRRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPPR1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S6PRRPPR1_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:02:26
% EndTime: 2019-03-08 21:02:38
% DurationCPUTime: 10.39s
% Computational Cost: add. (4928->532), mult. (11762->734), div. (0->0), fcn. (9474->18), ass. (0->240)
t533 = cos(qJ(3));
t640 = cos(pkin(11));
t571 = t640 * t533;
t502 = qJD(2) * t571;
t523 = sin(pkin(11));
t530 = sin(qJ(3));
t601 = qJD(2) * t530;
t470 = t523 * t601 - t502;
t462 = qJD(6) + t470;
t484 = t523 * t533 + t530 * t640;
t473 = t484 * qJD(2);
t522 = sin(pkin(12));
t526 = cos(pkin(12));
t449 = t526 * qJD(3) - t473 * t522;
t532 = cos(qJ(6));
t448 = qJD(3) * t522 + t473 * t526;
t529 = sin(qJ(6));
t635 = t448 * t529;
t657 = t449 * t532 - t635;
t660 = t462 * t657;
t531 = sin(qJ(2));
t525 = sin(pkin(6));
t604 = qJD(1) * t525;
t584 = t531 * t604;
t642 = qJD(3) * pkin(3);
t659 = t530 * t642 - t584;
t472 = t484 * qJD(3);
t547 = -t523 * t530 + t571;
t475 = t547 * qJD(3);
t658 = pkin(4) * t472 - qJ(5) * t475 - qJD(5) * t484 + t659;
t644 = qJ(4) + pkin(8);
t576 = qJD(3) * t644;
t463 = qJD(4) * t533 - t530 * t576;
t464 = -qJD(4) * t530 - t533 * t576;
t534 = cos(qJ(2));
t583 = t534 * t604;
t610 = t463 * t640 + t523 * t464 - t547 * t583;
t555 = -t448 * t532 - t449 * t529;
t655 = t462 * t555;
t485 = t522 * t532 + t526 * t529;
t477 = t485 * qJD(6);
t606 = t485 * t470 + t477;
t527 = cos(pkin(6));
t622 = t525 * t531;
t478 = t527 * t533 - t530 * t622;
t614 = -t610 * t522 + t526 * t658;
t613 = t522 * t658 + t610 * t526;
t611 = t463 * t523 - t640 * t464 - t484 * t583;
t641 = cos(pkin(10));
t573 = t641 * t534;
t524 = sin(pkin(10));
t624 = t524 * t531;
t466 = -t527 * t573 + t624;
t574 = t641 * t531;
t623 = t524 * t534;
t468 = t527 * t623 + t574;
t566 = g(1) * t468 + g(2) * t466;
t620 = t525 * t534;
t654 = -g(3) * t620 + t566;
t568 = qJD(2) * t644 + t584;
t603 = qJD(1) * t527;
t440 = t530 * t603 + t533 * t568;
t430 = t523 * t440;
t439 = -t530 * t568 + t533 * t603;
t383 = t439 * t640 - t430;
t589 = pkin(3) * t601;
t410 = pkin(4) * t473 + qJ(5) * t470 + t589;
t356 = t526 * t383 + t522 * t410;
t653 = qJD(5) * t526 - t356;
t355 = -t383 * t522 + t526 * t410;
t652 = -qJD(5) * t522 - t355;
t594 = qJD(2) * qJD(3);
t579 = t530 * t594;
t651 = pkin(3) * t579 + qJDD(4);
t593 = qJDD(2) * t530;
t423 = qJD(2) * t472 - qJDD(2) * t571 + t523 * t593;
t419 = qJDD(6) + t423;
t483 = t522 * t529 - t532 * t526;
t607 = t462 * t483;
t650 = -t419 * t485 + t462 * t607;
t535 = qJD(3) ^ 2;
t595 = qJD(1) * qJD(2);
t580 = t531 * t595;
t498 = t525 * t580;
t577 = qJDD(1) * t620;
t562 = t498 - t577;
t639 = qJDD(2) * pkin(2);
t649 = -pkin(8) * t535 + t525 * (-g(3) * t534 + t580) - t562 + t566 + 0.2e1 * t639;
t648 = pkin(3) * t523;
t647 = pkin(9) * t526;
t646 = g(3) * t525;
t507 = qJ(5) + t648;
t645 = pkin(9) + t507;
t643 = qJD(2) * pkin(2);
t638 = t657 * t473;
t637 = t555 * t473;
t634 = t470 * t522;
t633 = t475 * t522;
t632 = t484 * t522;
t631 = t484 * t526;
t630 = t507 * t522;
t629 = t507 * t526;
t518 = pkin(12) + qJ(6);
t514 = sin(t518);
t519 = qJ(3) + pkin(11);
t517 = cos(t519);
t628 = t514 * t517;
t516 = cos(t518);
t627 = t516 * t517;
t626 = t517 * t534;
t625 = t524 * t525;
t621 = t525 * t533;
t618 = t644 * t531;
t617 = qJDD(1) - g(3);
t591 = t527 * qJDD(1);
t501 = t533 * t591;
t454 = qJDD(2) * pkin(8) + (qJDD(1) * t531 + t534 * t595) * t525;
t538 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t603 + t454;
t553 = t568 * qJD(3);
t370 = qJDD(3) * pkin(3) - t530 * t538 - t533 * t553 + t501;
t371 = (-t553 + t591) * t530 + t538 * t533;
t342 = t523 * t370 + t640 * t371;
t339 = qJDD(3) * qJ(5) + qJD(3) * qJD(5) + t342;
t511 = pkin(3) * t533 + pkin(2);
t422 = -qJDD(2) * t511 + t562 + t651;
t424 = qJD(3) * t502 + qJDD(2) * t484 - t523 * t579;
t354 = pkin(4) * t423 - qJ(5) * t424 - qJD(5) * t473 + t422;
t335 = t526 * t339 + t522 * t354;
t616 = pkin(5) * t472 - t475 * t647 + t614;
t615 = pkin(9) * t633 - t613;
t434 = t439 + t642;
t572 = t640 * t440;
t379 = t523 * t434 + t572;
t374 = qJD(3) * qJ(5) + t379;
t461 = -qJD(2) * t511 + qJD(4) - t583;
t393 = pkin(4) * t470 - qJ(5) * t473 + t461;
t350 = t526 * t374 + t522 * t393;
t612 = pkin(5) * t633 + t611;
t420 = -pkin(4) * t547 - qJ(5) * t484 - t511;
t493 = t644 * t530;
t494 = t644 * t533;
t445 = -t523 * t493 + t494 * t640;
t376 = t522 * t420 + t526 * t445;
t467 = t527 * t574 + t623;
t609 = -t466 * t511 + t467 * t644;
t469 = -t527 * t624 + t573;
t608 = -t468 * t511 + t469 * t644;
t520 = t530 ^ 2;
t605 = -t533 ^ 2 + t520;
t602 = qJD(2) * t525;
t598 = qJD(6) * t529;
t597 = qJD(6) * t532;
t378 = t434 * t640 - t430;
t373 = -qJD(3) * pkin(4) + qJD(5) - t378;
t596 = -qJD(5) + t373;
t592 = qJDD(2) * t533;
t590 = g(3) * t622;
t406 = -t526 * qJDD(3) + t424 * t522;
t407 = qJDD(3) * t522 + t424 * t526;
t586 = -t529 * t406 + t532 * t407 + t449 * t597;
t585 = t640 * pkin(3);
t582 = t531 * t602;
t581 = t534 * t602;
t578 = t533 * t594;
t575 = t525 * t641;
t334 = -t339 * t522 + t526 * t354;
t330 = pkin(5) * t423 - pkin(9) * t407 + t334;
t331 = -pkin(9) * t406 + t335;
t570 = t532 * t330 - t529 * t331;
t349 = -t374 * t522 + t526 * t393;
t569 = t532 * t406 + t529 * t407;
t375 = t526 * t420 - t445 * t522;
t381 = t439 * t523 + t572;
t444 = t640 * t493 + t494 * t523;
t567 = (-t469 * t530 + t524 * t621) * pkin(3);
t510 = -t585 - pkin(4);
t565 = g(1) * t469 + g(2) * t467;
t564 = -t483 * t419 - t462 * t606;
t515 = sin(t519);
t563 = pkin(4) * t517 + qJ(5) * t515;
t341 = t370 * t640 - t523 * t371;
t561 = t529 * t330 + t532 * t331;
t560 = -t334 * t526 - t335 * t522;
t337 = pkin(5) * t470 - pkin(9) * t448 + t349;
t343 = pkin(9) * t449 + t350;
t332 = t337 * t532 - t343 * t529;
t333 = t337 * t529 + t343 * t532;
t359 = -pkin(5) * t547 - pkin(9) * t631 + t375;
t361 = -pkin(9) * t632 + t376;
t559 = t359 * t532 - t361 * t529;
t558 = t359 * t529 + t361 * t532;
t479 = t527 * t530 + t531 * t621;
t416 = t523 * t478 + t479 * t640;
t394 = -t416 * t522 - t526 * t620;
t395 = t416 * t526 - t522 * t620;
t557 = t394 * t532 - t395 * t529;
t556 = t394 * t529 + t395 * t532;
t554 = t478 * pkin(3);
t551 = -g(1) * t524 + g(2) * t641;
t481 = t645 * t526;
t549 = pkin(5) * t473 + qJD(6) * t481 + t470 * t647 - t652;
t480 = t645 * t522;
t548 = pkin(9) * t634 + qJD(6) * t480 - t653;
t344 = -t448 * t598 + t586;
t426 = t467 * t515 + t517 * t575;
t428 = t469 * t515 - t517 * t625;
t456 = t515 * t622 - t527 * t517;
t546 = g(1) * t428 + g(2) * t426 + g(3) * t456;
t340 = -qJDD(3) * pkin(4) + qJDD(5) - t341;
t544 = -t340 + t546;
t490 = -t583 - t643;
t543 = -qJD(2) * t490 - t454 + t565;
t542 = (-t467 * t530 - t533 * t575) * pkin(3);
t541 = t654 + t577;
t540 = t340 * t484 + t373 * t475 - t565;
t345 = -qJD(6) * t555 + t569;
t537 = -pkin(8) * qJDD(3) + (t490 + t583 - t643) * qJD(3);
t536 = qJD(2) ^ 2;
t491 = -t526 * pkin(5) + t510;
t487 = t511 * t620;
t465 = t470 ^ 2;
t457 = t515 * t527 + t517 * t622;
t438 = -qJD(3) * t479 - t530 * t581;
t437 = qJD(3) * t478 + t533 * t581;
t429 = t469 * t517 + t515 * t625;
t427 = t467 * t517 - t515 * t575;
t415 = -t478 * t640 + t479 * t523;
t414 = t483 * t484;
t413 = t485 * t484;
t405 = pkin(5) * t632 + t444;
t382 = t437 * t640 + t523 * t438;
t380 = t437 * t523 - t438 * t640;
t367 = t382 * t526 + t522 * t582;
t366 = -t382 * t522 + t526 * t582;
t364 = -pkin(5) * t634 + t381;
t363 = t475 * t485 + t597 * t631 - t598 * t632;
t362 = -t475 * t483 - t477 * t484;
t360 = -pkin(5) * t449 + t373;
t336 = t406 * pkin(5) + t340;
t1 = [t617 * MDP(1) + (qJD(3) * t438 + qJDD(3) * t478) * MDP(10) + (-qJD(3) * t437 - qJDD(3) * t479) * MDP(11) + (t380 * t473 - t382 * t470 + t415 * t424 - t416 * t423) * MDP(12) + (-t341 * t415 + t342 * t416 - t378 * t380 + t379 * t382 - g(3)) * MDP(13) + (t366 * t470 - t380 * t449 + t394 * t423 + t406 * t415) * MDP(14) + (-t367 * t470 + t380 * t448 - t395 * t423 + t407 * t415) * MDP(15) + (-t366 * t448 + t367 * t449 - t394 * t407 - t395 * t406) * MDP(16) + (t334 * t394 + t335 * t395 + t340 * t415 + t349 * t366 + t350 * t367 + t373 * t380 - g(3)) * MDP(17) + ((-qJD(6) * t556 + t366 * t532 - t367 * t529) * t462 + t557 * t419 - t380 * t657 + t415 * t345) * MDP(23) + (-(qJD(6) * t557 + t366 * t529 + t367 * t532) * t462 - t556 * t419 - t380 * t555 + t415 * t344) * MDP(24) + ((t461 * qJD(2) * MDP(13) - qJDD(2) * MDP(4) + (-t533 * MDP(10) + t530 * MDP(11) - MDP(3)) * t536) * t531 + (qJDD(2) * MDP(3) - t536 * MDP(4) + (-t579 + t592) * MDP(10) + (-t578 - t593) * MDP(11) - t422 * MDP(13)) * t534) * t525; qJDD(2) * MDP(2) + t541 * MDP(3) + (-t617 * t622 + t565) * MDP(4) + (qJDD(2) * t520 + 0.2e1 * t530 * t578) * MDP(5) + 0.2e1 * (t530 * t592 - t594 * t605) * MDP(6) + (qJDD(3) * t530 + t533 * t535) * MDP(7) + (qJDD(3) * t533 - t530 * t535) * MDP(8) + (t537 * t530 + t533 * t649) * MDP(10) + (-t530 * t649 + t537 * t533) * MDP(11) + (-t341 * t484 + t342 * t547 - t378 * t475 - t379 * t472 - t423 * t445 + t424 * t444 - t470 * t610 + t473 * t611 - t565 - t590) * MDP(12) + (t342 * t445 - t341 * t444 - t422 * t511 - g(1) * t608 - g(2) * t609 - g(3) * (t525 * t618 + t487) + t659 * t461 + t610 * t379 - t611 * t378) * MDP(13) + (-t334 * t547 + t349 * t472 + t375 * t423 + t444 * t406 + t654 * t526 * t517 + (t540 - t590) * t522 + t614 * t470 - t611 * t449) * MDP(14) + (t335 * t547 - t350 * t472 - t376 * t423 + t444 * t407 - t566 * t522 * t517 + t540 * t526 - (-t522 * t626 + t526 * t531) * t646 - t613 * t470 + t611 * t448) * MDP(15) + (-t375 * t407 - t376 * t406 + t560 * t484 + (-t349 * t526 - t350 * t522) * t475 - t614 * t448 + t613 * t449 + t654 * t515) * MDP(16) + (t335 * t376 + t334 * t375 + t340 * t444 - g(1) * (-t468 * t563 + t608) - g(2) * (-t466 * t563 + t609) - g(3) * t487 - (t534 * t563 + t618) * t646 + t611 * t373 + t613 * t350 + t614 * t349) * MDP(17) + (-t344 * t414 - t362 * t555) * MDP(18) + (-t344 * t413 + t345 * t414 + t362 * t657 + t363 * t555) * MDP(19) + (-t344 * t547 + t362 * t462 - t414 * t419 - t472 * t555) * MDP(20) + (t345 * t547 - t363 * t462 - t413 * t419 + t472 * t657) * MDP(21) + (-t419 * t547 + t462 * t472) * MDP(22) + (t559 * t419 - t570 * t547 + t332 * t472 + t405 * t345 + t336 * t413 + t360 * t363 - g(1) * (-t468 * t627 + t469 * t514) - g(2) * (-t466 * t627 + t467 * t514) - (t514 * t531 + t516 * t626) * t646 + (t529 * t615 + t532 * t616) * t462 - t612 * t657 + (t333 * t547 - t462 * t558) * qJD(6)) * MDP(23) + (-t558 * t419 + t561 * t547 - t333 * t472 + t405 * t344 - t336 * t414 + t360 * t362 - g(1) * (t468 * t628 + t469 * t516) - g(2) * (t466 * t628 + t467 * t516) - (-t514 * t626 + t516 * t531) * t646 + (-t529 * t616 + t532 * t615) * t462 - t612 * t555 + (t332 * t547 - t462 * t559) * qJD(6)) * MDP(24); MDP(7) * t593 + MDP(8) * t592 + qJDD(3) * MDP(9) + (-g(3) * t478 + t530 * t543 + t551 * t621 + t501) * MDP(10) + (g(3) * t479 + (-t525 * t551 - t591) * t530 + t543 * t533) * MDP(11) + ((t379 - t381) * t473 + (-t378 + t383) * t470 + (-t423 * t523 - t424 * t640) * pkin(3)) * MDP(12) + (-g(1) * t567 - g(2) * t542 - g(3) * t554 + t341 * t585 + t342 * t648 + t378 * t381 - t379 * t383 - t461 * t589) * MDP(13) + (-t423 * t630 - t349 * t473 + t381 * t449 + t406 * t510 + (t522 * t596 - t355) * t470 + t544 * t526) * MDP(14) + (-t423 * t629 + t350 * t473 - t381 * t448 + t407 * t510 + (t526 * t596 + t356) * t470 - t544 * t522) * MDP(15) + (-g(1) * t429 - g(2) * t427 - g(3) * t457 + t355 * t448 - t356 * t449 + (qJD(5) * t449 - t349 * t470 - t406 * t507 + t335) * t526 + (qJD(5) * t448 - t350 * t470 + t407 * t507 - t334) * t522) * MDP(16) + (t335 * t629 - t334 * t630 + t340 * t510 - t373 * t381 - g(1) * (-pkin(4) * t428 + qJ(5) * t429 + t567) - g(2) * (-t426 * pkin(4) + t427 * qJ(5) + t542) - g(3) * (-pkin(4) * t456 + qJ(5) * t457 + t554) + t653 * t350 + t652 * t349) * MDP(17) + (t344 * t485 + t555 * t607) * MDP(18) + (-t344 * t483 - t345 * t485 + t555 * t606 - t607 * t657) * MDP(19) + (t637 - t650) * MDP(20) + (t564 - t638) * MDP(21) - t462 * t473 * MDP(22) + ((-t480 * t532 - t481 * t529) * t419 + t491 * t345 + t336 * t483 - t332 * t473 + t364 * t657 + (t529 * t548 - t532 * t549) * t462 + t606 * t360 + t546 * t516) * MDP(23) + (-(-t480 * t529 + t481 * t532) * t419 + t491 * t344 + t336 * t485 + t333 * t473 + t364 * t555 + (t529 * t549 + t532 * t548) * t462 - t607 * t360 - t546 * t514) * MDP(24) + (-MDP(5) * t530 * t533 + MDP(6) * t605) * t536; (-t473 ^ 2 - t465) * MDP(12) + (-pkin(3) * t592 + t378 * t473 + t498 - t541 - t639 + t651) * MDP(13) + (t423 * t526 + t449 * t473) * MDP(14) + (-t423 * t522 - t448 * t473 - t465 * t526) * MDP(15) + (-t406 * t522 - t407 * t526) * MDP(16) + (-t373 * t473 - t560 - t654) * MDP(17) + (t564 + t638) * MDP(23) + (t637 + t650) * MDP(24) + (t379 * MDP(13) + (t448 * t522 + t449 * t526) * MDP(16) + (-t349 * t522 + t350 * t526) * MDP(17) - MDP(14) * t634) * t470; (t448 * t470 + t406) * MDP(14) + (t449 * t470 + t407) * MDP(15) + (-t448 ^ 2 - t449 ^ 2) * MDP(16) + (t349 * t448 - t350 * t449 - t544) * MDP(17) + (t345 - t655) * MDP(23) + (t344 + t660) * MDP(24); t555 * t657 * MDP(18) + (t555 ^ 2 - t657 ^ 2) * MDP(19) + (t586 - t660) * MDP(20) + (-t569 - t655) * MDP(21) + t419 * MDP(22) + (t333 * t462 + t360 * t555 - g(1) * (-t429 * t514 + t468 * t516) - g(2) * (-t427 * t514 + t466 * t516) - g(3) * (-t457 * t514 - t516 * t620) + t570) * MDP(23) + (t332 * t462 - t360 * t657 - g(1) * (-t429 * t516 - t468 * t514) - g(2) * (-t427 * t516 - t466 * t514) - g(3) * (-t457 * t516 + t514 * t620) - t561) * MDP(24) + (-MDP(20) * t635 + MDP(21) * t555 - MDP(23) * t333 - MDP(24) * t332) * qJD(6);];
tau  = t1;
