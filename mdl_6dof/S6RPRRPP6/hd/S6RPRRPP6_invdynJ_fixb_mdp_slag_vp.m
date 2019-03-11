% Calculate vector of inverse dynamics joint torques for
% S6RPRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRRPP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:48:42
% EndTime: 2019-03-09 04:48:52
% DurationCPUTime: 7.10s
% Computational Cost: add. (5684->562), mult. (11029->704), div. (0->0), fcn. (6973->10), ass. (0->234)
t530 = sin(qJ(4));
t534 = cos(qJ(3));
t533 = cos(qJ(4));
t603 = t533 * qJD(3);
t668 = qJD(1) * t530 - t534 * t603;
t613 = qJD(1) * t534;
t467 = t530 * t613 - t603;
t612 = qJD(3) * t530;
t469 = t533 * t613 + t612;
t527 = sin(pkin(9));
t528 = cos(pkin(9));
t413 = t528 * t467 + t469 * t527;
t531 = sin(qJ(3));
t615 = qJD(1) * t531;
t499 = qJD(4) + t615;
t667 = t413 * t499;
t465 = t527 * t533 + t528 * t530;
t446 = t465 * qJD(1);
t464 = t527 * t530 - t528 * t533;
t607 = qJD(4) * t531;
t610 = qJD(3) * t534;
t664 = -t464 * t610 - t465 * t607 - t446;
t445 = t465 * qJD(4);
t625 = t531 * t446 + t445;
t614 = qJD(1) * t533;
t580 = t528 * t614;
t608 = qJD(4) * t530;
t581 = t527 * t608;
t589 = t530 * t615;
t606 = qJD(4) * t533;
t624 = -t527 * t589 + t528 * t606 + t531 * t580 - t581;
t532 = sin(qJ(1));
t637 = t532 * t533;
t535 = cos(qJ(1));
t638 = t531 * t535;
t451 = t530 * t638 + t637;
t560 = -t467 * t527 + t528 * t469;
t666 = t560 ^ 2;
t662 = -t530 * t610 - t531 * t606;
t665 = t668 * t527 + t528 * t662 + t531 * t581 - t580;
t566 = pkin(3) * t534 + pkin(8) * t531;
t471 = t566 * qJD(1);
t454 = t533 * t471;
t536 = -pkin(1) - pkin(7);
t494 = qJD(1) * t536 + qJD(2);
t639 = t531 * t533;
t642 = t530 * t534;
t393 = -t494 * t642 + t454 + (pkin(4) * t534 + qJ(5) * t639) * qJD(1);
t635 = t533 * t534;
t623 = t530 * t471 + t494 * t635;
t402 = qJ(5) * t589 + t623;
t654 = qJ(5) + pkin(8);
t572 = qJD(4) * t654;
t604 = qJD(5) * t533;
t439 = -t530 * t572 + t604;
t545 = -qJD(5) * t530 - t533 * t572;
t628 = (-t393 + t545) * t528 + (t402 - t439) * t527;
t476 = pkin(3) * t531 - pkin(8) * t534 + qJ(2);
t621 = t530 * t476 + t536 * t639;
t523 = g(2) * t535;
t657 = g(1) * t532;
t663 = -t523 + t657;
t661 = MDP(21) + MDP(24);
t611 = qJD(3) * t531;
t571 = -qJDD(3) * pkin(3) + t494 * t611;
t490 = qJDD(1) * t536 + qJDD(2);
t645 = t490 * t534;
t422 = t571 - t645;
t656 = g(3) * t531;
t543 = -t534 * t663 + t656;
t660 = -qJD(4) * pkin(8) * t499 - t422 + t543;
t644 = t494 * t534;
t457 = -qJD(3) * pkin(3) - t644;
t420 = pkin(4) * t467 + qJD(5) + t457;
t363 = pkin(5) * t413 - qJ(6) * t560 + t420;
t524 = qJ(4) + pkin(9);
t513 = sin(t524);
t514 = cos(t524);
t632 = t535 * t514;
t640 = t531 * t532;
t430 = t513 * t640 - t632;
t432 = t513 * t638 + t514 * t532;
t448 = t476 * qJD(1);
t475 = t531 * t494;
t456 = qJD(3) * pkin(8) + t475;
t405 = t448 * t530 + t456 * t533;
t586 = t531 * t603;
t605 = qJD(4) * t534;
t546 = -t530 * t605 - t586;
t598 = qJDD(1) * t534;
t406 = qJD(1) * t546 + qJD(4) * t603 + t530 * qJDD(3) + t533 * t598;
t463 = qJD(3) * t566 + qJD(2);
t418 = qJD(1) * t463 + qJDD(1) * t476;
t409 = t533 * t418;
t423 = qJDD(3) * pkin(8) + t490 * t531 + t494 * t610;
t601 = qJD(1) * qJD(3);
t578 = t534 * t601;
t599 = qJDD(1) * t531;
t462 = qJDD(4) + t578 + t599;
t347 = pkin(4) * t462 - qJ(5) * t406 - qJD(4) * t405 - qJD(5) * t469 - t423 * t530 + t409;
t588 = t530 * t611;
t407 = -qJD(1) * t588 + qJD(4) * t469 - t533 * qJDD(3) + t530 * t598;
t592 = -t530 * t418 - t533 * t423 - t448 * t606;
t549 = -t456 * t608 - t592;
t354 = -qJ(5) * t407 - qJD(5) * t467 + t549;
t342 = t528 * t347 - t527 * t354;
t577 = -qJDD(6) + t342;
t655 = g(3) * t534;
t659 = g(1) * t430 - g(2) * t432 - t363 * t560 + t513 * t655 + t577;
t658 = pkin(5) * t462;
t653 = pkin(1) * qJDD(1);
t538 = qJD(1) ^ 2;
t652 = qJ(2) * t538;
t391 = -qJ(5) * t467 + t405;
t387 = t528 * t391;
t404 = t533 * t448 - t456 * t530;
t390 = -qJ(5) * t469 + t404;
t360 = t390 * t527 + t387;
t651 = t360 * t560;
t650 = t391 * t527;
t649 = t406 * t530;
t648 = t467 * t499;
t647 = t469 * t499;
t646 = t469 * t533;
t643 = t499 * t530;
t641 = t530 * t535;
t636 = t532 * t534;
t634 = t533 * t535;
t633 = t534 * t535;
t537 = qJD(3) ^ 2;
t631 = t536 * t537;
t343 = t527 * t347 + t528 * t354;
t441 = t533 * t463;
t576 = -t530 * t536 + pkin(4);
t371 = qJ(5) * t586 + t441 - t621 * qJD(4) + (qJ(5) * t608 + qJD(3) * t576 - t604) * t534;
t582 = t533 * t605;
t609 = qJD(3) * t536;
t584 = t534 * t609;
t591 = t530 * t463 + t476 * t606 + t533 * t584;
t376 = -qJ(5) * t582 + (-qJD(5) * t534 + (qJ(5) * qJD(3) - qJD(4) * t536) * t531) * t530 + t591;
t351 = t527 * t371 + t528 * t376;
t630 = -qJD(6) * t465 - t475 - t624 * qJ(6) + t625 * pkin(5) + (t589 + t608) * pkin(4);
t385 = pkin(4) * t499 + t390;
t359 = t527 * t385 + t387;
t368 = t527 * t393 + t528 * t402;
t627 = pkin(5) * t613 - t628;
t365 = qJ(6) * t613 + t368;
t400 = t528 * t439 + t527 * t545;
t626 = t400 - t365;
t461 = t533 * t476;
t411 = -qJ(5) * t635 + t531 * t576 + t461;
t421 = -qJ(5) * t642 + t621;
t381 = t527 * t411 + t528 * t421;
t622 = t451 * pkin(4);
t620 = g(1) * t633 + g(2) * t636;
t619 = t535 * pkin(1) + t532 * qJ(2);
t526 = t534 ^ 2;
t618 = t531 ^ 2 - t526;
t617 = -t537 - t538;
t361 = t390 * t528 - t650;
t602 = qJD(6) - t361;
t600 = qJDD(1) * qJ(2);
t596 = 0.2e1 * qJD(1) * qJD(2);
t595 = t530 * t640;
t593 = t462 * qJ(6) + t343;
t512 = pkin(4) * t533 + pkin(3);
t520 = t535 * qJ(2);
t590 = t512 * t638 - t633 * t654 + t520;
t579 = t654 * t530;
t574 = -t490 - t523;
t573 = -g(1) * t636 + t656;
t372 = t406 * t527 + t528 * t407;
t570 = t499 * t536 + t456;
t500 = pkin(4) * t642;
t569 = -t534 * t536 + t500;
t568 = -qJD(4) * t448 - t423;
t567 = qJDD(2) - t653;
t565 = g(1) * t535 + g(2) * t532;
t563 = pkin(5) * t514 + qJ(6) * t513;
t562 = t652 + t657;
t350 = t371 * t528 - t376 * t527;
t358 = t385 * t528 - t650;
t373 = t406 * t528 - t407 * t527;
t380 = t411 * t528 - t421 * t527;
t559 = g(1) * t640 + t655;
t558 = t596 + 0.2e1 * t600;
t557 = t512 + t563;
t555 = (-pkin(4) * t530 + t536) * t657;
t554 = pkin(4) * t641 + t535 * pkin(7) + t512 * t640 - t636 * t654 + t619;
t553 = t462 * t530 + t499 * t606;
t552 = t462 * t533 - t499 * t608;
t551 = t531 * t609 + (t582 - t588) * pkin(4);
t434 = t465 * t531;
t550 = 0.2e1 * qJ(2) * t601 + qJDD(3) * t536;
t548 = pkin(4) * t407 + qJDD(5) + t571;
t544 = -pkin(8) * t462 + t457 * t499;
t542 = t558 - t565;
t483 = t654 * t533;
t424 = t483 * t527 + t528 * t579;
t425 = t528 * t483 - t527 * t579;
t505 = g(2) * t638;
t540 = -t425 * t372 + t373 * t424 - t400 * t413 + t505 - t559;
t539 = pkin(5) * t372 - qJ(6) * t373 - qJD(6) * t560 + t548;
t517 = qJDD(3) * t534;
t510 = -pkin(4) * t528 - pkin(5);
t508 = pkin(4) * t527 + qJ(6);
t503 = pkin(4) * t634;
t488 = t654 * t638;
t473 = t512 * t636;
t452 = -t530 * t532 + t531 * t634;
t450 = t531 * t637 + t641;
t449 = -t595 + t634;
t437 = -t527 * t642 + t528 * t635;
t436 = t464 * t531;
t435 = t465 * t534;
t433 = -t513 * t532 + t531 * t632;
t431 = t513 * t535 + t514 * t640;
t410 = pkin(5) * t464 - qJ(6) * t465 - t512;
t398 = t445 * t534 - t527 * t588 + t528 * t586;
t396 = qJD(3) * t434 + t464 * t605;
t389 = pkin(5) * t435 - qJ(6) * t437 + t569;
t382 = t548 - t645;
t377 = -pkin(5) * t531 - t380;
t375 = qJ(6) * t531 + t381;
t369 = pkin(4) * t469 + pkin(5) * t560 + qJ(6) * t413;
t357 = qJ(6) * t499 + t359;
t356 = -pkin(5) * t499 + qJD(6) - t358;
t355 = -pkin(5) * t396 + qJ(6) * t398 - qJD(6) * t437 + t551;
t349 = -pkin(5) * t610 - t350;
t348 = qJ(6) * t610 + qJD(6) * t531 + t351;
t344 = t539 - t645;
t341 = -t577 - t658;
t340 = qJD(6) * t499 + t593;
t1 = [(-t550 * t531 + (t558 - t631) * t534 - t620) * MDP(13) + (t550 * t534 + (t542 - t631) * t531) * MDP(12) + 0.2e1 * (-t531 * t598 + t601 * t618) * MDP(8) + (-t567 * pkin(1) - g(1) * (-pkin(1) * t532 + t520) - g(2) * t619 + (t596 + t600) * qJ(2)) * MDP(6) + (-t591 * t499 - t621 * t462 + g(1) * t451 - g(2) * t449 + (t570 * t608 + (-t457 * t533 + t469 * t536) * qJD(3) + t592) * t531 + (-qJD(3) * t405 - t406 * t536 + t422 * t533 - t457 * t608) * t534) * MDP(20) + (-g(1) * t452 - g(2) * t450 + t441 * t499 + t461 * t462 + (t467 * t609 - t570 * t606 + t409) * t531 + (qJD(3) * t404 - t407 * t536 + t457 * t606) * t534 + ((-qJD(4) * t476 - t584) * t499 + t422 * t534 + (-qJD(3) * t457 - t462 * t536 + t568) * t531) * t530) * MDP(19) + (-g(1) * t433 - g(2) * t431 - t341 * t531 + t344 * t435 - t349 * t499 + t355 * t413 - t356 * t610 - t363 * t396 + t372 * t389 - t377 * t462) * MDP(23) + (t462 * t531 + t499 * t610) * MDP(18) + ((t499 * t612 - t407) * t531 + (-qJD(3) * t467 - t553) * t534) * MDP(17) + (-t342 * t437 - t343 * t435 - t350 * t560 - t351 * t413 + t358 * t398 + t359 * t396 - t372 * t381 - t373 * t380 + t620) * MDP(21) + (-t340 * t435 + t341 * t437 - t348 * t413 + t349 * t560 - t356 * t398 + t357 * t396 - t372 * t375 + t373 * t377 + t620) * MDP(24) + (-g(1) * t432 - g(2) * t430 + t340 * t531 - t344 * t437 + t348 * t499 - t355 * t560 + t357 * t610 + t363 * t398 - t373 * t389 + t375 * t462) * MDP(25) + ((-t499 * t603 + t406) * t531 + (qJD(3) * t469 + t552) * t534) * MDP(16) + (t406 * t635 + t469 * t546) * MDP(14) + (t340 * t375 + t357 * t348 + t344 * t389 + t363 * t355 + t341 * t377 + t356 * t349 - g(1) * (pkin(5) * t433 + qJ(6) * t432 + t590) - g(2) * (pkin(5) * t431 + qJ(6) * t430 + t554) - t555) * MDP(26) + (-g(1) * t590 - g(2) * t554 + t342 * t380 + t343 * t381 + t358 * t350 + t359 * t351 + t382 * t569 + t420 * t551 - t555) * MDP(22) + qJDD(1) * MDP(1) + (qJDD(1) * t526 - 0.2e1 * t531 * t578) * MDP(7) + ((t467 * t533 + t469 * t530) * t611 + (-t649 - t407 * t533 + (t467 * t530 - t646) * qJD(4)) * t534) * MDP(15) + t565 * MDP(3) + t542 * MDP(5) + (-qJDD(3) * t531 - t534 * t537) * MDP(10) + (-t531 * t537 + t517) * MDP(9) + t663 * MDP(2) + (qJDD(2) - t663 - 0.2e1 * t653) * MDP(4); qJDD(1) * MDP(4) - t538 * MDP(5) + (t523 - t562 + t567) * MDP(6) + t517 * MDP(12) + (-t342 * t434 - t343 * t436 + t358 * t665 + t359 * t664 - t663) * MDP(22) + (-t340 * t436 + t341 * t434 - t356 * t665 + t357 * t664 - t663) * MDP(26) + (-MDP(23) * t434 - MDP(25) * t436) * t462 + (MDP(13) * t617 - t407 * MDP(19) - t406 * MDP(20) - t382 * MDP(22) - t372 * MDP(23) + t373 * MDP(25) - t344 * MDP(26)) * t534 + (t617 * MDP(12) - qJDD(3) * MDP(13) + (-t530 * MDP(19) - t533 * MDP(20)) * t462 + (MDP(19) * t467 + MDP(20) * t469 + MDP(22) * t420 + MDP(23) * t413 - MDP(25) * t560 + MDP(26) * t363) * qJD(3)) * t531 + ((-t614 + t662) * MDP(19) + (t530 * t607 + t668) * MDP(20) + t664 * MDP(25)) * t499 + (t499 * MDP(23) - t560 * t661) * t665 + t661 * (t436 * t372 + t373 * t434 - t413 * t664); MDP(9) * t598 - MDP(10) * t599 + qJDD(3) * MDP(11) + ((-t574 - t652) * t534 + t573) * MDP(12) + (t655 - t505 + (-t490 + t562) * t531) * MDP(13) + (t499 * t646 + t649) * MDP(14) + ((t406 - t648) * t533 + (-t407 - t647) * t530) * MDP(15) + ((-t469 * t534 + t499 * t639) * qJD(1) + t553) * MDP(16) + ((t467 * t534 - t531 * t643) * qJD(1) + t552) * MDP(17) - t499 * MDP(18) * t613 + (-t404 * t613 - t467 * t475 - pkin(3) * t407 - t454 * t499 + (t499 * t644 + t544) * t530 + t660 * t533) * MDP(19) + (-pkin(3) * t406 + t405 * t613 - t469 * t475 + t623 * t499 - t530 * t660 + t544 * t533) * MDP(20) + (-t342 * t465 - t343 * t464 - t358 * t624 - t359 * t625 + t368 * t413 - t560 * t628 + t540) * MDP(21) + (t343 * t425 - t342 * t424 - t382 * t512 - g(1) * (t640 * t654 + t473) - g(2) * (-t512 * t633 - t488) - g(3) * (-t512 * t531 + t534 * t654) + (pkin(4) * t643 - t475) * t420 + (t400 - t368) * t359 + t628 * t358) * MDP(22) + (t344 * t464 + t356 * t613 + t363 * t625 + t372 * t410 + t413 * t630 - t424 * t462 - t499 * t627 + t514 * t543) * MDP(23) + (-t340 * t464 + t341 * t465 + t356 * t624 - t357 * t625 + t365 * t413 + t560 * t627 + t540) * MDP(24) + (-t344 * t465 - t357 * t613 - t363 * t624 - t373 * t410 + t425 * t462 + t499 * t626 + t513 * t543 - t560 * t630) * MDP(25) + (-g(1) * t473 + g(2) * t488 + t340 * t425 + t341 * t424 + t344 * t410 + t630 * t363 + t626 * t357 + t627 * t356 + (g(3) * t557 - t654 * t657) * t531 + (-g(3) * t654 + t523 * t557 - t563 * t657) * t534) * MDP(26) + (MDP(7) * t531 * t534 - MDP(8) * t618) * t538; t469 * t467 * MDP(14) + (-t467 ^ 2 + t469 ^ 2) * MDP(15) + (t406 + t648) * MDP(16) + (-t407 + t647) * MDP(17) + t462 * MDP(18) + (-t456 * t606 - g(1) * t449 - g(2) * t451 + t405 * t499 - t457 * t469 + t409 + (t568 + t655) * t530) * MDP(19) + (g(1) * t450 - g(2) * t452 + g(3) * t635 + t404 * t499 + t457 * t467 - t549) * MDP(20) + (t359 * t560 - t651 + (-t372 * t527 - t373 * t528) * pkin(4) + (-t358 + t361) * t413) * MDP(21) + (-t359 * t361 + t358 * t360 - g(1) * t503 - g(2) * t622 + (t342 * t528 + t343 * t527 - t420 * t469 + t530 * t559) * pkin(4)) * MDP(22) + (t360 * t499 - t369 * t413 + (pkin(5) - t510) * t462 + t659) * MDP(23) + (t357 * t560 - t372 * t508 + t373 * t510 - t651 + (t356 - t602) * t413) * MDP(24) + (-t514 * t655 - g(1) * t431 + g(2) * t433 - t363 * t413 + t369 * t560 + t462 * t508 + (0.2e1 * qJD(6) - t361) * t499 + t593) * MDP(25) + (t340 * t508 + t341 * t510 - t363 * t369 - t356 * t360 - g(1) * (-pkin(4) * t595 - pkin(5) * t430 + qJ(6) * t431 + t503) - g(2) * (pkin(5) * t432 - qJ(6) * t433 + t622) - g(3) * (-t500 + (-pkin(5) * t513 + qJ(6) * t514) * t534) + t602 * t357) * MDP(26); (t358 * t560 + t359 * t413 + t548 - t573) * MDP(22) + (t499 * t560 + t372) * MDP(23) + (-t373 + t667) * MDP(25) + (-t356 * t560 + t357 * t413 + t539 - t573) * MDP(26) + (MDP(22) + MDP(26)) * t534 * t574 + t661 * (-t413 ^ 2 - t666); (t413 * t560 - t462) * MDP(23) + (t373 + t667) * MDP(24) + (-t499 ^ 2 - t666) * MDP(25) + (-t357 * t499 - t658 - t659) * MDP(26);];
tau  = t1;
