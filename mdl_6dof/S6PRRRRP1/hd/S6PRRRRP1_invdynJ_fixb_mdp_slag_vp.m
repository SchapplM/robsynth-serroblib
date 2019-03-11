% Calculate vector of inverse dynamics joint torques for
% S6PRRRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRRP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6PRRRRP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:59:34
% EndTime: 2019-03-08 23:59:46
% DurationCPUTime: 9.00s
% Computational Cost: add. (5756->513), mult. (12885->684), div. (0->0), fcn. (10011->14), ass. (0->243)
t539 = sin(qJ(3));
t540 = sin(qJ(2));
t534 = sin(pkin(6));
t640 = qJD(1) * t534;
t614 = t540 * t640;
t675 = qJD(3) * pkin(3);
t687 = -t539 * t675 + t614;
t538 = sin(qJ(4));
t542 = cos(qJ(3));
t684 = cos(qJ(4));
t615 = t684 * t542;
t572 = -t538 * t539 + t615;
t626 = qJD(3) + qJD(4);
t453 = t626 * t572;
t654 = t538 * t542;
t494 = t539 * t684 + t654;
t454 = t626 * t494;
t701 = pkin(4) * t454 - pkin(10) * t453 - t687;
t698 = pkin(9) + pkin(8);
t616 = qJD(3) * t698;
t495 = t539 * t616;
t496 = t542 * t616;
t505 = t698 * t539;
t506 = t698 * t542;
t573 = -t505 * t684 - t538 * t506;
t408 = qJD(4) * t573 - t495 * t684 - t538 * t496;
t543 = cos(qJ(2));
t613 = t543 * t640;
t462 = t572 * t613;
t700 = t408 - t462;
t598 = qJD(2) * t698 + t614;
t535 = cos(pkin(6));
t639 = qJD(1) * t535;
t459 = t539 * t639 + t542 * t598;
t448 = t538 * t459;
t458 = -t598 * t539 + t542 * t639;
t406 = t458 * t684 - t448;
t607 = qJD(4) * t684;
t699 = -pkin(3) * t607 + t406;
t627 = t535 * qJDD(1);
t510 = t542 * t627;
t631 = qJD(1) * qJD(2);
t473 = qJDD(2) * pkin(8) + (qJDD(1) * t540 + t543 * t631) * t534;
t597 = pkin(9) * qJDD(2) + t473;
t392 = qJDD(3) * pkin(3) - qJD(3) * t459 - t597 * t539 + t510;
t400 = qJD(3) * t458 + t539 * t627 + t597 * t542;
t450 = t458 + t675;
t635 = qJD(4) * t538;
t591 = -t684 * t392 + t538 * t400 + t450 * t635 + t459 * t607;
t625 = qJDD(3) + qJDD(4);
t361 = -pkin(4) * t625 + t591;
t674 = cos(pkin(11));
t600 = t674 * t540;
t533 = sin(pkin(11));
t660 = t533 * t543;
t483 = t535 * t600 + t660;
t532 = qJ(3) + qJ(4);
t527 = sin(t532);
t528 = cos(t532);
t601 = t534 * t674;
t442 = t483 * t527 + t528 * t601;
t599 = t674 * t543;
t661 = t533 * t540;
t485 = -t535 * t661 + t599;
t662 = t533 * t534;
t444 = t485 * t527 - t528 * t662;
t659 = t534 * t540;
t474 = t527 * t659 - t535 * t528;
t568 = g(1) * t444 + g(2) * t442 + g(3) * t474;
t563 = -t361 + t568;
t537 = sin(qJ(5));
t541 = cos(qJ(5));
t697 = t462 * t537 + t701 * t541;
t608 = qJD(2) * t684;
t637 = qJD(2) * t539;
t488 = t538 * t637 - t542 * t608;
t667 = t488 * t537;
t696 = -qJ(6) * t667 + t541 * qJD(6);
t490 = -qJD(2) * t654 - t539 * t608;
t439 = -pkin(4) * t490 + pkin(10) * t488;
t423 = pkin(3) * t637 + t439;
t695 = t537 * t423 + t699 * t541;
t525 = pkin(3) * t542 + pkin(2);
t440 = -pkin(4) * t572 - pkin(10) * t494 - t525;
t632 = qJD(5) * t541;
t694 = t440 * t632 + t701 * t537 + t700 * t541;
t482 = -t535 * t599 + t661;
t484 = t535 * t660 + t600;
t589 = g(1) * t484 + g(2) * t482;
t657 = t534 * t543;
t564 = -g(3) * t657 + t589;
t693 = t564 * t527;
t469 = -t538 * t505 + t506 * t684;
t457 = t541 * t469;
t580 = -qJ(6) * t453 - qJD(6) * t494;
t692 = pkin(5) * t454 - t408 * t537 + t580 * t541 + (-t457 + (qJ(6) * t494 - t440) * t537) * qJD(5) + t697;
t609 = t494 * t632;
t691 = -qJ(6) * t609 + (-qJD(5) * t469 + t580) * t537 + t694;
t420 = t541 * t423;
t529 = t541 * qJ(6);
t587 = -t490 * pkin(5) + t488 * t529;
t522 = pkin(3) * t538 + pkin(10);
t653 = -qJ(6) - t522;
t595 = qJD(5) * t653;
t690 = -t541 * t595 + t420 + t587 + (qJD(6) - t699) * t537;
t449 = t684 * t459;
t405 = t538 * t458 + t449;
t590 = pkin(3) * t635 - t405;
t560 = t494 * t657;
t648 = -qJD(1) * t560 + qJD(4) * t469 - t538 * t495 + t496 * t684;
t689 = t537 * t595 - t695 + t696;
t633 = qJD(5) * t537;
t688 = (t633 + t667) * pkin(5);
t548 = t453 * qJD(2) + t494 * qJDD(2);
t606 = t540 * t631;
t586 = -qJDD(1) * t657 + t534 * t606;
t673 = qJDD(2) * pkin(2);
t472 = t586 - t673;
t545 = qJD(3) ^ 2;
t686 = -pkin(8) * t545 + t534 * (-g(3) * t543 + t606) - t472 + t589 + t673;
t566 = t541 * t490 - t537 * t626;
t386 = -qJD(5) * t566 + t537 * t548 - t541 * t625;
t685 = t566 ^ 2;
t678 = g(3) * t534;
t677 = t541 * pkin(5);
t536 = -qJ(6) - pkin(10);
t676 = qJD(2) * pkin(2);
t593 = t541 * t626;
t385 = -qJD(5) * t593 - t490 * t633 - t537 * t625 - t541 * t548;
t672 = t385 * t537;
t403 = t450 * t684 - t448;
t395 = -pkin(4) * t626 - t403;
t671 = t395 * t488;
t629 = qJDD(2) * t539;
t585 = -qJDD(2) * t615 + t538 * t629;
t417 = qJD(2) * t454 + t585;
t413 = qJDD(5) + t417;
t670 = t413 * t541;
t463 = -t490 * t537 - t593;
t481 = qJD(5) + t488;
t669 = t463 * t481;
t668 = t566 * t481;
t666 = t488 * t541;
t665 = t494 * t537;
t664 = t494 * t541;
t663 = t528 * t537;
t658 = t534 * t542;
t656 = t537 * t413;
t655 = t537 * t543;
t652 = qJDD(1) - g(3);
t404 = t538 * t450 + t449;
t396 = pkin(10) * t626 + t404;
t477 = -qJD(2) * t525 - t613;
t416 = pkin(4) * t488 + pkin(10) * t490 + t477;
t378 = -t396 * t537 + t541 * t416;
t368 = qJ(6) * t566 + t378;
t364 = pkin(5) * t481 + t368;
t651 = -t368 + t364;
t650 = t541 * t403 + t537 * t439;
t443 = t483 * t528 - t527 * t601;
t523 = pkin(4) + t677;
t647 = -t442 * t523 - t443 * t536;
t445 = t485 * t528 + t527 * t662;
t646 = -t444 * t523 - t445 * t536;
t645 = t537 * t440 + t457;
t475 = t527 * t535 + t528 * t659;
t644 = -t474 * t523 - t475 * t536;
t602 = qJD(5) * t536;
t643 = t537 * t602 - t650 + t696;
t430 = t541 * t439;
t642 = t541 * t602 - t430 - t587 + (-qJD(6) + t403) * t537;
t530 = t539 ^ 2;
t641 = -t542 ^ 2 + t530;
t638 = qJD(2) * t534;
t636 = qJD(2) * t540;
t634 = qJD(5) * t481;
t630 = qJD(2) * qJD(3);
t628 = qJDD(2) * t542;
t619 = t534 * t655;
t618 = t541 * t657;
t611 = t534 * t636;
t610 = t543 * t638;
t605 = t539 * t630;
t604 = t543 * t630;
t603 = pkin(5) * t537 + t698;
t596 = t378 * t490 + t395 * t633;
t594 = t481 * t541;
t524 = -pkin(3) * t684 - pkin(4);
t588 = g(1) * t485 + g(2) * t483;
t379 = t396 * t541 + t416 * t537;
t369 = -qJ(6) * t463 + t379;
t582 = -t364 * t541 - t369 * t537;
t581 = -t522 * t413 + t671;
t546 = qJD(2) ^ 2;
t579 = qJDD(2) * t543 - t540 * t546;
t578 = -g(1) * t533 + g(2) * t674;
t486 = t535 * t542 - t539 * t659;
t487 = t535 * t539 + t540 * t658;
t427 = t538 * t486 + t487 * t684;
t414 = -t427 * t537 - t618;
t576 = -t427 * t541 + t619;
t575 = t523 * t528 - t527 * t536 + t525;
t574 = t486 * t684 - t538 * t487;
t571 = t453 * t537 + t609;
t570 = t453 * t541 - t494 * t633;
t555 = t538 * t392 + t400 * t684 + t450 * t607 - t459 * t635;
t360 = pkin(10) * t625 + t555;
t518 = pkin(3) * t605;
t376 = -pkin(3) * t628 + t417 * pkin(4) - pkin(10) * t548 + t472 + t518;
t569 = -t541 * t360 - t537 * t376 + t396 * t633 - t416 * t632;
t567 = g(1) * t445 + g(2) * t443 + g(3) * t475;
t565 = -t379 * t490 + t395 * t632 - t537 * t563;
t498 = -t613 - t676;
t561 = -qJD(2) * t498 - t473 + t588;
t559 = t564 * t528;
t375 = t541 * t376;
t557 = -qJD(5) * t379 - t537 * t360 + t375;
t556 = -pkin(8) * qJDD(3) + (t498 + t613 - t676) * qJD(3);
t355 = t386 * pkin(5) + qJDD(6) + t361;
t554 = t477 * t490 + t568 - t591;
t553 = ((-t385 - t669) * t541 + (-t386 + t668) * t537) * MDP(20) + (-t566 * t594 - t672) * MDP(19) + (-t481 ^ 2 * t537 - t463 * t490 + t670) * MDP(22) + (t481 * t594 - t490 * t566 + t656) * MDP(21) + (t488 * t626 + t548) * MDP(14) + (-t585 + (-qJD(2) * t494 - t490) * t626) * MDP(15) + (-t488 ^ 2 + t490 ^ 2) * MDP(13) + t625 * MDP(16) + (-MDP(12) * t488 + MDP(23) * t481) * t490;
t351 = pkin(5) * t413 + qJ(6) * t385 + qJD(6) * t566 + t557;
t353 = -qJ(6) * t386 - qJD(6) * t463 - t569;
t551 = qJD(5) * t582 - t351 * t537 + t353 * t541 - t364 * t666 - t369 * t667 - t567;
t550 = t477 * t488 - t555 + t567;
t549 = -g(1) * (-t445 * t537 + t484 * t541) - g(2) * (-t443 * t537 + t482 * t541) - g(3) * (-t475 * t537 - t618);
t502 = pkin(10) * t541 + t529;
t501 = t536 * t537;
t492 = t522 * t541 + t529;
t491 = t653 * t537;
t460 = t463 ^ 2;
t456 = -qJD(3) * t487 - t539 * t610;
t455 = qJD(3) * t486 + t542 * t610;
t441 = -qJDD(2) * t525 + t518 + t586;
t434 = t541 * t440;
t387 = -qJ(6) * t665 + t645;
t384 = -pkin(5) * t572 - t469 * t537 - t494 * t529 + t434;
t383 = t463 * pkin(5) + qJD(6) + t395;
t382 = qJD(4) * t427 + t538 * t455 - t456 * t684;
t381 = qJD(4) * t574 + t455 * t684 + t538 * t456;
t366 = qJD(5) * t576 - t381 * t537 + t541 * t611;
t365 = qJD(5) * t414 + t381 * t541 + t537 * t611;
t1 = [t652 * MDP(1) + (qJD(3) * t456 + qJDD(3) * t486) * MDP(10) + (-qJD(3) * t455 - qJDD(3) * t487) * MDP(11) + (-t382 * t626 + t574 * t625) * MDP(17) + (-t381 * t626 - t427 * t625 - qJDD(2) * t560 + (-t453 * t543 - t540 * t490) * t638) * MDP(18) + (t366 * t481 + t382 * t463 - t386 * t574 + t413 * t414) * MDP(24) + (-t365 * t481 - t382 * t566 + t385 * t574 + t413 * t576) * MDP(25) + (-t365 * t463 + t366 * t566 + t385 * t414 + t386 * t576) * MDP(26) + (t351 * t414 - t353 * t576 - t355 * t574 + t364 * t366 + t365 * t369 + t382 * t383 - g(3)) * MDP(27) + (t579 * MDP(3) + (-qJDD(2) * t540 - t543 * t546) * MDP(4) + (-t539 * t604 + t542 * t579) * MDP(10) + (-t539 * t579 - t542 * t604) * MDP(11) + (-t417 * t543 + t488 * t636) * MDP(17)) * t534; qJDD(2) * MDP(2) + (t652 * t657 + t589) * MDP(3) + (-t652 * t659 + t588) * MDP(4) + (qJDD(2) * t530 + 0.2e1 * t542 * t605) * MDP(5) + 0.2e1 * (t539 * t628 - t630 * t641) * MDP(6) + (qJDD(3) * t539 + t542 * t545) * MDP(7) + (qJDD(3) * t542 - t539 * t545) * MDP(8) + (t556 * t539 + t542 * t686) * MDP(10) + (-t539 * t686 + t556 * t542) * MDP(11) + (-t490 * t453 + t494 * t548) * MDP(12) + (-t494 * t417 - t453 * t488 + t490 * t454 + t548 * t572) * MDP(13) + (t453 * t626 + t494 * t625) * MDP(14) + (-t454 * t626 + t572 * t625) * MDP(15) + (-t525 * t417 - t441 * t572 + t477 * t454 - t488 * t687 + t573 * t625 - t626 * t648 + t559) * MDP(17) + (t441 * t494 + t477 * t453 - t469 * t625 + t687 * t490 - t525 * t548 - t700 * t626 - t693) * MDP(18) + (-t385 * t664 - t566 * t570) * MDP(19) + ((-t463 * t541 + t537 * t566) * t453 + (t672 - t386 * t541 + (t463 * t537 + t541 * t566) * qJD(5)) * t494) * MDP(20) + (t385 * t572 + t413 * t664 - t454 * t566 + t481 * t570) * MDP(21) + (t386 * t572 - t454 * t463 - t481 * t571 - t494 * t656) * MDP(22) + (-t413 * t572 + t454 * t481) * MDP(23) + (-t375 * t572 + t378 * t454 - t573 * t386 + t434 * t413 + t697 * t481 + t648 * t463 + (t559 + (t395 * t494 + t396 * t572 - t469 * t481) * qJD(5)) * t541 + ((-qJD(5) * t440 - t408) * t481 - t469 * t413 - (-qJD(5) * t416 - t360) * t572 + t361 * t494 + t395 * t453 - g(3) * t659 - t588) * t537) * MDP(24) + (-t645 * t413 - t569 * t572 - t379 * t454 + t573 * t385 + t361 * t664 - g(1) * (t484 * t663 + t485 * t541) - g(2) * (t482 * t663 + t483 * t541) - (-t528 * t655 + t540 * t541) * t678 + (t469 * t633 - t694) * t481 - t648 * t566 + t570 * t395) * MDP(25) + (t384 * t385 - t386 * t387 + t692 * t566 - t691 * t463 + t582 * t453 + t693 + (-t351 * t541 - t353 * t537 + (t364 * t537 - t369 * t541) * qJD(5)) * t494) * MDP(26) + (t353 * t387 + t351 * t384 + t355 * (pkin(5) * t665 - t573) - g(1) * (-t484 * t575 + t485 * t603) - g(2) * (-t482 * t575 + t483 * t603) - (t540 * t603 + t543 * t575) * t678 + (pkin(5) * t571 + t648) * t383 + t691 * t369 + t692 * t364) * MDP(27); (t524 * t386 - t420 * t481 + t590 * t463 + (t481 * t699 + t581) * t537 + (-t522 * t634 + t563) * t541 + t596) * MDP(24) + (-t524 * t385 + t581 * t541 - t590 * t566 + (t522 * t633 + t695) * t481 + t565) * MDP(25) + (t353 * t492 + t351 * t491 + t355 * (t524 - t677) - g(1) * ((-t485 * t539 + t533 * t658) * pkin(3) + t646) - g(2) * ((-t483 * t539 - t542 * t601) * pkin(3) + t647) - g(3) * (pkin(3) * t486 + t644) + (t688 + t590) * t383 + t689 * t369 - t690 * t364) * MDP(27) + (t405 * t626 + (-t488 * t637 + t625 * t684 - t626 * t635) * pkin(3) + t554) * MDP(17) + (g(3) * t487 + (-t534 * t578 - t627) * t539 + t561 * t542) * MDP(11) + (t406 * t626 + (t490 * t637 - t538 * t625 - t607 * t626) * pkin(3) + t550) * MDP(18) + (t385 * t491 - t386 * t492 - t463 * t689 - t566 * t690 + t551) * MDP(26) + (-g(3) * t486 + t539 * t561 + t578 * t658 + t510) * MDP(10) + MDP(7) * t629 + MDP(8) * t628 + qJDD(3) * MDP(9) + t553 + (-MDP(5) * t539 * t542 + MDP(6) * t641) * t546; (pkin(4) * t385 + t650 * t481 + t404 * t566 + t395 * t666 + (t481 * t633 - t670) * pkin(10) + t565) * MDP(25) + (t353 * t502 + t351 * t501 - t355 * t523 - g(1) * t646 - g(2) * t647 - g(3) * t644 + (-t404 + t688) * t383 + t643 * t369 + t642 * t364) * MDP(27) + (-pkin(4) * t386 - t404 * t463 - t430 * t481 + (-pkin(10) * t413 + t403 * t481 + t671) * t537 + (-pkin(10) * t634 + t563) * t541 + t596) * MDP(24) + (t403 * t626 + t550) * MDP(18) + (t385 * t501 - t386 * t502 - t463 * t643 + t566 * t642 + t551) * MDP(26) + t553 + (t404 * t626 + t554) * MDP(17); -t566 * t463 * MDP(19) + (-t460 + t685) * MDP(20) + (-t385 + t669) * MDP(21) + (-t386 - t668) * MDP(22) + t413 * MDP(23) + (t379 * t481 + t395 * t566 + t549 + t557) * MDP(24) + (t378 * t481 + t395 * t463 - g(1) * (-t445 * t541 - t484 * t537) - g(2) * (-t443 * t541 - t482 * t537) - g(3) * (-t475 * t541 + t619) + t569) * MDP(25) + (pkin(5) * t385 - t463 * t651) * MDP(26) + (t651 * t369 + (t383 * t566 + t351 + t549) * pkin(5)) * MDP(27); (-t460 - t685) * MDP(26) + (-t364 * t566 + t369 * t463 + t355 - t568) * MDP(27);];
tau  = t1;
