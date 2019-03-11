% Calculate vector of inverse dynamics joint torques for
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRPPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PRRPPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:21:19
% EndTime: 2019-03-08 21:21:29
% DurationCPUTime: 8.37s
% Computational Cost: add. (3125->548), mult. (7025->739), div. (0->0), fcn. (5284->14), ass. (0->241)
t523 = sin(pkin(6));
t529 = sin(qJ(2));
t532 = cos(qJ(2));
t615 = qJD(1) * qJD(2);
t594 = t532 * t615;
t525 = cos(pkin(6));
t624 = qJD(3) * t525;
t662 = qJDD(2) * pkin(8);
t685 = t662 + (qJDD(1) * t529 + t594) * t523 + qJD(1) * t624;
t671 = pkin(4) + pkin(8);
t528 = sin(qJ(3));
t589 = -qJ(4) * t528 - pkin(2);
t629 = qJD(2) * t528;
t504 = qJD(6) + t629;
t524 = cos(pkin(11));
t521 = sin(pkin(11));
t626 = qJD(3) * t521;
t531 = cos(qJ(3));
t628 = qJD(2) * t531;
t464 = t524 * t628 + t626;
t625 = qJD(3) * t524;
t466 = -t521 * t628 + t625;
t527 = sin(qJ(6));
t530 = cos(qJ(6));
t563 = t464 * t527 - t466 * t530;
t684 = t504 * t563;
t469 = t521 * t530 + t524 * t527;
t453 = t469 * qJD(6);
t552 = t469 * t528;
t638 = -qJD(2) * t552 - t453;
t614 = qJD(2) * qJD(3);
t592 = t528 * t614;
t610 = qJDD(2) * t531;
t683 = t592 - t610;
t591 = t531 * t614;
t611 = qJDD(2) * t528;
t551 = t591 + t611;
t468 = qJDD(6) + t551;
t597 = t524 * t629;
t619 = qJD(6) * t530;
t620 = qJD(6) * t527;
t652 = t521 * t527;
t637 = -t521 * t620 + t524 * t619 + t530 * t597 - t629 * t652;
t682 = -t468 * t469 - t504 * t637;
t634 = qJD(1) * t523;
t603 = t529 * t634;
t478 = qJD(2) * pkin(8) + t603;
t633 = qJD(1) * t525;
t423 = t528 * t478 - t531 * t633;
t681 = qJD(4) + t423;
t623 = qJD(3) * t528;
t509 = pkin(3) * t623;
t664 = qJ(4) * t531;
t569 = qJ(5) * t528 - t664;
t621 = qJD(4) * t528;
t539 = qJD(3) * t569 - qJD(5) * t531 - t621;
t410 = t509 + t539;
t622 = qJD(3) * t531;
t475 = t671 * t622;
t645 = t528 * t532;
t640 = -t410 * t521 + t524 * t475 - (-t521 * t529 + t524 * t645) * t634;
t639 = t524 * t410 + t521 * t475 - (t521 * t645 + t524 * t529) * t634;
t424 = t531 * t478 + t528 * t633;
t518 = qJD(3) * qJ(4);
t414 = -t518 - t424;
t412 = -qJD(3) * pkin(3) + t681;
t617 = pkin(4) * t629 + t681;
t678 = -MDP(10) + MDP(13);
t677 = MDP(11) - MDP(14);
t632 = qJD(1) * t532;
t602 = t523 * t632;
t480 = -pkin(3) * t531 + t589;
t631 = qJD(2) * t480;
t425 = -t602 + t631;
t676 = t425 * t629 + qJDD(4);
t675 = t478 * t622 + t685 * t528;
t595 = t529 * t615;
t647 = t523 * t532;
t570 = -qJDD(1) * t647 + t523 * t595;
t667 = cos(pkin(10));
t586 = t667 * t532;
t522 = sin(pkin(10));
t651 = t522 * t529;
t449 = -t525 * t586 + t651;
t587 = t667 * t529;
t650 = t522 * t532;
t451 = t525 * t650 + t587;
t576 = g(1) * t451 + g(2) * t449;
t533 = qJD(3) ^ 2;
t670 = pkin(8) * t533;
t674 = 0.2e1 * qJDD(2) * pkin(2) + t523 * (-g(3) * t532 + t595) - t570 + t576 - t670;
t409 = pkin(4) * t628 + t424;
t390 = qJD(5) + t518 + t409;
t526 = -pkin(3) - qJ(5);
t673 = -t526 * t622 + t528 * (qJD(5) - t390);
t555 = -qJ(4) * t622 - t621;
t557 = pkin(3) * t592 + t570;
t612 = qJDD(2) * t480;
t373 = qJD(2) * t555 + t557 + t612;
t446 = t509 + t555;
t546 = g(3) * t647 - t576;
t672 = qJD(2) * (-t446 + t603) - t373 - t546 - t612 - t670;
t669 = -pkin(9) + t526;
t668 = qJD(2) * pkin(2);
t666 = pkin(8) * qJDD(3);
t661 = qJDD(3) * pkin(3);
t657 = t466 * t527;
t387 = t530 * t464 + t657;
t660 = t387 * t504;
t659 = t449 * t531;
t658 = t451 * t531;
t515 = pkin(11) + qJ(6);
t511 = sin(t515);
t655 = t511 * t528;
t512 = cos(t515);
t654 = t512 * t528;
t519 = t528 ^ 2;
t534 = qJD(2) ^ 2;
t653 = t519 * t534;
t649 = t523 * t529;
t648 = t523 * t531;
t646 = t524 * t531;
t644 = t531 * t532;
t643 = qJDD(1) - g(3);
t613 = qJDD(1) * t525;
t547 = -t531 * t613 + t675;
t544 = qJDD(4) + t547;
t359 = pkin(4) * t551 - qJD(3) * qJD(5) + qJDD(3) * t526 + t544;
t461 = t526 * t531 + t589;
t364 = qJD(2) * t539 + qJDD(2) * t461 + t557;
t351 = t521 * t359 + t524 * t364;
t558 = -pkin(9) * t521 * t528 + pkin(5) * t531;
t642 = qJD(3) * t558 + t640;
t641 = -pkin(9) * t524 * t623 - t639;
t386 = qJD(3) * t526 + t617;
t407 = qJD(2) * t461 - t602;
t363 = t521 * t386 + t524 * t407;
t510 = pkin(3) * t629;
t441 = qJD(2) * t569 + t510;
t372 = t521 * t409 + t524 * t441;
t489 = t671 * t528;
t395 = t524 * t461 + t521 * t489;
t490 = t671 * t531;
t520 = t531 ^ 2;
t636 = t519 - t520;
t635 = t519 + t520;
t630 = qJD(2) * t523;
t627 = qJD(3) * t424;
t604 = -pkin(5) * t524 - pkin(4);
t618 = -t604 * t629 + t681;
t609 = t528 * t649;
t608 = t528 * t531 * t534;
t426 = qJDD(3) * t521 - t524 * t683;
t579 = qJDD(3) * t524 - t521 * t610;
t427 = t521 * t592 + t579;
t607 = -t527 * t426 + t530 * t427 - t464 * t619;
t606 = -pkin(3) * t659 + t449 * t589;
t605 = -pkin(3) * t658 + t451 * t589;
t601 = t531 * t632;
t600 = t466 * t629;
t599 = t529 * t630;
t598 = t532 * t630;
t590 = t526 * t611;
t588 = t523 * t667;
t450 = t525 * t587 + t650;
t401 = t450 * t528 + t531 * t588;
t402 = t450 * t531 - t528 * t588;
t584 = -t401 * pkin(3) + qJ(4) * t402;
t452 = -t525 * t651 + t586;
t403 = t452 * t528 - t522 * t648;
t404 = t522 * t523 * t528 + t452 * t531;
t583 = -t403 * pkin(3) + qJ(4) * t404;
t457 = -t525 * t531 + t609;
t458 = t525 * t528 + t529 * t648;
t582 = -t457 * pkin(3) + qJ(4) * t458;
t350 = t524 * t359 - t364 * t521;
t346 = pkin(5) * t551 - pkin(9) * t427 + t350;
t349 = -pkin(9) * t426 + t351;
t581 = t530 * t346 - t527 * t349;
t362 = t524 * t386 - t407 * t521;
t371 = t524 * t409 - t441 * t521;
t580 = t530 * t426 + t527 * t427;
t578 = t478 * t623 - t528 * t613 - t685 * t531;
t577 = pkin(2) * t647 + pkin(8) * t649 + (pkin(3) * t644 + qJ(4) * t645) * t523;
t575 = g(1) * t452 + g(2) * t450;
t562 = -t524 * t530 + t652;
t574 = -t562 * t468 + t504 * t638;
t572 = t363 * t629 + t350;
t571 = -t362 * t629 + t351;
t568 = t527 * t346 + t530 * t349;
t353 = pkin(5) * t629 - pkin(9) * t466 + t362;
t356 = -pkin(9) * t464 + t363;
t347 = t353 * t530 - t356 * t527;
t348 = t353 * t527 + t356 * t530;
t472 = t524 * t489;
t379 = pkin(5) * t528 + t472 + (pkin(9) * t531 - t461) * t521;
t384 = -pkin(9) * t646 + t395;
t567 = t379 * t530 - t384 * t527;
t566 = t379 * t527 + t384 * t530;
t399 = t457 * t524 + t521 * t647;
t400 = t457 * t521 - t524 * t647;
t565 = t399 * t530 - t400 * t527;
t564 = t399 * t527 + t400 * t530;
t516 = qJDD(3) * qJ(4);
t517 = qJD(3) * qJD(4);
t366 = -t516 - t517 + t578;
t476 = t669 * t521;
t554 = qJD(2) * t558 + qJD(5) * t524 + qJD(6) * t476 + t371;
t477 = t669 * t524;
t553 = pkin(9) * t597 + qJD(5) * t521 - qJD(6) * t477 + t372;
t438 = t562 * t531;
t354 = -t466 * t620 + t607;
t549 = g(1) * t403 + g(2) * t401 + g(3) * t457;
t548 = -g(1) * t404 - g(2) * t402 - g(3) * t458;
t360 = -pkin(4) * t683 + qJDD(5) - t366;
t545 = t360 + t548;
t543 = g(3) * t649 - t360 * t531 + t575;
t355 = -qJD(6) * t563 + t580;
t479 = -t602 - t668;
t541 = -t666 + (t479 + t602 - t668) * qJD(3);
t540 = t666 + (-t425 - t602 - t631) * qJD(3);
t538 = qJD(3) * t423 + t548 - t578;
t537 = -t547 + t549;
t368 = t544 - t661;
t535 = -t366 * t531 + t368 * t528 + (t412 * t531 + t414 * t528) * qJD(3) - t575;
t505 = pkin(5) * t521 + qJ(4);
t474 = t671 * t623;
t473 = -qJ(4) * t628 + t510;
t448 = pkin(5) * t646 + t490;
t439 = t469 * t531;
t435 = (-pkin(8) + t604) * t623;
t406 = qJD(3) * t458 + t528 * t598;
t405 = -qJD(3) * t609 + (t598 + t624) * t531;
t394 = -t461 * t521 + t472;
t381 = t453 * t531 - t562 * t623;
t380 = qJD(3) * t552 + qJD(6) * t438;
t378 = t406 * t521 + t524 * t599;
t377 = t406 * t524 - t521 * t599;
t374 = pkin(5) * t464 + t390;
t352 = pkin(5) * t426 + t360;
t1 = [t643 * MDP(1) + (-t366 * t458 + t368 * t457 - t405 * t414 + t406 * t412 - g(3)) * MDP(15) + (t405 * t464 + t426 * t458) * MDP(16) + (t405 * t466 + t427 * t458) * MDP(17) + (-t377 * t466 - t378 * t464 - t399 * t427 - t400 * t426) * MDP(18) + (t350 * t399 + t351 * t400 + t360 * t458 + t362 * t377 + t363 * t378 + t390 * t405 - g(3)) * MDP(19) + ((-qJD(6) * t564 + t377 * t530 - t378 * t527) * t504 + t565 * t468 + t405 * t387 + t458 * t355) * MDP(25) + (-(qJD(6) * t565 + t377 * t527 + t378 * t530) * t504 - t564 * t468 - t405 * t563 + t458 * t354) * MDP(26) + (t458 * t531 * MDP(12) + (MDP(12) * t457 + MDP(16) * t399 - MDP(17) * t400) * t528) * qJDD(2) + ((t405 * t531 + t406 * t528 + t457 * t622 - t458 * t623) * MDP(12) + (t377 * t528 + t399 * t622) * MDP(16) + (-t378 * t528 - t400 * t622) * MDP(17)) * qJD(2) + ((qJD(2) * t425 * MDP(15) - qJDD(2) * MDP(4) + (t528 * t677 + t531 * t678 - MDP(3)) * t534) * t529 + (-t373 * MDP(15) + qJDD(2) * MDP(3) - t534 * MDP(4) - t677 * t551 + t678 * t683) * t532) * t523 - t677 * (qJD(3) * t405 + qJDD(3) * t458) + t678 * (qJD(3) * t406 + qJDD(3) * t457); qJDD(2) * MDP(2) + (t643 * t647 + t576) * MDP(3) + (-t643 * t649 + t575) * MDP(4) + (qJDD(2) * t519 + 0.2e1 * t528 * t591) * MDP(5) + 0.2e1 * (t528 * t610 - t614 * t636) * MDP(6) + (qJDD(3) * t528 + t531 * t533) * MDP(7) + (qJDD(3) * t531 - t528 * t533) * MDP(8) + (t541 * t528 + t531 * t674) * MDP(10) + (-t528 * t674 + t541 * t531) * MDP(11) + (t635 * t662 + (-g(3) * t529 - t594 * t635) * t523 + t535) * MDP(12) + (t540 * t528 - t531 * t672) * MDP(13) + (t528 * t672 + t540 * t531) * MDP(14) + (t373 * t480 + t425 * t446 - g(1) * t605 - g(2) * t606 - g(3) * t577 + (-t425 * t529 + (-t412 * t528 + t414 * t531) * t532) * t634 + t535 * pkin(8)) * MDP(15) + (t490 * t426 - t474 * t464 + (-t464 * t602 + (qJD(2) * t394 + t362) * qJD(3)) * t531 - t543 * t524 + (qJD(2) * t640 + t394 * qJDD(2) - t390 * t625 - t521 * t546 + t350) * t528) * MDP(16) + (t490 * t427 - t474 * t466 + (-t466 * t602 + (-qJD(2) * t395 - t363) * qJD(3)) * t531 + t543 * t521 + (-qJD(2) * t639 - t395 * qJDD(2) + t390 * t626 - t524 * t546 - t351) * t528) * MDP(17) + (-t394 * t427 - t395 * t426 - t640 * t466 - t639 * t464 + (-t362 * t521 + t363 * t524) * t623 + (t350 * t521 - t351 * t524 - t546) * t531) * MDP(18) + (t351 * t395 + t350 * t394 + t360 * t490 - g(1) * (-qJ(5) * t658 + t452 * t671 + t605) - g(2) * (-qJ(5) * t659 + t450 * t671 + t606) - g(3) * ((pkin(4) * t529 + qJ(5) * t644) * t523 + t577) + (-t523 * t601 - t474) * t390 + t639 * t363 + t640 * t362) * MDP(19) + (-t354 * t439 - t380 * t563) * MDP(20) + (t354 * t438 + t355 * t439 - t380 * t387 - t381 * t563) * MDP(21) + (t354 * t528 + t380 * t504 - t439 * t468 - t563 * t622) * MDP(22) + (-t355 * t528 + t381 * t504 - t387 * t622 + t438 * t468) * MDP(23) + (t468 * t528 + t504 * t622) * MDP(24) + (t567 * t468 + t581 * t528 + t347 * t622 + t435 * t387 + t448 * t355 - t352 * t438 - t374 * t381 - g(1) * (-t451 * t655 + t452 * t512) - g(2) * (-t449 * t655 + t450 * t512) + (t527 * t641 + t530 * t642) * t504 + (-t348 * t528 - t504 * t566) * qJD(6) + (-t387 * t601 - g(3) * (t511 * t645 + t512 * t529)) * t523) * MDP(25) + (-t566 * t468 - t568 * t528 - t348 * t622 - t435 * t563 + t448 * t354 - t352 * t439 + t374 * t380 - g(1) * (-t451 * t654 - t452 * t511) - g(2) * (-t449 * t654 - t450 * t511) + (-t527 * t642 + t530 * t641) * t504 + (-t347 * t528 - t504 * t567) * qJD(6) + (t563 * t601 - g(3) * (-t511 * t529 + t512 * t645)) * t523) * MDP(26); -MDP(5) * t608 + t636 * t534 * MDP(6) + MDP(7) * t611 + MDP(8) * t610 + qJDD(3) * MDP(9) + (-t479 * t629 + t537 + t627) * MDP(10) + (-t479 * t628 - t538) * MDP(11) + (-pkin(3) * t528 + t664) * qJDD(2) * MDP(12) + (-0.2e1 * t661 - t627 + (-qJD(2) * t473 - t613) * t531 - t549 + t675 + t676) * MDP(13) + (0.2e1 * t516 + 0.2e1 * t517 + (t425 * t531 + t473 * t528) * qJD(2) + t538) * MDP(14) + (-t368 * pkin(3) - g(1) * t583 - g(2) * t584 - g(3) * t582 - t366 * qJ(4) - t412 * t424 - t414 * t681 - t425 * t473) * MDP(15) + (t524 * t590 + qJ(4) * t426 + t617 * t464 + t545 * t521 + (-t362 * t531 - t371 * t528 - t524 * t673) * qJD(2)) * MDP(16) + (-t521 * t590 + qJ(4) * t427 + t617 * t466 + t545 * t524 + (t363 * t531 + t372 * t528 + t521 * t673) * qJD(2)) * MDP(17) + (t371 * t466 + t372 * t464 + (qJD(5) * t466 - t427 * t526 - t572) * t524 + (qJD(5) * t464 - t426 * t526 - t571) * t521 + t549) * MDP(18) + (t360 * qJ(4) - t363 * t372 - t362 * t371 - g(1) * (-qJ(5) * t403 + t583) - g(2) * (-qJ(5) * t401 + t584) - g(3) * (-qJ(5) * t457 + t582) + (t350 * t524 + t351 * t521) * t526 + t617 * t390 + (-t362 * t524 - t363 * t521) * qJD(5)) * MDP(19) + (-t354 * t562 - t563 * t638) * MDP(20) + (-t354 * t469 + t355 * t562 - t387 * t638 + t563 * t637) * MDP(21) + (t563 * t628 + t574) * MDP(22) + (t387 * t628 + t682) * MDP(23) - t504 * MDP(24) * t628 + ((-t476 * t527 + t477 * t530) * t468 + t505 * t355 + t352 * t469 - t347 * t628 + (t527 * t553 - t530 * t554) * t504 + t618 * t387 + t637 * t374 + t548 * t511) * MDP(25) + (-(t476 * t530 + t477 * t527) * t468 + t505 * t354 - t352 * t562 + t348 * t628 + (t527 * t554 + t530 * t553) * t504 - t618 * t563 + t638 * t374 + t548 * t512) * MDP(26); MDP(12) * t611 + (qJDD(3) + t608) * MDP(13) + (-t533 - t653) * MDP(14) + (-t537 - t661 + t676) * MDP(15) - t549 * MDP(19) + t574 * MDP(25) + t682 * MDP(26) + (t414 * MDP(15) - t464 * MDP(16) - t466 * MDP(17) - t390 * MDP(19) - t387 * MDP(25) + MDP(26) * t563) * qJD(3) + (t551 * MDP(16) - MDP(17) * t653 + (-t464 * t629 - t427) * MDP(18) + t572 * MDP(19)) * t524 + (-MDP(16) * t653 - t551 * MDP(17) + (-t426 + t600) * MDP(18) + t571 * MDP(19)) * t521; (t426 + t600) * MDP(16) + ((-t464 + t626) * t629 + t579) * MDP(17) + (-t464 ^ 2 - t466 ^ 2) * MDP(18) + (t362 * t466 + t363 * t464 + t545) * MDP(19) + (t355 - t684) * MDP(25) + (t354 - t660) * MDP(26); -t563 * t387 * MDP(20) + (-t387 ^ 2 + t563 ^ 2) * MDP(21) + (t607 + t660) * MDP(22) + (-t580 - t684) * MDP(23) + t468 * MDP(24) + (t348 * t504 + t374 * t563 - g(1) * (t403 * t512 - t451 * t511) - g(2) * (t401 * t512 - t449 * t511) - g(3) * (t457 * t512 + t511 * t647) + t581) * MDP(25) + (t347 * t504 + t374 * t387 - g(1) * (-t403 * t511 - t451 * t512) - g(2) * (-t401 * t511 - t449 * t512) - g(3) * (-t457 * t511 + t512 * t647) - t568) * MDP(26) + (-MDP(22) * t657 + MDP(23) * t563 - MDP(25) * t348 - MDP(26) * t347) * qJD(6);];
tau  = t1;
