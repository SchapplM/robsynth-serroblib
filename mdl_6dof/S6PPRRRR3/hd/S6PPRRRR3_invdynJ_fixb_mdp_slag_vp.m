% Calculate vector of inverse dynamics joint torques for
% S6PPRRRR3
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PPRRRR3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPRRRR3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_invdynJ_fixb_mdp_slag_vp: pkin has to be [14x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6PPRRRR3_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:12:00
% EndTime: 2019-03-08 19:12:14
% DurationCPUTime: 10.63s
% Computational Cost: add. (6884->601), mult. (18295->885), div. (0->0), fcn. (17674->18), ass. (0->268)
t565 = sin(pkin(6));
t561 = sin(pkin(14));
t574 = sin(qJ(3));
t578 = cos(qJ(3));
t566 = cos(pkin(14));
t569 = cos(pkin(7));
t698 = t566 * t569;
t615 = -t561 * t574 + t578 * t698;
t734 = t565 * t615;
t570 = cos(pkin(6));
t551 = qJD(1) * t570 + qJD(2);
t564 = sin(pkin(7));
t700 = t564 * t578;
t533 = t551 * t700;
t477 = qJD(1) * t734 + t533;
t701 = t564 * t574;
t614 = t561 * t578 + t574 * t698;
t733 = t614 * t565;
t478 = -qJD(1) * t733 - t551 * t701;
t577 = cos(qJ(4));
t568 = cos(pkin(8));
t573 = sin(qJ(4));
t696 = t568 * t573;
t563 = sin(pkin(8));
t705 = t563 * t573;
t552 = pkin(10) * t705;
t695 = t568 * t577;
t730 = pkin(3) * t695 - t552;
t736 = -t730 * qJD(4) + t477 * t577 + t478 * t696;
t550 = t570 * qJDD(1) + qJDD(2);
t473 = qJD(3) * pkin(3) + t477;
t720 = t473 * t568;
t728 = qJD(1) * t615;
t735 = -pkin(10) * qJDD(3) * t563 - (qJD(3) * t551 * t578 + t550 * t574) * t564 - (qJD(3) * t728 + qJDD(1) * t614) * t565 - qJD(4) * t720;
t485 = t570 * t701 + t733;
t678 = qJD(3) * t577;
t654 = t563 * t678;
t547 = -qJD(5) + t654;
t576 = cos(qJ(5));
t572 = sin(qJ(5));
t679 = qJD(3) * t573;
t655 = t563 * t679;
t633 = t572 * t655;
t680 = qJD(3) * t568;
t639 = qJD(4) + t680;
t510 = -t576 * t639 + t633;
t504 = qJD(6) + t510;
t667 = qJDD(3) * t568;
t549 = qJDD(4) + t667;
t620 = qJD(5) * t639;
t676 = qJD(4) * t577;
t650 = t572 * t676;
t666 = qJDD(3) * t573;
t672 = qJD(5) * t576;
t454 = -t576 * t549 + t563 * (qJD(3) * (t573 * t672 + t650) + t572 * t666) + t572 * t620;
t618 = t563 * (pkin(4) * t573 - pkin(11) * t577);
t521 = qJD(4) * t618;
t703 = t563 * t577;
t683 = pkin(3) * t696 + pkin(10) * t703;
t517 = pkin(11) * t568 + t683;
t631 = -pkin(4) * t577 - pkin(11) * t573;
t518 = (-pkin(3) + t631) * t563;
t685 = t576 * t517 + t572 * t518;
t704 = t563 * t576;
t732 = qJD(5) * t685 - t478 * t704 - t576 * t521 - t572 * t736;
t674 = qJD(5) * t572;
t706 = t563 * t572;
t731 = t478 * t706 - t517 * t674 + t518 * t672 + t572 * t521 - t576 * t736;
t684 = t683 * qJD(4) - t477 * t573 + t478 * t695;
t567 = cos(pkin(13));
t562 = sin(pkin(13));
t697 = t567 * t570;
t617 = -t561 * t562 + t566 * t697;
t702 = t564 * t565;
t729 = -t567 * t702 + t617 * t569;
t681 = qJD(3) * t563;
t472 = pkin(10) * t681 - t478;
t660 = t566 * t702;
t507 = -qJD(1) * t660 + t551 * t569;
t414 = -t573 * t472 + t577 * (t507 * t563 + t720);
t665 = qJDD(3) * t577;
t548 = t563 * t665;
t668 = qJD(3) * qJD(4);
t648 = t573 * t668;
t519 = t563 * t648 + qJDD(5) - t548;
t601 = qJDD(1) * t734 + t550 * t700;
t439 = qJDD(3) * pkin(3) + qJD(3) * t478 + t601;
t503 = -qJDD(1) * t660 + t550 * t569;
t651 = t563 * t676;
t677 = qJD(4) * t573;
t592 = -t439 * t696 + t472 * t677 - t503 * t705 - t507 * t651 + t577 * t735;
t378 = pkin(11) * t549 - t592;
t415 = t577 * t472 + t473 * t696 + t507 * t705;
t410 = pkin(11) * t639 + t415;
t498 = t568 * t507;
t435 = t498 + (qJD(3) * t631 - t473) * t563;
t389 = t410 * t576 + t435 * t572;
t497 = t568 * t503;
t609 = t648 - t665;
t647 = t577 * t668;
t610 = t647 + t666;
t405 = t497 + (pkin(4) * t609 - pkin(11) * t610 - t439) * t563;
t587 = -qJD(5) * t389 - t378 * t572 + t576 * t405;
t370 = -pkin(5) * t519 - t587;
t512 = t572 * t639 + t576 * t655;
t524 = t561 * t697 + t562 * t566;
t455 = -t524 * t574 + t578 * t729;
t456 = t524 * t578 + t574 * t729;
t699 = t565 * t569;
t586 = -t564 * t617 - t567 * t699;
t581 = t586 * t563;
t396 = t456 * t577 + (t455 * t568 + t581) * t573;
t707 = t562 * t570;
t525 = -t561 * t707 + t566 * t567;
t616 = -t561 * t567 - t566 * t707;
t584 = t562 * t702 + t569 * t616;
t457 = -t525 * t574 + t578 * t584;
t458 = t525 * t578 + t574 * t584;
t585 = t562 * t699 - t564 * t616;
t580 = t585 * t563;
t398 = t458 * t577 + (t457 * t568 + t580) * t573;
t484 = t570 * t700 + t734;
t613 = t569 * t570 - t660;
t596 = t613 * t563;
t588 = t484 * t568 + t596;
t427 = t485 * t577 + t573 * t588;
t459 = -t484 * t563 + t568 * t613;
t399 = t427 * t572 - t459 * t576;
t428 = -t455 * t563 + t568 * t586;
t429 = -t457 * t563 + t568 * t585;
t606 = g(1) * (-t398 * t572 + t429 * t576) + g(2) * (-t396 * t572 + t428 * t576) - g(3) * t399;
t727 = t504 * (pkin(5) * t512 + pkin(12) * t504) + t370 + t606;
t450 = qJDD(6) + t454;
t546 = -pkin(5) * t576 - pkin(12) * t572 - pkin(4);
t726 = t504 * (t415 + t547 * (pkin(5) * t572 - pkin(12) * t576)) - t546 * t450;
t579 = qJD(3) ^ 2;
t724 = pkin(11) * qJD(5);
t632 = t576 * t651;
t646 = t563 * t666;
t453 = qJD(3) * t632 - qJD(5) * t633 + t572 * t549 + (t620 + t646) * t576;
t571 = sin(qJ(6));
t575 = cos(qJ(6));
t670 = qJD(6) * t575;
t656 = t575 * t453 + t571 * t519 - t547 * t670;
t671 = qJD(6) * t571;
t416 = -t512 * t671 + t656;
t723 = t416 * t571;
t722 = t450 * t571;
t721 = t450 * t575;
t713 = t512 * t571;
t474 = t575 * t547 + t713;
t719 = t474 * t504;
t476 = t512 * t575 - t547 * t571;
t718 = t476 * t504;
t483 = t485 * qJD(3);
t557 = t563 ^ 2;
t717 = t483 * t557;
t716 = t485 * t573;
t715 = t510 * t547;
t714 = t512 * t547;
t526 = -t563 * t700 + t568 * t569;
t712 = t526 * t563;
t710 = t547 * t572;
t694 = t573 * t574;
t693 = t573 * t578;
t692 = t574 * t577;
t691 = t574 * t579;
t690 = t576 * t577;
t689 = t577 * t578;
t520 = qJD(3) * t618;
t687 = t576 * t414 + t572 * t520;
t652 = t563 * t677;
t686 = -pkin(5) * t652 + t732;
t559 = t573 ^ 2;
t682 = -t577 ^ 2 + t559;
t675 = qJD(5) * t571;
t673 = qJD(5) * t575;
t669 = t549 * MDP(10);
t661 = t571 * t703;
t649 = t563 * t568 * t579;
t608 = t576 * t378 + t572 * t405 - t410 * t674 + t435 * t672;
t369 = pkin(12) * t519 + t608;
t619 = t439 * t695 - t472 * t676 + t503 * t703 - t507 * t652 + t573 * t735;
t379 = -pkin(4) * t549 - t619;
t372 = pkin(5) * t454 - pkin(12) * t453 + t379;
t643 = -t571 * t369 + t575 * t372;
t641 = t453 * t571 - t575 * t519;
t640 = t504 * t575;
t638 = qJD(4) + 0.2e1 * t680;
t637 = t549 + t667;
t636 = t557 * t564 * t691;
t634 = t681 * t701;
t501 = (t571 * t573 + t575 * t690) * t681;
t629 = t575 * t672 - t501;
t516 = t552 + (-pkin(3) * t577 - pkin(4)) * t568;
t527 = -t576 * t568 + t572 * t705;
t528 = t568 * t572 + t573 * t704;
t462 = pkin(5) * t527 - pkin(12) * t528 + t516;
t628 = -pkin(12) * t652 - qJD(6) * t462 - t731;
t464 = -pkin(12) * t703 + t685;
t492 = -qJD(5) * t527 + t632;
t493 = qJD(5) * t528 + t563 * t650;
t627 = -pkin(5) * t493 + pkin(12) * t492 + qJD(6) * t464 - t684;
t626 = t575 * t369 + t571 * t372;
t387 = -pkin(12) * t547 + t389;
t409 = -pkin(4) * t639 - t414;
t394 = t510 * pkin(5) - t512 * pkin(12) + t409;
t374 = t387 * t575 + t394 * t571;
t625 = t387 * t571 - t394 * t575;
t400 = t427 * t576 + t459 * t572;
t426 = -t484 * t695 - t577 * t596 + t716;
t381 = t400 * t575 + t426 * t571;
t380 = -t400 * t571 + t426 * t575;
t388 = -t410 * t572 + t435 * t576;
t611 = t568 * t693 + t692;
t487 = t564 * t611 + t569 * t705;
t461 = t487 * t576 + t526 * t572;
t612 = t568 * t689 - t694;
t486 = -t564 * t612 - t569 * t703;
t624 = t461 * t575 + t486 * t571;
t623 = -t461 * t571 + t486 * t575;
t460 = t487 * t572 - t526 * t576;
t621 = -t517 * t572 + t518 * t576;
t494 = t528 * t571 + t575 * t703;
t395 = -t455 * t695 + t456 * t573 - t577 * t581;
t397 = -t457 * t695 + t458 * t573 - t577 * t580;
t605 = g(1) * t397 + g(2) * t395 + g(3) * t426;
t604 = -g(1) * t398 - g(2) * t396 - g(3) * t427;
t419 = t455 * t577 - t456 * t696;
t421 = t457 * t577 - t458 * t696;
t444 = t484 * t577 - t485 * t696;
t603 = g(1) * t421 + g(2) * t419 + g(3) * t444;
t602 = g(1) * t458 + g(2) * t456 + g(3) * t485;
t594 = -t379 + t605;
t386 = pkin(5) * t547 - t388;
t591 = -pkin(12) * t450 + (t386 + t388) * t504;
t590 = -pkin(11) * t519 - t409 * t547;
t589 = pkin(11) * qJD(6) * t504 - t605;
t582 = (pkin(12) * t655 - qJD(6) * t546 + t687) * t504 + t604;
t500 = t571 * t576 * t654 - t575 * t655;
t495 = t528 * t575 - t661;
t482 = t484 * qJD(3);
t463 = pkin(5) * t703 - t621;
t452 = t569 * t652 + (t611 * qJD(4) + (t568 * t692 + t693) * qJD(3)) * t564;
t451 = t569 * t651 + (t612 * qJD(4) + (-t568 * t694 + t689) * qJD(3)) * t564;
t445 = -t473 * t563 + t498;
t443 = t484 * t573 + t485 * t695;
t442 = -qJD(6) * t661 + t492 * t571 + t528 * t670 - t575 * t652;
t441 = -qJD(6) * t494 + t492 * t575 + t571 * t652;
t423 = t444 * t576 + t485 * t706;
t422 = -t439 * t563 + t497;
t420 = t457 * t573 + t458 * t695;
t418 = t455 * t573 + t456 * t695;
t417 = qJD(6) * t476 + t641;
t412 = qJD(5) * t461 + t451 * t572 - t576 * t634;
t411 = -qJD(5) * t460 + t451 * t576 + t572 * t634;
t401 = -pkin(5) * t655 + t414 * t572 - t520 * t576;
t393 = -t483 * t696 + t482 * t577 + (t577 * t588 - t716) * qJD(4);
t392 = qJD(4) * t427 + t482 * t573 + t483 * t695;
t391 = t421 * t576 + t458 * t706;
t390 = t419 * t576 + t456 * t706;
t385 = t398 * t576 + t429 * t572;
t383 = t396 * t576 + t428 * t572;
t376 = -qJD(5) * t399 + t393 * t576 + t483 * t706;
t375 = qJD(5) * t400 + t393 * t572 - t483 * t704;
t368 = -t374 * qJD(6) + t643;
t367 = -qJD(6) * t625 + t626;
t1 = [(qJDD(1) - g(3)) * MDP(1) + (t550 * t570 - g(3) + (t561 ^ 2 + t566 ^ 2) * t565 ^ 2 * qJDD(1)) * MDP(2) + (-qJD(3) * t483 + qJDD(3) * t484) * MDP(4) + (-qJD(3) * t482 - qJDD(3) * t485) * MDP(5) + (-t459 * t548 - t392 * qJD(4) - t426 * t549 + (-t392 * t568 + t459 * t652 - t577 * t717) * qJD(3)) * MDP(11) + (t459 * t646 - t393 * qJD(4) - t427 * t549 + (-t393 * t568 + t459 * t651 + t573 * t717) * qJD(3)) * MDP(12) + (t375 * t547 + t392 * t510 - t399 * t519 + t426 * t454) * MDP(18) + (t376 * t547 + t392 * t512 - t400 * t519 + t426 * t453) * MDP(19) + ((-qJD(6) * t381 - t376 * t571 + t392 * t575) * t504 + t380 * t450 + t375 * t474 + t399 * t417) * MDP(25) + (-(qJD(6) * t380 + t376 * t575 + t392 * t571) * t504 - t381 * t450 + t375 * t476 + t399 * t416) * MDP(26); (-g(3) * t570 + (-g(1) * t562 + g(2) * t567) * t565 + t550) * MDP(2) + (-t452 * t639 - t486 * t549 - t577 * t636 + t609 * t712) * MDP(11) + (-t451 * t639 - t487 * t549 + t573 * t636 + t610 * t712) * MDP(12) + (t412 * t547 + t452 * t510 + t454 * t486 - t460 * t519) * MDP(18) + (t411 * t547 + t452 * t512 + t453 * t486 - t461 * t519) * MDP(19) + ((-qJD(6) * t624 - t411 * t571 + t452 * t575) * t504 + t623 * t450 + t412 * t474 + t460 * t417) * MDP(25) + (-(qJD(6) * t623 + t411 * t575 + t452 * t571) * t504 - t624 * t450 + t412 * t476 + t460 * t416) * MDP(26) + ((qJDD(3) * t578 - t691) * MDP(4) + (-qJDD(3) * t574 - t578 * t579) * MDP(5)) * t564; (g(1) * t420 + g(2) * t418 + g(3) * t443 - t683 * t549 + t592 * t568 + t736 * t639) * MDP(12) + (t453 * t528 + t492 * t512) * MDP(13) + (-t453 * t527 - t454 * t528 - t492 * t510 - t493 * t512) * MDP(14) + (-t492 * t547 + t519 * t528) * MDP(15) + (t493 * t547 - t519 * t527) * MDP(16) + (-g(1) * t391 - g(2) * t390 - g(3) * t423 + t379 * t527 + t409 * t493 + t516 * t454 + t684 * t510 + t621 * t519 + t732 * t547) * MDP(18) + (t379 * t528 + t409 * t492 + t516 * t453 + t684 * t512 - t685 * t519 + t731 * t547 + t603 * t572) * MDP(19) + (t416 * t495 + t441 * t476) * MDP(20) + (-t416 * t494 - t417 * t495 - t441 * t474 - t442 * t476) * MDP(21) + (t416 * t527 + t441 * t504 + t450 * t495 + t476 * t493) * MDP(22) + (-t417 * t527 - t442 * t504 - t450 * t494 - t474 * t493) * MDP(23) + (t450 * t527 + t493 * t504) * MDP(24) + ((t462 * t575 - t464 * t571) * t450 + t368 * t527 - t625 * t493 + t463 * t417 + t370 * t494 + t386 * t442 - g(1) * (t391 * t575 + t420 * t571) - g(2) * (t390 * t575 + t418 * t571) - g(3) * (t423 * t575 + t443 * t571) + (t571 * t628 - t575 * t627) * t504 + t686 * t474) * MDP(25) + (-(t462 * t571 + t464 * t575) * t450 - t367 * t527 - t374 * t493 + t463 * t416 + t370 * t495 + t386 * t441 - g(1) * (-t391 * t571 + t420 * t575) - g(2) * (-t390 * t571 + t418 * t575) - g(3) * (-t423 * t571 + t443 * t575) + (t571 * t627 + t575 * t628) * t504 + t686 * t476) * MDP(26) + qJDD(3) * MDP(3) + (-g(1) * t457 - g(2) * t455 - g(3) * t484 + t601) * MDP(4) + (-t550 * t701 - qJDD(1) * t733 + (-t565 * t728 + t477 - t533) * qJD(3) + t602) * MDP(5) + (t730 * t549 + t619 * t568 - t684 * t639 - t603) * MDP(11) + t568 * t669 + ((-pkin(3) * t610 + t478 * t679) * MDP(12) + (qJDD(3) * t559 + 0.2e1 * t573 * t647) * MDP(6) + 0.2e1 * (t573 * t665 - t668 * t682) * MDP(7) + (-pkin(3) * t609 - t478 * t678) * MDP(11)) * t557 + ((t422 * t573 + t445 * t676) * MDP(12) + (-t453 * t577 + t512 * t677) * MDP(15) + (t454 * t577 - t510 * t677) * MDP(16) + (-t519 * t577 - t547 * t677) * MDP(17) + (t388 * t677 - t577 * t587) * MDP(18) + (-t389 * t677 - t576 * t602 + t577 * t608) * MDP(19) + (t573 * t637 + t638 * t676) * MDP(8) + (t577 * t637 - t638 * t677) * MDP(9) + (-t422 * t577 + t445 * t677) * MDP(11)) * t563; (-t577 * t649 + t646) * MDP(8) + (t573 * t649 + t548) * MDP(9) + t669 + (t415 * t639 - t445 * t655 + t605 + t619) * MDP(11) + (t414 * t639 - t445 * t654 + t592 - t604) * MDP(12) + (t453 * t572 - t576 * t714) * MDP(13) + ((t453 + t715) * t576 + (-t454 + t714) * t572) * MDP(14) + (-t547 * t672 + t519 * t572 + (-t512 * t573 + t547 * t690) * t681) * MDP(15) + (t547 * t674 + t519 * t576 + (t510 * t573 - t577 * t710) * t681) * MDP(16) + t547 * MDP(17) * t655 + (-t388 * t655 - pkin(4) * t454 - t415 * t510 + (-t414 * t547 + t590) * t572 + ((t520 + t724) * t547 + t594) * t576) * MDP(18) + (-pkin(4) * t453 - t687 * t547 + t389 * t655 - t415 * t512 + t590 * t576 + (-t547 * t724 - t594) * t572) * MDP(19) + (t416 * t572 * t575 + (-t572 * t671 + t629) * t476) * MDP(20) + (t474 * t501 + t476 * t500 + (-t474 * t575 - t476 * t571) * t672 + (-t723 - t417 * t575 + (t474 * t571 - t476 * t575) * qJD(6)) * t572) * MDP(21) + (-t416 * t576 + t629 * t504 + (-t476 * t547 - t504 * t671 + t721) * t572) * MDP(22) + (t417 * t576 + (-t571 * t672 + t500) * t504 + (t474 * t547 - t504 * t670 - t722) * t572) * MDP(23) + (-t450 * t576 - t504 * t710) * MDP(24) + (-t386 * t500 - t401 * t474 - t726 * t575 + t582 * t571 + (t386 * t675 - t368 + (qJD(5) * t474 - t722) * pkin(11) - t589 * t575) * t576 + (t386 * t670 + t370 * t571 + t547 * t625 + (t504 * t675 + t417) * pkin(11)) * t572) * MDP(25) + (-t386 * t501 - t401 * t476 + t726 * t571 + t582 * t575 + (t386 * t673 + t367 + (qJD(5) * t476 - t721) * pkin(11) + t589 * t571) * t576 + (-t386 * t671 + t370 * t575 + t547 * t374 + (t504 * t673 + t416) * pkin(11)) * t572) * MDP(26) + (-MDP(6) * t573 * t577 + MDP(7) * t682) * t557 * t579; -t510 ^ 2 * MDP(14) + (t453 - t715) * MDP(15) + (-t454 - t714) * MDP(16) + t519 * MDP(17) + (-t389 * t547 + t587 - t606) * MDP(18) + (g(1) * t385 + g(2) * t383 + g(3) * t400 - t388 * t547 + t409 * t510 - t608) * MDP(19) + (t476 * t640 + t723) * MDP(20) + ((t416 - t719) * t575 + (-t417 - t718) * t571) * MDP(21) + (t504 * t640 + t722) * MDP(22) + (-t504 ^ 2 * t571 + t721) * MDP(23) + (-pkin(5) * t417 - t389 * t474 + t591 * t571 - t575 * t727) * MDP(25) + (-pkin(5) * t416 - t389 * t476 + t571 * t727 + t591 * t575) * MDP(26) + (MDP(13) * t510 + MDP(14) * t512 - MDP(18) * t409 - MDP(22) * t476 + MDP(23) * t474 - MDP(24) * t504 + MDP(25) * t625 + MDP(26) * t374) * t512; t476 * t474 * MDP(20) + (-t474 ^ 2 + t476 ^ 2) * MDP(21) + (t656 + t719) * MDP(22) + (-t641 + t718) * MDP(23) + t450 * MDP(24) + (t374 * t504 - t386 * t476 - g(1) * (-t385 * t571 + t397 * t575) - g(2) * (-t383 * t571 + t395 * t575) - g(3) * t380 + t643) * MDP(25) + (-t625 * t504 + t386 * t474 - g(1) * (-t385 * t575 - t397 * t571) - g(2) * (-t383 * t575 - t395 * t571) + g(3) * t381 - t626) * MDP(26) + (-MDP(22) * t713 - MDP(23) * t476 - MDP(25) * t374 + MDP(26) * t625) * qJD(6);];
tau  = t1;
