% Calculate vector of inverse dynamics joint torques for
% S6RPRRRP3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRRRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:05:32
% EndTime: 2019-03-09 06:05:44
% DurationCPUTime: 8.24s
% Computational Cost: add. (7011->590), mult. (14291->734), div. (0->0), fcn. (9583->14), ass. (0->248)
t571 = sin(qJ(4));
t574 = cos(qJ(4));
t670 = qJD(3) * t574;
t572 = sin(qJ(3));
t676 = qJD(1) * t572;
t508 = -t571 * t676 + t670;
t672 = qJD(3) * t571;
t509 = t574 * t676 + t672;
t570 = sin(qJ(5));
t721 = cos(qJ(5));
t450 = -t721 * t508 + t509 * t570;
t608 = t570 * t508 + t509 * t721;
t750 = t450 * t608;
t568 = sin(pkin(10));
t550 = pkin(1) * t568 + pkin(7);
t528 = t550 * qJD(1);
t575 = cos(qJ(3));
t480 = qJD(2) * t575 - t572 * t528;
t626 = pkin(3) * t572 - pkin(8) * t575;
t513 = t626 * qJD(1);
t496 = t574 * t513;
t696 = t574 * t575;
t617 = pkin(4) * t572 - pkin(9) * t696;
t722 = pkin(9) + pkin(8);
t650 = qJD(4) * t722;
t749 = qJD(1) * t617 - t480 * t571 + t574 * t650 + t496;
t675 = qJD(1) * t575;
t648 = t571 * t675;
t683 = t574 * t480 + t571 * t513;
t748 = pkin(9) * t648 - t571 * t650 - t683;
t634 = qJD(4) + t675;
t657 = qJDD(1) * t572;
t595 = qJD(3) * t634 + t657;
t660 = qJD(1) * qJD(4);
t639 = t572 * t660;
t619 = -qJDD(3) + t639;
t604 = t619 * t571;
t583 = t595 * t574 - t604;
t616 = qJD(3) * qJD(4) + t657;
t661 = qJD(1) * qJD(3);
t640 = t575 * t661;
t593 = t616 + t640;
t631 = t571 * t593 + t574 * t639;
t600 = qJDD(3) * t574 - t631;
t641 = t721 * qJD(5);
t664 = qJD(5) * t570;
t392 = -t508 * t641 + t509 * t664 - t570 * t600 - t721 * t583;
t548 = -qJD(4) + t675;
t538 = -qJD(5) + t548;
t384 = -t450 * t538 - t392;
t393 = qJD(5) * t608 + t570 * t583 - t721 * t600;
t560 = t575 * qJDD(1);
t504 = t572 * t661 + qJDD(4) - t560;
t498 = qJDD(5) + t504;
t723 = t608 ^ 2;
t747 = t498 * MDP(23) + (-t538 * t608 - t393) * MDP(22) + MDP(19) * t750 + (-t450 ^ 2 + t723) * MDP(20) + t384 * MDP(21);
t531 = t722 * t571;
t532 = t722 * t574;
t460 = -t570 * t531 + t532 * t721;
t567 = qJ(4) + qJ(5);
t562 = sin(t567);
t564 = qJ(1) + pkin(10);
t558 = sin(t564);
t559 = cos(t564);
t623 = g(1) * t559 + g(2) * t558;
t613 = t623 * t572;
t714 = g(3) * t575;
t746 = t460 * t498 + (t613 - t714) * t562;
t470 = -qJD(3) * pkin(3) - t480;
t444 = -pkin(4) * t508 + t470;
t396 = pkin(5) * t450 - qJ(6) * t608 + t444;
t745 = t396 * t450;
t744 = t444 * t450;
t701 = t570 * t571;
t606 = t574 * t721 - t701;
t728 = qJD(4) + qJD(5);
t730 = t721 * qJD(4) + t641;
t685 = -t574 * t730 + t606 * t675 + t701 * t728;
t511 = t570 * t574 + t571 * t721;
t456 = t728 * t511;
t684 = -t511 * t675 + t456;
t523 = t550 * qJDD(1);
t743 = qJD(2) * qJD(3) + t523;
t669 = qJD(3) * t575;
t646 = t571 * t669;
t665 = qJD(4) * t574;
t740 = t572 * t665 + t646;
t738 = qJD(3) * t480;
t411 = pkin(5) * t608 + qJ(6) * t450;
t630 = t721 * t669;
t418 = t456 * t572 + t570 * t646 - t574 * t630;
t483 = t606 * t572;
t691 = t418 * t538 + t483 * t498;
t666 = qJD(4) * t572;
t645 = t571 * t666;
t700 = t571 * t572;
t419 = t571 * t630 - t570 * t645 - t664 * t700 + (t570 * t669 + t572 * t730) * t574;
t482 = t511 * t572;
t690 = t419 * t538 - t482 * t498;
t607 = -t531 * t721 - t570 * t532;
t734 = qJD(5) * t607 - t570 * t749 + t721 * t748;
t733 = qJD(5) * t460 + t570 * t748 + t721 * t749;
t481 = t572 * qJD(2) + t575 * t528;
t471 = qJD(3) * pkin(8) + t481;
t569 = cos(pkin(10));
t552 = -pkin(1) * t569 - pkin(2);
t499 = -pkin(3) * t575 - pkin(8) * t572 + t552;
t472 = t499 * qJD(1);
t426 = t574 * t471 + t571 * t472;
t417 = pkin(9) * t508 + t426;
t485 = t574 * t499;
t698 = t572 * t574;
t707 = t550 * t571;
t434 = -pkin(9) * t698 + t485 + (-pkin(4) - t707) * t575;
t512 = t550 * t696;
t680 = t571 * t499 + t512;
t443 = -pkin(9) * t700 + t680;
t732 = t570 * t434 + t721 * t443;
t671 = qJD(3) * t572;
t638 = t392 * t575 + t608 * t671;
t637 = -t393 * t575 + t450 * t671;
t667 = qJD(4) * t571;
t628 = -t481 + (-t648 + t667) * pkin(4);
t490 = t498 * qJ(6);
t520 = t538 * qJD(6);
t731 = t490 - t520;
t699 = t571 * t575;
t475 = t558 * t699 + t559 * t574;
t477 = t558 * t574 - t559 * t699;
t729 = -g(1) * t477 + g(2) * t475;
t491 = t498 * pkin(5);
t727 = t491 - qJDD(6);
t563 = cos(t567);
t705 = t562 * t575;
t465 = t558 * t705 + t559 * t563;
t467 = -t558 * t563 + t559 * t705;
t432 = qJDD(3) * pkin(8) + qJDD(2) * t572 + t523 * t575 + t738;
t516 = t626 * qJD(3);
t445 = qJD(1) * t516 + qJDD(1) * t499;
t440 = t574 * t445;
t655 = t571 * qJDD(3);
t656 = qJDD(1) * t574;
t379 = -t571 * t432 + t440 - (t572 * t656 + t574 * t640 + t655) * pkin(9) + t504 * pkin(4) - t417 * qJD(4);
t653 = t574 * t432 + t571 * t445 + t472 * t665;
t602 = -t471 * t667 + t653;
t382 = pkin(9) * t600 + t602;
t425 = -t471 * t571 + t574 * t472;
t416 = -pkin(9) * t509 + t425;
t409 = -pkin(4) * t548 + t416;
t632 = -t721 * t379 + t570 * t382 + t409 * t664 + t417 * t641;
t706 = t562 * t572;
t591 = g(1) * t467 + g(2) * t465 + g(3) * t706 - t632;
t584 = t396 * t608 - t591 - t727;
t726 = -t444 * t608 + t591;
t681 = t574 * t516 + t671 * t707;
t403 = t617 * qJD(3) + (-t512 + (pkin(9) * t572 - t499) * t571) * qJD(4) + t681;
t682 = t499 * t665 + t571 * t516;
t406 = (-t572 * t670 - t575 * t667) * t550 - t740 * pkin(9) + t682;
t725 = -qJD(5) * t732 + t403 * t721 - t570 * t406;
t715 = g(3) * t572;
t713 = qJDD(3) * pkin(3);
t649 = t721 * t417;
t388 = t570 * t409 + t649;
t712 = t388 * t538;
t708 = t509 * t548;
t704 = t563 * t572;
t703 = t563 * t575;
t702 = t570 * t417;
t697 = t572 * t722;
t695 = t575 * t548;
t694 = qJDD(2) - g(3);
t693 = -t483 * t393 + t418 * t450;
t692 = pkin(5) * t684 + qJ(6) * t685 - qJD(6) * t511 + t628;
t689 = -qJ(6) * t676 + t734;
t688 = pkin(5) * t676 + t733;
t549 = pkin(4) * t700;
t488 = t572 * t550 + t549;
t390 = t416 * t721 - t702;
t679 = pkin(4) * t641 + qJD(6) - t390;
t565 = t572 ^ 2;
t678 = -t575 ^ 2 + t565;
t529 = qJD(1) * t552;
t673 = qJD(3) * t508;
t668 = qJD(4) * t508;
t663 = t470 * qJD(4);
t387 = t409 * t721 - t702;
t662 = qJD(6) - t387;
t652 = t528 * t669 + t572 * t743;
t457 = pkin(4) * t740 + t550 * t669;
t557 = pkin(4) * t574 + pkin(3);
t651 = pkin(4) * t571 + pkin(7);
t647 = t548 * t672;
t644 = t548 * t666;
t636 = t548 * t550 + t471;
t635 = -qJD(4) * t472 - t432;
t633 = t570 * t379 + t721 * t382 + t409 * t641 - t417 * t664;
t629 = -pkin(5) * t706 + qJ(6) * t704;
t389 = t570 * t416 + t649;
t627 = pkin(4) * t664 - t389;
t625 = -g(1) * t465 + g(2) * t467;
t466 = t558 * t703 - t559 * t562;
t468 = t558 * t562 + t559 * t703;
t624 = g(1) * t466 - g(2) * t468;
t622 = g(1) * t558 - g(2) * t559;
t573 = sin(qJ(1));
t576 = cos(qJ(1));
t621 = g(1) * t573 - g(2) * t576;
t620 = -t548 + t634;
t618 = -t392 * t482 + t419 * t608;
t615 = t557 * t575 + pkin(2) + t697;
t614 = pkin(5) * t563 + qJ(6) * t562 + t557;
t611 = t434 * t721 - t570 * t443;
t605 = -t504 * t571 + t548 * t665;
t603 = t619 * t575;
t601 = t570 * t403 + t721 * t406 + t434 * t641 - t443 * t664;
t599 = -qJD(1) * t529 + t623;
t598 = -g(3) * t703 + t607 * t498 + t623 * t704;
t597 = -pkin(8) * t504 - t470 * t548;
t433 = -t575 * qJDD(2) + t652 - t713;
t596 = 0.2e1 * qJD(3) * t529 - qJDD(3) * t550;
t594 = -qJD(4) * pkin(8) * t548 + t433 - t713 + t714;
t592 = g(1) * t468 + g(2) * t466 + g(3) * t704 - t633;
t588 = -g(1) * (-t467 * pkin(5) + qJ(6) * t468) - g(2) * (-t465 * pkin(5) + qJ(6) * t466);
t578 = qJD(3) ^ 2;
t587 = -0.2e1 * qJDD(1) * t552 - t550 * t578 + t622;
t585 = -t387 * t538 + t592;
t405 = -pkin(4) * t600 + t433;
t556 = -pkin(4) * t721 - pkin(5);
t551 = pkin(4) * t570 + qJ(6);
t519 = qJDD(3) * t575 - t572 * t578;
t518 = qJDD(3) * t572 + t575 * t578;
t486 = t509 * t671;
t478 = t558 * t571 + t559 * t696;
t476 = -t558 * t696 + t559 * t571;
t446 = -pkin(5) * t606 - qJ(6) * t511 - t557;
t422 = pkin(5) * t482 - qJ(6) * t483 + t488;
t402 = pkin(4) * t509 + t411;
t399 = t575 * pkin(5) - t611;
t398 = -qJ(6) * t575 + t732;
t386 = -t538 * qJ(6) + t388;
t385 = t538 * pkin(5) + t662;
t383 = pkin(5) * t419 + qJ(6) * t418 - qJD(6) * t483 + t457;
t376 = -pkin(5) * t671 - t725;
t375 = qJ(6) * t671 - qJD(6) * t575 + t601;
t374 = t393 * pkin(5) + t392 * qJ(6) - qJD(6) * t608 + t405;
t373 = t632 - t727;
t372 = t633 + t731;
t1 = [(t387 * t671 + t488 * t393 + t405 * t482 + t444 * t419 + t457 * t450 + t611 * t498 - t538 * t725 + t632 * t575 + t624) * MDP(24) + (t638 + t691) * MDP(21) + (-t618 + t693) * MDP(20) + t518 * MDP(7) + t519 * MDP(8) + (t690 - t637) * MDP(22) + qJDD(1) * MDP(1) + (t682 * t548 - t680 * t504 - g(1) * t475 - g(2) * t477 + (-t636 * t667 + (t470 * t574 + t509 * t550) * qJD(3) + t653) * t575 + (-t571 * t663 + t433 * t574 + (t655 + (-t571 * t660 + t656) * t572) * t550 + (t550 * t574 * t620 - t426) * qJD(3)) * t572) * MDP(18) + (-(-t499 * t667 + t681) * t548 + t485 * t504 - g(1) * t476 - g(2) * t478 + (-t550 * t673 - t440 + t636 * t665 + (qJD(3) * t470 - t504 * t550 - t635) * t571) * t575 + (t425 * qJD(3) + t433 * t571 - t550 * t600 + t574 * t663) * t572) * MDP(17) + (-t388 * t671 - t488 * t392 + t405 * t483 - t444 * t418 + t457 * t608 - t498 * t732 + t538 * t601 + t575 * t633 + t625) * MDP(25) + ((t508 * t574 - t509 * t571) * t669 + ((-t509 * qJD(4) + t600) * t574 + (-t668 - t583) * t571) * t572) * MDP(13) + 0.2e1 * (t560 * t572 - t661 * t678) * MDP(6) + ((-t600 + t647) * t575 + (t605 + t673) * t572) * MDP(15) + (-t504 * t575 - t548 * t671) * MDP(16) + (-t498 * t575 - t538 * t671) * MDP(23) + (t373 * t575 + t374 * t482 + t376 * t538 + t383 * t450 - t385 * t671 + t393 * t422 + t396 * t419 - t399 * t498 + t624) * MDP(26) + (-t392 * t483 - t418 * t608) * MDP(19) + (-t372 * t575 - t374 * t483 - t375 * t538 - t383 * t608 + t386 * t671 + t392 * t422 + t396 * t418 + t398 * t498 - t625) * MDP(28) + (-t372 * t482 + t373 * t483 - t375 * t450 + t376 * t608 - t385 * t418 - t386 * t419 - t392 * t399 - t393 * t398 + t572 * t622) * MDP(27) + (t372 * t398 + t386 * t375 + t374 * t422 + t396 * t383 + t373 * t399 + t385 * t376 - g(1) * (-pkin(1) * t573 - pkin(5) * t466 - qJ(6) * t465) - g(2) * (pkin(1) * t576 + pkin(5) * t468 + qJ(6) * t467) + (-g(1) * t651 - g(2) * t615) * t559 + (g(1) * t615 - g(2) * t651) * t558) * MDP(29) + (qJDD(1) * t565 + 0.2e1 * t572 * t640) * MDP(5) + (t486 + (t603 + t644) * t571 + ((t504 - t560) * t572 + (-t548 - t634) * t669) * t574) * MDP(14) + (-t509 * t645 + (t509 * t669 - t572 * t604 + t593 * t698) * t574) * MDP(12) + t621 * MDP(2) + (t621 + (t568 ^ 2 + t569 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + (g(1) * t576 + g(2) * t573) * MDP(3) + (t572 * t596 + t575 * t587) * MDP(10) + (-t572 * t587 + t575 * t596) * MDP(11); t694 * MDP(4) + t519 * MDP(10) - t518 * MDP(11) + ((t600 + t647) * t575 + (t605 - t673) * t572) * MDP(17) + (t486 + (t603 - t644) * t571 + ((-t504 - t560) * t572 - t620 * t669) * t574) * MDP(18) + (t618 + t693) * MDP(27) + (t372 * t483 + t373 * t482 - t374 * t575 + t385 * t419 - t386 * t418 + t396 * t671 - g(3)) * MDP(29) + (MDP(24) + MDP(26)) * (t637 + t690) + (-MDP(25) + MDP(28)) * (t691 - t638); (-t619 * t571 ^ 2 + (t571 * t595 - t708) * t574) * MDP(12) + ((-t631 + t708) * t571 + (t668 + 0.2e1 * t655 + t616 * t574 + (-t645 + (-t508 + t670) * t575) * qJD(1)) * t574) * MDP(13) + (t738 + (qJD(3) * t528 - t694) * t572 + (t599 - t743) * t575) * MDP(11) + (qJD(3) * t481 + t572 * t599 + t575 * t694 - t652) * MDP(10) + (t548 * t667 + t574 * t504 + (-t508 * t572 - t571 * t695) * qJD(1)) * MDP(15) + ((-t509 * t572 + t574 * t695) * qJD(1) - t605) * MDP(14) + (-t374 * t606 + t393 * t446 + t396 * t684 + t450 * t692 + t538 * t688 + t598) * MDP(26) + qJDD(3) * MDP(9) + (-t392 * t606 - t393 * t511 + t450 * t685 - t608 * t684) * MDP(20) + (t498 * t511 + t538 * t685) * MDP(21) + (-t392 * t511 - t608 * t685) * MDP(19) + (-t683 * t548 - t481 * t509 + (-pkin(3) * t593 + t597) * t574 + ((pkin(3) * t660 - t623) * t572 + t594) * t571) * MDP(18) + (-t557 * t393 - t405 * t606 + t684 * t444 + t628 * t450 + t538 * t733 + t598) * MDP(24) + (t498 * t606 + t538 * t684) * MDP(22) + (-g(3) * t697 + t372 * t460 - t373 * t607 + t374 * t446 + t385 * t688 + t386 * t689 + t396 * t692 - t614 * t714 + t623 * (t572 * t614 - t575 * t722)) * MDP(29) + (-t374 * t511 + t392 * t446 + t685 * t396 - t689 * t538 - t692 * t608 + t746) * MDP(28) + (t557 * t392 + t405 * t511 - t685 * t444 + t734 * t538 + t628 * t608 - t746) * MDP(25) + (t372 * t606 + t373 * t511 - t385 * t685 - t386 * t684 + t392 * t607 - t393 * t460 - t450 * t689 - t575 * t623 + t608 * t688 - t715) * MDP(27) + MDP(8) * t560 + MDP(7) * t657 + (-pkin(3) * t631 + t496 * t548 + t481 * t508 + (-t480 * t548 + t597) * t571 + (t613 - t594) * t574) * MDP(17) + (t548 * MDP(16) - t425 * MDP(17) + t426 * MDP(18) - MDP(21) * t608 + t450 * MDP(22) + t538 * MDP(23) - t387 * MDP(24) + t388 * MDP(25) + t385 * MDP(26) - t386 * MDP(28)) * t676 + (-MDP(5) * t572 * t575 + MDP(6) * t678) * qJD(1) ^ 2; -t509 * t508 * MDP(12) + (-t508 ^ 2 + t509 ^ 2) * MDP(13) + (t508 * t548 + t583) * MDP(14) + (t600 - t708) * MDP(15) + t504 * MDP(16) + (-t471 * t665 - t426 * t548 - t470 * t509 + t440 + (t635 + t715) * t571 + t729) * MDP(17) + (g(1) * t478 - g(2) * t476 + g(3) * t698 - t425 * t548 - t470 * t508 - t602) * MDP(18) + (-t389 * t538 + (-t450 * t509 + t498 * t721 + t538 * t664) * pkin(4) + t726) * MDP(24) + (-t390 * t538 + t744 + (-t498 * t570 - t509 * t608 + t538 * t641) * pkin(4) + t592) * MDP(25) + (-t402 * t450 - t498 * t556 + t538 * t627 - t584) * MDP(26) + (-t392 * t556 - t393 * t551 + (t386 + t627) * t608 + (t385 - t679) * t450) * MDP(27) + (t402 * t608 + t498 * t551 - t538 * t679 - t592 + t731 - t745) * MDP(28) + (t372 * t551 + t373 * t556 - t396 * t402 - t385 * t389 - g(3) * (-t549 + t629) + t679 * t386 + (t385 * t664 + t729) * pkin(4) + t588) * MDP(29) + t747; (-t712 + t726) * MDP(24) + (t585 + t744) * MDP(25) + (-t411 * t450 + t491 - t584 - t712) * MDP(26) + (pkin(5) * t392 - qJ(6) * t393 + (t386 - t388) * t608 + (t385 - t662) * t450) * MDP(27) + (t411 * t608 + 0.2e1 * t490 - 0.2e1 * t520 - t585 - t745) * MDP(28) + (-t373 * pkin(5) - g(3) * t629 + t372 * qJ(6) - t385 * t388 + t386 * t662 - t396 * t411 + t588) * MDP(29) + t747; (-t498 + t750) * MDP(26) + t384 * MDP(27) + (-t538 ^ 2 - t723) * MDP(28) + (t386 * t538 + t584) * MDP(29);];
tau  = t1;
