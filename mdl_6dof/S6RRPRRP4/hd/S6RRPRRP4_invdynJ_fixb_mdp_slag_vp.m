% Calculate vector of inverse dynamics joint torques for
% S6RRPRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:56:00
% EndTime: 2019-03-09 11:56:16
% DurationCPUTime: 12.81s
% Computational Cost: add. (11366->633), mult. (25945->768), div. (0->0), fcn. (19439->14), ass. (0->276)
t621 = qJ(2) + pkin(10);
t617 = cos(t621);
t767 = g(3) * t617;
t616 = sin(t621);
t630 = sin(qJ(1));
t633 = cos(qJ(1));
t677 = g(1) * t633 + g(2) * t630;
t808 = t616 * t677;
t652 = -t767 + t808;
t625 = sin(pkin(10));
t629 = sin(qJ(2));
t632 = cos(qJ(2));
t760 = cos(pkin(10));
t584 = t625 * t632 + t629 * t760;
t572 = t584 * qJD(1);
t628 = sin(qJ(4));
t631 = cos(qJ(4));
t540 = t631 * qJD(2) - t628 * t572;
t541 = qJD(2) * t628 + t572 * t631;
t627 = sin(qJ(5));
t775 = cos(qJ(5));
t477 = -t775 * t540 + t541 * t627;
t662 = t627 * t540 + t541 * t775;
t807 = t477 * t662;
t685 = t760 * t632;
t603 = qJD(1) * t685;
t716 = qJD(1) * t629;
t570 = -t625 * t716 + t603;
t500 = pkin(2) * t716 + pkin(3) * t572 - pkin(8) * t570;
t493 = t631 * t500;
t626 = -qJ(3) - pkin(7);
t596 = t626 * t632;
t589 = qJD(1) * t596;
t575 = t625 * t589;
t595 = t626 * t629;
t588 = qJD(1) * t595;
t531 = t588 * t760 + t575;
t774 = pkin(2) * t625;
t608 = pkin(8) + t774;
t762 = pkin(9) + t608;
t687 = qJD(4) * t762;
t748 = t570 * t631;
t806 = pkin(4) * t572 - pkin(9) * t748 - t531 * t628 + t631 * t687 + t493;
t721 = t628 * t500 + t631 * t531;
t749 = t570 * t628;
t805 = -pkin(9) * t749 + t628 * t687 + t721;
t715 = qJD(4) * t628;
t804 = t715 - t749;
t711 = qJD(1) * qJD(2);
t691 = t629 * t711;
t647 = qJDD(1) * t584 - t625 * t691;
t641 = qJD(2) * t603 + t647;
t639 = t628 * qJDD(2) + t631 * t641;
t638 = qJD(4) * t540 + t639;
t710 = qJD(2) * qJD(4);
t714 = qJD(4) * t631;
t702 = t572 * t714 + (t641 + t710) * t628;
t668 = qJDD(2) * t631 - t702;
t692 = t775 * qJD(5);
t713 = qJD(5) * t627;
t429 = -t540 * t692 + t541 * t713 - t627 * t668 - t775 * t638;
t560 = qJD(4) - t570;
t557 = qJD(5) + t560;
t416 = t477 * t557 - t429;
t430 = t540 * t713 + t541 * t692 + t627 * t638 - t775 * t668;
t571 = t584 * qJD(2);
t709 = qJDD(1) * t629;
t673 = qJDD(1) * t685 - t625 * t709;
t521 = qJD(1) * t571 + qJDD(4) - t673;
t516 = qJDD(5) + t521;
t776 = t662 ^ 2;
t803 = t516 * MDP(24) + (t557 * t662 - t430) * MDP(23) + MDP(20) * t807 + (-t477 ^ 2 + t776) * MDP(21) + t416 * MDP(22);
t761 = qJD(2) * pkin(2);
t578 = t588 + t761;
t526 = t578 * t760 + t575;
t513 = -qJD(2) * pkin(3) - t526;
t475 = -t540 * pkin(4) + t513;
t433 = t477 * pkin(5) - qJ(6) * t662 + t475;
t802 = t433 * t477;
t801 = t475 * t477;
t700 = t775 * t628;
t587 = t627 * t631 + t700;
t495 = t587 * t570;
t783 = qJD(4) + qJD(5);
t534 = t783 * t587;
t723 = t495 - t534;
t737 = t627 * t628;
t660 = t775 * t631 - t737;
t496 = t660 * t570;
t786 = t775 * qJD(4) + t692;
t533 = -t631 * t786 + t737 * t783;
t722 = t496 + t533;
t658 = -t625 * t629 + t685;
t574 = t658 * qJD(2);
t800 = t574 * t628 + t584 * t714;
t579 = t762 * t628;
t580 = t762 * t631;
t523 = -t627 * t579 + t580 * t775;
t624 = qJ(4) + qJ(5);
t618 = sin(t624);
t797 = t523 * t516 + t618 * t652;
t449 = pkin(5) * t662 + qJ(6) * t477;
t754 = t521 * t628;
t763 = t632 * pkin(2);
t615 = pkin(1) + t763;
t525 = -pkin(3) * t658 - pkin(8) * t584 - t615;
t512 = t631 * t525;
t539 = t625 * t595 - t596 * t760;
t744 = t584 * t631;
t453 = -pkin(4) * t658 - pkin(9) * t744 - t539 * t628 + t512;
t532 = t631 * t539;
t720 = t628 * t525 + t532;
t745 = t584 * t628;
t462 = -pkin(9) * t745 + t720;
t792 = t627 * t453 + t775 * t462;
t688 = qJD(2) * t626;
t569 = -qJD(3) * t629 + t632 * t688;
t520 = qJDD(2) * pkin(2) + qJD(1) * t569 + qJDD(1) * t595;
t568 = qJD(3) * t632 + t629 * t688;
t529 = qJD(1) * t568 - qJDD(1) * t596;
t470 = t520 * t760 - t625 * t529;
t466 = -qJDD(2) * pkin(3) - t470;
t791 = qJD(4) * t608 * t560 + t466;
t661 = -t579 * t775 - t627 * t580;
t790 = -qJD(5) * t661 + t627 * t806 + t775 * t805;
t789 = -qJD(5) * t523 + t627 * t805 - t775 * t806;
t686 = t760 * t589;
t530 = t588 * t625 - t686;
t681 = pkin(4) * t804 - t530;
t505 = t516 * qJ(6);
t549 = t557 * qJD(6);
t788 = t505 + t549;
t787 = g(1) * t630 - g(2) * t633;
t697 = t584 * t715;
t746 = t574 * t631;
t785 = -t697 + t746;
t732 = t631 * t633;
t736 = t628 * t630;
t561 = t617 * t736 + t732;
t733 = t630 * t631;
t735 = t628 * t633;
t563 = -t617 * t735 + t733;
t784 = -g(1) * t563 + g(2) * t561;
t507 = t516 * pkin(5);
t782 = t507 - qJDD(6);
t619 = cos(t624);
t739 = t619 * t633;
t740 = t618 * t630;
t545 = t617 * t740 + t739;
t731 = t633 * t618;
t734 = t630 * t619;
t547 = t617 * t731 - t734;
t471 = t625 * t520 + t760 * t529;
t467 = qJDD(2) * pkin(8) + t471;
t593 = -qJD(1) * t615 + qJD(3);
t489 = -pkin(3) * t570 - pkin(8) * t572 + t593;
t528 = -qJD(2) * t572 + t673;
t707 = pkin(2) * t691 + qJDD(3);
t708 = qJDD(1) * t632;
t759 = qJDD(1) * pkin(1);
t459 = -pkin(2) * t708 - t528 * pkin(3) - pkin(8) * t641 + t707 - t759;
t527 = t625 * t578 - t686;
t514 = qJD(2) * pkin(8) + t527;
t674 = t631 * t459 - t514 * t714;
t410 = t521 * pkin(4) - pkin(9) * t638 - t628 * t467 - t489 * t715 + t674;
t657 = t628 * t459 + t631 * t467 + t489 * t714 - t514 * t715;
t414 = pkin(9) * t668 + t657;
t457 = t631 * t489 - t514 * t628;
t447 = -pkin(9) * t541 + t457;
t439 = pkin(4) * t560 + t447;
t458 = t628 * t489 + t631 * t514;
t448 = pkin(9) * t540 + t458;
t682 = -t775 * t410 + t627 * t414 + t439 * t713 + t448 * t692;
t743 = t616 * t618;
t648 = g(1) * t547 + g(2) * t545 + g(3) * t743 - t682;
t643 = t433 * t662 - t648 - t782;
t781 = -t475 * t662 + t648;
t780 = -t560 ^ 2 * t631 - t754;
t779 = t429 * t660 - t662 * t723;
t704 = t629 * t761;
t501 = pkin(3) * t571 - pkin(8) * t574 + t704;
t494 = t631 * t501;
t499 = t568 * t760 + t625 * t569;
t426 = -pkin(9) * t746 + pkin(4) * t571 - t499 * t628 + t494 + (-t532 + (pkin(9) * t584 - t525) * t628) * qJD(4);
t656 = t631 * t499 + t628 * t501 + t525 * t714 - t539 * t715;
t432 = -pkin(9) * t800 + t656;
t778 = -qJD(5) * t792 + t426 * t775 - t627 * t432;
t768 = g(3) * t616;
t766 = g(3) * t628;
t765 = g(3) * t632;
t764 = t631 * pkin(4);
t701 = t775 * t448;
t419 = t627 * t439 + t701;
t758 = t419 * t557;
t756 = t477 * t572;
t755 = t662 * t572;
t752 = t540 * t572;
t751 = t541 * t560;
t750 = t541 * t572;
t483 = t660 * t516;
t484 = t587 * t516;
t742 = t616 * t619;
t634 = -pkin(9) - pkin(8);
t741 = t616 * t634;
t738 = t627 * t448;
t506 = t631 * t521;
t728 = -pkin(5) * t723 + qJ(6) * t722 - qJD(6) * t587 + t681;
t727 = -qJ(6) * t572 - t790;
t726 = t572 * pkin(5) - t789;
t725 = -t534 * t557 + t483;
t724 = -t533 * t557 + t484;
t423 = t447 * t775 - t738;
t719 = pkin(4) * t692 + qJD(6) - t423;
t622 = t629 ^ 2;
t718 = -t632 ^ 2 + t622;
t418 = t439 * t775 - t738;
t712 = qJD(6) - t418;
t698 = t760 * pkin(2);
t694 = t513 * t714;
t690 = pkin(4) * t628 - t626;
t498 = t568 * t625 - t760 * t569;
t538 = -t760 * t595 - t596 * t625;
t684 = -qJD(4) * t489 - t467;
t683 = t627 * t410 + t775 * t414 + t439 * t692 - t448 * t713;
t610 = -t698 - pkin(3);
t422 = t627 * t447 + t701;
t680 = pkin(4) * t713 - t422;
t679 = -g(1) * t545 + g(2) * t547;
t546 = t617 * t734 - t731;
t548 = t617 * t739 + t740;
t678 = g(1) * t546 - g(2) * t548;
t675 = -t587 * t430 + t477 * t722;
t497 = pkin(4) * t745 + t538;
t614 = pkin(3) + t764;
t672 = t614 * t617 - t741;
t465 = pkin(4) * t800 + t498;
t670 = -t560 * t804 + t506;
t669 = pkin(5) * t619 + qJ(6) * t618 + t614;
t667 = -0.2e1 * pkin(1) * t711 - pkin(7) * qJDD(2);
t665 = t453 * t775 - t627 * t462;
t594 = t610 - t764;
t655 = t627 * t426 + t775 * t432 + t453 * t692 - t462 * t713;
t654 = -qJDD(1) * t615 + t707;
t653 = t661 * t516 + t619 * t652;
t635 = qJD(2) ^ 2;
t651 = -pkin(7) * t635 + 0.2e1 * t759 + t787;
t636 = qJD(1) ^ 2;
t650 = pkin(1) * t636 - pkin(7) * qJDD(1) + t677;
t649 = g(1) * t548 + g(2) * t546 + g(3) * t742 - t683;
t645 = t418 * t557 + t649;
t434 = -pkin(4) * t668 + t466;
t640 = -g(1) * (-t547 * pkin(5) + qJ(6) * t548) - g(2) * (-t545 * pkin(5) + qJ(6) * t546) - g(3) * (-pkin(5) * t743 + qJ(6) * t742);
t637 = t638 * t631;
t613 = -pkin(4) * t775 - pkin(5);
t609 = pkin(4) * t627 + qJ(6);
t598 = t633 * t615;
t564 = t617 * t732 + t736;
t562 = -t617 * t733 + t735;
t510 = -pkin(5) * t660 - t587 * qJ(6) + t594;
t509 = t660 * t584;
t508 = t587 * t584;
t442 = pkin(5) * t508 - qJ(6) * t509 + t497;
t441 = t574 * t700 - t627 * t697 - t713 * t745 + (t574 * t627 + t584 * t786) * t631;
t440 = t534 * t584 - t574 * t660;
t436 = pkin(4) * t541 + t449;
t428 = pkin(5) * t658 - t665;
t427 = -qJ(6) * t658 + t792;
t417 = t557 * qJ(6) + t419;
t415 = -t557 * pkin(5) + t712;
t411 = pkin(5) * t441 + qJ(6) * t440 - qJD(6) * t509 + t465;
t407 = t430 * pkin(5) + t429 * qJ(6) - qJD(6) * t662 + t434;
t406 = -t571 * pkin(5) - t778;
t405 = qJ(6) * t571 - qJD(6) * t658 + t655;
t404 = t682 - t782;
t403 = t683 + t788;
t1 = [(-t470 * t584 + t471 * t658 + t498 * t572 + t499 * t570 - t526 * t574 - t527 * t571 + t539 * t528 + t538 * t641 - t677) * MDP(11) + (-t521 * t745 + t540 * t571 - t560 * t800 - t658 * t668) * MDP(16) + (t540 * t785 - t541 * t800 - t638 * t745 + t668 * t744) * MDP(14) + (t506 * t584 + t541 * t571 + t560 * t785 - t638 * t658) * MDP(15) + (-g(1) * t561 - g(2) * t563 - t458 * t571 + t466 * t744 + t498 * t541 + t513 * t785 - t521 * t720 + t538 * t638 - t560 * t656 + t657 * t658) * MDP(19) + (t418 * t571 + t497 * t430 + t434 * t508 + t475 * t441 + t465 * t477 + t665 * t516 + t557 * t778 + t658 * t682 + t678) * MDP(25) + (-t516 * t658 + t557 * t571) * MDP(24) + (t430 * t658 - t441 * t557 - t477 * t571 - t508 * t516) * MDP(23) + (-t521 * t658 + t560 * t571) * MDP(17) + ((-t539 * t714 + t494) * t560 + t512 * t521 - t674 * t658 + t457 * t571 - t498 * t540 - t538 * t668 + t584 * t694 - g(1) * t562 - g(2) * t564 + ((-qJD(4) * t525 - t499) * t560 - t539 * t521 - t684 * t658 + t466 * t584 + t513 * t574) * t628) * MDP(18) + (t429 * t658 - t440 * t557 + t509 * t516 + t571 * t662) * MDP(22) + (-t403 * t658 + t405 * t557 - t407 * t509 - t411 * t662 + t417 * t571 + t427 * t516 + t429 * t442 + t433 * t440 - t679) * MDP(29) + (t404 * t658 - t406 * t557 + t407 * t508 + t411 * t477 - t415 * t571 - t428 * t516 + t430 * t442 + t433 * t441 + t678) * MDP(27) + (-t403 * t508 + t404 * t509 - t405 * t477 + t406 * t662 - t415 * t440 - t417 * t441 - t427 * t430 - t428 * t429 + t616 * t787) * MDP(28) + t787 * MDP(2) + (-t419 * t571 - t497 * t429 + t434 * t509 - t475 * t440 + t465 * t662 - t516 * t792 - t557 * t655 + t658 * t683 + t679) * MDP(26) + qJDD(1) * MDP(1) + 0.2e1 * (t629 * t708 - t711 * t718) * MDP(5) + (t471 * t539 + t527 * t499 - t470 * t538 - t526 * t498 - t654 * t615 + t593 * t704 - g(1) * (-t615 * t630 - t626 * t633) - g(2) * (-t626 * t630 + t598)) * MDP(12) + (t541 * t785 + t584 * t637) * MDP(13) + (qJDD(1) * t622 + 0.2e1 * t632 * t691) * MDP(4) + (t403 * t427 + t417 * t405 + t407 * t442 + t433 * t411 + t404 * t428 + t415 * t406 - g(1) * (-pkin(5) * t546 - qJ(6) * t545) - g(2) * (pkin(5) * t548 + qJ(6) * t547 + t598) + (-g(1) * t690 - g(2) * t672) * t633 + (-g(1) * (-t615 - t672) - g(2) * t690) * t630) * MDP(30) + (qJDD(2) * t629 + t632 * t635) * MDP(6) + (qJDD(2) * t632 - t629 * t635) * MDP(7) + (t429 * t508 - t430 * t509 + t440 * t477 - t441 * t662) * MDP(21) + (-t429 * t509 - t440 * t662) * MDP(20) + (t629 * t667 + t632 * t651) * MDP(9) + (-t629 * t651 + t632 * t667) * MDP(10) + t677 * MDP(3); (-MDP(4) * t629 * t632 + MDP(5) * t718) * t636 + (t528 * t774 - t641 * t698 - (-t527 + t530) * t572 + (t526 - t531) * t570) * MDP(11) + MDP(7) * t708 + MDP(6) * t709 + (t628 * t668 + t637 - t804 * t541 + (t714 - t748) * t540) * MDP(14) + (t458 * t572 - t506 * t608 - t513 * t748 - t530 * t541 + t560 * t721 + t610 * t638 + t617 * t766 + t694 + (t791 - t808) * t628) * MDP(19) + (t419 * t572 - t594 * t429 + t434 * t587 - t722 * t475 + t557 * t790 + t681 * t662 - t797) * MDP(26) + (-t407 * t587 - t417 * t572 + t429 * t510 + t433 * t722 + t557 * t727 - t662 * t728 + t797) * MDP(29) + (t403 * t523 + t407 * t510 - t404 * t661 - g(3) * (-t741 + t763) - t669 * t767 + t728 * t433 + t727 * t417 + t726 * t415 + t677 * (pkin(2) * t629 + t616 * t669 + t617 * t634)) * MDP(30) + (t403 * t660 + t404 * t587 - t415 * t722 + t417 * t723 + t429 * t661 - t430 * t523 - t477 * t727 - t617 * t677 + t662 * t726 - t768) * MDP(28) + (-t418 * t572 + t594 * t430 - t434 * t660 - t723 * t475 + t681 * t477 + t557 * t789 + t653) * MDP(25) + (-t407 * t660 + t415 * t572 + t430 * t510 - t433 * t723 + t477 * t728 - t557 * t726 + t653) * MDP(27) + qJDD(2) * MDP(8) + (t526 * t530 - t527 * t531 + (t760 * t470 - t765 + t471 * t625 + (-qJD(1) * t593 + t677) * t629) * pkin(2)) * MDP(12) + (t629 * t650 - t765) * MDP(9) + (t610 * t702 - t457 * t572 + t530 * t540 - t608 * t754 + (-t610 * qJDD(2) + t652 - t791) * t631 + (-t493 + (t513 + t531) * t628) * t560) * MDP(18) + (t675 - t779) * MDP(21) + (-t750 - t780) * MDP(15) - t557 * t572 * MDP(24) - t560 * t572 * MDP(17) + (-t429 * t587 - t662 * t722) * MDP(20) + ((-qJD(4) * t572 + qJDD(2)) * t628 ^ 2 + (((t603 + qJD(4)) * qJD(2) + t647) * t628 + t751) * t631) * MDP(13) + (g(3) * t629 + t632 * t650) * MDP(10) + (t670 - t752) * MDP(16) + (-t496 * t557 + t724 - t755) * MDP(22) + (t495 * t557 + t725 + t756) * MDP(23); (-t570 ^ 2 - t572 ^ 2) * MDP(11) + (t526 * t572 - t527 * t570 + t654 - t787) * MDP(12) + (t670 + t752) * MDP(18) + (-t750 + t780) * MDP(19) + (t725 - t756) * MDP(25) + (-t484 - t755) * MDP(26) + (t483 - t756) * MDP(27) + (t675 + t779) * MDP(28) + (t724 + t755) * MDP(29) + (t403 * t587 - t404 * t660 - t415 * t723 - t417 * t722 - t433 * t572 - t787) * MDP(30) + (t495 * MDP(25) + MDP(26) * t722 + MDP(27) * t723 - t496 * MDP(29)) * t557; -t541 * t540 * MDP(13) + (-t540 ^ 2 + t541 ^ 2) * MDP(14) + (-t540 * t560 - t572 * t715 + t631 * t710 + t639) * MDP(15) + (t668 + t751) * MDP(16) + t521 * MDP(17) + (t458 * t560 - t513 * t541 + (t684 + t768) * t628 + t674 + t784) * MDP(18) + (g(1) * t564 - g(2) * t562 + t457 * t560 - t513 * t540 + t631 * t768 - t657) * MDP(19) + (t422 * t557 + (-t477 * t541 + t516 * t775 - t557 * t713) * pkin(4) + t781) * MDP(25) + (t423 * t557 + t801 + (-t516 * t627 - t541 * t662 - t557 * t692) * pkin(4) + t649) * MDP(26) + (-t436 * t477 - t516 * t613 - t557 * t680 - t643) * MDP(27) + (-t429 * t613 - t430 * t609 + (t417 + t680) * t662 + (t415 - t719) * t477) * MDP(28) + (t436 * t662 + t516 * t609 + t557 * t719 - t649 + t788 - t802) * MDP(29) + (t403 * t609 + t404 * t613 - t433 * t436 - t415 * t422 + t719 * t417 + (t415 * t713 + t616 * t766 + t784) * pkin(4) + t640) * MDP(30) + t803; (t758 + t781) * MDP(25) + (t645 + t801) * MDP(26) + (-t449 * t477 + t507 - t643 + t758) * MDP(27) + (pkin(5) * t429 - qJ(6) * t430 + (t417 - t419) * t662 + (t415 - t712) * t477) * MDP(28) + (t449 * t662 + 0.2e1 * t505 + 0.2e1 * t549 - t645 - t802) * MDP(29) + (-t404 * pkin(5) + t403 * qJ(6) - t415 * t419 + t417 * t712 - t433 * t449 + t640) * MDP(30) + t803; (-qJDD(4) - qJDD(5) + t528 + t807) * MDP(27) + t416 * MDP(28) + (-t557 ^ 2 - t776) * MDP(29) + (-t417 * t557 + t643) * MDP(30);];
tau  = t1;
