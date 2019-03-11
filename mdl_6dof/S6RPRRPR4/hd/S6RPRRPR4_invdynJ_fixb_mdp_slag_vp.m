% Calculate vector of inverse dynamics joint torques for
% S6RPRRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPR4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPR4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RPRRPR4_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:10:36
% EndTime: 2019-03-09 05:10:51
% DurationCPUTime: 12.17s
% Computational Cost: add. (10748->555), mult. (26227->709), div. (0->0), fcn. (21246->18), ass. (0->257)
t642 = cos(pkin(10));
t648 = cos(qJ(3));
t731 = t648 * t642;
t640 = sin(pkin(10));
t645 = sin(qJ(3));
t732 = t640 * t645;
t578 = -t731 + t732;
t563 = t578 * qJD(1);
t580 = t640 * t648 + t642 * t645;
t564 = t580 * qJD(1);
t644 = sin(qJ(4));
t769 = cos(qJ(4));
t538 = t769 * t563 + t564 * t644;
t536 = qJD(6) + t538;
t639 = sin(pkin(11));
t641 = cos(pkin(11));
t643 = sin(qJ(6));
t647 = cos(qJ(6));
t577 = t639 * t643 - t647 * t641;
t810 = t536 * t577;
t716 = qJD(1) * qJD(3);
t705 = t645 * t716;
t714 = qJDD(1) * t648;
t793 = qJDD(1) * t645 + t648 * t716;
t709 = t640 * t714 + t642 * t793;
t545 = -t640 * t705 + t709;
t710 = -t640 * t793 - t642 * t705;
t666 = t642 * t714 + t710;
t672 = -t644 * t563 + t564 * t769;
t479 = qJD(4) * t672 + t644 * t545 - t769 * t666;
t475 = qJDD(6) + t479;
t579 = t639 * t647 + t641 * t643;
t813 = t579 * t475 - t536 * t810;
t638 = qJD(3) + qJD(4);
t526 = -t641 * t638 + t639 * t672;
t528 = t638 * t639 + t641 * t672;
t753 = t528 * t643;
t805 = -t647 * t526 - t753;
t812 = t536 * t805;
t566 = t579 * qJD(6);
t811 = t579 * t538 + t566;
t706 = qJD(4) * t769;
t722 = qJD(4) * t644;
t478 = t769 * t545 - t563 * t706 - t564 * t722 + t644 * t666;
t632 = qJDD(3) + qJDD(4);
t466 = t478 * t639 - t641 * t632;
t467 = t478 * t641 + t632 * t639;
t720 = qJD(6) * t647;
t712 = -t643 * t466 + t647 * t467 - t526 * t720;
t721 = qJD(6) * t643;
t421 = -t528 * t721 + t712;
t680 = t526 * t643 - t528 * t647;
t700 = t647 * t466 + t643 * t467;
t422 = -qJD(6) * t680 + t700;
t687 = -t577 * t475 - t536 * t811;
t750 = t538 * t638;
t752 = t672 * t638;
t755 = t680 * t672;
t756 = t805 * t672;
t809 = (t687 - t756) * MDP(29) + t632 * MDP(19) + (-t479 + t752) * MDP(18) - t538 ^ 2 * MDP(16) + (MDP(15) * t538 + MDP(16) * t672 - MDP(30) * t536) * t672 + (t478 + t750) * MDP(17) + (t755 + t813) * MDP(28) + t421 * t579 * MDP(26) + (-t421 * t577 - t579 * t422 - t810 * t805) * MDP(27) + (MDP(26) * t810 + MDP(27) * t811) * t680;
t799 = t538 * t639;
t808 = pkin(5) * t799;
t807 = pkin(9) * t799;
t761 = pkin(7) + qJ(2);
t598 = t761 * t640;
t583 = qJD(1) * t598;
t600 = t761 * t642;
t584 = qJD(1) * t600;
t678 = t583 * t645 - t584 * t648;
t531 = -pkin(8) * t563 - t678;
t524 = t769 * t531;
t781 = -t648 * t583 - t584 * t645;
t530 = -pkin(8) * t564 + t781;
t525 = qJD(3) * pkin(3) + t530;
t486 = t644 * t525 + t524;
t480 = qJ(5) * t638 + t486;
t708 = pkin(2) * t642 + pkin(1);
t588 = -qJD(1) * t708 + qJD(2);
t548 = pkin(3) * t563 + t588;
t487 = pkin(4) * t538 - qJ(5) * t672 + t548;
t440 = -t480 * t639 + t641 * t487;
t429 = pkin(5) * t538 - pkin(9) * t528 + t440;
t441 = t641 * t480 + t639 * t487;
t433 = -pkin(9) * t526 + t441;
t410 = t429 * t643 + t433 * t647;
t717 = qJD(1) * qJD(2);
t770 = qJDD(1) * t761 + t717;
t556 = t770 * t640;
t557 = t770 * t642;
t698 = -t648 * t556 - t645 * t557;
t473 = qJDD(3) * pkin(3) - pkin(8) * t545 + qJD(3) * t678 + t698;
t679 = -t645 * t556 + t648 * t557;
t477 = t666 * pkin(8) + qJD(3) * t781 + t679;
t694 = -t769 * t473 + t644 * t477 + t525 * t722 + t531 * t706;
t777 = -pkin(4) * t632 + qJDD(5);
t426 = t694 + t777;
t415 = pkin(5) * t466 + t426;
t523 = t644 * t531;
t485 = t525 * t769 - t523;
t476 = -t638 * pkin(4) + qJD(5) - t485;
t457 = t526 * pkin(5) + t476;
t637 = pkin(10) + qJ(3);
t629 = qJ(4) + t637;
t618 = cos(t629);
t636 = pkin(11) + qJ(6);
t625 = sin(t636);
t617 = sin(t629);
t646 = sin(qJ(1));
t649 = cos(qJ(1));
t689 = g(1) * t649 + g(2) * t646;
t675 = t617 * t689;
t764 = g(3) * t625;
t804 = t410 * t672 + t415 * t579 - t457 * t810 + t618 * t764 - t625 * t675;
t409 = t429 * t647 - t433 * t643;
t627 = cos(t636);
t763 = g(3) * t627;
t741 = t617 * t649;
t742 = t617 * t646;
t795 = g(1) * t741 + g(2) * t742;
t803 = -t409 * t672 + t415 * t577 + t457 * t811 - t618 * t763 + t627 * t795;
t801 = t440 * t538;
t800 = t536 * t680;
t628 = cos(t637);
t726 = t618 * pkin(4) + t617 * qJ(5);
t794 = pkin(3) * t628 + t726;
t759 = qJDD(1) * pkin(1);
t780 = g(1) * t646 - g(2) * t649;
t677 = -qJDD(2) + t759 + t780;
t791 = -g(3) * t618 + t795;
t508 = pkin(4) * t672 + qJ(5) * t538;
t654 = t644 * t473 + t477 * t769 + t525 * t706 - t531 * t722;
t424 = t632 * qJ(5) + t638 * qJD(5) + t654;
t532 = qJDD(2) - t710 * pkin(3) + (-pkin(1) + (-pkin(3) * t648 - pkin(2)) * t642) * qJDD(1);
t430 = t479 * pkin(4) - t478 * qJ(5) - qJD(5) * t672 + t532;
t407 = -t424 * t639 + t641 * t430;
t790 = -t441 * t538 - t407;
t738 = t618 * t649;
t739 = t618 * t646;
t711 = -g(1) * t738 - g(2) * t739 - g(3) * t617;
t789 = t548 * t538 - t654 - t711;
t785 = pkin(5) * t672;
t784 = t780 * t617;
t492 = t530 * t769 - t523;
t495 = pkin(3) * t564 + t508;
t446 = -t492 * t639 + t641 * t495;
t608 = pkin(3) * t706 + qJD(5);
t783 = -t608 * t639 - t446;
t447 = t641 * t492 + t639 * t495;
t782 = t608 * t641 - t447;
t491 = t644 * t530 + t524;
t691 = pkin(3) * t722 - t491;
t727 = -t645 * t598 + t648 * t600;
t779 = MDP(4) * t642 - t640 * MDP(5);
t778 = qJ(2) * qJDD(1);
t774 = -t440 * t672 + (-t426 + t791) * t641;
t740 = t618 * t639;
t773 = t441 * t672 + g(3) * t740 + (t426 - t675) * t639;
t661 = -t694 + t791;
t772 = -t548 * t672 + t661;
t568 = t580 * qJD(3);
t768 = pkin(3) * t568;
t762 = t641 * pkin(5);
t630 = t641 * pkin(9);
t408 = t641 * t424 + t639 * t430;
t406 = t408 * t641;
t758 = t479 * t639;
t757 = t479 * t641;
t567 = t578 * qJD(3);
t671 = -t578 * t769 - t644 * t580;
t511 = t671 * qJD(4) - t567 * t769 - t644 * t568;
t754 = t511 * t639;
t748 = t538 * t641;
t547 = -t644 * t578 + t580 * t769;
t747 = t547 * t639;
t746 = t547 * t641;
t737 = t625 * t646;
t736 = t625 * t649;
t735 = t627 * t646;
t734 = t627 * t649;
t575 = t648 * t598;
t660 = -qJD(3) * t575 + qJD(2) * t731 + (-qJD(2) * t640 - qJD(3) * t600) * t645;
t513 = -pkin(8) * t568 + t660;
t653 = -t580 * qJD(2) - qJD(3) * t727;
t514 = pkin(8) * t567 + t653;
t697 = -t600 * t645 - t575;
t533 = -pkin(8) * t580 + t697;
t534 = -pkin(8) * t578 + t727;
t673 = t533 * t769 - t644 * t534;
t444 = t673 * qJD(4) + t513 * t769 + t644 * t514;
t512 = t547 * qJD(4) - t644 * t567 + t568 * t769;
t450 = pkin(4) * t512 - qJ(5) * t511 - qJD(5) * t547 + t768;
t417 = t641 * t444 + t639 * t450;
t452 = t641 * t485 + t639 * t508;
t554 = pkin(3) * t578 - t708;
t499 = -pkin(4) * t671 - qJ(5) * t547 + t554;
t505 = t644 * t533 + t534 * t769;
t454 = t639 * t499 + t641 * t505;
t725 = t640 ^ 2 + t642 ^ 2;
t723 = qJD(3) * t564;
t718 = -qJD(5) + t476;
t707 = qJD(1) * t732;
t702 = t725 * qJD(1) ^ 2;
t404 = pkin(5) * t479 - pkin(9) * t467 + t407;
t405 = -pkin(9) * t466 + t408;
t701 = t647 * t404 - t405 * t643;
t416 = -t444 * t639 + t641 * t450;
t451 = -t485 * t639 + t641 * t508;
t453 = t641 * t499 - t505 * t639;
t696 = t406 + t711;
t693 = 0.2e1 * t725;
t623 = -pkin(3) * t769 - pkin(4);
t692 = t691 + t808;
t626 = sin(t637);
t690 = -pkin(3) * t626 - pkin(4) * t617;
t686 = t404 * t643 + t405 * t647;
t685 = -t407 * t641 - t408 * t639;
t684 = -t407 * t639 + t406;
t435 = -pkin(5) * t671 - pkin(9) * t746 + t453;
t443 = -pkin(9) * t747 + t454;
t683 = t435 * t647 - t443 * t643;
t682 = t435 * t643 + t443 * t647;
t681 = t440 * t639 - t441 * t641;
t676 = t708 + t794;
t674 = t780 * t618;
t619 = pkin(3) * t644 + qJ(5);
t572 = t619 * t641 + t630;
t670 = pkin(9) * t748 + qJD(6) * t572 - t783 + t785;
t571 = (-pkin(9) - t619) * t639;
t669 = -qJD(6) * t571 - t782 + t807;
t599 = qJ(5) * t641 + t630;
t668 = qJD(5) * t639 + qJD(6) * t599 + t538 * t630 + t451 + t785;
t597 = (-pkin(9) - qJ(5)) * t639;
t667 = -qJD(5) * t641 - qJD(6) * t597 + t452 + t807;
t663 = -t479 * t619 + (t476 - t608) * t538;
t658 = t426 * t547 + t476 * t511 - t689;
t656 = t693 * t717 - t689;
t445 = t505 * qJD(4) + t644 * t513 - t514 * t769;
t633 = -pkin(8) - t761;
t620 = -pkin(4) - t762;
t592 = t623 - t762;
t590 = qJ(5) * t738;
t589 = qJ(5) * t739;
t587 = -qJDD(1) * t708 + qJDD(2);
t552 = t618 * t734 + t737;
t551 = -t618 * t736 + t735;
t550 = -t618 * t735 + t736;
t549 = t618 * t737 + t734;
t510 = t577 * t547;
t509 = t579 * t547;
t465 = pkin(5) * t747 - t673;
t458 = t486 - t808;
t438 = t511 * t579 + t720 * t746 - t721 * t747;
t437 = -t511 * t577 - t547 * t566;
t431 = pkin(5) * t754 + t445;
t412 = -pkin(9) * t754 + t417;
t411 = pkin(5) * t512 - t511 * t630 + t416;
t1 = [(-t444 * t638 + t478 * t554 - t505 * t632 + t511 * t548 + t532 * t547 + t672 * t768 - t784) * MDP(21) + (-t416 * t528 - t417 * t526 - t453 * t467 - t454 * t466 + t784 + t685 * t547 + (-t440 * t641 - t441 * t639) * t511) * MDP(24) + (-t421 * t509 + t422 * t510 + t437 * t805 + t438 * t680) * MDP(27) + ((t411 * t647 - t412 * t643) * t536 + t683 * t475 - t701 * t671 + t409 * t512 - t431 * t805 + t465 * t422 + t415 * t509 + t457 * t438 - g(1) * t550 - g(2) * t552 + (t410 * t671 - t536 * t682) * qJD(6)) * MDP(31) + (t422 * t671 - t438 * t536 - t475 * t509 + t512 * t805) * MDP(29) + (-t421 * t510 - t437 * t680) * MDP(26) + (t407 * t453 + t408 * t454 + t440 * t416 + t441 * t417 - t426 * t673 + t476 * t445 + (g(1) * t633 - g(2) * t676) * t649 + (g(1) * t676 + g(2) * t633) * t646) * MDP(25) + (-t512 * t638 + t632 * t671) * MDP(18) + (t478 * t671 - t479 * t547 - t511 * t538 - t512 * t672) * MDP(16) + (-t475 * t671 + t512 * t536) * MDP(30) + (-(t411 * t643 + t412 * t647) * t536 - t682 * t475 + t686 * t671 - t410 * t512 - t431 * t680 + t465 * t421 - t415 * t510 + t457 * t437 - g(1) * t549 - g(2) * t551 + (t409 * t671 - t536 * t683) * qJD(6)) * MDP(32) + (-t421 * t671 + t437 * t536 - t475 * t510 - t512 * t680) * MDP(28) + (-t407 * t671 + t416 * t538 + t440 * t512 + t445 * t526 + t453 * t479 - t466 * t673 + t639 * t658 + t641 * t674) * MDP(22) + (-t445 * t638 + t479 * t554 + t512 * t548 - t532 * t671 + t538 * t768 + t632 * t673 + t674) * MDP(20) + (t408 * t671 - t417 * t538 - t441 * t512 + t445 * t528 - t454 * t479 - t467 * t673 + t641 * t658 - t740 * t780) * MDP(23) + (-qJD(3) * t660 - qJDD(3) * t727 - t545 * t708 - t588 * t567 + t587 * t580 - t626 * t780) * MDP(14) + (qJD(3) * t653 + qJDD(3) * t697 + t588 * t568 + t587 * t578 + t628 * t780 + t666 * t708) * MDP(13) + t780 * MDP(2) + t689 * MDP(3) + t779 * (t677 + t759) + (-qJD(3) * t567 + qJDD(3) * t580) * MDP(10) + (t545 * t580 - t564 * t567) * MDP(8) + (-t545 * t578 + t567 * t563 - t564 * t568 + t580 * t666) * MDP(9) + (t511 * t638 + t547 * t632) * MDP(17) + (t677 * pkin(1) + (t725 * t778 + t656) * qJ(2)) * MDP(7) + (t693 * t778 + t656) * MDP(6) + (-qJD(3) * t568 - qJDD(3) * t578) * MDP(11) + (t478 * t547 + t511 * t672) * MDP(15) + qJDD(1) * MDP(1); -MDP(6) * t702 + (-qJ(2) * t702 - t677) * MDP(7) + (-t666 + t723) * MDP(13) + ((-t563 - t707) * qJD(3) + t709) * MDP(14) + (t479 + t752) * MDP(20) + (t478 - t750) * MDP(21) + (-t526 * t672 - t538 * t799 + t757) * MDP(22) + (-t528 * t672 - t538 * t748 - t758) * MDP(23) + (-t466 * t639 - t467 * t641 - (t526 * t641 - t528 * t639) * t538) * MDP(24) + (-t476 * t672 - t538 * t681 - t685 - t780) * MDP(25) + (t687 + t756) * MDP(31) + (t755 - t813) * MDP(32) - t779 * qJDD(1); (t666 + t723) * MDP(11) + (t492 * t638 + (-t564 * t672 - t632 * t644 - t638 * t706) * pkin(3) + t789) * MDP(21) + ((t563 - t707) * qJD(3) + t709) * MDP(10) + (g(3) * t626 + t588 * t563 + t628 * t689 - t679) * MDP(14) + (t426 * t623 - g(1) * (t649 * t690 + t590) - g(2) * (t646 * t690 + t589) - g(3) * t794 + t684 * t619 + t691 * t476 + t782 * t441 + t783 * t440) * MDP(25) + (t491 * t638 + (-t538 * t564 + t632 * t769 - t638 * t722) * pkin(3) + t772) * MDP(20) + (t447 * t538 + t467 * t623 + t528 * t691 + t641 * t663 + t773) * MDP(23) + (-t446 * t538 + t466 * t623 + t526 * t691 + t639 * t663 + t774) * MDP(22) + (t446 * t528 + t447 * t526 + (-t466 * t619 - t526 * t608 - t801) * t641 + (t467 * t619 + t528 * t608 + t790) * t639 + t696) * MDP(24) + (-g(3) * t628 - t588 * t564 + t626 * t689 + t698) * MDP(13) + t564 * t563 * MDP(8) + (-t563 ^ 2 + t564 ^ 2) * MDP(9) + qJDD(3) * MDP(12) + (-(t571 * t643 + t572 * t647) * t475 + t592 * t421 + (t643 * t670 + t647 * t669) * t536 - t692 * t680 + t804) * MDP(32) + ((t571 * t647 - t572 * t643) * t475 + t592 * t422 + (t643 * t669 - t647 * t670) * t536 - t692 * t805 + t803) * MDP(31) + t809; (t486 * t638 + t772) * MDP(20) + (t485 * t638 + t789) * MDP(21) + (-qJ(5) * t758 - pkin(4) * t466 - t486 * t526 + (t639 * t718 - t451) * t538 + t774) * MDP(22) + (-qJ(5) * t757 - pkin(4) * t467 - t486 * t528 + (t641 * t718 + t452) * t538 + t773) * MDP(23) + (t451 * t528 + t452 * t526 + (-qJ(5) * t466 - qJD(5) * t526 - t801) * t641 + (qJ(5) * t467 + qJD(5) * t528 + t790) * t639 + t696) * MDP(24) + (-t426 * pkin(4) - t441 * t452 - t440 * t451 - t476 * t486 - g(1) * (-pkin(4) * t741 + t590) - g(2) * (-pkin(4) * t742 + t589) - g(3) * t726 - t681 * qJD(5) + t684 * qJ(5)) * MDP(25) + ((t597 * t647 - t599 * t643) * t475 + t620 * t422 + t458 * t805 + (t643 * t667 - t647 * t668) * t536 + t803) * MDP(31) + (-(t597 * t643 + t599 * t647) * t475 + t620 * t421 + t458 * t680 + (t643 * t668 + t647 * t667) * t536 + t804) * MDP(32) + t809; (t528 * t538 + t466) * MDP(22) + (-t526 * t538 + t467) * MDP(23) + (-t526 ^ 2 - t528 ^ 2) * MDP(24) + (t440 * t528 + t441 * t526 - t661 + t777) * MDP(25) + (t422 - t800) * MDP(31) + (t421 + t812) * MDP(32); t680 * t805 * MDP(26) + (t680 ^ 2 - t805 ^ 2) * MDP(27) + (t712 - t812) * MDP(28) + (-t700 - t800) * MDP(29) + t475 * MDP(30) + (-g(1) * t551 + g(2) * t549 + t410 * t536 + t457 * t680 + t617 * t764 + t701) * MDP(31) + (g(1) * t552 - g(2) * t550 + t409 * t536 - t457 * t805 + t617 * t763 - t686) * MDP(32) + (-MDP(28) * t753 + MDP(29) * t680 - MDP(31) * t410 - MDP(32) * t409) * qJD(6);];
tau  = t1;
