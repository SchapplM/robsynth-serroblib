% Calculate vector of inverse dynamics joint torques for
% S6PRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6PRRRRR2_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRRR2_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR2_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR2_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR2_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR2_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR2_invdynJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6PRRRRR2_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:45:44
% EndTime: 2019-03-09 00:45:59
% DurationCPUTime: 11.84s
% Computational Cost: add. (7167->605), mult. (15908->823), div. (0->0), fcn. (12880->18), ass. (0->275)
t648 = qJD(3) + qJD(4);
t658 = sin(qJ(4));
t663 = cos(qJ(3));
t817 = cos(qJ(4));
t733 = qJD(2) * t817;
t659 = sin(qJ(3));
t761 = qJD(2) * t659;
t834 = -t658 * t761 + t663 * t733;
t839 = t648 * t834;
t647 = qJDD(3) + qJDD(4);
t660 = sin(qJ(2));
t654 = sin(pkin(6));
t764 = qJD(1) * t654;
t740 = t660 * t764;
t819 = pkin(8) + pkin(9);
t723 = qJD(2) * t819 + t740;
t655 = cos(pkin(6));
t763 = qJD(1) * t655;
t567 = t659 * t763 + t663 * t723;
t750 = t655 * qJDD(1);
t623 = t663 * t750;
t664 = cos(qJ(2));
t754 = qJD(1) * qJD(2);
t580 = qJDD(2) * pkin(8) + (qJDD(1) * t660 + t664 * t754) * t654;
t722 = pkin(9) * qJDD(2) + t580;
t480 = qJDD(3) * pkin(3) - qJD(3) * t567 - t722 * t659 + t623;
t566 = -t723 * t659 + t663 * t763;
t488 = qJD(3) * t566 + t659 * t750 + t722 * t663;
t809 = qJD(3) * pkin(3);
t555 = t566 + t809;
t732 = qJD(4) * t817;
t760 = qJD(4) * t658;
t714 = -t817 * t480 + t658 * t488 + t555 * t760 + t567 * t732;
t440 = -pkin(4) * t647 + t714;
t808 = cos(pkin(12));
t725 = t808 * t660;
t653 = sin(pkin(12));
t786 = t653 * t664;
t591 = t655 * t725 + t786;
t724 = t808 * t664;
t787 = t653 * t660;
t593 = -t655 * t787 + t724;
t652 = qJ(3) + qJ(4);
t643 = sin(t652);
t645 = cos(t652);
t726 = t654 * t808;
t785 = t654 * t660;
t788 = t653 * t654;
t829 = g(3) * (-t643 * t785 + t645 * t655) + g(2) * (-t591 * t643 - t645 * t726) + g(1) * (-t593 * t643 + t645 * t788);
t838 = t440 + t829;
t821 = qJD(5) + qJD(6);
t837 = -t834 + t821;
t683 = t659 * t809 - t740;
t689 = -t658 * t659 + t663 * t817;
t559 = t648 * t689;
t778 = t658 * t663;
t604 = t659 * t817 + t778;
t560 = t648 * t604;
t836 = pkin(4) * t560 - pkin(10) * t559 + t683;
t742 = qJD(3) * t819;
t605 = t659 * t742;
t606 = t663 * t742;
t616 = t819 * t659;
t618 = t819 * t663;
t690 = -t616 * t817 - t658 * t618;
t505 = qJD(4) * t690 - t605 * t817 - t658 * t606;
t739 = t664 * t764;
t569 = t689 * t739;
t835 = -t505 + t569;
t553 = t658 * t567;
t494 = t566 * t817 - t553;
t830 = -pkin(3) * t732 + t494;
t598 = -qJD(2) * t778 - t659 * t733;
t657 = sin(qJ(5));
t662 = cos(qJ(5));
t570 = -t598 * t657 - t662 * t648;
t661 = cos(qJ(6));
t656 = sin(qJ(6));
t695 = t598 * t662 - t648 * t657;
t797 = t695 * t656;
t501 = t661 * t570 - t797;
t589 = qJD(5) - t834;
t586 = qJD(6) + t589;
t833 = t501 * t586;
t696 = t570 * t656 + t661 * t695;
t832 = t586 * t696;
t781 = t656 * t662;
t603 = t657 * t661 + t781;
t831 = t837 * t603;
t601 = t656 * t657 - t661 * t662;
t769 = t837 * t601;
t554 = t817 * t567;
t492 = t658 * t555 + t554;
t484 = pkin(10) * t648 + t492;
t641 = -pkin(3) * t663 - pkin(2);
t587 = qJD(2) * t641 - t739;
t515 = -pkin(4) * t834 + pkin(10) * t598 + t587;
t461 = t484 * t662 + t515 * t657;
t447 = -pkin(11) * t570 + t461;
t756 = qJD(6) * t656;
t445 = t447 * t756;
t491 = t555 * t817 - t553;
t483 = -t648 * pkin(4) - t491;
t466 = t570 * pkin(5) + t483;
t548 = t591 * t645 - t643 * t726;
t550 = t593 * t645 + t643 * t788;
t584 = t643 * t655 + t645 * t785;
t590 = -t655 * t724 + t787;
t592 = t655 * t786 + t725;
t651 = qJ(5) + qJ(6);
t642 = sin(t651);
t644 = cos(t651);
t783 = t654 * t664;
t828 = t466 * t501 - g(1) * (-t550 * t644 - t592 * t642) - g(2) * (-t548 * t644 - t590 * t642) - g(3) * (-t584 * t644 + t642 * t783) + t445;
t727 = qJDD(2) * t817;
t751 = qJDD(2) * t663;
t516 = t658 * t751 + t659 * t727 + t839;
t757 = qJD(5) * t662;
t758 = qJD(5) * t657;
t471 = t662 * t516 + t598 * t758 + t657 * t647 + t648 * t757;
t752 = qJDD(2) * t659;
t703 = t658 * t752 - t663 * t727;
t517 = qJD(2) * t560 + t703;
t512 = qJDD(5) + t517;
t671 = t658 * t480 + t488 * t817 + t555 * t732 - t567 * t760;
t439 = t647 * pkin(10) + t671;
t731 = t660 * t754;
t704 = -qJDD(1) * t783 + t654 * t731;
t753 = qJD(2) * qJD(3);
t730 = t659 * t753;
t546 = pkin(3) * t730 + qJDD(2) * t641 + t704;
t456 = pkin(4) * t517 - pkin(10) * t516 + t546;
t455 = t662 * t456;
t673 = -qJD(5) * t461 - t657 * t439 + t455;
t422 = pkin(5) * t512 - pkin(11) * t471 + t673;
t472 = -qJD(5) * t695 + t516 * t657 - t662 * t647;
t686 = -t662 * t439 - t657 * t456 + t484 * t758 - t515 * t757;
t423 = -pkin(11) * t472 - t686;
t721 = t661 * t422 - t656 * t423;
t827 = t466 * t696 - g(1) * (-t550 * t642 + t592 * t644) - g(2) * (-t548 * t642 + t590 * t644) - g(3) * (-t584 * t642 - t644 * t783) + t721;
t511 = qJDD(6) + t512;
t826 = t511 * MDP(30) + (-t501 ^ 2 + t696 ^ 2) * MDP(27) - t501 * MDP(26) * t696;
t534 = t603 * t604;
t825 = t569 * t657 + t662 * t836;
t493 = t658 * t566 + t554;
t712 = pkin(3) * t760 - t493;
t575 = -t658 * t616 + t618 * t817;
t771 = qJD(4) * t575 - t604 * t739 - t658 * t605 + t606 * t817;
t545 = -pkin(4) * t689 - pkin(10) * t604 + t641;
t824 = -t545 * t757 + t575 * t758 - t657 * t836 + t835 * t662;
t796 = t834 * t657;
t823 = (t758 - t796) * pkin(5);
t544 = -pkin(4) * t598 - pkin(10) * t834;
t526 = pkin(3) * t761 + t544;
t822 = t657 * t526 + t662 * t830;
t665 = qJD(3) ^ 2;
t710 = g(1) * t592 + g(2) * t590;
t820 = 0.2e1 * qJDD(2) * pkin(2) - pkin(8) * t665 + t654 * (-g(3) * t664 + t731) - t704 + t710;
t719 = t471 * t656 + t661 * t472;
t436 = -qJD(6) * t696 + t719;
t818 = -pkin(10) - pkin(11);
t813 = g(3) * t654;
t812 = t662 * pkin(5);
t638 = pkin(3) * t658 + pkin(10);
t811 = -pkin(11) - t638;
t810 = qJD(2) * pkin(2);
t460 = -t484 * t657 + t662 * t515;
t446 = pkin(11) * t695 + t460;
t441 = pkin(5) * t589 + t446;
t806 = t441 * t661;
t805 = t447 * t661;
t804 = t471 * t657;
t803 = t483 * t834;
t802 = t512 * t662;
t801 = t559 * t657;
t800 = t559 * t662;
t799 = t570 * t589;
t798 = t695 * t589;
t795 = t834 * t662;
t794 = t604 * t657;
t793 = t604 * t662;
t792 = t642 * t645;
t791 = t644 * t645;
t790 = t645 * t657;
t789 = t645 * t664;
t784 = t654 * t663;
t782 = t656 * t422;
t780 = t657 * t512;
t779 = t657 * t664;
t565 = t662 * t575;
t777 = qJDD(1) - g(3);
t776 = -pkin(11) * t800 + pkin(5) * t560 - t505 * t657 + (-t565 + (pkin(11) * t604 - t545) * t657) * qJD(5) + t825;
t688 = t604 * t757 + t801;
t775 = pkin(11) * t688 + t824;
t774 = pkin(5) * t688 + t771;
t773 = t662 * t491 + t657 * t544;
t768 = t657 * t545 + t565;
t767 = t823 + t712;
t649 = t659 ^ 2;
t766 = -t663 ^ 2 + t649;
t762 = qJD(2) * t654;
t759 = qJD(5) * t589;
t755 = qJD(6) * t661;
t749 = pkin(11) * t796;
t745 = t654 * t779;
t744 = t662 * t783;
t743 = t661 * t471 - t656 * t472 - t570 * t755;
t741 = qJD(5) * t818;
t737 = t660 * t762;
t736 = t664 * t762;
t734 = t604 * t758;
t729 = t663 * t753;
t728 = qJD(5) * t811;
t720 = t460 * t598 + t483 * t758;
t717 = t589 * t662;
t716 = qJD(6) * t441 + t423;
t639 = -pkin(3) * t817 - pkin(4);
t713 = -t598 * pkin(5) - pkin(11) * t795;
t711 = -t492 + t823;
t709 = g(1) * t593 + g(2) * t591;
t520 = t662 * t526;
t646 = t662 * pkin(11);
t600 = t638 * t662 + t646;
t708 = qJD(6) * t600 - t657 * t830 - t662 * t728 + t520 + t713;
t533 = t662 * t544;
t617 = pkin(10) * t662 + t646;
t707 = qJD(6) * t617 - t491 * t657 - t662 * t741 + t533 + t713;
t599 = t811 * t657;
t706 = -qJD(6) * t599 - t657 * t728 - t749 + t822;
t615 = t818 * t657;
t705 = -qJD(6) * t615 - t657 * t741 - t749 + t773;
t426 = t441 * t656 + t805;
t539 = t662 * t545;
t467 = -pkin(5) * t689 - pkin(11) * t793 - t575 * t657 + t539;
t475 = -pkin(11) * t794 + t768;
t700 = t467 * t656 + t475 * t661;
t699 = -t638 * t512 - t803;
t594 = t655 * t663 - t659 * t785;
t595 = t655 * t659 + t660 * t784;
t530 = t658 * t594 + t595 * t817;
t513 = -t530 * t657 - t744;
t692 = -t530 * t662 + t745;
t698 = t513 * t661 + t656 * t692;
t697 = t513 * t656 - t661 * t692;
t694 = -g(1) * t653 + g(2) * t808;
t691 = t594 * t817 - t658 * t595;
t687 = -t734 + t800;
t435 = t695 * t756 + t743;
t682 = -t461 * t598 + t483 * t757 + t657 * t838;
t681 = g(3) * t783 - t710;
t611 = -t739 - t810;
t678 = -qJD(2) * t611 - t580 + t709;
t677 = t681 * t645;
t430 = pkin(5) * t472 + t440;
t676 = -t426 * t598 + t430 * t603 - t466 * t769 + t642 * t829;
t672 = -pkin(8) * qJDD(3) + (t611 + t739 - t810) * qJD(3);
t670 = t587 * t598 - t714 - t829;
t425 = -t447 * t656 + t806;
t669 = t425 * t598 + t430 * t601 + t466 * t831 - t644 * t829;
t668 = g(1) * t550 + g(2) * t548 + g(3) * t584 - t587 * t834 - t671;
t667 = (-t435 * t601 - t436 * t603 + t501 * t769 + t696 * t831) * MDP(27) + (t435 * t603 + t696 * t769) * MDP(26) + ((t471 - t799) * t662 + (-t472 + t798) * t657) * MDP(20) + (t511 * t603 - t586 * t769 - t598 * t696) * MDP(28) + (-t501 * t598 - t511 * t601 - t586 * t831) * MDP(29) + (-t695 * t717 + t804) * MDP(19) + (-t589 ^ 2 * t657 - t570 * t598 + t802) * MDP(22) + (t589 * t717 - t598 * t695 + t780) * MDP(21) + (t516 - t839) * MDP(14) + (-t703 + (-qJD(2) * t604 - t598) * t648) * MDP(15) + (t598 ^ 2 - t834 ^ 2) * MDP(13) + t647 * MDP(16) + (MDP(12) * t834 + MDP(23) * t589 + MDP(30) * t586) * t598;
t666 = qJD(2) ^ 2;
t640 = -pkin(4) - t812;
t613 = t639 - t812;
t564 = -qJD(3) * t595 - t659 * t736;
t563 = qJD(3) * t594 + t663 * t736;
t535 = t601 * t604;
t523 = pkin(5) * t794 - t690;
t465 = qJD(4) * t530 + t658 * t563 - t564 * t817;
t464 = qJD(4) * t691 + t563 * t817 + t658 * t564;
t451 = t559 * t781 - t656 * t734 - t756 * t794 + (t793 * t821 + t801) * t661;
t450 = -t534 * t821 - t601 * t559;
t443 = qJD(5) * t692 - t464 * t657 + t662 * t737;
t442 = qJD(5) * t513 + t464 * t662 + t657 * t737;
t1 = [t777 * MDP(1) + (qJD(3) * t564 + qJDD(3) * t594) * MDP(10) + (-qJD(3) * t563 - qJDD(3) * t595) * MDP(11) + (-t465 * t648 + t647 * t691) * MDP(17) + (-t464 * t648 - t530 * t647) * MDP(18) + (t443 * t589 + t465 * t570 - t472 * t691 + t512 * t513) * MDP(24) + (-t442 * t589 - t465 * t695 - t471 * t691 + t512 * t692) * MDP(25) + ((-qJD(6) * t697 - t442 * t656 + t443 * t661) * t586 + t698 * t511 + t465 * t501 - t691 * t436) * MDP(31) + (-(qJD(6) * t698 + t442 * t661 + t443 * t656) * t586 - t697 * t511 - t465 * t696 - t691 * t435) * MDP(32) + ((-qJDD(2) * MDP(4) + (-MDP(17) * t834 - t598 * MDP(18)) * qJD(2) + (-MDP(10) * t663 + MDP(11) * t659 - MDP(3)) * t666) * t660 + (qJDD(2) * MDP(3) - t666 * MDP(4) + (-t730 + t751) * MDP(10) + (-t729 - t752) * MDP(11) - t517 * MDP(17) - t516 * MDP(18)) * t664) * t654; (t516 * t689 - t517 * t604 + t559 * t834 + t560 * t598) * MDP(13) + (t517 * t641 - t546 * t689 + t560 * t587 + t647 * t690 - t648 * t771 - t683 * t834 - t677) * MDP(17) + (t516 * t641 + t546 * t604 + t559 * t587 - t575 * t647 - t683 * t598 + t681 * t643 + t648 * t835) * MDP(18) + qJDD(2) * MDP(2) + (t672 * t659 + t663 * t820) * MDP(10) + (-t659 * t820 + t672 * t663) * MDP(11) + (-t768 * t512 - t686 * t689 - t461 * t560 - t690 * t471 + t440 * t793 - g(1) * (t592 * t790 + t593 * t662) - g(2) * (t590 * t790 + t591 * t662) - (-t645 * t779 + t660 * t662) * t813 + t824 * t589 - t771 * t695 + t687 * t483) * MDP(25) + (-t435 * t535 - t450 * t696) * MDP(26) + (-t435 * t534 + t436 * t535 - t450 * t501 + t451 * t696) * MDP(27) + ((-t570 * t662 + t657 * t695) * t559 + (-t804 - t472 * t662 + (t570 * t657 + t662 * t695) * qJD(5)) * t604) * MDP(20) + (t471 * t793 - t687 * t695) * MDP(19) + ((t467 * t661 - t475 * t656) * t511 - t721 * t689 + t425 * t560 + t523 * t436 + t430 * t534 + t466 * t451 - g(1) * (-t592 * t791 + t593 * t642) - g(2) * (-t590 * t791 + t591 * t642) - (t642 * t660 + t644 * t789) * t813 + (t656 * t775 + t661 * t776) * t586 + t774 * t501 + (t426 * t689 - t586 * t700) * qJD(6)) * MDP(31) + (t472 * t689 - t560 * t570 - t589 * t688 - t604 * t780) * MDP(22) + (-t560 * t648 + t647 * t689) * MDP(15) + (-t511 * t689 + t560 * t586) * MDP(30) + (t436 * t689 - t451 * t586 - t501 * t560 - t511 * t534) * MDP(29) + (-t512 * t689 + t560 * t589) * MDP(23) + (-t700 * t511 + (t716 * t661 - t445 + t782) * t689 - t426 * t560 + t523 * t435 - t430 * t535 + t466 * t450 - g(1) * (t592 * t792 + t593 * t644) - g(2) * (t590 * t792 + t591 * t644) - (-t642 * t789 + t644 * t660) * t813 + ((-qJD(6) * t467 + t775) * t661 + (qJD(6) * t475 - t776) * t656) * t586 - t774 * t696) * MDP(32) + (-t435 * t689 + t450 * t586 - t511 * t535 - t560 * t696) * MDP(28) + (-t471 * t689 + t512 * t793 - t560 * t695 + t589 * t687) * MDP(21) + (-t455 * t689 + t460 * t560 - t690 * t472 + t539 * t512 + t825 * t589 + t771 * t570 + (-t677 + (t483 * t604 + t484 * t689 - t575 * t589) * qJD(5)) * t662 + ((-qJD(5) * t545 - t505) * t589 - t575 * t512 - (-qJD(5) * t515 - t439) * t689 + t440 * t604 + t483 * t559 - g(3) * t785 - t709) * t657) * MDP(24) + (t777 * t783 + t710) * MDP(3) + (-t777 * t785 + t709) * MDP(4) + 0.2e1 * (t659 * t751 - t753 * t766) * MDP(6) + (qJDD(2) * t649 + 0.2e1 * t659 * t729) * MDP(5) + (qJDD(3) * t663 - t659 * t665) * MDP(8) + (qJDD(3) * t659 + t663 * t665) * MDP(7) + (t559 * t648 + t604 * t647) * MDP(14) + (t516 * t604 - t559 * t598) * MDP(12); (t639 * t472 - t520 * t589 + t712 * t570 + (t589 * t830 + t699) * t657 + (-t638 * t759 - t838) * t662 + t720) * MDP(24) + (t639 * t471 + t699 * t662 - t712 * t695 + (t638 * t758 + t822) * t589 + t682) * MDP(25) + MDP(7) * t752 + MDP(8) * t751 + (t494 * t648 + (t598 * t761 - t647 * t658 - t648 * t732) * pkin(3) + t668) * MDP(18) + (-g(3) * t594 + t659 * t678 + t694 * t784 + t623) * MDP(10) + (-(t599 * t656 + t600 * t661) * t511 + t613 * t435 + (t656 * t708 + t661 * t706) * t586 - t767 * t696 + t676) * MDP(32) + ((t599 * t661 - t600 * t656) * t511 + t613 * t436 + (t656 * t706 - t661 * t708) * t586 + t767 * t501 + t669) * MDP(31) + qJDD(3) * MDP(9) + (g(3) * t595 + (-t654 * t694 - t750) * t659 + t678 * t663) * MDP(11) + (t493 * t648 + (t647 * t817 - t648 * t760 + t761 * t834) * pkin(3) + t670) * MDP(17) + t667 + (-MDP(5) * t659 * t663 + MDP(6) * t766) * t666; (-(t615 * t656 + t617 * t661) * t511 + t640 * t435 + (t656 * t707 + t661 * t705) * t586 - t711 * t696 + t676) * MDP(32) + (-pkin(4) * t472 - t492 * t570 - t533 * t589 + (-pkin(10) * t512 + t491 * t589 - t803) * t657 + (-pkin(10) * t759 - t838) * t662 + t720) * MDP(24) + (-pkin(4) * t471 + t773 * t589 + t492 * t695 - t483 * t795 + (t589 * t758 - t802) * pkin(10) + t682) * MDP(25) + ((t615 * t661 - t617 * t656) * t511 + t640 * t436 + (t656 * t705 - t661 * t707) * t586 + t711 * t501 + t669) * MDP(31) + (t491 * t648 + t668) * MDP(18) + (t492 * t648 + t670) * MDP(17) + t667; -t695 * t570 * MDP(19) + (-t570 ^ 2 + t695 ^ 2) * MDP(20) + (t471 + t799) * MDP(21) + (-t472 - t798) * MDP(22) + t512 * MDP(23) + (t461 * t589 + t483 * t695 - g(1) * (-t550 * t657 + t592 * t662) - g(2) * (-t548 * t657 + t590 * t662) - g(3) * (-t584 * t657 - t744) + t673) * MDP(24) + (t460 * t589 + t483 * t570 - g(1) * (-t550 * t662 - t592 * t657) - g(2) * (-t548 * t662 - t590 * t657) - g(3) * (-t584 * t662 + t745) + t686) * MDP(25) + (t435 + t833) * MDP(28) + (-t436 - t832) * MDP(29) + (-(-t446 * t656 - t805) * t586 - t426 * qJD(6) + (t501 * t695 + t511 * t661 - t586 * t756) * pkin(5) + t827) * MDP(31) + ((-t447 * t586 - t422) * t656 + (t446 * t586 - t716) * t661 + (-t511 * t656 - t586 * t755 - t695 * t696) * pkin(5) + t828) * MDP(32) + t826; (t743 + t833) * MDP(28) + (-t719 - t832) * MDP(29) + (t426 * t586 + t827) * MDP(31) + (-t661 * t423 + t425 * t586 - t782 + t828) * MDP(32) + (MDP(28) * t797 + MDP(29) * t696 - MDP(31) * t426 - MDP(32) * t806) * qJD(6) + t826;];
tau  = t1;
