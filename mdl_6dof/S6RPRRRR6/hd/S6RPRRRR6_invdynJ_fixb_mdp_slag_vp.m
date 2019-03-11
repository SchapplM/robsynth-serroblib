% Calculate vector of inverse dynamics joint torques for
% S6RPRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRRR6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_invdynJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR6_invdynJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:14:56
% EndTime: 2019-03-09 07:15:17
% DurationCPUTime: 15.01s
% Computational Cost: add. (11139->625), mult. (26473->803), div. (0->0), fcn. (21458->18), ass. (0->273)
t685 = sin(pkin(11));
t686 = cos(pkin(11));
t690 = sin(qJ(3));
t695 = cos(qJ(3));
t637 = t685 * t695 + t686 * t690;
t625 = t637 * qJD(1);
t689 = sin(qJ(4));
t694 = cos(qJ(4));
t759 = t694 * qJD(3);
t591 = t625 * t689 - t759;
t593 = qJD(3) * t689 + t625 * t694;
t688 = sin(qJ(5));
t693 = cos(qJ(5));
t532 = t693 * t591 + t593 * t688;
t692 = cos(qJ(6));
t687 = sin(qJ(6));
t722 = t591 * t688 - t693 * t593;
t817 = t722 * t687;
t480 = -t692 * t532 + t817;
t627 = t637 * qJD(3);
t756 = qJDD(1) * t695;
t660 = t686 * t756;
t757 = qJDD(1) * t690;
t847 = -qJD(1) * t627 - t685 * t757 + t660;
t570 = qJDD(4) - t847;
t567 = qJDD(5) + t570;
t556 = qJDD(6) + t567;
t723 = t532 * t687 + t692 * t722;
t862 = t556 * MDP(33) + (-t480 ^ 2 + t723 ^ 2) * MDP(30) + t480 * MDP(29) * t723;
t768 = qJD(1) * t695;
t663 = t686 * t768;
t769 = qJD(1) * t690;
t750 = t685 * t769;
t624 = t663 - t750;
t833 = qJD(4) + qJD(5);
t861 = t624 - t833;
t752 = qJD(3) * t663 + t685 * t756 + t686 * t757;
t576 = -qJD(3) * t750 + t752;
t765 = qJD(4) * t689;
t524 = qJD(4) * t759 + t689 * qJDD(3) + t694 * t576 - t625 * t765;
t525 = qJD(4) * t593 - t694 * qJDD(3) + t576 * t689;
t762 = qJD(5) * t693;
t763 = qJD(5) * t688;
t463 = t693 * t524 - t688 * t525 - t591 * t762 - t593 * t763;
t701 = qJD(5) * t722 - t524 * t688 - t693 * t525;
t760 = qJD(6) * t692;
t753 = t692 * t463 - t532 * t760 + t687 * t701;
t761 = qJD(6) * t687;
t434 = t722 * t761 + t753;
t615 = qJD(4) - t624;
t611 = qJD(5) + t615;
t744 = t463 * t687 - t692 * t701;
t702 = qJD(6) * t723 - t744;
t605 = qJD(6) + t611;
t850 = t605 * t723;
t851 = t480 * t605;
t860 = t567 * MDP(26) + (-t532 ^ 2 + t722 ^ 2) * MDP(23) + (t532 * t611 + t463) * MDP(24) + (-t611 * t722 + t701) * MDP(25) - t532 * MDP(22) * t722 + (t702 - t850) * MDP(32) + (t434 - t851) * MDP(31) + t862;
t665 = -pkin(2) * t686 - pkin(1);
t650 = qJD(1) * t665 + qJD(2);
t543 = -pkin(3) * t624 - pkin(8) * t625 + t650;
t825 = pkin(7) + qJ(2);
t653 = t825 * t685;
t638 = qJD(1) * t653;
t654 = t825 * t686;
t639 = qJD(1) * t654;
t580 = -t690 * t638 + t695 * t639;
t572 = qJD(3) * pkin(8) + t580;
t506 = t694 * t543 - t572 * t689;
t489 = -pkin(9) * t593 + t506;
t470 = pkin(4) * t615 + t489;
t507 = t543 * t689 + t572 * t694;
t490 = -pkin(9) * t591 + t507;
t486 = t693 * t490;
t453 = t470 * t688 + t486;
t853 = pkin(10) * t532;
t447 = t453 - t853;
t445 = t447 * t761;
t838 = -t638 * t695 - t690 * t639;
t571 = -qJD(3) * pkin(3) - t838;
t526 = pkin(4) * t591 + t571;
t482 = pkin(5) * t532 + t526;
t684 = qJ(4) + qJ(5);
t680 = qJ(6) + t684;
t668 = sin(t680);
t669 = cos(t680);
t696 = cos(qJ(1));
t683 = pkin(11) + qJ(3);
t674 = cos(t683);
t691 = sin(qJ(1));
t801 = t674 * t691;
t595 = t668 * t696 - t669 * t801;
t800 = t674 * t696;
t597 = t668 * t691 + t669 * t800;
t673 = sin(t683);
t827 = g(3) * t673;
t845 = g(1) * t597 - g(2) * t595 - t482 * t480 + t669 * t827 + t445;
t640 = t688 * t689 - t693 * t694;
t775 = t861 * t640;
t791 = t688 * t694;
t641 = t689 * t693 + t791;
t774 = t861 * t641;
t594 = t668 * t801 + t669 * t696;
t596 = -t668 * t800 + t669 * t691;
t649 = qJDD(1) * t665 + qJDD(2);
t519 = -pkin(3) * t847 - pkin(8) * t576 + t649;
t511 = t694 * t519;
t758 = qJD(1) * qJD(2);
t829 = qJDD(1) * t825 + t758;
t609 = t829 * t685;
t610 = t829 * t686;
t721 = -t609 * t690 + t610 * t695;
t516 = qJDD(3) * pkin(8) + qJD(3) * t838 + t721;
t438 = pkin(4) * t570 - pkin(9) * t524 - qJD(4) * t507 - t516 * t689 + t511;
t764 = qJD(4) * t694;
t712 = t694 * t516 + t689 * t519 + t543 * t764 - t572 * t765;
t443 = -pkin(9) * t525 + t712;
t746 = t693 * t438 - t688 * t443;
t703 = -qJD(5) * t453 + t746;
t428 = pkin(5) * t567 - pkin(10) * t463 + t703;
t735 = -t688 * t438 - t693 * t443 - t470 * t762 + t490 * t763;
t429 = pkin(10) * t701 - t735;
t747 = t692 * t428 - t687 * t429;
t844 = -g(1) * t596 + g(2) * t594 + t482 * t723 + t668 * t827 + t747;
t573 = pkin(3) * t625 - pkin(8) * t624;
t559 = t694 * t573;
t828 = pkin(8) + pkin(9);
t751 = qJD(4) * t828;
t857 = pkin(4) * t625 - t689 * t838 + t559 + (-pkin(9) * t624 + t751) * t694;
t777 = t689 * t573 + t694 * t838;
t809 = t624 * t689;
t856 = -pkin(9) * t809 + t689 * t751 + t777;
t766 = qJD(3) * t695;
t767 = qJD(3) * t690;
t718 = -t609 * t695 - t690 * t610 + t638 * t767 - t639 * t766;
t517 = -qJDD(3) * pkin(3) - t718;
t732 = g(1) * t696 + g(2) * t691;
t707 = -g(3) * t674 + t673 * t732;
t855 = -qJD(4) * pkin(8) * t615 - t517 + t707;
t854 = t765 - t809;
t852 = pkin(10) * t722;
t582 = -t640 * t687 + t641 * t692;
t782 = qJD(6) * t582 + t687 * t775 - t692 * t774;
t748 = t637 * t764;
t636 = t685 * t690 - t695 * t686;
t626 = t636 * qJD(3);
t808 = t626 * t689;
t846 = t748 - t808;
t679 = cos(t684);
t797 = t679 * t691;
t678 = sin(t684);
t798 = t678 * t696;
t600 = -t674 * t797 + t798;
t796 = t679 * t696;
t799 = t678 * t691;
t602 = t674 * t796 + t799;
t843 = g(1) * t602 - g(2) * t600 + t526 * t532 + t679 * t827 + t735;
t599 = t674 * t799 + t796;
t601 = -t674 * t798 + t797;
t842 = -g(1) * t601 + g(2) * t599 + t526 * t722 + t678 * t827 + t703;
t561 = t641 * t637;
t575 = pkin(3) * t636 - pkin(8) * t637 + t665;
t564 = t694 * t575;
t590 = -t653 * t690 + t654 * t695;
t805 = t637 * t694;
t495 = pkin(4) * t636 - pkin(9) * t805 - t590 * t689 + t564;
t583 = t694 * t590;
t776 = t689 * t575 + t583;
t806 = t637 * t689;
t509 = -pkin(9) * t806 + t776;
t780 = t688 * t495 + t693 * t509;
t733 = pkin(4) * t854 - t580;
t837 = t857 * t693;
t589 = t653 * t695 + t690 * t654;
t655 = t828 * t689;
t656 = t828 * t694;
t773 = -t688 * t655 + t693 * t656;
t836 = t655 * t762 + t656 * t763 + t688 * t857 + t856 * t693;
t834 = qJ(2) * qJDD(1);
t731 = g(1) * t691 - g(2) * t696;
t832 = qJDD(2) - t731;
t581 = t692 * t640 + t641 * t687;
t783 = -qJD(6) * t581 + t687 * t774 + t692 * t775;
t831 = -t556 * t582 - t605 * t783;
t830 = -t567 * t641 - t611 * t775;
t824 = qJDD(1) * pkin(1);
t484 = t688 * t490;
t452 = t693 * t470 - t484;
t446 = t452 + t852;
t444 = pkin(5) * t611 + t446;
t823 = t444 * t692;
t822 = t480 * t625;
t821 = t723 * t625;
t820 = t524 * t689;
t819 = t532 * t625;
t818 = t722 * t625;
t815 = t556 * t688;
t813 = t591 * t615;
t812 = t591 * t625;
t811 = t593 * t615;
t810 = t593 * t625;
t807 = t626 * t694;
t670 = pkin(4) * t693 + pkin(5);
t802 = t670 * t556;
t794 = t686 * MDP(4);
t792 = t688 * t692;
t790 = t689 * t570;
t789 = t689 * t691;
t788 = t689 * t696;
t787 = t691 * t694;
t786 = t692 * t447;
t555 = t694 * t570;
t785 = t694 * t696;
t784 = t693 * t489 - t484;
t778 = -pkin(5) * t774 + t733;
t772 = t685 ^ 2 + t686 ^ 2;
t671 = -pkin(4) * t694 - pkin(3);
t749 = t637 * t765;
t544 = -t636 * qJD(2) - qJD(3) * t589;
t574 = pkin(3) * t627 + pkin(8) * t626;
t560 = t694 * t574;
t458 = pkin(9) * t807 + pkin(4) * t627 - t544 * t689 + t560 + (-t583 + (pkin(9) * t637 - t575) * t689) * qJD(4);
t711 = t694 * t544 + t689 * t574 + t575 * t764 - t590 * t765;
t466 = -pkin(9) * t846 + t711;
t745 = t693 * t458 - t466 * t688;
t743 = -t489 * t688 - t486;
t742 = t693 * t495 - t509 * t688;
t739 = -t693 * t655 - t656 * t688;
t738 = t615 * t694;
t737 = -qJD(4) * t543 - t516;
t736 = qJD(6) * t444 + t429;
t545 = qJD(2) * t637 - t653 * t767 + t654 * t766;
t734 = 0.2e1 * t772;
t730 = -t581 * t556 - t605 * t782;
t729 = -t640 * t567 + t611 * t774;
t546 = pkin(4) * t806 + t589;
t728 = -t572 * t764 + t511;
t554 = -pkin(10) * t640 + t773;
t727 = pkin(5) * t625 + pkin(10) * t775 + qJD(5) * t773 + qJD(6) * t554 - t688 * t856 + t837;
t553 = -pkin(10) * t641 + t739;
t726 = -pkin(10) * t774 - qJD(6) * t553 + t836;
t433 = t687 * t444 + t786;
t562 = t640 * t637;
t513 = t692 * t561 - t562 * t687;
t514 = -t561 * t687 - t562 * t692;
t719 = t824 - t832;
t717 = -t615 * t854 + t555;
t518 = pkin(4) * t846 + t545;
t715 = -t749 - t807;
t714 = -pkin(8) * t570 + t571 * t615;
t710 = t688 * t458 + t693 * t466 + t495 * t762 - t509 * t763;
t467 = pkin(4) * t525 + t517;
t700 = t734 * t758 - t732;
t619 = t674 * t785 + t789;
t618 = -t674 * t788 + t787;
t617 = -t674 * t787 + t788;
t616 = t674 * t789 + t785;
t608 = pkin(5) * t640 + t671;
t515 = pkin(4) * t593 - pkin(5) * t722;
t512 = pkin(5) * t561 + t546;
t473 = -t626 * t791 - t688 * t749 - t763 * t806 + (t805 * t833 - t808) * t693;
t472 = -t561 * t833 + t640 * t626;
t459 = pkin(5) * t473 + t518;
t455 = -pkin(10) * t561 + t780;
t454 = pkin(5) * t636 + pkin(10) * t562 + t742;
t449 = t784 + t852;
t448 = t743 + t853;
t441 = qJD(6) * t514 + t472 * t687 + t692 * t473;
t440 = -qJD(6) * t513 + t472 * t692 - t473 * t687;
t439 = -pkin(5) * t701 + t467;
t432 = -t447 * t687 + t823;
t431 = -pkin(10) * t473 + t710;
t430 = pkin(5) * t627 - pkin(10) * t472 - qJD(5) * t780 + t745;
t1 = [(-t434 * t513 + t440 * t480 + t441 * t723 + t514 * t702) * MDP(30) + (-t441 * t605 + t480 * t627 - t513 * t556 + t636 * t702) * MDP(32) + ((t430 * t692 - t431 * t687) * t605 + (t454 * t692 - t455 * t687) * t556 + t747 * t636 + t432 * t627 - t459 * t480 - t512 * t702 + t439 * t513 + t482 * t441 - g(1) * t595 - g(2) * t597 + ((-t454 * t687 - t455 * t692) * t605 - t433 * t636) * qJD(6)) * MDP(34) + (-MDP(5) * t685 + t794) * (t719 + t824) + (-t463 * t561 - t472 * t532 + t473 * t722 - t562 * t701) * MDP(23) + (-t473 * t611 - t532 * t627 - t561 * t567 + t636 * t701) * MDP(25) + (t745 * t611 + t742 * t567 + t746 * t636 + t452 * t627 + t518 * t532 - t546 * t701 + t467 * t561 + t526 * t473 - g(1) * t600 - g(2) * t602 + (-t453 * t636 - t611 * t780) * qJD(5)) * MDP(27) + (-qJD(3) * t627 - qJDD(3) * t636) * MDP(11) + (-t576 * t636 - t624 * t626 - t625 * t627 + t637 * t847) * MDP(9) + (-qJD(3) * t545 - qJDD(3) * t589 + t627 * t650 + t636 * t649 - t665 * t847 + t674 * t731) * MDP(13) + (-t525 * t636 - t591 * t627 - t615 * t846 - t637 * t790) * MDP(18) + (t576 * t637 - t625 * t626) * MDP(8) + (-qJD(3) * t626 + qJDD(3) * t637) * MDP(10) + (-(-t591 * t694 - t593 * t689) * t626 + (-t820 - t525 * t694 + (t591 * t689 - t593 * t694) * qJD(4)) * t637) * MDP(16) + ((-t590 * t764 + t560) * t615 + t564 * t570 + t728 * t636 + t506 * t627 + t545 * t591 + t589 * t525 + t571 * t748 - g(1) * t617 - g(2) * t619 + ((-qJD(4) * t575 - t544) * t615 - t590 * t570 + t737 * t636 + t517 * t637 - t571 * t626) * t689) * MDP(20) + (-qJD(3) * t544 - qJDD(3) * t590 + t576 * t665 - t626 * t650 + t637 * t649 - t673 * t731) * MDP(14) + (t434 * t636 + t440 * t605 + t514 * t556 - t627 * t723) * MDP(31) + (t434 * t514 - t440 * t723) * MDP(29) + (-g(1) * t594 - g(2) * t596 - t433 * t627 + t512 * t434 + t439 * t514 + t482 * t440 + t445 * t636 - t459 * t723 + (-(-qJD(6) * t455 + t430) * t605 - t454 * t556 - t428 * t636) * t687 + (-(qJD(6) * t454 + t431) * t605 - t455 * t556 - t736 * t636) * t692) * MDP(35) + (t463 * t636 + t472 * t611 - t562 * t567 - t627 * t722) * MDP(24) + (-t463 * t562 - t472 * t722) * MDP(22) + (-g(1) * t599 - g(2) * t601 - t453 * t627 + t546 * t463 - t467 * t562 + t526 * t472 - t518 * t722 - t567 * t780 - t611 * t710 + t636 * t735) * MDP(28) + (t556 * t636 + t605 * t627) * MDP(33) + (t567 * t636 + t611 * t627) * MDP(26) + (t570 * t636 + t615 * t627) * MDP(19) + qJDD(1) * MDP(1) + (pkin(1) * t719 + (t772 * t834 + t700) * qJ(2)) * MDP(7) + (t734 * t834 + t700) * MDP(6) + (t524 * t805 + t593 * t715) * MDP(15) + (-g(1) * t616 - g(2) * t618 - t507 * t627 + t517 * t805 + t589 * t524 + t545 * t593 - t570 * t776 + t571 * t715 - t615 * t711 - t636 * t712) * MDP(21) + (t524 * t636 + t555 * t637 + t593 * t627 + t615 * t715) * MDP(17) + t731 * MDP(2) + t732 * MDP(3); t832 * MDP(7) - t660 * MDP(13) + t752 * MDP(14) + (t717 - t812) * MDP(20) + (-t615 ^ 2 * t694 - t790 - t810) * MDP(21) + (t729 - t819) * MDP(27) + (t818 + t830) * MDP(28) + (t730 + t822) * MDP(34) + (t821 + t831) * MDP(35) + (-t794 - pkin(1) * MDP(7) + (MDP(13) * t690 + MDP(5)) * t685) * qJDD(1) + ((t685 * t768 + t686 * t769 + t625) * MDP(13) + (t624 - t750) * MDP(14)) * qJD(3) + (-MDP(7) * qJ(2) - MDP(6)) * qJD(1) ^ 2 * t772; -t624 ^ 2 * MDP(9) + (-t434 * t581 + t480 * t783 + t582 * t702 + t723 * t782) * MDP(30) + ((t553 * t692 - t554 * t687) * t556 - t608 * t702 + t439 * t581 + (t687 * t726 - t692 * t727) * t605 + t782 * t482 - t778 * t480 + t707 * t669) * MDP(34) + (-(t553 * t687 + t554 * t692) * t556 + t608 * t434 + t439 * t582 + (t687 * t727 + t692 * t726) * t605 + t783 * t482 - t778 * t723 - t707 * t668) * MDP(35) + (t671 * t463 + t467 * t641 + t775 * t526 - t773 * t567 + t611 * t836 - t678 * t707 - t722 * t733) * MDP(28) + (t730 - t822) * MDP(32) + t847 * MDP(11) + (-pkin(3) * t525 - t559 * t615 - t580 * t591 + (t615 * t838 + t714) * t689 + t855 * t694) * MDP(20) + (-pkin(3) * t524 - t580 * t593 + t777 * t615 - t689 * t855 + t714 * t694) * MDP(21) + (-t624 * t650 + t732 * t674 - t721 + t827) * MDP(14) + (t434 * t582 - t723 * t783) * MDP(29) + (t463 * t641 - t722 * t775) * MDP(22) + (t818 - t830) * MDP(24) + (t821 - t831) * MDP(31) + (qJD(3) * t580 + t707 + t718) * MDP(13) + (qJD(3) * MDP(11) - t650 * MDP(13) - t615 * MDP(19) - t506 * MDP(20) + t507 * MDP(21) - t611 * MDP(26) - t452 * MDP(27) + t453 * MDP(28) - t605 * MDP(33) - t432 * MDP(34) + t433 * MDP(35) - MDP(8) * t624 + MDP(9) * t625) * t625 + qJDD(3) * MDP(12) + (t593 * t738 + t820) * MDP(15) + (t729 + t819) * MDP(25) + (t717 + t812) * MDP(18) + ((t524 - t813) * t694 + (-t525 - t811) * t689) * MDP(16) + (t615 * t738 + t790 - t810) * MDP(17) + (t739 * t567 - t671 * t701 + t467 * t640 + (-t656 * t762 + (qJD(5) * t655 + t856) * t688 - t837) * t611 + t733 * t532 - t774 * t526 + t707 * t679) * MDP(27) + (-t463 * t640 - t532 * t775 + t641 * t701 - t722 * t774) * MDP(23) + ((-t624 - t750) * qJD(3) + t752) * MDP(10); (-t525 + t811) * MDP(18) + (t515 * t723 + (-t802 - t428 + (t448 - (-qJD(5) - qJD(6)) * t688 * pkin(4)) * t605) * t687 + (-pkin(4) * t815 + (-pkin(4) * t762 - qJD(6) * t670 + t449) * t605 - t736) * t692 + t845) * MDP(35) + (-t743 * t611 + (-t532 * t593 + t567 * t693 - t611 * t763) * pkin(4) + t842) * MDP(27) + (t784 * t611 + (-t567 * t688 + t593 * t722 - t611 * t762) * pkin(4) + t843) * MDP(28) + t570 * MDP(19) + (g(1) * t619 - g(2) * t617 + t506 * t615 + t571 * t591 + t694 * t827 - t712) * MDP(21) + (-g(1) * t618 + g(2) * t616 + t507 * t615 - t571 * t593 + (t737 + t827) * t689 + t728) * MDP(20) + (t524 + t813) * MDP(17) + (t692 * t802 - (t448 * t692 - t449 * t687) * t605 + t515 * t480 + (-t687 * t815 + (-t687 * t693 - t792) * t605 * qJD(5)) * pkin(4) + ((-pkin(4) * t792 - t670 * t687) * t605 - t433) * qJD(6) + t844) * MDP(34) + t593 * t591 * MDP(15) + (-t591 ^ 2 + t593 ^ 2) * MDP(16) + t860; (t453 * t611 + t842) * MDP(27) + (t452 * t611 + t843) * MDP(28) + (-(-t446 * t687 - t786) * t605 - t433 * qJD(6) + (-t480 * t722 + t556 * t692 - t605 * t761) * pkin(5) + t844) * MDP(34) + ((-t447 * t605 - t428) * t687 + (t446 * t605 - t736) * t692 + (-t556 * t687 - t605 * t760 - t722 * t723) * pkin(5) + t845) * MDP(35) + t860; (t753 - t851) * MDP(31) + (-t744 - t850) * MDP(32) + (t433 * t605 + t844) * MDP(34) + (-t687 * t428 - t692 * t429 + t432 * t605 + t845) * MDP(35) + (MDP(31) * t817 + MDP(32) * t723 - MDP(34) * t433 - MDP(35) * t823) * qJD(6) + t862;];
tau  = t1;
