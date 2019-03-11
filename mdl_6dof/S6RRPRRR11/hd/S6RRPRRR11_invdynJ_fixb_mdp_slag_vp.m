% Calculate vector of inverse dynamics joint torques for
% S6RRPRRR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRR11_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRRR11_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR11_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR11_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR11_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RRPRRR11_invdynJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:33:11
% EndTime: 2019-03-09 14:33:31
% DurationCPUTime: 14.63s
% Computational Cost: add. (7508->650), mult. (15626->833), div. (0->0), fcn. (10826->14), ass. (0->291)
t677 = cos(qJ(4));
t672 = sin(qJ(4));
t774 = qJD(2) * t672;
t678 = cos(qJ(2));
t777 = qJD(1) * t678;
t608 = t677 * t777 + t774;
t747 = t672 * t777;
t772 = qJD(2) * t677;
t610 = -t747 + t772;
t671 = sin(qJ(5));
t676 = cos(qJ(5));
t534 = t676 * t608 + t610 * t671;
t675 = cos(qJ(6));
t670 = sin(qJ(6));
t716 = t608 * t671 - t676 * t610;
t816 = t716 * t670;
t477 = -t675 * t534 + t816;
t759 = qJD(1) * qJD(2);
t745 = t678 * t759;
t673 = sin(qJ(2));
t758 = qJDD(1) * t673;
t704 = t745 + t758;
t607 = qJDD(4) + t704;
t599 = qJDD(5) + t607;
t582 = qJDD(6) + t599;
t718 = t534 * t670 + t675 * t716;
t858 = t582 * MDP(33) + (-t477 ^ 2 + t718 ^ 2) * MDP(30) + t477 * MDP(29) * t718;
t746 = t673 * t759;
t757 = qJDD(1) * t678;
t852 = t746 - t757;
t522 = -qJD(4) * t608 + t677 * qJDD(2) + t672 * t852;
t523 = -qJD(4) * t747 + qJDD(2) * t672 + (qJD(2) * qJD(4) - t852) * t677;
t763 = qJD(5) * t676;
t764 = qJD(5) * t671;
t456 = t676 * t522 - t671 * t523 - t608 * t763 - t610 * t764;
t688 = qJD(5) * t716 - t522 * t671 - t676 * t523;
t761 = qJD(6) * t675;
t754 = t675 * t456 - t534 * t761 + t670 * t688;
t762 = qJD(6) * t670;
t433 = t716 * t762 + t754;
t778 = qJD(1) * t673;
t641 = qJD(4) + t778;
t634 = qJD(5) + t641;
t739 = t456 * t670 - t675 * t688;
t689 = qJD(6) * t718 - t739;
t626 = qJD(6) + t634;
t848 = t626 * t718;
t849 = t477 * t626;
t857 = t599 * MDP(26) + (-t534 ^ 2 + t716 ^ 2) * MDP(23) + (t534 * t634 + t456) * MDP(24) + (-t634 * t716 + t688) * MDP(25) - t534 * MDP(22) * t716 + (t689 - t848) * MDP(32) + (t433 - t849) * MDP(31) + t858;
t680 = -pkin(2) - pkin(8);
t821 = qJ(3) * t673;
t600 = t678 * t680 - pkin(1) - t821;
t557 = t600 * qJD(1);
t650 = pkin(7) * t778;
t845 = qJD(3) + t650;
t760 = pkin(3) * t778 + t845;
t564 = qJD(2) * t680 + t760;
t500 = -t557 * t672 + t677 * t564;
t487 = -pkin(9) * t610 + t500;
t480 = pkin(4) * t641 + t487;
t501 = t557 * t677 + t564 * t672;
t488 = -pkin(9) * t608 + t501;
t484 = t676 * t488;
t450 = t480 * t671 + t484;
t851 = pkin(10) * t534;
t444 = t450 - t851;
t441 = t444 * t762;
t651 = pkin(7) * t777;
t618 = pkin(3) * t777 + t651;
t666 = qJD(2) * qJ(3);
t589 = t666 + t618;
t542 = pkin(4) * t608 + t589;
t485 = pkin(5) * t534 + t542;
t669 = qJ(4) + qJ(5);
t659 = qJ(6) + t669;
t644 = sin(t659);
t645 = cos(t659);
t674 = sin(qJ(1));
t679 = cos(qJ(1));
t802 = t673 * t679;
t552 = t644 * t802 + t645 * t674;
t804 = t673 * t674;
t554 = -t644 * t804 + t645 * t679;
t663 = g(3) * t678;
t843 = g(1) * t552 - g(2) * t554 - t485 * t477 - t644 * t663 + t441;
t551 = -t644 * t674 + t645 * t802;
t553 = t644 * t679 + t645 * t804;
t640 = pkin(2) * t746;
t820 = qJ(3) * t678;
t721 = pkin(8) * t673 - t820;
t770 = qJD(3) * t673;
t693 = qJD(2) * t721 - t770;
t496 = qJD(1) * t693 + qJDD(1) * t600 + t640;
t639 = pkin(7) * t745;
t647 = pkin(7) * t758;
t756 = qJDD(3) + t647;
t744 = t639 + t756;
t529 = t704 * pkin(3) + qJDD(2) * t680 + t744;
t736 = -t496 * t672 + t677 * t529;
t691 = -qJD(4) * t501 + t736;
t437 = pkin(4) * t607 - pkin(9) * t522 + t691;
t767 = qJD(4) * t677;
t753 = -t677 * t496 - t672 * t529 - t564 * t767;
t768 = qJD(4) * t672;
t442 = -pkin(9) * t523 - t557 * t768 - t753;
t740 = t676 * t437 - t671 * t442;
t690 = -qJD(5) * t450 + t740;
t427 = pkin(5) * t599 - pkin(10) * t456 + t690;
t730 = -t671 * t437 - t676 * t442 - t480 * t763 + t488 * t764;
t428 = pkin(10) * t688 - t730;
t741 = t675 * t427 - t670 * t428;
t842 = -g(1) * t551 - g(2) * t553 + t485 * t718 + t645 * t663 + t741;
t655 = pkin(2) * t778;
t571 = qJD(1) * t721 + t655;
t597 = t677 * t618;
t806 = t672 * t673;
t712 = pkin(4) * t678 - pkin(9) * t806;
t823 = pkin(9) - t680;
t854 = qJD(1) * t712 - t571 * t672 - t823 * t768 + t597;
t622 = t823 * t677;
t750 = t677 * t778;
t786 = t677 * t571 + t672 * t618;
t853 = pkin(9) * t750 + qJD(4) * t622 + t786;
t611 = t671 * t677 + t672 * t676;
t706 = t611 * t673;
t833 = qJD(4) + qJD(5);
t788 = -qJD(1) * t706 - t833 * t611;
t830 = pkin(3) + pkin(7);
t850 = pkin(10) * t716;
t798 = t676 * t677;
t809 = t671 * t672;
t787 = -t671 * t768 - t672 * t764 + t676 * t750 - t778 * t809 + t798 * t833;
t723 = pkin(2) * t678 + t821;
t711 = pkin(1) + t723;
t590 = t711 * qJD(1);
t844 = qJDD(1) * t711;
t657 = sin(t669);
t658 = cos(t669);
t566 = t657 * t802 + t658 * t674;
t568 = -t657 * t804 + t658 * t679;
t841 = g(1) * t566 - g(2) * t568 + t534 * t542 - t657 * t663 + t730;
t565 = -t657 * t674 + t658 * t802;
t567 = t657 * t679 + t658 * t804;
t840 = -g(1) * t565 - g(2) * t567 + t542 * t716 + t658 * t663 + t690;
t628 = t830 * t673;
t614 = t677 * t628;
t743 = pkin(9) * t678 - t600;
t511 = pkin(4) * t673 + t672 * t743 + t614;
t613 = t672 * t628;
t784 = t677 * t600 + t613;
t797 = t677 * t678;
t518 = -pkin(9) * t797 + t784;
t790 = t671 * t511 + t676 * t518;
t837 = t788 * t675;
t836 = t854 * t676;
t714 = -t798 + t809;
t539 = -t611 * t670 - t675 * t714;
t621 = t823 * t672;
t783 = -t676 * t621 - t671 * t622;
t751 = -pkin(4) * t677 - pkin(3);
t782 = pkin(4) * t767 - t751 * t778 + t845;
t835 = -t621 * t764 + t622 * t763 + t671 * t854 + t853 * t676;
t581 = t677 * t607;
t834 = -t641 * t768 + t581;
t728 = g(1) * t674 - g(2) * t679;
t729 = g(1) * t679 + g(2) * t674;
t832 = t589 * t641 + t680 * t607;
t831 = t678 * t833;
t829 = pkin(4) * t671;
t824 = g(3) * t673;
t822 = pkin(7) * qJDD(2);
t819 = qJDD(2) * pkin(2);
t482 = t671 * t488;
t449 = t676 * t480 - t482;
t443 = t449 + t850;
t439 = pkin(5) * t634 + t443;
t818 = t439 * t675;
t817 = t522 * t677;
t538 = t675 * t611 - t670 * t714;
t815 = t538 * t582;
t814 = t539 * t582;
t813 = t608 * t641;
t812 = t610 * t641;
t810 = t670 * t582;
t808 = t671 * t675;
t807 = t672 * t607;
t805 = t672 * t678;
t803 = t673 * t677;
t801 = t674 * t677;
t800 = t675 * t444;
t799 = t675 * t582;
t796 = t677 * t679;
t643 = t672 * pkin(4) + qJ(3);
t794 = -qJD(6) * t538 - t670 * t787 + t837;
t793 = qJD(6) * t539 + t670 * t788 + t675 * t787;
t792 = t676 * t487 - t482;
t789 = pkin(5) * t787 + t782;
t629 = t830 * t678;
t667 = t673 ^ 2;
t668 = t678 ^ 2;
t781 = t667 - t668;
t776 = qJD(2) * t608;
t775 = qJD(2) * t610;
t773 = qJD(2) * t673;
t771 = qJD(2) * t678;
t769 = qJD(4) * t557;
t766 = qJD(4) * t678;
t765 = qJD(4) * t680;
t682 = qJD(1) ^ 2;
t755 = t673 * t678 * t682;
t648 = pkin(7) * t757;
t664 = qJDD(2) * qJ(3);
t665 = qJD(2) * qJD(3);
t752 = t648 + t664 + t665;
t585 = pkin(4) * t797 + t629;
t748 = t672 * t766;
t742 = -qJD(2) * pkin(2) + qJD(3);
t654 = pkin(2) * t773;
t549 = t654 + t693;
t619 = t830 * t771;
t733 = -t549 * t672 + t677 * t619;
t464 = t712 * qJD(2) + (t677 * t743 - t613) * qJD(4) + t733;
t703 = t677 * t549 - t600 * t768 + t672 * t619 + t628 * t767;
t468 = (t673 * t772 + t748) * pkin(9) + t703;
t738 = t676 * t464 - t468 * t671;
t737 = -t487 * t671 - t484;
t735 = t676 * t511 - t518 * t671;
t732 = t621 * t671 - t676 * t622;
t731 = qJD(6) * t439 + t428;
t617 = t830 * t773;
t727 = -t714 * t599 + t634 * t788;
t509 = -pkin(10) * t611 + t783;
t726 = pkin(5) * t777 + pkin(10) * t788 + qJD(5) * t783 + qJD(6) * t509 - t671 * t853 + t836;
t508 = pkin(10) * t714 + t732;
t725 = pkin(10) * t787 - qJD(6) * t508 + t835;
t724 = qJD(6) * t714 - t787;
t722 = pkin(2) * t673 - t820;
t430 = t670 * t439 + t800;
t572 = t714 * t678;
t573 = t611 * t678;
t717 = t675 * t572 + t573 * t670;
t507 = t572 * t670 - t573 * t675;
t620 = t650 + t742;
t627 = -t651 - t666;
t715 = t620 * t678 + t627 * t673;
t713 = t641 * t672;
t709 = -0.2e1 * pkin(1) * t759 - t822;
t708 = -t641 * t767 - t807;
t707 = -qJ(3) * t771 - t770;
t705 = -t599 * t611 - t634 * t787;
t702 = t671 * t464 + t676 * t468 + t511 * t763 - t518 * t764;
t700 = pkin(1) * t682 + t729;
t681 = qJD(2) ^ 2;
t699 = pkin(7) * t681 - t728;
t698 = 0.2e1 * qJD(2) * t590 + t822;
t696 = -t678 * t729 - t824;
t694 = 0.2e1 * qJDD(1) * pkin(1) - t699;
t532 = pkin(3) * t757 - qJD(1) * t617 + t752;
t692 = t532 + t696;
t687 = -t590 * t778 - t673 * t729 + t663 + t756;
t543 = -pkin(4) * t748 + (-pkin(7) + t751) * t773;
t479 = pkin(4) * t523 + t532;
t524 = qJD(1) * t707 + t640 - t844;
t576 = t654 + t707;
t685 = qJD(1) * t576 + t524 + t699 - t844;
t560 = pkin(7) * t746 - t752;
t570 = t744 - t819;
t683 = qJD(2) * t715 - t560 * t678 + t570 * t673 - t729;
t646 = pkin(4) * t676 + pkin(5);
t615 = -qJ(3) * t777 + t655;
t594 = -t672 * t804 + t796;
t593 = t672 * t679 + t673 * t801;
t592 = t672 * t802 + t801;
t591 = -t672 * t674 + t673 * t796;
t561 = pkin(5) * t611 + t643;
t528 = -pkin(5) * t572 + t585;
t495 = pkin(4) * t610 - pkin(5) * t716;
t490 = t611 * t831 - t714 * t773;
t489 = qJD(2) * t706 + t714 * t831;
t469 = -pkin(5) * t490 + t543;
t459 = pkin(10) * t572 + t790;
t458 = pkin(5) * t673 + pkin(10) * t573 + t735;
t448 = t792 + t850;
t447 = t737 + t851;
t446 = qJD(6) * t507 + t489 * t670 - t675 * t490;
t445 = qJD(6) * t717 + t489 * t675 + t490 * t670;
t438 = -pkin(5) * t688 + t479;
t432 = pkin(10) * t490 + t702;
t431 = pkin(5) * t771 - pkin(10) * t489 - qJD(5) * t790 + t738;
t429 = -t444 * t670 + t818;
t1 = [(t490 * t634 - t534 * t771 + t572 * t599 + t673 * t688) * MDP(25) + (t456 * t572 - t489 * t534 - t490 * t716 - t573 * t688) * MDP(23) + (t738 * t634 + t735 * t599 + t740 * t673 + t449 * t771 + t543 * t534 - t585 * t688 - t479 * t572 - t542 * t490 - g(1) * t568 - g(2) * t566 + (-t450 * t673 - t634 * t790) * qJD(5)) * MDP(27) + (pkin(7) * t683 - t590 * t576 + (-t524 + t728) * t711) * MDP(14) + ((t667 + t668) * qJDD(1) * pkin(7) + t683) * MDP(11) + (qJDD(1) * t667 + 0.2e1 * t673 * t745) * MDP(4) + ((-t608 * t672 + t610 * t677) * t773 + (-t817 + t523 * t672 + (t608 * t677 + t610 * t672) * qJD(4)) * t678) * MDP(16) + (t673 * t698 + t678 * t685) * MDP(12) + (-t673 * t685 + t678 * t698) * MDP(13) + ((t641 * t774 + t522) * t673 + (t708 + t775) * t678) * MDP(17) + 0.2e1 * (t673 * t757 - t759 * t781) * MDP(5) + t728 * MDP(2) + t729 * MDP(3) + (t673 * t709 + t678 * t694) * MDP(9) + (-t673 * t694 + t678 * t709) * MDP(10) + ((t641 * t772 - t523) * t673 + (-t776 - t834) * t678) * MDP(18) + (t733 * t641 + (-t600 * t672 + t614) * t607 + t736 * t673 - t617 * t608 + t629 * t523 + t532 * t797 - g(1) * t594 - g(2) * t592 + (t500 * t678 - t589 * t803) * qJD(2) + (-t501 * t673 - t589 * t805 - t641 * t784) * qJD(4)) * MDP(20) + (-t522 * t805 + (t672 * t773 - t677 * t766) * t610) * MDP(15) + (t433 * t717 + t445 * t477 + t446 * t718 + t507 * t689) * MDP(30) + (-t446 * t626 + t477 * t771 + t582 * t717 + t673 * t689) * MDP(32) + ((t431 * t675 - t432 * t670) * t626 + (t458 * t675 - t459 * t670) * t582 + t741 * t673 + t429 * t771 - t469 * t477 - t528 * t689 - t438 * t717 + t485 * t446 - g(1) * t554 - g(2) * t552 + ((-t458 * t670 - t459 * t675) * t626 - t430 * t673) * qJD(6)) * MDP(34) + (t607 * t673 + t641 * t771) * MDP(19) + (t599 * t673 + t634 * t771) * MDP(26) + (t582 * t673 + t626 * t771) * MDP(33) + qJDD(1) * MDP(1) + (t433 * t673 + t445 * t626 + t507 * t582 - t718 * t771) * MDP(31) + (-t430 * t771 + g(1) * t553 - g(2) * t551 + t528 * t433 + t438 * t507 + t441 * t673 + t485 * t445 - t469 * t718 + (-(-qJD(6) * t459 + t431) * t626 - t458 * t582 - t427 * t673) * t670 + (-(qJD(6) * t458 + t432) * t626 - t459 * t582 - t731 * t673) * t675) * MDP(35) + (t433 * t507 - t445 * t718) * MDP(29) + (t456 * t673 + t489 * t634 - t573 * t599 - t716 * t771) * MDP(24) + (g(1) * t567 - g(2) * t565 - t450 * t771 + t585 * t456 - t479 * t573 + t542 * t489 - t543 * t716 - t599 * t790 - t634 * t702 + t673 * t730) * MDP(28) + (-t456 * t573 - t489 * t716) * MDP(22) + (-t703 * t641 - t784 * t607 - t617 * t610 + t629 * t522 + g(1) * t593 - g(2) * t591 + ((qJD(2) * t589 + t769) * t672 + t753) * t673 + (-qJD(2) * t501 - t532 * t672 - t589 * t767) * t678) * MDP(21) + (qJDD(2) * t673 + t678 * t681) * MDP(6) + (qJDD(2) * t678 - t673 * t681) * MDP(7); (-t456 * t611 - t534 * t788 - t688 * t714 + t716 * t787) * MDP(23) + ((-t523 - t812) * t677 + (-t522 + t813) * t672) * MDP(16) + (-t722 * qJDD(1) + ((-t627 - t666) * t673 + (-t620 + t742) * t678) * qJD(1)) * MDP(11) + (t687 - 0.2e1 * t819) * MDP(12) + (qJ(3) * t523 - t597 * t641 + t760 * t608 + t832 * t677 + ((t571 - t765) * t641 + t692) * t672) * MDP(20) + (qJ(3) * t522 + t786 * t641 + t760 * t610 - t832 * t672 + (-t641 * t765 + t692) * t677) * MDP(21) - MDP(4) * t755 + (t673 * t700 - t647 - t663) * MDP(9) + (-pkin(7) * qJD(1) * t715 - t570 * pkin(2) - g(3) * t723 - t560 * qJ(3) - t627 * qJD(3) + t590 * t615 + t722 * t729) * MDP(14) + ((-t610 * t678 - t641 * t806) * qJD(1) + t834) * MDP(17) + ((t608 * t678 - t641 * t803) * qJD(1) + t708) * MDP(18) + ((t508 * t675 - t509 * t670) * t582 - t561 * t689 + t438 * t538 + (t670 * t725 - t675 * t726) * t626 + t793 * t485 - t789 * t477 + t696 * t644) * MDP(34) + (-t433 * t538 + t477 * t794 + t539 * t689 + t718 * t793) * MDP(30) + (-MDP(12) * t615 - t641 * MDP(19) - t500 * MDP(20) + t501 * MDP(21) + MDP(24) * t716 + t534 * MDP(25) - t634 * MDP(26) - t449 * MDP(27) + t450 * MDP(28) + MDP(31) * t718 - MDP(32) * t477 - t626 * MDP(33) - t429 * MDP(34) + t430 * MDP(35)) * t777 + (t732 * t599 - t643 * t688 + t479 * t611 + (t621 * t763 + (qJD(5) * t622 + t853) * t671 - t836) * t634 + t787 * t542 + t782 * t534 + t696 * t657) * MDP(27) + (t678 * t700 - t648 + t824) * MDP(10) + t727 * MDP(24) + t781 * MDP(5) * t682 + (-t626 * t793 - t815) * MDP(32) + (t626 * t794 + t814) * MDP(31) + qJDD(2) * MDP(8) + (t433 * t539 - t718 * t794) * MDP(29) + (-(t508 * t670 + t509 * t675) * t582 + t561 * t433 + t438 * t539 + (t670 * t726 + t675 * t725) * t626 + t794 * t485 - t789 * t718 + t696 * t645) * MDP(35) + (t643 * t456 - t479 * t714 + t788 * t542 - t783 * t599 + t634 * t835 + t696 * t658 - t716 * t782) * MDP(28) + (-t456 * t714 - t716 * t788) * MDP(22) + t705 * MDP(25) + (t648 + 0.2e1 * t664 + 0.2e1 * t665 + (qJD(1) * t615 - g(3)) * t673 + (-qJD(1) * t590 - t729) * t678) * MDP(13) + (-t610 * t713 + t817) * MDP(15) + MDP(7) * t757 + MDP(6) * t758; MDP(11) * t758 + (qJDD(2) + t755) * MDP(12) + (-t667 * t682 - t681) * MDP(13) + (qJD(2) * t627 + t639 + t687 - t819) * MDP(14) + (-t641 * t713 + t581 - t776) * MDP(20) + (-t641 ^ 2 * t677 - t775 - t807) * MDP(21) + (-qJD(2) * t534 + t727) * MDP(27) + (qJD(2) * t716 + t705) * MDP(28) + (t814 + qJD(2) * t477 + (-t611 * t761 + t670 * t724 + t837) * t626) * MDP(34) + (-t815 + qJD(2) * t718 + (t724 * t675 + (qJD(6) * t611 - t788) * t670) * t626) * MDP(35); (-g(1) * t591 - g(2) * t593 + g(3) * t797 + t501 * t641 - t589 * t610 + t691) * MDP(20) + (-t523 + t812) * MDP(18) + (t646 * t799 - (t447 * t675 - t448 * t670) * t626 + t495 * t477 + (-t671 * t810 + (-t670 * t676 - t808) * t626 * qJD(5)) * pkin(4) + ((-pkin(4) * t808 - t646 * t670) * t626 - t430) * qJD(6) + t842) * MDP(34) + (t522 + t813) * MDP(17) + (g(1) * t592 - g(2) * t594 + t500 * t641 + t589 * t608 + (t769 - t663) * t672 + t753) * MDP(21) + t610 * t608 * MDP(15) + (t495 * t718 + (-t646 * t582 - t427 + (t447 - (-qJD(5) - qJD(6)) * t829) * t626) * t670 + (-t582 * t829 + (-pkin(4) * t763 - qJD(6) * t646 + t448) * t626 - t731) * t675 + t843) * MDP(35) + (-t737 * t634 + (-t534 * t610 + t599 * t676 - t634 * t764) * pkin(4) + t840) * MDP(27) + (t792 * t634 + (-t599 * t671 + t610 * t716 - t634 * t763) * pkin(4) + t841) * MDP(28) + t607 * MDP(19) + (-t608 ^ 2 + t610 ^ 2) * MDP(16) + t857; (t450 * t634 + t840) * MDP(27) + (t449 * t634 + t841) * MDP(28) + (-(-t443 * t670 - t800) * t626 - t430 * qJD(6) + (-t477 * t716 - t626 * t762 + t799) * pkin(5) + t842) * MDP(34) + ((-t444 * t626 - t427) * t670 + (t443 * t626 - t731) * t675 + (-t626 * t761 - t716 * t718 - t810) * pkin(5) + t843) * MDP(35) + t857; (t754 - t849) * MDP(31) + (-t739 - t848) * MDP(32) + (t430 * t626 + t842) * MDP(34) + (-t670 * t427 - t675 * t428 + t429 * t626 + t843) * MDP(35) + (MDP(31) * t816 + MDP(32) * t718 - MDP(34) * t430 - MDP(35) * t818) * qJD(6) + t858;];
tau  = t1;
