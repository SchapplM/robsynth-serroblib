% Calculate vector of inverse dynamics joint torques for
% S6RRRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP6_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPP6_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP6_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP6_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP6_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRRPP6_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:15:25
% EndTime: 2019-03-09 21:15:43
% DurationCPUTime: 13.73s
% Computational Cost: add. (9781->747), mult. (21008->872), div. (0->0), fcn. (14453->10), ass. (0->312)
t697 = cos(qJ(3));
t695 = sin(qJ(2));
t698 = cos(qJ(2));
t845 = t697 * t698;
t750 = pkin(3) * t695 - pkin(9) * t845;
t700 = -pkin(9) - pkin(8);
t796 = qJD(3) * t700;
t765 = pkin(2) * t695 - pkin(8) * t698;
t615 = t765 * qJD(1);
t694 = sin(qJ(3));
t824 = qJD(1) * t695;
t788 = t694 * t824;
t830 = pkin(7) * t788 + t697 * t615;
t914 = qJD(1) * t750 - t697 * t796 + t830;
t598 = t694 * t615;
t850 = t695 * t697;
t853 = t694 * t698;
t913 = -t598 - (-pkin(7) * t850 - pkin(9) * t853) * qJD(1) + t694 * t796;
t819 = qJD(2) * t697;
t611 = -t788 + t819;
t821 = qJD(2) * t694;
t612 = t697 * t824 + t821;
t693 = sin(qJ(4));
t883 = cos(qJ(4));
t547 = -t883 * t611 + t612 * t693;
t741 = t693 * t611 + t612 * t883;
t863 = t547 * t741;
t823 = qJD(1) * t698;
t794 = t694 * t823;
t795 = t883 * t697;
t856 = t693 * t694;
t889 = qJD(3) + qJD(4);
t784 = t883 * qJD(4);
t892 = t883 * qJD(3) + t784;
t834 = -t693 * t794 - t697 * t892 + t795 * t823 + t856 * t889;
t775 = qJD(3) + t823;
t804 = qJDD(1) * t695;
t722 = qJD(2) * t775 + t804;
t806 = qJD(1) * qJD(3);
t781 = t695 * t806;
t753 = -qJDD(2) + t781;
t736 = t753 * t694;
t707 = t722 * t697 - t736;
t748 = qJD(2) * qJD(3) + t804;
t807 = qJD(1) * qJD(2);
t782 = t698 * t807;
t909 = t748 + t782;
t772 = t694 * t909 + t697 * t781;
t730 = qJDD(2) * t697 - t772;
t814 = qJD(4) * t693;
t486 = -t611 * t784 + t612 * t814 - t693 * t730 - t883 * t707;
t664 = -qJD(3) + t823;
t647 = -qJD(4) + t664;
t468 = -t547 * t647 - t486;
t545 = t741 ^ 2;
t681 = t698 * qJDD(1);
t893 = -t695 * t807 + t681;
t608 = qJDD(3) - t893;
t603 = qJDD(4) + t608;
t487 = qJD(4) * t741 + t693 * t707 - t883 * t730;
t703 = -t647 * t741 - t487;
t884 = t547 ^ 2;
t912 = MDP(18) * t863 + t703 * MDP(21) + t603 * MDP(22) + t468 * MDP(20) + (t545 - t884) * MDP(19);
t867 = qJ(5) * t547;
t869 = qJD(2) * pkin(2);
t632 = pkin(7) * t824 - t869;
t567 = -pkin(3) * t611 + t632;
t724 = -qJ(5) * t741 + t567;
t496 = pkin(4) * t547 + t724;
t911 = t496 * t547;
t634 = t700 * t694;
t635 = t700 * t697;
t910 = t634 * t784 + t635 * t814 - t914 * t693 + t913 * t883;
t614 = t693 * t697 + t694 * t883;
t562 = t889 * t614;
t833 = -t614 * t823 + t562;
t678 = pkin(7) * t823;
t816 = qJD(3) * t694;
t767 = -t678 + (-t794 + t816) * pkin(3);
t818 = qJD(2) * t698;
t793 = t694 * t818;
t815 = qJD(3) * t697;
t908 = t695 * t815 + t793;
t699 = cos(qJ(1));
t849 = t695 * t699;
t696 = sin(qJ(1));
t851 = t695 * t696;
t907 = g(1) * t849 + g(2) * t851;
t627 = -pkin(2) * t698 - pkin(8) * t695 - pkin(1);
t604 = t627 * qJD(1);
t633 = qJD(2) * pkin(8) + t678;
t557 = t697 * t604 - t633 * t694;
t530 = -pkin(9) * t612 + t557;
t524 = -pkin(3) * t664 + t530;
t558 = t694 * t604 + t697 * t633;
t531 = pkin(9) * t611 + t558;
t857 = t693 * t531;
t489 = -t883 * t524 + t857;
t901 = pkin(5) * t741;
t745 = t489 + t901;
t809 = qJD(5) + t745;
t890 = pkin(5) * t547 - qJD(6);
t871 = pkin(4) + qJ(6);
t480 = t547 * t871 + t724;
t888 = -pkin(5) * t487 + qJDD(6);
t905 = -t480 * t547 + t888;
t691 = qJ(3) + qJ(4);
t683 = cos(t691);
t682 = sin(t691);
t843 = t699 * t682;
t846 = t696 * t698;
t578 = t683 * t846 - t843;
t844 = t698 * t699;
t580 = t682 * t696 + t683 * t844;
t618 = t765 * qJD(2);
t563 = qJD(1) * t618 + qJDD(1) * t627;
t553 = t697 * t563;
t588 = pkin(7) * t893 + qJDD(2) * pkin(8);
t802 = t694 * qJDD(2);
t803 = qJDD(1) * t697;
t474 = -t694 * t588 + t553 - (t695 * t803 + t697 * t782 + t802) * pkin(9) + t608 * pkin(3) - t531 * qJD(3);
t733 = t694 * t563 + t697 * t588 + t604 * t815 - t633 * t816;
t481 = pkin(9) * t730 + t733;
t774 = -t693 * t474 - t883 * t481 - t524 * t784 + t531 * t814;
t859 = t683 * t695;
t719 = g(1) * t580 + g(2) * t578 + g(3) * t859 + t774;
t904 = t547 * t567 + t719;
t902 = pkin(4) * t741;
t900 = t496 * t741;
t875 = g(2) * t696;
t762 = g(1) * t699 + t875;
t899 = t698 * t762;
t898 = t741 * t871;
t528 = t883 * t531;
t493 = t693 * t530 + t528;
t766 = pkin(3) * t814 - t493;
t838 = qJ(5) * t824 - t910;
t566 = t693 * t634 - t635 * t883;
t897 = qJD(4) * t566 + t913 * t693 + t914 * t883;
t587 = t603 * qJ(5);
t626 = qJD(5) * t647;
t895 = t626 - t587;
t494 = t530 * t883 - t857;
t828 = -pkin(3) * t784 - qJD(5) + t494;
t894 = qJ(5) * t834 - qJD(5) * t614 + t767;
t887 = t486 * pkin(5) + t603 * qJ(6) - qJD(6) * t647;
t590 = t603 * pkin(4);
t577 = t682 * t846 + t683 * t699;
t579 = -t696 * t683 + t698 * t843;
t773 = -t883 * t474 + t693 * t481 + t524 * t814 + t531 * t784;
t861 = t682 * t695;
t718 = g(1) * t579 + g(2) * t577 + g(3) * t861 - t773;
t715 = -qJDD(5) + t718;
t713 = t590 + t715;
t708 = -t480 * t741 + t713 + t887;
t886 = -t567 * t741 + t718;
t885 = -0.2e1 * pkin(1);
t639 = t647 ^ 2;
t880 = pkin(7) * t694;
t878 = g(1) * t696;
t593 = t694 * t846 + t697 * t699;
t876 = g(2) * t593;
t874 = g(2) * t699;
t873 = g(3) * t695;
t872 = g(3) * t698;
t684 = t695 * pkin(7);
t870 = pkin(5) - t700;
t868 = qJ(5) * t487;
t866 = qJDD(2) * pkin(2);
t671 = pkin(3) * t693 + qJ(5);
t865 = t487 * t671;
t490 = t693 * t524 + t528;
t864 = t490 * t647;
t862 = t612 * t664;
t860 = t682 * t698;
t858 = t683 * t698;
t855 = t694 * t695;
t854 = t694 * t696;
t852 = t694 * t699;
t848 = t695 * t700;
t847 = t696 * t697;
t675 = pkin(3) * t697 + pkin(2);
t638 = t698 * t675;
t613 = -t795 + t856;
t842 = qJD(6) * t613 + t833 * t871 + t894;
t841 = -pkin(5) * t833 - t838;
t787 = t871 * t695;
t840 = -pkin(5) * t834 + qJD(1) * t787 + t897;
t839 = pkin(4) * t833 + t894;
t837 = pkin(4) * t824 + t897;
t610 = t697 * t627;
t556 = -pkin(9) * t850 + t610 + (-pkin(3) - t880) * t698;
t668 = pkin(7) * t845;
t826 = t694 * t627 + t668;
t564 = -pkin(9) * t855 + t826;
t835 = t693 * t556 + t883 * t564;
t832 = t694 * t618 + t627 * t815;
t820 = qJD(2) * t695;
t831 = t697 * t618 + t820 * t880;
t829 = t766 + t890;
t827 = -t828 + t901;
t665 = pkin(3) * t855;
t619 = t684 + t665;
t689 = t695 ^ 2;
t825 = -t698 ^ 2 + t689;
t822 = qJD(2) * t611;
t817 = qJD(3) * t611;
t812 = t612 * qJD(2);
t811 = t632 * qJD(3);
t810 = -qJD(5) - t489;
t808 = -t490 + t890;
t798 = t693 * t855;
t797 = t694 * t844;
t568 = pkin(3) * t908 + pkin(7) * t818;
t792 = t664 * t815;
t791 = t664 * t816;
t790 = t695 * t816;
t786 = t870 * t699;
t780 = -pkin(1) - t638;
t779 = -t577 * pkin(4) + qJ(5) * t578;
t778 = -t579 * pkin(4) + qJ(5) * t580;
t777 = -qJ(5) * t682 - t675;
t776 = -qJD(3) * t604 - t588;
t674 = -pkin(3) * t883 - pkin(4);
t771 = g(3) * (pkin(4) * t858 + qJ(5) * t860 + t638);
t770 = t883 * t818;
t636 = qJ(5) * t859;
t769 = -pkin(4) * t861 + t636;
t670 = g(1) * t851;
t768 = -g(2) * t849 + t670;
t764 = g(1) * t577 - g(2) * t579;
t763 = g(1) * t578 - g(2) * t580;
t676 = pkin(7) * t804;
t589 = pkin(7) * t782 + t676 - t866;
t484 = qJ(5) * t647 - t490;
t506 = qJ(5) * t698 - t835;
t583 = t695 * t795 - t798;
t760 = -qJ(5) * t583 + t619;
t758 = t556 * t883 - t693 * t564;
t757 = t633 * t815 - t553;
t756 = g(2) * t779;
t755 = pkin(3) * t612 + t867;
t754 = -pkin(8) * t608 + t811;
t752 = -g(3) * t860 + t682 * t907;
t751 = -g(3) * t858 + t683 * t907;
t749 = -qJ(5) * t614 - t675;
t507 = t698 * pkin(4) - t758;
t747 = -t682 * t787 + t636;
t458 = t774 + t895;
t744 = -pkin(7) * qJDD(2) + t807 * t885;
t743 = -qJ(6) * t577 + t779;
t742 = -qJ(6) * t579 + t778;
t565 = -t634 * t883 - t693 * t635;
t459 = -t590 + qJDD(5) + t773;
t740 = t608 * t694 - t792;
t739 = t697 * t608 + t791;
t737 = pkin(3) * t852 - t578 * pkin(4) + t699 * pkin(7) - qJ(5) * t577 + t696 * t848;
t510 = t750 * qJD(2) + (-t668 + (pkin(9) * t695 - t627) * t694) * qJD(3) + t831;
t513 = -t908 * pkin(9) + (-t695 * t819 - t698 * t816) * pkin(7) + t832;
t735 = -t510 * t883 + t693 * t513 + t556 * t814 + t564 * t784;
t734 = -t694 * t806 + t803;
t732 = t693 * t510 + t883 * t513 + t556 * t784 - t564 * t814;
t702 = qJD(1) ^ 2;
t731 = pkin(1) * t702 + t762;
t701 = qJD(2) ^ 2;
t729 = pkin(7) * t701 + qJDD(1) * t885 + t874;
t727 = t699 * pkin(1) + pkin(3) * t854 + t580 * pkin(4) + t696 * pkin(7) + qJ(5) * t579 + t675 * t844;
t726 = -t565 * t603 + t751;
t725 = t566 * t603 + t752;
t514 = t562 * t695 + t693 * t793 - t697 * t770;
t723 = qJ(5) * t514 - qJD(5) * t583 + t568;
t721 = -t873 - t899;
t720 = -qJD(3) * pkin(8) * t664 + t589 - t866 + t872;
t714 = -t719 - t895;
t464 = -qJ(5) * t820 + qJD(5) * t698 - t732;
t523 = -pkin(3) * t730 + t589;
t712 = -t715 + t900;
t711 = t603 * t671 + t714;
t709 = -t719 + t905;
t706 = t486 * qJ(5) - qJD(5) * t741 + t523;
t667 = pkin(3) * t847;
t662 = -qJ(6) + t674;
t596 = t697 * t844 + t854;
t595 = -t797 + t847;
t594 = -t696 * t845 + t852;
t582 = t614 * t695;
t544 = pkin(4) * t613 + t749;
t537 = -t613 * pkin(5) + t566;
t536 = t614 * pkin(5) + t565;
t529 = t613 * t871 + t749;
t525 = pkin(4) * t582 + t760;
t515 = t694 * t770 - t693 * t790 - qJD(4) * t798 + (t693 * t818 + t695 * t892) * t697;
t509 = t867 + t902;
t508 = t582 * t871 + t760;
t501 = t755 + t902;
t497 = -pkin(5) * t582 - t506;
t495 = t583 * pkin(5) + t698 * qJ(6) + t507;
t492 = t867 + t898;
t483 = pkin(4) * t647 - t810;
t482 = t755 + t898;
t469 = -t484 - t890;
t467 = t647 * t871 + t809;
t466 = pkin(4) * t515 + t723;
t465 = -pkin(4) * t820 + t735;
t463 = qJD(6) * t582 + t515 * t871 + t723;
t462 = -pkin(5) * t515 - t464;
t461 = -t514 * pkin(5) - qJD(2) * t787 + t698 * qJD(6) + t735;
t460 = t487 * pkin(4) + t706;
t457 = t547 * qJD(6) + t487 * t871 + t706;
t456 = -t458 + t888;
t455 = t459 - t887;
t1 = [(t460 * t525 + t496 * t466 + t458 * t506 + t484 * t464 + t459 * t507 + t483 * t465 - g(1) * (t696 * t780 + t737) - g(2) * (-t699 * t848 + t727)) * MDP(28) + (-t612 * t790 + (-t695 * t736 + t698 * t812 + t850 * t909) * t697) * MDP(11) + ((t611 * t697 - t612 * t694) * t818 + ((-t612 * qJD(3) + t730) * t697 + (-t817 - t707) * t694) * t695) * MDP(12) + (t457 * t508 + t480 * t463 + t455 * t495 + t467 * t461 + t456 * t497 + t469 * t462 - g(1) * (-qJ(6) * t578 + t737) - g(2) * (qJ(6) * t580 + t695 * t786 + t727) - (-pkin(5) * t695 + t780) * t878) * MDP(32) + (-t874 + t878) * MDP(2) + (t744 * t695 + (-t729 + t878) * t698) * MDP(9) + ((-t802 + (-t664 - t775) * t819) * t698 + (-t698 * t734 + t739 + t812) * t695) * MDP(13) + (-(-t627 * t816 + t831) * t664 + t610 * t608 - g(1) * t594 - g(2) * t596 + ((t792 - t822) * pkin(7) + (-pkin(7) * t608 + qJD(2) * t632 - t776) * t694 + t757) * t698 + (-pkin(7) * t730 + t557 * qJD(2) + t589 * t694 + t697 * t811) * t695) * MDP(16) + (t832 * t664 - t826 * t608 - g(1) * t593 - g(2) * t595 + (t632 * t819 + (-t791 + t812) * pkin(7) + t733) * t698 + (-t694 * t811 - t558 * qJD(2) + t589 * t697 + (t802 + t734 * t695 + (-t664 + t775) * t819) * pkin(7)) * t695) * MDP(17) + (t619 * t487 - t489 * t820 + t567 * t515 + t523 * t582 + t568 * t547 + t603 * t758 + t647 * t735 + t698 * t773 + t763) * MDP(23) + (-t608 * t698 - t664 * t820) * MDP(15) + (-t603 * t698 - t647 * t820) * MDP(22) + (t487 * t698 + t515 * t647 - t547 * t820 - t582 * t603) * MDP(21) + (t455 * t698 + t457 * t582 + t461 * t647 + t463 * t547 - t467 * t820 + t480 * t515 + t487 * t508 - t495 * t603 + t763) * MDP(31) + (-t459 * t698 - t460 * t582 - t465 * t647 - t466 * t547 + t483 * t820 - t487 * t525 - t496 * t515 + t507 * t603 - t763) * MDP(26) + (qJDD(1) * t689 + 0.2e1 * t695 * t782) * MDP(4) + (t695 * t729 + t698 * t744 - t670) * MDP(10) + ((t664 * t821 - t730) * t698 + (-t740 + t822) * t695) * MDP(14) + 0.2e1 * (t681 * t695 - t807 * t825) * MDP(5) + (qJDD(2) * t695 + t698 * t701) * MDP(6) + (qJDD(2) * t698 - t695 * t701) * MDP(7) + t762 * MDP(3) + (t458 * t582 + t459 * t583 + t464 * t547 + t465 * t741 - t483 * t514 + t484 * t515 - t486 * t507 + t487 * t506 + t768) * MDP(25) + (t455 * t583 - t456 * t582 + t461 * t741 - t462 * t547 - t467 * t514 - t469 * t515 - t486 * t495 - t487 * t497 + t768) * MDP(29) + (t458 * t698 - t460 * t583 + t464 * t647 - t466 * t741 - t484 * t820 + t486 * t525 + t496 * t514 - t506 * t603 + t764) * MDP(27) + (t486 * t698 + t514 * t647 + t583 * t603 + t741 * t820) * MDP(20) + (-t456 * t698 - t457 * t583 - t462 * t647 - t463 * t741 + t469 * t820 + t480 * t514 + t486 * t508 + t497 * t603 + t764) * MDP(30) + (-t619 * t486 - t490 * t820 - t567 * t514 + t523 * t583 + t568 * t741 - t603 * t835 + t647 * t732 - t698 * t774 - t764) * MDP(24) + (t486 * t582 - t487 * t583 + t514 * t547 - t515 * t741) * MDP(19) + (-t486 * t583 - t514 * t741) * MDP(18) + qJDD(1) * MDP(1); ((-t611 * t695 - t664 * t853) * qJD(1) + t739) * MDP(14) + (t675 * t486 + t523 * t614 - t834 * t567 + t647 * t910 + t767 * t741 - t725) * MDP(24) + (-t753 * t694 ^ 2 + (t694 * t722 - t862) * t697) * MDP(11) + ((-t772 + t862) * t694 + (t817 + 0.2e1 * t802 + t748 * t697 + (-t790 + (-t611 + t819) * t698) * qJD(1)) * t697) * MDP(12) + (-MDP(4) * t695 * t698 + MDP(5) * t825) * t702 + (t603 * t614 + t647 * t834) * MDP(20) + (-t603 * t613 + t647 * t833) * MDP(21) + (-pkin(2) * t772 + t830 * t664 + t754 * t694 + (-t557 * t695 + (pkin(7) * t611 - t632 * t694) * t698) * qJD(1) + (t695 * t762 - t720) * t697) * MDP(16) + (-t598 * t664 + (-t698 * pkin(7) * t612 + t558 * t695) * qJD(1) + (-pkin(2) * t748 + (t664 * t684 + (-t632 - t869) * t698) * qJD(1) + t754) * t697 + ((pkin(2) * t806 - t762) * t695 + t720) * t694) * MDP(17) + (t695 * t731 - t676 - t872) * MDP(9) + (t873 + (-pkin(7) * qJDD(1) + t731) * t698) * MDP(10) + ((-t612 * t695 + t664 * t845) * qJD(1) + t740) * MDP(13) + (-t460 * t613 - t487 * t544 - t496 * t833 - t547 * t839 - t647 * t837 - t726) * MDP(26) + (t457 * t613 + t480 * t833 + t487 * t529 - t536 * t603 + t547 * t842 + t647 * t840 + t751) * MDP(31) + MDP(7) * t681 + MDP(6) * t804 + (t664 * MDP(15) - MDP(20) * t741 + t547 * MDP(21) + t647 * MDP(22) + t489 * MDP(23) + t490 * MDP(24) - t483 * MDP(26) + t484 * MDP(27) - t469 * MDP(30) + t467 * MDP(31)) * t824 + (t458 * t613 + t459 * t614 - t483 * t834 + t484 * t833 - t486 * t565 - t487 * t566 + t547 * t838 + t741 * t837 + t721) * MDP(25) + (t486 * t613 - t487 * t614 + t547 * t834 - t741 * t833) * MDP(19) + (-t486 * t614 - t741 * t834) * MDP(18) + (t455 * t614 - t456 * t613 - t467 * t834 - t469 * t833 - t486 * t536 - t487 * t537 - t547 * t841 + t741 * t840 + t721) * MDP(29) + (-t460 * t614 + t486 * t544 + t496 * t834 + t647 * t838 - t741 * t839 + t725) * MDP(27) + (-t457 * t614 + t480 * t834 + t486 * t529 + t537 * t603 - t647 * t841 - t741 * t842 + t752) * MDP(30) + qJDD(2) * MDP(8) + (t457 * t529 + t455 * t536 + t456 * t537 - t771 + t842 * t480 + t841 * t469 + t840 * t467 + (-g(3) * qJ(6) * t683 - g(1) * t786 - t870 * t875) * t698 + (-g(3) * t870 + t762 * (t683 * t871 - t777)) * t695) * MDP(32) + (-t675 * t487 + t523 * t613 + t767 * t547 + t833 * t567 + t647 * t897 + t726) * MDP(23) + (t460 * t544 - t458 * t566 + t459 * t565 - t771 + t700 * t899 + t839 * t496 + t838 * t484 + t837 * t483 + (g(3) * t700 + t762 * (pkin(4) * t683 - t777)) * t695) * MDP(28); (-t494 * t647 + (-t603 * t693 - t612 * t741 + t647 * t784) * pkin(3) + t904) * MDP(24) + (-t486 * t674 - t865 + (-t484 + t766) * t741 + (t483 + t828) * t547) * MDP(25) + (t501 * t741 + t647 * t828 + t711 - t911) * MDP(27) + (-t486 * t662 - t865 + (t469 + t829) * t741 + (t467 - t827) * t547) * MDP(29) + (t482 * t741 - t647 * t827 + t711 + t905) * MDP(30) + (g(1) * t596 - g(2) * t594 + g(3) * t850 - t557 * t664 - t611 * t632 - t733) * MDP(17) + (t730 - t862) * MDP(14) - t612 * t611 * MDP(11) + (-g(1) * t595 + t876 - t558 * t664 - t612 * t632 + (t776 + t873) * t694 - t757) * MDP(16) + (-t458 * t671 + t459 * t674 - t496 * t501 - t483 * t493 - g(1) * (t667 + t778) - t756 - g(3) * (-t665 + t769) + t828 * t484 + (g(1) * t797 + t483 * t814 + t876) * pkin(3)) * MDP(28) + (t455 * t662 + t456 * t671 - t480 * t482 - g(1) * (-pkin(3) * t797 + t667 + t742) - g(2) * (-pkin(3) * t593 + t743) - g(3) * (-t665 + t747) + t827 * t469 + t829 * t467) * MDP(32) + (-t611 ^ 2 + t612 ^ 2) * MDP(12) + t608 * MDP(15) + (t611 * t664 + t707) * MDP(13) + (t501 * t547 + t603 * t674 - t647 * t766 - t713 + t900) * MDP(26) + (-t493 * t647 + (-t547 * t612 + t603 * t883 + t647 * t814) * pkin(3) + t886) * MDP(23) + (-t482 * t547 - t603 * t662 + t647 * t829 + t708) * MDP(31) + t912; (-t864 + t886) * MDP(23) + (t489 * t647 + t904) * MDP(24) + (pkin(4) * t486 - t868 + (-t484 - t490) * t741 + (t483 + t810) * t547) * MDP(25) + (t509 * t547 - 0.2e1 * t590 + t712 + t864) * MDP(26) + (t509 * t741 + t647 * t810 + t587 + t714 - t911) * MDP(27) + (-t459 * pkin(4) - g(1) * t778 - g(3) * t769 - t458 * qJ(5) - t483 * t490 + t484 * t810 - t496 * t509 - t756) * MDP(28) + (-t868 + t486 * t871 + (t469 + t808) * t741 + (t467 - t809) * t547) * MDP(29) + (t492 * t741 - t647 * t745 + 0.2e1 * t587 - 0.2e1 * t626 + t709) * MDP(30) + (-t492 * t547 + t603 * t871 + t647 * t808 + t708) * MDP(31) + (-g(1) * t742 - g(2) * t743 - g(3) * t747 + t456 * qJ(5) - t455 * t871 + t467 * t808 + t469 * t809 - t480 * t492) * MDP(32) + t912; (-t484 * t647 - t590 + t712) * MDP(28) + (t469 * t647 - t708) * MDP(32) + (MDP(26) - MDP(31)) * (t603 - t863) + (MDP(25) + MDP(29)) * t468 + (MDP(27) + MDP(30)) * (-t545 - t639); t703 * MDP(29) + (t603 + t863) * MDP(30) + (-t639 - t884) * MDP(31) + (-t467 * t647 + t709 - t895) * MDP(32);];
tau  = t1;
