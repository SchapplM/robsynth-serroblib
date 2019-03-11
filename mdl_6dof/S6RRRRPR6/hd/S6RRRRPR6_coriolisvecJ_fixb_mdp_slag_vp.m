% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRPR6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:23:36
% EndTime: 2019-03-09 22:23:56
% DurationCPUTime: 13.35s
% Computational Cost: add. (13052->592), mult. (32383->801), div. (0->0), fcn. (23927->10), ass. (0->269)
t728 = cos(qJ(3));
t725 = sin(qJ(2));
t729 = cos(qJ(2));
t829 = t728 * t729;
t745 = pkin(3) * t725 - pkin(9) * t829;
t852 = pkin(8) + pkin(9);
t788 = qJD(3) * t852;
t760 = pkin(2) * t725 - pkin(8) * t729;
t679 = t760 * qJD(1);
t724 = sin(qJ(3));
t802 = qJD(1) * t725;
t787 = t724 * t802;
t808 = pkin(7) * t787 + t728 * t679;
t890 = qJD(1) * t745 + t728 * t788 + t808;
t658 = t724 * t679;
t833 = t725 * t728;
t834 = t724 * t729;
t887 = t658 + (-pkin(7) * t833 - pkin(9) * t834) * qJD(1) + t724 * t788;
t799 = qJD(2) * t725;
t776 = MDP(31) * t799;
t722 = sin(qJ(6));
t726 = cos(qJ(6));
t792 = t728 * qJD(2);
t671 = t787 - t792;
t800 = qJD(2) * t724;
t673 = t728 * t802 + t800;
t723 = sin(qJ(4));
t727 = cos(qJ(4));
t622 = t727 * t671 + t673 * t723;
t720 = sin(pkin(11));
t721 = cos(pkin(11));
t749 = t671 * t723 - t727 * t673;
t861 = t622 * t720 + t721 * t749;
t862 = -t622 * t721 + t720 * t749;
t873 = t722 * t862 - t726 * t861;
t830 = t726 * t862;
t846 = t861 * t722;
t880 = t830 + t846;
t889 = qJD(1) * t776 + (t873 ^ 2 - t880 ^ 2) * MDP(28) - t880 * MDP(27) * t873;
t791 = qJD(1) * qJD(2);
t780 = t729 * t791;
t797 = qJD(3) * t724;
t784 = t725 * t797;
t790 = qJD(2) * qJD(3);
t637 = -qJD(1) * t784 + (t780 + t790) * t728;
t796 = qJD(3) * t728;
t782 = t725 * t796;
t798 = qJD(2) * t729;
t786 = t724 * t798;
t737 = t782 + t786;
t638 = qJD(1) * t737 + t724 * t790;
t794 = qJD(4) * t727;
t795 = qJD(4) * t723;
t556 = t727 * t637 - t723 * t638 - t671 * t794 - t673 * t795;
t733 = qJD(4) * t749 - t637 * t723 - t727 * t638;
t499 = -t556 * t720 + t721 * t733;
t500 = t556 * t721 + t720 * t733;
t789 = qJD(6) * t830 + t722 * t499 + t726 * t500;
t793 = qJD(6) * t722;
t455 = t793 * t861 + t789;
t773 = -t726 * t499 + t500 * t722;
t456 = qJD(6) * t873 + t773;
t801 = qJD(1) * t729;
t704 = -qJD(3) + t801;
t695 = -qJD(4) + t704;
t778 = MDP(22) * t802;
t687 = -qJD(6) + t695;
t847 = t880 * t687;
t848 = t873 * t687;
t888 = qJD(2) * t778 + (-t622 ^ 2 + t749 ^ 2) * MDP(19) + (-t622 * t695 + t556) * MDP(20) + (t695 * t749 + t733) * MDP(21) - t622 * MDP(18) * t749 + (-t456 - t848) * MDP(30) + t889 + (t455 + t847) * MDP(29);
t674 = t723 * t724 - t727 * t728;
t740 = t674 * t729;
t853 = qJD(3) + qJD(4);
t882 = -qJD(1) * t740 + t853 * t674;
t675 = t723 * t728 + t724 * t727;
t815 = (-t801 + t853) * t675;
t684 = -pkin(2) * t729 - pkin(8) * t725 - pkin(1);
t664 = t684 * qJD(1);
t713 = pkin(7) * t801;
t691 = qJD(2) * pkin(8) + t713;
t628 = t728 * t664 - t691 * t724;
t601 = -pkin(9) * t673 + t628;
t592 = -pkin(3) * t704 + t601;
t836 = t724 * t664;
t629 = t691 * t728 + t836;
t602 = -pkin(9) * t671 + t629;
t599 = t727 * t602;
t540 = t592 * t723 + t599;
t682 = t760 * qJD(2);
t665 = qJD(1) * t682;
t781 = t725 * t791;
t764 = pkin(7) * t781;
t814 = -t728 * t665 - t724 * t764;
t736 = -qJD(3) * t629 - t814;
t550 = pkin(3) * t781 - pkin(9) * t637 + t736;
t744 = t664 * t796 + t724 * t665 - t691 * t797;
t732 = -t728 * t764 + t744;
t562 = -pkin(9) * t638 + t732;
t772 = t727 * t550 - t723 * t562;
t734 = -qJD(4) * t540 + t772;
t459 = pkin(4) * t781 - qJ(5) * t556 + qJD(5) * t749 + t734;
t763 = -t723 * t550 - t727 * t562 - t592 * t794 + t602 * t795;
t461 = qJ(5) * t733 - qJD(5) * t622 - t763;
t447 = t721 * t459 - t461 * t720;
t445 = pkin(5) * t781 - pkin(10) * t500 + t447;
t448 = t720 * t459 + t721 * t461;
t446 = pkin(10) * t499 + t448;
t597 = t723 * t602;
t539 = t727 * t592 - t597;
t868 = qJ(5) * t749;
t517 = t539 + t868;
t512 = -pkin(4) * t695 + t517;
t869 = qJ(5) * t622;
t518 = t540 - t869;
t838 = t721 * t518;
t473 = t720 * t512 + t838;
t858 = pkin(10) * t862;
t464 = t473 + t858;
t463 = t464 * t793;
t690 = -qJD(2) * pkin(2) + pkin(7) * t802;
t639 = pkin(3) * t671 + t690;
t588 = pkin(4) * t622 + qJD(5) + t639;
t522 = -pkin(5) * t862 + t588;
t886 = -t722 * t445 - t726 * t446 - t522 * t880 + t463;
t883 = t890 * t727;
t692 = t852 * t724;
t693 = t852 * t728;
t881 = -t692 * t794 - t693 * t795 - t723 * t890 - t887 * t727;
t774 = t726 * t445 - t722 * t446;
t874 = -t522 * t873 + t774;
t878 = pkin(5) * t861;
t870 = pkin(10) * t861;
t809 = -t723 * t692 + t727 * t693;
t877 = -pkin(4) * t802 + qJ(5) * t882 - qJD(4) * t809 - qJD(5) * t675 + t723 * t887 - t883;
t876 = -qJ(5) * t815 - qJD(5) * t674 + t881;
t761 = -t713 + (-t724 * t801 + t797) * pkin(3);
t875 = pkin(4) * t815 + t761;
t872 = pkin(4) * t749;
t867 = t720 * t882 - t721 * t815;
t866 = -t720 * t815 - t721 * t882;
t864 = t622 * t639 + t763;
t863 = t639 * t749 + t734;
t859 = -0.2e1 * t791;
t857 = MDP(4) * t725;
t718 = t725 ^ 2;
t856 = MDP(5) * (-t729 ^ 2 + t718);
t650 = t675 * t725;
t822 = -t720 * t876 + t721 * t877;
t821 = t720 * t877 + t721 * t876;
t670 = t728 * t684;
t850 = pkin(7) * t724;
t627 = -pkin(9) * t833 + t670 + (-pkin(3) - t850) * t729;
t706 = pkin(7) * t829;
t806 = t724 * t684 + t706;
t835 = t724 * t725;
t633 = -pkin(9) * t835 + t806;
t817 = t723 * t627 + t727 * t633;
t770 = -t601 * t723 - t599;
t523 = t770 + t869;
t819 = t727 * t601 - t597;
t524 = t819 + t868;
t837 = t721 * t723;
t849 = pkin(3) * qJD(4);
t813 = -t721 * t523 + t524 * t720 + (-t720 * t727 - t837) * t849;
t839 = t720 * t723;
t812 = -t720 * t523 - t721 * t524 + (t721 * t727 - t839) * t849;
t854 = t729 * t792 - t784;
t851 = pkin(4) * t720;
t845 = t637 * t724;
t844 = t671 * t704;
t843 = t673 * t704;
t842 = t690 * t724;
t841 = t690 * t728;
t840 = t704 * t728;
t513 = t720 * t518;
t730 = qJD(2) ^ 2;
t832 = t725 * t730;
t472 = t721 * t512 - t513;
t462 = -pkin(5) * t695 + t472 + t870;
t831 = t726 * t462;
t828 = t729 * t730;
t731 = qJD(1) ^ 2;
t827 = t729 * t731;
t826 = -t813 - t858;
t825 = t812 - t870;
t586 = -qJD(2) * t740 - t650 * t853;
t651 = t674 * t725;
t810 = t728 * t682 + t799 * t850;
t581 = t745 * qJD(2) + (-t706 + (pkin(9) * t725 - t684) * t724) * qJD(3) + t810;
t783 = t729 * t797;
t811 = t724 * t682 + t684 * t796;
t585 = -t737 * pkin(9) + (-t725 * t792 - t783) * pkin(7) + t811;
t771 = t727 * t581 - t585 * t723;
t476 = pkin(4) * t799 - qJ(5) * t586 - qJD(4) * t817 + qJD(5) * t651 + t771;
t587 = -t795 * t835 + (t833 * t853 + t786) * t727 + t854 * t723;
t739 = t723 * t581 + t727 * t585 + t627 * t794 - t633 * t795;
t480 = -qJ(5) * t587 - qJD(5) * t650 + t739;
t454 = t720 * t476 + t721 * t480;
t619 = -t674 * t721 - t675 * t720;
t620 = -t674 * t720 + t675 * t721;
t754 = t726 * t619 - t620 * t722;
t824 = qJD(6) * t754 + t722 * t867 + t726 * t866;
t577 = t619 * t722 + t620 * t726;
t823 = qJD(6) * t577 + t722 * t866 - t726 * t867;
t479 = t721 * t517 - t513;
t820 = -pkin(5) * t867 + t875;
t769 = t727 * t627 - t633 * t723;
t558 = -pkin(4) * t729 + qJ(5) * t651 + t769;
t564 = -qJ(5) * t650 + t817;
t502 = t720 * t558 + t721 * t564;
t767 = -t727 * t692 - t693 * t723;
t610 = -qJ(5) * t675 + t767;
t611 = -qJ(5) * t674 + t809;
t561 = t720 * t610 + t721 * t611;
t683 = pkin(3) * t835 + t725 * pkin(7);
t640 = pkin(3) * t737 + pkin(7) * t798;
t711 = -pkin(3) * t728 - pkin(2);
t779 = MDP(15) * t802;
t618 = pkin(3) * t638 + pkin(7) * t780;
t775 = pkin(1) * t859;
t453 = t721 * t476 - t480 * t720;
t478 = -t517 * t720 - t838;
t501 = t721 * t558 - t564 * t720;
t560 = t721 * t610 - t611 * t720;
t766 = t671 + t792;
t765 = -t673 + t800;
t710 = pkin(3) * t727 + pkin(4);
t653 = -pkin(3) * t839 + t721 * t710;
t762 = pkin(4) * t650 + t683;
t759 = pkin(3) * t673 - t872;
t531 = pkin(10) * t619 + t561;
t758 = pkin(5) * t802 + pkin(10) * t866 + qJD(6) * t531 - t822;
t530 = -pkin(10) * t620 + t560;
t757 = pkin(10) * t867 + qJD(6) * t530 + t821;
t450 = t722 * t462 + t726 * t464;
t604 = -t650 * t720 - t651 * t721;
t487 = -pkin(5) * t729 - pkin(10) * t604 + t501;
t603 = -t650 * t721 + t651 * t720;
t488 = pkin(10) * t603 + t502;
t756 = t487 * t722 + t488 * t726;
t755 = t726 * t603 - t604 * t722;
t552 = t603 * t722 + t604 * t726;
t644 = pkin(5) + t653;
t654 = pkin(3) * t837 + t710 * t720;
t751 = t644 * t726 - t654 * t722;
t750 = t644 * t722 + t654 * t726;
t748 = qJD(1) * t718 - t704 * t729;
t747 = pkin(4) * t674 + t711;
t746 = pkin(4) * t587 + t640;
t528 = -pkin(4) * t733 + t618;
t707 = pkin(4) * t721 + pkin(5);
t743 = t707 * t722 + t726 * t851;
t742 = t707 * t726 - t722 * t851;
t591 = -pkin(5) * t619 + t747;
t569 = -pkin(5) * t603 + t762;
t532 = -t872 - t878;
t529 = t759 - t878;
t527 = t586 * t721 - t587 * t720;
t526 = -t586 * t720 - t587 * t721;
t493 = -pkin(5) * t526 + t746;
t471 = -pkin(5) * t499 + t528;
t468 = qJD(6) * t552 - t726 * t526 + t527 * t722;
t467 = qJD(6) * t755 + t526 * t722 + t527 * t726;
t466 = t479 + t870;
t465 = t478 - t858;
t452 = pkin(10) * t526 + t454;
t451 = pkin(5) * t799 - pkin(10) * t527 + t453;
t449 = -t464 * t722 + t831;
t1 = [(-pkin(7) * t828 + t725 * t775) * MDP(9) + 0.2e1 * t780 * t857 + ((-t671 * t728 - t673 * t724) * t798 + (-t845 - t638 * t728 + (t671 * t724 - t673 * t728) * qJD(3)) * t725) * MDP(12) + ((-pkin(7) * t783 + t811) * t704 + t744 * t729 + (pkin(7) * t637 - t690 * t797) * t725 + ((pkin(7) * t673 + t841) * t729 + (-pkin(7) * t840 - qJD(1) * t806 - t629) * t725) * qJD(2)) * MDP(17) + (-(-t684 * t797 + t810) * t704 + (t690 * t796 + pkin(7) * t638 + (qJD(1) * t670 + t628) * qJD(2)) * t725 + ((pkin(7) * t671 + t842) * qJD(2) + (t836 + (pkin(7) * t704 + t691) * t728) * qJD(3) + t814) * t729) * MDP(16) + (t456 * t729 + t468 * t687) * MDP(30) + (-t455 * t729 - t467 * t687) * MDP(29) + (-t556 * t729 - t586 * t695) * MDP(20) + (t704 * t782 + t638 * t729 + (-t671 * t725 - t724 * t748) * qJD(2)) * MDP(14) + (-t771 * t695 - t772 * t729 + t640 * t622 - t683 * t733 + t618 * t650 + t639 * t587 + (t540 * t729 + t695 * t817) * qJD(4)) * MDP(23) + (-t447 * t604 + t448 * t603 + t453 * t861 + t454 * t862 - t472 * t527 + t473 * t526 + t499 * t502 - t500 * t501) * MDP(25) + (t587 * t695 - t729 * t733) * MDP(21) + (-t556 * t650 - t586 * t622 + t587 * t749 - t651 * t733) * MDP(19) + t856 * t859 + (t569 * t455 - t463 * t729 + t522 * t467 + t471 * t552 + t493 * t873 + ((-qJD(6) * t488 + t451) * t687 + t445 * t729) * t722 + ((qJD(6) * t487 + t452) * t687 + (qJD(6) * t462 + t446) * t729) * t726) * MDP(33) + (t455 * t552 + t467 * t873) * MDP(27) + (t447 * t501 + t448 * t502 + t472 * t453 + t473 * t454 + t528 * t762 + t588 * t746) * MDP(26) + (-t556 * t651 - t586 * t749) * MDP(18) + (t683 * t556 + t639 * t586 - t618 * t651 - t640 * t749 + t739 * t695 - t763 * t729) * MDP(24) + (t455 * t755 - t456 * t552 + t467 * t880 - t468 * t873) * MDP(28) + (((t487 * t726 - t488 * t722) * qJD(1) + t449) * MDP(32) + (-qJD(1) * t651 - t749) * MDP(20) + (-qJD(1) * t650 - t622) * MDP(21) + (qJD(1) * t769 + t539) * MDP(23) + (qJD(1) * t552 + t873) * MDP(29) + (qJD(1) * t755 + t880) * MDP(30) + (-qJD(1) * t756 - t450) * MDP(33) + (-qJD(1) * t817 - t540) * MDP(24) + (-t704 - t801) * MDP(15) + (-t695 - t801) * MDP(22)) * t799 + (-(t451 * t726 - t452 * t722) * t687 - t774 * t729 - t493 * t880 + t569 * t456 - t471 * t755 + t522 * t468 + (t450 * t729 + t687 * t756) * qJD(6)) * MDP(32) + (t704 * t784 - t637 * t729 + (t673 * t725 + t728 * t748) * qJD(2)) * MDP(13) + (-t687 - t801) * t776 - MDP(7) * t832 + (pkin(7) * t832 + t729 * t775) * MDP(10) + MDP(6) * t828 + (t637 * t833 + t673 * t854) * MDP(11); (-t704 * t796 + (t704 * t829 + t725 * t765) * qJD(1)) * MDP(13) + t731 * t856 - t827 * t857 + ((t637 + t844) * t728 + (-t638 + t843) * t724) * MDP(12) + (-t673 * t840 + t845) * MDP(11) + (-pkin(2) * t637 - t658 * t704 + (-pkin(8) * t704 * t724 + t841) * qJD(3) + (-t690 * t829 + (-pkin(8) * t792 + t629) * t725 + (t704 * t833 + t729 * t765) * pkin(7)) * qJD(1)) * MDP(17) + (-pkin(2) * t638 + t808 * t704 + (pkin(8) * t840 + t842) * qJD(3) + ((-pkin(8) * t800 - t628) * t725 + (-pkin(7) * t766 - t842) * t729) * qJD(1)) * MDP(16) + (-t447 * t620 + t448 * t619 - t472 * t866 + t473 * t867 + t499 * t561 - t500 * t560 + t821 * t862 + t822 * t861) * MDP(25) + (-t824 * t687 + (qJD(2) * t577 - t873) * t802) * MDP(29) + (t591 * t455 + t471 * t577 + (-t722 * t758 + t726 * t757) * t687 + t824 * t522 + t820 * t873 + (-(t530 * t722 + t531 * t726) * qJD(2) + t450) * t802) * MDP(33) + (t455 * t577 + t824 * t873) * MDP(27) + (t447 * t560 + t448 * t561 + t822 * t472 + t821 * t473 + t528 * t747 + t588 * t875) * MDP(26) + (t556 * t675 + t749 * t882) * MDP(18) + (t882 * t695 + (qJD(2) * t675 + t749) * t802) * MDP(20) + (t711 * t556 + t618 * t675 + t881 * t695 - t882 * t639 - t761 * t749 + (-qJD(2) * t809 + t540) * t802) * MDP(24) + (-t556 * t674 + t622 * t882 + t675 * t733 + t749 * t815) * MDP(19) + t687 * MDP(31) * t802 + (t455 * t754 - t456 * t577 - t823 * t873 + t824 * t880) * MDP(28) + (t591 * t456 - t471 * t754 + (t722 * t757 + t726 * t758) * t687 + t823 * t522 - t820 * t880 + ((t530 * t726 - t531 * t722) * qJD(2) - t449) * t802) * MDP(32) + (t823 * t687 + (qJD(2) * t754 - t880) * t802) * MDP(30) + (t815 * t695 + (-qJD(2) * t674 + t622) * t802) * MDP(21) + (t704 * t797 + (-t704 * t834 + t725 * t766) * qJD(1)) * MDP(14) + (MDP(9) * t725 * t731 + MDP(10) * t827) * pkin(1) + t695 * t778 + t704 * t779 + (-t711 * t733 + t618 * t674 + (t693 * t794 + (-qJD(4) * t692 - t887) * t723 + t883) * t695 + t815 * t639 + t761 * t622 + (qJD(2) * t767 - t539) * t802) * MDP(23); (t751 * t781 + t529 * t880 + (t722 * t825 + t726 * t826) * t687 + (t687 * t750 - t450) * qJD(6) + t874) * MDP(32) + (t637 - t844) * MDP(13) + (-t750 * t781 - t529 * t873 + (-t722 * t826 + t726 * t825) * t687 + (t687 * t751 - t831) * qJD(6) + t886) * MDP(33) + (t499 * t654 - t500 * t653 + (t472 + t812) * t862 + (-t473 + t813) * t861) * MDP(25) + (t770 * t695 + (-t622 * t673 + t695 * t795 + t727 * t781) * pkin(3) + t863) * MDP(23) + (-t819 * t695 + (t673 * t749 + t695 * t794 - t723 * t781) * pkin(3) + t864) * MDP(24) + (-t628 * t704 + t671 * t690 - t732) * MDP(17) + (t447 * t653 + t448 * t654 + t472 * t813 + t473 * t812 - t588 * t759) * MDP(26) + (-t629 * t704 - t673 * t690 + t736) * MDP(16) + t673 * t671 * MDP(11) + qJD(2) * t779 + (-t671 ^ 2 + t673 ^ 2) * MDP(12) + (-t638 - t843) * MDP(14) + t888; (-t540 * t695 + t863) * MDP(23) + (-t539 * t695 + t864) * MDP(24) + ((t499 * t720 - t500 * t721) * pkin(4) + (t472 - t479) * t862 + (-t473 - t478) * t861) * MDP(25) + (-t472 * t478 - t473 * t479 + (t447 * t721 + t448 * t720 + t588 * t749) * pkin(4)) * MDP(26) + (t742 * t781 + (t465 * t726 - t466 * t722) * t687 + t532 * t880 + (t687 * t743 - t450) * qJD(6) + t874) * MDP(32) + (-t743 * t781 - (t465 * t722 + t466 * t726) * t687 - t532 * t873 + (t687 * t742 - t831) * qJD(6) + t886) * MDP(33) + t888; (-t861 ^ 2 - t862 ^ 2) * MDP(25) + (-t472 * t861 - t473 * t862 + t528) * MDP(26) + (t456 - t848) * MDP(32) + (t455 - t847) * MDP(33); (t789 + t847) * MDP(29) + (-t773 - t848) * MDP(30) + (-t450 * t687 + t874) * MDP(32) + (-t449 * t687 + t886) * MDP(33) + (MDP(29) * t846 - MDP(30) * t873 - MDP(32) * t450 - MDP(33) * t831) * qJD(6) + t889;];
tauc  = t1;
