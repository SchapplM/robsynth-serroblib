% Calculate vector of inverse dynamics joint torques for
% S6RRRRPP3
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
%   see S6RRRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRRPP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:56:51
% EndTime: 2019-03-09 20:57:06
% DurationCPUTime: 11.97s
% Computational Cost: add. (10019->697), mult. (21810->793), div. (0->0), fcn. (15371->10), ass. (0->300)
t701 = sin(qJ(2));
t704 = cos(qJ(2));
t884 = sin(qJ(3));
t793 = qJDD(1) * t884;
t885 = cos(qJ(3));
t794 = qJDD(1) * t885;
t758 = t701 * t793 - t704 * t794;
t742 = -t701 * t885 - t704 * t884;
t813 = qJD(2) + qJD(3);
t573 = t813 * t742;
t888 = qJD(1) * t573;
t538 = t758 - t888;
t703 = cos(qJ(4));
t820 = qJD(4) * t703;
t801 = qJD(1) * t884;
t802 = qJD(1) * t885;
t616 = t701 * t801 - t704 * t802;
t851 = t616 * t703;
t905 = t820 + t851;
t700 = sin(qJ(4));
t821 = qJD(4) * t700;
t852 = t616 * t700;
t904 = t821 + t852;
t537 = qJDD(4) + t538;
t871 = pkin(4) + qJ(6);
t798 = t871 * t537;
t698 = qJ(2) + qJ(3);
t692 = sin(t698);
t705 = cos(qJ(1));
t845 = t692 * t705;
t702 = sin(qJ(1));
t847 = t692 * t702;
t906 = -g(1) * t845 - g(2) * t847;
t706 = -pkin(8) - pkin(7);
t648 = t706 * t704;
t635 = qJD(1) * t648;
t622 = t885 * t635;
t646 = t706 * t701;
t633 = qJD(1) * t646;
t569 = t884 * t633 - t622;
t799 = qJD(3) * t884;
t764 = pkin(2) * t799 - t569;
t630 = t701 * t884 - t704 * t885;
t572 = t813 * t630;
t858 = t572 * t700;
t751 = -t742 * t820 - t858;
t618 = -t701 * t802 - t704 * t801;
t745 = t703 * t618 - t700 * t813;
t903 = qJD(4) * t745;
t872 = t704 * pkin(2);
t688 = pkin(1) + t872;
t644 = t688 * qJD(1);
t549 = pkin(3) * t616 + pkin(9) * t618 - t644;
t870 = qJD(2) * pkin(2);
t624 = t633 + t870;
t566 = t884 * t624 - t622;
t553 = pkin(9) * t813 + t566;
t517 = -t703 * t549 + t553 * t700;
t759 = -pkin(5) * t745 + t517;
t819 = qJD(5) + t759;
t693 = cos(t698);
t902 = g(3) * t693 + t906;
t901 = pkin(9) * t742 - t688;
t610 = qJD(4) + t616;
t859 = t537 * t703;
t900 = pkin(9) * (t610 * t821 - t859);
t860 = t537 * t700;
t899 = pkin(9) * (t610 * t820 + t860);
t535 = t537 * qJ(5);
t595 = t610 * qJD(5);
t898 = -t535 - t595;
t786 = pkin(4) * t821 - qJD(5) * t700;
t897 = (-qJ(5) * qJD(4) - qJD(6)) * t703 + t786 + t904 * qJ(6);
t585 = -t885 * t646 - t884 * t648;
t895 = -t693 * pkin(3) - t692 * pkin(9);
t894 = pkin(5) * t905 - t618 * t871;
t773 = -pkin(4) * t852 + qJ(5) * t851;
t893 = -t773 + t764;
t713 = t572 * qJD(1);
t710 = -t742 * qJDD(1) - t713;
t818 = qJD(5) + t517;
t500 = -pkin(4) * t610 + t818;
t518 = t700 * t549 + t703 * t553;
t501 = -qJ(5) * t610 - t518;
t892 = -t500 * t700 + t501 * t703;
t775 = g(1) * t705 + g(2) * t702;
t800 = qJD(3) * t885;
t812 = qJDD(2) + qJDD(3);
t674 = t703 * t812;
t516 = t700 * t710 - t674 - t903;
t890 = -pkin(5) * t516 + qJDD(6);
t889 = t692 * t775;
t581 = -t618 * t700 - t703 * t813;
t887 = t581 ^ 2;
t579 = t745 ^ 2;
t603 = t610 ^ 2;
t886 = 0.2e1 * t535;
t883 = pkin(2) * t701;
t882 = pkin(4) * t537;
t881 = pkin(4) * t618;
t785 = qJD(4) * t813;
t515 = -t618 * t821 - t700 * t812 + (-t710 - t785) * t703;
t880 = pkin(5) * t515;
t878 = pkin(5) * t581;
t875 = g(2) * t706;
t873 = g(3) * t704;
t678 = t692 * pkin(5);
t869 = qJ(5) * t516;
t868 = qJ(5) * t581;
t816 = qJD(1) * qJD(2);
t797 = t701 * t816;
t675 = pkin(2) * t797;
t490 = t538 * pkin(3) + pkin(9) * t713 + qJDD(1) * t901 + t675;
t796 = t704 * t816;
t815 = qJDD(1) * t701;
t575 = qJDD(2) * pkin(2) - t706 * (-t796 - t815);
t814 = qJDD(1) * t704;
t580 = t706 * (-t797 + t814);
t722 = t575 * t884 - t580 * t885 + t624 * t800 + t635 * t799;
t505 = pkin(9) * t812 + t722;
t783 = -t700 * t490 - t703 * t505 - t549 * t820 + t553 * t821;
t463 = t783 + t898;
t867 = t463 * t703;
t782 = -t703 * t490 + t700 * t505 + t549 * t821 + t553 * t820;
t766 = -qJDD(5) - t782;
t465 = -t766 - t882;
t464 = t465 * t700;
t864 = t515 * t700;
t863 = t516 * t703;
t862 = t518 * t610;
t686 = pkin(2) * t884 + pkin(9);
t861 = t537 * t686;
t857 = t581 * t745;
t856 = t581 * t610;
t855 = t581 * t700;
t854 = t745 * t610;
t853 = t745 * t703;
t850 = t742 * t700;
t849 = t742 * t703;
t848 = t692 * t700;
t846 = t692 * t703;
t844 = t693 * t702;
t843 = t693 * t703;
t842 = t693 * t705;
t841 = t700 * qJ(5);
t840 = t700 * t702;
t839 = t702 * t703;
t838 = t703 * t705;
t837 = t705 * t700;
t560 = -pkin(3) * t618 + pkin(9) * t616;
t621 = t884 * t635;
t565 = t624 * t885 + t621;
t836 = t700 * t560 + t703 * t565;
t822 = qJD(1) * t701;
t551 = pkin(2) * t822 + t560;
t570 = t633 * t885 + t621;
t835 = t700 * t551 + t703 * t570;
t834 = t893 + t897;
t564 = pkin(3) * t630 + t901;
t586 = t646 * t884 - t648 * t885;
t833 = t700 * t564 + t703 * t586;
t523 = t773 + t566;
t832 = -t523 + t897;
t784 = pkin(2) * t800;
t770 = t703 * t784;
t604 = t618 * qJ(5);
t779 = pkin(5) * t852 - t604;
t831 = t770 + (-pkin(5) - t686) * t821 - t779 - t835;
t771 = t700 * t784;
t736 = t686 * t820 + t771;
t561 = t700 * t570;
t790 = -t551 * t703 + t561;
t830 = t736 - t790 + t894;
t609 = -qJ(5) * t820 + t786;
t829 = t609 + t893;
t828 = t609 - t523;
t827 = (-pkin(5) - pkin(9)) * t821 - t779 - t836;
t789 = t560 * t703 - t700 * t565;
t826 = pkin(9) * t820 + t789 + t894;
t653 = pkin(9) * t844;
t825 = pkin(5) * t844 + t653;
t657 = pkin(9) * t842;
t824 = pkin(5) * t842 + t657;
t696 = t701 ^ 2;
t823 = -t704 ^ 2 + t696;
t492 = t518 - t878;
t817 = -qJD(6) - t492;
t811 = t885 * pkin(2);
t810 = t701 * t870;
t809 = t902 * t700;
t808 = g(3) * t843 + t703 * t906;
t807 = g(1) * t842 + g(2) * t844 + g(3) * t692;
t806 = qJD(2) * t706;
t804 = t742 * t821;
t795 = -pkin(3) - t841;
t605 = t693 * t840 + t838;
t606 = t693 * t839 - t837;
t792 = -t605 * pkin(4) + qJ(5) * t606;
t607 = t693 * t837 - t839;
t608 = t693 * t838 + t840;
t791 = -t607 * pkin(4) + qJ(5) * t608;
t576 = t700 * t586;
t788 = t564 * t703 - t576;
t787 = t610 * t703;
t781 = -t885 * t575 - t884 * t580 + t624 * t799 - t635 * t800;
t780 = pkin(4) * t843 + t693 * t841 - t895;
t778 = g(1) * t847 - g(2) * t845;
t777 = g(1) * t605 - g(2) * t607;
t776 = g(1) * t606 - g(2) * t608;
t774 = g(1) * t702 - g(2) * t705;
t520 = -qJ(5) * t630 - t833;
t769 = -qJ(5) * t703 + qJ(6) * t700;
t768 = t500 * t703 + t501 * t700;
t552 = -pkin(3) * t813 - t565;
t767 = t552 * t616 - t861;
t727 = qJ(5) * t745 + t552;
t510 = t581 * pkin(4) + t727;
t765 = -qJD(4) * t510 + t861;
t763 = t703 * pkin(4) - t795;
t762 = -t688 + t895;
t528 = -pkin(3) * t573 + pkin(9) * t572 + t810;
t634 = t701 * t806;
t636 = t704 * t806;
t532 = -qJD(3) * t585 + t885 * t634 + t884 * t636;
t761 = t528 * t703 - t700 * t532 - t564 * t821 - t586 * t820;
t760 = -pkin(4) * t850 + t585;
t757 = t500 * t820 + t501 * t821 + t464 - t807;
t756 = qJ(6) * t843 + t678 + t780;
t755 = -0.2e1 * pkin(1) * t816 - pkin(7) * qJDD(2);
t754 = -t606 * pkin(4) - qJ(5) * t605 - t705 * t706;
t750 = t572 * t703 - t804;
t749 = pkin(3) * t842 + t608 * pkin(4) + pkin(9) * t845 + qJ(5) * t607 + t705 * t688;
t533 = t884 * t634 - t636 * t885 + t646 * t799 - t648 * t800;
t746 = t700 * t528 + t703 * t532 + t564 * t820 - t586 * t821;
t744 = t703 * t871 - t795;
t743 = -qJDD(4) - t758;
t506 = -pkin(3) * t812 + t781;
t720 = t515 * qJ(5) + qJD(5) * t745 + t506;
t469 = t516 * pkin(4) + t720;
t741 = t469 * t703 + t500 * t618 - t510 * t852 + t808;
t740 = t506 * t700 - t518 * t618 + t552 * t820 + t809;
t739 = -t469 * t700 - t501 * t618 - t510 * t851 - t809;
t738 = -t506 * t703 - t517 * t618 + t552 * t821 - t808;
t735 = -t686 * t821 + t770;
t483 = -t515 + t856;
t707 = qJD(2) ^ 2;
t734 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t707 + t774;
t708 = qJD(1) ^ 2;
t733 = pkin(1) * t708 - pkin(7) * qJDD(1) + t775;
t732 = g(1) * t607 + g(2) * t605 + g(3) * t848 - t782;
t731 = pkin(4) * t751 - qJ(5) * t804 + t533;
t729 = -qJDD(5) + t732;
t728 = -t644 * t618 - t781 - t902;
t726 = g(1) * t608 + g(2) * t606 + g(3) * t846 + t783;
t472 = qJ(5) * t573 - qJD(5) * t630 - t746;
t725 = qJD(4) * t768 + t464 - t867;
t462 = t581 * qJD(6) + t516 * t871 + t720;
t480 = -t610 * t871 + t819;
t489 = t581 * t871 + t727;
t724 = -t462 * t703 - t480 * t618 + t489 * t904 - t808;
t482 = qJD(6) - t501 - t878;
t723 = -t462 * t700 + t482 * t618 - t489 * t905 - t809;
t721 = -t510 * t745 - t729;
t459 = -qJD(6) * t610 - t766 - t798 - t880;
t461 = -t463 + t890;
t719 = t459 * t700 + t461 * t703 + t480 * t905 - t482 * t904 - t807;
t718 = ((-t515 - t856) * t703 + (-t516 + t854) * t700) * MDP(19) + (-t745 * t787 - t864) * MDP(18) + (-t581 * t618 - t603 * t700 + t859) * MDP(21) + (t610 * t787 - t618 * t745 + t860) * MDP(20) + (t616 * t813 + t710) * MDP(13) + (-t618 * t813 - t538) * MDP(14) + (-t616 ^ 2 + t618 ^ 2) * MDP(12) + t812 * MDP(15) + (-MDP(11) * t616 + t610 * MDP(22)) * t618;
t717 = -t489 * t745 - t729 - t880;
t716 = -t489 * t581 - t726 + t890;
t715 = -t644 * t616 - t722 + t807;
t712 = t744 * t889;
t711 = -g(1) * t657 - g(2) * t653 - g(3) * t780 + t763 * t889;
t695 = t703 * pkin(5);
t694 = t700 * pkin(5);
t687 = -t811 - pkin(3);
t649 = qJ(5) * t846;
t647 = pkin(9) * t703 + t695;
t645 = pkin(9) * t700 + t694;
t627 = t686 * t703 + t695;
t626 = t686 * t700 + t694;
t625 = -t811 - t763;
t613 = -t688 * qJDD(1) + t675;
t594 = -t811 - t744;
t534 = -pkin(4) * t745 + t868;
t531 = qJ(5) * t849 + t760;
t522 = -t742 * t769 + t760;
t521 = -pkin(4) * t630 - t788;
t519 = -t745 * t871 + t868;
t513 = -t789 + t881;
t512 = t604 - t836;
t509 = t790 + t881;
t508 = t604 - t835;
t507 = pkin(5) * t850 - t520;
t499 = t576 + (-pkin(5) * t742 - t564) * t703 - t871 * t630;
t474 = (qJ(5) * t572 + qJD(5) * t742) * t703 + t731;
t473 = pkin(4) * t573 - t761;
t471 = -t769 * t572 - (qJD(6) * t700 + (qJ(6) * qJD(4) - qJD(5)) * t703) * t742 + t731;
t470 = -pkin(5) * t751 - t472;
t467 = -pkin(5) * t750 - qJD(6) * t630 + t573 * t871 - t761;
t1 = [(-(-t581 * t703 + t700 * t745) * t572 - (t864 - t863 + (t853 + t855) * qJD(4)) * t742) * MDP(19) + (-t467 * t745 - t470 * t581 - t499 * t515 - t507 * t516 - (t480 * t703 - t482 * t700) * t572 - (t459 * t703 - t461 * t700 + (-t480 * t700 - t482 * t703) * qJD(4)) * t742 + t778) * MDP(29) + (t472 * t581 - t473 * t745 - t515 * t521 + t516 * t520 - t768 * t572 - (qJD(4) * t892 + t463 * t700 + t465 * t703) * t742 + t778) * MDP(25) + (-t572 * t813 - t742 * t812) * MDP(13) + (-t532 * t813 + t644 * t572 - t586 * t812 - t613 * t742 - t618 * t810 - t688 * t710 - t778) * MDP(17) + (t618 * t572 - t710 * t742) * MDP(11) + (t515 * t849 + t745 * t750) * MDP(18) + (qJDD(1) * t696 + 0.2e1 * t701 * t796) * MDP(4) + (t462 * t522 + t489 * t471 + t459 * t499 + t480 * t467 + t461 * t507 + t482 * t470 - g(1) * (-qJ(6) * t606 + t754) - g(2) * (pkin(5) * t845 + qJ(6) * t608 + t749) + (-g(1) * (t762 - t678) + t875) * t702) * MDP(32) + (t469 * t531 + t510 * t474 + t463 * t520 + t501 * t472 + t465 * t521 + t500 * t473 - g(1) * t754 - g(2) * t749 + (-g(1) * t762 + t875) * t702) * MDP(28) + 0.2e1 * (t701 * t814 - t816 * t823) * MDP(5) + t774 * MDP(2) + t775 * MDP(3) + (t701 * t755 + t704 * t734) * MDP(9) + (-t701 * t734 + t704 * t755) * MDP(10) + (t538 * t742 + t572 * t616 - t618 * t573 - t630 * t710) * MDP(12) + (-t463 * t630 + t469 * t849 - t472 * t610 + t474 * t745 + t501 * t573 + t510 * t750 + t515 * t531 - t520 * t537 + t777) * MDP(27) + (t461 * t630 + t462 * t849 + t470 * t610 + t471 * t745 - t482 * t573 + t489 * t750 + t507 * t537 + t515 * t522 + t777) * MDP(30) + (-t515 * t630 - t537 * t849 + t573 * t745 - t610 * t750) * MDP(20) + (-t506 * t849 - t585 * t515 + t518 * t573 - t533 * t745 - t537 * t833 - t552 * t750 - t610 * t746 + t630 * t783 - t777) * MDP(24) + (t465 * t630 + t469 * t850 + t473 * t610 - t474 * t581 - t500 * t573 - t510 * t751 - t516 * t531 + t521 * t537 - t776) * MDP(26) + (-t506 * t850 + t585 * t516 + t517 * t573 + t533 * t581 + t537 * t788 + t552 * t751 + t610 * t761 - t630 * t782 + t776) * MDP(23) + (-t459 * t630 - t462 * t850 - t467 * t610 + t471 * t581 + t480 * t573 + t489 * t751 - t499 * t537 + t516 * t522 + t776) * MDP(31) + (-t516 * t630 + t537 * t850 + t573 * t581 - t610 * t751) * MDP(21) + (t573 * t813 - t630 * t812) * MDP(14) + (t537 * t630 - t573 * t610) * MDP(22) + (-t533 * t813 - t688 * t538 + t644 * t573 - t585 * t812 + t613 * t630 + t616 * t810 + t693 * t774) * MDP(16) + (qJDD(2) * t701 + t704 * t707) * MDP(6) + (qJDD(2) * t704 - t701 * t707) * MDP(7) + qJDD(1) * MDP(1); (-t515 * t626 - t516 * t627 - t581 * t831 - t745 * t830 + t719) * MDP(29) + (t515 * t594 + t537 * t627 + t610 * t831 + t745 * t834 + t723) * MDP(30) + (-t625 * t516 + t765 * t700 - t829 * t581 + (-t509 + t736) * t610 + t741) * MDP(26) + (t625 * t515 + t765 * t703 + t829 * t745 + (t508 + t735) * t610 + t739) * MDP(27) + (t701 * t733 - t873) * MDP(9) + (g(3) * t701 + t704 * t733) * MDP(10) + (-t508 * t581 + t509 * t745 + (-t745 * t784 + t501 * t616 + (qJD(4) * t581 - t515) * t686) * t700 + (-t581 * t784 + t500 * t616 - t463 + (-t516 - t903) * t686) * t703 + t757) * MDP(25) + (-t687 * t515 + t767 * t703 - t764 * t745 + (-t735 + t835) * t610 + t740) * MDP(24) + (t469 * t625 - t501 * t508 - t500 * t509 + t829 * t510 + (t775 * t701 - t800 * t892 - t873) * pkin(2) + t725 * t686 + t711) * MDP(28) + (t570 * t813 + (t618 * t822 - t800 * t813 - t812 * t884) * pkin(2) + t715) * MDP(17) + (t516 * t594 - t537 * t626 + t581 * t834 - t610 * t830 + t724) * MDP(31) + t718 + (t687 * t516 + t767 * t700 + t764 * t581 + (-t771 + t561 + (-qJD(4) * t686 - t551) * t703) * t610 + t738) * MDP(23) + (t462 * t594 + t459 * t626 + t461 * t627 - g(1) * (-t705 * t883 + t824) - g(2) * (-t702 * t883 + t825) - g(3) * (t756 + t872) + t834 * t489 + t831 * t482 + t830 * t480 + t712) * MDP(32) + MDP(7) * t814 + MDP(6) * t815 + qJDD(2) * MDP(8) + (t569 * t813 + (-t616 * t822 - t799 * t813 + t812 * t885) * pkin(2) + t728) * MDP(16) + (-MDP(4) * t701 * t704 + MDP(5) * t823) * t708; (-t515 * t645 - t516 * t647 - t581 * t827 - t745 * t826 + t719) * MDP(29) + (-t515 * t744 + t537 * t647 + t610 * t827 + t745 * t832 + t723) * MDP(30) + (-t510 * t821 - t513 * t610 + t516 * t763 - t581 * t828 + t741 + t899) * MDP(26) + (-t510 * t820 + t512 * t610 - t515 * t763 + t745 * t828 + t739 - t900) * MDP(27) + (-t516 * t744 - t537 * t645 + t581 * t832 - t610 * t826 + t724) * MDP(31) + (pkin(9) * t725 - t469 * t763 - t500 * t513 - t501 * t512 + t510 * t828 + t711) * MDP(28) + (t565 * t813 + t715) * MDP(17) + t718 + (t566 * t813 + t728) * MDP(16) + (-g(1) * t824 - g(2) * t825 - g(3) * t756 + t459 * t645 + t461 * t647 - t462 * t744 + t480 * t826 + t482 * t827 + t489 * t832 + t712) * MDP(32) + (-t867 - t512 * t581 + t513 * t745 + t768 * t616 + (-t864 - t863 + (-t853 + t855) * qJD(4)) * pkin(9) + t757) * MDP(25) + (-pkin(3) * t516 + t552 * t852 - t566 * t581 - t610 * t789 + t738 - t899) * MDP(23) + (pkin(3) * t515 + t552 * t851 + t566 * t745 + t610 * t836 + t740 + t900) * MDP(24); -MDP(18) * t857 + (t579 - t887) * MDP(19) + t483 * MDP(20) + (-t516 - t854) * MDP(21) + t537 * MDP(22) + (t552 * t745 + t732 + t862) * MDP(23) + (-t517 * t610 + t552 * t581 + t726) * MDP(24) + (pkin(4) * t515 - t869 - (-t501 - t518) * t745 + (t500 - t818) * t581) * MDP(25) + (t534 * t581 + t721 - t862 - 0.2e1 * t882) * MDP(26) + (-t510 * t581 - t534 * t745 + t610 * t818 + t595 - t726 + t886) * MDP(27) + (-t463 * qJ(5) - t465 * pkin(4) - t510 * t534 - t500 * t518 - g(1) * t791 - g(2) * t792 - g(3) * (-pkin(4) * t848 + t649) - t818 * t501) * MDP(28) + (-t869 + t515 * t871 - (t482 + t817) * t745 + (t480 - t819) * t581) * MDP(29) + (-t519 * t745 + t610 * t759 + 0.2e1 * t595 + t716 + t886) * MDP(30) + (-t519 * t581 + (0.2e1 * qJD(6) + t492) * t610 + 0.2e1 * t798 - t717) * MDP(31) + (-t459 * t871 + t461 * qJ(5) - t489 * t519 - g(1) * (-qJ(6) * t607 + t791) - g(2) * (-qJ(6) * t605 + t792) - g(3) * (-t848 * t871 + t649) + t819 * t482 + t817 * t480) * MDP(32); (t501 * t610 + t721 - t882) * MDP(28) + ((-qJD(6) - t482) * t610 - t798 + t717) * MDP(32) + (MDP(25) + MDP(29)) * t483 + (MDP(27) + MDP(30)) * (-t579 - t603) + (-t743 + t857 - t888) * (MDP(26) - MDP(31)); (t618 * t820 + t674 - t854 + (-t701 * t794 - t704 * t793 - t785) * t700) * MDP(29) + (-t743 - t857) * MDP(30) + (-t603 - t887) * MDP(31) + (t480 * t610 + t716 - t898) * MDP(32) + (MDP(29) * t858 - MDP(30) * t573) * qJD(1);];
tau  = t1;
