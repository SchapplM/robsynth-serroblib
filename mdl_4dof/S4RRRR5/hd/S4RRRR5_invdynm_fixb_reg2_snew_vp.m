% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S4RRRR5
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
%
% Output:
% m_new_reg [(3*5)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S4RRRR5_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_invdynm_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_invdynm_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR5_invdynm_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR5_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_invdynm_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:28:26
% EndTime: 2019-12-31 17:28:33
% DurationCPUTime: 6.80s
% Computational Cost: add. (33469->493), mult. (67765->668), div. (0->0), fcn. (45833->8), ass. (0->333)
t844 = sin(qJ(4));
t845 = sin(qJ(3));
t849 = cos(qJ(3));
t846 = sin(qJ(2));
t887 = qJD(1) * t846;
t806 = -t849 * qJD(2) + t845 * t887;
t807 = t845 * qJD(2) + t849 * t887;
t848 = cos(qJ(4));
t768 = t848 * t806 + t844 * t807;
t770 = -t844 * t806 + t848 * t807;
t725 = t770 * t768;
t884 = qJD(1) * qJD(2);
t835 = t846 * t884;
t850 = cos(qJ(2));
t881 = t850 * qJDD(1);
t812 = -t835 + t881;
t803 = -qJDD(3) + t812;
t800 = -qJDD(4) + t803;
t919 = -t725 - t800;
t923 = t844 * t919;
t922 = t848 * t919;
t910 = t807 * t806;
t859 = -t803 - t910;
t921 = t845 * t859;
t920 = t849 * t859;
t836 = t850 * t884;
t883 = t846 * qJDD(1);
t811 = t836 + t883;
t864 = -t845 * qJDD(2) - t849 * t811;
t764 = -t806 * qJD(3) - t864;
t865 = t849 * qJDD(2) - t845 * t811;
t858 = t807 * qJD(3) - t865;
t697 = -t768 * qJD(4) + t848 * t764 - t844 * t858;
t886 = t850 * qJD(1);
t832 = -qJD(3) + t886;
t824 = -qJD(4) + t832;
t756 = t768 * t824;
t918 = t756 + t697;
t793 = t806 * t832;
t739 = t764 - t793;
t875 = t844 * t764 + t848 * t858;
t668 = (qJD(4) + t824) * t770 + t875;
t766 = t768 ^ 2;
t767 = t770 ^ 2;
t917 = t806 ^ 2;
t802 = t807 ^ 2;
t823 = t824 ^ 2;
t830 = t832 ^ 2;
t916 = qJD(2) ^ 2;
t915 = pkin(2) * t846;
t914 = pkin(2) * t850;
t847 = sin(qJ(1));
t851 = cos(qJ(1));
t821 = t847 * g(1) - t851 * g(2);
t852 = qJD(1) ^ 2;
t798 = qJDD(1) * pkin(1) + t852 * pkin(5) + t821;
t867 = -t812 + t835;
t868 = t811 + t836;
t734 = pkin(2) * t867 - pkin(6) * t868 - t798;
t822 = t851 * g(1) + t847 * g(2);
t799 = -t852 * pkin(1) + qJDD(1) * pkin(5) - t822;
t787 = -t846 * g(3) + t850 * t799;
t870 = -pkin(6) * t846 - t914;
t809 = t870 * qJD(1);
t751 = -t916 * pkin(2) + qJDD(2) * pkin(6) + t809 * t886 + t787;
t699 = -t849 * t734 + t845 * t751;
t651 = t859 * pkin(3) - t739 * pkin(7) - t699;
t700 = t845 * t734 + t849 * t751;
t788 = -t832 * pkin(3) - t807 * pkin(7);
t655 = -t917 * pkin(3) - pkin(7) * t858 + t832 * t788 + t700;
t619 = -t848 * t651 + t844 * t655;
t620 = t844 * t651 + t848 * t655;
t589 = -t848 * t619 + t844 * t620;
t913 = pkin(3) * t589;
t671 = -t756 + t697;
t635 = -t668 * t844 - t848 * t671;
t912 = pkin(3) * t635;
t911 = t850 * g(3);
t909 = t824 * t770;
t908 = t824 * t844;
t907 = t824 * t848;
t906 = t832 * t845;
t905 = t832 * t849;
t840 = t846 ^ 2;
t904 = t840 * t852;
t750 = t911 - qJDD(2) * pkin(2) - t916 * pkin(6) + (qJD(1) * t809 + t799) * t846;
t681 = pkin(3) * t858 - t917 * pkin(7) + t807 * t788 + t750;
t903 = t844 * t681;
t708 = -t725 + t800;
t902 = t844 * t708;
t901 = t845 * t589;
t900 = t845 * t750;
t759 = t803 - t910;
t899 = t845 * t759;
t898 = t846 * t798;
t831 = t850 * t852 * t846;
t819 = qJDD(2) + t831;
t897 = t846 * t819;
t820 = qJDD(2) - t831;
t896 = t846 * t820;
t895 = t848 * t681;
t894 = t848 * t708;
t893 = t849 * t589;
t892 = t849 * t750;
t891 = t849 * t759;
t890 = t850 * t798;
t889 = t850 * t820;
t841 = t850 ^ 2;
t888 = t840 + t841;
t882 = t847 * qJDD(1);
t880 = t851 * qJDD(1);
t879 = t846 * t725;
t878 = t846 * t910;
t877 = t850 * t725;
t876 = t850 * t910;
t590 = t844 * t619 + t848 * t620;
t649 = t845 * t699 + t849 * t700;
t786 = t846 * t799 + t911;
t743 = t846 * t786 + t850 * t787;
t874 = -t847 * t821 - t851 * t822;
t873 = t847 * t831;
t872 = t851 * t831;
t871 = -pkin(2) * t750 + pkin(6) * t649;
t816 = -t847 * t852 + t880;
t869 = -pkin(4) * t816 - t847 * g(3);
t648 = -t849 * t699 + t845 * t700;
t742 = t850 * t786 - t846 * t787;
t866 = t851 * t821 - t847 * t822;
t714 = -t823 - t766;
t661 = t844 * t714 + t922;
t863 = pkin(3) * t661 - t619;
t765 = -t830 - t917;
t712 = t849 * t765 - t921;
t794 = t832 * t807;
t736 = t794 - t858;
t862 = pkin(2) * t736 + pkin(6) * t712 - t892;
t773 = -t802 - t830;
t716 = -t845 * t773 + t891;
t740 = (qJD(3) - t832) * t806 + t864;
t861 = pkin(2) * t740 + pkin(6) * t716 + t900;
t744 = -t767 - t823;
t677 = t848 * t744 + t902;
t860 = pkin(3) * t677 - t620;
t737 = (-qJD(3) - t832) * t807 + t865;
t695 = t849 * t737 + t845 * t739;
t758 = t802 + t917;
t857 = pkin(2) * t758 + pkin(6) * t695 + t649;
t637 = -t668 * t848 + t844 * t671;
t701 = -t766 - t767;
t580 = -pkin(3) * t701 + pkin(7) * t637 + t590;
t581 = -pkin(7) * t635 - t589;
t600 = -t845 * t635 + t849 * t637;
t856 = -pkin(2) * t701 + pkin(6) * t600 + t849 * t580 + t845 * t581;
t662 = t848 * t714 - t923;
t667 = (qJD(4) - t824) * t770 + t875;
t609 = -pkin(3) * t667 + pkin(7) * t662 - t895;
t627 = -t845 * t661 + t849 * t662;
t638 = -pkin(7) * t661 + t903;
t855 = -pkin(2) * t667 + pkin(6) * t627 + t849 * t609 + t845 * t638;
t678 = -t844 * t744 + t894;
t613 = -pkin(3) * t918 + pkin(7) * t678 + t903;
t640 = -t845 * t677 + t849 * t678;
t641 = -pkin(7) * t677 + t895;
t854 = -pkin(2) * t918 + pkin(6) * t640 + t849 * t613 + t845 * t641;
t578 = t849 * t590 - t901;
t586 = -pkin(3) * t681 + pkin(7) * t590;
t853 = -pkin(2) * t681 + pkin(6) * t578 - pkin(7) * t901 + t849 * t586;
t838 = t841 * t852;
t829 = -t838 - t916;
t828 = t838 - t916;
t827 = -t904 - t916;
t826 = -t904 + t916;
t818 = -t838 + t904;
t817 = t838 + t904;
t815 = t851 * t852 + t882;
t814 = t888 * qJDD(1);
t813 = -0.2e1 * t835 + t881;
t810 = 0.2e1 * t836 + t883;
t805 = t850 * t819;
t804 = t888 * t884;
t796 = -pkin(4) * t815 + t851 * g(3);
t792 = -t802 + t830;
t791 = -t830 + t917;
t790 = t850 * t811 - t840 * t884;
t789 = -t846 * t812 - t841 * t884;
t785 = -t846 * t827 - t889;
t784 = -t846 * t826 + t805;
t783 = t850 * t829 - t897;
t782 = t850 * t828 - t896;
t781 = t850 * t827 - t896;
t780 = t850 * t826 + t897;
t779 = t846 * t829 + t805;
t778 = t846 * t828 + t889;
t777 = t868 * t846;
t776 = t867 * t850;
t774 = t802 - t917;
t772 = -t846 * t810 + t850 * t813;
t771 = t850 * t810 + t846 * t813;
t755 = -t767 + t823;
t754 = t766 - t823;
t753 = -pkin(5) * t781 - t890;
t752 = -pkin(5) * t779 - t898;
t749 = (t806 * t849 - t807 * t845) * t832;
t748 = (t806 * t845 + t807 * t849) * t832;
t746 = -pkin(1) * t781 + t787;
t745 = -pkin(1) * t779 + t786;
t738 = t764 + t793;
t735 = t794 + t858;
t733 = pkin(1) * t813 + pkin(5) * t783 + t890;
t732 = -pkin(1) * t810 + pkin(5) * t785 - t898;
t729 = t849 * t764 + t807 * t906;
t728 = t845 * t764 - t807 * t905;
t727 = -t806 * t905 + t845 * t858;
t726 = t806 * t906 + t849 * t858;
t724 = t850 * t749 - t846 * t803;
t723 = t846 * t749 + t850 * t803;
t722 = t767 - t766;
t721 = t849 * t791 + t899;
t720 = -t845 * t792 + t920;
t719 = t845 * t791 - t891;
t718 = t849 * t792 + t921;
t717 = pkin(1) * t798 + pkin(5) * t743;
t715 = t849 * t773 + t899;
t713 = pkin(1) * t817 + pkin(5) * t814 + t743;
t711 = t845 * t765 + t920;
t707 = (t768 * t848 - t770 * t844) * t824;
t706 = (t768 * t844 + t770 * t848) * t824;
t705 = t850 * t729 + t878;
t704 = t850 * t727 - t878;
t703 = t846 * t729 - t876;
t702 = t846 * t727 + t876;
t696 = -t770 * qJD(4) - t875;
t694 = t849 * t736 - t845 * t738;
t693 = t845 * t737 - t849 * t739;
t692 = t845 * t736 + t849 * t738;
t691 = -pkin(6) * t715 + t892;
t690 = t850 * t721 - t846 * t735;
t689 = t850 * t720 + t846 * t739;
t688 = t846 * t721 + t850 * t735;
t687 = t846 * t720 - t850 * t739;
t686 = t848 * t754 + t902;
t685 = -t844 * t755 + t922;
t684 = t844 * t754 - t894;
t683 = t848 * t755 + t923;
t682 = -pkin(6) * t711 + t900;
t680 = t850 * t716 - t846 * t740;
t679 = t846 * t716 + t850 * t740;
t676 = t850 * t712 - t846 * t736;
t675 = t846 * t712 + t850 * t736;
t674 = t850 * t694 + t846 * t774;
t673 = t846 * t694 - t850 * t774;
t666 = t848 * t697 + t770 * t908;
t665 = t844 * t697 - t770 * t907;
t664 = -t844 * t696 - t768 * t907;
t663 = t848 * t696 - t768 * t908;
t660 = t850 * t695 - t846 * t758;
t659 = t846 * t695 + t850 * t758;
t658 = -pkin(2) * t715 + t700;
t657 = -t845 * t706 + t849 * t707;
t656 = t849 * t706 + t845 * t707;
t654 = -pkin(2) * t711 + t699;
t653 = t850 * t657 - t846 * t800;
t652 = t846 * t657 + t850 * t800;
t647 = -t845 * t684 + t849 * t686;
t646 = -t845 * t683 + t849 * t685;
t645 = t849 * t684 + t845 * t686;
t644 = t849 * t683 + t845 * t685;
t643 = t850 * t649 + t846 * t750;
t642 = t846 * t649 - t850 * t750;
t639 = t849 * t677 + t845 * t678;
t636 = -t848 * t667 - t844 * t918;
t634 = -t844 * t667 + t848 * t918;
t633 = -pkin(1) * t679 - t861;
t632 = -t845 * t665 + t849 * t666;
t631 = -t845 * t663 + t849 * t664;
t630 = t849 * t665 + t845 * t666;
t629 = t849 * t663 + t845 * t664;
t628 = -pkin(1) * t675 - t862;
t626 = t849 * t661 + t845 * t662;
t625 = -pkin(6) * t693 - t648;
t624 = t850 * t632 + t879;
t623 = t850 * t631 - t879;
t622 = t846 * t632 - t877;
t621 = t846 * t631 + t877;
t617 = t850 * t647 - t846 * t668;
t616 = t850 * t646 + t846 * t671;
t615 = t846 * t647 + t850 * t668;
t614 = t846 * t646 - t850 * t671;
t612 = t850 * t640 + t846 * t918;
t611 = t846 * t640 - t850 * t918;
t610 = -pkin(5) * t679 - t846 * t658 + t850 * t691;
t608 = -pkin(5) * t675 - t846 * t654 + t850 * t682;
t607 = t850 * t627 + t846 * t667;
t606 = t846 * t627 - t850 * t667;
t605 = -pkin(1) * t715 + pkin(5) * t680 + t850 * t658 + t846 * t691;
t604 = -pkin(1) * t659 - t857;
t603 = -pkin(1) * t642 - t871;
t602 = -pkin(1) * t711 + pkin(5) * t676 + t850 * t654 + t846 * t682;
t601 = -pkin(5) * t659 + t850 * t625 + t693 * t915;
t599 = -t845 * t634 + t849 * t636;
t598 = t849 * t635 + t845 * t637;
t597 = t849 * t634 + t845 * t636;
t596 = t850 * t599 + t846 * t722;
t595 = t846 * t599 - t850 * t722;
t594 = -pkin(5) * t642 + (-pkin(6) * t850 + t915) * t648;
t593 = t850 * t600 + t846 * t701;
t592 = t846 * t600 - t850 * t701;
t591 = pkin(5) * t660 + t846 * t625 + (-pkin(1) - t914) * t693;
t588 = -pkin(2) * t639 - t860;
t587 = -pkin(2) * t626 - t863;
t585 = -pkin(2) * t598 - t912;
t584 = pkin(5) * t643 + (-pkin(1) + t870) * t648;
t583 = -pkin(6) * t639 - t845 * t613 + t849 * t641;
t582 = -pkin(6) * t626 - t845 * t609 + t849 * t638;
t579 = -pkin(1) * t611 - t854;
t577 = t845 * t590 + t893;
t576 = -pkin(1) * t606 - t855;
t575 = t850 * t578 + t846 * t681;
t574 = t846 * t578 - t850 * t681;
t573 = -pkin(5) * t611 + t850 * t583 - t846 * t588;
t572 = -pkin(2) * t577 - t913;
t571 = -pkin(5) * t606 + t850 * t582 - t846 * t587;
t570 = -pkin(1) * t639 + pkin(5) * t612 + t846 * t583 + t850 * t588;
t569 = -pkin(1) * t626 + pkin(5) * t607 + t846 * t582 + t850 * t587;
t568 = -pkin(6) * t598 - t845 * t580 + t849 * t581;
t567 = -pkin(6) * t577 - pkin(7) * t893 - t845 * t586;
t566 = -pkin(1) * t592 - t856;
t565 = -pkin(5) * t592 + t850 * t568 - t846 * t585;
t564 = -pkin(1) * t574 - t853;
t563 = -pkin(1) * t598 + pkin(5) * t593 + t846 * t568 + t850 * t585;
t562 = -pkin(5) * t574 + t850 * t567 - t846 * t572;
t561 = -pkin(1) * t577 + pkin(5) * t575 + t846 * t567 + t850 * t572;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t816, 0, -t815, 0, t869, -t796, -t866, -pkin(4) * t866, t851 * t790 - t873, t851 * t772 + t847 * t818, t851 * t784 + t846 * t882, t851 * t789 + t873, t851 * t782 + t847 * t881, t847 * qJDD(2) + t851 * t804, t851 * t752 - t847 * t745 - pkin(4) * (t847 * t783 + t851 * t813), t851 * t753 - t847 * t746 - pkin(4) * (t847 * t785 - t851 * t810), t851 * t742 - pkin(4) * (t847 * t814 + t851 * t817), -pkin(4) * (t847 * t743 + t851 * t798) - (t847 * pkin(1) - t851 * pkin(5)) * t742, t851 * t705 + t847 * t728, t851 * t674 + t847 * t692, t851 * t689 + t847 * t718, t851 * t704 - t847 * t726, t851 * t690 + t847 * t719, t851 * t724 + t847 * t748, t851 * t608 - t847 * t628 - pkin(4) * (t847 * t676 - t851 * t711), t851 * t610 - t847 * t633 - pkin(4) * (t847 * t680 - t851 * t715), t851 * t601 - t847 * t604 - pkin(4) * (t847 * t660 - t851 * t693), t851 * t594 - t847 * t603 - pkin(4) * (t847 * t643 - t851 * t648), t851 * t624 + t847 * t630, t851 * t596 + t847 * t597, t851 * t616 + t847 * t644, t851 * t623 + t847 * t629, t851 * t617 + t847 * t645, t851 * t653 + t847 * t656, t851 * t571 - t847 * t576 - pkin(4) * (t847 * t607 - t851 * t626), t851 * t573 - t847 * t579 - pkin(4) * (t847 * t612 - t851 * t639), t851 * t565 - t847 * t566 - pkin(4) * (t847 * t593 - t851 * t598), t851 * t562 - t847 * t564 - pkin(4) * (t847 * t575 - t851 * t577); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t815, 0, t816, 0, t796, t869, t874, pkin(4) * t874, t847 * t790 + t872, t847 * t772 - t851 * t818, t847 * t784 - t846 * t880, t847 * t789 - t872, t847 * t782 - t850 * t880, -t851 * qJDD(2) + t847 * t804, t847 * t752 + t851 * t745 + pkin(4) * (t851 * t783 - t847 * t813), t847 * t753 + t851 * t746 + pkin(4) * (t851 * t785 + t847 * t810), t847 * t742 + pkin(4) * (t851 * t814 - t847 * t817), pkin(4) * (t851 * t743 - t847 * t798) - (-t851 * pkin(1) - t847 * pkin(5)) * t742, t847 * t705 - t851 * t728, t847 * t674 - t851 * t692, t847 * t689 - t851 * t718, t847 * t704 + t851 * t726, t847 * t690 - t851 * t719, t847 * t724 - t851 * t748, t847 * t608 + t851 * t628 + pkin(4) * (t851 * t676 + t847 * t711), t847 * t610 + t851 * t633 + pkin(4) * (t851 * t680 + t847 * t715), t847 * t601 + t851 * t604 + pkin(4) * (t851 * t660 + t847 * t693), t847 * t594 + t851 * t603 + pkin(4) * (t851 * t643 + t847 * t648), t847 * t624 - t851 * t630, t847 * t596 - t851 * t597, t847 * t616 - t851 * t644, t847 * t623 - t851 * t629, t847 * t617 - t851 * t645, t847 * t653 - t851 * t656, t847 * t571 + t851 * t576 + pkin(4) * (t851 * t607 + t847 * t626), t847 * t573 + t851 * t579 + pkin(4) * (t851 * t612 + t847 * t639), t847 * t565 + t851 * t566 + pkin(4) * (t851 * t593 + t847 * t598), t847 * t562 + t851 * t564 + pkin(4) * (t851 * t575 + t847 * t577); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t821, t822, 0, 0, t777, t771, t780, -t776, t778, 0, t733, t732, t713, t717, t703, t673, t687, t702, t688, t723, t602, t605, t591, t584, t622, t595, t614, t621, t615, t652, t569, t570, t563, t561; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t852, 0, 0, -g(3), -t821, 0, t790, t772, t784, t789, t782, t804, t752, t753, t742, pkin(5) * t742, t705, t674, t689, t704, t690, t724, t608, t610, t601, t594, t624, t596, t616, t623, t617, t653, t571, t573, t565, t562; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t852, 0, qJDD(1), 0, g(3), 0, -t822, 0, t831, -t818, -t883, -t831, -t881, -qJDD(2), t745, t746, 0, pkin(1) * t742, -t728, -t692, -t718, t726, -t719, -t748, t628, t633, t604, t603, -t630, -t597, -t644, -t629, -t645, -t656, t576, t579, t566, t564; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t821, t822, 0, 0, t777, t771, t780, -t776, t778, 0, t733, t732, t713, t717, t703, t673, t687, t702, t688, t723, t602, t605, t591, t584, t622, t595, t614, t621, t615, t652, t569, t570, t563, t561; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t811, t813, t819, -t836, t828, t836, 0, -t798, t786, 0, t729, t694, t720, t727, t721, t749, t682, t691, t625, -pkin(6) * t648, t632, t599, t646, t631, t647, t657, t582, t583, t568, t567; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t835, t810, t826, t812, t820, -t835, t798, 0, t787, 0, -t910, -t774, -t739, t910, t735, t803, t654, t658, -pkin(2) * t693, -pkin(2) * t648, -t725, -t722, -t671, t725, t668, t800, t587, t588, t585, t572; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t831, t818, t883, t831, t881, qJDD(2), -t786, -t787, 0, 0, t728, t692, t718, -t726, t719, t748, t862, t861, t857, t871, t630, t597, t644, t629, t645, t656, t855, t854, t856, t853; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t764, t736, t859, -t793, t791, t793, 0, t750, t699, 0, t666, t636, t685, t664, t686, t707, t638, t641, t581, -pkin(7) * t589; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t794, t738, t792, -t858, -t759, t794, -t750, 0, t700, 0, t665, t634, t683, t663, t684, t706, t609, t613, t580, t586; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t910, t774, t739, -t910, -t735, -t803, -t699, -t700, 0, 0, t725, t722, t671, -t725, -t668, -t800, t863, t860, t912, t913; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t697, -t667, t919, -t756, t754, t756, 0, t681, t619, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t909, t918, t755, t696, -t708, t909, -t681, 0, t620, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t725, t722, t671, -t725, -t668, -t800, -t619, -t620, 0, 0;];
m_new_reg = t1;