% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RRPRR1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR1_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_invdynm_fixb_reg2_snew_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:22:42
% EndTime: 2019-07-18 17:22:49
% DurationCPUTime: 7.43s
% Computational Cost: add. (32476->466), mult. (71438->583), div. (0->0), fcn. (45788->8), ass. (0->346)
t906 = cos(qJ(2));
t902 = sin(qJ(2));
t909 = qJD(1) ^ 2;
t963 = t902 * t909;
t884 = t906 * t963;
t874 = qJDD(2) + t884;
t994 = pkin(1) * t874;
t900 = sin(qJ(5));
t901 = sin(qJ(4));
t905 = cos(qJ(4));
t858 = (t901 * t906 + t902 * t905) * qJD(1);
t904 = cos(qJ(5));
t954 = qJD(2) + qJD(4);
t831 = t900 * t858 - t904 * t954;
t833 = t904 * t858 + t900 * t954;
t797 = t833 * t831;
t956 = qJD(1) * qJD(2);
t890 = t906 * t956;
t892 = t902 * qJDD(1);
t868 = t892 + t890;
t889 = t902 * t956;
t893 = t906 * qJDD(1);
t869 = t893 - t889;
t934 = -t901 * t868 + t905 * t869;
t917 = t858 * qJD(4) - t934;
t798 = qJDD(5) + t917;
t989 = -t797 + t798;
t993 = t900 * t989;
t957 = qJD(1) * t902;
t856 = -t905 * t906 * qJD(1) + t901 * t957;
t826 = t858 * t856;
t896 = qJDD(2) + qJDD(4);
t988 = -t826 + t896;
t992 = t901 * t988;
t991 = t904 * t989;
t990 = t905 * t988;
t952 = t954 ^ 2;
t986 = t858 ^ 2;
t841 = -t952 - t986;
t898 = t906 ^ 2;
t894 = t898 * t909;
t903 = sin(qJ(1));
t907 = cos(qJ(1));
t877 = g(1) * t907 + g(2) * t903;
t852 = -t902 * g(3) - t906 * t877;
t929 = -t869 * qJ(3) - t852;
t955 = qJD(1) * qJD(3);
t951 = 0.2e1 * t955;
t920 = t906 * t951 - t929;
t978 = pkin(3) + qJ(3);
t985 = pkin(1) + pkin(2);
t760 = t869 * pkin(3) - t985 * t894 + (-qJD(2) * t985 + t957 * t978) * qJD(2) + t920;
t942 = -qJ(3) * t890 - t994;
t964 = t902 * t877;
t987 = -0.2e1 * t902 * t955 + t964;
t910 = qJDD(2) * pkin(2) - t978 * t868 + (pkin(2) * t963 + pkin(3) * t956 - g(3)) * t906 - t942 + t987;
t705 = t901 * t760 - t905 * t910;
t706 = t905 * t760 + t901 * t910;
t936 = t905 * t705 - t901 * t706;
t664 = t901 * t705 + t905 * t706;
t921 = t905 * t868 + t901 * t869;
t800 = -t856 * qJD(4) + t921;
t848 = t954 * t856;
t780 = t800 - t848;
t823 = t826 + t896;
t908 = qJD(2) ^ 2;
t881 = t894 - t908;
t876 = t903 * g(1) - t907 * g(2);
t918 = qJD(2) * pkin(1) - qJ(3) * t957;
t919 = -t918 * t957 - qJDD(3) + t876;
t768 = t985 * t869 + t978 * t894 - (qJD(2) * pkin(2) - pkin(3) * t957) * t957 + t919;
t829 = t831 ^ 2;
t830 = t833 ^ 2;
t854 = qJD(5) + t856;
t850 = t854 ^ 2;
t855 = t856 ^ 2;
t980 = t906 * g(3);
t912 = t868 * qJ(3) + t942 + t980;
t932 = -t877 + t951;
t793 = t902 * t932 + t912;
t984 = pkin(1) * t793;
t694 = pkin(4) * t823 + t706;
t911 = -pkin(4) * t780 - t768;
t660 = t900 * t694 - t904 * t911;
t661 = t904 * t694 + t900 * t911;
t632 = -t660 * t904 + t661 * t900;
t983 = pkin(4) * t632;
t633 = t900 * t660 + t904 * t661;
t982 = pkin(4) * t633;
t981 = t903 * g(3);
t979 = t907 * g(3);
t977 = qJ(3) * t793;
t976 = qJ(3) * t902;
t744 = t797 + t798;
t975 = t744 * t900;
t974 = t744 * t904;
t973 = t823 * t901;
t972 = t823 * t905;
t971 = t854 * t900;
t970 = t854 * t904;
t969 = t874 * t902;
t968 = t874 * t906;
t897 = t902 ^ 2;
t967 = t897 * t909;
t695 = pkin(4) * t841 + t705;
t966 = t900 * t695;
t763 = t901 * t768;
t864 = t902 * t876;
t962 = t904 * t695;
t961 = t905 * t768;
t865 = t906 * t876;
t935 = -t900 * t800 + t904 * t896;
t739 = (-qJD(5) + t854) * t833 + t935;
t922 = -t904 * t800 - t900 * t896;
t755 = -qJD(5) * t831 - t922;
t809 = t854 * t831;
t741 = t755 + t809;
t689 = t739 * t900 - t741 * t904;
t625 = -pkin(4) * t689 - t632;
t691 = t739 * t904 + t741 * t900;
t761 = t829 + t830;
t667 = t691 * t905 - t761 * t901;
t960 = pkin(3) * t667 + t901 * t625;
t813 = -t952 - t855;
t767 = t813 * t905 - t992;
t959 = pkin(3) * t767 + t961;
t783 = -t841 * t901 - t972;
t958 = pkin(3) * t783 - t763;
t953 = t905 * t983;
t950 = pkin(1) * t892;
t949 = t901 * t797;
t948 = t905 * t797;
t947 = t903 * t826;
t946 = t907 * t826;
t774 = -t850 - t829;
t703 = t774 * t900 + t991;
t653 = -pkin(4) * t703 + t966;
t704 = t774 * t904 - t993;
t810 = t854 * t833;
t916 = qJD(5) * t833 - t935;
t738 = -t810 - t916;
t674 = t704 * t905 - t738 * t901;
t945 = pkin(3) * t674 + t901 * t653 + t905 * t660;
t792 = -t830 - t850;
t707 = t792 * t904 - t975;
t657 = -pkin(4) * t707 + t962;
t708 = -t792 * t900 - t974;
t742 = (qJD(5) + t854) * t831 + t922;
t678 = t708 * t905 - t742 * t901;
t944 = pkin(3) * t678 + t901 * t657 + t905 * t661;
t778 = qJD(2) * t858 + t934;
t781 = t800 + t848;
t726 = t778 * t905 + t781 * t901;
t943 = pkin(3) * t726 + t664;
t941 = -pkin(4) * t901 - pkin(2);
t940 = t903 * t892;
t939 = t907 * t892;
t938 = t905 * t653 - t901 * t660;
t937 = t905 * t657 - t901 * t661;
t933 = t967 - t881;
t849 = t954 * t858;
t931 = t903 * t884;
t930 = t907 * t884;
t928 = t901 * t848;
t927 = t901 * t849;
t926 = t905 * t848;
t925 = t905 * t849;
t924 = pkin(4) * t704 - t962;
t923 = pkin(4) * t708 + t966;
t914 = pkin(4) * t691 + t633;
t913 = -t869 * pkin(1) - t919;
t879 = t908 - t967;
t875 = qJDD(2) - t884;
t873 = -t894 + t967;
t872 = qJDD(1) * t907 - t903 * t909;
t871 = qJDD(1) * t903 + t907 * t909;
t870 = t893 - 0.2e1 * t889;
t867 = t892 + 0.2e1 * t890;
t863 = (t897 + t898) * t956;
t851 = -t964 + t980;
t847 = t952 - t986;
t846 = t855 - t952;
t845 = qJDD(2) * t903 + t863 * t907;
t844 = t868 * t906 - t897 * t956;
t843 = -qJDD(2) * t907 + t863 * t903;
t842 = -t869 * t902 - t898 * t956;
t840 = -t879 * t902 + t968;
t839 = -t875 * t902 + t881 * t906;
t838 = t879 * t906 + t969;
t837 = t875 * t906 + t881 * t902;
t836 = (t868 + t890) * t902;
t835 = (t869 - t889) * t906;
t834 = -pkin(1) * t867 - qJ(3) * t875;
t828 = -t867 * t902 + t870 * t906;
t827 = t867 * t906 + t870 * t902;
t825 = -t855 + t986;
t821 = t844 * t907 - t931;
t820 = t842 * t907 + t931;
t819 = t844 * t903 + t930;
t818 = t842 * t903 - t930;
t817 = t840 * t907 + t940;
t816 = t839 * t907 + t893 * t903;
t815 = t840 * t903 - t939;
t814 = t839 * t903 - t893 * t907;
t812 = t851 * t906 - t852 * t902;
t811 = t851 * t902 + t852 * t906;
t808 = qJ(3) * t894 - t913;
t807 = -t830 + t850;
t806 = t829 - t850;
t805 = -t926 + t927;
t804 = -t928 - t925;
t803 = -t855 - t986;
t802 = t828 * t907 + t873 * t903;
t801 = t828 * t903 - t873 * t907;
t796 = qJ(3) * t933 + t913;
t795 = -pkin(1) * t894 - qJD(2) * t918 + t920;
t794 = t830 - t829;
t791 = (qJ(3) * qJDD(1) + t932) * t902 + t912;
t790 = (-qJD(2) * t976 - 0.2e1 * qJD(3) * t906) * qJD(1) + (-t933 + t908) * pkin(1) + t929;
t789 = t846 * t905 - t973;
t788 = -t847 * t901 + t990;
t787 = t846 * t901 + t972;
t786 = t847 * t905 + t992;
t785 = -t912 + t987 + t994;
t784 = -t908 * qJ(3) + (t869 + t870) * pkin(1) + t919;
t782 = t841 * t905 - t973;
t779 = (-0.2e1 * qJD(4) - qJD(2)) * t856 + t921;
t777 = -t849 + t917;
t776 = t849 + t917;
t773 = (t889 + t893) * qJ(3) + ((t897 - t898) * t909 + t881) * pkin(1) + t920;
t772 = t905 * t800 - t927;
t771 = t901 * t800 + t925;
t770 = t901 * t917 + t926;
t769 = -t905 * t917 + t928;
t766 = t813 * t901 + t990;
t759 = -qJ(3) * t968 - t784 * t902;
t758 = -qJ(3) * t969 + t784 * t906;
t757 = t796 * t906 - t834 * t902;
t756 = t796 * t902 + t834 * t906;
t750 = (-t831 * t904 + t833 * t900) * t854;
t749 = (-t831 * t900 - t833 * t904) * t854;
t748 = -t804 * t902 + t805 * t906;
t747 = t804 * t906 + t805 * t902;
t746 = pkin(1) * t808 + qJ(3) * t795;
t740 = t755 - t809;
t737 = -t810 + t916;
t736 = -t787 * t902 + t789 * t906;
t735 = -t786 * t902 + t788 * t906;
t734 = t787 * t906 + t789 * t902;
t733 = t786 * t906 + t788 * t902;
t732 = t755 * t904 - t833 * t971;
t731 = t755 * t900 + t833 * t970;
t730 = t831 * t970 + t900 * t916;
t729 = -t831 * t971 + t904 * t916;
t728 = -t773 * t902 + t791 * t906;
t727 = t773 * t906 + t791 * t902;
t725 = -t776 * t905 - t780 * t901;
t724 = t778 * t901 - t781 * t905;
t723 = -t776 * t901 + t780 * t905;
t721 = t750 * t905 + t798 * t901;
t720 = t750 * t901 - t798 * t905;
t718 = t806 * t904 - t975;
t717 = -t807 * t900 + t991;
t716 = t806 * t900 + t974;
t715 = t807 * t904 + t993;
t714 = -t771 * t902 + t772 * t906;
t713 = -t769 * t902 + t770 * t906;
t712 = t771 * t906 + t772 * t902;
t711 = t769 * t906 + t770 * t902;
t710 = -t746 * t902 + t906 * t977;
t709 = t746 * t906 + t793 * t976;
t699 = t732 * t905 + t949;
t698 = t730 * t905 - t949;
t697 = t732 * t901 - t948;
t696 = t730 * t901 + t948;
t693 = -t782 * t978 - t961;
t692 = -t766 * t978 - t763;
t690 = t738 * t904 - t740 * t900;
t688 = t738 * t900 + t740 * t904;
t687 = t985 * t724;
t686 = -t723 * t902 + t725 * t906;
t685 = t723 * t906 + t725 * t902;
t684 = t718 * t905 - t737 * t901;
t683 = t717 * t905 + t741 * t901;
t682 = t718 * t901 + t737 * t905;
t681 = t717 * t901 - t741 * t905;
t680 = -t720 * t902 + t721 * t906;
t679 = t720 * t906 + t721 * t902;
t677 = t708 * t901 + t742 * t905;
t675 = t782 * t985 - t706;
t673 = t704 * t901 + t738 * t905;
t671 = t690 * t905 + t794 * t901;
t670 = t690 * t901 - t794 * t905;
t669 = t766 * t985 - t705;
t668 = qJ(3) * t783 - t779 * t985 + t958;
t666 = t691 * t901 + t761 * t905;
t662 = pkin(3) * t664;
t656 = qJ(3) * t767 - t776 * t985 + t959;
t650 = -t697 * t902 + t699 * t906;
t649 = -t696 * t902 + t698 * t906;
t648 = t697 * t906 + t699 * t902;
t647 = t696 * t906 + t698 * t902;
t646 = -t668 * t902 + t693 * t906;
t645 = t668 * t906 + t693 * t902;
t644 = -t682 * t902 + t684 * t906;
t643 = -t681 * t902 + t683 * t906;
t642 = t682 * t906 + t684 * t902;
t641 = t681 * t906 + t683 * t902;
t640 = -t656 * t902 + t692 * t906;
t639 = t656 * t906 + t692 * t902;
t638 = -t724 * t978 + t936;
t637 = -t670 * t902 + t671 * t906;
t636 = t670 * t906 + t671 * t902;
t635 = t985 * t936;
t634 = t978 * t936;
t631 = qJ(3) * t726 - t803 * t985 + t943;
t630 = qJ(3) * t664 + t768 * t985 + t662;
t629 = t633 * t905 + t695 * t901;
t628 = t633 * t901 - t695 * t905;
t627 = pkin(3) * t629;
t626 = t677 * t985 + t923;
t624 = t905 * t625;
t622 = t673 * t985 + t924;
t621 = -t677 * t978 + t937;
t620 = -t673 * t978 + t938;
t619 = -t631 * t902 + t638 * t906;
t618 = t631 * t906 + t638 * t902;
t617 = -t630 * t902 + t634 * t906;
t616 = t630 * t906 + t634 * t902;
t615 = qJ(3) * t678 - t707 * t985 + t944;
t614 = t666 * t985 + t914;
t613 = -t666 * t978 + t624;
t612 = qJ(3) * t674 - t703 * t985 + t945;
t611 = qJ(3) * t667 - t689 * t985 + t960;
t610 = t628 * t985 + t982;
t609 = -t628 * t978 - t953;
t608 = -t615 * t902 + t621 * t906;
t607 = t615 * t906 + t621 * t902;
t606 = -t612 * t902 + t620 * t906;
t605 = t612 * t906 + t620 * t902;
t604 = -t611 * t902 + t613 * t906;
t603 = t611 * t906 + t613 * t902;
t602 = qJ(3) * t629 + t627 + (-pkin(1) + t941) * t632;
t601 = -t602 * t902 + t609 * t906;
t600 = t602 * t906 + t609 * t902;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t872, 0, -t871, 0, -t981, -t979, -t876 * t907 + t877 * t903, 0, t821, t802, t817, t820, t816, t845, -t851 * t903 - t864 * t907, -t852 * t903 - t865 * t907, t907 * t812, 0, t821, t802, t817, t820, t816, t845, t759 * t907 + t785 * t903, t757 * t907 + t790 * t903, -pkin(1) * t940 + t728 * t907, t710 * t907 - t903 * t984, t714 * t907 + t947, t686 * t907 + t825 * t903, t735 * t907 + t781 * t903, t713 * t907 - t947, t736 * t907 - t777 * t903, t748 * t907 + t896 * t903, t640 * t907 + t669 * t903, t646 * t907 + t675 * t903, t619 * t907 + t687 * t903, t617 * t907 - t635 * t903, t650 * t907 + t731 * t903, t637 * t907 + t688 * t903, t643 * t907 + t715 * t903, t649 * t907 - t729 * t903, t644 * t907 + t716 * t903, t680 * t907 + t749 * t903, t606 * t907 + t622 * t903, t608 * t907 + t626 * t903, t604 * t907 + t614 * t903, t601 * t907 + t610 * t903; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t871, 0, t872, 0, t979, -t981, -t876 * t903 - t877 * t907, 0, t819, t801, t815, t818, t814, t843, t851 * t907 - t864 * t903, t852 * t907 - t865 * t903, t903 * t812, 0, t819, t801, t815, t818, t814, t843, t759 * t903 - t785 * t907, t757 * t903 - t790 * t907, pkin(1) * t939 + t728 * t903, t710 * t903 + t907 * t984, t714 * t903 - t946, t686 * t903 - t825 * t907, t735 * t903 - t781 * t907, t713 * t903 + t946, t736 * t903 + t777 * t907, t748 * t903 - t896 * t907, t640 * t903 - t669 * t907, t646 * t903 - t675 * t907, t619 * t903 - t687 * t907, t617 * t903 + t635 * t907, t650 * t903 - t731 * t907, t637 * t903 - t688 * t907, t643 * t903 - t715 * t907, t649 * t903 + t729 * t907, t644 * t903 - t716 * t907, t680 * t903 - t749 * t907, t606 * t903 - t622 * t907, t608 * t903 - t626 * t907, t604 * t903 - t614 * t907, t601 * t903 - t610 * t907; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t876, t877, 0, 0, t836, t827, t838, t835, t837, 0, t865, -t864, t811, 0, t836, t827, t838, t835, t837, 0, t758, t756, t727, t709, t712, t685, t733, t711, t734, t747, t639, t645, t618, t616, t648, t636, t641, t647, t642, t679, t605, t607, t603, t600; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t909, 0, 0, -g(3), -t876, 0, t844, t828, t840, t842, t839, t863, -t864, -t865, t812, 0, t844, t828, t840, t842, t839, t863, t759, t757, t728, t710, t714, t686, t735, t713, t736, t748, t640, t646, t619, t617, t650, t637, t643, t649, t644, t680, t606, t608, t604, t601; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t909, 0, qJDD(1), 0, g(3), 0, -t877, 0, t884, -t873, -t892, -t884, -t893, -qJDD(2), t851, t852, 0, 0, t884, -t873, -t892, -t884, -t893, -qJDD(2), -t785, -t790, t950, t984, -t826, -t825, -t781, t826, t777, -t896, -t669, -t675, -t687, t635, -t731, -t688, -t715, t729, -t716, -t749, -t622, -t626, -t614, -t610; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t876, t877, 0, 0, t836, t827, t838, t835, t837, 0, t865, -t864, t811, 0, t836, t827, t838, t835, t837, 0, t758, t756, t727, t709, t712, t685, t733, t711, t734, t747, t639, t645, t618, t616, t648, t636, t641, t647, t642, t679, t605, t607, t603, t600; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t868, t870, t874, -t890, t881, t890, 0, -t876, t851, 0, t868, t870, t874, -t890, t881, t890, -qJ(3) * t874, t796, t791, t977, t772, t725, t788, t770, t789, t805, t692, t693, t638, t634, t699, t671, t683, t698, t684, t721, t620, t621, t613, t609; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t889, t867, t879, t869, t875, -t889, t876, 0, t852, 0, t889, t867, t879, t869, t875, -t889, t784, t834, t773, t746, t771, t723, t786, t769, t787, t804, t656, t668, t631, t630, t697, t670, t681, t696, t682, t720, t612, t615, t611, t602; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t884, t873, t892, t884, t893, qJDD(2), -t851, -t852, 0, 0, -t884, t873, t892, t884, t893, qJDD(2), t785, t790, -t950, -t984, t826, t825, t781, -t826, -t777, t896, t669, t675, t687, -t635, t731, t688, t715, -t729, t716, t749, t622, t626, t614, t610; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t868, t870, t874, -t890, t881, t890, 0, -t808, t793, 0, t772, t725, t788, t770, t789, t805, -pkin(3) * t766 - t763, -pkin(3) * t782 - t961, -pkin(3) * t724 + t936, pkin(3) * t936, t699, t671, t683, t698, t684, t721, -pkin(3) * t673 + t938, -pkin(3) * t677 + t937, -pkin(3) * t666 + t624, -pkin(3) * t628 - t953; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t889, t867, t879, t869, t875, -t889, t808, 0, t795, 0, t771, t723, t786, t769, t787, t804, -pkin(2) * t776 + t959, -pkin(2) * t779 + t958, -pkin(2) * t803 + t943, pkin(2) * t768 + t662, t697, t670, t681, t696, t682, t720, -pkin(2) * t703 + t945, -pkin(2) * t707 + t944, -pkin(2) * t689 + t960, t632 * t941 + t627; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t884, t873, t892, t884, t893, qJDD(2), -t793, -t795, 0, 0, t826, t825, t781, -t826, -t777, t896, pkin(2) * t766 - t705, pkin(2) * t782 - t706, pkin(2) * t724, -pkin(2) * t936, t731, t688, t715, -t729, t716, t749, pkin(2) * t673 + t924, pkin(2) * t677 + t923, pkin(2) * t666 + t914, pkin(2) * t628 + t982; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t800, -t776, t988, t848, t846, -t848, 0, -t768, t705, 0, t732, t690, t717, t730, t718, t750, t653, t657, t625, -t983; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t849, t780, t847, -t917, t823, -t849, t768, 0, t706, 0, -t797, -t794, -t741, t797, t737, -t798, t660, t661, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t826, t825, t781, -t826, -t777, t896, -t705, -t706, 0, 0, t731, t688, t715, -t729, t716, t749, t924, t923, t914, t982; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t755, t738, t989, t809, t806, -t809, 0, t695, t660, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t810, t740, t807, -t916, t744, -t810, -t695, 0, t661, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t797, t794, t741, -t797, -t737, t798, -t660, -t661, 0, 0;];
m_new_reg  = t1;