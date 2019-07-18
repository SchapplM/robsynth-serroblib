% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RPRRPR2
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
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 22:08
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RPRRPR2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR2_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_invdynB_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:06:55
% EndTime: 2019-05-05 22:07:32
% DurationCPUTime: 35.43s
% Computational Cost: add. (262018->877), mult. (531738->1353), div. (0->0), fcn. (371592->12), ass. (0->599)
t956 = sin(qJ(1));
t960 = cos(qJ(1));
t921 = g(1) * t960 + g(2) * t956;
t961 = qJD(1) ^ 2;
t905 = -pkin(1) * t961 - t921;
t948 = sin(pkin(10));
t950 = cos(pkin(10));
t920 = g(1) * t956 - t960 * g(2);
t967 = qJDD(1) * pkin(1) + t920;
t854 = t948 * t905 - t950 * t967;
t855 = t950 * t905 + t948 * t967;
t785 = t854 * t950 - t855 * t948;
t1018 = t785 * t956;
t976 = t854 * t948 + t950 * t855;
t707 = t960 * t976 + t1018;
t1017 = t785 * t960;
t1057 = -t956 * t976 + t1017;
t912 = qJDD(1) * t950 - t948 * t961;
t994 = g(3) - qJDD(2);
t1040 = -qJ(2) * t912 - t948 * t994;
t911 = qJDD(1) * t948 + t950 * t961;
t1041 = t960 * t911 + t912 * t956;
t882 = -qJ(2) * t911 + t950 * t994;
t1056 = -pkin(6) * t1041 + t1040 * t956 + t882 * t960;
t857 = -t911 * t956 + t960 * t912;
t1055 = -pkin(6) * t857 + t1040 * t960 - t882 * t956;
t947 = sin(pkin(11));
t954 = sin(qJ(4));
t958 = cos(qJ(4));
t955 = sin(qJ(3));
t992 = qJD(1) * t955;
t902 = -t958 * qJD(3) + t954 * t992;
t903 = qJD(3) * t954 + t958 * t992;
t949 = cos(pkin(11));
t851 = t949 * t902 + t903 * t947;
t853 = -t902 * t947 + t903 * t949;
t1010 = t851 * t853;
t988 = qJD(1) * qJD(3);
t936 = t955 * t988;
t959 = cos(qJ(3));
t986 = qJDD(1) * t959;
t909 = -t936 + t986;
t899 = -qJDD(4) + t909;
t966 = -t899 - t1010;
t1051 = t947 * t966;
t1050 = t949 * t966;
t953 = sin(qJ(6));
t957 = cos(qJ(6));
t780 = t957 * t851 + t853 * t953;
t782 = -t851 * t953 + t853 * t957;
t1019 = t780 * t782;
t892 = -qJDD(6) + t899;
t962 = -t892 - t1019;
t1049 = t953 * t962;
t1007 = t902 * t903;
t965 = -t899 - t1007;
t1048 = t954 * t965;
t1046 = t957 * t962;
t1045 = t958 * t965;
t991 = qJD(1) * t959;
t931 = -qJD(4) + t991;
t923 = -qJD(6) + t931;
t757 = t780 * t923;
t979 = t959 * t988;
t987 = qJDD(1) * t955;
t908 = t979 + t987;
t847 = -t902 * qJD(4) + t954 * qJDD(3) + t958 * t908;
t968 = t958 * qJDD(3) - t954 * t908;
t963 = qJD(4) * t903 - t968;
t771 = t949 * t847 - t947 * t963;
t977 = t847 * t947 + t949 * t963;
t964 = qJD(6) * t780 - t771 * t957 + t953 * t977;
t1043 = t757 - t964;
t826 = t851 * t931;
t736 = t826 - t771;
t1042 = t826 + t771;
t884 = t902 * t931;
t807 = t884 - t847;
t805 = t884 + t847;
t978 = t953 * t771 + t957 * t977;
t632 = (qJD(6) + t923) * t782 + t978;
t808 = (qJD(4) + t931) * t903 - t968;
t775 = t780 ^ 2;
t776 = t782 ^ 2;
t1039 = t851 ^ 2;
t850 = t853 ^ 2;
t1038 = t902 ^ 2;
t898 = t903 ^ 2;
t922 = t923 ^ 2;
t1037 = t931 ^ 2;
t1036 = qJD(3) ^ 2;
t1035 = pkin(3) * t955;
t1034 = pkin(3) * t959;
t832 = -qJDD(1) * pkin(2) - t961 * pkin(7) + t854;
t969 = -t909 + t936;
t970 = t908 + t979;
t769 = pkin(3) * t969 - pkin(8) * t970 + t832;
t833 = -pkin(2) * t961 + qJDD(1) * pkin(7) + t855;
t814 = t959 * t833 - t955 * t994;
t972 = -pkin(8) * t955 - t1034;
t906 = t972 * qJD(1);
t788 = -pkin(3) * t1036 + qJDD(3) * pkin(8) + t906 * t991 + t814;
t703 = -t958 * t769 + t954 * t788;
t656 = t965 * pkin(4) + qJ(5) * t807 - t703;
t704 = t954 * t769 + t958 * t788;
t873 = -pkin(4) * t931 - qJ(5) * t903;
t668 = -pkin(4) * t1038 - qJ(5) * t963 + t931 * t873 + t704;
t579 = 0.2e1 * qJD(5) * t853 - t949 * t656 + t947 * t668;
t540 = t966 * pkin(5) + pkin(9) * t736 - t579;
t580 = -0.2e1 * qJD(5) * t851 + t947 * t656 + t949 * t668;
t821 = -pkin(5) * t931 - pkin(9) * t853;
t548 = -pkin(5) * t1039 - pkin(9) * t977 + t821 * t931 + t580;
t488 = -t957 * t540 + t548 * t953;
t489 = t953 * t540 + t957 * t548;
t443 = -t488 * t957 + t489 * t953;
t1031 = t443 * t947;
t1030 = t443 * t949;
t518 = -t579 * t949 + t580 * t947;
t1029 = t518 * t954;
t1028 = t518 * t958;
t935 = t959 * t994;
t787 = t935 - qJDD(3) * pkin(3) - t1036 * pkin(8) + (qJD(1) * t906 + t833) * t955;
t700 = t963 * pkin(4) - t1038 * qJ(5) + t903 * t873 + qJDD(5) + t787;
t612 = pkin(5) * t977 - pkin(9) * t1039 + t853 * t821 + t700;
t1027 = t612 * t953;
t1026 = t612 * t957;
t697 = t892 - t1019;
t1025 = t697 * t953;
t1024 = t697 * t957;
t1023 = t700 * t947;
t1022 = t700 * t949;
t764 = t899 - t1010;
t1021 = t764 * t947;
t1020 = t764 * t949;
t1016 = t787 * t954;
t1015 = t787 * t958;
t828 = t899 - t1007;
t1014 = t828 * t954;
t1013 = t828 * t958;
t1012 = t832 * t955;
t1011 = t832 * t959;
t1009 = t853 * t931;
t1008 = t899 * t955;
t930 = t959 * t961 * t955;
t918 = -t930 + qJDD(3);
t1004 = t918 * t955;
t1003 = t918 * t959;
t919 = qJDD(3) + t930;
t1002 = t919 * t955;
t1001 = t923 * t953;
t1000 = t923 * t957;
t999 = t931 * t947;
t998 = t931 * t949;
t997 = t931 * t954;
t996 = t931 * t958;
t943 = t955 ^ 2;
t995 = t943 * t961;
t944 = t959 ^ 2;
t993 = t943 + t944;
t985 = t955 * t1019;
t984 = t959 * t1019;
t983 = t955 * t1010;
t982 = t959 * t1010;
t981 = t955 * t1007;
t980 = t959 * t1007;
t444 = t488 * t953 + t957 * t489;
t519 = t579 * t947 + t949 * t580;
t812 = t833 * t955 + t935;
t744 = t812 * t955 + t959 * t814;
t866 = -t920 * t956 - t960 * t921;
t974 = t948 * t930;
t973 = t950 * t930;
t915 = qJDD(1) * t960 - t956 * t961;
t971 = -pkin(6) * t915 - g(3) * t956;
t633 = -t703 * t958 + t704 * t954;
t634 = t703 * t954 + t704 * t958;
t743 = t812 * t959 - t814 * t955;
t865 = t920 * t960 - t921 * t956;
t737 = t977 + t1009;
t941 = t944 * t961;
t928 = -t941 - t1036;
t927 = t941 - t1036;
t926 = -t995 - t1036;
t925 = -t995 + t1036;
t917 = t941 - t995;
t916 = t941 + t995;
t914 = qJDD(1) * t956 + t960 * t961;
t913 = t993 * qJDD(1);
t910 = -0.2e1 * t936 + t986;
t907 = 0.2e1 * t979 + t987;
t901 = t959 * t919;
t900 = t993 * t988;
t886 = t959 * t899;
t885 = -pkin(6) * t914 + g(3) * t960;
t879 = -t898 + t1037;
t878 = -t1037 + t1038;
t877 = t908 * t959 - t943 * t988;
t876 = -t909 * t955 - t944 * t988;
t875 = qJDD(3) * t948 + t900 * t950;
t874 = -qJDD(3) * t950 + t900 * t948;
t872 = -t926 * t955 - t1003;
t871 = -t925 * t955 + t901;
t870 = t928 * t959 - t1002;
t869 = t927 * t959 - t1004;
t868 = t926 * t959 - t1004;
t867 = t928 * t955 + t901;
t864 = -t898 + t1038;
t863 = -t898 - t1037;
t862 = t913 * t950 - t916 * t948;
t861 = t913 * t948 + t916 * t950;
t856 = -t907 * t955 + t910 * t959;
t849 = -t1037 - t1038;
t843 = t877 * t950 - t974;
t842 = t876 * t950 + t974;
t841 = t877 * t948 + t973;
t840 = t876 * t948 - t973;
t839 = t871 * t950 + t948 * t987;
t838 = t869 * t950 + t948 * t986;
t837 = t871 * t948 - t950 * t987;
t836 = t869 * t948 - t950 * t986;
t827 = t898 + t1038;
t823 = -t850 + t1037;
t822 = -t1037 + t1039;
t820 = t872 * t950 + t907 * t948;
t819 = t870 * t950 - t910 * t948;
t818 = t872 * t948 - t907 * t950;
t817 = t870 * t948 + t910 * t950;
t816 = (t902 * t958 - t903 * t954) * t931;
t815 = (-t902 * t954 - t903 * t958) * t931;
t813 = t856 * t950 - t917 * t948;
t811 = t856 * t948 + t917 * t950;
t810 = -t850 - t1037;
t803 = (-qJD(4) + t931) * t903 + t968;
t802 = t847 * t958 + t903 * t997;
t801 = -t847 * t954 + t903 * t996;
t800 = -t902 * t996 + t954 * t963;
t799 = t902 * t997 + t958 * t963;
t798 = -t861 * t956 + t862 * t960;
t797 = t861 * t960 + t862 * t956;
t796 = t816 * t959 - t1008;
t795 = -pkin(7) * t868 + t1011;
t794 = t878 * t958 + t1014;
t793 = -t879 * t954 + t1045;
t792 = -pkin(7) * t867 + t1012;
t791 = -t878 * t954 + t1013;
t790 = -t879 * t958 - t1048;
t789 = -t850 + t1039;
t778 = -pkin(2) * t868 + t814;
t777 = -pkin(2) * t867 + t812;
t774 = -t863 * t954 + t1013;
t773 = t863 * t958 + t1014;
t772 = -t1037 - t1039;
t768 = pkin(1) * t994 + qJ(2) * t976;
t763 = t849 * t958 - t1048;
t762 = t849 * t954 + t1045;
t756 = -t776 + t922;
t755 = t775 - t922;
t754 = (t851 * t949 - t853 * t947) * t931;
t753 = (t851 * t947 + t853 * t949) * t931;
t752 = -t776 - t922;
t751 = t802 * t959 + t981;
t750 = t800 * t959 - t981;
t749 = -t818 * t956 + t820 * t960;
t748 = -t817 * t956 + t819 * t960;
t747 = t818 * t960 + t820 * t956;
t746 = t817 * t960 + t819 * t956;
t745 = -t850 - t1039;
t741 = -t807 * t954 - t808 * t958;
t740 = t803 * t958 - t805 * t954;
t739 = t807 * t958 - t808 * t954;
t738 = -t803 * t954 - t805 * t958;
t732 = t977 - t1009;
t731 = t796 * t950 - t815 * t948;
t730 = t796 * t948 + t815 * t950;
t729 = t771 * t949 + t853 * t999;
t728 = t771 * t947 - t853 * t998;
t727 = -t851 * t998 + t947 * t977;
t726 = -t851 * t999 - t949 * t977;
t725 = t794 * t959 - t808 * t955;
t724 = t793 * t959 - t807 * t955;
t723 = t822 * t949 + t1021;
t722 = -t823 * t947 + t1050;
t721 = t822 * t947 - t1020;
t720 = t823 * t949 + t1051;
t719 = t774 * t959 + t805 * t955;
t718 = -t810 * t947 + t1020;
t717 = t774 * t955 - t805 * t959;
t716 = t810 * t949 + t1021;
t715 = -qJ(2) * t861 + t743 * t950;
t714 = qJ(2) * t862 + t743 * t948;
t713 = t763 * t959 - t803 * t955;
t712 = t763 * t955 + t803 * t959;
t711 = -t776 + t775;
t710 = t744 * t950 + t832 * t948;
t709 = t744 * t948 - t832 * t950;
t708 = -pkin(8) * t773 + t1015;
t705 = t740 * t959 - t864 * t955;
t702 = -pkin(8) * t762 + t1016;
t701 = -t922 - t775;
t696 = t772 * t949 - t1051;
t695 = t772 * t947 + t1050;
t694 = t751 * t950 - t801 * t948;
t693 = t750 * t950 - t799 * t948;
t692 = t751 * t948 + t801 * t950;
t691 = t750 * t948 + t799 * t950;
t690 = t741 * t959 - t827 * t955;
t689 = t741 * t955 + t827 * t959;
t688 = (t780 * t957 - t782 * t953) * t923;
t687 = (t780 * t953 + t782 * t957) * t923;
t686 = -t753 * t954 + t754 * t958;
t685 = -t753 * t958 - t754 * t954;
t684 = -qJ(2) * t818 - t778 * t948 + t795 * t950;
t683 = -qJ(2) * t817 - t777 * t948 + t792 * t950;
t682 = t686 * t959 - t1008;
t681 = t725 * t950 - t791 * t948;
t680 = t724 * t950 - t790 * t948;
t679 = t725 * t948 + t791 * t950;
t678 = t724 * t948 + t790 * t950;
t677 = t719 * t950 + t773 * t948;
t676 = t719 * t948 - t773 * t950;
t675 = -pkin(1) * t868 + qJ(2) * t820 + t778 * t950 + t795 * t948;
t674 = -pkin(1) * t867 + qJ(2) * t819 + t777 * t950 + t792 * t948;
t673 = -t775 - t776;
t672 = t713 * t950 + t762 * t948;
t671 = t713 * t948 - t762 * t950;
t670 = -pkin(3) * t773 + t704;
t669 = -pkin(3) * t762 + t703;
t666 = -qJD(6) * t782 - t978;
t664 = -t736 * t947 - t737 * t949;
t663 = -t1042 * t947 - t732 * t949;
t662 = t736 * t949 - t737 * t947;
t661 = t1042 * t949 - t732 * t947;
t660 = t755 * t957 + t1025;
t659 = -t756 * t953 + t1046;
t658 = t755 * t953 - t1024;
t657 = t756 * t957 + t1049;
t655 = -t728 * t954 + t729 * t958;
t654 = -t726 * t954 + t727 * t958;
t653 = -t728 * t958 - t729 * t954;
t652 = -t726 * t958 - t727 * t954;
t649 = -t752 * t953 + t1024;
t648 = t752 * t957 + t1025;
t647 = -t721 * t954 + t723 * t958;
t646 = -t720 * t954 + t722 * t958;
t645 = -t721 * t958 - t723 * t954;
t644 = -t720 * t958 - t722 * t954;
t643 = t705 * t950 - t738 * t948;
t642 = t705 * t948 + t738 * t950;
t641 = -t716 * t954 + t718 * t958;
t640 = t716 * t958 + t718 * t954;
t639 = t690 * t950 + t739 * t948;
t638 = t690 * t948 - t739 * t950;
t637 = -qJ(5) * t716 + t1022;
t636 = -t709 * t956 + t710 * t960;
t635 = t709 * t960 + t710 * t956;
t631 = t757 + t964;
t627 = (qJD(6) - t923) * t782 + t978;
t626 = t701 * t957 - t1049;
t625 = t701 * t953 + t1046;
t624 = t1001 * t782 - t957 * t964;
t623 = -t1000 * t782 - t953 * t964;
t622 = -t1000 * t780 - t666 * t953;
t621 = -t1001 * t780 + t666 * t957;
t620 = -qJ(5) * t695 + t1023;
t619 = t655 * t959 + t983;
t618 = t654 * t959 - t983;
t617 = -t695 * t954 + t696 * t958;
t616 = t695 * t958 + t696 * t954;
t615 = -pkin(2) * t717 + pkin(3) * t805 - pkin(8) * t774 - t1016;
t614 = -t687 * t947 + t688 * t949;
t613 = t687 * t949 + t688 * t947;
t611 = -pkin(2) * t712 - pkin(3) * t803 - pkin(8) * t763 + t1015;
t610 = -qJ(2) * t709 - (pkin(2) * t948 - pkin(7) * t950) * t743;
t609 = t647 * t959 - t737 * t955;
t608 = t646 * t959 - t736 * t955;
t607 = t682 * t950 - t685 * t948;
t606 = t682 * t948 + t685 * t950;
t605 = t634 * t959 + t787 * t955;
t604 = t634 * t955 - t787 * t959;
t603 = t1042 * t955 + t641 * t959;
t602 = -t1042 * t959 + t641 * t955;
t601 = -pkin(4) * t1042 + qJ(5) * t718 + t1023;
t600 = -pkin(8) * t739 - t633;
t599 = -t676 * t956 + t677 * t960;
t598 = t676 * t960 + t677 * t956;
t597 = -pkin(4) * t732 + qJ(5) * t696 - t1022;
t596 = t617 * t959 + t732 * t955;
t595 = t617 * t955 - t732 * t959;
t594 = -t671 * t956 + t672 * t960;
t593 = t671 * t960 + t672 * t956;
t592 = qJ(2) * t710 - (-pkin(2) * t950 - pkin(7) * t948 - pkin(1)) * t743;
t591 = -t662 * t954 + t664 * t958;
t590 = -t661 * t954 + t663 * t958;
t589 = t662 * t958 + t664 * t954;
t588 = -t661 * t958 - t663 * t954;
t587 = -t658 * t947 + t660 * t949;
t586 = -t657 * t947 + t659 * t949;
t585 = t658 * t949 + t660 * t947;
t584 = t657 * t949 + t659 * t947;
t583 = -pkin(7) * t717 - t670 * t955 + t708 * t959;
t582 = -t648 * t947 + t649 * t949;
t581 = t648 * t949 + t649 * t947;
t577 = -pkin(7) * t712 - t669 * t955 + t702 * t959;
t576 = t590 * t959 - t789 * t955;
t575 = -t638 * t956 + t639 * t960;
t574 = t638 * t960 + t639 * t956;
t573 = t619 * t950 - t653 * t948;
t572 = t618 * t950 - t652 * t948;
t571 = t619 * t948 + t653 * t950;
t570 = t618 * t948 + t652 * t950;
t569 = t591 * t959 + t745 * t955;
t568 = t591 * t955 - t745 * t959;
t567 = -pkin(9) * t648 + t1026;
t566 = -pkin(2) * t689 - pkin(3) * t827 - pkin(8) * t741 - t634;
t565 = -t631 * t953 - t632 * t957;
t564 = -t1043 * t953 - t627 * t957;
t563 = t631 * t957 - t632 * t953;
t562 = t1043 * t957 - t627 * t953;
t561 = -t625 * t947 + t626 * t949;
t560 = t625 * t949 + t626 * t947;
t559 = -t623 * t947 + t624 * t949;
t558 = -t621 * t947 + t622 * t949;
t557 = t623 * t949 + t624 * t947;
t556 = t621 * t949 + t622 * t947;
t555 = t609 * t950 - t645 * t948;
t554 = t608 * t950 - t644 * t948;
t553 = t609 * t948 + t645 * t950;
t552 = t608 * t948 + t644 * t950;
t551 = -pkin(9) * t625 + t1027;
t550 = t603 * t950 + t640 * t948;
t549 = t603 * t948 - t640 * t950;
t547 = -t613 * t954 + t614 * t958;
t546 = -t613 * t958 - t614 * t954;
t544 = t605 * t950 + t633 * t948;
t543 = t605 * t948 - t633 * t950;
t542 = t547 * t959 - t892 * t955;
t541 = -pkin(7) * t689 + t1035 * t739 + t600 * t959;
t537 = -pkin(3) * t589 - pkin(4) * t662;
t536 = -pkin(2) * t604 + pkin(3) * t787 - pkin(8) * t634;
t535 = t596 * t950 + t616 * t948;
t534 = t596 * t948 - t616 * t950;
t533 = -pkin(5) * t1043 + pkin(9) * t649 + t1027;
t532 = -pkin(3) * t640 - pkin(4) * t716 + t580;
t531 = -pkin(5) * t627 + pkin(9) * t626 - t1026;
t530 = -pkin(8) * t640 - t601 * t954 + t637 * t958;
t529 = -pkin(3) * t616 - pkin(4) * t695 + t579;
t528 = -t585 * t954 + t587 * t958;
t527 = -t584 * t954 + t586 * t958;
t526 = -t585 * t958 - t587 * t954;
t525 = -t584 * t958 - t586 * t954;
t524 = -pkin(7) * t604 + (-pkin(8) * t959 + t1035) * t633;
t523 = -t581 * t954 + t582 * t958;
t522 = t581 * t958 + t582 * t954;
t521 = t576 * t950 - t588 * t948;
t520 = t576 * t948 + t588 * t950;
t517 = -qJ(2) * t676 + t583 * t950 - t615 * t948;
t516 = t569 * t950 + t589 * t948;
t515 = t569 * t948 - t589 * t950;
t514 = -pkin(8) * t616 - t597 * t954 + t620 * t958;
t513 = -qJ(2) * t671 + t577 * t950 - t611 * t948;
t512 = -pkin(1) * t717 + qJ(2) * t677 + t583 * t948 + t615 * t950;
t511 = -pkin(4) * t700 + qJ(5) * t519;
t510 = -pkin(1) * t712 + qJ(2) * t672 + t577 * t948 + t611 * t950;
t509 = -t563 * t947 + t565 * t949;
t508 = -t562 * t947 + t564 * t949;
t507 = t563 * t949 + t565 * t947;
t506 = t562 * t949 + t564 * t947;
t505 = -qJ(5) * t662 - t518;
t504 = -t560 * t954 + t561 * t958;
t503 = t560 * t958 + t561 * t954;
t502 = -t557 * t954 + t559 * t958;
t501 = -t556 * t954 + t558 * t958;
t500 = -t557 * t958 - t559 * t954;
t499 = -t556 * t958 - t558 * t954;
t498 = t528 * t959 - t632 * t955;
t497 = t527 * t959 - t631 * t955;
t496 = -t549 * t956 + t550 * t960;
t495 = t549 * t960 + t550 * t956;
t494 = t1043 * t955 + t523 * t959;
t493 = -t1043 * t959 + t523 * t955;
t492 = -pkin(4) * t745 + qJ(5) * t664 + t519;
t491 = t542 * t950 - t546 * t948;
t490 = t542 * t948 + t546 * t950;
t486 = -t543 * t956 + t544 * t960;
t485 = t543 * t960 + t544 * t956;
t484 = t502 * t959 + t985;
t483 = t501 * t959 - t985;
t482 = -pkin(2) * t602 + pkin(3) * t1042 - pkin(8) * t641 - t601 * t958 - t637 * t954;
t481 = -t534 * t956 + t535 * t960;
t480 = t534 * t960 + t535 * t956;
t479 = t504 * t959 + t627 * t955;
t478 = t504 * t955 - t627 * t959;
t477 = -qJ(2) * t638 + t541 * t950 - t566 * t948;
t476 = -pkin(2) * t595 + pkin(3) * t732 - pkin(8) * t617 - t597 * t958 - t620 * t954;
t475 = -pkin(1) * t689 + qJ(2) * t639 + t541 * t948 + t566 * t950;
t474 = -qJ(5) * t581 - t533 * t947 + t567 * t949;
t473 = t519 * t958 - t1029;
t472 = t519 * t954 + t1028;
t471 = -t515 * t956 + t516 * t960;
t470 = t515 * t960 + t516 * t956;
t469 = t473 * t959 + t700 * t955;
t468 = t473 * t955 - t700 * t959;
t467 = -qJ(5) * t560 - t531 * t947 + t551 * t949;
t466 = -pkin(7) * t602 + t530 * t959 - t532 * t955;
t465 = t498 * t950 - t526 * t948;
t464 = t497 * t950 - t525 * t948;
t463 = t498 * t948 + t526 * t950;
t462 = t497 * t948 + t525 * t950;
t461 = -pkin(4) * t1043 + qJ(5) * t582 + t533 * t949 + t567 * t947;
t460 = t494 * t950 + t522 * t948;
t459 = t494 * t948 - t522 * t950;
t458 = -qJ(2) * t543 + t524 * t950 - t536 * t948;
t457 = -t507 * t954 + t509 * t958;
t456 = -t506 * t954 + t508 * t958;
t455 = t507 * t958 + t509 * t954;
t454 = -t506 * t958 - t508 * t954;
t453 = -pkin(7) * t595 + t514 * t959 - t529 * t955;
t452 = -pkin(4) * t627 + qJ(5) * t561 + t531 * t949 + t551 * t947;
t451 = t456 * t959 - t711 * t955;
t450 = t484 * t950 - t500 * t948;
t449 = t483 * t950 - t499 * t948;
t448 = t484 * t948 + t500 * t950;
t447 = t483 * t948 + t499 * t950;
t446 = t457 * t959 + t673 * t955;
t445 = t457 * t955 - t673 * t959;
t442 = -pkin(1) * t604 + qJ(2) * t544 + t524 * t948 + t536 * t950;
t441 = t479 * t950 + t503 * t948;
t440 = t479 * t948 - t503 * t950;
t439 = -pkin(3) * t472 - pkin(4) * t518;
t438 = -pkin(5) * t612 + pkin(9) * t444;
t437 = -pkin(8) * t589 - t492 * t954 + t505 * t958;
t436 = -pkin(3) * t522 - pkin(4) * t581 - pkin(5) * t648 + t489;
t435 = -pkin(9) * t563 - t443;
t434 = -pkin(3) * t503 - pkin(4) * t560 - pkin(5) * t625 + t488;
t433 = -pkin(5) * t673 + pkin(9) * t565 + t444;
t432 = -pkin(2) * t568 + pkin(3) * t745 - pkin(8) * t591 - t492 * t958 - t505 * t954;
t431 = -pkin(8) * t472 - qJ(5) * t1028 - t511 * t954;
t430 = t469 * t950 + t472 * t948;
t429 = t469 * t948 - t472 * t950;
t428 = -pkin(3) * t455 - pkin(4) * t507 - pkin(5) * t563;
t427 = -pkin(7) * t568 + t437 * t959 - t537 * t955;
t426 = -qJ(2) * t549 + t466 * t950 - t482 * t948;
t425 = -t459 * t956 + t460 * t960;
t424 = t459 * t960 + t460 * t956;
t423 = -pkin(1) * t602 + qJ(2) * t550 + t466 * t948 + t482 * t950;
t422 = -qJ(2) * t534 + t453 * t950 - t476 * t948;
t421 = t451 * t950 - t454 * t948;
t420 = t451 * t948 + t454 * t950;
t419 = t446 * t950 + t455 * t948;
t418 = t446 * t948 - t455 * t950;
t417 = -pkin(1) * t595 + qJ(2) * t535 + t453 * t948 + t476 * t950;
t416 = -pkin(8) * t522 - t461 * t954 + t474 * t958;
t415 = t444 * t949 - t1031;
t414 = t444 * t947 + t1030;
t413 = -t440 * t956 + t441 * t960;
t412 = t440 * t960 + t441 * t956;
t411 = -pkin(8) * t503 - t452 * t954 + t467 * t958;
t410 = -pkin(2) * t468 + pkin(3) * t700 - pkin(8) * t473 + qJ(5) * t1029 - t511 * t958;
t409 = -pkin(2) * t493 + pkin(3) * t1043 - pkin(8) * t523 - t461 * t958 - t474 * t954;
t408 = -pkin(2) * t478 + pkin(3) * t627 - pkin(8) * t504 - t452 * t958 - t467 * t954;
t407 = -qJ(5) * t507 - t433 * t947 + t435 * t949;
t406 = -pkin(4) * t673 + qJ(5) * t509 + t433 * t949 + t435 * t947;
t405 = -t429 * t956 + t430 * t960;
t404 = t429 * t960 + t430 * t956;
t403 = -qJ(2) * t515 + t427 * t950 - t432 * t948;
t402 = -pkin(7) * t468 + t431 * t959 - t439 * t955;
t401 = -pkin(7) * t493 + t416 * t959 - t436 * t955;
t400 = -pkin(1) * t568 + qJ(2) * t516 + t427 * t948 + t432 * t950;
t399 = -t418 * t956 + t419 * t960;
t398 = t418 * t960 + t419 * t956;
t397 = -t414 * t954 + t415 * t958;
t396 = t414 * t958 + t415 * t954;
t395 = -pkin(7) * t478 + t411 * t959 - t434 * t955;
t394 = -pkin(9) * t1030 - qJ(5) * t414 - t438 * t947;
t393 = t397 * t959 + t612 * t955;
t392 = t397 * t955 - t612 * t959;
t391 = -pkin(4) * t612 - pkin(9) * t1031 + qJ(5) * t415 + t438 * t949;
t390 = -qJ(2) * t459 + t401 * t950 - t409 * t948;
t389 = -pkin(8) * t455 - t406 * t954 + t407 * t958;
t388 = -pkin(3) * t396 - pkin(4) * t414 - pkin(5) * t443;
t387 = -pkin(1) * t493 + qJ(2) * t460 + t401 * t948 + t409 * t950;
t386 = -qJ(2) * t429 + t402 * t950 - t410 * t948;
t385 = -qJ(2) * t440 + t395 * t950 - t408 * t948;
t384 = -pkin(2) * t445 + pkin(3) * t673 - pkin(8) * t457 - t406 * t958 - t407 * t954;
t383 = t393 * t950 + t396 * t948;
t382 = t393 * t948 - t396 * t950;
t381 = -pkin(1) * t468 + qJ(2) * t430 + t402 * t948 + t410 * t950;
t380 = -pkin(1) * t478 + qJ(2) * t441 + t395 * t948 + t408 * t950;
t379 = -pkin(7) * t445 + t389 * t959 - t428 * t955;
t378 = -pkin(8) * t396 - t391 * t954 + t394 * t958;
t377 = -t382 * t956 + t383 * t960;
t376 = t382 * t960 + t383 * t956;
t375 = -pkin(2) * t392 + pkin(3) * t612 - pkin(8) * t397 - t391 * t958 - t394 * t954;
t374 = -qJ(2) * t418 + t379 * t950 - t384 * t948;
t373 = -pkin(1) * t445 + qJ(2) * t419 + t379 * t948 + t384 * t950;
t372 = -pkin(7) * t392 + t378 * t959 - t388 * t955;
t371 = -qJ(2) * t382 + t372 * t950 - t375 * t948;
t370 = -pkin(1) * t392 + qJ(2) * t383 + t372 * t948 + t375 * t950;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t914, -t915, 0, t866, 0, 0, 0, 0, 0, 0, -t1041, -t857, 0, t707, 0, 0, 0, 0, 0, 0, t748, t749, t798, t636, 0, 0, 0, 0, 0, 0, t594, t599, t575, t486, 0, 0, 0, 0, 0, 0, t481, t496, t471, t405, 0, 0, 0, 0, 0, 0, t413, t425, t399, t377; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t915, -t914, 0, t865, 0, 0, 0, 0, 0, 0, t857, -t1041, 0, -t1057, 0, 0, 0, 0, 0, 0, t746, t747, t797, t635, 0, 0, 0, 0, 0, 0, t593, t598, t574, t485, 0, 0, 0, 0, 0, 0, t480, t495, t470, t404, 0, 0, 0, 0, 0, 0, t412, t424, t398, t376; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t994, 0, 0, 0, 0, 0, 0, t867, t868, 0, -t743, 0, 0, 0, 0, 0, 0, t712, t717, t689, t604, 0, 0, 0, 0, 0, 0, t595, t602, t568, t468, 0, 0, 0, 0, 0, 0, t478, t493, t445, t392; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t915, 0, -t914, 0, t971, -t885, -t865, -pkin(6) * t865, 0, 0, t857, 0, -t1041, 0, t1055, -t1056, t1057, pkin(6) * t1057 + qJ(2) * t1017 - t768 * t956, -t841 * t956 + t843 * t960, -t811 * t956 + t813 * t960, -t837 * t956 + t839 * t960, -t840 * t956 + t842 * t960, -t836 * t956 + t838 * t960, -t874 * t956 + t875 * t960, -pkin(6) * t746 - t674 * t956 + t683 * t960, -pkin(6) * t747 - t675 * t956 + t684 * t960, -pkin(6) * t797 - t714 * t956 + t715 * t960, -pkin(6) * t635 - t592 * t956 + t610 * t960, -t692 * t956 + t694 * t960, -t642 * t956 + t643 * t960, -t678 * t956 + t680 * t960, -t691 * t956 + t693 * t960, -t679 * t956 + t681 * t960, -t730 * t956 + t731 * t960, -pkin(6) * t593 - t510 * t956 + t513 * t960, -pkin(6) * t598 - t512 * t956 + t517 * t960, -pkin(6) * t574 - t475 * t956 + t477 * t960, -pkin(6) * t485 - t442 * t956 + t458 * t960, -t571 * t956 + t573 * t960, -t520 * t956 + t521 * t960, -t552 * t956 + t554 * t960, -t570 * t956 + t572 * t960, -t553 * t956 + t555 * t960, -t606 * t956 + t607 * t960, -pkin(6) * t480 - t417 * t956 + t422 * t960, -pkin(6) * t495 - t423 * t956 + t426 * t960, -pkin(6) * t470 - t400 * t956 + t403 * t960, -pkin(6) * t404 - t381 * t956 + t386 * t960, -t448 * t956 + t450 * t960, -t420 * t956 + t421 * t960, -t462 * t956 + t464 * t960, -t447 * t956 + t449 * t960, -t463 * t956 + t465 * t960, -t490 * t956 + t491 * t960, -pkin(6) * t412 - t380 * t956 + t385 * t960, -pkin(6) * t424 - t387 * t956 + t390 * t960, -pkin(6) * t398 - t373 * t956 + t374 * t960, -pkin(6) * t376 - t370 * t956 + t371 * t960; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t914, 0, t915, 0, t885, t971, t866, pkin(6) * t866, 0, 0, t1041, 0, t857, 0, t1056, t1055, t707, pkin(6) * t707 + qJ(2) * t1018 + t768 * t960, t841 * t960 + t843 * t956, t811 * t960 + t813 * t956, t837 * t960 + t839 * t956, t840 * t960 + t842 * t956, t836 * t960 + t838 * t956, t874 * t960 + t875 * t956, pkin(6) * t748 + t674 * t960 + t683 * t956, pkin(6) * t749 + t675 * t960 + t684 * t956, pkin(6) * t798 + t714 * t960 + t715 * t956, pkin(6) * t636 + t592 * t960 + t610 * t956, t692 * t960 + t694 * t956, t642 * t960 + t643 * t956, t678 * t960 + t680 * t956, t691 * t960 + t693 * t956, t679 * t960 + t681 * t956, t730 * t960 + t731 * t956, pkin(6) * t594 + t510 * t960 + t513 * t956, pkin(6) * t599 + t512 * t960 + t517 * t956, pkin(6) * t575 + t475 * t960 + t477 * t956, pkin(6) * t486 + t442 * t960 + t458 * t956, t571 * t960 + t573 * t956, t520 * t960 + t521 * t956, t552 * t960 + t554 * t956, t570 * t960 + t572 * t956, t553 * t960 + t555 * t956, t606 * t960 + t607 * t956, pkin(6) * t481 + t417 * t960 + t422 * t956, pkin(6) * t496 + t423 * t960 + t426 * t956, pkin(6) * t471 + t400 * t960 + t403 * t956, pkin(6) * t405 + t381 * t960 + t386 * t956, t448 * t960 + t450 * t956, t420 * t960 + t421 * t956, t462 * t960 + t464 * t956, t447 * t960 + t449 * t956, t463 * t960 + t465 * t956, t490 * t960 + t491 * t956, pkin(6) * t413 + t380 * t960 + t385 * t956, pkin(6) * t425 + t387 * t960 + t390 * t956, pkin(6) * t399 + t373 * t960 + t374 * t956, pkin(6) * t377 + t370 * t960 + t371 * t956; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t920, t921, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * t912 - t854, -pkin(1) * t911 - t855, 0, -pkin(1) * t785, t970 * t955, t907 * t959 + t910 * t955, t925 * t959 + t1002, -t969 * t959, t927 * t955 + t1003, 0, pkin(1) * t817 + pkin(2) * t910 + pkin(7) * t870 - t1011, pkin(1) * t818 - pkin(2) * t907 + pkin(7) * t872 + t1012, pkin(1) * t861 + pkin(2) * t916 + pkin(7) * t913 + t744, pkin(1) * t709 - pkin(2) * t832 + pkin(7) * t744, t802 * t955 - t980, t740 * t955 + t864 * t959, t793 * t955 + t807 * t959, t800 * t955 + t980, t794 * t955 + t808 * t959, t816 * t955 + t886, pkin(1) * t671 - pkin(2) * t762 + pkin(7) * t713 + t669 * t959 + t702 * t955, pkin(1) * t676 - pkin(2) * t773 + pkin(7) * t719 + t670 * t959 + t708 * t955, pkin(1) * t638 + pkin(7) * t690 + t955 * t600 + (-pkin(2) - t1034) * t739, pkin(1) * t543 + pkin(7) * t605 + (-pkin(2) + t972) * t633, t655 * t955 - t982, t590 * t955 + t789 * t959, t646 * t955 + t736 * t959, t654 * t955 + t982, t647 * t955 + t737 * t959, t686 * t955 + t886, pkin(1) * t534 - pkin(2) * t616 + pkin(7) * t596 + t514 * t955 + t529 * t959, pkin(1) * t549 - pkin(2) * t640 + pkin(7) * t603 + t530 * t955 + t532 * t959, pkin(1) * t515 - pkin(2) * t589 + pkin(7) * t569 + t437 * t955 + t537 * t959, pkin(1) * t429 - pkin(2) * t472 + pkin(7) * t469 + t431 * t955 + t439 * t959, t502 * t955 - t984, t456 * t955 + t711 * t959, t527 * t955 + t631 * t959, t501 * t955 + t984, t528 * t955 + t632 * t959, t547 * t955 + t892 * t959, pkin(1) * t440 - pkin(2) * t503 + pkin(7) * t479 + t411 * t955 + t434 * t959, pkin(1) * t459 - pkin(2) * t522 + pkin(7) * t494 + t416 * t955 + t436 * t959, pkin(1) * t418 - pkin(2) * t455 + pkin(7) * t446 + t389 * t955 + t428 * t959, pkin(1) * t382 - pkin(2) * t396 + pkin(7) * t393 + t378 * t955 + t388 * t959;];
tauB_reg  = t1;