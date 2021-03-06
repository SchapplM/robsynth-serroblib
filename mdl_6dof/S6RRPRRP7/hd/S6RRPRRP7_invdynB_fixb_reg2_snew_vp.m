% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RRPRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 18:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RRPRRP7_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP7_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP7_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP7_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP7_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP7_invdynB_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:16:40
% EndTime: 2019-05-06 18:17:13
% DurationCPUTime: 25.22s
% Computational Cost: add. (56303->684), mult. (116091->911), div. (0->0), fcn. (75403->8), ass. (0->496)
t873 = sin(qJ(4));
t874 = sin(qJ(2));
t877 = cos(qJ(4));
t878 = cos(qJ(2));
t821 = (t873 * t874 + t877 * t878) * qJD(1);
t817 = qJD(5) + t821;
t1010 = t817 ^ 2;
t950 = qJD(1) * t874;
t823 = -t873 * t878 * qJD(1) + t877 * t950;
t872 = sin(qJ(5));
t876 = cos(qJ(5));
t939 = qJD(2) - qJD(4);
t784 = t872 * t823 + t876 * t939;
t1011 = t784 ^ 2;
t748 = t1011 - t1010;
t786 = t876 * t823 - t872 * t939;
t731 = t786 * t784;
t945 = qJD(1) * qJD(2);
t928 = t878 * t945;
t943 = qJDD(1) * t874;
t837 = t928 + t943;
t929 = t874 * t945;
t941 = qJDD(1) * t878;
t838 = -t929 + t941;
t923 = t873 * t837 + t877 * t838;
t739 = -t823 * qJD(4) - t923;
t734 = qJDD(5) - t739;
t1023 = t731 + t734;
t990 = t1023 * t872;
t635 = t748 * t876 - t990;
t740 = -t821 * qJD(4) + t877 * t837 - t873 * t838;
t938 = qJDD(2) - qJDD(4);
t924 = t872 * t740 + t876 * t938;
t657 = (qJD(5) - t817) * t786 + t924;
t572 = t635 * t873 + t657 * t877;
t577 = t635 * t877 - t657 * t873;
t528 = t572 * t874 + t577 * t878;
t989 = t1023 * t876;
t631 = t748 * t872 + t989;
t875 = sin(qJ(1));
t879 = cos(qJ(1));
t1122 = t528 * t875 + t631 * t879;
t1121 = t528 * t879 - t631 * t875;
t946 = qJD(5) + t817;
t656 = t786 * t946 + t924;
t891 = -t876 * t740 + t872 * t938;
t889 = -t784 * qJD(5) - t891;
t981 = t784 * t817;
t1020 = -t981 + t889;
t993 = t1020 * t872;
t585 = t656 * t876 + t993;
t782 = t786 ^ 2;
t725 = t782 - t1011;
t555 = t585 * t873 + t725 * t877;
t556 = t585 * t877 - t725 * t873;
t509 = t555 * t874 + t556 * t878;
t583 = t1020 * t876 - t656 * t872;
t1120 = t509 * t875 - t583 * t879;
t1119 = t509 * t879 + t583 * t875;
t1021 = -t981 - t889;
t1061 = t1021 * t876 - t657 * t872;
t1019 = t782 + t1011;
t1060 = -t1021 * t872 - t657 * t876;
t1080 = -t1019 * t873 + t1060 * t877;
t1081 = t1019 * t877 + t1060 * t873;
t1091 = t1080 * t878 + t1081 * t874;
t1111 = -t1061 * t875 + t1091 * t879;
t1118 = pkin(6) * t1111;
t1112 = t1061 * t879 + t1091 * t875;
t1117 = pkin(6) * t1112;
t1116 = t572 * t878 - t577 * t874;
t1090 = t1080 * t874 - t1081 * t878;
t1115 = pkin(7) * t1090;
t1008 = pkin(2) + pkin(3);
t1114 = -pkin(1) * t1090 + pkin(4) * t1019 + pkin(9) * t1060 - qJ(3) * t1080 + t1008 * t1081;
t1113 = pkin(1) * t1061 + pkin(7) * t1091;
t1110 = t555 * t878 - t556 * t874;
t722 = -t782 - t1010;
t618 = t722 * t876 - t990;
t1109 = pkin(1) * t618;
t1108 = pkin(4) * t618;
t880 = qJD(2) ^ 2;
t870 = t874 ^ 2;
t881 = qJD(1) ^ 2;
t963 = t870 * t881;
t851 = t880 + t963;
t856 = t878 * t881 * t874;
t847 = qJDD(2) - t856;
t957 = t878 * t847;
t800 = -t851 * t874 + t957;
t836 = 0.2e1 * t928 + t943;
t756 = t800 * t875 + t836 * t879;
t1107 = pkin(6) * t756;
t759 = t800 * t879 - t836 * t875;
t1106 = pkin(6) * t759;
t614 = t722 * t872 + t989;
t1105 = pkin(9) * t614;
t1104 = pkin(9) * t618;
t1103 = pkin(8) * t1080;
t1102 = pkin(8) * t1081;
t1101 = qJ(3) * t618;
t1100 = t614 * t873;
t1099 = t614 * t877;
t1098 = t618 * t875;
t1097 = t618 * t879;
t1094 = t1008 * t618;
t749 = -t782 + t1010;
t1024 = -t731 + t734;
t988 = t1024 * t872;
t1049 = t876 * t749 + t988;
t987 = t1024 * t876;
t1047 = -t749 * t872 + t987;
t1054 = -t1021 * t873 + t1047 * t877;
t1055 = -t1021 * t877 - t1047 * t873;
t1079 = t1054 * t878 - t1055 * t874;
t1093 = t1049 * t879 + t1079 * t875;
t1092 = -t1049 * t875 + t1079 * t879;
t1087 = pkin(9) * t1061;
t1078 = t1054 * t874 + t1055 * t878;
t1077 = pkin(7) * t800;
t1018 = -t1010 - t1011;
t1036 = t1018 * t872 + t987;
t1076 = pkin(1) * t1036;
t1075 = pkin(4) * t1036;
t1035 = t1018 * t876 - t988;
t1074 = pkin(9) * t1035;
t1073 = pkin(9) * t1036;
t1070 = qJ(3) * t1036;
t1069 = qJ(6) * t1020;
t1068 = t1035 * t873;
t1067 = t1035 * t877;
t1066 = t1036 * t875;
t1065 = t1036 * t879;
t1062 = t1008 * t1036;
t686 = -qJD(5) * t786 - t924;
t978 = t817 * t876;
t934 = t784 * t978;
t898 = -t686 * t872 + t934;
t936 = t873 * t731;
t1015 = t877 * t898 - t936;
t935 = t877 * t731;
t1017 = -t873 * t898 - t935;
t1032 = t1015 * t878 - t1017 * t874;
t979 = t817 * t872;
t913 = t876 * t686 + t784 * t979;
t1059 = t1032 * t875 + t879 * t913;
t742 = t786 * t979;
t912 = t742 - t934;
t1014 = t734 * t873 + t877 * t912;
t1016 = t734 * t877 - t873 * t912;
t1034 = t1014 * t878 - t1016 * t874;
t892 = (-t784 * t872 - t786 * t876) * t817;
t1058 = t1034 * t875 + t879 * t892;
t1057 = t1032 * t879 - t875 * t913;
t1056 = t1034 * t879 - t875 * t892;
t1009 = 2 * qJD(3);
t1053 = 2 * qJD(6);
t966 = t847 * t874;
t794 = t851 * t878 + t966;
t1052 = pkin(1) * t794;
t1051 = pkin(7) * t794;
t921 = t821 * t939;
t1048 = t740 + t921;
t839 = -0.2e1 * t929 + t941;
t970 = t839 * t878;
t974 = t836 * t874;
t775 = -t970 + t974;
t871 = t878 ^ 2;
t845 = (t870 - t871) * t881;
t1042 = t775 * t875 + t845 * t879;
t1041 = t775 * t879 - t845 * t875;
t962 = t871 * t881;
t853 = -t880 + t962;
t798 = -t853 * t878 + t966;
t940 = qJDD(1) * t879;
t1040 = t798 * t875 + t878 * t940;
t1039 = t798 * t879 - t875 * t941;
t653 = t876 * t889 - t742;
t899 = -t653 * t873 + t935;
t914 = t877 * t653 + t936;
t1012 = -t874 * t899 + t878 * t914;
t652 = t786 * t978 + t872 * t889;
t1038 = t1012 * t875 + t879 * t652;
t1037 = t1012 * t879 - t652 * t875;
t1033 = t1014 * t874 + t1016 * t878;
t1031 = t1015 * t874 + t1017 * t878;
t937 = t939 ^ 2;
t977 = t821 * t823;
t901 = -t938 - t977;
t1030 = t873 * t901;
t1027 = t877 * t901;
t1006 = pkin(2) * t878;
t907 = -qJ(3) * t874 - t1006;
t894 = t881 * t907;
t850 = g(1) * t879 + g(2) * t875;
t825 = -pkin(1) * t881 + qJDD(1) * pkin(7) - t850;
t925 = t874 * g(3) - t878 * t825;
t895 = qJDD(2) * qJ(3) + qJD(2) * t1009 + t878 * t894 - t925;
t1001 = t838 * pkin(2);
t848 = -qJD(2) * pkin(3) - pkin(8) * t950;
t849 = t875 * g(1) - t879 * g(2);
t824 = qJDD(1) * pkin(1) + t881 * pkin(7) + t849;
t890 = -pkin(2) * t929 + t824;
t949 = qJD(2) * t878;
t677 = t1001 + t837 * qJ(3) + t838 * pkin(3) - pkin(8) * t962 + (qJ(3) * t949 + (t1009 + t848) * t874) * qJD(1) + t890;
t920 = t939 * t823;
t596 = -t1048 * pkin(9) + (-t739 - t920) * pkin(4) + t677;
t733 = -pkin(2) * t880 + t895;
t694 = -pkin(3) * t962 - pkin(8) * t838 + qJD(2) * t848 + t733;
t846 = qJDD(2) + t856;
t803 = t878 * g(3) + t874 * t825;
t897 = qJDD(3) + t803;
t956 = t880 * qJ(3);
t999 = qJDD(2) * pkin(2);
t882 = -t999 - t837 * pkin(8) - t956 - t846 * pkin(3) + (pkin(8) * t949 + t907 * t950) * qJD(1) + t897;
t623 = t877 * t694 + t873 * t882;
t771 = pkin(4) * t821 - pkin(9) * t823;
t600 = -pkin(4) * t937 - pkin(9) * t938 - t821 * t771 + t623;
t538 = t872 * t596 + t876 * t600;
t724 = pkin(5) * t784 - qJ(6) * t786;
t904 = t734 * qJ(6) + t1053 * t817 - t784 * t724 + t538;
t1022 = t853 * t874 + t957;
t1013 = t874 * t914 + t878 * t899;
t819 = t821 ^ 2;
t820 = t823 ^ 2;
t830 = t878 * t846;
t854 = -t880 - t962;
t792 = t854 * t874 + t830;
t1007 = pkin(1) * t792;
t1005 = pkin(5) * t876;
t967 = t846 * t874;
t797 = t854 * t878 - t967;
t755 = t797 * t875 + t839 * t879;
t1004 = pkin(6) * t755;
t951 = t870 + t871;
t841 = t951 * qJDD(1);
t844 = t951 * t881;
t778 = t841 * t875 + t844 * t879;
t1003 = pkin(6) * t778;
t1002 = pkin(7) * t792;
t1000 = qJ(6) * t876;
t768 = t938 - t977;
t983 = t768 * t873;
t982 = t768 * t877;
t980 = t817 * t786;
t976 = t824 * t874;
t975 = t824 * t878;
t971 = t839 * t874;
t622 = t873 * t694 - t877 * t882;
t599 = t938 * pkin(4) - t937 * pkin(9) + t823 * t771 + t622;
t961 = t872 * t599;
t960 = t873 * t677;
t959 = t876 * t599;
t958 = t877 * t677;
t537 = -t876 * t596 + t872 * t600;
t955 = t1019 - t1010;
t954 = pkin(1) * t839 + pkin(7) * t797;
t953 = pkin(1) * t844 + pkin(7) * t841;
t952 = t844 - t880;
t942 = qJDD(1) * t875;
t933 = t875 * t977;
t932 = t879 * t977;
t927 = pkin(4) * t873 + qJ(3);
t926 = -qJ(6) * t872 - pkin(4);
t487 = t872 * t537 + t876 * t538;
t729 = t803 * t874 - t878 * t925;
t789 = -t849 * t875 - t879 * t850;
t922 = pkin(4) * t877 + t1008;
t919 = t875 * t856;
t918 = t879 * t856;
t916 = t786 * t724 + qJDD(6) + t537;
t843 = -t875 * t881 + t940;
t915 = -pkin(6) * t843 - g(3) * t875;
t911 = t873 * t921;
t910 = t873 * t920;
t909 = t877 * t921;
t908 = t877 * t920;
t906 = pkin(2) * t874 - qJ(3) * t878;
t905 = t837 + t928;
t486 = -t537 * t876 + t538 * t872;
t552 = -t622 * t877 + t623 * t873;
t553 = t873 * t622 + t877 * t623;
t728 = t803 * t878 + t874 * t925;
t903 = t836 * t878 + t971;
t788 = t849 * t879 - t850 * t875;
t896 = -qJD(2) * t823 - t923;
t893 = -t734 * pkin(5) + t916;
t888 = -t686 * pkin(5) - t1069 + t599;
t887 = t1053 * t786 - t888;
t886 = t1009 * t950 + t890;
t885 = t874 * t894 + t897;
t884 = t885 - t999;
t883 = qJ(3) * t905 + t886;
t852 = t880 - t963;
t842 = t879 * t881 + t942;
t833 = t906 * qJDD(1);
t829 = t951 * t945;
t818 = -pkin(6) * t842 + g(3) * t879;
t810 = -t820 + t937;
t809 = t819 - t937;
t808 = qJDD(2) * t875 + t829 * t879;
t807 = t837 * t878 - t870 * t945;
t806 = -qJDD(2) * t879 + t829 * t875;
t805 = -t838 * t874 - t871 * t945;
t802 = -t820 - t937;
t799 = -t852 * t874 + t830;
t793 = t852 * t878 + t967;
t791 = t905 * t874;
t790 = (t838 - t929) * t878;
t779 = t841 * t879 - t844 * t875;
t776 = pkin(6) * t779;
t773 = t820 - t819;
t767 = t807 * t879 - t919;
t766 = t805 * t879 + t919;
t765 = t807 * t875 + t918;
t764 = t805 * t875 - t918;
t763 = t799 * t879 + t874 * t942;
t762 = t799 * t875 - t874 * t940;
t761 = -t937 - t819;
t758 = t797 * t879 - t839 * t875;
t750 = pkin(6) * t758;
t746 = t909 - t910;
t745 = -t911 - t908;
t744 = -t975 + t1051;
t743 = -t976 - t1002;
t741 = -t819 - t820;
t737 = -t884 + t956;
t736 = -t925 + t1052;
t735 = t803 - t1007;
t723 = qJ(3) * t952 + t884;
t720 = pkin(2) * t952 + t895;
t719 = t883 + t1001;
t718 = t809 * t877 + t983;
t717 = -t810 * t873 + t1027;
t716 = -t809 * t873 + t982;
t715 = -t810 * t877 - t1030;
t714 = -t802 * t873 + t982;
t713 = t802 * t877 + t983;
t712 = t740 - t921;
t707 = (0.2e1 * qJD(4) - qJD(2)) * t823 + t923;
t704 = t877 * t740 + t910;
t703 = -t873 * t740 + t908;
t702 = -t873 * t739 - t909;
t701 = -t877 * t739 + t911;
t700 = (t838 + t839) * pkin(2) + t883;
t699 = t1001 + (t836 + t905) * qJ(3) + t886;
t698 = t729 * t879 - t824 * t875;
t697 = t729 * t875 + t824 * t879;
t696 = t761 * t877 - t1030;
t695 = t761 * t873 + t1027;
t678 = -t1007 + (-t854 - t880) * qJ(3) + (-qJDD(2) - t846) * pkin(2) + t885;
t676 = -t1052 - qJ(3) * t847 + (-t851 + t880) * pkin(2) - t895;
t675 = -t745 * t874 + t746 * t878;
t669 = t733 * t878 - t737 * t874;
t668 = t733 * t874 + t737 * t878;
t667 = -pkin(2) * t974 + t699 * t878 - t1051;
t666 = qJ(3) * t970 - t700 * t874 - t1002;
t665 = -t720 * t874 + t723 * t878;
t664 = t784 * t946 + t891;
t658 = -t686 + t980;
t655 = -t716 * t874 + t718 * t878;
t654 = -t715 * t874 + t717 * t878;
t645 = t713 * t874 + t714 * t878;
t644 = -t713 * t878 + t714 * t874;
t643 = t712 * t873 + t877 * t896;
t642 = -t1048 * t873 - t707 * t877;
t641 = -t712 * t877 + t873 * t896;
t640 = -t1048 * t877 + t707 * t873;
t627 = -t703 * t874 + t704 * t878;
t626 = -t701 * t874 + t702 * t878;
t625 = t695 * t874 + t696 * t878;
t624 = -t695 * t878 + t696 * t874;
t617 = t669 * t879 - t719 * t875;
t616 = t669 * t875 + t719 * t879;
t603 = -pkin(1) * t668 - pkin(2) * t737 - qJ(3) * t733;
t602 = -t1048 * t875 + t645 * t879;
t601 = t1048 * t879 + t645 * t875;
t595 = -pkin(8) * t713 + qJ(3) * t1048 + t958;
t592 = t625 * t879 - t707 * t875;
t591 = t625 * t875 + t707 * t879;
t590 = -pkin(7) * t668 - t719 * t906;
t589 = -pkin(8) * t695 + qJ(3) * t707 + t960;
t580 = t641 * t874 + t643 * t878;
t579 = -t640 * t874 + t642 * t878;
t578 = -t641 * t878 + t643 * t874;
t569 = -pkin(8) * t714 + t1008 * t1048 - t960;
t566 = t658 * t873 + t1067;
t565 = -t664 * t873 - t1099;
t564 = -t658 * t877 + t1068;
t563 = t664 * t877 - t1100;
t562 = -pkin(8) * t696 + t1008 * t707 + t958;
t561 = t656 * t873 + t1067;
t560 = -t1020 * t873 + t1099;
t559 = -t656 * t877 + t1068;
t558 = t1020 * t877 + t1100;
t551 = t580 * t879 - t741 * t875;
t550 = t580 * t875 + t741 * t879;
t541 = t959 - t1104;
t540 = t961 - t1073;
t539 = -pkin(8) * t552 + qJ(3) * t677;
t535 = -pkin(1) * t644 - qJ(3) * t714 + t1008 * t713 - t623;
t534 = (pkin(5) * t817 - (2 * qJD(6))) * t786 + t888;
t533 = -pkin(8) * t641 + qJ(3) * t741 - t552;
t532 = -pkin(4) * t1061 - pkin(5) * t1021 + qJ(6) * t657;
t531 = -pkin(8) * t553 + t1008 * t677;
t530 = -pkin(1) * t624 - qJ(3) * t696 + t1008 * t695 - t622;
t529 = -pkin(8) * t643 + t1008 * t741 - t553;
t524 = qJ(6) * t1010 - t893;
t523 = t564 * t874 + t566 * t878;
t522 = t563 * t874 + t565 * t878;
t521 = -t564 * t878 + t566 * t874;
t520 = -t563 * t878 + t565 * t874;
t519 = -pkin(5) * t1010 + t904;
t518 = t538 - t1108;
t517 = t559 * t874 + t561 * t878;
t516 = t558 * t874 + t560 * t878;
t515 = -t559 * t878 + t561 * t874;
t514 = -t558 * t878 + t560 * t874;
t513 = t537 - t1075;
t512 = (-t658 - t980) * pkin(5) + t887;
t511 = -pkin(5) * t980 + t1069 + t887;
t508 = t552 * t874 + t553 * t878;
t507 = -t552 * t878 + t553 * t874;
t506 = -pkin(7) * t644 - t569 * t874 + t595 * t878;
t505 = qJ(6) * t955 + t893;
t500 = pkin(5) * t955 + t904;
t499 = -pkin(1) * t578 - qJ(3) * t643 + t1008 * t641;
t498 = -pkin(7) * t624 - t562 * t874 + t589 * t878;
t497 = t508 * t879 - t677 * t875;
t496 = t508 * t875 + t677 * t879;
t495 = t523 * t879 - t1066;
t494 = t522 * t879 - t1098;
t493 = t523 * t875 + t1065;
t492 = t522 * t875 + t1097;
t491 = t517 * t879 - t1066;
t490 = t516 * t879 + t1098;
t489 = t517 * t875 + t1065;
t488 = t516 * t875 - t1097;
t485 = -t1075 + (-t1018 - t1010) * qJ(6) + (-t1024 - t734) * pkin(5) + t916;
t484 = -t1000 * t658 - t512 * t872 - t1073;
t483 = -pkin(5) * t993 + t511 * t876 + t1104;
t478 = t1108 - qJ(6) * t1023 + (t722 + t1010) * pkin(5) - t904;
t477 = t487 * t877 + t599 * t873;
t476 = t487 * t873 - t599 * t877;
t475 = -t486 - t1087;
t474 = t519 * t876 - t524 * t872;
t473 = t519 * t872 + t524 * t876;
t472 = -pkin(7) * t578 - t529 * t874 + t533 * t878;
t471 = -pkin(8) * t563 - t518 * t873 + t541 * t877 + t1101;
t470 = -pkin(8) * t559 - t513 * t873 + t540 * t877 + t1070;
t469 = -pkin(7) * t507 - t531 * t874 + t539 * t878;
t468 = -t500 * t872 + t505 * t876 - t1087;
t467 = -pkin(8) * t565 - t877 * t518 - t873 * t541 + t1094;
t466 = -pkin(1) * t507 - qJ(3) * t553 + t1008 * t552;
t465 = -pkin(8) * t561 - t877 * t513 - t873 * t540 + t1062;
t464 = t474 * t877 + t534 * t873;
t463 = t474 * t873 - t534 * t877;
t462 = -pkin(1) * t520 + pkin(4) * t664 - qJ(3) * t565 + t1008 * t563 - t1105 + t961;
t461 = t1061 * t927 + t877 * t475 - t1102;
t460 = -pkin(1) * t515 - pkin(4) * t656 - qJ(3) * t561 + t1008 * t559 + t1074 - t959;
t459 = -pkin(9) * t473 + (pkin(5) * t872 - t1000) * t534;
t458 = t476 * t874 + t477 * t878;
t457 = -t476 * t878 + t477 * t874;
t456 = t1061 * t922 - t873 * t475 - t1103;
t455 = -pkin(8) * t564 + t484 * t877 - t485 * t873 + t1070;
t454 = -pkin(8) * t558 - t478 * t873 + t483 * t877 - t1101;
t453 = -pkin(4) * t473 - pkin(5) * t524 - qJ(6) * t519;
t452 = -pkin(8) * t566 - t873 * t484 - t877 * t485 + t1062;
t451 = -pkin(1) * t521 - qJ(3) * t566 + t1008 * t564 + t876 * t512 + t658 * t926 + t1074;
t450 = -pkin(8) * t560 - t877 * t478 - t873 * t483 - t1094;
t449 = -pkin(1) * t514 + t1105 - qJ(3) * t560 + t872 * t511 + (pkin(4) + t1005) * t1020 + t1008 * t558;
t448 = qJ(3) * t1061 + t468 * t877 - t532 * t873 - t1102;
t447 = t1114 + t487;
t446 = t1008 * t1061 - t873 * t468 - t877 * t532 - t1103;
t445 = t458 * t879 - t486 * t875;
t444 = t458 * t875 + t486 * t879;
t443 = t876 * t500 + t872 * t505 + t1114;
t442 = t463 * t874 + t464 * t878;
t441 = -t463 * t878 + t464 * t874;
t440 = -pkin(7) * t520 - t467 * t874 + t471 * t878;
t439 = -pkin(7) * t515 - t465 * t874 + t470 * t878;
t438 = -pkin(8) * t476 + (-pkin(9) * t877 + t927) * t486;
t437 = t442 * t879 - t473 * t875;
t436 = t442 * t875 + t473 * t879;
t435 = -pkin(8) * t477 + (pkin(9) * t873 + t922) * t486;
t434 = -t456 * t874 + t461 * t878 - t1115;
t433 = -pkin(7) * t521 - t452 * t874 + t455 * t878;
t432 = -pkin(7) * t514 - t450 * t874 + t454 * t878;
t431 = -t446 * t874 + t448 * t878 - t1115;
t430 = -pkin(1) * t457 - pkin(4) * t599 + pkin(9) * t487 - qJ(3) * t477 + t1008 * t476;
t429 = -pkin(8) * t463 + qJ(3) * t473 - t453 * t873 + t459 * t877;
t428 = -pkin(8) * t464 + t1008 * t473 - t877 * t453 - t873 * t459;
t427 = -pkin(7) * t457 - t435 * t874 + t438 * t878;
t426 = -pkin(1) * t441 + pkin(9) * t474 - qJ(3) * t464 + t1008 * t463 + (t926 - t1005) * t534;
t425 = -pkin(7) * t441 - t428 * t874 + t429 * t878;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t842, -t843, 0, t789, 0, 0, 0, 0, 0, 0, t758, -t759, t779, t698, 0, 0, 0, 0, 0, 0, t758, t779, t759, t617, 0, 0, 0, 0, 0, 0, t592, t602, t551, t497, 0, 0, 0, 0, 0, 0, t491, t494, t1111, t445, 0, 0, 0, 0, 0, 0, t495, t1111, t490, t437; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t843, -t842, 0, t788, 0, 0, 0, 0, 0, 0, t755, -t756, t778, t697, 0, 0, 0, 0, 0, 0, t755, t778, t756, t616, 0, 0, 0, 0, 0, 0, t591, t601, t550, t496, 0, 0, 0, 0, 0, 0, t489, t492, t1112, t444, 0, 0, 0, 0, 0, 0, t493, t1112, t488, t436; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t792, -t794, 0, -t728, 0, 0, 0, 0, 0, 0, t792, 0, t794, t668, 0, 0, 0, 0, 0, 0, t624, t644, t578, t507, 0, 0, 0, 0, 0, 0, t515, t520, t1090, t457, 0, 0, 0, 0, 0, 0, t521, t1090, t514, t441; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t843, 0, -t842, 0, t915, -t818, -t788, -pkin(6) * t788, t767, -t1041, t763, t766, -t1039, t808, -t735 * t875 + t743 * t879 - t1004, -t736 * t875 + t744 * t879 + t1107, t728 * t879 - t1003, -pkin(6) * t697 - (pkin(1) * t875 - pkin(7) * t879) * t728, t767, t763, t1041, t808, t1039, t766, t666 * t879 - t678 * t875 - t1004, t665 * t879 - t833 * t875 - t1003, t667 * t879 - t676 * t875 - t1107, -pkin(6) * t616 + t590 * t879 - t603 * t875, t627 * t879 - t933, t579 * t879 - t773 * t875, t654 * t879 - t712 * t875, t626 * t879 + t933, t655 * t879 - t875 * t896, t879 * t675 + t875 * t938, -pkin(6) * t591 + t498 * t879 - t530 * t875, -pkin(6) * t601 + t506 * t879 - t535 * t875, -pkin(6) * t550 + t472 * t879 - t499 * t875, -pkin(6) * t496 - t466 * t875 + t469 * t879, t1037, -t1119, t1092, t1057, t1121, t1056, -pkin(6) * t489 + t439 * t879 - t460 * t875, -pkin(6) * t492 + t440 * t879 - t462 * t875, t434 * t879 - t447 * t875 - t1117, -pkin(6) * t444 + t427 * t879 - t430 * t875, t1037, t1092, t1119, t1056, -t1121, t1057, -pkin(6) * t493 + t433 * t879 - t451 * t875, t431 * t879 - t443 * t875 - t1117, -pkin(6) * t488 + t432 * t879 - t449 * t875, -pkin(6) * t436 + t425 * t879 - t426 * t875; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t842, 0, t843, 0, t818, t915, t789, pkin(6) * t789, t765, -t1042, t762, t764, -t1040, t806, t735 * t879 + t743 * t875 + t750, t736 * t879 + t744 * t875 - t1106, t728 * t875 + t776, pkin(6) * t698 - (-pkin(1) * t879 - pkin(7) * t875) * t728, t765, t762, t1042, t806, t1040, t764, t666 * t875 + t678 * t879 + t750, t665 * t875 + t833 * t879 + t776, t667 * t875 + t676 * t879 + t1106, pkin(6) * t617 + t590 * t875 + t603 * t879, t627 * t875 + t932, t579 * t875 + t773 * t879, t654 * t875 + t712 * t879, t626 * t875 - t932, t655 * t875 + t879 * t896, t875 * t675 - t879 * t938, pkin(6) * t592 + t498 * t875 + t530 * t879, pkin(6) * t602 + t506 * t875 + t535 * t879, pkin(6) * t551 + t472 * t875 + t499 * t879, pkin(6) * t497 + t466 * t879 + t469 * t875, t1038, -t1120, t1093, t1059, t1122, t1058, pkin(6) * t491 + t439 * t875 + t460 * t879, pkin(6) * t494 + t440 * t875 + t462 * t879, t434 * t875 + t447 * t879 + t1118, pkin(6) * t445 + t427 * t875 + t430 * t879, t1038, t1093, t1120, t1058, -t1122, t1059, pkin(6) * t495 + t433 * t875 + t451 * t879, t431 * t875 + t443 * t879 + t1118, pkin(6) * t490 + t432 * t875 + t449 * t879, pkin(6) * t437 + t425 * t875 + t426 * t879; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t849, t850, 0, 0, t791, t903, t793, t790, t1022, 0, t954 + t975, -pkin(1) * t836 - t1077 - t976, t729 + t953, pkin(1) * t824 + pkin(7) * t729, t791, t793, -t903, 0, -t1022, t790, qJ(3) * t971 + t700 * t878 + t954, t720 * t878 + t723 * t874 + t953, t1077 + t874 * t699 + (pkin(1) + t1006) * t836, pkin(7) * t669 + (pkin(1) - t907) * t719, t703 * t878 + t704 * t874, t640 * t878 + t642 * t874, t715 * t878 + t717 * t874, t701 * t878 + t702 * t874, t716 * t878 + t718 * t874, t745 * t878 + t746 * t874, pkin(1) * t707 + pkin(7) * t625 + t562 * t878 + t589 * t874, pkin(1) * t1048 + pkin(7) * t645 + t569 * t878 + t595 * t874, pkin(1) * t741 + pkin(7) * t580 + t529 * t878 + t533 * t874, pkin(1) * t677 + pkin(7) * t508 + t531 * t878 + t539 * t874, t1013, t1110, t1078, t1031, -t1116, t1033, pkin(7) * t517 + t465 * t878 + t470 * t874 + t1076, pkin(7) * t522 + t467 * t878 + t471 * t874 + t1109, t456 * t878 + t461 * t874 + t1113, pkin(1) * t486 + pkin(7) * t458 + t435 * t878 + t438 * t874, t1013, t1078, -t1110, t1033, t1116, t1031, pkin(7) * t523 + t452 * t878 + t455 * t874 + t1076, t446 * t878 + t448 * t874 + t1113, pkin(7) * t516 + t450 * t878 + t454 * t874 - t1109, pkin(1) * t473 + pkin(7) * t442 + t428 * t878 + t429 * t874;];
tauB_reg  = t1;
