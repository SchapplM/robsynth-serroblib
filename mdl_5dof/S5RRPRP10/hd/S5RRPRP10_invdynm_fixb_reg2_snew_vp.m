% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RRPRP10
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RRPRP10_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP10_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_invdynm_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:11:20
% EndTime: 2019-12-31 20:11:28
% DurationCPUTime: 8.70s
% Computational Cost: add. (23894->565), mult. (50099->568), div. (0->0), fcn. (28480->6), ass. (0->377)
t986 = cos(qJ(2));
t988 = qJD(1) ^ 2;
t1068 = t986 * t988;
t983 = sin(qJ(2));
t964 = t983 * t1068;
t955 = qJDD(2) - t964;
t1069 = t986 * t955;
t1122 = qJD(2) ^ 2;
t979 = t983 ^ 2;
t972 = t979 * t988;
t960 = -t972 - t1122;
t903 = t960 * t983 + t1069;
t1054 = qJD(1) * qJD(2);
t1038 = t986 * t1054;
t1051 = t983 * qJDD(1);
t945 = 0.2e1 * t1038 + t1051;
t984 = sin(qJ(1));
t987 = cos(qJ(1));
t1159 = pkin(5) * (t903 * t987 - t945 * t984);
t1158 = pkin(5) * (t903 * t984 + t945 * t987);
t1057 = qJD(1) * t986;
t982 = sin(qJ(4));
t985 = cos(qJ(4));
t938 = qJD(2) * t982 + t985 * t1057;
t940 = t985 * qJD(2) - t982 * t1057;
t885 = t940 * t938;
t946 = t1038 + t1051;
t931 = qJDD(4) + t946;
t1141 = t885 - t931;
t1157 = t1141 * pkin(4);
t954 = t964 + qJDD(2);
t1070 = t986 * t954;
t1123 = t986 ^ 2;
t974 = t1123 * t988;
t962 = -t974 - t1122;
t893 = t962 * t983 + t1070;
t1156 = pkin(1) * t893;
t1155 = pkin(6) * t893;
t1154 = pkin(6) * t903;
t1080 = t954 * t983;
t902 = -t962 * t986 + t1080;
t1039 = t983 * t1054;
t971 = t986 * qJDD(1);
t947 = t971 - 0.2e1 * t1039;
t1153 = pkin(5) * (t902 * t987 + t947 * t984);
t1152 = pkin(5) * (t902 * t984 - t947 * t987);
t1089 = t1141 * t985;
t929 = t938 ^ 2;
t1058 = qJD(1) * t983;
t966 = qJD(4) + t1058;
t963 = t966 ^ 2;
t873 = -t963 - t929;
t797 = t873 * t982 - t1089;
t794 = pkin(3) * t797;
t1018 = -t971 + t1039;
t872 = -t938 * qJD(4) + t985 * qJDD(2) + t982 * t1018;
t1011 = t1018 * pkin(2);
t956 = pkin(3) * t1058 - qJD(2) * pkin(7);
t957 = t984 * g(1) - t987 * g(2);
t1013 = -qJDD(1) * pkin(1) - t957;
t1145 = 2 * qJD(3);
t1130 = -pkin(2) * t1039 + t1058 * t1145;
t1019 = t946 + t1038;
t1131 = t1019 * qJ(3);
t990 = t1013 - t1130 - t1131;
t783 = t1011 + t1018 * pkin(7) - t956 * t1058 + (-pkin(3) * t1123 - pkin(6)) * t988 + t990;
t1102 = t986 * g(3);
t1015 = -qJDD(2) * pkin(2) - t1122 * qJ(3) + qJDD(3) + t1102;
t1096 = qJ(3) * t983;
t1112 = pkin(2) * t986;
t1021 = -t1096 - t1112;
t958 = g(1) * t987 + g(2) * t984;
t924 = -pkin(1) * t988 + qJDD(1) * pkin(6) - t958;
t1034 = t988 * t1021 + t924;
t800 = -qJDD(2) * pkin(7) + (t946 - t1038) * pkin(3) + (-pkin(7) * t1068 + t1034) * t983 + t1015;
t739 = t982 * t783 - t985 * t800;
t917 = t966 * t938;
t993 = qJ(5) * t917 + 0.2e1 * qJD(5) * t940 + t1157 + t739;
t716 = qJ(5) * t872 + t993;
t991 = -t716 - t1157;
t1151 = -t794 - t991;
t976 = t983 * g(3);
t999 = -t1122 * pkin(2) + t1034 * t986 - t976;
t992 = qJD(2) * t1145 + t999;
t1150 = (qJDD(2) + t955) * qJ(3) - pkin(2) * t960 + t992;
t959 = -t972 + t1122;
t897 = -t959 * t983 + t1070;
t1149 = -t987 * t1051 + t897 * t984;
t1148 = t984 * t1051 + t897 * t987;
t1079 = t955 * t983;
t895 = -t960 * t986 + t1079;
t1144 = pkin(1) * t895;
t1143 = pkin(6) * t895;
t1142 = pkin(6) * t902;
t1071 = t982 * t1141;
t961 = t974 - t1122;
t900 = -t961 * t986 + t1079;
t1140 = t900 * t984 + t987 * t971;
t1139 = t900 * t987 - t984 * t971;
t1033 = -t982 * qJDD(2) + t985 * t1018;
t1009 = qJD(4) * t940 - t1033;
t918 = t966 * t940;
t834 = t1009 + t918;
t1136 = -pkin(2) * t797 + qJ(3) * t834;
t1132 = -t917 + t872;
t860 = t885 + t931;
t1091 = t860 * t982;
t930 = t940 ^ 2;
t880 = -t930 - t963;
t809 = t880 * t985 - t1091;
t1135 = -pkin(2) * t809 + qJ(3) * t1132;
t836 = (-qJD(4) + t966) * t940 + t1033;
t839 = t917 + t872;
t775 = t836 * t982 - t839 * t985;
t858 = -t929 - t930;
t1134 = -pkin(2) * t775 + qJ(3) * t858;
t740 = t985 * t783 + t982 * t800;
t705 = -t985 * t739 + t740 * t982;
t891 = t961 * t983 + t1069;
t1052 = qJDD(2) * qJ(3);
t1129 = -t1018 * pkin(3) - pkin(7) * t974 + t1052;
t846 = t1034 * t983 + t1015;
t1128 = -pkin(2) * t954 - qJ(3) * t962 + t846;
t1094 = t716 * t985;
t905 = pkin(4) * t966 - qJ(5) * t940;
t1017 = t929 * pkin(4) + qJ(5) * t1009 + t966 * t905 - t740;
t1056 = qJD(5) * t938;
t920 = -0.2e1 * t1056;
t720 = t920 - t1017;
t1124 = pkin(4) * t1009 - t929 * qJ(5) + t940 * t905 + qJDD(5);
t799 = (t1145 + t956) * qJD(2) + t999 + t1129;
t741 = t799 + t1124;
t701 = -pkin(4) * t741 + qJ(5) * t720;
t1010 = qJ(5) * t1094 - t982 * t701;
t1121 = pkin(2) + pkin(7);
t693 = t720 * t982 - t1094;
t1127 = qJ(3) * t741 - t1121 * t693 + t1010;
t1125 = qJ(3) * t799 - t1121 * t705;
t743 = -t775 * t986 + t858 * t983;
t1120 = pkin(1) * t743;
t753 = -t797 * t986 + t834 * t983;
t1119 = pkin(1) * t753;
t756 = t1132 * t983 - t809 * t986;
t1118 = pkin(1) * t756;
t777 = t836 * t985 + t982 * t839;
t1117 = pkin(2) * t777;
t798 = t873 * t985 + t1071;
t1116 = pkin(2) * t798;
t1090 = t860 * t985;
t810 = -t982 * t880 - t1090;
t1115 = pkin(2) * t810;
t1111 = pkin(3) * t705;
t1110 = pkin(3) * t799;
t744 = t775 * t983 + t858 * t986;
t1109 = pkin(5) * (t744 * t984 - t777 * t987);
t754 = t797 * t983 + t834 * t986;
t1108 = pkin(5) * (t754 * t984 - t798 * t987);
t757 = t1132 * t986 + t809 * t983;
t1107 = pkin(5) * (t757 * t984 - t810 * t987);
t1049 = t979 + t1123;
t949 = t1049 * qJDD(1);
t952 = t972 + t974;
t1106 = pkin(5) * (t949 * t984 + t952 * t987);
t1105 = pkin(6) * t743;
t1104 = pkin(6) * t753;
t1103 = pkin(6) * t756;
t1101 = t988 * pkin(6);
t1100 = qJ(3) * t798;
t1098 = qJ(3) * t810;
t1092 = t799 * t982;
t923 = -t1013 + t1101;
t1088 = t923 * t983;
t1087 = t923 * t986;
t1085 = t945 * t983;
t887 = t947 * t986;
t1074 = t966 * t982;
t1073 = t966 * t985;
t1072 = t982 * t716;
t789 = t985 * t799;
t1067 = -pkin(1) * t777 + pkin(6) * t744;
t1066 = -pkin(1) * t798 + pkin(6) * t754;
t1065 = -pkin(1) * t810 + pkin(6) * t757;
t1064 = -pkin(3) * t858 + pkin(7) * t777;
t802 = pkin(7) * t809;
t1063 = t789 - t802;
t1062 = pkin(3) * t834 - pkin(7) * t798;
t1061 = pkin(3) * t1132 - pkin(7) * t810;
t1060 = pkin(1) * t952 + pkin(6) * t949;
t1048 = t983 * t885;
t1047 = t986 * t885;
t1046 = t940 * t1074;
t804 = pkin(3) * t809;
t1045 = -t804 + t740;
t714 = pkin(4) * t716;
t1037 = -pkin(3) * t693 + t714;
t771 = pkin(3) * t775;
t725 = -qJ(3) * t777 + t771;
t792 = pkin(7) * t797;
t1036 = -t792 + t1092;
t906 = t983 * t924 + t1102;
t907 = t986 * t924 - t976;
t842 = t906 * t983 + t986 * t907;
t1035 = -t957 * t984 - t987 * t958;
t1032 = t984 * t964;
t1031 = t987 * t964;
t1030 = t1063 + t1135;
t951 = qJDD(1) * t987 - t984 * t988;
t1029 = -pkin(5) * t951 - g(3) * t984;
t704 = -pkin(4) * t858 + qJ(5) * t836 + t720;
t708 = (t839 + t872) * qJ(5) + t993;
t769 = pkin(7) * t775;
t1028 = -t704 * t982 + t985 * t708 - t769;
t736 = -qJ(5) * t880 + t741;
t780 = -pkin(4) * t1132 - qJ(5) * t860;
t1027 = t985 * t736 - t780 * t982 - t802;
t1026 = -t769 - t705;
t1025 = t739 - t794;
t1024 = -t1062 - t789;
t1023 = -t1061 + t1092;
t843 = t992 + t1052;
t1022 = -pkin(2) * t846 + qJ(3) * t843;
t1020 = pkin(2) * t983 - qJ(3) * t986;
t706 = t982 * t739 + t740 * t985;
t841 = t906 * t986 - t907 * t983;
t888 = t959 * t986 + t1080;
t1016 = t957 * t987 - t958 * t984;
t1014 = t1036 + t1136;
t1012 = pkin(4) * t880 + t1017;
t1008 = t1028 + t1134;
t1007 = t1027 + t1135;
t1006 = t1026 + t1134;
t1004 = -t804 - t1012;
t723 = -pkin(4) * t834 + qJ(5) * t873 - qJD(2) * t956 - t1124 - t1129 - t992;
t1003 = qJ(5) * t1089 - t723 * t982 - t792;
t1002 = t704 * t985 + t708 * t982 + t1064;
t1001 = t736 * t982 + t780 * t985 - t1061;
t1000 = t706 + t1064;
t998 = qJ(5) * t1071 + t723 * t985 - t1062;
t997 = pkin(3) * t741 - qJ(5) * t1072 - t985 * t701;
t994 = t1003 + t1136;
t989 = -t1011 + t923 + t1130;
t953 = t972 - t974;
t950 = qJDD(1) * t984 + t987 * t988;
t941 = t1020 * qJDD(1);
t933 = t1049 * t1054;
t921 = 0.2e1 * t1056;
t919 = -pkin(5) * t950 + g(3) * t987;
t914 = -t930 + t963;
t913 = t929 - t963;
t912 = qJDD(2) * t984 + t933 * t987;
t911 = -t979 * t1054 + t946 * t986;
t910 = -qJDD(2) * t987 + t933 * t984;
t909 = t983 * t1018 - t1123 * t1054;
t908 = t938 * t1073;
t890 = t1019 * t983;
t882 = t930 - t929;
t879 = pkin(5) * (t949 * t987 - t952 * t984);
t877 = t887 - t1085;
t876 = t945 * t986 + t947 * t983;
t870 = t911 * t987 - t1032;
t869 = t909 * t987 + t1032;
t868 = t911 * t984 + t1031;
t867 = t909 * t984 - t1031;
t852 = -t1087 + t1143;
t851 = -t1088 - t1155;
t850 = t908 - t1046;
t849 = (-t938 * t982 - t940 * t985) * t966;
t848 = t877 * t987 + t953 * t984;
t847 = t877 * t984 - t953 * t987;
t845 = t907 + t1144;
t844 = t906 - t1156;
t835 = t1009 - t918;
t831 = pkin(4) * t839;
t828 = pkin(1) * t947 + t1087 - t1142;
t827 = -pkin(1) * t945 - t1088 - t1154;
t824 = qJ(3) * t952 + t846;
t823 = pkin(2) * t952 + t843;
t822 = t940 * t1073 + t872 * t982;
t821 = t1009 * t982 + t908;
t820 = -t985 * t872 + t1046;
t819 = t1009 * t985 - t938 * t1074;
t818 = t989 + t1131;
t817 = t849 * t983 + t931 * t986;
t816 = -t849 * t986 + t931 * t983;
t815 = t914 * t982 + t1089;
t814 = t913 * t982 + t1090;
t813 = -t985 * t913 + t1091;
t812 = t914 * t985 - t1071;
t811 = pkin(1) * t923 + pkin(6) * t842;
t808 = -t1101 + (t1018 - t947) * pkin(2) + t990;
t807 = (t1019 + t945) * qJ(3) + t989;
t801 = t842 + t1060;
t788 = -t1128 + t1156;
t787 = t822 * t983 + t1047;
t786 = -t819 * t983 - t1047;
t785 = -t822 * t986 + t1048;
t784 = t819 * t986 - t1048;
t782 = -t1144 - t1150;
t779 = t843 * t986 + t846 * t983;
t778 = t843 * t983 - t846 * t986;
t776 = t1132 * t982 + t985 * t834;
t774 = t1132 * t985 - t834 * t982;
t768 = t817 * t987 - t850 * t984;
t767 = t817 * t984 + t850 * t987;
t766 = -pkin(2) * t1085 + t807 * t986 - t1143;
t765 = -qJ(3) * t887 - t808 * t983 + t1155;
t764 = -t823 * t983 + t824 * t986;
t763 = t812 * t983 + t839 * t986;
t762 = t814 * t983 - t835 * t986;
t761 = -t812 * t986 + t839 * t983;
t760 = -t814 * t986 - t835 * t983;
t759 = t1154 + t983 * t807 + (pkin(1) + t1112) * t945;
t758 = t1142 + t986 * t808 + (-pkin(1) - t1096) * t947;
t751 = t823 * t986 + t824 * t983 + t1060;
t750 = t774 * t983 + t882 * t986;
t749 = -t774 * t986 + t882 * t983;
t748 = t787 * t987 - t820 * t984;
t747 = t786 * t987 + t821 * t984;
t746 = t787 * t984 + t820 * t987;
t745 = t786 * t984 - t821 * t987;
t734 = t763 * t987 - t815 * t984;
t733 = t762 * t987 - t813 * t984;
t732 = t763 * t984 + t815 * t987;
t731 = t762 * t984 + t813 * t987;
t730 = -pkin(1) * t778 - t1022;
t728 = pkin(5) * (t757 * t987 + t810 * t984);
t726 = pkin(5) * (t754 * t987 + t798 * t984);
t724 = -pkin(6) * t778 - t1020 * t818;
t722 = t750 * t987 - t776 * t984;
t721 = t750 * t984 + t776 * t987;
t718 = pkin(5) * (t744 * t987 + t777 * t984);
t717 = -t1023 - t1115;
t715 = t725 - t831;
t712 = pkin(6) * t779 + (pkin(1) - t1021) * t818;
t711 = -t1024 - t1116;
t710 = -t1045 - t1098;
t709 = -t1025 - t1100;
t703 = -t1030 - t1118;
t702 = -t1014 - t1119;
t700 = t705 * t983 + t799 * t986;
t699 = -t705 * t986 + t799 * t983;
t698 = -t1004 + t921 - t1098;
t697 = -t1100 - t1151;
t696 = -t1001 - t1115;
t695 = -t998 - t1116;
t694 = t720 * t985 + t1072;
t692 = -t1000 - t1117;
t691 = -t1007 - t1118;
t690 = -t994 - t1119;
t689 = -qJ(3) * t706 + t1111;
t688 = t693 * t983 + t741 * t986;
t687 = -t693 * t986 + t741 * t983;
t686 = t710 * t986 - t717 * t983 - t1103;
t685 = t709 * t986 - t711 * t983 - t1104;
t684 = -t1006 - t1120;
t683 = -t1121 * t706 + t1110;
t682 = t710 * t983 + t717 * t986 + t1065;
t681 = t709 * t983 + t711 * t986 + t1066;
t680 = -t1002 - t1117;
t679 = -t692 * t983 + t725 * t986 - t1105;
t678 = -t696 * t983 + t698 * t986 - t1103;
t677 = t692 * t986 + t725 * t983 + t1067;
t676 = -t695 * t983 + t697 * t986 - t1104;
t675 = -t1008 - t1120;
t674 = t696 * t986 + t698 * t983 + t1065;
t673 = t695 * t986 + t697 * t983 + t1066;
t672 = -pkin(1) * t699 - t1125;
t671 = -qJ(3) * t694 - t1037;
t670 = -t680 * t983 + t715 * t986 - t1105;
t669 = t680 * t986 + t715 * t983 + t1067;
t668 = -pkin(6) * t699 - t683 * t983 + t689 * t986;
t667 = -t1121 * t694 + t997;
t666 = -pkin(1) * t706 + pkin(6) * t700 + t683 * t986 + t689 * t983;
t665 = -pkin(1) * t687 - t1127;
t664 = -pkin(6) * t687 - t667 * t983 + t671 * t986;
t663 = -pkin(1) * t694 + pkin(6) * t688 + t667 * t986 + t671 * t983;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t951, 0, -t950, 0, t1029, -t919, -t1016, -pkin(5) * t1016, t870, t848, t1148, t869, -t1139, t912, -t984 * t844 + t987 * t851 + t1152, -t984 * t845 + t987 * t852 + t1158, t841 * t987 - t1106, -pkin(5) * (t842 * t984 + t923 * t987) - (pkin(1) * t984 - pkin(6) * t987) * t841, t912, -t1148, t1139, t870, t848, t869, t764 * t987 - t941 * t984 - t1106, t987 * t765 - t984 * t788 - t1152, t987 * t766 - t984 * t782 - t1158, t987 * t724 - t984 * t730 - pkin(5) * (t779 * t984 + t818 * t987), t748, t722, t734, t747, t733, t768, t685 * t987 - t702 * t984 - t1108, t686 * t987 - t703 * t984 - t1107, t679 * t987 - t684 * t984 - t1109, t987 * t668 - t984 * t672 - pkin(5) * (t700 * t984 - t706 * t987), t748, t722, t734, t747, t733, t768, t676 * t987 - t690 * t984 - t1108, t678 * t987 - t691 * t984 - t1107, t670 * t987 - t675 * t984 - t1109, t987 * t664 - t984 * t665 - pkin(5) * (t688 * t984 - t694 * t987); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t950, 0, t951, 0, t919, t1029, t1035, pkin(5) * t1035, t868, t847, t1149, t867, -t1140, t910, t987 * t844 + t984 * t851 - t1153, t987 * t845 + t984 * t852 - t1159, t841 * t984 + t879, pkin(5) * (t842 * t987 - t923 * t984) - (-pkin(1) * t987 - pkin(6) * t984) * t841, t910, -t1149, t1140, t868, t847, t867, t764 * t984 + t941 * t987 + t879, t984 * t765 + t987 * t788 + t1153, t984 * t766 + t987 * t782 + t1159, t984 * t724 + t987 * t730 + pkin(5) * (t779 * t987 - t818 * t984), t746, t721, t732, t745, t731, t767, t685 * t984 + t702 * t987 + t726, t686 * t984 + t703 * t987 + t728, t679 * t984 + t684 * t987 + t718, t984 * t668 + t987 * t672 + pkin(5) * (t700 * t987 + t706 * t984), t746, t721, t732, t745, t731, t767, t676 * t984 + t690 * t987 + t726, t678 * t984 + t691 * t987 + t728, t670 * t984 + t675 * t987 + t718, t984 * t664 + t987 * t665 + pkin(5) * (t688 * t987 + t694 * t984); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t957, t958, 0, 0, t890, t876, t888, t887, t891, 0, t828, t827, t801, t811, 0, -t888, -t891, t890, t876, t887, t751, t758, t759, t712, t785, t749, t761, t784, t760, t816, t681, t682, t677, t666, t785, t749, t761, t784, t760, t816, t673, t674, t669, t663; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t988, 0, 0, -g(3), -t957, 0, t911, t877, t897, t909, -t900, t933, t851, t852, t841, pkin(6) * t841, t933, -t897, t900, t911, t877, t909, t764, t765, t766, t724, t787, t750, t763, t786, t762, t817, t685, t686, t679, t668, t787, t750, t763, t786, t762, t817, t676, t678, t670, t664; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t988, 0, qJDD(1), 0, g(3), 0, -t958, 0, t964, -t953, -t1051, -t964, -t971, -qJDD(2), t844, t845, 0, pkin(1) * t841, -qJDD(2), t1051, t971, t964, -t953, -t964, t941, t788, t782, t730, t820, t776, t815, -t821, t813, t850, t702, t703, t684, t672, t820, t776, t815, -t821, t813, t850, t690, t691, t675, t665; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t957, t958, 0, 0, t890, t876, t888, t887, t891, 0, t828, t827, t801, t811, 0, -t888, -t891, t890, t876, t887, t751, t758, t759, t712, t785, t749, t761, t784, t760, t816, t681, t682, t677, t666, t785, t749, t761, t784, t760, t816, t673, t674, t669, t663; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t946, t947, t954, -t1038, t961, t1038, 0, -t923, t906, 0, t1038, -t954, -t961, t946, t947, -t1038, t824, -qJ(3) * t947, t807, qJ(3) * t818, t885, t882, t839, -t885, -t835, t931, t709, t710, t725, t689, t885, t882, t839, -t885, -t835, t931, t697, t698, t715, t671; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1039, t945, t959, -t1018, t955, -t1039, t923, 0, t907, 0, -t1039, -t959, -t955, t1039, t945, -t1018, t823, t808, pkin(2) * t945, pkin(2) * t818, -t822, -t774, -t812, t819, -t814, -t849, t711, t717, t692, t683, -t822, -t774, -t812, t819, -t814, -t849, t695, t696, t680, t667; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t964, t953, t1051, t964, t971, qJDD(2), -t906, -t907, 0, 0, qJDD(2), -t1051, -t971, -t964, t953, t964, -t941, t1128, t1150, t1022, -t820, -t776, -t815, t821, -t813, -t850, t1014, t1030, t1006, t1125, -t820, -t776, -t815, t821, -t813, -t850, t994, t1007, t1008, t1127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t1051, -t971, -t964, t953, t964, 0, t846, t843, 0, -t820, -t776, -t815, t821, -t813, -t850, t1036, t1063, t1026, -pkin(7) * t705, -t820, -t776, -t815, t821, -t813, -t850, t1003, t1027, t1028, -pkin(7) * t693 + t1010; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1038, t954, t961, -t946, -t947, t1038, -t846, 0, -t818, 0, -t885, -t882, -t839, t885, t835, -t931, t1025, t1045, -t771, -t1111, -t885, -t882, -t839, t885, t835, -t931, t1151, t920 + t1004, t831 - t771, t1037; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1039, t959, t955, -t1039, -t945, t1018, -t843, t818, 0, 0, t822, t774, t812, -t819, t814, t849, t1024, t1023, t1000, pkin(7) * t706 - t1110, t822, t774, t812, -t819, t814, t849, t998, t1001, t1002, pkin(7) * t694 - t997; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t872, -t834, -t1141, t917, t913, -t917, 0, t799, t739, 0, t872, -t834, -t1141, t917, t913, -t917, qJ(5) * t1141, t736, t708, qJ(5) * t716; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t918, t1132, t914, -t1009, t860, -t918, -t799, 0, t740, 0, t918, t1132, t914, -t1009, t860, -t918, t723, t780, t704, t701; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t885, t882, t839, -t885, -t835, t931, -t739, -t740, 0, 0, t885, t882, t839, -t885, -t835, t931, t991, t921 + t1012, -t831, -t714; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t872, -t834, -t1141, t917, t913, -t917, 0, t741, t716, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t918, t1132, t914, -t1009, t860, -t918, -t741, 0, t720, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t885, t882, t839, -t885, -t835, t931, -t716, -t720, 0, 0;];
m_new_reg = t1;
