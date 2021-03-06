% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5PRPRR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5PRPRR4_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRR4_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRR4_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_invdynm_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:51:56
% EndTime: 2019-12-05 15:52:13
% DurationCPUTime: 17.44s
% Computational Cost: add. (71718->653), mult. (134283->984), div. (0->0), fcn. (99159->12), ass. (0->457)
t1057 = sin(qJ(2));
t1047 = sin(pkin(10));
t1050 = cos(pkin(10));
t1048 = sin(pkin(9));
t1051 = cos(pkin(9));
t1019 = t1051 * g(1) + t1048 * g(2);
t1060 = cos(qJ(2));
t1018 = t1048 * g(1) - t1051 * g(2);
t1049 = sin(pkin(5));
t1052 = cos(pkin(5));
t1154 = g(3) - qJDD(1);
t1174 = -t1052 * t1018 + t1049 * t1154;
t934 = -t1057 * t1019 + t1174 * t1060;
t1076 = qJDD(2) * pkin(2) - t934;
t1061 = qJD(2) ^ 2;
t935 = -t1060 * t1019 - t1174 * t1057;
t925 = -t1061 * pkin(2) + t935;
t874 = t1047 * t925 - t1050 * t1076;
t875 = t1047 * t1076 + t1050 * t925;
t1123 = t1047 * t874 + t1050 * t875;
t793 = t1047 * t875 - t1050 * t874;
t1146 = t1060 * t793;
t1099 = t1057 * t1123 + t1146;
t1150 = t1057 * t793;
t1175 = t1060 * t1123 - t1150;
t1012 = t1050 * qJDD(2) - t1047 * t1061;
t985 = t1049 * t1018 + t1052 * t1154;
t978 = -qJDD(3) + t985;
t1109 = -qJ(3) * t1012 - t1047 * t978;
t1011 = t1047 * qJDD(2) + t1050 * t1061;
t930 = qJ(3) * t1011 - t1050 * t978;
t1173 = t1057 * t930 + t1060 * t1109;
t1172 = -t1057 * t1109 + t1060 * t930;
t1083 = t1060 * t1011 + t1057 * t1012;
t1160 = t1083 * t1052;
t958 = t1057 * t1011 - t1060 * t1012;
t1171 = t1048 * t958 - t1051 * t1160;
t1170 = t1048 * t1160 + t1051 * t958;
t1056 = sin(qJ(4));
t1059 = cos(qJ(4));
t866 = -t1061 * pkin(3) + qJDD(2) * pkin(7) + t875;
t967 = t1059 * t978;
t835 = t1056 * t866 + t967;
t836 = -t1056 * t978 + t1059 * t866;
t774 = t1056 * t835 + t1059 * t836;
t1169 = t1048 * t1154;
t1168 = t1051 * t1154;
t1055 = sin(qJ(5));
t1058 = cos(qJ(5));
t1143 = qJD(2) * t1056;
t1001 = -t1058 * qJD(4) + t1055 * t1143;
t1003 = t1055 * qJD(4) + t1058 * t1143;
t1142 = t1003 * t1001;
t1034 = qJD(4) * t1143;
t1130 = t1059 * qJDD(2);
t1008 = -t1034 + t1130;
t996 = -qJDD(5) + t1008;
t1065 = -t996 - t1142;
t1167 = t1055 * t1065;
t1164 = t1058 * t1065;
t1161 = t1083 * t1049;
t1159 = qJD(4) ^ 2;
t994 = t1001 ^ 2;
t995 = t1003 ^ 2;
t1135 = t1059 * qJD(2);
t1030 = -qJD(5) + t1135;
t1028 = t1030 ^ 2;
t1158 = pkin(2) * t793;
t879 = t1057 * t934 + t1060 * t935;
t1157 = pkin(6) * t879;
t1156 = pkin(4) * t1056;
t1155 = pkin(4) * t1059;
t865 = -qJDD(2) * pkin(3) - t1061 * pkin(7) + t874;
t1153 = -pkin(3) * t865 + pkin(7) * t774;
t1114 = -pkin(8) * t1056 - t1155;
t1005 = t1114 * qJD(2);
t811 = -qJDD(4) * pkin(4) - t1159 * pkin(8) + t967 + (qJD(2) * t1005 + t866) * t1056;
t1152 = t1055 * t811;
t931 = t996 - t1142;
t1151 = t1055 * t931;
t861 = t1056 * t865;
t1149 = t1057 * t985;
t1148 = t1058 * t811;
t1147 = t1058 * t931;
t862 = t1059 * t865;
t1145 = t1060 * t985;
t1144 = qJD(2) * qJD(4);
t1141 = t1030 * t1055;
t1140 = t1030 * t1058;
t1043 = t1056 ^ 2;
t1139 = t1043 * t1061;
t1029 = t1056 * t1061 * t1059;
t1020 = qJDD(4) + t1029;
t1137 = t1056 * t1020;
t1021 = qJDD(4) - t1029;
t1136 = t1056 * t1021;
t1134 = t1059 * t1021;
t1044 = t1059 ^ 2;
t1133 = t1043 + t1044;
t1132 = t1049 * qJDD(2);
t1131 = t1056 * qJDD(2);
t1035 = qJD(4) * t1135;
t1006 = 0.2e1 * t1035 + t1131;
t1025 = -t1139 - t1159;
t977 = -t1056 * t1025 - t1134;
t1129 = -pkin(3) * t1006 + pkin(7) * t977 + t861;
t1009 = -0.2e1 * t1034 + t1130;
t1041 = t1044 * t1061;
t1027 = -t1041 - t1159;
t975 = t1059 * t1027 - t1137;
t1128 = pkin(3) * t1009 + pkin(7) * t975 - t862;
t1127 = t1056 * t1142;
t1126 = t1059 * t1142;
t1125 = t1051 * t1132;
t812 = -t1159 * pkin(4) + qJDD(4) * pkin(8) + t1005 * t1135 + t836;
t1110 = -t1008 + t1034;
t1007 = t1035 + t1131;
t1111 = t1007 + t1035;
t820 = pkin(4) * t1110 - pkin(8) * t1111 + t865;
t761 = t1055 * t812 - t1058 * t820;
t762 = t1055 * t820 + t1058 * t812;
t718 = t1055 * t761 + t1058 * t762;
t1122 = -t1048 * t1018 - t1051 * t1019;
t954 = -t1028 - t994;
t881 = t1055 * t954 + t1164;
t744 = -pkin(4) * t881 + t761;
t778 = -pkin(8) * t881 + t1152;
t882 = t1058 * t954 - t1167;
t1119 = t1058 * qJDD(4) - t1055 * t1007;
t944 = -t1003 * qJD(5) + t1119;
t986 = t1003 * t1030;
t910 = t944 + t986;
t817 = -t1056 * t910 + t1059 * t882;
t1121 = -pkin(3) * t881 + pkin(7) * t817 + t1056 * t778 + t1059 * t744;
t962 = -t995 - t1028;
t888 = t1058 * t962 + t1151;
t746 = -pkin(4) * t888 + t762;
t782 = -pkin(8) * t888 + t1148;
t889 = -t1055 * t962 + t1147;
t1080 = -t1055 * qJDD(4) - t1058 * t1007;
t914 = (qJD(5) - t1030) * t1001 + t1080;
t822 = -t1056 * t914 + t1059 * t889;
t1120 = -pkin(3) * t888 + pkin(7) * t822 + t1056 * t782 + t1059 * t746;
t1013 = t1133 * qJDD(2);
t1016 = t1041 + t1139;
t1118 = pkin(3) * t1016 + pkin(7) * t1013 + t774;
t1117 = t1047 * t1029;
t1116 = t1050 * t1029;
t1115 = -pkin(4) * t811 + pkin(8) * t718;
t1014 = t1060 * qJDD(2) - t1057 * t1061;
t1113 = -pkin(6) * t1014 - t1149;
t1081 = t1057 * qJDD(2) + t1060 * t1061;
t1112 = -pkin(6) * t1081 + t1145;
t717 = t1055 * t762 - t1058 * t761;
t772 = t1056 * t836 - t1059 * t835;
t705 = t1056 * t811 + t1059 * t718;
t681 = t1047 * t705 - t1050 * t717;
t682 = t1047 * t717 + t1050 * t705;
t1108 = t1057 * t682 + t1060 * t681;
t747 = t1047 * t774 - t1050 * t865;
t748 = t1047 * t865 + t1050 * t774;
t1107 = t1057 * t748 + t1060 * t747;
t911 = (-qJD(5) - t1030) * t1003 + t1119;
t945 = -t1001 * qJD(5) - t1080;
t987 = t1001 * t1030;
t913 = t945 - t987;
t843 = t1055 * t913 + t1058 * t911;
t926 = t994 + t995;
t802 = -t1056 * t926 + t1059 * t843;
t841 = t1055 * t911 - t1058 * t913;
t763 = t1047 * t802 - t1050 * t841;
t764 = t1047 * t841 + t1050 * t802;
t1106 = t1057 * t764 + t1060 * t763;
t912 = t945 + t987;
t842 = -t1055 * t912 + t1058 * t910;
t963 = t995 - t994;
t810 = t1056 * t963 + t1059 * t842;
t840 = t1055 * t910 + t1058 * t912;
t765 = t1047 * t810 - t1050 * t840;
t766 = t1047 * t840 + t1050 * t810;
t1105 = t1057 * t766 + t1060 * t765;
t780 = t1047 * t817 - t1050 * t881;
t781 = t1047 * t881 + t1050 * t817;
t1104 = t1057 * t781 + t1060 * t780;
t783 = t1047 * t822 - t1050 * t888;
t784 = t1047 * t888 + t1050 * t822;
t1103 = t1057 * t784 + t1060 * t783;
t984 = -t995 + t1028;
t896 = -t1055 * t984 + t1164;
t827 = t1056 * t913 + t1059 * t896;
t894 = t1058 * t984 + t1167;
t788 = t1047 * t827 - t1050 * t894;
t790 = t1047 * t894 + t1050 * t827;
t1102 = t1057 * t790 + t1060 * t788;
t983 = t994 - t1028;
t897 = t1058 * t983 + t1151;
t909 = -t944 + t986;
t828 = -t1056 * t909 + t1059 * t897;
t895 = t1055 * t983 - t1147;
t789 = t1047 * t828 - t1050 * t895;
t791 = t1047 * t895 + t1050 * t828;
t1101 = t1057 * t791 + t1060 * t789;
t904 = -t1001 * t1140 - t1055 * t944;
t872 = t1059 * t904 - t1127;
t903 = t1001 * t1141 - t1058 * t944;
t803 = t1047 * t872 + t1050 * t903;
t805 = -t1047 * t903 + t1050 * t872;
t1098 = t1057 * t805 + t1060 * t803;
t906 = t1003 * t1141 + t1058 * t945;
t873 = t1059 * t906 + t1127;
t905 = -t1003 * t1140 + t1055 * t945;
t804 = t1047 * t873 - t1050 * t905;
t806 = t1047 * t905 + t1050 * t873;
t1097 = t1057 * t806 + t1060 * t804;
t919 = (t1001 * t1058 - t1003 * t1055) * t1030;
t899 = -t1056 * t996 + t1059 * t919;
t918 = (t1001 * t1055 + t1003 * t1058) * t1030;
t837 = t1047 * t899 - t1050 * t918;
t838 = t1047 * t918 + t1050 * t899;
t1096 = t1057 * t838 + t1060 * t837;
t1017 = -t1041 + t1139;
t956 = -t1056 * t1006 + t1059 * t1009;
t916 = -t1050 * t1017 + t1047 * t956;
t917 = t1047 * t1017 + t1050 * t956;
t1095 = t1057 * t917 + t1060 * t916;
t921 = t1050 * t1009 + t1047 * t975;
t923 = -t1047 * t1009 + t1050 * t975;
t1094 = t1057 * t923 + t1060 * t921;
t922 = -t1050 * t1006 + t1047 * t977;
t924 = t1047 * t1006 + t1050 * t977;
t1093 = t1057 * t924 + t1060 * t922;
t1092 = t1057 * t935 - t1060 * t934;
t1026 = t1041 - t1159;
t974 = t1059 * t1026 - t1136;
t936 = t1047 * t974 - t1050 * t1130;
t938 = t1047 * t1130 + t1050 * t974;
t1091 = t1057 * t938 + t1060 * t936;
t1000 = t1059 * t1020;
t1024 = -t1139 + t1159;
t976 = -t1056 * t1024 + t1000;
t937 = t1047 * t976 - t1050 * t1131;
t939 = t1047 * t1131 + t1050 * t976;
t1090 = t1057 * t939 + t1060 * t937;
t981 = -t1056 * t1008 - t1044 * t1144;
t940 = t1047 * t981 - t1116;
t942 = t1050 * t981 + t1117;
t1089 = t1057 * t942 + t1060 * t940;
t982 = t1059 * t1007 - t1043 * t1144;
t941 = t1047 * t982 + t1116;
t943 = t1050 * t982 - t1117;
t1088 = t1057 * t943 + t1060 * t941;
t960 = t1047 * t1013 + t1050 * t1016;
t961 = t1050 * t1013 - t1047 * t1016;
t1087 = t1057 * t961 + t1060 * t960;
t999 = t1133 * t1144;
t979 = -t1050 * qJDD(4) + t1047 * t999;
t980 = t1047 * qJDD(4) + t1050 * t999;
t1086 = t1057 * t980 + t1060 * t979;
t990 = t1081 * t1052;
t1085 = t1051 * t1014 - t1048 * t990;
t1084 = t1048 * t1014 + t1051 * t990;
t1082 = t1051 * t1018 - t1048 * t1019;
t713 = -pkin(8) * t841 - t717;
t1079 = pkin(7) * t802 + t1056 * t713 + (-pkin(3) - t1155) * t841;
t1078 = pkin(4) * t914 + pkin(8) * t889 + t1152;
t1077 = pkin(4) * t910 + pkin(8) * t882 - t1148;
t704 = t1056 * t718 - t1059 * t811;
t672 = -pkin(7) * t704 + (-pkin(8) * t1059 + t1156) * t717;
t675 = -pkin(3) * t704 - t1115;
t655 = -pkin(2) * t704 + qJ(3) * t682 + t1047 * t672 + t1050 * t675;
t660 = -qJ(3) * t681 - t1047 * t675 + t1050 * t672;
t669 = -t1057 * t681 + t1060 * t682;
t1075 = pkin(6) * t669 + t1057 * t660 + t1060 * t655;
t801 = t1056 * t843 + t1059 * t926;
t695 = -pkin(7) * t801 + t1059 * t713 + t841 * t1156;
t1064 = pkin(4) * t926 + pkin(8) * t843 + t718;
t698 = -pkin(3) * t801 - t1064;
t670 = -pkin(2) * t801 + qJ(3) * t764 + t1047 * t695 + t1050 * t698;
t671 = -qJ(3) * t763 - t1047 * t698 + t1050 * t695;
t721 = -t1057 * t763 + t1060 * t764;
t1074 = pkin(6) * t721 + t1057 * t671 + t1060 * t670;
t816 = t1056 * t882 + t1059 * t910;
t700 = -pkin(7) * t816 - t1056 * t744 + t1059 * t778;
t733 = -pkin(3) * t816 - t1077;
t673 = -pkin(2) * t816 + qJ(3) * t781 + t1047 * t700 + t1050 * t733;
t676 = -qJ(3) * t780 - t1047 * t733 + t1050 * t700;
t727 = -t1057 * t780 + t1060 * t781;
t1073 = pkin(6) * t727 + t1057 * t676 + t1060 * t673;
t821 = t1056 * t889 + t1059 * t914;
t707 = -pkin(7) * t821 - t1056 * t746 + t1059 * t782;
t734 = -pkin(3) * t821 - t1078;
t674 = -pkin(2) * t821 + qJ(3) * t784 + t1047 * t707 + t1050 * t734;
t680 = -qJ(3) * t783 - t1047 * t734 + t1050 * t707;
t728 = -t1057 * t783 + t1060 * t784;
t1072 = pkin(6) * t728 + t1057 * t680 + t1060 * t674;
t684 = qJ(3) * t748 + (-pkin(3) * t1050 - pkin(7) * t1047 - pkin(2)) * t772;
t699 = -qJ(3) * t747 + (pkin(3) * t1047 - pkin(7) * t1050) * t772;
t708 = -t1057 * t747 + t1060 * t748;
t1071 = pkin(6) * t708 + t1057 * t699 + t1060 * t684;
t971 = t1056 * t1027 + t1000;
t807 = -pkin(3) * t971 + t835;
t833 = -pkin(7) * t971 + t861;
t749 = -pkin(2) * t971 + qJ(3) * t923 + t1047 * t833 + t1050 * t807;
t753 = -qJ(3) * t921 - t1047 * t807 + t1050 * t833;
t868 = -t1057 * t921 + t1060 * t923;
t1070 = pkin(6) * t868 + t1057 * t753 + t1060 * t749;
t973 = t1059 * t1025 - t1136;
t808 = -pkin(3) * t973 + t836;
t834 = -pkin(7) * t973 + t862;
t750 = -pkin(2) * t973 + qJ(3) * t924 + t1047 * t834 + t1050 * t808;
t754 = -qJ(3) * t922 - t1047 * t808 + t1050 * t834;
t869 = -t1057 * t922 + t1060 * t924;
t1069 = pkin(6) * t869 + t1057 * t754 + t1060 * t750;
t759 = qJ(3) * t961 - t1047 * t772;
t760 = -qJ(3) * t960 - t1050 * t772;
t902 = -t1057 * t960 + t1060 * t961;
t1068 = pkin(6) * t902 + t1057 * t760 + t1060 * t759;
t1067 = -pkin(6) * t1083 - t1172;
t1066 = pkin(6) * t958 + t1173;
t785 = pkin(2) * t978 + qJ(3) * t1123;
t1063 = pkin(6) * t1175 - qJ(3) * t1150 + t1060 * t785;
t1062 = pkin(7) * t705 + (-pkin(3) + t1114) * t717;
t1036 = t1052 * qJDD(2);
t1023 = t1048 * t1132;
t991 = t1014 * t1052;
t989 = t1014 * t1049;
t988 = t1081 * t1049;
t972 = t1059 * t1024 + t1137;
t970 = t1056 * t1026 + t1134;
t969 = t1111 * t1056;
t968 = t1110 * t1059;
t955 = t1059 * t1006 + t1056 * t1009;
t953 = -t1048 * t991 - t1051 * t1081;
t952 = -t1048 * t1081 + t1051 * t991;
t951 = t958 * t1052;
t948 = t958 * t1049;
t915 = -t1057 * t979 + t1060 * t980;
t908 = t1086 * t1052;
t907 = t1086 * t1049;
t901 = -t1145 + (t1049 * t988 + t1052 * t990) * pkin(6);
t900 = -t1149 + (-t1049 * t989 - t1052 * t991) * pkin(6);
t898 = t1056 * t919 + t1059 * t996;
t893 = t1087 * t1052;
t892 = t1087 * t1049;
t891 = t1048 * t951 - t1051 * t1083;
t890 = -t1048 * t1083 - t1051 * t951;
t886 = -t1057 * t941 + t1060 * t943;
t885 = -t1057 * t940 + t1060 * t942;
t884 = -t1057 * t937 + t1060 * t939;
t883 = -t1057 * t936 + t1060 * t938;
t877 = t879 * t1052;
t876 = t879 * t1049;
t871 = t1056 * t906 - t1126;
t870 = t1056 * t904 + t1126;
t860 = -pkin(1) * t989 + t1049 * t934 + t1052 * t1112;
t859 = pkin(1) * t988 + t1049 * t935 + t1052 * t1113;
t858 = pkin(1) * t991 + t1049 * t1112 - t1052 * t934;
t857 = -pkin(1) * t990 + t1049 * t1113 - t1052 * t935;
t856 = -t1057 * t916 + t1060 * t917;
t855 = -pkin(2) * t1011 - t875;
t854 = pkin(2) * t1012 - t874;
t853 = -t1049 * t969 + t1052 * t1088;
t852 = t1049 * t968 + t1052 * t1089;
t851 = t1049 * t1088 + t1052 * t969;
t850 = t1049 * t1089 - t1052 * t968;
t849 = -t1049 * t972 + t1052 * t1090;
t848 = -t1049 * t970 + t1052 * t1091;
t847 = t1049 * t1090 + t1052 * t972;
t846 = t1049 * t1091 + t1052 * t970;
t845 = t1049 * t985 + t1052 * t1092;
t844 = t1049 * t1092 - t1052 * t985;
t832 = -t1049 * t973 + t1052 * t1093;
t831 = -t1049 * t971 + t1052 * t1094;
t830 = t1049 * t1093 + t1052 * t973;
t829 = t1049 * t1094 + t1052 * t971;
t826 = t1056 * t897 + t1059 * t909;
t825 = t1056 * t896 - t1059 * t913;
t815 = -t1049 * t955 + t1052 * t1095;
t814 = t1049 * t1095 + t1052 * t955;
t809 = t1056 * t842 - t1059 * t963;
t800 = (t1049 * t1161 + t1052 * t1160) * pkin(6) + t1172;
t799 = (t1049 * t948 + t1052 * t951) * pkin(6) + t1173;
t797 = pkin(2) * t921 + t1128;
t796 = pkin(2) * t922 + t1129;
t787 = pkin(1) * t845 + t1049 * t1157;
t786 = -pkin(1) * t844 + t1052 * t1157;
t776 = (-t1049 * t844 - t1052 * t845) * pkin(6);
t775 = -t1057 * t837 + t1060 * t838;
t770 = pkin(1) * t1161 - t1049 * t855 + t1052 * t1066;
t769 = pkin(1) * t948 - t1049 * t854 + t1052 * t1067;
t768 = -pkin(1) * t1160 + t1049 * t1066 + t1052 * t855;
t767 = -pkin(1) * t951 + t1049 * t1067 + t1052 * t854;
t757 = pkin(2) * t960 + t1118;
t756 = -t1049 * t898 + t1052 * t1096;
t755 = t1049 * t1096 + t1052 * t898;
t752 = -t1057 * t804 + t1060 * t806;
t751 = -t1057 * t803 + t1060 * t805;
t740 = t1175 * t1052;
t739 = t1175 * t1049;
t738 = -t1057 * t789 + t1060 * t791;
t737 = -t1057 * t788 + t1060 * t790;
t736 = t1049 * t978 + t1052 * t1099;
t735 = t1049 * t1099 - t1052 * t978;
t732 = -t1049 * t871 + t1052 * t1097;
t731 = -t1049 * t870 + t1052 * t1098;
t730 = t1049 * t1097 + t1052 * t871;
t729 = t1049 * t1098 + t1052 * t870;
t726 = -t1049 * t826 + t1052 * t1101;
t725 = -t1049 * t825 + t1052 * t1102;
t724 = t1049 * t1101 + t1052 * t826;
t723 = t1049 * t1102 + t1052 * t825;
t722 = -t1057 * t765 + t1060 * t766;
t720 = -t1049 * t821 + t1052 * t1103;
t719 = t1049 * t1103 + t1052 * t821;
t715 = -t1049 * t816 + t1052 * t1104;
t714 = t1049 * t1104 + t1052 * t816;
t711 = -t1049 * t809 + t1052 * t1105;
t710 = t1049 * t1105 + t1052 * t809;
t709 = pkin(2) * t747 + t1153;
t706 = -t1057 * t759 + t1060 * t760 + (-t1049 * t892 - t1052 * t893) * pkin(6);
t702 = -t1049 * t801 + t1052 * t1106;
t701 = t1049 * t1106 + t1052 * t801;
t697 = -t1057 * t750 + t1060 * t754 + (-t1049 * t830 - t1052 * t832) * pkin(6);
t696 = -t1057 * t749 + t1060 * t753 + (-t1049 * t829 - t1052 * t831) * pkin(6);
t694 = -t1049 * t772 + t1052 * t1107;
t693 = t1049 * t1107 + t1052 * t772;
t692 = pkin(2) * t783 + t1120;
t691 = -pkin(1) * t892 - t1049 * t757 + t1052 * t1068;
t690 = pkin(1) * t893 + t1049 * t1068 + t1052 * t757;
t689 = pkin(2) * t780 + t1121;
t688 = -pkin(1) * t830 - t1049 * t796 + t1052 * t1069;
t687 = -pkin(1) * t829 - t1049 * t797 + t1052 * t1070;
t686 = pkin(1) * t832 + t1049 * t1069 + t1052 * t796;
t685 = pkin(1) * t831 + t1049 * t1070 + t1052 * t797;
t683 = -qJ(3) * t1146 - t1057 * t785 + (-t1049 * t735 - t1052 * t736) * pkin(6);
t679 = pkin(2) * t763 + t1079;
t678 = -pkin(1) * t735 - t1049 * t1158 + t1052 * t1063;
t677 = pkin(1) * t736 + t1049 * t1063 + t1052 * t1158;
t668 = -t1049 * t704 + t1052 * t1108;
t667 = t1049 * t1108 + t1052 * t704;
t666 = pkin(2) * t681 + t1062;
t665 = -t1057 * t684 + t1060 * t699 + (-t1049 * t693 - t1052 * t694) * pkin(6);
t664 = -t1057 * t674 + t1060 * t680 + (-t1049 * t719 - t1052 * t720) * pkin(6);
t663 = -t1057 * t673 + t1060 * t676 + (-t1049 * t714 - t1052 * t715) * pkin(6);
t662 = -pkin(1) * t693 - t1049 * t709 + t1052 * t1071;
t661 = pkin(1) * t694 + t1049 * t1071 + t1052 * t709;
t659 = -pkin(1) * t719 - t1049 * t692 + t1052 * t1072;
t658 = pkin(1) * t720 + t1049 * t1072 + t1052 * t692;
t657 = -pkin(1) * t714 - t1049 * t689 + t1052 * t1073;
t656 = pkin(1) * t715 + t1049 * t1073 + t1052 * t689;
t654 = -t1057 * t670 + t1060 * t671 + (-t1049 * t701 - t1052 * t702) * pkin(6);
t653 = -pkin(1) * t701 - t1049 * t679 + t1052 * t1074;
t652 = pkin(1) * t702 + t1049 * t1074 + t1052 * t679;
t651 = -t1057 * t655 + t1060 * t660 + (-t1049 * t667 - t1052 * t668) * pkin(6);
t650 = -pkin(1) * t667 - t1049 * t666 + t1052 * t1075;
t649 = pkin(1) * t668 + t1049 * t1075 + t1052 * t666;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t1169, -t1168, -t1082, -qJ(1) * t1082, 0, 0, t1085, 0, t953, t1023, -qJ(1) * t952 - t1048 * t860 + t1051 * t900, qJ(1) * t1084 - t1048 * t859 + t1051 * t901, -t1048 * t877 - t1051 * t1092, t1051 * t776 - t1048 * t786 - qJ(1) * (t1048 * t879 + t1051 * t845), 0, 0, -t1170, 0, t891, t1023, -qJ(1) * t890 - t1048 * t769 + t1051 * t799, -qJ(1) * t1171 - t1048 * t770 + t1051 * t800, -t1048 * t740 - t1051 * t1099, t1051 * t683 - t1048 * t678 - qJ(1) * (t1048 * t1175 + t1051 * t736), -t1048 * t853 + t1051 * t886, -t1048 * t815 + t1051 * t856, -t1048 * t849 + t1051 * t884, -t1048 * t852 + t1051 * t885, -t1048 * t848 + t1051 * t883, -t1048 * t908 + t1051 * t915, t1051 * t696 - t1048 * t687 - qJ(1) * (t1048 * t868 + t1051 * t831), t1051 * t697 - t1048 * t688 - qJ(1) * (t1048 * t869 + t1051 * t832), t1051 * t706 - t1048 * t691 - qJ(1) * (t1048 * t902 + t1051 * t893), t1051 * t665 - t1048 * t662 - qJ(1) * (t1048 * t708 + t1051 * t694), -t1048 * t732 + t1051 * t752, -t1048 * t711 + t1051 * t722, -t1048 * t725 + t1051 * t737, -t1048 * t731 + t1051 * t751, -t1048 * t726 + t1051 * t738, -t1048 * t756 + t1051 * t775, t1051 * t663 - t1048 * t657 - qJ(1) * (t1048 * t727 + t1051 * t715), t1051 * t664 - t1048 * t659 - qJ(1) * (t1048 * t728 + t1051 * t720), t1051 * t654 - t1048 * t653 - qJ(1) * (t1048 * t721 + t1051 * t702), t1051 * t651 - t1048 * t650 - qJ(1) * (t1048 * t669 + t1051 * t668); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t1168, -t1169, t1122, qJ(1) * t1122, 0, 0, t1084, 0, t952, -t1125, qJ(1) * t953 + t1048 * t900 + t1051 * t860, -qJ(1) * t1085 + t1048 * t901 + t1051 * t859, -t1048 * t1092 + t1051 * t877, t1048 * t776 + t1051 * t786 + qJ(1) * (-t1048 * t845 + t1051 * t879), 0, 0, -t1171, 0, t890, -t1125, qJ(1) * t891 + t1048 * t799 + t1051 * t769, qJ(1) * t1170 + t1048 * t800 + t1051 * t770, -t1048 * t1099 + t1051 * t740, t1048 * t683 + t1051 * t678 + qJ(1) * (-t1048 * t736 + t1051 * t1175), t1048 * t886 + t1051 * t853, t1048 * t856 + t1051 * t815, t1048 * t884 + t1051 * t849, t1048 * t885 + t1051 * t852, t1048 * t883 + t1051 * t848, t1048 * t915 + t1051 * t908, t1048 * t696 + t1051 * t687 + qJ(1) * (-t1048 * t831 + t1051 * t868), t1048 * t697 + t1051 * t688 + qJ(1) * (-t1048 * t832 + t1051 * t869), t1048 * t706 + t1051 * t691 + qJ(1) * (-t1048 * t893 + t1051 * t902), t1048 * t665 + t1051 * t662 + qJ(1) * (-t1048 * t694 + t1051 * t708), t1048 * t752 + t1051 * t732, t1048 * t722 + t1051 * t711, t1048 * t737 + t1051 * t725, t1048 * t751 + t1051 * t731, t1048 * t738 + t1051 * t726, t1048 * t775 + t1051 * t756, t1048 * t663 + t1051 * t657 + qJ(1) * (-t1048 * t715 + t1051 * t727), t1048 * t664 + t1051 * t659 + qJ(1) * (-t1048 * t720 + t1051 * t728), t1048 * t654 + t1051 * t653 + qJ(1) * (-t1048 * t702 + t1051 * t721), t1048 * t651 + t1051 * t650 + qJ(1) * (-t1048 * t668 + t1051 * t669); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t1018, t1019, 0, 0, 0, 0, t988, 0, t989, t1036, t858, t857, t876, t787, 0, 0, t1161, 0, -t948, t1036, t767, t768, t739, t677, t851, t814, t847, t850, t846, t907, t685, t686, t690, t661, t730, t710, t723, t729, t724, t755, t656, t658, t652, t649; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1154, -t1018, 0, 0, 0, t1014, 0, -t1081, 0, t900, t901, -t1092, t776, 0, 0, -t958, 0, -t1083, 0, t799, t800, -t1099, t683, t886, t856, t884, t885, t883, t915, t696, t697, t706, t665, t752, t722, t737, t751, t738, t775, t663, t664, t654, t651; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1154, 0, -t1019, 0, 0, 0, t990, 0, t991, -t1132, t860, t859, t877, t786, 0, 0, t1160, 0, -t951, -t1132, t769, t770, t740, t678, t853, t815, t849, t852, t848, t908, t687, t688, t691, t662, t732, t711, t725, t731, t726, t756, t657, t659, t653, t650; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1018, t1019, 0, 0, 0, 0, t988, 0, t989, t1036, t858, t857, t876, t787, 0, 0, t1161, 0, -t948, t1036, t767, t768, t739, t677, t851, t814, t847, t850, t846, t907, t685, t686, t690, t661, t730, t710, t723, t729, t724, t755, t656, t658, t652, t649; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t1061, 0, 0, -t985, t934, 0, 0, 0, t1012, 0, -t1011, 0, t1109, t930, -t793, -qJ(3) * t793, t943, t917, t939, t942, t938, t980, t753, t754, t760, t699, t806, t766, t790, t805, t791, t838, t676, t680, t671, t660; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1061, 0, qJDD(2), 0, t985, 0, t935, 0, 0, 0, t1011, 0, t1012, 0, -t930, t1109, t1123, t785, t941, t916, t937, t940, t936, t979, t749, t750, t759, t684, t804, t765, t788, t803, t789, t837, t673, t674, t670, t655; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t934, -t935, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t854, t855, 0, t1158, t969, t955, t972, -t968, t970, 0, t797, t796, t757, t709, t871, t809, t825, t870, t826, t898, t689, t692, t679, t666; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t1061, 0, 0, -t978, t874, 0, t982, t956, t976, t981, t974, t999, t833, t834, -t772, -pkin(7) * t772, t873, t810, t827, t872, t828, t899, t700, t707, t695, t672; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1061, 0, qJDD(2), 0, t978, 0, t875, 0, t1029, -t1017, -t1131, -t1029, -t1130, -qJDD(4), t807, t808, 0, -pkin(3) * t772, -t905, -t840, -t894, t903, -t895, -t918, t733, t734, t698, t675; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t874, -t875, 0, 0, t969, t955, t972, -t968, t970, 0, t1128, t1129, t1118, t1153, t871, t809, t825, t870, t826, t898, t1121, t1120, t1079, t1062; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1007, t1009, t1020, -t1035, t1026, t1035, 0, t865, t835, 0, t906, t842, t896, t904, t897, t919, t778, t782, t713, -pkin(8) * t717; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1034, t1006, t1024, t1008, t1021, -t1034, -t865, 0, t836, 0, -t1142, -t963, -t913, t1142, t909, t996, t744, t746, -pkin(4) * t841, -pkin(4) * t717; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1029, t1017, t1131, t1029, t1130, qJDD(4), -t835, -t836, 0, 0, t905, t840, t894, -t903, t895, t918, t1077, t1078, t1064, t1115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t945, t910, t1065, -t987, t983, t987, 0, t811, t761, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t986, t912, t984, t944, -t931, t986, -t811, 0, t762, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1142, t963, t913, -t1142, -t909, -t996, -t761, -t762, 0, 0;];
m_new_reg = t1;
