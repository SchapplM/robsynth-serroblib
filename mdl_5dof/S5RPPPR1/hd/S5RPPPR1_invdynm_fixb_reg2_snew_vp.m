% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RPPPR1
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RPPPR1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR1_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:55
% EndTime: 2020-01-03 11:21:09
% DurationCPUTime: 14.41s
% Computational Cost: add. (53224->619), mult. (128162->868), div. (0->0), fcn. (82460->10), ass. (0->434)
t1063 = sin(qJ(5));
t1061 = cos(pkin(8));
t1127 = t1061 * qJDD(1);
t1036 = -qJDD(5) + t1127;
t1058 = sin(pkin(8));
t1057 = sin(pkin(9));
t1060 = cos(pkin(9));
t1065 = cos(qJ(5));
t1089 = t1057 * t1065 + t1060 * t1063;
t1082 = t1089 * t1058;
t981 = qJD(1) * t1082;
t1138 = t1058 * t1060;
t1140 = t1057 * t1058;
t983 = (-t1063 * t1140 + t1065 * t1138) * qJD(1);
t1166 = t983 * t981;
t1184 = -t1036 - t1166;
t1193 = t1063 * t1184;
t1192 = t1065 * t1184;
t1064 = sin(qJ(1));
t1059 = sin(pkin(7));
t1062 = cos(pkin(7));
t1066 = cos(qJ(1));
t1029 = t1064 * g(2) - t1066 * g(3);
t1067 = qJD(1) ^ 2;
t1013 = -t1067 * pkin(1) - t1029;
t1030 = t1066 * g(2) + t1064 * g(3);
t1081 = qJDD(1) * pkin(1) - t1030;
t954 = t1059 * t1013 - t1062 * t1081;
t955 = t1062 * t1013 + t1059 * t1081;
t1111 = t1059 * t954 + t1062 * t955;
t885 = t1059 * t955 - t1062 * t954;
t1145 = t1066 * t885;
t1191 = t1064 * t1111 + t1145;
t1148 = t1064 * t885;
t1190 = -t1066 * t1111 + t1148;
t1126 = t1062 * qJDD(1);
t1017 = -t1059 * t1067 + t1126;
t1165 = g(1) - qJDD(2);
t1175 = -qJ(2) * t1017 - t1059 * t1165;
t1128 = t1059 * qJDD(1);
t1016 = t1062 * t1067 + t1128;
t1177 = t1064 * t1016 - t1066 * t1017;
t986 = -qJ(2) * t1016 + t1062 * t1165;
t1189 = -pkin(5) * t1177 + t1064 * t986 - t1066 * t1175;
t1176 = -t1066 * t1016 - t1064 * t1017;
t1188 = pkin(5) * t1176 + t1064 * t1175 + t1066 * t986;
t1168 = pkin(3) * t1061;
t1095 = -qJ(4) * t1058 - t1168;
t936 = -qJDD(1) * pkin(2) - t1067 * qJ(3) + qJDD(3) + t954;
t1187 = -0.2e1 * qJD(1) * qJD(4) * t1058 + t1095 * qJDD(1) + t936;
t1136 = t1061 * qJD(1);
t1037 = -qJD(5) + t1136;
t1144 = t981 * t1037;
t1043 = t1058 * qJDD(1);
t1115 = t1060 * t1043;
t1118 = t1057 * t1043;
t926 = -t981 * qJD(5) - t1063 * t1118 + t1065 * t1115;
t1186 = t926 + t1144;
t1052 = t1058 ^ 2;
t1054 = t1061 ^ 2;
t1135 = t1061 * t1067;
t1008 = (t1052 + t1054) * t1135;
t1114 = t1061 * t1126;
t968 = -t1059 * t1008 + t1114;
t1116 = t1059 * t1127;
t970 = t1062 * t1008 + t1116;
t1182 = t1064 * t970 - t1066 * t968;
t1181 = t1064 * t968 + t1066 * t970;
t1040 = t1061 * t1165;
t1173 = 2 * qJD(3);
t939 = -t1067 * pkin(2) + qJDD(1) * qJ(3) + t955;
t1110 = qJD(1) * t1173 + t939;
t913 = t1110 * t1058 + t1040;
t1112 = t1058 * t1165;
t914 = t1110 * t1061 - t1112;
t846 = t1058 * t913 + t1061 * t914;
t1180 = (qJD(5) + t1037) * t983;
t979 = t981 ^ 2;
t980 = t983 ^ 2;
t1035 = t1037 ^ 2;
t1174 = t1057 ^ 2;
t1085 = t1187 * t1060;
t1087 = -pkin(4) * t1061 - pkin(6) * t1138;
t1137 = t1058 * t1061;
t1167 = pkin(4) * t1052;
t1132 = t1095 * qJD(1) + t1173;
t1092 = t1132 * qJD(1) + t939;
t880 = t1092 * t1061 - t1112;
t809 = t1087 * qJDD(1) + (-t880 + (pkin(6) * t1137 - t1060 * t1167) * t1067) * t1057 + t1085;
t1000 = t1087 * qJD(1);
t1139 = t1057 * t1067;
t828 = t1187 * t1057 + t1060 * t880;
t810 = t1000 * t1136 + (-pkin(6) * t1043 - t1139 * t1167) * t1057 + t828;
t759 = t1063 * t810 - t1065 * t809;
t760 = t1063 * t809 + t1065 * t810;
t711 = t1063 * t760 - t1065 * t759;
t1172 = pkin(4) * t711;
t1076 = qJDD(1) * t1082;
t874 = -t1076 - t1180;
t877 = t926 - t1144;
t819 = t1063 * t874 - t1065 * t877;
t1171 = pkin(4) * t819;
t1169 = pkin(3) * t1058;
t1164 = -pkin(2) * t936 + qJ(3) * t846;
t1161 = t1037 * t983;
t1160 = t1057 * t711;
t1125 = qJDD(4) + t1040;
t878 = t1058 * t1092 + t1125;
t1159 = t1057 * t878;
t1120 = t1060 * t1139;
t1018 = t1052 * t1120;
t992 = -t1018 + t1127;
t1158 = t1057 * t992;
t993 = -t1018 - t1127;
t1157 = t1057 * t993;
t931 = t1058 * t936;
t1156 = t1059 * t936;
t1155 = t1060 * t711;
t1154 = t1060 * t878;
t1153 = t1060 * t992;
t1152 = t1060 * t993;
t932 = t1061 * t936;
t1151 = t1062 * t936;
t1141 = t1052 * t1067;
t1033 = t1174 * t1141;
t1130 = qJDD(1) * t1057;
t847 = -pkin(6) * t1033 + (pkin(4) * t1130 + t939 + (t1000 * t1060 + t1132) * qJD(1)) * t1058 + t1125;
t1150 = t1063 * t847;
t917 = t1036 - t1166;
t1149 = t1063 * t917;
t1147 = t1065 * t847;
t1146 = t1065 * t917;
t1143 = t1037 * t1063;
t1142 = t1037 * t1065;
t1129 = qJDD(1) * t1060;
t1124 = t1058 * t1166;
t1123 = t1061 * t1166;
t1122 = pkin(2) * t1127 - qJ(3) * t1008 - t932;
t1053 = t1060 ^ 2;
t1121 = t1053 * t1141;
t1119 = t1060 * t1135;
t1034 = t1058 * t1135;
t1117 = t1058 * t1127;
t1025 = -t1064 * qJDD(1) - t1066 * t1067;
t1113 = pkin(5) * t1025 + t1066 * g(1);
t827 = t1057 * t880 - t1085;
t778 = t1057 * t827 + t1060 * t828;
t712 = t1063 * t759 + t1065 * t760;
t696 = t1057 * t712 + t1155;
t705 = -pkin(4) * t847 + pkin(6) * t712;
t678 = -pkin(6) * t1155 - qJ(4) * t696 - t1057 * t705;
t686 = -pkin(3) * t696 - t1172;
t697 = t1060 * t712 - t1160;
t694 = t1058 * t847 + t1061 * t697;
t1108 = -pkin(2) * t696 + qJ(3) * t694 + t1058 * t678 + t1061 * t686;
t821 = t1063 * t877 + t1065 * t874;
t896 = -t979 - t980;
t703 = -pkin(4) * t896 + pkin(6) * t821 + t712;
t704 = -pkin(6) * t819 - t711;
t766 = t1057 * t821 + t1060 * t819;
t688 = -qJ(4) * t766 - t1057 * t703 + t1060 * t704;
t736 = -pkin(3) * t766 - t1171;
t768 = -t1057 * t819 + t1060 * t821;
t750 = t1058 * t896 + t1061 * t768;
t1107 = -pkin(2) * t766 + qJ(3) * t750 + t1058 * t688 + t1061 * t736;
t915 = -t1035 - t979;
t857 = t1065 * t915 - t1193;
t1075 = t1089 * t1043;
t872 = (qJD(5) - t1037) * t983 + t1075;
t772 = -pkin(4) * t872 + pkin(6) * t857 - t1147;
t856 = t1063 * t915 + t1192;
t786 = -pkin(6) * t856 + t1150;
t788 = t1057 * t857 + t1060 * t856;
t707 = -qJ(4) * t788 - t1057 * t772 + t1060 * t786;
t1086 = pkin(4) * t856 - t759;
t718 = -pkin(3) * t788 - t1086;
t789 = -t1057 * t856 + t1060 * t857;
t771 = t1058 * t872 + t1061 * t789;
t1106 = -pkin(2) * t788 + qJ(3) * t771 + t1058 * t707 + t1061 * t718;
t948 = -t980 - t1035;
t859 = -t1063 * t948 + t1146;
t776 = -pkin(4) * t1186 + pkin(6) * t859 + t1150;
t858 = t1065 * t948 + t1149;
t794 = -pkin(6) * t858 + t1147;
t798 = t1057 * t859 + t1060 * t858;
t716 = -qJ(4) * t798 - t1057 * t776 + t1060 * t794;
t1078 = pkin(4) * t858 - t760;
t726 = -pkin(3) * t798 - t1078;
t799 = -t1057 * t858 + t1060 * t859;
t781 = t1058 * t1186 + t1061 * t799;
t1105 = -pkin(2) * t798 + qJ(3) * t781 + t1058 * t716 + t1061 * t726;
t1046 = t1054 * t1067;
t1004 = -t1033 - t1046;
t941 = t1057 * t1004 + t1152;
t807 = -pkin(3) * t941 + t827;
t842 = -qJ(4) * t941 + t1159;
t945 = t1060 * t1004 - t1157;
t1083 = t1119 - t1130;
t989 = t1083 * t1058;
t902 = -t1058 * t989 + t1061 * t945;
t1104 = -pkin(2) * t941 + qJ(3) * t902 + t1058 * t842 + t1061 * t807;
t1007 = -t1046 - t1121;
t942 = t1060 * t1007 + t1158;
t808 = -pkin(3) * t942 + t828;
t843 = -qJ(4) * t942 + t1154;
t946 = -t1057 * t1007 + t1153;
t1020 = t1057 * t1034;
t990 = t1020 + t1115;
t903 = t1058 * t990 + t1061 * t946;
t1103 = -pkin(2) * t942 + qJ(3) * t903 + t1058 * t843 + t1061 * t808;
t1102 = -t1064 * t1029 - t1066 * t1030;
t1042 = t1052 * qJDD(1);
t1044 = t1054 * qJDD(1);
t1014 = t1044 + t1042;
t1021 = t1046 + t1141;
t1101 = pkin(2) * t1021 + qJ(3) * t1014 + t846;
t1051 = t1058 * t1052;
t1100 = t1051 * t1120;
t1099 = t1057 * t1119;
t1098 = t1060 * t1034;
t1097 = t1058 * t1114;
t1096 = -pkin(3) * t878 + qJ(4) * t778;
t777 = t1057 * t828 - t1060 * t827;
t844 = t1058 * t914 - t1061 * t913;
t975 = t1016 * t1137;
t976 = -t1059 * t1034 + t1097;
t1094 = t1064 * t976 + t1066 * t975;
t1093 = t1064 * t975 - t1066 * t976;
t1091 = t1052 * t1099;
t1090 = t1066 * t1029 - t1064 * t1030;
t1005 = (t1054 * t1058 + t1051) * t1067;
t1088 = -pkin(2) * t1043 + qJ(3) * t1005 + t931;
t988 = (t1119 + t1130) * t1058;
t991 = -t1020 + t1115;
t922 = -t1057 * t988 - t1060 * t991;
t762 = -qJ(4) * t922 - t777;
t924 = t1057 * t991 - t1060 * t988;
t994 = t1033 + t1121;
t887 = -t1058 * t994 + t1061 * t924;
t1084 = qJ(3) * t887 + t1058 * t762 + (-pkin(2) - t1168) * t922;
t1080 = -pkin(3) * t990 + qJ(4) * t946 + t1159;
t1079 = pkin(3) * t989 + qJ(4) * t945 - t1154;
t1077 = pkin(3) * t994 + qJ(4) * t924 + t778;
t757 = t1058 * t878 + t1061 * t778;
t1074 = qJ(3) * t757 + (-pkin(2) + t1095) * t777;
t1073 = -pkin(3) * t896 + qJ(4) * t768 + t1057 * t704 + t1060 * t703;
t1072 = -pkin(3) * t872 + qJ(4) * t789 + t1057 * t786 + t1060 * t772;
t1071 = -pkin(3) * t1186 + qJ(4) * t799 + t1057 * t794 + t1060 * t776;
t1070 = -pkin(3) * t847 - pkin(6) * t1160 + qJ(4) * t697 + t1060 * t705;
t1031 = 0.2e1 * t1117;
t1026 = t1066 * qJDD(1) - t1064 * t1067;
t1022 = -t1046 + t1141;
t1015 = t1044 - t1042;
t1006 = t1046 - t1121;
t1003 = t1033 - t1046;
t996 = pkin(5) * t1026 + t1064 * g(1);
t995 = -t1033 + t1121;
t978 = (t1053 + t1174) * t1034;
t977 = (qJDD(1) * t1053 + t1099) * t1058;
t974 = (-t1053 * t1135 + t1057 * t1129) * t1058;
t973 = (t1057 * t1135 + t1129) * t1140;
t972 = t1083 * t1140;
t969 = t1062 * t1005 + t1058 * t1128;
t966 = t1059 * t1005 - t1058 * t1126;
t961 = t1062 * t1015 + t1059 * t1022;
t960 = t1062 * t1014 - t1059 * t1021;
t959 = t1059 * t1015 - t1062 * t1022;
t958 = t1059 * t1014 + t1062 * t1021;
t957 = -t980 + t1035;
t956 = t979 - t1035;
t952 = t1061 * t977 + t1100;
t951 = -t1061 * t972 - t1100;
t950 = t1059 * t978 - t1097;
t949 = -t1058 * t1116 - t1062 * t978;
t947 = -t1057 * t1006 + t1152;
t944 = t1060 * t1003 + t1158;
t943 = t1060 * t1006 + t1157;
t940 = t1057 * t1003 - t1153;
t935 = t1058 * t977 - t1091;
t934 = -t1058 * t972 + t1091;
t929 = -pkin(1) * t1016 - t955;
t928 = pkin(1) * t1017 - t954;
t927 = t980 - t979;
t925 = -t983 * qJD(5) - t1076;
t923 = -t1057 * t990 + t1060 * t989;
t921 = t1057 * t989 + t1060 * t990;
t908 = t1064 * t966 - t1066 * t969;
t907 = t1064 * t969 + t1066 * t966;
t906 = (-t1063 * t983 + t1065 * t981) * t1037;
t905 = (t1063 * t981 + t1065 * t983) * t1037;
t904 = t1058 * t991 + t1061 * t947;
t901 = -t1058 * t988 + t1061 * t944;
t900 = t1058 * t947 - t1061 * t991;
t899 = t1058 * t946 - t1061 * t990;
t898 = t1058 * t945 + t1061 * t989;
t897 = t1058 * t944 + t1061 * t988;
t893 = t1059 * t974 + t1062 * t952;
t892 = -t1059 * t973 + t1062 * t951;
t891 = t1059 * t952 - t1062 * t974;
t890 = t1059 * t951 + t1062 * t973;
t886 = t1058 * t995 + t1061 * t923;
t883 = t1058 * t924 + t1061 * t994;
t882 = t1058 * t923 - t1061 * t995;
t881 = pkin(1) * t885;
t873 = t1075 + t1180;
t870 = pkin(1) * t1165 + qJ(2) * t1111;
t869 = t1065 * t926 + t983 * t1143;
t868 = t1063 * t926 - t983 * t1142;
t867 = -t1063 * t925 - t981 * t1142;
t866 = t1065 * t925 - t981 * t1143;
t865 = t1065 * t956 + t1149;
t864 = -t1063 * t957 + t1192;
t863 = t1063 * t956 - t1146;
t862 = t1065 * t957 + t1193;
t861 = pkin(1) * t968 + t1122;
t860 = pkin(1) * t966 + t1088;
t855 = t1059 * t943 + t1062 * t904;
t854 = t1059 * t942 + t1062 * t903;
t853 = t1059 * t941 + t1062 * t902;
t852 = t1059 * t940 + t1062 * t901;
t851 = t1059 * t904 - t1062 * t943;
t850 = t1059 * t903 - t1062 * t942;
t849 = t1059 * t902 - t1062 * t941;
t848 = t1059 * t901 - t1062 * t940;
t838 = t1059 * t922 + t1062 * t887;
t837 = t1059 * t921 + t1062 * t886;
t836 = t1059 * t887 - t1062 * t922;
t835 = t1059 * t886 - t1062 * t921;
t834 = -qJ(2) * t966 - t1059 * t914 + t1061 * t1151;
t833 = -qJ(2) * t968 + t1058 * t1151 - t1059 * t913;
t832 = qJ(2) * t969 + t1059 * t932 + t1062 * t914;
t831 = -qJ(2) * t970 + t1058 * t1156 + t1062 * t913;
t830 = -t1057 * t905 + t1060 * t906;
t829 = t1057 * t906 + t1060 * t905;
t825 = -t1058 * t1036 + t1061 * t830;
t824 = t1061 * t1036 + t1058 * t830;
t823 = -qJ(2) * t958 - t1062 * t844;
t822 = qJ(2) * t960 - t1059 * t844;
t820 = -t1063 * t1186 - t1065 * t872;
t818 = -t1063 * t872 + t1065 * t1186;
t817 = t1062 * t846 + t1156;
t816 = t1059 * t846 - t1151;
t815 = -t1057 * t868 + t1060 * t869;
t814 = -t1057 * t866 + t1060 * t867;
t813 = t1057 * t869 + t1060 * t868;
t812 = t1057 * t867 + t1060 * t866;
t811 = pkin(1) * t958 + t1101;
t803 = -t1057 * t863 + t1060 * t865;
t802 = -t1057 * t862 + t1060 * t864;
t801 = t1057 * t865 + t1060 * t863;
t800 = t1057 * t864 + t1060 * t862;
t796 = -pkin(2) * t898 - t1079;
t795 = -pkin(2) * t899 - t1080;
t793 = t1061 * t815 + t1124;
t792 = t1061 * t814 - t1124;
t791 = t1058 * t815 - t1123;
t790 = t1058 * t814 + t1123;
t785 = -t1058 * t873 + t1061 * t803;
t784 = t1058 * t877 + t1061 * t802;
t783 = t1058 * t803 + t1061 * t873;
t782 = t1058 * t802 - t1061 * t877;
t780 = t1058 * t799 - t1061 * t1186;
t774 = t1059 * t829 + t1062 * t825;
t773 = t1059 * t825 - t1062 * t829;
t770 = t1058 * t789 - t1061 * t872;
t767 = -t1057 * t818 + t1060 * t820;
t765 = t1057 * t820 + t1060 * t818;
t763 = pkin(1) * t816 + t1164;
t756 = t1058 * t778 - t1061 * t878;
t755 = t1058 * t927 + t1061 * t767;
t754 = t1058 * t767 - t1061 * t927;
t753 = -qJ(3) * t899 - t1058 * t808 + t1061 * t843;
t752 = -qJ(3) * t898 - t1058 * t807 + t1061 * t842;
t749 = t1058 * t768 - t1061 * t896;
t747 = t1059 * t813 + t1062 * t793;
t746 = t1059 * t812 + t1062 * t792;
t745 = t1059 * t793 - t1062 * t813;
t744 = t1059 * t792 - t1062 * t812;
t743 = -qJ(2) * t816 + (pkin(2) * t1059 - qJ(3) * t1062) * t844;
t742 = -pkin(2) * t883 - t1077;
t741 = t1059 * t801 + t1062 * t785;
t740 = t1059 * t800 + t1062 * t784;
t739 = t1059 * t785 - t1062 * t801;
t738 = t1059 * t784 - t1062 * t800;
t737 = -qJ(3) * t883 + t1061 * t762 + t922 * t1169;
t734 = t1059 * t798 + t1062 * t781;
t733 = t1059 * t781 - t1062 * t798;
t732 = pkin(1) * t850 + t1103;
t731 = pkin(1) * t849 + t1104;
t730 = qJ(2) * t817 + (-pkin(2) * t1062 - qJ(3) * t1059 - pkin(1)) * t844;
t729 = t1059 * t788 + t1062 * t771;
t728 = t1059 * t771 - t1062 * t788;
t727 = pkin(1) * t836 + t1084;
t724 = t1059 * t777 + t1062 * t757;
t723 = t1059 * t757 - t1062 * t777;
t722 = -qJ(2) * t850 - t1059 * t795 + t1062 * t753;
t721 = -qJ(2) * t849 - t1059 * t796 + t1062 * t752;
t720 = t1059 * t765 + t1062 * t755;
t719 = t1059 * t755 - t1062 * t765;
t714 = t1059 * t766 + t1062 * t750;
t713 = t1059 * t750 - t1062 * t766;
t710 = -pkin(2) * t756 - t1096;
t709 = -pkin(1) * t899 + qJ(2) * t854 + t1059 * t753 + t1062 * t795;
t708 = -pkin(1) * t898 + qJ(2) * t853 + t1059 * t752 + t1062 * t796;
t702 = -qJ(3) * t756 + (-qJ(4) * t1061 + t1169) * t777;
t701 = -qJ(2) * t836 - t1059 * t742 + t1062 * t737;
t700 = -pkin(2) * t780 - t1071;
t699 = -pkin(1) * t883 + qJ(2) * t838 + t1059 * t737 + t1062 * t742;
t698 = -pkin(2) * t770 - t1072;
t693 = t1058 * t697 - t1061 * t847;
t691 = -qJ(3) * t780 - t1058 * t726 + t1061 * t716;
t690 = -qJ(3) * t770 - t1058 * t718 + t1061 * t707;
t689 = pkin(1) * t723 + t1074;
t684 = pkin(1) * t733 + t1105;
t683 = -qJ(2) * t723 - t1059 * t710 + t1062 * t702;
t682 = pkin(1) * t728 + t1106;
t681 = -pkin(2) * t749 - t1073;
t680 = -pkin(1) * t756 + qJ(2) * t724 + t1059 * t702 + t1062 * t710;
t679 = -qJ(3) * t749 - t1058 * t736 + t1061 * t688;
t676 = t1059 * t696 + t1062 * t694;
t675 = t1059 * t694 - t1062 * t696;
t674 = -qJ(2) * t733 - t1059 * t700 + t1062 * t691;
t673 = -qJ(2) * t728 - t1059 * t698 + t1062 * t690;
t672 = -pkin(1) * t780 + qJ(2) * t734 + t1059 * t691 + t1062 * t700;
t671 = -pkin(1) * t770 + qJ(2) * t729 + t1059 * t690 + t1062 * t698;
t670 = pkin(1) * t713 + t1107;
t669 = -pkin(2) * t693 - t1070;
t668 = -qJ(2) * t713 - t1059 * t681 + t1062 * t679;
t667 = -pkin(1) * t749 + qJ(2) * t714 + t1059 * t679 + t1062 * t681;
t666 = -qJ(3) * t693 - t1058 * t686 + t1061 * t678;
t665 = pkin(1) * t675 + t1108;
t664 = -qJ(2) * t675 - t1059 * t669 + t1062 * t666;
t663 = -pkin(1) * t693 + qJ(2) * t676 + t1059 * t666 + t1062 * t669;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, qJDD(1), -t1030, t1029, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t928, t929, 0, t881, t1042, t1031, 0, t1044, 0, 0, t861, t860, t811, t763, t935, t882, t900, t934, t897, t1044, t731, t732, t727, t689, t791, t754, t782, t790, t783, t824, t682, t684, t670, t665; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, -t1025, 0, t1026, 0, t1113, -t996, -t1090, -pkin(5) * t1090, 0, 0, -t1176, 0, -t1177, 0, t1188, -t1189, -t1190, -pkin(5) * t1190 - qJ(2) * t1148 + t1066 * t870, t1094, t1064 * t961 + t1066 * t959, t907, -t1094, t1182, 0, -pkin(5) * t1181 + t1064 * t833 + t1066 * t831, -pkin(5) * t908 + t1064 * t834 + t1066 * t832, t1064 * t823 + t1066 * t822 - pkin(5) * (t1064 * t958 - t1066 * t960), t1064 * t743 + t1066 * t730 - pkin(5) * (t1064 * t816 - t1066 * t817), t1064 * t893 + t1066 * t891, t1064 * t837 + t1066 * t835, t1064 * t855 + t1066 * t851, t1064 * t892 + t1066 * t890, t1064 * t852 + t1066 * t848, t1064 * t950 + t1066 * t949, t1064 * t721 + t1066 * t708 - pkin(5) * (t1064 * t849 - t1066 * t853), t1064 * t722 + t1066 * t709 - pkin(5) * (t1064 * t850 - t1066 * t854), t1064 * t701 + t1066 * t699 - pkin(5) * (t1064 * t836 - t1066 * t838), t1064 * t683 + t1066 * t680 - pkin(5) * (t1064 * t723 - t1066 * t724), t1064 * t747 + t1066 * t745, t1064 * t720 + t1066 * t719, t1064 * t740 + t1066 * t738, t1064 * t746 + t1066 * t744, t1064 * t741 + t1066 * t739, t1064 * t774 + t1066 * t773, t1064 * t673 + t1066 * t671 - pkin(5) * (t1064 * t728 - t1066 * t729), t1064 * t674 + t1066 * t672 - pkin(5) * (t1064 * t733 - t1066 * t734), t1064 * t668 + t1066 * t667 - pkin(5) * (t1064 * t713 - t1066 * t714), t1064 * t664 + t1066 * t663 - pkin(5) * (t1064 * t675 - t1066 * t676); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, -t1026, 0, -t1025, 0, t996, t1113, t1102, pkin(5) * t1102, 0, 0, t1177, 0, -t1176, 0, t1189, t1188, t1191, pkin(5) * t1191 + qJ(2) * t1145 + t1064 * t870, t1093, t1064 * t959 - t1066 * t961, t908, -t1093, -t1181, 0, -pkin(5) * t1182 + t1064 * t831 - t1066 * t833, pkin(5) * t907 + t1064 * t832 - t1066 * t834, -t1066 * t823 + t1064 * t822 + pkin(5) * (t1064 * t960 + t1066 * t958), -t1066 * t743 + t1064 * t730 + pkin(5) * (t1064 * t817 + t1066 * t816), t1064 * t891 - t1066 * t893, t1064 * t835 - t1066 * t837, t1064 * t851 - t1066 * t855, t1064 * t890 - t1066 * t892, t1064 * t848 - t1066 * t852, t1064 * t949 - t1066 * t950, -t1066 * t721 + t1064 * t708 + pkin(5) * (t1064 * t853 + t1066 * t849), -t1066 * t722 + t1064 * t709 + pkin(5) * (t1064 * t854 + t1066 * t850), -t1066 * t701 + t1064 * t699 + pkin(5) * (t1064 * t838 + t1066 * t836), -t1066 * t683 + t1064 * t680 + pkin(5) * (t1064 * t724 + t1066 * t723), t1064 * t745 - t1066 * t747, t1064 * t719 - t1066 * t720, t1064 * t738 - t1066 * t740, t1064 * t744 - t1066 * t746, t1064 * t739 - t1066 * t741, t1064 * t773 - t1066 * t774, -t1066 * t673 + t1064 * t671 + pkin(5) * (t1064 * t729 + t1066 * t728), -t1066 * t674 + t1064 * t672 + pkin(5) * (t1064 * t734 + t1066 * t733), -t1066 * t668 + t1064 * t667 + pkin(5) * (t1064 * t714 + t1066 * t713), -t1066 * t664 + t1064 * t663 + pkin(5) * (t1064 * t676 + t1066 * t675); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1067, 0, 0, -g(1), t1030, 0, 0, 0, t1017, 0, -t1016, 0, t1175, -t986, -t885, -qJ(2) * t885, t976, t961, t969, -t976, t970, 0, t833, t834, t823, t743, t893, t837, t855, t892, t852, t950, t721, t722, t701, t683, t747, t720, t740, t746, t741, t774, t673, t674, t668, t664; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1067, 0, qJDD(1), 0, g(1), 0, -t1029, 0, 0, 0, t1016, 0, t1017, 0, t986, t1175, t1111, t870, t975, t959, t966, -t975, -t968, 0, t831, t832, t822, t730, t891, t835, t851, t890, t848, t949, t708, t709, t699, t680, t745, t719, t738, t744, t739, t773, t671, t672, t667, t663; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1030, t1029, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t928, t929, 0, t881, t1042, t1031, 0, t1044, 0, 0, t861, t860, t811, t763, t935, t882, t900, t934, t897, t1044, t731, t732, t727, t689, t791, t754, t782, t790, t783, t824, t682, t684, t670, t665; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1067, 0, 0, -t1165, t954, 0, t1117, t1015, t1005, -t1117, t1008, 0, t931, t932, -t844, -qJ(3) * t844, t952, t886, t904, t951, t901, -t1117, t752, t753, t737, t702, t793, t755, t784, t792, t785, t825, t690, t691, t679, t666; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1067, 0, qJDD(1), 0, t1165, 0, t955, 0, t1034, -t1022, -t1043, -t1034, -t1127, 0, t913, t914, 0, -pkin(2) * t844, -t974, -t921, -t943, t973, -t940, -t978, t796, t795, t742, t710, -t813, -t765, -t800, -t812, -t801, -t829, t698, t700, t681, t669; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t954, -t955, 0, 0, t1042, t1031, 0, t1044, 0, 0, t1122, t1088, t1101, t1164, t935, t882, t900, t934, t897, t1044, t1104, t1103, t1084, t1074, t791, t754, t782, t790, t783, t824, t1106, t1105, t1107, t1108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1043, t1127, t1034, 0, t1046, 0, 0, t936, t913, 0, t977, t923, t947, -t972, t944, 0, t842, t843, t762, -qJ(4) * t777, t815, t767, t802, t814, t803, t830, t707, t716, t688, t678; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1043, -t1141, t1127, -t1034, 0, -t936, 0, t914, 0, -t1018, -t995, -t991, t1018, t988, t1127, t807, t808, -pkin(3) * t922, -pkin(3) * t777, -t1166, -t927, -t877, t1166, t873, t1036, t718, t726, t736, t686; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1034, t1022, t1043, t1034, t1127, 0, -t913, -t914, 0, 0, t974, t921, t943, -t973, t940, t978, t1079, t1080, t1077, t1096, t813, t765, t800, t812, t801, t829, t1072, t1071, t1073, t1070; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1115, t989, t993, -t1020, t1003, t1020, 0, t878, t827, 0, t869, t820, t864, t867, t865, t906, t786, t794, t704, -pkin(6) * t711; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1098, t990, t1006, -t1118, -t992, t1098, -t878, 0, t828, 0, t868, t818, t862, t866, t863, t905, t772, t776, t703, t705; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1018, t995, t991, -t1018, -t988, -t1127, -t827, -t828, 0, 0, t1166, t927, t877, -t1166, -t873, -t1036, t1086, t1078, t1171, t1172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t926, -t872, t1184, -t1144, t956, t1144, 0, t847, t759, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1161, t1186, t957, t925, -t917, t1161, -t847, 0, t760, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1166, t927, t877, -t1166, -t873, -t1036, -t759, -t760, 0, 0;];
m_new_reg = t1;
