% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RRRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 20:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RRRRPR6_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR6_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR6_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR6_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_invdynB_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:52:41
% EndTime: 2019-05-07 20:54:09
% DurationCPUTime: 64.34s
% Computational Cost: add. (625694->945), mult. (1312322->1460), div. (0->0), fcn. (987166->12), ass. (0->652)
t1127 = sin(qJ(3));
t1128 = sin(qJ(2));
t1192 = qJD(1) * t1128;
t1114 = qJD(2) * t1192;
t1133 = cos(qJ(2));
t1163 = t1133 * qJDD(1);
t1087 = -t1114 + t1163;
t1079 = -qJDD(3) + t1087;
t1132 = cos(qJ(3));
t1082 = qJD(2) * t1127 + t1132 * t1192;
t1153 = qJD(2) * t1132 - t1127 * t1192;
t1144 = t1153 * t1082;
t1231 = -t1079 + t1144;
t1233 = t1127 * t1231;
t1232 = t1132 * t1231;
t1151 = t1153 ^ 2;
t1125 = sin(qJ(6));
t1130 = cos(qJ(6));
t1122 = sin(pkin(11));
t1123 = cos(pkin(11));
t1126 = sin(qJ(4));
t1131 = cos(qJ(4));
t1043 = t1131 * t1082 + t1126 * t1153;
t1191 = qJD(1) * t1133;
t1156 = qJD(2) * t1191;
t1166 = qJDD(1) * t1128;
t1086 = t1156 + t1166;
t1038 = t1153 * qJD(3) + t1127 * qJDD(2) + t1132 * t1086;
t1143 = t1132 * qJDD(2) - t1127 * t1086;
t1139 = qJD(3) * t1082 - t1143;
t1152 = t1126 * t1038 + t1131 * t1139;
t1140 = qJD(4) * t1043 + t1152;
t1042 = t1082 * t1126 - t1131 * t1153;
t958 = -t1042 * qJD(4) + t1131 * t1038 - t1126 * t1139;
t1155 = t1122 * t958 + t1123 * t1140;
t887 = -t1122 * t1140 + t1123 * t958;
t987 = t1123 * t1042 + t1043 * t1122;
t989 = -t1042 * t1122 + t1043 * t1123;
t932 = t1125 * t989 + t1130 * t987;
t796 = -qJD(6) * t932 - t1125 * t1155 + t1130 * t887;
t1111 = -qJD(3) + t1191;
t1102 = -qJD(4) + t1111;
t1095 = -qJD(6) + t1102;
t912 = t932 * t1095;
t1230 = t796 + t912;
t970 = t987 * t1102;
t864 = t970 - t887;
t1229 = t970 + t887;
t1029 = t1042 * t1102;
t931 = t1029 - t958;
t1228 = t1029 + t958;
t1075 = -qJDD(4) + t1079;
t1213 = t987 * t989;
t1138 = -t1075 - t1213;
t1227 = t1122 * t1138;
t1226 = t1123 * t1138;
t1070 = -qJDD(6) + t1075;
t934 = -t1125 * t987 + t1130 * t989;
t1214 = t932 * t934;
t1136 = -t1070 - t1214;
t1225 = t1125 * t1136;
t1186 = t1042 * t1043;
t1137 = -t1075 - t1186;
t1224 = t1126 * t1137;
t1223 = t1130 * t1136;
t1222 = t1131 * t1137;
t1065 = t1153 * t1111;
t1008 = -t1038 - t1065;
t1154 = t1125 * t887 + t1130 * t1155;
t757 = (qJD(6) + t1095) * t934 + t1154;
t927 = (qJD(4) + t1102) * t1043 + t1152;
t1003 = (qJD(3) + t1111) * t1082 - t1143;
t924 = t932 ^ 2;
t925 = t934 ^ 2;
t1221 = t987 ^ 2;
t986 = t989 ^ 2;
t1220 = qJD(2) ^ 2;
t1219 = t1042 ^ 2;
t1041 = t1043 ^ 2;
t1078 = t1082 ^ 2;
t1092 = t1095 ^ 2;
t1218 = t1102 ^ 2;
t1108 = t1111 ^ 2;
t1217 = pkin(2) * t1128;
t1216 = pkin(2) * t1133;
t1215 = t1133 * g(3);
t1129 = sin(qJ(1));
t1134 = cos(qJ(1));
t1099 = t1129 * g(1) - t1134 * g(2);
t1135 = qJD(1) ^ 2;
t1072 = qJDD(1) * pkin(1) + t1135 * pkin(7) + t1099;
t1145 = -t1087 + t1114;
t1146 = t1086 + t1156;
t1002 = pkin(2) * t1145 - pkin(8) * t1146 - t1072;
t1100 = g(1) * t1134 + g(2) * t1129;
t1073 = -pkin(1) * t1135 + qJDD(1) * pkin(7) - t1100;
t1059 = -g(3) * t1128 + t1133 * t1073;
t1148 = -pkin(8) * t1128 - t1216;
t1084 = t1148 * qJD(1);
t1019 = -pkin(2) * t1220 + qJDD(2) * pkin(8) + t1084 * t1191 + t1059;
t959 = -t1132 * t1002 + t1127 * t1019;
t897 = pkin(3) * t1231 + pkin(9) * t1008 - t959;
t1060 = -pkin(3) * t1111 - pkin(9) * t1082;
t960 = t1127 * t1002 + t1132 * t1019;
t901 = -pkin(3) * t1151 - pkin(9) * t1139 + t1111 * t1060 + t960;
t827 = t1126 * t901 - t1131 * t897;
t784 = t1137 * pkin(4) + qJ(5) * t931 - t827;
t1022 = -pkin(4) * t1102 - qJ(5) * t1043;
t828 = t1126 * t897 + t1131 * t901;
t798 = -pkin(4) * t1219 - qJ(5) * t1140 + t1102 * t1022 + t828;
t711 = 0.2e1 * qJD(5) * t989 + t1122 * t798 - t1123 * t784;
t679 = t1138 * pkin(5) + pkin(10) * t864 - t711;
t712 = -0.2e1 * qJD(5) * t987 + t1122 * t784 + t1123 * t798;
t967 = -pkin(5) * t1102 - pkin(10) * t989;
t684 = -pkin(5) * t1221 - pkin(10) * t1155 + t1102 * t967 + t712;
t627 = t1125 * t679 + t1130 * t684;
t1212 = t1102 * t989;
t626 = t1125 * t684 - t1130 * t679;
t584 = t1125 * t627 - t1130 * t626;
t1211 = t1122 * t584;
t1018 = t1215 - qJDD(2) * pkin(2) - t1220 * pkin(8) + (qJD(1) * t1084 + t1073) * t1128;
t943 = t1139 * pkin(3) - t1151 * pkin(9) + t1082 * t1060 + t1018;
t831 = t1140 * pkin(4) - t1219 * qJ(5) + t1043 * t1022 + qJDD(5) + t943;
t1210 = t1122 * t831;
t916 = t1075 - t1213;
t1209 = t1122 * t916;
t1208 = t1123 * t584;
t1207 = t1123 * t831;
t1206 = t1123 * t916;
t749 = pkin(5) * t1155 - pkin(10) * t1221 + t989 * t967 + t831;
t1205 = t1125 * t749;
t846 = t1070 - t1214;
t1204 = t1125 * t846;
t653 = t1122 * t712 - t1123 * t711;
t1203 = t1126 * t653;
t1202 = t1126 * t943;
t972 = t1075 - t1186;
t1201 = t1126 * t972;
t750 = t1126 * t828 - t1131 * t827;
t1200 = t1127 * t750;
t1199 = t1130 * t749;
t1198 = t1130 * t846;
t1197 = t1131 * t653;
t1196 = t1131 * t943;
t1195 = t1131 * t972;
t1194 = t1132 * t750;
t1193 = qJD(1) * qJD(2);
t1190 = t1018 * t1127;
t1189 = t1018 * t1132;
t1032 = t1079 + t1144;
t1188 = t1032 * t1127;
t1187 = t1032 * t1132;
t1185 = t1072 * t1128;
t1184 = t1072 * t1133;
t1183 = t1075 * t1128;
t1182 = t1095 * t1125;
t1181 = t1095 * t1130;
t1110 = t1133 * t1135 * t1128;
t1097 = -t1110 + qJDD(2);
t1180 = t1097 * t1128;
t1179 = t1097 * t1133;
t1098 = qJDD(2) + t1110;
t1178 = t1098 * t1128;
t1177 = t1102 * t1122;
t1176 = t1102 * t1123;
t1175 = t1102 * t1126;
t1174 = t1102 * t1131;
t1119 = t1128 ^ 2;
t1173 = t1119 * t1135;
t1172 = t1127 * t1082;
t1171 = t1132 * t1082;
t1120 = t1133 ^ 2;
t1167 = t1119 + t1120;
t1165 = qJDD(1) * t1129;
t1164 = qJDD(1) * t1134;
t1162 = t1128 * t1214;
t1161 = t1128 * t1213;
t1160 = t1133 * t1214;
t1159 = t1133 * t1213;
t1158 = t1128 * t1186;
t1157 = t1133 * t1186;
t654 = t1122 * t711 + t1123 * t712;
t585 = t1125 * t626 + t1130 * t627;
t751 = t1126 * t827 + t1131 * t828;
t1058 = t1073 * t1128 + t1215;
t1011 = t1058 * t1128 + t1133 * t1059;
t1051 = -t1099 * t1129 - t1134 * t1100;
t1150 = t1129 * t1110;
t1149 = t1134 * t1110;
t1006 = -t1065 + t1038;
t1091 = -t1129 * t1135 + t1164;
t1147 = -pkin(6) * t1091 - g(3) * t1129;
t890 = t1127 * t960 - t1132 * t959;
t891 = t1127 * t959 + t1132 * t960;
t1010 = t1058 * t1133 - t1059 * t1128;
t1050 = t1099 * t1134 - t1100 * t1129;
t860 = t1155 + t1212;
t1142 = t1128 * t1144;
t1141 = t1133 * t1144;
t1117 = t1120 * t1135;
t1107 = -t1117 - t1220;
t1106 = t1117 - t1220;
t1105 = -t1173 - t1220;
t1104 = -t1173 + t1220;
t1094 = t1117 - t1173;
t1093 = t1117 + t1173;
t1090 = t1134 * t1135 + t1165;
t1089 = t1167 * qJDD(1);
t1088 = -0.2e1 * t1114 + t1163;
t1085 = 0.2e1 * t1156 + t1166;
t1081 = t1133 * t1098;
t1080 = t1167 * t1193;
t1069 = -pkin(6) * t1090 + g(3) * t1134;
t1067 = t1133 * t1075;
t1064 = -t1078 + t1108;
t1063 = t1151 - t1108;
t1062 = t1086 * t1133 - t1119 * t1193;
t1061 = -t1087 * t1128 - t1120 * t1193;
t1057 = -t1105 * t1128 - t1179;
t1056 = -t1104 * t1128 + t1081;
t1055 = t1107 * t1133 - t1178;
t1054 = t1106 * t1133 - t1180;
t1053 = t1105 * t1133 - t1180;
t1052 = t1107 * t1128 + t1081;
t1048 = -t1078 + t1151;
t1047 = t1089 * t1134 - t1093 * t1129;
t1046 = t1089 * t1129 + t1093 * t1134;
t1045 = -t1078 - t1108;
t1044 = -t1085 * t1128 + t1088 * t1133;
t1040 = -t1108 - t1151;
t1031 = t1151 + t1078;
t1028 = -t1041 + t1218;
t1027 = -t1218 + t1219;
t1026 = t1057 * t1134 + t1085 * t1129;
t1025 = t1055 * t1134 - t1088 * t1129;
t1024 = t1057 * t1129 - t1085 * t1134;
t1023 = t1055 * t1129 + t1088 * t1134;
t1021 = -pkin(7) * t1053 - t1184;
t1020 = -pkin(7) * t1052 - t1185;
t1017 = (-t1132 * t1153 - t1172) * t1111;
t1016 = (t1127 * t1153 - t1171) * t1111;
t1014 = -pkin(1) * t1053 + t1059;
t1013 = -pkin(1) * t1052 + t1058;
t1012 = -t1041 - t1218;
t1004 = (-qJD(3) + t1111) * t1082 + t1143;
t999 = t1038 * t1132 + t1111 * t1172;
t998 = -t1038 * t1127 + t1111 * t1171;
t997 = t1065 * t1132 + t1127 * t1139;
t996 = -t1065 * t1127 + t1132 * t1139;
t995 = t1017 * t1133 - t1079 * t1128;
t994 = -t1041 + t1219;
t993 = t1063 * t1132 + t1188;
t992 = -t1064 * t1127 + t1232;
t991 = -t1063 * t1127 + t1187;
t990 = -t1064 * t1132 - t1233;
t984 = -t1045 * t1127 + t1187;
t983 = t1045 * t1132 + t1188;
t982 = -t1218 - t1219;
t979 = t1011 * t1134 - t1072 * t1129;
t978 = t1011 * t1129 + t1072 * t1134;
t976 = t1040 * t1132 - t1233;
t975 = t1040 * t1127 + t1232;
t969 = -t986 + t1218;
t968 = -t1218 + t1221;
t966 = (t1042 * t1131 - t1043 * t1126) * t1102;
t965 = (t1042 * t1126 + t1043 * t1131) * t1102;
t964 = -t986 - t1218;
t963 = t1133 * t999 - t1142;
t962 = t1133 * t997 + t1142;
t961 = -t1041 - t1219;
t955 = -t1003 * t1132 - t1008 * t1127;
t954 = t1004 * t1132 - t1006 * t1127;
t953 = -t1003 * t1127 + t1008 * t1132;
t952 = -t1004 * t1127 - t1006 * t1132;
t951 = -pkin(8) * t983 + t1189;
t950 = -t1003 * t1128 + t1133 * t993;
t949 = -t1008 * t1128 + t1133 * t992;
t948 = t1027 * t1131 + t1201;
t947 = -t1028 * t1126 + t1222;
t946 = t1027 * t1126 - t1195;
t945 = t1028 * t1131 + t1224;
t944 = -pkin(8) * t975 + t1190;
t942 = t1006 * t1128 + t1133 * t984;
t941 = -t1006 * t1133 + t1128 * t984;
t940 = -t1012 * t1126 + t1195;
t939 = t1012 * t1131 + t1201;
t938 = -t1004 * t1128 + t1133 * t976;
t937 = t1004 * t1133 + t1128 * t976;
t936 = -t986 + t1221;
t935 = -t1048 * t1128 + t1133 * t954;
t926 = (qJD(4) - t1102) * t1043 + t1152;
t923 = -t1218 - t1221;
t922 = t1043 * t1175 + t1131 * t958;
t921 = -t1043 * t1174 + t1126 * t958;
t920 = -t1042 * t1174 + t1126 * t1140;
t919 = -t1042 * t1175 - t1131 * t1140;
t915 = t1131 * t982 - t1224;
t914 = t1126 * t982 + t1222;
t911 = -t925 + t1092;
t910 = t924 - t1092;
t909 = -t1031 * t1128 + t1133 * t955;
t908 = t1031 * t1133 + t1128 * t955;
t907 = (-t1122 * t989 + t1123 * t987) * t1102;
t906 = (t1122 * t987 + t1123 * t989) * t1102;
t905 = -t925 - t1092;
t904 = -pkin(2) * t983 + t960;
t903 = -t1127 * t965 + t1132 * t966;
t902 = -t1127 * t966 - t1132 * t965;
t900 = -pkin(2) * t975 + t959;
t898 = t1133 * t903 - t1183;
t894 = t1129 * t983 + t1134 * t942;
t893 = t1129 * t942 - t1134 * t983;
t892 = -t986 - t1221;
t889 = t1129 * t975 + t1134 * t938;
t888 = t1129 * t938 - t1134 * t975;
t883 = t1123 * t968 + t1209;
t882 = -t1122 * t969 + t1226;
t881 = t1122 * t968 - t1206;
t880 = t1123 * t969 + t1227;
t879 = -t1122 * t964 + t1206;
t878 = t1123 * t964 + t1209;
t877 = -t1127 * t946 + t1132 * t948;
t876 = -t1127 * t945 + t1132 * t947;
t875 = -t1127 * t948 - t1132 * t946;
t874 = -t1127 * t947 - t1132 * t945;
t873 = t1018 * t1128 + t1133 * t891;
t872 = -t1018 * t1133 + t1128 * t891;
t871 = -pkin(9) * t939 + t1196;
t870 = -t1127 * t939 + t1132 * t940;
t869 = t1127 * t940 + t1132 * t939;
t868 = t1129 * t953 + t1134 * t909;
t867 = t1129 * t909 - t1134 * t953;
t866 = -pkin(9) * t914 + t1202;
t865 = -t925 + t924;
t859 = t1155 - t1212;
t858 = -t1126 * t931 - t1131 * t927;
t857 = -t1126 * t1228 - t1131 * t926;
t856 = -t1126 * t927 + t1131 * t931;
t855 = -t1126 * t926 + t1131 * t1228;
t854 = -t1092 - t924;
t853 = -pkin(1) * t941 + pkin(2) * t1006 - pkin(8) * t984 - t1190;
t852 = t1123 * t887 + t1177 * t989;
t851 = t1122 * t887 - t1176 * t989;
t850 = t1122 * t1155 - t1176 * t987;
t849 = -t1123 * t1155 - t1177 * t987;
t845 = t1123 * t923 - t1227;
t844 = t1122 * t923 + t1226;
t843 = -t1127 * t921 + t1132 * t922;
t842 = -t1127 * t919 + t1132 * t920;
t841 = -t1127 * t922 - t1132 * t921;
t840 = -t1127 * t920 - t1132 * t919;
t839 = (-t1125 * t934 + t1130 * t932) * t1095;
t838 = (t1125 * t932 + t1130 * t934) * t1095;
t837 = -pkin(1) * t937 - pkin(2) * t1004 - pkin(8) * t976 + t1189;
t836 = -t1127 * t914 + t1132 * t915;
t835 = t1127 * t915 + t1132 * t914;
t834 = -pkin(8) * t953 - t890;
t833 = -t1126 * t906 + t1131 * t907;
t832 = t1126 * t907 + t1131 * t906;
t830 = t1133 * t843 + t1158;
t829 = t1133 * t842 - t1158;
t825 = -t1128 * t927 + t1133 * t877;
t824 = -t1128 * t931 + t1133 * t876;
t823 = -pkin(3) * t1228 + pkin(9) * t940 + t1202;
t822 = t1128 * t1228 + t1133 * t870;
t821 = t1128 * t870 - t1133 * t1228;
t820 = -pkin(7) * t941 - t1128 * t904 + t1133 * t951;
t819 = -pkin(3) * t926 + pkin(9) * t915 - t1196;
t818 = -t924 - t925;
t817 = -pkin(7) * t937 - t1128 * t900 + t1133 * t944;
t816 = t1128 * t926 + t1133 * t836;
t815 = t1128 * t836 - t1133 * t926;
t814 = t1130 * t910 + t1204;
t813 = -t1125 * t911 + t1223;
t812 = t1125 * t910 - t1198;
t811 = t1130 * t911 + t1225;
t810 = -t1125 * t905 + t1198;
t809 = t1130 * t905 + t1204;
t808 = t1129 * t890 + t1134 * t873;
t807 = t1129 * t873 - t1134 * t890;
t806 = -t1126 * t881 + t1131 * t883;
t805 = -t1126 * t880 + t1131 * t882;
t804 = t1126 * t883 + t1131 * t881;
t803 = t1126 * t882 + t1131 * t880;
t802 = -t1126 * t878 + t1131 * t879;
t801 = t1126 * t879 + t1131 * t878;
t800 = -pkin(1) * t908 - pkin(2) * t1031 - pkin(8) * t955 - t891;
t799 = -pkin(1) * t872 + pkin(2) * t1018 - pkin(8) * t891;
t795 = -qJD(6) * t934 - t1154;
t794 = -qJ(5) * t878 + t1207;
t793 = -pkin(7) * t908 + t1133 * t834 + t1217 * t953;
t792 = -t1127 * t856 + t1132 * t858;
t791 = -t1127 * t855 + t1132 * t857;
t790 = t1127 * t858 + t1132 * t856;
t789 = -t1127 * t857 - t1132 * t855;
t788 = -t1122 * t864 - t1123 * t860;
t787 = -t1122 * t1229 - t1123 * t859;
t786 = -t1122 * t860 + t1123 * t864;
t785 = -t1122 * t859 + t1123 * t1229;
t781 = t1130 * t854 - t1225;
t780 = t1125 * t854 + t1223;
t779 = -t1126 * t851 + t1131 * t852;
t778 = -t1126 * t849 + t1131 * t850;
t777 = t1126 * t852 + t1131 * t851;
t776 = t1126 * t850 + t1131 * t849;
t775 = -t1126 * t844 + t1131 * t845;
t774 = t1126 * t845 + t1131 * t844;
t773 = -t1122 * t838 + t1123 * t839;
t772 = t1122 * t839 + t1123 * t838;
t771 = -qJ(5) * t844 + t1210;
t770 = -t1127 * t832 + t1132 * t833;
t769 = -t1127 * t833 - t1132 * t832;
t768 = -t1128 * t994 + t1133 * t791;
t767 = t1133 * t770 - t1183;
t766 = -pkin(7) * t872 + (-pkin(8) * t1133 + t1217) * t890;
t765 = t1129 * t869 + t1134 * t822;
t764 = t1129 * t822 - t1134 * t869;
t763 = t1128 * t961 + t1133 * t792;
t762 = t1128 * t792 - t1133 * t961;
t761 = -t796 + t912;
t756 = (qJD(6) - t1095) * t934 + t1154;
t755 = t1130 * t796 + t1182 * t934;
t754 = t1125 * t796 - t1181 * t934;
t753 = -t1125 * t795 - t1181 * t932;
t752 = t1130 * t795 - t1182 * t932;
t748 = -pkin(2) * t869 - pkin(3) * t939 + t828;
t747 = t1129 * t835 + t1134 * t816;
t746 = t1129 * t816 - t1134 * t835;
t745 = -pkin(2) * t835 - pkin(3) * t914 + t827;
t744 = -pkin(4) * t1229 + qJ(5) * t879 + t1210;
t743 = -pkin(3) * t943 + pkin(9) * t751;
t742 = -pkin(2) * t790 - pkin(3) * t856;
t741 = -t1122 * t812 + t1123 * t814;
t740 = -t1122 * t811 + t1123 * t813;
t739 = t1122 * t814 + t1123 * t812;
t738 = t1122 * t813 + t1123 * t811;
t737 = -pkin(8) * t869 - t1127 * t823 + t1132 * t871;
t736 = -pkin(4) * t859 + qJ(5) * t845 - t1207;
t735 = -t1122 * t809 + t1123 * t810;
t734 = t1122 * t810 + t1123 * t809;
t733 = -t1127 * t804 + t1132 * t806;
t732 = -t1127 * t803 + t1132 * t805;
t731 = -t1127 * t806 - t1132 * t804;
t730 = -t1127 * t805 - t1132 * t803;
t729 = -t1127 * t801 + t1132 * t802;
t728 = t1127 * t802 + t1132 * t801;
t727 = -pkin(8) * t835 - t1127 * t819 + t1132 * t866;
t726 = -pkin(9) * t856 - t750;
t725 = -pkin(3) * t961 + pkin(9) * t858 + t751;
t724 = -t1126 * t786 + t1131 * t788;
t723 = -t1126 * t785 + t1131 * t787;
t722 = t1126 * t788 + t1131 * t786;
t721 = t1126 * t787 + t1131 * t785;
t720 = -t1122 * t780 + t1123 * t781;
t719 = t1122 * t781 + t1123 * t780;
t718 = -t1127 * t777 + t1132 * t779;
t717 = -t1127 * t776 + t1132 * t778;
t716 = -t1127 * t779 - t1132 * t777;
t715 = -t1127 * t778 - t1132 * t776;
t714 = -t1127 * t774 + t1132 * t775;
t713 = t1127 * t775 + t1132 * t774;
t709 = -pkin(10) * t809 + t1199;
t708 = -t1128 * t860 + t1133 * t733;
t707 = -t1128 * t864 + t1133 * t732;
t706 = -t1126 * t772 + t1131 * t773;
t705 = t1126 * t773 + t1131 * t772;
t704 = t1129 * t790 + t1134 * t763;
t703 = t1129 * t763 - t1134 * t790;
t702 = t1128 * t1229 + t1133 * t729;
t701 = t1128 * t729 - t1133 * t1229;
t700 = t1133 * t718 + t1161;
t699 = t1133 * t717 - t1161;
t698 = -pkin(10) * t780 + t1205;
t697 = -t1125 * t761 - t1130 * t757;
t696 = -t1125 * t1230 - t1130 * t756;
t695 = -t1125 * t757 + t1130 * t761;
t694 = -t1125 * t756 + t1130 * t1230;
t693 = -t1122 * t754 + t1123 * t755;
t692 = -t1122 * t752 + t1123 * t753;
t691 = t1122 * t755 + t1123 * t754;
t690 = t1122 * t753 + t1123 * t752;
t689 = -pkin(1) * t821 + pkin(2) * t1228 - pkin(8) * t870 - t1127 * t871 - t1132 * t823;
t688 = t1132 * t751 - t1200;
t687 = t1127 * t751 + t1194;
t686 = t1128 * t859 + t1133 * t714;
t685 = t1128 * t714 - t1133 * t859;
t682 = -pkin(1) * t815 + pkin(2) * t926 - pkin(8) * t836 - t1127 * t866 - t1132 * t819;
t681 = t1128 * t943 + t1133 * t688;
t680 = t1128 * t688 - t1133 * t943;
t676 = -t1126 * t739 + t1131 * t741;
t675 = -t1126 * t738 + t1131 * t740;
t674 = t1126 * t741 + t1131 * t739;
t673 = t1126 * t740 + t1131 * t738;
t672 = -t1126 * t734 + t1131 * t735;
t671 = t1126 * t735 + t1131 * t734;
t670 = -pkin(9) * t801 - t1126 * t744 + t1131 * t794;
t669 = -pkin(5) * t1230 + pkin(10) * t810 + t1205;
t668 = -pkin(5) * t756 + pkin(10) * t781 - t1199;
t667 = -pkin(7) * t821 - t1128 * t748 + t1133 * t737;
t666 = -pkin(9) * t774 - t1126 * t736 + t1131 * t771;
t665 = -pkin(3) * t1229 + pkin(9) * t802 + t1126 * t794 + t1131 * t744;
t664 = -pkin(2) * t687 - pkin(3) * t750;
t663 = t1129 * t728 + t1134 * t702;
t662 = t1129 * t702 - t1134 * t728;
t661 = -pkin(7) * t815 - t1128 * t745 + t1133 * t727;
t660 = -t1127 * t722 + t1132 * t724;
t659 = -t1127 * t721 + t1132 * t723;
t658 = t1127 * t724 + t1132 * t722;
t657 = -t1127 * t723 - t1132 * t721;
t656 = -t1126 * t719 + t1131 * t720;
t655 = t1126 * t720 + t1131 * t719;
t652 = -pkin(3) * t859 + pkin(9) * t775 + t1126 * t771 + t1131 * t736;
t651 = -t1127 * t705 + t1132 * t706;
t650 = -t1127 * t706 - t1132 * t705;
t649 = -t1128 * t936 + t1133 * t659;
t648 = -t1070 * t1128 + t1133 * t651;
t647 = t1128 * t892 + t1133 * t660;
t646 = t1128 * t660 - t1133 * t892;
t645 = t1129 * t713 + t1134 * t686;
t644 = t1129 * t686 - t1134 * t713;
t643 = -pkin(8) * t790 - t1127 * t725 + t1132 * t726;
t642 = -pkin(4) * t831 + qJ(5) * t654;
t641 = -t1122 * t695 + t1123 * t697;
t640 = -t1122 * t694 + t1123 * t696;
t639 = t1122 * t697 + t1123 * t695;
t638 = t1122 * t696 + t1123 * t694;
t637 = -pkin(2) * t728 - pkin(3) * t801 - pkin(4) * t878 + t712;
t636 = -t1126 * t691 + t1131 * t693;
t635 = -t1126 * t690 + t1131 * t692;
t634 = t1126 * t693 + t1131 * t691;
t633 = t1126 * t692 + t1131 * t690;
t632 = -qJ(5) * t786 - t653;
t631 = -pkin(8) * t687 - pkin(9) * t1194 - t1127 * t743;
t630 = t1129 * t687 + t1134 * t681;
t629 = t1129 * t681 - t1134 * t687;
t628 = -pkin(4) * t892 + qJ(5) * t788 + t654;
t624 = -pkin(2) * t713 - pkin(3) * t774 - pkin(4) * t844 + t711;
t623 = -pkin(1) * t762 + pkin(2) * t961 - pkin(8) * t792 - t1127 * t726 - t1132 * t725;
t622 = -t1127 * t674 + t1132 * t676;
t621 = -t1127 * t673 + t1132 * t675;
t620 = -t1127 * t676 - t1132 * t674;
t619 = -t1127 * t675 - t1132 * t673;
t618 = -t1127 * t671 + t1132 * t672;
t617 = t1127 * t672 + t1132 * t671;
t616 = -pkin(2) * t658 - pkin(3) * t722 - pkin(4) * t786;
t615 = -qJ(5) * t734 - t1122 * t669 + t1123 * t709;
t614 = -pkin(7) * t762 - t1128 * t742 + t1133 * t643;
t613 = -t1128 * t757 + t1133 * t622;
t612 = -t1128 * t761 + t1133 * t621;
t611 = t1128 * t1230 + t1133 * t618;
t610 = t1128 * t618 - t1133 * t1230;
t609 = -qJ(5) * t719 - t1122 * t668 + t1123 * t698;
t608 = -t1127 * t655 + t1132 * t656;
t607 = t1127 * t656 + t1132 * t655;
t606 = -pkin(4) * t1230 + qJ(5) * t735 + t1122 * t709 + t1123 * t669;
t605 = -pkin(1) * t680 + pkin(2) * t943 - pkin(8) * t688 + pkin(9) * t1200 - t1132 * t743;
t604 = t1131 * t654 - t1203;
t603 = t1126 * t654 + t1197;
t602 = t1129 * t658 + t1134 * t647;
t601 = t1129 * t647 - t1134 * t658;
t600 = -pkin(8) * t728 - t1127 * t665 + t1132 * t670;
t599 = -pkin(4) * t756 + qJ(5) * t720 + t1122 * t698 + t1123 * t668;
t598 = t1128 * t756 + t1133 * t608;
t597 = t1128 * t608 - t1133 * t756;
t596 = -t1126 * t639 + t1131 * t641;
t595 = -t1126 * t638 + t1131 * t640;
t594 = t1126 * t641 + t1131 * t639;
t593 = t1126 * t640 + t1131 * t638;
t592 = -pkin(8) * t713 - t1127 * t652 + t1132 * t666;
t591 = -t1127 * t634 + t1132 * t636;
t590 = -t1127 * t633 + t1132 * t635;
t589 = -t1127 * t636 - t1132 * t634;
t588 = -t1127 * t635 - t1132 * t633;
t587 = t1133 * t591 + t1162;
t586 = t1133 * t590 - t1162;
t583 = -pkin(1) * t701 + pkin(2) * t1229 - pkin(8) * t729 - t1127 * t670 - t1132 * t665;
t582 = -pkin(7) * t680 - t1128 * t664 + t1133 * t631;
t581 = -pkin(5) * t749 + pkin(10) * t585;
t580 = -pkin(9) * t722 - t1126 * t628 + t1131 * t632;
t579 = -pkin(1) * t685 + pkin(2) * t859 - pkin(8) * t714 - t1127 * t666 - t1132 * t652;
t578 = -pkin(3) * t892 + pkin(9) * t724 + t1126 * t632 + t1131 * t628;
t577 = -pkin(10) * t695 - t584;
t576 = t1129 * t617 + t1134 * t611;
t575 = t1129 * t611 - t1134 * t617;
t574 = -pkin(5) * t818 + pkin(10) * t697 + t585;
t573 = -pkin(7) * t701 - t1128 * t637 + t1133 * t600;
t572 = -t1127 * t603 + t1132 * t604;
t571 = t1127 * t604 + t1132 * t603;
t570 = -pkin(9) * t603 - qJ(5) * t1197 - t1126 * t642;
t569 = t1128 * t831 + t1133 * t572;
t568 = t1128 * t572 - t1133 * t831;
t567 = -pkin(2) * t617 - pkin(3) * t671 - pkin(4) * t734 - pkin(5) * t809 + t627;
t566 = t1129 * t607 + t1134 * t598;
t565 = t1129 * t598 - t1134 * t607;
t564 = -pkin(3) * t831 + pkin(9) * t604 - qJ(5) * t1203 + t1131 * t642;
t563 = -pkin(7) * t685 - t1128 * t624 + t1133 * t592;
t562 = -pkin(9) * t671 - t1126 * t606 + t1131 * t615;
t561 = -pkin(3) * t1230 + pkin(9) * t672 + t1126 * t615 + t1131 * t606;
t560 = -t1127 * t594 + t1132 * t596;
t559 = -t1127 * t593 + t1132 * t595;
t558 = t1127 * t596 + t1132 * t594;
t557 = -t1127 * t595 - t1132 * t593;
t556 = -t1128 * t865 + t1133 * t559;
t555 = t1128 * t818 + t1133 * t560;
t554 = t1128 * t560 - t1133 * t818;
t553 = -pkin(2) * t607 - pkin(3) * t655 - pkin(4) * t719 - pkin(5) * t780 + t626;
t552 = -pkin(9) * t655 - t1126 * t599 + t1131 * t609;
t551 = t1123 * t585 - t1211;
t550 = t1122 * t585 + t1208;
t549 = -pkin(3) * t756 + pkin(9) * t656 + t1126 * t609 + t1131 * t599;
t548 = -pkin(2) * t571 - pkin(3) * t603 - pkin(4) * t653;
t547 = -pkin(8) * t658 - t1127 * t578 + t1132 * t580;
t546 = -qJ(5) * t639 - t1122 * t574 + t1123 * t577;
t545 = -pkin(4) * t818 + qJ(5) * t641 + t1122 * t577 + t1123 * t574;
t544 = t1129 * t571 + t1134 * t569;
t543 = t1129 * t569 - t1134 * t571;
t542 = -pkin(1) * t646 + pkin(2) * t892 - pkin(8) * t660 - t1127 * t580 - t1132 * t578;
t541 = -pkin(2) * t558 - pkin(3) * t594 - pkin(4) * t639 - pkin(5) * t695;
t540 = t1129 * t558 + t1134 * t555;
t539 = t1129 * t555 - t1134 * t558;
t538 = -pkin(7) * t646 - t1128 * t616 + t1133 * t547;
t537 = -t1126 * t550 + t1131 * t551;
t536 = t1126 * t551 + t1131 * t550;
t535 = -pkin(8) * t617 - t1127 * t561 + t1132 * t562;
t534 = -pkin(10) * t1208 - qJ(5) * t550 - t1122 * t581;
t533 = -pkin(4) * t749 - pkin(10) * t1211 + qJ(5) * t551 + t1123 * t581;
t532 = -pkin(8) * t607 - t1127 * t549 + t1132 * t552;
t531 = -pkin(8) * t571 - t1127 * t564 + t1132 * t570;
t530 = -pkin(1) * t610 + pkin(2) * t1230 - pkin(8) * t618 - t1127 * t562 - t1132 * t561;
t529 = -pkin(1) * t597 + pkin(2) * t756 - pkin(8) * t608 - t1127 * t552 - t1132 * t549;
t528 = -pkin(1) * t568 + pkin(2) * t831 - pkin(8) * t572 - t1127 * t570 - t1132 * t564;
t527 = -pkin(7) * t610 - t1128 * t567 + t1133 * t535;
t526 = -pkin(9) * t594 - t1126 * t545 + t1131 * t546;
t525 = -pkin(3) * t818 + pkin(9) * t596 + t1126 * t546 + t1131 * t545;
t524 = -pkin(7) * t597 - t1128 * t553 + t1133 * t532;
t523 = -t1127 * t536 + t1132 * t537;
t522 = t1127 * t537 + t1132 * t536;
t521 = t1128 * t749 + t1133 * t523;
t520 = t1128 * t523 - t1133 * t749;
t519 = -pkin(7) * t568 - t1128 * t548 + t1133 * t531;
t518 = -pkin(9) * t536 - t1126 * t533 + t1131 * t534;
t517 = -pkin(3) * t749 + pkin(9) * t537 + t1126 * t534 + t1131 * t533;
t516 = -pkin(2) * t522 - pkin(3) * t536 - pkin(4) * t550 - pkin(5) * t584;
t515 = -pkin(8) * t558 - t1127 * t525 + t1132 * t526;
t514 = t1129 * t522 + t1134 * t521;
t513 = t1129 * t521 - t1134 * t522;
t512 = -pkin(1) * t554 + pkin(2) * t818 - pkin(8) * t560 - t1127 * t526 - t1132 * t525;
t511 = -pkin(7) * t554 - t1128 * t541 + t1133 * t515;
t510 = -pkin(8) * t522 - t1127 * t517 + t1132 * t518;
t509 = -pkin(1) * t520 + pkin(2) * t749 - pkin(8) * t523 - t1127 * t518 - t1132 * t517;
t508 = -pkin(7) * t520 - t1128 * t516 + t1133 * t510;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t1090, -t1091, 0, t1051, 0, 0, 0, 0, 0, 0, t1025, t1026, t1047, t979, 0, 0, 0, 0, 0, 0, t889, t894, t868, t808, 0, 0, 0, 0, 0, 0, t747, t765, t704, t630, 0, 0, 0, 0, 0, 0, t645, t663, t602, t544, 0, 0, 0, 0, 0, 0, t566, t576, t540, t514; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1091, -t1090, 0, t1050, 0, 0, 0, 0, 0, 0, t1023, t1024, t1046, t978, 0, 0, 0, 0, 0, 0, t888, t893, t867, t807, 0, 0, 0, 0, 0, 0, t746, t764, t703, t629, 0, 0, 0, 0, 0, 0, t644, t662, t601, t543, 0, 0, 0, 0, 0, 0, t565, t575, t539, t513; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1052, t1053, 0, -t1010, 0, 0, 0, 0, 0, 0, t937, t941, t908, t872, 0, 0, 0, 0, 0, 0, t815, t821, t762, t680, 0, 0, 0, 0, 0, 0, t685, t701, t646, t568, 0, 0, 0, 0, 0, 0, t597, t610, t554, t520; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1091, 0, -t1090, 0, t1147, -t1069, -t1050, -pkin(6) * t1050, t1062 * t1134 - t1150, t1044 * t1134 - t1094 * t1129, t1056 * t1134 + t1128 * t1165, t1061 * t1134 + t1150, t1054 * t1134 + t1129 * t1163, qJDD(2) * t1129 + t1080 * t1134, -pkin(6) * t1023 - t1013 * t1129 + t1020 * t1134, -pkin(6) * t1024 - t1014 * t1129 + t1021 * t1134, -pkin(6) * t1046 + t1010 * t1134, -pkin(6) * t978 - (pkin(1) * t1129 - pkin(7) * t1134) * t1010, -t1129 * t998 + t1134 * t963, -t1129 * t952 + t1134 * t935, -t1129 * t990 + t1134 * t949, -t1129 * t996 + t1134 * t962, -t1129 * t991 + t1134 * t950, -t1016 * t1129 + t1134 * t995, -pkin(6) * t888 - t1129 * t837 + t1134 * t817, -pkin(6) * t893 - t1129 * t853 + t1134 * t820, -pkin(6) * t867 - t1129 * t800 + t1134 * t793, -pkin(6) * t807 - t1129 * t799 + t1134 * t766, -t1129 * t841 + t1134 * t830, -t1129 * t789 + t1134 * t768, -t1129 * t874 + t1134 * t824, -t1129 * t840 + t1134 * t829, -t1129 * t875 + t1134 * t825, -t1129 * t902 + t1134 * t898, -pkin(6) * t746 - t1129 * t682 + t1134 * t661, -pkin(6) * t764 - t1129 * t689 + t1134 * t667, -pkin(6) * t703 - t1129 * t623 + t1134 * t614, -pkin(6) * t629 - t1129 * t605 + t1134 * t582, -t1129 * t716 + t1134 * t700, -t1129 * t657 + t1134 * t649, -t1129 * t730 + t1134 * t707, -t1129 * t715 + t1134 * t699, -t1129 * t731 + t1134 * t708, -t1129 * t769 + t1134 * t767, -pkin(6) * t644 - t1129 * t579 + t1134 * t563, -pkin(6) * t662 - t1129 * t583 + t1134 * t573, -pkin(6) * t601 - t1129 * t542 + t1134 * t538, -pkin(6) * t543 - t1129 * t528 + t1134 * t519, -t1129 * t589 + t1134 * t587, -t1129 * t557 + t1134 * t556, -t1129 * t619 + t1134 * t612, -t1129 * t588 + t1134 * t586, -t1129 * t620 + t1134 * t613, -t1129 * t650 + t1134 * t648, -pkin(6) * t565 - t1129 * t529 + t1134 * t524, -pkin(6) * t575 - t1129 * t530 + t1134 * t527, -pkin(6) * t539 - t1129 * t512 + t1134 * t511, -pkin(6) * t513 - t1129 * t509 + t1134 * t508; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1090, 0, t1091, 0, t1069, t1147, t1051, pkin(6) * t1051, t1062 * t1129 + t1149, t1044 * t1129 + t1094 * t1134, t1056 * t1129 - t1128 * t1164, t1061 * t1129 - t1149, t1054 * t1129 - t1134 * t1163, -qJDD(2) * t1134 + t1080 * t1129, pkin(6) * t1025 + t1013 * t1134 + t1020 * t1129, pkin(6) * t1026 + t1014 * t1134 + t1021 * t1129, pkin(6) * t1047 + t1010 * t1129, pkin(6) * t979 - (-pkin(1) * t1134 - pkin(7) * t1129) * t1010, t1129 * t963 + t1134 * t998, t1129 * t935 + t1134 * t952, t1129 * t949 + t1134 * t990, t1129 * t962 + t1134 * t996, t1129 * t950 + t1134 * t991, t1016 * t1134 + t1129 * t995, pkin(6) * t889 + t1129 * t817 + t1134 * t837, pkin(6) * t894 + t1129 * t820 + t1134 * t853, pkin(6) * t868 + t1129 * t793 + t1134 * t800, pkin(6) * t808 + t1129 * t766 + t1134 * t799, t1129 * t830 + t1134 * t841, t1129 * t768 + t1134 * t789, t1129 * t824 + t1134 * t874, t1129 * t829 + t1134 * t840, t1129 * t825 + t1134 * t875, t1129 * t898 + t1134 * t902, pkin(6) * t747 + t1129 * t661 + t1134 * t682, pkin(6) * t765 + t1129 * t667 + t1134 * t689, pkin(6) * t704 + t1129 * t614 + t1134 * t623, pkin(6) * t630 + t1129 * t582 + t1134 * t605, t1129 * t700 + t1134 * t716, t1129 * t649 + t1134 * t657, t1129 * t707 + t1134 * t730, t1129 * t699 + t1134 * t715, t1129 * t708 + t1134 * t731, t1129 * t767 + t1134 * t769, pkin(6) * t645 + t1129 * t563 + t1134 * t579, pkin(6) * t663 + t1129 * t573 + t1134 * t583, pkin(6) * t602 + t1129 * t538 + t1134 * t542, pkin(6) * t544 + t1129 * t519 + t1134 * t528, t1129 * t587 + t1134 * t589, t1129 * t556 + t1134 * t557, t1129 * t612 + t1134 * t619, t1129 * t586 + t1134 * t588, t1129 * t613 + t1134 * t620, t1129 * t648 + t1134 * t650, pkin(6) * t566 + t1129 * t524 + t1134 * t529, pkin(6) * t576 + t1129 * t527 + t1134 * t530, pkin(6) * t540 + t1129 * t511 + t1134 * t512, pkin(6) * t514 + t1129 * t508 + t1134 * t509; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1099, t1100, 0, 0, t1146 * t1128, t1085 * t1133 + t1088 * t1128, t1104 * t1133 + t1178, -t1145 * t1133, t1106 * t1128 + t1179, 0, pkin(1) * t1088 + pkin(7) * t1055 + t1184, -pkin(1) * t1085 + pkin(7) * t1057 - t1185, pkin(1) * t1093 + pkin(7) * t1089 + t1011, pkin(1) * t1072 + pkin(7) * t1011, t1128 * t999 + t1141, t1048 * t1133 + t1128 * t954, t1008 * t1133 + t1128 * t992, t1128 * t997 - t1141, t1003 * t1133 + t1128 * t993, t1017 * t1128 + t1079 * t1133, -pkin(1) * t975 + pkin(7) * t938 + t1128 * t944 + t1133 * t900, -pkin(1) * t983 + pkin(7) * t942 + t1128 * t951 + t1133 * t904, pkin(7) * t909 + t1128 * t834 + (-pkin(1) - t1216) * t953, pkin(7) * t873 + (-pkin(1) + t1148) * t890, t1128 * t843 - t1157, t1128 * t791 + t1133 * t994, t1128 * t876 + t1133 * t931, t1128 * t842 + t1157, t1128 * t877 + t1133 * t927, t1128 * t903 + t1067, -pkin(1) * t835 + pkin(7) * t816 + t1128 * t727 + t1133 * t745, -pkin(1) * t869 + pkin(7) * t822 + t1128 * t737 + t1133 * t748, -pkin(1) * t790 + pkin(7) * t763 + t1128 * t643 + t1133 * t742, -pkin(1) * t687 + pkin(7) * t681 + t1128 * t631 + t1133 * t664, t1128 * t718 - t1159, t1128 * t659 + t1133 * t936, t1128 * t732 + t1133 * t864, t1128 * t717 + t1159, t1128 * t733 + t1133 * t860, t1128 * t770 + t1067, -pkin(1) * t713 + pkin(7) * t686 + t1128 * t592 + t1133 * t624, -pkin(1) * t728 + pkin(7) * t702 + t1128 * t600 + t1133 * t637, -pkin(1) * t658 + pkin(7) * t647 + t1128 * t547 + t1133 * t616, -pkin(1) * t571 + pkin(7) * t569 + t1128 * t531 + t1133 * t548, t1128 * t591 - t1160, t1128 * t559 + t1133 * t865, t1128 * t621 + t1133 * t761, t1128 * t590 + t1160, t1128 * t622 + t1133 * t757, t1070 * t1133 + t1128 * t651, -pkin(1) * t607 + pkin(7) * t598 + t1128 * t532 + t1133 * t553, -pkin(1) * t617 + pkin(7) * t611 + t1128 * t535 + t1133 * t567, -pkin(1) * t558 + pkin(7) * t555 + t1128 * t515 + t1133 * t541, -pkin(1) * t522 + pkin(7) * t521 + t1128 * t510 + t1133 * t516;];
tauB_reg  = t1;