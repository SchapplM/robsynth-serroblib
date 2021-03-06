% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RRRPPR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 04:48
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RRRPPR4_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR4_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR4_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_invdynB_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:46:49
% EndTime: 2019-05-07 04:47:44
% DurationCPUTime: 50.80s
% Computational Cost: add. (153933->801), mult. (322412->1153), div. (0->0), fcn. (231960->10), ass. (0->573)
t1068 = sin(qJ(1));
t1072 = cos(qJ(1));
t1067 = sin(qJ(2));
t1071 = cos(qJ(2));
t1062 = sin(pkin(10));
t1063 = cos(pkin(10));
t1066 = sin(qJ(3));
t1070 = cos(qJ(3));
t1137 = qJD(1) * t1067;
t1026 = qJD(2) * t1066 + t1070 * t1137;
t1136 = qJD(1) * t1071;
t1106 = qJD(2) * t1136;
t1117 = qJDD(1) * t1067;
t1083 = t1106 + t1117;
t1076 = t1070 * qJDD(2) - t1066 * t1083;
t1075 = t1026 * qJD(3) - t1076;
t1025 = -t1070 * qJD(2) + t1066 * t1137;
t1105 = t1025 * qJD(3) - t1066 * qJDD(2) - t1070 * t1083;
t1167 = -t1062 * t1075 - t1063 * t1105;
t1050 = -qJD(3) + t1136;
t984 = t1063 * t1025 + t1026 * t1062;
t965 = t984 * t1050;
t1194 = t965 + t1167;
t1053 = qJD(2) * t1137;
t1115 = qJDD(1) * t1071;
t1030 = -t1053 + t1115;
t1022 = -qJDD(3) + t1030;
t986 = -t1025 * t1062 + t1026 * t1063;
t1170 = t984 * t986;
t1080 = t1022 - t1170;
t1154 = t1063 * t1080;
t1176 = t1050 ^ 2;
t983 = t986 ^ 2;
t946 = -t983 - t1176;
t821 = t1062 * t946 - t1154;
t1159 = t1062 * t1080;
t832 = t1063 * t946 + t1159;
t731 = t1066 * t832 + t1070 * t821;
t709 = -t1067 * t1194 + t1071 * t731;
t729 = t1066 * t821 - t1070 * t832;
t647 = t1068 * t709 - t1072 * t729;
t1306 = pkin(6) * t647;
t649 = t1068 * t729 + t1072 * t709;
t1305 = pkin(6) * t649;
t707 = t1067 * t731 + t1071 * t1194;
t1304 = pkin(7) * t707;
t1303 = -pkin(1) * t707 - pkin(2) * t1194 - pkin(8) * t731;
t1302 = -pkin(1) * t729 + pkin(7) * t709;
t1179 = t984 ^ 2;
t959 = t1179 - t1176;
t842 = t1062 * t959 - t1154;
t846 = t1063 * t959 + t1159;
t768 = t1066 * t842 - t1070 * t846;
t1140 = t986 * t1050;
t902 = -t1062 * t1105 + t1063 * t1075;
t860 = t902 + t1140;
t725 = t1067 * t860 + t1071 * t768;
t764 = t1066 * t846 + t1070 * t842;
t1301 = t1068 * t725 + t1072 * t764;
t1086 = t902 - t1140;
t786 = -t1062 * t1086 + t1063 * t1194;
t1161 = t1062 * t1194;
t788 = t1063 * t1086 + t1161;
t702 = t1066 * t786 + t1070 * t788;
t919 = t983 - t1179;
t683 = -t1067 * t919 + t1071 * t702;
t700 = t1066 * t788 - t1070 * t786;
t1300 = t1068 * t683 - t1072 * t700;
t1299 = -t1068 * t764 + t1072 * t725;
t1298 = t1068 * t700 + t1072 * t683;
t1296 = pkin(8) * t729;
t1289 = -pkin(2) * t729 + pkin(3) * t832;
t1287 = t1067 * t768 - t1071 * t860;
t1286 = t1067 * t702 + t1071 * t919;
t1103 = t1022 + t1170;
t1158 = t1062 * t1103;
t1187 = -t1176 - t1179;
t1203 = t1063 * t1187 + t1158;
t1193 = t1063 * t1103;
t1206 = t1062 * t1187 - t1193;
t1224 = t1066 * t1203 + t1070 * t1206;
t1225 = -t1066 * t1206 + t1070 * t1203;
t1245 = t1067 * t1086 + t1071 * t1225;
t1263 = t1068 * t1224 + t1072 * t1245;
t1285 = pkin(6) * t1263;
t1195 = t965 - t1167;
t1205 = -t1062 * t1195 - t1063 * t860;
t1207 = -t1062 * t860 + t1063 * t1195;
t1222 = t1066 * t1205 + t1070 * t1207;
t1223 = -t1066 * t1207 + t1070 * t1205;
t876 = t983 + t1179;
t1244 = -t1067 * t876 + t1071 * t1223;
t1264 = t1068 * t1222 + t1072 * t1244;
t1284 = pkin(6) * t1264;
t1266 = t1068 * t1245 - t1072 * t1224;
t1283 = pkin(6) * t1266;
t1267 = t1068 * t1244 - t1072 * t1222;
t1282 = pkin(6) * t1267;
t960 = -t983 + t1176;
t1196 = t1063 * t960 - t1158;
t1231 = -t1062 * t960 - t1193;
t1249 = -t1066 * t1231 - t1070 * t1196;
t1246 = -t1066 * t1196 + t1070 * t1231;
t1265 = -t1067 * t1195 + t1071 * t1246;
t1281 = t1068 * t1265 + t1072 * t1249;
t1280 = -t1068 * t1249 + t1072 * t1265;
t1247 = t1067 * t1223 + t1071 * t876;
t1278 = pkin(7) * t1247;
t1248 = t1067 * t1225 - t1071 * t1086;
t1277 = pkin(7) * t1248;
t1276 = qJ(4) * t821;
t1275 = qJ(4) * t832;
t1272 = -pkin(1) * t1247 - pkin(2) * t876 - pkin(8) * t1223;
t1271 = -pkin(1) * t1248 + pkin(2) * t1086 - pkin(8) * t1225;
t1270 = -pkin(1) * t1222 + pkin(7) * t1244;
t1269 = -pkin(1) * t1224 + pkin(7) * t1245;
t1268 = t1067 * t1246 + t1071 * t1195;
t1260 = pkin(8) * t1222;
t1259 = pkin(8) * t1224;
t654 = -pkin(2) * t1222 - pkin(3) * t1207;
t1252 = -pkin(2) * t1224 - pkin(3) * t1206;
t1240 = qJ(4) * t1203;
t1239 = qJ(4) * t1206;
t1238 = qJ(4) * t1207;
t1237 = qJ(5) * t1194;
t1001 = -pkin(3) * t1050 - qJ(4) * t1026;
t1177 = t1025 ^ 2;
t1041 = g(1) * t1072 + g(2) * t1068;
t1073 = qJD(1) ^ 2;
t1014 = -pkin(1) * t1073 + qJDD(1) * pkin(7) - t1041;
t1173 = pkin(2) * t1071;
t1093 = -pkin(8) * t1067 - t1173;
t1028 = t1093 * qJD(1);
t1171 = t1071 * g(3);
t1178 = qJD(2) ^ 2;
t949 = (qJD(1) * t1028 + t1014) * t1067 - qJDD(2) * pkin(2) - t1178 * pkin(8) + t1171;
t838 = t1075 * pkin(3) - t1177 * qJ(4) + t1026 * t1001 + qJDD(4) + t949;
t1232 = t902 * pkin(4) - t1237 + t838;
t1082 = (t1062 * t984 + t1063 * t986) * t1050;
t1123 = t1050 * t1063;
t1109 = t984 * t1123;
t1124 = t1050 * t1062;
t954 = t986 * t1124;
t1090 = -t954 + t1109;
t1183 = -t1066 * t1090 - t1070 * t1082;
t1132 = t1022 * t1067;
t1182 = -t1066 * t1082 + t1070 * t1090;
t1200 = t1071 * t1182 - t1132;
t1230 = t1068 * t1200 + t1072 * t1183;
t1085 = t1062 * t902 - t1109;
t1091 = -t1063 * t902 - t1124 * t984;
t1180 = -t1066 * t1085 - t1070 * t1091;
t1112 = t1067 * t1170;
t1181 = -t1066 * t1091 + t1070 * t1085;
t1201 = t1071 * t1181 - t1112;
t1229 = t1068 * t1201 + t1072 * t1180;
t1228 = -t1068 * t1183 + t1072 * t1200;
t1227 = -t1068 * t1180 + t1072 * t1201;
t1226 = pkin(3) * t876 + qJ(4) * t1205;
t1221 = -2 * qJD(5);
t1065 = sin(qJ(6));
t1017 = qJDD(6) + t1022;
t1069 = cos(qJ(6));
t914 = t1065 * t986 - t1069 * t984;
t916 = t1065 * t984 + t1069 * t986;
t829 = t916 * t914;
t1197 = t1017 - t829;
t1219 = t1065 * t1197;
t1213 = t1069 * t1197;
t1008 = t1025 * t1050;
t938 = t1008 + t1105;
t1011 = t1071 * t1022;
t1204 = t1067 * t1182 + t1011;
t1110 = t1071 * t1170;
t1202 = t1067 * t1181 + t1110;
t792 = -t916 * qJD(6) - t1065 * t1167 + t1069 * t902;
t1043 = qJD(6) + t1050;
t888 = t1043 * t916;
t1198 = t792 + t888;
t793 = -t914 * qJD(6) + t1065 * t902 + t1069 * t1167;
t887 = t1043 * t914;
t747 = -t887 + t793;
t1131 = t1025 * t1026;
t1079 = -t1022 - t1131;
t1192 = t1066 * t1079;
t1190 = t1070 * t1079;
t1040 = t1068 * g(1) - t1072 * g(2);
t1013 = qJDD(1) * pkin(1) + t1073 * pkin(7) + t1040;
t1029 = 0.2e1 * t1106 + t1117;
t1089 = -t1030 + t1053;
t932 = pkin(2) * t1089 - pkin(8) * t1029 - t1013;
t1003 = -g(3) * t1067 + t1071 * t1014;
t950 = -pkin(2) * t1178 + qJDD(2) * pkin(8) + t1028 * t1136 + t1003;
t871 = t1066 * t950 - t1070 * t932;
t803 = t1079 * pkin(3) + qJ(4) * t938 - t871;
t872 = t1066 * t932 + t1070 * t950;
t808 = -pkin(3) * t1177 - qJ(4) * t1075 + t1050 * t1001 + t872;
t1169 = t1062 * t803 + t1063 * t808;
t917 = pkin(4) * t984 - qJ(5) * t986;
t1188 = -t1022 * qJ(5) + t1050 * t1221 - t984 * t917 + t1169;
t937 = t1105 - t1008;
t854 = t1062 * t1167 - t1123 * t986;
t855 = t1063 * t1167 + t954;
t779 = -t1066 * t854 + t1070 * t855;
t1096 = t1071 * t779 + t1112;
t776 = -t1066 * t855 - t1070 * t854;
t1185 = t1068 * t1096 + t1072 * t776;
t933 = (qJD(3) + t1050) * t1026 - t1076;
t1184 = -t1068 * t776 + t1072 * t1096;
t912 = t914 ^ 2;
t913 = t916 ^ 2;
t1021 = t1026 ^ 2;
t1042 = t1043 ^ 2;
t1175 = pkin(4) + pkin(5);
t1174 = pkin(2) * t1067;
t1172 = pkin(4) * t1063;
t1168 = t1062 * t808 - t1063 * t803;
t1166 = qJ(5) * t1063;
t1165 = qJD(4) * t984;
t1164 = qJD(4) * t986;
t1163 = t1062 * t838;
t1156 = t1063 * t838;
t1097 = pkin(5) * t1050 - pkin(9) * t986;
t716 = (-pkin(4) * t1050 + t1221) * t986 + t1232;
t681 = t902 * pkin(5) + pkin(9) * t1179 - t1097 * t986 + t716;
t1153 = t1065 * t681;
t824 = t1017 + t829;
t1152 = t1065 * t824;
t978 = 0.2e1 * t1164;
t713 = t978 + t1168;
t976 = -0.2e1 * t1165;
t714 = t976 + t1169;
t643 = t1062 * t714 - t1063 * t713;
t1151 = t1066 * t643;
t1150 = t1066 * t949;
t967 = t1022 - t1131;
t1149 = t1066 * t967;
t1146 = t1069 * t681;
t1145 = t1069 * t824;
t1144 = t1070 * t643;
t1143 = t1070 * t949;
t1142 = t1070 * t967;
t1139 = -t1176 + t876;
t1138 = qJD(1) * qJD(2);
t1134 = t1013 * t1067;
t1133 = t1013 * t1071;
t1130 = t1029 * t1067;
t1049 = t1071 * t1073 * t1067;
t1037 = -t1049 + qJDD(2);
t1129 = t1037 * t1067;
t1128 = t1037 * t1071;
t1038 = qJDD(2) + t1049;
t1127 = t1038 * t1067;
t1126 = t1043 * t1065;
t1125 = t1043 * t1069;
t1122 = t1050 * t1066;
t1121 = t1050 * t1070;
t1059 = t1067 ^ 2;
t1120 = t1059 * t1073;
t1060 = t1071 ^ 2;
t1118 = t1059 + t1060;
t1116 = qJDD(1) * t1068;
t1114 = qJDD(1) * t1072;
t1113 = t1067 * t829;
t1111 = t1071 * t829;
t1108 = t1067 * t1131;
t1107 = t1071 * t1131;
t1104 = -qJ(5) * t1062 - pkin(3);
t644 = t1062 * t713 + t1063 * t714;
t1100 = t986 * t917 + qJDD(5) + t1168;
t1084 = t1022 * pkin(4) + t1100;
t1078 = -qJ(5) * t1176 + t1084;
t653 = pkin(5) * t1103 + pkin(9) * t1195 + t1078 + t978;
t1087 = t976 + t1188;
t682 = -pkin(4) * t1176 + t1087;
t655 = -pkin(5) * t1179 + t902 * pkin(9) - t1050 * t1097 + t682;
t598 = t1065 * t655 - t1069 * t653;
t1002 = t1067 * t1014 + t1171;
t941 = t1002 * t1067 + t1071 * t1003;
t994 = -t1040 * t1068 - t1072 * t1041;
t1099 = t1068 * t1049;
t1098 = t1072 * t1049;
t1095 = t1067 * t779 - t1110;
t1034 = -t1068 * t1073 + t1114;
t1092 = -pkin(6) * t1034 - g(3) * t1068;
t599 = t1065 * t653 + t1069 * t655;
t561 = t1065 * t599 - t1069 * t598;
t562 = t1065 * t598 + t1069 * t599;
t796 = t1066 * t872 - t1070 * t871;
t797 = t1066 * t871 + t1070 * t872;
t940 = t1002 * t1071 - t1003 * t1067;
t993 = t1040 * t1072 - t1041 * t1068;
t1074 = 0.2e1 * qJD(5) * t986 - t1232;
t1057 = t1060 * t1073;
t1047 = -t1057 - t1178;
t1046 = t1057 - t1178;
t1045 = -t1120 - t1178;
t1044 = -t1120 + t1178;
t1036 = t1057 - t1120;
t1035 = t1057 + t1120;
t1033 = t1072 * t1073 + t1116;
t1032 = t1118 * qJDD(1);
t1031 = -0.2e1 * t1053 + t1115;
t1024 = t1071 * t1038;
t1023 = t1118 * t1138;
t1010 = -pkin(6) * t1033 + g(3) * t1072;
t1007 = -t1021 + t1176;
t1006 = -t1176 + t1177;
t1005 = -t1059 * t1138 + t1071 * t1083;
t1004 = -t1030 * t1067 - t1060 * t1138;
t1000 = -t1045 * t1067 - t1128;
t999 = -t1044 * t1067 + t1024;
t998 = t1047 * t1071 - t1127;
t997 = t1046 * t1071 - t1129;
t996 = t1045 * t1071 - t1129;
t995 = t1047 * t1067 + t1024;
t991 = -t1021 + t1177;
t990 = t1032 * t1072 - t1035 * t1068;
t989 = t1032 * t1068 + t1035 * t1072;
t988 = -t1021 - t1176;
t987 = t1031 * t1071 - t1130;
t980 = -t1176 - t1177;
t966 = t1021 + t1177;
t958 = t1000 * t1072 + t1029 * t1068;
t957 = -t1031 * t1068 + t1072 * t998;
t956 = t1000 * t1068 - t1029 * t1072;
t955 = t1031 * t1072 + t1068 * t998;
t953 = -pkin(7) * t996 - t1133;
t952 = -pkin(7) * t995 - t1134;
t948 = (t1025 * t1070 - t1026 * t1066) * t1050;
t947 = (-t1025 * t1066 - t1026 * t1070) * t1050;
t943 = -pkin(1) * t996 + t1003;
t942 = -pkin(1) * t995 + t1002;
t934 = (-qJD(3) + t1050) * t1026 + t1076;
t929 = t1026 * t1122 - t1070 * t1105;
t928 = t1026 * t1121 + t1066 * t1105;
t927 = -t1025 * t1121 + t1066 * t1075;
t926 = t1025 * t1122 + t1070 * t1075;
t925 = t1071 * t948 - t1132;
t923 = t1006 * t1070 + t1149;
t922 = -t1007 * t1066 + t1190;
t921 = -t1006 * t1066 + t1142;
t920 = -t1007 * t1070 - t1192;
t911 = -t1066 * t988 + t1142;
t910 = t1070 * t988 + t1149;
t905 = -t1013 * t1068 + t1072 * t941;
t904 = t1013 * t1072 + t1068 * t941;
t894 = t1070 * t980 - t1192;
t893 = t1066 * t980 + t1190;
t885 = -t913 + t1042;
t884 = t912 - t1042;
t879 = -t913 - t1042;
t878 = t1071 * t929 + t1108;
t877 = t1071 * t927 - t1108;
t870 = -t1066 * t938 - t1070 * t933;
t869 = t1066 * t937 + t1070 * t934;
t868 = -t1066 * t933 + t1070 * t938;
t867 = -t1066 * t934 + t1070 * t937;
t856 = -pkin(8) * t910 + t1143;
t849 = -t1067 * t933 + t1071 * t923;
t848 = -t1067 * t938 + t1071 * t922;
t839 = -pkin(8) * t893 + t1150;
t837 = -t1067 * t937 + t1071 * t911;
t834 = t1067 * t911 + t1071 * t937;
t831 = -t1067 * t934 + t1071 * t894;
t830 = t1067 * t894 + t1071 * t934;
t828 = t913 - t912;
t827 = -t1067 * t991 + t1071 * t869;
t826 = -t1042 - t912;
t818 = -t1067 * t966 + t1071 * t870;
t817 = t1067 * t870 + t1071 * t966;
t816 = (t1065 * t916 - t1069 * t914) * t1043;
t815 = (t1065 * t914 + t1069 * t916) * t1043;
t810 = -pkin(2) * t910 + t872;
t809 = -pkin(2) * t893 + t871;
t800 = t1068 * t910 + t1072 * t837;
t799 = t1068 * t837 - t1072 * t910;
t798 = -t912 - t913;
t795 = t1068 * t893 + t1072 * t831;
t794 = t1068 * t831 - t1072 * t893;
t783 = t1069 * t884 - t1152;
t782 = -t1065 * t885 + t1213;
t781 = -t1065 * t884 - t1145;
t780 = -t1069 * t885 - t1219;
t773 = -t1065 * t879 - t1145;
t772 = t1069 * t879 - t1152;
t760 = t1067 * t949 + t1071 * t797;
t759 = t1067 * t797 - t1071 * t949;
t758 = t1156 - t1275;
t753 = t1068 * t868 + t1072 * t818;
t752 = t1068 * t818 - t1072 * t868;
t751 = t1163 - t1239;
t750 = -pkin(1) * t834 - pkin(2) * t937 - pkin(8) * t911 - t1150;
t748 = t887 + t793;
t745 = t792 - t888;
t743 = t1069 * t826 - t1219;
t742 = t1065 * t826 + t1213;
t741 = t1069 * t793 - t1126 * t916;
t740 = -t1065 * t793 - t1125 * t916;
t739 = -t1065 * t792 + t1125 * t914;
t738 = -t1069 * t792 - t1126 * t914;
t737 = -pkin(1) * t830 - pkin(2) * t934 - pkin(8) * t894 + t1143;
t728 = -pkin(8) * t868 - t796;
t727 = -t1062 * t815 + t1063 * t816;
t726 = t1062 * t816 + t1063 * t815;
t721 = -pkin(3) * t1194 + t1163 - t1276;
t715 = -pkin(3) * t1086 - t1156 + t1240;
t711 = -pkin(7) * t834 - t1067 * t810 + t1071 * t856;
t706 = -pkin(7) * t830 - t1067 * t809 + t1071 * t839;
t697 = t1068 * t796 + t1072 * t760;
t696 = t1068 * t760 - t1072 * t796;
t695 = -t1062 * t781 + t1063 * t783;
t694 = -t1062 * t780 + t1063 * t782;
t693 = t1062 * t783 + t1063 * t781;
t692 = t1062 * t782 + t1063 * t780;
t691 = t1062 * t772 + t1063 * t773;
t690 = t1062 * t773 - t1063 * t772;
t689 = (-t1086 + t1140) * pkin(4) + t1074;
t688 = pkin(4) * t1140 + t1074 + t1237;
t687 = -pkin(1) * t817 - pkin(2) * t966 - pkin(8) * t870 - t797;
t686 = -pkin(1) * t759 + pkin(2) * t949 - pkin(8) * t797;
t685 = -t1078 - 0.2e1 * t1164;
t676 = -pkin(7) * t817 + t1071 * t728 + t1174 * t868;
t675 = t1065 * t748 + t1069 * t1198;
t674 = -t1065 * t747 + t1069 * t745;
t673 = t1065 * t1198 - t1069 * t748;
t672 = -t1065 * t745 - t1069 * t747;
t671 = t1062 * t742 + t1063 * t743;
t670 = t1062 * t743 - t1063 * t742;
t669 = -t1062 * t740 + t1063 * t741;
t668 = -t1062 * t738 + t1063 * t739;
t667 = t1062 * t741 + t1063 * t740;
t666 = t1062 * t739 + t1063 * t738;
t665 = qJ(5) * t1139 + t1084 + t978;
t664 = pkin(4) * t1139 + t1087;
t659 = -t1066 * t726 + t1070 * t727;
t658 = -t1066 * t727 - t1070 * t726;
t657 = -t1017 * t1067 + t1071 * t659;
t656 = -pkin(7) * t759 + (-pkin(8) * t1071 + t1174) * t796;
t651 = -t1062 * t689 - t1086 * t1166 - t1239;
t646 = -pkin(4) * t1161 + t1063 * t688 + t1275;
t645 = -t1289 + t714;
t642 = t1063 * t689 + t1086 * t1104 + t1240;
t641 = t1252 + t713;
t640 = t1276 + t1062 * t688 - (-pkin(3) - t1172) * t1194;
t639 = -t1066 * t721 + t1070 * t758 + t1296;
t638 = -pkin(4) * t1195 + qJ(5) * t860 + t654;
t637 = -pkin(3) * t838 + qJ(4) * t644;
t636 = -t1066 * t693 + t1070 * t695;
t635 = -t1066 * t692 + t1070 * t694;
t634 = -t1066 * t695 - t1070 * t693;
t633 = -t1066 * t694 - t1070 * t692;
t632 = -t1066 * t690 + t1070 * t691;
t631 = t1066 * t691 + t1070 * t690;
t630 = -t1066 * t715 + t1070 * t751 - t1259;
t625 = t978 + (-t1176 - t1187) * qJ(5) + (t1022 + t1103) * pkin(4) + t1100 + t1252;
t624 = -t643 - t1238;
t623 = -t1062 * t685 + t1063 * t682;
t622 = t1062 * t682 + t1063 * t685;
t621 = -pkin(9) * t772 + qJ(5) * t747 - t1146;
t620 = qJ(5) * t1080 + 0.2e1 * t1165 + (t1176 + t946) * pkin(4) - t1188 + t1289;
t619 = t1226 + t644;
t618 = -pkin(9) * t742 - qJ(5) * t745 - t1153;
t617 = t1062 * t673 + t1063 * t675;
t616 = -t1062 * t672 + t1063 * t674;
t615 = t1062 * t675 - t1063 * t673;
t614 = t1062 * t674 + t1063 * t672;
t613 = -t1066 * t670 + t1070 * t671;
t612 = t1066 * t671 + t1070 * t670;
t611 = -t1066 * t667 + t1070 * t669;
t610 = -t1066 * t666 + t1070 * t668;
t609 = -t1066 * t669 - t1070 * t667;
t608 = -t1066 * t668 - t1070 * t666;
t607 = -t1067 * t1198 + t1071 * t636;
t606 = -t1067 * t748 + t1071 * t635;
t605 = -t1067 * t747 + t1071 * t632;
t604 = t1067 * t632 + t1071 * t747;
t603 = t1071 * t611 - t1113;
t602 = t1071 * t610 + t1113;
t601 = -t1066 * t758 - t1070 * t721 - t1303;
t600 = -pkin(9) * t773 + t1175 * t747 + t1153;
t597 = -pkin(9) * t743 - t1175 * t745 - t1146;
t596 = -t1062 * t664 + t1063 * t665 - t1238;
t595 = -t1066 * t751 - t1070 * t715 + t1271;
t594 = t1067 * t745 + t1071 * t613;
t593 = t1067 * t613 - t1071 * t745;
t592 = t1062 * t665 + t1063 * t664 + t1226;
t591 = t1070 * t644 - t1151;
t590 = t1066 * t644 + t1144;
t589 = t1067 * t838 + t1071 * t591;
t588 = t1067 * t591 - t1071 * t838;
t587 = -qJ(4) * t622 + (pkin(4) * t1062 - t1166) * t716;
t586 = -t1066 * t642 + t1070 * t651 - t1259;
t585 = -t1066 * t640 + t1070 * t646 - t1296;
t584 = -t1067 * t645 + t1071 * t639 + t1304;
t583 = -t1066 * t622 + t1070 * t623;
t582 = t1066 * t623 + t1070 * t622;
t581 = qJ(4) * t623 + (t1104 - t1172) * t716;
t580 = t1068 * t631 + t1072 * t605;
t579 = t1068 * t605 - t1072 * t631;
t578 = -t1067 * t641 + t1071 * t630 - t1277;
t577 = -t1066 * t651 - t1070 * t642 + t1271;
t576 = -t1066 * t615 + t1070 * t617;
t575 = -t1066 * t614 + t1070 * t616;
t574 = t1066 * t617 + t1070 * t615;
t573 = -t1066 * t616 - t1070 * t614;
t572 = -pkin(2) * t590 - pkin(3) * t643;
t571 = -t1066 * t646 - t1070 * t640 + t1303;
t570 = t1067 * t716 + t1071 * t583;
t569 = t1067 * t583 - t1071 * t716;
t568 = -t1067 * t828 + t1071 * t575;
t567 = -t1067 * t798 + t1071 * t576;
t566 = t1067 * t576 + t1071 * t798;
t565 = -t1066 * t619 + t1070 * t624 - t1260;
t564 = t1068 * t612 + t1072 * t594;
t563 = t1068 * t594 - t1072 * t612;
t560 = -qJ(4) * t690 - t1062 * t600 + t1063 * t621;
t559 = -qJ(4) * t670 - t1062 * t597 + t1063 * t618;
t558 = -t1067 * t625 + t1071 * t586 - t1277;
t557 = pkin(3) * t747 + qJ(4) * t691 + t1062 * t621 + t1063 * t600;
t556 = -t1066 * t624 - t1070 * t619 + t1272;
t555 = -pkin(9) * t561 - qJ(5) * t681;
t554 = -pkin(8) * t590 - qJ(4) * t1144 - t1066 * t637;
t553 = t1068 * t590 + t1072 * t589;
t552 = t1068 * t589 - t1072 * t590;
t551 = -t1067 * t620 + t1071 * t585 - t1304;
t550 = -t1066 * t592 + t1070 * t596 - t1260;
t549 = -pkin(3) * t745 + qJ(4) * t671 + t1062 * t618 + t1063 * t597;
t548 = -pkin(2) * t631 - pkin(3) * t690 - qJ(5) * t773 + t1175 * t772 - t599;
t547 = -pkin(9) * t673 + qJ(5) * t798 - t561;
t546 = -pkin(9) * t675 + t1175 * t798 - t562;
t545 = -t1067 * t654 + t1071 * t565 - t1278;
t544 = -pkin(9) * t562 - t1175 * t681;
t543 = -pkin(2) * t612 - pkin(3) * t670 - qJ(5) * t743 + t1175 * t742 - t598;
t542 = -t1066 * t596 - t1070 * t592 + t1272;
t541 = -pkin(2) * t582 - pkin(3) * t622 - pkin(4) * t685 - qJ(5) * t682;
t540 = t1068 * t582 + t1072 * t570;
t539 = t1068 * t570 - t1072 * t582;
t538 = -t1067 * t638 + t1071 * t550 - t1278;
t537 = t1068 * t574 + t1072 * t567;
t536 = t1068 * t567 - t1072 * t574;
t535 = -pkin(1) * t588 + pkin(2) * t838 - pkin(8) * t591 + qJ(4) * t1151 - t1070 * t637;
t534 = -pkin(2) * t574 - pkin(3) * t615 - qJ(5) * t675 + t1175 * t673;
t533 = t1062 * t561 + t1063 * t562;
t532 = t1062 * t562 - t1063 * t561;
t531 = -pkin(8) * t582 - t1066 * t581 + t1070 * t587;
t530 = -pkin(8) * t631 - t1066 * t557 + t1070 * t560;
t529 = -pkin(7) * t588 - t1067 * t572 + t1071 * t554;
t528 = -pkin(8) * t612 - t1066 * t549 + t1070 * t559;
t527 = -qJ(4) * t615 - t1062 * t546 + t1063 * t547;
t526 = pkin(3) * t798 + qJ(4) * t617 + t1062 * t547 + t1063 * t546;
t525 = -pkin(1) * t604 - pkin(2) * t747 - pkin(8) * t632 - t1066 * t560 - t1070 * t557;
t524 = -pkin(1) * t569 + pkin(2) * t716 - pkin(8) * t583 - t1066 * t587 - t1070 * t581;
t523 = -pkin(1) * t593 + pkin(2) * t745 - pkin(8) * t613 - t1066 * t559 - t1070 * t549;
t522 = -t1066 * t532 + t1070 * t533;
t521 = t1066 * t533 + t1070 * t532;
t520 = t1067 * t681 + t1071 * t522;
t519 = t1067 * t522 - t1071 * t681;
t518 = -pkin(7) * t604 - t1067 * t548 + t1071 * t530;
t517 = -qJ(4) * t532 - t1062 * t544 + t1063 * t555;
t516 = -pkin(7) * t569 - t1067 * t541 + t1071 * t531;
t515 = -pkin(3) * t681 + qJ(4) * t533 + t1062 * t555 + t1063 * t544;
t514 = -pkin(7) * t593 - t1067 * t543 + t1071 * t528;
t513 = -pkin(8) * t574 - t1066 * t526 + t1070 * t527;
t512 = -pkin(1) * t566 - pkin(2) * t798 - pkin(8) * t576 - t1066 * t527 - t1070 * t526;
t511 = t1068 * t521 + t1072 * t520;
t510 = t1068 * t520 - t1072 * t521;
t509 = -pkin(2) * t521 - pkin(3) * t532 - qJ(5) * t562 + t1175 * t561;
t508 = -pkin(7) * t566 - t1067 * t534 + t1071 * t513;
t507 = -pkin(8) * t521 - t1066 * t515 + t1070 * t517;
t506 = -pkin(1) * t519 + pkin(2) * t681 - pkin(8) * t522 - t1066 * t517 - t1070 * t515;
t505 = -pkin(7) * t519 - t1067 * t509 + t1071 * t507;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t1033, -t1034, 0, t994, 0, 0, 0, 0, 0, 0, t957, t958, t990, t905, 0, 0, 0, 0, 0, 0, t795, t800, t753, t697, 0, 0, 0, 0, 0, 0, t1263, -t649, t1264, t553, 0, 0, 0, 0, 0, 0, t1263, t1264, t649, t540, 0, 0, 0, 0, 0, 0, t564, t580, t537, t511; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1034, -t1033, 0, t993, 0, 0, 0, 0, 0, 0, t955, t956, t989, t904, 0, 0, 0, 0, 0, 0, t794, t799, t752, t696, 0, 0, 0, 0, 0, 0, t1266, -t647, t1267, t552, 0, 0, 0, 0, 0, 0, t1266, t1267, t647, t539, 0, 0, 0, 0, 0, 0, t563, t579, t536, t510; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t995, t996, 0, -t940, 0, 0, 0, 0, 0, 0, t830, t834, t817, t759, 0, 0, 0, 0, 0, 0, t1248, -t707, t1247, t588, 0, 0, 0, 0, 0, 0, t1248, t1247, t707, t569, 0, 0, 0, 0, 0, 0, t593, t604, t566, t519; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1034, 0, -t1033, 0, t1092, -t1010, -t993, -pkin(6) * t993, t1005 * t1072 - t1099, -t1036 * t1068 + t1072 * t987, t1067 * t1116 + t1072 * t999, t1004 * t1072 + t1099, t1068 * t1115 + t1072 * t997, qJDD(2) * t1068 + t1023 * t1072, -pkin(6) * t955 - t1068 * t942 + t1072 * t952, -pkin(6) * t956 - t1068 * t943 + t1072 * t953, -pkin(6) * t989 + t1072 * t940, -pkin(6) * t904 - (pkin(1) * t1068 - pkin(7) * t1072) * t940, -t1068 * t928 + t1072 * t878, -t1068 * t867 + t1072 * t827, -t1068 * t920 + t1072 * t848, -t1068 * t926 + t1072 * t877, -t1068 * t921 + t1072 * t849, -t1068 * t947 + t1072 * t925, -pkin(6) * t794 - t1068 * t737 + t1072 * t706, -pkin(6) * t799 - t1068 * t750 + t1072 * t711, -pkin(6) * t752 - t1068 * t687 + t1072 * t676, -pkin(6) * t696 - t1068 * t686 + t1072 * t656, t1184, -t1298, t1280, t1227, -t1299, t1228, -t1068 * t595 + t1072 * t578 - t1283, -t1068 * t601 + t1072 * t584 + t1306, -t1068 * t556 + t1072 * t545 - t1282, -pkin(6) * t552 - t1068 * t535 + t1072 * t529, t1184, t1280, t1298, t1228, t1299, t1227, -t1068 * t577 + t1072 * t558 - t1283, -t1068 * t542 + t1072 * t538 - t1282, -t1068 * t571 + t1072 * t551 - t1306, -pkin(6) * t539 - t1068 * t524 + t1072 * t516, -t1068 * t609 + t1072 * t603, -t1068 * t573 + t1072 * t568, -t1068 * t633 + t1072 * t606, -t1068 * t608 + t1072 * t602, -t1068 * t634 + t1072 * t607, -t1068 * t658 + t1072 * t657, -pkin(6) * t563 - t1068 * t523 + t1072 * t514, -pkin(6) * t579 - t1068 * t525 + t1072 * t518, -pkin(6) * t536 - t1068 * t512 + t1072 * t508, -pkin(6) * t510 - t1068 * t506 + t1072 * t505; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1033, 0, t1034, 0, t1010, t1092, t994, pkin(6) * t994, t1005 * t1068 + t1098, t1036 * t1072 + t1068 * t987, -t1067 * t1114 + t1068 * t999, t1004 * t1068 - t1098, t1068 * t997 - t1071 * t1114, -qJDD(2) * t1072 + t1023 * t1068, pkin(6) * t957 + t1068 * t952 + t1072 * t942, pkin(6) * t958 + t1068 * t953 + t1072 * t943, pkin(6) * t990 + t1068 * t940, pkin(6) * t905 - (-pkin(1) * t1072 - pkin(7) * t1068) * t940, t1068 * t878 + t1072 * t928, t1068 * t827 + t1072 * t867, t1068 * t848 + t1072 * t920, t1068 * t877 + t1072 * t926, t1068 * t849 + t1072 * t921, t1068 * t925 + t1072 * t947, pkin(6) * t795 + t1068 * t706 + t1072 * t737, pkin(6) * t800 + t1068 * t711 + t1072 * t750, pkin(6) * t753 + t1068 * t676 + t1072 * t687, pkin(6) * t697 + t1068 * t656 + t1072 * t686, t1185, -t1300, t1281, t1229, -t1301, t1230, t1068 * t578 + t1072 * t595 + t1285, t1068 * t584 + t1072 * t601 - t1305, t1068 * t545 + t1072 * t556 + t1284, pkin(6) * t553 + t1068 * t529 + t1072 * t535, t1185, t1281, t1300, t1230, t1301, t1229, t1068 * t558 + t1072 * t577 + t1285, t1068 * t538 + t1072 * t542 + t1284, t1068 * t551 + t1072 * t571 + t1305, pkin(6) * t540 + t1068 * t516 + t1072 * t524, t1068 * t603 + t1072 * t609, t1068 * t568 + t1072 * t573, t1068 * t606 + t1072 * t633, t1068 * t602 + t1072 * t608, t1068 * t607 + t1072 * t634, t1068 * t657 + t1072 * t658, pkin(6) * t564 + t1068 * t514 + t1072 * t523, pkin(6) * t580 + t1068 * t518 + t1072 * t525, pkin(6) * t537 + t1068 * t508 + t1072 * t512, pkin(6) * t511 + t1068 * t505 + t1072 * t506; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1040, t1041, 0, 0, t1130, t1029 * t1071 + t1031 * t1067, t1044 * t1071 + t1127, -t1089 * t1071, t1046 * t1067 + t1128, 0, pkin(1) * t1031 + pkin(7) * t998 + t1133, -pkin(1) * t1029 + pkin(7) * t1000 - t1134, pkin(1) * t1035 + pkin(7) * t1032 + t941, pkin(1) * t1013 + pkin(7) * t941, t1067 * t929 - t1107, t1067 * t869 + t1071 * t991, t1067 * t922 + t1071 * t938, t1067 * t927 + t1107, t1067 * t923 + t1071 * t933, t1067 * t948 + t1011, -pkin(1) * t893 + pkin(7) * t831 + t1067 * t839 + t1071 * t809, -pkin(1) * t910 + pkin(7) * t837 + t1067 * t856 + t1071 * t810, pkin(7) * t818 + t1067 * t728 + (-pkin(1) - t1173) * t868, pkin(7) * t760 + (-pkin(1) + t1093) * t796, t1095, -t1286, t1268, t1202, -t1287, t1204, t1067 * t630 + t1071 * t641 + t1269, t1067 * t639 + t1071 * t645 - t1302, t1067 * t565 + t1071 * t654 + t1270, -pkin(1) * t590 + pkin(7) * t589 + t1067 * t554 + t1071 * t572, t1095, t1268, t1286, t1204, t1287, t1202, t1067 * t586 + t1071 * t625 + t1269, t1067 * t550 + t1071 * t638 + t1270, t1067 * t585 + t1071 * t620 + t1302, -pkin(1) * t582 + pkin(7) * t570 + t1067 * t531 + t1071 * t541, t1067 * t611 + t1111, t1067 * t575 + t1071 * t828, t1067 * t635 + t1071 * t748, t1067 * t610 - t1111, t1067 * t636 + t1071 * t1198, t1017 * t1071 + t1067 * t659, -pkin(1) * t612 + pkin(7) * t594 + t1067 * t528 + t1071 * t543, -pkin(1) * t631 + pkin(7) * t605 + t1067 * t530 + t1071 * t548, -pkin(1) * t574 + pkin(7) * t567 + t1067 * t513 + t1071 * t534, -pkin(1) * t521 + pkin(7) * t520 + t1067 * t507 + t1071 * t509;];
tauB_reg  = t1;
