% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6RPPRRP6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 15:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6RPPRRP6_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP6_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_invdynm_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:01:41
% EndTime: 2019-05-05 15:01:57
% DurationCPUTime: 17.12s
% Computational Cost: add. (31402->659), mult. (59500->639), div. (0->0), fcn. (33547->6), ass. (0->400)
t1078 = sin(qJ(1));
t1081 = cos(qJ(1));
t1077 = sin(qJ(4));
t1080 = cos(qJ(4));
t1079 = cos(qJ(5));
t1076 = sin(qJ(5));
t1203 = qJD(1) * t1080;
t1058 = qJD(4) * t1203;
t1186 = t1077 * qJDD(1);
t1145 = t1058 + t1186;
t1018 = qJDD(5) + t1145;
t1024 = -t1079 * qJD(4) + t1076 * t1203;
t1026 = t1076 * qJD(4) + t1079 * t1203;
t1199 = t1026 * t1024;
t1265 = t1018 + t1199;
t1219 = t1076 * t1265;
t1194 = t1077 * qJD(1);
t1055 = qJD(5) + t1194;
t1244 = t1055 ^ 2;
t1245 = t1024 ^ 2;
t988 = t1245 - t1244;
t1126 = t1079 * t988 - t1219;
t1172 = qJD(4) * t1194;
t1185 = t1080 * qJDD(1);
t1030 = -t1172 + t1185;
t1148 = t1079 * qJDD(4) - t1076 * t1030;
t941 = qJD(5) * t1026 - t1148;
t996 = t1055 * t1026;
t915 = t941 - t996;
t822 = t1077 * t1126 + t1080 * t915;
t1210 = t1079 * t1265;
t882 = t1076 * t988 + t1210;
t1338 = t1078 * t822 + t1081 * t882;
t1200 = t1024 * t1055;
t1117 = -t1076 * qJDD(4) - t1079 * t1030;
t942 = -qJD(5) * t1024 - t1117;
t1273 = t942 - t1200;
t1221 = t1076 * t1273;
t1275 = t941 + t996;
t1129 = t1079 * t1275 + t1221;
t1017 = t1026 ^ 2;
t1268 = t1017 - t1245;
t792 = t1077 * t1129 + t1080 * t1268;
t1212 = t1079 * t1273;
t837 = -t1076 * t1275 + t1212;
t1337 = t1078 * t792 - t1081 * t837;
t1336 = -t1078 * t882 + t1081 * t822;
t1335 = t1078 * t837 + t1081 * t792;
t946 = t1244 + t1017;
t858 = t1079 * t946 + t1219;
t1334 = pkin(2) * t858;
t1333 = pkin(3) * t858;
t1332 = pkin(4) * t858;
t1331 = pkin(8) * t858;
t872 = t1076 * t946 - t1210;
t1330 = pkin(8) * t872;
t1329 = qJ(3) * t858;
t1328 = t1077 * t872;
t1326 = t1078 * t858;
t1324 = t1080 * t872;
t1322 = t1081 * t858;
t1237 = pkin(1) + qJ(3);
t1320 = t1237 * t858;
t1266 = t1018 - t1199;
t1218 = t1076 * t1266;
t989 = t1017 - t1244;
t1294 = -t1079 * t989 + t1218;
t1272 = t942 + t1200;
t1209 = t1079 * t1266;
t1293 = t1076 * t989 + t1209;
t1308 = -t1077 * t1293 + t1080 * t1272;
t1319 = -t1078 * t1308 + t1081 * t1294;
t825 = -t1077 * t915 + t1080 * t1126;
t793 = -t1077 * t1268 + t1080 * t1129;
t1318 = -t1078 * t1294 - t1081 * t1308;
t1264 = -t1244 - t1245;
t1277 = t1079 * t1264 - t1218;
t1292 = t1077 * t1277 - t1080 * t1275;
t1317 = pkin(3) * t1292;
t1316 = pkin(7) * t1292;
t1291 = t1077 * t1275 + t1080 * t1277;
t1315 = qJ(2) * t1291;
t1314 = qJ(3) * t1291;
t1311 = t1237 * t1291;
t1242 = pkin(2) + pkin(3);
t1310 = t1242 * t1292;
t1279 = t1076 * t1264 + t1209;
t1309 = -pkin(3) * t1279 + pkin(7) * t1291;
t1307 = t1077 * t1272 + t1080 * t1293;
t1236 = qJ(2) - pkin(7);
t1306 = t1236 * t1292 + t1237 * t1279;
t1305 = pkin(6) * (t1078 * t1292 + t1081 * t1279);
t1304 = pkin(6) * (-t1078 * t1279 + t1081 * t1292);
t1303 = pkin(2) * t1279;
t1301 = pkin(4) * t1279;
t1300 = pkin(8) * t1277;
t1299 = pkin(8) * t1279;
t1298 = qJ(3) * t1279;
t1267 = t1017 + t1245;
t1290 = pkin(4) * t1267;
t1289 = qJ(6) * t1076 + pkin(4);
t1288 = t1273 * qJ(6);
t1287 = t1077 * t1267;
t1283 = t1080 * t1267;
t1198 = t1055 * t1076;
t1173 = t1024 * t1198;
t1099 = -t1079 * t941 + t1173;
t1197 = t1055 * t1079;
t1100 = t1024 * t1197 + t1076 * t941;
t1174 = t1080 * t1199;
t1250 = t1077 * t1100 + t1174;
t1278 = t1078 * t1250 + t1081 * t1099;
t1276 = -t1078 * t1099 + t1081 * t1250;
t1084 = qJD(1) ^ 2;
t1074 = t1084 * pkin(7);
t1045 = t1078 * g(1) - t1081 * g(2);
t1134 = -qJDD(2) + t1045;
t1205 = t1084 * qJ(2);
t1095 = t1134 + t1205;
t1233 = pkin(4) * t1080;
t1139 = pkin(8) * t1077 + t1233;
t1234 = pkin(4) * t1077;
t1243 = 2 * qJD(3);
t884 = t1058 * pkin(4) - t1030 * pkin(8) - t1074 + (t1234 + t1237) * qJDD(1) + (qJD(4) * t1139 + t1243) * qJD(1) + t1095;
t1083 = qJD(4) ^ 2;
t1230 = pkin(8) * t1080;
t1102 = t1084 * (-t1230 + t1234);
t1070 = qJDD(1) * qJ(2);
t1046 = t1081 * g(1) + t1078 * g(2);
t1110 = 0.2e1 * qJD(2) * qJD(1) - t1046;
t1097 = qJDD(3) + t1110;
t1093 = t1070 + t1097;
t1269 = t1084 * t1237;
t964 = -qJDD(1) * pkin(7) + t1093 - t1269;
t944 = t1080 * g(3) - t1077 * t964;
t908 = -t1083 * pkin(4) + qJDD(4) * pkin(8) - t1077 * t1102 - t944;
t817 = t1076 * t908 - t1079 * t884;
t818 = t1076 * t884 + t1079 * t908;
t762 = t1076 * t817 + t1079 * t818;
t943 = t1077 * g(3) + t1080 * t964;
t865 = -t1077 * t943 - t1080 * t944;
t864 = -t1077 * t944 + t1080 * t943;
t1274 = t1242 * t864;
t1071 = t1077 ^ 2;
t1196 = t1071 * t1084;
t1050 = -t1083 - t1196;
t1053 = t1077 * t1084 * t1080;
t1042 = qJDD(4) - t1053;
t1191 = t1080 * t1042;
t974 = t1050 * t1077 + t1191;
t1263 = -t1242 * t974 - t943;
t1072 = t1080 ^ 2;
t1195 = t1072 * t1084;
t1051 = -t1083 - t1195;
t1041 = qJDD(4) + t1053;
t1193 = t1077 * t1041;
t976 = t1051 * t1080 - t1193;
t1262 = -t1242 * t976 - t944;
t896 = -t1026 * t1198 + t1079 * t942;
t1101 = t1077 * t896 - t1174;
t984 = t1026 * t1197;
t895 = -t1076 * t942 - t984;
t1261 = t1078 * t1101 - t1081 * t895;
t907 = qJDD(4) * pkin(4) + t1083 * pkin(8) - t1080 * t1102 + t943;
t756 = -t1077 * t907 + t1080 * t762;
t761 = t1076 * t818 - t1079 * t817;
t1227 = pkin(3) * t761 - pkin(7) * t756;
t1260 = (pkin(2) + t1139) * t761 + t1227;
t1086 = -pkin(5) * t996 + 0.2e1 * qJD(6) * t1026 + t907;
t1085 = t1086 + t1288;
t768 = (-t1275 - t941) * pkin(5) + t1085;
t1150 = -t1079 * t768 + t1275 * t1289 - t1300;
t1259 = t1150 - t1310;
t1220 = t1076 * t1272;
t916 = (-qJD(5) + t1055) * t1026 + t1148;
t840 = t1079 * t916 + t1220;
t1151 = pkin(8) * t840 + t1290 + t762;
t788 = t1077 * t840 + t1283;
t1258 = -t1242 * t788 - t1151;
t1238 = t941 * pkin(5);
t784 = t1085 - t1238;
t1240 = pkin(5) * t784;
t1201 = qJD(6) * t1055;
t1043 = 0.2e1 * t1201;
t956 = pkin(5) * t1024 - qJ(6) * t1026;
t1130 = -pkin(5) * t1244 + t1018 * qJ(6) - t1024 * t956 + t818;
t779 = t1043 + t1130;
t780 = -t1018 * pkin(5) - qJ(6) * t1244 + t1026 * t956 + qJDD(6) + t817;
t740 = t1076 * t780 + t1079 * t779;
t1154 = pkin(8) * t740 + t1079 * t1240 + t1289 * t784;
t735 = t1077 * t740 + t1080 * t784;
t1257 = -t1242 * t735 - t1154;
t765 = pkin(5) * t1267 + t779;
t770 = qJ(6) * t1267 + t780;
t839 = -t1079 * t915 + t1220;
t1155 = pkin(8) * t839 + t1076 * t770 + t1079 * t765 + t1290;
t787 = t1077 * t839 + t1283;
t1256 = -t1242 * t787 - t1155;
t767 = t1086 - t1238 + 0.2e1 * t1288;
t1156 = pkin(4) * t1273 + pkin(5) * t1212 + t1076 * t767 - t1330;
t800 = t1080 * t1273 - t1328;
t1255 = -t1242 * t800 - t1156;
t1211 = t1079 * t1272;
t838 = t1076 * t916 - t1211;
t746 = -pkin(8) * t838 - t761;
t790 = t1080 * t840 - t1287;
t1184 = pkin(3) * t838 - pkin(7) * t790 - t1077 * t746;
t1254 = (pkin(2) + t1233) * t838 + t1184;
t898 = t1079 * t907;
t1182 = -pkin(4) * t1275 + t1300 + t898;
t1253 = -t1182 - t1310;
t897 = t1076 * t907;
t921 = (qJD(5) + t1055) * t1024 + t1117;
t1183 = -pkin(4) * t921 - t1330 + t897;
t811 = t1080 * t921 + t1328;
t1252 = -t1242 * t811 + t1183;
t1226 = pkin(4) * t907 + pkin(8) * t762;
t755 = t1077 * t762 + t1080 * t907;
t1251 = -t1242 * t755 - t1226;
t1175 = t1077 * t1199;
t1249 = t1080 * t1100 - t1175;
t1248 = t1078 * t895 + t1081 * t1101;
t1096 = (-t1024 * t1079 + t1026 * t1076) * t1055;
t1090 = t1080 * t1018 - t1077 * t1096;
t1131 = -t984 - t1173;
t1247 = -t1078 * t1131 - t1081 * t1090;
t1246 = -t1078 * t1090 + t1081 * t1131;
t1094 = qJD(1) * t1243 + t1134;
t1089 = t1094 + t1205;
t982 = qJDD(1) * t1237 + t1089;
t1241 = pkin(2) * t982;
t1239 = pkin(7) * t864;
t1235 = pkin(1) * t1084;
t1187 = qJDD(1) * t1081;
t1037 = -t1078 * t1084 + t1187;
t1232 = pkin(6) * t1037;
t1188 = qJDD(1) * t1078;
t1038 = t1081 * t1084 + t1188;
t1231 = pkin(6) * t1038;
t1229 = qJDD(1) * pkin(1);
t1228 = t1080 * t746 + t838 * t1234;
t963 = -t1074 + t982;
t1225 = -pkin(3) * t963 + pkin(7) * t865;
t1223 = qJ(6) * t1079;
t952 = t1077 * t963;
t1206 = t1080 * t963;
t1204 = qJD(1) * qJD(4);
t1192 = t1080 * t1041;
t1189 = t1071 + t1072;
t1039 = t1189 * t1084;
t1181 = pkin(3) * t1039 + t865;
t1031 = -0.2e1 * t1172 + t1185;
t981 = -t1051 * t1077 - t1192;
t1180 = -pkin(3) * t1031 + pkin(7) * t981 + t952;
t1029 = 0.2e1 * t1058 + t1186;
t1022 = t1077 * t1042;
t978 = t1050 * t1080 - t1022;
t1179 = -pkin(3) * t1029 + pkin(7) * t978 - t1206;
t758 = t761 * t1234;
t1177 = -pkin(7) * t755 + t758;
t1176 = -pkin(7) * t976 + t1206;
t1171 = qJ(3) - t1230;
t1138 = -pkin(5) * t780 + qJ(6) * t779;
t739 = t1076 * t779 - t1079 * t780;
t721 = -pkin(4) * t739 - t1138;
t726 = -pkin(8) * t739 + (-pkin(5) * t1076 + t1223) * t784;
t1170 = -t1077 * t721 + t1080 * t726;
t836 = -t1076 * t915 - t1211;
t733 = -pkin(8) * t836 - t1076 * t765 + t1079 * t770;
t1137 = -pkin(5) * t1272 - qJ(6) * t915;
t772 = -pkin(4) * t836 - t1137;
t1169 = -t1077 * t772 + t1080 * t733;
t745 = -pkin(5) * t1221 + t1079 * t767 - t1331;
t1091 = pkin(5) * t946 + qJ(6) * t1265 + t1130;
t751 = -t1091 - 0.2e1 * t1201 - t1332;
t1168 = -t1077 * t751 + t1080 * t745;
t749 = -t1076 * t768 - t1223 * t1275 - t1299;
t1087 = pkin(5) * t1266 + qJ(6) * t1264 - t780;
t753 = -t1087 - t1301;
t1167 = -t1077 * t753 + t1080 * t749;
t776 = t817 - t1301;
t799 = -t897 - t1299;
t1166 = -t1077 * t776 + t1080 * t799;
t778 = t818 + t1332;
t810 = -t898 + t1331;
t1165 = -t1077 * t778 + t1080 * t810;
t1092 = -t1070 - t1110;
t983 = -qJDD(3) + t1092 + t1269;
t1163 = -t1078 * t982 - t1081 * t983;
t997 = t1092 + t1235;
t999 = t1095 + t1229;
t1162 = -t1078 * t999 - t1081 * t997;
t736 = -t1077 * t784 + t1080 * t740;
t1160 = -pkin(3) * t739 + pkin(7) * t736 + t1077 * t726 + t1080 * t721;
t789 = t1080 * t839 - t1287;
t1159 = -pkin(3) * t836 + pkin(7) * t789 + t1077 * t733 + t1080 * t772;
t802 = -t1077 * t1273 - t1324;
t1158 = pkin(7) * t802 + t1077 * t745 + t1080 * t751 - t1333;
t1157 = t1077 * t749 + t1080 * t753 + t1309;
t1153 = t1077 * t799 + t1080 * t776 + t1309;
t813 = -t1077 * t921 + t1324;
t1152 = pkin(7) * t813 + t1077 * t810 + t1080 * t778 + t1333;
t1149 = -t1045 * t1078 - t1081 * t1046;
t1147 = -t1083 + t1196;
t1146 = -t1083 + t1195;
t1144 = pkin(2) * t963 - t1225;
t1143 = -pkin(7) * t788 + t1228;
t1142 = t1078 * t1053;
t1141 = t1081 * t1053;
t1140 = pkin(2) * t983 - qJ(3) * g(3);
t1136 = g(3) * t1078 + t1232;
t1135 = g(3) * t1081 - t1231;
t1133 = -pkin(7) * t974 + t952;
t1132 = t1080 * t896 + t1175;
t1123 = t1078 * t983 - t1081 * t982;
t1122 = t1078 * t997 - t1081 * t999;
t1120 = pkin(2) * t1029 - t1179;
t1119 = pkin(2) * t1031 - t1180;
t1118 = t1045 * t1081 - t1046 * t1078;
t1116 = -pkin(7) * t735 + t1170;
t1115 = -pkin(7) * t787 + t1169;
t1114 = -pkin(7) * t800 + t1168;
t1113 = t1167 - t1316;
t1112 = t1166 - t1316;
t1111 = -pkin(7) * t811 + t1165;
t1036 = t1189 * qJDD(1);
t1109 = pkin(7) * t1036 - t864;
t1108 = pkin(2) * t739 - t1160;
t1107 = pkin(2) * t836 - t1159;
t1106 = -t1158 + t1334;
t1105 = -t1157 + t1303;
t1104 = -t1153 + t1303;
t1103 = -t1152 - t1334;
t1064 = pkin(2) * t1084 - g(3);
t1098 = -pkin(2) * t1188 - t1081 * t1064 - t1231;
t1088 = t1077 * t1018 + t1080 * t1096;
t1082 = pkin(1) * g(3);
t1075 = qJ(2) * g(3);
t1073 = pkin(2) * qJDD(1);
t1066 = 0.2e1 * t1070;
t1065 = 0.2e1 * qJDD(1) * qJ(3);
t1040 = (-t1071 + t1072) * t1084;
t1028 = pkin(3) * t1036;
t1021 = t1189 * t1204;
t1014 = -t1134 - 0.2e1 * t1229;
t1006 = t1066 + t1110;
t998 = t1066 + t1097;
t992 = t1065 + t1094 + 0.2e1 * t1229;
t987 = t1030 * t1077 + t1072 * t1204;
t986 = -t1071 * t1204 + t1080 * t1145;
t980 = t1077 * t1146 + t1191;
t979 = (t1030 - t1172) * t1080;
t977 = t1080 * t1147 - t1193;
t975 = t1080 * t1146 - t1022;
t973 = -t1077 * t1147 - t1192;
t972 = (t1145 + t1058) * t1077;
t971 = pkin(2) * t1036 + t1028;
t965 = -pkin(2) * t1187 + t1064 * t1078 - t1232;
t957 = t1075 - t1241;
t951 = -t1029 * t1080 - t1031 * t1077;
t950 = t1077 * t1029 - t1080 * t1031;
t947 = t1082 - t1140;
t930 = pkin(1) * t999 - qJ(2) * t997;
t863 = -qJ(2) * t983 + t1237 * t982;
t845 = -pkin(2) * t1039 - t1181;
t844 = -qJ(2) * t981 - t1119;
t843 = -qJ(2) * t978 - t1120;
t842 = t1031 * t1237 + t1236 * t976 + t1206;
t841 = t1029 * t1237 + t1236 * t974 + t952;
t827 = -t1036 * t1236 - t1039 * t1237 - t864;
t809 = -t1237 * t981 - t1262;
t808 = -t1237 * t978 - t1263;
t775 = -qJ(2) * t865 - t1144;
t774 = t1236 * t864 + t1237 * t963;
t757 = -t1237 * t865 + t1274;
t730 = -t1237 * t813 - t1252;
t729 = -t1253 - t1311;
t728 = -qJ(2) * t813 - t1103;
t727 = t1236 * t811 + t1165 - t1320;
t723 = -t1104 - t1315;
t722 = t1166 + t1306;
t719 = -t1259 - t1311;
t718 = -t1237 * t802 - t1255;
t717 = -qJ(2) * t790 - t1254;
t716 = t1236 * t788 + t1237 * t838 + t1228;
t715 = -t1237 * t790 - t1258;
t714 = -t1105 - t1315;
t713 = t1167 + t1306;
t712 = -qJ(2) * t802 - t1106;
t711 = t1236 * t800 + t1168 + t1320;
t710 = -t1237 * t789 - t1256;
t709 = -qJ(2) * t789 - t1107;
t708 = t1236 * t787 + t1237 * t836 + t1169;
t707 = -t1237 * t756 - t1251;
t706 = -qJ(2) * t756 - t1260;
t705 = t758 + t1236 * t755 + (pkin(1) + t1171) * t761;
t704 = -t1237 * t736 - t1257;
t703 = -qJ(2) * t736 - t1108;
t702 = t1236 * t735 + t1237 * t739 + t1170;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1037, 0, -t1038, 0, -t1136, -t1135, -t1118, -pkin(6) * t1118, 0, -t1037, t1038, 0, 0, 0, t1122, t1136, t1135, pkin(6) * t1122 + (-pkin(1) * t1078 + qJ(2) * t1081) * g(3), 0, t1038, t1037, 0, 0, 0, t1123, t1098, t965, pkin(6) * t1123 - t1078 * t947 + t1081 * t957, t1081 * t987 - t1142, -t1040 * t1078 - t1081 * t950, -t1078 * t1185 - t1081 * t975, -t1081 * t986 + t1142, t1078 * t1186 - t1081 * t973, -qJDD(4) * t1078 - t1021 * t1081, t1081 * t843 - t1078 * t808 - pkin(6) * (t1029 * t1081 + t1078 * t974), t1081 * t844 - t1078 * t809 - pkin(6) * (t1031 * t1081 + t1078 * t976), -t1081 * t845 + t1078 * t971 - pkin(6) * (-t1036 * t1078 - t1039 * t1081), t1081 * t775 - t1078 * t757 - pkin(6) * (t1078 * t864 + t1081 * t963), t1248, -t1335, t1318, t1276, t1336, t1247, -t1078 * t729 + t1081 * t723 - t1305, t1081 * t728 - t1078 * t730 - pkin(6) * (t1078 * t811 - t1322), t1081 * t717 - t1078 * t715 - pkin(6) * (t1078 * t788 + t1081 * t838), t1081 * t706 - t1078 * t707 - pkin(6) * (t1078 * t755 + t1081 * t761), t1248, t1318, t1335, t1247, -t1336, t1276, -t1078 * t719 + t1081 * t714 - t1305, t1081 * t709 - t1078 * t710 - pkin(6) * (t1078 * t787 + t1081 * t836), t1081 * t712 - t1078 * t718 - pkin(6) * (t1078 * t800 + t1322), t1081 * t703 - t1078 * t704 - pkin(6) * (t1078 * t735 + t1081 * t739); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1038, 0, t1037, 0, t1135, -t1136, t1149, pkin(6) * t1149, 0, -t1038, -t1037, 0, 0, 0, t1162, -t1135, t1136, pkin(6) * t1162 + (pkin(1) * t1081 + qJ(2) * t1078) * g(3), 0, -t1037, t1038, 0, 0, 0, t1163, -t965, t1098, pkin(6) * t1163 + t1078 * t957 + t1081 * t947, t1078 * t987 + t1141, t1040 * t1081 - t1078 * t950, -t1078 * t975 + t1081 * t1185, -t1078 * t986 - t1141, -t1078 * t973 - t1081 * t1186, qJDD(4) * t1081 - t1021 * t1078, t1078 * t843 + t1081 * t808 + pkin(6) * (-t1029 * t1078 + t1081 * t974), t1078 * t844 + t1081 * t809 + pkin(6) * (-t1031 * t1078 + t1081 * t976), -t1078 * t845 - t1081 * t971 + pkin(6) * (-t1036 * t1081 + t1039 * t1078), t1078 * t775 + t1081 * t757 + pkin(6) * (-t1078 * t963 + t1081 * t864), t1261, -t1337, t1319, t1278, t1338, t1246, t1078 * t723 + t1081 * t729 + t1304, t1078 * t728 + t1081 * t730 + pkin(6) * (t1081 * t811 + t1326), t1078 * t717 + t1081 * t715 + pkin(6) * (-t1078 * t838 + t1081 * t788), t1078 * t706 + t1081 * t707 + pkin(6) * (-t1078 * t761 + t1081 * t755), t1261, t1319, t1337, t1246, -t1338, t1278, t1078 * t714 + t1081 * t719 + t1304, t1078 * t709 + t1081 * t710 + pkin(6) * (-t1078 * t836 + t1081 * t787), t1078 * t712 + t1081 * t718 + pkin(6) * (t1081 * t800 - t1326), t1078 * t703 + t1081 * t704 + pkin(6) * (-t1078 * t739 + t1081 * t735); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1045, t1046, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t1014, t1006, t930, qJDD(1), 0, 0, 0, 0, 0, 0, t998, t992, t863, t979, t951, t980, t972, t977, 0, t841, t842, t827, t774, t1132, -t793, t1307, t1249, t825, t1088, t722, t727, t716, t705, t1132, t1307, t793, t1088, -t825, t1249, t713, t708, t711, t702; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1084, 0, 0, -g(3), -t1045, 0, 0, -qJDD(1), t1084, 0, 0, 0, -t999, 0, g(3), t1075, 0, t1084, qJDD(1), 0, 0, 0, -t982, -t1064, -t1073, t957, t987, -t950, -t975, -t986, -t973, -t1021, t843, t844, -t845, t775, t1101, -t792, -t1308, t1250, t822, -t1090, t723, t728, t717, t706, t1101, -t1308, t792, -t1090, -t822, t1250, t714, t709, t712, t703; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1084, 0, qJDD(1), 0, g(3), 0, -t1046, 0, 0, -t1084, -qJDD(1), 0, 0, 0, -t997, -g(3), 0, t1082, 0, -qJDD(1), t1084, 0, 0, 0, -t983, t1073, -t1064, t947, t1053, t1040, t1185, -t1053, -t1186, qJDD(4), t808, t809, -t971, t757, -t895, t837, t1294, t1099, t882, t1131, t729, t730, t715, t707, -t895, t1294, -t837, t1131, -t882, t1099, t719, t710, t718, t704; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1045, t1046, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t1014, t1006, t930, qJDD(1), 0, 0, 0, 0, 0, 0, t998, t992, t863, t979, t951, t980, t972, t977, 0, t841, t842, t827, t774, t1132, -t793, t1307, t1249, t825, t1088, t722, t727, t716, t705, t1132, t1307, t793, t1088, -t825, t1249, t713, t708, t711, t702; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t999, -t997, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t1093 - t1235, t1065 + t1089 + t1229, qJ(3) * t982, t979, t951, t980, t972, t977, 0, qJ(3) * t1029 + t1133, qJ(3) * t1031 + t1176, -qJ(3) * t1039 + t1109, qJ(3) * t963 - t1239, t1132, -t793, t1307, t1249, t825, t1088, t1112 + t1298, t1111 - t1329, qJ(3) * t838 + t1143, t1171 * t761 + t1177, t1132, t1307, t793, t1088, -t825, t1249, t1113 + t1298, qJ(3) * t836 + t1115, t1114 + t1329, qJ(3) * t739 + t1116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1084, 0, 0, 0, t999, 0, -g(3), 0, 0, -t1084, -qJDD(1), 0, 0, 0, t982, t1064, t1073, t1241, -t987, t950, t975, t986, t973, t1021, t1120, t1119, t845, t1144, -t1101, t792, t1308, -t1250, -t822, t1090, t1104, t1103, t1254, t1260, -t1101, t1308, -t792, t1090, t822, -t1250, t1105, t1107, t1106, t1108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1084, qJDD(1), 0, 0, 0, t997, g(3), 0, 0, 0, qJDD(1), -t1084, 0, 0, 0, t983, -t1073, t1064, t1140, -t1053, -t1040, -t1185, t1053, t1186, -qJDD(4), qJ(3) * t978 + t1263, qJ(3) * t981 + t1262, t971, qJ(3) * t865 - t1274, t895, -t837, -t1294, -t1099, -t882, -t1131, t1253 + t1314, qJ(3) * t813 + t1252, qJ(3) * t790 + t1258, qJ(3) * t756 + t1251, t895, -t1294, t837, -t1131, t882, -t1099, t1259 + t1314, qJ(3) * t789 + t1256, qJ(3) * t802 + t1255, qJ(3) * t736 + t1257; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t983, t982, 0, t979, t951, t980, t972, t977, 0, t1133, t1176, t1109, -t1239, t1132, -t793, t1307, t1249, t825, t1088, t1112, t1111, t1143, -t1230 * t761 + t1177, t1132, t1307, t793, t1088, -t825, t1249, t1113, t1115, t1114, t1116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1084, 0, 0, 0, t983, 0, -g(3), 0, -t1053, -t1040, -t1185, t1053, t1186, -qJDD(4), -pkin(3) * t974 - t943, -pkin(3) * t976 - t944, t1028, -pkin(3) * t864, t895, -t837, -t1294, -t1099, -t882, -t1131, -t1182 - t1317, -pkin(3) * t811 + t1183, -pkin(3) * t788 - t1151, -pkin(3) * t755 - t1226, t895, -t1294, t837, -t1131, t882, -t1099, t1150 - t1317, -pkin(3) * t787 - t1155, -pkin(3) * t800 - t1156, -pkin(3) * t735 - t1154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1084, qJDD(1), 0, 0, 0, -t982, g(3), 0, 0, t987, -t950, -t975, -t986, -t973, -t1021, t1179, t1180, t1181, t1225, t1101, -t792, -t1308, t1250, t822, -t1090, t1153, t1152, -t1233 * t838 - t1184, -t1139 * t761 - t1227, t1101, -t1308, t792, -t1090, -t822, t1250, t1157, t1159, t1158, t1160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1030, -t1029, t1042, t1172, t1147, -t1172, 0, t963, -t943, 0, t896, -t1129, t1293, t1100, t1126, t1096, t799, t810, t746, -pkin(8) * t761, t896, t1293, t1129, t1096, -t1126, t1100, t749, t733, t745, t726; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1058, t1031, -t1146, -t1145, t1041, -t1058, -t963, 0, -t944, 0, -t1199, -t1268, -t1272, t1199, t915, -t1018, t776, t778, -pkin(4) * t838, -pkin(4) * t761, -t1199, -t1272, t1268, -t1018, -t915, t1199, t753, t772, t751, t721; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1053, t1040, t1185, -t1053, -t1186, qJDD(4), t943, t944, 0, 0, -t895, t837, t1294, t1099, t882, t1131, t1182, -t1183, t1151, t1226, -t895, t1294, -t837, t1131, -t882, t1099, -t1150, t1155, t1156, t1154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t942, -t1275, t1266, t1200, t988, -t1200, 0, -t907, t817, 0, t942, t1266, t1275, -t1200, -t988, t1200, -qJ(6) * t1275, t770, t767, qJ(6) * t784; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t996, t1273, -t989, -t941, t1265, -t996, t907, 0, t818, 0, t996, -t989, -t1273, -t996, -t1265, -t941, t768, t765, pkin(5) * t1273, t1240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1199, t1268, t1272, -t1199, -t915, t1018, -t817, -t818, 0, 0, t1199, t1272, -t1268, t1018, t915, -t1199, t1087, t1137, t1043 + t1091, t1138; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t942, t1266, t1275, -t1200, -t988, t1200, 0, t780, t784, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1199, t1272, -t1268, t1018, t915, -t1199, -t780, 0, t779, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t996, t989, t1273, t996, t1265, t941, -t784, -t779, 0, 0;];
m_new_reg  = t1;
