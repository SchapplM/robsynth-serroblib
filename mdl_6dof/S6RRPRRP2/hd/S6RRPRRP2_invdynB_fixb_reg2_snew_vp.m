% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6RRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 17:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6RRPRRP2_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRP2_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP2_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_invdynB_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:25:44
% EndTime: 2019-05-06 17:26:37
% DurationCPUTime: 52.39s
% Computational Cost: add. (187738->850), mult. (437921->1222), div. (0->0), fcn. (329476->10), ass. (0->611)
t1133 = sin(qJ(1));
t1137 = cos(qJ(1));
t1132 = sin(qJ(2));
t1136 = cos(qJ(2));
t1128 = sin(pkin(10));
t1129 = cos(pkin(10));
t1131 = sin(qJ(4));
t1135 = cos(qJ(4));
t1134 = cos(qJ(5));
t1130 = sin(qJ(5));
t1225 = qJD(1) * t1136;
t1226 = qJD(1) * t1132;
t1088 = -t1128 * t1226 + t1129 * t1225;
t1089 = (t1128 * t1136 + t1129 * t1132) * qJD(1);
t1053 = t1088 * t1131 + t1089 * t1135;
t1200 = qJD(2) + qJD(4);
t1022 = t1130 * t1053 - t1134 * t1200;
t1192 = qJDD(2) + qJDD(4);
t1051 = -t1135 * t1088 + t1089 * t1131;
t1179 = qJD(2) * t1225;
t1196 = qJDD(1) * t1132;
t1098 = t1179 + t1196;
t1122 = t1136 * qJDD(1);
t1188 = qJD(2) * t1226;
t1099 = t1122 - t1188;
t1063 = t1129 * t1098 + t1128 * t1099;
t1172 = t1098 * t1128 - t1129 * t1099;
t968 = -t1051 * qJD(4) + t1135 * t1063 - t1131 * t1172;
t1151 = -t1130 * t1192 - t1134 * t968;
t1148 = -t1022 * qJD(5) - t1151;
t1048 = qJD(5) + t1051;
t1219 = t1022 * t1048;
t1269 = -t1219 + t1148;
t1249 = t1130 * t1269;
t1024 = t1134 * t1053 + t1130 * t1200;
t1176 = t1130 * t968 - t1134 * t1192;
t1198 = qJD(5) + t1048;
t887 = t1024 * t1198 + t1176;
t807 = t1134 * t887 + t1249;
t1020 = t1024 ^ 2;
t1260 = t1022 ^ 2;
t974 = t1020 - t1260;
t763 = t1131 * t807 + t1135 * t974;
t765 = -t1131 * t974 + t1135 * t807;
t702 = t1128 * t765 + t1129 * t763;
t705 = t1128 * t763 - t1129 * t765;
t657 = t1132 * t702 + t1136 * t705;
t803 = -t1130 * t887 + t1134 * t1269;
t1379 = t1133 * t657 - t1137 * t803;
t1173 = t1131 * t1063 + t1135 * t1172;
t967 = -t1053 * qJD(4) - t1173;
t966 = qJDD(5) - t967;
t977 = t1024 * t1022;
t1278 = t966 + t977;
t1247 = t1130 * t1278;
t1259 = t1048 ^ 2;
t984 = t1260 - t1259;
t862 = t1134 * t984 - t1247;
t888 = (qJD(5) - t1048) * t1024 + t1176;
t788 = t1131 * t862 + t1135 * t888;
t792 = -t1131 * t888 + t1135 * t862;
t724 = t1128 * t792 + t1129 * t788;
t727 = t1128 * t788 - t1129 * t792;
t669 = t1132 * t724 + t1136 * t727;
t1235 = t1134 * t1278;
t857 = t1130 * t984 + t1235;
t1378 = t1133 * t669 + t1137 * t857;
t1377 = t1133 * t803 + t1137 * t657;
t1376 = -t1133 * t857 + t1137 * t669;
t1375 = t1132 * t705 - t1136 * t702;
t1374 = t1132 * t727 - t1136 * t724;
t1270 = -t1219 - t1148;
t1306 = -t1130 * t888 + t1134 * t1270;
t1271 = t1020 + t1260;
t1304 = -t1130 * t1270 - t1134 * t888;
t1324 = -t1131 * t1271 + t1135 * t1304;
t1328 = t1131 * t1304 + t1135 * t1271;
t1339 = t1128 * t1324 + t1129 * t1328;
t1340 = -t1128 * t1328 + t1129 * t1324;
t1360 = -t1132 * t1339 + t1136 * t1340;
t1365 = t1133 * t1306 + t1137 * t1360;
t1373 = pkin(6) * t1365;
t1367 = t1133 * t1360 - t1137 * t1306;
t1372 = pkin(6) * t1367;
t1361 = t1132 * t1340 + t1136 * t1339;
t1371 = pkin(7) * t1361;
t1370 = -pkin(1) * t1361 - pkin(2) * t1339 - pkin(3) * t1328 - pkin(4) * t1271 - pkin(9) * t1304;
t1369 = -pkin(1) * t1306 + pkin(7) * t1360;
t1279 = t966 - t977;
t1246 = t1130 * t1279;
t985 = -t1020 + t1259;
t1308 = -t1134 * t985 - t1246;
t1234 = t1134 * t1279;
t1307 = -t1130 * t985 + t1234;
t1323 = -t1131 * t1270 + t1135 * t1307;
t1327 = t1131 * t1307 + t1135 * t1270;
t1341 = t1128 * t1323 + t1129 * t1327;
t1342 = -t1128 * t1327 + t1129 * t1323;
t1358 = -t1132 * t1341 + t1136 * t1342;
t1368 = t1133 * t1358 + t1137 * t1308;
t1366 = -t1133 * t1308 + t1137 * t1358;
t1364 = qJ(3) * t1339;
t1362 = -pkin(2) * t1306 + qJ(3) * t1340;
t1359 = t1132 * t1342 + t1136 * t1341;
t942 = t1259 + t1020;
t833 = t1134 * t942 + t1247;
t1357 = pkin(1) * t833;
t1356 = pkin(2) * t833;
t1355 = pkin(3) * t833;
t1354 = pkin(4) * t833;
t1353 = pkin(9) * t833;
t841 = t1130 * t942 - t1235;
t1352 = pkin(9) * t841;
t1351 = pkin(8) * t1324;
t1350 = pkin(8) * t1328;
t1349 = t1131 * t841;
t1348 = t1133 * t833;
t1346 = t1135 * t841;
t1345 = t1137 * t833;
t1336 = pkin(9) * t1306;
t1217 = t1048 * t1130;
t1185 = t1022 * t1217;
t1216 = t1048 * t1134;
t979 = t1024 * t1216;
t1164 = t979 + t1185;
t1184 = t1022 * t1216;
t978 = t1024 * t1217;
t1163 = t978 - t1184;
t1263 = t1131 * t966 + t1135 * t1163;
t1266 = t1131 * t1163 - t1135 * t966;
t1284 = t1128 * t1263 + t1129 * t1266;
t1285 = -t1128 * t1266 + t1129 * t1263;
t1299 = -t1132 * t1284 + t1136 * t1285;
t1326 = t1133 * t1299 + t1137 * t1164;
t921 = -qJD(5) * t1024 - t1176;
t1152 = -t1134 * t921 - t1185;
t1153 = -t1130 * t921 + t1184;
t1187 = t1131 * t977;
t1264 = t1135 * t1153 - t1187;
t1186 = t1135 * t977;
t1265 = t1131 * t1153 + t1186;
t1282 = t1128 * t1264 + t1129 * t1265;
t1283 = -t1128 * t1265 + t1129 * t1264;
t1301 = -t1132 * t1282 + t1136 * t1283;
t1325 = t1133 * t1301 + t1137 * t1152;
t1322 = -t1133 * t1164 + t1137 * t1299;
t1321 = -t1133 * t1152 + t1137 * t1301;
t1268 = -t1259 - t1260;
t1287 = t1130 * t1268 + t1234;
t1320 = pkin(1) * t1287;
t1319 = pkin(2) * t1287;
t1318 = pkin(3) * t1287;
t1317 = pkin(4) * t1287;
t1286 = t1134 * t1268 - t1246;
t1316 = pkin(9) * t1286;
t1315 = pkin(9) * t1287;
t1312 = t1131 * t1286;
t1311 = t1133 * t1287;
t1310 = t1135 * t1286;
t1309 = t1137 * t1287;
t1127 = t1136 ^ 2;
t1139 = qJD(1) ^ 2;
t1155 = qJD(2) * pkin(2) - qJ(3) * t1226;
t1108 = t1133 * g(1) - t1137 * g(2);
t1156 = qJDD(1) * pkin(1) + t1108;
t1021 = t1099 * pkin(2) - t1155 * t1226 - qJDD(3) + t1156 + (qJ(3) * t1127 + pkin(7)) * t1139;
t886 = t1134 * t1148 - t978;
t1165 = t1131 * t886 - t1186;
t1166 = t1135 * t886 + t1187;
t1261 = -t1128 * t1165 + t1129 * t1166;
t1262 = t1128 * t1166 + t1129 * t1165;
t1281 = -t1132 * t1262 + t1136 * t1261;
t885 = -t1130 * t1148 - t979;
t1305 = t1133 * t1281 + t1137 * t885;
t1303 = -t1133 * t885 + t1137 * t1281;
t1302 = t1132 * t1283 + t1136 * t1282;
t1300 = t1132 * t1285 + t1136 * t1284;
t1298 = 2 * qJD(6);
t1296 = qJ(6) * t1269;
t1061 = t1088 * t1089;
t1267 = qJDD(2) + t1061;
t1295 = t1128 * t1267;
t1294 = t1129 * t1267;
t997 = t1053 * t1051;
t1276 = -t997 + t1192;
t1291 = t1131 * t1276;
t1288 = t1135 * t1276;
t1086 = t1088 ^ 2;
t1168 = qJD(2) * pkin(3) - pkin(8) * t1089;
t946 = -pkin(3) * t1172 + t1086 * pkin(8) - t1089 * t1168 + t1021;
t1280 = t1132 * t1261 + t1136 * t1262;
t1191 = t1200 ^ 2;
t1040 = t1200 * t1051;
t1277 = -t968 + t1040;
t1124 = t1127 * t1139;
t1138 = qJD(2) ^ 2;
t1114 = -t1124 - t1138;
t1224 = qJD(2) * t1088;
t1029 = -t1063 + t1224;
t1109 = g(1) * t1137 + g(2) * t1133;
t1092 = -pkin(1) * t1139 + qJDD(1) * pkin(7) - t1109;
t1255 = t1132 * g(3);
t1142 = -pkin(2) * t1124 + t1099 * qJ(3) - qJD(2) * t1155 - t1255;
t1202 = t1132 * t1139;
t1227 = qJD(1) * qJD(2);
t1149 = pkin(2) * t1202 + qJ(3) * t1227 - g(3);
t1203 = t1132 * t1092;
t1150 = qJDD(2) * pkin(2) - t1098 * qJ(3) - t1203;
t1258 = 2 * qJD(3);
t1190 = t1089 * t1258;
t947 = t1128 * t1142 - t1129 * t1150 + t1136 * (t1128 * t1092 - t1129 * t1149) + t1190;
t1140 = pkin(3) * t1267 + pkin(8) * t1029 - t947;
t1201 = t1136 * t1092;
t948 = t1129 * (t1142 + t1201) + t1128 * (t1136 * t1149 + t1150) + t1088 * t1258;
t909 = -t1086 * pkin(3) - pkin(8) * t1172 - qJD(2) * t1168 + t948;
t818 = t1131 * t1140 + t1135 * t909;
t994 = pkin(4) * t1051 - pkin(9) * t1053;
t798 = -pkin(4) * t1191 + pkin(9) * t1192 - t1051 * t994 + t818;
t1175 = t1200 * t1053;
t815 = t1277 * pkin(9) + (-t967 + t1175) * pkin(4) - t946;
t735 = t1130 * t815 + t1134 * t798;
t973 = pkin(5) * t1022 - qJ(6) * t1024;
t1157 = t966 * qJ(6) - t1022 * t973 + t1048 * t1298 + t735;
t1049 = t1051 ^ 2;
t1050 = t1053 ^ 2;
t1087 = t1089 ^ 2;
t1257 = pkin(4) * t1131;
t1256 = pkin(5) * t1134;
t734 = t1130 * t798 - t1134 * t815;
t1254 = qJ(6) * t1134;
t817 = t1131 * t909 - t1135 * t1140;
t740 = t1131 * t818 - t1135 * t817;
t1253 = t1128 * t740;
t1252 = t1129 * t740;
t797 = -t1192 * pkin(4) - t1191 * pkin(9) + t1053 * t994 + t817;
t1251 = t1130 * t797;
t1243 = t1131 * t946;
t992 = t997 + t1192;
t1241 = t1131 * t992;
t870 = t1128 * t948 - t1129 * t947;
t1240 = t1132 * t870;
t1238 = t1134 * t797;
t1231 = t1135 * t946;
t1230 = t1135 * t992;
t1229 = t1136 * t870;
t1228 = -t1259 + t1271;
t1223 = qJD(2) * t1089;
t1221 = t1021 * t1128;
t1220 = t1021 * t1129;
t1218 = t1048 * t1024;
t1057 = qJDD(2) - t1061;
t1215 = t1057 * t1128;
t1214 = t1057 * t1129;
t1213 = t1088 * t1128;
t1212 = t1088 * t1129;
t1211 = t1089 * t1128;
t1210 = t1089 * t1129;
t1091 = t1139 * pkin(7) + t1156;
t1209 = t1091 * t1132;
t1208 = t1091 * t1136;
t1115 = t1136 * t1202;
t1106 = qJDD(2) + t1115;
t1207 = t1106 * t1132;
t1107 = qJDD(2) - t1115;
t1206 = t1107 * t1132;
t1205 = t1107 * t1136;
t1126 = t1132 ^ 2;
t1204 = t1126 * t1139;
t1197 = t1126 + t1127;
t1195 = qJDD(1) * t1133;
t1194 = qJDD(1) * t1137;
t1193 = qJDD(2) * t1137;
t1189 = -pkin(4) * t1135 - pkin(3);
t1183 = t1133 * t997;
t1182 = t1137 * t997;
t1181 = t1133 * t1061;
t1180 = t1137 * t1061;
t1178 = qJ(6) * t1130 + pkin(4);
t871 = t1128 * t947 + t1129 * t948;
t741 = t1131 * t817 + t1135 * t818;
t1076 = t1136 * g(3) + t1203;
t1077 = t1201 - t1255;
t1016 = t1076 * t1132 + t1136 * t1077;
t1069 = -t1108 * t1133 - t1137 * t1109;
t1171 = t1024 * t973 + qJDD(6) + t734;
t1170 = t1133 * t1115;
t1169 = t1137 * t1115;
t1103 = -t1133 * t1139 + t1194;
t1167 = -pkin(6) * t1103 - g(3) * t1133;
t1161 = t1131 * t1040;
t1160 = t1131 * t1175;
t1159 = t1135 * t1040;
t1158 = t1135 * t1175;
t673 = t1130 * t735 - t1134 * t734;
t674 = t1130 * t734 + t1134 * t735;
t1015 = t1076 * t1136 - t1077 * t1132;
t1068 = t1108 * t1137 - t1109 * t1133;
t1154 = -pkin(5) * t966 + t1171;
t934 = qJD(2) * t1053 - t1173;
t1027 = -t1172 + t1223;
t1146 = -t921 * pkin(5) - t1296 + t797;
t1145 = t1024 * t1298 - t1146;
t1121 = t1133 * qJDD(2);
t1113 = t1124 - t1138;
t1112 = -t1138 - t1204;
t1111 = t1138 - t1204;
t1105 = t1124 - t1204;
t1104 = t1124 + t1204;
t1102 = t1137 * t1139 + t1195;
t1101 = t1197 * qJDD(1);
t1100 = t1122 - 0.2e1 * t1188;
t1097 = 0.2e1 * t1179 + t1196;
t1095 = t1136 * t1106;
t1094 = t1197 * t1227;
t1085 = -pkin(6) * t1102 + g(3) * t1137;
t1082 = -t1087 - t1138;
t1081 = -t1087 + t1138;
t1080 = t1086 - t1138;
t1079 = t1098 * t1136 - t1126 * t1227;
t1078 = -t1099 * t1132 - t1127 * t1227;
t1075 = -t1112 * t1132 - t1205;
t1074 = -t1111 * t1132 + t1095;
t1073 = t1114 * t1136 - t1207;
t1072 = t1113 * t1136 - t1206;
t1071 = t1112 * t1136 - t1206;
t1070 = t1114 * t1132 + t1095;
t1066 = t1101 * t1137 - t1104 * t1133;
t1065 = t1101 * t1133 + t1104 * t1137;
t1064 = -t1097 * t1132 + t1100 * t1136;
t1060 = -t1087 + t1086;
t1055 = -t1138 - t1086;
t1046 = (t1211 + t1212) * qJD(2);
t1045 = (-t1210 + t1213) * qJD(2);
t1044 = t1075 * t1137 + t1097 * t1133;
t1043 = t1073 * t1137 - t1100 * t1133;
t1042 = t1075 * t1133 - t1097 * t1137;
t1041 = t1073 * t1133 + t1100 * t1137;
t1037 = -t1050 + t1191;
t1036 = t1049 - t1191;
t1035 = -pkin(7) * t1071 - t1208;
t1034 = -pkin(7) * t1070 - t1209;
t1033 = -t1050 - t1191;
t1032 = -pkin(1) * t1071 + t1077;
t1031 = -pkin(1) * t1070 + t1076;
t1028 = t1063 + t1224;
t1026 = t1172 + t1223;
t1025 = -t1086 - t1087;
t1013 = -qJD(2) * t1211 + t1063 * t1129;
t1012 = qJD(2) * t1210 + t1063 * t1128;
t1011 = -qJD(2) * t1212 + t1128 * t1172;
t1010 = -qJD(2) * t1213 - t1129 * t1172;
t1007 = -t1082 * t1128 - t1214;
t1006 = -t1081 * t1128 + t1294;
t1005 = t1080 * t1129 - t1215;
t1004 = t1082 * t1129 - t1215;
t1003 = t1081 * t1129 + t1295;
t1002 = t1080 * t1128 + t1214;
t1001 = t1016 * t1137 - t1091 * t1133;
t1000 = t1016 * t1133 + t1091 * t1137;
t999 = t1055 * t1129 - t1295;
t998 = t1055 * t1128 + t1294;
t996 = t1049 - t1050;
t990 = -t1191 - t1049;
t982 = -t1045 * t1132 + t1046 * t1136;
t981 = -t1159 + t1160;
t980 = -t1161 - t1158;
t972 = t1027 * t1129 - t1029 * t1128;
t971 = -t1026 * t1129 - t1028 * t1128;
t970 = t1027 * t1128 + t1029 * t1129;
t969 = -t1026 * t1128 + t1028 * t1129;
t963 = -qJ(3) * t1004 - t1220;
t962 = -t1049 - t1050;
t961 = -t1012 * t1132 + t1013 * t1136;
t960 = -t1010 * t1132 + t1011 * t1136;
t959 = -t1004 * t1132 + t1007 * t1136;
t958 = -t1003 * t1132 + t1006 * t1136;
t957 = -t1002 * t1132 + t1005 * t1136;
t956 = t1004 * t1136 + t1007 * t1132;
t953 = -qJ(3) * t998 - t1221;
t952 = t1036 * t1135 - t1241;
t951 = -t1037 * t1131 + t1288;
t950 = t1036 * t1131 + t1230;
t949 = t1037 * t1135 + t1291;
t945 = -t1033 * t1131 - t1230;
t944 = t1033 * t1135 - t1241;
t937 = -t1040 - t968;
t933 = (0.2e1 * qJD(4) + qJD(2)) * t1053 + t1173;
t930 = t1135 * t968 - t1160;
t929 = t1131 * t968 + t1158;
t928 = -t1131 * t967 + t1159;
t927 = t1135 * t967 + t1161;
t926 = -t1132 * t998 + t1136 * t999;
t925 = t1132 * t999 + t1136 * t998;
t924 = -pkin(2) * t1028 + qJ(3) * t1007 - t1221;
t923 = t1135 * t990 - t1291;
t922 = t1131 * t990 + t1288;
t918 = -pkin(2) * t1026 + qJ(3) * t999 + t1220;
t917 = t1028 * t1133 + t1137 * t959;
t916 = -t1028 * t1137 + t1133 * t959;
t911 = -t1128 * t980 + t1129 * t981;
t910 = t1128 * t981 + t1129 * t980;
t902 = t1026 * t1133 + t1137 * t926;
t901 = -t1132 * t970 + t1136 * t972;
t900 = -t1132 * t969 + t1136 * t971;
t899 = -t1026 * t1137 + t1133 * t926;
t898 = t1132 * t972 + t1136 * t970;
t895 = t1022 * t1198 + t1151;
t889 = -t921 + t1218;
t877 = t1025 * t1133 + t1137 * t901;
t876 = -t1025 * t1137 + t1133 * t901;
t875 = -t1128 * t950 + t1129 * t952;
t874 = -t1128 * t949 + t1129 * t951;
t873 = t1128 * t952 + t1129 * t950;
t872 = t1128 * t951 + t1129 * t949;
t869 = -pkin(8) * t944 - t1231;
t868 = -t1128 * t944 + t1129 * t945;
t867 = t1128 * t945 + t1129 * t944;
t854 = -t1131 * t937 + t1135 * t934;
t853 = t1131 * t1277 - t1135 * t933;
t852 = t1131 * t934 + t1135 * t937;
t851 = -t1131 * t933 - t1135 * t1277;
t850 = -pkin(8) * t922 - t1243;
t849 = -t1128 * t929 + t1129 * t930;
t848 = -t1128 * t927 + t1129 * t928;
t847 = t1128 * t930 + t1129 * t929;
t846 = t1128 * t928 + t1129 * t927;
t845 = -pkin(1) * t898 - pkin(2) * t970;
t844 = -pkin(1) * t956 - pkin(2) * t1004 + t948;
t843 = pkin(2) * t1021 + qJ(3) * t871;
t838 = -t1128 * t922 + t1129 * t923;
t837 = t1128 * t923 + t1129 * t922;
t824 = -pkin(1) * t925 + t1128 * t1077 + t1129 * t1076 + t1190 + (t1128 * (t1099 + t1188) - t1129 * (-t1098 + t1179)) * qJ(3) + (-t1129 * t1106 + t1114 * t1128 - t998) * pkin(2);
t823 = -qJ(3) * t970 - t870;
t822 = -t1132 * t910 + t1136 * t911;
t821 = -pkin(7) * t956 - t1132 * t924 + t1136 * t963;
t820 = pkin(3) * t1277 + pkin(8) * t945 - t1243;
t819 = -pkin(2) * t1025 + qJ(3) * t972 + t871;
t812 = -pkin(3) * t933 + pkin(8) * t923 + t1231;
t811 = -pkin(7) * t925 - t1132 * t918 + t1136 * t953;
t802 = -t1132 * t873 + t1136 * t875;
t801 = -t1132 * t872 + t1136 * t874;
t800 = t1136 * t871 - t1240;
t799 = t1132 * t871 + t1229;
t794 = -t1132 * t867 + t1136 * t868;
t793 = t1132 * t868 + t1136 * t867;
t784 = -t1021 * t1133 + t1137 * t800;
t783 = t1021 * t1137 + t1133 * t800;
t778 = t1131 * t889 + t1310;
t777 = -t1131 * t895 + t1346;
t776 = -t1135 * t889 + t1312;
t775 = t1135 * t895 + t1349;
t774 = t1131 * t887 + t1310;
t773 = -t1131 * t1269 - t1346;
t772 = -t1135 * t887 + t1312;
t771 = t1135 * t1269 - t1349;
t770 = -t1128 * t852 + t1129 * t854;
t769 = -t1128 * t851 + t1129 * t853;
t768 = t1128 * t854 + t1129 * t852;
t767 = t1128 * t853 + t1129 * t851;
t762 = -t1132 * t847 + t1136 * t849;
t761 = -t1132 * t846 + t1136 * t848;
t756 = -t1132 * t837 + t1136 * t838;
t755 = t1132 * t838 + t1136 * t837;
t754 = -t1133 * t1277 + t1137 * t794;
t753 = t1133 * t794 + t1137 * t1277;
t744 = -pkin(1) * t799 - pkin(2) * t870;
t743 = t1133 * t933 + t1137 * t756;
t742 = t1133 * t756 - t1137 * t933;
t739 = t1238 + t1353;
t738 = t1251 - t1315;
t737 = -qJ(3) * t867 - t1128 * t820 + t1129 * t869;
t736 = pkin(3) * t946 + pkin(8) * t741;
t733 = -pkin(7) * t898 - t1132 * t819 + t1136 * t823;
t732 = -pkin(4) * t1306 - pkin(5) * t1270 + qJ(6) * t888;
t731 = -qJ(3) * t837 - t1128 * t812 + t1129 * t850;
t730 = (pkin(5) * t1048 - (2 * qJD(6))) * t1024 + t1146;
t729 = pkin(2) * t1277 + qJ(3) * t868 + t1128 * t869 + t1129 * t820;
t720 = -pkin(7) * t799 - qJ(3) * t1229 - t1132 * t843;
t719 = -pkin(8) * t852 - t740;
t716 = -t1128 * t776 + t1129 * t778;
t715 = -t1128 * t775 + t1129 * t777;
t714 = t1128 * t778 + t1129 * t776;
t713 = t1128 * t777 + t1129 * t775;
t712 = -t1128 * t772 + t1129 * t774;
t711 = -t1128 * t771 + t1129 * t773;
t710 = t1128 * t774 + t1129 * t772;
t709 = t1128 * t773 + t1129 * t771;
t708 = -t1132 * t768 + t1136 * t770;
t707 = -t1132 * t767 + t1136 * t769;
t706 = t1132 * t770 + t1136 * t768;
t701 = -pkin(2) * t933 + qJ(3) * t838 + t1128 * t850 + t1129 * t812;
t700 = -pkin(3) * t962 + pkin(8) * t854 + t741;
t699 = qJ(6) * t1259 - t1154;
t698 = -pkin(5) * t1259 + t1157;
t697 = t735 + t1354;
t696 = t734 - t1317;
t691 = t1133 * t962 + t1137 * t708;
t690 = t1133 * t708 - t1137 * t962;
t689 = (-t889 - t1218) * pkin(5) + t1145;
t688 = -pkin(5) * t1218 + t1145 + t1296;
t687 = -pkin(1) * t793 - pkin(2) * t867 - pkin(3) * t944 + t818;
t686 = qJ(6) * t1228 + t1154;
t685 = pkin(5) * t1228 + t1157;
t680 = -pkin(1) * t755 - pkin(2) * t837 - pkin(3) * t922 + t817;
t679 = t1129 * t741 - t1253;
t678 = t1128 * t741 + t1252;
t677 = -t1317 + (-t1259 - t1268) * qJ(6) + (-t1279 - t966) * pkin(5) + t1171;
t676 = -t1130 * t689 - t1254 * t889 - t1315;
t675 = -pkin(5) * t1249 + t1134 * t688 - t1353;
t672 = -t1354 - qJ(6) * t1278 + (t1259 - t942) * pkin(5) - t1157;
t671 = -pkin(1) * t706 - pkin(2) * t768 - pkin(3) * t852;
t666 = -t1132 * t714 + t1136 * t716;
t665 = -t1132 * t713 + t1136 * t715;
t664 = t1132 * t716 + t1136 * t714;
t663 = t1132 * t715 + t1136 * t713;
t662 = -t1132 * t710 + t1136 * t712;
t661 = -t1132 * t709 + t1136 * t711;
t660 = t1132 * t712 + t1136 * t710;
t659 = t1132 * t711 + t1136 * t709;
t658 = -t673 - t1336;
t655 = t1131 * t797 + t1135 * t674;
t654 = t1131 * t674 - t1135 * t797;
t653 = -t1130 * t699 + t1134 * t698;
t652 = t1130 * t698 + t1134 * t699;
t651 = -pkin(7) * t793 - t1132 * t729 + t1136 * t737;
t646 = t1137 * t666 + t1311;
t645 = t1137 * t665 - t1348;
t644 = t1133 * t666 - t1309;
t643 = t1133 * t665 + t1345;
t642 = t1137 * t662 + t1311;
t641 = t1137 * t661 + t1348;
t640 = t1133 * t662 - t1309;
t639 = t1133 * t661 - t1345;
t638 = -pkin(8) * t775 - t1131 * t697 + t1135 * t739;
t637 = -pkin(8) * t772 - t1131 * t696 + t1135 * t738;
t636 = -qJ(3) * t768 - t1128 * t700 + t1129 * t719;
t635 = -pkin(7) * t755 - t1132 * t701 + t1136 * t731;
t634 = -pkin(2) * t962 + qJ(3) * t770 + t1128 * t719 + t1129 * t700;
t633 = pkin(8) * t777 + t1131 * t739 + t1135 * t697 + t1355;
t632 = pkin(8) * t774 + t1131 * t738 + t1135 * t696 - t1318;
t627 = -t1132 * t678 + t1136 * t679;
t626 = t1132 * t679 + t1136 * t678;
t625 = -t1130 * t685 + t1134 * t686 - t1336;
t624 = -pkin(8) * t1252 - qJ(3) * t678 - t1128 * t736;
t623 = -t1133 * t946 + t1137 * t627;
t622 = t1133 * t627 + t1137 * t946;
t621 = pkin(2) * t946 - pkin(8) * t1253 + qJ(3) * t679 + t1129 * t736;
t620 = t1135 * t658 + t1257 * t1306 - t1350;
t619 = t1131 * t730 + t1135 * t653;
t618 = t1131 * t653 - t1135 * t730;
t617 = t1131 * t658 + t1189 * t1306 + t1351;
t616 = -pkin(8) * t776 - t1131 * t677 + t1135 * t676;
t615 = -pkin(8) * t771 - t1131 * t672 + t1135 * t675;
t614 = -pkin(9) * t652 + (pkin(5) * t1130 - t1254) * t730;
t613 = pkin(8) * t778 + t1131 * t676 + t1135 * t677 - t1318;
t612 = pkin(8) * t773 + t1131 * t675 + t1135 * t672 - t1355;
t611 = -t1128 * t654 + t1129 * t655;
t610 = t1128 * t655 + t1129 * t654;
t609 = -t1131 * t732 + t1135 * t625 - t1350;
t608 = -pkin(4) * t652 - pkin(5) * t699 - qJ(6) * t698;
t607 = -pkin(1) * t663 - pkin(2) * t713 - pkin(3) * t775 - pkin(4) * t895 - t1251 - t1352;
t606 = -pkin(1) * t660 - pkin(2) * t710 - pkin(3) * t772 + pkin(4) * t887 + t1238 - t1316;
t605 = -pkin(3) * t1306 + t1131 * t625 + t1135 * t732 + t1351;
t604 = -pkin(1) * t626 - pkin(2) * t678 - pkin(3) * t740;
t603 = -pkin(1) * t664 - pkin(2) * t714 - pkin(3) * t776 - t1134 * t689 + t1178 * t889 - t1316;
t602 = -pkin(8) * t654 + (-pkin(9) * t1135 + t1257) * t673;
t601 = -pkin(1) * t659 - pkin(2) * t709 - pkin(3) * t771 + t1352 - t1130 * t688 + (-pkin(4) - t1256) * t1269;
t600 = -qJ(3) * t713 - t1128 * t633 + t1129 * t638;
t599 = -qJ(3) * t710 - t1128 * t632 + t1129 * t637;
t598 = t1370 - t674;
t597 = -pkin(7) * t706 - t1132 * t634 + t1136 * t636;
t596 = -t1128 * t618 + t1129 * t619;
t595 = t1128 * t619 + t1129 * t618;
t594 = qJ(3) * t715 + t1128 * t638 + t1129 * t633 + t1356;
t593 = qJ(3) * t712 + t1128 * t637 + t1129 * t632 - t1319;
t592 = -t1130 * t686 - t1134 * t685 + t1370;
t591 = pkin(8) * t655 + (-pkin(9) * t1131 + t1189) * t673;
t590 = -t1128 * t617 + t1129 * t620 - t1364;
t589 = t1128 * t620 + t1129 * t617 + t1362;
t588 = -t1132 * t610 + t1136 * t611;
t587 = t1132 * t611 + t1136 * t610;
t586 = -qJ(3) * t714 - t1128 * t613 + t1129 * t616;
t585 = -pkin(7) * t626 - t1132 * t621 + t1136 * t624;
t584 = -qJ(3) * t709 - t1128 * t612 + t1129 * t615;
t583 = qJ(3) * t716 + t1128 * t616 + t1129 * t613 - t1319;
t582 = qJ(3) * t711 + t1128 * t615 + t1129 * t612 - t1356;
t581 = -t1128 * t605 + t1129 * t609 - t1364;
t580 = t1133 * t673 + t1137 * t588;
t579 = t1133 * t588 - t1137 * t673;
t578 = t1128 * t609 + t1129 * t605 + t1362;
t577 = -pkin(8) * t618 - t1131 * t608 + t1135 * t614;
t576 = -t1132 * t595 + t1136 * t596;
t575 = t1132 * t596 + t1136 * t595;
t574 = -pkin(3) * t652 + pkin(8) * t619 + t1131 * t614 + t1135 * t608;
t573 = -pkin(7) * t663 - t1132 * t594 + t1136 * t600;
t572 = -pkin(7) * t660 - t1132 * t593 + t1136 * t599;
t571 = t1133 * t652 + t1137 * t576;
t570 = t1133 * t576 - t1137 * t652;
t569 = -qJ(3) * t610 - t1128 * t591 + t1129 * t602;
t568 = -pkin(1) * t587 - pkin(2) * t610 - pkin(3) * t654 + pkin(4) * t797 - pkin(9) * t674;
t567 = -t1132 * t589 + t1136 * t590 - t1371;
t566 = -pkin(2) * t673 + qJ(3) * t611 + t1128 * t602 + t1129 * t591;
t565 = -pkin(7) * t664 - t1132 * t583 + t1136 * t586;
t564 = -pkin(7) * t659 - t1132 * t582 + t1136 * t584;
t563 = -t1132 * t578 + t1136 * t581 - t1371;
t562 = -pkin(1) * t575 - pkin(2) * t595 - pkin(3) * t618 - pkin(9) * t653 + (t1178 + t1256) * t730;
t561 = -qJ(3) * t595 - t1128 * t574 + t1129 * t577;
t560 = -pkin(2) * t652 + qJ(3) * t596 + t1128 * t577 + t1129 * t574;
t559 = -pkin(7) * t587 - t1132 * t566 + t1136 * t569;
t558 = -pkin(7) * t575 - t1132 * t560 + t1136 * t561;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, -t1102, -t1103, 0, t1069, 0, 0, 0, 0, 0, 0, t1043, t1044, t1066, t1001, 0, 0, 0, 0, 0, 0, t902, t917, t877, t784, 0, 0, 0, 0, 0, 0, t743, t754, t691, t623, 0, 0, 0, 0, 0, 0, t642, t645, t1365, t580, 0, 0, 0, 0, 0, 0, t646, t1365, t641, t571; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t1103, -t1102, 0, t1068, 0, 0, 0, 0, 0, 0, t1041, t1042, t1065, t1000, 0, 0, 0, 0, 0, 0, t899, t916, t876, t783, 0, 0, 0, 0, 0, 0, t742, t753, t690, t622, 0, 0, 0, 0, 0, 0, t640, t643, t1367, t579, 0, 0, 0, 0, 0, 0, t644, t1367, t639, t570; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t1070, t1071, 0, -t1015, 0, 0, 0, 0, 0, 0, t925, t956, t898, t799, 0, 0, 0, 0, 0, 0, t755, t793, t706, t626, 0, 0, 0, 0, 0, 0, t660, t663, t1361, t587, 0, 0, 0, 0, 0, 0, t664, t1361, t659, t575; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1103, 0, -t1102, 0, t1167, -t1085, -t1068, -pkin(6) * t1068, t1079 * t1137 - t1170, t1064 * t1137 - t1105 * t1133, t1074 * t1137 + t1132 * t1195, t1078 * t1137 + t1170, t1072 * t1137 + t1122 * t1133, t1094 * t1137 + t1121, -pkin(6) * t1041 - t1031 * t1133 + t1034 * t1137, -pkin(6) * t1042 - t1032 * t1133 + t1035 * t1137, -pkin(6) * t1065 + t1015 * t1137, -pkin(6) * t1000 - (pkin(1) * t1133 - pkin(7) * t1137) * t1015, t1137 * t961 - t1181, -t1060 * t1133 + t1137 * t900, -t1029 * t1133 + t1137 * t958, t1137 * t960 + t1181, t1027 * t1133 + t1137 * t957, t1137 * t982 + t1121, -pkin(6) * t899 - t1133 * t824 + t1137 * t811, -pkin(6) * t916 - t1133 * t844 + t1137 * t821, -pkin(6) * t876 - t1133 * t845 + t1137 * t733, -pkin(6) * t783 - t1133 * t744 + t1137 * t720, t1137 * t762 + t1183, -t1133 * t996 + t1137 * t707, -t1133 * t937 + t1137 * t801, t1137 * t761 - t1183, t1133 * t934 + t1137 * t802, t1133 * t1192 + t1137 * t822, -pkin(6) * t742 - t1133 * t680 + t1137 * t635, -pkin(6) * t753 - t1133 * t687 + t1137 * t651, -pkin(6) * t690 - t1133 * t671 + t1137 * t597, -pkin(6) * t622 - t1133 * t604 + t1137 * t585, t1303, t1377, t1366, t1321, -t1376, t1322, -pkin(6) * t640 - t1133 * t606 + t1137 * t572, -pkin(6) * t643 - t1133 * t607 + t1137 * t573, -t1133 * t598 + t1137 * t567 - t1372, -pkin(6) * t579 - t1133 * t568 + t1137 * t559, t1303, t1366, -t1377, t1322, t1376, t1321, -pkin(6) * t644 - t1133 * t603 + t1137 * t565, -t1133 * t592 + t1137 * t563 - t1372, -pkin(6) * t639 - t1133 * t601 + t1137 * t564, -pkin(6) * t570 - t1133 * t562 + t1137 * t558; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1102, 0, t1103, 0, t1085, t1167, t1069, pkin(6) * t1069, t1079 * t1133 + t1169, t1064 * t1133 + t1105 * t1137, t1074 * t1133 - t1132 * t1194, t1078 * t1133 - t1169, t1072 * t1133 - t1122 * t1137, t1094 * t1133 - t1193, pkin(6) * t1043 + t1031 * t1137 + t1034 * t1133, pkin(6) * t1044 + t1032 * t1137 + t1035 * t1133, pkin(6) * t1066 + t1015 * t1133, pkin(6) * t1001 - (-pkin(1) * t1137 - pkin(7) * t1133) * t1015, t1133 * t961 + t1180, t1060 * t1137 + t1133 * t900, t1029 * t1137 + t1133 * t958, t1133 * t960 - t1180, -t1027 * t1137 + t1133 * t957, t1133 * t982 - t1193, pkin(6) * t902 + t1133 * t811 + t1137 * t824, pkin(6) * t917 + t1133 * t821 + t1137 * t844, pkin(6) * t877 + t1133 * t733 + t1137 * t845, pkin(6) * t784 + t1133 * t720 + t1137 * t744, t1133 * t762 - t1182, t1133 * t707 + t1137 * t996, t1133 * t801 + t1137 * t937, t1133 * t761 + t1182, t1133 * t802 - t1137 * t934, t1133 * t822 - t1137 * t1192, pkin(6) * t743 + t1133 * t635 + t1137 * t680, pkin(6) * t754 + t1133 * t651 + t1137 * t687, pkin(6) * t691 + t1133 * t597 + t1137 * t671, pkin(6) * t623 + t1133 * t585 + t1137 * t604, t1305, t1379, t1368, t1325, -t1378, t1326, pkin(6) * t642 + t1133 * t572 + t1137 * t606, pkin(6) * t645 + t1133 * t573 + t1137 * t607, t1133 * t567 + t1137 * t598 + t1373, pkin(6) * t580 + t1133 * t559 + t1137 * t568, t1305, t1368, -t1379, t1326, t1378, t1325, pkin(6) * t646 + t1133 * t565 + t1137 * t603, t1133 * t563 + t1137 * t592 + t1373, pkin(6) * t641 + t1133 * t564 + t1137 * t601, pkin(6) * t571 + t1133 * t558 + t1137 * t562; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1108, t1109, 0, 0, (t1098 + t1179) * t1132, t1097 * t1136 + t1100 * t1132, t1111 * t1136 + t1207, (t1099 - t1188) * t1136, t1113 * t1132 + t1205, 0, pkin(1) * t1100 + pkin(7) * t1073 + t1208, -pkin(1) * t1097 + pkin(7) * t1075 - t1209, pkin(1) * t1104 + pkin(7) * t1101 + t1016, pkin(1) * t1091 + pkin(7) * t1016, t1012 * t1136 + t1013 * t1132, t1132 * t971 + t1136 * t969, t1003 * t1136 + t1006 * t1132, t1010 * t1136 + t1011 * t1132, t1002 * t1136 + t1005 * t1132, t1045 * t1136 + t1046 * t1132, -pkin(1) * t1026 + pkin(7) * t926 + t1132 * t953 + t1136 * t918, -pkin(1) * t1028 + pkin(7) * t959 + t1132 * t963 + t1136 * t924, -pkin(1) * t1025 + pkin(7) * t901 + t1132 * t823 + t1136 * t819, pkin(1) * t1021 + pkin(7) * t800 - qJ(3) * t1240 + t1136 * t843, t1132 * t849 + t1136 * t847, t1132 * t769 + t1136 * t767, t1132 * t874 + t1136 * t872, t1132 * t848 + t1136 * t846, t1132 * t875 + t1136 * t873, t1132 * t911 + t1136 * t910, -pkin(1) * t933 + pkin(7) * t756 + t1132 * t731 + t1136 * t701, pkin(1) * t1277 + pkin(7) * t794 + t1132 * t737 + t1136 * t729, -pkin(1) * t962 + pkin(7) * t708 + t1132 * t636 + t1136 * t634, pkin(1) * t946 + pkin(7) * t627 + t1132 * t624 + t1136 * t621, t1280, t1375, t1359, t1302, -t1374, t1300, pkin(7) * t662 + t1132 * t599 + t1136 * t593 - t1320, pkin(7) * t665 + t1132 * t600 + t1136 * t594 + t1357, t1132 * t590 + t1136 * t589 + t1369, -pkin(1) * t673 + pkin(7) * t588 + t1132 * t569 + t1136 * t566, t1280, t1359, -t1375, t1300, t1374, t1302, pkin(7) * t666 + t1132 * t586 + t1136 * t583 - t1320, t1132 * t581 + t1136 * t578 + t1369, pkin(7) * t661 + t1132 * t584 + t1136 * t582 - t1357, -pkin(1) * t652 + pkin(7) * t576 + t1132 * t561 + t1136 * t560;];
tauB_reg  = t1;