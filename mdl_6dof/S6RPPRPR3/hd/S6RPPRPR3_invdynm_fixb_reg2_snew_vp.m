% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6RPPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:10
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6RPPRPR3_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR3_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR3_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_invdynm_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:10:19
% EndTime: 2019-05-05 14:10:43
% DurationCPUTime: 25.94s
% Computational Cost: add. (101247->778), mult. (210699->1012), div. (0->0), fcn. (135785->10), ass. (0->522)
t1246 = sin(pkin(9));
t1248 = cos(pkin(9));
t1259 = qJD(1) ^ 2;
t1206 = -qJDD(1) * t1248 + t1246 * t1259;
t1242 = g(3) - qJDD(2);
t1166 = qJ(2) * t1206 - t1242 * t1246;
t1254 = sin(qJ(1));
t1257 = cos(qJ(1));
t1205 = qJDD(1) * t1246 + t1248 * t1259;
t1322 = t1205 * t1257 - t1206 * t1254;
t1329 = -qJ(2) * t1205 + t1242 * t1248;
t1428 = -pkin(6) * t1322 + t1166 * t1254 + t1257 * t1329;
t1245 = sin(pkin(10));
t1247 = cos(pkin(10));
t1253 = sin(qJ(4));
t1256 = cos(qJ(4));
t1184 = (t1245 * t1256 + t1247 * t1253) * qJD(1);
t1390 = qJD(1) * t1256;
t1391 = qJD(1) * t1253;
t1186 = -t1245 * t1391 + t1247 * t1390;
t1368 = t1186 * t1184;
t1423 = qJDD(4) - t1368;
t1430 = t1245 * t1423;
t1429 = t1247 * t1423;
t1215 = g(1) * t1254 - g(2) * t1257;
t1198 = qJDD(1) * pkin(1) + t1215;
t1216 = g(1) * t1257 + g(2) * t1254;
t1199 = -pkin(1) * t1259 - t1216;
t1128 = -t1198 * t1248 + t1199 * t1246;
t1129 = t1246 * t1198 + t1248 * t1199;
t1325 = t1128 * t1246 + t1129 * t1248;
t1054 = t1128 * t1248 - t1129 * t1246;
t1381 = t1054 * t1254;
t1425 = t1257 * t1325 + t1381;
t1233 = t1256 * qJDD(1);
t1334 = qJD(4) * t1391;
t1203 = t1233 - t1334;
t1333 = qJD(4) * t1390;
t1346 = t1253 * qJDD(1);
t1285 = -t1333 - t1346;
t1131 = t1203 * t1247 + t1245 * t1285;
t1389 = qJD(4) * t1184;
t1085 = -t1389 + t1131;
t1380 = t1054 * t1257;
t1424 = -t1254 * t1325 + t1380;
t1320 = t1205 * t1254 + t1206 * t1257;
t1406 = pkin(6) * t1320 + t1166 * t1257 - t1254 * t1329;
t1252 = sin(qJ(6));
t1255 = cos(qJ(6));
t1148 = -qJD(4) * t1255 + t1186 * t1252;
t1150 = qJD(4) * t1252 + t1186 * t1255;
t1080 = t1150 * t1148;
t1323 = t1203 * t1245 - t1247 * t1285;
t1126 = qJDD(6) + t1323;
t1410 = -t1080 + t1126;
t1420 = t1252 * t1410;
t1419 = t1255 * t1410;
t1178 = qJD(4) * t1186;
t1082 = t1178 + t1323;
t1386 = qJDD(1) * qJ(3);
t1263 = -pkin(2) * t1259 + t1129 + t1386;
t1405 = 2 * qJD(3);
t1345 = qJD(1) * t1405;
t1108 = t1263 + t1345;
t1244 = qJDD(1) * pkin(2);
t1304 = qJDD(3) + t1128;
t1110 = -qJ(3) * t1259 - t1244 + t1304;
t1033 = t1108 * t1246 - t1110 * t1248;
t1326 = t1108 * t1248 + t1110 * t1246;
t1418 = -t1033 * t1254 + t1257 * t1326;
t1417 = t1033 * t1257 + t1254 * t1326;
t1119 = pkin(5) * t1184 - pkin(8) * t1186;
t1258 = qJD(4) ^ 2;
t1388 = qJD(5) * t1184;
t1174 = -0.2e1 * t1388;
t1103 = -qJDD(1) * pkin(7) + t1110;
t1073 = t1103 * t1253 - t1242 * t1256;
t1212 = qJD(4) * pkin(4) - qJ(5) * t1390;
t1240 = t1253 ^ 2;
t1357 = t1240 * t1259;
t1039 = -pkin(4) * t1357 + qJ(5) * t1285 - qJD(4) * t1212 + t1073;
t1352 = t1256 * t1259;
t1353 = t1256 * t1103;
t1392 = qJD(1) * qJD(4);
t1260 = qJDD(4) * pkin(4) - t1203 * qJ(5) + t1353 + (-pkin(4) * t1352 - qJ(5) * t1392 + t1242) * t1253;
t1349 = t1039 * t1247 + t1245 * t1260;
t942 = t1174 + t1349;
t919 = -pkin(5) * t1258 + qJDD(4) * pkin(8) - t1119 * t1184 + t942;
t1250 = t1259 * pkin(7);
t1044 = qJDD(5) - t1285 * pkin(4) - qJ(5) * t1357 - t1250 + (t1212 * t1256 + t1405) * qJD(1) + t1263;
t950 = pkin(5) * t1082 - pkin(8) * t1085 + t1044;
t880 = t1252 * t919 - t1255 * t950;
t881 = t1252 * t950 + t1255 * t919;
t826 = t1252 * t880 + t1255 * t881;
t1072 = -t1253 * t1242 - t1353;
t1005 = -t1072 * t1256 + t1073 * t1253;
t1084 = -t1323 + t1178;
t1146 = t1148 ^ 2;
t1147 = t1150 ^ 2;
t1180 = qJD(6) + t1184;
t1179 = t1180 ^ 2;
t1182 = t1184 ^ 2;
t1183 = t1186 ^ 2;
t1404 = -pkin(2) - pkin(7);
t1403 = pkin(1) * t1205;
t1402 = pkin(1) * t1206;
t1401 = pkin(3) * t1005;
t1093 = t1108 - t1250;
t1400 = pkin(3) * t1093;
t1241 = t1256 ^ 2;
t1347 = t1240 + t1241;
t1207 = t1347 * qJDD(1);
t1399 = pkin(3) * t1207;
t1398 = pkin(5) * t1245;
t1397 = pkin(7) * t1005;
t1327 = t1245 * t1039 - t1247 * t1260;
t918 = -qJDD(4) * pkin(5) - t1258 * pkin(8) + (0.2e1 * qJD(5) + t1119) * t1186 + t1327;
t1396 = -pkin(5) * t918 + pkin(8) * t826;
t915 = t1252 * t918;
t1387 = qJD(5) * t1186;
t1344 = 0.2e1 * t1387;
t941 = t1327 + t1344;
t885 = t1245 * t942 - t1247 * t941;
t1394 = t1253 * t885;
t916 = t1255 * t918;
t1393 = t1256 * t885;
t1041 = t1080 + t1126;
t1385 = t1041 * t1252;
t1384 = t1041 * t1255;
t1383 = t1044 * t1245;
t1382 = t1044 * t1247;
t1123 = qJDD(4) + t1368;
t1378 = t1123 * t1245;
t1377 = t1123 * t1247;
t1372 = t1180 * t1252;
t1371 = t1180 * t1255;
t1370 = t1184 * t1245;
t1369 = t1184 * t1247;
t1367 = t1186 * t1245;
t1366 = t1186 * t1247;
t1202 = 0.2e1 * t1333 + t1346;
t1151 = t1202 * t1253;
t1363 = t1207 * t1246;
t1362 = t1207 * t1248;
t1222 = t1253 * t1352;
t1213 = qJDD(4) + t1222;
t1361 = t1213 * t1253;
t1360 = t1213 * t1256;
t1214 = qJDD(4) - t1222;
t1359 = t1214 * t1253;
t1358 = t1214 * t1256;
t1356 = t1241 * t1259;
t1354 = t1253 * t1093;
t1089 = t1256 * t1093;
t1348 = -pkin(2) * t1110 + qJ(3) * t1108;
t819 = t1245 * t826 - t1247 * t918;
t1343 = pkin(4) * t819 + t1396;
t1299 = -t1252 * qJDD(4) - t1255 * t1131;
t1017 = (qJD(6) + t1180) * t1148 + t1299;
t1066 = -t1147 - t1179;
t980 = -t1066 * t1252 - t1384;
t1342 = pkin(5) * t1017 + pkin(8) * t980 + t915;
t1107 = t1180 * t1150;
t1324 = qJDD(4) * t1255 - t1252 * t1131;
t1286 = qJD(6) * t1150 - t1324;
t1013 = -t1107 - t1286;
t1057 = -t1179 - t1146;
t975 = t1057 * t1255 - t1420;
t1341 = pkin(5) * t1013 + pkin(8) * t975 - t916;
t1340 = -pkin(5) * t1247 - pkin(4);
t886 = t1245 * t941 + t1247 * t942;
t835 = t1253 * t886 + t1393;
t884 = pkin(4) * t885;
t1339 = -pkin(3) * t835 - t884;
t1338 = t1245 * t1080;
t1337 = t1247 * t1080;
t1336 = t1246 * t1368;
t1335 = t1248 * t1368;
t1173 = -t1183 - t1258;
t1062 = t1173 * t1247 - t1378;
t1332 = pkin(4) * t1062 - t1349;
t1086 = t1389 + t1131;
t1021 = t1084 * t1245 - t1086 * t1247;
t1019 = pkin(4) * t1021;
t1023 = t1084 * t1247 + t1086 * t1245;
t935 = t1021 * t1256 + t1023 * t1253;
t1331 = -pkin(3) * t935 - t1019;
t1221 = -t1258 - t1356;
t1155 = t1221 * t1256 - t1361;
t1330 = -pkin(7) * t1155 + t1089;
t1318 = -t1215 * t1254 - t1216 * t1257;
t1050 = t1146 + t1147;
t1014 = (-qJD(6) + t1180) * t1150 + t1324;
t1058 = -qJD(6) * t1148 - t1299;
t1106 = t1180 * t1148;
t1016 = t1058 + t1106;
t930 = t1014 * t1255 + t1016 * t1252;
t1317 = pkin(5) * t1050 + pkin(8) * t930 + t826;
t911 = t1017 * t1247 + t1245 * t980;
t1316 = pkin(4) * t911 + t1342;
t904 = t1013 * t1247 + t1245 * t975;
t1315 = pkin(4) * t904 + t1341;
t1314 = t1246 * t1222;
t1313 = t1248 * t1222;
t1312 = -pkin(3) * t1155 + t1073;
t1311 = -pkin(2) * t1005 + qJ(3) * t1093 - t1397;
t1209 = qJDD(1) * t1257 - t1254 * t1259;
t1310 = -pkin(6) * t1209 - g(3) * t1254;
t894 = t1050 * t1247 + t1245 * t930;
t1309 = pkin(4) * t894 + t1317;
t1121 = -t1258 - t1182;
t1047 = t1121 * t1245 + t1429;
t1308 = pkin(4) * t1047 - t1327;
t1307 = pkin(3) * t1202 + t1089;
t1204 = t1233 - 0.2e1 * t1334;
t1306 = pkin(3) * t1204 - t1354;
t1219 = -t1258 - t1357;
t1153 = t1219 * t1253 + t1358;
t1305 = -pkin(7) * t1153 + t1354;
t820 = t1245 * t918 + t1247 * t826;
t798 = t1253 * t820 + t1256 * t819;
t1303 = -pkin(3) * t798 - t1343;
t825 = t1252 * t881 - t1255 * t880;
t1006 = t1072 * t1253 + t1073 * t1256;
t1302 = t1215 * t1257 - t1216 * t1254;
t1065 = -t1173 * t1245 - t1377;
t991 = t1062 * t1256 + t1065 * t1253;
t1300 = -pkin(3) * t991 - t1332;
t1298 = -pkin(2) * t1155 + qJ(3) * t1204 + t1330;
t787 = qJ(5) * t820 + (-pkin(8) * t1245 + t1340) * t825;
t792 = -qJ(5) * t819 + (-pkin(8) * t1247 + t1398) * t825;
t1297 = -pkin(7) * t798 - t1253 * t787 + t1256 * t792;
t928 = t1014 * t1252 - t1016 * t1255;
t821 = -pkin(8) * t928 - t825;
t895 = -t1050 * t1245 + t1247 * t930;
t802 = qJ(5) * t895 + t1245 * t821 + t1340 * t928;
t807 = -qJ(5) * t894 + t1247 * t821 + t1398 * t928;
t850 = t1253 * t895 + t1256 * t894;
t1296 = -pkin(7) * t850 - t1253 * t802 + t1256 * t807;
t974 = t1057 * t1252 + t1419;
t853 = -pkin(5) * t974 + t880;
t882 = -pkin(8) * t974 + t915;
t905 = -t1013 * t1245 + t1247 * t975;
t811 = -pkin(4) * t974 + qJ(5) * t905 + t1245 * t882 + t1247 * t853;
t815 = -qJ(5) * t904 - t1245 * t853 + t1247 * t882;
t859 = t1253 * t905 + t1256 * t904;
t1295 = -pkin(7) * t859 - t1253 * t811 + t1256 * t815;
t979 = t1066 * t1255 - t1385;
t854 = -pkin(5) * t979 + t881;
t883 = -pkin(8) * t979 + t916;
t912 = -t1017 * t1245 + t1247 * t980;
t812 = -pkin(4) * t979 + qJ(5) * t912 + t1245 * t883 + t1247 * t854;
t817 = -qJ(5) * t911 - t1245 * t854 + t1247 * t883;
t863 = t1253 * t912 + t1256 * t911;
t1294 = -pkin(7) * t863 - t1253 * t812 + t1256 * t817;
t1079 = -t1182 - t1183;
t858 = -pkin(4) * t1079 + qJ(5) * t1023 + t886;
t871 = -qJ(5) * t1021 - t885;
t1293 = -pkin(7) * t935 - t1253 * t858 + t1256 * t871;
t1048 = t1121 * t1247 - t1430;
t938 = -pkin(4) * t1082 + qJ(5) * t1048 - t1382;
t963 = -qJ(5) * t1047 + t1383;
t964 = t1047 * t1256 + t1048 * t1253;
t1292 = -pkin(7) * t964 - t1253 * t938 + t1256 * t963;
t940 = -pkin(4) * t1085 + qJ(5) * t1065 + t1383;
t981 = -qJ(5) * t1062 + t1382;
t1291 = -pkin(7) * t991 - t1253 * t940 + t1256 * t981;
t1290 = pkin(7) * t1207 - t1005;
t1289 = -0.2e1 * t1244 + t1304;
t1288 = -pkin(3) * t859 - t1315;
t1287 = -pkin(3) * t863 - t1316;
t1284 = -pkin(3) * t964 - t1308;
t1283 = -pkin(3) * t850 - t1309;
t1282 = -pkin(2) * t1153 + qJ(3) * t1202 + t1305;
t1281 = pkin(3) * t825 - t1253 * t792 - t1256 * t787;
t1280 = pkin(3) * t928 - t1253 * t807 - t1256 * t802;
t1279 = pkin(3) * t974 - t1253 * t815 - t1256 * t811;
t1278 = pkin(3) * t979 - t1253 * t817 - t1256 * t812;
t1277 = pkin(3) * t1079 - t1253 * t871 - t1256 * t858;
t1276 = pkin(3) * t1082 - t1253 * t963 - t1256 * t938;
t1275 = pkin(3) * t1085 - t1253 * t981 - t1256 * t940;
t1274 = -pkin(2) * t798 + qJ(3) * t825 + t1297;
t1273 = -pkin(2) * t850 + qJ(3) * t928 + t1296;
t1272 = -pkin(2) * t859 + qJ(3) * t974 + t1295;
t1271 = -pkin(2) * t863 + qJ(3) * t979 + t1294;
t1270 = -pkin(3) * t1153 + t1072;
t1269 = -pkin(2) * t935 + qJ(3) * t1079 + t1293;
t1268 = -pkin(2) * t964 + qJ(3) * t1082 + t1292;
t1267 = -pkin(2) * t991 + qJ(3) * t1085 + t1291;
t872 = -pkin(4) * t1044 + qJ(5) * t886;
t1266 = -pkin(7) * t835 - qJ(5) * t1393 - t1253 * t872;
t1210 = t1347 * t1259;
t1265 = pkin(2) * t1207 - qJ(3) * t1210 + t1290;
t1264 = pkin(3) * t1044 + qJ(5) * t1394 - t1256 * t872;
t1262 = t1129 + 0.2e1 * t1386 + t1345;
t1261 = -pkin(2) * t835 + qJ(3) * t1044 + t1266;
t1230 = t1248 * qJDD(4);
t1229 = t1246 * qJDD(4);
t1220 = t1258 - t1356;
t1218 = -t1258 + t1357;
t1211 = (-t1240 + t1241) * t1259;
t1208 = qJDD(1) * t1254 + t1257 * t1259;
t1196 = t1347 * t1392;
t1181 = -pkin(6) * t1208 + g(3) * t1257;
t1176 = -0.2e1 * t1387;
t1175 = 0.2e1 * t1388;
t1172 = -t1183 + t1258;
t1171 = t1182 - t1258;
t1164 = t1203 * t1253 + t1241 * t1392;
t1163 = t1240 * t1392 + t1256 * t1285;
t1162 = -t1196 * t1246 + t1230;
t1161 = t1196 * t1248 + t1229;
t1160 = -t1221 * t1253 - t1360;
t1159 = -t1220 * t1253 + t1358;
t1158 = (t1203 - t1334) * t1256;
t1157 = t1219 * t1256 - t1359;
t1156 = t1218 * t1256 - t1361;
t1154 = t1220 * t1256 + t1359;
t1152 = t1218 * t1253 + t1360;
t1143 = -t1210 * t1248 - t1363;
t1142 = -t1210 * t1246 + t1362;
t1133 = -t1202 * t1256 - t1204 * t1253;
t1132 = t1204 * t1256 - t1151;
t1125 = t1183 - t1182;
t1118 = t1163 * t1246 - t1313;
t1117 = t1164 * t1246 + t1313;
t1116 = -t1163 * t1248 - t1314;
t1115 = -t1164 * t1248 + t1314;
t1114 = t1154 * t1246 + t1233 * t1248;
t1113 = t1152 * t1246 - t1248 * t1346;
t1112 = -t1154 * t1248 + t1233 * t1246;
t1111 = -t1152 * t1248 - t1246 * t1346;
t1105 = (t1367 - t1369) * qJD(4);
t1104 = (-t1366 - t1370) * qJD(4);
t1101 = t1155 * t1246 + t1204 * t1248;
t1100 = t1153 * t1246 + t1202 * t1248;
t1099 = -t1155 * t1248 + t1204 * t1246;
t1098 = -t1153 * t1248 + t1202 * t1246;
t1097 = -t1147 + t1179;
t1096 = t1146 - t1179;
t1095 = -t1128 - t1402;
t1094 = -t1129 - t1403;
t1088 = t1132 * t1246 + t1211 * t1248;
t1087 = -t1132 * t1248 + t1211 * t1246;
t1081 = t1289 + t1402;
t1076 = t1147 - t1146;
t1075 = t1262 + t1403;
t1071 = -qJD(4) * t1367 + t1131 * t1247;
t1070 = qJD(4) * t1366 + t1131 * t1245;
t1069 = qJD(4) * t1369 + t1245 * t1323;
t1068 = qJD(4) * t1370 - t1247 * t1323;
t1064 = -t1172 * t1245 + t1429;
t1063 = t1171 * t1247 - t1378;
t1061 = t1172 * t1247 + t1430;
t1060 = t1171 * t1245 + t1377;
t1051 = pkin(1) * t1054;
t1046 = pkin(1) * t1242 + qJ(2) * t1325;
t1031 = (-t1148 * t1255 + t1150 * t1252) * t1180;
t1030 = (-t1148 * t1252 - t1150 * t1255) * t1180;
t1029 = -t1104 * t1253 + t1105 * t1256;
t1028 = t1104 * t1256 + t1105 * t1253;
t1027 = t1028 * t1246 + t1230;
t1026 = -t1028 * t1248 + t1229;
t1025 = -qJ(3) * t1160 - t1312;
t1024 = -qJ(3) * t1157 - t1270;
t1022 = -t1082 * t1247 - t1085 * t1245;
t1020 = -t1082 * t1245 + t1085 * t1247;
t1018 = -qJ(2) * t1033 + (-pkin(2) * t1246 + qJ(3) * t1248) * t1242;
t1015 = t1058 - t1106;
t1012 = -t1107 + t1286;
t1009 = qJ(2) * t1326 + (pkin(2) * t1248 + qJ(3) * t1246 + pkin(1)) * t1242;
t1008 = t1157 * t1404 + t1307;
t1007 = t1160 * t1404 + t1306;
t1003 = t1058 * t1255 - t1150 * t1372;
t1002 = t1058 * t1252 + t1150 * t1371;
t1001 = t1148 * t1371 + t1252 * t1286;
t1000 = -t1148 * t1372 + t1255 * t1286;
t999 = -t1070 * t1253 + t1071 * t1256;
t998 = -t1068 * t1253 + t1069 * t1256;
t997 = t1070 * t1256 + t1071 * t1253;
t996 = t1068 * t1256 + t1069 * t1253;
t995 = pkin(3) * t1210 + t1006;
t994 = -t1062 * t1253 + t1065 * t1256;
t993 = -t1061 * t1253 + t1064 * t1256;
t992 = -t1060 * t1253 + t1063 * t1256;
t990 = t1061 * t1256 + t1064 * t1253;
t989 = t1060 * t1256 + t1063 * t1253;
t988 = t1031 * t1247 + t1126 * t1245;
t987 = t1031 * t1245 - t1126 * t1247;
t985 = t1096 * t1255 - t1385;
t984 = -t1097 * t1252 + t1419;
t983 = t1096 * t1252 + t1384;
t982 = t1097 * t1255 + t1420;
t973 = pkin(1) * t1099 + t1298;
t972 = pkin(1) * t1098 + t1282;
t970 = t1246 * t997 + t1335;
t969 = t1246 * t996 - t1335;
t968 = -t1248 * t997 + t1336;
t967 = -t1248 * t996 - t1336;
t965 = -t1047 * t1253 + t1048 * t1256;
t961 = t1005 * t1246 + t1093 * t1248;
t960 = -t1005 * t1248 + t1093 * t1246;
t958 = t1003 * t1247 + t1338;
t957 = t1001 * t1247 - t1338;
t956 = t1003 * t1245 - t1337;
t955 = t1001 * t1245 + t1337;
t954 = -pkin(3) * t1362 - qJ(2) * t1142 + t1246 * t995;
t953 = -pkin(3) * t1363 + qJ(2) * t1143 - t1248 * t995;
t952 = pkin(1) * t1033 + t1348;
t951 = pkin(1) * t1142 + t1265;
t948 = t1086 * t1248 + t1246 * t990;
t947 = t1085 * t1248 + t1246 * t991;
t946 = t1084 * t1248 + t1246 * t989;
t945 = t1086 * t1246 - t1248 * t990;
t944 = t1085 * t1246 - t1248 * t991;
t943 = t1084 * t1246 - t1248 * t989;
t937 = -t1021 * t1253 + t1023 * t1256;
t936 = -t1020 * t1253 + t1022 * t1256;
t934 = t1020 * t1256 + t1022 * t1253;
t933 = t1082 * t1248 + t1246 * t964;
t932 = t1082 * t1246 - t1248 * t964;
t929 = t1013 * t1255 - t1015 * t1252;
t927 = t1013 * t1252 + t1015 * t1255;
t924 = -qJ(3) * t1006 + t1401;
t923 = -t1012 * t1245 + t1247 * t985;
t922 = t1016 * t1245 + t1247 * t984;
t921 = t1012 * t1247 + t1245 * t985;
t920 = -t1016 * t1247 + t1245 * t984;
t914 = t1125 * t1248 + t1246 * t934;
t913 = t1125 * t1246 - t1248 * t934;
t909 = -t1253 * t987 + t1256 * t988;
t908 = t1253 * t988 + t1256 * t987;
t907 = -qJ(2) * t1099 - t1007 * t1246 + t1025 * t1248;
t906 = -qJ(2) * t1098 - t1008 * t1246 + t1024 * t1248;
t902 = t1079 * t1248 + t1246 * t935;
t901 = t1079 * t1246 - t1248 * t935;
t900 = t1006 * t1404 + t1400;
t899 = t1076 * t1245 + t1247 * t929;
t898 = -t1076 * t1247 + t1245 * t929;
t897 = -pkin(1) * t1160 + qJ(2) * t1101 + t1007 * t1248 + t1025 * t1246;
t896 = -pkin(1) * t1157 + qJ(2) * t1100 + t1008 * t1248 + t1024 * t1246;
t892 = -t1253 * t956 + t1256 * t958;
t891 = -t1253 * t955 + t1256 * t957;
t890 = t1253 * t958 + t1256 * t956;
t889 = t1253 * t957 + t1256 * t955;
t888 = t1030 * t1248 + t1246 * t908;
t887 = t1030 * t1246 - t1248 * t908;
t877 = pkin(1) * t960 + t1311;
t876 = t1002 * t1248 + t1246 * t890;
t875 = -t1000 * t1248 + t1246 * t889;
t874 = t1002 * t1246 - t1248 * t890;
t873 = -t1000 * t1246 - t1248 * t889;
t870 = -t1253 * t921 + t1256 * t923;
t869 = -t1253 * t920 + t1256 * t922;
t868 = t1253 * t923 + t1256 * t921;
t867 = t1253 * t922 + t1256 * t920;
t865 = -qJ(3) * t994 + t1175 - t1300;
t864 = -t1253 * t911 + t1256 * t912;
t861 = -qJ(3) * t937 - t1331;
t860 = -t1253 * t904 + t1256 * t905;
t856 = -t1253 * t898 + t1256 * t899;
t855 = t1253 * t899 + t1256 * t898;
t852 = -qJ(3) * t965 + t1176 - t1284;
t851 = -t1253 * t894 + t1256 * t895;
t848 = t1246 * t868 + t1248 * t983;
t847 = t1246 * t867 + t1248 * t982;
t846 = t1246 * t983 - t1248 * t868;
t845 = t1246 * t982 - t1248 * t867;
t844 = t1404 * t994 + t1275;
t843 = -qJ(2) * t960 - t1246 * t900 + t1248 * t924;
t842 = t1246 * t863 + t1248 * t979;
t841 = t1246 * t979 - t1248 * t863;
t840 = t1246 * t859 + t1248 * t974;
t839 = t1246 * t974 - t1248 * t859;
t838 = t1404 * t965 + t1276;
t837 = -pkin(1) * t1006 + qJ(2) * t961 + t1246 * t924 + t1248 * t900;
t836 = t1256 * t886 - t1394;
t833 = t1246 * t855 + t1248 * t927;
t832 = t1246 * t927 - t1248 * t855;
t831 = t1246 * t850 + t1248 * t928;
t830 = t1246 * t928 - t1248 * t850;
t829 = pkin(1) * t944 + t1267;
t828 = t1044 * t1248 + t1246 * t835;
t827 = t1044 * t1246 - t1248 * t835;
t822 = pkin(1) * t932 + t1268;
t813 = -qJ(2) * t944 - t1246 * t844 + t1248 * t865;
t810 = -pkin(1) * t994 + qJ(2) * t947 + t1246 * t865 + t1248 * t844;
t809 = t1404 * t937 + t1277;
t808 = -qJ(2) * t932 - t1246 * t838 + t1248 * t852;
t805 = -pkin(1) * t965 + qJ(2) * t933 + t1246 * t852 + t1248 * t838;
t804 = pkin(1) * t901 + t1269;
t803 = -qJ(3) * t864 - t1287;
t801 = -qJ(3) * t860 - t1288;
t800 = -qJ(3) * t836 - t1339;
t799 = -t1253 * t819 + t1256 * t820;
t796 = -qJ(2) * t901 - t1246 * t809 + t1248 * t861;
t795 = t1404 * t836 + t1264;
t794 = -pkin(1) * t937 + qJ(2) * t902 + t1246 * t861 + t1248 * t809;
t793 = -qJ(3) * t851 - t1283;
t790 = t1246 * t798 + t1248 * t825;
t789 = t1246 * t825 - t1248 * t798;
t788 = pkin(1) * t827 + t1261;
t786 = t1404 * t864 + t1278;
t785 = t1404 * t860 + t1279;
t784 = pkin(1) * t841 + t1271;
t783 = pkin(1) * t839 + t1272;
t782 = t1404 * t851 + t1280;
t781 = -qJ(2) * t827 - t1246 * t795 + t1248 * t800;
t780 = pkin(1) * t830 + t1273;
t779 = -pkin(1) * t836 + qJ(2) * t828 + t1246 * t800 + t1248 * t795;
t778 = -qJ(2) * t841 - t1246 * t786 + t1248 * t803;
t777 = -qJ(2) * t839 - t1246 * t785 + t1248 * t801;
t776 = -qJ(3) * t799 - t1303;
t775 = -pkin(1) * t864 + qJ(2) * t842 + t1246 * t803 + t1248 * t786;
t774 = -pkin(1) * t860 + qJ(2) * t840 + t1246 * t801 + t1248 * t785;
t773 = -qJ(2) * t830 - t1246 * t782 + t1248 * t793;
t772 = -pkin(1) * t851 + qJ(2) * t831 + t1246 * t793 + t1248 * t782;
t771 = t1404 * t799 + t1281;
t770 = pkin(1) * t789 + t1274;
t769 = -qJ(2) * t789 - t1246 * t771 + t1248 * t776;
t768 = -pkin(1) * t799 + qJ(2) * t790 + t1246 * t776 + t1248 * t771;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1209, 0, -t1208, 0, t1310, -t1181, -t1302, -pkin(6) * t1302, 0, 0, -t1320, 0, -t1322, 0, t1406, -t1428, t1424, pkin(6) * t1424 + qJ(2) * t1380 - t1254 * t1046, 0, t1320, t1322, 0, 0, 0, -t1417, -t1406, t1428, -pkin(6) * t1417 - t1254 * t1009 + t1257 * t1018, -t1115 * t1254 + t1117 * t1257, -t1087 * t1254 + t1088 * t1257, -t1112 * t1254 + t1114 * t1257, -t1116 * t1254 + t1118 * t1257, -t1111 * t1254 + t1113 * t1257, -t1161 * t1254 + t1162 * t1257, t1257 * t906 - t1254 * t896 - pkin(6) * (t1098 * t1257 + t1100 * t1254), t1257 * t907 - t1254 * t897 - pkin(6) * (t1099 * t1257 + t1101 * t1254), t1257 * t954 - t1254 * t953 - pkin(6) * (t1142 * t1257 + t1143 * t1254), t1257 * t843 - t1254 * t837 - pkin(6) * (t1254 * t961 + t1257 * t960), -t1254 * t968 + t1257 * t970, -t1254 * t913 + t1257 * t914, -t1254 * t945 + t1257 * t948, -t1254 * t967 + t1257 * t969, -t1254 * t943 + t1257 * t946, -t1026 * t1254 + t1027 * t1257, t1257 * t808 - t1254 * t805 - pkin(6) * (t1254 * t933 + t1257 * t932), t1257 * t813 - t1254 * t810 - pkin(6) * (t1254 * t947 + t1257 * t944), t1257 * t796 - t1254 * t794 - pkin(6) * (t1254 * t902 + t1257 * t901), t1257 * t781 - t1254 * t779 - pkin(6) * (t1254 * t828 + t1257 * t827), -t1254 * t874 + t1257 * t876, -t1254 * t832 + t1257 * t833, -t1254 * t845 + t1257 * t847, -t1254 * t873 + t1257 * t875, -t1254 * t846 + t1257 * t848, -t1254 * t887 + t1257 * t888, t1257 * t777 - t1254 * t774 - pkin(6) * (t1254 * t840 + t1257 * t839), t1257 * t778 - t1254 * t775 - pkin(6) * (t1254 * t842 + t1257 * t841), t1257 * t773 - t1254 * t772 - pkin(6) * (t1254 * t831 + t1257 * t830), t1257 * t769 - t1254 * t768 - pkin(6) * (t1254 * t790 + t1257 * t789); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1208, 0, t1209, 0, t1181, t1310, t1318, pkin(6) * t1318, 0, 0, t1322, 0, -t1320, 0, t1428, t1406, t1425, pkin(6) * t1425 + qJ(2) * t1381 + t1257 * t1046, 0, -t1322, t1320, 0, 0, 0, t1418, -t1428, -t1406, pkin(6) * t1418 + t1257 * t1009 + t1254 * t1018, t1115 * t1257 + t1117 * t1254, t1087 * t1257 + t1088 * t1254, t1112 * t1257 + t1114 * t1254, t1116 * t1257 + t1118 * t1254, t1111 * t1257 + t1113 * t1254, t1161 * t1257 + t1162 * t1254, t1254 * t906 + t1257 * t896 + pkin(6) * (-t1098 * t1254 + t1100 * t1257), t1254 * t907 + t1257 * t897 + pkin(6) * (-t1099 * t1254 + t1101 * t1257), t1254 * t954 + t1257 * t953 + pkin(6) * (-t1142 * t1254 + t1143 * t1257), t1254 * t843 + t1257 * t837 + pkin(6) * (-t1254 * t960 + t1257 * t961), t1254 * t970 + t1257 * t968, t1254 * t914 + t1257 * t913, t1254 * t948 + t1257 * t945, t1254 * t969 + t1257 * t967, t1254 * t946 + t1257 * t943, t1026 * t1257 + t1027 * t1254, t1254 * t808 + t1257 * t805 + pkin(6) * (-t1254 * t932 + t1257 * t933), t1254 * t813 + t1257 * t810 + pkin(6) * (-t1254 * t944 + t1257 * t947), t1254 * t796 + t1257 * t794 + pkin(6) * (-t1254 * t901 + t1257 * t902), t1254 * t781 + t1257 * t779 + pkin(6) * (-t1254 * t827 + t1257 * t828), t1254 * t876 + t1257 * t874, t1254 * t833 + t1257 * t832, t1254 * t847 + t1257 * t845, t1254 * t875 + t1257 * t873, t1254 * t848 + t1257 * t846, t1254 * t888 + t1257 * t887, t1254 * t777 + t1257 * t774 + pkin(6) * (-t1254 * t839 + t1257 * t840), t1254 * t778 + t1257 * t775 + pkin(6) * (-t1254 * t841 + t1257 * t842), t1254 * t773 + t1257 * t772 + pkin(6) * (-t1254 * t830 + t1257 * t831), t1254 * t769 + t1257 * t768 + pkin(6) * (-t1254 * t789 + t1257 * t790); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1215, t1216, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1095, t1094, 0, -t1051, qJDD(1), 0, 0, 0, 0, 0, 0, t1081, t1075, t952, t1158, t1133, t1159, t1151, t1156, 0, t972, t973, t951, t877, t999, t936, t993, t998, t992, t1029, t822, t829, t804, t788, t892, t856, t869, t891, t870, t909, t783, t784, t780, t770; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1259, 0, 0, -g(3), -t1215, 0, 0, 0, -t1206, 0, -t1205, 0, t1166, -t1329, t1054, qJ(2) * t1054, 0, t1206, t1205, 0, 0, 0, -t1033, -t1166, t1329, t1018, t1117, t1088, t1114, t1118, t1113, t1162, t906, t907, t954, t843, t970, t914, t948, t969, t946, t1027, t808, t813, t796, t781, t876, t833, t847, t875, t848, t888, t777, t778, t773, t769; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1259, 0, qJDD(1), 0, g(3), 0, -t1216, 0, 0, 0, t1205, 0, -t1206, 0, t1329, t1166, t1325, t1046, 0, -t1205, t1206, 0, 0, 0, t1326, -t1329, -t1166, t1009, t1115, t1087, t1112, t1116, t1111, t1161, t896, t897, t953, t837, t968, t913, t945, t967, t943, t1026, t805, t810, t794, t779, t874, t832, t845, t873, t846, t887, t774, t775, t772, t768; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1215, t1216, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1095, t1094, 0, -t1051, qJDD(1), 0, 0, 0, 0, 0, 0, t1081, t1075, t952, t1158, t1133, t1159, t1151, t1156, 0, t972, t973, t951, t877, t999, t936, t993, t998, t992, t1029, t822, t829, t804, t788, t892, t856, t869, t891, t870, t909, t783, t784, t780, t770; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1259, 0, 0, -t1242, t1128, 0, 0, -qJDD(1), t1259, 0, 0, 0, t1110, 0, t1242, qJ(3) * t1242, t1222, t1211, t1233, -t1222, -t1346, qJDD(4), t1024, t1025, -t1399, t924, t1368, t1125, t1086, -t1368, t1084, qJDD(4), t852, t865, t861, t800, t1002, t927, t982, -t1000, t983, t1030, t801, t803, t793, t776; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1259, 0, qJDD(1), 0, t1242, 0, t1129, 0, 0, -t1259, -qJDD(1), 0, 0, 0, t1108, -t1242, 0, pkin(2) * t1242, -t1164, -t1132, -t1154, -t1163, -t1152, t1196, t1008, t1007, -t995, t900, -t997, -t934, -t990, -t996, -t989, -t1028, t838, t844, t809, t795, -t890, -t855, -t867, -t889, -t868, -t908, t785, t786, t782, t771; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1128, -t1129, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t1289, t1262, t1348, t1158, t1133, t1159, t1151, t1156, 0, t1282, t1298, t1265, t1311, t999, t936, t993, t998, t992, t1029, t1268, t1267, t1269, t1261, t892, t856, t869, t891, t870, t909, t1272, t1271, t1273, t1274; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t1110, t1108, 0, t1158, t1133, t1159, t1151, t1156, 0, t1305, t1330, t1290, -t1397, t999, t936, t993, t998, t992, t1029, t1292, t1291, t1293, t1266, t892, t856, t869, t891, t870, t909, t1295, t1294, t1296, t1297; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1259, 0, 0, 0, -t1110, 0, -t1242, 0, -t1222, -t1211, -t1233, t1222, t1346, -qJDD(4), t1270, t1312, t1399, -t1401, -t1368, -t1125, -t1086, t1368, -t1084, -qJDD(4), t1284 + t1344, t1174 + t1300, t1331, t1339, -t1002, -t927, -t982, t1000, -t983, -t1030, t1288, t1287, t1283, t1303; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1259, qJDD(1), 0, 0, 0, -t1108, t1242, 0, 0, t1164, t1132, t1154, t1163, t1152, -t1196, pkin(7) * t1157 - t1307, pkin(7) * t1160 - t1306, t995, pkin(7) * t1006 - t1400, t997, t934, t990, t996, t989, t1028, pkin(7) * t965 - t1276, pkin(7) * t994 - t1275, pkin(7) * t937 - t1277, pkin(7) * t836 - t1264, t890, t855, t867, t889, t868, t908, pkin(7) * t860 - t1279, pkin(7) * t864 - t1278, pkin(7) * t851 - t1280, pkin(7) * t799 - t1281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1203, -t1202, t1214, t1334, t1218, -t1334, 0, t1093, t1072, 0, t1071, t1022, t1064, t1069, t1063, t1105, t963, t981, t871, -qJ(5) * t885, t958, t899, t922, t957, t923, t988, t815, t817, t807, t792; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1333, t1204, t1220, t1285, t1213, -t1333, -t1093, 0, t1073, 0, t1070, t1020, t1061, t1068, t1060, t1104, t938, t940, t858, t872, t956, t898, t920, t955, t921, t987, t811, t812, t802, t787; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1222, t1211, t1233, -t1222, -t1346, qJDD(4), -t1072, -t1073, 0, 0, t1368, t1125, t1086, -t1368, t1084, qJDD(4), t1176 + t1308, t1175 + t1332, t1019, t884, t1002, t927, t982, -t1000, t983, t1030, t1315, t1316, t1309, t1343; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1131, -t1082, t1423, t1389, t1171, -t1389, 0, t1044, t941, 0, t1003, t929, t984, t1001, t985, t1031, t882, t883, t821, -pkin(8) * t825; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1178, t1085, t1172, -t1323, t1123, -t1178, -t1044, 0, t942, 0, -t1080, -t1076, -t1016, t1080, t1012, -t1126, t853, t854, -pkin(5) * t928, -pkin(5) * t825; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1368, t1125, t1086, -t1368, t1084, qJDD(4), -t941, -t942, 0, 0, t1002, t927, t982, -t1000, t983, t1030, t1341, t1342, t1317, t1396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1058, t1013, t1410, t1106, t1096, -t1106, 0, t918, t880, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1107, t1015, t1097, -t1286, t1041, -t1107, -t918, 0, t881, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1080, t1076, t1016, -t1080, -t1012, t1126, -t880, -t881, 0, 0;];
m_new_reg  = t1;