% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6RRPPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 08:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6RRPPPR3_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR3_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR3_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:35:57
% EndTime: 2019-05-06 08:36:18
% DurationCPUTime: 22.36s
% Computational Cost: add. (84889->834), mult. (188264->914), div. (0->0), fcn. (105548->8), ass. (0->493)
t1275 = sin(qJ(2));
t1268 = t1275 ^ 2;
t1280 = qJD(1) ^ 2;
t1388 = t1268 * t1280;
t1462 = qJD(2) ^ 2;
t1244 = t1388 + t1462;
t1278 = cos(qJ(2));
t1483 = t1278 * t1280;
t1362 = t1275 * t1483;
t1239 = -qJDD(2) + t1362;
t1384 = t1278 * t1239;
t1179 = t1244 * t1275 + t1384;
t1416 = qJD(2) * t1278;
t1257 = qJD(1) * t1416;
t1375 = t1275 * qJDD(1);
t1226 = 0.2e1 * t1257 + t1375;
t1276 = sin(qJ(1));
t1279 = cos(qJ(1));
t1114 = pkin(6) * (t1179 * t1279 + t1226 * t1276);
t1451 = pkin(6) * (t1179 * t1276 - t1226 * t1279);
t1238 = qJDD(2) + t1362;
t1219 = t1278 * t1238;
t1459 = t1278 ^ 2;
t1383 = t1459 * t1280;
t1246 = t1383 + t1462;
t1177 = -t1246 * t1275 + t1219;
t1455 = pkin(1) * t1177;
t1395 = t1238 * t1275;
t1181 = t1246 * t1278 + t1395;
t1419 = qJD(1) * t1275;
t1256 = qJD(2) * t1419;
t1374 = t1278 * qJDD(1);
t1229 = -0.2e1 * t1256 + t1374;
t1115 = pkin(6) * (t1181 * t1279 + t1229 * t1276);
t1450 = pkin(6) * (t1181 * t1276 - t1229 * t1279);
t1447 = pkin(7) * t1177;
t1271 = sin(pkin(9));
t1272 = cos(pkin(9));
t1418 = qJD(1) * t1278;
t1211 = -qJD(2) * t1272 + t1271 * t1418;
t1461 = t1211 ^ 2;
t1154 = -t1388 - t1461;
t1212 = qJD(2) * t1271 + t1272 * t1418;
t1162 = t1211 * t1212;
t1227 = t1257 + t1375;
t1479 = -t1162 + t1227;
t1494 = t1272 * t1479;
t1074 = t1271 * t1154 + t1494;
t1228 = -t1256 + t1374;
t1240 = t1276 * g(1) - t1279 * g(2);
t1205 = qJDD(1) * pkin(1) + t1280 * pkin(7) + t1240;
t1292 = -pkin(2) * t1256 + t1205;
t1474 = t1228 * pkin(3) - qJ(4) * t1383 + qJDD(4);
t1288 = t1292 + t1474;
t1237 = -qJD(2) * pkin(3) - qJ(4) * t1419;
t1382 = (2 * qJD(3)) + t1237;
t1350 = t1382 * t1275;
t1456 = qJ(3) + pkin(4);
t1000 = (pkin(2) + qJ(5)) * t1228 + t1456 * t1227 + (t1350 + (-qJ(5) * t1275 + t1278 * t1456) * qJD(2)) * qJD(1) + t1288;
t1445 = t1278 * g(3);
t1318 = -qJDD(2) * pkin(2) - t1462 * qJ(3) + qJDD(3) + t1445;
t1241 = g(1) * t1279 + g(2) * t1276;
t1206 = -pkin(1) * t1280 + qJDD(1) * pkin(7) - t1241;
t1385 = t1275 * t1206;
t1496 = pkin(3) * t1238;
t1289 = -t1227 * qJ(4) + t1318 + t1385 - t1496;
t1330 = pkin(4) * t1275 + qJ(5) * t1278;
t1370 = qJ(4) * t1416;
t1331 = -pkin(2) * t1278 - qJ(3) * t1275;
t1224 = t1331 * qJD(1);
t1381 = -(2 * qJD(4)) + t1224;
t1044 = -t1462 * pkin(4) - qJDD(2) * qJ(5) + (t1370 + (-qJD(1) * t1330 + t1381) * t1275) * qJD(1) + t1289;
t1498 = 2 * qJD(5);
t1369 = t1272 * t1000 - t1271 * t1044 + t1212 * t1498;
t1501 = -pkin(4) * t1074 - t1369;
t1380 = pkin(1) * t1226 - pkin(7) * t1179;
t1394 = t1239 * t1275;
t1171 = -t1244 * t1278 + t1394;
t1497 = pkin(1) * t1171;
t1448 = pkin(7) * t1171;
t1495 = t1271 * t1479;
t1274 = sin(qJ(6));
t1277 = cos(qJ(6));
t1151 = -t1277 * t1211 - t1212 * t1274;
t1153 = t1211 * t1274 - t1212 * t1277;
t1082 = t1153 * t1151;
t1214 = qJDD(6) + t1227;
t1481 = -t1082 + t1214;
t1493 = t1274 * t1481;
t1492 = t1277 * t1481;
t1264 = t1275 * g(3);
t1376 = t1462 * pkin(2) + t1264;
t1412 = qJDD(2) * qJ(3);
t1313 = -0.2e1 * qJD(3) * qJD(2) + t1376 - t1412;
t1347 = qJD(1) * t1224 + t1206;
t1101 = t1278 * t1347 - t1313;
t1478 = pkin(2) * t1244 - qJ(3) * t1239;
t1491 = t1101 + t1478;
t1105 = t1275 * t1347 + t1318;
t1439 = qJ(3) * t1246;
t1490 = pkin(2) * t1238 - t1105 - t1439;
t1379 = pkin(1) * t1229 - pkin(7) * t1181;
t1399 = t1229 * t1278;
t1405 = t1226 * t1275;
t1158 = -t1399 + t1405;
t1236 = (t1268 - t1459) * t1280;
t1489 = t1158 * t1276 + t1236 * t1279;
t1488 = t1158 * t1279 - t1236 * t1276;
t1247 = t1383 - t1462;
t1182 = t1247 * t1278 + t1394;
t1142 = t1182 * t1276 - t1279 * t1374;
t1144 = t1182 * t1279 + t1276 * t1374;
t1371 = t1268 + t1459;
t1232 = t1371 * qJDD(1);
t1235 = t1371 * t1280;
t1159 = pkin(6) * (t1232 * t1279 - t1235 * t1276);
t1460 = t1212 ^ 2;
t1339 = -t1388 - t1460;
t1139 = t1162 + t1227;
t1387 = t1271 * t1139;
t1086 = t1272 * t1339 - t1387;
t962 = t1271 * t1000 + t1272 * t1044 + t1211 * t1498;
t1485 = -pkin(4) * t1086 + t962;
t902 = -t1271 * t1369 + t1272 * t962;
t1484 = t1271 * t962 + t1272 * t1369;
t1194 = -t1271 * qJDD(2) - t1228 * t1272;
t1319 = qJDD(2) * t1272 - t1228 * t1271;
t1064 = -t1151 * qJD(6) + t1277 * t1194 - t1274 * t1319;
t1250 = qJD(6) + t1419;
t1129 = t1250 * t1151;
t1480 = -t1129 + t1064;
t1364 = t1211 * t1419;
t1137 = t1194 - t1364;
t1378 = pkin(1) * t1235 + pkin(7) * t1232;
t1477 = pkin(3) * t1383 + t1228 * qJ(4);
t1176 = -t1247 * t1275 + t1384;
t1476 = pkin(3) * t1229 + qJ(4) * t1246;
t1329 = t1227 + t1257;
t1372 = 0.2e1 * t1419;
t1349 = qJD(3) * t1372;
t1284 = qJ(3) * t1329 + t1292 + t1349;
t1079 = (t1228 + t1229) * pkin(2) + t1284;
t1149 = t1151 ^ 2;
t1150 = t1153 ^ 2;
t1049 = -t1149 - t1150;
t931 = pkin(5) * t1479 - pkin(8) * t1137 + t1369;
t1195 = pkin(5) * t1419 + pkin(8) * t1212;
t940 = -pkin(5) * t1461 - pkin(8) * t1319 - t1195 * t1419 + t962;
t889 = t1274 * t940 - t1277 * t931;
t890 = t1274 * t931 + t1277 * t940;
t878 = t1274 * t889 + t1277 * t890;
t1346 = -t1274 * t1194 - t1277 * t1319;
t1019 = (-qJD(6) + t1250) * t1153 + t1346;
t1022 = t1129 + t1064;
t969 = t1019 * t1277 + t1022 * t1274;
t866 = -pkin(5) * t1049 + pkin(8) * t969 + t878;
t877 = t1274 * t890 - t1277 * t889;
t967 = t1019 * t1274 - t1022 * t1277;
t870 = -pkin(8) * t967 - t877;
t908 = -t1271 * t967 + t1272 * t969;
t1344 = -pkin(4) * t1049 + qJ(5) * t908 + t1271 * t870 + t1272 * t866;
t1458 = pkin(2) + pkin(3);
t1473 = qJ(3) * t1049 - t1458 * t908 - t1344;
t1290 = qJD(2) * t1382 - t1376 - t1477;
t1326 = qJD(1) * t1381 + t1206;
t1286 = t1278 * t1326 + t1290;
t1043 = -t1462 * qJ(5) + qJDD(2) * t1456 - t1330 * t1483 + qJDD(5) + t1286;
t984 = pkin(5) * t1319 - pkin(8) * t1461 - t1195 * t1212 + t1043;
t1426 = t1274 * t984;
t1248 = t1250 ^ 2;
t1102 = -t1150 - t1248;
t1066 = t1082 + t1214;
t1410 = t1066 * t1277;
t999 = -t1102 * t1274 - t1410;
t917 = -pkin(5) * t1480 + pkin(8) * t999 + t1426;
t1425 = t1277 * t984;
t1411 = t1066 * t1274;
t998 = t1102 * t1277 - t1411;
t943 = -pkin(8) * t998 + t1425;
t949 = -t1271 * t998 + t1272 * t999;
t1342 = -pkin(4) * t1480 + qJ(5) * t949 + t1271 * t943 + t1272 * t917;
t1472 = qJ(3) * t1480 - t1458 * t949 - t1342;
t1130 = t1250 * t1153;
t1304 = qJD(6) * t1153 - t1346;
t1017 = t1130 + t1304;
t1077 = -t1248 - t1149;
t990 = t1077 * t1277 - t1493;
t912 = -pkin(5) * t1017 + pkin(8) * t990 - t1425;
t989 = t1077 * t1274 + t1492;
t930 = -pkin(8) * t989 + t1426;
t939 = -t1271 * t989 + t1272 * t990;
t1343 = -pkin(4) * t1017 + qJ(5) * t939 + t1271 * t930 + t1272 * t912;
t1471 = qJ(3) * t1017 - t1458 * t939 - t1343;
t1430 = t1271 * t877;
t856 = t1272 * t878 - t1430;
t873 = -pkin(5) * t984 + pkin(8) * t878;
t1348 = pkin(4) * t984 + pkin(8) * t1430 - qJ(5) * t856 - t1272 * t873;
t1470 = qJ(3) * t984 - t1458 * t856 + t1348;
t1422 = pkin(4) * t1043 - qJ(5) * t902;
t1469 = qJ(3) * t1043 - t1458 * t902 + t1422;
t1406 = t1212 * t1275;
t1363 = qJD(1) * t1406;
t1134 = t1319 + t1363;
t1060 = -t1134 * t1272 + t1271 * t1137;
t1127 = -t1460 - t1461;
t1317 = -pkin(4) * t1127 + qJ(5) * t1060 + t902;
t1468 = qJ(3) * t1127 - t1060 * t1458 - t1317;
t1057 = t1286 + t1412;
t1062 = (t1275 * t1381 + t1370) * qJD(1) + t1289;
t1467 = qJ(3) * t1057 - t1062 * t1458;
t1075 = t1154 * t1272 - t1495;
t1133 = t1319 - t1363;
t1386 = t1272 * t1043;
t1360 = pkin(4) * t1133 - qJ(5) * t1075 + t1386;
t1466 = qJ(3) * t1133 - t1075 * t1458 + t1360;
t1409 = t1139 * t1272;
t1089 = -t1271 * t1339 - t1409;
t1136 = t1194 + t1364;
t1031 = t1271 * t1043;
t1361 = -pkin(4) * t1136 + qJ(5) * t1089 + t1031;
t1465 = qJ(3) * t1136 - t1089 * t1458 - t1361;
t1464 = -t1238 * t1458 + t1439;
t1463 = -qJD(2) * t1237 + 0.2e1 * qJD(4) * t1418 + t1477;
t1457 = -pkin(3) - qJ(5);
t1452 = pkin(3) * t1244;
t1449 = pkin(6) * (t1232 * t1276 + t1235 * t1279);
t1446 = t1228 * pkin(2);
t1428 = t1272 * t877;
t855 = t1271 * t878 + t1428;
t876 = pkin(5) * t877;
t1444 = -pkin(4) * t855 - t876;
t906 = t1271 * t969 + t1272 * t967;
t966 = pkin(5) * t967;
t1443 = -pkin(4) * t906 - t966;
t1441 = qJ(3) * t1229;
t1440 = qJ(3) * t1235;
t1438 = qJ(4) * t1043;
t1437 = qJ(4) * t1057;
t1436 = qJ(4) * t1062;
t1435 = qJ(4) * t1238;
t1434 = qJ(4) * t1239;
t1433 = qJ(4) * t1244;
t1421 = qJ(4) * qJDD(1);
t1420 = qJD(1) * qJD(2);
t1408 = t1205 * t1275;
t1407 = t1205 * t1278;
t1403 = t1226 * t1278;
t1401 = t1229 * t1275;
t1390 = t1250 * t1274;
t1389 = t1250 * t1277;
t1373 = pkin(2) - t1457;
t1368 = t1275 * t1082;
t1367 = t1278 * t1082;
t1366 = t1211 * t1406;
t1365 = t1278 * t1162;
t900 = pkin(4) * t1484;
t1357 = -qJ(4) * t902 + t900;
t1356 = t1271 * t866 - t1272 * t870;
t1355 = t1271 * t912 - t1272 * t930;
t1354 = t1271 * t917 - t1272 * t943;
t1352 = -qJ(4) * t1136 - t1386;
t1058 = -t1134 * t1271 - t1137 * t1272;
t1056 = pkin(4) * t1058;
t1351 = -qJ(4) * t1060 + t1056;
t1185 = t1385 + t1445;
t1186 = t1278 * t1206 - t1264;
t1100 = t1185 * t1275 + t1278 * t1186;
t1345 = -t1240 * t1276 - t1279 * t1241;
t1341 = t1276 * t1362;
t1340 = t1279 * t1362;
t1338 = t1388 - t1460;
t1337 = -qJ(4) * t856 - t1444;
t1336 = -qJ(4) * t908 - t1443;
t1335 = pkin(5) * t989 - t889;
t1234 = qJDD(1) * t1279 - t1276 * t1280;
t1334 = -pkin(6) * t1234 - g(3) * t1276;
t1332 = -pkin(2) * t1105 + qJ(3) * t1101;
t1327 = -qJ(4) * t1133 - t1031;
t1325 = -t1017 * t1274 + t1277 * t1480;
t1120 = -t1150 + t1248;
t1324 = t1120 * t1277 + t1493;
t1119 = t1149 - t1248;
t1323 = t1119 * t1274 + t1410;
t1099 = t1185 * t1278 - t1186 * t1275;
t1155 = t1401 + t1403;
t1321 = t1240 * t1279 - t1241 * t1276;
t938 = t1271 * t990 + t1272 * t989;
t1320 = -pkin(4) * t938 - t1335;
t1316 = pkin(8) * t1428 + t1271 * t873;
t1312 = -qJ(4) * t1017 + t1355;
t1311 = -qJ(4) * t1480 + t1354;
t1310 = -qJ(4) * t1049 + t1356;
t1309 = -qJ(4) * t1127 + t1484;
t1308 = pkin(5) * t998 - t890;
t1306 = t1151 * t1390 - t1277 * t1304;
t1305 = t1064 * t1274 + t1153 * t1389;
t1299 = -qJ(4) * t1075 - t1501;
t1298 = -qJ(4) * t1089 - t1485;
t1297 = (-t1151 * t1274 - t1153 * t1277) * t1250;
t948 = t1271 * t999 + t1272 * t998;
t1296 = -pkin(4) * t948 - t1308;
t1295 = -qJ(4) * t939 - t1320;
t1294 = qJ(4) * t984 - t1316;
t1291 = -qJ(4) * t949 - t1296;
t1051 = -t1491 + t1497;
t1285 = t1446 + (t1226 + t1329) * qJ(3) + t1292;
t1283 = qJD(4) * t1372 + (-t1275 * t1224 - t1370) * qJD(1) - t1289;
t1282 = qJ(4) * t1375 + t1283;
t1281 = t1057 + t1452;
t1050 = t1446 + t1227 * qJ(3) + (qJ(3) * t1416 + t1350) * qJD(1) + t1288;
t1259 = qJ(3) * t1374;
t1245 = t1388 - t1462;
t1233 = qJDD(1) * t1276 + t1279 * t1280;
t1222 = pkin(2) * t1375 - t1259;
t1218 = t1371 * t1420;
t1208 = t1278 * t1227;
t1207 = t1275 * t1227;
t1199 = -pkin(6) * t1233 + g(3) * t1279;
t1197 = t1375 * t1458 - t1259;
t1196 = -t1388 + t1461;
t1193 = qJDD(2) * t1276 + t1218 * t1279;
t1192 = -t1268 * t1420 + t1208;
t1191 = -qJDD(2) * t1279 + t1218 * t1276;
t1190 = -t1228 * t1275 - t1420 * t1459;
t1178 = t1245 * t1275 + t1219;
t1173 = t1257 * t1275 + t1207;
t1170 = -t1245 * t1278 + t1395;
t1169 = (t1228 - t1256) * t1278;
t1166 = -t1435 - t1441;
t1160 = t1460 - t1461;
t1148 = t1192 * t1279 - t1341;
t1147 = t1190 * t1279 + t1341;
t1146 = t1192 * t1276 + t1340;
t1145 = t1190 * t1276 - t1340;
t1143 = t1178 * t1279 + t1276 * t1375;
t1141 = t1178 * t1276 - t1279 * t1375;
t1140 = t1226 * t1458 + t1434;
t1124 = (t1211 * t1272 - t1212 * t1271) * t1419;
t1123 = (t1211 * t1271 + t1212 * t1272) * t1419;
t1113 = -t1407 - t1448;
t1112 = -t1408 - t1447;
t1111 = t1194 * t1272 + t1271 * t1363;
t1110 = t1194 * t1271 - t1272 * t1363;
t1109 = -t1271 * t1319 + t1272 * t1364;
t1108 = t1271 * t1364 + t1272 * t1319;
t1104 = t1185 - t1455;
t1103 = t1186 - t1497;
t1097 = t1379 + t1407;
t1096 = -t1380 - t1408;
t1095 = t1105 + t1440;
t1094 = t1124 * t1275 + t1208;
t1093 = -t1124 * t1278 + t1207;
t1092 = pkin(2) * t1235 + t1101;
t1091 = t1284 + t1446;
t1090 = t1196 * t1272 - t1387;
t1088 = -t1271 * t1338 + t1494;
t1087 = t1196 * t1271 + t1409;
t1085 = t1272 * t1338 + t1495;
t1081 = t1150 - t1149;
t1080 = pkin(1) * t1205 + pkin(7) * t1100;
t1078 = t1285 + t1349;
t1076 = t1100 + t1378;
t1073 = t1111 * t1275 + t1365;
t1072 = -t1109 * t1275 - t1365;
t1071 = -t1111 * t1278 + t1366;
t1070 = t1109 * t1278 - t1366;
t1061 = -t1133 * t1272 - t1271 * t1136;
t1059 = -t1133 * t1271 + t1136 * t1272;
t1054 = (-t1151 * t1277 + t1153 * t1274) * t1250;
t1052 = -t1455 - t1490;
t1048 = t1282 - t1440;
t1047 = t1101 * t1278 + t1105 * t1275;
t1046 = t1101 * t1275 - t1105 * t1278;
t1042 = t1088 * t1275 + t1137 * t1278;
t1041 = t1089 * t1275 + t1136 * t1278;
t1040 = t1090 * t1275 - t1134 * t1278;
t1039 = -t1088 * t1278 + t1137 * t1275;
t1038 = -t1089 * t1278 + t1136 * t1275;
t1037 = -t1090 * t1278 - t1134 * t1275;
t1030 = -t1458 * t1235 + (-t1347 + t1421) * t1278 + t1313 + t1463;
t1029 = qJD(1) * t1350 + t1285 + t1433 + t1474;
t1028 = -pkin(2) * t1405 + t1078 * t1278 + t1448;
t1027 = qJ(3) * t1399 - t1079 * t1275 - t1447;
t1026 = t1061 * t1275 + t1160 * t1278;
t1025 = -t1061 * t1278 + t1160 * t1275;
t1024 = t1075 * t1275 + t1133 * t1278;
t1023 = -t1075 * t1278 + t1133 * t1275;
t1018 = -t1130 + t1304;
t1014 = -t1092 * t1275 + t1095 * t1278;
t1013 = t1119 * t1277 - t1411;
t1012 = -t1120 * t1274 + t1492;
t1011 = -t1237 * t1419 - t1079 - t1474 - t1476;
t1010 = t1064 * t1277 - t1153 * t1390;
t1009 = t1151 * t1389 + t1274 * t1304;
t1008 = pkin(2) * t1403 + t1078 * t1275 + t1380;
t1007 = qJ(3) * t1401 + t1079 * t1278 + t1379;
t1004 = t1060 * t1275 + t1127 * t1278;
t1003 = -t1060 * t1278 + t1127 * t1275;
t994 = t1051 - t1452 + t1463;
t992 = t1283 - t1464 + t1455;
t991 = t1092 * t1278 + t1095 * t1275 + t1378;
t986 = t1057 * t1278 + t1062 * t1275;
t985 = t1057 * t1275 - t1062 * t1278;
t983 = t1054 * t1272 - t1271 * t1297;
t982 = t1271 * t1054 + t1272 * t1297;
t980 = qJ(3) * t1050 - t1436;
t979 = -t1011 * t1275 + t1166 * t1278 + t1447;
t978 = t1029 * t1278 - t1140 * t1275 + t1448;
t977 = t1011 * t1278 + t1166 * t1275 - t1379;
t976 = t1029 * t1275 + t1140 * t1278 + t1380;
t975 = t1214 * t1278 + t1275 * t983;
t974 = t1214 * t1275 - t1278 * t983;
t973 = -pkin(1) * t1046 - t1332;
t972 = -t1030 * t1275 + t1048 * t1278;
t971 = -pkin(7) * t1046 + (-pkin(2) * t1275 + qJ(3) * t1278) * t1091;
t970 = t1030 * t1278 + t1048 * t1275 - t1378;
t968 = -t1017 * t1277 - t1274 * t1480;
t963 = qJ(3) * t1058 + t1351;
t960 = t1013 * t1272 - t1271 * t1323;
t959 = t1012 * t1272 - t1271 * t1324;
t958 = t1271 * t1013 + t1272 * t1323;
t957 = t1271 * t1012 + t1272 * t1324;
t953 = t1010 * t1272 - t1271 * t1305;
t952 = t1009 * t1272 - t1271 * t1306;
t951 = t1271 * t1010 + t1272 * t1305;
t950 = t1271 * t1009 + t1272 * t1306;
t947 = t1050 * t1458 - t1437;
t944 = pkin(7) * t1047 + (pkin(1) - t1331) * t1091;
t935 = t1275 * t953 + t1367;
t934 = t1275 * t952 - t1367;
t933 = -t1278 * t953 + t1368;
t932 = -t1278 * t952 - t1368;
t926 = t1086 * t1373 + t1352;
t925 = t1074 * t1373 + t1327;
t924 = t1022 * t1278 + t1275 * t959;
t923 = -t1018 * t1278 + t1275 * t960;
t922 = t1022 * t1275 - t1278 * t959;
t921 = -t1018 * t1275 - t1278 * t960;
t920 = t1275 * t949 + t1278 * t1480;
t919 = t1275 * t1480 - t1278 * t949;
t918 = -pkin(1) * t985 - t1467;
t915 = t1017 * t1278 + t1275 * t939;
t914 = t1017 * t1275 - t1278 * t939;
t913 = qJ(3) * t1086 + t1298;
t910 = qJ(3) * t1074 + t1299;
t909 = -pkin(1) * t1038 - t1465;
t907 = -t1271 * t1325 + t1272 * t968;
t905 = t1271 * t968 + t1272 * t1325;
t898 = -pkin(1) * t1023 - t1466;
t897 = t1081 * t1278 + t1275 * t907;
t896 = t1081 * t1275 - t1278 * t907;
t895 = t1049 * t1278 + t1275 * t908;
t894 = t1049 * t1275 - t1278 * t908;
t893 = -pkin(7) * t985 - t1275 * t947 + t1278 * t980;
t892 = t1043 * t1278 + t1275 * t902;
t891 = t1043 * t1275 - t1278 * t902;
t887 = pkin(1) * t1050 + pkin(7) * t986 + t1275 * t980 + t1278 * t947;
t886 = t1058 * t1373 + t1309;
t885 = -pkin(7) * t1038 - t1275 * t926 + t1278 * t913;
t884 = pkin(1) * t1086 + pkin(7) * t1041 + t1275 * t913 + t1278 * t926;
t883 = -pkin(1) * t1003 - t1468;
t882 = -pkin(7) * t1023 - t1275 * t925 + t1278 * t910;
t881 = pkin(1) * t1074 + pkin(7) * t1024 + t1275 * t910 + t1278 * t925;
t880 = -pkin(7) * t1003 - t1275 * t886 + t1278 * t963;
t879 = pkin(1) * t1058 + pkin(7) * t1004 + t1275 * t963 + t1278 * t886;
t871 = qJ(3) * t1484 + t1357;
t867 = qJ(3) * t948 + t1291;
t865 = t1373 * t1484 - t1438;
t863 = t1373 * t948 + t1311;
t862 = qJ(3) * t906 + t1336;
t861 = qJ(3) * t938 + t1295;
t860 = t1373 * t938 + t1312;
t859 = -pkin(1) * t919 - t1472;
t858 = -pkin(1) * t891 - t1469;
t857 = -pkin(1) * t914 - t1471;
t852 = t1275 * t856 + t1278 * t984;
t851 = t1275 * t984 - t1278 * t856;
t850 = -pkin(7) * t919 - t1275 * t863 + t1278 * t867;
t849 = -pkin(7) * t891 - t1275 * t865 + t1278 * t871;
t848 = pkin(1) * t948 + pkin(7) * t920 + t1275 * t867 + t1278 * t863;
t847 = -pkin(7) * t914 - t1275 * t860 + t1278 * t861;
t846 = pkin(1) * t1484 + pkin(7) * t892 + t1275 * t871 + t1278 * t865;
t845 = pkin(1) * t938 + pkin(7) * t915 + t1275 * t861 + t1278 * t860;
t844 = t1373 * t906 + t1310;
t843 = -pkin(1) * t894 - t1473;
t842 = -pkin(7) * t894 - t1275 * t844 + t1278 * t862;
t841 = pkin(1) * t906 + pkin(7) * t895 + t1275 * t862 + t1278 * t844;
t840 = qJ(3) * t855 + t1337;
t839 = t1373 * t855 - t1294;
t838 = -pkin(1) * t851 - t1470;
t837 = -pkin(7) * t851 - t1275 * t839 + t1278 * t840;
t836 = pkin(1) * t855 + pkin(7) * t852 + t1275 * t840 + t1278 * t839;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1234, 0, -t1233, 0, t1334, -t1199, -t1321, -pkin(6) * t1321, t1148, -t1488, t1143, t1147, t1144, t1193, -t1104 * t1276 + t1112 * t1279 + t1450, -t1276 * t1103 + t1279 * t1113 - t1451, t1099 * t1279 - t1449, -pkin(6) * (t1100 * t1276 + t1205 * t1279) - (pkin(1) * t1276 - pkin(7) * t1279) * t1099, t1148, t1143, t1488, t1193, -t1144, t1147, t1027 * t1279 - t1052 * t1276 + t1450, t1014 * t1279 - t1222 * t1276 - t1449, t1028 * t1279 - t1051 * t1276 + t1451, t1279 * t971 - t1276 * t973 - pkin(6) * (t1047 * t1276 + t1091 * t1279), t1147, -t1488, t1144, t1148, t1143, t1193, -t1276 * t994 + t1279 * t978 + t1451, -t1276 * t992 + t1279 * t979 - t1450, t1276 * t1197 + t1279 * t972 + t1449, t1279 * t893 - t1276 * t918 - pkin(6) * (t1050 * t1279 + t1276 * t986), t1073 * t1279 - t1110 * t1276, t1026 * t1279 - t1059 * t1276, t1042 * t1279 - t1085 * t1276, t1072 * t1279 + t1108 * t1276, t1040 * t1279 - t1087 * t1276, t1094 * t1279 - t1123 * t1276, t1279 * t882 - t1276 * t898 - pkin(6) * (t1024 * t1276 + t1074 * t1279), t1279 * t885 - t1276 * t909 - pkin(6) * (t1041 * t1276 + t1086 * t1279), t1279 * t880 - t1276 * t883 - pkin(6) * (t1004 * t1276 + t1058 * t1279), t1279 * t849 - t1276 * t858 - pkin(6) * (t1276 * t892 + t1279 * t1484), -t1276 * t951 + t1279 * t935, -t1276 * t905 + t1279 * t897, -t1276 * t957 + t1279 * t924, -t1276 * t950 + t1279 * t934, -t1276 * t958 + t1279 * t923, -t1276 * t982 + t1279 * t975, t1279 * t847 - t1276 * t857 - pkin(6) * (t1276 * t915 + t1279 * t938), t1279 * t850 - t1276 * t859 - pkin(6) * (t1276 * t920 + t1279 * t948), t1279 * t842 - t1276 * t843 - pkin(6) * (t1276 * t895 + t1279 * t906), t1279 * t837 - t1276 * t838 - pkin(6) * (t1276 * t852 + t1279 * t855); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1233, 0, t1234, 0, t1199, t1334, t1345, pkin(6) * t1345, t1146, -t1489, t1141, t1145, t1142, t1191, t1104 * t1279 + t1112 * t1276 - t1115, t1279 * t1103 + t1276 * t1113 + t1114, t1099 * t1276 + t1159, pkin(6) * (t1100 * t1279 - t1205 * t1276) - (-pkin(1) * t1279 - pkin(7) * t1276) * t1099, t1146, t1141, t1489, t1191, -t1142, t1145, t1027 * t1276 + t1052 * t1279 - t1115, t1014 * t1276 + t1222 * t1279 + t1159, t1028 * t1276 + t1051 * t1279 - t1114, t1276 * t971 + t1279 * t973 + pkin(6) * (t1047 * t1279 - t1091 * t1276), t1145, -t1489, t1142, t1146, t1141, t1191, t1276 * t978 + t1279 * t994 - t1114, t1276 * t979 + t1279 * t992 + t1115, -t1279 * t1197 + t1276 * t972 - t1159, t1276 * t893 + t1279 * t918 + pkin(6) * (-t1050 * t1276 + t1279 * t986), t1073 * t1276 + t1110 * t1279, t1026 * t1276 + t1059 * t1279, t1042 * t1276 + t1085 * t1279, t1072 * t1276 - t1108 * t1279, t1040 * t1276 + t1087 * t1279, t1094 * t1276 + t1123 * t1279, t1276 * t882 + t1279 * t898 + pkin(6) * (t1024 * t1279 - t1074 * t1276), t1276 * t885 + t1279 * t909 + pkin(6) * (t1041 * t1279 - t1086 * t1276), t1276 * t880 + t1279 * t883 + pkin(6) * (t1004 * t1279 - t1058 * t1276), t1276 * t849 + t1279 * t858 + pkin(6) * (-t1276 * t1484 + t1279 * t892), t1276 * t935 + t1279 * t951, t1276 * t897 + t1279 * t905, t1276 * t924 + t1279 * t957, t1276 * t934 + t1279 * t950, t1276 * t923 + t1279 * t958, t1276 * t975 + t1279 * t982, t1276 * t847 + t1279 * t857 + pkin(6) * (-t1276 * t938 + t1279 * t915), t1276 * t850 + t1279 * t859 + pkin(6) * (-t1276 * t948 + t1279 * t920), t1276 * t842 + t1279 * t843 + pkin(6) * (-t1276 * t906 + t1279 * t895), t1276 * t837 + t1279 * t838 + pkin(6) * (-t1276 * t855 + t1279 * t852); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1240, t1241, 0, 0, t1173, t1155, t1170, t1169, -t1176, 0, t1097, t1096, t1076, t1080, t1173, t1170, -t1155, 0, t1176, t1169, t1007, t991, t1008, t944, t1169, t1155, -t1176, t1173, t1170, 0, t976, t977, t970, t887, t1071, t1025, t1039, t1070, t1037, t1093, t881, t884, t879, t846, t933, t896, t922, t932, t921, t974, t845, t848, t841, t836; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1280, 0, 0, -g(3), -t1240, 0, t1192, -t1158, t1178, t1190, t1182, t1218, t1112, t1113, t1099, pkin(7) * t1099, t1192, t1178, t1158, t1218, -t1182, t1190, t1027, t1014, t1028, t971, t1190, -t1158, t1182, t1192, t1178, t1218, t978, t979, t972, t893, t1073, t1026, t1042, t1072, t1040, t1094, t882, t885, t880, t849, t935, t897, t924, t934, t923, t975, t847, t850, t842, t837; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1280, 0, qJDD(1), 0, g(3), 0, -t1241, 0, t1362, -t1236, -t1375, -t1362, -t1374, -qJDD(2), t1104, t1103, 0, pkin(1) * t1099, t1362, -t1375, t1236, -qJDD(2), t1374, -t1362, t1052, t1222, t1051, t973, -t1362, -t1236, -t1374, t1362, -t1375, -qJDD(2), t994, t992, -t1197, t918, t1110, t1059, t1085, -t1108, t1087, t1123, t898, t909, t883, t858, t951, t905, t957, t950, t958, t982, t857, t859, t843, t838; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1240, t1241, 0, 0, t1173, t1155, t1170, t1169, -t1176, 0, t1097, t1096, t1076, t1080, t1173, t1170, -t1155, 0, t1176, t1169, t1007, t991, t1008, t944, t1169, t1155, -t1176, t1173, t1170, 0, t976, t977, t970, t887, t1071, t1025, t1039, t1070, t1037, t1093, t881, t884, t879, t846, t933, t896, t922, t932, t921, t974, t845, t848, t841, t836; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1227, t1229, t1238, -t1257, t1247, t1257, 0, -t1205, t1185, 0, t1227, t1238, -t1229, t1257, -t1247, -t1257, t1441, t1095, t1078, qJ(3) * t1091, -t1257, t1229, t1247, t1227, t1238, t1257, t1029, t1166, t1048, t980, t1162, t1160, t1137, -t1162, -t1134, t1227, t910, t913, t963, t871, t1082, t1081, t1022, -t1082, -t1018, t1214, t861, t867, t862, t840; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1256, t1226, -t1245, t1228, -t1239, -t1256, t1205, 0, t1186, 0, t1256, -t1245, -t1226, -t1256, t1239, t1228, t1079, t1092, pkin(2) * t1226, pkin(2) * t1091, t1228, t1226, -t1239, t1256, -t1245, -t1256, t1140, t1011, t1030, t947, -t1111, -t1061, -t1088, t1109, -t1090, -t1124, t925, t926, t886, t865, -t953, -t907, -t959, -t952, -t960, -t983, t860, t863, t844, t839; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1362, t1236, t1375, t1362, t1374, qJDD(2), -t1185, -t1186, 0, 0, -t1362, t1375, -t1236, qJDD(2), -t1374, t1362, t1490, -t1222, t1491, t1332, t1362, t1236, t1374, -t1362, t1375, qJDD(2), t1281 + t1478, t1062 + t1464, t1197, t1467, -t1110, -t1059, -t1085, t1108, -t1087, -t1123, t1466, t1465, t1468, t1469, -t951, -t905, -t957, -t950, -t958, -t982, t1471, t1472, t1473, t1470; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1227, t1238, -t1229, t1257, -t1247, -t1257, 0, t1105, t1091, 0, -t1257, t1229, t1247, t1227, t1238, t1257, t1050 + t1433, -t1435, t1282, -t1436, t1162, t1160, t1137, -t1162, -t1134, t1227, t1299, t1298, t1351, t1357, t1082, t1081, t1022, -t1082, -t1018, t1214, t1295, t1291, t1336, t1337; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1362, t1375, -t1236, qJDD(2), -t1374, t1362, -t1105, 0, t1101, 0, t1362, t1236, t1374, -t1362, t1375, qJDD(2), t1281, t1062 - t1496, pkin(3) * t1375, -pkin(3) * t1062, -t1110, -t1059, -t1085, t1108, -t1087, -t1123, -pkin(3) * t1075 + t1360, -pkin(3) * t1089 - t1361, -pkin(3) * t1060 - t1317, -pkin(3) * t902 + t1422, -t951, -t905, -t957, -t950, -t958, -t982, -pkin(3) * t939 - t1343, -pkin(3) * t949 - t1342, -pkin(3) * t908 - t1344, -pkin(3) * t856 + t1348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1256, t1245, t1226, t1256, -t1239, -t1228, -t1091, -t1101, 0, 0, -t1228, -t1226, t1239, -t1256, t1245, t1256, -pkin(3) * t1226 - t1434, t1050 + t1476, pkin(3) * t1235 + t1412 + (t1326 - t1421) * t1278 + t1290, -pkin(3) * t1050 + t1437, t1111, t1061, t1088, -t1109, t1090, t1124, t1074 * t1457 - t1327, t1086 * t1457 - t1352, t1058 * t1457 - t1309, t1457 * t1484 + t1438, t953, t907, t959, t952, t960, t983, t1457 * t938 - t1312, t1457 * t948 - t1311, t1457 * t906 - t1310, t1457 * t855 + t1294; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1228, -t1226, t1239, -t1256, t1245, t1256, 0, t1050, t1057, 0, t1111, t1061, t1088, -t1109, t1090, t1124, -qJ(5) * t1074 + t1031, -qJ(5) * t1086 + t1386, -qJ(5) * t1058 - t1484, -qJ(5) * t1484, t953, t907, t959, t952, t960, t983, -qJ(5) * t938 - t1355, -qJ(5) * t948 - t1354, -qJ(5) * t906 - t1356, -qJ(5) * t855 - t1316; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1257, -t1229, -t1247, -t1227, -t1238, -t1257, -t1050, 0, t1062, 0, -t1162, -t1160, -t1137, t1162, t1134, -t1227, t1501, t1485, -t1056, -t900, -t1082, -t1081, -t1022, t1082, t1018, -t1214, t1320, t1296, t1443, t1444; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1362, -t1236, -t1374, t1362, -t1375, -qJDD(2), -t1057, -t1062, 0, 0, t1110, t1059, t1085, -t1108, t1087, t1123, -t1360, t1361, t1317, -t1422, t951, t905, t957, t950, t958, t982, t1343, t1342, t1344, -t1348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1194, -t1133, t1479, -t1364, t1196, t1364, 0, t1043, -t1369, 0, t1010, t968, t1012, t1009, t1013, t1054, t930, t943, t870, -pkin(8) * t877; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1363, t1136, t1338, -t1319, t1139, t1363, -t1043, 0, t962, 0, t1305, t1325, t1324, t1306, t1323, t1297, t912, t917, t866, t873; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1162, t1160, t1137, -t1162, -t1134, t1227, t1369, -t962, 0, 0, t1082, t1081, t1022, -t1082, -t1018, t1214, t1335, t1308, t966, t876; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1064, -t1017, t1481, t1129, t1119, -t1129, 0, t984, t889, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1130, t1480, t1120, -t1304, t1066, -t1130, -t984, 0, t890, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1082, t1081, t1022, -t1082, -t1018, t1214, -t889, -t890, 0, 0;];
m_new_reg  = t1;
