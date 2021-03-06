% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RRPRP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RRPRP5_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP5_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_invdynm_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:55:27
% EndTime: 2019-12-31 19:55:50
% DurationCPUTime: 23.90s
% Computational Cost: add. (79514->652), mult. (188602->816), div. (0->0), fcn. (133086->8), ass. (0->449)
t1359 = sin(pkin(8));
t1360 = cos(pkin(8));
t1365 = cos(qJ(2));
t1462 = qJD(1) * t1365;
t1362 = sin(qJ(2));
t1463 = qJD(1) * t1362;
t1311 = t1359 * t1463 - t1360 * t1462;
t1313 = (t1359 * t1365 + t1360 * t1362) * qJD(1);
t1361 = sin(qJ(4));
t1364 = cos(qJ(4));
t1270 = -t1311 * t1361 + t1313 * t1364;
t1267 = t1270 ^ 2;
t1356 = qJD(2) + qJD(4);
t1471 = t1356 ^ 2;
t1201 = t1471 + t1267;
t1355 = qJDD(2) + qJDD(4);
t1268 = t1311 * t1364 + t1313 * t1361;
t1457 = t1270 * t1268;
t1481 = t1355 + t1457;
t1439 = t1361 * t1481;
t1125 = t1201 * t1364 + t1439;
t1426 = t1364 * t1481;
t1155 = t1201 * t1361 - t1426;
t1070 = t1125 * t1360 - t1155 * t1359;
t1096 = t1125 * t1359 + t1155 * t1360;
t1039 = t1070 * t1362 + t1096 * t1365;
t1363 = sin(qJ(1));
t1366 = cos(qJ(1));
t1344 = qJD(2) * t1462;
t1417 = t1362 * qJDD(1);
t1321 = t1344 + t1417;
t1350 = t1365 * qJDD(1);
t1412 = qJD(2) * t1463;
t1322 = t1350 - t1412;
t1279 = t1321 * t1360 + t1322 * t1359;
t1385 = t1321 * t1359 - t1322 * t1360;
t1185 = -qJD(4) * t1268 + t1279 * t1364 - t1361 * t1385;
t1452 = t1356 * t1268;
t1483 = t1185 - t1452;
t1556 = pkin(5) * (t1039 * t1366 + t1363 * t1483);
t1555 = pkin(5) * (t1039 * t1363 - t1366 * t1483);
t1011 = t1070 * t1365 - t1096 * t1362;
t1554 = pkin(1) * t1011;
t1553 = pkin(6) * t1011;
t1552 = -pkin(1) * t1483 + pkin(6) * t1039;
t1472 = t1268 ^ 2;
t1482 = t1267 - t1472;
t1400 = -t1361 * t1279 - t1364 * t1385;
t1184 = qJD(4) * t1270 - t1400;
t1261 = t1356 * t1270;
t1484 = t1184 + t1261;
t1088 = -t1361 * t1484 + t1364 * t1483;
t1442 = t1361 * t1483;
t1090 = t1364 * t1484 + t1442;
t1027 = -t1088 * t1360 + t1090 * t1359;
t1031 = t1088 * t1359 + t1090 * t1360;
t989 = t1027 * t1362 - t1031 * t1365;
t1551 = t1363 * t989 - t1366 * t1482;
t1550 = t1363 * t1482 + t1366 * t989;
t1253 = t1472 - t1471;
t1162 = t1253 * t1361 + t1426;
t1166 = t1253 * t1364 - t1439;
t1103 = t1162 * t1360 + t1166 * t1359;
t1108 = t1162 * t1359 - t1166 * t1360;
t1052 = t1103 * t1362 + t1108 * t1365;
t1143 = t1184 - t1261;
t1547 = t1052 * t1363 - t1366 * t1143;
t1546 = t1052 * t1366 + t1363 * t1143;
t1545 = pkin(2) * t1070;
t1544 = qJ(3) * t1070;
t1543 = pkin(2) * t1483 - qJ(3) * t1096;
t983 = t1027 * t1365 + t1031 * t1362;
t1047 = t1103 * t1365 - t1108 * t1362;
t1480 = t1452 + t1185;
t1254 = t1267 - t1471;
t1479 = -t1457 + t1355;
t1438 = t1361 * t1479;
t1513 = -t1254 * t1364 + t1438;
t1196 = t1364 * t1479;
t1514 = t1254 * t1361 + t1196;
t1525 = -t1359 * t1513 + t1360 * t1514;
t1526 = t1359 * t1514 + t1360 * t1513;
t1533 = -t1362 * t1526 + t1365 * t1525;
t1542 = t1363 * t1533 - t1366 * t1480;
t1541 = t1363 * t1480 + t1366 * t1533;
t1478 = -t1471 - t1472;
t1489 = t1364 * t1478 - t1438;
t1490 = t1361 * t1478 + t1196;
t1506 = t1359 * t1489 + t1360 * t1490;
t1507 = -t1359 * t1490 + t1360 * t1489;
t1524 = t1362 * t1507 + t1365 * t1506;
t1540 = pkin(1) * t1524;
t1123 = pkin(3) * t1125;
t1539 = pkin(6) * t1524;
t1538 = pkin(7) * t1125;
t1537 = pkin(7) * t1155;
t1523 = -t1362 * t1506 + t1365 * t1507;
t1536 = -pkin(1) * t1484 + pkin(6) * t1523;
t1535 = pkin(5) * (t1363 * t1523 - t1366 * t1484);
t1534 = pkin(5) * (t1363 * t1484 + t1366 * t1523);
t1532 = t1362 * t1525 + t1365 * t1526;
t1531 = pkin(2) * t1506;
t1530 = qJ(3) * t1506;
t1527 = -pkin(2) * t1484 + qJ(3) * t1507;
t1180 = -t1472 - t1267;
t1522 = pkin(1) * t1180;
t1521 = pkin(2) * t1180;
t1520 = pkin(3) * t1180;
t1519 = pkin(3) * t1490;
t1518 = pkin(7) * t1489;
t1517 = pkin(7) * t1490;
t1516 = t1363 * t1180;
t1515 = t1366 * t1180;
t1309 = t1311 ^ 2;
t1367 = qJD(2) ^ 2;
t1272 = -t1367 - t1309;
t1278 = t1313 * t1311;
t1477 = qJDD(2) - t1278;
t1499 = t1360 * t1477;
t1213 = t1272 * t1359 + t1499;
t1368 = qJD(1) ^ 2;
t1431 = t1362 * t1368;
t1332 = g(1) * t1366 + g(2) * t1363;
t1315 = -pkin(1) * t1368 + qJDD(1) * pkin(6) - t1332;
t1434 = t1362 * t1315;
t1464 = qJD(1) * qJD(2);
t1232 = qJDD(2) * pkin(2) - t1321 * qJ(3) - t1434 + (pkin(2) * t1431 + qJ(3) * t1464 - g(3)) * t1365;
t1294 = -t1362 * g(3) + t1315 * t1365;
t1358 = t1365 ^ 2;
t1352 = t1358 * t1368;
t1380 = qJD(2) * pkin(2) - qJ(3) * t1463;
t1233 = -pkin(2) * t1352 + t1322 * qJ(3) - qJD(2) * t1380 + t1294;
t1388 = -0.2e1 * qJD(3) * t1313 + t1232 * t1360 - t1359 * t1233;
t1512 = pkin(2) * t1213 + t1388;
t1378 = (-t1268 * t1361 - t1270 * t1364) * t1356;
t1451 = t1356 * t1361;
t1250 = t1270 * t1451;
t1450 = t1356 * t1364;
t1409 = t1268 * t1450;
t1390 = t1250 - t1409;
t1475 = -t1359 * t1378 + t1360 * t1390;
t1476 = t1359 * t1390 + t1360 * t1378;
t1486 = -t1362 * t1476 + t1365 * t1475;
t1511 = -t1355 * t1366 + t1363 * t1486;
t1408 = t1366 * t1457;
t1379 = t1184 * t1361 + t1409;
t1391 = -t1184 * t1364 + t1268 * t1451;
t1473 = t1359 * t1379 + t1360 * t1391;
t1474 = -t1359 * t1391 + t1360 * t1379;
t1487 = -t1362 * t1473 + t1365 * t1474;
t1510 = t1363 * t1487 + t1408;
t1411 = t1363 * t1457;
t1509 = t1366 * t1487 - t1411;
t1508 = t1355 * t1363 + t1366 * t1486;
t1501 = t1483 * qJ(5);
t1500 = t1359 * t1477;
t1488 = t1362 * t1474 + t1365 * t1473;
t1485 = t1362 * t1475 + t1365 * t1476;
t1307 = qJD(2) * t1311;
t1244 = t1279 + t1307;
t1310 = t1313 ^ 2;
t1460 = qJD(3) * t1311;
t1303 = -0.2e1 * t1460;
t1419 = t1232 * t1359 + t1233 * t1360;
t1159 = t1303 + t1419;
t1099 = t1159 * t1359 + t1360 * t1388;
t1470 = pkin(2) * t1099;
t1461 = qJD(2) * t1313;
t1242 = -t1385 + t1461;
t1187 = t1242 * t1359 - t1244 * t1360;
t1469 = pkin(2) * t1187;
t1467 = pkin(4) * t1364;
t1466 = t1184 * pkin(4);
t1465 = qJ(5) * t1364;
t1458 = qJD(5) * t1356;
t1456 = t1311 * t1359;
t1455 = t1311 * t1360;
t1454 = t1313 * t1359;
t1453 = t1313 * t1360;
t1111 = pkin(3) * t1477 - pkin(7) * t1244 + t1388;
t1395 = qJD(2) * pkin(3) - pkin(7) * t1313;
t1117 = -t1309 * pkin(3) - pkin(7) * t1385 - qJD(2) * t1395 + t1159;
t1059 = -t1111 * t1364 + t1117 * t1361;
t1060 = t1111 * t1361 + t1117 * t1364;
t1008 = -t1059 * t1364 + t1060 * t1361;
t1449 = t1359 * t1008;
t1331 = t1363 * g(1) - g(2) * t1366;
t1381 = qJDD(1) * pkin(1) + t1331;
t1238 = t1322 * pkin(2) - qJDD(3) - t1380 * t1463 + (qJ(3) * t1358 + pkin(6)) * t1368 + t1381;
t1448 = t1359 * t1238;
t1274 = qJDD(2) + t1278;
t1447 = t1359 * t1274;
t1446 = t1360 * t1008;
t1445 = t1360 * t1238;
t1444 = t1360 * t1274;
t1441 = t1361 * t1480;
t1157 = -pkin(3) * t1385 + t1309 * pkin(7) - t1313 * t1395 + t1238;
t1440 = t1361 * t1157;
t1436 = t1362 * t1099;
t1314 = t1368 * pkin(6) + t1381;
t1435 = t1362 * t1314;
t1339 = t1365 * t1431;
t1329 = qJDD(2) + t1339;
t1433 = t1362 * t1329;
t1330 = qJDD(2) - t1339;
t1432 = t1362 * t1330;
t1427 = t1364 * t1157;
t1425 = t1365 * t1099;
t1424 = t1365 * t1314;
t1423 = t1365 * t1330;
t1342 = 0.2e1 * t1458;
t1207 = pkin(4) * t1268 - qJ(5) * t1270;
t1386 = -pkin(4) * t1471 + qJ(5) * t1355 - t1207 * t1268 + t1060;
t1036 = t1342 + t1386;
t1042 = -pkin(4) * t1355 - qJ(5) * t1471 + t1207 * t1270 + qJDD(5) + t1059;
t1421 = -pkin(4) * t1042 + qJ(5) * t1036;
t1420 = -pkin(4) * t1480 - qJ(5) * t1143;
t1357 = t1362 ^ 2;
t1418 = t1357 + t1358;
t1416 = t1363 * qJDD(1);
t1415 = t1366 * qJDD(1);
t1414 = t1366 * qJDD(2);
t996 = t1036 * t1361 - t1042 * t1364;
t1413 = pkin(3) * t996 + t1421;
t1410 = t1363 * t1278;
t1407 = t1366 * t1278;
t1406 = -t1060 - t1123;
t1138 = t1364 * t1480;
t1087 = -t1143 * t1361 - t1138;
t1405 = pkin(3) * t1087 + t1420;
t1007 = pkin(3) * t1008;
t1009 = t1059 * t1361 + t1060 * t1364;
t975 = t1009 * t1359 + t1446;
t1404 = pkin(2) * t975 + t1007;
t1403 = -qJ(5) * t1361 - pkin(3);
t1144 = (-qJD(4) + t1356) * t1270 + t1400;
t1089 = t1144 * t1361 - t1138;
t1093 = t1144 * t1364 + t1441;
t1030 = t1089 * t1360 + t1093 * t1359;
t1085 = pkin(3) * t1089;
t1402 = pkin(2) * t1030 + t1085;
t1100 = t1159 * t1360 - t1359 * t1388;
t1293 = t1365 * g(3) + t1434;
t1236 = t1293 * t1362 + t1294 * t1365;
t1399 = -t1331 * t1363 - t1332 * t1366;
t1398 = t1363 * t1339;
t1397 = t1366 * t1339;
t1299 = -t1310 - t1367;
t1219 = t1299 * t1360 - t1447;
t1396 = pkin(2) * t1219 - t1419;
t1326 = -t1363 * t1368 + t1415;
t1394 = -pkin(5) * t1326 - g(3) * t1363;
t1136 = t1185 * t1361 + t1270 * t1450;
t1137 = t1185 * t1364 - t1250;
t1079 = t1136 * t1360 + t1137 * t1359;
t1082 = -t1136 * t1359 + t1137 * t1360;
t1024 = -t1079 * t1362 + t1082 * t1365;
t1393 = t1024 * t1363 - t1408;
t1392 = t1024 * t1366 + t1411;
t1389 = -t1059 + t1519;
t997 = t1036 * t1364 + t1042 * t1361;
t969 = t1359 * t997 + t1360 * t996;
t1387 = pkin(2) * t969 + t1413;
t1235 = t1293 * t1365 - t1294 * t1362;
t1384 = t1331 * t1366 - t1332 * t1363;
t1091 = -t1143 * t1364 + t1441;
t1028 = t1087 * t1360 + t1091 * t1359;
t1383 = pkin(2) * t1028 + t1405;
t1382 = t1406 - t1545;
t1377 = t1389 + t1531;
t1376 = pkin(4) * t1201 + qJ(5) * t1481 + t1386;
t1375 = t1342 + t1376;
t1374 = t1123 + t1376 + t1545;
t1373 = pkin(4) * t1479 + qJ(5) * t1478 - t1042;
t1372 = t1373 + t1519;
t1371 = t1372 + t1531;
t1370 = -pkin(4) * t1261 + 0.2e1 * qJD(5) * t1270 + t1157;
t1369 = t1370 + t1501;
t1351 = t1357 * t1368;
t1349 = t1363 * qJDD(2);
t1338 = -t1352 - t1367;
t1337 = t1352 - t1367;
t1336 = -t1351 - t1367;
t1335 = -t1351 + t1367;
t1328 = t1352 - t1351;
t1327 = t1351 + t1352;
t1325 = t1366 * t1368 + t1416;
t1324 = t1418 * qJDD(1);
t1323 = t1350 - 0.2e1 * t1412;
t1320 = 0.2e1 * t1344 + t1417;
t1318 = t1365 * t1329;
t1317 = t1418 * t1464;
t1308 = -pkin(5) * t1325 + g(3) * t1366;
t1298 = -t1310 + t1367;
t1297 = t1309 - t1367;
t1296 = t1321 * t1365 - t1357 * t1464;
t1295 = -t1322 * t1362 - t1358 * t1464;
t1292 = -t1336 * t1362 - t1423;
t1291 = -t1335 * t1362 + t1318;
t1290 = t1338 * t1365 - t1433;
t1289 = t1337 * t1365 - t1432;
t1288 = t1336 * t1365 - t1432;
t1287 = t1335 * t1365 + t1433;
t1286 = t1338 * t1362 + t1318;
t1285 = t1337 * t1362 + t1423;
t1284 = (t1321 + t1344) * t1362;
t1283 = (t1322 - t1412) * t1365;
t1281 = -t1320 * t1362 + t1323 * t1365;
t1280 = t1320 * t1365 + t1323 * t1362;
t1277 = t1310 - t1309;
t1263 = (t1454 - t1455) * qJD(2);
t1262 = (-t1453 - t1456) * qJD(2);
t1252 = -pkin(6) * t1288 - t1424;
t1251 = -pkin(6) * t1286 - t1435;
t1246 = -pkin(1) * t1288 + t1294;
t1245 = -pkin(1) * t1286 + t1293;
t1243 = t1279 - t1307;
t1240 = t1385 + t1461;
t1239 = -t1309 - t1310;
t1231 = -qJD(2) * t1454 + t1279 * t1360;
t1230 = qJD(2) * t1453 + t1279 * t1359;
t1229 = qJD(2) * t1455 + t1359 * t1385;
t1228 = qJD(2) * t1456 - t1360 * t1385;
t1226 = pkin(1) * t1323 + pkin(6) * t1290 + t1424;
t1225 = -pkin(1) * t1320 + pkin(6) * t1292 - t1435;
t1222 = -t1299 * t1359 - t1444;
t1221 = -t1298 * t1359 + t1499;
t1220 = t1297 * t1360 - t1447;
t1218 = t1298 * t1360 + t1500;
t1217 = t1297 * t1359 + t1444;
t1216 = pkin(1) * t1314 + pkin(6) * t1236;
t1215 = pkin(1) * t1327 + pkin(6) * t1324 + t1236;
t1214 = t1272 * t1360 - t1500;
t1195 = -t1262 * t1362 + t1263 * t1365;
t1194 = t1262 * t1365 + t1263 * t1362;
t1189 = t1242 * t1360 + t1244 * t1359;
t1188 = -t1240 * t1360 - t1243 * t1359;
t1186 = -t1240 * t1359 + t1243 * t1360;
t1182 = -qJ(3) * t1219 - t1445;
t1179 = -t1230 * t1362 + t1231 * t1365;
t1178 = -t1228 * t1362 + t1229 * t1365;
t1177 = t1230 * t1365 + t1231 * t1362;
t1176 = t1228 * t1365 + t1229 * t1362;
t1174 = -t1219 * t1362 + t1222 * t1365;
t1173 = -t1218 * t1362 + t1221 * t1365;
t1172 = -t1217 * t1362 + t1220 * t1365;
t1171 = t1219 * t1365 + t1222 * t1362;
t1170 = t1218 * t1365 + t1221 * t1362;
t1169 = t1217 * t1365 + t1220 * t1362;
t1168 = -qJ(3) * t1213 - t1448;
t1131 = -t1213 * t1362 + t1214 * t1365;
t1130 = t1213 * t1365 + t1214 * t1362;
t1129 = -pkin(2) * t1243 + qJ(3) * t1222 - t1448;
t1122 = -pkin(2) * t1240 + qJ(3) * t1214 + t1445;
t1115 = -t1187 * t1362 + t1189 * t1365;
t1114 = -t1186 * t1362 + t1188 * t1365;
t1113 = t1187 * t1365 + t1189 * t1362;
t1112 = t1186 * t1365 + t1188 * t1362;
t1098 = -t1427 + t1538;
t1083 = -t1440 - t1517;
t1076 = -pkin(1) * t1113 - t1469;
t1075 = -pkin(1) * t1171 + t1303 - t1396;
t1074 = pkin(2) * t1238 + qJ(3) * t1100;
t1069 = -pkin(1) * t1130 - t1512;
t1068 = -qJ(3) * t1187 - t1099;
t1063 = -pkin(6) * t1171 - t1129 * t1362 + t1182 * t1365;
t1062 = -pkin(3) * t1483 - t1440 + t1537;
t1061 = -pkin(2) * t1239 + qJ(3) * t1189 + t1100;
t1057 = -pkin(3) * t1484 + t1427 + t1518;
t1056 = -pkin(1) * t1243 + pkin(6) * t1174 + t1129 * t1365 + t1182 * t1362;
t1055 = -pkin(6) * t1130 - t1122 * t1362 + t1168 * t1365;
t1054 = t1369 - t1466;
t1053 = -pkin(1) * t1240 + pkin(6) * t1131 + t1122 * t1365 + t1168 * t1362;
t1044 = t1100 * t1365 - t1436;
t1043 = t1100 * t1362 + t1425;
t1034 = -t1089 * t1359 + t1093 * t1360;
t1032 = -t1087 * t1359 + t1091 * t1360;
t1026 = t1369 + (-t1184 - t1484) * pkin(4);
t1025 = t1370 - t1466 + 0.2e1 * t1501;
t1021 = t1079 * t1365 + t1082 * t1362;
t1016 = -qJ(5) * t1180 + t1042;
t1015 = -pkin(4) * t1180 + t1036;
t1010 = -pkin(1) * t1043 - t1470;
t1006 = -t1026 * t1361 - t1465 * t1484 - t1517;
t1005 = -pkin(4) * t1442 + t1025 * t1364 - t1538;
t1004 = -t1062 * t1359 + t1098 * t1360 + t1544;
t1003 = pkin(3) * t1157 + pkin(7) * t1009;
t1002 = -pkin(6) * t1113 - t1061 * t1362 + t1068 * t1365;
t1001 = -pkin(1) * t1239 + pkin(6) * t1115 + t1061 * t1365 + t1068 * t1362;
t1000 = t1364 * t1026 + t1403 * t1484 + t1518;
t999 = -t1537 + t1361 * t1025 + (pkin(3) + t1467) * t1483;
t998 = -t1057 * t1359 + t1083 * t1360 - t1530;
t994 = t1062 * t1360 + t1098 * t1359 - t1543;
t993 = -pkin(6) * t1043 - qJ(3) * t1425 - t1074 * t1362;
t992 = -pkin(7) * t1089 - t1008;
t991 = pkin(1) * t1238 + pkin(6) * t1044 - qJ(3) * t1436 + t1074 * t1365;
t990 = -t1030 * t1362 + t1034 * t1365;
t988 = -t1028 * t1362 + t1032 * t1365;
t986 = t1030 * t1365 + t1034 * t1362;
t984 = t1028 * t1365 + t1032 * t1362;
t982 = t1057 * t1360 + t1083 * t1359 + t1527;
t981 = pkin(7) * t1093 + t1009 - t1520;
t980 = -t1382 + t1554;
t979 = -t1377 - t1540;
t978 = -t1371 - t1540;
t977 = -pkin(7) * t1087 - t1015 * t1361 + t1016 * t1364;
t976 = t1009 * t1360 - t1449;
t974 = pkin(7) * t1091 + t1015 * t1364 + t1016 * t1361 - t1520;
t973 = -t1374 - 0.2e1 * t1458 - t1554;
t972 = -pkin(7) * t996 + (-pkin(4) * t1361 + t1465) * t1054;
t971 = -pkin(1) * t986 - t1402;
t970 = -t1359 * t996 + t1360 * t997;
t968 = -t1000 * t1359 + t1006 * t1360 - t1530;
t967 = t1005 * t1360 - t1359 * t999 - t1544;
t966 = t1000 * t1360 + t1006 * t1359 + t1527;
t965 = -pkin(1) * t984 - t1383;
t964 = pkin(7) * t997 + (-t1403 + t1467) * t1054;
t963 = t1005 * t1359 + t1360 * t999 + t1543;
t962 = t1004 * t1365 - t1362 * t994 + t1553;
t961 = t1004 * t1362 + t1365 * t994 + t1552;
t960 = -qJ(3) * t1030 - t1359 * t981 + t1360 * t992;
t959 = -t1362 * t982 + t1365 * t998 - t1539;
t958 = qJ(3) * t1034 + t1359 * t992 + t1360 * t981 - t1521;
t957 = t1362 * t998 + t1365 * t982 + t1536;
t956 = -t1362 * t975 + t1365 * t976;
t955 = t1362 * t976 + t1365 * t975;
t954 = -pkin(7) * t1446 - qJ(3) * t975 - t1003 * t1359;
t953 = pkin(2) * t1157 - pkin(7) * t1449 + qJ(3) * t976 + t1003 * t1360;
t952 = -qJ(3) * t1028 - t1359 * t974 + t1360 * t977;
t951 = qJ(3) * t1032 + t1359 * t977 + t1360 * t974 - t1521;
t950 = -t1362 * t969 + t1365 * t970;
t949 = t1362 * t970 + t1365 * t969;
t948 = -t1362 * t966 + t1365 * t968 - t1539;
t947 = t1362 * t968 + t1365 * t966 + t1536;
t946 = -t1362 * t963 + t1365 * t967 - t1553;
t945 = t1362 * t967 + t1365 * t963 - t1552;
t944 = -pkin(1) * t955 - t1404;
t943 = -qJ(3) * t969 - t1359 * t964 + t1360 * t972;
t942 = pkin(2) * t1054 + qJ(3) * t970 + t1359 * t972 + t1360 * t964;
t941 = -pkin(6) * t986 - t1362 * t958 + t1365 * t960;
t940 = pkin(6) * t990 + t1362 * t960 + t1365 * t958 - t1522;
t939 = -pkin(1) * t949 - t1387;
t938 = -pkin(6) * t984 - t1362 * t951 + t1365 * t952;
t937 = pkin(6) * t988 + t1362 * t952 + t1365 * t951 - t1522;
t936 = -pkin(6) * t955 - t1362 * t953 + t1365 * t954;
t935 = pkin(1) * t1157 + pkin(6) * t956 + t1362 * t954 + t1365 * t953;
t934 = -pkin(6) * t949 - t1362 * t942 + t1365 * t943;
t933 = pkin(1) * t1054 + pkin(6) * t950 + t1362 * t943 + t1365 * t942;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1326, 0, -t1325, 0, t1394, -t1308, -t1384, -pkin(5) * t1384, t1296 * t1366 - t1398, t1281 * t1366 - t1328 * t1363, t1291 * t1366 + t1362 * t1416, t1295 * t1366 + t1398, t1289 * t1366 + t1350 * t1363, t1317 * t1366 + t1349, t1366 * t1251 - t1363 * t1245 - pkin(5) * (t1290 * t1363 + t1323 * t1366), t1366 * t1252 - t1363 * t1246 - pkin(5) * (t1292 * t1363 - t1320 * t1366), t1366 * t1235 - pkin(5) * (t1324 * t1363 + t1327 * t1366), -pkin(5) * (t1236 * t1363 + t1314 * t1366) - (pkin(1) * t1363 - pkin(6) * t1366) * t1235, t1179 * t1366 + t1410, t1114 * t1366 + t1277 * t1363, t1173 * t1366 + t1244 * t1363, t1178 * t1366 - t1410, t1172 * t1366 + t1242 * t1363, t1195 * t1366 + t1349, t1366 * t1055 - t1363 * t1069 - pkin(5) * (t1131 * t1363 - t1240 * t1366), t1366 * t1063 - t1363 * t1075 - pkin(5) * (t1174 * t1363 - t1243 * t1366), t1366 * t1002 - t1363 * t1076 - pkin(5) * (t1115 * t1363 - t1239 * t1366), t1366 * t993 - t1363 * t1010 - pkin(5) * (t1044 * t1363 + t1238 * t1366), t1392, t1550, t1541, t1509, -t1546, t1508, -t1363 * t979 + t1366 * t959 - t1535, -t1363 * t980 + t1366 * t962 - t1555, t1366 * t941 - t1363 * t971 - pkin(5) * (t1363 * t990 - t1515), t1366 * t936 - t1363 * t944 - pkin(5) * (t1157 * t1366 + t1363 * t956), t1392, t1541, -t1550, t1508, t1546, t1509, -t1363 * t978 + t1366 * t948 - t1535, t1366 * t938 - t1363 * t965 - pkin(5) * (t1363 * t988 - t1515), -t1363 * t973 + t1366 * t946 + t1555, t1366 * t934 - t1363 * t939 - pkin(5) * (t1054 * t1366 + t1363 * t950); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1325, 0, t1326, 0, t1308, t1394, t1399, pkin(5) * t1399, t1296 * t1363 + t1397, t1281 * t1363 + t1328 * t1366, t1291 * t1363 - t1362 * t1415, t1295 * t1363 - t1397, t1289 * t1363 - t1365 * t1415, t1317 * t1363 - t1414, t1363 * t1251 + t1366 * t1245 + pkin(5) * (t1290 * t1366 - t1323 * t1363), t1363 * t1252 + t1366 * t1246 + pkin(5) * (t1292 * t1366 + t1320 * t1363), t1363 * t1235 + pkin(5) * (t1324 * t1366 - t1327 * t1363), pkin(5) * (t1236 * t1366 - t1314 * t1363) - (-pkin(1) * t1366 - pkin(6) * t1363) * t1235, t1179 * t1363 - t1407, t1114 * t1363 - t1277 * t1366, t1173 * t1363 - t1244 * t1366, t1178 * t1363 + t1407, t1172 * t1363 - t1242 * t1366, t1195 * t1363 - t1414, t1363 * t1055 + t1366 * t1069 + pkin(5) * (t1131 * t1366 + t1240 * t1363), t1363 * t1063 + t1366 * t1075 + pkin(5) * (t1174 * t1366 + t1243 * t1363), t1363 * t1002 + t1366 * t1076 + pkin(5) * (t1115 * t1366 + t1239 * t1363), t1363 * t993 + t1366 * t1010 + pkin(5) * (t1044 * t1366 - t1238 * t1363), t1393, t1551, t1542, t1510, -t1547, t1511, t1363 * t959 + t1366 * t979 + t1534, t1363 * t962 + t1366 * t980 + t1556, t1363 * t941 + t1366 * t971 + pkin(5) * (t1366 * t990 + t1516), t1363 * t936 + t1366 * t944 + pkin(5) * (-t1157 * t1363 + t1366 * t956), t1393, t1542, -t1551, t1511, t1547, t1510, t1363 * t948 + t1366 * t978 + t1534, t1363 * t938 + t1366 * t965 + pkin(5) * (t1366 * t988 + t1516), t1363 * t946 + t1366 * t973 - t1556, t1363 * t934 + t1366 * t939 + pkin(5) * (-t1054 * t1363 + t1366 * t950); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1331, t1332, 0, 0, t1284, t1280, t1287, t1283, t1285, 0, t1226, t1225, t1215, t1216, t1177, t1112, t1170, t1176, t1169, t1194, t1053, t1056, t1001, t991, t1021, -t983, t1532, t1488, t1047, t1485, t957, t961, t940, t935, t1021, t1532, t983, t1485, -t1047, t1488, t947, t937, t945, t933; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1368, 0, 0, -g(3), -t1331, 0, t1296, t1281, t1291, t1295, t1289, t1317, t1251, t1252, t1235, pkin(6) * t1235, t1179, t1114, t1173, t1178, t1172, t1195, t1055, t1063, t1002, t993, t1024, t989, t1533, t1487, -t1052, t1486, t959, t962, t941, t936, t1024, t1533, -t989, t1486, t1052, t1487, t948, t938, t946, t934; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1368, 0, qJDD(1), 0, g(3), 0, -t1332, 0, t1339, t1328, -t1417, -t1339, -t1350, -qJDD(2), t1245, t1246, 0, pkin(1) * t1235, -t1278, -t1277, -t1244, t1278, -t1242, -qJDD(2), t1069, t1075, t1076, t1010, -t1457, -t1482, -t1480, t1457, t1143, -t1355, t979, t980, t971, t944, -t1457, -t1480, t1482, -t1355, -t1143, t1457, t978, t965, t973, t939; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1331, t1332, 0, 0, t1284, t1280, t1287, t1283, t1285, 0, t1226, t1225, t1215, t1216, t1177, t1112, t1170, t1176, t1169, t1194, t1053, t1056, t1001, t991, t1021, -t983, t1532, t1488, t1047, t1485, t957, t961, t940, t935, t1021, t1532, t983, t1485, -t1047, t1488, t947, t937, t945, t933; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1321, t1323, t1329, -t1344, t1337, t1344, 0, -t1314, t1293, 0, t1231, t1188, t1221, t1229, t1220, t1263, t1168, t1182, t1068, -qJ(3) * t1099, t1082, -t1031, t1525, t1474, -t1108, t1475, t998, t1004, t960, t954, t1082, t1525, t1031, t1475, t1108, t1474, t968, t952, t967, t943; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1412, t1320, t1335, t1322, t1330, -t1412, t1314, 0, t1294, 0, t1230, t1186, t1218, t1228, t1217, t1262, t1122, t1129, t1061, t1074, t1079, -t1027, t1526, t1473, t1103, t1476, t982, t994, t958, t953, t1079, t1526, t1027, t1476, -t1103, t1473, t966, t951, t963, t942; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1339, -t1328, t1417, t1339, t1350, qJDD(2), -t1293, -t1294, 0, 0, t1278, t1277, t1244, -t1278, t1242, qJDD(2), t1512, t1396 + 0.2e1 * t1460, t1469, t1470, t1457, t1482, t1480, -t1457, -t1143, t1355, t1377, t1382, t1402, t1404, t1457, t1480, -t1482, t1355, t1143, -t1457, t1371, t1383, t1342 + t1374, t1387; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1279, -t1240, t1477, t1307, t1297, -t1307, 0, -t1238, -t1388, 0, t1137, -t1090, t1514, t1379, t1166, t1390, t1083, t1098, t992, -pkin(7) * t1008, t1137, t1514, t1090, t1390, -t1166, t1379, t1006, t977, t1005, t972; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1461, t1243, t1298, -t1385, t1274, -t1461, t1238, 0, t1159, 0, t1136, t1088, t1513, t1391, t1162, t1378, t1057, t1062, t981, t1003, t1136, t1513, -t1088, t1378, -t1162, t1391, t1000, t974, t999, t964; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1278, t1277, t1244, -t1278, t1242, qJDD(2), t1388, -t1159, 0, 0, t1457, t1482, t1480, -t1457, -t1143, t1355, t1389, t1406, t1085, t1007, t1457, t1480, -t1482, t1355, t1143, -t1457, t1372, t1405, t1123 + t1375, t1413; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1185, -t1484, t1479, t1452, t1253, -t1452, 0, -t1157, t1059, 0, t1185, t1479, t1484, -t1452, -t1253, t1452, -qJ(5) * t1484, t1016, t1025, qJ(5) * t1054; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1261, t1483, -t1254, -t1184, t1481, -t1261, t1157, 0, t1060, 0, t1261, -t1254, -t1483, -t1261, -t1481, -t1184, t1026, t1015, pkin(4) * t1483, pkin(4) * t1054; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1457, t1482, t1480, -t1457, -t1143, t1355, -t1059, -t1060, 0, 0, t1457, t1480, -t1482, t1355, t1143, -t1457, t1373, t1420, t1375, t1421; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1185, t1479, t1484, -t1452, -t1253, t1452, 0, t1042, t1054, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1457, t1480, -t1482, t1355, t1143, -t1457, -t1042, 0, t1036, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1261, t1254, t1483, t1261, t1481, t1184, -t1054, -t1036, 0, 0;];
m_new_reg = t1;
