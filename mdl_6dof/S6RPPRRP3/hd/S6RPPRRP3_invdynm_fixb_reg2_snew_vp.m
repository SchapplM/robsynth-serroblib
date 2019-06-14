% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6RPPRRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6RPPRRP3_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP3_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:52:22
% EndTime: 2019-05-05 14:52:45
% DurationCPUTime: 24.93s
% Computational Cost: add. (50210->711), mult. (94020->806), div. (0->0), fcn. (56422->8), ass. (0->459)
t1332 = sin(qJ(1));
t1335 = cos(qJ(1));
t1333 = cos(qJ(5));
t1330 = sin(qJ(5));
t1334 = cos(qJ(4));
t1476 = qJD(1) * t1334;
t1272 = -t1333 * qJD(4) + t1330 * t1476;
t1313 = t1334 * qJDD(1);
t1331 = sin(qJ(4));
t1477 = qJD(1) * t1331;
t1425 = qJD(4) * t1477;
t1282 = t1313 - t1425;
t1385 = -t1330 * qJDD(4) - t1333 * t1282;
t1200 = -qJD(5) * t1272 - t1385;
t1305 = qJD(5) + t1477;
t1450 = t1272 * t1305;
t1518 = t1200 - t1450;
t1467 = t1518 * t1330;
t1274 = qJD(4) * t1330 + t1333 * t1476;
t1254 = t1305 * t1274;
t1409 = t1333 * qJDD(4) - t1330 * t1282;
t1367 = qJD(5) * t1274 - t1409;
t1514 = t1254 + t1367;
t1063 = -t1514 * t1333 - t1467;
t1269 = t1274 ^ 2;
t1498 = t1272 ^ 2;
t1513 = t1269 - t1498;
t1024 = t1063 * t1331 - t1334 * t1513;
t1466 = t1518 * t1333;
t1059 = -t1514 * t1330 + t1466;
t1325 = sin(pkin(9));
t1326 = cos(pkin(9));
t958 = t1024 * t1326 - t1059 * t1325;
t961 = t1024 * t1325 + t1059 * t1326;
t1614 = t1332 * t961 - t1335 * t958;
t1613 = t1332 * t958 + t1335 * t961;
t1497 = t1305 ^ 2;
t1241 = t1498 - t1497;
t1223 = t1274 * t1272;
t1423 = qJD(4) * t1476;
t1430 = t1331 * qJDD(1);
t1281 = -t1423 - t1430;
t1270 = qJDD(5) - t1281;
t1515 = t1223 + t1270;
t1542 = t1515 * t1330;
t1123 = t1241 * t1333 - t1542;
t1145 = -t1254 + t1367;
t1041 = t1123 * t1331 + t1145 * t1334;
t1541 = t1515 * t1333;
t1119 = t1241 * t1330 + t1541;
t988 = t1041 * t1326 - t1119 * t1325;
t991 = t1041 * t1325 + t1119 * t1326;
t1612 = t1332 * t991 - t1335 * t988;
t1611 = t1332 * t988 + t1335 * t991;
t1242 = t1269 - t1497;
t1516 = -t1223 + t1270;
t1460 = t1516 * t1330;
t1564 = -t1242 * t1333 + t1460;
t1517 = t1200 + t1450;
t1459 = t1516 * t1333;
t1563 = t1242 * t1330 + t1459;
t1575 = t1331 * t1563 - t1334 * t1517;
t1591 = t1325 * t1564 - t1326 * t1575;
t1592 = t1325 * t1575 + t1326 * t1564;
t1610 = t1332 * t1592 + t1335 * t1591;
t1609 = -t1332 * t1591 + t1335 * t1592;
t1511 = -t1497 - t1498;
t1536 = t1330 * t1511 + t1459;
t1535 = t1333 * t1511 - t1460;
t1560 = t1331 * t1535 - t1334 * t1514;
t1576 = t1325 * t1536 - t1326 * t1560;
t1608 = pkin(1) * t1576;
t1202 = t1497 + t1269;
t1093 = t1202 * t1333 + t1542;
t1607 = pkin(3) * t1093;
t1606 = pkin(4) * t1093;
t1605 = pkin(8) * t1093;
t1110 = t1202 * t1330 - t1541;
t1604 = pkin(8) * t1110;
t1603 = qJ(2) * t1576;
t1602 = qJ(3) * t1093;
t1599 = t1093 * t1325;
t1598 = t1093 * t1326;
t1597 = t1110 * t1331;
t1596 = t1110 * t1334;
t1336 = qJD(1) ^ 2;
t1285 = -t1326 * qJDD(1) + t1325 * t1336;
t1322 = g(3) - qJDD(2);
t1246 = qJ(2) * t1285 - t1322 * t1325;
t1284 = qJDD(1) * t1325 + t1326 * t1336;
t1414 = t1335 * t1284 - t1285 * t1332;
t1419 = -qJ(2) * t1284 + t1326 * t1322;
t1580 = -pkin(6) * t1414 + t1246 * t1332 + t1335 * t1419;
t1559 = t1331 * t1514 + t1334 * t1535;
t1577 = t1325 * t1560 + t1326 * t1536;
t1593 = -pkin(1) * t1559 + qJ(2) * t1577;
t1026 = t1063 * t1334 + t1331 * t1513;
t1045 = t1123 * t1334 - t1145 * t1331;
t1590 = pkin(6) * (t1332 * t1577 + t1335 * t1576);
t1589 = pkin(6) * (-t1332 * t1576 + t1335 * t1577);
t1587 = pkin(3) * t1560;
t1586 = pkin(7) * t1559;
t1585 = pkin(7) * t1560;
t1584 = qJ(3) * t1559;
t1495 = -pkin(2) - pkin(7);
t1581 = t1495 * t1559;
t1578 = -pkin(2) * t1560 + qJ(3) * t1536;
t1574 = t1331 * t1517 + t1334 * t1563;
t1573 = pkin(3) * t1536;
t1572 = pkin(4) * t1536;
t1570 = pkin(8) * t1535;
t1569 = pkin(8) * t1536;
t1568 = -qJ(6) * t1330 - pkin(4);
t1296 = g(1) * t1335 + g(2) * t1332;
t1276 = -pkin(1) * t1336 - t1296;
t1295 = g(1) * t1332 - t1335 * g(2);
t1380 = qJDD(1) * pkin(1) + t1295;
t1203 = t1276 * t1325 - t1326 * t1380;
t1204 = t1326 * t1276 + t1325 * t1380;
t1415 = t1203 * t1325 + t1326 * t1204;
t1115 = t1203 * t1326 - t1204 * t1325;
t1473 = t1115 * t1332;
t1562 = t1335 * t1415 + t1473;
t1472 = t1115 * t1335;
t1561 = -t1332 * t1415 + t1472;
t1441 = t1305 * t1330;
t1236 = t1274 * t1441;
t1440 = t1305 * t1333;
t1424 = t1272 * t1440;
t1394 = t1236 - t1424;
t1506 = -t1270 * t1334 + t1331 * t1394;
t1519 = (t1272 * t1330 + t1274 * t1333) * t1305;
t1531 = t1325 * t1506 - t1326 * t1519;
t1532 = -t1325 * t1519 - t1326 * t1506;
t1558 = -t1332 * t1532 + t1335 * t1531;
t1557 = t1332 * t1531 + t1335 * t1532;
t1132 = t1272 * t1441 - t1333 * t1367;
t1342 = t1330 * t1367 + t1424;
t1426 = t1334 * t1223;
t1507 = t1331 * t1342 + t1426;
t1533 = t1132 * t1325 - t1326 * t1507;
t1534 = t1132 * t1326 + t1325 * t1507;
t1556 = t1332 * t1534 + t1335 * t1533;
t1555 = -t1332 * t1533 + t1335 * t1534;
t1412 = t1284 * t1332 + t1335 * t1285;
t1500 = pkin(6) * t1412 + t1246 * t1335 - t1332 * t1419;
t1512 = t1269 + t1498;
t1554 = pkin(4) * t1512;
t1553 = qJ(6) * t1518;
t1550 = t1331 * t1512;
t1546 = t1334 * t1512;
t1475 = qJD(6) * t1305;
t1293 = 0.2e1 * t1475;
t1328 = t1336 * pkin(7);
t1474 = qJDD(1) * qJ(3);
t1339 = -t1336 * pkin(2) + t1204 + t1474;
t1487 = pkin(4) * t1334;
t1400 = pkin(8) * t1331 + t1487;
t1496 = 2 * qJD(3);
t1088 = -t1328 - t1282 * pkin(8) - t1281 * pkin(4) + (qJD(4) * t1400 + t1496) * qJD(1) + t1339;
t1488 = pkin(4) * t1331;
t1539 = -pkin(8) * t1334 + t1488;
t1279 = t1539 * qJD(1);
t1324 = qJDD(1) * pkin(2);
t1390 = qJDD(3) + t1203;
t1180 = -t1336 * qJ(3) - t1324 + t1390;
t1345 = -qJDD(1) * pkin(7) + t1180;
t1432 = t1334 * t1322 - t1331 * t1345;
t1499 = qJD(4) ^ 2;
t1107 = -pkin(4) * t1499 + qJDD(4) * pkin(8) - t1279 * t1477 - t1432;
t1000 = t1330 * t1088 + t1333 * t1107;
t1207 = pkin(5) * t1272 - qJ(6) * t1274;
t1388 = -pkin(5) * t1497 + t1270 * qJ(6) - t1272 * t1207 + t1000;
t971 = t1293 + t1388;
t999 = -t1333 * t1088 + t1107 * t1330;
t972 = -t1270 * pkin(5) - qJ(6) * t1497 + t1207 * t1274 + qJDD(6) + t999;
t930 = t1330 * t972 + t1333 * t971;
t1417 = t1331 * t1322 + t1334 * t1345;
t1106 = -qJDD(4) * pkin(4) - pkin(8) * t1499 + t1279 * t1476 - t1417;
t1341 = t1367 * pkin(5) + t1106 - t1553;
t984 = (pkin(5) * t1305 - 0.2e1 * qJD(6)) * t1274 + t1341;
t1348 = pkin(8) * t930 + (-pkin(5) * t1333 + t1568) * t984;
t920 = t1331 * t930 - t1334 * t984;
t1540 = -pkin(3) * t920 - t1348;
t1429 = qJD(1) * t1496;
t1176 = t1339 + t1429;
t1083 = t1176 * t1325 - t1180 * t1326;
t1416 = t1326 * t1176 + t1180 * t1325;
t1538 = -t1083 * t1332 + t1335 * t1416;
t1537 = t1083 * t1335 + t1332 * t1416;
t1136 = t1200 * t1330 + t1274 * t1440;
t1137 = t1200 * t1333 - t1236;
t1368 = t1137 * t1331 - t1426;
t1508 = t1326 * t1136 + t1325 * t1368;
t1509 = t1325 * t1136 - t1326 * t1368;
t1530 = -t1332 * t1509 + t1335 * t1508;
t1529 = t1332 * t1508 + t1335 * t1509;
t947 = t1333 * t1000 + t1330 * t999;
t946 = t1000 * t1330 - t1333 * t999;
t1526 = (pkin(3) + t1400) * t946;
t1066 = -t1331 * t1432 + t1334 * t1417;
t1146 = (-qJD(5) + t1305) * t1274 + t1409;
t1464 = t1517 * t1333;
t1060 = t1146 * t1330 - t1464;
t932 = -pkin(8) * t1060 - t946;
t1510 = -t1060 * (pkin(3) + t1487) + t1331 * t932;
t1427 = t1331 * t1223;
t1505 = t1334 * t1342 - t1427;
t1504 = t1270 * t1331 + t1334 * t1394;
t1493 = pkin(1) * t1284;
t1492 = pkin(1) * t1285;
t1491 = pkin(3) * t1066;
t1167 = -t1328 + t1176;
t1490 = pkin(3) * t1167;
t1320 = t1331 ^ 2;
t1321 = t1334 ^ 2;
t1431 = t1320 + t1321;
t1286 = t1431 * qJDD(1);
t1489 = pkin(3) * t1286;
t1486 = pkin(7) * t1066;
t1481 = qJ(6) * t1333;
t1479 = -pkin(4) * t1106 + pkin(8) * t947;
t1478 = qJD(1) * qJD(4);
t1465 = t1517 * t1330;
t1447 = t1286 * t1325;
t1446 = t1286 * t1326;
t1304 = t1331 * t1336 * t1334;
t1291 = qJDD(4) + t1304;
t1445 = t1291 * t1331;
t1444 = t1291 * t1334;
t1292 = qJDD(4) - t1304;
t1443 = t1292 * t1331;
t1442 = t1292 * t1334;
t1439 = t1320 * t1336;
t1438 = t1321 * t1336;
t1098 = t1330 * t1106;
t1436 = t1331 * t1167;
t1099 = t1333 * t1106;
t1158 = t1334 * t1167;
t1433 = -pkin(2) * t1180 + qJ(3) * t1176;
t1151 = (qJD(5) + t1305) * t1272 + t1385;
t1422 = pkin(4) * t1151 + t1098 + t1604;
t1421 = -pkin(4) * t1514 - t1099 + t1570;
t1302 = -t1438 - t1499;
t1230 = t1302 * t1334 - t1445;
t1420 = -pkin(7) * t1230 + t1158;
t1410 = -t1295 * t1332 - t1335 * t1296;
t1062 = -t1145 * t1333 + t1465;
t955 = pkin(5) * t1512 + t971;
t962 = qJ(6) * t1512 + t972;
t1408 = pkin(8) * t1062 + t1330 * t962 + t1333 * t955 + t1554;
t1064 = t1146 * t1333 + t1465;
t1407 = pkin(8) * t1064 + t1554 + t947;
t935 = -t1106 * t1334 + t1331 * t947;
t1406 = -pkin(3) * t935 - t1479;
t1405 = t1325 * t1304;
t1404 = t1326 * t1304;
t1007 = t1064 * t1331 + t1546;
t1403 = -pkin(7) * t1007 + t1060 * t1488 + t1334 * t932;
t1402 = -pkin(3) * t1230 - t1432;
t1401 = -pkin(2) * t1066 + qJ(3) * t1167 - t1486;
t1399 = -pkin(5) * t972 + qJ(6) * t971;
t1288 = qJDD(1) * t1335 - t1332 * t1336;
t1398 = -pkin(6) * t1288 - g(3) * t1332;
t1397 = -pkin(5) * t1517 - qJ(6) * t1145;
t1395 = t1334 * t1137 + t1427;
t1280 = 0.2e1 * t1423 + t1430;
t1393 = pkin(3) * t1280 + t1158;
t1283 = t1313 - 0.2e1 * t1425;
t1392 = pkin(3) * t1283 - t1436;
t1300 = -t1439 - t1499;
t1228 = t1300 * t1331 + t1442;
t1391 = -pkin(7) * t1228 + t1436;
t1067 = -t1331 * t1417 - t1334 * t1432;
t1386 = t1295 * t1335 - t1296 * t1332;
t1384 = -t1421 - t1587;
t1035 = t1151 * t1334 + t1597;
t1383 = -pkin(3) * t1035 - t1422;
t1382 = -pkin(2) * t1230 + qJ(3) * t1283 + t1420;
t929 = t1330 * t971 - t1333 * t972;
t908 = -pkin(4) * t929 - t1399;
t910 = -pkin(8) * t929 + (pkin(5) * t1330 - t1481) * t984;
t1381 = -pkin(7) * t920 - t1331 * t908 + t1334 * t910;
t1006 = t1062 * t1331 + t1546;
t1058 = -t1145 * t1330 - t1464;
t923 = -pkin(8) * t1058 - t1330 * t955 + t1333 * t962;
t973 = -pkin(4) * t1058 - t1397;
t1379 = -pkin(7) * t1006 - t1331 * t973 + t1334 * t923;
t1029 = t1334 * t1518 - t1597;
t1340 = 0.2e1 * qJD(6) * t1274 - t1341;
t967 = -pkin(5) * t1254 + t1340 + t1553;
t938 = -pkin(5) * t1467 + t1333 * t967 - t1605;
t1344 = pkin(5) * t1202 + qJ(6) * t1515 + t1388;
t941 = -t1344 - 0.2e1 * t1475 - t1606;
t1378 = -pkin(7) * t1029 - t1331 * t941 + t1334 * t938;
t968 = (-t1514 - t1254) * pkin(5) + t1340;
t940 = -t1330 * t968 - t1481 * t1514 - t1569;
t1337 = pkin(5) * t1516 + qJ(6) * t1511 - t972;
t942 = -t1337 - t1572;
t1377 = -t1331 * t942 + t1334 * t940 - t1585;
t1018 = t1098 - t1569;
t969 = t999 - t1572;
t1376 = t1334 * t1018 - t1331 * t969 - t1585;
t1022 = t1099 + t1605;
t970 = t1000 + t1606;
t1375 = -pkin(7) * t1035 + t1334 * t1022 - t1331 * t970;
t1374 = -pkin(3) * t1228 - t1417;
t1373 = pkin(7) * t1286 - t1066;
t1372 = -0.2e1 * t1324 + t1390;
t1371 = -pkin(7) * t935 + t1539 * t946;
t1370 = -pkin(3) * t1006 - t1408;
t1369 = -pkin(3) * t1007 - t1407;
t1366 = pkin(4) * t1518 + pkin(5) * t1466 + t1330 * t967 - t1604;
t1365 = -pkin(2) * t1007 + qJ(3) * t1060 + t1403;
t1364 = t1333 * t968 + t1514 * t1568 + t1570;
t1363 = -pkin(2) * t1228 + qJ(3) * t1280 + t1391;
t1362 = pkin(3) * t929 - t1331 * t910 - t1334 * t908;
t1361 = pkin(3) * t1058 - t1331 * t923 - t1334 * t973;
t1360 = -t1331 * t938 - t1334 * t941 + t1607;
t1359 = -t1331 * t940 - t1334 * t942 + t1573;
t1358 = -t1331 * t1018 - t1334 * t969 + t1573;
t1357 = -t1331 * t1022 - t1334 * t970 - t1607;
t1356 = -pkin(2) * t920 + qJ(3) * t929 + t1381;
t1355 = -pkin(2) * t1006 + qJ(3) * t1058 + t1379;
t1354 = -pkin(2) * t1029 + t1378 + t1602;
t1353 = t1377 + t1578;
t1352 = t1376 + t1578;
t1351 = -pkin(2) * t1035 + t1375 - t1602;
t1289 = t1431 * t1336;
t1350 = pkin(2) * t1286 - qJ(3) * t1289 + t1373;
t1349 = -pkin(2) * t935 + qJ(3) * t946 + t1371;
t1347 = -pkin(3) * t1029 - t1366;
t1346 = -t1364 - t1587;
t1338 = t1204 + 0.2e1 * t1474 + t1429;
t1301 = -t1438 + t1499;
t1299 = t1439 - t1499;
t1290 = (-t1320 + t1321) * t1336;
t1287 = qJDD(1) * t1332 + t1335 * t1336;
t1271 = t1431 * t1478;
t1255 = -pkin(6) * t1287 + g(3) * t1335;
t1240 = t1282 * t1331 + t1321 * t1478;
t1239 = t1281 * t1334 + t1320 * t1478;
t1238 = qJDD(4) * t1326 - t1271 * t1325;
t1237 = qJDD(4) * t1325 + t1271 * t1326;
t1235 = -t1302 * t1331 - t1444;
t1234 = -t1301 * t1331 + t1442;
t1233 = (t1282 - t1425) * t1334;
t1232 = t1300 * t1334 - t1443;
t1231 = t1299 * t1334 - t1445;
t1229 = t1301 * t1334 + t1443;
t1227 = t1299 * t1331 + t1444;
t1226 = (-t1281 + t1423) * t1331;
t1217 = -t1289 * t1326 - t1447;
t1216 = -t1289 * t1325 + t1446;
t1206 = -t1280 * t1334 - t1283 * t1331;
t1205 = -t1280 * t1331 + t1283 * t1334;
t1198 = t1239 * t1325 - t1404;
t1197 = t1240 * t1325 + t1404;
t1196 = -t1239 * t1326 - t1405;
t1195 = -t1240 * t1326 + t1405;
t1194 = t1229 * t1325 + t1313 * t1326;
t1193 = t1227 * t1325 - t1326 * t1430;
t1192 = -t1229 * t1326 + t1313 * t1325;
t1191 = -t1227 * t1326 - t1325 * t1430;
t1173 = t1230 * t1325 + t1283 * t1326;
t1172 = t1228 * t1325 + t1280 * t1326;
t1171 = -t1230 * t1326 + t1283 * t1325;
t1170 = -t1228 * t1326 + t1280 * t1325;
t1169 = -t1203 - t1492;
t1168 = -t1204 - t1493;
t1157 = t1205 * t1325 + t1290 * t1326;
t1156 = -t1205 * t1326 + t1290 * t1325;
t1155 = t1372 + t1492;
t1154 = t1338 + t1493;
t1112 = pkin(1) * t1115;
t1097 = pkin(1) * t1322 + qJ(2) * t1415;
t1073 = -qJ(3) * t1235 - t1402;
t1072 = -qJ(3) * t1232 - t1374;
t1071 = -qJ(2) * t1083 + (-pkin(2) * t1325 + qJ(3) * t1326) * t1322;
t1070 = qJ(2) * t1416 + (pkin(2) * t1326 + qJ(3) * t1325 + pkin(1)) * t1322;
t1069 = t1232 * t1495 + t1393;
t1068 = t1235 * t1495 + t1392;
t1047 = pkin(3) * t1289 + t1067;
t1037 = -t1151 * t1331 + t1596;
t1031 = -t1331 * t1518 - t1596;
t1021 = pkin(1) * t1171 + t1382;
t1020 = pkin(1) * t1170 + t1363;
t1009 = t1064 * t1334 - t1550;
t1008 = t1062 * t1334 - t1550;
t1002 = t1066 * t1325 + t1167 * t1326;
t1001 = -t1066 * t1326 + t1167 * t1325;
t996 = -pkin(3) * t1446 - qJ(2) * t1216 + t1047 * t1325;
t995 = -pkin(3) * t1447 + qJ(2) * t1217 - t1047 * t1326;
t994 = pkin(1) * t1083 + t1433;
t993 = pkin(1) * t1216 + t1350;
t981 = t1035 * t1325 - t1598;
t979 = -t1035 * t1326 - t1599;
t977 = t1029 * t1325 + t1598;
t975 = -t1029 * t1326 + t1599;
t974 = -qJ(3) * t1067 + t1491;
t964 = -qJ(2) * t1171 - t1068 * t1325 + t1073 * t1326;
t963 = -qJ(2) * t1170 - t1069 * t1325 + t1072 * t1326;
t957 = t1067 * t1495 + t1490;
t953 = t1007 * t1325 + t1060 * t1326;
t952 = t1006 * t1325 + t1058 * t1326;
t951 = -t1007 * t1326 + t1060 * t1325;
t950 = -t1006 * t1326 + t1058 * t1325;
t949 = -pkin(1) * t1235 + qJ(2) * t1173 + t1068 * t1326 + t1073 * t1325;
t948 = -pkin(1) * t1232 + qJ(2) * t1172 + t1069 * t1326 + t1072 * t1325;
t936 = t1106 * t1331 + t1334 * t947;
t933 = pkin(1) * t1001 + t1401;
t926 = -qJ(3) * t1037 - t1383;
t925 = -t1384 - t1584;
t924 = -qJ(2) * t1001 - t1325 * t957 + t1326 * t974;
t921 = t1331 * t984 + t1334 * t930;
t918 = -t1346 - t1584;
t917 = -qJ(3) * t1031 - t1347;
t916 = t1037 * t1495 + t1357;
t915 = t1325 * t935 + t1326 * t946;
t914 = t1325 * t946 - t1326 * t935;
t913 = t1358 + t1581;
t912 = -pkin(1) * t1067 + qJ(2) * t1002 + t1325 * t974 + t1326 * t957;
t911 = -qJ(3) * t1009 - t1369;
t907 = pkin(1) * t979 + t1351;
t906 = t1352 + t1608;
t905 = t1009 * t1495 - t1510;
t904 = -qJ(3) * t1008 - t1370;
t903 = t1359 + t1581;
t902 = t1031 * t1495 + t1360;
t901 = t1325 * t920 + t1326 * t929;
t900 = t1325 * t929 - t1326 * t920;
t899 = -qJ(3) * t936 - t1406;
t898 = pkin(1) * t951 + t1365;
t897 = t1353 + t1608;
t896 = t1008 * t1495 + t1361;
t895 = pkin(1) * t975 + t1354;
t894 = -qJ(2) * t979 - t1325 * t916 + t1326 * t926;
t893 = -t1325 * t913 + t1326 * t925 - t1603;
t892 = -pkin(1) * t1037 + qJ(2) * t981 + t1325 * t926 + t1326 * t916;
t891 = pkin(1) * t950 + t1355;
t890 = t1325 * t925 + t1326 * t913 + t1593;
t889 = t1495 * t936 + t1526;
t888 = -t1325 * t903 + t1326 * t918 - t1603;
t887 = -qJ(2) * t951 - t1325 * t905 + t1326 * t911;
t886 = -qJ(2) * t975 - t1325 * t902 + t1326 * t917;
t885 = t1325 * t918 + t1326 * t903 + t1593;
t884 = -pkin(1) * t1031 + qJ(2) * t977 + t1325 * t917 + t1326 * t902;
t883 = -pkin(1) * t1009 + qJ(2) * t953 + t1325 * t911 + t1326 * t905;
t882 = -qJ(3) * t921 - t1540;
t881 = pkin(1) * t914 + t1349;
t880 = -qJ(2) * t950 - t1325 * t896 + t1326 * t904;
t879 = -pkin(1) * t1008 + qJ(2) * t952 + t1325 * t904 + t1326 * t896;
t878 = t1495 * t921 + t1362;
t877 = -qJ(2) * t914 - t1325 * t889 + t1326 * t899;
t876 = -pkin(1) * t936 + qJ(2) * t915 + t1325 * t899 + t1326 * t889;
t875 = pkin(1) * t900 + t1356;
t874 = -qJ(2) * t900 - t1325 * t878 + t1326 * t882;
t873 = -pkin(1) * t921 + qJ(2) * t901 + t1325 * t882 + t1326 * t878;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1288, 0, -t1287, 0, t1398, -t1255, -t1386, -pkin(6) * t1386, 0, 0, -t1412, 0, -t1414, 0, t1500, -t1580, t1561, pkin(6) * t1561 + qJ(2) * t1472 - t1332 * t1097, 0, t1412, t1414, 0, 0, 0, -t1537, -t1500, t1580, -pkin(6) * t1537 - t1332 * t1070 + t1335 * t1071, -t1195 * t1332 + t1197 * t1335, -t1156 * t1332 + t1157 * t1335, -t1192 * t1332 + t1194 * t1335, -t1196 * t1332 + t1198 * t1335, -t1191 * t1332 + t1193 * t1335, -t1237 * t1332 + t1238 * t1335, t1335 * t963 - t1332 * t948 - pkin(6) * (t1170 * t1335 + t1172 * t1332), t1335 * t964 - t1332 * t949 - pkin(6) * (t1171 * t1335 + t1173 * t1332), t1335 * t996 - t1332 * t995 - pkin(6) * (t1216 * t1335 + t1217 * t1332), t1335 * t924 - t1332 * t912 - pkin(6) * (t1001 * t1335 + t1002 * t1332), t1530, t1613, t1609, t1555, t1611, t1558, -t1332 * t890 + t1335 * t893 - t1590, t1335 * t894 - t1332 * t892 - pkin(6) * (t1332 * t981 + t1335 * t979), t1335 * t887 - t1332 * t883 - pkin(6) * (t1332 * t953 + t1335 * t951), t1335 * t877 - t1332 * t876 - pkin(6) * (t1332 * t915 + t1335 * t914), t1530, t1609, -t1613, t1558, -t1611, t1555, -t1332 * t885 + t1335 * t888 - t1590, t1335 * t880 - t1332 * t879 - pkin(6) * (t1332 * t952 + t1335 * t950), t1335 * t886 - t1332 * t884 - pkin(6) * (t1332 * t977 + t1335 * t975), t1335 * t874 - t1332 * t873 - pkin(6) * (t1332 * t901 + t1335 * t900); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1287, 0, t1288, 0, t1255, t1398, t1410, pkin(6) * t1410, 0, 0, t1414, 0, -t1412, 0, t1580, t1500, t1562, pkin(6) * t1562 + qJ(2) * t1473 + t1335 * t1097, 0, -t1414, t1412, 0, 0, 0, t1538, -t1580, -t1500, pkin(6) * t1538 + t1335 * t1070 + t1332 * t1071, t1195 * t1335 + t1197 * t1332, t1156 * t1335 + t1157 * t1332, t1192 * t1335 + t1194 * t1332, t1196 * t1335 + t1198 * t1332, t1191 * t1335 + t1193 * t1332, t1237 * t1335 + t1238 * t1332, t1332 * t963 + t1335 * t948 + pkin(6) * (-t1170 * t1332 + t1172 * t1335), t1332 * t964 + t1335 * t949 + pkin(6) * (-t1171 * t1332 + t1173 * t1335), t1332 * t996 + t1335 * t995 + pkin(6) * (-t1216 * t1332 + t1217 * t1335), t1332 * t924 + t1335 * t912 + pkin(6) * (-t1001 * t1332 + t1002 * t1335), t1529, t1614, t1610, t1556, t1612, t1557, t1332 * t893 + t1335 * t890 + t1589, t1332 * t894 + t1335 * t892 + pkin(6) * (-t1332 * t979 + t1335 * t981), t1332 * t887 + t1335 * t883 + pkin(6) * (-t1332 * t951 + t1335 * t953), t1332 * t877 + t1335 * t876 + pkin(6) * (-t1332 * t914 + t1335 * t915), t1529, t1610, -t1614, t1557, -t1612, t1556, t1332 * t888 + t1335 * t885 + t1589, t1332 * t880 + t1335 * t879 + pkin(6) * (-t1332 * t950 + t1335 * t952), t1332 * t886 + t1335 * t884 + pkin(6) * (-t1332 * t975 + t1335 * t977), t1332 * t874 + t1335 * t873 + pkin(6) * (-t1332 * t900 + t1335 * t901); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1295, t1296, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1169, t1168, 0, -t1112, qJDD(1), 0, 0, 0, 0, 0, 0, t1155, t1154, t994, t1233, t1206, t1234, t1226, t1231, 0, t1020, t1021, t993, t933, t1395, t1026, t1574, t1505, t1045, t1504, t906, t907, t898, t881, t1395, t1574, -t1026, t1504, -t1045, t1505, t897, t891, t895, t875; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1336, 0, 0, -g(3), -t1295, 0, 0, 0, -t1285, 0, -t1284, 0, t1246, -t1419, t1115, qJ(2) * t1115, 0, t1285, t1284, 0, 0, 0, -t1083, -t1246, t1419, t1071, t1197, t1157, t1194, t1198, t1193, t1238, t963, t964, t996, t924, t1508, t961, t1592, t1534, t991, t1531, t893, t894, t887, t877, t1508, t1592, -t961, t1531, -t991, t1534, t888, t880, t886, t874; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1336, 0, qJDD(1), 0, g(3), 0, -t1296, 0, 0, 0, t1284, 0, -t1285, 0, t1419, t1246, t1415, t1097, 0, -t1284, t1285, 0, 0, 0, t1416, -t1419, -t1246, t1070, t1195, t1156, t1192, t1196, t1191, t1237, t948, t949, t995, t912, t1509, -t958, t1591, t1533, -t988, t1532, t890, t892, t883, t876, t1509, t1591, t958, t1532, t988, t1533, t885, t879, t884, t873; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1295, t1296, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1169, t1168, 0, -t1112, qJDD(1), 0, 0, 0, 0, 0, 0, t1155, t1154, t994, t1233, t1206, t1234, t1226, t1231, 0, t1020, t1021, t993, t933, t1395, t1026, t1574, t1505, t1045, t1504, t906, t907, t898, t881, t1395, t1574, -t1026, t1504, -t1045, t1505, t897, t891, t895, t875; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1336, 0, 0, -t1322, t1203, 0, 0, -qJDD(1), t1336, 0, 0, 0, t1180, 0, t1322, qJ(3) * t1322, t1304, t1290, t1313, -t1304, -t1430, qJDD(4), t1072, t1073, -t1489, t974, t1136, t1059, t1564, t1132, t1119, -t1519, t925, t926, t911, t899, t1136, t1564, -t1059, -t1519, -t1119, t1132, t918, t904, t917, t882; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1336, 0, qJDD(1), 0, t1322, 0, t1204, 0, 0, -t1336, -qJDD(1), 0, 0, 0, t1176, -t1322, 0, pkin(2) * t1322, -t1240, -t1205, -t1229, -t1239, -t1227, t1271, t1069, t1068, -t1047, t957, -t1368, -t1024, -t1575, -t1507, -t1041, -t1506, t913, t916, t905, t889, -t1368, -t1575, t1024, -t1506, t1041, -t1507, t903, t896, t902, t878; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1203, -t1204, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t1372, t1338, t1433, t1233, t1206, t1234, t1226, t1231, 0, t1363, t1382, t1350, t1401, t1395, t1026, t1574, t1505, t1045, t1504, t1352, t1351, t1365, t1349, t1395, t1574, -t1026, t1504, -t1045, t1505, t1353, t1355, t1354, t1356; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t1180, t1176, 0, t1233, t1206, t1234, t1226, t1231, 0, t1391, t1420, t1373, -t1486, t1395, t1026, t1574, t1505, t1045, t1504, t1376, t1375, t1403, t1371, t1395, t1574, -t1026, t1504, -t1045, t1505, t1377, t1379, t1378, t1381; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1336, 0, 0, 0, -t1180, 0, -t1322, 0, -t1304, -t1290, -t1313, t1304, t1430, -qJDD(4), t1374, t1402, t1489, -t1491, -t1136, -t1059, -t1564, -t1132, -t1119, t1519, t1384, t1383, t1369, t1406, -t1136, -t1564, t1059, t1519, t1119, -t1132, t1346, t1370, t1347, t1540; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1336, qJDD(1), 0, 0, 0, -t1176, t1322, 0, 0, t1240, t1205, t1229, t1239, t1227, -t1271, pkin(7) * t1232 - t1393, pkin(7) * t1235 - t1392, t1047, pkin(7) * t1067 - t1490, t1368, t1024, t1575, t1507, t1041, t1506, -t1358 + t1586, pkin(7) * t1037 - t1357, pkin(7) * t1009 + t1510, pkin(7) * t936 - t1526, t1368, t1575, -t1024, t1506, -t1041, t1507, -t1359 + t1586, pkin(7) * t1008 - t1361, pkin(7) * t1031 - t1360, pkin(7) * t921 - t1362; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1282, -t1280, t1292, t1425, t1299, -t1425, 0, t1167, -t1417, 0, t1137, t1063, t1563, t1342, t1123, t1394, t1018, t1022, t932, -pkin(8) * t946, t1137, t1563, -t1063, t1394, -t1123, t1342, t940, t923, t938, t910; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1423, t1283, t1301, t1281, t1291, -t1423, -t1167, 0, -t1432, 0, -t1223, -t1513, -t1517, t1223, t1145, -t1270, t969, t970, -pkin(4) * t1060, -pkin(4) * t946, -t1223, -t1517, t1513, -t1270, -t1145, t1223, t942, t973, t941, t908; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1304, t1290, t1313, -t1304, -t1430, qJDD(4), t1417, t1432, 0, 0, t1136, t1059, t1564, t1132, t1119, -t1519, t1421, t1422, t1407, t1479, t1136, t1564, -t1059, -t1519, -t1119, t1132, t1364, t1408, t1366, t1348; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1200, -t1514, t1516, t1450, t1241, -t1450, 0, t1106, t999, 0, t1200, t1516, t1514, -t1450, -t1241, t1450, -qJ(6) * t1514, t962, t967, -qJ(6) * t984; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1254, t1518, -t1242, -t1367, t1515, -t1254, -t1106, 0, t1000, 0, t1254, -t1242, -t1518, -t1254, -t1515, -t1367, t968, t955, pkin(5) * t1518, -pkin(5) * t984; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1223, t1513, t1517, -t1223, -t1145, t1270, -t999, -t1000, 0, 0, t1223, t1517, -t1513, t1270, t1145, -t1223, t1337, t1397, t1293 + t1344, t1399; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1200, t1516, t1514, -t1450, -t1241, t1450, 0, t972, -t984, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1223, t1517, -t1513, t1270, t1145, -t1223, -t972, 0, t971, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1254, t1242, t1518, t1254, t1515, t1367, t984, -t971, 0, 0;];
m_new_reg  = t1;