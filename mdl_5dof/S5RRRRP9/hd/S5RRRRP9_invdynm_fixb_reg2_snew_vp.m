% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RRRRP9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RRRRP9_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP9_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP9_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_invdynm_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:06:45
% EndTime: 2019-12-31 22:07:09
% DurationCPUTime: 25.64s
% Computational Cost: add. (78026->672), mult. (156952->816), div. (0->0), fcn. (109516->8), ass. (0->443)
t1413 = sin(qJ(3));
t1417 = cos(qJ(3));
t1414 = sin(qJ(2));
t1512 = qJD(1) * t1414;
t1372 = -t1417 * qJD(2) + t1413 * t1512;
t1373 = qJD(2) * t1413 + t1417 * t1512;
t1412 = sin(qJ(4));
t1416 = cos(qJ(4));
t1332 = -t1372 * t1412 + t1373 * t1416;
t1329 = t1332 ^ 2;
t1418 = cos(qJ(2));
t1511 = qJD(1) * t1418;
t1401 = -qJD(3) + t1511;
t1393 = -qJD(4) + t1401;
t1521 = t1393 ^ 2;
t1261 = t1521 + t1329;
t1404 = qJD(2) * t1512;
t1471 = t1418 * qJDD(1);
t1377 = -t1404 + t1471;
t1369 = -qJDD(3) + t1377;
t1366 = -qJDD(4) + t1369;
t1330 = t1416 * t1372 + t1373 * t1412;
t1489 = t1332 * t1330;
t1535 = -t1489 + t1366;
t1498 = t1535 * t1412;
t1204 = -t1261 * t1416 + t1498;
t1497 = t1535 * t1416;
t1206 = t1261 * t1412 + t1497;
t1122 = t1204 * t1413 - t1206 * t1417;
t1465 = qJD(2) * t1511;
t1472 = t1414 * qJDD(1);
t1442 = t1465 + t1472;
t1426 = -t1413 * qJDD(2) - t1417 * t1442;
t1325 = -t1372 * qJD(3) - t1426;
t1425 = t1417 * qJDD(2) - t1413 * t1442;
t1423 = t1373 * qJD(3) - t1425;
t1233 = -t1330 * qJD(4) + t1416 * t1325 - t1412 * t1423;
t1490 = t1330 * t1393;
t1538 = t1233 + t1490;
t1089 = t1122 * t1418 - t1414 * t1538;
t1144 = t1204 * t1417 + t1206 * t1413;
t1415 = sin(qJ(1));
t1419 = cos(qJ(1));
t1627 = pkin(5) * (t1089 * t1419 - t1144 * t1415);
t1626 = pkin(5) * (t1089 * t1415 + t1144 * t1419);
t1087 = t1122 * t1414 + t1418 * t1538;
t1625 = pkin(1) * t1087;
t1624 = pkin(6) * t1087;
t1623 = pkin(1) * t1144 + pkin(6) * t1089;
t1462 = -t1325 * t1412 - t1416 * t1423;
t1232 = qJD(4) * t1332 - t1462;
t1314 = t1332 * t1393;
t1539 = t1232 - t1314;
t1137 = -t1539 * t1412 + t1416 * t1538;
t1501 = t1538 * t1412;
t1139 = t1539 * t1416 + t1501;
t1072 = t1137 * t1413 + t1139 * t1417;
t1523 = t1330 ^ 2;
t1537 = t1329 - t1523;
t1064 = t1072 * t1418 - t1414 * t1537;
t1068 = -t1137 * t1417 + t1139 * t1413;
t1620 = t1064 * t1415 - t1068 * t1419;
t1619 = t1064 * t1419 + t1068 * t1415;
t1310 = t1523 - t1521;
t1215 = t1310 * t1412 - t1497;
t1219 = t1310 * t1416 + t1498;
t1157 = t1215 * t1413 - t1219 * t1417;
t1192 = t1232 + t1314;
t1106 = t1157 * t1418 + t1192 * t1414;
t1154 = t1215 * t1417 + t1219 * t1413;
t1618 = t1106 * t1415 + t1154 * t1419;
t1617 = t1106 * t1419 - t1154 * t1415;
t1615 = pkin(2) * t1144;
t1614 = pkin(7) * t1144;
t1611 = pkin(2) * t1538 + pkin(7) * t1122;
t1062 = t1072 * t1414 + t1418 * t1537;
t1102 = t1157 * t1414 - t1192 * t1418;
t1606 = pkin(3) * t1204;
t1605 = pkin(8) * t1204;
t1604 = pkin(8) * t1206;
t1311 = t1329 - t1521;
t1536 = -t1489 - t1366;
t1496 = t1536 * t1412;
t1568 = -t1416 * t1311 + t1496;
t1248 = t1416 * t1536;
t1569 = t1311 * t1412 + t1248;
t1582 = t1413 * t1569 + t1417 * t1568;
t1534 = -t1490 + t1233;
t1581 = -t1413 * t1568 + t1417 * t1569;
t1595 = t1414 * t1534 + t1418 * t1581;
t1603 = t1415 * t1595 - t1419 * t1582;
t1602 = t1415 * t1582 + t1419 * t1595;
t1533 = -t1521 - t1523;
t1544 = t1416 * t1533 - t1496;
t1547 = t1412 * t1533 + t1248;
t1563 = -t1413 * t1547 + t1417 * t1544;
t1584 = t1414 * t1563 - t1418 * t1539;
t1601 = pkin(1) * t1584;
t1600 = pkin(6) * t1584;
t1562 = t1413 * t1544 + t1417 * t1547;
t1583 = t1414 * t1539 + t1418 * t1563;
t1597 = -pkin(1) * t1562 + pkin(6) * t1583;
t1596 = t1414 * t1581 - t1418 * t1534;
t1594 = pkin(5) * (t1415 * t1583 - t1419 * t1562);
t1593 = pkin(5) * (t1415 * t1562 + t1419 * t1583);
t1591 = pkin(2) * t1562;
t1590 = pkin(7) * t1562;
t1585 = -pkin(2) * t1539 + pkin(7) * t1563;
t1238 = -t1523 - t1329;
t1580 = pkin(2) * t1238;
t1579 = pkin(3) * t1238;
t1578 = pkin(3) * t1547;
t1577 = pkin(8) * t1544;
t1576 = pkin(8) * t1547;
t1575 = qJ(5) * t1538;
t1574 = t1238 * t1414;
t1573 = t1238 * t1418;
t1391 = g(1) * t1419 + g(2) * t1415;
t1420 = qJD(1) ^ 2;
t1365 = -pkin(1) * t1420 + qJDD(1) * pkin(6) - t1391;
t1519 = pkin(2) * t1418;
t1453 = -pkin(7) * t1414 - t1519;
t1375 = t1453 * qJD(1);
t1515 = g(3) * t1418;
t1524 = qJD(2) ^ 2;
t1304 = -qJDD(2) * pkin(2) - t1524 * pkin(7) + (qJD(1) * t1375 + t1365) * t1414 + t1515;
t1350 = -pkin(3) * t1401 - pkin(8) * t1373;
t1522 = t1372 ^ 2;
t1210 = t1423 * pkin(3) - t1522 * pkin(8) + t1350 * t1373 + t1304;
t1570 = t1232 * pkin(4) + t1210 - t1575;
t1441 = (t1330 * t1412 + t1332 * t1416) * t1393;
t1480 = t1393 * t1412;
t1308 = t1332 * t1480;
t1479 = t1393 * t1416;
t1468 = t1330 * t1479;
t1447 = -t1308 + t1468;
t1528 = t1413 * t1447 + t1417 * t1441;
t1527 = -t1413 * t1441 + t1417 * t1447;
t1542 = -t1366 * t1414 + t1418 * t1527;
t1567 = t1415 * t1542 - t1419 * t1528;
t1443 = t1232 * t1412 - t1468;
t1448 = -t1416 * t1232 - t1330 * t1480;
t1525 = t1413 * t1443 + t1417 * t1448;
t1470 = t1414 * t1489;
t1526 = -t1413 * t1448 + t1417 * t1443;
t1543 = t1418 * t1526 - t1470;
t1566 = t1415 * t1543 - t1419 * t1525;
t1565 = t1415 * t1528 + t1419 * t1542;
t1564 = t1415 * t1525 + t1419 * t1543;
t1546 = t1418 * t1366 + t1414 * t1527;
t1469 = t1418 * t1489;
t1545 = t1414 * t1526 + t1469;
t1484 = t1373 * t1372;
t1437 = -t1369 - t1484;
t1541 = t1413 * t1437;
t1540 = t1417 * t1437;
t1355 = t1372 * t1401;
t1290 = -t1355 + t1325;
t1186 = t1233 * t1412 - t1332 * t1479;
t1187 = t1233 * t1416 + t1308;
t1130 = t1186 * t1417 + t1187 * t1413;
t1133 = -t1186 * t1413 + t1187 * t1417;
t1449 = t1418 * t1133 + t1470;
t1530 = -t1419 * t1130 + t1415 * t1449;
t1529 = t1130 * t1415 + t1419 * t1449;
t1368 = t1373 ^ 2;
t1399 = t1401 ^ 2;
t1520 = pkin(2) * t1414;
t1390 = t1415 * g(1) - t1419 * g(2);
t1364 = qJDD(1) * pkin(1) + t1420 * pkin(6) + t1390;
t1376 = 0.2e1 * t1465 + t1472;
t1451 = -t1377 + t1404;
t1285 = pkin(2) * t1451 - pkin(7) * t1376 - t1364;
t1349 = -g(3) * t1414 + t1418 * t1365;
t1305 = -t1524 * pkin(2) + qJDD(2) * pkin(7) + t1375 * t1511 + t1349;
t1235 = -t1417 * t1285 + t1305 * t1413;
t1163 = t1437 * pkin(3) - pkin(8) * t1290 - t1235;
t1236 = t1413 * t1285 + t1417 * t1305;
t1170 = -t1522 * pkin(3) - pkin(8) * t1423 + t1401 * t1350 + t1236;
t1109 = -t1416 * t1163 + t1170 * t1412;
t1110 = t1412 * t1163 + t1416 * t1170;
t1054 = -t1109 * t1416 + t1110 * t1412;
t1518 = pkin(3) * t1054;
t1188 = t1416 * t1534;
t1194 = (-qJD(4) - t1393) * t1332 + t1462;
t1138 = t1194 * t1412 - t1188;
t1517 = pkin(3) * t1138;
t1516 = pkin(4) * t1416;
t1514 = qJ(5) * t1416;
t1513 = qJD(1) * qJD(2);
t1510 = qJD(5) * t1393;
t1509 = t1054 * t1413;
t1508 = t1054 * t1417;
t1502 = t1534 * t1412;
t1500 = t1210 * t1412;
t1499 = t1210 * t1416;
t1495 = t1304 * t1413;
t1494 = t1304 * t1417;
t1319 = t1369 - t1484;
t1492 = t1319 * t1413;
t1491 = t1319 * t1417;
t1487 = t1364 * t1414;
t1486 = t1364 * t1418;
t1339 = t1376 * t1414;
t1400 = t1418 * t1420 * t1414;
t1388 = -t1400 + qJDD(2);
t1483 = t1388 * t1414;
t1482 = t1388 * t1418;
t1389 = qJDD(2) + t1400;
t1481 = t1389 * t1414;
t1478 = t1401 * t1413;
t1477 = t1401 * t1417;
t1408 = t1414 ^ 2;
t1476 = t1408 * t1420;
t1379 = -0.2e1 * t1510;
t1265 = pkin(4) * t1330 - qJ(5) * t1332;
t1446 = -pkin(4) * t1521 - t1366 * qJ(5) - t1330 * t1265 + t1110;
t1083 = t1379 + t1446;
t1085 = t1366 * pkin(4) - qJ(5) * t1521 + t1265 * t1332 + qJDD(5) + t1109;
t1475 = -pkin(4) * t1085 + qJ(5) * t1083;
t1474 = -pkin(4) * t1534 - qJ(5) * t1192;
t1409 = t1418 ^ 2;
t1473 = t1408 + t1409;
t1467 = t1414 * t1484;
t1466 = t1418 * t1484;
t1464 = -qJ(5) * t1412 - pkin(3);
t1055 = t1109 * t1412 + t1416 * t1110;
t1160 = t1235 * t1413 + t1417 * t1236;
t1348 = t1365 * t1414 + t1515;
t1294 = t1348 * t1414 + t1418 * t1349;
t1461 = -t1390 * t1415 - t1419 * t1391;
t1460 = t1415 * t1400;
t1459 = t1419 * t1400;
t1044 = t1083 * t1412 - t1085 * t1416;
t1458 = pkin(3) * t1044 + t1475;
t1136 = -t1192 * t1412 - t1188;
t1457 = pkin(3) * t1136 + t1474;
t1456 = -t1110 + t1606;
t1454 = -pkin(2) * t1304 + pkin(7) * t1160;
t1383 = qJDD(1) * t1419 - t1415 * t1420;
t1452 = -pkin(5) * t1383 - g(3) * t1415;
t1450 = t1414 * t1133 - t1469;
t1159 = -t1235 * t1417 + t1236 * t1413;
t1293 = t1348 * t1418 - t1349 * t1414;
t1445 = t1390 * t1419 - t1391 * t1415;
t1444 = -t1109 + t1578;
t1326 = -t1399 - t1522;
t1257 = t1326 * t1417 - t1541;
t1356 = t1401 * t1373;
t1287 = t1356 - t1423;
t1440 = pkin(2) * t1287 + pkin(7) * t1257 - t1494;
t1335 = -t1368 - t1399;
t1263 = -t1335 * t1413 + t1491;
t1291 = (qJD(3) - t1401) * t1372 + t1426;
t1439 = pkin(2) * t1291 + pkin(7) * t1263 + t1495;
t1438 = pkin(4) * t1261 - qJ(5) * t1535 + t1446;
t1288 = (-qJD(3) - t1401) * t1373 + t1425;
t1228 = t1288 * t1417 + t1290 * t1413;
t1318 = t1368 + t1522;
t1436 = pkin(2) * t1318 + pkin(7) * t1228 + t1160;
t1435 = t1438 - t1606;
t1045 = t1083 * t1416 + t1085 * t1412;
t1091 = (-pkin(4) * t1393 - 0.2e1 * qJD(5)) * t1332 + t1570;
t1019 = pkin(8) * t1045 + (t1464 - t1516) * t1091;
t1026 = -t1044 * t1413 + t1045 * t1417;
t1027 = -pkin(8) * t1044 + (pkin(4) * t1412 - t1514) * t1091;
t1434 = -pkin(2) * t1091 + pkin(7) * t1026 + t1019 * t1417 + t1027 * t1413;
t1066 = -pkin(4) * t1238 + t1083;
t1067 = -qJ(5) * t1238 + t1085;
t1140 = -t1192 * t1416 + t1502;
t1031 = pkin(8) * t1140 + t1066 * t1416 + t1067 * t1412 - t1579;
t1035 = -pkin(8) * t1136 - t1066 * t1412 + t1067 * t1416;
t1073 = -t1136 * t1413 + t1140 * t1417;
t1433 = pkin(7) * t1073 + t1031 * t1417 + t1035 * t1413 - t1580;
t1142 = t1194 * t1416 + t1502;
t1037 = pkin(8) * t1142 + t1055 - t1579;
t1039 = -pkin(8) * t1138 - t1054;
t1075 = -t1138 * t1413 + t1142 * t1417;
t1432 = pkin(7) * t1075 + t1037 * t1417 + t1039 * t1413 - t1580;
t1421 = 0.2e1 * qJD(5) * t1332 - t1570;
t1076 = pkin(4) * t1314 + t1421 + t1575;
t1040 = -t1604 + t1076 * t1412 + (pkin(3) + t1516) * t1538;
t1047 = -pkin(4) * t1501 + t1076 * t1416 + t1605;
t1431 = t1040 * t1417 + t1047 * t1413 + t1611;
t1077 = (-t1539 + t1314) * pkin(4) + t1421;
t1043 = t1077 * t1416 + t1464 * t1539 + t1577;
t1051 = -t1077 * t1412 - t1514 * t1539 - t1576;
t1430 = t1043 * t1417 + t1051 * t1413 + t1585;
t1093 = -pkin(3) * t1539 - t1499 + t1577;
t1143 = t1500 - t1576;
t1429 = t1093 * t1417 + t1143 * t1413 + t1585;
t1099 = -pkin(3) * t1538 + t1500 + t1604;
t1148 = t1499 - t1605;
t1428 = t1099 * t1417 + t1148 * t1413 - t1611;
t1033 = t1055 * t1417 - t1509;
t1050 = -pkin(3) * t1210 + pkin(8) * t1055;
t1427 = -pkin(2) * t1210 + pkin(7) * t1033 - pkin(8) * t1509 + t1050 * t1417;
t1424 = pkin(4) * t1536 + qJ(5) * t1533 - t1085;
t1422 = t1424 + t1578;
t1406 = t1409 * t1420;
t1398 = -t1406 - t1524;
t1397 = t1406 - t1524;
t1396 = -t1476 - t1524;
t1395 = -t1476 + t1524;
t1385 = -t1406 + t1476;
t1384 = t1406 + t1476;
t1382 = qJDD(1) * t1415 + t1419 * t1420;
t1381 = t1473 * qJDD(1);
t1378 = -0.2e1 * t1404 + t1471;
t1371 = t1418 * t1389;
t1370 = t1473 * t1513;
t1360 = -pkin(5) * t1382 + g(3) * t1419;
t1354 = -t1368 + t1399;
t1353 = -t1399 + t1522;
t1352 = -t1408 * t1513 + t1418 * t1442;
t1351 = -t1377 * t1414 - t1409 * t1513;
t1347 = -t1396 * t1414 - t1482;
t1346 = -t1395 * t1414 + t1371;
t1345 = t1398 * t1418 - t1481;
t1344 = t1397 * t1418 - t1483;
t1343 = t1396 * t1418 - t1483;
t1342 = t1395 * t1418 + t1481;
t1341 = t1398 * t1414 + t1371;
t1340 = t1397 * t1414 + t1482;
t1338 = t1451 * t1418;
t1336 = t1368 - t1522;
t1334 = t1378 * t1418 - t1339;
t1333 = t1376 * t1418 + t1378 * t1414;
t1307 = -pkin(6) * t1343 - t1486;
t1306 = -pkin(6) * t1341 - t1487;
t1302 = (t1372 * t1417 - t1373 * t1413) * t1401;
t1301 = (t1372 * t1413 + t1373 * t1417) * t1401;
t1299 = -pkin(1) * t1343 + t1349;
t1298 = -pkin(1) * t1341 + t1348;
t1289 = t1355 + t1325;
t1286 = t1356 + t1423;
t1284 = pkin(1) * t1378 + pkin(6) * t1345 + t1486;
t1283 = -pkin(1) * t1376 + pkin(6) * t1347 - t1487;
t1280 = t1325 * t1417 + t1373 * t1478;
t1279 = t1325 * t1413 - t1373 * t1477;
t1278 = -t1372 * t1477 + t1413 * t1423;
t1277 = t1372 * t1478 + t1417 * t1423;
t1273 = t1302 * t1418 - t1369 * t1414;
t1272 = t1302 * t1414 + t1369 * t1418;
t1269 = t1353 * t1417 + t1492;
t1268 = -t1354 * t1413 + t1540;
t1267 = t1353 * t1413 - t1491;
t1266 = t1354 * t1417 + t1541;
t1264 = pkin(1) * t1364 + pkin(6) * t1294;
t1262 = t1335 * t1417 + t1492;
t1258 = pkin(1) * t1384 + pkin(6) * t1381 + t1294;
t1256 = t1326 * t1413 + t1540;
t1243 = t1280 * t1418 + t1467;
t1242 = t1278 * t1418 - t1467;
t1241 = t1280 * t1414 - t1466;
t1240 = t1278 * t1414 + t1466;
t1227 = t1287 * t1417 - t1289 * t1413;
t1226 = t1288 * t1413 - t1290 * t1417;
t1225 = t1287 * t1413 + t1289 * t1417;
t1224 = -pkin(7) * t1262 + t1494;
t1223 = t1269 * t1418 - t1286 * t1414;
t1222 = t1268 * t1418 + t1290 * t1414;
t1221 = t1269 * t1414 + t1286 * t1418;
t1220 = t1268 * t1414 - t1290 * t1418;
t1211 = -pkin(7) * t1256 + t1495;
t1209 = t1263 * t1418 - t1291 * t1414;
t1208 = t1263 * t1414 + t1291 * t1418;
t1203 = t1257 * t1418 - t1287 * t1414;
t1202 = t1257 * t1414 + t1287 * t1418;
t1201 = t1227 * t1418 + t1336 * t1414;
t1200 = t1227 * t1414 - t1336 * t1418;
t1177 = t1228 * t1418 - t1318 * t1414;
t1176 = t1228 * t1414 + t1318 * t1418;
t1175 = -pkin(2) * t1262 + t1236;
t1169 = -pkin(2) * t1256 + t1235;
t1150 = t1160 * t1418 + t1304 * t1414;
t1149 = t1160 * t1414 - t1304 * t1418;
t1134 = -pkin(1) * t1208 - t1439;
t1124 = -pkin(1) * t1202 - t1440;
t1119 = -pkin(7) * t1226 - t1159;
t1094 = -pkin(6) * t1208 - t1175 * t1414 + t1224 * t1418;
t1092 = -pkin(6) * t1202 - t1169 * t1414 + t1211 * t1418;
t1086 = -pkin(1) * t1262 + pkin(6) * t1209 + t1175 * t1418 + t1224 * t1414;
t1081 = -pkin(1) * t1176 - t1436;
t1080 = -pkin(1) * t1149 - t1454;
t1079 = -pkin(1) * t1256 + pkin(6) * t1203 + t1169 * t1418 + t1211 * t1414;
t1078 = -pkin(6) * t1176 + t1119 * t1418 + t1226 * t1520;
t1071 = t1138 * t1417 + t1142 * t1413;
t1069 = t1136 * t1417 + t1140 * t1413;
t1061 = -pkin(6) * t1149 + (-pkin(7) * t1418 + t1520) * t1159;
t1060 = t1075 * t1418 + t1574;
t1059 = t1073 * t1418 + t1574;
t1058 = t1075 * t1414 - t1573;
t1057 = t1073 * t1414 - t1573;
t1056 = pkin(6) * t1177 + t1119 * t1414 + (-pkin(1) - t1519) * t1226;
t1053 = -t1456 - t1615;
t1052 = -t1444 - t1591;
t1049 = -pkin(2) * t1071 - t1517;
t1048 = pkin(6) * t1150 + (-pkin(1) + t1453) * t1159;
t1046 = -t1099 * t1413 + t1148 * t1417 - t1614;
t1042 = -t1422 - t1591;
t1041 = -t1093 * t1413 + t1143 * t1417 - t1590;
t1038 = -t1435 + 0.2e1 * t1510 + t1615;
t1036 = -pkin(2) * t1069 - t1457;
t1034 = -t1428 + t1625;
t1032 = t1055 * t1413 + t1508;
t1030 = -t1429 - t1601;
t1029 = t1033 * t1418 + t1210 * t1414;
t1028 = t1033 * t1414 - t1210 * t1418;
t1025 = t1044 * t1417 + t1045 * t1413;
t1024 = t1046 * t1418 - t1053 * t1414 + t1624;
t1023 = -t1043 * t1413 + t1051 * t1417 - t1590;
t1022 = -t1040 * t1413 + t1047 * t1417 + t1614;
t1021 = -pkin(2) * t1032 - t1518;
t1020 = t1041 * t1418 - t1052 * t1414 - t1600;
t1018 = t1046 * t1414 + t1053 * t1418 - t1623;
t1017 = t1026 * t1418 + t1091 * t1414;
t1016 = t1026 * t1414 - t1091 * t1418;
t1015 = t1041 * t1414 + t1052 * t1418 + t1597;
t1014 = -t1430 - t1601;
t1013 = -pkin(7) * t1071 - t1037 * t1413 + t1039 * t1417;
t1012 = -t1431 - t1625;
t1011 = -pkin(7) * t1032 - pkin(8) * t1508 - t1050 * t1413;
t1010 = -pkin(1) * t1058 - t1432;
t1009 = t1023 * t1418 - t1042 * t1414 - t1600;
t1008 = -pkin(7) * t1069 - t1031 * t1413 + t1035 * t1417;
t1007 = t1022 * t1418 - t1038 * t1414 - t1624;
t1006 = t1023 * t1414 + t1042 * t1418 + t1597;
t1005 = -pkin(2) * t1025 - t1458;
t1004 = t1022 * t1414 + t1038 * t1418 + t1623;
t1003 = -pkin(6) * t1058 + t1013 * t1418 - t1049 * t1414;
t1002 = -pkin(1) * t1057 - t1433;
t1001 = -pkin(1) * t1028 - t1427;
t1000 = -pkin(1) * t1071 + pkin(6) * t1060 + t1013 * t1414 + t1049 * t1418;
t999 = -pkin(6) * t1057 + t1008 * t1418 - t1036 * t1414;
t998 = -pkin(1) * t1069 + pkin(6) * t1059 + t1008 * t1414 + t1036 * t1418;
t997 = -pkin(7) * t1025 - t1019 * t1413 + t1027 * t1417;
t996 = -pkin(6) * t1028 + t1011 * t1418 - t1021 * t1414;
t995 = -pkin(1) * t1032 + pkin(6) * t1029 + t1011 * t1414 + t1021 * t1418;
t994 = -pkin(1) * t1016 - t1434;
t993 = -pkin(6) * t1016 - t1005 * t1414 + t1418 * t997;
t992 = -pkin(1) * t1025 + pkin(6) * t1017 + t1005 * t1418 + t1414 * t997;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1383, 0, -t1382, 0, t1452, -t1360, -t1445, -pkin(5) * t1445, t1352 * t1419 - t1460, t1334 * t1419 + t1385 * t1415, t1346 * t1419 + t1415 * t1472, t1351 * t1419 + t1460, t1344 * t1419 + t1415 * t1471, qJDD(2) * t1415 + t1370 * t1419, t1419 * t1306 - t1415 * t1298 - pkin(5) * (t1345 * t1415 + t1378 * t1419), t1419 * t1307 - t1415 * t1299 - pkin(5) * (t1347 * t1415 - t1376 * t1419), t1419 * t1293 - pkin(5) * (t1381 * t1415 + t1384 * t1419), -pkin(5) * (t1294 * t1415 + t1364 * t1419) - (pkin(1) * t1415 - pkin(6) * t1419) * t1293, t1243 * t1419 + t1279 * t1415, t1201 * t1419 + t1225 * t1415, t1222 * t1419 + t1266 * t1415, t1242 * t1419 - t1277 * t1415, t1223 * t1419 + t1267 * t1415, t1273 * t1419 + t1301 * t1415, t1419 * t1092 - t1415 * t1124 - pkin(5) * (t1203 * t1415 - t1256 * t1419), t1419 * t1094 - t1415 * t1134 - pkin(5) * (t1209 * t1415 - t1262 * t1419), t1419 * t1078 - t1415 * t1081 - pkin(5) * (t1177 * t1415 - t1226 * t1419), t1419 * t1061 - t1415 * t1080 - pkin(5) * (t1150 * t1415 - t1159 * t1419), t1529, -t1619, t1602, t1564, -t1617, t1565, t1419 * t1020 - t1415 * t1030 - t1594, t1419 * t1024 - t1415 * t1034 + t1626, t1419 * t1003 - t1415 * t1010 - pkin(5) * (t1060 * t1415 - t1071 * t1419), t1419 * t996 - t1415 * t1001 - pkin(5) * (t1029 * t1415 - t1032 * t1419), t1529, t1602, t1619, t1565, t1617, t1564, t1419 * t1009 - t1415 * t1014 - t1594, t1419 * t999 - t1415 * t1002 - pkin(5) * (t1059 * t1415 - t1069 * t1419), t1419 * t1007 - t1415 * t1012 - t1626, t1419 * t993 - t1415 * t994 - pkin(5) * (t1017 * t1415 - t1025 * t1419); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1382, 0, t1383, 0, t1360, t1452, t1461, pkin(5) * t1461, t1352 * t1415 + t1459, t1334 * t1415 - t1385 * t1419, t1346 * t1415 - t1419 * t1472, t1351 * t1415 - t1459, t1344 * t1415 - t1419 * t1471, -qJDD(2) * t1419 + t1370 * t1415, t1415 * t1306 + t1419 * t1298 + pkin(5) * (t1345 * t1419 - t1378 * t1415), t1415 * t1307 + t1419 * t1299 + pkin(5) * (t1347 * t1419 + t1376 * t1415), t1415 * t1293 + pkin(5) * (t1381 * t1419 - t1384 * t1415), pkin(5) * (t1294 * t1419 - t1364 * t1415) - (-pkin(1) * t1419 - pkin(6) * t1415) * t1293, t1243 * t1415 - t1279 * t1419, t1201 * t1415 - t1225 * t1419, t1222 * t1415 - t1266 * t1419, t1242 * t1415 + t1277 * t1419, t1223 * t1415 - t1267 * t1419, t1273 * t1415 - t1301 * t1419, t1415 * t1092 + t1419 * t1124 + pkin(5) * (t1203 * t1419 + t1256 * t1415), t1415 * t1094 + t1419 * t1134 + pkin(5) * (t1209 * t1419 + t1262 * t1415), t1415 * t1078 + t1419 * t1081 + pkin(5) * (t1177 * t1419 + t1226 * t1415), t1415 * t1061 + t1419 * t1080 + pkin(5) * (t1150 * t1419 + t1159 * t1415), t1530, -t1620, t1603, t1566, -t1618, t1567, t1415 * t1020 + t1419 * t1030 + t1593, t1415 * t1024 + t1419 * t1034 - t1627, t1415 * t1003 + t1419 * t1010 + pkin(5) * (t1060 * t1419 + t1071 * t1415), t1415 * t996 + t1419 * t1001 + pkin(5) * (t1029 * t1419 + t1032 * t1415), t1530, t1603, t1620, t1567, t1618, t1566, t1415 * t1009 + t1419 * t1014 + t1593, t1415 * t999 + t1419 * t1002 + pkin(5) * (t1059 * t1419 + t1069 * t1415), t1415 * t1007 + t1419 * t1012 + t1627, t1415 * t993 + t1419 * t994 + pkin(5) * (t1017 * t1419 + t1025 * t1415); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1390, t1391, 0, 0, t1339, t1333, t1342, -t1338, t1340, 0, t1284, t1283, t1258, t1264, t1241, t1200, t1220, t1240, t1221, t1272, t1079, t1086, t1056, t1048, t1450, -t1062, t1596, t1545, -t1102, t1546, t1015, t1018, t1000, t995, t1450, t1596, t1062, t1546, t1102, t1545, t1006, t998, t1004, t992; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1420, 0, 0, -g(3), -t1390, 0, t1352, t1334, t1346, t1351, t1344, t1370, t1306, t1307, t1293, pkin(6) * t1293, t1243, t1201, t1222, t1242, t1223, t1273, t1092, t1094, t1078, t1061, t1449, -t1064, t1595, t1543, -t1106, t1542, t1020, t1024, t1003, t996, t1449, t1595, t1064, t1542, t1106, t1543, t1009, t999, t1007, t993; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1420, 0, qJDD(1), 0, g(3), 0, -t1391, 0, t1400, -t1385, -t1472, -t1400, -t1471, -qJDD(2), t1298, t1299, 0, pkin(1) * t1293, -t1279, -t1225, -t1266, t1277, -t1267, -t1301, t1124, t1134, t1081, t1080, -t1130, t1068, -t1582, -t1525, -t1154, -t1528, t1030, t1034, t1010, t1001, -t1130, -t1582, -t1068, -t1528, t1154, -t1525, t1014, t1002, t1012, t994; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1390, t1391, 0, 0, t1339, t1333, t1342, -t1338, t1340, 0, t1284, t1283, t1258, t1264, t1241, t1200, t1220, t1240, t1221, t1272, t1079, t1086, t1056, t1048, t1450, -t1062, t1596, t1545, -t1102, t1546, t1015, t1018, t1000, t995, t1450, t1596, t1062, t1546, t1102, t1545, t1006, t998, t1004, t992; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1442, t1378, t1389, -t1465, t1397, t1465, 0, -t1364, t1348, 0, t1280, t1227, t1268, t1278, t1269, t1302, t1211, t1224, t1119, -pkin(7) * t1159, t1133, -t1072, t1581, t1526, -t1157, t1527, t1041, t1046, t1013, t1011, t1133, t1581, t1072, t1527, t1157, t1526, t1023, t1008, t1022, t997; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1404, t1376, t1395, t1377, t1388, -t1404, t1364, 0, t1349, 0, -t1484, -t1336, -t1290, t1484, t1286, t1369, t1169, t1175, -pkin(2) * t1226, -pkin(2) * t1159, -t1489, -t1537, -t1534, t1489, t1192, t1366, t1052, t1053, t1049, t1021, -t1489, -t1534, t1537, t1366, -t1192, t1489, t1042, t1036, t1038, t1005; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1400, t1385, t1472, t1400, t1471, qJDD(2), -t1348, -t1349, 0, 0, t1279, t1225, t1266, -t1277, t1267, t1301, t1440, t1439, t1436, t1454, t1130, -t1068, t1582, t1525, t1154, t1528, t1429, t1428, t1432, t1427, t1130, t1582, t1068, t1528, -t1154, t1525, t1430, t1433, t1431, t1434; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1325, t1287, t1437, -t1355, t1353, t1355, 0, t1304, t1235, 0, t1187, -t1139, t1569, t1443, t1219, t1447, t1143, t1148, t1039, -pkin(8) * t1054, t1187, t1569, t1139, t1447, -t1219, t1443, t1051, t1035, t1047, t1027; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1356, t1289, t1354, -t1423, -t1319, t1356, -t1304, 0, t1236, 0, t1186, t1137, t1568, t1448, t1215, t1441, t1093, t1099, t1037, t1050, t1186, t1568, -t1137, t1441, -t1215, t1448, t1043, t1031, t1040, t1019; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1484, t1336, t1290, -t1484, -t1286, -t1369, -t1235, -t1236, 0, 0, t1489, t1537, t1534, -t1489, -t1192, -t1366, t1444, t1456, t1517, t1518, t1489, t1534, -t1537, -t1366, t1192, -t1489, t1422, t1457, t1379 + t1435, t1458; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1233, -t1539, t1536, -t1490, t1310, t1490, 0, t1210, t1109, 0, t1233, t1536, t1539, t1490, -t1310, -t1490, -qJ(5) * t1539, t1067, t1076, -qJ(5) * t1091; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1314, t1538, -t1311, -t1232, -t1535, t1314, -t1210, 0, t1110, 0, -t1314, -t1311, -t1538, t1314, t1535, -t1232, t1077, t1066, pkin(4) * t1538, -pkin(4) * t1091; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1489, t1537, t1534, -t1489, -t1192, -t1366, -t1109, -t1110, 0, 0, t1489, t1534, -t1537, -t1366, t1192, -t1489, t1424, t1474, t1379 + t1438, t1475; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1233, t1536, t1539, t1490, -t1310, -t1490, 0, t1085, -t1091, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1489, t1534, -t1537, -t1366, t1192, -t1489, -t1085, 0, t1083, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1314, t1311, t1538, -t1314, -t1535, t1232, t1091, -t1083, 0, 0;];
m_new_reg = t1;
