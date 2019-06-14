% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6RPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 01:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6RPRRRP2_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP2_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP2_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_invdynm_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:15:50
% EndTime: 2019-05-06 01:16:16
% DurationCPUTime: 27.53s
% Computational Cost: add. (199780->831), mult. (392361->1074), div. (0->0), fcn. (269695->10), ass. (0->580)
t1525 = sin(qJ(4));
t1529 = cos(qJ(4));
t1526 = sin(qJ(3));
t1643 = qJD(1) * t1526;
t1476 = -t1529 * qJD(3) + t1525 * t1643;
t1477 = qJD(3) * t1525 + t1529 * t1643;
t1524 = sin(qJ(5));
t1528 = cos(qJ(5));
t1422 = t1528 * t1476 + t1477 * t1524;
t1424 = -t1524 * t1476 + t1528 * t1477;
t1358 = t1424 * t1422;
t1511 = qJD(3) * t1643;
t1530 = cos(qJ(3));
t1600 = t1530 * qJDD(1);
t1483 = -t1511 + t1600;
t1471 = -qJDD(4) + t1483;
t1466 = -qJDD(5) + t1471;
t1687 = t1358 + t1466;
t1693 = pkin(5) * t1687;
t1527 = sin(qJ(1));
t1531 = cos(qJ(1));
t1496 = g(1) * t1531 + g(2) * t1527;
t1532 = qJD(1) ^ 2;
t1479 = -pkin(1) * t1532 - t1496;
t1520 = sin(pkin(10));
t1521 = cos(pkin(10));
t1495 = g(1) * t1527 - t1531 * g(2);
t1558 = qJDD(1) * pkin(1) + t1495;
t1426 = t1521 * t1479 + t1520 * t1558;
t1402 = -t1532 * pkin(2) + qJDD(1) * pkin(7) + t1426;
t1666 = pkin(3) * t1530;
t1565 = -pkin(8) * t1526 - t1666;
t1481 = t1565 * qJD(1);
t1651 = g(3) - qJDD(2);
t1510 = t1530 * t1651;
t1674 = qJD(3) ^ 2;
t1346 = (qJD(1) * t1481 + t1402) * t1526 - qJDD(3) * pkin(3) - pkin(8) * t1674 + t1510;
t1642 = qJD(1) * t1530;
t1506 = -qJD(4) + t1642;
t1447 = -pkin(4) * t1506 - pkin(9) * t1477;
t1595 = qJD(3) * t1642;
t1601 = t1526 * qJDD(1);
t1556 = t1595 + t1601;
t1537 = t1529 * qJDD(3) - t1525 * t1556;
t1534 = t1477 * qJD(4) - t1537;
t1673 = t1476 ^ 2;
t1258 = t1534 * pkin(4) - t1673 * pkin(9) + t1477 * t1447 + t1346;
t1425 = t1520 * t1479 - t1521 * t1558;
t1581 = t1425 * t1520 + t1521 * t1426;
t1344 = t1425 * t1521 - t1426 * t1520;
t1629 = t1344 * t1527;
t1692 = t1531 * t1581 + t1629;
t1628 = t1344 * t1531;
t1691 = -t1527 * t1581 + t1628;
t1486 = qJDD(1) * t1520 + t1521 * t1532;
t1487 = qJDD(1) * t1521 - t1520 * t1532;
t1429 = -t1486 * t1527 + t1531 * t1487;
t1456 = -qJ(2) * t1486 + t1521 * t1651;
t1675 = -qJ(2) * t1487 - t1520 * t1651;
t1690 = -pkin(6) * t1429 - t1456 * t1527 + t1531 * t1675;
t1677 = t1531 * t1486 + t1487 * t1527;
t1689 = -pkin(6) * t1677 + t1456 * t1531 + t1527 * t1675;
t1631 = t1687 * t1524;
t1630 = t1687 * t1528;
t1420 = t1422 ^ 2;
t1498 = -qJD(5) + t1506;
t1497 = t1498 ^ 2;
t1335 = -t1497 - t1420;
t1250 = t1335 * t1524 - t1630;
t1245 = pkin(4) * t1250;
t1538 = -t1525 * qJDD(3) - t1529 * t1556;
t1416 = -t1476 * qJD(4) - t1538;
t1310 = -t1422 * qJD(5) + t1528 * t1416 - t1524 * t1534;
t1401 = -qJDD(1) * pkin(2) - t1532 * pkin(7) + t1425;
t1482 = 0.2e1 * t1595 + t1601;
t1563 = -t1483 + t1511;
t1333 = pkin(3) * t1563 - pkin(8) * t1482 + t1401;
t1377 = t1530 * t1402 - t1526 * t1651;
t1347 = -pkin(3) * t1674 + qJDD(3) * pkin(8) + t1481 * t1642 + t1377;
t1262 = t1525 * t1333 + t1529 * t1347;
t1222 = -pkin(4) * t1673 - pkin(9) * t1534 + t1506 * t1447 + t1262;
t1261 = -t1529 * t1333 + t1525 * t1347;
t1458 = t1476 * t1506;
t1367 = -t1458 + t1416;
t1622 = t1477 * t1476;
t1544 = -t1471 - t1622;
t1533 = t1544 * pkin(4) - pkin(9) * t1367 - t1261;
t1153 = t1524 * t1222 - t1528 * t1533;
t1393 = t1422 * t1498;
t1541 = -qJ(6) * t1393 + 0.2e1 * qJD(6) * t1424 + t1153 + t1693;
t1120 = qJ(6) * t1310 + t1541;
t1536 = -t1120 - t1693;
t1688 = t1245 + t1536;
t1582 = t1524 * t1416 + t1528 * t1534;
t1309 = -qJD(5) * t1424 - t1582;
t1380 = -pkin(5) * t1498 - qJ(6) * t1424;
t1180 = -t1309 * pkin(5) - t1420 * qJ(6) + t1424 * t1380 + qJDD(6) + t1258;
t1683 = t1525 * t1544;
t1681 = t1529 * t1544;
t1375 = t1402 * t1526 + t1510;
t1314 = t1526 * t1375 + t1530 * t1377;
t1678 = t1393 + t1310;
t1268 = (qJD(5) + t1498) * t1424 + t1582;
t1421 = t1424 ^ 2;
t1470 = t1477 ^ 2;
t1504 = t1506 ^ 2;
t1271 = -t1393 + t1310;
t1199 = -t1268 * t1524 - t1271 * t1528;
t1201 = -t1268 * t1528 + t1271 * t1524;
t1148 = -t1199 * t1525 + t1201 * t1529;
t1316 = -t1420 - t1421;
t1122 = t1148 * t1526 - t1316 * t1530;
t1672 = pkin(2) * t1122;
t1251 = t1335 * t1528 + t1631;
t1187 = -t1250 * t1525 + t1251 * t1529;
t1267 = (qJD(5) - t1498) * t1424 + t1582;
t1156 = t1187 * t1526 - t1267 * t1530;
t1671 = pkin(2) * t1156;
t1373 = -t1421 - t1497;
t1327 = -t1358 + t1466;
t1633 = t1327 * t1524;
t1285 = t1373 * t1528 + t1633;
t1632 = t1327 * t1528;
t1286 = -t1373 * t1524 + t1632;
t1209 = -t1285 * t1525 + t1286 * t1529;
t1162 = t1209 * t1526 - t1530 * t1678;
t1670 = pkin(2) * t1162;
t1186 = t1250 * t1529 + t1251 * t1525;
t1669 = pkin(3) * t1186;
t1208 = t1285 * t1529 + t1286 * t1525;
t1668 = pkin(3) * t1208;
t1667 = pkin(3) * t1526;
t1154 = t1528 * t1222 + t1524 * t1533;
t1097 = -t1153 * t1528 + t1154 * t1524;
t1665 = pkin(4) * t1097;
t1123 = t1148 * t1530 + t1316 * t1526;
t1146 = t1199 * t1529 + t1201 * t1525;
t1073 = t1123 * t1520 - t1146 * t1521;
t1074 = t1123 * t1521 + t1146 * t1520;
t1664 = pkin(6) * (t1073 * t1531 + t1074 * t1527);
t1157 = t1187 * t1530 + t1267 * t1526;
t1112 = t1157 * t1520 - t1186 * t1521;
t1113 = t1157 * t1521 + t1186 * t1520;
t1663 = pkin(6) * (t1112 * t1531 + t1113 * t1527);
t1163 = t1209 * t1530 + t1526 * t1678;
t1126 = t1163 * t1520 - t1208 * t1521;
t1127 = t1163 * t1521 + t1208 * t1520;
t1662 = pkin(6) * (t1126 * t1531 + t1127 * t1527);
t1660 = pkin(7) * t1122;
t1659 = pkin(7) * t1156;
t1658 = pkin(7) * t1162;
t1657 = pkin(8) * t1146;
t1656 = pkin(8) * t1186;
t1655 = pkin(8) * t1208;
t1654 = pkin(9) * t1199;
t1653 = pkin(9) * t1250;
t1652 = pkin(9) * t1285;
t1650 = qJ(2) * t1073;
t1649 = qJ(2) * t1112;
t1648 = qJ(2) * t1126;
t1644 = qJD(1) * qJD(3);
t1641 = qJD(6) * t1422;
t1639 = t1097 * t1525;
t1638 = t1097 * t1529;
t1637 = t1120 * t1524;
t1636 = t1120 * t1528;
t1635 = t1258 * t1524;
t1634 = t1258 * t1528;
t1627 = t1346 * t1525;
t1626 = t1346 * t1529;
t1397 = t1471 - t1622;
t1625 = t1397 * t1525;
t1624 = t1397 * t1529;
t1438 = t1482 * t1526;
t1505 = t1530 * t1532 * t1526;
t1493 = -t1505 + qJDD(3);
t1619 = t1493 * t1526;
t1618 = t1493 * t1530;
t1494 = qJDD(3) + t1505;
t1617 = t1494 * t1526;
t1616 = t1498 * t1424;
t1615 = t1498 * t1524;
t1614 = t1498 * t1528;
t1613 = t1506 * t1525;
t1612 = t1506 * t1529;
t1516 = t1526 ^ 2;
t1611 = t1516 * t1532;
t1390 = t1526 * t1401;
t1391 = t1530 * t1401;
t1609 = -pkin(2) * t1146 + pkin(7) * t1123;
t1608 = -pkin(3) * t1316 + pkin(8) * t1148;
t1607 = -pkin(2) * t1186 + pkin(7) * t1157;
t1606 = -pkin(2) * t1208 + pkin(7) * t1163;
t1605 = -pkin(3) * t1267 + pkin(8) * t1187;
t1604 = -pkin(3) * t1678 + pkin(8) * t1209;
t1603 = -pkin(2) * t1401 + pkin(7) * t1314;
t1517 = t1530 ^ 2;
t1602 = t1516 + t1517;
t1599 = t1526 * t1358;
t1598 = t1530 * t1358;
t1597 = t1526 * t1622;
t1596 = t1530 * t1622;
t1282 = pkin(4) * t1285;
t1594 = t1282 - t1154;
t1501 = -t1611 - t1674;
t1446 = -t1501 * t1526 - t1618;
t1593 = -pkin(2) * t1482 + pkin(7) * t1446 + t1390;
t1514 = t1517 * t1532;
t1503 = -t1514 - t1674;
t1444 = t1503 * t1530 - t1617;
t1484 = -0.2e1 * t1511 + t1600;
t1592 = pkin(2) * t1484 + pkin(7) * t1444 - t1391;
t1590 = -pkin(1) * t1122 + qJ(2) * t1074;
t1589 = -pkin(1) * t1156 + qJ(2) * t1113;
t1588 = -pkin(1) * t1162 + qJ(2) * t1127;
t1197 = pkin(4) * t1199;
t1104 = -pkin(3) * t1146 - t1197;
t1413 = -0.2e1 * t1641;
t1560 = t1420 * pkin(5) - t1309 * qJ(6) - t1498 * t1380 - t1154;
t1128 = t1413 - t1560;
t1058 = t1128 * t1524 - t1636;
t1119 = pkin(5) * t1120;
t1587 = pkin(4) * t1058 - t1119;
t1586 = -pkin(4) * t1316 + pkin(9) * t1201;
t1585 = -pkin(4) * t1267 + pkin(9) * t1251;
t1584 = -pkin(4) * t1678 + pkin(9) * t1286;
t1098 = t1153 * t1524 + t1528 * t1154;
t1195 = t1261 * t1525 + t1529 * t1262;
t1580 = -t1495 * t1527 - t1531 * t1496;
t1059 = t1128 * t1528 + t1637;
t1026 = -t1058 * t1525 + t1059 * t1529;
t1020 = t1026 * t1530 + t1180 * t1526;
t1025 = t1058 * t1529 + t1059 * t1525;
t1089 = -pkin(5) * t1180 + qJ(6) * t1128;
t1015 = -pkin(4) * t1180 + pkin(9) * t1059 + qJ(6) * t1637 + t1089 * t1528;
t1021 = -pkin(9) * t1058 + qJ(6) * t1636 - t1089 * t1524;
t984 = -pkin(8) * t1025 - t1015 * t1525 + t1021 * t1529;
t999 = -pkin(3) * t1025 - t1587;
t1579 = -pkin(2) * t1025 + pkin(7) * t1020 + t1526 * t984 + t1530 * t999;
t1578 = t1520 * t1505;
t1577 = t1521 * t1505;
t1101 = -pkin(5) * t1316 - qJ(6) * t1268 + t1128;
t1105 = (t1271 + t1310) * qJ(6) + t1541;
t1039 = t1101 * t1528 + t1105 * t1524 + t1586;
t1040 = -t1101 * t1524 + t1105 * t1528 - t1654;
t1004 = -t1039 * t1525 + t1040 * t1529 - t1657;
t1263 = pkin(5) * t1271;
t1091 = t1104 + t1263;
t1576 = t1526 * t1004 + t1530 * t1091 + t1609;
t1050 = t1098 * t1525 + t1638;
t1084 = -pkin(4) * t1258 + pkin(9) * t1098;
t1012 = -pkin(8) * t1050 - pkin(9) * t1638 - t1084 * t1525;
t1028 = -pkin(3) * t1050 - t1665;
t1051 = t1098 * t1529 - t1639;
t1045 = t1051 * t1530 + t1258 * t1526;
t1575 = -pkin(2) * t1050 + pkin(7) * t1045 + t1526 * t1012 + t1530 * t1028;
t1066 = t1098 + t1586;
t1079 = -t1097 - t1654;
t1023 = -t1066 * t1525 + t1079 * t1529 - t1657;
t1574 = t1526 * t1023 + t1530 * t1104 + t1609;
t1149 = -pkin(5) * t1267 + qJ(6) * t1335 - t1180;
t1096 = qJ(6) * t1631 + t1149 * t1528 + t1585;
t1108 = qJ(6) * t1630 - t1149 * t1524 - t1653;
t1042 = -t1096 * t1525 + t1108 * t1529 - t1656;
t1076 = -t1669 - t1688;
t1573 = t1526 * t1042 + t1530 * t1076 + t1607;
t1166 = -qJ(6) * t1373 + t1180;
t1221 = -pkin(5) * t1678 + qJ(6) * t1327;
t1102 = t1166 * t1524 + t1221 * t1528 + t1584;
t1115 = t1166 * t1528 - t1221 * t1524 - t1652;
t1048 = -t1102 * t1525 + t1115 * t1529 - t1655;
t1557 = pkin(5) * t1373 + t1560;
t1551 = t1282 + t1557;
t1083 = t1413 - t1551 - t1668;
t1572 = t1526 * t1048 + t1530 * t1083 + t1606;
t1158 = t1585 - t1634;
t1192 = t1635 - t1653;
t1086 = -t1158 * t1525 + t1192 * t1529 - t1656;
t1562 = t1153 - t1245;
t1100 = t1562 - t1669;
t1571 = t1526 * t1086 + t1530 * t1100 + t1607;
t1160 = t1584 + t1635;
t1202 = t1634 - t1652;
t1094 = -t1160 * t1525 + t1202 * t1529 - t1655;
t1107 = -t1594 - t1668;
t1570 = t1526 * t1094 + t1530 * t1107 + t1606;
t1418 = -t1504 - t1673;
t1330 = t1418 * t1525 + t1681;
t1220 = -pkin(3) * t1330 + t1261;
t1260 = -pkin(8) * t1330 + t1627;
t1331 = t1418 * t1529 - t1683;
t1459 = t1506 * t1477;
t1364 = t1459 - t1534;
t1280 = t1331 * t1530 - t1364 * t1526;
t1569 = -pkin(2) * t1330 + pkin(7) * t1280 + t1530 * t1220 + t1526 * t1260;
t1433 = -t1470 - t1504;
t1336 = t1433 * t1529 + t1625;
t1224 = -pkin(3) * t1336 + t1262;
t1275 = -pkin(8) * t1336 + t1626;
t1337 = -t1433 * t1525 + t1624;
t1368 = (qJD(4) - t1506) * t1476 + t1538;
t1289 = t1337 * t1530 - t1368 * t1526;
t1568 = -pkin(2) * t1336 + pkin(7) * t1289 + t1530 * t1224 + t1526 * t1275;
t1488 = t1602 * qJDD(1);
t1491 = t1514 + t1611;
t1567 = pkin(2) * t1491 + pkin(7) * t1488 + t1314;
t1566 = -pkin(3) * t1346 + pkin(8) * t1195;
t1490 = qJDD(1) * t1531 - t1527 * t1532;
t1564 = -pkin(6) * t1490 - g(3) * t1527;
t1194 = -t1261 * t1529 + t1262 * t1525;
t1313 = t1375 * t1530 - t1377 * t1526;
t1559 = t1495 * t1531 - t1496 * t1527;
t1365 = (-qJD(4) - t1506) * t1477 + t1537;
t1304 = t1365 * t1525 - t1367 * t1529;
t1165 = -pkin(8) * t1304 - t1194;
t1306 = t1365 * t1529 + t1367 * t1525;
t1395 = t1470 + t1673;
t1243 = t1306 * t1530 - t1395 * t1526;
t1555 = pkin(7) * t1243 + t1526 * t1165 + (-pkin(2) - t1666) * t1304;
t1553 = pkin(3) * t1364 + pkin(8) * t1331 - t1626;
t1552 = pkin(3) * t1368 + pkin(8) * t1337 + t1627;
t1550 = t1039 * t1529 + t1040 * t1525 + t1608;
t1549 = t1066 * t1529 + t1079 * t1525 + t1608;
t1548 = t1096 * t1529 + t1108 * t1525 + t1605;
t1547 = t1102 * t1529 + t1115 * t1525 + t1604;
t1546 = t1158 * t1529 + t1192 * t1525 + t1605;
t1545 = t1160 * t1529 + t1202 * t1525 + t1604;
t1543 = pkin(3) * t1395 + pkin(8) * t1306 + t1195;
t1179 = t1195 * t1530 + t1346 * t1526;
t1542 = pkin(7) * t1179 + (-pkin(2) + t1565) * t1194;
t1540 = -pkin(3) * t1180 + pkin(8) * t1026 + t1015 * t1529 + t1021 * t1525;
t1539 = -pkin(3) * t1258 + pkin(8) * t1051 - pkin(9) * t1639 + t1084 * t1529;
t1502 = t1514 - t1674;
t1500 = -t1611 + t1674;
t1492 = -t1514 + t1611;
t1489 = qJDD(1) * t1527 + t1531 * t1532;
t1475 = t1530 * t1494;
t1474 = t1602 * t1644;
t1460 = -pkin(6) * t1489 + g(3) * t1531;
t1453 = -t1470 + t1504;
t1452 = -t1504 + t1673;
t1451 = -t1516 * t1644 + t1530 * t1556;
t1450 = -t1483 * t1526 - t1517 * t1644;
t1449 = qJDD(3) * t1520 + t1474 * t1521;
t1448 = -qJDD(3) * t1521 + t1474 * t1520;
t1445 = -t1500 * t1526 + t1475;
t1443 = t1502 * t1530 - t1619;
t1442 = t1501 * t1530 - t1619;
t1441 = t1500 * t1530 + t1617;
t1440 = t1503 * t1526 + t1475;
t1439 = t1502 * t1526 + t1618;
t1437 = t1563 * t1530;
t1434 = t1470 - t1673;
t1432 = t1488 * t1521 - t1491 * t1520;
t1431 = t1488 * t1520 + t1491 * t1521;
t1428 = t1484 * t1530 - t1438;
t1427 = t1482 * t1530 + t1484 * t1526;
t1414 = 0.2e1 * t1641;
t1411 = t1451 * t1521 - t1578;
t1410 = t1450 * t1521 + t1578;
t1409 = t1451 * t1520 + t1577;
t1408 = t1450 * t1520 - t1577;
t1407 = t1445 * t1521 + t1520 * t1601;
t1406 = t1443 * t1521 + t1520 * t1600;
t1405 = t1445 * t1520 - t1521 * t1601;
t1404 = t1443 * t1520 - t1521 * t1600;
t1388 = -t1421 + t1497;
t1387 = t1420 - t1497;
t1386 = t1446 * t1521 + t1482 * t1520;
t1385 = t1444 * t1521 - t1484 * t1520;
t1384 = t1446 * t1520 - t1482 * t1521;
t1383 = t1444 * t1520 + t1484 * t1521;
t1382 = -pkin(1) * t1486 - t1426;
t1381 = pkin(1) * t1487 - t1425;
t1379 = (t1476 * t1529 - t1477 * t1525) * t1506;
t1378 = (t1476 * t1525 + t1477 * t1529) * t1506;
t1376 = t1428 * t1521 + t1492 * t1520;
t1374 = t1428 * t1520 - t1492 * t1521;
t1366 = t1458 + t1416;
t1363 = t1459 + t1534;
t1362 = t1416 * t1529 + t1477 * t1613;
t1361 = t1416 * t1525 - t1477 * t1612;
t1360 = -t1476 * t1612 + t1525 * t1534;
t1359 = t1476 * t1613 + t1529 * t1534;
t1357 = t1379 * t1530 - t1471 * t1526;
t1356 = t1379 * t1526 + t1471 * t1530;
t1354 = t1421 - t1420;
t1353 = -pkin(7) * t1442 + t1391;
t1352 = t1452 * t1529 + t1625;
t1351 = -t1453 * t1525 + t1681;
t1350 = -pkin(7) * t1440 + t1390;
t1349 = t1452 * t1525 - t1624;
t1348 = t1453 * t1529 + t1683;
t1341 = pkin(1) * t1344;
t1339 = -pkin(2) * t1442 + t1377;
t1338 = -pkin(2) * t1440 + t1375;
t1332 = pkin(1) * t1651 + qJ(2) * t1581;
t1322 = (t1422 * t1528 - t1424 * t1524) * t1498;
t1321 = (t1422 * t1524 + t1424 * t1528) * t1498;
t1320 = t1362 * t1530 + t1597;
t1319 = t1360 * t1530 - t1597;
t1318 = t1362 * t1526 - t1596;
t1317 = t1360 * t1526 + t1596;
t1305 = t1364 * t1529 - t1366 * t1525;
t1303 = t1364 * t1525 + t1366 * t1529;
t1301 = t1357 * t1521 + t1378 * t1520;
t1300 = t1357 * t1520 - t1378 * t1521;
t1299 = pkin(1) * t1383 + t1592;
t1298 = pkin(1) * t1384 + t1593;
t1297 = t1352 * t1530 - t1363 * t1526;
t1296 = t1351 * t1530 + t1367 * t1526;
t1295 = t1352 * t1526 + t1363 * t1530;
t1294 = t1351 * t1526 - t1367 * t1530;
t1293 = t1387 * t1528 + t1633;
t1292 = -t1388 * t1524 - t1630;
t1291 = t1387 * t1524 - t1632;
t1290 = t1388 * t1528 - t1631;
t1288 = t1337 * t1526 + t1368 * t1530;
t1284 = -qJ(2) * t1431 + t1313 * t1521;
t1283 = qJ(2) * t1432 + t1313 * t1520;
t1279 = t1331 * t1526 + t1364 * t1530;
t1277 = t1314 * t1521 + t1401 * t1520;
t1276 = t1314 * t1520 - t1401 * t1521;
t1274 = t1305 * t1530 + t1434 * t1526;
t1273 = t1305 * t1526 - t1434 * t1530;
t1259 = pkin(1) * t1431 + t1567;
t1256 = t1310 * t1528 + t1424 * t1615;
t1255 = t1310 * t1524 - t1424 * t1614;
t1254 = -t1309 * t1524 - t1422 * t1614;
t1253 = t1309 * t1528 - t1422 * t1615;
t1249 = t1320 * t1521 + t1361 * t1520;
t1248 = t1319 * t1521 - t1359 * t1520;
t1247 = t1320 * t1520 - t1361 * t1521;
t1246 = t1319 * t1520 + t1359 * t1521;
t1242 = t1306 * t1526 + t1395 * t1530;
t1240 = -t1321 * t1525 + t1322 * t1529;
t1239 = t1321 * t1529 + t1322 * t1525;
t1238 = -qJ(2) * t1384 - t1339 * t1520 + t1353 * t1521;
t1237 = -qJ(2) * t1383 - t1338 * t1520 + t1350 * t1521;
t1236 = t1297 * t1521 + t1349 * t1520;
t1235 = t1296 * t1521 + t1348 * t1520;
t1234 = t1297 * t1520 - t1349 * t1521;
t1233 = t1296 * t1520 - t1348 * t1521;
t1232 = t1240 * t1530 - t1466 * t1526;
t1231 = t1240 * t1526 + t1466 * t1530;
t1230 = t1289 * t1521 + t1336 * t1520;
t1229 = t1289 * t1520 - t1336 * t1521;
t1228 = -pkin(1) * t1442 + qJ(2) * t1386 + t1339 * t1521 + t1353 * t1520;
t1227 = -pkin(1) * t1440 + qJ(2) * t1385 + t1338 * t1521 + t1350 * t1520;
t1226 = t1280 * t1521 + t1330 * t1520;
t1225 = t1280 * t1520 - t1330 * t1521;
t1215 = -t1291 * t1525 + t1293 * t1529;
t1214 = -t1290 * t1525 + t1292 * t1529;
t1213 = t1291 * t1529 + t1293 * t1525;
t1212 = t1290 * t1529 + t1292 * t1525;
t1211 = t1274 * t1521 + t1303 * t1520;
t1210 = t1274 * t1520 - t1303 * t1521;
t1205 = t1243 * t1521 + t1304 * t1520;
t1204 = t1243 * t1520 - t1304 * t1521;
t1203 = pkin(1) * t1276 + t1603;
t1200 = -t1267 * t1528 - t1524 * t1678;
t1198 = -t1267 * t1524 + t1528 * t1678;
t1191 = -t1255 * t1525 + t1256 * t1529;
t1190 = -t1253 * t1525 + t1254 * t1529;
t1189 = t1255 * t1529 + t1256 * t1525;
t1188 = t1253 * t1529 + t1254 * t1525;
t1183 = -pkin(2) * t1288 - t1552;
t1182 = -pkin(2) * t1279 - t1553;
t1181 = -qJ(2) * t1276 - (pkin(2) * t1520 - pkin(7) * t1521) * t1313;
t1178 = t1195 * t1526 - t1346 * t1530;
t1177 = t1232 * t1521 + t1239 * t1520;
t1176 = t1232 * t1520 - t1239 * t1521;
t1175 = t1191 * t1530 + t1599;
t1174 = t1190 * t1530 - t1599;
t1173 = t1191 * t1526 - t1598;
t1172 = t1190 * t1526 + t1598;
t1170 = t1215 * t1530 - t1268 * t1526;
t1169 = t1214 * t1530 + t1271 * t1526;
t1168 = t1215 * t1526 + t1268 * t1530;
t1167 = t1214 * t1526 - t1271 * t1530;
t1159 = qJ(2) * t1277 - (-pkin(2) * t1521 - pkin(7) * t1520 - pkin(1)) * t1313;
t1151 = -pkin(7) * t1288 - t1224 * t1526 + t1275 * t1530;
t1150 = -pkin(7) * t1279 - t1220 * t1526 + t1260 * t1530;
t1147 = -t1198 * t1525 + t1200 * t1529;
t1145 = t1198 * t1529 + t1200 * t1525;
t1142 = -pkin(2) * t1242 - t1543;
t1141 = t1170 * t1521 + t1213 * t1520;
t1140 = t1169 * t1521 + t1212 * t1520;
t1139 = t1170 * t1520 - t1213 * t1521;
t1138 = t1169 * t1520 - t1212 * t1521;
t1137 = t1147 * t1530 + t1354 * t1526;
t1136 = t1147 * t1526 - t1354 * t1530;
t1135 = t1179 * t1521 + t1194 * t1520;
t1134 = t1179 * t1520 - t1194 * t1521;
t1133 = -pkin(7) * t1242 + t1165 * t1530 + t1304 * t1667;
t1132 = t1175 * t1521 + t1189 * t1520;
t1131 = t1174 * t1521 + t1188 * t1520;
t1130 = t1175 * t1520 - t1189 * t1521;
t1129 = t1174 * t1520 - t1188 * t1521;
t1125 = pkin(1) * t1126;
t1118 = -pkin(2) * t1178 - t1566;
t1117 = -t1176 * t1527 + t1177 * t1531;
t1116 = t1176 * t1531 + t1177 * t1527;
t1114 = pkin(1) * t1229 + t1568;
t1111 = pkin(1) * t1112;
t1109 = pkin(1) * t1225 + t1569;
t1095 = -pkin(7) * t1178 + (-pkin(8) * t1530 + t1667) * t1194;
t1092 = -qJ(2) * t1229 + t1151 * t1521 - t1183 * t1520;
t1088 = pkin(1) * t1204 + t1555;
t1087 = -qJ(2) * t1225 + t1150 * t1521 - t1182 * t1520;
t1081 = -pkin(1) * t1288 + qJ(2) * t1230 + t1151 * t1520 + t1183 * t1521;
t1080 = -pkin(1) * t1279 + qJ(2) * t1226 + t1150 * t1520 + t1182 * t1521;
t1078 = t1137 * t1521 + t1145 * t1520;
t1077 = t1137 * t1520 - t1145 * t1521;
t1072 = pkin(1) * t1073;
t1070 = -t1139 * t1527 + t1141 * t1531;
t1069 = -t1138 * t1527 + t1140 * t1531;
t1068 = t1139 * t1531 + t1141 * t1527;
t1067 = t1138 * t1531 + t1140 * t1527;
t1065 = -t1130 * t1527 + t1132 * t1531;
t1064 = -t1129 * t1527 + t1131 * t1531;
t1063 = t1130 * t1531 + t1132 * t1527;
t1062 = t1129 * t1531 + t1131 * t1527;
t1060 = pkin(6) * (-t1126 * t1527 + t1127 * t1531);
t1057 = -t1545 - t1670;
t1055 = pkin(6) * (-t1112 * t1527 + t1113 * t1531);
t1054 = -qJ(2) * t1204 + t1133 * t1521 - t1142 * t1520;
t1053 = -t1546 - t1671;
t1052 = -pkin(1) * t1242 + qJ(2) * t1205 + t1133 * t1520 + t1142 * t1521;
t1047 = pkin(1) * t1134 + t1542;
t1044 = t1051 * t1526 - t1258 * t1530;
t1038 = t1094 * t1530 - t1107 * t1526 - t1658;
t1037 = -qJ(2) * t1134 + t1095 * t1521 - t1118 * t1520;
t1036 = t1086 * t1530 - t1100 * t1526 - t1659;
t1035 = -t1547 - t1670;
t1034 = -t1077 * t1527 + t1078 * t1531;
t1033 = t1077 * t1531 + t1078 * t1527;
t1031 = pkin(6) * (-t1073 * t1527 + t1074 * t1531);
t1030 = -t1548 - t1671;
t1029 = -pkin(1) * t1178 + qJ(2) * t1135 + t1095 * t1520 + t1118 * t1521;
t1019 = t1026 * t1526 - t1180 * t1530;
t1017 = t1125 + t1570;
t1016 = t1048 * t1530 - t1083 * t1526 - t1658;
t1014 = t1111 + t1571;
t1013 = t1042 * t1530 - t1076 * t1526 - t1659;
t1010 = t1045 * t1521 + t1050 * t1520;
t1009 = t1045 * t1520 - t1050 * t1521;
t1008 = -t1549 - t1672;
t1007 = t1125 + t1572;
t1006 = t1038 * t1521 - t1057 * t1520 - t1648;
t1005 = t1023 * t1530 - t1104 * t1526 - t1660;
t1002 = t1111 + t1573;
t1001 = t1038 * t1520 + t1057 * t1521 + t1588;
t1000 = t1036 * t1521 - t1053 * t1520 - t1649;
t997 = t1036 * t1520 + t1053 * t1521 + t1589;
t996 = -t1550 - t1672;
t995 = t1020 * t1521 + t1025 * t1520;
t994 = t1020 * t1520 - t1025 * t1521;
t993 = -pkin(2) * t1044 - t1539;
t992 = t1004 * t1530 - t1091 * t1526 - t1660;
t991 = t1072 + t1574;
t990 = t1016 * t1521 - t1035 * t1520 - t1648;
t989 = t1016 * t1520 + t1035 * t1521 + t1588;
t988 = t1013 * t1521 - t1030 * t1520 - t1649;
t987 = t1013 * t1520 + t1030 * t1521 + t1589;
t986 = t1072 + t1576;
t985 = -pkin(7) * t1044 + t1012 * t1530 - t1028 * t1526;
t982 = t1005 * t1521 - t1008 * t1520 - t1650;
t981 = t1005 * t1520 + t1008 * t1521 + t1590;
t980 = -pkin(2) * t1019 - t1540;
t979 = -t1520 * t996 + t1521 * t992 - t1650;
t978 = t1520 * t992 + t1521 * t996 + t1590;
t977 = pkin(1) * t1009 + t1575;
t976 = -pkin(7) * t1019 - t1526 * t999 + t1530 * t984;
t975 = -qJ(2) * t1009 - t1520 * t993 + t1521 * t985;
t974 = -pkin(1) * t1044 + qJ(2) * t1010 + t1520 * t985 + t1521 * t993;
t973 = pkin(1) * t994 + t1579;
t972 = -qJ(2) * t994 - t1520 * t980 + t1521 * t976;
t971 = -pkin(1) * t1019 + qJ(2) * t995 + t1520 * t976 + t1521 * t980;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1490, 0, -t1489, 0, t1564, -t1460, -t1559, -pkin(6) * t1559, 0, 0, t1429, 0, -t1677, 0, t1690, -t1689, t1691, pkin(6) * t1691 + qJ(2) * t1628 - t1527 * t1332, -t1409 * t1527 + t1411 * t1531, -t1374 * t1527 + t1376 * t1531, -t1405 * t1527 + t1407 * t1531, -t1408 * t1527 + t1410 * t1531, -t1404 * t1527 + t1406 * t1531, -t1448 * t1527 + t1449 * t1531, t1531 * t1237 - t1527 * t1227 - pkin(6) * (t1383 * t1531 + t1385 * t1527), t1531 * t1238 - t1527 * t1228 - pkin(6) * (t1384 * t1531 + t1386 * t1527), t1531 * t1284 - t1527 * t1283 - pkin(6) * (t1431 * t1531 + t1432 * t1527), t1531 * t1181 - t1527 * t1159 - pkin(6) * (t1276 * t1531 + t1277 * t1527), -t1247 * t1527 + t1249 * t1531, -t1210 * t1527 + t1211 * t1531, -t1233 * t1527 + t1235 * t1531, -t1246 * t1527 + t1248 * t1531, -t1234 * t1527 + t1236 * t1531, -t1300 * t1527 + t1301 * t1531, t1531 * t1087 - t1527 * t1080 - pkin(6) * (t1225 * t1531 + t1226 * t1527), t1531 * t1092 - t1527 * t1081 - pkin(6) * (t1229 * t1531 + t1230 * t1527), t1531 * t1054 - t1527 * t1052 - pkin(6) * (t1204 * t1531 + t1205 * t1527), t1531 * t1037 - t1527 * t1029 - pkin(6) * (t1134 * t1531 + t1135 * t1527), t1065, t1034, t1069, t1064, t1070, t1117, t1000 * t1531 - t1527 * t997 - t1663, -t1001 * t1527 + t1006 * t1531 - t1662, -t1527 * t981 + t1531 * t982 - t1664, t1531 * t975 - t1527 * t974 - pkin(6) * (t1009 * t1531 + t1010 * t1527), t1065, t1034, t1069, t1064, t1070, t1117, -t1527 * t987 + t1531 * t988 - t1663, -t1527 * t989 + t1531 * t990 - t1662, -t1527 * t978 + t1531 * t979 - t1664, t1531 * t972 - t1527 * t971 - pkin(6) * (t1527 * t995 + t1531 * t994); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1489, 0, t1490, 0, t1460, t1564, t1580, pkin(6) * t1580, 0, 0, t1677, 0, t1429, 0, t1689, t1690, t1692, pkin(6) * t1692 + qJ(2) * t1629 + t1531 * t1332, t1409 * t1531 + t1411 * t1527, t1374 * t1531 + t1376 * t1527, t1405 * t1531 + t1407 * t1527, t1408 * t1531 + t1410 * t1527, t1404 * t1531 + t1406 * t1527, t1448 * t1531 + t1449 * t1527, t1527 * t1237 + t1531 * t1227 + pkin(6) * (-t1383 * t1527 + t1385 * t1531), t1527 * t1238 + t1531 * t1228 + pkin(6) * (-t1384 * t1527 + t1386 * t1531), t1527 * t1284 + t1531 * t1283 + pkin(6) * (-t1431 * t1527 + t1432 * t1531), t1527 * t1181 + t1531 * t1159 + pkin(6) * (-t1276 * t1527 + t1277 * t1531), t1247 * t1531 + t1249 * t1527, t1210 * t1531 + t1211 * t1527, t1233 * t1531 + t1235 * t1527, t1246 * t1531 + t1248 * t1527, t1234 * t1531 + t1236 * t1527, t1300 * t1531 + t1301 * t1527, t1527 * t1087 + t1531 * t1080 + pkin(6) * (-t1225 * t1527 + t1226 * t1531), t1527 * t1092 + t1531 * t1081 + pkin(6) * (-t1229 * t1527 + t1230 * t1531), t1527 * t1054 + t1531 * t1052 + pkin(6) * (-t1204 * t1527 + t1205 * t1531), t1527 * t1037 + t1531 * t1029 + pkin(6) * (-t1134 * t1527 + t1135 * t1531), t1063, t1033, t1067, t1062, t1068, t1116, t1000 * t1527 + t1531 * t997 + t1055, t1001 * t1531 + t1006 * t1527 + t1060, t1527 * t982 + t1531 * t981 + t1031, t1527 * t975 + t1531 * t974 + pkin(6) * (-t1009 * t1527 + t1010 * t1531), t1063, t1033, t1067, t1062, t1068, t1116, t1527 * t988 + t1531 * t987 + t1055, t1527 * t990 + t1531 * t989 + t1060, t1527 * t979 + t1531 * t978 + t1031, t1527 * t972 + t1531 * t971 + pkin(6) * (-t1527 * t994 + t1531 * t995); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1495, t1496, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1381, t1382, 0, -t1341, t1438, t1427, t1441, -t1437, t1439, 0, t1299, t1298, t1259, t1203, t1318, t1273, t1294, t1317, t1295, t1356, t1109, t1114, t1088, t1047, t1173, t1136, t1167, t1172, t1168, t1231, t1014, t1017, t991, t977, t1173, t1136, t1167, t1172, t1168, t1231, t1002, t1007, t986, t973; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1532, 0, 0, -g(3), -t1495, 0, 0, 0, t1487, 0, -t1486, 0, t1675, -t1456, t1344, qJ(2) * t1344, t1411, t1376, t1407, t1410, t1406, t1449, t1237, t1238, t1284, t1181, t1249, t1211, t1235, t1248, t1236, t1301, t1087, t1092, t1054, t1037, t1132, t1078, t1140, t1131, t1141, t1177, t1000, t1006, t982, t975, t1132, t1078, t1140, t1131, t1141, t1177, t988, t990, t979, t972; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1532, 0, qJDD(1), 0, g(3), 0, -t1496, 0, 0, 0, t1486, 0, t1487, 0, t1456, t1675, t1581, t1332, t1409, t1374, t1405, t1408, t1404, t1448, t1227, t1228, t1283, t1159, t1247, t1210, t1233, t1246, t1234, t1300, t1080, t1081, t1052, t1029, t1130, t1077, t1138, t1129, t1139, t1176, t997, t1001, t981, t974, t1130, t1077, t1138, t1129, t1139, t1176, t987, t989, t978, t971; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1495, t1496, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1381, t1382, 0, -t1341, t1438, t1427, t1441, -t1437, t1439, 0, t1299, t1298, t1259, t1203, t1318, t1273, t1294, t1317, t1295, t1356, t1109, t1114, t1088, t1047, t1173, t1136, t1167, t1172, t1168, t1231, t1014, t1017, t991, t977, t1173, t1136, t1167, t1172, t1168, t1231, t1002, t1007, t986, t973; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1532, 0, 0, -t1651, t1425, 0, t1451, t1428, t1445, t1450, t1443, t1474, t1350, t1353, t1313, pkin(7) * t1313, t1320, t1274, t1296, t1319, t1297, t1357, t1150, t1151, t1133, t1095, t1175, t1137, t1169, t1174, t1170, t1232, t1036, t1038, t1005, t985, t1175, t1137, t1169, t1174, t1170, t1232, t1013, t1016, t992, t976; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1532, 0, qJDD(1), 0, t1651, 0, t1426, 0, t1505, -t1492, -t1601, -t1505, -t1600, -qJDD(3), t1338, t1339, 0, pkin(2) * t1313, -t1361, -t1303, -t1348, t1359, -t1349, -t1378, t1182, t1183, t1142, t1118, -t1189, -t1145, -t1212, -t1188, -t1213, -t1239, t1053, t1057, t1008, t993, -t1189, -t1145, -t1212, -t1188, -t1213, -t1239, t1030, t1035, t996, t980; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1425, -t1426, 0, 0, t1438, t1427, t1441, -t1437, t1439, 0, t1592, t1593, t1567, t1603, t1318, t1273, t1294, t1317, t1295, t1356, t1569, t1568, t1555, t1542, t1173, t1136, t1167, t1172, t1168, t1231, t1571, t1570, t1574, t1575, t1173, t1136, t1167, t1172, t1168, t1231, t1573, t1572, t1576, t1579; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1556, t1484, t1494, -t1595, t1502, t1595, 0, t1401, t1375, 0, t1362, t1305, t1351, t1360, t1352, t1379, t1260, t1275, t1165, -pkin(8) * t1194, t1191, t1147, t1214, t1190, t1215, t1240, t1086, t1094, t1023, t1012, t1191, t1147, t1214, t1190, t1215, t1240, t1042, t1048, t1004, t984; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1511, t1482, t1500, t1483, t1493, -t1511, -t1401, 0, t1377, 0, -t1622, -t1434, -t1367, t1622, t1363, t1471, t1220, t1224, -pkin(3) * t1304, -pkin(3) * t1194, -t1358, -t1354, -t1271, t1358, t1268, t1466, t1100, t1107, t1104, t1028, -t1358, -t1354, -t1271, t1358, t1268, t1466, t1076, t1083, t1091, t999; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1505, t1492, t1601, t1505, t1600, qJDD(3), -t1375, -t1377, 0, 0, t1361, t1303, t1348, -t1359, t1349, t1378, t1553, t1552, t1543, t1566, t1189, t1145, t1212, t1188, t1213, t1239, t1546, t1545, t1549, t1539, t1189, t1145, t1212, t1188, t1213, t1239, t1548, t1547, t1550, t1540; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1416, t1364, t1544, -t1458, t1452, t1458, 0, t1346, t1261, 0, t1256, t1200, t1292, t1254, t1293, t1322, t1192, t1202, t1079, -pkin(9) * t1097, t1256, t1200, t1292, t1254, t1293, t1322, t1108, t1115, t1040, t1021; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1459, t1366, t1453, -t1534, -t1397, t1459, -t1346, 0, t1262, 0, t1255, t1198, t1290, t1253, t1291, t1321, t1158, t1160, t1066, t1084, t1255, t1198, t1290, t1253, t1291, t1321, t1096, t1102, t1039, t1015; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1622, t1434, t1367, -t1622, -t1363, -t1471, -t1261, -t1262, 0, 0, t1358, t1354, t1271, -t1358, -t1268, -t1466, -t1562, t1594, t1197, t1665, t1358, t1354, t1271, -t1358, -t1268, -t1466, t1688, t1414 + t1551, -t1263 + t1197, t1587; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1310, -t1267, -t1687, -t1393, t1387, t1393, 0, t1258, t1153, 0, t1310, -t1267, -t1687, -t1393, t1387, t1393, qJ(6) * t1687, t1166, t1105, qJ(6) * t1120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1616, t1678, t1388, t1309, -t1327, t1616, -t1258, 0, t1154, 0, -t1616, t1678, t1388, t1309, -t1327, t1616, t1149, t1221, t1101, t1089; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1358, t1354, t1271, -t1358, -t1268, -t1466, -t1153, -t1154, 0, 0, t1358, t1354, t1271, -t1358, -t1268, -t1466, t1536, t1414 + t1557, -t1263, -t1119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1310, -t1267, -t1687, -t1393, t1387, t1393, 0, t1180, t1120, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1616, t1678, t1388, t1309, -t1327, t1616, -t1180, 0, t1128, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1358, t1354, t1271, -t1358, -t1268, -t1466, -t1120, -t1128, 0, 0;];
m_new_reg  = t1;