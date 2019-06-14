% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6RRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 22:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6RRPRRR7_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR7_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR7_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_invdynm_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 22:21:55
% EndTime: 2019-05-06 22:22:29
% DurationCPUTime: 37.03s
% Computational Cost: add. (263279->932), mult. (536702->1167), div. (0->0), fcn. (361056->10), ass. (0->593)
t1538 = sin(qJ(2));
t1531 = t1538 ^ 2;
t1545 = qJD(1) ^ 2;
t1642 = t1531 * t1545;
t1693 = qJD(2) ^ 2;
t1509 = t1642 + t1693;
t1543 = cos(qJ(2));
t1709 = t1543 * t1545;
t1514 = t1538 * t1709;
t1505 = qJDD(2) - t1514;
t1639 = t1543 * t1505;
t1459 = -t1509 * t1538 + t1639;
t1673 = qJD(2) * t1543;
t1520 = qJD(1) * t1673;
t1632 = t1538 * qJDD(1);
t1495 = 0.2e1 * t1520 + t1632;
t1539 = sin(qJ(1));
t1544 = cos(qJ(1));
t1728 = pkin(6) * (t1459 * t1544 - t1495 * t1539);
t1727 = pkin(6) * (t1459 * t1539 + t1495 * t1544);
t1537 = sin(qJ(4));
t1631 = qJDD(2) - qJDD(4);
t1542 = cos(qJ(4));
t1481 = (t1537 * t1538 + t1542 * t1543) * qJD(1);
t1674 = qJD(1) * t1538;
t1483 = -t1537 * t1543 * qJD(1) + t1542 * t1674;
t1655 = t1483 * t1481;
t1723 = -t1631 - t1655;
t1726 = t1537 * t1723;
t1725 = t1542 * t1723;
t1724 = pkin(7) * t1459;
t1645 = t1505 * t1538;
t1453 = t1509 * t1543 + t1645;
t1720 = pkin(1) * t1453;
t1719 = pkin(7) * t1453;
t1535 = sin(qJ(6));
t1536 = sin(qJ(5));
t1541 = cos(qJ(5));
t1638 = qJD(2) - qJD(4);
t1444 = t1483 * t1536 + t1541 * t1638;
t1445 = t1541 * t1483 - t1536 * t1638;
t1540 = cos(qJ(6));
t1391 = t1540 * t1444 + t1445 * t1535;
t1393 = -t1444 * t1535 + t1445 * t1540;
t1333 = t1393 * t1391;
t1496 = t1520 + t1632;
t1524 = t1543 * qJDD(1);
t1622 = qJD(2) * t1674;
t1591 = -t1524 + t1622;
t1612 = t1537 * t1496 - t1542 * t1591;
t1407 = -t1483 * qJD(4) - t1612;
t1403 = qJDD(5) - t1407;
t1401 = qJDD(6) + t1403;
t1708 = -t1333 + t1401;
t1718 = t1535 * t1708;
t1400 = t1445 * t1444;
t1706 = -t1400 + t1403;
t1717 = t1536 * t1706;
t1716 = t1540 * t1708;
t1715 = t1541 * t1706;
t1408 = -t1481 * qJD(4) + t1542 * t1496 + t1537 * t1591;
t1615 = t1481 * t1638;
t1714 = t1408 + t1615;
t1497 = t1524 - 0.2e1 * t1622;
t1447 = t1497 * t1543;
t1652 = t1495 * t1538;
t1437 = -t1447 + t1652;
t1532 = t1543 ^ 2;
t1503 = (-t1531 + t1532) * t1545;
t1713 = t1437 * t1539 - t1503 * t1544;
t1712 = t1437 * t1544 + t1503 * t1539;
t1641 = t1532 * t1545;
t1511 = t1641 - t1693;
t1457 = -t1511 * t1543 + t1645;
t1711 = t1457 * t1539 + t1544 * t1524;
t1710 = t1457 * t1544 - t1539 * t1524;
t1630 = t1638 ^ 2;
t1506 = -qJD(2) * pkin(3) - pkin(8) * t1674;
t1507 = t1539 * g(1) - t1544 * g(2);
t1484 = qJDD(1) * pkin(1) + t1545 * pkin(7) + t1507;
t1554 = -pkin(2) * t1622 + t1484;
t1691 = 2 * qJD(3);
t1344 = t1524 * pkin(2) + t1496 * qJ(3) - t1591 * pkin(3) - pkin(8) * t1641 + (qJ(3) * t1673 + (-pkin(2) * qJD(2) + t1506 + t1691) * t1538) * qJD(1) + t1554;
t1614 = t1638 * t1483;
t1262 = -t1714 * pkin(9) + (-t1407 - t1614) * pkin(4) + t1344;
t1672 = qJD(3) * qJD(2);
t1525 = 0.2e1 * t1672;
t1508 = g(1) * t1544 + g(2) * t1539;
t1553 = -pkin(1) * t1545 + qJDD(1) * pkin(7) - t1508;
t1463 = -t1538 * g(3) + t1543 * t1553;
t1688 = pkin(2) * t1543;
t1594 = -qJ(3) * t1538 - t1688;
t1585 = -t1693 * pkin(2) + qJDD(2) * qJ(3) + t1594 * t1709 + t1463;
t1402 = t1525 + t1585;
t1357 = -pkin(3) * t1641 + pkin(8) * t1591 + qJD(2) * t1506 + t1402;
t1504 = qJDD(2) + t1514;
t1462 = t1543 * g(3) + t1538 * t1553;
t1549 = qJDD(2) * pkin(2) + t1693 * qJ(3) - qJDD(3) - t1462;
t1569 = t1538 * t1594;
t1546 = -t1496 * pkin(8) - t1504 * pkin(3) + (pkin(8) * t1673 + qJD(1) * t1569) * qJD(1) - t1549;
t1288 = t1542 * t1357 + t1537 * t1546;
t1431 = pkin(4) * t1481 - pkin(9) * t1483;
t1268 = -pkin(4) * t1630 - pkin(9) * t1631 - t1481 * t1431 + t1288;
t1189 = -t1541 * t1262 + t1536 * t1268;
t1190 = t1536 * t1262 + t1541 * t1268;
t1123 = t1536 * t1189 + t1541 * t1190;
t1287 = t1537 * t1357 - t1542 * t1546;
t1206 = -t1542 * t1287 + t1288 * t1537;
t1563 = -t1541 * t1408 + t1536 * t1631;
t1351 = -t1444 * qJD(5) - t1563;
t1562 = -t1536 * t1408 - t1541 * t1631;
t1552 = t1445 * qJD(5) - t1562;
t1257 = -t1391 * qJD(6) + t1540 * t1351 - t1535 * t1552;
t1477 = qJD(5) + t1481;
t1474 = qJD(6) + t1477;
t1362 = t1474 * t1391;
t1707 = -t1362 + t1257;
t1418 = t1477 * t1444;
t1327 = t1418 + t1351;
t1449 = t1511 * t1538 + t1639;
t1145 = pkin(5) * t1706 - pkin(10) * t1327 - t1189;
t1410 = pkin(5) * t1477 - pkin(10) * t1445;
t1692 = t1444 ^ 2;
t1162 = -pkin(5) * t1692 - pkin(10) * t1552 - t1477 * t1410 + t1190;
t1103 = -t1540 * t1145 + t1162 * t1535;
t1104 = t1145 * t1535 + t1162 * t1540;
t1063 = t1103 * t1535 + t1540 * t1104;
t1062 = -t1103 * t1540 + t1104 * t1535;
t1671 = t1062 * t1536;
t1042 = t1063 * t1541 - t1671;
t1267 = t1631 * pkin(4) - t1630 * pkin(9) + t1431 * t1483 + t1287;
t1193 = pkin(5) * t1552 - pkin(10) * t1692 + t1410 * t1445 + t1267;
t1037 = t1042 * t1537 - t1193 * t1542;
t1038 = t1042 * t1542 + t1193 * t1537;
t1055 = -pkin(5) * t1193 + pkin(10) * t1063;
t1607 = pkin(4) * t1193 - pkin(9) * t1042 + pkin(10) * t1671 - t1541 * t1055;
t1690 = pkin(2) + pkin(3);
t1705 = qJ(3) * t1038 - t1037 * t1690 + t1607;
t1613 = t1535 * t1351 + t1540 * t1552;
t1219 = (qJD(6) - t1474) * t1393 + t1613;
t1222 = t1362 + t1257;
t1159 = -t1219 * t1535 - t1222 * t1540;
t1161 = -t1219 * t1540 + t1222 * t1535;
t1109 = -t1159 * t1536 + t1161 * t1541;
t1389 = t1391 ^ 2;
t1390 = t1393 ^ 2;
t1278 = -t1389 - t1390;
t1096 = t1109 * t1537 - t1278 * t1542;
t1097 = t1109 * t1542 + t1278 * t1537;
t1051 = -pkin(5) * t1278 + pkin(10) * t1161 + t1063;
t1053 = -pkin(10) * t1159 - t1062;
t1608 = -pkin(4) * t1278 + pkin(9) * t1109 + t1541 * t1051 + t1536 * t1053;
t1704 = qJ(3) * t1097 - t1096 * t1690 - t1608;
t1112 = t1123 * t1537 - t1267 * t1542;
t1113 = t1123 * t1542 + t1267 * t1537;
t1636 = pkin(4) * t1267 - pkin(9) * t1123;
t1703 = qJ(3) * t1113 - t1112 * t1690 + t1636;
t1472 = t1474 ^ 2;
t1310 = -t1472 - t1389;
t1227 = t1310 * t1535 + t1716;
t1228 = t1310 * t1540 - t1718;
t1168 = -t1227 * t1536 + t1228 * t1541;
t1218 = (qJD(6) + t1474) * t1393 + t1613;
t1118 = t1168 * t1537 - t1218 * t1542;
t1119 = t1168 * t1542 + t1218 * t1537;
t1668 = t1193 * t1540;
t1115 = -pkin(5) * t1218 + pkin(10) * t1228 - t1668;
t1669 = t1193 * t1535;
t1139 = -pkin(10) * t1227 + t1669;
t1606 = -pkin(4) * t1218 + pkin(9) * t1168 + t1541 * t1115 + t1536 * t1139;
t1702 = qJ(3) * t1119 - t1118 * t1690 - t1606;
t1340 = -t1390 - t1472;
t1283 = t1333 + t1401;
t1667 = t1283 * t1535;
t1243 = t1340 * t1540 - t1667;
t1666 = t1283 * t1540;
t1244 = -t1340 * t1535 - t1666;
t1179 = -t1243 * t1536 + t1244 * t1541;
t1124 = t1179 * t1537 - t1542 * t1707;
t1125 = t1179 * t1542 + t1537 * t1707;
t1117 = -pkin(5) * t1707 + pkin(10) * t1244 + t1669;
t1147 = -pkin(10) * t1243 + t1668;
t1605 = -pkin(4) * t1707 + pkin(9) * t1179 + t1541 * t1117 + t1536 * t1147;
t1701 = qJ(3) * t1125 - t1124 * t1690 - t1605;
t1325 = (-qJD(5) + t1477) * t1445 + t1562;
t1250 = t1325 * t1541 + t1327 * t1536;
t1442 = t1445 ^ 2;
t1356 = t1442 + t1692;
t1204 = t1250 * t1537 + t1356 * t1542;
t1205 = t1250 * t1542 - t1356 * t1537;
t1604 = pkin(4) * t1356 + pkin(9) * t1250 + t1123;
t1700 = qJ(3) * t1205 - t1204 * t1690 - t1604;
t1207 = t1537 * t1287 + t1542 * t1288;
t1699 = qJ(3) * t1207 - t1206 * t1690;
t1475 = t1477 ^ 2;
t1372 = -t1475 - t1692;
t1277 = t1372 * t1541 - t1717;
t1419 = t1477 * t1445;
t1324 = -t1419 - t1552;
t1224 = t1277 * t1537 + t1324 * t1542;
t1225 = t1277 * t1542 - t1324 * t1537;
t1264 = t1541 * t1267;
t1620 = -pkin(4) * t1324 - pkin(9) * t1277 + t1264;
t1698 = qJ(3) * t1225 - t1224 * t1690 + t1620;
t1387 = -t1442 - t1475;
t1337 = t1400 + t1403;
t1663 = t1337 * t1541;
t1281 = -t1387 * t1536 - t1663;
t1328 = (qJD(5) + t1477) * t1444 + t1563;
t1229 = t1281 * t1537 + t1328 * t1542;
t1230 = t1281 * t1542 - t1328 * t1537;
t1263 = t1536 * t1267;
t1621 = pkin(4) * t1328 + pkin(9) * t1281 + t1263;
t1697 = qJ(3) * t1230 - t1229 * t1690 - t1621;
t1374 = qJD(2) * t1483 + t1612;
t1378 = -t1615 + t1408;
t1305 = -t1374 * t1537 - t1378 * t1542;
t1307 = -t1374 * t1542 + t1378 * t1537;
t1696 = qJ(3) * t1307 - t1305 * t1690;
t1479 = t1481 ^ 2;
t1421 = -t1630 - t1479;
t1360 = t1421 * t1537 + t1725;
t1361 = t1421 * t1542 - t1726;
t1695 = qJ(3) * t1361 - t1360 * t1690 + t1287;
t1480 = t1483 ^ 2;
t1461 = -t1480 - t1630;
t1428 = -t1655 + t1631;
t1662 = t1428 * t1537;
t1379 = t1461 * t1542 + t1662;
t1661 = t1428 * t1542;
t1380 = -t1461 * t1537 + t1661;
t1694 = qJ(3) * t1380 - t1379 * t1690 + t1288;
t1489 = t1543 * t1504;
t1512 = -t1641 - t1693;
t1450 = t1512 * t1538 + t1489;
t1689 = pkin(1) * t1450;
t1687 = pkin(4) * t1537;
t1686 = pkin(5) * t1062;
t1685 = pkin(5) * t1159;
t1646 = t1504 * t1538;
t1456 = t1512 * t1543 - t1646;
t1684 = pkin(6) * (t1456 * t1539 + t1497 * t1544);
t1633 = t1531 + t1532;
t1499 = t1633 * qJDD(1);
t1502 = t1633 * t1545;
t1683 = pkin(6) * (t1499 * t1539 + t1502 * t1544);
t1682 = pkin(7) * t1450;
t1681 = pkin(8) * t1113;
t1680 = pkin(8) * t1206;
t1679 = pkin(8) * t1207;
t1678 = pkin(9) * t1542;
t1675 = qJD(1) * qJD(2);
t1670 = t1062 * t1541;
t1664 = t1337 * t1536;
t1660 = t1474 * t1393;
t1659 = t1474 * t1535;
t1658 = t1474 * t1540;
t1657 = t1477 * t1536;
t1656 = t1477 * t1541;
t1654 = t1484 * t1538;
t1653 = t1484 * t1543;
t1649 = t1497 * t1538;
t1640 = t1537 * t1344;
t1339 = t1542 * t1344;
t1635 = pkin(1) * t1497 + pkin(7) * t1456;
t1634 = pkin(1) * t1502 + pkin(7) * t1499;
t1629 = pkin(4) * t1542 + pkin(3);
t1628 = t1537 * t1333;
t1627 = t1542 * t1333;
t1626 = t1537 * t1400;
t1625 = t1542 * t1400;
t1624 = t1539 * t1655;
t1623 = t1544 * t1655;
t1122 = -t1189 * t1541 + t1190 * t1536;
t1617 = -pkin(8) * t1112 + t1122 * t1687;
t1616 = -pkin(8) * t1379 + t1339;
t1399 = t1462 * t1538 + t1543 * t1463;
t1611 = -t1507 * t1539 - t1544 * t1508;
t1610 = t1539 * t1514;
t1609 = t1544 * t1514;
t1248 = t1325 * t1536 - t1327 * t1541;
t1111 = -pkin(9) * t1248 - t1122;
t1603 = -pkin(8) * t1204 + t1542 * t1111 + t1248 * t1687;
t1501 = qJDD(1) * t1544 - t1539 * t1545;
t1602 = -pkin(6) * t1501 - g(3) * t1539;
t1599 = t1537 * t1615;
t1598 = t1537 * t1614;
t1597 = t1542 * t1615;
t1596 = t1542 * t1614;
t1406 = t1545 * t1569 - t1549;
t1595 = -pkin(2) * t1406 + qJ(3) * t1402;
t1593 = pkin(2) * t1538 - qJ(3) * t1543;
t1592 = t1496 + t1520;
t1590 = -pkin(8) * t1205 - t1537 * t1111;
t1589 = -pkin(8) * t1360 + t1640;
t1588 = -pkin(8) * t1361 + t1339;
t1587 = -pkin(8) * t1380 - t1640;
t1586 = -pkin(9) * t1537 - t1629;
t1398 = t1462 * t1543 - t1463 * t1538;
t1434 = t1495 * t1543 + t1649;
t1584 = t1507 * t1544 - t1508 * t1539;
t1580 = t1591 * pkin(2);
t1579 = pkin(5) * t1227 - t1103;
t1041 = t1063 * t1536 + t1670;
t1026 = -pkin(9) * t1041 - pkin(10) * t1670 - t1055 * t1536;
t1029 = -pkin(4) * t1041 - t1686;
t1578 = -pkin(8) * t1037 + t1542 * t1026 - t1029 * t1537;
t1107 = t1159 * t1541 + t1161 * t1536;
t1028 = -pkin(9) * t1107 - t1051 * t1536 + t1053 * t1541;
t1076 = -pkin(4) * t1107 - t1685;
t1577 = -pkin(8) * t1096 + t1542 * t1028 - t1076 * t1537;
t1167 = t1227 * t1541 + t1228 * t1536;
t1065 = -pkin(9) * t1167 - t1115 * t1536 + t1139 * t1541;
t1069 = -pkin(4) * t1167 - t1579;
t1576 = -pkin(8) * t1118 + t1542 * t1065 - t1069 * t1537;
t1178 = t1243 * t1541 + t1244 * t1536;
t1071 = -pkin(9) * t1178 - t1117 * t1536 + t1147 * t1541;
t1561 = pkin(5) * t1243 - t1104;
t1072 = -pkin(4) * t1178 - t1561;
t1575 = -pkin(8) * t1124 + t1542 * t1071 - t1072 * t1537;
t1276 = t1372 * t1536 + t1715;
t1157 = -pkin(4) * t1276 + t1189;
t1195 = -pkin(9) * t1276 + t1263;
t1574 = -pkin(8) * t1224 - t1157 * t1537 + t1542 * t1195;
t1280 = t1387 * t1541 - t1664;
t1166 = -pkin(4) * t1280 + t1190;
t1197 = -pkin(9) * t1280 + t1264;
t1573 = -pkin(8) * t1229 - t1166 * t1537 + t1542 * t1197;
t1572 = -pkin(8) * t1305 - t1206;
t1571 = -pkin(8) * t1307 - t1207;
t1560 = -pkin(8) * t1038 - t1537 * t1026 - t1542 * t1029;
t1559 = -pkin(8) * t1097 - t1537 * t1028 - t1542 * t1076;
t1558 = -pkin(8) * t1119 - t1537 * t1065 - t1542 * t1069;
t1557 = -pkin(8) * t1125 - t1537 * t1071 - t1542 * t1072;
t1556 = -pkin(8) * t1225 - t1542 * t1157 - t1537 * t1195;
t1555 = -pkin(8) * t1230 - t1542 * t1166 - t1537 * t1197;
t1551 = pkin(2) * t1509 + qJ(3) * t1505 + t1585;
t1550 = t1674 * t1691 + t1554;
t1548 = qJ(3) * t1592 + t1550;
t1547 = pkin(2) * t1504 + qJ(3) * t1512 - t1406;
t1510 = t1642 - t1693;
t1500 = qJDD(1) * t1539 + t1544 * t1545;
t1492 = t1593 * qJDD(1);
t1488 = t1633 * t1675;
t1478 = -pkin(6) * t1500 + g(3) * t1544;
t1469 = -t1480 + t1630;
t1468 = t1479 - t1630;
t1467 = qJDD(2) * t1539 + t1488 * t1544;
t1466 = t1496 * t1543 - t1531 * t1675;
t1465 = -qJDD(2) * t1544 + t1488 * t1539;
t1464 = -t1532 * t1675 + t1538 * t1591;
t1458 = t1510 * t1538 + t1489;
t1452 = -t1510 * t1543 + t1646;
t1448 = t1592 * t1538;
t1439 = pkin(6) * (t1499 * t1544 - t1502 * t1539);
t1432 = t1480 - t1479;
t1427 = t1466 * t1544 - t1610;
t1426 = t1464 * t1544 + t1610;
t1425 = t1466 * t1539 + t1609;
t1424 = t1464 * t1539 - t1609;
t1423 = t1458 * t1544 + t1539 * t1632;
t1422 = t1458 * t1539 - t1544 * t1632;
t1417 = pkin(6) * (t1456 * t1544 - t1497 * t1539);
t1416 = -t1442 + t1475;
t1415 = -t1475 + t1692;
t1414 = t1597 - t1598;
t1413 = t1599 + t1596;
t1412 = -t1653 + t1719;
t1411 = -t1654 - t1682;
t1409 = -t1479 - t1480;
t1405 = t1463 + t1720;
t1404 = t1462 - t1689;
t1396 = t1442 - t1692;
t1395 = t1635 + t1653;
t1394 = -pkin(1) * t1495 - t1654 - t1724;
t1388 = qJ(3) * t1502 + t1406;
t1386 = pkin(2) * t1502 + t1402;
t1385 = -t1580 + t1548;
t1384 = t1468 * t1542 + t1662;
t1383 = -t1469 * t1537 + t1725;
t1382 = t1468 * t1537 - t1661;
t1381 = t1469 * t1542 + t1726;
t1373 = (0.2e1 * qJD(4) - qJD(2)) * t1483 + t1612;
t1371 = pkin(1) * t1484 + pkin(7) * t1399;
t1370 = t1542 * t1408 + t1598;
t1369 = t1537 * t1408 - t1596;
t1368 = -t1537 * t1407 - t1597;
t1367 = -t1542 * t1407 + t1599;
t1366 = (t1497 - t1591) * pkin(2) + t1548;
t1365 = -t1580 + (t1495 + t1592) * qJ(3) + t1550;
t1364 = t1399 + t1634;
t1359 = -t1390 + t1472;
t1358 = t1389 - t1472;
t1347 = (-t1444 * t1541 + t1445 * t1536) * t1477;
t1346 = (-t1444 * t1536 - t1445 * t1541) * t1477;
t1345 = -t1547 - t1689;
t1343 = -t1551 - 0.2e1 * t1672 - t1720;
t1342 = t1413 * t1538 + t1414 * t1543;
t1341 = -t1413 * t1543 + t1414 * t1538;
t1335 = t1402 * t1543 + t1406 * t1538;
t1334 = t1402 * t1538 - t1406 * t1543;
t1332 = t1390 - t1389;
t1331 = -pkin(2) * t1652 + t1365 * t1543 - t1719;
t1330 = qJ(3) * t1447 - t1366 * t1538 - t1682;
t1329 = -t1386 * t1538 + t1388 * t1543;
t1326 = -t1418 + t1351;
t1323 = -t1419 + t1552;
t1320 = t1382 * t1538 + t1384 * t1543;
t1319 = t1381 * t1538 + t1383 * t1543;
t1318 = -t1382 * t1543 + t1384 * t1538;
t1317 = -t1381 * t1543 + t1383 * t1538;
t1316 = t1724 + t1538 * t1365 + (pkin(1) + t1688) * t1495;
t1315 = qJ(3) * t1649 + t1366 * t1543 + t1635;
t1314 = t1351 * t1541 - t1445 * t1657;
t1313 = t1351 * t1536 + t1445 * t1656;
t1312 = t1444 * t1656 + t1536 * t1552;
t1311 = -t1444 * t1657 + t1541 * t1552;
t1309 = t1379 * t1538 + t1380 * t1543;
t1308 = -t1379 * t1543 + t1380 * t1538;
t1306 = -t1373 * t1542 - t1537 * t1714;
t1304 = -t1373 * t1537 + t1542 * t1714;
t1303 = t1347 * t1542 + t1403 * t1537;
t1302 = t1347 * t1537 - t1403 * t1542;
t1301 = t1386 * t1543 + t1388 * t1538 + t1634;
t1300 = t1415 * t1541 - t1664;
t1299 = -t1416 * t1536 + t1715;
t1298 = t1415 * t1536 + t1663;
t1297 = t1416 * t1541 + t1717;
t1296 = t1369 * t1538 + t1370 * t1543;
t1295 = -t1367 * t1538 + t1368 * t1543;
t1294 = -t1369 * t1543 + t1370 * t1538;
t1293 = t1367 * t1543 + t1368 * t1538;
t1292 = (-t1391 * t1540 + t1393 * t1535) * t1474;
t1291 = (-t1391 * t1535 - t1393 * t1540) * t1474;
t1290 = t1360 * t1538 + t1361 * t1543;
t1289 = -t1360 * t1543 + t1361 * t1538;
t1273 = t1314 * t1542 + t1626;
t1272 = t1312 * t1542 - t1626;
t1271 = t1314 * t1537 - t1625;
t1270 = t1312 * t1537 + t1625;
t1269 = -pkin(1) * t1334 - t1595;
t1261 = qJ(3) * t1714 + t1616;
t1258 = -pkin(7) * t1334 - t1385 * t1593;
t1256 = -qJD(6) * t1393 - t1613;
t1255 = qJ(3) * t1373 + t1589;
t1254 = t1358 * t1540 - t1667;
t1253 = -t1359 * t1535 + t1716;
t1252 = t1358 * t1535 + t1666;
t1251 = t1359 * t1540 + t1718;
t1249 = t1324 * t1541 - t1326 * t1536;
t1247 = t1324 * t1536 + t1326 * t1541;
t1242 = t1305 * t1538 + t1307 * t1543;
t1241 = t1304 * t1538 + t1306 * t1543;
t1240 = -t1305 * t1543 + t1307 * t1538;
t1239 = -t1304 * t1543 + t1306 * t1538;
t1238 = t1300 * t1542 - t1323 * t1537;
t1237 = t1299 * t1542 + t1327 * t1537;
t1236 = t1300 * t1537 + t1323 * t1542;
t1235 = t1299 * t1537 - t1327 * t1542;
t1234 = pkin(7) * t1335 + (pkin(1) - t1594) * t1385;
t1233 = t1690 * t1714 + t1587;
t1232 = t1302 * t1538 + t1303 * t1543;
t1231 = -t1302 * t1543 + t1303 * t1538;
t1226 = t1373 * t1690 + t1588;
t1215 = t1257 * t1540 - t1393 * t1659;
t1214 = t1257 * t1535 + t1393 * t1658;
t1213 = -t1256 * t1535 + t1391 * t1658;
t1212 = t1256 * t1540 + t1391 * t1659;
t1211 = t1249 * t1542 + t1396 * t1537;
t1210 = t1249 * t1537 - t1396 * t1542;
t1209 = -t1291 * t1536 + t1292 * t1541;
t1208 = t1291 * t1541 + t1292 * t1536;
t1203 = t1209 * t1542 + t1401 * t1537;
t1202 = t1209 * t1537 - t1401 * t1542;
t1201 = t1271 * t1538 + t1273 * t1543;
t1200 = t1270 * t1538 + t1272 * t1543;
t1199 = -t1271 * t1543 + t1273 * t1538;
t1198 = -t1270 * t1543 + t1272 * t1538;
t1191 = qJ(3) * t1344 - t1680;
t1186 = -pkin(1) * t1308 - t1694;
t1185 = -t1252 * t1536 + t1254 * t1541;
t1184 = -t1251 * t1536 + t1253 * t1541;
t1183 = t1252 * t1541 + t1254 * t1536;
t1182 = t1251 * t1541 + t1253 * t1536;
t1181 = qJ(3) * t1409 + t1572;
t1180 = t1344 * t1690 - t1679;
t1176 = -pkin(1) * t1289 - t1695;
t1175 = t1409 * t1690 + t1571;
t1174 = t1236 * t1538 + t1238 * t1543;
t1173 = t1235 * t1538 + t1237 * t1543;
t1172 = -t1236 * t1543 + t1238 * t1538;
t1171 = -t1235 * t1543 + t1237 * t1538;
t1170 = t1229 * t1538 + t1230 * t1543;
t1169 = -t1229 * t1543 + t1230 * t1538;
t1164 = t1224 * t1538 + t1225 * t1543;
t1163 = -t1224 * t1543 + t1225 * t1538;
t1160 = -t1218 * t1540 - t1535 * t1707;
t1158 = -t1218 * t1535 + t1540 * t1707;
t1156 = -t1214 * t1536 + t1215 * t1541;
t1155 = -t1212 * t1536 + t1213 * t1541;
t1154 = t1214 * t1541 + t1215 * t1536;
t1153 = t1212 * t1541 + t1213 * t1536;
t1152 = t1210 * t1538 + t1211 * t1543;
t1151 = -t1210 * t1543 + t1211 * t1538;
t1150 = t1206 * t1538 + t1207 * t1543;
t1149 = -t1206 * t1543 + t1207 * t1538;
t1148 = -pkin(7) * t1308 - t1233 * t1538 + t1261 * t1543;
t1143 = t1204 * t1538 + t1205 * t1543;
t1142 = -t1204 * t1543 + t1205 * t1538;
t1141 = -pkin(1) * t1240 - t1696;
t1140 = -pkin(7) * t1289 - t1226 * t1538 + t1255 * t1543;
t1137 = pkin(1) * t1714 + pkin(7) * t1309 + t1233 * t1543 + t1261 * t1538;
t1136 = t1202 * t1538 + t1203 * t1543;
t1135 = -t1202 * t1543 + t1203 * t1538;
t1134 = t1156 * t1542 + t1628;
t1133 = t1155 * t1542 - t1628;
t1132 = t1156 * t1537 - t1627;
t1131 = t1155 * t1537 + t1627;
t1130 = pkin(1) * t1373 + pkin(7) * t1290 + t1226 * t1543 + t1255 * t1538;
t1129 = t1185 * t1542 - t1219 * t1537;
t1128 = t1184 * t1542 + t1222 * t1537;
t1127 = t1185 * t1537 + t1219 * t1542;
t1126 = t1184 * t1537 - t1222 * t1542;
t1108 = -t1158 * t1536 + t1160 * t1541;
t1106 = t1158 * t1541 + t1160 * t1536;
t1101 = -pkin(7) * t1240 - t1175 * t1538 + t1181 * t1543;
t1100 = t1108 * t1542 + t1332 * t1537;
t1099 = t1108 * t1537 - t1332 * t1542;
t1098 = pkin(1) * t1409 + pkin(7) * t1242 + t1175 * t1543 + t1181 * t1538;
t1095 = qJ(3) * t1280 + t1573;
t1094 = qJ(3) * t1276 + t1574;
t1093 = t1132 * t1538 + t1134 * t1543;
t1092 = t1131 * t1538 + t1133 * t1543;
t1091 = -t1132 * t1543 + t1134 * t1538;
t1090 = -t1131 * t1543 + t1133 * t1538;
t1089 = -pkin(7) * t1149 - t1180 * t1538 + t1191 * t1543;
t1088 = t1127 * t1538 + t1129 * t1543;
t1087 = t1126 * t1538 + t1128 * t1543;
t1086 = -t1127 * t1543 + t1129 * t1538;
t1085 = -t1126 * t1543 + t1128 * t1538;
t1084 = t1280 * t1690 + t1555;
t1083 = -pkin(1) * t1149 - t1699;
t1082 = t1276 * t1690 + t1556;
t1081 = pkin(1) * t1344 + pkin(7) * t1150 + t1180 * t1543 + t1191 * t1538;
t1080 = t1124 * t1538 + t1125 * t1543;
t1079 = -t1124 * t1543 + t1125 * t1538;
t1078 = t1118 * t1538 + t1119 * t1543;
t1077 = -t1118 * t1543 + t1119 * t1538;
t1075 = -pkin(1) * t1169 - t1697;
t1074 = qJ(3) * t1248 + t1603;
t1073 = -pkin(1) * t1163 - t1698;
t1068 = t1112 * t1538 + t1113 * t1543;
t1067 = -t1112 * t1543 + t1113 * t1538;
t1066 = (pkin(2) + t1629) * t1248 + t1590;
t1060 = t1099 * t1538 + t1100 * t1543;
t1059 = -t1099 * t1543 + t1100 * t1538;
t1058 = t1096 * t1538 + t1097 * t1543;
t1057 = -t1096 * t1543 + t1097 * t1538;
t1056 = -pkin(1) * t1142 - t1700;
t1049 = -pkin(7) * t1169 - t1084 * t1538 + t1095 * t1543;
t1048 = -pkin(7) * t1163 - t1082 * t1538 + t1094 * t1543;
t1047 = (qJ(3) - t1678) * t1122 + t1617;
t1046 = pkin(1) * t1280 + pkin(7) * t1170 + t1084 * t1543 + t1095 * t1538;
t1045 = pkin(1) * t1276 + pkin(7) * t1164 + t1082 * t1543 + t1094 * t1538;
t1044 = -t1681 + (pkin(2) - t1586) * t1122;
t1043 = -pkin(7) * t1142 - t1066 * t1538 + t1074 * t1543;
t1039 = pkin(1) * t1248 + pkin(7) * t1143 + t1066 * t1543 + t1074 * t1538;
t1036 = qJ(3) * t1178 + t1575;
t1035 = qJ(3) * t1167 + t1576;
t1034 = t1178 * t1690 + t1557;
t1033 = t1167 * t1690 + t1558;
t1032 = -pkin(1) * t1079 - t1701;
t1031 = -pkin(1) * t1067 - t1703;
t1030 = -pkin(1) * t1077 - t1702;
t1024 = -pkin(7) * t1067 - t1044 * t1538 + t1047 * t1543;
t1023 = t1037 * t1538 + t1038 * t1543;
t1022 = -t1037 * t1543 + t1038 * t1538;
t1021 = pkin(1) * t1122 + pkin(7) * t1068 + t1044 * t1543 + t1047 * t1538;
t1020 = qJ(3) * t1107 + t1577;
t1019 = -pkin(7) * t1079 - t1034 * t1538 + t1036 * t1543;
t1018 = t1107 * t1690 + t1559;
t1017 = pkin(1) * t1178 + pkin(7) * t1080 + t1034 * t1543 + t1036 * t1538;
t1016 = -pkin(7) * t1077 - t1033 * t1538 + t1035 * t1543;
t1015 = -pkin(1) * t1057 - t1704;
t1014 = pkin(1) * t1167 + pkin(7) * t1078 + t1033 * t1543 + t1035 * t1538;
t1013 = -pkin(7) * t1057 - t1018 * t1538 + t1020 * t1543;
t1012 = qJ(3) * t1041 + t1578;
t1011 = pkin(1) * t1107 + pkin(7) * t1058 + t1018 * t1543 + t1020 * t1538;
t1010 = t1041 * t1690 + t1560;
t1009 = -pkin(1) * t1022 - t1705;
t1008 = -pkin(7) * t1022 - t1010 * t1538 + t1012 * t1543;
t1007 = pkin(1) * t1041 + pkin(7) * t1023 + t1010 * t1543 + t1012 * t1538;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1501, 0, -t1500, 0, t1602, -t1478, -t1584, -pkin(6) * t1584, t1427, -t1712, t1423, t1426, -t1710, t1467, -t1404 * t1539 + t1411 * t1544 - t1684, -t1539 * t1405 + t1544 * t1412 + t1727, t1398 * t1544 - t1683, -pkin(6) * (t1399 * t1539 + t1484 * t1544) - (pkin(1) * t1539 - pkin(7) * t1544) * t1398, t1427, t1423, t1712, t1467, t1710, t1426, t1330 * t1544 - t1345 * t1539 - t1684, t1329 * t1544 - t1492 * t1539 - t1683, t1544 * t1331 - t1539 * t1343 - t1727, t1544 * t1258 - t1539 * t1269 - pkin(6) * (t1335 * t1539 + t1385 * t1544), t1296 * t1544 - t1624, t1241 * t1544 - t1432 * t1539, t1319 * t1544 - t1378 * t1539, t1295 * t1544 + t1624, t1320 * t1544 + t1374 * t1539, t1544 * t1342 + t1539 * t1631, t1544 * t1140 - t1539 * t1176 - pkin(6) * (t1290 * t1539 + t1373 * t1544), t1544 * t1148 - t1539 * t1186 - pkin(6) * (t1309 * t1539 + t1544 * t1714), t1544 * t1101 - t1539 * t1141 - pkin(6) * (t1242 * t1539 + t1409 * t1544), t1544 * t1089 - t1539 * t1083 - pkin(6) * (t1150 * t1539 + t1344 * t1544), t1201 * t1544 - t1313 * t1539, t1152 * t1544 - t1247 * t1539, t1173 * t1544 - t1297 * t1539, t1200 * t1544 + t1311 * t1539, t1174 * t1544 - t1298 * t1539, t1232 * t1544 - t1346 * t1539, t1544 * t1048 - t1539 * t1073 - pkin(6) * (t1164 * t1539 + t1276 * t1544), t1544 * t1049 - t1539 * t1075 - pkin(6) * (t1170 * t1539 + t1280 * t1544), t1544 * t1043 - t1539 * t1056 - pkin(6) * (t1143 * t1539 + t1248 * t1544), t1544 * t1024 - t1539 * t1031 - pkin(6) * (t1068 * t1539 + t1122 * t1544), t1093 * t1544 - t1154 * t1539, t1060 * t1544 - t1106 * t1539, t1087 * t1544 - t1182 * t1539, t1092 * t1544 - t1153 * t1539, t1088 * t1544 - t1183 * t1539, t1136 * t1544 - t1208 * t1539, t1544 * t1016 - t1539 * t1030 - pkin(6) * (t1078 * t1539 + t1167 * t1544), t1544 * t1019 - t1539 * t1032 - pkin(6) * (t1080 * t1539 + t1178 * t1544), t1544 * t1013 - t1539 * t1015 - pkin(6) * (t1058 * t1539 + t1107 * t1544), t1544 * t1008 - t1539 * t1009 - pkin(6) * (t1023 * t1539 + t1041 * t1544); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1500, 0, t1501, 0, t1478, t1602, t1611, pkin(6) * t1611, t1425, -t1713, t1422, t1424, -t1711, t1465, t1404 * t1544 + t1411 * t1539 + t1417, t1544 * t1405 + t1539 * t1412 - t1728, t1398 * t1539 + t1439, pkin(6) * (t1399 * t1544 - t1484 * t1539) - (-pkin(1) * t1544 - pkin(7) * t1539) * t1398, t1425, t1422, t1713, t1465, t1711, t1424, t1330 * t1539 + t1345 * t1544 + t1417, t1329 * t1539 + t1492 * t1544 + t1439, t1539 * t1331 + t1544 * t1343 + t1728, t1539 * t1258 + t1544 * t1269 + pkin(6) * (t1335 * t1544 - t1385 * t1539), t1296 * t1539 + t1623, t1241 * t1539 + t1432 * t1544, t1319 * t1539 + t1378 * t1544, t1295 * t1539 - t1623, t1320 * t1539 - t1374 * t1544, t1539 * t1342 - t1544 * t1631, t1539 * t1140 + t1544 * t1176 + pkin(6) * (t1290 * t1544 - t1373 * t1539), t1539 * t1148 + t1544 * t1186 + pkin(6) * (t1309 * t1544 - t1539 * t1714), t1539 * t1101 + t1544 * t1141 + pkin(6) * (t1242 * t1544 - t1409 * t1539), t1539 * t1089 + t1544 * t1083 + pkin(6) * (t1150 * t1544 - t1344 * t1539), t1201 * t1539 + t1313 * t1544, t1152 * t1539 + t1247 * t1544, t1173 * t1539 + t1297 * t1544, t1200 * t1539 - t1311 * t1544, t1174 * t1539 + t1298 * t1544, t1232 * t1539 + t1346 * t1544, t1539 * t1048 + t1544 * t1073 + pkin(6) * (t1164 * t1544 - t1276 * t1539), t1539 * t1049 + t1544 * t1075 + pkin(6) * (t1170 * t1544 - t1280 * t1539), t1539 * t1043 + t1544 * t1056 + pkin(6) * (t1143 * t1544 - t1248 * t1539), t1539 * t1024 + t1544 * t1031 + pkin(6) * (t1068 * t1544 - t1122 * t1539), t1093 * t1539 + t1154 * t1544, t1060 * t1539 + t1106 * t1544, t1087 * t1539 + t1182 * t1544, t1092 * t1539 + t1153 * t1544, t1088 * t1539 + t1183 * t1544, t1136 * t1539 + t1208 * t1544, t1539 * t1016 + t1544 * t1030 + pkin(6) * (t1078 * t1544 - t1167 * t1539), t1539 * t1019 + t1544 * t1032 + pkin(6) * (t1080 * t1544 - t1178 * t1539), t1539 * t1013 + t1544 * t1015 + pkin(6) * (t1058 * t1544 - t1107 * t1539), t1539 * t1008 + t1544 * t1009 + pkin(6) * (t1023 * t1544 - t1041 * t1539); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1507, t1508, 0, 0, t1448, t1434, t1452, t1447, t1449, 0, t1395, t1394, t1364, t1371, t1448, t1452, -t1434, 0, -t1449, t1447, t1315, t1301, t1316, t1234, t1294, t1239, t1317, t1293, t1318, t1341, t1130, t1137, t1098, t1081, t1199, t1151, t1171, t1198, t1172, t1231, t1045, t1046, t1039, t1021, t1091, t1059, t1085, t1090, t1086, t1135, t1014, t1017, t1011, t1007; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1545, 0, 0, -g(3), -t1507, 0, t1466, -t1437, t1458, t1464, -t1457, t1488, t1411, t1412, t1398, pkin(7) * t1398, t1466, t1458, t1437, t1488, t1457, t1464, t1330, t1329, t1331, t1258, t1296, t1241, t1319, t1295, t1320, t1342, t1140, t1148, t1101, t1089, t1201, t1152, t1173, t1200, t1174, t1232, t1048, t1049, t1043, t1024, t1093, t1060, t1087, t1092, t1088, t1136, t1016, t1019, t1013, t1008; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1545, 0, qJDD(1), 0, g(3), 0, -t1508, 0, t1514, t1503, -t1632, -t1514, -t1524, -qJDD(2), t1404, t1405, 0, pkin(1) * t1398, t1514, -t1632, -t1503, -qJDD(2), t1524, -t1514, t1345, t1492, t1343, t1269, t1655, t1432, t1378, -t1655, -t1374, -t1631, t1176, t1186, t1141, t1083, t1313, t1247, t1297, -t1311, t1298, t1346, t1073, t1075, t1056, t1031, t1154, t1106, t1182, t1153, t1183, t1208, t1030, t1032, t1015, t1009; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1507, t1508, 0, 0, t1448, t1434, t1452, t1447, t1449, 0, t1395, t1394, t1364, t1371, t1448, t1452, -t1434, 0, -t1449, t1447, t1315, t1301, t1316, t1234, t1294, t1239, t1317, t1293, t1318, t1341, t1130, t1137, t1098, t1081, t1199, t1151, t1171, t1198, t1172, t1231, t1045, t1046, t1039, t1021, t1091, t1059, t1085, t1090, t1086, t1135, t1014, t1017, t1011, t1007; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1496, t1497, t1504, -t1520, t1511, t1520, 0, -t1484, t1462, 0, t1496, t1504, -t1497, t1520, -t1511, -t1520, qJ(3) * t1497, t1388, t1365, qJ(3) * t1385, t1370, t1306, t1383, t1368, t1384, t1414, t1255, t1261, t1181, t1191, t1273, t1211, t1237, t1272, t1238, t1303, t1094, t1095, t1074, t1047, t1134, t1100, t1128, t1133, t1129, t1203, t1035, t1036, t1020, t1012; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1622, t1495, -t1510, -t1591, t1505, -t1622, t1484, 0, t1463, 0, t1622, -t1510, -t1495, -t1622, -t1505, -t1591, t1366, t1386, pkin(2) * t1495, pkin(2) * t1385, -t1369, -t1304, -t1381, t1367, -t1382, -t1413, t1226, t1233, t1175, t1180, -t1271, -t1210, -t1235, -t1270, -t1236, -t1302, t1082, t1084, t1066, t1044, -t1132, -t1099, -t1126, -t1131, -t1127, -t1202, t1033, t1034, t1018, t1010; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1514, -t1503, t1632, t1514, t1524, qJDD(2), -t1462, -t1463, 0, 0, -t1514, t1632, t1503, qJDD(2), -t1524, t1514, t1547, -t1492, t1525 + t1551, t1595, -t1655, -t1432, -t1378, t1655, t1374, t1631, t1695, t1694, t1696, t1699, -t1313, -t1247, -t1297, t1311, -t1298, -t1346, t1698, t1697, t1700, t1703, -t1154, -t1106, -t1182, -t1153, -t1183, -t1208, t1702, t1701, t1704, t1705; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1496, t1504, -t1497, t1520, -t1511, -t1520, 0, t1406, t1385, 0, t1370, t1306, t1383, t1368, t1384, t1414, t1589, t1616, t1572, -t1680, t1273, t1211, t1237, t1272, t1238, t1303, t1574, t1573, t1603, -t1122 * t1678 + t1617, t1134, t1100, t1128, t1133, t1129, t1203, t1576, t1575, t1577, t1578; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1514, t1632, t1503, qJDD(2), -t1524, t1514, -t1406, 0, t1402, 0, -t1655, -t1432, -t1378, t1655, t1374, t1631, -pkin(3) * t1360 + t1287, -pkin(3) * t1379 + t1288, -pkin(3) * t1305, -pkin(3) * t1206, -t1313, -t1247, -t1297, t1311, -t1298, -t1346, -pkin(3) * t1224 + t1620, -pkin(3) * t1229 - t1621, -pkin(3) * t1204 - t1604, -pkin(3) * t1112 + t1636, -t1154, -t1106, -t1182, -t1153, -t1183, -t1208, -pkin(3) * t1118 - t1606, -pkin(3) * t1124 - t1605, -pkin(3) * t1096 - t1608, -pkin(3) * t1037 + t1607; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1622, t1510, t1495, t1622, t1505, t1591, -t1385, -t1402, 0, 0, t1369, t1304, t1381, -t1367, t1382, t1413, -pkin(3) * t1373 - t1588, -pkin(3) * t1714 - t1587, -pkin(3) * t1409 - t1571, -pkin(3) * t1344 + t1679, t1271, t1210, t1235, t1270, t1236, t1302, -pkin(3) * t1276 - t1556, -pkin(3) * t1280 - t1555, -t1248 * t1629 - t1590, t1122 * t1586 + t1681, t1132, t1099, t1126, t1131, t1127, t1202, -pkin(3) * t1167 - t1558, -pkin(3) * t1178 - t1557, -pkin(3) * t1107 - t1559, -pkin(3) * t1041 - t1560; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1408, -t1373, t1723, -t1615, t1468, t1615, 0, t1344, t1287, 0, t1314, t1249, t1299, t1312, t1300, t1347, t1195, t1197, t1111, -pkin(9) * t1122, t1156, t1108, t1184, t1155, t1185, t1209, t1065, t1071, t1028, t1026; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1614, t1714, t1469, t1407, -t1428, t1614, -t1344, 0, t1288, 0, -t1400, -t1396, -t1327, t1400, t1323, -t1403, t1157, t1166, -pkin(4) * t1248, -pkin(4) * t1122, -t1333, -t1332, -t1222, t1333, t1219, -t1401, t1069, t1072, t1076, t1029; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1655, t1432, t1378, -t1655, -t1374, -t1631, -t1287, -t1288, 0, 0, t1313, t1247, t1297, -t1311, t1298, t1346, -t1620, t1621, t1604, -t1636, t1154, t1106, t1182, t1153, t1183, t1208, t1606, t1605, t1608, -t1607; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1351, t1324, t1706, t1418, t1415, -t1418, 0, t1267, t1189, 0, t1215, t1160, t1253, t1213, t1254, t1292, t1139, t1147, t1053, -pkin(10) * t1062; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1419, t1326, t1416, -t1552, t1337, -t1419, -t1267, 0, t1190, 0, t1214, t1158, t1251, t1212, t1252, t1291, t1115, t1117, t1051, t1055; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1400, t1396, t1327, -t1400, -t1323, t1403, -t1189, -t1190, 0, 0, t1333, t1332, t1222, -t1333, -t1219, t1401, t1579, t1561, t1685, t1686; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1257, -t1218, t1708, t1362, t1358, -t1362, 0, t1193, t1103, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1660, t1707, t1359, t1256, t1283, -t1660, -t1193, 0, t1104, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1333, t1332, t1222, -t1333, -t1219, t1401, -t1103, -t1104, 0, 0;];
m_new_reg  = t1;