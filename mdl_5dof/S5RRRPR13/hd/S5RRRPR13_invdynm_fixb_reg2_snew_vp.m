% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S5RRRPR13
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
%
% Output:
% m_new_reg [(3*6)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S5RRRPR13_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_invdynm_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_invdynm_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR13_invdynm_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR13_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:47:31
% EndTime: 2019-12-31 21:48:08
% DurationCPUTime: 38.75s
% Computational Cost: add. (119385->789), mult. (261193->1088), div. (0->0), fcn. (200513->10), ass. (0->546)
t1522 = sin(pkin(5));
t1530 = cos(qJ(2));
t1660 = t1522 * t1530;
t1509 = qJD(1) * t1660 - qJD(3);
t1703 = t1509 ^ 2;
t1525 = sin(qJ(3));
t1529 = cos(qJ(3));
t1523 = cos(pkin(5));
t1688 = qJD(1) * t1523;
t1634 = qJD(2) + t1688;
t1526 = sin(qJ(2));
t1661 = t1522 * t1526;
t1645 = qJD(1) * t1661;
t1481 = t1525 * t1634 + t1529 * t1645;
t1704 = t1481 ^ 2;
t1440 = t1704 + t1703;
t1648 = qJDD(1) * t1522;
t1490 = -qJD(2) * t1645 + t1530 * t1648;
t1484 = -qJDD(3) + t1490;
t1479 = t1525 * t1645 - t1529 * t1634;
t1671 = t1481 * t1479;
t1550 = t1484 - t1671;
t1725 = t1550 * t1525;
t1352 = t1440 * t1529 - t1725;
t1724 = t1550 * t1529;
t1353 = t1440 * t1525 + t1724;
t1647 = t1526 * qJDD(1);
t1687 = qJD(1) * t1530;
t1489 = (qJD(2) * t1687 + t1647) * t1522;
t1628 = qJDD(1) * t1523 + qJDD(2);
t1551 = t1529 * t1489 + t1525 * t1628;
t1394 = (qJD(3) - t1509) * t1479 - t1551;
t1595 = t1353 * t1526 + t1394 * t1530;
t1240 = t1522 * t1352 + t1523 * t1595;
t1292 = t1353 * t1530 - t1394 * t1526;
t1527 = sin(qJ(1));
t1531 = cos(qJ(1));
t1807 = pkin(6) * (t1240 * t1531 + t1292 * t1527);
t1806 = pkin(6) * (t1240 * t1527 - t1292 * t1531);
t1476 = t1479 ^ 2;
t1433 = -t1703 - t1476;
t1549 = t1484 + t1671;
t1726 = t1529 * t1549;
t1343 = t1433 * t1525 - t1726;
t1727 = t1525 * t1549;
t1345 = t1433 * t1529 + t1727;
t1630 = t1489 * t1525 - t1529 * t1628;
t1569 = qJD(3) * t1481 + t1630;
t1665 = t1509 * t1481;
t1388 = t1665 - t1569;
t1597 = t1345 * t1526 + t1388 * t1530;
t1235 = -t1522 * t1343 + t1523 * t1597;
t1289 = t1345 * t1530 - t1388 * t1526;
t1805 = pkin(6) * (t1235 * t1531 + t1289 * t1527);
t1804 = pkin(6) * (t1235 * t1527 - t1289 * t1531);
t1238 = -t1523 * t1352 + t1522 * t1595;
t1803 = pkin(7) * (t1238 * t1522 + t1240 * t1523);
t1233 = t1523 * t1343 + t1522 * t1597;
t1802 = pkin(7) * (t1233 * t1522 + t1235 * t1523);
t1801 = pkin(1) * t1233;
t1800 = pkin(1) * t1235;
t1799 = pkin(1) * t1238;
t1798 = pkin(1) * t1240;
t1454 = t1703 - t1476;
t1358 = t1454 * t1525 + t1724;
t1362 = t1454 * t1529 - t1725;
t1721 = t1665 + t1569;
t1590 = t1362 * t1526 - t1530 * t1721;
t1255 = -t1522 * t1358 + t1523 * t1590;
t1305 = t1362 * t1530 + t1526 * t1721;
t1793 = t1255 * t1527 - t1305 * t1531;
t1792 = t1255 * t1531 + t1305 * t1527;
t1789 = pkin(7) * t1289;
t1788 = pkin(7) * t1292;
t1251 = t1523 * t1358 + t1522 * t1590;
t1781 = pkin(2) * t1343;
t1780 = pkin(2) * t1352;
t1779 = pkin(8) * t1343;
t1778 = pkin(8) * t1345;
t1777 = pkin(8) * t1352;
t1776 = pkin(8) * t1353;
t1720 = -t1704 + t1703;
t1359 = t1525 * t1720 + t1726;
t1775 = t1359 * t1526;
t1774 = t1359 * t1530;
t1356 = t1529 * t1720 - t1727;
t1771 = t1522 * t1356;
t1767 = t1523 * t1356;
t1702 = -2 * qJD(4);
t1410 = -t1704 - t1476;
t1765 = pkin(2) * t1410;
t1764 = t1388 * t1525;
t1682 = t1388 * t1529;
t1759 = t1410 * t1526;
t1758 = t1410 * t1530;
t1428 = -t1479 * qJD(3) + t1551;
t1462 = t1479 * t1509;
t1392 = t1462 + t1428;
t1699 = pkin(2) * t1530;
t1626 = -pkin(8) * t1526 - t1699;
t1689 = qJD(1) * t1522;
t1488 = t1626 * t1689;
t1627 = t1634 ^ 2;
t1512 = g(1) * t1531 + g(2) * t1527;
t1532 = qJD(1) ^ 2;
t1485 = -pkin(1) * t1532 + pkin(7) * t1648 - t1512;
t1511 = t1527 * g(1) - t1531 * g(2);
t1693 = pkin(7) * t1522;
t1546 = qJDD(1) * pkin(1) + t1532 * t1693 + t1511;
t1545 = t1523 * t1546;
t1631 = t1526 * t1485 - t1530 * t1545;
t1369 = (qJD(1) * t1488 * t1526 + g(3) * t1530) * t1522 - t1628 * pkin(2) - t1627 * pkin(8) + t1631;
t1534 = t1569 * pkin(3) - t1392 * qJ(4) + t1369;
t1757 = t1481 * t1702 + t1534;
t1663 = t1509 * t1529;
t1541 = -t1479 * t1663 + t1525 * t1569;
t1642 = t1526 * t1671;
t1710 = t1530 * t1541 - t1642;
t1664 = t1509 * t1525;
t1621 = -t1479 * t1664 - t1529 * t1569;
t1641 = t1530 * t1671;
t1712 = t1526 * t1541 + t1641;
t1735 = -t1522 * t1621 + t1523 * t1712;
t1756 = -t1527 * t1735 + t1531 * t1710;
t1620 = t1529 * t1428 + t1481 * t1664;
t1708 = t1530 * t1620 + t1642;
t1622 = t1525 * t1428 - t1481 * t1663;
t1711 = t1526 * t1620 - t1641;
t1736 = -t1522 * t1622 + t1523 * t1711;
t1755 = -t1527 * t1736 + t1531 * t1708;
t1754 = t1527 * t1710 + t1531 * t1735;
t1753 = t1527 * t1708 + t1531 * t1736;
t1524 = sin(qJ(5));
t1528 = cos(qJ(5));
t1446 = -t1528 * t1479 - t1509 * t1524;
t1448 = t1479 * t1524 - t1509 * t1528;
t1402 = t1448 * t1446;
t1425 = qJDD(5) + t1428;
t1722 = -t1402 + t1425;
t1752 = t1524 * t1722;
t1750 = t1525 * t1721;
t1719 = t1704 - t1476;
t1749 = t1526 * t1719;
t1745 = t1528 * t1722;
t1743 = t1529 * t1721;
t1742 = t1530 * t1719;
t1738 = t1522 * t1711 + t1523 * t1622;
t1737 = t1522 * t1712 + t1523 * t1621;
t1567 = (t1479 * t1525 + t1481 * t1529) * t1509;
t1568 = (t1479 * t1529 - t1481 * t1525) * t1509;
t1706 = -t1522 * t1567 + (t1484 * t1530 + t1526 * t1568) * t1523;
t1709 = -t1484 * t1526 + t1530 * t1568;
t1734 = t1527 * t1709 + t1531 * t1706;
t1733 = -t1527 * t1706 + t1531 * t1709;
t1519 = t1522 ^ 2;
t1662 = t1519 * t1532;
t1614 = t1634 * qJD(1);
t1732 = t1519 * (-t1523 * t1532 + t1614);
t1539 = -g(3) * t1661 + t1526 * t1545;
t1370 = t1628 * pkin(8) - t1627 * pkin(2) + (t1488 * t1689 + t1485) * t1530 + t1539;
t1574 = t1526 * t1614;
t1575 = t1530 * t1614;
t1692 = t1523 * g(3);
t1371 = -t1490 * pkin(2) - t1489 * pkin(8) - t1692 + (pkin(2) * t1574 - pkin(8) * t1575 - t1546) * t1522;
t1300 = t1525 * t1370 - t1529 * t1371;
t1437 = pkin(3) * t1479 - qJ(4) * t1481;
t1257 = t1484 * pkin(3) - t1703 * qJ(4) + t1481 * t1437 + qJDD(4) + t1300;
t1393 = -t1462 + t1428;
t1201 = t1393 * pkin(4) + t1549 * pkin(9) + t1257;
t1451 = pkin(4) * t1481 + pkin(9) * t1509;
t1636 = -pkin(3) * t1509 + t1702;
t1209 = t1630 * pkin(9) - t1476 * pkin(4) + (pkin(9) * qJD(3) - t1451 + t1636) * t1481 + t1534;
t1168 = -t1201 * t1528 + t1524 * t1209;
t1169 = t1201 * t1524 + t1209 * t1528;
t1723 = -t1528 * t1168 + t1524 * t1169;
t1301 = t1529 * t1370 + t1525 * t1371;
t1547 = -t1703 * pkin(3) - t1479 * t1437 + t1301;
t1542 = 0.2e1 * qJD(4) * t1509 - t1547;
t1717 = pkin(3) * t1440 - qJ(4) * (t1550 + t1484) - t1542;
t1690 = t1484 * qJ(4);
t1208 = -t1690 - t1569 * pkin(4) - t1476 * pkin(9) + (t1702 - t1451) * t1509 + t1547;
t1701 = pkin(3) + pkin(9);
t1716 = qJ(4) * t1208 - t1701 * t1723;
t1632 = t1524 * t1484 + t1528 * t1569;
t1475 = qJD(5) + t1481;
t1652 = qJD(5) - t1475;
t1315 = t1652 * t1448 - t1632;
t1543 = -t1528 * t1484 + t1524 * t1569;
t1317 = -t1652 * t1446 + t1543;
t1228 = -t1315 * t1524 - t1317 * t1528;
t1445 = t1446 ^ 2;
t1705 = t1448 ^ 2;
t1349 = -t1445 - t1705;
t1715 = qJ(4) * t1349 - t1701 * t1228 - t1723;
t1472 = t1475 ^ 2;
t1363 = -t1472 - t1445;
t1277 = t1524 * t1363 + t1745;
t1651 = qJD(5) + t1475;
t1314 = t1448 * t1651 - t1632;
t1657 = t1524 * t1208;
t1714 = qJ(4) * t1314 - t1701 * t1277 + t1657;
t1206 = t1528 * t1208;
t1638 = -t1472 - t1705;
t1340 = t1402 + t1425;
t1656 = t1524 * t1340;
t1283 = t1528 * t1638 - t1656;
t1348 = -t1446 * qJD(5) + t1543;
t1675 = t1475 * t1446;
t1318 = t1348 - t1675;
t1713 = qJ(4) * t1318 - t1701 * t1283 + t1206;
t1707 = t1484 * t1660 + t1523 * t1567 + t1568 * t1661;
t1700 = pkin(2) * t1526;
t1698 = pkin(3) * t1525;
t1697 = pkin(3) * t1529;
t1696 = pkin(4) * t1723;
t1695 = pkin(4) * t1208;
t1694 = pkin(4) * t1228;
t1685 = t1340 * t1528;
t1684 = t1369 * t1525;
t1683 = t1369 * t1529;
t1674 = t1475 * t1448;
t1673 = t1475 * t1524;
t1672 = t1475 * t1528;
t1508 = t1530 * t1526 * t1662;
t1486 = -t1508 + t1628;
t1669 = t1486 * t1526;
t1668 = t1486 * t1530;
t1487 = t1508 + t1628;
t1667 = t1487 * t1526;
t1666 = t1487 * t1530;
t1466 = t1522 * t1546 + t1692;
t1655 = t1526 * t1466;
t1654 = t1530 * t1466;
t1520 = t1526 ^ 2;
t1521 = t1530 ^ 2;
t1649 = t1520 + t1521;
t1646 = t1472 - t1705;
t1644 = t1525 * t1402;
t1643 = t1529 * t1402;
t1640 = t1520 * t1662;
t1639 = t1521 * t1662;
t1637 = qJ(4) * t1525 + pkin(2);
t1219 = t1300 * t1525 + t1529 * t1301;
t1629 = -t1511 * t1527 - t1531 * t1512;
t1507 = qJDD(1) * t1531 - t1527 * t1532;
t1625 = -pkin(6) * t1507 - g(3) * t1527;
t1246 = -t1542 - t1690;
t1624 = -pkin(3) * t1257 + qJ(4) * t1246;
t1390 = (-qJD(3) - t1509) * t1479 + t1551;
t1623 = -pkin(3) * t1390 - qJ(4) * t1721;
t1619 = pkin(4) * t1314 + t1206;
t1618 = pkin(4) * t1318 - t1657;
t1473 = -t1640 - t1627;
t1436 = -t1473 * t1526 - t1668;
t1617 = pkin(7) * t1436 - t1655;
t1494 = -t1627 - t1639;
t1444 = t1494 * t1530 - t1667;
t1616 = pkin(7) * t1444 + t1654;
t1123 = t1208 * t1529 + t1525 * t1723;
t1129 = t1524 * t1168 + t1528 * t1169;
t1613 = t1123 * t1526 - t1129 * t1530;
t1183 = t1246 * t1529 + t1257 * t1525;
t1258 = t1481 * t1636 + t1534;
t1612 = t1183 * t1526 - t1258 * t1530;
t1196 = t1228 * t1525 + t1349 * t1529;
t1230 = -t1315 * t1528 + t1524 * t1317;
t1611 = t1196 * t1526 - t1230 * t1530;
t1535 = -t1446 * t1651 + t1543;
t1227 = -t1524 * t1314 + t1528 * t1535;
t1401 = -t1445 + t1705;
t1200 = t1227 * t1525 + t1401 * t1529;
t1229 = -t1314 * t1528 - t1524 * t1535;
t1610 = t1200 * t1526 - t1229 * t1530;
t1211 = t1277 * t1525 + t1314 * t1529;
t1278 = t1363 * t1528 - t1752;
t1609 = t1211 * t1526 - t1278 * t1530;
t1217 = t1283 * t1525 + t1318 * t1529;
t1284 = -t1524 * t1638 - t1685;
t1608 = t1217 * t1526 - t1284 * t1530;
t1607 = t1219 * t1526 - t1369 * t1530;
t1412 = t1445 - t1472;
t1296 = t1412 * t1524 + t1685;
t1222 = t1296 * t1525 - t1315 * t1529;
t1298 = t1412 * t1528 - t1656;
t1606 = t1222 * t1526 - t1298 * t1530;
t1295 = t1528 * t1646 + t1752;
t1223 = t1295 * t1525 + t1317 * t1529;
t1297 = -t1524 * t1646 + t1745;
t1605 = t1223 * t1526 - t1297 * t1530;
t1347 = -qJD(5) * t1448 + t1632;
t1308 = -t1347 * t1528 - t1446 * t1673;
t1263 = -t1308 * t1525 - t1643;
t1309 = t1347 * t1524 - t1446 * t1672;
t1604 = t1263 * t1526 + t1309 * t1530;
t1310 = t1348 * t1524 + t1448 * t1672;
t1264 = t1310 * t1525 + t1643;
t1311 = t1348 * t1528 - t1448 * t1673;
t1603 = t1264 * t1526 - t1311 * t1530;
t1218 = -t1300 * t1529 + t1301 * t1525;
t1341 = (-t1446 * t1524 - t1448 * t1528) * t1475;
t1307 = t1341 * t1525 + t1425 * t1529;
t1342 = (-t1446 * t1528 + t1448 * t1524) * t1475;
t1602 = t1307 * t1526 - t1342 * t1530;
t1323 = t1390 * t1525 - t1743;
t1601 = t1323 * t1526 - t1758;
t1324 = -t1392 * t1525 + t1682;
t1600 = t1324 * t1526 - t1742;
t1325 = t1394 * t1525 + t1682;
t1599 = t1325 * t1526 - t1742;
t1326 = t1393 * t1525 - t1743;
t1598 = t1326 * t1526 - t1758;
t1593 = t1390 * t1530 + t1775;
t1592 = -t1393 * t1530 - t1775;
t1434 = g(3) * t1660 + t1631;
t1435 = t1530 * t1485 + t1539;
t1589 = -t1530 * t1434 + t1526 * t1435;
t1350 = t1434 * t1526 + t1435 * t1530;
t1498 = t1522 * t1575;
t1458 = t1498 + t1489;
t1497 = t1522 * t1574;
t1461 = t1490 - t1497;
t1588 = t1458 * t1530 + t1461 * t1526;
t1459 = -t1498 + t1489;
t1460 = t1490 + t1497;
t1587 = -t1459 * t1530 + t1460 * t1526;
t1586 = t1473 * t1530 - t1669;
t1493 = -t1627 + t1639;
t1585 = t1493 * t1526 + t1668;
t1492 = t1627 - t1640;
t1584 = t1492 * t1530 + t1667;
t1583 = t1494 * t1526 + t1666;
t1582 = t1511 * t1531 - t1512 * t1527;
t1581 = t1522 * t1628;
t1578 = -pkin(4) * t1277 + t1168;
t1576 = t1522 * t1614;
t1566 = pkin(4) * t1349 - t1129;
t1565 = -pkin(4) * t1283 + t1169;
t1104 = -t1701 * t1129 + t1695;
t1106 = -qJ(4) * t1129 + t1696;
t1122 = t1208 * t1525 - t1529 * t1723;
t1088 = -pkin(8) * t1122 - t1104 * t1525 + t1106 * t1529;
t1100 = -pkin(2) * t1122 - t1716;
t1103 = t1123 * t1530 + t1129 * t1526;
t1564 = pkin(7) * t1103 + t1088 * t1526 + t1100 * t1530;
t1117 = -t1701 * t1230 + t1566;
t1176 = -qJ(4) * t1230 + t1694;
t1195 = -t1228 * t1529 + t1349 * t1525;
t1105 = -pkin(8) * t1195 - t1117 * t1525 + t1176 * t1529;
t1108 = -pkin(2) * t1195 - t1715;
t1170 = t1196 * t1530 + t1230 * t1526;
t1563 = pkin(7) * t1170 + t1105 * t1526 + t1108 * t1530;
t1141 = -qJ(4) * t1278 - t1578;
t1153 = -t1701 * t1278 + t1619;
t1210 = -t1277 * t1529 + t1314 * t1525;
t1112 = -pkin(8) * t1210 + t1141 * t1529 - t1153 * t1525;
t1134 = -pkin(2) * t1210 - t1714;
t1179 = t1211 * t1530 + t1278 * t1526;
t1562 = pkin(7) * t1179 + t1112 * t1526 + t1134 * t1530;
t1144 = -qJ(4) * t1284 - t1565;
t1154 = -t1701 * t1284 + t1618;
t1216 = -t1283 * t1529 + t1318 * t1525;
t1113 = -pkin(8) * t1216 + t1144 * t1529 - t1154 * t1525;
t1136 = -pkin(2) * t1216 - t1713;
t1181 = t1217 * t1530 + t1284 * t1526;
t1561 = pkin(7) * t1181 + t1113 * t1526 + t1136 * t1530;
t1182 = t1246 * t1525 - t1257 * t1529;
t1148 = -pkin(2) * t1182 - t1624;
t1149 = -pkin(8) * t1182 + (-qJ(4) * t1529 + t1698) * t1258;
t1163 = t1183 * t1530 + t1258 * t1526;
t1560 = pkin(7) * t1163 + t1148 * t1530 + t1149 * t1526;
t1226 = -pkin(3) * t1410 + t1246;
t1231 = -qJ(4) * t1410 + t1257;
t1319 = -t1390 * t1529 - t1750;
t1162 = -pkin(8) * t1319 - t1226 * t1525 + t1231 * t1529;
t1232 = -pkin(2) * t1319 - t1623;
t1267 = t1323 * t1530 + t1759;
t1559 = pkin(7) * t1267 + t1162 * t1526 + t1232 * t1530;
t1225 = (-t1388 - t1665) * pkin(3) + t1757;
t1184 = -qJ(4) * t1682 - t1225 * t1525 + t1779;
t1537 = pkin(3) * t1549 - qJ(4) * t1433 + t1257;
t1189 = -t1537 + t1781;
t1558 = t1184 * t1526 + t1189 * t1530 - t1789;
t1224 = pkin(3) * t1665 - qJ(4) * t1394 - t1757;
t1187 = t1224 * t1529 + t1394 * t1698 - t1777;
t1190 = -t1717 - t1780;
t1557 = t1187 * t1526 + t1190 * t1530 - t1788;
t1242 = t1300 - t1781;
t1285 = t1684 - t1779;
t1556 = t1242 * t1530 + t1285 * t1526 + t1789;
t1245 = t1301 + t1780;
t1291 = t1683 + t1777;
t1555 = t1245 * t1530 + t1291 * t1526 + t1788;
t1403 = t1459 * t1526 + t1460 * t1530;
t1554 = pkin(7) * t1403 + t1350;
t1322 = -t1393 * t1529 - t1750;
t1188 = -pkin(8) * t1322 - t1218;
t1268 = t1326 * t1530 + t1759;
t1553 = pkin(7) * t1268 + t1188 * t1526 - t1322 * t1699;
t1193 = t1219 * t1530 + t1369 * t1526;
t1548 = pkin(7) * t1193 + t1218 * t1626;
t1544 = t1522 * t1662 + t1523 * t1576;
t1506 = qJDD(1) * t1527 + t1531 * t1532;
t1503 = t1523 * t1628;
t1496 = t1649 * t1662;
t1495 = (t1520 - t1521) * t1662;
t1491 = -pkin(6) * t1506 + g(3) * t1531;
t1465 = t1634 * t1649 * t1689;
t1457 = (t1647 + (0.2e1 * qJD(2) + t1688) * t1687) * t1522;
t1450 = t1530 * t1489 - t1520 * t1576;
t1449 = -t1526 * t1490 - t1521 * t1576;
t1443 = t1493 * t1530 - t1669;
t1442 = -t1492 * t1526 + t1666;
t1432 = (t1523 * t1489 + t1530 * t1544) * t1526;
t1431 = (t1522 * t1489 + t1530 * t1732) * t1526;
t1430 = (t1522 * t1490 - t1526 * t1732) * t1530;
t1429 = (t1523 * t1490 - t1526 * t1544) * t1530;
t1404 = -t1458 * t1526 + t1461 * t1530;
t1400 = t1522 * t1461 + t1523 * t1583;
t1399 = -t1522 * t1460 + t1523 * t1585;
t1398 = -t1522 * t1459 + t1523 * t1584;
t1397 = -t1523 * t1461 + t1522 * t1583;
t1396 = t1523 * t1460 + t1522 * t1585;
t1395 = t1523 * t1459 + t1522 * t1584;
t1383 = -t1522 * t1457 + t1523 * t1586;
t1382 = t1523 * t1457 + t1522 * t1586;
t1368 = -t1522 * t1495 + t1523 * t1588;
t1367 = t1522 * t1496 + t1523 * t1587;
t1366 = t1523 * t1495 + t1522 * t1588;
t1365 = -t1523 * t1496 + t1522 * t1587;
t1332 = t1522 * t1466 + t1523 * t1589;
t1331 = -t1523 * t1466 + t1522 * t1589;
t1321 = -t1394 * t1529 + t1764;
t1320 = t1392 * t1529 + t1764;
t1306 = -t1341 * t1529 + t1425 * t1525;
t1303 = t1393 * t1526 - t1774;
t1302 = -t1390 * t1526 + t1774;
t1294 = -t1655 + (-t1397 * t1522 - t1400 * t1523) * pkin(7);
t1288 = -t1654 + (-t1382 * t1522 - t1383 * t1523) * pkin(7);
t1287 = -pkin(1) * t1397 + t1522 * t1434 + t1523 * t1616;
t1286 = pkin(1) * t1400 - t1523 * t1434 + t1522 * t1616;
t1282 = t1325 * t1530 + t1749;
t1281 = t1324 * t1530 + t1749;
t1280 = -pkin(1) * t1382 + t1522 * t1435 + t1523 * t1617;
t1279 = pkin(1) * t1383 - t1523 * t1435 + t1522 * t1617;
t1266 = pkin(1) * t1332 + t1350 * t1693;
t1265 = pkin(7) * t1350 * t1523 - pkin(1) * t1331;
t1262 = -t1310 * t1529 + t1644;
t1261 = t1308 * t1529 - t1644;
t1260 = -pkin(1) * t1365 + t1523 * t1554;
t1259 = pkin(1) * t1367 + t1522 * t1554;
t1256 = pkin(2) * t1394 + t1684 + t1776;
t1253 = t1523 * t1592 - t1771;
t1252 = t1523 * t1593 + t1771;
t1249 = t1522 * t1592 + t1767;
t1248 = t1522 * t1593 - t1767;
t1247 = (-t1331 * t1522 - t1332 * t1523) * pkin(7);
t1244 = pkin(2) * t1388 - t1683 + t1778;
t1243 = (-t1365 * t1522 - t1367 * t1523) * pkin(7) - t1589;
t1237 = t1307 * t1530 + t1342 * t1526;
t1221 = -t1295 * t1529 + t1317 * t1525;
t1220 = -t1296 * t1529 - t1315 * t1525;
t1215 = -t1522 * t1321 + t1523 * t1599;
t1214 = -t1522 * t1320 + t1523 * t1600;
t1213 = t1523 * t1321 + t1522 * t1599;
t1212 = t1523 * t1320 + t1522 * t1600;
t1205 = -t1522 * t1322 + t1523 * t1598;
t1204 = -t1522 * t1319 + t1523 * t1601;
t1203 = t1523 * t1322 + t1522 * t1598;
t1202 = t1523 * t1319 + t1522 * t1601;
t1199 = -t1227 * t1529 + t1401 * t1525;
t1198 = t1264 * t1530 + t1311 * t1526;
t1197 = t1263 * t1530 - t1309 * t1526;
t1194 = -pkin(2) * t1369 + pkin(8) * t1219;
t1192 = -t1522 * t1306 + t1523 * t1602;
t1191 = t1523 * t1306 + t1522 * t1602;
t1186 = t1223 * t1530 + t1297 * t1526;
t1185 = t1222 * t1530 + t1298 * t1526;
t1180 = pkin(8) * t1326 + t1219 - t1765;
t1178 = -t1776 + t1525 * t1224 - (pkin(2) + t1697) * t1394;
t1177 = t1529 * t1225 - t1388 * t1637 - t1778;
t1175 = -t1522 * t1262 + t1523 * t1603;
t1174 = -t1522 * t1261 + t1523 * t1604;
t1173 = t1523 * t1262 + t1522 * t1603;
t1172 = t1523 * t1261 + t1522 * t1604;
t1171 = t1200 * t1530 + t1229 * t1526;
t1165 = -t1522 * t1218 + t1523 * t1607;
t1164 = t1523 * t1218 + t1522 * t1607;
t1161 = pkin(8) * t1323 + t1226 * t1529 + t1231 * t1525 - t1765;
t1160 = -t1522 * t1221 + t1523 * t1605;
t1159 = -t1522 * t1220 + t1523 * t1606;
t1158 = t1523 * t1221 + t1522 * t1605;
t1157 = t1523 * t1220 + t1522 * t1606;
t1156 = -t1522 * t1216 + t1523 * t1608;
t1155 = t1523 * t1216 + t1522 * t1608;
t1152 = -t1522 * t1210 + t1523 * t1609;
t1151 = t1523 * t1210 + t1522 * t1609;
t1150 = -t1526 * t1245 + t1530 * t1291 - t1803;
t1147 = -t1526 * t1242 + t1530 * t1285 - t1802;
t1146 = -t1522 * t1199 + t1523 * t1610;
t1145 = t1523 * t1199 + t1522 * t1610;
t1143 = -t1522 * t1256 + t1523 * t1555 - t1799;
t1142 = t1523 * t1256 + t1522 * t1555 + t1798;
t1140 = -t1522 * t1195 + t1523 * t1611;
t1139 = t1523 * t1195 + t1522 * t1611;
t1138 = -t1522 * t1244 + t1523 * t1556 - t1801;
t1137 = t1523 * t1244 + t1522 * t1556 + t1800;
t1135 = pkin(8) * t1183 + (-t1637 - t1697) * t1258;
t1133 = -t1522 * t1182 + t1523 * t1612;
t1132 = t1523 * t1182 + t1522 * t1612;
t1131 = t1322 * t1700 + t1530 * t1188 + (-t1203 * t1522 - t1205 * t1523) * pkin(7);
t1130 = t1530 * t1187 - t1526 * t1190 + t1803;
t1127 = t1530 * t1184 - t1526 * t1189 + t1802;
t1126 = -pkin(1) * t1203 - t1522 * t1180 + t1523 * t1553;
t1125 = pkin(1) * t1205 + t1523 * t1180 + t1522 * t1553;
t1124 = t1530 * t1162 - t1526 * t1232 + (-t1202 * t1522 - t1204 * t1523) * pkin(7);
t1121 = -t1522 * t1178 + t1523 * t1557 + t1799;
t1120 = t1523 * t1178 + t1522 * t1557 - t1798;
t1119 = -t1522 * t1177 + t1523 * t1558 + t1801;
t1118 = t1523 * t1177 + t1522 * t1558 - t1800;
t1116 = (-pkin(8) * t1530 + t1700) * t1218 + (-t1164 * t1522 - t1165 * t1523) * pkin(7);
t1115 = -pkin(1) * t1164 - t1522 * t1194 + t1523 * t1548;
t1114 = pkin(1) * t1165 + t1523 * t1194 + t1522 * t1548;
t1111 = -pkin(1) * t1202 - t1522 * t1161 + t1523 * t1559;
t1110 = pkin(1) * t1204 + t1523 * t1161 + t1522 * t1559;
t1109 = -pkin(2) * t1284 + pkin(8) * t1217 + t1144 * t1525 + t1154 * t1529;
t1107 = -pkin(2) * t1278 + pkin(8) * t1211 + t1141 * t1525 + t1153 * t1529;
t1102 = -pkin(2) * t1230 + pkin(8) * t1196 + t1117 * t1529 + t1176 * t1525;
t1101 = -t1526 * t1148 + t1530 * t1149 + (-t1132 * t1522 - t1133 * t1523) * pkin(7);
t1099 = -t1522 * t1122 + t1523 * t1613;
t1098 = t1523 * t1122 + t1522 * t1613;
t1097 = -pkin(1) * t1132 - t1522 * t1135 + t1523 * t1560;
t1096 = pkin(1) * t1133 + t1523 * t1135 + t1522 * t1560;
t1095 = t1530 * t1113 - t1526 * t1136 + (-t1155 * t1522 - t1156 * t1523) * pkin(7);
t1094 = t1530 * t1112 - t1526 * t1134 + (-t1151 * t1522 - t1152 * t1523) * pkin(7);
t1093 = -pkin(1) * t1155 - t1522 * t1109 + t1523 * t1561;
t1092 = pkin(1) * t1156 + t1523 * t1109 + t1522 * t1561;
t1091 = -pkin(1) * t1151 - t1522 * t1107 + t1523 * t1562;
t1090 = pkin(1) * t1152 + t1523 * t1107 + t1522 * t1562;
t1089 = t1530 * t1105 - t1526 * t1108 + (-t1139 * t1522 - t1140 * t1523) * pkin(7);
t1087 = -pkin(2) * t1129 + pkin(8) * t1123 + t1104 * t1529 + t1106 * t1525;
t1086 = -pkin(1) * t1139 - t1522 * t1102 + t1523 * t1563;
t1085 = pkin(1) * t1140 + t1523 * t1102 + t1522 * t1563;
t1084 = t1530 * t1088 - t1526 * t1100 + (-t1098 * t1522 - t1099 * t1523) * pkin(7);
t1083 = -pkin(1) * t1098 - t1522 * t1087 + t1523 * t1564;
t1082 = pkin(1) * t1099 + t1523 * t1087 + t1522 * t1564;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1507, 0, -t1506, 0, t1625, -t1491, -t1582, -pkin(6) * t1582, -t1432 * t1527 + t1450 * t1531, -t1368 * t1527 + t1404 * t1531, -t1398 * t1527 + t1442 * t1531, -t1429 * t1527 + t1449 * t1531, -t1399 * t1527 + t1443 * t1531, t1531 * t1465 + t1527 * t1581, t1531 * t1294 - t1527 * t1287 - pkin(6) * (t1400 * t1531 + t1444 * t1527), t1531 * t1288 - t1527 * t1280 - pkin(6) * (t1383 * t1531 + t1436 * t1527), t1531 * t1243 - t1527 * t1260 - pkin(6) * (t1367 * t1531 + t1403 * t1527), t1531 * t1247 - t1527 * t1265 - pkin(6) * (t1332 * t1531 + t1350 * t1527), t1755, -t1214 * t1527 + t1281 * t1531, -t1253 * t1527 + t1303 * t1531, t1756, t1793, t1733, -t1527 * t1138 + t1531 * t1147 - t1805, -t1527 * t1143 + t1531 * t1150 - t1807, t1531 * t1131 - t1527 * t1126 - pkin(6) * (t1205 * t1531 + t1268 * t1527), t1531 * t1116 - t1527 * t1115 - pkin(6) * (t1165 * t1531 + t1193 * t1527), t1733, -t1252 * t1527 + t1302 * t1531, -t1793, t1755, -t1215 * t1527 + t1282 * t1531, t1756, t1531 * t1124 - t1527 * t1111 - pkin(6) * (t1204 * t1531 + t1267 * t1527), -t1527 * t1119 + t1531 * t1127 + t1805, -t1527 * t1121 + t1531 * t1130 + t1807, t1531 * t1101 - t1527 * t1097 - pkin(6) * (t1133 * t1531 + t1163 * t1527), -t1175 * t1527 + t1198 * t1531, -t1146 * t1527 + t1171 * t1531, -t1160 * t1527 + t1186 * t1531, -t1174 * t1527 + t1197 * t1531, -t1159 * t1527 + t1185 * t1531, -t1192 * t1527 + t1237 * t1531, t1531 * t1094 - t1527 * t1091 - pkin(6) * (t1152 * t1531 + t1179 * t1527), t1531 * t1095 - t1527 * t1093 - pkin(6) * (t1156 * t1531 + t1181 * t1527), t1531 * t1089 - t1527 * t1086 - pkin(6) * (t1140 * t1531 + t1170 * t1527), t1531 * t1084 - t1527 * t1083 - pkin(6) * (t1099 * t1531 + t1103 * t1527); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1506, 0, t1507, 0, t1491, t1625, t1629, pkin(6) * t1629, t1432 * t1531 + t1450 * t1527, t1368 * t1531 + t1404 * t1527, t1398 * t1531 + t1442 * t1527, t1429 * t1531 + t1449 * t1527, t1399 * t1531 + t1443 * t1527, t1527 * t1465 - t1531 * t1581, t1527 * t1294 + t1531 * t1287 + pkin(6) * (-t1400 * t1527 + t1444 * t1531), t1527 * t1288 + t1531 * t1280 + pkin(6) * (-t1383 * t1527 + t1436 * t1531), t1527 * t1243 + t1531 * t1260 + pkin(6) * (-t1367 * t1527 + t1403 * t1531), t1527 * t1247 + t1531 * t1265 + pkin(6) * (-t1332 * t1527 + t1350 * t1531), t1753, t1214 * t1531 + t1281 * t1527, t1253 * t1531 + t1303 * t1527, t1754, -t1792, t1734, t1531 * t1138 + t1527 * t1147 - t1804, t1531 * t1143 + t1527 * t1150 - t1806, t1527 * t1131 + t1531 * t1126 + pkin(6) * (-t1205 * t1527 + t1268 * t1531), t1527 * t1116 + t1531 * t1115 + pkin(6) * (-t1165 * t1527 + t1193 * t1531), t1734, t1252 * t1531 + t1302 * t1527, t1792, t1753, t1215 * t1531 + t1282 * t1527, t1754, t1527 * t1124 + t1531 * t1111 + pkin(6) * (-t1204 * t1527 + t1267 * t1531), t1531 * t1119 + t1527 * t1127 + t1804, t1531 * t1121 + t1527 * t1130 + t1806, t1527 * t1101 + t1531 * t1097 + pkin(6) * (-t1133 * t1527 + t1163 * t1531), t1175 * t1531 + t1198 * t1527, t1146 * t1531 + t1171 * t1527, t1160 * t1531 + t1186 * t1527, t1174 * t1531 + t1197 * t1527, t1159 * t1531 + t1185 * t1527, t1192 * t1531 + t1237 * t1527, t1527 * t1094 + t1531 * t1091 + pkin(6) * (-t1152 * t1527 + t1179 * t1531), t1527 * t1095 + t1531 * t1093 + pkin(6) * (-t1156 * t1527 + t1181 * t1531), t1527 * t1089 + t1531 * t1086 + pkin(6) * (-t1140 * t1527 + t1170 * t1531), t1527 * t1084 + t1531 * t1083 + pkin(6) * (-t1099 * t1527 + t1103 * t1531); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1511, t1512, 0, 0, t1431, t1366, t1395, t1430, t1396, t1503, t1286, t1279, t1259, t1266, t1738, t1212, t1249, t1737, -t1251, t1707, t1137, t1142, t1125, t1114, t1707, t1248, t1251, t1738, t1213, t1737, t1110, t1118, t1120, t1096, t1173, t1145, t1158, t1172, t1157, t1191, t1090, t1092, t1085, t1082; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1532, 0, 0, -g(3), -t1511, 0, t1450, t1404, t1442, t1449, t1443, t1465, t1294, t1288, t1243, t1247, t1708, t1281, t1303, t1710, -t1305, t1709, t1147, t1150, t1131, t1116, t1709, t1302, t1305, t1708, t1282, t1710, t1124, t1127, t1130, t1101, t1198, t1171, t1186, t1197, t1185, t1237, t1094, t1095, t1089, t1084; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1532, 0, qJDD(1), 0, g(3), 0, -t1512, 0, t1432, t1368, t1398, t1429, t1399, -t1581, t1287, t1280, t1260, t1265, t1736, t1214, t1253, t1735, -t1255, t1706, t1138, t1143, t1126, t1115, t1706, t1252, t1255, t1736, t1215, t1735, t1111, t1119, t1121, t1097, t1175, t1146, t1160, t1174, t1159, t1192, t1091, t1093, t1086, t1083; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1511, t1512, 0, 0, t1431, t1366, t1395, t1430, t1396, t1503, t1286, t1279, t1259, t1266, t1738, t1212, t1249, t1737, -t1251, t1707, t1137, t1142, t1125, t1114, t1707, t1248, t1251, t1738, t1213, t1737, t1110, t1118, t1120, t1096, t1173, t1145, t1158, t1172, t1157, t1191, t1090, t1092, t1085, t1082; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1489, t1461, t1487, -t1498, t1493, t1498, 0, -t1466, t1434, 0, t1620, t1324, -t1359, t1541, -t1362, t1568, t1285, t1291, t1188, -pkin(8) * t1218, t1568, t1359, t1362, t1620, t1325, t1541, t1162, t1184, t1187, t1149, t1264, t1200, t1223, t1263, t1222, t1307, t1112, t1113, t1105, t1088; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1497, t1458, t1492, t1490, t1486, -t1497, t1466, 0, t1435, 0, -t1671, -t1719, -t1393, t1671, t1721, t1484, t1242, t1245, -pkin(2) * t1322, -pkin(2) * t1218, t1484, t1390, -t1721, -t1671, -t1719, t1671, t1232, t1189, t1190, t1148, -t1311, -t1229, -t1297, t1309, -t1298, -t1342, t1134, t1136, t1108, t1100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1508, t1495, t1459, t1508, t1460, t1628, -t1434, -t1435, 0, 0, t1622, t1320, t1356, t1621, -t1358, t1567, t1244, t1256, t1180, t1194, t1567, -t1356, t1358, t1622, t1321, t1621, t1161, t1177, t1178, t1135, t1262, t1199, t1221, t1261, t1220, t1306, t1107, t1109, t1102, t1087; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1428, t1388, -t1549, -t1462, -t1454, t1462, 0, t1369, t1300, 0, t1462, t1549, t1454, t1428, t1388, -t1462, t1231, -qJ(4) * t1388, t1224, -qJ(4) * t1258, t1402, t1401, t1317, -t1402, -t1315, t1425, t1141, t1144, t1176, t1106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1665, t1392, t1720, -t1569, -t1550, t1665, -t1369, 0, t1301, 0, t1665, -t1720, t1550, -t1665, -t1394, -t1569, t1226, t1225, -pkin(3) * t1394, -pkin(3) * t1258, -t1310, -t1227, -t1295, t1308, -t1296, -t1341, t1153, t1154, t1117, t1104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1671, t1719, t1393, -t1671, -t1721, -t1484, -t1300, -t1301, 0, 0, -t1484, -t1390, t1721, t1671, t1719, -t1671, t1623, t1537, t1717, t1624, t1311, t1229, t1297, -t1309, t1298, t1342, t1714, t1713, t1715, t1716; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1484, -t1390, t1721, t1671, t1719, -t1671, 0, t1257, t1246, 0, t1311, t1229, t1297, -t1309, t1298, t1342, -pkin(9) * t1277 + t1657, -pkin(9) * t1283 + t1206, -pkin(9) * t1228 - t1723, -pkin(9) * t1723; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1462, -t1549, -t1454, -t1428, -t1388, t1462, -t1257, 0, t1258, 0, -t1402, -t1401, -t1317, t1402, t1315, -t1425, t1578, t1565, -t1694, -t1696; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1665, t1720, -t1550, t1665, t1394, t1569, -t1246, -t1258, 0, 0, t1310, t1227, t1295, -t1308, t1296, t1341, pkin(9) * t1278 - t1619, pkin(9) * t1284 - t1618, pkin(9) * t1230 - t1566, pkin(9) * t1129 - t1695; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1348, -t1314, t1722, t1675, t1412, -t1675, 0, t1208, t1168, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1674, t1535, t1646, t1347, t1340, -t1674, -t1208, 0, t1169, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1402, t1401, t1317, -t1402, -t1315, t1425, -t1168, -t1169, 0, 0;];
m_new_reg = t1;
