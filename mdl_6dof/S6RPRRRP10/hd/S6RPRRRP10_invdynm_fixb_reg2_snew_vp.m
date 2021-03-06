% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6RPRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 01:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6RPRRRP10_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP10_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP10_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP10_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:58:56
% EndTime: 2019-05-06 01:59:29
% DurationCPUTime: 34.81s
% Computational Cost: add. (119826->795), mult. (234766->894), div. (0->0), fcn. (157590->8), ass. (0->516)
t1515 = sin(qJ(4));
t1519 = cos(qJ(4));
t1520 = cos(qJ(3));
t1664 = qJD(1) * t1520;
t1475 = -t1519 * qJD(3) + t1515 * t1664;
t1477 = qJD(3) * t1515 + t1519 * t1664;
t1514 = sin(qJ(5));
t1518 = cos(qJ(5));
t1427 = -t1475 * t1514 + t1477 * t1518;
t1424 = t1427 ^ 2;
t1516 = sin(qJ(3));
t1506 = t1516 * qJD(1);
t1503 = t1506 + qJD(4);
t1496 = qJD(5) + t1503;
t1678 = t1496 ^ 2;
t1356 = t1678 + t1424;
t1612 = qJD(3) * t1664;
t1621 = t1516 * qJDD(1);
t1480 = -t1612 - t1621;
t1473 = qJDD(4) - t1480;
t1469 = qJDD(5) + t1473;
t1425 = t1518 * t1475 + t1477 * t1514;
t1644 = t1427 * t1425;
t1690 = t1469 + t1644;
t1651 = t1690 * t1514;
t1251 = t1356 * t1518 + t1651;
t1650 = t1690 * t1518;
t1287 = t1356 * t1514 - t1650;
t1219 = t1251 * t1515 + t1287 * t1519;
t1613 = qJD(3) * t1506;
t1620 = t1520 * qJDD(1);
t1558 = t1613 - t1620;
t1530 = -t1515 * qJDD(3) + t1519 * t1558;
t1420 = -t1475 * qJD(4) - t1530;
t1529 = t1519 * qJDD(3) + t1515 * t1558;
t1527 = -t1477 * qJD(4) + t1529;
t1320 = -t1425 * qJD(5) + t1518 * t1420 + t1514 * t1527;
t1645 = t1425 * t1496;
t1693 = t1320 - t1645;
t1160 = t1219 * t1516 - t1520 * t1693;
t1189 = t1251 * t1519 - t1287 * t1515;
t1517 = sin(qJ(1));
t1521 = cos(qJ(1));
t1789 = pkin(6) * (t1160 * t1521 + t1189 * t1517);
t1788 = pkin(6) * (t1160 * t1517 - t1189 * t1521);
t1787 = pkin(2) * t1160;
t1786 = pkin(7) * t1160;
t1162 = t1219 * t1520 + t1516 * t1693;
t1785 = pkin(7) * t1162;
t1784 = qJ(2) * t1162;
t1677 = pkin(7) + pkin(1);
t1783 = t1162 * t1677;
t1598 = -t1514 * t1420 + t1518 * t1527;
t1319 = qJD(5) * t1427 - t1598;
t1407 = t1496 * t1427;
t1694 = t1319 + t1407;
t1204 = -t1694 * t1514 + t1518 * t1693;
t1654 = t1693 * t1514;
t1206 = t1694 * t1518 + t1654;
t1142 = t1204 * t1515 + t1206 * t1519;
t1679 = t1425 ^ 2;
t1692 = t1424 - t1679;
t1126 = t1142 * t1516 + t1520 * t1692;
t1138 = -t1204 * t1519 + t1206 * t1515;
t1782 = t1126 * t1517 + t1138 * t1521;
t1781 = t1126 * t1521 - t1138 * t1517;
t1778 = -qJ(2) * t1189 - t1160 * t1677;
t1402 = t1679 - t1678;
t1299 = t1402 * t1514 + t1650;
t1303 = t1402 * t1518 - t1651;
t1228 = t1299 * t1515 - t1303 * t1519;
t1271 = t1319 - t1407;
t1169 = t1228 * t1516 - t1271 * t1520;
t1225 = t1299 * t1519 + t1303 * t1515;
t1777 = t1169 * t1517 - t1225 * t1521;
t1776 = t1169 * t1521 + t1225 * t1517;
t1775 = pkin(2) * t1189;
t1774 = pkin(3) * t1189;
t1773 = pkin(8) * t1189;
t1765 = pkin(3) * t1693 - pkin(8) * t1219;
t1128 = t1142 * t1520 - t1516 * t1692;
t1173 = t1228 * t1520 + t1271 * t1516;
t1403 = t1424 - t1678;
t1689 = -t1644 + t1469;
t1649 = t1689 * t1514;
t1726 = -t1518 * t1403 + t1649;
t1339 = t1518 * t1689;
t1727 = t1403 * t1514 + t1339;
t1736 = t1515 * t1727 + t1519 * t1726;
t1688 = t1645 + t1320;
t1735 = -t1515 * t1726 + t1519 * t1727;
t1752 = t1516 * t1735 - t1520 * t1688;
t1764 = t1517 * t1752 + t1521 * t1736;
t1763 = t1517 * t1736 - t1521 * t1752;
t1687 = -t1678 - t1679;
t1699 = t1518 * t1687 - t1649;
t1702 = t1514 * t1687 + t1339;
t1721 = -t1515 * t1702 + t1519 * t1699;
t1738 = t1516 * t1721 - t1520 * t1694;
t1762 = pkin(2) * t1738;
t1761 = pkin(4) * t1251;
t1737 = t1516 * t1694 + t1520 * t1721;
t1760 = pkin(7) * t1737;
t1759 = pkin(7) * t1738;
t1758 = pkin(9) * t1251;
t1757 = pkin(9) * t1287;
t1756 = qJ(2) * t1737;
t1753 = t1677 * t1737;
t1751 = t1516 * t1688 + t1520 * t1735;
t1720 = t1515 * t1699 + t1519 * t1702;
t1750 = qJ(2) * t1720 - t1677 * t1738;
t1749 = pkin(6) * (t1517 * t1738 + t1521 * t1720);
t1748 = pkin(6) * (t1517 * t1720 - t1521 * t1738);
t1747 = pkin(2) * t1720;
t1746 = pkin(3) * t1720;
t1745 = pkin(8) * t1720;
t1739 = -pkin(3) * t1694 + pkin(8) * t1721;
t1324 = -t1679 - t1424;
t1734 = pkin(3) * t1324;
t1733 = pkin(4) * t1324;
t1732 = pkin(4) * t1702;
t1731 = pkin(9) * t1699;
t1730 = pkin(9) * t1702;
t1729 = t1324 * t1516;
t1728 = t1324 * t1520;
t1547 = (-t1425 * t1514 - t1427 * t1518) * t1496;
t1635 = t1496 * t1514;
t1400 = t1427 * t1635;
t1634 = t1496 * t1518;
t1616 = t1425 * t1634;
t1573 = t1400 - t1616;
t1683 = t1515 * t1573 + t1519 * t1547;
t1682 = -t1515 * t1547 + t1519 * t1573;
t1700 = -t1469 * t1520 + t1516 * t1682;
t1725 = t1517 * t1700 + t1521 * t1683;
t1559 = t1319 * t1514 + t1616;
t1574 = -t1518 * t1319 + t1425 * t1635;
t1680 = t1515 * t1559 + t1519 * t1574;
t1617 = t1520 * t1644;
t1681 = -t1515 * t1574 + t1519 * t1559;
t1701 = t1516 * t1681 + t1617;
t1724 = t1517 * t1701 + t1521 * t1680;
t1723 = t1517 * t1683 - t1521 * t1700;
t1722 = t1517 * t1680 - t1521 * t1701;
t1717 = qJ(6) * t1693;
t1435 = t1477 * t1475;
t1691 = -t1435 + t1473;
t1716 = t1515 * t1691;
t1709 = t1519 * t1691;
t1557 = 0.2e1 * t1613 - t1620;
t1444 = t1557 * t1520;
t1618 = t1516 * t1644;
t1698 = t1520 * t1681 - t1618;
t1697 = t1469 * t1516 + t1520 * t1682;
t1663 = qJD(2) * qJD(1);
t1509 = 0.2e1 * t1663;
t1494 = t1521 * g(1) + t1517 * g(2);
t1511 = qJDD(1) * qJ(2);
t1563 = t1494 - t1511;
t1576 = -t1480 + t1612;
t1523 = qJD(1) ^ 2;
t1695 = t1523 * t1677;
t1375 = pkin(3) * t1576 + pkin(8) * t1557 + t1509 - t1563 - t1695;
t1493 = t1517 * g(1) - t1521 * g(2);
t1577 = qJDD(2) - t1493;
t1546 = -t1523 * qJ(2) + t1577;
t1454 = -qJDD(1) * t1677 + t1546;
t1432 = t1520 * g(3) - t1516 * t1454;
t1522 = qJD(3) ^ 2;
t1669 = pkin(8) * t1520;
t1673 = pkin(3) * t1516;
t1561 = t1523 * (-t1669 + t1673);
t1396 = -t1522 * pkin(3) + qJDD(3) * pkin(8) - t1516 * t1561 - t1432;
t1311 = -t1519 * t1375 + t1515 * t1396;
t1312 = t1515 * t1375 + t1519 * t1396;
t1232 = -t1311 * t1519 + t1312 * t1515;
t1619 = pkin(3) * t1520 + pkin(2);
t1696 = t1232 * (pkin(8) * t1516 + t1619);
t1233 = t1515 * t1311 + t1519 * t1312;
t1431 = t1516 * g(3) + t1520 * t1454;
t1373 = t1520 * t1431 - t1516 * t1432;
t1455 = t1503 * t1475;
t1386 = t1455 + t1420;
t1259 = t1320 * t1514 + t1427 * t1634;
t1260 = t1320 * t1518 - t1400;
t1198 = t1259 * t1519 + t1260 * t1515;
t1201 = -t1259 * t1515 + t1260 * t1519;
t1560 = t1201 * t1516 - t1617;
t1686 = t1517 * t1198 - t1521 * t1560;
t1685 = t1521 * t1198 + t1517 * t1560;
t1384 = (-qJD(4) + t1503) * t1477 + t1529;
t1315 = t1384 * t1515 - t1386 * t1519;
t1186 = -pkin(8) * t1315 - t1232;
t1684 = t1516 * t1186 - t1315 * t1619;
t1471 = t1475 ^ 2;
t1472 = t1477 ^ 2;
t1501 = t1503 ^ 2;
t1676 = pkin(2) * t1373;
t1533 = t1563 - 0.2e1 * t1663;
t1449 = t1533 + t1695;
t1675 = pkin(2) * t1449;
t1512 = t1516 ^ 2;
t1513 = t1520 ^ 2;
t1622 = t1512 + t1513;
t1484 = t1622 * qJDD(1);
t1674 = pkin(2) * t1484;
t1236 = pkin(4) * t1691 - pkin(9) * t1386 - t1311;
t1580 = pkin(4) * t1503 - pkin(9) * t1477;
t1243 = -t1471 * pkin(4) + pkin(9) * t1527 - t1503 * t1580 + t1312;
t1175 = -t1518 * t1236 + t1243 * t1514;
t1176 = t1514 * t1236 + t1518 * t1243;
t1116 = -t1175 * t1518 + t1176 * t1514;
t1672 = pkin(4) * t1116;
t1261 = t1518 * t1688;
t1273 = (-qJD(5) + t1496) * t1427 + t1598;
t1205 = t1273 * t1514 - t1261;
t1671 = pkin(4) * t1205;
t1670 = pkin(5) * t1518;
t1668 = t1319 * pkin(5);
t1667 = qJDD(1) * pkin(1);
t1666 = qJ(6) * t1518;
t1665 = qJD(1) * qJD(3);
t1662 = qJD(6) * t1496;
t1661 = t1116 * t1515;
t1660 = t1116 * t1519;
t1655 = t1688 * t1514;
t1395 = qJDD(3) * pkin(3) + t1522 * pkin(8) - t1520 * t1561 + t1431;
t1289 = pkin(4) * t1527 + t1471 * pkin(9) - t1477 * t1580 + t1395;
t1653 = t1289 * t1514;
t1652 = t1289 * t1518;
t1414 = t1435 + t1473;
t1647 = t1414 * t1515;
t1646 = t1414 * t1519;
t1641 = t1484 * t1517;
t1640 = t1484 * t1521;
t1502 = t1516 * t1523 * t1520;
t1491 = t1502 + qJDD(3);
t1639 = t1491 * t1516;
t1638 = t1491 * t1520;
t1492 = qJDD(3) - t1502;
t1637 = t1492 * t1516;
t1636 = t1492 * t1520;
t1633 = t1503 * t1515;
t1632 = t1503 * t1519;
t1631 = t1512 * t1523;
t1630 = t1513 * t1523;
t1388 = t1515 * t1395;
t1627 = t1516 * t1449;
t1389 = t1519 * t1395;
t1436 = t1520 * t1449;
t1482 = 0.2e1 * t1662;
t1359 = pkin(5) * t1425 - qJ(6) * t1427;
t1569 = -pkin(5) * t1678 + t1469 * qJ(6) - t1425 * t1359 + t1176;
t1147 = t1482 + t1569;
t1149 = -t1469 * pkin(5) - qJ(6) * t1678 + t1359 * t1427 + qJDD(6) + t1175;
t1626 = -pkin(5) * t1149 + qJ(6) * t1147;
t1625 = t1520 * t1186 + t1315 * t1673;
t1624 = -pkin(5) * t1688 - qJ(6) * t1271;
t1623 = pkin(3) * t1395 + pkin(8) * t1233;
t1615 = t1516 * t1435;
t1614 = t1520 * t1435;
t1433 = -t1472 - t1501;
t1358 = -t1433 * t1515 - t1646;
t1387 = (qJD(4) + t1503) * t1475 + t1530;
t1611 = pkin(3) * t1387 + pkin(8) * t1358 - t1388;
t1421 = -t1501 - t1471;
t1349 = t1421 * t1519 - t1716;
t1456 = t1503 * t1477;
t1383 = -t1456 + t1527;
t1610 = pkin(3) * t1383 + pkin(8) * t1349 + t1389;
t1609 = -qJ(6) * t1514 - pkin(4);
t1095 = t1147 * t1518 + t1149 * t1514;
t1525 = -pkin(5) * t1407 + 0.2e1 * qJD(6) * t1427 + t1289;
t1524 = t1525 + t1717;
t1155 = t1524 - t1668;
t1064 = pkin(9) * t1095 + (-t1609 + t1670) * t1155;
t1094 = t1147 * t1514 - t1149 * t1518;
t1071 = t1094 * t1519 + t1095 * t1515;
t1075 = -pkin(9) * t1094 + (-pkin(5) * t1514 + t1666) * t1155;
t1039 = -pkin(8) * t1071 - t1064 * t1515 + t1075 * t1519;
t1583 = pkin(4) * t1094 + t1626;
t1048 = -pkin(3) * t1071 - t1583;
t1608 = t1520 * t1039 - t1516 * t1048;
t1132 = -pkin(5) * t1324 + t1147;
t1133 = -qJ(6) * t1324 + t1149;
t1207 = -t1271 * t1518 + t1655;
t1080 = pkin(9) * t1207 + t1132 * t1518 + t1133 * t1514 - t1733;
t1203 = -t1271 * t1514 - t1261;
t1085 = -pkin(9) * t1203 - t1132 * t1514 + t1133 * t1518;
t1139 = t1203 * t1519 + t1207 * t1515;
t1050 = -pkin(8) * t1139 - t1080 * t1515 + t1085 * t1519;
t1582 = pkin(4) * t1203 + t1624;
t1090 = -pkin(3) * t1139 - t1582;
t1607 = t1520 * t1050 - t1516 * t1090;
t1117 = t1175 * t1514 + t1518 * t1176;
t1082 = t1117 * t1515 + t1660;
t1106 = pkin(4) * t1289 + pkin(9) * t1117;
t1052 = -pkin(8) * t1082 - pkin(9) * t1660 - t1106 * t1515;
t1065 = -pkin(3) * t1082 - t1672;
t1606 = t1520 * t1052 - t1516 * t1065;
t1209 = t1273 * t1518 + t1655;
t1089 = pkin(9) * t1209 + t1117 - t1733;
t1093 = -pkin(9) * t1205 - t1116;
t1141 = t1205 * t1519 + t1209 * t1515;
t1058 = -pkin(8) * t1141 - t1089 * t1515 + t1093 * t1519;
t1109 = -pkin(3) * t1141 - t1671;
t1605 = t1520 * t1058 - t1516 * t1109;
t1134 = t1525 - t1668 + 0.2e1 * t1717;
t1097 = -t1757 + t1514 * t1134 + (pkin(4) + t1670) * t1693;
t1108 = -pkin(5) * t1654 + t1134 * t1518 - t1758;
t1067 = -t1097 * t1515 + t1108 * t1519 - t1773;
t1544 = pkin(5) * t1356 + qJ(6) * t1690 + t1569;
t1531 = t1544 + t1761;
t1091 = -t1531 - 0.2e1 * t1662 - t1774;
t1604 = t1520 * t1067 - t1516 * t1091;
t1135 = (-t1319 - t1694) * pkin(5) + t1524;
t1102 = t1518 * t1135 + t1609 * t1694 + t1731;
t1113 = -t1135 * t1514 - t1666 * t1694 - t1730;
t1069 = -t1102 * t1515 + t1113 * t1519 - t1745;
t1528 = pkin(5) * t1689 + qJ(6) * t1687 - t1149;
t1526 = t1528 + t1732;
t1098 = -t1526 - t1746;
t1603 = t1520 * t1069 - t1516 * t1098;
t1159 = -pkin(4) * t1694 + t1652 + t1731;
t1211 = -t1653 - t1730;
t1100 = -t1159 * t1515 + t1211 * t1519 - t1745;
t1562 = -t1175 + t1732;
t1114 = -t1562 - t1746;
t1602 = t1520 * t1100 - t1516 * t1114;
t1165 = -pkin(4) * t1693 - t1653 + t1757;
t1221 = -t1652 + t1758;
t1104 = -t1165 * t1515 + t1221 * t1519 + t1773;
t1581 = -t1176 - t1761;
t1118 = -t1581 + t1774;
t1601 = t1520 * t1104 - t1516 * t1118;
t1348 = t1421 * t1515 + t1709;
t1242 = -pkin(3) * t1348 + t1311;
t1292 = -pkin(8) * t1348 - t1388;
t1600 = -t1516 * t1242 + t1520 * t1292;
t1357 = t1433 * t1519 - t1647;
t1244 = -pkin(3) * t1357 + t1312;
t1295 = -pkin(8) * t1357 - t1389;
t1599 = -t1516 * t1244 + t1520 * t1295;
t1457 = t1523 * pkin(1) + t1533;
t1458 = -t1546 + t1667;
t1596 = -t1521 * t1457 - t1458 * t1517;
t1595 = -t1493 * t1517 - t1521 * t1494;
t1594 = t1517 * t1502;
t1593 = t1521 * t1502;
t1072 = -t1094 * t1515 + t1095 * t1519;
t1592 = pkin(3) * t1155 + pkin(8) * t1072 + t1519 * t1064 + t1515 * t1075;
t1143 = -t1203 * t1515 + t1207 * t1519;
t1591 = pkin(8) * t1143 + t1519 * t1080 + t1515 * t1085 - t1734;
t1145 = -t1205 * t1515 + t1209 * t1519;
t1590 = pkin(8) * t1145 + t1519 * t1089 + t1515 * t1093 - t1734;
t1589 = t1519 * t1097 + t1515 * t1108 + t1765;
t1588 = t1519 * t1102 + t1515 * t1113 + t1739;
t1587 = t1519 * t1159 + t1515 * t1211 + t1739;
t1586 = t1519 * t1165 + t1515 * t1221 - t1765;
t1317 = t1384 * t1519 + t1386 * t1515;
t1411 = t1471 + t1472;
t1585 = pkin(3) * t1411 + pkin(8) * t1317 + t1233;
t1212 = t1233 * t1516 + t1395 * t1520;
t1584 = -pkin(2) * t1212 - t1623;
t1485 = qJDD(1) * t1521 - t1517 * t1523;
t1579 = pkin(6) * t1485 + g(3) * t1517;
t1486 = qJDD(1) * t1517 + t1521 * t1523;
t1578 = -pkin(6) * t1486 + g(3) * t1521;
t1575 = t1520 * t1201 + t1618;
t1479 = 0.2e1 * t1612 + t1621;
t1572 = pkin(2) * t1479 - t1436;
t1571 = -pkin(2) * t1557 + t1627;
t1374 = -t1431 * t1516 - t1432 * t1520;
t1568 = t1457 * t1517 - t1458 * t1521;
t1567 = t1493 * t1521 - t1494 * t1517;
t1281 = t1349 * t1516 + t1383 * t1520;
t1566 = -pkin(2) * t1281 - t1610;
t1290 = t1358 * t1516 + t1387 * t1520;
t1565 = -pkin(2) * t1290 - t1611;
t1500 = -t1522 - t1630;
t1441 = t1500 * t1520 - t1639;
t1564 = -pkin(2) * t1441 - t1432;
t1083 = t1117 * t1519 - t1661;
t1556 = pkin(3) * t1289 + pkin(8) * t1083 - pkin(9) * t1661 + t1519 * t1106;
t1061 = t1072 * t1516 + t1155 * t1520;
t1555 = -pkin(2) * t1061 - t1592;
t1119 = t1143 * t1516 - t1728;
t1554 = -pkin(2) * t1119 - t1591;
t1120 = t1145 * t1516 - t1728;
t1553 = -pkin(2) * t1120 - t1590;
t1552 = -t1589 + t1787;
t1551 = -t1588 - t1762;
t1550 = -t1587 - t1762;
t1549 = -t1586 - t1787;
t1249 = t1317 * t1516 + t1411 * t1520;
t1548 = -pkin(2) * t1249 - t1585;
t1498 = -t1522 - t1631;
t1439 = t1498 * t1516 + t1636;
t1545 = -pkin(2) * t1439 - t1431;
t1543 = pkin(2) * t1071 - t1516 * t1039 - t1520 * t1048;
t1542 = pkin(2) * t1139 - t1516 * t1050 - t1520 * t1090;
t1541 = pkin(2) * t1082 - t1516 * t1052 - t1520 * t1065;
t1540 = pkin(2) * t1141 - t1516 * t1058 - t1520 * t1109;
t1539 = -t1516 * t1067 - t1520 * t1091 + t1775;
t1538 = -t1516 * t1069 - t1520 * t1098 + t1747;
t1537 = -t1516 * t1100 - t1520 * t1114 + t1747;
t1536 = -t1516 * t1104 - t1520 * t1118 - t1775;
t1535 = pkin(2) * t1348 - t1520 * t1242 - t1516 * t1292;
t1534 = pkin(2) * t1357 - t1520 * t1244 - t1516 * t1295;
t1077 = t1083 * t1516 + t1289 * t1520;
t1532 = -pkin(2) * t1077 - t1556;
t1499 = t1522 - t1630;
t1497 = -t1522 + t1631;
t1488 = (-t1512 + t1513) * t1523;
t1487 = t1622 * t1523;
t1474 = t1622 * t1665;
t1470 = t1577 - 0.2e1 * t1667;
t1467 = -t1494 + t1509 + 0.2e1 * t1511;
t1453 = -t1472 + t1501;
t1452 = t1471 - t1501;
t1451 = t1513 * t1665 - t1516 * t1558;
t1450 = t1480 * t1520 + t1512 * t1665;
t1446 = -t1500 * t1516 - t1638;
t1445 = -t1499 * t1516 + t1636;
t1443 = t1498 * t1520 - t1637;
t1442 = t1497 * t1520 - t1639;
t1440 = t1499 * t1520 + t1637;
t1438 = t1497 * t1516 + t1638;
t1437 = t1576 * t1516;
t1434 = t1472 - t1471;
t1430 = -t1479 * t1520 + t1516 * t1557;
t1429 = -t1479 * t1516 - t1444;
t1412 = pkin(1) * t1458 - qJ(2) * t1457;
t1398 = (-t1475 * t1519 + t1477 * t1515) * t1503;
t1397 = (-t1475 * t1515 - t1477 * t1519) * t1503;
t1385 = -t1455 + t1420;
t1382 = -t1456 - t1527;
t1379 = t1420 * t1519 - t1477 * t1633;
t1378 = t1420 * t1515 + t1477 * t1632;
t1377 = t1475 * t1632 - t1515 * t1527;
t1376 = -t1475 * t1633 - t1519 * t1527;
t1367 = t1398 * t1520 + t1473 * t1516;
t1366 = t1398 * t1516 - t1473 * t1520;
t1363 = t1452 * t1519 - t1647;
t1362 = -t1453 * t1515 + t1709;
t1361 = t1452 * t1515 + t1646;
t1360 = t1453 * t1519 + t1716;
t1354 = pkin(2) * t1487 + t1374;
t1351 = -qJ(2) * t1446 - t1564;
t1350 = -qJ(2) * t1443 - t1545;
t1334 = -t1443 * t1677 + t1572;
t1333 = -t1446 * t1677 + t1571;
t1332 = -qJ(2) * t1557 - t1441 * t1677 - t1436;
t1331 = qJ(2) * t1479 - t1439 * t1677 - t1627;
t1330 = t1379 * t1520 + t1615;
t1329 = t1377 * t1520 - t1615;
t1328 = t1379 * t1516 - t1614;
t1327 = t1377 * t1516 + t1614;
t1326 = -qJ(2) * t1487 + t1484 * t1677 - t1373;
t1316 = t1383 * t1519 - t1385 * t1515;
t1314 = t1383 * t1515 + t1385 * t1519;
t1307 = t1363 * t1520 - t1382 * t1516;
t1306 = t1362 * t1520 + t1386 * t1516;
t1305 = t1363 * t1516 + t1382 * t1520;
t1304 = t1362 * t1516 - t1386 * t1520;
t1294 = -qJ(2) * t1374 + t1676;
t1291 = t1358 * t1520 - t1387 * t1516;
t1282 = t1349 * t1520 - t1383 * t1516;
t1280 = t1316 * t1520 + t1434 * t1516;
t1279 = t1316 * t1516 - t1434 * t1520;
t1269 = -t1374 * t1677 - t1675;
t1268 = -qJ(2) * t1449 - t1373 * t1677;
t1250 = t1317 * t1520 - t1411 * t1516;
t1230 = t1232 * t1673;
t1213 = t1233 * t1520 - t1395 * t1516;
t1157 = -qJ(2) * t1291 - t1565;
t1156 = -qJ(2) * t1282 - t1566;
t1131 = -t1291 * t1677 + t1534;
t1130 = qJ(2) * t1357 - t1290 * t1677 + t1599;
t1125 = -t1282 * t1677 + t1535;
t1124 = qJ(2) * t1348 - t1281 * t1677 + t1600;
t1123 = -qJ(2) * t1250 - t1548;
t1122 = t1145 * t1520 + t1729;
t1121 = t1143 * t1520 + t1729;
t1115 = -qJ(2) * t1213 - t1584;
t1112 = -t1250 * t1677 - t1684;
t1111 = qJ(2) * t1315 - t1249 * t1677 + t1625;
t1087 = -t1213 * t1677 + t1696;
t1086 = t1230 + (qJ(2) - t1669) * t1232 - t1677 * t1212;
t1078 = t1083 * t1520 - t1289 * t1516;
t1076 = -t1549 - t1784;
t1073 = -t1550 - t1756;
t1062 = t1072 * t1520 - t1155 * t1516;
t1060 = t1536 - t1783;
t1059 = t1601 + t1778;
t1056 = t1537 - t1753;
t1055 = t1602 + t1750;
t1054 = -t1551 - t1756;
t1053 = -t1552 + t1784;
t1047 = t1538 - t1753;
t1046 = t1603 + t1750;
t1045 = -qJ(2) * t1122 - t1553;
t1044 = t1539 + t1783;
t1043 = t1604 - t1778;
t1042 = -qJ(2) * t1121 - t1554;
t1041 = -t1122 * t1677 + t1540;
t1040 = qJ(2) * t1141 - t1120 * t1677 + t1605;
t1037 = -t1121 * t1677 + t1542;
t1036 = qJ(2) * t1139 - t1119 * t1677 + t1607;
t1035 = -qJ(2) * t1078 - t1532;
t1034 = -t1078 * t1677 + t1541;
t1033 = qJ(2) * t1082 - t1077 * t1677 + t1606;
t1032 = -qJ(2) * t1062 - t1555;
t1031 = -t1062 * t1677 + t1543;
t1030 = qJ(2) * t1071 - t1061 * t1677 + t1608;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1485, 0, -t1486, 0, -t1579, -t1578, -t1567, -pkin(6) * t1567, 0, -t1485, t1486, 0, 0, 0, t1568, t1579, t1578, pkin(6) * t1568 + (-pkin(1) * t1517 + qJ(2) * t1521) * g(3), t1451 * t1517 + t1593, t1429 * t1517 + t1488 * t1521, t1440 * t1517 + t1521 * t1620, t1450 * t1517 - t1593, t1438 * t1517 - t1521 * t1621, qJDD(3) * t1521 - t1474 * t1517, t1521 * t1350 - t1517 * t1334 - pkin(6) * (-t1439 * t1521 + t1479 * t1517), t1521 * t1351 - t1517 * t1333 - pkin(6) * (-t1441 * t1521 - t1517 * t1557), -pkin(2) * t1640 + t1517 * t1354 - pkin(6) * (-t1487 * t1517 + t1640), t1521 * t1294 - t1517 * t1269 - pkin(6) * (-t1373 * t1521 - t1449 * t1517), t1328 * t1517 + t1378 * t1521, t1279 * t1517 + t1314 * t1521, t1304 * t1517 + t1360 * t1521, t1327 * t1517 - t1376 * t1521, t1305 * t1517 + t1361 * t1521, t1366 * t1517 + t1397 * t1521, t1521 * t1156 - t1517 * t1125 - pkin(6) * (-t1281 * t1521 + t1348 * t1517), t1521 * t1157 - t1517 * t1131 - pkin(6) * (-t1290 * t1521 + t1357 * t1517), t1521 * t1123 - t1517 * t1112 - pkin(6) * (-t1249 * t1521 + t1315 * t1517), t1521 * t1115 - t1517 * t1087 - pkin(6) * (-t1212 * t1521 + t1232 * t1517), t1685, -t1782, t1764, t1724, -t1777, t1725, -t1517 * t1056 + t1521 * t1073 - t1748, -t1517 * t1060 + t1521 * t1076 + t1789, t1521 * t1045 - t1517 * t1041 - pkin(6) * (-t1120 * t1521 + t1141 * t1517), t1521 * t1035 - t1517 * t1034 - pkin(6) * (-t1077 * t1521 + t1082 * t1517), t1685, t1764, t1782, t1725, t1777, t1724, -t1517 * t1047 + t1521 * t1054 - t1748, t1521 * t1042 - t1517 * t1037 - pkin(6) * (-t1119 * t1521 + t1139 * t1517), -t1517 * t1044 + t1521 * t1053 - t1789, t1521 * t1032 - t1517 * t1031 - pkin(6) * (-t1061 * t1521 + t1071 * t1517); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1486, 0, t1485, 0, t1578, -t1579, t1595, pkin(6) * t1595, 0, -t1486, -t1485, 0, 0, 0, t1596, -t1578, t1579, pkin(6) * t1596 + (pkin(1) * t1521 + qJ(2) * t1517) * g(3), -t1451 * t1521 + t1594, -t1429 * t1521 + t1488 * t1517, -t1440 * t1521 + t1517 * t1620, -t1450 * t1521 - t1594, -t1438 * t1521 - t1517 * t1621, qJDD(3) * t1517 + t1474 * t1521, t1517 * t1350 + t1521 * t1334 + pkin(6) * (t1439 * t1517 + t1479 * t1521), t1517 * t1351 + t1521 * t1333 + pkin(6) * (t1441 * t1517 - t1521 * t1557), -pkin(2) * t1641 - t1521 * t1354 + pkin(6) * (-t1487 * t1521 - t1641), t1517 * t1294 + t1521 * t1269 + pkin(6) * (t1373 * t1517 - t1449 * t1521), -t1328 * t1521 + t1378 * t1517, -t1279 * t1521 + t1314 * t1517, -t1304 * t1521 + t1360 * t1517, -t1327 * t1521 - t1376 * t1517, -t1305 * t1521 + t1361 * t1517, -t1366 * t1521 + t1397 * t1517, t1517 * t1156 + t1521 * t1125 + pkin(6) * (t1281 * t1517 + t1348 * t1521), t1517 * t1157 + t1521 * t1131 + pkin(6) * (t1290 * t1517 + t1357 * t1521), t1517 * t1123 + t1521 * t1112 + pkin(6) * (t1249 * t1517 + t1315 * t1521), t1517 * t1115 + t1521 * t1087 + pkin(6) * (t1212 * t1517 + t1232 * t1521), t1686, t1781, t1763, t1722, t1776, t1723, t1521 * t1056 + t1517 * t1073 + t1749, t1521 * t1060 + t1517 * t1076 + t1788, t1517 * t1045 + t1521 * t1041 + pkin(6) * (t1120 * t1517 + t1141 * t1521), t1517 * t1035 + t1521 * t1034 + pkin(6) * (t1077 * t1517 + t1082 * t1521), t1686, t1763, -t1781, t1723, -t1776, t1722, t1521 * t1047 + t1517 * t1054 + t1749, t1517 * t1042 + t1521 * t1037 + pkin(6) * (t1119 * t1517 + t1139 * t1521), t1521 * t1044 + t1517 * t1053 - t1788, t1517 * t1032 + t1521 * t1031 + pkin(6) * (t1061 * t1517 + t1071 * t1521); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1493, t1494, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t1470, t1467, t1412, -t1444, t1430, t1445, t1437, t1442, 0, t1331, t1332, t1326, t1268, t1330, t1280, t1306, t1329, t1307, t1367, t1124, t1130, t1111, t1086, t1575, -t1128, t1751, t1698, -t1173, t1697, t1055, t1059, t1040, t1033, t1575, t1751, t1128, t1697, t1173, t1698, t1046, t1036, t1043, t1030; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1523, 0, 0, -g(3), -t1493, 0, 0, -qJDD(1), t1523, 0, 0, 0, -t1458, 0, g(3), qJ(2) * g(3), t1502, t1488, t1620, -t1502, -t1621, qJDD(3), t1350, t1351, -t1674, t1294, t1378, t1314, t1360, -t1376, t1361, t1397, t1156, t1157, t1123, t1115, t1198, -t1138, t1736, t1680, t1225, t1683, t1073, t1076, t1045, t1035, t1198, t1736, t1138, t1683, -t1225, t1680, t1054, t1042, t1053, t1032; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1523, 0, qJDD(1), 0, g(3), 0, -t1494, 0, 0, -t1523, -qJDD(1), 0, 0, 0, -t1457, -g(3), 0, pkin(1) * g(3), -t1451, -t1429, -t1440, -t1450, -t1438, t1474, t1334, t1333, -t1354, t1269, -t1328, -t1279, -t1304, -t1327, -t1305, -t1366, t1125, t1131, t1112, t1087, -t1560, t1126, -t1752, -t1701, t1169, -t1700, t1056, t1060, t1041, t1034, -t1560, -t1752, -t1126, -t1700, -t1169, -t1701, t1047, t1037, t1044, t1031; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1493, t1494, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t1470, t1467, t1412, -t1444, t1430, t1445, t1437, t1442, 0, t1331, t1332, t1326, t1268, t1330, t1280, t1306, t1329, t1307, t1367, t1124, t1130, t1111, t1086, t1575, -t1128, t1751, t1698, -t1173, t1697, t1055, t1059, t1040, t1033, t1575, t1751, t1128, t1697, t1173, t1698, t1046, t1036, t1043, t1030; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, -t1458, -t1457, 0, -t1444, t1430, t1445, t1437, t1442, 0, -pkin(7) * t1439 - t1627, -pkin(7) * t1441 - t1436, pkin(7) * t1484 - t1373, -pkin(7) * t1373, t1330, t1280, t1306, t1329, t1307, t1367, -pkin(7) * t1281 + t1600, -pkin(7) * t1290 + t1599, -pkin(7) * t1249 + t1625, -pkin(7) * t1212 - t1232 * t1669 + t1230, t1575, -t1128, t1751, t1698, -t1173, t1697, t1602 - t1759, t1601 - t1786, -pkin(7) * t1120 + t1605, -pkin(7) * t1077 + t1606, t1575, t1751, t1128, t1697, t1173, t1698, t1603 - t1759, -pkin(7) * t1119 + t1607, t1604 + t1786, -pkin(7) * t1061 + t1608; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1523, 0, 0, 0, t1458, 0, -g(3), 0, -t1502, -t1488, -t1620, t1502, t1621, -qJDD(3), t1545, t1564, t1674, -t1676, -t1378, -t1314, -t1360, t1376, -t1361, -t1397, t1566, t1565, t1548, t1584, -t1198, t1138, -t1736, -t1680, -t1225, -t1683, t1550, t1549, t1553, t1532, -t1198, -t1736, -t1138, -t1683, t1225, -t1680, t1551, t1554, t1552, t1555; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1523, qJDD(1), 0, 0, 0, t1457, g(3), 0, 0, t1451, t1429, t1440, t1450, t1438, -t1474, pkin(7) * t1443 - t1572, pkin(7) * t1446 - t1571, t1354, pkin(7) * t1374 + t1675, t1328, t1279, t1304, t1327, t1305, t1366, pkin(7) * t1282 - t1535, pkin(7) * t1291 - t1534, pkin(7) * t1250 + t1684, pkin(7) * t1213 - t1696, t1560, -t1126, t1752, t1701, -t1169, t1700, -t1537 + t1760, -t1536 + t1785, pkin(7) * t1122 - t1540, pkin(7) * t1078 - t1541, t1560, t1752, t1126, t1700, t1169, t1701, -t1538 + t1760, pkin(7) * t1121 - t1542, -t1539 - t1785, pkin(7) * t1062 - t1543; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1558, -t1479, t1492, t1613, t1497, -t1613, 0, -t1449, -t1431, 0, t1379, t1316, t1362, t1377, t1363, t1398, t1292, t1295, t1186, -pkin(8) * t1232, t1201, -t1142, t1735, t1681, -t1228, t1682, t1100, t1104, t1058, t1052, t1201, t1735, t1142, t1682, t1228, t1681, t1069, t1050, t1067, t1039; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1612, -t1557, t1499, t1480, t1491, -t1612, t1449, 0, -t1432, 0, -t1435, -t1434, -t1386, t1435, t1382, -t1473, t1242, t1244, -pkin(3) * t1315, -pkin(3) * t1232, -t1644, -t1692, -t1688, t1644, t1271, -t1469, t1114, t1118, t1109, t1065, -t1644, -t1688, t1692, -t1469, -t1271, t1644, t1098, t1090, t1091, t1048; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1502, t1488, t1620, -t1502, -t1621, qJDD(3), t1431, t1432, 0, 0, t1378, t1314, t1360, -t1376, t1361, t1397, t1610, t1611, t1585, t1623, t1198, -t1138, t1736, t1680, t1225, t1683, t1587, t1586, t1590, t1556, t1198, t1736, t1138, t1683, -t1225, t1680, t1588, t1591, t1589, t1592; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1420, t1383, t1691, t1455, t1452, -t1455, 0, -t1395, t1311, 0, t1260, -t1206, t1727, t1559, t1303, t1573, t1211, t1221, t1093, -pkin(9) * t1116, t1260, t1727, t1206, t1573, -t1303, t1559, t1113, t1085, t1108, t1075; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1456, t1385, t1453, t1527, t1414, -t1456, t1395, 0, t1312, 0, t1259, t1204, t1726, t1574, t1299, t1547, t1159, t1165, t1089, t1106, t1259, t1726, -t1204, t1547, -t1299, t1574, t1102, t1080, t1097, t1064; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1435, t1434, t1386, -t1435, -t1382, t1473, -t1311, -t1312, 0, 0, t1644, t1692, t1688, -t1644, -t1271, t1469, t1562, t1581, t1671, t1672, t1644, t1688, -t1692, t1469, t1271, -t1644, t1526, t1582, t1482 + t1531, t1583; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1320, -t1694, t1689, t1645, t1402, -t1645, 0, -t1289, t1175, 0, t1320, t1689, t1694, -t1645, -t1402, t1645, -qJ(6) * t1694, t1133, t1134, qJ(6) * t1155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1407, t1693, -t1403, -t1319, t1690, -t1407, t1289, 0, t1176, 0, t1407, -t1403, -t1693, -t1407, -t1690, -t1319, t1135, t1132, pkin(5) * t1693, pkin(5) * t1155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1644, t1692, t1688, -t1644, -t1271, t1469, -t1175, -t1176, 0, 0, t1644, t1688, -t1692, t1469, t1271, -t1644, t1528, t1624, t1482 + t1544, t1626; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1320, t1689, t1694, -t1645, -t1402, t1645, 0, t1149, t1155, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1644, t1688, -t1692, t1469, t1271, -t1644, -t1149, 0, t1147, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1407, t1403, t1693, t1407, t1690, t1319, -t1155, -t1147, 0, 0;];
m_new_reg  = t1;
