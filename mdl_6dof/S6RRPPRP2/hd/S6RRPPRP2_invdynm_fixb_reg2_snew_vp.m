% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6RRPPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 09:13
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6RRPPRP2_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP2_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP2_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:12:31
% EndTime: 2019-05-06 09:13:04
% DurationCPUTime: 35.56s
% Computational Cost: add. (104989->823), mult. (245766->922), div. (0->0), fcn. (168247->8), ass. (0->564)
t1617 = qJD(2) ^ 2;
t1608 = sin(pkin(9));
t1609 = cos(pkin(9));
t1612 = sin(qJ(2));
t1615 = cos(qJ(2));
t1570 = (t1608 * t1615 + t1609 * t1612) * qJD(1);
t1790 = t1570 ^ 2;
t1548 = t1790 + t1617;
t1758 = qJD(1) * t1615;
t1759 = qJD(1) * t1612;
t1568 = t1608 * t1759 - t1609 * t1758;
t1727 = t1570 * t1568;
t1838 = qJDD(2) + t1727;
t1850 = t1838 * t1608;
t1435 = t1548 * t1609 + t1850;
t1849 = t1838 * t1609;
t1441 = -t1548 * t1608 + t1849;
t1366 = t1435 * t1612 - t1441 * t1615;
t1613 = sin(qJ(1));
t1616 = cos(qJ(1));
t1598 = qJD(2) * t1758;
t1701 = t1612 * qJDD(1);
t1580 = t1598 + t1701;
t1695 = qJD(2) * t1759;
t1700 = t1615 * qJDD(1);
t1656 = t1695 - t1700;
t1514 = t1609 * t1580 - t1608 * t1656;
t1730 = t1568 * qJD(2);
t1800 = -t1730 + t1514;
t1899 = pkin(6) * (t1366 * t1616 + t1613 * t1800);
t1898 = pkin(6) * (t1366 * t1613 - t1616 * t1800);
t1359 = t1435 * t1615 + t1441 * t1612;
t1897 = pkin(7) * t1359;
t1888 = pkin(2) * t1435;
t1896 = pkin(1) * t1359 + t1888;
t1895 = pkin(1) * t1800 - pkin(7) * t1366;
t1799 = -t1790 + t1617;
t1497 = t1727 - qJDD(2);
t1868 = t1497 * t1608;
t1431 = -t1609 * t1799 + t1868;
t1867 = t1497 * t1609;
t1437 = t1608 * t1799 + t1867;
t1362 = t1431 * t1612 - t1437 * t1615;
t1803 = t1514 + t1730;
t1894 = t1362 * t1613 - t1616 * t1803;
t1893 = t1362 * t1616 + t1613 * t1803;
t1565 = t1568 ^ 2;
t1547 = t1617 - t1565;
t1433 = -t1547 * t1608 + t1849;
t1440 = t1547 * t1609 + t1850;
t1364 = t1433 * t1612 + t1440 * t1615;
t1684 = t1580 * t1608 + t1609 * t1656;
t1757 = qJD(2) * t1570;
t1840 = -t1757 + t1684;
t1892 = t1364 * t1613 - t1616 * t1840;
t1891 = t1364 * t1616 + t1613 * t1840;
t1496 = -t1617 - t1565;
t1413 = t1496 * t1608 - t1867;
t1416 = -t1496 * t1609 - t1868;
t1330 = t1413 * t1612 + t1416 * t1615;
t1839 = t1757 + t1684;
t1887 = pkin(6) * (t1330 * t1616 - t1613 * t1839);
t1886 = pkin(6) * (t1330 * t1613 + t1616 * t1839);
t1885 = qJ(3) * t1435;
t1884 = qJ(3) * t1441;
t1798 = t1790 - t1565;
t1847 = -t1608 * t1839 + t1609 * t1800;
t1740 = t1800 * t1608;
t1741 = t1839 * t1609;
t1863 = -t1740 - t1741;
t1872 = -t1612 * t1847 + t1615 * t1863;
t1883 = t1613 * t1872 - t1616 * t1798;
t1882 = t1613 * t1798 + t1616 * t1872;
t1355 = t1431 * t1615 + t1437 * t1612;
t1357 = t1433 * t1615 - t1440 * t1612;
t1327 = t1413 * t1615 - t1416 * t1612;
t1881 = pkin(1) * t1327;
t1846 = t1608 * t1803 - t1609 * t1840;
t1848 = -t1608 * t1840 - t1609 * t1803;
t1860 = t1612 * t1846 + t1615 * t1848;
t1880 = pkin(1) * t1860;
t1879 = pkin(7) * t1327;
t1878 = pkin(7) * t1860;
t1877 = pkin(1) * t1839 + pkin(7) * t1330;
t1801 = -t1565 - t1790;
t1861 = -t1612 * t1848 + t1615 * t1846;
t1876 = -pkin(1) * t1801 + pkin(7) * t1861;
t1873 = t1612 * t1863 + t1615 * t1847;
t1871 = pkin(6) * (t1613 * t1861 - t1616 * t1801);
t1870 = pkin(6) * (t1613 * t1801 + t1616 * t1861);
t1784 = pkin(2) * t1848;
t1869 = qJ(3) * t1848;
t1755 = qJD(3) * t1570;
t1858 = pkin(2) * t1413;
t1864 = -0.2e1 * t1755 + t1858;
t1862 = -pkin(2) * t1801 + qJ(3) * t1846;
t1611 = sin(qJ(5));
t1614 = cos(qJ(5));
t1525 = qJD(2) * t1611 - t1614 * t1568;
t1527 = t1614 * qJD(2) + t1611 * t1568;
t1467 = t1527 * t1525;
t1510 = qJDD(5) + t1514;
t1819 = t1467 - t1510;
t1857 = pkin(5) * t1819;
t1856 = qJ(3) * t1413;
t1855 = qJ(3) * t1416;
t1812 = qJ(4) * t1800;
t1523 = t1525 ^ 2;
t1563 = qJD(5) + t1570;
t1561 = t1563 ^ 2;
t1423 = -t1561 - t1523;
t1744 = t1819 * t1614;
t1336 = t1423 * t1611 - t1744;
t1333 = pkin(4) * t1336;
t1427 = -t1525 * qJD(5) + t1614 * qJDD(2) + t1611 * t1684;
t1542 = pkin(4) * t1570 - qJD(2) * pkin(8);
t1606 = t1615 ^ 2;
t1618 = qJD(1) ^ 2;
t1589 = t1613 * g(1) - t1616 * g(2);
t1659 = qJDD(1) * pkin(1) + t1589;
t1464 = (qJ(3) * t1606 + pkin(7)) * t1618 - pkin(2) * t1656 - qJDD(3) - (qJD(2) * pkin(2) - qJ(3) * t1759) * t1759 + t1659;
t1837 = 2 * qJD(4);
t1820 = pkin(3) * t1757 - t1570 * t1837 - t1464;
t1620 = -t1812 + t1820;
t1665 = t1684 * pkin(3);
t1289 = -pkin(4) * t1565 + pkin(8) * t1684 - t1570 * t1542 + t1620 + t1665;
t1714 = t1612 * t1618;
t1590 = g(1) * t1616 + g(2) * t1613;
t1574 = -pkin(1) * t1618 + qJDD(1) * pkin(7) - t1590;
t1715 = t1612 * t1574;
t1760 = qJD(1) * qJD(2);
t1455 = qJDD(2) * pkin(2) - t1580 * qJ(3) - t1715 + (pkin(2) * t1714 + qJ(3) * t1760 - g(3)) * t1615;
t1539 = -t1612 * g(3) + t1615 * t1574;
t1601 = t1606 * t1618;
t1595 = -t1601 - t1617;
t1456 = pkin(2) * t1595 + qJ(3) * t1700 + t1539;
t1685 = -t1609 * t1455 + t1608 * t1456;
t1646 = -qJDD(2) * pkin(3) - t1617 * qJ(4) + qJDD(4) + t1685;
t1495 = pkin(3) * t1568 - qJ(4) * t1570;
t1713 = 0.2e1 * qJD(3) + t1495;
t1621 = -qJDD(2) * pkin(8) + t1803 * pkin(4) + (pkin(8) * t1568 + t1713) * t1570 + t1646;
t1222 = t1611 * t1289 - t1614 * t1621;
t1489 = t1563 * t1525;
t1631 = qJ(6) * t1489 + 0.2e1 * qJD(6) * t1527 + t1222 + t1857;
t1193 = qJ(6) * t1427 + t1631;
t1628 = -t1193 - t1857;
t1841 = -t1333 - t1628;
t1716 = t1611 * t1819;
t1599 = t1613 * qJDD(2);
t1725 = t1570 * t1609;
t1729 = t1568 * t1608;
t1652 = (-t1725 - t1729) * qJD(2);
t1726 = t1570 * t1608;
t1728 = t1568 * t1609;
t1653 = (t1726 - t1728) * qJD(2);
t1793 = -t1612 * t1652 + t1615 * t1653;
t1818 = t1616 * t1793 + t1599;
t1696 = t1616 * t1727;
t1672 = -qJD(2) * t1726 + t1609 * t1514;
t1674 = qJD(2) * t1725 + t1608 * t1514;
t1792 = -t1612 * t1674 + t1615 * t1672;
t1817 = t1613 * t1792 - t1696;
t1632 = qJD(2) * t1728 + t1608 * t1684;
t1673 = qJD(2) * t1729 - t1609 * t1684;
t1796 = -t1612 * t1673 + t1615 * t1632;
t1816 = t1613 * t1796 + t1696;
t1702 = qJDD(2) * t1616;
t1815 = t1613 * t1793 - t1702;
t1697 = t1613 * t1727;
t1814 = t1616 * t1792 + t1697;
t1813 = t1616 * t1796 - t1697;
t1223 = t1614 * t1289 + t1611 * t1621;
t1169 = -t1614 * t1222 + t1223 * t1611;
t1490 = t1563 * t1527;
t1682 = -t1611 * qJDD(2) + t1614 * t1684;
t1655 = qJD(5) * t1527 - t1682;
t1387 = t1490 + t1655;
t1811 = -pkin(3) * t1336 + qJ(4) * t1387;
t1524 = t1527 ^ 2;
t1443 = -t1524 - t1561;
t1411 = t1467 + t1510;
t1717 = t1611 * t1411;
t1343 = t1443 * t1614 - t1717;
t1807 = -t1489 + t1427;
t1810 = -pkin(3) * t1343 + qJ(4) * t1807;
t1389 = (-qJD(5) + t1563) * t1527 + t1682;
t1391 = t1489 + t1427;
t1302 = t1389 * t1611 - t1391 * t1614;
t1420 = -t1523 - t1524;
t1809 = -pkin(3) * t1302 + qJ(4) * t1420;
t1806 = pkin(3) * t1548 + qJ(4) * t1838;
t1805 = pkin(3) * t1497 - qJ(4) * t1496;
t1804 = -t1684 * pkin(4) - t1565 * pkin(8);
t1479 = pkin(5) * t1563 - qJ(6) * t1527;
t1797 = pkin(5) * t1655 - t1523 * qJ(6) + t1527 * t1479 + qJDD(6);
t1795 = t1612 * t1632 + t1615 * t1673;
t1794 = t1612 * t1653 + t1615 * t1652;
t1791 = t1612 * t1672 + t1615 * t1674;
t1789 = -pkin(3) - pkin(8);
t1267 = -t1302 * t1609 + t1420 * t1608;
t1268 = t1302 * t1608 + t1420 * t1609;
t1208 = t1267 * t1615 + t1268 * t1612;
t1788 = pkin(1) * t1208;
t1277 = -t1336 * t1609 + t1387 * t1608;
t1278 = t1336 * t1608 + t1387 * t1609;
t1218 = t1277 * t1615 + t1278 * t1612;
t1787 = pkin(1) * t1218;
t1286 = -t1343 * t1609 + t1608 * t1807;
t1287 = t1343 * t1608 + t1609 * t1807;
t1225 = t1286 * t1615 + t1287 * t1612;
t1786 = pkin(1) * t1225;
t1346 = t1685 + 0.2e1 * t1755;
t1662 = t1608 * t1455 + t1609 * t1456;
t1756 = qJD(3) * t1568;
t1347 = t1662 - 0.2e1 * t1756;
t1273 = -t1346 * t1609 + t1347 * t1608;
t1785 = pkin(2) * t1273;
t1304 = t1389 * t1614 + t1611 * t1391;
t1779 = pkin(3) * t1304;
t1337 = t1423 * t1614 + t1716;
t1778 = pkin(3) * t1337;
t1745 = t1411 * t1614;
t1344 = -t1611 * t1443 - t1745;
t1777 = pkin(3) * t1344;
t1776 = pkin(3) * t1609;
t1775 = pkin(4) * t1169;
t1634 = -t1617 * pkin(3) + qJDD(2) * qJ(4) + t1662;
t1626 = -t1568 * t1713 + t1634;
t1281 = (t1837 + t1542) * qJD(2) + t1626 + t1804;
t1774 = pkin(4) * t1281;
t1209 = -t1267 * t1612 + t1268 * t1615;
t1773 = pkin(6) * (t1209 * t1613 - t1304 * t1616);
t1219 = -t1277 * t1612 + t1278 * t1615;
t1772 = pkin(6) * (t1219 * t1613 - t1337 * t1616);
t1226 = -t1286 * t1612 + t1287 * t1615;
t1771 = pkin(6) * (t1226 * t1613 - t1344 * t1616);
t1770 = pkin(7) * t1208;
t1769 = pkin(7) * t1218;
t1768 = pkin(7) * t1225;
t1767 = pkin(8) * t1169;
t1766 = qJ(3) * t1267;
t1765 = qJ(3) * t1277;
t1764 = qJ(3) * t1286;
t1763 = qJ(4) * t1337;
t1762 = qJ(4) * t1344;
t1754 = qJD(4) * qJD(2);
t1752 = qJD(6) * t1525;
t1750 = t1193 * t1614;
t1748 = t1273 * t1612;
t1747 = t1273 * t1615;
t1746 = t1281 * t1611;
t1743 = t1464 * t1608;
t1742 = t1464 * t1609;
t1732 = t1563 * t1611;
t1731 = t1563 * t1614;
t1573 = t1618 * pkin(7) + t1659;
t1724 = t1573 * t1612;
t1723 = t1573 * t1615;
t1581 = -0.2e1 * t1695 + t1700;
t1528 = t1581 * t1615;
t1596 = t1615 * t1714;
t1587 = qJDD(2) + t1596;
t1722 = t1587 * t1612;
t1588 = qJDD(2) - t1596;
t1721 = t1588 * t1612;
t1720 = t1588 * t1615;
t1605 = t1612 ^ 2;
t1719 = t1605 * t1618;
t1718 = t1611 * t1193;
t1279 = t1614 * t1281;
t1712 = -pkin(1) * t1304 + pkin(7) * t1209;
t1711 = -pkin(1) * t1337 + pkin(7) * t1219;
t1710 = -pkin(1) * t1344 + pkin(7) * t1226;
t1338 = pkin(8) * t1343;
t1709 = t1279 - t1338;
t1708 = -pkin(4) * t1420 + pkin(8) * t1304;
t1316 = t1626 + 0.2e1 * t1754;
t1319 = t1713 * t1570 + t1646;
t1707 = -pkin(3) * t1319 + qJ(4) * t1316;
t1706 = pkin(4) * t1387 - pkin(8) * t1337;
t1705 = pkin(4) * t1807 - pkin(8) * t1344;
t1704 = -pkin(3) * t1803 - qJ(4) * t1840;
t1703 = t1605 + t1606;
t1699 = t1608 * t1467;
t1698 = t1609 * t1467;
t1340 = pkin(4) * t1343;
t1694 = -t1340 + t1223;
t1693 = qJ(4) * t1608 + pkin(2);
t1691 = -pkin(2) * t1304 + qJ(3) * t1268;
t1690 = -pkin(2) * t1337 + qJ(3) * t1278;
t1689 = -pkin(2) * t1344 + qJ(3) * t1287;
t1517 = -0.2e1 * t1752;
t1663 = t1523 * pkin(5) + qJ(6) * t1655 + t1563 * t1479 - t1223;
t1195 = t1517 - t1663;
t1156 = t1195 * t1611 - t1750;
t1192 = pkin(5) * t1193;
t1688 = -pkin(4) * t1156 + t1192;
t1298 = pkin(4) * t1302;
t1241 = -qJ(4) * t1304 + t1298;
t1331 = pkin(8) * t1336;
t1686 = -t1331 + t1746;
t1274 = t1346 * t1608 + t1609 * t1347;
t1538 = t1615 * g(3) + t1715;
t1459 = t1538 * t1612 + t1615 * t1539;
t1683 = -t1589 * t1613 - t1616 * t1590;
t1681 = t1613 * t1596;
t1680 = t1616 * t1596;
t1679 = t1709 + t1810;
t1249 = t1316 * t1608 - t1319 * t1609;
t1678 = pkin(2) * t1249 + t1707;
t1677 = t1704 + t1784;
t1676 = -pkin(3) * t1169 + qJ(4) * t1281 - t1767;
t1584 = qJDD(1) * t1616 - t1613 * t1618;
t1675 = -pkin(6) * t1584 - g(3) * t1613;
t1178 = -pkin(5) * t1420 + qJ(6) * t1389 + t1195;
t1180 = (t1391 + t1427) * qJ(6) + t1631;
t1296 = pkin(8) * t1302;
t1671 = -t1178 * t1611 + t1614 * t1180 - t1296;
t1670 = -t1296 - t1169;
t1244 = t1281 + t1797;
t1234 = -qJ(6) * t1443 + t1244;
t1318 = -pkin(5) * t1807 - qJ(6) * t1411;
t1669 = t1614 * t1234 - t1318 * t1611 - t1338;
t1668 = t1222 - t1333;
t1667 = -t1706 - t1279;
t1666 = -t1705 + t1746;
t1285 = pkin(2) * t1286;
t1664 = t1285 + t1679;
t1170 = t1611 * t1222 + t1223 * t1614;
t1458 = t1538 * t1615 - t1539 * t1612;
t1661 = t1589 * t1616 - t1590 * t1613;
t1658 = t1686 + t1811;
t1657 = pkin(5) * t1443 + t1663;
t1160 = -t1169 * t1609 + t1281 * t1608;
t1654 = pkin(2) * t1160 + t1676;
t1651 = t1671 + t1809;
t1650 = t1670 + t1809;
t1649 = t1669 + t1810;
t1276 = pkin(2) * t1277;
t1647 = t1276 + t1658;
t1554 = 0.2e1 * t1756;
t1624 = t1568 * t1495 + t1554 - t1634 - 0.2e1 * t1754;
t1210 = -pkin(5) * t1387 + qJ(6) * t1423 - qJD(2) * t1542 + t1624 - t1797 - t1804;
t1645 = qJ(6) * t1744 - t1210 * t1611 - t1331;
t1644 = 0.2e1 * t1752 + t1657;
t1643 = t1178 * t1614 + t1180 * t1611 + t1708;
t1642 = t1170 + t1708;
t1641 = t1234 * t1611 + t1318 * t1614 - t1705;
t1266 = pkin(2) * t1267;
t1640 = t1266 + t1651;
t1639 = t1266 + t1650;
t1638 = t1285 + t1649;
t1637 = qJ(6) * t1716 + t1210 * t1614 - t1706;
t1165 = -pkin(5) * t1244 + qJ(6) * t1195;
t1636 = pkin(4) * t1244 - qJ(6) * t1718 - t1614 * t1165;
t1635 = -pkin(8) * t1156 + qJ(6) * t1750 - t1165 * t1611;
t1633 = t1645 + t1811;
t1630 = t1276 + t1633;
t1629 = -pkin(3) * t1156 + qJ(4) * t1244 + t1635;
t1146 = -t1156 * t1609 + t1244 * t1608;
t1627 = pkin(2) * t1146 + t1629;
t1625 = t1319 + t1805;
t1622 = t1316 + t1806;
t1619 = -t1665 - t1820;
t1594 = t1601 - t1617;
t1593 = -t1617 - t1719;
t1592 = t1617 - t1719;
t1586 = -t1601 + t1719;
t1585 = t1601 + t1719;
t1583 = qJDD(1) * t1613 + t1616 * t1618;
t1582 = t1703 * qJDD(1);
t1579 = 0.2e1 * t1598 + t1701;
t1577 = t1615 * t1587;
t1576 = t1703 * t1760;
t1564 = -pkin(6) * t1583 + g(3) * t1616;
t1541 = t1580 * t1615 - t1605 * t1760;
t1540 = -t1606 * t1760 + t1612 * t1656;
t1537 = -t1593 * t1612 - t1720;
t1536 = -t1592 * t1612 + t1577;
t1535 = t1595 * t1615 - t1722;
t1534 = t1594 * t1615 - t1721;
t1533 = t1593 * t1615 - t1721;
t1532 = t1592 * t1615 + t1722;
t1531 = t1595 * t1612 + t1577;
t1530 = t1594 * t1612 + t1720;
t1529 = (t1580 + t1598) * t1612;
t1516 = -t1579 * t1612 + t1528;
t1515 = t1579 * t1615 + t1581 * t1612;
t1483 = -t1524 + t1561;
t1482 = t1523 - t1561;
t1481 = -pkin(7) * t1533 - t1723;
t1480 = -pkin(7) * t1531 - t1724;
t1477 = -pkin(1) * t1533 + t1539;
t1476 = -pkin(1) * t1531 + t1538;
t1460 = t1524 - t1523;
t1446 = pkin(1) * t1581 + pkin(7) * t1535 + t1723;
t1445 = -pkin(1) * t1579 + pkin(7) * t1537 - t1724;
t1421 = pkin(1) * t1573 + pkin(7) * t1459;
t1419 = pkin(1) * t1585 + pkin(7) * t1582 + t1459;
t1406 = (-t1525 * t1614 + t1527 * t1611) * t1563;
t1405 = (-t1525 * t1611 - t1527 * t1614) * t1563;
t1388 = -t1490 + t1655;
t1384 = pkin(5) * t1391;
t1379 = t1427 * t1614 - t1527 * t1732;
t1378 = t1427 * t1611 + t1527 * t1731;
t1377 = -t1525 * t1731 - t1611 * t1655;
t1376 = -t1525 * t1732 + t1614 * t1655;
t1375 = -t1742 + t1885;
t1354 = t1405 * t1608 + t1510 * t1609;
t1353 = -t1405 * t1609 + t1510 * t1608;
t1352 = t1482 * t1614 - t1717;
t1351 = -t1611 * t1483 - t1744;
t1350 = t1482 * t1611 + t1745;
t1349 = t1483 * t1614 - t1716;
t1348 = -t1743 - t1856;
t1326 = -pkin(2) * t1800 - t1743 - t1884;
t1325 = t1619 + t1812;
t1324 = t1378 * t1608 + t1698;
t1323 = -t1376 * t1608 - t1698;
t1322 = -t1378 * t1609 + t1699;
t1321 = t1376 * t1609 - t1699;
t1320 = -pkin(2) * t1839 + t1742 - t1855;
t1306 = (t1839 + t1684) * pkin(3) + t1620;
t1305 = t1619 + 0.2e1 * t1812;
t1303 = -t1387 * t1614 - t1611 * t1807;
t1301 = -t1387 * t1611 + t1614 * t1807;
t1295 = -qJ(4) * t1801 + t1319;
t1294 = -pkin(3) * t1801 + t1316;
t1293 = t1349 * t1608 + t1391 * t1609;
t1292 = t1350 * t1608 - t1388 * t1609;
t1291 = -t1349 * t1609 + t1391 * t1608;
t1290 = -t1350 * t1609 - t1388 * t1608;
t1283 = -t1353 * t1612 + t1354 * t1615;
t1282 = t1353 * t1615 + t1354 * t1612;
t1272 = t1301 * t1608 + t1460 * t1609;
t1271 = -t1301 * t1609 + t1460 * t1608;
t1264 = -t1784 - t1880;
t1263 = t1347 + t1896;
t1262 = pkin(2) * t1464 + qJ(3) * t1274;
t1261 = -pkin(3) * t1740 + t1305 * t1609 - t1885;
t1260 = -t1322 * t1612 + t1324 * t1615;
t1259 = -t1321 * t1612 + t1323 * t1615;
t1258 = t1322 * t1615 + t1324 * t1612;
t1257 = t1321 * t1615 + t1323 * t1612;
t1256 = qJ(4) * t1741 - t1306 * t1608 + t1856;
t1255 = t1283 * t1616 + t1406 * t1613;
t1254 = t1283 * t1613 - t1406 * t1616;
t1253 = t1346 - t1858 - t1881;
t1252 = -t1273 - t1869;
t1251 = t1884 + t1608 * t1305 + (pkin(2) + t1776) * t1800;
t1250 = t1316 * t1609 + t1319 * t1608;
t1248 = t1609 * t1306 + t1693 * t1839 + t1855;
t1247 = -t1326 * t1612 + t1375 * t1615 + t1897;
t1246 = -t1677 - t1880;
t1245 = t1274 + t1862;
t1242 = t1624 - t1806 - t1896;
t1240 = t1260 * t1616 + t1379 * t1613;
t1239 = t1259 * t1616 - t1377 * t1613;
t1238 = t1260 * t1613 - t1379 * t1616;
t1237 = t1259 * t1613 + t1377 * t1616;
t1236 = t1326 * t1615 + t1375 * t1612 - t1895;
t1235 = -t1320 * t1612 + t1348 * t1615 - t1879;
t1232 = -t1291 * t1612 + t1293 * t1615;
t1231 = -t1290 * t1612 + t1292 * t1615;
t1230 = t1291 * t1615 + t1293 * t1612;
t1229 = t1290 * t1615 + t1292 * t1612;
t1228 = -t1570 * t1495 - t1646 - t1805 + t1864 + t1881;
t1227 = t1320 * t1615 + t1348 * t1612 - t1877;
t1216 = t1274 * t1615 - t1748;
t1215 = t1274 * t1612 + t1747;
t1214 = t1241 - t1384;
t1213 = -t1271 * t1612 + t1272 * t1615;
t1212 = t1271 * t1615 + t1272 * t1612;
t1211 = -t1294 * t1608 + t1295 * t1609 - t1869;
t1206 = t1232 * t1616 + t1351 * t1613;
t1205 = t1231 * t1616 + t1352 * t1613;
t1204 = t1232 * t1613 - t1351 * t1616;
t1203 = t1231 * t1613 - t1352 * t1616;
t1202 = t1294 * t1609 + t1295 * t1608 + t1862;
t1201 = -t1666 - t1777;
t1199 = pkin(6) * (t1226 * t1616 + t1344 * t1613);
t1198 = -t1667 - t1778;
t1196 = pkin(6) * (t1219 * t1616 + t1337 * t1613);
t1194 = -qJ(3) * t1249 + (-pkin(3) * t1608 + qJ(4) * t1609) * t1325;
t1190 = -t1249 * t1612 + t1250 * t1615;
t1189 = t1249 * t1615 + t1250 * t1612;
t1188 = t1213 * t1616 + t1303 * t1613;
t1187 = t1213 * t1613 - t1303 * t1616;
t1186 = -pkin(1) * t1215 - t1785;
t1185 = -t1694 - t1762;
t1183 = pkin(6) * (t1209 * t1616 + t1304 * t1613);
t1182 = -t1668 - t1763;
t1181 = -t1251 * t1612 + t1261 * t1615 - t1897;
t1177 = qJ(3) * t1250 + (t1693 + t1776) * t1325;
t1176 = t1251 * t1615 + t1261 * t1612 + t1895;
t1175 = -t1248 * t1612 + t1256 * t1615 + t1879;
t1174 = t1248 * t1615 + t1256 * t1612 + t1877;
t1173 = -t1245 * t1612 + t1252 * t1615 - t1878;
t1172 = t1245 * t1615 + t1252 * t1612 + t1876;
t1171 = t1340 + t1644 - t1762;
t1167 = -t1641 - t1777;
t1166 = -t1637 - t1778;
t1164 = -t1763 - t1841;
t1163 = -pkin(7) * t1215 - qJ(3) * t1747 - t1262 * t1612;
t1162 = pkin(1) * t1464 + pkin(7) * t1216 - qJ(3) * t1748 + t1262 * t1615;
t1161 = t1169 * t1608 + t1281 * t1609;
t1159 = -t1202 * t1612 + t1211 * t1615 - t1878;
t1158 = -t1664 - t1786;
t1157 = t1195 * t1614 + t1718;
t1154 = -t1647 - t1787;
t1153 = t1202 * t1615 + t1211 * t1612 + t1876;
t1152 = -pkin(1) * t1189 - t1678;
t1151 = -t1642 - t1779;
t1150 = t1185 * t1609 - t1201 * t1608 - t1764;
t1149 = t1182 * t1609 - t1198 * t1608 - t1765;
t1148 = -t1638 - t1786;
t1147 = t1156 * t1608 + t1244 * t1609;
t1145 = t1185 * t1608 + t1201 * t1609 + t1689;
t1144 = -t1630 - t1787;
t1143 = t1182 * t1608 + t1198 * t1609 + t1690;
t1142 = -qJ(4) * t1170 + t1775;
t1141 = -t1643 - t1779;
t1140 = -t1151 * t1608 + t1241 * t1609 - t1766;
t1139 = t1170 * t1789 + t1774;
t1138 = -t1167 * t1608 + t1171 * t1609 - t1764;
t1137 = -pkin(7) * t1189 - t1177 * t1612 + t1194 * t1615;
t1136 = -t1639 - t1788;
t1135 = t1164 * t1609 - t1166 * t1608 - t1765;
t1134 = t1151 * t1609 + t1241 * t1608 + t1691;
t1133 = t1167 * t1609 + t1171 * t1608 + t1689;
t1132 = pkin(1) * t1325 + pkin(7) * t1190 + t1177 * t1615 + t1194 * t1612;
t1131 = t1164 * t1608 + t1166 * t1609 + t1690;
t1130 = -t1160 * t1612 + t1161 * t1615;
t1129 = t1160 * t1615 + t1161 * t1612;
t1128 = -t1141 * t1608 + t1214 * t1609 - t1766;
t1127 = -t1640 - t1788;
t1126 = t1141 * t1609 + t1214 * t1608 + t1691;
t1125 = -qJ(4) * t1157 - t1688;
t1124 = -t1146 * t1612 + t1147 * t1615;
t1123 = t1146 * t1615 + t1147 * t1612;
t1122 = -t1145 * t1612 + t1150 * t1615 - t1768;
t1121 = -t1143 * t1612 + t1149 * t1615 - t1769;
t1120 = t1145 * t1615 + t1150 * t1612 + t1710;
t1119 = t1143 * t1615 + t1149 * t1612 + t1711;
t1118 = t1157 * t1789 + t1636;
t1117 = -t1134 * t1612 + t1140 * t1615 - t1770;
t1116 = -t1133 * t1612 + t1138 * t1615 - t1768;
t1115 = t1133 * t1615 + t1138 * t1612 + t1710;
t1114 = t1134 * t1615 + t1140 * t1612 + t1712;
t1113 = -t1131 * t1612 + t1135 * t1615 - t1769;
t1112 = t1131 * t1615 + t1135 * t1612 + t1711;
t1111 = -qJ(3) * t1160 - t1139 * t1608 + t1142 * t1609;
t1110 = -pkin(1) * t1129 - t1654;
t1109 = -pkin(2) * t1170 + qJ(3) * t1161 + t1139 * t1609 + t1142 * t1608;
t1108 = -t1126 * t1612 + t1128 * t1615 - t1770;
t1107 = t1126 * t1615 + t1128 * t1612 + t1712;
t1106 = -qJ(3) * t1146 - t1118 * t1608 + t1125 * t1609;
t1105 = -pkin(1) * t1123 - t1627;
t1104 = -pkin(2) * t1157 + qJ(3) * t1147 + t1118 * t1609 + t1125 * t1608;
t1103 = -pkin(7) * t1129 - t1109 * t1612 + t1111 * t1615;
t1102 = -pkin(1) * t1170 + pkin(7) * t1130 + t1109 * t1615 + t1111 * t1612;
t1101 = -pkin(7) * t1123 - t1104 * t1612 + t1106 * t1615;
t1100 = -pkin(1) * t1157 + pkin(7) * t1124 + t1104 * t1615 + t1106 * t1612;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1584, 0, -t1583, 0, t1675, -t1564, -t1661, -pkin(6) * t1661, t1541 * t1616 - t1681, t1516 * t1616 + t1586 * t1613, t1536 * t1616 + t1613 * t1701, t1540 * t1616 + t1681, t1534 * t1616 + t1613 * t1700, t1576 * t1616 + t1599, t1616 * t1480 - t1613 * t1476 - pkin(6) * (t1535 * t1613 + t1581 * t1616), t1616 * t1481 - t1613 * t1477 - pkin(6) * (t1537 * t1613 - t1579 * t1616), t1616 * t1458 - pkin(6) * (t1582 * t1613 + t1585 * t1616), -pkin(6) * (t1459 * t1613 + t1573 * t1616) - (pkin(1) * t1613 - pkin(7) * t1616) * t1458, t1814, t1882, t1893, t1813, -t1891, t1818, t1616 * t1235 - t1613 * t1253 + t1886, t1616 * t1247 - t1613 * t1263 - t1898, t1616 * t1173 - t1613 * t1264 - t1871, t1616 * t1163 - t1613 * t1186 - pkin(6) * (t1216 * t1613 + t1464 * t1616), t1818, -t1893, t1891, t1814, t1882, t1813, t1616 * t1159 - t1613 * t1246 - t1871, t1616 * t1175 - t1613 * t1228 - t1886, t1616 * t1181 - t1613 * t1242 + t1898, t1616 * t1137 - t1613 * t1152 - pkin(6) * (t1190 * t1613 + t1325 * t1616), t1240, t1188, t1206, t1239, t1205, t1255, t1121 * t1616 - t1154 * t1613 - t1772, t1122 * t1616 - t1158 * t1613 - t1771, t1117 * t1616 - t1136 * t1613 - t1773, t1616 * t1103 - t1613 * t1110 - pkin(6) * (t1130 * t1613 - t1170 * t1616), t1240, t1188, t1206, t1239, t1205, t1255, t1113 * t1616 - t1144 * t1613 - t1772, t1116 * t1616 - t1148 * t1613 - t1771, t1108 * t1616 - t1127 * t1613 - t1773, t1616 * t1101 - t1613 * t1105 - pkin(6) * (t1124 * t1613 - t1157 * t1616); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1583, 0, t1584, 0, t1564, t1675, t1683, pkin(6) * t1683, t1541 * t1613 + t1680, t1516 * t1613 - t1586 * t1616, t1536 * t1613 - t1616 * t1701, t1540 * t1613 - t1680, t1534 * t1613 - t1616 * t1700, t1576 * t1613 - t1702, t1613 * t1480 + t1616 * t1476 + pkin(6) * (t1535 * t1616 - t1581 * t1613), t1613 * t1481 + t1616 * t1477 + pkin(6) * (t1537 * t1616 + t1579 * t1613), t1613 * t1458 + pkin(6) * (t1582 * t1616 - t1585 * t1613), pkin(6) * (t1459 * t1616 - t1573 * t1613) - (-pkin(1) * t1616 - pkin(7) * t1613) * t1458, t1817, t1883, t1894, t1816, -t1892, t1815, t1613 * t1235 + t1616 * t1253 - t1887, t1613 * t1247 + t1616 * t1263 + t1899, t1613 * t1173 + t1616 * t1264 + t1870, t1613 * t1163 + t1616 * t1186 + pkin(6) * (t1216 * t1616 - t1464 * t1613), t1815, -t1894, t1892, t1817, t1883, t1816, t1613 * t1159 + t1616 * t1246 + t1870, t1613 * t1175 + t1616 * t1228 + t1887, t1613 * t1181 + t1616 * t1242 - t1899, t1613 * t1137 + t1616 * t1152 + pkin(6) * (t1190 * t1616 - t1325 * t1613), t1238, t1187, t1204, t1237, t1203, t1254, t1121 * t1613 + t1154 * t1616 + t1196, t1122 * t1613 + t1158 * t1616 + t1199, t1117 * t1613 + t1136 * t1616 + t1183, t1613 * t1103 + t1616 * t1110 + pkin(6) * (t1130 * t1616 + t1170 * t1613), t1238, t1187, t1204, t1237, t1203, t1254, t1113 * t1613 + t1144 * t1616 + t1196, t1116 * t1613 + t1148 * t1616 + t1199, t1108 * t1613 + t1127 * t1616 + t1183, t1613 * t1101 + t1616 * t1105 + pkin(6) * (t1124 * t1616 + t1157 * t1613); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1589, t1590, 0, 0, t1529, t1515, t1532, t1528, t1530, 0, t1446, t1445, t1419, t1421, t1791, t1873, -t1355, t1795, t1357, t1794, t1227, t1236, t1172, t1162, t1794, t1355, -t1357, t1791, t1873, t1795, t1153, t1174, t1176, t1132, t1258, t1212, t1230, t1257, t1229, t1282, t1119, t1120, t1114, t1102, t1258, t1212, t1230, t1257, t1229, t1282, t1112, t1115, t1107, t1100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1618, 0, 0, -g(3), -t1589, 0, t1541, t1516, t1536, t1540, t1534, t1576, t1480, t1481, t1458, pkin(7) * t1458, t1792, t1872, t1362, t1796, -t1364, t1793, t1235, t1247, t1173, t1163, t1793, -t1362, t1364, t1792, t1872, t1796, t1159, t1175, t1181, t1137, t1260, t1213, t1232, t1259, t1231, t1283, t1121, t1122, t1117, t1103, t1260, t1213, t1232, t1259, t1231, t1283, t1113, t1116, t1108, t1101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1618, 0, qJDD(1), 0, g(3), 0, -t1590, 0, t1596, -t1586, -t1701, -t1596, -t1700, -qJDD(2), t1476, t1477, 0, pkin(1) * t1458, -t1727, -t1798, -t1803, t1727, t1840, -qJDD(2), t1253, t1263, t1264, t1186, -qJDD(2), t1803, -t1840, -t1727, -t1798, t1727, t1246, t1228, t1242, t1152, -t1379, -t1303, -t1351, t1377, -t1352, -t1406, t1154, t1158, t1136, t1110, -t1379, -t1303, -t1351, t1377, -t1352, -t1406, t1144, t1148, t1127, t1105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1589, t1590, 0, 0, t1529, t1515, t1532, t1528, t1530, 0, t1446, t1445, t1419, t1421, t1791, t1873, -t1355, t1795, t1357, t1794, t1227, t1236, t1172, t1162, t1794, t1355, -t1357, t1791, t1873, t1795, t1153, t1174, t1176, t1132, t1258, t1212, t1230, t1257, t1229, t1282, t1119, t1120, t1114, t1102, t1258, t1212, t1230, t1257, t1229, t1282, t1112, t1115, t1107, t1100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1580, t1581, t1587, -t1598, t1594, t1598, 0, -t1573, t1538, 0, t1672, t1863, -t1437, t1632, -t1440, t1653, t1348, t1375, t1252, -qJ(3) * t1273, t1653, t1437, t1440, t1672, t1863, t1632, t1211, t1256, t1261, t1194, t1324, t1272, t1293, t1323, t1292, t1354, t1149, t1150, t1140, t1111, t1324, t1272, t1293, t1323, t1292, t1354, t1135, t1138, t1128, t1106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1695, t1579, t1592, -t1656, t1588, -t1695, t1573, 0, t1539, 0, t1674, t1847, -t1431, t1673, t1433, t1652, t1320, t1326, t1245, t1262, t1652, t1431, -t1433, t1674, t1847, t1673, t1202, t1248, t1251, t1177, t1322, t1271, t1291, t1321, t1290, t1353, t1143, t1145, t1134, t1109, t1322, t1271, t1291, t1321, t1290, t1353, t1131, t1133, t1126, t1104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1596, t1586, t1701, t1596, t1700, qJDD(2), -t1538, -t1539, 0, 0, t1727, t1798, t1803, -t1727, -t1840, qJDD(2), -t1685 + t1864, t1554 - t1662 - t1888, t1784, t1785, qJDD(2), -t1803, t1840, t1727, t1798, -t1727, t1677, t1625 - t1858, t1622 + t1888, t1678, t1379, t1303, t1351, -t1377, t1352, t1406, t1647, t1664, t1639, t1654, t1379, t1303, t1351, -t1377, t1352, t1406, t1630, t1638, t1640, t1627; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1514, -t1839, -t1497, t1730, -t1547, -t1730, 0, -t1464, t1346, 0, -t1730, t1497, t1547, t1514, -t1839, t1730, t1295, qJ(4) * t1839, t1305, qJ(4) * t1325, t1467, t1460, t1391, -t1467, -t1388, t1510, t1182, t1185, t1241, t1142, t1467, t1460, t1391, -t1467, -t1388, t1510, t1164, t1171, t1214, t1125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1757, t1800, t1799, -t1684, t1838, -t1757, t1464, 0, t1347, 0, -t1757, -t1799, -t1838, t1757, t1800, -t1684, t1294, t1306, pkin(3) * t1800, pkin(3) * t1325, -t1378, -t1301, -t1349, t1376, -t1350, -t1405, t1198, t1201, t1151, t1139, -t1378, -t1301, -t1349, t1376, -t1350, -t1405, t1166, t1167, t1141, t1118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1727, t1798, t1803, -t1727, -t1840, qJDD(2), -t1346, -t1347, 0, 0, qJDD(2), -t1803, t1840, t1727, t1798, -t1727, t1704, t1625, t1622, t1707, t1379, t1303, t1351, -t1377, t1352, t1406, t1658, t1679, t1650, t1676, t1379, t1303, t1351, -t1377, t1352, t1406, t1633, t1649, t1651, t1629; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t1803, t1840, t1727, t1798, -t1727, 0, t1319, t1316, 0, t1379, t1303, t1351, -t1377, t1352, t1406, t1686, t1709, t1670, -t1767, t1379, t1303, t1351, -t1377, t1352, t1406, t1645, t1669, t1671, t1635; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1730, -t1497, -t1547, -t1514, t1839, -t1730, -t1319, 0, -t1325, 0, -t1467, -t1460, -t1391, t1467, t1388, -t1510, t1668, t1694, -t1298, -t1775, -t1467, -t1460, -t1391, t1467, t1388, -t1510, t1841, t1517 - t1340 - t1657, t1384 - t1298, t1688; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1757, t1799, t1838, -t1757, -t1800, t1684, -t1316, t1325, 0, 0, t1378, t1301, t1349, -t1376, t1350, t1405, t1667, t1666, t1642, pkin(8) * t1170 - t1774, t1378, t1301, t1349, -t1376, t1350, t1405, t1637, t1641, t1643, pkin(8) * t1157 - t1636; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1427, -t1387, -t1819, t1489, t1482, -t1489, 0, t1281, t1222, 0, t1427, -t1387, -t1819, t1489, t1482, -t1489, qJ(6) * t1819, t1234, t1180, qJ(6) * t1193; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1490, t1807, t1483, -t1655, t1411, -t1490, -t1281, 0, t1223, 0, t1490, t1807, t1483, -t1655, t1411, -t1490, t1210, t1318, t1178, t1165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1467, t1460, t1391, -t1467, -t1388, t1510, -t1222, -t1223, 0, 0, t1467, t1460, t1391, -t1467, -t1388, t1510, t1628, t1644, -t1384, -t1192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1427, -t1387, -t1819, t1489, t1482, -t1489, 0, t1244, t1193, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1490, t1807, t1483, -t1655, t1411, -t1490, -t1244, 0, t1195, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1467, t1460, t1391, -t1467, -t1388, t1510, -t1193, -t1195, 0, 0;];
m_new_reg  = t1;