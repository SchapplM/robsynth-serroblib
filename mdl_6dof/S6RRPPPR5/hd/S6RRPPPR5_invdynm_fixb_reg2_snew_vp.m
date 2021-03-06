% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6RRPPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 08:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6RRPPPR5_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR5_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR5_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 08:55:04
% EndTime: 2019-05-06 08:55:37
% DurationCPUTime: 35.49s
% Computational Cost: add. (76959->871), mult. (172643->927), div. (0->0), fcn. (111167->8), ass. (0->524)
t1629 = sin(qJ(1));
t1632 = cos(qJ(1));
t1626 = cos(pkin(9));
t1625 = sin(pkin(9));
t1628 = sin(qJ(2));
t1792 = qJD(1) * t1628;
t1584 = qJD(2) * t1625 + t1626 * t1792;
t1579 = t1584 ^ 2;
t1631 = cos(qJ(2));
t1623 = t1631 ^ 2;
t1633 = qJD(1) ^ 2;
t1733 = t1623 * t1633;
t1845 = t1579 + t1733;
t1616 = qJD(2) * t1792;
t1721 = t1631 * qJDD(1);
t1593 = -t1616 + t1721;
t1582 = -t1626 * qJD(2) + t1625 * t1792;
t1743 = t1584 * t1582;
t1841 = -t1743 + t1593;
t1882 = t1841 * t1625;
t1891 = t1845 * t1626 - t1882;
t1881 = t1841 * t1626;
t1433 = -t1845 * t1625 - t1881;
t1746 = t1582 * t1631;
t1563 = qJD(1) * t1746;
t1791 = qJD(1) * t1631;
t1713 = qJD(2) * t1791;
t1722 = t1628 * qJDD(1);
t1665 = t1713 + t1722;
t1644 = t1625 * qJDD(2) + t1626 * t1665;
t1850 = t1563 + t1644;
t1922 = t1433 * t1631 - t1850 * t1628;
t1926 = pkin(6) * (t1629 * t1891 + t1632 * t1922);
t1851 = -t1563 + t1644;
t1742 = t1584 * t1631;
t1562 = qJD(1) * t1742;
t1695 = qJDD(2) * t1626 - t1625 * t1665;
t1861 = t1695 - t1562;
t1383 = -t1625 * t1861 + t1851 * t1626;
t1578 = t1582 ^ 2;
t1488 = -t1579 - t1578;
t1879 = t1851 * t1625 + t1626 * t1861;
t1886 = t1488 * t1628 + t1631 * t1879;
t1894 = pkin(6) * (-t1383 * t1629 + t1632 * t1886);
t1927 = pkin(6) * (t1629 * t1922 - t1632 * t1891);
t1895 = pkin(6) * (t1383 * t1632 + t1629 * t1886);
t1514 = t1733 + t1578;
t1844 = t1593 + t1743;
t1864 = t1844 * t1625;
t1413 = t1626 * t1514 - t1864;
t1852 = t1562 + t1695;
t1344 = t1413 * t1631 + t1852 * t1628;
t1863 = t1844 * t1626;
t1410 = t1514 * t1625 + t1863;
t1928 = pkin(6) * (t1344 * t1632 + t1410 * t1629);
t1929 = pkin(6) * (t1344 * t1629 - t1410 * t1632);
t1914 = pkin(7) * t1886;
t1946 = -pkin(1) * t1383 - t1914;
t1553 = t1733 - t1578;
t1436 = t1553 * t1626 - t1882;
t1372 = t1436 * t1631 - t1861 * t1628;
t1427 = t1553 * t1625 + t1881;
t1945 = t1372 * t1629 - t1427 * t1632;
t1944 = t1372 * t1632 + t1427 * t1629;
t1341 = t1413 * t1628 - t1631 * t1852;
t1943 = pkin(1) * t1341;
t1888 = -t1488 * t1631 + t1628 * t1879;
t1919 = pkin(1) * t1888;
t1923 = t1433 * t1628 + t1631 * t1850;
t1942 = pkin(1) * t1923;
t1941 = pkin(7) * t1341;
t1913 = pkin(7) * t1888;
t1940 = pkin(7) * t1923;
t1939 = -pkin(1) * t1410 + pkin(7) * t1344;
t1938 = -pkin(1) * t1891 + pkin(7) * t1922;
t1774 = t1850 * t1625;
t1862 = t1852 * t1626;
t1387 = t1774 - t1862;
t1518 = -t1579 + t1578;
t1354 = t1387 * t1631 + t1518 * t1628;
t1773 = t1850 * t1626;
t1878 = t1852 * t1625 + t1773;
t1935 = t1354 * t1629 + t1632 * t1878;
t1934 = t1354 * t1632 - t1629 * t1878;
t1846 = -t1579 + t1733;
t1425 = t1846 * t1626 - t1864;
t1880 = t1846 * t1625 + t1863;
t1887 = -t1851 * t1628 + t1631 * t1880;
t1897 = t1425 * t1632 + t1629 * t1887;
t1896 = -t1425 * t1629 + t1632 * t1887;
t1898 = -pkin(2) * t1850 - qJ(3) * t1433;
t1363 = t1436 * t1628 + t1631 * t1861;
t1917 = pkin(2) * t1383;
t1916 = pkin(2) * t1410;
t1915 = pkin(2) * t1891;
t1912 = qJ(3) * t1383;
t1911 = qJ(3) * t1410;
t1910 = qJ(3) * t1413;
t1909 = qJ(3) * t1891;
t1719 = qJD(4) * t1791;
t1790 = qJD(3) * t1582;
t1848 = 0.2e1 * t1790 + 0.2e1 * t1719;
t1900 = t1848 - t1915;
t1899 = pkin(2) * t1852 - t1910;
t1890 = -pkin(2) * t1488 + qJ(3) * t1879;
t1351 = t1387 * t1628 - t1518 * t1631;
t1889 = t1628 * t1880 + t1631 * t1851;
t1876 = pkin(3) * t1852;
t1885 = qJ(4) * t1841;
t1877 = 2 * qJD(5);
t1805 = pkin(4) * t1851;
t1875 = qJ(4) * t1852;
t1874 = qJ(4) * t1861;
t1873 = qJ(5) * t1844;
t1551 = t1626 * t1562;
t1694 = t1625 * t1563;
t1479 = t1551 + t1694;
t1872 = t1479 * t1629;
t1871 = t1479 * t1632;
t1627 = sin(qJ(6));
t1630 = cos(qJ(6));
t1511 = t1582 * t1627 - t1630 * t1584;
t1513 = t1582 * t1630 + t1584 * t1627;
t1420 = t1513 * t1511;
t1586 = -qJDD(6) + t1593;
t1855 = -t1420 - t1586;
t1870 = t1627 * t1855;
t1552 = t1626 * t1563;
t1693 = t1625 * t1562;
t1670 = t1552 - t1693;
t1741 = t1593 * t1628;
t1819 = t1631 * t1670 - t1741;
t1869 = t1629 * t1819;
t1868 = t1630 * t1855;
t1867 = t1632 * t1819;
t1466 = -t1626 * t1695 + t1694;
t1698 = -t1625 * t1695 - t1552;
t1715 = t1628 * t1743;
t1817 = t1631 * t1698 - t1715;
t1828 = t1632 * t1466 + t1629 * t1817;
t1815 = -t1466 * t1629 + t1632 * t1817;
t1671 = t1626 * t1644 + t1693;
t1818 = t1631 * t1671 + t1715;
t1843 = t1625 * t1644 - t1551;
t1816 = t1629 * t1843 + t1632 * t1818;
t1827 = t1629 * t1818 - t1632 * t1843;
t1812 = pkin(4) + pkin(8);
t1550 = -pkin(5) * t1791 - pkin(8) * t1582;
t1602 = t1629 * g(1) - t1632 * g(2);
t1573 = qJDD(1) * pkin(1) + t1633 * pkin(7) + t1602;
t1592 = 0.2e1 * t1713 + t1722;
t1449 = (-t1593 + t1616) * pkin(2) - t1592 * qJ(3) - t1573;
t1603 = g(1) * t1632 + g(2) * t1629;
t1574 = -pkin(1) * t1633 + qJDD(1) * pkin(7) - t1603;
t1540 = -g(3) * t1628 + t1631 * t1574;
t1808 = pkin(2) * t1631;
t1687 = -qJ(3) * t1628 - t1808;
t1591 = t1687 * qJD(1);
t1814 = qJD(2) ^ 2;
t1464 = -pkin(2) * t1814 + qJDD(2) * qJ(3) + t1591 * t1791 + t1540;
t1700 = -t1626 * t1449 + t1625 * t1464;
t1663 = t1593 * pkin(3) - qJ(4) * t1733 + qJDD(4) + t1700;
t1638 = t1791 * t1877 + t1663 + t1805 + t1873;
t1515 = pkin(3) * t1582 - qJ(4) * t1584;
t1728 = 0.2e1 * qJD(3) + t1515;
t1255 = t1644 * pkin(8) + t1550 * t1791 + (-pkin(5) * t1584 + t1728) * t1584 + t1638;
t1548 = pkin(4) * t1584 + qJ(5) * t1791;
t1725 = t1625 * t1449 + t1626 * t1464;
t1679 = pkin(3) * t1733 + t1593 * qJ(4) + t1582 * t1515 - t1725;
t1654 = -pkin(4) * t1695 + t1578 * qJ(5) + t1548 * t1791 - qJDD(5) + t1679;
t1568 = -0.2e1 * t1790;
t1849 = t1568 - 0.2e1 * t1719;
t1285 = -t1654 + t1849;
t1256 = -pkin(5) * t1841 + pkin(8) * t1861 + t1285;
t1200 = t1627 * t1255 - t1630 * t1256;
t1201 = t1630 * t1255 + t1627 * t1256;
t1856 = -t1630 * t1200 + t1627 * t1201;
t1860 = t1856 * t1812;
t1857 = (qJD(1) * t1591 + t1574) * t1628;
t1169 = t1627 * t1200 + t1201 * t1630;
t1392 = -t1511 * qJD(6) + t1627 * t1644 - t1630 * t1695;
t1610 = -qJD(6) + t1791;
t1491 = t1511 * t1610;
t1854 = t1491 + t1392;
t1789 = qJD(3) * t1584;
t1847 = -t1584 * t1515 - 0.2e1 * t1789;
t1576 = t1578 * pkin(4);
t1842 = (-0.2e1 * qJD(4) - t1548) * t1584 - t1576;
t1840 = -t1644 * qJ(4) - t1876;
t1787 = qJD(5) * t1582;
t1839 = t1584 * t1548 + t1576 - 0.2e1 * t1787;
t1838 = -qJ(5) * (t1852 + t1695) + pkin(4) * t1514;
t1800 = t1631 * g(3);
t1674 = -qJDD(2) * pkin(2) - qJ(3) * t1814 + qJDD(3) + t1800;
t1658 = t1674 + t1840;
t1730 = t1628 * t1574;
t1636 = -qJD(1) * (qJ(4) * t1746 - t1591 * t1628) + t1658 + t1730;
t1167 = pkin(5) * t1856;
t1811 = pkin(3) + qJ(5);
t1837 = qJ(4) * t1856 - t1169 * t1811 + t1167;
t1463 = t1674 + t1857;
t1794 = t1695 * qJ(5);
t1267 = -t1794 - t1644 * pkin(5) - t1579 * pkin(8) + (-qJ(4) * t1791 + t1550 + t1877) * t1582 + t1463 + t1840 + t1842;
t1263 = t1627 * t1267;
t1510 = t1511 ^ 2;
t1608 = t1610 ^ 2;
t1417 = -t1608 - t1510;
t1309 = t1627 * t1417 + t1868;
t1836 = -t1309 * t1812 + t1263;
t1704 = t1728 * t1584;
t1273 = t1704 + t1638;
t1835 = qJ(4) * t1285 - t1273 * t1811;
t1699 = -t1627 * t1695 - t1630 * t1644;
t1337 = (qJD(6) + t1610) * t1513 + t1699;
t1340 = -t1491 + t1392;
t1287 = -t1337 * t1627 - t1340 * t1630;
t1834 = -t1287 * t1812 - t1856;
t1284 = pkin(5) * t1287;
t1289 = -t1337 * t1630 + t1627 * t1340;
t1833 = qJ(4) * t1287 - t1289 * t1811 + t1284;
t1310 = t1417 * t1630 - t1870;
t1711 = pkin(5) * t1309 - t1200;
t1832 = qJ(4) * t1309 - t1310 * t1811 + t1711;
t1299 = -t1679 + t1849;
t1301 = t1704 + t1663;
t1241 = t1299 * t1626 + t1301 * t1625;
t1788 = qJD(4) * t1584;
t1311 = t1636 - 0.2e1 * t1788;
t1707 = qJ(4) * t1625 + pkin(2);
t1831 = -t1311 * (pkin(3) * t1626 + t1707) + qJ(3) * t1241;
t1813 = t1513 ^ 2;
t1708 = -t1608 - t1813;
t1396 = -t1420 + t1586;
t1731 = t1627 * t1396;
t1314 = t1630 * t1708 + t1731;
t1729 = t1630 * t1267;
t1830 = -t1314 * t1812 + t1729;
t1785 = t1396 * t1630;
t1315 = -t1627 * t1708 + t1785;
t1712 = -pkin(5) * t1314 + t1201;
t1829 = qJ(4) * t1314 - t1315 * t1811 - t1712;
t1303 = t1311 - t1876;
t1826 = t1626 * t1303 - t1707 * t1852 + t1910;
t1825 = t1811 * t1851 - t1874;
t1795 = qJ(4) * t1514;
t1824 = -t1811 * t1844 - t1795;
t1823 = t1811 * t1845 - t1654 - t1885;
t1575 = t1631 * t1593;
t1822 = t1628 * t1670 + t1575;
t1714 = t1582 * t1742;
t1821 = t1628 * t1671 - t1714;
t1820 = t1628 * t1698 + t1714;
t1809 = pkin(2) * t1628;
t1807 = pkin(4) * t1273;
t1806 = pkin(4) * t1285;
t1804 = pkin(4) * t1844;
t1803 = pkin(4) * t1841;
t1801 = pkin(4) * t1845;
t1797 = qJ(4) * t1488;
t1793 = qJD(1) * qJD(2);
t1784 = t1463 * t1625;
t1783 = t1463 * t1626;
t1748 = t1573 * t1628;
t1747 = t1573 * t1631;
t1530 = t1592 * t1628;
t1609 = t1631 * t1633 * t1628;
t1600 = -t1609 + qJDD(2);
t1740 = t1600 * t1628;
t1739 = t1600 * t1631;
t1601 = qJDD(2) + t1609;
t1738 = t1601 * t1628;
t1737 = t1610 * t1513;
t1736 = t1610 * t1627;
t1735 = t1610 * t1630;
t1622 = t1628 ^ 2;
t1734 = t1622 * t1633;
t1726 = -pkin(5) * t1267 + pkin(8) * t1169;
t1723 = t1622 + t1623;
t1720 = t1608 - t1813;
t1717 = t1628 * t1420;
t1716 = t1631 * t1420;
t1710 = -pkin(5) * t1854 + pkin(8) * t1315 + t1263;
t1336 = (qJD(6) - t1610) * t1513 + t1699;
t1709 = -pkin(5) * t1336 + pkin(8) * t1310 - t1729;
t1348 = t1568 + t1725;
t1347 = t1700 + 0.2e1 * t1789;
t1291 = t1347 * t1625 + t1626 * t1348;
t1539 = t1730 + t1800;
t1452 = t1539 * t1628 + t1631 * t1540;
t1696 = -t1602 * t1629 - t1632 * t1603;
t1692 = t1629 * t1609;
t1691 = t1632 * t1609;
t1690 = -pkin(4) * t1169 - t1726;
t1597 = qJDD(1) * t1632 - t1629 * t1633;
t1689 = -pkin(6) * t1597 - g(3) * t1629;
t1688 = -pkin(2) * t1463 + qJ(3) * t1291;
t1686 = -pkin(3) * t1301 + qJ(4) * t1299;
t1685 = -pkin(3) * t1851 + t1874;
t1290 = -t1347 * t1626 + t1348 * t1625;
t1451 = t1539 * t1631 - t1540 * t1628;
t1678 = t1602 * t1632 - t1603 * t1629;
t1676 = -pkin(4) * t1315 - t1710;
t1675 = -pkin(4) * t1310 - t1709;
t1376 = -t1510 - t1813;
t1669 = -pkin(5) * t1376 + pkin(8) * t1289 + t1169;
t1667 = -t1783 + t1899;
t1666 = -qJD(6) * t1513 - t1699;
t1662 = t1784 + t1898;
t1660 = -pkin(4) * t1289 - t1669;
t1659 = t1291 + t1890;
t1567 = 0.2e1 * t1788;
t1244 = t1567 - t1636 - t1838 + t1839 + t1876;
t1395 = t1804 + t1875;
t1657 = t1244 * t1626 + t1395 * t1625 + t1899;
t1302 = t1567 - t1857 + (t1850 + t1563) * qJ(4) - t1658;
t1270 = t1302 + t1794 - t1801 + t1839;
t1349 = t1811 * t1850 - t1803;
t1656 = t1270 * t1625 + t1349 * t1626 - t1898;
t1655 = pkin(3) * t1773 + t1302 * t1625 - t1898;
t1146 = -qJ(4) * t1267 - t1690;
t1147 = -t1267 * t1811 + t1860;
t1152 = t1169 * t1625 + t1626 * t1856;
t1653 = -pkin(2) * t1267 + qJ(3) * t1152 + t1146 * t1625 + t1147 * t1626;
t1157 = -qJ(4) * t1376 - t1660;
t1158 = -t1376 * t1811 - t1834;
t1220 = t1287 * t1626 + t1289 * t1625;
t1652 = -pkin(2) * t1376 + qJ(3) * t1220 + t1157 * t1625 + t1158 * t1626;
t1180 = -qJ(4) * t1336 - t1675;
t1181 = -t1336 * t1811 - t1836;
t1258 = t1309 * t1626 + t1310 * t1625;
t1651 = -pkin(2) * t1336 + qJ(3) * t1258 + t1180 * t1625 + t1181 * t1626;
t1185 = -qJ(4) * t1854 - t1676;
t1186 = -t1811 * t1854 - t1830;
t1272 = t1314 * t1626 + t1315 * t1625;
t1650 = -pkin(2) * t1854 + qJ(3) * t1272 + t1185 * t1625 + t1186 * t1626;
t1634 = t1636 + 0.2e1 * t1787 + t1842;
t1293 = t1634 - t1794;
t1195 = -t1293 * t1811 + t1806;
t1214 = t1273 * t1625 + t1285 * t1626;
t1223 = -qJ(4) * t1293 + t1807;
t1649 = -pkin(2) * t1293 + qJ(3) * t1214 + t1195 * t1626 + t1223 * t1625;
t1641 = -pkin(4) * t1861 + t1654;
t1242 = t1488 * t1811 + t1641 + t1848;
t1637 = -t1638 + t1847;
t1243 = t1637 + t1797 - t1805;
t1648 = t1242 * t1626 + t1243 * t1625 - t1890;
t1292 = -pkin(3) * t1488 + t1299;
t1294 = t1301 - t1797;
t1647 = t1292 * t1626 + t1294 * t1625 + t1890;
t1646 = pkin(3) * t1845 - t1679 - t1885;
t1642 = pkin(3) * t1844 + t1663 + t1795;
t1607 = -t1733 - t1814;
t1606 = t1733 - t1814;
t1605 = -t1734 - t1814;
t1604 = -t1734 + t1814;
t1599 = (t1622 - t1623) * t1633;
t1598 = t1723 * t1633;
t1596 = qJDD(1) * t1629 + t1632 * t1633;
t1595 = t1723 * qJDD(1);
t1594 = -0.2e1 * t1616 + t1721;
t1589 = t1631 * t1601;
t1588 = t1723 * t1793;
t1564 = -pkin(6) * t1596 + g(3) * t1632;
t1547 = -t1622 * t1793 + t1631 * t1665;
t1546 = -t1623 * t1793 - t1741;
t1538 = -t1605 * t1628 - t1739;
t1537 = -t1604 * t1628 + t1589;
t1536 = t1607 * t1631 - t1738;
t1535 = t1606 * t1631 - t1740;
t1534 = t1605 * t1631 - t1740;
t1533 = t1604 * t1631 + t1738;
t1532 = t1607 * t1628 + t1589;
t1531 = t1606 * t1628 + t1739;
t1529 = -t1628 * t1713 + t1575;
t1517 = t1594 * t1631 - t1530;
t1516 = t1592 * t1631 + t1594 * t1628;
t1481 = (t1582 * t1626 - t1584 * t1625) * t1791;
t1480 = (t1582 * t1625 + t1584 * t1626) * t1791;
t1477 = t1510 - t1608;
t1475 = -pkin(7) * t1534 - t1747;
t1474 = -pkin(7) * t1532 - t1748;
t1462 = -pkin(1) * t1534 + t1540;
t1461 = -pkin(1) * t1532 + t1539;
t1448 = pkin(1) * t1594 + pkin(7) * t1536 + t1747;
t1447 = -pkin(1) * t1592 + pkin(7) * t1538 - t1748;
t1442 = t1481 * t1631 - t1741;
t1439 = t1481 * t1628 + t1575;
t1419 = -t1510 + t1813;
t1418 = pkin(1) * t1573 + pkin(7) * t1452;
t1416 = pkin(1) * t1598 + pkin(7) * t1595 + t1452;
t1379 = (t1511 * t1630 - t1513 * t1627) * t1610;
t1378 = (t1511 * t1627 + t1513 * t1630) * t1610;
t1356 = t1783 + t1909;
t1330 = t1477 * t1630 + t1731;
t1329 = -t1627 * t1720 + t1868;
t1328 = t1477 * t1627 - t1785;
t1327 = t1630 * t1720 + t1870;
t1326 = t1392 * t1630 + t1513 * t1736;
t1325 = t1392 * t1627 - t1513 * t1735;
t1324 = t1511 * t1735 + t1627 * t1666;
t1323 = t1511 * t1736 - t1630 * t1666;
t1322 = t1784 + t1911;
t1306 = -t1685 + t1917;
t1305 = t1378 * t1626 + t1379 * t1625;
t1304 = t1378 * t1625 - t1379 * t1626;
t1300 = t1348 + t1915;
t1298 = t1305 * t1631 - t1586 * t1628;
t1297 = t1305 * t1628 + t1586 * t1631;
t1296 = t1347 + t1916;
t1295 = -t1825 - t1917;
t1288 = -t1336 * t1630 - t1627 * t1854;
t1286 = -t1627 * t1336 + t1630 * t1854;
t1282 = t1328 * t1626 + t1330 * t1625;
t1281 = t1327 * t1626 + t1329 * t1625;
t1280 = t1328 * t1625 - t1330 * t1626;
t1279 = t1327 * t1625 - t1329 * t1626;
t1278 = t1325 * t1626 + t1326 * t1625;
t1277 = -t1323 * t1626 - t1324 * t1625;
t1276 = t1325 * t1625 - t1326 * t1626;
t1275 = -t1323 * t1625 + t1324 * t1626;
t1274 = -t1662 + t1942;
t1271 = t1314 * t1625 - t1315 * t1626;
t1269 = -t1667 + t1943;
t1268 = -pkin(3) * t1774 + t1302 * t1626 - t1909;
t1262 = t1291 * t1631 + t1463 * t1628;
t1261 = t1291 * t1628 - t1463 * t1631;
t1260 = -qJ(4) * t1862 - t1303 * t1625 - t1911;
t1259 = -t1646 + t1900;
t1257 = t1309 * t1625 - t1310 * t1626;
t1250 = -t1290 + t1912;
t1249 = -t1642 + t1847 - t1916;
t1248 = t1278 * t1631 + t1717;
t1247 = t1277 * t1631 - t1717;
t1246 = t1278 * t1628 - t1716;
t1245 = t1277 * t1628 + t1716;
t1240 = t1299 * t1625 - t1301 * t1626;
t1239 = t1282 * t1631 - t1337 * t1628;
t1238 = t1281 * t1631 + t1340 * t1628;
t1237 = t1282 * t1628 + t1337 * t1631;
t1236 = t1281 * t1628 - t1340 * t1631;
t1235 = -t1823 + t1900;
t1234 = t1272 * t1631 + t1628 * t1854;
t1233 = t1272 * t1628 - t1631 * t1854;
t1232 = -t1300 * t1628 + t1356 * t1631 + t1940;
t1231 = -t1655 - t1942;
t1230 = t1258 * t1631 + t1336 * t1628;
t1229 = t1258 * t1628 - t1336 * t1631;
t1228 = -t1296 * t1628 + t1322 * t1631 + t1941;
t1227 = t1273 - t1824 + t1916;
t1226 = -t1826 - t1943;
t1225 = t1270 * t1626 - t1349 * t1625 - t1909;
t1224 = t1300 * t1631 + t1356 * t1628 - t1938;
t1222 = -t1244 * t1625 + t1395 * t1626 + t1911;
t1221 = -t1659 - t1919;
t1219 = t1286 * t1626 + t1288 * t1625;
t1218 = t1287 * t1625 - t1289 * t1626;
t1217 = t1286 * t1625 - t1288 * t1626;
t1216 = t1241 * t1631 + t1311 * t1628;
t1215 = t1241 * t1628 - t1311 * t1631;
t1213 = -t1273 * t1626 + t1285 * t1625;
t1212 = t1296 * t1631 + t1322 * t1628 - t1939;
t1211 = -t1292 * t1625 + t1294 * t1626 + t1912;
t1210 = t1250 * t1631 - t1383 * t1809 - t1913;
t1209 = t1219 * t1631 + t1419 * t1628;
t1208 = t1219 * t1628 - t1419 * t1631;
t1207 = -pkin(1) * t1261 - t1688;
t1206 = t1220 * t1631 + t1376 * t1628;
t1205 = t1220 * t1628 - t1376 * t1631;
t1204 = -t1656 - t1942;
t1203 = t1914 + t1628 * t1250 - (-pkin(1) - t1808) * t1383;
t1202 = -t1657 + t1943;
t1196 = -qJ(3) * t1240 + (pkin(3) * t1625 - qJ(4) * t1626) * t1311;
t1194 = -t1259 * t1628 + t1268 * t1631 - t1940;
t1193 = -t1647 - t1919;
t1192 = -pkin(2) * t1240 - t1686;
t1191 = t1214 * t1631 + t1293 * t1628;
t1190 = t1214 * t1628 - t1293 * t1631;
t1189 = -t1249 * t1628 + t1260 * t1631 - t1941;
t1188 = t1259 * t1631 + t1268 * t1628 + t1938;
t1187 = -pkin(7) * t1261 + (-qJ(3) * t1631 + t1809) * t1290;
t1184 = -t1242 * t1625 + t1243 * t1626 - t1912;
t1183 = t1249 * t1631 + t1260 * t1628 + t1939;
t1182 = t1211 * t1631 - t1306 * t1628 - t1913;
t1179 = t1211 * t1628 + t1306 * t1631 - t1946;
t1178 = t1225 * t1631 - t1235 * t1628 - t1940;
t1177 = -t1648 + t1919;
t1176 = t1225 * t1628 + t1235 * t1631 + t1938;
t1175 = pkin(7) * t1262 + (-pkin(1) + t1687) * t1290;
t1174 = t1222 * t1631 - t1227 * t1628 + t1941;
t1173 = t1184 * t1631 - t1295 * t1628 + t1913;
t1172 = t1222 * t1628 + t1227 * t1631 - t1939;
t1171 = t1184 * t1628 + t1295 * t1631 + t1946;
t1170 = -pkin(2) * t1213 - t1835;
t1165 = -pkin(1) * t1215 - t1831;
t1164 = -pkin(2) * t1271 - t1829;
t1163 = -pkin(2) * t1218 - t1833;
t1162 = -pkin(2) * t1257 - t1832;
t1161 = -qJ(3) * t1213 - t1195 * t1625 + t1223 * t1626;
t1160 = -qJ(3) * t1271 + t1185 * t1626 - t1186 * t1625;
t1159 = -pkin(7) * t1215 - t1192 * t1628 + t1196 * t1631;
t1156 = -qJ(3) * t1257 + t1180 * t1626 - t1181 * t1625;
t1155 = -pkin(1) * t1240 + pkin(7) * t1216 + t1192 * t1631 + t1196 * t1628;
t1154 = -pkin(1) * t1233 - t1650;
t1153 = -pkin(1) * t1229 - t1651;
t1151 = -t1169 * t1626 + t1625 * t1856;
t1150 = -pkin(1) * t1190 - t1649;
t1149 = t1152 * t1631 + t1267 * t1628;
t1148 = t1152 * t1628 - t1267 * t1631;
t1145 = -pkin(7) * t1233 + t1160 * t1631 - t1164 * t1628;
t1144 = -pkin(7) * t1190 + t1161 * t1631 - t1170 * t1628;
t1143 = -pkin(1) * t1271 + pkin(7) * t1234 + t1160 * t1628 + t1164 * t1631;
t1142 = -pkin(7) * t1229 + t1156 * t1631 - t1162 * t1628;
t1141 = -pkin(1) * t1213 + pkin(7) * t1191 + t1161 * t1628 + t1170 * t1631;
t1140 = -pkin(1) * t1257 + pkin(7) * t1230 + t1156 * t1628 + t1162 * t1631;
t1139 = -qJ(3) * t1218 + t1157 * t1626 - t1158 * t1625;
t1138 = -pkin(1) * t1205 - t1652;
t1137 = -pkin(7) * t1205 + t1139 * t1631 - t1163 * t1628;
t1136 = -pkin(1) * t1218 + pkin(7) * t1206 + t1139 * t1628 + t1163 * t1631;
t1135 = -pkin(2) * t1151 - t1837;
t1134 = -qJ(3) * t1151 + t1146 * t1626 - t1147 * t1625;
t1133 = -pkin(1) * t1148 - t1653;
t1132 = -pkin(7) * t1148 + t1134 * t1631 - t1135 * t1628;
t1131 = -pkin(1) * t1151 + pkin(7) * t1149 + t1134 * t1628 + t1135 * t1631;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1597, 0, -t1596, 0, t1689, -t1564, -t1678, -pkin(6) * t1678, t1547 * t1632 - t1692, t1517 * t1632 + t1599 * t1629, t1537 * t1632 + t1629 * t1722, t1546 * t1632 + t1692, t1535 * t1632 + t1629 * t1721, qJDD(2) * t1629 + t1588 * t1632, t1632 * t1474 - t1629 * t1461 - pkin(6) * (t1536 * t1629 + t1594 * t1632), t1632 * t1475 - t1629 * t1462 - pkin(6) * (t1538 * t1629 - t1592 * t1632), t1632 * t1451 - pkin(6) * (t1595 * t1629 + t1598 * t1632), -pkin(6) * (t1452 * t1629 + t1573 * t1632) - (pkin(1) * t1629 - pkin(7) * t1632) * t1451, t1816, -t1934, -t1896, t1815, -t1944, t1442 * t1632 + t1872, t1632 * t1228 - t1629 * t1269 + t1929, t1632 * t1232 - t1629 * t1274 + t1927, t1632 * t1210 - t1629 * t1221 - t1895, t1632 * t1187 - t1629 * t1207 - pkin(6) * (t1262 * t1629 - t1290 * t1632), t1480 * t1629 + t1867, t1896, t1944, t1816, -t1934, t1815, t1632 * t1182 - t1629 * t1193 - t1895, t1632 * t1189 - t1629 * t1226 - t1929, t1632 * t1194 - t1629 * t1231 - t1927, t1632 * t1159 - t1629 * t1165 - pkin(6) * (t1216 * t1629 - t1240 * t1632), t1815, -t1944, t1934, t1867 + t1872, t1896, t1816, t1632 * t1178 - t1629 * t1204 - t1927, t1632 * t1173 - t1629 * t1177 + t1895, t1632 * t1174 - t1629 * t1202 + t1929, t1632 * t1144 - t1629 * t1150 - pkin(6) * (t1191 * t1629 - t1213 * t1632), t1248 * t1632 + t1276 * t1629, t1209 * t1632 + t1217 * t1629, t1238 * t1632 + t1279 * t1629, t1247 * t1632 + t1275 * t1629, t1239 * t1632 + t1280 * t1629, t1298 * t1632 + t1304 * t1629, t1632 * t1142 - t1629 * t1153 - pkin(6) * (t1230 * t1629 - t1257 * t1632), t1632 * t1145 - t1629 * t1154 - pkin(6) * (t1234 * t1629 - t1271 * t1632), t1632 * t1137 - t1629 * t1138 - pkin(6) * (t1206 * t1629 - t1218 * t1632), t1632 * t1132 - t1629 * t1133 - pkin(6) * (t1149 * t1629 - t1151 * t1632); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1596, 0, t1597, 0, t1564, t1689, t1696, pkin(6) * t1696, t1547 * t1629 + t1691, t1517 * t1629 - t1599 * t1632, t1537 * t1629 - t1632 * t1722, t1546 * t1629 - t1691, t1535 * t1629 - t1632 * t1721, -qJDD(2) * t1632 + t1588 * t1629, t1629 * t1474 + t1632 * t1461 + pkin(6) * (t1536 * t1632 - t1594 * t1629), t1629 * t1475 + t1632 * t1462 + pkin(6) * (t1538 * t1632 + t1592 * t1629), t1629 * t1451 + pkin(6) * (t1595 * t1632 - t1598 * t1629), pkin(6) * (t1452 * t1632 - t1573 * t1629) - (-pkin(1) * t1632 - pkin(7) * t1629) * t1451, t1827, -t1935, -t1897, t1828, -t1945, t1442 * t1629 - t1871, t1629 * t1228 + t1632 * t1269 - t1928, t1629 * t1232 + t1632 * t1274 - t1926, t1629 * t1210 + t1632 * t1221 + t1894, t1629 * t1187 + t1632 * t1207 + pkin(6) * (t1262 * t1632 + t1290 * t1629), -t1480 * t1632 + t1869, t1897, t1945, t1827, -t1935, t1828, t1629 * t1182 + t1632 * t1193 + t1894, t1629 * t1189 + t1632 * t1226 + t1928, t1629 * t1194 + t1632 * t1231 + t1926, t1629 * t1159 + t1632 * t1165 + pkin(6) * (t1216 * t1632 + t1240 * t1629), t1828, -t1945, t1935, t1869 - t1871, t1897, t1827, t1629 * t1178 + t1632 * t1204 + t1926, t1629 * t1173 + t1632 * t1177 - t1894, t1629 * t1174 + t1632 * t1202 - t1928, t1629 * t1144 + t1632 * t1150 + pkin(6) * (t1191 * t1632 + t1213 * t1629), t1248 * t1629 - t1276 * t1632, t1209 * t1629 - t1217 * t1632, t1238 * t1629 - t1279 * t1632, t1247 * t1629 - t1275 * t1632, t1239 * t1629 - t1280 * t1632, t1298 * t1629 - t1304 * t1632, t1629 * t1142 + t1632 * t1153 + pkin(6) * (t1230 * t1632 + t1257 * t1629), t1629 * t1145 + t1632 * t1154 + pkin(6) * (t1234 * t1632 + t1271 * t1629), t1629 * t1137 + t1632 * t1138 + pkin(6) * (t1206 * t1632 + t1218 * t1629), t1629 * t1132 + t1632 * t1133 + pkin(6) * (t1149 * t1632 + t1151 * t1629); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1602, t1603, 0, 0, t1530, t1516, t1533, t1529, t1531, 0, t1448, t1447, t1416, t1418, t1821, -t1351, -t1889, t1820, -t1363, t1439, t1212, t1224, t1203, t1175, t1822, t1889, t1363, t1821, -t1351, t1820, t1179, t1183, t1188, t1155, t1820, -t1363, t1351, t1822, t1889, t1821, t1176, t1171, t1172, t1141, t1246, t1208, t1236, t1245, t1237, t1297, t1140, t1143, t1136, t1131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1633, 0, 0, -g(3), -t1602, 0, t1547, t1517, t1537, t1546, t1535, t1588, t1474, t1475, t1451, pkin(7) * t1451, t1818, -t1354, -t1887, t1817, -t1372, t1442, t1228, t1232, t1210, t1187, t1819, t1887, t1372, t1818, -t1354, t1817, t1182, t1189, t1194, t1159, t1817, -t1372, t1354, t1819, t1887, t1818, t1178, t1173, t1174, t1144, t1248, t1209, t1238, t1247, t1239, t1298, t1142, t1145, t1137, t1132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1633, 0, qJDD(1), 0, g(3), 0, -t1603, 0, t1609, -t1599, -t1722, -t1609, -t1721, -qJDD(2), t1461, t1462, 0, pkin(1) * t1451, -t1843, -t1878, -t1425, t1466, t1427, -t1479, t1269, t1274, t1221, t1207, -t1480, t1425, -t1427, -t1843, -t1878, t1466, t1193, t1226, t1231, t1165, t1466, t1427, t1878, -t1479, t1425, -t1843, t1204, t1177, t1202, t1150, -t1276, -t1217, -t1279, -t1275, -t1280, -t1304, t1153, t1154, t1138, t1133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1602, t1603, 0, 0, t1530, t1516, t1533, t1529, t1531, 0, t1448, t1447, t1416, t1418, t1821, -t1351, -t1889, t1820, -t1363, t1439, t1212, t1224, t1203, t1175, t1822, t1889, t1363, t1821, -t1351, t1820, t1179, t1183, t1188, t1155, t1820, -t1363, t1351, t1822, t1889, t1821, t1176, t1171, t1172, t1141, t1246, t1208, t1236, t1245, t1237, t1297, t1140, t1143, t1136, t1131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1665, t1594, t1601, -t1713, t1606, t1713, 0, -t1573, t1539, 0, t1671, -t1387, -t1880, t1698, -t1436, t1481, t1322, t1356, t1250, -qJ(3) * t1290, t1670, t1880, t1436, t1671, -t1387, t1698, t1211, t1260, t1268, t1196, t1698, -t1436, t1387, t1670, t1880, t1671, t1225, t1184, t1222, t1161, t1278, t1219, t1281, t1277, t1282, t1305, t1156, t1160, t1139, t1134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1616, t1592, t1604, t1593, t1600, -t1616, t1573, 0, t1540, 0, -t1743, t1518, -t1851, t1743, -t1861, t1593, t1296, t1300, t1917, -pkin(2) * t1290, t1593, t1851, t1861, -t1743, t1518, t1743, t1306, t1249, t1259, t1192, t1743, -t1861, -t1518, t1593, t1851, -t1743, t1235, t1295, t1227, t1170, -t1420, -t1419, -t1340, t1420, t1337, t1586, t1162, t1164, t1163, t1135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1609, t1599, t1722, t1609, t1721, qJDD(2), -t1539, -t1540, 0, 0, t1843, t1878, t1425, -t1466, -t1427, t1479, t1667, t1662, t1659, t1688, t1480, -t1425, t1427, t1843, t1878, -t1466, t1647, t1826, t1655, t1831, -t1466, -t1427, -t1878, t1479, -t1425, t1843, t1656, t1648, t1657, t1649, t1276, t1217, t1279, t1275, t1280, t1304, t1651, t1650, t1652, t1653; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1644, t1852, -t1844, -t1563, -t1553, t1563, 0, t1463, t1347, 0, t1563, t1844, t1553, t1644, t1852, -t1563, t1294, -t1875, t1302, -qJ(4) * t1311, -t1563, -t1553, -t1852, t1563, t1844, t1644, t1270, t1243, t1395, t1223, t1325, t1286, t1327, -t1323, t1328, t1378, t1180, t1185, t1157, t1146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1562, t1850, t1846, t1695, -t1841, t1562, -t1463, 0, t1348, 0, t1562, -t1846, t1841, -t1562, t1850, t1695, t1292, t1303, pkin(3) * t1850, -pkin(3) * t1311, t1695, -t1841, -t1850, t1562, -t1846, -t1562, t1349, t1242, t1244, t1195, -t1326, -t1288, -t1329, t1324, -t1330, -t1379, t1181, t1186, t1158, t1147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1743, -t1518, t1851, -t1743, t1861, -t1593, -t1347, -t1348, 0, 0, -t1593, -t1851, -t1861, t1743, -t1518, -t1743, t1685, t1704 + t1642, t1646 + t1849, t1686, -t1743, t1861, t1518, -t1593, -t1851, t1743, t1823 + t1849, t1825, t1637 + t1824, t1835, t1420, t1419, t1340, -t1420, -t1337, -t1586, t1832, t1829, t1833, t1837; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1593, -t1851, -t1861, t1743, -t1518, -t1743, 0, t1301, t1299, 0, -t1743, t1861, t1518, -t1593, -t1851, t1743, qJ(5) * t1845 + t1285, qJ(5) * t1851, t1637 - t1873, -qJ(5) * t1273, t1420, t1419, t1340, -t1420, -t1337, -t1586, -qJ(5) * t1310 + t1711, -qJ(5) * t1315 - t1712, -qJ(5) * t1289 + t1284, -qJ(5) * t1169 + t1167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1563, -t1844, -t1553, -t1644, -t1852, t1563, -t1301, 0, t1311, 0, t1563, t1553, t1852, -t1563, -t1844, -t1644, t1293 + t1801, t1273 + t1805, -t1804, -t1807, -t1325, -t1286, -t1327, t1323, -t1328, -t1378, t1675, t1676, t1660, t1690; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1562, t1846, -t1841, t1562, -t1850, -t1695, -t1299, -t1311, 0, 0, -t1695, t1841, t1850, -t1562, t1846, t1562, -qJ(5) * t1850 + t1803, -qJ(5) * t1488 - t1641 + t1849, t1634 + t1838, qJ(5) * t1293 - t1806, t1326, t1288, t1329, -t1324, t1330, t1379, qJ(5) * t1336 + t1836, qJ(5) * t1854 + t1830, qJ(5) * t1376 + t1834, qJ(5) * t1267 - t1860; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1695, t1841, t1850, -t1562, t1846, t1562, 0, t1285, t1293, 0, t1326, t1288, t1329, -t1324, t1330, t1379, -pkin(8) * t1309 + t1263, -pkin(8) * t1314 + t1729, -pkin(8) * t1287 - t1856, -pkin(8) * t1856; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1743, -t1861, -t1518, t1593, t1851, -t1743, -t1285, 0, t1273, 0, -t1420, -t1419, -t1340, t1420, t1337, t1586, -t1711, t1712, -t1284, -t1167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1563, -t1553, -t1852, t1563, t1844, t1644, -t1293, -t1273, 0, 0, t1325, t1286, t1327, -t1323, t1328, t1378, t1709, t1710, t1669, t1726; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1392, -t1336, t1855, -t1491, t1477, t1491, 0, t1267, t1200, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1737, t1854, t1720, t1666, -t1396, t1737, -t1267, 0, t1201, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1420, t1419, t1340, -t1420, -t1337, -t1586, -t1200, -t1201, 0, 0;];
m_new_reg  = t1;
