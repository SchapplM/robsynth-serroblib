% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6RPRRPR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 22:53
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6RPRRPR6_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR6_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR6_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_invdynm_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 22:51:58
% EndTime: 2019-05-05 22:52:58
% DurationCPUTime: 63.45s
% Computational Cost: add. (719978->1007), mult. (1743907->1437), div. (0->0), fcn. (1364392->12), ass. (0->688)
t1798 = sin(pkin(10));
t1800 = cos(pkin(10));
t1804 = sin(qJ(3));
t1808 = cos(qJ(3));
t1769 = (t1798 * t1804 - t1800 * t1808) * qJD(1);
t1762 = qJD(4) + t1769;
t1930 = t1762 ^ 2;
t1836 = t1798 * t1808 + t1800 * t1804;
t1771 = t1836 * qJD(1);
t1803 = sin(qJ(4));
t1807 = cos(qJ(4));
t1749 = -t1807 * qJD(3) + t1771 * t1803;
t1750 = qJD(3) * t1803 + t1771 * t1807;
t1797 = sin(pkin(11));
t1799 = cos(pkin(11));
t1708 = t1749 * t1799 + t1750 * t1797;
t1932 = t1708 ^ 2;
t1625 = -t1930 - t1932;
t1709 = -t1749 * t1797 + t1750 * t1799;
t1647 = t1709 * t1708;
t1788 = t1798 * qJDD(1);
t1790 = t1800 * qJDD(1);
t1767 = t1788 * t1804 - t1808 * t1790;
t1885 = t1771 * qJD(3);
t1738 = -t1767 - t1885;
t1729 = qJDD(4) - t1738;
t1941 = -t1647 + t1729;
t1951 = t1799 * t1941;
t1526 = t1625 * t1797 + t1951;
t1795 = t1800 ^ 2;
t1921 = t1798 * g(3);
t1929 = 2 * qJD(2);
t1805 = sin(qJ(1));
t1809 = cos(qJ(1));
t1782 = t1809 * g(1) + t1805 * g(2);
t1935 = (pkin(7) + qJ(2)) * qJDD(1) - t1782;
t1722 = -t1921 + t1935 * t1800 + (t1800 * t1929 + (-pkin(1) * t1800 - pkin(2) * t1795) * qJD(1)) * qJD(1);
t1872 = pkin(2) * t1800 + pkin(1);
t1920 = t1800 * g(3);
t1812 = -t1920 + ((qJD(1) * t1872 - (2 * qJD(2))) * qJD(1) - t1935) * t1798;
t1659 = t1808 * t1722 + t1804 * t1812;
t1730 = pkin(3) * t1769 - pkin(8) * t1771;
t1933 = qJD(3) ^ 2;
t1603 = -pkin(3) * t1933 + qJDD(3) * pkin(8) - t1730 * t1769 + t1659;
t1760 = t1769 * qJD(3);
t1781 = t1805 * g(1) - t1809 * g(2);
t1840 = -qJDD(2) + t1781;
t1794 = t1798 ^ 2;
t1874 = t1794 + t1795;
t1934 = qJD(1) ^ 2;
t1944 = -(pkin(7) * t1874 + qJ(2)) * t1934 - t1840;
t1613 = 0.2e1 * t1760 * pkin(8) + (-t1738 + t1885) * pkin(3) + (-pkin(8) * t1836 - t1872) * qJDD(1) + t1944;
t1524 = t1803 * t1603 - t1807 * t1613;
t1768 = t1836 * qJDD(1);
t1819 = t1760 - t1768;
t1815 = -t1803 * qJDD(3) + t1807 * t1819;
t1686 = -t1749 * qJD(4) - t1815;
t1723 = t1762 * t1749;
t1654 = t1723 + t1686;
t1713 = t1750 * t1749;
t1939 = -t1713 + t1729;
t1472 = pkin(4) * t1939 - qJ(5) * t1654 - t1524;
t1525 = t1807 * t1603 + t1803 * t1613;
t1716 = pkin(4) * t1762 - qJ(5) * t1750;
t1814 = t1807 * qJDD(3) + t1803 * t1819;
t1813 = t1750 * qJD(4) - t1814;
t1931 = t1749 ^ 2;
t1480 = -pkin(4) * t1931 - qJ(5) * t1813 - t1762 * t1716 + t1525;
t1838 = -0.2e1 * qJD(5) * t1709 + t1799 * t1472 - t1797 * t1480;
t1953 = pkin(4) * t1526 + t1838;
t1952 = t1797 * t1941;
t1802 = sin(qJ(6));
t1806 = cos(qJ(6));
t1638 = t1806 * t1708 + t1709 * t1802;
t1640 = -t1708 * t1802 + t1709 * t1806;
t1551 = t1640 * t1638;
t1728 = qJDD(6) + t1729;
t1943 = -t1551 + t1728;
t1950 = t1802 * t1943;
t1949 = t1803 * t1939;
t1742 = t1771 * t1769;
t1938 = -t1742 + qJDD(3);
t1948 = t1804 * t1938;
t1947 = t1806 * t1943;
t1946 = t1807 * t1939;
t1945 = t1808 * t1938;
t1429 = t1803 * t1524 + t1807 * t1525;
t1610 = t1799 * t1686 - t1797 * t1813;
t1854 = t1686 * t1797 + t1799 * t1813;
t1499 = -t1638 * qJD(6) + t1806 * t1610 - t1802 * t1854;
t1757 = qJD(6) + t1762;
t1616 = t1757 * t1638;
t1942 = -t1616 + t1499;
t1681 = t1762 * t1708;
t1940 = -t1681 + t1610;
t1573 = t1681 + t1610;
t1936 = t1934 * t1874;
t1855 = t1802 * t1610 + t1806 * t1854;
t1453 = (qJD(6) - t1757) * t1640 + t1855;
t1636 = t1638 ^ 2;
t1637 = t1640 ^ 2;
t1707 = t1709 ^ 2;
t1748 = t1750 ^ 2;
t1756 = t1757 ^ 2;
t1765 = t1769 ^ 2;
t1766 = t1771 ^ 2;
t1658 = t1722 * t1804 - t1808 * t1812;
t1575 = -t1658 * t1808 + t1659 * t1804;
t1927 = pkin(2) * t1575;
t1676 = -t1767 * t1804 - t1768 * t1808;
t1926 = pkin(2) * t1676;
t1925 = pkin(3) * t1804;
t1918 = qJD(5) * t1708;
t1699 = -0.2e1 * t1918;
t1876 = t1797 * t1472 + t1799 * t1480;
t1386 = t1699 + t1876;
t1320 = t1386 * t1797 + t1799 * t1838;
t1924 = pkin(4) * t1320;
t1892 = t1762 * t1709;
t1570 = t1854 - t1892;
t1487 = -t1570 * t1797 - t1573 * t1799;
t1923 = pkin(4) * t1487;
t1919 = qJDD(1) * pkin(1);
t1351 = pkin(5) * t1941 - pkin(9) * t1573 + t1838;
t1670 = pkin(5) * t1762 - pkin(9) * t1709;
t1354 = -pkin(5) * t1932 - pkin(9) * t1854 - t1762 * t1670 + t1386;
t1286 = -t1806 * t1351 + t1354 * t1802;
t1287 = t1351 * t1802 + t1354 * t1806;
t1236 = -t1286 * t1806 + t1287 * t1802;
t1916 = t1236 * t1797;
t1915 = t1236 * t1799;
t1914 = t1320 * t1803;
t1913 = t1320 * t1807;
t1602 = -qJDD(3) * pkin(3) - t1933 * pkin(8) + t1730 * t1771 + t1658;
t1515 = t1813 * pkin(4) - t1931 * qJ(5) + t1716 * t1750 + qJDD(5) + t1602;
t1422 = pkin(5) * t1854 - pkin(9) * t1932 + t1670 * t1709 + t1515;
t1912 = t1422 * t1802;
t1911 = t1422 * t1806;
t1910 = t1515 * t1797;
t1909 = t1515 * t1799;
t1534 = t1551 + t1728;
t1908 = t1534 * t1802;
t1907 = t1534 * t1806;
t1906 = t1575 * t1798;
t1905 = t1575 * t1800;
t1606 = t1647 + t1729;
t1904 = t1606 * t1797;
t1903 = t1606 * t1799;
t1668 = t1713 + t1729;
t1902 = t1668 * t1803;
t1901 = t1668 * t1807;
t1900 = t1729 * t1804;
t1733 = t1872 * qJDD(1) - t1944;
t1899 = t1733 * t1804;
t1898 = t1733 * t1808;
t1734 = qJDD(3) + t1742;
t1897 = t1734 * t1804;
t1896 = t1734 * t1808;
t1895 = t1757 * t1640;
t1894 = t1757 * t1802;
t1893 = t1757 * t1806;
t1891 = t1762 * t1797;
t1890 = t1762 * t1799;
t1889 = t1762 * t1803;
t1888 = t1762 * t1807;
t1887 = t1769 * t1804;
t1886 = t1769 * t1808;
t1884 = t1771 * t1804;
t1883 = t1771 * t1808;
t1882 = t1794 * t1934;
t1763 = qJ(2) * t1934 + t1840 + t1919;
t1881 = t1798 * t1763;
t1880 = t1798 * t1800;
t1879 = t1800 * t1763;
t1593 = t1803 * t1602;
t1878 = t1805 * t1763;
t1594 = t1807 * t1602;
t1875 = -pkin(3) * t1602 + pkin(8) * t1429;
t1873 = t1805 * qJDD(1);
t1871 = -pkin(3) * t1808 - pkin(2);
t1870 = t1804 * t1551;
t1869 = t1808 * t1551;
t1868 = t1804 * t1647;
t1867 = t1808 * t1647;
t1866 = t1804 * t1713;
t1865 = t1808 * t1713;
t1864 = t1809 * t1742;
t1863 = t1805 * t1742;
t1697 = -t1748 - t1930;
t1601 = -t1697 * t1803 - t1901;
t1655 = (qJD(4) + t1762) * t1749 + t1815;
t1862 = pkin(3) * t1655 + pkin(8) * t1601 + t1593;
t1690 = -t1930 - t1931;
t1587 = t1690 * t1807 - t1949;
t1724 = t1762 * t1750;
t1651 = -t1724 - t1813;
t1861 = pkin(3) * t1651 + pkin(8) * t1587 - t1594;
t1860 = t1798 * t1790;
t1237 = t1286 * t1802 + t1806 * t1287;
t1208 = t1237 * t1797 + t1915;
t1235 = pkin(5) * t1236;
t1859 = pkin(4) * t1208 + t1235;
t1456 = t1616 + t1499;
t1374 = -t1453 * t1802 - t1456 * t1806;
t1376 = -t1453 * t1806 + t1456 * t1802;
t1312 = t1374 * t1799 + t1376 * t1797;
t1372 = pkin(5) * t1374;
t1858 = pkin(4) * t1312 + t1372;
t1857 = t1763 + t1919;
t1321 = t1799 * t1386 - t1797 * t1838;
t1576 = t1658 * t1804 + t1808 * t1659;
t1816 = qJDD(1) * qJ(2) + (-pkin(1) * qJD(1) + t1929) * qJD(1) - t1782;
t1743 = t1798 * t1816 + t1920;
t1744 = t1800 * t1816 - t1921;
t1689 = t1743 * t1798 + t1800 * t1744;
t1853 = -t1781 * t1805 - t1809 * t1782;
t1209 = t1237 * t1799 - t1916;
t1230 = -pkin(5) * t1422 + pkin(9) * t1237;
t1183 = -pkin(4) * t1422 - pkin(9) * t1916 + qJ(5) * t1209 + t1230 * t1799;
t1187 = -pkin(9) * t1915 - qJ(5) * t1208 - t1230 * t1797;
t1192 = -t1208 * t1803 + t1209 * t1807;
t1852 = -pkin(3) * t1422 + pkin(8) * t1192 + t1807 * t1183 + t1803 * t1187;
t1505 = -t1636 - t1637;
t1224 = -pkin(5) * t1505 + pkin(9) * t1376 + t1237;
t1226 = -pkin(9) * t1374 - t1236;
t1314 = -t1374 * t1797 + t1376 * t1799;
t1202 = -pkin(4) * t1505 + qJ(5) * t1314 + t1224 * t1799 + t1226 * t1797;
t1204 = -qJ(5) * t1312 - t1224 * t1797 + t1226 * t1799;
t1262 = -t1312 * t1803 + t1314 * t1807;
t1851 = -pkin(3) * t1505 + pkin(8) * t1262 + t1807 * t1202 + t1803 * t1204;
t1452 = (qJD(6) + t1757) * t1640 + t1855;
t1543 = -t1756 - t1636;
t1462 = t1543 * t1806 - t1950;
t1333 = -pkin(5) * t1452 + pkin(9) * t1462 - t1911;
t1461 = t1543 * t1802 + t1947;
t1358 = -pkin(9) * t1461 + t1912;
t1379 = -t1461 * t1797 + t1462 * t1799;
t1253 = -pkin(4) * t1452 + qJ(5) * t1379 + t1333 * t1799 + t1358 * t1797;
t1378 = t1461 * t1799 + t1462 * t1797;
t1273 = -qJ(5) * t1378 - t1333 * t1797 + t1358 * t1799;
t1317 = -t1378 * t1803 + t1379 * t1807;
t1850 = -pkin(3) * t1452 + pkin(8) * t1317 + t1807 * t1253 + t1803 * t1273;
t1588 = -t1637 - t1756;
t1492 = -t1588 * t1802 - t1907;
t1338 = -pkin(5) * t1942 + pkin(9) * t1492 + t1912;
t1491 = t1588 * t1806 - t1908;
t1377 = -pkin(9) * t1491 + t1911;
t1395 = -t1491 * t1797 + t1492 * t1799;
t1268 = -pkin(4) * t1942 + qJ(5) * t1395 + t1338 * t1799 + t1377 * t1797;
t1394 = t1491 * t1799 + t1492 * t1797;
t1276 = -qJ(5) * t1394 - t1338 * t1797 + t1377 * t1799;
t1331 = -t1394 * t1803 + t1395 * t1807;
t1849 = -pkin(3) * t1942 + pkin(8) * t1331 + t1807 * t1268 + t1803 * t1276;
t1489 = -t1570 * t1799 + t1573 * t1797;
t1582 = -t1707 - t1932;
t1289 = -pkin(4) * t1582 + qJ(5) * t1489 + t1321;
t1296 = -qJ(5) * t1487 - t1320;
t1393 = -t1487 * t1803 + t1489 * t1807;
t1848 = -pkin(3) * t1582 + pkin(8) * t1393 + t1807 * t1289 + t1803 * t1296;
t1527 = t1625 * t1799 - t1952;
t1569 = t1854 + t1892;
t1404 = -pkin(4) * t1569 + qJ(5) * t1527 - t1909;
t1426 = -qJ(5) * t1526 + t1910;
t1433 = -t1526 * t1803 + t1527 * t1807;
t1847 = -pkin(3) * t1569 + pkin(8) * t1433 + t1807 * t1404 + t1803 * t1426;
t1661 = -t1707 - t1930;
t1549 = -t1661 * t1797 - t1903;
t1409 = -pkin(4) * t1940 + qJ(5) * t1549 + t1910;
t1548 = t1661 * t1799 - t1904;
t1443 = -qJ(5) * t1548 + t1909;
t1469 = -t1548 * t1803 + t1549 * t1807;
t1846 = -pkin(3) * t1940 + pkin(8) * t1469 + t1807 * t1409 + t1803 * t1443;
t1652 = (-qJD(4) + t1762) * t1750 + t1814;
t1562 = t1652 * t1807 + t1654 * t1803;
t1680 = t1748 + t1931;
t1845 = pkin(3) * t1680 + pkin(8) * t1562 + t1429;
t1410 = t1429 * t1804 - t1602 * t1808;
t1844 = pkin(2) * t1410 + t1875;
t1754 = -t1766 - t1933;
t1693 = t1754 * t1808 - t1897;
t1843 = pkin(2) * t1693 - t1659;
t1842 = pkin(4) * t1548 - t1876;
t1780 = qJDD(1) * t1809 - t1805 * t1934;
t1841 = -pkin(6) * t1780 - g(3) * t1805;
t1839 = pkin(5) * t1461 - t1286;
t1428 = -t1524 * t1807 + t1525 * t1803;
t1688 = t1743 * t1800 - t1744 * t1798;
t1837 = t1781 * t1809 - t1805 * t1782;
t1779 = t1809 * t1934 + t1873;
t1536 = t1587 * t1804 + t1651 * t1808;
t1835 = pkin(2) * t1536 + t1861;
t1541 = t1601 * t1804 + t1655 * t1808;
t1834 = pkin(2) * t1541 + t1862;
t1731 = -t1765 - t1933;
t1671 = t1731 * t1804 + t1945;
t1833 = pkin(2) * t1671 - t1658;
t1832 = pkin(5) * t1491 - t1287;
t1773 = t1800 * t1936;
t1831 = -t1805 * t1773 + t1790 * t1809;
t1830 = t1773 * t1809 + t1800 * t1873;
t1184 = t1192 * t1804 - t1422 * t1808;
t1829 = pkin(2) * t1184 + t1852;
t1244 = t1262 * t1804 - t1505 * t1808;
t1828 = pkin(2) * t1244 + t1851;
t1283 = t1317 * t1804 - t1452 * t1808;
t1827 = pkin(2) * t1283 + t1850;
t1301 = t1331 * t1804 - t1808 * t1942;
t1826 = pkin(2) * t1301 + t1849;
t1380 = t1393 * t1804 - t1582 * t1808;
t1825 = pkin(2) * t1380 + t1848;
t1406 = t1433 * t1804 - t1569 * t1808;
t1824 = pkin(2) * t1406 + t1847;
t1413 = t1469 * t1804 - t1808 * t1940;
t1823 = pkin(2) * t1413 + t1846;
t1519 = t1562 * t1804 + t1680 * t1808;
t1822 = pkin(2) * t1519 + t1845;
t1271 = t1321 * t1807 - t1914;
t1300 = -pkin(4) * t1515 + qJ(5) * t1321;
t1821 = -pkin(3) * t1515 + pkin(8) * t1271 - qJ(5) * t1914 + t1807 * t1300;
t1820 = pkin(4) * t1378 + t1839;
t1818 = pkin(4) * t1394 + t1832;
t1263 = t1271 * t1804 - t1515 * t1808;
t1817 = pkin(2) * t1263 + t1821;
t1791 = t1795 * t1934;
t1789 = t1795 * qJDD(1);
t1787 = t1794 * qJDD(1);
t1784 = t1934 * t1880;
t1783 = 0.2e1 * t1860;
t1778 = -t1791 + t1882;
t1777 = t1791 + t1882;
t1776 = t1789 - t1787;
t1775 = t1789 + t1787;
t1772 = t1798 * t1936;
t1764 = -pkin(6) * t1779 + g(3) * t1809;
t1755 = -t1766 + t1933;
t1753 = t1765 - t1933;
t1752 = t1780 * t1880;
t1751 = t1779 * t1880;
t1747 = t1772 * t1809 + t1798 * t1873;
t1746 = t1805 * t1772 - t1788 * t1809;
t1741 = t1766 - t1765;
t1739 = -0.2e1 * t1760 + t1768;
t1737 = t1767 + 0.2e1 * t1885;
t1727 = t1808 * t1729;
t1726 = (t1884 - t1886) * qJD(3);
t1725 = (-t1883 - t1887) * qJD(3);
t1721 = -qJ(2) * t1773 + t1800 * t1857;
t1720 = qJ(2) * t1772 - t1798 * t1857;
t1719 = -t1748 + t1930;
t1718 = -t1930 + t1931;
t1712 = -t1765 - t1766;
t1711 = t1748 - t1931;
t1705 = -qJD(3) * t1884 - t1808 * t1819;
t1704 = qJD(3) * t1883 - t1804 * t1819;
t1703 = qJD(3) * t1886 - t1738 * t1804;
t1702 = qJD(3) * t1887 + t1738 * t1808;
t1696 = -t1754 * t1804 - t1896;
t1695 = -t1755 * t1804 + t1945;
t1694 = t1753 * t1808 - t1897;
t1692 = t1755 * t1808 + t1948;
t1691 = t1753 * t1804 + t1896;
t1679 = -t1737 * t1808 - t1739 * t1804;
t1678 = -t1767 * t1808 + t1768 * t1804;
t1677 = -t1737 * t1804 + t1739 * t1808;
t1675 = -t1707 + t1930;
t1674 = -t1930 + t1932;
t1672 = t1731 * t1808 - t1948;
t1666 = pkin(1) * t1763 + qJ(2) * t1689;
t1665 = (-t1749 * t1807 + t1750 * t1803) * t1762;
t1664 = (-t1749 * t1803 - t1750 * t1807) * t1762;
t1663 = -t1725 * t1798 + t1726 * t1800;
t1662 = t1725 * t1800 + t1726 * t1798;
t1660 = pkin(1) * t1777 + qJ(2) * t1775 + t1689;
t1656 = -pkin(7) * t1693 - t1898;
t1653 = -t1723 + t1686;
t1650 = -t1724 + t1813;
t1646 = t1707 - t1932;
t1645 = -pkin(7) * t1671 - t1899;
t1644 = t1686 * t1807 - t1750 * t1889;
t1643 = t1686 * t1803 + t1750 * t1888;
t1642 = t1749 * t1888 + t1803 * t1813;
t1641 = -t1749 * t1889 + t1807 * t1813;
t1635 = -t1704 * t1798 + t1705 * t1800;
t1634 = -t1702 * t1798 + t1703 * t1800;
t1633 = t1704 * t1800 + t1705 * t1798;
t1632 = t1702 * t1800 + t1703 * t1798;
t1631 = -t1693 * t1798 + t1696 * t1800;
t1630 = -t1692 * t1798 + t1695 * t1800;
t1629 = -t1691 * t1798 + t1694 * t1800;
t1628 = t1693 * t1800 + t1696 * t1798;
t1627 = t1692 * t1800 + t1695 * t1798;
t1626 = t1691 * t1800 + t1694 * t1798;
t1624 = t1665 * t1808 + t1900;
t1623 = t1665 * t1804 - t1727;
t1622 = t1718 * t1807 - t1902;
t1621 = -t1719 * t1803 + t1946;
t1620 = t1718 * t1803 + t1901;
t1619 = t1719 * t1807 + t1949;
t1618 = -pkin(2) * t1739 + pkin(7) * t1696 - t1899;
t1615 = -t1637 + t1756;
t1614 = t1636 - t1756;
t1604 = -pkin(2) * t1737 + pkin(7) * t1672 + t1898;
t1600 = t1697 * t1807 - t1902;
t1597 = (-t1708 * t1799 + t1709 * t1797) * t1762;
t1596 = (-t1708 * t1797 - t1709 * t1799) * t1762;
t1592 = -t1677 * t1798 + t1679 * t1800;
t1591 = -t1676 * t1798 + t1678 * t1800;
t1590 = t1677 * t1800 + t1679 * t1798;
t1589 = t1676 * t1800 + t1678 * t1798;
t1586 = t1690 * t1803 + t1946;
t1585 = -t1671 * t1798 + t1672 * t1800;
t1584 = t1671 * t1800 + t1672 * t1798;
t1580 = t1644 * t1808 + t1866;
t1579 = t1642 * t1808 - t1866;
t1578 = t1644 * t1804 - t1865;
t1577 = t1642 * t1804 + t1865;
t1566 = t1610 * t1799 - t1709 * t1891;
t1565 = t1610 * t1797 + t1709 * t1890;
t1564 = t1708 * t1890 + t1797 * t1854;
t1563 = t1708 * t1891 - t1799 * t1854;
t1561 = t1651 * t1807 - t1653 * t1803;
t1560 = t1652 * t1803 - t1654 * t1807;
t1559 = t1651 * t1803 + t1653 * t1807;
t1557 = t1674 * t1799 - t1904;
t1556 = -t1675 * t1797 + t1951;
t1555 = t1674 * t1797 + t1903;
t1554 = t1675 * t1799 + t1952;
t1553 = pkin(2) * t1733 + pkin(7) * t1576;
t1552 = -pkin(1) * t1589 - t1926;
t1550 = t1637 - t1636;
t1547 = t1622 * t1808 - t1650 * t1804;
t1546 = t1621 * t1808 + t1654 * t1804;
t1545 = t1622 * t1804 + t1650 * t1808;
t1544 = t1621 * t1804 - t1654 * t1808;
t1542 = t1601 * t1808 - t1655 * t1804;
t1540 = -pkin(7) * t1676 - t1575;
t1539 = -t1623 * t1798 + t1624 * t1800;
t1538 = t1623 * t1800 + t1624 * t1798;
t1537 = t1587 * t1808 - t1651 * t1804;
t1532 = (-t1638 * t1806 + t1640 * t1802) * t1757;
t1531 = (-t1638 * t1802 - t1640 * t1806) * t1757;
t1530 = t1561 * t1808 + t1711 * t1804;
t1529 = t1561 * t1804 - t1711 * t1808;
t1528 = -pkin(1) * t1628 - t1843;
t1521 = -pkin(8) * t1600 + t1594;
t1520 = t1562 * t1808 - t1680 * t1804;
t1518 = -pkin(8) * t1586 + t1593;
t1517 = -t1596 * t1803 + t1597 * t1807;
t1516 = t1596 * t1807 + t1597 * t1803;
t1513 = -pkin(2) * t1712 + pkin(7) * t1678 + t1576;
t1512 = -pkin(1) * t1584 - t1833;
t1511 = -t1578 * t1798 + t1580 * t1800;
t1510 = -t1577 * t1798 + t1579 * t1800;
t1509 = t1578 * t1800 + t1580 * t1798;
t1508 = t1577 * t1800 + t1579 * t1798;
t1507 = t1517 * t1808 + t1900;
t1506 = t1517 * t1804 - t1727;
t1503 = -qJ(2) * t1628 - t1618 * t1798 + t1656 * t1800;
t1502 = t1576 * t1800 - t1906;
t1501 = t1576 * t1798 + t1905;
t1500 = -pkin(1) * t1739 + qJ(2) * t1631 + t1618 * t1800 + t1656 * t1798;
t1498 = -qJD(6) * t1640 - t1855;
t1497 = t1614 * t1806 - t1908;
t1496 = -t1615 * t1802 + t1947;
t1495 = t1614 * t1802 + t1907;
t1494 = t1615 * t1806 + t1950;
t1493 = -qJ(2) * t1584 - t1604 * t1798 + t1645 * t1800;
t1488 = -t1569 * t1799 - t1797 * t1940;
t1486 = -t1569 * t1797 + t1799 * t1940;
t1485 = -pkin(3) * t1600 + t1525;
t1484 = -t1565 * t1803 + t1566 * t1807;
t1483 = -t1563 * t1803 + t1564 * t1807;
t1482 = t1565 * t1807 + t1566 * t1803;
t1481 = t1563 * t1807 + t1564 * t1803;
t1479 = -pkin(3) * t1586 + t1524;
t1477 = -t1555 * t1803 + t1557 * t1807;
t1476 = -t1554 * t1803 + t1556 * t1807;
t1475 = t1555 * t1807 + t1557 * t1803;
t1474 = t1554 * t1807 + t1556 * t1803;
t1473 = -pkin(1) * t1737 + qJ(2) * t1585 + t1604 * t1800 + t1645 * t1798;
t1468 = t1548 * t1807 + t1549 * t1803;
t1466 = -t1545 * t1798 + t1547 * t1800;
t1465 = -t1544 * t1798 + t1546 * t1800;
t1464 = t1545 * t1800 + t1547 * t1798;
t1463 = t1544 * t1800 + t1546 * t1798;
t1459 = -t1541 * t1798 + t1542 * t1800;
t1458 = t1541 * t1800 + t1542 * t1798;
t1449 = -t1536 * t1798 + t1537 * t1800;
t1448 = t1536 * t1800 + t1537 * t1798;
t1447 = t1499 * t1806 - t1640 * t1894;
t1446 = t1499 * t1802 + t1640 * t1893;
t1445 = -t1498 * t1802 + t1638 * t1893;
t1444 = t1498 * t1806 + t1638 * t1894;
t1442 = -t1531 * t1797 + t1532 * t1799;
t1441 = t1531 * t1799 + t1532 * t1797;
t1439 = t1484 * t1808 + t1868;
t1438 = t1483 * t1808 - t1868;
t1437 = t1484 * t1804 - t1867;
t1436 = t1483 * t1804 + t1867;
t1435 = -t1529 * t1798 + t1530 * t1800;
t1434 = t1529 * t1800 + t1530 * t1798;
t1432 = t1526 * t1807 + t1527 * t1803;
t1430 = -pkin(1) * t1501 - t1927;
t1424 = -t1519 * t1798 + t1520 * t1800;
t1423 = t1519 * t1800 + t1520 * t1798;
t1420 = t1477 * t1808 - t1570 * t1804;
t1419 = t1476 * t1808 + t1573 * t1804;
t1418 = t1477 * t1804 + t1570 * t1808;
t1417 = t1476 * t1804 - t1573 * t1808;
t1416 = -t1506 * t1798 + t1507 * t1800;
t1415 = t1506 * t1800 + t1507 * t1798;
t1414 = t1469 * t1808 + t1804 * t1940;
t1412 = -qJ(2) * t1589 - t1513 * t1798 + t1540 * t1800;
t1411 = t1429 * t1808 + t1602 * t1804;
t1407 = t1433 * t1808 + t1569 * t1804;
t1405 = -pkin(1) * t1712 + qJ(2) * t1591 + t1513 * t1800 + t1540 * t1798;
t1403 = -pkin(7) * t1905 - qJ(2) * t1501 - t1553 * t1798;
t1401 = -pkin(8) * t1560 - t1428;
t1400 = -t1495 * t1797 + t1497 * t1799;
t1399 = -t1494 * t1797 + t1496 * t1799;
t1398 = t1495 * t1799 + t1497 * t1797;
t1397 = t1494 * t1799 + t1496 * t1797;
t1396 = pkin(1) * t1733 - pkin(7) * t1906 + qJ(2) * t1502 + t1553 * t1800;
t1392 = -t1486 * t1803 + t1488 * t1807;
t1391 = t1487 * t1807 + t1489 * t1803;
t1390 = t1486 * t1807 + t1488 * t1803;
t1388 = t1392 * t1808 + t1646 * t1804;
t1387 = t1392 * t1804 - t1646 * t1808;
t1383 = -pkin(7) * t1541 - t1485 * t1804 + t1521 * t1808;
t1382 = -pkin(7) * t1536 - t1479 * t1804 + t1518 * t1808;
t1381 = t1393 * t1808 + t1582 * t1804;
t1375 = -t1452 * t1806 - t1802 * t1942;
t1373 = -t1452 * t1802 + t1806 * t1942;
t1371 = -t1446 * t1797 + t1447 * t1799;
t1370 = -t1444 * t1797 + t1445 * t1799;
t1369 = t1446 * t1799 + t1447 * t1797;
t1368 = t1444 * t1799 + t1445 * t1797;
t1367 = -t1441 * t1803 + t1442 * t1807;
t1366 = t1441 * t1807 + t1442 * t1803;
t1365 = -t1437 * t1798 + t1439 * t1800;
t1364 = -t1436 * t1798 + t1438 * t1800;
t1363 = t1437 * t1800 + t1439 * t1798;
t1362 = t1436 * t1800 + t1438 * t1798;
t1361 = -pkin(2) * t1600 + pkin(7) * t1542 + t1485 * t1808 + t1521 * t1804;
t1360 = t1367 * t1808 + t1728 * t1804;
t1359 = t1367 * t1804 - t1728 * t1808;
t1357 = -pkin(2) * t1586 + pkin(7) * t1537 + t1479 * t1808 + t1518 * t1804;
t1356 = -pkin(1) * t1458 - t1834;
t1355 = -pkin(1) * t1448 - t1835;
t1353 = -pkin(3) * t1391 - t1923;
t1352 = -pkin(7) * t1519 + t1401 * t1808 + t1560 * t1925;
t1350 = -t1418 * t1798 + t1420 * t1800;
t1349 = -t1417 * t1798 + t1419 * t1800;
t1348 = t1418 * t1800 + t1420 * t1798;
t1347 = t1417 * t1800 + t1419 * t1798;
t1345 = -t1413 * t1798 + t1414 * t1800;
t1344 = t1413 * t1800 + t1414 * t1798;
t1343 = -t1410 * t1798 + t1411 * t1800;
t1342 = t1410 * t1800 + t1411 * t1798;
t1341 = pkin(7) * t1520 + t1804 * t1401 + t1560 * t1871;
t1340 = -t1406 * t1798 + t1407 * t1800;
t1339 = t1406 * t1800 + t1407 * t1798;
t1337 = -t1398 * t1803 + t1400 * t1807;
t1336 = -t1397 * t1803 + t1399 * t1807;
t1335 = t1398 * t1807 + t1400 * t1803;
t1334 = t1397 * t1807 + t1399 * t1803;
t1332 = -pkin(3) * t1468 + t1699 - t1842;
t1330 = t1394 * t1807 + t1395 * t1803;
t1328 = -pkin(1) * t1423 - t1822;
t1327 = -pkin(8) * t1468 - t1409 * t1803 + t1443 * t1807;
t1326 = -pkin(3) * t1432 - t1953;
t1325 = -pkin(7) * t1410 + (-pkin(8) * t1808 + t1925) * t1428;
t1324 = -pkin(8) * t1432 - t1404 * t1803 + t1426 * t1807;
t1323 = -t1387 * t1798 + t1388 * t1800;
t1322 = t1387 * t1800 + t1388 * t1798;
t1319 = -t1380 * t1798 + t1381 * t1800;
t1318 = t1380 * t1800 + t1381 * t1798;
t1316 = t1378 * t1807 + t1379 * t1803;
t1313 = -t1373 * t1797 + t1375 * t1799;
t1311 = t1373 * t1799 + t1375 * t1797;
t1310 = -t1369 * t1803 + t1371 * t1807;
t1309 = -t1368 * t1803 + t1370 * t1807;
t1308 = t1369 * t1807 + t1371 * t1803;
t1307 = t1368 * t1807 + t1370 * t1803;
t1306 = t1337 * t1808 - t1453 * t1804;
t1305 = t1336 * t1808 + t1456 * t1804;
t1304 = t1337 * t1804 + t1453 * t1808;
t1303 = t1336 * t1804 - t1456 * t1808;
t1302 = t1331 * t1808 + t1804 * t1942;
t1298 = -t1359 * t1798 + t1360 * t1800;
t1297 = t1359 * t1800 + t1360 * t1798;
t1294 = t1310 * t1808 + t1870;
t1293 = t1309 * t1808 - t1870;
t1292 = t1310 * t1804 - t1869;
t1291 = t1309 * t1804 + t1869;
t1290 = pkin(7) * t1411 + (-pkin(8) * t1804 + t1871) * t1428;
t1284 = t1317 * t1808 + t1452 * t1804;
t1282 = -qJ(2) * t1458 - t1361 * t1798 + t1383 * t1800;
t1281 = -qJ(2) * t1448 - t1357 * t1798 + t1382 * t1800;
t1280 = -pkin(1) * t1600 + qJ(2) * t1459 + t1361 * t1800 + t1383 * t1798;
t1279 = -pkin(1) * t1586 + qJ(2) * t1449 + t1357 * t1800 + t1382 * t1798;
t1278 = -pkin(1) * t1342 - t1844;
t1277 = -qJ(2) * t1423 - t1341 * t1798 + t1352 * t1800;
t1274 = -pkin(1) * t1560 + qJ(2) * t1424 + t1341 * t1800 + t1352 * t1798;
t1270 = t1321 * t1803 + t1913;
t1266 = -pkin(1) * t1344 - t1823;
t1265 = -pkin(7) * t1413 + t1327 * t1808 - t1332 * t1804;
t1264 = t1271 * t1808 + t1515 * t1804;
t1261 = -t1311 * t1803 + t1313 * t1807;
t1260 = t1312 * t1807 + t1314 * t1803;
t1259 = t1311 * t1807 + t1313 * t1803;
t1257 = -t1304 * t1798 + t1306 * t1800;
t1256 = -t1303 * t1798 + t1305 * t1800;
t1255 = t1304 * t1800 + t1306 * t1798;
t1254 = t1303 * t1800 + t1305 * t1798;
t1251 = -pkin(1) * t1339 - t1824;
t1250 = -t1301 * t1798 + t1302 * t1800;
t1249 = t1301 * t1800 + t1302 * t1798;
t1248 = -pkin(7) * t1406 + t1324 * t1808 - t1326 * t1804;
t1247 = t1261 * t1808 + t1550 * t1804;
t1246 = t1261 * t1804 - t1550 * t1808;
t1245 = t1262 * t1808 + t1505 * t1804;
t1243 = -pkin(2) * t1468 + pkin(7) * t1414 + t1327 * t1804 + t1332 * t1808;
t1242 = -t1292 * t1798 + t1294 * t1800;
t1241 = -t1291 * t1798 + t1293 * t1800;
t1240 = t1292 * t1800 + t1294 * t1798;
t1239 = t1291 * t1800 + t1293 * t1798;
t1238 = -pkin(2) * t1432 + pkin(7) * t1407 + t1324 * t1804 + t1326 * t1808;
t1234 = -t1283 * t1798 + t1284 * t1800;
t1233 = t1283 * t1800 + t1284 * t1798;
t1232 = -pkin(3) * t1330 - t1818;
t1231 = -pkin(8) * t1391 - t1289 * t1803 + t1296 * t1807;
t1229 = -qJ(2) * t1342 - t1290 * t1798 + t1325 * t1800;
t1228 = -pkin(3) * t1270 - t1924;
t1227 = -pkin(3) * t1316 - t1820;
t1225 = -pkin(1) * t1428 + qJ(2) * t1343 + t1290 * t1800 + t1325 * t1798;
t1223 = -pkin(3) * t1260 - t1858;
t1222 = -pkin(7) * t1380 + t1231 * t1808 - t1353 * t1804;
t1221 = -pkin(8) * t1270 - qJ(5) * t1913 - t1300 * t1803;
t1220 = -t1263 * t1798 + t1264 * t1800;
t1219 = t1263 * t1800 + t1264 * t1798;
t1218 = -pkin(2) * t1391 + pkin(7) * t1381 + t1231 * t1804 + t1353 * t1808;
t1217 = -pkin(1) * t1318 - t1825;
t1216 = -t1246 * t1798 + t1247 * t1800;
t1215 = t1246 * t1800 + t1247 * t1798;
t1214 = -pkin(8) * t1330 - t1268 * t1803 + t1276 * t1807;
t1213 = -t1244 * t1798 + t1245 * t1800;
t1212 = t1244 * t1800 + t1245 * t1798;
t1211 = -qJ(2) * t1344 - t1243 * t1798 + t1265 * t1800;
t1210 = -pkin(8) * t1316 - t1253 * t1803 + t1273 * t1807;
t1207 = -pkin(1) * t1468 + qJ(2) * t1345 + t1243 * t1800 + t1265 * t1798;
t1206 = -qJ(2) * t1339 - t1238 * t1798 + t1248 * t1800;
t1205 = -pkin(1) * t1432 + qJ(2) * t1340 + t1238 * t1800 + t1248 * t1798;
t1200 = -pkin(1) * t1249 - t1826;
t1199 = -pkin(7) * t1301 + t1214 * t1808 - t1232 * t1804;
t1198 = -pkin(7) * t1263 + t1221 * t1808 - t1228 * t1804;
t1197 = -pkin(2) * t1330 + pkin(7) * t1302 + t1214 * t1804 + t1232 * t1808;
t1196 = -pkin(1) * t1233 - t1827;
t1195 = -qJ(2) * t1318 - t1218 * t1798 + t1222 * t1800;
t1194 = -pkin(7) * t1283 + t1210 * t1808 - t1227 * t1804;
t1193 = -pkin(1) * t1391 + qJ(2) * t1319 + t1218 * t1800 + t1222 * t1798;
t1191 = t1208 * t1807 + t1209 * t1803;
t1189 = -pkin(1) * t1219 - t1817;
t1188 = -pkin(2) * t1316 + pkin(7) * t1284 + t1210 * t1804 + t1227 * t1808;
t1185 = t1192 * t1808 + t1422 * t1804;
t1181 = -pkin(2) * t1270 + pkin(7) * t1264 + t1221 * t1804 + t1228 * t1808;
t1180 = -pkin(8) * t1260 - t1202 * t1803 + t1204 * t1807;
t1179 = -pkin(3) * t1191 - t1859;
t1178 = -qJ(2) * t1249 - t1197 * t1798 + t1199 * t1800;
t1177 = -pkin(1) * t1330 + qJ(2) * t1250 + t1197 * t1800 + t1199 * t1798;
t1176 = -t1184 * t1798 + t1185 * t1800;
t1175 = t1184 * t1800 + t1185 * t1798;
t1174 = -qJ(2) * t1233 - t1188 * t1798 + t1194 * t1800;
t1173 = -pkin(7) * t1244 + t1180 * t1808 - t1223 * t1804;
t1172 = -pkin(1) * t1316 + qJ(2) * t1234 + t1188 * t1800 + t1194 * t1798;
t1171 = -pkin(1) * t1212 - t1828;
t1170 = -qJ(2) * t1219 - t1181 * t1798 + t1198 * t1800;
t1169 = -pkin(2) * t1260 + pkin(7) * t1245 + t1180 * t1804 + t1223 * t1808;
t1168 = -pkin(1) * t1270 + qJ(2) * t1220 + t1181 * t1800 + t1198 * t1798;
t1167 = -pkin(8) * t1191 - t1183 * t1803 + t1187 * t1807;
t1166 = -qJ(2) * t1212 - t1169 * t1798 + t1173 * t1800;
t1165 = -pkin(1) * t1260 + qJ(2) * t1213 + t1169 * t1800 + t1173 * t1798;
t1164 = -pkin(7) * t1184 + t1167 * t1808 - t1179 * t1804;
t1163 = -pkin(1) * t1175 - t1829;
t1162 = -pkin(2) * t1191 + pkin(7) * t1185 + t1167 * t1804 + t1179 * t1808;
t1161 = -qJ(2) * t1175 - t1162 * t1798 + t1164 * t1800;
t1160 = -pkin(1) * t1191 + qJ(2) * t1176 + t1162 * t1800 + t1164 * t1798;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1780, 0, -t1779, 0, t1841, -t1764, -t1837, -pkin(6) * t1837, t1752, t1776 * t1809 + t1805 * t1778, t1747, -t1752, t1830, 0, -pkin(6) * t1831 - t1805 * t1743 - t1809 * t1881, -pkin(6) * t1746 - t1805 * t1744 - t1809 * t1879, t1809 * t1688 - pkin(6) * (t1805 * t1775 + t1777 * t1809), -pkin(6) * (t1805 * t1689 + t1763 * t1809) - (pkin(1) * t1805 - qJ(2) * t1809) * t1688, t1635 * t1809 + t1863, t1592 * t1809 + t1805 * t1741, t1630 * t1809 + t1805 * t1768, t1634 * t1809 - t1863, t1629 * t1809 - t1805 * t1767, t1805 * qJDD(3) + t1663 * t1809, t1809 * t1493 - t1805 * t1512 - pkin(6) * (t1805 * t1585 - t1737 * t1809), t1809 * t1503 - t1805 * t1528 - pkin(6) * (t1805 * t1631 - t1739 * t1809), t1809 * t1412 - t1805 * t1552 - pkin(6) * (t1805 * t1591 - t1712 * t1809), t1809 * t1403 - t1805 * t1430 - pkin(6) * (t1805 * t1502 + t1733 * t1809), t1511 * t1809 + t1805 * t1643, t1435 * t1809 + t1805 * t1559, t1465 * t1809 + t1805 * t1619, t1510 * t1809 - t1805 * t1641, t1466 * t1809 + t1805 * t1620, t1539 * t1809 + t1805 * t1664, t1809 * t1281 - t1805 * t1355 - pkin(6) * (t1805 * t1449 - t1586 * t1809), t1809 * t1282 - t1805 * t1356 - pkin(6) * (t1805 * t1459 - t1600 * t1809), t1809 * t1277 - t1805 * t1328 - pkin(6) * (t1805 * t1424 - t1560 * t1809), t1809 * t1229 - t1805 * t1278 - pkin(6) * (t1805 * t1343 - t1428 * t1809), t1365 * t1809 + t1805 * t1482, t1323 * t1809 + t1805 * t1390, t1349 * t1809 + t1805 * t1474, t1364 * t1809 + t1805 * t1481, t1350 * t1809 + t1805 * t1475, t1416 * t1809 + t1805 * t1516, t1809 * t1206 - t1805 * t1251 - pkin(6) * (t1805 * t1340 - t1432 * t1809), t1809 * t1211 - t1805 * t1266 - pkin(6) * (t1805 * t1345 - t1468 * t1809), t1809 * t1195 - t1805 * t1217 - pkin(6) * (t1805 * t1319 - t1391 * t1809), t1809 * t1170 - t1805 * t1189 - pkin(6) * (t1805 * t1220 - t1270 * t1809), t1242 * t1809 + t1805 * t1308, t1216 * t1809 + t1805 * t1259, t1256 * t1809 + t1805 * t1334, t1241 * t1809 + t1805 * t1307, t1257 * t1809 + t1805 * t1335, t1298 * t1809 + t1805 * t1366, t1809 * t1174 - t1805 * t1196 - pkin(6) * (t1805 * t1234 - t1316 * t1809), t1809 * t1178 - t1805 * t1200 - pkin(6) * (t1805 * t1250 - t1330 * t1809), t1809 * t1166 - t1805 * t1171 - pkin(6) * (t1805 * t1213 - t1260 * t1809), t1809 * t1161 - t1805 * t1163 - pkin(6) * (t1805 * t1176 - t1191 * t1809); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1779, 0, t1780, 0, t1764, t1841, t1853, pkin(6) * t1853, t1751, t1805 * t1776 - t1778 * t1809, t1746, -t1751, -t1831, 0, -pkin(6) * t1830 + t1809 * t1743 - t1798 * t1878, pkin(6) * t1747 + t1744 * t1809 - t1800 * t1878, t1805 * t1688 + pkin(6) * (t1775 * t1809 - t1805 * t1777), pkin(6) * (t1689 * t1809 - t1878) - (-pkin(1) * t1809 - qJ(2) * t1805) * t1688, t1805 * t1635 - t1864, t1805 * t1592 - t1741 * t1809, t1805 * t1630 - t1768 * t1809, t1805 * t1634 + t1864, t1805 * t1629 + t1767 * t1809, -qJDD(3) * t1809 + t1805 * t1663, t1805 * t1493 + t1809 * t1512 + pkin(6) * (t1585 * t1809 + t1805 * t1737), t1805 * t1503 + t1809 * t1528 + pkin(6) * (t1631 * t1809 + t1805 * t1739), t1805 * t1412 + t1809 * t1552 + pkin(6) * (t1591 * t1809 + t1805 * t1712), t1805 * t1403 + t1809 * t1430 + pkin(6) * (t1502 * t1809 - t1805 * t1733), t1805 * t1511 - t1643 * t1809, t1805 * t1435 - t1559 * t1809, t1805 * t1465 - t1619 * t1809, t1805 * t1510 + t1641 * t1809, t1805 * t1466 - t1620 * t1809, t1805 * t1539 - t1664 * t1809, t1805 * t1281 + t1809 * t1355 + pkin(6) * (t1449 * t1809 + t1805 * t1586), t1805 * t1282 + t1809 * t1356 + pkin(6) * (t1459 * t1809 + t1805 * t1600), t1805 * t1277 + t1809 * t1328 + pkin(6) * (t1424 * t1809 + t1805 * t1560), t1805 * t1229 + t1809 * t1278 + pkin(6) * (t1343 * t1809 + t1805 * t1428), t1805 * t1365 - t1482 * t1809, t1805 * t1323 - t1390 * t1809, t1805 * t1349 - t1474 * t1809, t1805 * t1364 - t1481 * t1809, t1805 * t1350 - t1475 * t1809, t1805 * t1416 - t1516 * t1809, t1805 * t1206 + t1809 * t1251 + pkin(6) * (t1340 * t1809 + t1805 * t1432), t1805 * t1211 + t1809 * t1266 + pkin(6) * (t1345 * t1809 + t1805 * t1468), t1805 * t1195 + t1809 * t1217 + pkin(6) * (t1319 * t1809 + t1805 * t1391), t1805 * t1170 + t1809 * t1189 + pkin(6) * (t1220 * t1809 + t1805 * t1270), t1805 * t1242 - t1308 * t1809, t1805 * t1216 - t1259 * t1809, t1805 * t1256 - t1334 * t1809, t1805 * t1241 - t1307 * t1809, t1805 * t1257 - t1335 * t1809, t1805 * t1298 - t1366 * t1809, t1805 * t1174 + t1809 * t1196 + pkin(6) * (t1234 * t1809 + t1805 * t1316), t1805 * t1178 + t1809 * t1200 + pkin(6) * (t1250 * t1809 + t1805 * t1330), t1805 * t1166 + t1809 * t1171 + pkin(6) * (t1213 * t1809 + t1805 * t1260), t1805 * t1161 + t1809 * t1163 + pkin(6) * (t1176 * t1809 + t1805 * t1191); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1781, t1782, 0, 0, t1787, t1783, 0, t1789, 0, 0, t1721, t1720, t1660, t1666, t1633, t1590, t1627, t1632, t1626, t1662, t1473, t1500, t1405, t1396, t1509, t1434, t1463, t1508, t1464, t1538, t1279, t1280, t1274, t1225, t1363, t1322, t1347, t1362, t1348, t1415, t1205, t1207, t1193, t1168, t1240, t1215, t1254, t1239, t1255, t1297, t1172, t1177, t1165, t1160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1934, 0, 0, -g(3), -t1781, 0, t1860, t1776, t1772, -t1860, t1773, 0, -t1881, -t1879, t1688, qJ(2) * t1688, t1635, t1592, t1630, t1634, t1629, t1663, t1493, t1503, t1412, t1403, t1511, t1435, t1465, t1510, t1466, t1539, t1281, t1282, t1277, t1229, t1365, t1323, t1349, t1364, t1350, t1416, t1206, t1211, t1195, t1170, t1242, t1216, t1256, t1241, t1257, t1298, t1174, t1178, t1166, t1161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1934, 0, qJDD(1), 0, g(3), 0, -t1782, 0, t1784, -t1778, -t1788, -t1784, -t1790, 0, t1743, t1744, 0, pkin(1) * t1688, -t1742, -t1741, -t1768, t1742, t1767, -qJDD(3), t1512, t1528, t1552, t1430, -t1643, -t1559, -t1619, t1641, -t1620, -t1664, t1355, t1356, t1328, t1278, -t1482, -t1390, -t1474, -t1481, -t1475, -t1516, t1251, t1266, t1217, t1189, -t1308, -t1259, -t1334, -t1307, -t1335, -t1366, t1196, t1200, t1171, t1163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1781, t1782, 0, 0, t1787, t1783, 0, t1789, 0, 0, t1721, t1720, t1660, t1666, t1633, t1590, t1627, t1632, t1626, t1662, t1473, t1500, t1405, t1396, t1509, t1434, t1463, t1508, t1464, t1538, t1279, t1280, t1274, t1225, t1363, t1322, t1347, t1362, t1348, t1415, t1205, t1207, t1193, t1168, t1240, t1215, t1254, t1239, t1255, t1297, t1172, t1177, t1165, t1160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1788, t1790, t1784, 0, t1791, 0, 0, -t1763, t1743, 0, t1705, t1679, t1695, t1703, t1694, t1726, t1645, t1656, t1540, -pkin(7) * t1575, t1580, t1530, t1546, t1579, t1547, t1624, t1382, t1383, t1352, t1325, t1439, t1388, t1419, t1438, t1420, t1507, t1248, t1265, t1222, t1198, t1294, t1247, t1305, t1293, t1306, t1360, t1194, t1199, t1173, t1164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1788, -t1882, t1790, -t1784, 0, t1763, 0, t1744, 0, t1704, t1677, t1692, t1702, t1691, t1725, t1604, t1618, t1513, t1553, t1578, t1529, t1544, t1577, t1545, t1623, t1357, t1361, t1341, t1290, t1437, t1387, t1417, t1436, t1418, t1506, t1238, t1243, t1218, t1181, t1292, t1246, t1303, t1291, t1304, t1359, t1188, t1197, t1169, t1162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1784, t1778, t1788, t1784, t1790, 0, -t1743, -t1744, 0, 0, t1742, t1741, t1768, -t1742, -t1767, qJDD(3), t1833, t1843, t1926, t1927, t1643, t1559, t1619, -t1641, t1620, t1664, t1835, t1834, t1822, t1844, t1482, t1390, t1474, t1481, t1475, t1516, t1824, t1823, t1825, t1817, t1308, t1259, t1334, t1307, t1335, t1366, t1827, t1826, t1828, t1829; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1819, -t1737, t1938, t1760, t1753, -t1760, 0, -t1733, t1658, 0, t1644, t1561, t1621, t1642, t1622, t1665, t1518, t1521, t1401, -pkin(8) * t1428, t1484, t1392, t1476, t1483, t1477, t1517, t1324, t1327, t1231, t1221, t1310, t1261, t1336, t1309, t1337, t1367, t1210, t1214, t1180, t1167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1885, t1739, t1755, t1738, t1734, -t1885, t1733, 0, t1659, 0, -t1713, -t1711, -t1654, t1713, t1650, -t1729, t1479, t1485, -pkin(3) * t1560, -pkin(3) * t1428, -t1647, -t1646, -t1573, t1647, t1570, -t1729, t1326, t1332, t1353, t1228, -t1551, -t1550, -t1456, t1551, t1453, -t1728, t1227, t1232, t1223, t1179; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1742, t1741, t1768, -t1742, -t1767, qJDD(3), -t1658, -t1659, 0, 0, t1643, t1559, t1619, -t1641, t1620, t1664, t1861, t1862, t1845, t1875, t1482, t1390, t1474, t1481, t1475, t1516, t1847, t1846, t1848, t1821, t1308, t1259, t1334, t1307, t1335, t1366, t1850, t1849, t1851, t1852; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1686, t1651, t1939, t1723, t1718, -t1723, 0, t1602, t1524, 0, t1566, t1488, t1556, t1564, t1557, t1597, t1426, t1443, t1296, -qJ(5) * t1320, t1371, t1313, t1399, t1370, t1400, t1442, t1273, t1276, t1204, t1187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1724, t1653, t1719, -t1813, t1668, -t1724, -t1602, 0, t1525, 0, t1565, t1486, t1554, t1563, t1555, t1596, t1404, t1409, t1289, t1300, t1369, t1311, t1397, t1368, t1398, t1441, t1253, t1268, t1202, t1183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1713, t1711, t1654, -t1713, -t1650, t1729, -t1524, -t1525, 0, 0, t1647, t1646, t1573, -t1647, -t1570, t1729, t1953, t1842 + 0.2e1 * t1918, t1923, t1924, t1551, t1550, t1456, -t1551, -t1453, t1728, t1820, t1818, t1858, t1859; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1610, -t1569, t1941, t1681, t1674, -t1681, 0, t1515, -t1838, 0, t1447, t1375, t1496, t1445, t1497, t1532, t1358, t1377, t1226, -pkin(9) * t1236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1892, t1940, t1675, -t1854, t1606, -t1892, -t1515, 0, t1386, 0, t1446, t1373, t1494, t1444, t1495, t1531, t1333, t1338, t1224, t1230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1647, t1646, t1573, -t1647, -t1570, t1729, t1838, -t1386, 0, 0, t1551, t1550, t1456, -t1551, -t1453, t1728, t1839, t1832, t1372, t1235; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1499, -t1452, t1943, t1616, t1614, -t1616, 0, t1422, t1286, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1895, t1942, t1615, t1498, t1534, -t1895, -t1422, 0, t1287, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1551, t1550, t1456, -t1551, -t1453, t1728, -t1286, -t1287, 0, 0;];
m_new_reg  = t1;