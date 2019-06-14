% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6RRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 04:34
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6RRRPPR3_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR3_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_invdynm_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:33:34
% EndTime: 2019-05-07 04:34:14
% DurationCPUTime: 42.67s
% Computational Cost: add. (97609->876), mult. (206756->978), div. (0->0), fcn. (140344->8), ass. (0->564)
t1691 = sin(qJ(1));
t1689 = sin(qJ(3));
t1690 = sin(qJ(2));
t1693 = cos(qJ(3));
t1694 = cos(qJ(2));
t1641 = (t1689 * t1694 + t1690 * t1693) * qJD(1);
t1852 = qJD(1) * t1694;
t1674 = qJD(2) * t1852;
t1794 = t1690 * qJDD(1);
t1651 = t1674 + t1794;
t1681 = t1694 * qJDD(1);
t1853 = qJD(1) * t1690;
t1784 = qJD(2) * t1853;
t1758 = -t1681 + t1784;
t1769 = t1689 * t1651 + t1693 * t1758;
t1740 = -qJD(3) * t1641 - t1769;
t1685 = qJD(2) + qJD(3);
t1814 = t1685 * t1641;
t1922 = -t1740 + t1814;
t1639 = t1689 * t1853 - t1693 * t1852;
t1636 = t1639 ^ 2;
t1880 = t1685 ^ 2;
t1778 = -t1880 - t1636;
t1684 = qJDD(2) + qJDD(3);
t1822 = t1641 * t1639;
t1921 = -t1822 + t1684;
t1945 = t1921 * t1689;
t1906 = t1778 * t1693 - t1945;
t1944 = t1921 * t1693;
t1907 = t1778 * t1689 + t1944;
t1416 = t1690 * t1907 - t1694 * t1906;
t1695 = cos(qJ(1));
t1949 = t1416 * t1695;
t2003 = pkin(6) * (-t1691 * t1922 + t1949);
t1951 = t1416 * t1691;
t2002 = pkin(6) * (t1695 * t1922 + t1951);
t1804 = qJD(3) - t1685;
t1496 = -t1641 * t1804 - t1769;
t1891 = t1636 - t1880;
t1920 = t1822 + t1684;
t1962 = t1920 * t1689;
t1524 = -t1693 * t1891 + t1962;
t1961 = t1920 * t1693;
t1958 = t1689 * t1891 + t1961;
t1434 = t1524 * t1694 + t1690 * t1958;
t1976 = t1434 * t1691;
t2001 = t1496 * t1695 + t1976;
t1974 = t1434 * t1695;
t2000 = t1496 * t1691 - t1974;
t1725 = t1693 * t1651 - t1689 * t1758;
t1501 = t1639 * t1804 - t1725;
t1486 = t1693 * t1501;
t1497 = t1814 + t1740;
t1402 = t1497 * t1689 + t1486;
t1841 = t1501 * t1689;
t1408 = t1497 * t1693 - t1841;
t1334 = t1402 * t1690 - t1408 * t1694;
t1881 = t1641 ^ 2;
t1549 = t1881 + t1636;
t1927 = t1549 * t1691;
t1999 = pkin(6) * (t1334 * t1695 + t1927);
t1926 = t1549 * t1695;
t1998 = pkin(6) * (t1334 * t1691 - t1926);
t1953 = pkin(7) * t1416;
t1995 = pkin(1) * t1922 + t1953;
t1330 = t1402 * t1694 + t1408 * t1690;
t1994 = pkin(1) * t1330;
t1917 = t1690 * t1906 + t1694 * t1907;
t1956 = pkin(1) * t1917;
t1607 = -t1881 - t1880;
t1968 = -t1607 * t1693 + t1962;
t1969 = t1607 * t1689 + t1961;
t1970 = t1690 * t1969 + t1694 * t1968;
t1983 = pkin(1) * t1970;
t1993 = pkin(7) * t1330;
t1414 = t1690 * t1968 - t1694 * t1969;
t1992 = pkin(7) * t1414;
t1954 = pkin(7) * t1917;
t1982 = pkin(7) * t1970;
t1991 = t1414 * t1691;
t1990 = t1414 * t1695;
t1939 = pkin(1) * t1549;
t1989 = pkin(7) * t1334 - t1939;
t1546 = -t1639 * qJD(3) + t1725;
t1815 = t1685 * t1639;
t1500 = -t1815 + t1546;
t1988 = pkin(1) * t1500 - t1992;
t1985 = pkin(6) * (t1500 * t1695 - t1991);
t1984 = pkin(6) * (-t1500 * t1691 - t1990);
t1428 = t1524 * t1690 - t1694 * t1958;
t1875 = pkin(2) * t1402;
t1981 = pkin(2) * t1907;
t1873 = pkin(2) * t1968;
t1980 = pkin(8) * t1402;
t1935 = pkin(8) * t1906;
t1934 = pkin(8) * t1907;
t1979 = pkin(8) * t1968;
t1978 = pkin(8) * t1969;
t1888 = t1881 - t1880;
t1959 = t1888 * t1689 + t1944;
t1960 = -t1693 * t1888 + t1945;
t1965 = -t1690 * t1960 + t1694 * t1959;
t1977 = t1691 * t1965;
t1975 = t1695 * t1965;
t1938 = pkin(2) * t1549;
t1973 = pkin(8) * t1408 + t1938;
t1503 = t1815 + t1546;
t1972 = -t1503 * t1695 + t1977;
t1971 = t1503 * t1691 + t1975;
t1964 = t1690 * t1959 + t1694 * t1960;
t1963 = pkin(3) * (t1740 - t1922);
t1893 = -pkin(3) * t1607 + qJ(4) * t1920;
t1787 = t1695 * t1822;
t1812 = t1685 * t1693;
t1755 = t1689 * t1546 + t1641 * t1812;
t1813 = t1685 * t1689;
t1785 = t1641 * t1813;
t1925 = t1693 * t1546 - t1785;
t1941 = -t1690 * t1755 + t1694 * t1925;
t1757 = t1691 * t1941 - t1787;
t1788 = t1691 * t1822;
t1756 = t1695 * t1941 + t1788;
t1952 = qJ(4) * t1922;
t1892 = -pkin(3) * t1921 - qJ(4) * t1778;
t1942 = t1690 * t1925 + t1694 * t1755;
t1936 = pkin(3) * t1501;
t1752 = qJ(4) * t1500;
t1933 = qJ(4) * t1549;
t1842 = t1500 * t1689;
t1931 = t1500 * t1693;
t1889 = t1881 - t1636;
t1911 = t1691 * t1889;
t1908 = t1695 * t1889;
t1687 = t1694 ^ 2;
t1697 = qJD(1) ^ 2;
t1660 = t1691 * g(1) - t1695 * g(2);
t1746 = qJDD(1) * pkin(1) + t1660;
t1550 = (pkin(8) * t1687 + pkin(7)) * t1697 - pkin(2) * t1758 - (qJD(2) * pkin(2) - pkin(8) * t1853) * t1853 + t1746;
t1705 = -pkin(3) * t1814 + t1550;
t1701 = t1752 + t1705;
t1887 = -pkin(4) * t1740 + t1636 * qJ(5) - qJDD(5);
t1924 = t1701 - t1887;
t1578 = pkin(3) * t1639 - qJ(4) * t1641;
t1807 = t1690 * t1697;
t1661 = g(1) * t1695 + g(2) * t1691;
t1644 = -pkin(1) * t1697 + qJDD(1) * pkin(7) - t1661;
t1808 = t1690 * t1644;
t1854 = qJD(1) * qJD(2);
t1528 = qJDD(2) * pkin(2) - t1651 * pkin(8) - t1808 + (pkin(2) * t1807 + pkin(8) * t1854 - g(3)) * t1694;
t1610 = -t1690 * g(3) + t1694 * t1644;
t1682 = t1687 * t1697;
t1696 = qJD(2) ^ 2;
t1665 = -t1682 - t1696;
t1533 = pkin(2) * t1665 + t1681 * pkin(8) + t1610;
t1448 = t1689 * t1528 + t1693 * t1533;
t1727 = -t1880 * pkin(3) + t1684 * qJ(4) + t1448;
t1716 = -0.2e1 * qJD(4) * t1685 + t1639 * t1578 - t1727;
t1710 = -t1716 + t1893;
t1923 = t1710 + t1873;
t1919 = 0.2e1 * t1752;
t1688 = sin(qJ(6));
t1692 = cos(qJ(6));
t1593 = t1639 * t1688 + t1692 * t1685;
t1595 = t1639 * t1692 - t1685 * t1688;
t1534 = t1595 * t1593;
t1540 = qJDD(6) + t1546;
t1895 = -t1534 + t1540;
t1913 = t1688 * t1895;
t1910 = t1692 * t1895;
t1670 = t1695 * t1684;
t1786 = t1639 * t1812;
t1753 = t1785 - t1786;
t1734 = (-t1639 * t1689 - t1641 * t1693) * t1685;
t1900 = t1690 * t1734;
t1885 = t1694 * t1753 - t1900;
t1905 = t1691 * t1885 - t1670;
t1741 = -t1689 * t1740 + t1786;
t1754 = t1639 * t1813 + t1693 * t1740;
t1884 = -t1690 * t1754 + t1694 * t1741;
t1904 = t1691 * t1884 + t1787;
t1903 = t1695 * t1884 - t1788;
t1816 = t1684 * t1691;
t1902 = t1695 * t1885 + t1816;
t1901 = qJ(5) * (-t1497 - t1740);
t1898 = t1694 * t1734;
t1611 = -pkin(4) * t1685 - qJ(5) * t1641;
t1803 = 0.2e1 * qJD(4) + t1611;
t1283 = -(-pkin(3) - pkin(9)) * t1740 + t1500 * pkin(5) + (-pkin(9) * t1685 + t1803) * t1641 + t1924;
t1579 = pkin(5) * t1641 - pkin(9) * t1639;
t1447 = -t1693 * t1528 + t1689 * t1533;
t1733 = -t1684 * pkin(3) - t1880 * qJ(4) + qJDD(4) + t1447;
t1721 = -t1684 * pkin(4) - qJ(5) * t1503 + t1733;
t1776 = -pkin(4) * t1639 - t1578;
t1879 = -2 * qJD(5);
t1763 = t1879 - t1776;
t1299 = -t1880 * pkin(5) - t1684 * pkin(9) + (-t1579 + t1763) * t1641 + t1721;
t1239 = -t1692 * t1283 + t1688 * t1299;
t1240 = t1283 * t1688 + t1299 * t1692;
t1211 = t1688 * t1239 + t1240 * t1692;
t1896 = -t1692 * t1239 + t1688 * t1240;
t1890 = -pkin(4) * t1922 - qJ(5) * t1778;
t1634 = qJD(6) + t1641;
t1771 = t1692 * t1684 - t1688 * t1740;
t1441 = (qJD(6) - t1634) * t1595 + t1771;
t1886 = t1690 * t1753 + t1898;
t1883 = t1690 * t1741 + t1694 * t1754;
t1589 = t1593 ^ 2;
t1882 = t1595 ^ 2;
t1631 = t1634 ^ 2;
t1878 = pkin(3) + pkin(4);
t1877 = -pkin(4) - pkin(9);
t1349 = -t1447 * t1693 + t1448 * t1689;
t1876 = pkin(2) * t1349;
t1871 = pkin(3) * t1693;
t1320 = t1641 * t1763 + t1721;
t1870 = pkin(4) * t1320;
t1868 = pkin(4) * t1501;
t1867 = pkin(4) * t1607;
t1866 = t1740 * pkin(3);
t1864 = qJ(4) * t1693;
t1633 = t1636 * pkin(4);
t1712 = t1685 * t1803 - t1633 + t1727;
t1802 = (2 * qJD(5)) - t1578;
t1856 = t1740 * qJ(5);
t1308 = -t1856 + t1684 * pkin(5) - t1880 * pkin(9) + (t1579 + t1802) * t1639 + t1712;
t1863 = qJ(5) * t1308;
t1862 = qJ(5) * t1320;
t1708 = t1639 * t1802 + t1712;
t1322 = t1708 - t1856;
t1861 = qJ(5) * t1322;
t1860 = qJ(5) * t1607;
t1859 = qJ(5) * t1921;
t1858 = qJ(5) * t1920;
t1855 = pkin(3) - t1877;
t1851 = qJD(4) * t1641;
t1848 = t1349 * t1690;
t1847 = t1349 * t1694;
t1796 = -t1685 - qJD(3);
t1492 = -t1641 * t1796 + t1769;
t1846 = t1492 * t1689;
t1845 = t1492 * t1693;
t1840 = t1550 * t1689;
t1839 = t1550 * t1693;
t1827 = t1634 * t1593;
t1826 = t1634 * t1595;
t1825 = t1634 * t1688;
t1824 = t1634 * t1692;
t1823 = t1641 * t1578;
t1643 = t1697 * pkin(7) + t1746;
t1821 = t1643 * t1690;
t1820 = t1643 * t1694;
t1652 = t1681 - 0.2e1 * t1784;
t1596 = t1652 * t1694;
t1667 = t1694 * t1807;
t1658 = qJDD(2) + t1667;
t1819 = t1658 * t1690;
t1659 = qJDD(2) - t1667;
t1818 = t1659 * t1690;
t1817 = t1659 * t1694;
t1686 = t1690 ^ 2;
t1811 = t1686 * t1697;
t1301 = t1688 * t1308;
t1453 = t1534 + t1540;
t1809 = t1688 * t1453;
t1806 = t1692 * t1308;
t1805 = t1692 * t1453;
t1800 = -pkin(5) * t1308 + pkin(9) * t1211;
t1368 = t1733 + t1823;
t1799 = -pkin(3) * t1368 - qJ(4) * t1716;
t1798 = qJ(4) * t1496 + t1936;
t1795 = t1686 + t1687;
t1793 = 0.2e1 * t1851;
t1792 = t1631 - t1882;
t1790 = t1689 * t1534;
t1789 = t1693 * t1534;
t1485 = -t1631 - t1589;
t1365 = t1688 * t1485 + t1910;
t1783 = -pkin(5) * t1365 + t1239;
t1779 = -t1631 - t1882;
t1374 = -t1688 * t1779 - t1805;
t1780 = t1593 * qJD(6) + t1688 * t1684 + t1692 * t1740;
t1443 = -t1780 - t1827;
t1782 = -pkin(5) * t1443 + pkin(9) * t1374 + t1301;
t1366 = t1692 * t1485 - t1913;
t1440 = (qJD(6) + t1634) * t1595 + t1771;
t1781 = -pkin(5) * t1440 + pkin(9) * t1366 - t1806;
t1777 = -qJ(4) * t1689 - pkin(2);
t1208 = pkin(5) * t1896;
t1775 = -qJ(5) * t1211 + t1208;
t1774 = -qJ(5) * t1443 - t1806;
t1444 = t1780 - t1827;
t1346 = -t1441 * t1688 + t1444 * t1692;
t1343 = pkin(5) * t1346;
t1348 = -t1441 * t1692 - t1688 * t1444;
t1773 = -qJ(5) * t1348 + t1343;
t1350 = t1447 * t1689 + t1693 * t1448;
t1609 = t1694 * g(3) + t1808;
t1532 = t1609 * t1690 + t1694 * t1610;
t1768 = -t1660 * t1691 - t1695 * t1661;
t1767 = t1691 * t1667;
t1766 = t1695 * t1667;
t1296 = -t1368 * t1693 - t1689 * t1716;
t1765 = pkin(2) * t1296 + t1799;
t1401 = t1496 * t1689 + t1486;
t1764 = pkin(2) * t1401 + t1798;
t1762 = -pkin(4) * t1211 - t1800;
t1761 = -pkin(3) * t1320 + qJ(4) * t1322 - t1870;
t1760 = -qJ(4) * t1497 - t1868 - t1936;
t1655 = qJDD(1) * t1695 - t1691 * t1697;
t1759 = -pkin(6) * t1655 - g(3) * t1691;
t1751 = -qJ(5) * t1440 - t1301;
t1531 = t1609 * t1694 - t1610 * t1690;
t1750 = t1660 * t1695 - t1661 * t1691;
t1748 = -pkin(4) * t1374 - t1782;
t1747 = -pkin(4) * t1366 - t1781;
t1745 = -qJ(5) * t1366 - t1783;
t1466 = -t1589 - t1882;
t1744 = -pkin(5) * t1466 + pkin(9) * t1348 + t1211;
t1743 = -t1447 + t1981;
t1742 = -qJ(5) * t1466 + t1896;
t1739 = -qJD(6) * t1595 - t1771;
t1373 = t1692 * t1779 - t1809;
t1738 = -pkin(5) * t1373 + t1240;
t1737 = -pkin(3) * t1211 + qJ(4) * t1308 + t1762;
t1264 = -t1320 * t1693 + t1322 * t1689;
t1736 = pkin(2) * t1264 + t1761;
t1735 = t1760 - t1875;
t1732 = -t1448 - t1873;
t1731 = -pkin(3) * t1374 + qJ(4) * t1443 + t1748;
t1730 = -pkin(3) * t1366 + qJ(4) * t1440 + t1747;
t1729 = -pkin(4) * t1348 - t1744;
t1204 = -t1211 * t1693 + t1308 * t1689;
t1728 = pkin(2) * t1204 + t1737;
t1726 = -qJ(5) * t1374 - t1738;
t1315 = -t1374 * t1693 + t1443 * t1689;
t1724 = pkin(2) * t1315 + t1731;
t1312 = -t1366 * t1693 + t1440 * t1689;
t1723 = pkin(2) * t1312 + t1730;
t1722 = -pkin(3) * t1348 + qJ(4) * t1466 + t1729;
t1720 = -t1641 * t1879 - t1721;
t1718 = -t1368 - t1892;
t1305 = -t1348 * t1693 + t1466 * t1689;
t1717 = pkin(2) * t1305 + t1722;
t1714 = t1639 * t1796 + t1725;
t1713 = t1718 + t1981;
t1709 = -t1685 * t1611 + t1639 * t1879 + t1633 + t1716;
t1707 = -pkin(4) * t1921 + t1320;
t1706 = -qJ(5) * t1501 + t1641 * t1776 + t1720;
t1704 = t1707 + t1892;
t1703 = t1705 + t1866;
t1702 = t1322 - t1867;
t1700 = t1702 + t1893;
t1699 = t1701 + t1793;
t1698 = t1803 * t1641 + t1703 - t1887;
t1314 = t1752 + t1698;
t1664 = t1682 - t1696;
t1663 = -t1696 - t1811;
t1662 = t1696 - t1811;
t1657 = -t1682 + t1811;
t1656 = t1682 + t1811;
t1654 = qJDD(1) * t1691 + t1695 * t1697;
t1653 = t1795 * qJDD(1);
t1650 = 0.2e1 * t1674 + t1794;
t1648 = t1694 * t1658;
t1647 = t1795 * t1854;
t1635 = -pkin(6) * t1654 + g(3) * t1695;
t1614 = t1651 * t1694 - t1686 * t1854;
t1613 = -t1687 * t1854 + t1690 * t1758;
t1605 = -t1663 * t1690 - t1817;
t1604 = -t1662 * t1690 + t1648;
t1603 = t1665 * t1694 - t1819;
t1602 = t1664 * t1694 - t1818;
t1601 = t1663 * t1694 - t1818;
t1600 = t1662 * t1694 + t1819;
t1599 = t1665 * t1690 + t1648;
t1598 = t1664 * t1690 + t1817;
t1597 = (t1651 + t1674) * t1690;
t1587 = -t1650 * t1690 + t1596;
t1586 = t1650 * t1694 + t1652 * t1690;
t1559 = t1589 - t1631;
t1558 = (-t1639 * t1693 + t1641 * t1689) * t1685;
t1552 = -pkin(7) * t1601 - t1820;
t1551 = -pkin(7) * t1599 - t1821;
t1543 = -pkin(1) * t1601 + t1610;
t1542 = -pkin(1) * t1599 + t1609;
t1529 = -t1589 + t1882;
t1526 = pkin(1) * t1652 + pkin(7) * t1603 + t1820;
t1525 = -pkin(1) * t1650 + pkin(7) * t1605 - t1821;
t1484 = pkin(1) * t1643 + pkin(7) * t1532;
t1473 = pkin(1) * t1656 + pkin(7) * t1653 + t1532;
t1461 = (-t1593 * t1692 + t1595 * t1688) * t1634;
t1460 = (-t1593 * t1688 - t1595 * t1692) * t1634;
t1459 = t1558 * t1694 - t1900;
t1456 = t1558 * t1690 + t1898;
t1449 = -t1859 + t1952;
t1445 = -t1839 + t1979;
t1439 = -t1840 - t1934;
t1420 = -t1595 * t1825 - t1692 * t1780;
t1419 = t1595 * t1824 - t1688 * t1780;
t1418 = -t1593 * t1824 + t1688 * t1739;
t1417 = -t1593 * t1825 - t1692 * t1739;
t1410 = -t1693 * t1922 - t1842;
t1409 = -t1689 * t1714 - t1845;
t1407 = t1496 * t1693 - t1841;
t1405 = t1842 + t1845;
t1404 = -t1689 * t1922 + t1931;
t1403 = t1693 * t1714 - t1846;
t1399 = t1846 - t1931;
t1398 = t1461 * t1689 + t1540 * t1693;
t1397 = -t1461 * t1693 + t1540 * t1689;
t1396 = t1559 * t1692 - t1809;
t1395 = -t1688 * t1792 + t1910;
t1394 = t1559 * t1688 + t1805;
t1393 = t1692 * t1792 + t1913;
t1369 = t1500 * t1878 - t1858;
t1359 = -pkin(2) * t1714 - t1840 - t1978;
t1358 = t1420 * t1689 + t1789;
t1357 = -t1418 * t1689 - t1789;
t1356 = -t1420 * t1693 + t1790;
t1355 = t1418 * t1693 - t1790;
t1354 = -pkin(2) * t1492 + t1839 + t1935;
t1353 = t1699 + t1866;
t1352 = t1368 + t1933;
t1351 = pkin(3) * t1549 - t1716;
t1347 = -t1692 * t1440 - t1443 * t1688;
t1345 = -t1688 * t1440 + t1443 * t1692;
t1341 = t1699 + t1963;
t1340 = t1703 + t1793 + t1919;
t1339 = pkin(2) * t1550 + pkin(8) * t1350;
t1338 = -t1404 * t1690 + t1410 * t1694;
t1337 = -t1403 * t1690 + t1409 * t1694;
t1335 = -t1401 * t1690 + t1407 * t1694;
t1333 = -t1399 * t1690 + t1405 * t1694;
t1332 = t1404 * t1694 + t1410 * t1690;
t1331 = t1403 * t1694 + t1409 * t1690;
t1329 = t1401 * t1694 + t1407 * t1690;
t1327 = t1399 * t1694 + t1405 * t1690;
t1326 = t1396 * t1689 - t1441 * t1693;
t1325 = t1395 * t1689 - t1444 * t1693;
t1324 = -t1396 * t1693 - t1441 * t1689;
t1323 = -t1395 * t1693 - t1444 * t1689;
t1318 = -t1397 * t1690 + t1398 * t1694;
t1317 = t1397 * t1694 + t1398 * t1690;
t1316 = t1374 * t1689 + t1443 * t1693;
t1313 = t1366 * t1689 + t1440 * t1693;
t1311 = -t1732 + t1983;
t1310 = t1347 * t1689 + t1529 * t1693;
t1309 = -t1347 * t1693 + t1529 * t1689;
t1306 = t1348 * t1689 + t1466 * t1693;
t1300 = -t1743 - t1956;
t1297 = t1368 * t1689 - t1693 * t1716;
t1295 = -t1349 - t1980;
t1294 = -t1356 * t1690 + t1358 * t1694;
t1293 = -t1355 * t1690 + t1357 * t1694;
t1292 = t1356 * t1694 + t1358 * t1690;
t1291 = t1355 * t1694 + t1357 * t1690;
t1290 = -t1341 * t1689 - t1864 * t1922 - t1934;
t1289 = t1698 - t1860 + t1919;
t1288 = t1706 - t1933;
t1287 = -pkin(3) * t1842 + t1340 * t1693 - t1979;
t1286 = t1350 + t1973;
t1285 = -t1875 - t1994;
t1284 = -t1713 - t1956;
t1281 = -t1359 * t1690 + t1445 * t1694 + t1982;
t1280 = -t1549 * t1878 + t1709 - t1901;
t1279 = -t1641 * t1611 - 0.2e1 * t1851 - t1890 - t1924 - t1963;
t1278 = t1350 * t1694 - t1848;
t1277 = t1350 * t1690 + t1847;
t1276 = t1693 * t1341 + t1777 * t1922 + t1935;
t1275 = t1978 + t1689 * t1340 + (pkin(2) + t1871) * t1500;
t1274 = -t1923 - t1983;
t1273 = -t1354 * t1690 + t1439 * t1694 - t1954;
t1272 = -pkin(1) * t1714 + t1359 * t1694 + t1445 * t1690 + t1992;
t1271 = -pkin(1) * t1492 + t1354 * t1694 + t1439 * t1690 - t1953;
t1270 = -pkin(1) * t1329 - t1764;
t1269 = -t1324 * t1690 + t1326 * t1694;
t1268 = -t1323 * t1690 + t1325 * t1694;
t1267 = t1324 * t1694 + t1326 * t1690;
t1266 = t1323 * t1694 + t1325 * t1690;
t1265 = t1320 * t1689 + t1322 * t1693;
t1263 = -pkin(8) * t1401 - t1351 * t1689 + t1352 * t1693;
t1262 = (t1921 - t1822) * pkin(4) + t1981 + t1720 - t1823 + t1956 - t1892;
t1261 = qJ(4) * t1314 - t1862;
t1260 = -t1315 * t1690 + t1316 * t1694;
t1259 = t1315 * t1694 + t1316 * t1690;
t1258 = -t1312 * t1690 + t1313 * t1694;
t1257 = t1312 * t1694 + t1313 * t1690;
t1256 = -t1735 + t1994;
t1255 = t1709 + t1856 + t1867 - t1873 - t1893 - t1983;
t1254 = -t1279 * t1689 + t1449 * t1693 + t1934;
t1253 = pkin(8) * t1407 + t1351 * t1693 + t1352 * t1689 + t1938;
t1252 = -t1309 * t1690 + t1310 * t1694;
t1251 = t1309 * t1694 + t1310 * t1690;
t1250 = t1289 * t1693 - t1369 * t1689 - t1979;
t1249 = -t1305 * t1690 + t1306 * t1694;
t1248 = t1305 * t1694 + t1306 * t1690;
t1247 = -pkin(1) * t1277 - t1876;
t1246 = pkin(2) * t1922 + t1279 * t1693 + t1449 * t1689 - t1935;
t1245 = -t1296 * t1690 + t1297 * t1694;
t1244 = t1296 * t1694 + t1297 * t1690;
t1243 = pkin(2) * t1500 + t1289 * t1689 + t1369 * t1693 + t1978;
t1242 = -pkin(8) * t1296 + (-pkin(3) * t1689 + t1864) * t1353;
t1241 = qJ(4) * t1346 + t1773;
t1235 = t1314 * t1878 - t1861;
t1234 = pkin(8) * t1297 + (-t1777 + t1871) * t1353;
t1233 = t1373 * t1855 + t1774;
t1232 = -pkin(7) * t1277 - pkin(8) * t1847 - t1339 * t1690;
t1231 = -t1280 * t1689 + t1288 * t1693 + t1980;
t1230 = t1365 * t1855 + t1751;
t1229 = -t1276 * t1690 + t1290 * t1694 - t1954;
t1228 = pkin(1) * t1550 + pkin(7) * t1278 - pkin(8) * t1848 + t1339 * t1694;
t1227 = t1280 * t1693 + t1288 * t1689 - t1973;
t1226 = -t1275 * t1690 + t1287 * t1694 - t1982;
t1225 = -t1286 * t1690 + t1295 * t1694 - t1993;
t1224 = t1276 * t1694 + t1290 * t1690 - t1995;
t1223 = t1286 * t1694 + t1295 * t1690 - t1989;
t1222 = t1275 * t1694 + t1287 * t1690 + t1988;
t1221 = -t1264 * t1690 + t1265 * t1694;
t1220 = t1264 * t1694 + t1265 * t1690;
t1219 = qJ(4) * t1373 + t1726;
t1218 = qJ(4) * t1365 + t1745;
t1217 = -pkin(1) * t1244 - t1765;
t1216 = -t1246 * t1690 + t1254 * t1694 + t1954;
t1215 = -pkin(7) * t1329 - t1253 * t1690 + t1263 * t1694;
t1214 = t1246 * t1694 + t1254 * t1690 + t1995;
t1213 = -t1243 * t1690 + t1250 * t1694 - t1982;
t1212 = pkin(7) * t1335 + t1253 * t1694 + t1263 * t1690 + t1939;
t1206 = t1243 * t1694 + t1250 * t1690 + t1988;
t1205 = t1211 * t1689 + t1308 * t1693;
t1203 = -pkin(1) * t1259 - t1724;
t1202 = -pkin(1) * t1257 - t1723;
t1201 = -pkin(8) * t1264 - t1235 * t1689 + t1261 * t1693;
t1200 = -t1227 * t1690 + t1231 * t1694 + t1993;
t1199 = t1227 * t1694 + t1231 * t1690 + t1989;
t1198 = pkin(2) * t1314 + pkin(8) * t1265 + t1235 * t1693 + t1261 * t1689;
t1197 = -pkin(7) * t1244 - t1234 * t1690 + t1242 * t1694;
t1196 = t1346 * t1855 + t1742;
t1195 = -pkin(8) * t1315 + t1219 * t1693 - t1233 * t1689;
t1194 = -pkin(1) * t1220 - t1736;
t1193 = -pkin(8) * t1312 + t1218 * t1693 - t1230 * t1689;
t1192 = pkin(1) * t1353 + pkin(7) * t1245 + t1234 * t1694 + t1242 * t1690;
t1191 = pkin(2) * t1373 + pkin(8) * t1316 + t1219 * t1689 + t1233 * t1693;
t1190 = pkin(2) * t1365 + pkin(8) * t1313 + t1218 * t1689 + t1230 * t1693;
t1189 = -pkin(1) * t1248 - t1717;
t1188 = -pkin(8) * t1305 - t1196 * t1689 + t1241 * t1693;
t1187 = pkin(2) * t1346 + pkin(8) * t1306 + t1196 * t1693 + t1241 * t1689;
t1186 = -t1204 * t1690 + t1205 * t1694;
t1185 = t1204 * t1694 + t1205 * t1690;
t1184 = qJ(4) * t1896 + t1775;
t1183 = t1855 * t1896 - t1863;
t1182 = -pkin(7) * t1220 - t1198 * t1690 + t1201 * t1694;
t1181 = -pkin(7) * t1259 - t1191 * t1690 + t1195 * t1694;
t1180 = pkin(1) * t1314 + pkin(7) * t1221 + t1198 * t1694 + t1201 * t1690;
t1179 = -pkin(7) * t1257 - t1190 * t1690 + t1193 * t1694;
t1178 = pkin(1) * t1373 + pkin(7) * t1260 + t1191 * t1694 + t1195 * t1690;
t1177 = pkin(1) * t1365 + pkin(7) * t1258 + t1190 * t1694 + t1193 * t1690;
t1176 = -pkin(7) * t1248 - t1187 * t1690 + t1188 * t1694;
t1175 = pkin(1) * t1346 + pkin(7) * t1249 + t1187 * t1694 + t1188 * t1690;
t1174 = -pkin(1) * t1185 - t1728;
t1173 = -pkin(8) * t1204 - t1183 * t1689 + t1184 * t1693;
t1172 = pkin(2) * t1896 + pkin(8) * t1205 + t1183 * t1693 + t1184 * t1689;
t1171 = -pkin(7) * t1185 - t1172 * t1690 + t1173 * t1694;
t1170 = pkin(1) * t1896 + pkin(7) * t1186 + t1172 * t1694 + t1173 * t1690;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1655, 0, -t1654, 0, t1759, -t1635, -t1750, -pkin(6) * t1750, t1614 * t1695 - t1767, t1587 * t1695 + t1657 * t1691, t1604 * t1695 + t1691 * t1794, t1613 * t1695 + t1767, t1602 * t1695 + t1681 * t1691, qJDD(2) * t1691 + t1647 * t1695, t1695 * t1551 - t1691 * t1542 - pkin(6) * (t1603 * t1691 + t1652 * t1695), t1695 * t1552 - t1691 * t1543 - pkin(6) * (t1605 * t1691 - t1650 * t1695), t1695 * t1531 - pkin(6) * (t1653 * t1691 + t1656 * t1695), -pkin(6) * (t1532 * t1691 + t1643 * t1695) - (pkin(1) * t1691 - pkin(7) * t1695) * t1531, t1756, t1337 * t1695 + t1911, -t1501 * t1691 + t1975, t1903, t2000, t1902, t1695 * t1273 - t1691 * t1300 - pkin(6) * (-t1492 * t1695 - t1951), t1695 * t1281 - t1691 * t1311 - pkin(6) * (-t1695 * t1714 + t1991), t1695 * t1225 - t1691 * t1285 + t1998, t1695 * t1232 - t1691 * t1247 - pkin(6) * (t1278 * t1691 + t1550 * t1695), t1756, t1971, t1333 * t1695 - t1911, t1902, -t2000, t1903, t1695 * t1229 - t1691 * t1284 + t2002, t1695 * t1215 - t1691 * t1270 - pkin(6) * (t1335 * t1691 + t1926), t1695 * t1226 - t1691 * t1274 - t1985, t1695 * t1197 - t1691 * t1217 - pkin(6) * (t1245 * t1691 + t1353 * t1695), t1903, t1338 * t1695 + t1911, t1497 * t1691 - t1974, t1756, t1971, t1459 * t1695 + t1816, t1695 * t1213 - t1691 * t1255 - t1985, t1695 * t1216 - t1691 * t1262 - t2002, t1695 * t1200 - t1691 * t1256 - t1998, t1695 * t1182 - t1691 * t1194 - pkin(6) * (t1221 * t1691 + t1314 * t1695), t1294 * t1695 - t1419 * t1691, t1252 * t1695 - t1345 * t1691, t1268 * t1695 - t1393 * t1691, t1293 * t1695 + t1417 * t1691, t1269 * t1695 - t1394 * t1691, t1318 * t1695 - t1460 * t1691, t1695 * t1179 - t1691 * t1202 - pkin(6) * (t1258 * t1691 + t1365 * t1695), t1695 * t1181 - t1691 * t1203 - pkin(6) * (t1260 * t1691 + t1373 * t1695), t1695 * t1176 - t1691 * t1189 - pkin(6) * (t1249 * t1691 + t1346 * t1695), t1695 * t1171 - t1691 * t1174 - pkin(6) * (t1186 * t1691 + t1695 * t1896); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1654, 0, t1655, 0, t1635, t1759, t1768, pkin(6) * t1768, t1614 * t1691 + t1766, t1587 * t1691 - t1657 * t1695, t1604 * t1691 - t1695 * t1794, t1613 * t1691 - t1766, t1602 * t1691 - t1681 * t1695, -qJDD(2) * t1695 + t1647 * t1691, t1691 * t1551 + t1695 * t1542 + pkin(6) * (t1603 * t1695 - t1652 * t1691), t1691 * t1552 + t1695 * t1543 + pkin(6) * (t1605 * t1695 + t1650 * t1691), t1691 * t1531 + pkin(6) * (t1653 * t1695 - t1656 * t1691), pkin(6) * (t1532 * t1695 - t1643 * t1691) - (-pkin(1) * t1695 - pkin(7) * t1691) * t1531, t1757, t1337 * t1691 - t1908, t1501 * t1695 + t1977, t1904, -t2001, t1905, t1691 * t1273 + t1695 * t1300 + pkin(6) * (t1492 * t1691 - t1949), t1691 * t1281 + t1695 * t1311 + pkin(6) * (t1691 * t1714 + t1990), t1691 * t1225 + t1695 * t1285 - t1999, t1691 * t1232 + t1695 * t1247 + pkin(6) * (t1278 * t1695 - t1550 * t1691), t1757, t1972, t1333 * t1691 + t1908, t1905, t2001, t1904, t1691 * t1229 + t1695 * t1284 - t2003, t1691 * t1215 + t1695 * t1270 + pkin(6) * (t1335 * t1695 - t1927), t1691 * t1226 + t1695 * t1274 + t1984, t1691 * t1197 + t1695 * t1217 + pkin(6) * (t1245 * t1695 - t1353 * t1691), t1904, t1338 * t1691 - t1908, -t1497 * t1695 - t1976, t1757, t1972, t1459 * t1691 - t1670, t1691 * t1213 + t1695 * t1255 + t1984, t1691 * t1216 + t1695 * t1262 + t2003, t1691 * t1200 + t1695 * t1256 + t1999, t1691 * t1182 + t1695 * t1194 + pkin(6) * (t1221 * t1695 - t1314 * t1691), t1294 * t1691 + t1419 * t1695, t1252 * t1691 + t1345 * t1695, t1268 * t1691 + t1393 * t1695, t1293 * t1691 - t1417 * t1695, t1269 * t1691 + t1394 * t1695, t1318 * t1691 + t1460 * t1695, t1691 * t1179 + t1695 * t1202 + pkin(6) * (t1258 * t1695 - t1365 * t1691), t1691 * t1181 + t1695 * t1203 + pkin(6) * (t1260 * t1695 - t1373 * t1691), t1691 * t1176 + t1695 * t1189 + pkin(6) * (t1249 * t1695 - t1346 * t1691), t1691 * t1171 + t1695 * t1174 + pkin(6) * (t1186 * t1695 - t1691 * t1896); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1660, t1661, 0, 0, t1597, t1586, t1600, t1596, t1598, 0, t1526, t1525, t1473, t1484, t1942, t1331, t1964, t1883, -t1428, t1886, t1271, t1272, t1223, t1228, t1942, t1964, t1327, t1886, t1428, t1883, t1224, t1212, t1222, t1192, t1883, t1332, -t1428, t1942, t1964, t1456, t1206, t1214, t1199, t1180, t1292, t1251, t1266, t1291, t1267, t1317, t1177, t1178, t1175, t1170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1697, 0, 0, -g(3), -t1660, 0, t1614, t1587, t1604, t1613, t1602, t1647, t1551, t1552, t1531, pkin(7) * t1531, t1941, t1337, t1965, t1884, -t1434, t1885, t1273, t1281, t1225, t1232, t1941, t1965, t1333, t1885, t1434, t1884, t1229, t1215, t1226, t1197, t1884, t1338, -t1434, t1941, t1965, t1459, t1213, t1216, t1200, t1182, t1294, t1252, t1268, t1293, t1269, t1318, t1179, t1181, t1176, t1171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1697, 0, qJDD(1), 0, g(3), 0, -t1661, 0, t1667, -t1657, -t1794, -t1667, -t1681, -qJDD(2), t1542, t1543, 0, pkin(1) * t1531, -t1822, -t1889, t1501, t1822, -t1496, -t1684, t1300, t1311, t1285, t1247, -t1822, -t1503, t1889, -t1684, t1496, t1822, t1284, t1270, t1274, t1217, t1822, -t1889, -t1497, -t1822, -t1503, -t1684, t1255, t1262, t1256, t1194, t1419, t1345, t1393, -t1417, t1394, t1460, t1202, t1203, t1189, t1174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1660, t1661, 0, 0, t1597, t1586, t1600, t1596, t1598, 0, t1526, t1525, t1473, t1484, t1942, t1331, t1964, t1883, -t1428, t1886, t1271, t1272, t1223, t1228, t1942, t1964, t1327, t1886, t1428, t1883, t1224, t1212, t1222, t1192, t1883, t1332, -t1428, t1942, t1964, t1456, t1206, t1214, t1199, t1180, t1292, t1251, t1266, t1291, t1267, t1317, t1177, t1178, t1175, t1170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1651, t1652, t1658, -t1674, t1664, t1674, 0, -t1643, t1609, 0, t1925, t1409, t1959, t1741, -t1524, t1753, t1439, t1445, t1295, -pkin(8) * t1349, t1925, t1959, t1405, t1753, t1524, t1741, t1290, t1263, t1287, t1242, t1741, t1410, -t1524, t1925, t1959, t1558, t1250, t1254, t1231, t1201, t1358, t1310, t1325, t1357, t1326, t1398, t1193, t1195, t1188, t1173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1784, t1650, t1662, -t1758, t1659, -t1784, t1643, 0, t1610, 0, t1755, t1403, t1960, t1754, t1958, t1734, t1354, t1359, t1286, t1339, t1755, t1960, t1399, t1734, -t1958, t1754, t1276, t1253, t1275, t1234, t1754, t1404, t1958, t1755, t1960, t1734, t1243, t1246, t1227, t1198, t1356, t1309, t1323, t1355, t1324, t1397, t1190, t1191, t1187, t1172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1667, t1657, t1794, t1667, t1681, qJDD(2), -t1609, -t1610, 0, 0, t1822, t1889, -t1501, -t1822, t1496, t1684, t1743, t1732, t1875, t1876, t1822, t1503, -t1889, t1684, -t1496, -t1822, t1713, t1764, t1923, t1765, -t1822, t1889, t1497, t1822, t1503, t1684, t1700 + t1873, t1704 - t1981, t1735, t1736, -t1419, -t1345, -t1393, t1417, -t1394, -t1460, t1723, t1724, t1717, t1728; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1546, -t1492, t1921, t1815, t1891, -t1815, 0, -t1550, t1447, 0, t1546, t1921, t1492, -t1815, -t1891, t1815, -t1952, t1352, t1340, qJ(4) * t1353, t1815, -t1922, t1891, t1546, t1921, -t1815, t1289, t1449, t1288, t1261, t1534, t1529, -t1444, -t1534, -t1441, t1540, t1218, t1219, t1241, t1184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1814, t1714, -t1888, t1740, t1920, -t1814, t1550, 0, t1448, 0, t1814, -t1888, -t1500, -t1814, -t1920, t1740, t1341, t1351, pkin(3) * t1500, pkin(3) * t1353, t1740, t1500, t1920, t1814, -t1888, -t1814, t1369, t1279, t1280, t1235, -t1420, -t1347, -t1395, t1418, -t1396, -t1461, t1230, t1233, t1196, t1183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1822, t1889, -t1501, -t1822, t1496, t1684, -t1447, -t1448, 0, 0, t1822, t1503, -t1889, t1684, -t1496, -t1822, t1718, t1798, t1710, t1799, -t1822, t1889, t1497, t1822, t1503, t1684, t1700, t1704, t1760, t1761, -t1419, -t1345, -t1393, t1417, -t1394, -t1460, t1730, t1731, t1722, t1737; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1546, t1921, t1492, -t1815, -t1891, t1815, 0, t1368, t1353, 0, t1815, -t1922, t1891, t1546, t1921, -t1815, t1314 - t1860, -t1859, t1706, -t1862, t1534, t1529, -t1444, -t1534, -t1441, t1540, t1745, t1726, t1773, t1775; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1822, t1503, -t1889, t1684, -t1496, -t1822, -t1368, 0, -t1716, 0, -t1822, t1889, t1497, t1822, t1503, t1684, t1702, t1707, -t1868, -t1870, -t1419, -t1345, -t1393, t1417, -t1394, -t1460, t1747, t1748, t1729, t1762; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1814, t1888, t1500, t1814, t1920, -t1740, -t1353, t1716, 0, 0, -t1740, -t1500, -t1920, -t1814, t1888, t1814, -pkin(4) * t1500 + t1858, t1314 + t1890, pkin(4) * t1549 + t1708 + t1901, -pkin(4) * t1314 + t1861, t1420, t1347, t1395, -t1418, t1396, t1461, t1365 * t1877 - t1751, t1373 * t1877 - t1774, t1346 * t1877 - t1742, t1877 * t1896 + t1863; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1740, -t1500, -t1920, -t1814, t1888, t1814, 0, t1314, t1322, 0, t1420, t1347, t1395, -t1418, t1396, t1461, -pkin(9) * t1365 + t1301, -pkin(9) * t1373 + t1806, -pkin(9) * t1346 - t1896, -pkin(9) * t1896; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1815, t1922, -t1891, -t1546, -t1921, t1815, -t1314, 0, t1320, 0, -t1534, -t1529, t1444, t1534, t1441, -t1540, t1783, t1738, -t1343, -t1208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1822, -t1889, -t1497, -t1822, -t1503, -t1684, -t1322, -t1320, 0, 0, t1419, t1345, t1393, -t1417, t1394, t1460, t1781, t1782, t1744, t1800; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1780, -t1440, t1895, t1827, t1559, -t1827, 0, t1308, t1239, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1826, t1443, t1792, t1739, t1453, -t1826, -t1308, 0, t1240, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1534, t1529, -t1444, -t1534, -t1441, t1540, -t1239, -t1240, 0, 0;];
m_new_reg  = t1;