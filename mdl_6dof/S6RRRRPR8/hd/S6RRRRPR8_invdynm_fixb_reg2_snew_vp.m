% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6RRRRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 21:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6RRRRPR8_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR8_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR8_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR8_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_invdynm_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 21:56:13
% EndTime: 2019-05-07 21:57:12
% DurationCPUTime: 63.05s
% Computational Cost: add. (304060->959), mult. (614016->1198), div. (0->0), fcn. (448348->10), ass. (0->626)
t1891 = sin(qJ(3));
t1896 = cos(qJ(3));
t1892 = sin(qJ(2));
t2026 = qJD(1) * t1892;
t1848 = qJD(2) * t1891 + t1896 * t2026;
t1890 = sin(qJ(4));
t1895 = cos(qJ(4));
t1967 = qJD(2) * t1896 - t1891 * t2026;
t1805 = t1895 * t1848 + t1890 * t1967;
t1802 = t1805 ^ 2;
t1897 = cos(qJ(2));
t2025 = qJD(1) * t1897;
t1879 = -qJD(3) + t2025;
t1870 = -qJD(4) + t1879;
t2040 = t1870 ^ 2;
t1728 = t2040 + t1802;
t1882 = qJD(2) * t2026;
t1978 = t1897 * qJDD(1);
t1852 = -t1882 + t1978;
t1845 = -qJDD(3) + t1852;
t1842 = -qJDD(4) + t1845;
t1803 = t1848 * t1890 - t1895 * t1967;
t2000 = t1805 * t1803;
t2053 = -t2000 + t1842;
t2009 = t2053 * t1890;
t1665 = -t1728 * t1895 + t2009;
t2008 = t2053 * t1895;
t1667 = t1728 * t1890 + t2008;
t1560 = t1665 * t1891 - t1667 * t1896;
t1972 = qJD(2) * t2025;
t1979 = t1892 * qJDD(1);
t1931 = t1972 + t1979;
t1798 = qJD(3) * t1967 + t1891 * qJDD(2) + t1896 * t1931;
t1906 = t1896 * qJDD(2) - t1891 * t1931;
t1903 = -t1848 * qJD(3) + t1906;
t1694 = -t1803 * qJD(4) + t1895 * t1798 + t1890 * t1903;
t2001 = t1803 * t1870;
t2058 = t1694 + t2001;
t1515 = t1560 * t1897 - t1892 * t2058;
t1588 = t1665 * t1896 + t1667 * t1891;
t1893 = sin(qJ(1));
t1898 = cos(qJ(1));
t2149 = pkin(6) * (t1515 * t1898 - t1588 * t1893);
t2148 = pkin(6) * (t1515 * t1893 + t1588 * t1898);
t1513 = t1560 * t1892 + t1897 * t2058;
t2147 = pkin(1) * t1513;
t2146 = pkin(7) * t1513;
t2145 = pkin(1) * t1588 + pkin(7) * t1515;
t1786 = t1805 * t1870;
t1965 = t1890 * t1798 - t1895 * t1903;
t1932 = qJD(4) * t1805 + t1965;
t2056 = -t1786 + t1932;
t1579 = -t2056 * t1890 + t1895 * t2058;
t2012 = t2058 * t1890;
t1583 = -t2056 * t1895 - t2012;
t1494 = -t1579 * t1891 + t1583 * t1896;
t2041 = t1803 ^ 2;
t2055 = t1802 - t2041;
t1479 = t1494 * t1897 + t1892 * t2055;
t1490 = t1579 * t1896 + t1583 * t1891;
t2142 = t1479 * t1893 - t1490 * t1898;
t2141 = t1479 * t1898 + t1490 * t1893;
t1782 = t2041 - t2040;
t1676 = t1782 * t1890 - t2008;
t1680 = t1782 * t1895 + t2009;
t1601 = t1676 * t1891 - t1680 * t1896;
t1651 = t1786 + t1932;
t1532 = t1601 * t1897 + t1651 * t1892;
t1598 = t1676 * t1896 + t1680 * t1891;
t2140 = t1532 * t1893 + t1598 * t1898;
t2139 = t1532 * t1898 - t1598 * t1893;
t2137 = pkin(2) * t1588;
t2136 = pkin(8) * t1588;
t2133 = pkin(2) * t2058 + pkin(8) * t1560;
t1477 = t1494 * t1892 - t1897 * t2055;
t1528 = t1601 * t1892 - t1651 * t1897;
t2128 = pkin(3) * t1665;
t2127 = pkin(9) * t1665;
t2126 = pkin(9) * t1667;
t1783 = t1802 - t2040;
t2054 = -t2000 - t1842;
t2007 = t2054 * t1890;
t2089 = -t1895 * t1783 + t2007;
t1712 = t1895 * t2054;
t2090 = t1783 * t1890 + t1712;
t2104 = t1891 * t2090 + t1896 * t2089;
t2052 = -t2001 + t1694;
t2103 = -t1891 * t2089 + t1896 * t2090;
t2117 = t1892 * t2052 + t1897 * t2103;
t2125 = t1893 * t2117 - t1898 * t2104;
t2124 = t1893 * t2104 + t1898 * t2117;
t2051 = -t2040 - t2041;
t2062 = t1895 * t2051 - t2007;
t2065 = t1890 * t2051 + t1712;
t2084 = -t1891 * t2065 + t1896 * t2062;
t2106 = t1892 * t2084 - t1897 * t2056;
t2123 = pkin(1) * t2106;
t2122 = pkin(7) * t2106;
t2083 = t1891 * t2062 + t1896 * t2065;
t2105 = t1892 * t2056 + t1897 * t2084;
t2119 = -pkin(1) * t2083 + pkin(7) * t2105;
t2118 = t1892 * t2103 - t1897 * t2052;
t2116 = pkin(6) * (t1893 * t2105 - t1898 * t2083);
t2115 = pkin(6) * (t1893 * t2083 + t1898 * t2105);
t2113 = pkin(2) * t2083;
t2112 = pkin(8) * t2083;
t2107 = -pkin(2) * t2056 + pkin(8) * t2084;
t1699 = -t2041 - t1802;
t2102 = pkin(2) * t1699;
t2101 = pkin(3) * t1699;
t2100 = pkin(3) * t2065;
t2099 = pkin(9) * t2062;
t2098 = pkin(9) * t2065;
t2097 = qJ(5) * t2058;
t2096 = t1699 * t1892;
t2095 = t1699 * t1897;
t1941 = t1967 * t1848;
t2066 = -t1845 + t1941;
t2092 = t1891 * t2066;
t2091 = t1896 * t2066;
t1927 = (t1803 * t1890 + t1805 * t1895) * t1870;
t1989 = t1870 * t1890;
t1780 = t1805 * t1989;
t1988 = t1870 * t1895;
t1973 = t1803 * t1988;
t1945 = -t1780 + t1973;
t2044 = t1891 * t1945 + t1896 * t1927;
t2043 = -t1891 * t1927 + t1896 * t1945;
t2060 = -t1842 * t1892 + t1897 * t2043;
t2088 = t1893 * t2060 - t1898 * t2044;
t1904 = t1890 * t1932 - t1973;
t1946 = -t1803 * t1989 - t1895 * t1932;
t2045 = t1891 * t1904 + t1896 * t1946;
t1975 = t1892 * t2000;
t2046 = -t1891 * t1946 + t1896 * t1904;
t2061 = t1897 * t2046 - t1975;
t2087 = t1893 * t2061 - t1898 * t2045;
t2086 = t1893 * t2044 + t1898 * t2060;
t2085 = t1893 * t2045 + t1898 * t2061;
t1889 = sin(qJ(6));
t1894 = cos(qJ(6));
t1734 = -t1894 * t1803 + t1805 * t1889;
t1736 = t1803 * t1889 + t1805 * t1894;
t1662 = t1736 * t1734;
t1835 = qJDD(6) + t1842;
t2059 = -t1662 + t1835;
t2080 = t1889 * t2059;
t2073 = t1894 * t2059;
t1868 = g(1) * t1898 + g(2) * t1893;
t1899 = qJD(1) ^ 2;
t1841 = -pkin(1) * t1899 + qJDD(1) * pkin(7) - t1868;
t2037 = pkin(2) * t1897;
t1952 = -pkin(8) * t1892 - t2037;
t1850 = t1952 * qJD(1);
t2029 = t1897 * g(3);
t2042 = qJD(2) ^ 2;
t1776 = -qJDD(2) * pkin(2) - pkin(8) * t2042 + (qJD(1) * t1850 + t1841) * t1892 + t2029;
t1823 = -pkin(3) * t1879 - pkin(9) * t1848;
t1963 = t1967 ^ 2;
t1671 = -t1903 * pkin(3) - t1963 * pkin(9) + t1848 * t1823 + t1776;
t1901 = t1932 * pkin(4) + t1671 - t2097;
t2064 = t1897 * t1842 + t1892 * t2043;
t1974 = t1897 * t2000;
t2063 = t1892 * t2046 + t1974;
t1867 = t1893 * g(1) - t1898 * g(2);
t1840 = qJDD(1) * pkin(1) + t1899 * pkin(7) + t1867;
t1851 = 0.2e1 * t1972 + t1979;
t1949 = -t1852 + t1882;
t1757 = pkin(2) * t1949 - pkin(8) * t1851 - t1840;
t1822 = -g(3) * t1892 + t1897 * t1841;
t1777 = -pkin(2) * t2042 + qJDD(2) * pkin(8) + t1850 * t2025 + t1822;
t1696 = -t1896 * t1757 + t1891 * t1777;
t1828 = t1967 * t1879;
t1762 = t1828 + t1798;
t1616 = pkin(3) * t2066 - pkin(9) * t1762 - t1696;
t1697 = t1891 * t1757 + t1896 * t1777;
t1623 = -pkin(3) * t1963 + pkin(9) * t1903 + t1879 * t1823 + t1697;
t1535 = -t1895 * t1616 + t1890 * t1623;
t1737 = pkin(4) * t1803 - qJ(5) * t1805;
t1505 = t1842 * pkin(4) - qJ(5) * t2040 + t1805 * t1737 + qJDD(5) + t1535;
t1465 = -pkin(5) * t2054 - pkin(10) * t2052 + t1505;
t2024 = qJD(5) * t1870;
t1854 = -0.2e1 * t2024;
t1536 = t1890 * t1616 + t1895 * t1623;
t1940 = -pkin(4) * t2040 - t1842 * qJ(5) - t1803 * t1737 + t1536;
t1503 = t1854 + t1940;
t1951 = pkin(5) * t1870 - pkin(10) * t1805;
t1468 = -pkin(5) * t2041 + pkin(10) * t1932 - t1870 * t1951 + t1503;
t1415 = -t1894 * t1465 + t1468 * t1889;
t1416 = t1889 * t1465 + t1894 * t1468;
t1372 = -t1894 * t1415 + t1416 * t1889;
t1586 = -t1734 * qJD(6) + t1894 * t1694 + t1889 * t1932;
t1864 = qJD(6) + t1870;
t1713 = t1864 * t1734;
t2057 = -t1713 + t1586;
t1763 = t1828 - t1798;
t1644 = t1694 * t1890 - t1805 * t1988;
t1645 = t1694 * t1895 + t1780;
t1568 = t1644 * t1896 + t1645 * t1891;
t1571 = -t1644 * t1891 + t1645 * t1896;
t1947 = t1897 * t1571 + t1975;
t2048 = -t1898 * t1568 + t1893 * t1947;
t1966 = t1889 * t1694 - t1894 * t1932;
t1550 = (qJD(6) - t1864) * t1736 + t1966;
t2047 = t1568 * t1893 + t1898 * t1947;
t1731 = t1734 ^ 2;
t1732 = t1736 ^ 2;
t1844 = t1848 ^ 2;
t1859 = t1864 ^ 2;
t1876 = t1879 ^ 2;
t2039 = pkin(4) + pkin(5);
t2038 = pkin(2) * t1892;
t1451 = -t1535 * t1895 + t1536 * t1890;
t2036 = pkin(3) * t1451;
t1646 = t1895 * t2052;
t1653 = (-qJD(4) - t1870) * t1805 - t1965;
t1580 = t1653 * t1890 - t1646;
t2035 = pkin(3) * t1580;
t2034 = pkin(4) * t1895;
t2033 = pkin(5) * t1372;
t1553 = t1713 + t1586;
t1460 = -t1550 * t1889 - t1553 * t1894;
t2032 = pkin(5) * t1460;
t2031 = pkin(10) * t1372;
t1373 = t1889 * t1415 + t1894 * t1416;
t2030 = pkin(10) * t1373;
t2028 = qJ(5) * t1895;
t2027 = qJD(1) * qJD(2);
t2022 = t1451 * t1891;
t2021 = t1451 * t1896;
t1638 = t1662 + t1835;
t2019 = t1638 * t1889;
t2018 = t1638 * t1894;
t2013 = t2052 * t1890;
t2011 = t1671 * t1890;
t2010 = t1671 * t1895;
t2006 = t1776 * t1891;
t2005 = t1776 * t1896;
t1791 = t1845 + t1941;
t2003 = t1791 * t1891;
t2002 = t1791 * t1896;
t1999 = t1840 * t1892;
t1998 = t1840 * t1897;
t1812 = t1851 * t1892;
t1995 = t1864 * t1736;
t1994 = t1864 * t1889;
t1993 = t1864 * t1894;
t1878 = t1897 * t1899 * t1892;
t1865 = -t1878 + qJDD(2);
t1992 = t1865 * t1892;
t1991 = t1865 * t1897;
t1866 = qJDD(2) + t1878;
t1990 = t1866 * t1892;
t1885 = t1892 ^ 2;
t1987 = t1885 * t1899;
t1970 = -pkin(4) * t1870 - 0.2e1 * qJD(5);
t1481 = t1965 * pkin(5) + t2041 * pkin(10) + t1901 + (pkin(5) * qJD(4) - t1951 + t1970) * t1805;
t1986 = t1889 * t1481;
t1985 = t1891 * t1848;
t1480 = t1894 * t1481;
t1984 = t1896 * t1848;
t1982 = -pkin(4) * t1505 + qJ(5) * t1503;
t1981 = -pkin(4) * t2052 - qJ(5) * t1651;
t1886 = t1897 ^ 2;
t1980 = t1885 + t1886;
t1977 = t1892 * t1662;
t1976 = t1897 * t1662;
t1971 = -qJ(5) * t1890 - pkin(3);
t1705 = -t1732 - t1859;
t1605 = t1705 * t1894 - t2019;
t1969 = -pkin(10) * t1605 - t1480;
t1452 = t1535 * t1890 + t1895 * t1536;
t1612 = t1696 * t1891 + t1896 * t1697;
t1821 = t1892 * t1841 + t2029;
t1766 = t1821 * t1892 + t1897 * t1822;
t1964 = -t1867 * t1893 - t1898 * t1868;
t1962 = t1893 * t1878;
t1961 = t1898 * t1878;
t1435 = t1503 * t1890 - t1505 * t1895;
t1960 = pkin(3) * t1435 + t1982;
t1578 = -t1651 * t1890 - t1646;
t1959 = pkin(3) * t1578 + t1981;
t1958 = -t1536 + t2128;
t1957 = -pkin(4) * t1372 + qJ(5) * t1373 - t2033;
t1462 = -t1550 * t1894 + t1553 * t1889;
t1956 = -pkin(4) * t1460 + qJ(5) * t1462 - t2032;
t1955 = -pkin(5) * t1605 + t1416;
t1953 = -pkin(2) * t1776 + pkin(8) * t1612;
t1858 = qJDD(1) * t1898 - t1893 * t1899;
t1950 = -pkin(6) * t1858 - g(3) * t1893;
t1948 = t1892 * t1571 - t1974;
t1647 = -t1859 - t1731;
t1574 = t1647 * t1889 + t2073;
t1944 = -pkin(10) * t1574 - t1986;
t1575 = t1647 * t1894 - t2080;
t1943 = -pkin(10) * t1575 - t1480;
t1606 = -t1705 * t1889 - t2018;
t1942 = -pkin(10) * t1606 + t1986;
t1611 = -t1696 * t1896 + t1697 * t1891;
t1765 = t1821 * t1897 - t1822 * t1892;
t1939 = t1867 * t1898 - t1868 * t1893;
t1938 = -t1535 + t2100;
t1937 = -pkin(5) * t1574 + t1415;
t1936 = -pkin(10) * t1460 - t1372;
t1935 = -pkin(10) * t1462 - t1373;
t1934 = t1892 * t1941;
t1933 = t1897 * t1941;
t1930 = -pkin(4) * t1605 + qJ(5) * t1606 + t1955;
t1336 = -t1372 * t1895 + t1373 * t1890;
t1929 = pkin(3) * t1336 + t1957;
t1410 = -t1460 * t1895 + t1462 * t1890;
t1928 = pkin(3) * t1410 + t1956;
t1799 = -t1876 - t1963;
t1723 = t1799 * t1896 - t2092;
t1829 = t1879 * t1848;
t1759 = t1829 + t1903;
t1926 = pkin(2) * t1759 + pkin(8) * t1723 - t2005;
t1808 = -t1844 - t1876;
t1730 = -t1808 * t1891 + t2002;
t1925 = pkin(2) * t1763 + pkin(8) * t1730 + t2006;
t1924 = pkin(4) * t1728 - qJ(5) * t2053 + t1940;
t1923 = -pkin(4) * t1574 + qJ(5) * t1575 + t1937;
t1760 = (-qJD(3) - t1879) * t1848 + t1906;
t1689 = t1760 * t1896 + t1762 * t1891;
t1790 = t1963 + t1844;
t1922 = pkin(2) * t1790 + pkin(8) * t1689 + t1612;
t1507 = -t1605 * t1895 + t1606 * t1890;
t1921 = pkin(3) * t1507 + t1930;
t1920 = t1924 - t2128;
t1337 = t1372 * t1890 + t1373 * t1895;
t1344 = -t1481 * t2039 - t2030;
t1358 = -qJ(5) * t1481 - t2031;
t1316 = -pkin(3) * t1481 + pkin(9) * t1337 + t1344 * t1895 + t1358 * t1890;
t1317 = -pkin(9) * t1336 - t1344 * t1890 + t1358 * t1895;
t1323 = -t1336 * t1891 + t1337 * t1896;
t1919 = -pkin(2) * t1481 + pkin(8) * t1323 + t1316 * t1896 + t1317 * t1891;
t1613 = -t1731 - t1732;
t1348 = t1613 * t2039 + t1935;
t1351 = qJ(5) * t1613 + t1936;
t1412 = t1460 * t1890 + t1462 * t1895;
t1326 = pkin(3) * t1613 + pkin(9) * t1412 + t1348 * t1895 + t1351 * t1890;
t1327 = -pkin(9) * t1410 - t1348 * t1890 + t1351 * t1895;
t1368 = -t1410 * t1891 + t1412 * t1896;
t1918 = pkin(2) * t1613 + pkin(8) * t1368 + t1326 * t1896 + t1327 * t1891;
t1549 = (qJD(6) + t1864) * t1736 + t1966;
t1389 = t1549 * t2039 + t1943;
t1423 = qJ(5) * t1549 + t1944;
t1485 = t1574 * t1890 + t1575 * t1895;
t1347 = pkin(3) * t1549 + pkin(9) * t1485 + t1389 * t1895 + t1423 * t1890;
t1484 = -t1574 * t1895 + t1575 * t1890;
t1355 = -pkin(9) * t1484 - t1389 * t1890 + t1423 * t1895;
t1426 = -t1484 * t1891 + t1485 * t1896;
t1917 = pkin(2) * t1549 + pkin(8) * t1426 + t1347 * t1896 + t1355 * t1891;
t1396 = t2039 * t2057 + t1942;
t1424 = qJ(5) * t2057 + t1969;
t1508 = t1605 * t1890 + t1606 * t1895;
t1352 = pkin(3) * t2057 + pkin(9) * t1508 + t1396 * t1895 + t1424 * t1890;
t1359 = -pkin(9) * t1507 - t1396 * t1890 + t1424 * t1895;
t1438 = -t1507 * t1891 + t1508 * t1896;
t1916 = pkin(2) * t2057 + pkin(8) * t1438 + t1352 * t1896 + t1359 * t1891;
t1436 = t1503 * t1895 + t1505 * t1890;
t1517 = t1970 * t1805 + t1901;
t1380 = pkin(9) * t1436 + (t1971 - t2034) * t1517;
t1387 = -t1435 * t1891 + t1436 * t1896;
t1388 = -pkin(9) * t1435 + (pkin(4) * t1890 - t2028) * t1517;
t1915 = -pkin(2) * t1517 + pkin(8) * t1387 + t1380 * t1896 + t1388 * t1891;
t1486 = -pkin(4) * t1699 + t1503;
t1487 = -qJ(5) * t1699 + t1505;
t1582 = -t1651 * t1895 + t2013;
t1400 = pkin(9) * t1582 + t1486 * t1895 + t1487 * t1890 - t2101;
t1404 = -pkin(9) * t1578 - t1486 * t1890 + t1487 * t1895;
t1493 = -t1578 * t1891 + t1582 * t1896;
t1914 = pkin(8) * t1493 + t1400 * t1896 + t1404 * t1891 - t2102;
t1584 = t1653 * t1895 + t2013;
t1428 = pkin(9) * t1584 + t1452 - t2101;
t1430 = -pkin(9) * t1580 - t1451;
t1495 = -t1580 * t1891 + t1584 * t1896;
t1913 = pkin(8) * t1495 + t1428 * t1896 + t1430 * t1891 - t2102;
t1900 = 0.2e1 * qJD(5) * t1805 - t1901;
t1496 = pkin(4) * t1786 + t1900 + t2097;
t1431 = -t2126 + t1890 * t1496 + (pkin(3) + t2034) * t2058;
t1444 = -pkin(4) * t2012 + t1496 * t1895 + t2127;
t1912 = t1431 * t1896 + t1444 * t1891 + t2133;
t1497 = (-t2056 + t1786) * pkin(4) + t1900;
t1434 = t1895 * t1497 + t1971 * t2056 + t2099;
t1448 = -t1497 * t1890 - t2028 * t2056 - t2098;
t1911 = t1434 * t1896 + t1448 * t1891 + t2107;
t1519 = -pkin(3) * t2056 - t2010 + t2099;
t1587 = t2011 - t2098;
t1910 = t1519 * t1896 + t1587 * t1891 + t2107;
t1525 = -pkin(3) * t2058 + t2011 + t2126;
t1592 = t2010 - t2127;
t1909 = t1525 * t1896 + t1592 * t1891 - t2133;
t1908 = pkin(3) * t1484 + t1923;
t1402 = t1452 * t1896 - t2022;
t1447 = -pkin(3) * t1671 + pkin(9) * t1452;
t1907 = -pkin(2) * t1671 + pkin(8) * t1402 - pkin(9) * t2022 + t1447 * t1896;
t1905 = pkin(4) * t2054 + qJ(5) * t2051 - t1505;
t1902 = t1905 + t2100;
t1883 = t1886 * t1899;
t1875 = -t1883 - t2042;
t1874 = t1883 - t2042;
t1873 = -t1987 - t2042;
t1872 = -t1987 + t2042;
t1861 = -t1883 + t1987;
t1860 = t1883 + t1987;
t1857 = qJDD(1) * t1893 + t1898 * t1899;
t1856 = t1980 * qJDD(1);
t1853 = -0.2e1 * t1882 + t1978;
t1847 = t1897 * t1866;
t1846 = t1980 * t2027;
t1833 = -pkin(6) * t1857 + g(3) * t1898;
t1827 = -t1844 + t1876;
t1826 = t1963 - t1876;
t1825 = -t1885 * t2027 + t1897 * t1931;
t1824 = -t1852 * t1892 - t1886 * t2027;
t1820 = -t1873 * t1892 - t1991;
t1819 = -t1872 * t1892 + t1847;
t1818 = t1875 * t1897 - t1990;
t1817 = t1874 * t1897 - t1992;
t1816 = t1873 * t1897 - t1992;
t1815 = t1872 * t1897 + t1990;
t1814 = t1875 * t1892 + t1847;
t1813 = t1874 * t1892 + t1991;
t1811 = t1949 * t1897;
t1809 = t1844 - t1963;
t1807 = t1853 * t1897 - t1812;
t1806 = t1851 * t1897 + t1853 * t1892;
t1779 = -pkin(7) * t1816 - t1998;
t1778 = -pkin(7) * t1814 - t1999;
t1774 = (-t1896 * t1967 - t1985) * t1879;
t1773 = (-t1891 * t1967 + t1984) * t1879;
t1771 = -pkin(1) * t1816 + t1822;
t1770 = -pkin(1) * t1814 + t1821;
t1758 = t1829 - t1903;
t1756 = pkin(1) * t1853 + pkin(7) * t1818 + t1998;
t1755 = -pkin(1) * t1851 + pkin(7) * t1820 - t1999;
t1752 = t1798 * t1896 + t1879 * t1985;
t1751 = t1798 * t1891 - t1879 * t1984;
t1750 = t1828 * t1896 - t1891 * t1903;
t1749 = -t1828 * t1891 - t1896 * t1903;
t1745 = t1774 * t1897 - t1845 * t1892;
t1744 = t1774 * t1892 + t1845 * t1897;
t1741 = t1826 * t1896 + t2003;
t1740 = -t1827 * t1891 + t2091;
t1739 = t1826 * t1891 - t2002;
t1738 = t1827 * t1896 + t2092;
t1733 = pkin(1) * t1840 + pkin(7) * t1766;
t1729 = t1808 * t1896 + t2003;
t1725 = pkin(1) * t1860 + pkin(7) * t1856 + t1766;
t1722 = t1799 * t1891 + t2091;
t1711 = -t1732 + t1859;
t1710 = t1731 - t1859;
t1704 = t1897 * t1752 - t1934;
t1703 = t1897 * t1750 + t1934;
t1702 = t1892 * t1752 + t1933;
t1701 = t1892 * t1750 - t1933;
t1688 = t1759 * t1896 + t1763 * t1891;
t1687 = t1760 * t1891 - t1762 * t1896;
t1686 = t1759 * t1891 - t1763 * t1896;
t1685 = -pkin(8) * t1729 + t2005;
t1684 = t1741 * t1897 - t1758 * t1892;
t1683 = t1740 * t1897 + t1762 * t1892;
t1682 = t1741 * t1892 + t1758 * t1897;
t1681 = t1740 * t1892 - t1762 * t1897;
t1672 = -pkin(8) * t1722 + t2006;
t1670 = t1730 * t1897 - t1763 * t1892;
t1669 = t1730 * t1892 + t1763 * t1897;
t1664 = t1723 * t1897 - t1759 * t1892;
t1663 = t1723 * t1892 + t1759 * t1897;
t1661 = t1732 - t1731;
t1660 = t1688 * t1897 + t1809 * t1892;
t1659 = t1688 * t1892 - t1809 * t1897;
t1632 = t1689 * t1897 - t1790 * t1892;
t1631 = t1689 * t1892 + t1790 * t1897;
t1630 = (-t1734 * t1894 + t1736 * t1889) * t1864;
t1629 = (-t1734 * t1889 - t1736 * t1894) * t1864;
t1628 = -pkin(2) * t1729 + t1697;
t1622 = -pkin(2) * t1722 + t1696;
t1610 = t1710 * t1894 - t2019;
t1609 = -t1711 * t1889 + t2073;
t1608 = t1710 * t1889 + t2018;
t1607 = t1711 * t1894 + t2080;
t1594 = t1612 * t1897 + t1776 * t1892;
t1593 = t1612 * t1892 - t1776 * t1897;
t1585 = -qJD(6) * t1736 - t1966;
t1576 = -pkin(1) * t1669 - t1925;
t1562 = -pkin(1) * t1663 - t1926;
t1557 = -pkin(8) * t1687 - t1611;
t1556 = t1629 * t1890 + t1630 * t1895;
t1555 = -t1629 * t1895 + t1630 * t1890;
t1548 = t1586 * t1894 - t1736 * t1994;
t1547 = t1586 * t1889 + t1736 * t1993;
t1546 = -t1585 * t1889 + t1734 * t1993;
t1545 = -t1585 * t1894 - t1734 * t1994;
t1520 = -pkin(7) * t1669 - t1628 * t1892 + t1685 * t1897;
t1518 = -pkin(7) * t1663 - t1622 * t1892 + t1672 * t1897;
t1512 = t1608 * t1890 + t1610 * t1895;
t1511 = t1607 * t1890 + t1609 * t1895;
t1510 = -t1608 * t1895 + t1610 * t1890;
t1509 = -t1607 * t1895 + t1609 * t1890;
t1506 = -pkin(1) * t1729 + pkin(7) * t1670 + t1628 * t1897 + t1685 * t1892;
t1501 = -pkin(1) * t1631 - t1922;
t1500 = -pkin(1) * t1593 - t1953;
t1499 = -pkin(1) * t1722 + pkin(7) * t1664 + t1622 * t1897 + t1672 * t1892;
t1498 = -pkin(7) * t1631 + t1557 * t1897 + t1687 * t2038;
t1491 = t1580 * t1896 + t1584 * t1891;
t1489 = t1578 * t1896 + t1582 * t1891;
t1483 = -t1555 * t1891 + t1556 * t1896;
t1482 = t1555 * t1896 + t1556 * t1891;
t1475 = t1483 * t1897 - t1835 * t1892;
t1474 = t1483 * t1892 + t1835 * t1897;
t1473 = -pkin(7) * t1593 + (-pkin(8) * t1897 + t2038) * t1611;
t1472 = t1495 * t1897 + t2096;
t1471 = t1493 * t1897 + t2096;
t1470 = t1495 * t1892 - t2095;
t1469 = t1493 * t1892 - t2095;
t1466 = pkin(7) * t1632 + t1892 * t1557 + (-pkin(1) - t2037) * t1687;
t1461 = -t1549 * t1894 - t1889 * t2057;
t1459 = -t1549 * t1889 + t1894 * t2057;
t1456 = t1547 * t1890 + t1548 * t1895;
t1455 = -t1545 * t1890 + t1546 * t1895;
t1454 = -t1547 * t1895 + t1548 * t1890;
t1453 = t1545 * t1895 + t1546 * t1890;
t1450 = -t1958 - t2137;
t1449 = -t1938 - t2113;
t1446 = -pkin(2) * t1491 - t2035;
t1445 = pkin(7) * t1594 + (-pkin(1) + t1952) * t1611;
t1443 = -t1525 * t1891 + t1592 * t1896 - t2136;
t1442 = -t1510 * t1891 + t1512 * t1896;
t1441 = -t1509 * t1891 + t1511 * t1896;
t1440 = t1510 * t1896 + t1512 * t1891;
t1439 = t1509 * t1896 + t1511 * t1891;
t1437 = t1507 * t1896 + t1508 * t1891;
t1433 = -t1902 - t2113;
t1432 = -t1519 * t1891 + t1587 * t1896 - t2112;
t1429 = -t1920 + 0.2e1 * t2024 + t2137;
t1427 = -pkin(2) * t1489 - t1959;
t1425 = t1484 * t1896 + t1485 * t1891;
t1422 = t1442 * t1897 + t1550 * t1892;
t1421 = t1441 * t1897 - t1553 * t1892;
t1420 = t1442 * t1892 - t1550 * t1897;
t1419 = t1441 * t1892 + t1553 * t1897;
t1418 = t1438 * t1897 - t1892 * t2057;
t1417 = t1438 * t1892 + t1897 * t2057;
t1411 = t1459 * t1890 + t1461 * t1895;
t1409 = -t1459 * t1895 + t1461 * t1890;
t1408 = -t1454 * t1891 + t1456 * t1896;
t1407 = -t1453 * t1891 + t1455 * t1896;
t1406 = t1454 * t1896 + t1456 * t1891;
t1405 = t1453 * t1896 + t1455 * t1891;
t1403 = -t1909 + t2147;
t1401 = t1452 * t1891 + t2021;
t1399 = t1426 * t1897 - t1549 * t1892;
t1398 = t1426 * t1892 + t1549 * t1897;
t1397 = -t1910 - t2123;
t1395 = t1408 * t1897 - t1977;
t1394 = t1407 * t1897 + t1977;
t1393 = t1408 * t1892 + t1976;
t1392 = t1407 * t1892 - t1976;
t1391 = t1402 * t1897 + t1671 * t1892;
t1390 = t1402 * t1892 - t1671 * t1897;
t1386 = t1435 * t1896 + t1436 * t1891;
t1385 = t1443 * t1897 - t1450 * t1892 + t2146;
t1384 = -t1434 * t1891 + t1448 * t1896 - t2112;
t1383 = -t1431 * t1891 + t1444 * t1896 + t2136;
t1382 = -pkin(2) * t1401 - t2036;
t1381 = t1432 * t1897 - t1449 * t1892 - t2122;
t1379 = t1443 * t1892 + t1450 * t1897 - t2145;
t1378 = t1387 * t1897 + t1517 * t1892;
t1377 = t1387 * t1892 - t1517 * t1897;
t1376 = t1432 * t1892 + t1449 * t1897 + t2119;
t1375 = -t1911 - t2123;
t1374 = -pkin(8) * t1491 - t1428 * t1891 + t1430 * t1896;
t1369 = -t1912 - t2147;
t1367 = -t1409 * t1891 + t1411 * t1896;
t1366 = t1410 * t1896 + t1412 * t1891;
t1365 = t1409 * t1896 + t1411 * t1891;
t1364 = t1367 * t1897 - t1661 * t1892;
t1363 = t1367 * t1892 + t1661 * t1897;
t1362 = t1368 * t1897 - t1613 * t1892;
t1361 = t1368 * t1892 + t1613 * t1897;
t1360 = -pkin(8) * t1401 - pkin(9) * t2021 - t1447 * t1891;
t1357 = -pkin(1) * t1470 - t1913;
t1356 = t1384 * t1897 - t1433 * t1892 - t2122;
t1354 = -pkin(2) * t1437 - t1921;
t1353 = -pkin(8) * t1489 - t1400 * t1891 + t1404 * t1896;
t1350 = t1383 * t1897 - t1429 * t1892 - t2146;
t1349 = t1384 * t1892 + t1433 * t1897 + t2119;
t1346 = -pkin(2) * t1386 - t1960;
t1345 = t1383 * t1892 + t1429 * t1897 + t2145;
t1343 = -pkin(2) * t1425 - t1908;
t1342 = -pkin(7) * t1470 + t1374 * t1897 - t1446 * t1892;
t1341 = -pkin(1) * t1469 - t1914;
t1340 = -pkin(1) * t1390 - t1907;
t1339 = -pkin(1) * t1491 + pkin(7) * t1472 + t1374 * t1892 + t1446 * t1897;
t1338 = -pkin(7) * t1469 + t1353 * t1897 - t1427 * t1892;
t1335 = -pkin(1) * t1489 + pkin(7) * t1471 + t1353 * t1892 + t1427 * t1897;
t1334 = -pkin(8) * t1386 - t1380 * t1891 + t1388 * t1896;
t1333 = -pkin(2) * t1366 - t1928;
t1332 = -pkin(7) * t1390 + t1360 * t1897 - t1382 * t1892;
t1331 = -pkin(8) * t1437 - t1352 * t1891 + t1359 * t1896;
t1330 = -pkin(1) * t1401 + pkin(7) * t1391 + t1360 * t1892 + t1382 * t1897;
t1329 = -pkin(8) * t1425 - t1347 * t1891 + t1355 * t1896;
t1328 = -pkin(1) * t1377 - t1915;
t1325 = -pkin(1) * t1417 - t1916;
t1324 = -pkin(1) * t1398 - t1917;
t1322 = t1336 * t1896 + t1337 * t1891;
t1321 = t1323 * t1897 + t1481 * t1892;
t1320 = t1323 * t1892 - t1481 * t1897;
t1319 = -pkin(7) * t1377 + t1334 * t1897 - t1346 * t1892;
t1318 = -pkin(7) * t1417 + t1331 * t1897 - t1354 * t1892;
t1315 = -pkin(1) * t1437 + pkin(7) * t1418 + t1331 * t1892 + t1354 * t1897;
t1314 = -pkin(7) * t1398 + t1329 * t1897 - t1343 * t1892;
t1313 = -pkin(1) * t1386 + pkin(7) * t1378 + t1334 * t1892 + t1346 * t1897;
t1312 = -pkin(1) * t1425 + pkin(7) * t1399 + t1329 * t1892 + t1343 * t1897;
t1311 = -pkin(8) * t1366 - t1326 * t1891 + t1327 * t1896;
t1310 = -pkin(1) * t1361 - t1918;
t1309 = -pkin(2) * t1322 - t1929;
t1308 = -pkin(7) * t1361 + t1311 * t1897 - t1333 * t1892;
t1307 = -pkin(1) * t1366 + pkin(7) * t1362 + t1311 * t1892 + t1333 * t1897;
t1306 = -pkin(8) * t1322 - t1316 * t1891 + t1317 * t1896;
t1305 = -pkin(1) * t1320 - t1919;
t1304 = -pkin(7) * t1320 + t1306 * t1897 - t1309 * t1892;
t1303 = -pkin(1) * t1322 + pkin(7) * t1321 + t1306 * t1892 + t1309 * t1897;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1858, 0, -t1857, 0, t1950, -t1833, -t1939, -pkin(6) * t1939, t1825 * t1898 - t1962, t1807 * t1898 + t1861 * t1893, t1819 * t1898 + t1893 * t1979, t1824 * t1898 + t1962, t1817 * t1898 + t1893 * t1978, qJDD(2) * t1893 + t1846 * t1898, t1898 * t1778 - t1893 * t1770 - pkin(6) * (t1818 * t1893 + t1853 * t1898), t1898 * t1779 - t1893 * t1771 - pkin(6) * (t1820 * t1893 - t1851 * t1898), t1898 * t1765 - pkin(6) * (t1856 * t1893 + t1860 * t1898), -pkin(6) * (t1766 * t1893 + t1840 * t1898) - (pkin(1) * t1893 - pkin(7) * t1898) * t1765, t1704 * t1898 + t1751 * t1893, t1660 * t1898 + t1686 * t1893, t1683 * t1898 + t1738 * t1893, t1703 * t1898 - t1749 * t1893, t1684 * t1898 + t1739 * t1893, t1745 * t1898 + t1773 * t1893, t1898 * t1518 - t1893 * t1562 - pkin(6) * (t1664 * t1893 - t1722 * t1898), t1898 * t1520 - t1893 * t1576 - pkin(6) * (t1670 * t1893 - t1729 * t1898), t1898 * t1498 - t1893 * t1501 - pkin(6) * (t1632 * t1893 - t1687 * t1898), t1898 * t1473 - t1893 * t1500 - pkin(6) * (t1594 * t1893 - t1611 * t1898), t2047, t2141, t2124, t2085, -t2139, t2086, t1898 * t1381 - t1893 * t1397 - t2116, t1898 * t1385 - t1893 * t1403 + t2148, t1898 * t1342 - t1893 * t1357 - pkin(6) * (t1472 * t1893 - t1491 * t1898), t1898 * t1332 - t1893 * t1340 - pkin(6) * (t1391 * t1893 - t1401 * t1898), t2047, t2124, -t2141, t2086, t2139, t2085, t1898 * t1356 - t1893 * t1375 - t2116, t1898 * t1338 - t1893 * t1341 - pkin(6) * (t1471 * t1893 - t1489 * t1898), t1898 * t1350 - t1893 * t1369 - t2148, t1898 * t1319 - t1893 * t1328 - pkin(6) * (t1378 * t1893 - t1386 * t1898), t1395 * t1898 + t1406 * t1893, t1364 * t1898 + t1365 * t1893, t1421 * t1898 + t1439 * t1893, t1394 * t1898 + t1405 * t1893, t1422 * t1898 + t1440 * t1893, t1475 * t1898 + t1482 * t1893, t1898 * t1314 - t1893 * t1324 - pkin(6) * (t1399 * t1893 - t1425 * t1898), t1898 * t1318 - t1893 * t1325 - pkin(6) * (t1418 * t1893 - t1437 * t1898), t1898 * t1308 - t1893 * t1310 - pkin(6) * (t1362 * t1893 - t1366 * t1898), t1898 * t1304 - t1893 * t1305 - pkin(6) * (t1321 * t1893 - t1322 * t1898); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1857, 0, t1858, 0, t1833, t1950, t1964, pkin(6) * t1964, t1825 * t1893 + t1961, t1807 * t1893 - t1861 * t1898, t1819 * t1893 - t1898 * t1979, t1824 * t1893 - t1961, t1817 * t1893 - t1898 * t1978, -qJDD(2) * t1898 + t1846 * t1893, t1893 * t1778 + t1898 * t1770 + pkin(6) * (t1818 * t1898 - t1853 * t1893), t1893 * t1779 + t1898 * t1771 + pkin(6) * (t1820 * t1898 + t1851 * t1893), t1893 * t1765 + pkin(6) * (t1856 * t1898 - t1860 * t1893), pkin(6) * (t1766 * t1898 - t1840 * t1893) - (-pkin(1) * t1898 - pkin(7) * t1893) * t1765, t1704 * t1893 - t1751 * t1898, t1660 * t1893 - t1686 * t1898, t1683 * t1893 - t1738 * t1898, t1703 * t1893 + t1749 * t1898, t1684 * t1893 - t1739 * t1898, t1745 * t1893 - t1773 * t1898, t1893 * t1518 + t1898 * t1562 + pkin(6) * (t1664 * t1898 + t1722 * t1893), t1893 * t1520 + t1898 * t1576 + pkin(6) * (t1670 * t1898 + t1729 * t1893), t1893 * t1498 + t1898 * t1501 + pkin(6) * (t1632 * t1898 + t1687 * t1893), t1893 * t1473 + t1898 * t1500 + pkin(6) * (t1594 * t1898 + t1611 * t1893), t2048, t2142, t2125, t2087, -t2140, t2088, t1893 * t1381 + t1898 * t1397 + t2115, t1893 * t1385 + t1898 * t1403 - t2149, t1893 * t1342 + t1898 * t1357 + pkin(6) * (t1472 * t1898 + t1491 * t1893), t1893 * t1332 + t1898 * t1340 + pkin(6) * (t1391 * t1898 + t1401 * t1893), t2048, t2125, -t2142, t2088, t2140, t2087, t1893 * t1356 + t1898 * t1375 + t2115, t1893 * t1338 + t1898 * t1341 + pkin(6) * (t1471 * t1898 + t1489 * t1893), t1893 * t1350 + t1898 * t1369 + t2149, t1893 * t1319 + t1898 * t1328 + pkin(6) * (t1378 * t1898 + t1386 * t1893), t1395 * t1893 - t1406 * t1898, t1364 * t1893 - t1365 * t1898, t1421 * t1893 - t1439 * t1898, t1394 * t1893 - t1405 * t1898, t1422 * t1893 - t1440 * t1898, t1475 * t1893 - t1482 * t1898, t1893 * t1314 + t1898 * t1324 + pkin(6) * (t1399 * t1898 + t1425 * t1893), t1893 * t1318 + t1898 * t1325 + pkin(6) * (t1418 * t1898 + t1437 * t1893), t1893 * t1308 + t1898 * t1310 + pkin(6) * (t1362 * t1898 + t1366 * t1893), t1893 * t1304 + t1898 * t1305 + pkin(6) * (t1321 * t1898 + t1322 * t1893); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1867, t1868, 0, 0, t1812, t1806, t1815, -t1811, t1813, 0, t1756, t1755, t1725, t1733, t1702, t1659, t1681, t1701, t1682, t1744, t1499, t1506, t1466, t1445, t1948, t1477, t2118, t2063, -t1528, t2064, t1376, t1379, t1339, t1330, t1948, t2118, -t1477, t2064, t1528, t2063, t1349, t1335, t1345, t1313, t1393, t1363, t1419, t1392, t1420, t1474, t1312, t1315, t1307, t1303; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1899, 0, 0, -g(3), -t1867, 0, t1825, t1807, t1819, t1824, t1817, t1846, t1778, t1779, t1765, pkin(7) * t1765, t1704, t1660, t1683, t1703, t1684, t1745, t1518, t1520, t1498, t1473, t1947, t1479, t2117, t2061, -t1532, t2060, t1381, t1385, t1342, t1332, t1947, t2117, -t1479, t2060, t1532, t2061, t1356, t1338, t1350, t1319, t1395, t1364, t1421, t1394, t1422, t1475, t1314, t1318, t1308, t1304; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1899, 0, qJDD(1), 0, g(3), 0, -t1868, 0, t1878, -t1861, -t1979, -t1878, -t1978, -qJDD(2), t1770, t1771, 0, pkin(1) * t1765, -t1751, -t1686, -t1738, t1749, -t1739, -t1773, t1562, t1576, t1501, t1500, -t1568, -t1490, -t2104, -t2045, -t1598, -t2044, t1397, t1403, t1357, t1340, -t1568, -t2104, t1490, -t2044, t1598, -t2045, t1375, t1341, t1369, t1328, -t1406, -t1365, -t1439, -t1405, -t1440, -t1482, t1324, t1325, t1310, t1305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1867, t1868, 0, 0, t1812, t1806, t1815, -t1811, t1813, 0, t1756, t1755, t1725, t1733, t1702, t1659, t1681, t1701, t1682, t1744, t1499, t1506, t1466, t1445, t1948, t1477, t2118, t2063, -t1528, t2064, t1376, t1379, t1339, t1330, t1948, t2118, -t1477, t2064, t1528, t2063, t1349, t1335, t1345, t1313, t1393, t1363, t1419, t1392, t1420, t1474, t1312, t1315, t1307, t1303; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1931, t1853, t1866, -t1972, t1874, t1972, 0, -t1840, t1821, 0, t1752, t1688, t1740, t1750, t1741, t1774, t1672, t1685, t1557, -pkin(8) * t1611, t1571, t1494, t2103, t2046, -t1601, t2043, t1432, t1443, t1374, t1360, t1571, t2103, -t1494, t2043, t1601, t2046, t1384, t1353, t1383, t1334, t1408, t1367, t1441, t1407, t1442, t1483, t1329, t1331, t1311, t1306; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1882, t1851, t1872, t1852, t1865, -t1882, t1840, 0, t1822, 0, t1941, -t1809, -t1762, -t1941, t1758, t1845, t1622, t1628, -pkin(2) * t1687, -pkin(2) * t1611, -t2000, -t2055, -t2052, t2000, t1651, t1842, t1449, t1450, t1446, t1382, -t2000, -t2052, t2055, t1842, -t1651, t2000, t1433, t1427, t1429, t1346, t1662, t1661, t1553, -t1662, -t1550, t1835, t1343, t1354, t1333, t1309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1878, t1861, t1979, t1878, t1978, qJDD(2), -t1821, -t1822, 0, 0, t1751, t1686, t1738, -t1749, t1739, t1773, t1926, t1925, t1922, t1953, t1568, t1490, t2104, t2045, t1598, t2044, t1910, t1909, t1913, t1907, t1568, t2104, -t1490, t2044, -t1598, t2045, t1911, t1914, t1912, t1915, t1406, t1365, t1439, t1405, t1440, t1482, t1917, t1916, t1918, t1919; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1798, t1759, t2066, t1828, t1826, -t1828, 0, t1776, t1696, 0, t1645, t1583, t2090, t1904, t1680, t1945, t1587, t1592, t1430, -pkin(9) * t1451, t1645, t2090, -t1583, t1945, -t1680, t1904, t1448, t1404, t1444, t1388, t1456, t1411, t1511, t1455, t1512, t1556, t1355, t1359, t1327, t1317; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1829, -t1763, t1827, t1903, -t1791, t1829, -t1776, 0, t1697, 0, t1644, t1579, t2089, t1946, t1676, t1927, t1519, t1525, t1428, t1447, t1644, t2089, -t1579, t1927, -t1676, t1946, t1434, t1400, t1431, t1380, t1454, t1409, t1509, t1453, t1510, t1555, t1347, t1352, t1326, t1316; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1941, t1809, t1762, t1941, -t1758, -t1845, -t1696, -t1697, 0, 0, t2000, t2055, t2052, -t2000, -t1651, -t1842, t1938, t1958, t2035, t2036, t2000, t2052, -t2055, -t1842, t1651, -t2000, t1902, t1959, t1854 + t1920, t1960, -t1662, -t1661, -t1553, t1662, t1550, -t1835, t1908, t1921, t1928, t1929; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1694, -t2056, t2054, -t2001, t1782, t2001, 0, t1671, t1535, 0, t1694, t2054, t2056, t2001, -t1782, -t2001, -qJ(5) * t2056, t1487, t1496, -qJ(5) * t1517, t1548, t1461, t1609, t1546, t1610, t1630, t1423, t1424, t1351, t1358; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1786, t2058, -t1783, -t1932, -t2053, t1786, -t1671, 0, t1536, 0, -t1786, -t1783, -t2058, t1786, t2053, -t1932, t1497, t1486, pkin(4) * t2058, -pkin(4) * t1517, -t1547, -t1459, -t1607, t1545, -t1608, -t1629, t1389, t1396, t1348, t1344; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2000, t2055, t2052, -t2000, -t1651, -t1842, -t1535, -t1536, 0, 0, t2000, t2052, -t2055, -t1842, t1651, -t2000, t1905, t1981, t1854 + t1924, t1982, -t1662, -t1661, -t1553, t1662, t1550, -t1835, t1923, t1930, t1956, t1957; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1694, t2054, t2056, t2001, -t1782, -t2001, 0, t1505, -t1517, 0, t1548, t1461, t1609, t1546, t1610, t1630, t1944, t1969, t1936, -t2031; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2000, t2052, -t2055, -t1842, t1651, -t2000, -t1505, 0, t1503, 0, -t1662, -t1661, -t1553, t1662, t1550, -t1835, t1937, t1955, -t2032, -t2033; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1786, t1783, t2058, -t1786, -t2053, t1932, t1517, -t1503, 0, 0, t1547, t1459, t1607, -t1545, t1608, t1629, -pkin(5) * t1549 - t1943, -pkin(5) * t2057 - t1942, -pkin(5) * t1613 - t1935, pkin(5) * t1481 + t2030; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1586, -t1549, t2059, t1713, t1710, -t1713, 0, -t1481, t1415, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1995, t2057, t1711, t1585, t1638, -t1995, t1481, 0, t1416, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1662, t1661, t1553, -t1662, -t1550, t1835, -t1415, -t1416, 0, 0;];
m_new_reg  = t1;