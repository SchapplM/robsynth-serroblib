% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6PRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 03:15
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6PRRPPR4_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR4_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR4_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_invdynm_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:13:51
% EndTime: 2019-05-05 03:14:55
% DurationCPUTime: 69.51s
% Computational Cost: add. (183357->983), mult. (390195->1330), div. (0->0), fcn. (273389->12), ass. (0->671)
t1843 = sin(pkin(11));
t1846 = cos(pkin(11));
t1851 = sin(qJ(3));
t2022 = qJD(2) * t1851;
t1794 = t1843 * qJD(3) + t1846 * t2022;
t1790 = t1794 ^ 2;
t1854 = cos(qJ(3));
t1840 = t1854 ^ 2;
t1856 = qJD(2) ^ 2;
t1835 = t1840 * t1856;
t1739 = t1835 + t1790;
t1792 = -t1846 * qJD(3) + t1843 * t2022;
t1749 = t1794 * t1792;
t1826 = qJD(3) * t2022;
t1975 = t1854 * qJDD(2);
t1804 = -t1826 + t1975;
t2056 = -t1749 + t1804;
t2000 = t2056 * t1846;
t1668 = t1739 * t1843 + t2000;
t1982 = t1854 * qJD(2);
t1827 = qJD(3) * t1982;
t1976 = t1851 * qJDD(2);
t1803 = t1827 + t1976;
t1768 = qJDD(3) * t1843 + t1803 * t1846;
t1968 = t1792 * t1982;
t1867 = t1768 + t1968;
t1580 = t1668 * t1851 - t1854 * t1867;
t1845 = sin(pkin(6));
t1848 = cos(pkin(6));
t1586 = t1668 * t1854 + t1851 * t1867;
t2001 = t2056 * t1843;
t1638 = t1739 * t1846 - t2001;
t1852 = sin(qJ(2));
t1855 = cos(qJ(2));
t1919 = t1586 * t1852 + t1638 * t1855;
t1465 = t1848 * t1580 + t1845 * t1919;
t1471 = -t1845 * t1580 + t1848 * t1919;
t2207 = pkin(7) * (t1465 * t1845 + t1471 * t1848);
t1525 = t1586 * t1855 - t1638 * t1852;
t1844 = sin(pkin(10));
t1847 = cos(pkin(10));
t2206 = qJ(1) * (t1471 * t1847 + t1525 * t1844);
t2205 = qJ(1) * (t1471 * t1844 - t1525 * t1847);
t2204 = pkin(1) * t1465;
t2203 = pkin(1) * t1471;
t2034 = t1792 ^ 2;
t1772 = t2034 - t1835;
t1673 = t1772 * t1846 + t2001;
t1776 = t1794 * t1982;
t1960 = -t1846 * qJDD(3) + t1803 * t1843;
t1712 = -t1960 - t1776;
t1585 = t1673 * t1851 - t1712 * t1854;
t1591 = t1673 * t1854 + t1712 * t1851;
t1667 = t1772 * t1843 - t2000;
t1914 = t1591 * t1852 - t1667 * t1855;
t1476 = -t1845 * t1585 + t1848 * t1914;
t1530 = t1591 * t1855 + t1667 * t1852;
t2200 = t1476 * t1844 - t1530 * t1847;
t2199 = t1476 * t1847 + t1530 * t1844;
t2004 = t1867 * t1843;
t2076 = t1776 - t1960;
t2006 = t2076 * t1846;
t1620 = t2006 - t2004;
t2054 = t1790 - t2034;
t1575 = t1620 * t1851 - t1854 * t2054;
t1577 = t1620 * t1854 + t1851 * t2054;
t2007 = t2076 * t1843;
t1616 = t1846 * t1867 + t2007;
t1920 = t1577 * t1852 - t1616 * t1855;
t1455 = -t1845 * t1575 + t1848 * t1920;
t1512 = t1577 * t1855 + t1616 * t1852;
t2197 = t1455 * t1844 - t1512 * t1847;
t2196 = t1455 * t1847 + t1512 * t1844;
t2195 = pkin(7) * t1525;
t1470 = t1848 * t1585 + t1845 * t1914;
t1453 = t1848 * t1575 + t1845 * t1920;
t2053 = t1790 + t2034;
t1949 = -t1768 + t1968;
t2094 = t1712 * t1846 - t1949 * t1843;
t2116 = t1851 * t2094 + t1854 * t2053;
t2095 = t1712 * t1843 + t1949 * t1846;
t2113 = -t1851 * t2053 + t1854 * t2094;
t2138 = t1852 * t2113 - t1855 * t2095;
t2162 = -t1845 * t2116 + t1848 * t2138;
t2188 = pkin(1) * t2162;
t2164 = t1845 * t2138 + t1848 * t2116;
t2187 = pkin(1) * t2164;
t2186 = pkin(2) * t1580;
t2185 = pkin(8) * t1580;
t2180 = -pkin(2) * t1638 - pkin(8) * t1586;
t1773 = t1790 - t1835;
t2057 = -t1749 - t1804;
t1999 = t2057 * t1843;
t2100 = -t1773 * t1846 + t1999;
t1998 = t2057 * t1846;
t2099 = t1773 * t1843 + t1998;
t2114 = -t1851 * t1949 + t1854 * t2099;
t2135 = t1852 * t2100 + t1855 * t2114;
t2117 = t1851 * t2099 + t1854 * t1949;
t2137 = t1852 * t2114 - t1855 * t2100;
t2163 = -t1845 * t2117 + t1848 * t2137;
t2179 = -t1844 * t2163 + t1847 * t2135;
t2178 = t1844 * t2135 + t1847 * t2163;
t2136 = t1852 * t2095 + t1855 * t2113;
t2177 = qJ(1) * (-t1844 * t2162 + t1847 * t2136);
t2176 = qJ(1) * (t1844 * t2136 + t1847 * t2162);
t2175 = (-t1845 * t2164 - t1848 * t2162) * pkin(7);
t2052 = -t2034 - t1835;
t2073 = t1846 * t2052 - t1999;
t2096 = t1851 * t2073 + t1854 * t2076;
t2074 = t1843 * t2052 + t1998;
t2089 = -t1851 * t2076 + t1854 * t2073;
t2115 = t1852 * t2089 - t1855 * t2074;
t2139 = -t1845 * t2096 + t1848 * t2115;
t2174 = pkin(1) * t2139;
t2140 = t1845 * t2115 + t1848 * t2096;
t2173 = pkin(1) * t2140;
t2172 = pkin(7) * t2136;
t2165 = t1845 * t2137 + t1848 * t2117;
t2112 = t1852 * t2074 + t1855 * t2089;
t2161 = qJ(1) * (-t1844 * t2139 + t1847 * t2112);
t2160 = qJ(1) * (t1844 * t2112 + t1847 * t2139);
t2159 = (-t1845 * t2140 - t1848 * t2139) * pkin(7);
t2157 = pkin(2) * t2116;
t2156 = pkin(3) * t1638;
t2155 = pkin(7) * t2112;
t2154 = pkin(8) * t2113;
t2153 = pkin(8) * t2116;
t2152 = qJ(4) * t1638;
t2151 = qJ(4) * t1668;
t2134 = pkin(2) * t2096;
t2133 = pkin(3) * t2095;
t2132 = pkin(4) * t2076;
t2131 = pkin(8) * t2096;
t2130 = qJ(4) * t2095;
t2123 = -pkin(2) * t2074 + pkin(8) * t2089;
t2122 = pkin(3) * t2053 + qJ(4) * t2094;
t1959 = t1843 * t1968;
t1687 = t1846 * t1960 + t1959;
t1858 = t1843 * t1960 - t1846 * t1968;
t1970 = t1851 * t1749;
t2041 = t1854 * t1858 - t1970;
t2071 = -t1687 * t1852 + t1855 * t2041;
t1969 = t1854 * t1749;
t2042 = t1851 * t1858 + t1969;
t2072 = t1687 * t1855 + t1852 * t2041;
t2090 = -t1845 * t2042 + t1848 * t2072;
t2121 = -t1844 * t2090 + t1847 * t2071;
t1769 = t1846 * t1776;
t1700 = t1769 + t1959;
t1871 = (t1792 * t1846 - t1794 * t1843) * t1982;
t1994 = t1804 * t1851;
t2040 = t1854 * t1871 - t1994;
t2069 = t1700 * t1852 + t1855 * t2040;
t1787 = t1854 * t1804;
t2045 = t1851 * t1871 + t1787;
t2070 = -t1700 * t1855 + t1852 * t2040;
t2091 = -t1845 * t2045 + t1848 * t2070;
t2120 = -t1844 * t2091 + t1847 * t2069;
t2119 = t1844 * t2071 + t1847 * t2090;
t2118 = t1844 * t2069 + t1847 * t2091;
t2110 = pkin(3) * t2074;
t2109 = qJ(4) * t2074;
t1691 = t1768 * t1843 - t1769;
t1692 = t1768 * t1846 + t1776 * t1843;
t1947 = t1854 * t1692 + t1970;
t2098 = -t1691 * t1855 + t1852 * t1947;
t2097 = pkin(3) * t2076 + qJ(4) * t2073;
t2093 = t1845 * t2070 + t1848 * t2045;
t2092 = t1845 * t2072 + t1848 * t2042;
t1850 = sin(qJ(6));
t1853 = cos(qJ(6));
t1733 = -t1853 * t1792 + t1794 * t1850;
t1735 = t1792 * t1850 + t1794 * t1853;
t1661 = t1735 * t1733;
t1797 = qJDD(6) + t1804;
t2059 = -t1661 + t1797;
t2083 = t1850 * t2059;
t2079 = t1853 * t2059;
t2027 = g(3) - qJDD(1);
t1965 = t1845 * t2027;
t1966 = g(1) * t1844 - t1847 * g(2);
t2075 = t1848 * t1966 - t1965;
t1948 = t1851 * t1692 - t1969;
t2036 = -t1845 * t1948 + t1848 * t2098;
t2039 = t1691 * t1852 + t1855 * t1947;
t2068 = t1844 * t2039 + t1847 * t2036;
t2067 = -t1844 * t2036 + t1847 * t2039;
t1811 = g(1) * t1847 + g(2) * t1844;
t1717 = -t1852 * t1811 - t1855 * t2075;
t1723 = -t1855 * t1811 + t1852 * t2075;
t1626 = t1717 * t1852 + t1723 * t1855;
t1604 = t1626 * t1845;
t2066 = t1844 * t2027;
t2065 = t1847 * t2027;
t1795 = t1845 * t1966;
t1941 = t1848 * t2027 + t1795;
t1892 = t1852 * t1941;
t1891 = t1855 * t1941;
t1694 = -t1856 * pkin(2) + qJDD(2) * pkin(8) + t1723;
t2031 = pkin(3) * t1854;
t1953 = -qJ(4) * t1851 - t2031;
t1801 = t1953 * qJD(2);
t2060 = (qJD(2) * t1801 + t1694) * t1851;
t1651 = t1854 * t1694 - t1851 * t1941;
t2035 = qJD(3) ^ 2;
t1603 = -pkin(3) * t2035 + qJDD(3) * qJ(4) + t1801 * t1982 + t1651;
t1693 = -qJDD(2) * pkin(2) - t1856 * pkin(8) + t1717;
t1950 = t1803 + t1827;
t1622 = -t1950 * qJ(4) + (-t1804 + t1826) * pkin(3) + t1693;
t1962 = t1843 * t1603 - t1846 * t1622;
t1890 = t1804 * pkin(4) - qJ(5) * t1835 + qJDD(5) + t1962;
t1740 = pkin(4) * t1792 - qJ(5) * t1794;
t1981 = (2 * qJD(4)) + t1740;
t1442 = t1804 * pkin(5) + t1949 * pkin(9) + (pkin(5) * t1792 + t1981) * t1794 + t1890;
t1979 = t1846 * t1603 + t1843 * t1622;
t1939 = -pkin(4) * t1835 - t1804 * qJ(5) - t1792 * t1740 + t1979;
t2021 = qJD(4) * t1792;
t1783 = -0.2e1 * t2021;
t1973 = qJD(5) * t1982;
t2055 = t1783 - 0.2e1 * t1973;
t1480 = t1939 + t2055;
t1900 = pkin(5) * t1982 - pkin(9) * t1794;
t1461 = -pkin(5) * t2034 + pkin(9) * t1960 - t1900 * t1982 + t1480;
t1381 = -t1853 * t1442 + t1850 * t1461;
t1382 = t1850 * t1442 + t1853 * t1461;
t1338 = -t1853 * t1381 + t1382 * t1850;
t1624 = -t1733 * qJD(6) + t1853 * t1768 + t1850 * t1960;
t1822 = qJD(6) + t1982;
t1708 = t1822 * t1733;
t2058 = -t1708 + t1624;
t1762 = t1854 * t1941;
t1956 = -qJDD(3) * pkin(3) - qJ(4) * t2035 + qJDD(4) + t1762;
t1886 = t1768 * qJ(5) - t1956 + t2132;
t1983 = t1851 * t1694;
t2051 = -(qJ(5) * t1792 * t1854 - t1801 * t1851) * qJD(2) - t1886 + t1983;
t1339 = t1850 * t1381 + t1853 * t1382;
t2033 = pkin(4) + pkin(5);
t2050 = qJ(5) * t1339 - t1338 * t2033;
t1961 = t1850 * t1768 - t1853 * t1960;
t1565 = (qJD(6) - t1822) * t1735 + t1961;
t1569 = t1708 + t1624;
t1494 = -t1565 * t1850 - t1569 * t1853;
t1496 = -t1565 * t1853 + t1569 * t1850;
t2049 = qJ(5) * t1496 - t1494 * t2033;
t1964 = t1981 * t1794;
t1491 = t1964 + t1890;
t1402 = t1480 * t1846 + t1491 * t1843;
t2019 = qJD(5) * t1794;
t1524 = -0.2e1 * t2019 + t2051;
t1974 = pkin(4) * t1846 + pkin(3);
t2048 = -t1524 * (qJ(5) * t1843 + t1974) + qJ(4) * t1402;
t1729 = t1733 ^ 2;
t1820 = t1822 ^ 2;
t1646 = -t1820 - t1729;
t1541 = t1646 * t1850 + t2079;
t1542 = t1646 * t1853 - t2083;
t2047 = qJ(5) * t1542 - t1541 * t2033 + t1381;
t1730 = t1735 ^ 2;
t1684 = -t1730 - t1820;
t1628 = t1661 + t1797;
t2015 = t1628 * t1850;
t1546 = t1684 * t1853 - t2015;
t2014 = t1628 * t1853;
t1547 = -t1684 * t1850 - t2014;
t2046 = qJ(5) * t1547 - t1546 * t2033 + t1382;
t2044 = -t1847 * t1811 - t1844 * t1966;
t2043 = -t1844 * t1811 + t1847 * t1966;
t1782 = 0.2e1 * t2019;
t1498 = t1782 - t2060 + (t1867 + t1968) * qJ(5) + t1886;
t2038 = t1843 * t1498 + t1867 * t1974 - t2151;
t2037 = t1845 * t2098 + t1848 * t1948;
t2032 = pkin(3) * t1851;
t2030 = pkin(7) * t1848;
t2029 = pkin(9) * t1338;
t2028 = pkin(9) * t1339;
t2023 = qJD(2) * qJD(3);
t2020 = qJD(4) * t1794;
t1602 = t1956 + t2060;
t2017 = t1602 * t1843;
t2016 = t1602 * t1846;
t2011 = t1693 * t1851;
t2010 = t1693 * t1854;
t1997 = t1735 * t1822;
t1821 = t1851 * t1856 * t1854;
t1812 = -t1821 + qJDD(3);
t1993 = t1812 * t1851;
t1992 = t1812 * t1854;
t1813 = qJDD(3) + t1821;
t1991 = t1813 * t1851;
t1990 = t1822 * t1850;
t1989 = t1822 * t1853;
t1839 = t1851 ^ 2;
t1988 = t1839 * t1856;
t1492 = pkin(5) * t1960 + pkin(9) * t2034 - t1794 * t1900 + t1524;
t1984 = t1850 * t1492;
t1486 = t1853 * t1492;
t1978 = t1839 + t1840;
t1977 = t1845 * qJDD(2);
t1972 = t1851 * t1661;
t1971 = t1854 * t1661;
t1522 = t1783 + t1979;
t1967 = -pkin(9) * t1546 - t1486;
t1521 = t1962 + 0.2e1 * t2020;
t1419 = t1521 * t1843 + t1846 * t1522;
t1650 = t1762 + t1983;
t1545 = t1650 * t1851 + t1854 * t1651;
t1958 = t1852 * t1821;
t1957 = t1855 * t1821;
t1954 = -pkin(3) * t1602 + qJ(4) * t1419;
t1952 = -pkin(4) * t1491 + qJ(5) * t1480;
t1951 = pkin(4) * t1949 + qJ(5) * t1712;
t1544 = t1650 * t1854 - t1651 * t1851;
t1806 = t1978 * qJDD(2);
t1809 = t1835 + t1988;
t1745 = t1806 * t1855 - t1809 * t1852;
t1945 = pkin(7) * t1745 + t1544 * t1852;
t1944 = -pkin(9) * t1541 - t1984;
t1943 = -pkin(9) * t1542 - t1486;
t1942 = -pkin(9) * t1547 + t1984;
t1308 = t1338 * t1843 + t1339 * t1846;
t1299 = t1308 * t1854 + t1492 * t1851;
t1307 = -t1338 * t1846 + t1339 * t1843;
t1938 = t1299 * t1852 - t1307 * t1855;
t1375 = t1402 * t1854 + t1524 * t1851;
t1401 = t1480 * t1843 - t1491 * t1846;
t1937 = t1375 * t1852 - t1401 * t1855;
t1408 = t1494 * t1843 + t1496 * t1846;
t1594 = -t1729 - t1730;
t1390 = t1408 * t1854 - t1594 * t1851;
t1406 = -t1494 * t1846 + t1496 * t1843;
t1936 = t1390 * t1852 - t1406 * t1855;
t1564 = (qJD(6) + t1822) * t1735 + t1961;
t1493 = -t1564 * t1850 + t1853 * t2058;
t1495 = -t1564 * t1853 - t1850 * t2058;
t1407 = t1493 * t1843 + t1495 * t1846;
t1652 = t1730 - t1729;
t1397 = t1407 * t1854 - t1652 * t1851;
t1405 = -t1493 * t1846 + t1495 * t1843;
t1935 = t1397 * t1852 - t1405 * t1855;
t1410 = t1419 * t1854 + t1602 * t1851;
t1418 = -t1521 * t1846 + t1522 * t1843;
t1934 = t1410 * t1852 - t1418 * t1855;
t1460 = t1541 * t1843 + t1542 * t1846;
t1414 = t1460 * t1854 - t1564 * t1851;
t1459 = -t1541 * t1846 + t1542 * t1843;
t1933 = t1414 * t1852 - t1459 * t1855;
t1479 = t1546 * t1843 + t1547 * t1846;
t1417 = t1479 * t1854 - t1851 * t2058;
t1478 = -t1546 * t1846 + t1547 * t1843;
t1932 = t1417 * t1852 - t1478 * t1855;
t1699 = -t1730 + t1820;
t1560 = t1699 * t1853 + t2083;
t1562 = -t1699 * t1850 + t2079;
t1489 = t1560 * t1843 + t1562 * t1846;
t1424 = t1489 * t1854 - t1569 * t1851;
t1487 = -t1560 * t1846 + t1562 * t1843;
t1931 = t1424 * t1852 - t1487 * t1855;
t1698 = t1729 - t1820;
t1561 = t1698 * t1850 + t2014;
t1563 = t1698 * t1853 - t2015;
t1490 = t1561 * t1843 + t1563 * t1846;
t1425 = t1490 * t1854 + t1565 * t1851;
t1488 = -t1561 * t1846 + t1563 * t1843;
t1930 = t1425 * t1852 - t1488 * t1855;
t1623 = -qJD(6) * t1735 - t1961;
t1556 = -t1623 * t1853 - t1733 * t1990;
t1557 = -t1623 * t1850 + t1733 * t1989;
t1483 = -t1556 * t1843 + t1557 * t1846;
t1446 = t1483 * t1854 + t1972;
t1481 = t1556 * t1846 + t1557 * t1843;
t1929 = t1446 * t1852 - t1481 * t1855;
t1558 = t1624 * t1850 + t1735 * t1989;
t1559 = t1624 * t1853 - t1735 * t1990;
t1484 = t1558 * t1843 + t1559 * t1846;
t1447 = t1484 * t1854 - t1972;
t1482 = -t1558 * t1846 + t1559 * t1843;
t1928 = t1447 * t1852 - t1482 * t1855;
t1608 = (-t1733 * t1850 - t1735 * t1853) * t1822;
t1609 = (-t1733 * t1853 + t1735 * t1850) * t1822;
t1532 = t1608 * t1843 + t1609 * t1846;
t1519 = t1532 * t1854 - t1797 * t1851;
t1531 = -t1608 * t1846 + t1609 * t1843;
t1927 = t1519 * t1852 - t1531 * t1855;
t1926 = t1545 * t1852 - t1693 * t1855;
t1625 = t1855 * t1717 - t1852 * t1723;
t1802 = 0.2e1 * t1827 + t1976;
t1805 = -0.2e1 * t1826 + t1975;
t1742 = -t1802 * t1851 + t1805 * t1854;
t1810 = -t1835 + t1988;
t1909 = t1742 * t1852 - t1810 * t1855;
t1819 = -t1835 - t2035;
t1758 = t1819 * t1854 - t1991;
t1908 = t1758 * t1852 + t1805 * t1855;
t1817 = -t1988 - t2035;
t1760 = -t1817 * t1851 - t1992;
t1907 = t1760 * t1852 - t1802 * t1855;
t1903 = qJDD(2) * t1852 + t1855 * t1856;
t1780 = t1903 * t1848;
t1807 = qJDD(2) * t1855 - t1852 * t1856;
t1906 = t1780 * t1847 + t1807 * t1844;
t1905 = t1780 * t1844 - t1807 * t1847;
t1904 = t1806 * t1852 + t1809 * t1855;
t1799 = t1978 * t2023;
t1902 = -qJDD(3) * t1855 + t1799 * t1852;
t1899 = -pkin(9) * t1494 - t1338;
t1898 = -pkin(9) * t1496 - t1339;
t1896 = -t2016 + t2097;
t1818 = t1835 - t2035;
t1757 = t1818 * t1854 - t1993;
t1894 = t1757 * t1852 - t1855 * t1975;
t1800 = t1854 * t1813;
t1816 = -t1988 + t2035;
t1759 = -t1816 * t1851 + t1800;
t1893 = t1759 * t1852 - t1855 * t1976;
t1765 = -t1840 * t2023 - t1994;
t1889 = t1765 * t1852 - t1957;
t1766 = t1803 * t1854 - t1839 * t2023;
t1888 = t1766 * t1852 + t1957;
t1887 = -pkin(3) * t1867 + t2017 + t2151;
t1281 = -pkin(3) * t1307 - t2050;
t1321 = -t1492 * t2033 - t2028;
t1330 = -qJ(5) * t1492 - t2029;
t1285 = -qJ(4) * t1307 - t1321 * t1843 + t1330 * t1846;
t1298 = t1308 * t1851 - t1492 * t1854;
t1267 = -pkin(8) * t1298 - t1281 * t1851 + t1285 * t1854;
t1865 = -pkin(3) * t1492 + qJ(4) * t1308 + t1321 * t1846 + t1330 * t1843;
t1275 = -pkin(2) * t1298 - t1865;
t1287 = t1299 * t1855 + t1307 * t1852;
t1885 = pkin(7) * t1287 + t1267 * t1852 + t1275 * t1855;
t1323 = t1594 * t2033 + t1898;
t1327 = qJ(5) * t1594 + t1899;
t1290 = -qJ(4) * t1406 - t1323 * t1843 + t1327 * t1846;
t1342 = -pkin(3) * t1406 - t2049;
t1389 = t1408 * t1851 + t1594 * t1854;
t1286 = -pkin(8) * t1389 + t1290 * t1854 - t1342 * t1851;
t1864 = pkin(3) * t1594 + qJ(4) * t1408 + t1323 * t1846 + t1327 * t1843;
t1289 = -pkin(2) * t1389 - t1864;
t1346 = t1390 * t1855 + t1406 * t1852;
t1884 = pkin(7) * t1346 + t1286 * t1852 + t1289 * t1855;
t1383 = t1564 * t2033 + t1943;
t1399 = qJ(5) * t1564 + t1944;
t1332 = -qJ(4) * t1459 - t1383 * t1843 + t1399 * t1846;
t1333 = -pkin(3) * t1459 - t2047;
t1413 = t1460 * t1851 + t1564 * t1854;
t1294 = -pkin(8) * t1413 + t1332 * t1854 - t1333 * t1851;
t1863 = pkin(3) * t1564 + qJ(4) * t1460 + t1383 * t1846 + t1399 * t1843;
t1312 = -pkin(2) * t1413 - t1863;
t1373 = t1414 * t1855 + t1459 * t1852;
t1883 = pkin(7) * t1373 + t1294 * t1852 + t1312 * t1855;
t1388 = t2033 * t2058 + t1942;
t1403 = qJ(5) * t2058 + t1967;
t1335 = -qJ(4) * t1478 - t1388 * t1843 + t1403 * t1846;
t1340 = -pkin(3) * t1478 - t2046;
t1416 = t1479 * t1851 + t1854 * t2058;
t1297 = -pkin(8) * t1416 + t1335 * t1854 - t1340 * t1851;
t1862 = pkin(3) * t2058 + qJ(4) * t1479 + t1388 * t1846 + t1403 * t1843;
t1317 = -pkin(2) * t1416 - t1862;
t1378 = t1417 * t1855 + t1478 * t1852;
t1882 = pkin(7) * t1378 + t1297 * t1852 + t1317 * t1855;
t1353 = -pkin(3) * t1401 - t1952;
t1366 = -qJ(4) * t1401 + (pkin(4) * t1843 - qJ(5) * t1846) * t1524;
t1374 = t1402 * t1851 - t1524 * t1854;
t1306 = -pkin(8) * t1374 - t1353 * t1851 + t1366 * t1854;
t1322 = -pkin(2) * t1374 - t2048;
t1343 = t1375 * t1855 + t1401 * t1852;
t1881 = pkin(7) * t1343 + t1306 * t1852 + t1322 * t1855;
t1409 = t1419 * t1851 - t1602 * t1854;
t1341 = -pkin(8) * t1409 + (-qJ(4) * t1854 + t2032) * t1418;
t1359 = -pkin(2) * t1409 - t1954;
t1367 = t1410 * t1855 + t1418 * t1852;
t1880 = pkin(7) * t1367 + t1341 * t1852 + t1359 * t1855;
t1463 = pkin(4) * t2053 + t1480;
t1464 = qJ(5) * t2053 + t1491;
t1377 = -t1463 * t1843 + t1464 * t1846 - t2130;
t1540 = -t1951 - t2133;
t1360 = t1377 * t1854 - t1540 * t1851 - t2153;
t1861 = t1463 * t1846 + t1464 * t1843 + t2122;
t1370 = -t1861 - t2157;
t1879 = t1360 * t1852 + t1370 * t1855 + t2172;
t1860 = pkin(4) * t1739 - qJ(5) * t2056 + t1939;
t1415 = -t1860 + 0.2e1 * t1973 + 0.2e1 * t2021 - t2156;
t1438 = -pkin(4) * t2004 + t1498 * t1846 - t2152;
t1369 = -t1415 * t1851 + t1438 * t1854 + t2185;
t1404 = -t2038 + t2186;
t1878 = t1369 * t1852 + t1404 * t1855 - t2195;
t1857 = -pkin(4) * t2057 - qJ(5) * t2052 + t1890;
t1430 = t1857 + t1964 - t2110;
t1499 = t1782 - t2051 + t2132;
t1441 = qJ(5) * t2006 - t1499 * t1843 - t2109;
t1371 = -t1430 * t1851 + t1441 * t1854 - t2131;
t1866 = qJ(5) * t2007 + t1499 * t1846 + t2097;
t1411 = -t1866 - t2134;
t1877 = t1371 * t1852 + t1411 * t1855 + t2155;
t1412 = -t1418 - t2130;
t1376 = t1412 * t1854 + t2032 * t2095 - t2153;
t1870 = t1419 + t2122;
t1384 = -t1870 - t2157;
t1876 = t1376 * t1852 + t1384 * t1855 + t2172;
t1477 = t1521 - t2110;
t1533 = t2017 - t2109;
t1398 = -t1477 * t1851 + t1533 * t1854 - t2131;
t1443 = -t1896 - t2134;
t1875 = t1398 * t1852 + t1443 * t1855 + t2155;
t1485 = t1522 + t2156;
t1539 = t2016 + t2152;
t1400 = -t1485 * t1851 + t1539 * t1854 - t2185;
t1462 = -t1887 - t2186;
t1874 = t1400 * t1852 + t1462 * t1855 + t2195;
t1754 = t1819 * t1851 + t1800;
t1600 = -pkin(2) * t1754 + t1650;
t1648 = -pkin(8) * t1754 + t2011;
t1696 = t1758 * t1855 - t1805 * t1852;
t1873 = pkin(7) * t1696 + t1600 * t1855 + t1648 * t1852;
t1756 = t1817 * t1854 - t1993;
t1601 = -pkin(2) * t1756 + t1651;
t1649 = -pkin(8) * t1756 + t2010;
t1697 = t1760 * t1855 + t1802 * t1852;
t1872 = pkin(7) * t1697 + t1601 * t1855 + t1649 * t1852;
t1520 = t1545 * t1855 + t1693 * t1852;
t1869 = pkin(7) * t1520 - (-pkin(2) * t1855 - pkin(8) * t1852) * t1544;
t1831 = t1848 * qJDD(2);
t1781 = t1807 * t1848;
t1779 = t1807 * t1845;
t1778 = t1903 * t1845;
t1767 = qJDD(3) * t1852 + t1799 * t1855;
t1755 = t1816 * t1854 + t1991;
t1753 = t1818 * t1851 + t1992;
t1752 = t1950 * t1851;
t1751 = -t1827 * t1851 + t1787;
t1747 = t1902 * t1848;
t1746 = t1902 * t1845;
t1741 = t1802 * t1854 + t1805 * t1851;
t1737 = t1904 * t1848;
t1736 = t1904 * t1845;
t1732 = -t1781 * t1844 - t1847 * t1903;
t1731 = t1781 * t1847 - t1844 * t1903;
t1727 = t1766 * t1855 - t1958;
t1726 = t1765 * t1855 + t1958;
t1725 = t1759 * t1855 + t1852 * t1976;
t1724 = t1757 * t1855 + t1852 * t1975;
t1686 = t1742 * t1855 + t1810 * t1852;
t1675 = -t1891 + (t1778 * t1845 + t1780 * t1848) * pkin(7);
t1674 = -t1892 + (-t1779 * t1845 - t1781 * t1848) * pkin(7);
t1660 = -t1845 * t1752 + t1848 * t1888;
t1659 = -t1845 * t1751 + t1848 * t1889;
t1658 = t1848 * t1752 + t1845 * t1888;
t1657 = t1848 * t1751 + t1845 * t1889;
t1656 = -t1845 * t1755 + t1848 * t1893;
t1655 = -t1845 * t1753 + t1848 * t1894;
t1654 = t1848 * t1755 + t1845 * t1893;
t1653 = t1848 * t1753 + t1845 * t1894;
t1645 = -t1845 * t1756 + t1848 * t1907;
t1644 = -t1845 * t1754 + t1848 * t1908;
t1643 = t1848 * t1756 + t1845 * t1907;
t1642 = t1848 * t1754 + t1845 * t1908;
t1613 = -t1845 * t1741 + t1848 * t1909;
t1612 = t1848 * t1741 + t1845 * t1909;
t1611 = pkin(2) * t1805 + pkin(8) * t1758 - t2010;
t1610 = -pkin(2) * t1802 + pkin(8) * t1760 + t2011;
t1605 = t1626 * t1848;
t1598 = -pkin(1) * t1779 + t1845 * t1717 + t1848 * t1891 - t1903 * t2030;
t1597 = pkin(1) * t1778 + t1845 * t1723 - t1807 * t2030 - t1848 * t1892;
t1596 = pkin(1) * t1781 - t1848 * t1717 + (-pkin(7) * t1903 + t1891) * t1845;
t1595 = -pkin(1) * t1780 - t1848 * t1723 + (-pkin(7) * t1807 - t1892) * t1845;
t1579 = t1845 * t1795 + (t1965 - t1625) * t1848;
t1578 = -t1625 * t1845 - t1848 * t1941;
t1538 = pkin(2) * t1809 + pkin(8) * t1806 + t1545;
t1523 = -pkin(2) * t1693 + pkin(8) * t1545;
t1518 = t1532 * t1851 + t1797 * t1854;
t1517 = pkin(1) * t1579 + pkin(7) * t1604;
t1516 = -pkin(1) * t1578 + t1626 * t2030;
t1500 = t1855 * t1544 + (-t1736 * t1845 - t1737 * t1848) * pkin(7);
t1497 = (-t1578 * t1845 - t1579 * t1848) * pkin(7);
t1457 = -t1852 * t1601 + t1855 * t1649 + (-t1643 * t1845 - t1645 * t1848) * pkin(7);
t1456 = -t1852 * t1600 + t1855 * t1648 + (-t1642 * t1845 - t1644 * t1848) * pkin(7);
t1445 = t1484 * t1851 + t1971;
t1444 = t1483 * t1851 - t1971;
t1437 = t1544 * t1845 + t1848 * t1926;
t1436 = -t1544 * t1848 + t1845 * t1926;
t1431 = t1519 * t1855 + t1531 * t1852;
t1429 = -pkin(1) * t1643 - t1845 * t1610 + t1848 * t1872;
t1428 = -pkin(1) * t1642 - t1845 * t1611 + t1848 * t1873;
t1427 = pkin(1) * t1645 + t1848 * t1610 + t1845 * t1872;
t1426 = pkin(1) * t1644 + t1848 * t1611 + t1845 * t1873;
t1423 = t1490 * t1851 - t1565 * t1854;
t1422 = t1489 * t1851 + t1569 * t1854;
t1421 = -pkin(1) * t1736 - t1845 * t1538 + t1848 * t1945;
t1420 = pkin(1) * t1737 + t1848 * t1538 + t1845 * t1945;
t1396 = t1407 * t1851 + t1652 * t1854;
t1395 = -t1845 * t1518 + t1848 * t1927;
t1394 = t1848 * t1518 + t1845 * t1927;
t1393 = t1447 * t1855 + t1482 * t1852;
t1392 = t1446 * t1855 + t1481 * t1852;
t1391 = t1485 * t1854 + t1539 * t1851 - t2180;
t1387 = t1425 * t1855 + t1488 * t1852;
t1386 = t1424 * t1855 + t1487 * t1852;
t1385 = t1477 * t1854 + t1533 * t1851 + t2123;
t1372 = t2154 + t1851 * t1412 + (-pkin(2) - t2031) * t2095;
t1368 = t1430 * t1854 + t1441 * t1851 + t2123;
t1365 = -t1845 * t1445 + t1848 * t1928;
t1364 = -t1845 * t1444 + t1848 * t1929;
t1363 = t1848 * t1445 + t1845 * t1928;
t1362 = t1848 * t1444 + t1845 * t1929;
t1361 = t1415 * t1854 + t1438 * t1851 + t2180;
t1358 = -(pkin(2) * t1852 - pkin(8) * t1855) * t1544 + (-t1436 * t1845 - t1437 * t1848) * pkin(7);
t1357 = -t1845 * t1423 + t1848 * t1930;
t1356 = -t1845 * t1422 + t1848 * t1931;
t1355 = t1848 * t1423 + t1845 * t1930;
t1354 = t1848 * t1422 + t1845 * t1931;
t1352 = -pkin(1) * t1436 - t1845 * t1523 + t1848 * t1869;
t1351 = pkin(1) * t1437 + t1848 * t1523 + t1845 * t1869;
t1350 = -pkin(2) * t2095 + t1377 * t1851 + t1540 * t1854 + t2154;
t1349 = -t1845 * t1416 + t1848 * t1932;
t1348 = t1848 * t1416 + t1845 * t1932;
t1347 = t1397 * t1855 + t1405 * t1852;
t1345 = -t1845 * t1413 + t1848 * t1933;
t1344 = t1848 * t1413 + t1845 * t1933;
t1337 = -t1845 * t1409 + t1848 * t1934;
t1336 = t1848 * t1409 + t1845 * t1934;
t1334 = t1855 * t1400 - t1852 * t1462 - t2207;
t1331 = t1855 * t1398 - t1852 * t1443 + t2159;
t1329 = -t1845 * t1396 + t1848 * t1935;
t1328 = t1848 * t1396 + t1845 * t1935;
t1326 = pkin(8) * t1410 + (-pkin(2) + t1953) * t1418;
t1325 = -t1845 * t1389 + t1848 * t1936;
t1324 = t1848 * t1389 + t1845 * t1936;
t1320 = t1855 * t1371 - t1852 * t1411 + t2159;
t1319 = -t1845 * t1391 + t1848 * t1874 - t2204;
t1318 = t1848 * t1391 + t1845 * t1874 + t2203;
t1316 = -t1845 * t1374 + t1848 * t1937;
t1315 = t1848 * t1374 + t1845 * t1937;
t1314 = t1855 * t1369 - t1852 * t1404 + t2207;
t1313 = t1855 * t1376 - t1852 * t1384 + t2175;
t1311 = -t1845 * t1385 + t1848 * t1875 - t2173;
t1310 = t1848 * t1385 + t1845 * t1875 + t2174;
t1309 = t1855 * t1360 - t1852 * t1370 + t2175;
t1305 = -t1845 * t1368 + t1848 * t1877 - t2173;
t1304 = t1848 * t1368 + t1845 * t1877 + t2174;
t1303 = -t1845 * t1372 + t1848 * t1876 - t2187;
t1302 = t1848 * t1372 + t1845 * t1876 + t2188;
t1301 = -t1845 * t1361 + t1848 * t1878 + t2204;
t1300 = t1848 * t1361 + t1845 * t1878 - t2203;
t1296 = -pkin(2) * t1401 + pkin(8) * t1375 + t1353 * t1854 + t1366 * t1851;
t1295 = -pkin(2) * t1478 + pkin(8) * t1417 + t1335 * t1851 + t1340 * t1854;
t1293 = -t1845 * t1350 + t1848 * t1879 - t2187;
t1292 = t1848 * t1350 + t1845 * t1879 + t2188;
t1291 = -pkin(2) * t1459 + pkin(8) * t1414 + t1332 * t1851 + t1333 * t1854;
t1288 = t1855 * t1341 - t1852 * t1359 + (-t1336 * t1845 - t1337 * t1848) * pkin(7);
t1284 = -pkin(1) * t1336 - t1845 * t1326 + t1848 * t1880;
t1283 = pkin(1) * t1337 + t1848 * t1326 + t1845 * t1880;
t1282 = -pkin(2) * t1406 + pkin(8) * t1390 + t1290 * t1851 + t1342 * t1854;
t1280 = t1855 * t1297 - t1852 * t1317 + (-t1348 * t1845 - t1349 * t1848) * pkin(7);
t1279 = t1855 * t1294 - t1852 * t1312 + (-t1344 * t1845 - t1345 * t1848) * pkin(7);
t1278 = t1855 * t1306 - t1852 * t1322 + (-t1315 * t1845 - t1316 * t1848) * pkin(7);
t1277 = -t1845 * t1298 + t1848 * t1938;
t1276 = t1848 * t1298 + t1845 * t1938;
t1274 = -pkin(1) * t1348 - t1845 * t1295 + t1848 * t1882;
t1273 = pkin(1) * t1349 + t1848 * t1295 + t1845 * t1882;
t1272 = -pkin(1) * t1315 - t1845 * t1296 + t1848 * t1881;
t1271 = pkin(1) * t1316 + t1848 * t1296 + t1845 * t1881;
t1270 = -pkin(1) * t1344 - t1845 * t1291 + t1848 * t1883;
t1269 = pkin(1) * t1345 + t1848 * t1291 + t1845 * t1883;
t1268 = t1855 * t1286 - t1852 * t1289 + (-t1324 * t1845 - t1325 * t1848) * pkin(7);
t1266 = -pkin(2) * t1307 + pkin(8) * t1299 + t1281 * t1854 + t1285 * t1851;
t1265 = -pkin(1) * t1324 - t1845 * t1282 + t1848 * t1884;
t1264 = pkin(1) * t1325 + t1848 * t1282 + t1845 * t1884;
t1263 = t1855 * t1267 - t1852 * t1275 + (-t1276 * t1845 - t1277 * t1848) * pkin(7);
t1262 = -pkin(1) * t1276 - t1845 * t1266 + t1848 * t1885;
t1261 = pkin(1) * t1277 + t1848 * t1266 + t1845 * t1885;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t2066, -t2065, -t2043, -qJ(1) * t2043, 0, 0, -t1905, 0, t1732, t1844 * t1977, -qJ(1) * t1731 - t1598 * t1844 + t1674 * t1847, qJ(1) * t1906 - t1844 * t1597 + t1847 * t1675, -t1605 * t1844 + t1625 * t1847, t1847 * t1497 - t1844 * t1516 - qJ(1) * (t1579 * t1847 + t1626 * t1844), -t1660 * t1844 + t1727 * t1847, -t1613 * t1844 + t1686 * t1847, -t1656 * t1844 + t1725 * t1847, -t1659 * t1844 + t1726 * t1847, -t1655 * t1844 + t1724 * t1847, -t1747 * t1844 + t1767 * t1847, t1847 * t1456 - t1844 * t1428 - qJ(1) * (t1644 * t1847 + t1696 * t1844), t1847 * t1457 - t1844 * t1429 - qJ(1) * (t1645 * t1847 + t1697 * t1844), t1847 * t1500 - t1844 * t1421 - qJ(1) * (t1737 * t1847 + t1745 * t1844), t1847 * t1358 - t1844 * t1352 - qJ(1) * (t1437 * t1847 + t1520 * t1844), t2067, -t2197, t2179, t2121, -t2200, t2120, -t1844 * t1311 + t1847 * t1331 - t2160, -t1844 * t1319 + t1847 * t1334 - t2206, -t1844 * t1303 + t1847 * t1313 - t2176, t1847 * t1288 - t1844 * t1284 - qJ(1) * (t1337 * t1847 + t1367 * t1844), t2067, t2179, t2197, t2120, t2200, t2121, -t1844 * t1305 + t1847 * t1320 - t2160, -t1844 * t1293 + t1847 * t1309 - t2176, -t1844 * t1301 + t1847 * t1314 + t2206, t1847 * t1278 - t1844 * t1272 - qJ(1) * (t1316 * t1847 + t1343 * t1844), -t1365 * t1844 + t1393 * t1847, -t1329 * t1844 + t1347 * t1847, -t1356 * t1844 + t1386 * t1847, -t1364 * t1844 + t1392 * t1847, -t1357 * t1844 + t1387 * t1847, -t1395 * t1844 + t1431 * t1847, t1847 * t1279 - t1844 * t1270 - qJ(1) * (t1345 * t1847 + t1373 * t1844), t1847 * t1280 - t1844 * t1274 - qJ(1) * (t1349 * t1847 + t1378 * t1844), t1847 * t1268 - t1844 * t1265 - qJ(1) * (t1325 * t1847 + t1346 * t1844), t1847 * t1263 - t1844 * t1262 - qJ(1) * (t1277 * t1847 + t1287 * t1844); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t2065, -t2066, t2044, qJ(1) * t2044, 0, 0, t1906, 0, t1731, -t1847 * t1977, qJ(1) * t1732 + t1598 * t1847 + t1674 * t1844, qJ(1) * t1905 + t1847 * t1597 + t1844 * t1675, t1605 * t1847 + t1625 * t1844, t1844 * t1497 + t1847 * t1516 + qJ(1) * (-t1579 * t1844 + t1626 * t1847), t1660 * t1847 + t1727 * t1844, t1613 * t1847 + t1686 * t1844, t1656 * t1847 + t1725 * t1844, t1659 * t1847 + t1726 * t1844, t1655 * t1847 + t1724 * t1844, t1747 * t1847 + t1767 * t1844, t1844 * t1456 + t1847 * t1428 + qJ(1) * (-t1644 * t1844 + t1696 * t1847), t1844 * t1457 + t1847 * t1429 + qJ(1) * (-t1645 * t1844 + t1697 * t1847), t1844 * t1500 + t1847 * t1421 + qJ(1) * (-t1737 * t1844 + t1745 * t1847), t1844 * t1358 + t1847 * t1352 + qJ(1) * (-t1437 * t1844 + t1520 * t1847), t2068, t2196, t2178, t2119, t2199, t2118, t1847 * t1311 + t1844 * t1331 + t2161, t1847 * t1319 + t1844 * t1334 - t2205, t1847 * t1303 + t1844 * t1313 + t2177, t1844 * t1288 + t1847 * t1284 + qJ(1) * (-t1337 * t1844 + t1367 * t1847), t2068, t2178, -t2196, t2118, -t2199, t2119, t1847 * t1305 + t1844 * t1320 + t2161, t1847 * t1293 + t1844 * t1309 + t2177, t1847 * t1301 + t1844 * t1314 + t2205, t1844 * t1278 + t1847 * t1272 + qJ(1) * (-t1316 * t1844 + t1343 * t1847), t1365 * t1847 + t1393 * t1844, t1329 * t1847 + t1347 * t1844, t1356 * t1847 + t1386 * t1844, t1364 * t1847 + t1392 * t1844, t1357 * t1847 + t1387 * t1844, t1395 * t1847 + t1431 * t1844, t1844 * t1279 + t1847 * t1270 + qJ(1) * (-t1345 * t1844 + t1373 * t1847), t1844 * t1280 + t1847 * t1274 + qJ(1) * (-t1349 * t1844 + t1378 * t1847), t1844 * t1268 + t1847 * t1265 + qJ(1) * (-t1325 * t1844 + t1346 * t1847), t1844 * t1263 + t1847 * t1262 + qJ(1) * (-t1277 * t1844 + t1287 * t1847); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t1966, t1811, 0, 0, 0, 0, t1778, 0, t1779, t1831, t1596, t1595, t1604, t1517, t1658, t1612, t1654, t1657, t1653, t1746, t1426, t1427, t1420, t1351, t2037, t1453, t2165, t2092, t1470, t2093, t1310, t1318, t1302, t1283, t2037, t2165, -t1453, t2093, -t1470, t2092, t1304, t1292, t1300, t1271, t1363, t1328, t1354, t1362, t1355, t1394, t1269, t1273, t1264, t1261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2027, -t1966, 0, 0, 0, t1807, 0, -t1903, 0, t1674, t1675, t1625, t1497, t1727, t1686, t1725, t1726, t1724, t1767, t1456, t1457, t1500, t1358, t2039, t1512, t2135, t2071, t1530, t2069, t1331, t1334, t1313, t1288, t2039, t2135, -t1512, t2069, -t1530, t2071, t1320, t1309, t1314, t1278, t1393, t1347, t1386, t1392, t1387, t1431, t1279, t1280, t1268, t1263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2027, 0, -t1811, 0, 0, 0, t1780, 0, t1781, -t1977, t1598, t1597, t1605, t1516, t1660, t1613, t1656, t1659, t1655, t1747, t1428, t1429, t1421, t1352, t2036, t1455, t2163, t2090, t1476, t2091, t1311, t1319, t1303, t1284, t2036, t2163, -t1455, t2091, -t1476, t2090, t1305, t1293, t1301, t1272, t1365, t1329, t1356, t1364, t1357, t1395, t1270, t1274, t1265, t1262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1966, t1811, 0, 0, 0, 0, t1778, 0, t1779, t1831, t1596, t1595, t1604, t1517, t1658, t1612, t1654, t1657, t1653, t1746, t1426, t1427, t1420, t1351, t2037, t1453, t2165, t2092, t1470, t2093, t1310, t1318, t1302, t1283, t2037, t2165, -t1453, t2093, -t1470, t2092, t1304, t1292, t1300, t1271, t1363, t1328, t1354, t1362, t1355, t1394, t1269, t1273, t1264, t1261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t1856, 0, 0, -t1941, t1717, 0, t1766, t1742, t1759, t1765, t1757, t1799, t1648, t1649, t1544, pkin(8) * t1544, t1947, t1577, t2114, t2041, t1591, t2040, t1398, t1400, t1376, t1341, t1947, t2114, -t1577, t2040, -t1591, t2041, t1371, t1360, t1369, t1306, t1447, t1397, t1424, t1446, t1425, t1519, t1294, t1297, t1286, t1267; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1856, 0, qJDD(2), 0, t1941, 0, t1723, 0, t1821, -t1810, -t1976, -t1821, -t1975, -qJDD(3), t1600, t1601, 0, pkin(2) * t1544, -t1691, -t1616, -t2100, t1687, -t1667, -t1700, t1443, t1462, t1384, t1359, -t1691, -t2100, t1616, -t1700, t1667, t1687, t1411, t1370, t1404, t1322, -t1482, -t1405, -t1487, -t1481, -t1488, -t1531, t1312, t1317, t1289, t1275; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t1717, -t1723, 0, 0, t1752, t1741, t1755, t1751, t1753, 0, t1611, t1610, t1538, t1523, t1948, t1575, t2117, t2042, t1585, t2045, t1385, t1391, t1372, t1326, t1948, t2117, -t1575, t2045, -t1585, t2042, t1368, t1350, t1361, t1296, t1445, t1396, t1422, t1444, t1423, t1518, t1291, t1295, t1282, t1266; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1803, t1805, t1813, -t1827, t1818, t1827, 0, t1693, t1650, 0, t1692, t1620, t2099, t1858, t1673, t1871, t1533, t1539, t1412, -qJ(4) * t1418, t1692, t2099, -t1620, t1871, -t1673, t1858, t1441, t1377, t1438, t1366, t1484, t1407, t1489, t1483, t1490, t1532, t1332, t1335, t1290, t1285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1826, t1802, t1816, t1804, t1812, -t1826, -t1693, 0, t1651, 0, -t1749, -t2054, t1949, t1749, -t1712, t1804, t1477, t1485, -t2133, -pkin(3) * t1418, -t1749, t1949, t2054, t1804, t1712, t1749, t1430, t1540, t1415, t1353, t1661, t1652, t1569, -t1661, -t1565, t1797, t1333, t1340, t1342, t1281; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1821, t1810, t1976, t1821, t1975, qJDD(3), -t1650, -t1651, 0, 0, t1691, t1616, t2100, -t1687, t1667, t1700, t1896, t1887, t1870, t1954, t1691, t2100, -t1616, t1700, -t1667, -t1687, t1866, t1861, t2038, t2048, t1482, t1405, t1487, t1481, t1488, t1531, t1863, t1862, t1864, t1865; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1768, t2076, t2057, -t1968, t1772, t1968, 0, t1602, t1521, 0, t1768, t2057, -t2076, t1968, -t1772, -t1968, qJ(5) * t2076, t1464, t1498, -qJ(5) * t1524, t1559, t1495, t1562, t1557, t1563, t1609, t1399, t1403, t1327, t1330; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1776, t1867, -t1773, -t1960, -t2056, t1776, -t1602, 0, t1522, 0, -t1776, -t1773, -t1867, t1776, t2056, -t1960, t1499, t1463, pkin(4) * t1867, -pkin(4) * t1524, -t1558, -t1493, -t1560, t1556, -t1561, -t1608, t1383, t1388, t1323, t1321; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1749, t2054, -t1949, -t1749, t1712, -t1804, -t1521, -t1522, 0, 0, t1749, -t1949, -t2054, -t1804, -t1712, -t1749, -t1794 * t1740 - t1857 - 0.2e1 * t2020, t1951, t1860 + t2055, t1952, -t1661, -t1652, -t1569, t1661, t1565, -t1797, t2047, t2046, t2049, t2050; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1768, t2057, -t2076, t1968, -t1772, -t1968, 0, t1491, -t1524, 0, t1559, t1495, t1562, t1557, t1563, t1609, t1944, t1967, t1899, -t2029; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1749, -t1949, -t2054, -t1804, -t1712, -t1749, -t1491, 0, t1480, 0, -t1661, -t1652, -t1569, t1661, t1565, -t1797, -pkin(5) * t1541 + t1381, -pkin(5) * t1546 + t1382, -pkin(5) * t1494, -pkin(5) * t1338; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1776, t1773, t1867, -t1776, -t2056, t1960, t1524, -t1480, 0, 0, t1558, t1493, t1560, -t1556, t1561, t1608, -pkin(5) * t1564 - t1943, -pkin(5) * t2058 - t1942, -pkin(5) * t1594 - t1898, pkin(5) * t1492 + t2028; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1624, -t1564, t2059, t1708, t1698, -t1708, 0, -t1492, t1381, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1997, t2058, t1699, t1623, t1628, -t1997, t1492, 0, t1382, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1661, t1652, t1569, -t1661, -t1565, t1797, -t1381, -t1382, 0, 0;];
m_new_reg  = t1;