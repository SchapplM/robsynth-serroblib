% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6RPRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 21:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6RPRRPP1_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP1_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP1_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_invdynm_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:18:27
% EndTime: 2019-05-05 21:19:17
% DurationCPUTime: 53.58s
% Computational Cost: add. (179782->862), mult. (361623->1076), div. (0->0), fcn. (245008->10), ass. (0->589)
t1859 = sin(qJ(4));
t1862 = cos(qJ(4));
t1860 = sin(qJ(3));
t1979 = qJD(1) * t1860;
t1808 = -t1862 * qJD(3) + t1859 * t1979;
t1809 = qJD(3) * t1859 + t1862 * t1979;
t1853 = sin(pkin(10));
t1855 = cos(pkin(10));
t1754 = -t1853 * t1808 + t1855 * t1809;
t1751 = t1754 ^ 2;
t1863 = cos(qJ(3));
t1978 = qJD(1) * t1863;
t1839 = -qJD(4) + t1978;
t1993 = t1839 ^ 2;
t1657 = t1993 + t1751;
t1844 = qJD(3) * t1979;
t1933 = t1863 * qJDD(1);
t1815 = -t1844 + t1933;
t1803 = -qJDD(4) + t1815;
t1752 = t1855 * t1808 + t1809 * t1853;
t1952 = t1754 * t1752;
t2007 = -t1952 + t1803;
t1963 = t2007 * t1853;
t1574 = -t1657 * t1855 + t1963;
t1962 = t2007 * t1855;
t1577 = t1657 * t1853 + t1962;
t1464 = t1574 * t1859 - t1577 * t1862;
t1927 = qJD(3) * t1978;
t1934 = t1860 * qJDD(1);
t1889 = t1927 + t1934;
t1870 = -t1859 * qJDD(3) - t1862 * t1889;
t1747 = -t1808 * qJD(4) - t1870;
t1869 = t1862 * qJDD(3) - t1859 * t1889;
t1867 = t1809 * qJD(4) - t1869;
t1653 = t1855 * t1747 - t1853 * t1867;
t1953 = t1752 * t1839;
t2012 = t1653 + t1953;
t1429 = t1464 * t1863 - t1860 * t2012;
t1484 = t1574 * t1862 + t1577 * t1859;
t1854 = sin(pkin(9));
t1856 = cos(pkin(9));
t1360 = t1429 * t1854 + t1484 * t1856;
t1362 = t1429 * t1856 - t1484 * t1854;
t1861 = sin(qJ(1));
t1864 = cos(qJ(1));
t2155 = pkin(6) * (t1360 * t1861 - t1362 * t1864);
t2154 = pkin(6) * (t1360 * t1864 + t1362 * t1861);
t2153 = pkin(1) * t1360;
t2152 = qJ(2) * t1360;
t1427 = t1464 * t1860 + t1863 * t2012;
t2151 = -pkin(1) * t1427 + qJ(2) * t1362;
t1652 = t1747 * t1853 + t1855 * t1867;
t1724 = t1754 * t1839;
t2013 = t1652 - t1724;
t1512 = -t2013 * t1853 + t1855 * t2012;
t1965 = t2012 * t1853;
t1514 = t2013 * t1855 + t1965;
t1420 = t1512 * t1859 + t1514 * t1862;
t1995 = t1752 ^ 2;
t2010 = t1751 - t1995;
t1405 = t1420 * t1863 - t1860 * t2010;
t1416 = -t1512 * t1862 + t1514 * t1859;
t1341 = t1405 * t1854 - t1416 * t1856;
t1343 = t1405 * t1856 + t1416 * t1854;
t2148 = t1341 * t1864 + t1343 * t1861;
t2147 = t1341 * t1861 - t1343 * t1864;
t1713 = t1995 - t1993;
t1582 = t1713 * t1853 - t1962;
t1588 = t1713 * t1855 + t1963;
t1497 = t1582 * t1859 - t1588 * t1862;
t1608 = t1652 + t1724;
t1456 = t1497 * t1863 + t1608 * t1860;
t1492 = t1582 * t1862 + t1588 * t1859;
t1381 = t1456 * t1854 + t1492 * t1856;
t1385 = t1456 * t1856 - t1492 * t1854;
t2146 = t1381 * t1864 + t1385 * t1861;
t2145 = t1381 * t1861 - t1385 * t1864;
t2143 = pkin(2) * t1427;
t2142 = pkin(7) * t1427;
t2141 = pkin(2) * t1484 + pkin(7) * t1429;
t2139 = pkin(3) * t1484;
t2138 = pkin(8) * t1484;
t2135 = pkin(3) * t2012 + pkin(8) * t1464;
t2006 = -t1953 + t1653;
t2053 = -t1608 * t1855 + t2006 * t1853;
t2054 = -t1608 * t1853 - t1855 * t2006;
t2076 = t1859 * t2053 + t1862 * t2054;
t1624 = -t1995 - t1751;
t2077 = -t1859 * t2054 + t1862 * t2053;
t2090 = t1624 * t1860 + t1863 * t2077;
t2113 = t1854 * t2090 - t1856 * t2076;
t2134 = pkin(1) * t2113;
t2133 = qJ(2) * t2113;
t2092 = -t1624 * t1863 + t1860 * t2077;
t2111 = t1854 * t2076 + t1856 * t2090;
t2128 = -pkin(1) * t2092 + qJ(2) * t2111;
t1403 = t1420 * t1860 + t1863 * t2010;
t1452 = t1497 * t1860 - t1608 * t1863;
t1714 = t1751 - t1993;
t2008 = -t1952 - t1803;
t1961 = t2008 * t1853;
t2057 = -t1855 * t1714 + t1961;
t1634 = t1855 * t2008;
t2058 = t1714 * t1853 + t1634;
t2075 = t1859 * t2058 + t1862 * t2057;
t2074 = -t1859 * t2057 + t1862 * t2058;
t2091 = t1860 * t2006 + t1863 * t2074;
t2110 = t1854 * t2075 + t1856 * t2091;
t2112 = t1854 * t2091 - t1856 * t2075;
t2127 = -t1861 * t2112 + t1864 * t2110;
t2126 = t1861 * t2110 + t1864 * t2112;
t2125 = pkin(6) * (-t1861 * t2113 + t1864 * t2111);
t2124 = pkin(6) * (t1861 * t2111 + t1864 * t2113);
t2005 = -t1993 - t1995;
t2028 = t1855 * t2005 - t1961;
t2030 = t1853 * t2005 + t1634;
t2045 = t1859 * t2028 + t1862 * t2030;
t2046 = -t1859 * t2030 + t1862 * t2028;
t2078 = t1860 * t2013 + t1863 * t2046;
t2095 = t1854 * t2078 - t1856 * t2045;
t2122 = pkin(1) * t2095;
t2121 = pkin(2) * t2092;
t2120 = pkin(4) * t1574;
t2119 = pkin(7) * t2092;
t2118 = qJ(2) * t2095;
t2117 = qJ(5) * t1574;
t2116 = qJ(5) * t1577;
t2115 = -pkin(2) * t2076 + pkin(7) * t2090;
t2079 = t1860 * t2046 - t1863 * t2013;
t2094 = t1854 * t2045 + t1856 * t2078;
t2114 = -pkin(1) * t2079 + qJ(2) * t2094;
t2109 = pkin(6) * (-t1861 * t2095 + t1864 * t2094);
t2108 = pkin(6) * (t1861 * t2094 + t1864 * t2095);
t2105 = pkin(2) * t2079;
t2104 = pkin(3) * t2076;
t2103 = pkin(7) * t2079;
t2102 = pkin(8) * t2076;
t2097 = -pkin(2) * t2045 + pkin(7) * t2078;
t2096 = -pkin(3) * t1624 + pkin(8) * t2077;
t2093 = t1860 * t2074 - t1863 * t2006;
t1989 = pkin(4) * t2054;
t2088 = pkin(8) * t2045;
t2087 = qJ(5) * t2054;
t2067 = pkin(4) * t2030;
t2082 = -pkin(3) * t2045 - t2067;
t2081 = -pkin(3) * t2013 + pkin(8) * t2046;
t2080 = -pkin(4) * t1624 + qJ(5) * t2053;
t1943 = t1839 * t1855;
t1930 = t1752 * t1943;
t1890 = t1652 * t1853 - t1930;
t1944 = t1839 * t1853;
t1896 = -t1855 * t1652 - t1752 * t1944;
t1997 = t1859 * t1890 + t1862 * t1896;
t1932 = t1860 * t1952;
t1998 = -t1859 * t1896 + t1862 * t1890;
t2025 = t1863 * t1998 - t1932;
t2047 = t1854 * t1997 + t1856 * t2025;
t2049 = t1854 * t2025 - t1856 * t1997;
t2073 = -t1861 * t2049 + t1864 * t2047;
t2072 = t1861 * t2047 + t1864 * t2049;
t1887 = (t1752 * t1853 + t1754 * t1855) * t1839;
t1712 = t1754 * t1944;
t1895 = -t1712 + t1930;
t2000 = t1859 * t1895 + t1862 * t1887;
t1951 = t1803 * t1860;
t1999 = -t1859 * t1887 + t1862 * t1895;
t2024 = t1863 * t1999 - t1951;
t2048 = t1854 * t2000 + t1856 * t2024;
t2050 = t1854 * t2024 - t1856 * t2000;
t2071 = -t1861 * t2050 + t1864 * t2048;
t2070 = t1861 * t2048 + t1864 * t2050;
t2066 = qJ(5) * t2028;
t2065 = qJ(5) * t2030;
t2064 = qJ(6) * t2012;
t1830 = g(1) * t1864 + g(2) * t1861;
t1865 = qJD(1) ^ 2;
t1811 = -pkin(1) * t1865 - t1830;
t1829 = g(1) * t1861 - t1864 * g(2);
t1891 = qJDD(1) * pkin(1) + t1829;
t1756 = t1856 * t1811 + t1854 * t1891;
t1732 = -t1865 * pkin(2) + qJDD(1) * pkin(7) + t1756;
t1991 = pkin(3) * t1863;
t1901 = -pkin(8) * t1860 - t1991;
t1813 = t1901 * qJD(1);
t1984 = g(3) - qJDD(2);
t1843 = t1863 * t1984;
t1996 = qJD(3) ^ 2;
t1668 = (qJD(1) * t1813 + t1732) * t1860 - qJDD(3) * pkin(3) - pkin(8) * t1996 + t1843;
t1777 = -pkin(4) * t1839 - qJ(5) * t1809;
t1994 = t1808 ^ 2;
t1555 = t1867 * pkin(4) - t1994 * qJ(5) + t1809 * t1777 + qJDD(5) + t1668;
t2059 = t1652 * pkin(5) + t1555 - t2064;
t1755 = t1854 * t1811 - t1856 * t1891;
t1920 = t1755 * t1854 + t1856 * t1756;
t1666 = t1755 * t1856 - t1756 * t1854;
t1960 = t1666 * t1861;
t2056 = t1864 * t1920 + t1960;
t1959 = t1666 * t1864;
t2055 = -t1861 * t1920 + t1959;
t1818 = qJDD(1) * t1854 + t1856 * t1865;
t1819 = qJDD(1) * t1856 - t1854 * t1865;
t1759 = -t1818 * t1861 + t1864 * t1819;
t1786 = -qJ(2) * t1818 + t1856 * t1984;
t2002 = -qJ(2) * t1819 - t1854 * t1984;
t2052 = -pkin(6) * t1759 - t1786 * t1861 + t1864 * t2002;
t2009 = t1864 * t1818 + t1819 * t1861;
t2051 = -pkin(6) * t2009 + t1786 * t1864 + t1861 * t2002;
t1791 = t1863 * t1803;
t2029 = t1860 * t1999 + t1791;
t1931 = t1863 * t1952;
t2027 = t1860 * t1998 + t1931;
t1598 = t1653 * t1853 - t1754 * t1943;
t1599 = t1653 * t1855 + t1712;
t1505 = t1598 * t1862 + t1599 * t1859;
t1508 = -t1598 * t1859 + t1599 * t1862;
t1897 = t1863 * t1508 + t1932;
t2001 = t1505 * t1854 + t1856 * t1897;
t2003 = -t1856 * t1505 + t1854 * t1897;
t2022 = -t1861 * t2003 + t1864 * t2001;
t2021 = t1861 * t2001 + t1864 * t2003;
t1950 = t1809 * t1808;
t1883 = -t1803 - t1950;
t2019 = t1859 * t1883;
t2017 = t1862 * t1883;
t2015 = pkin(5) * t1657 - qJ(6) * t2007;
t2014 = -pkin(5) * t2008 - qJ(6) * t2005;
t1700 = t1860 * t1732 + t1843;
t1702 = t1863 * t1732 - t1860 * t1984;
t1623 = t1860 * t1700 + t1863 * t1702;
t1977 = qJD(5) * t1752;
t1743 = -0.2e1 * t1977;
t1975 = qJD(6) * t1839;
t2011 = t1743 - 0.2e1 * t1975;
t1788 = t1808 * t1839;
t1692 = -t1788 + t1747;
t1802 = t1809 ^ 2;
t1992 = pkin(3) * t1860;
t1731 = -qJDD(1) * pkin(2) - t1865 * pkin(7) + t1755;
t1814 = 0.2e1 * t1927 + t1934;
t1899 = -t1815 + t1844;
t1651 = pkin(3) * t1899 - pkin(8) * t1814 + t1731;
t1669 = -pkin(3) * t1996 + qJDD(3) * pkin(8) + t1813 * t1978 + t1702;
t1560 = -t1862 * t1651 + t1859 * t1669;
t1509 = t1883 * pkin(4) - qJ(5) * t1692 - t1560;
t1561 = t1859 * t1651 + t1862 * t1669;
t1520 = -pkin(4) * t1994 - qJ(5) * t1867 + t1839 * t1777 + t1561;
t1921 = -t1855 * t1509 + t1853 * t1520;
t1976 = qJD(5) * t1754;
t1411 = t1921 + 0.2e1 * t1976;
t1938 = t1853 * t1509 + t1855 * t1520;
t1412 = t1743 + t1938;
t1339 = -t1411 * t1855 + t1412 * t1853;
t1990 = pkin(4) * t1339;
t1986 = pkin(5) * t1855;
t1981 = qJ(6) * t1855;
t1980 = qJD(1) * qJD(3);
t1974 = t1339 * t1859;
t1973 = t1339 * t1862;
t1971 = t1555 * t1853;
t1970 = t1555 * t1855;
t1958 = t1668 * t1859;
t1957 = t1668 * t1862;
t1727 = t1803 - t1950;
t1955 = t1727 * t1859;
t1954 = t1727 * t1862;
t1768 = t1814 * t1860;
t1838 = t1863 * t1865 * t1860;
t1825 = -t1838 + qJDD(3);
t1947 = t1825 * t1860;
t1946 = t1825 * t1863;
t1826 = qJDD(3) + t1838;
t1945 = t1826 * t1860;
t1942 = t1839 * t1859;
t1941 = t1839 * t1862;
t1849 = t1860 ^ 2;
t1940 = t1849 * t1865;
t1717 = t1860 * t1731;
t1718 = t1863 * t1731;
t1670 = pkin(5) * t1752 - qJ(6) * t1754;
t1893 = -pkin(5) * t1993 - t1803 * qJ(6) - t1752 * t1670 + t1938;
t1377 = t1893 + t2011;
t1886 = t1803 * pkin(5) - qJ(6) * t1993 + qJDD(6) + t1921;
t1387 = (0.2e1 * qJD(5) + t1670) * t1754 + t1886;
t1939 = -pkin(5) * t1387 + qJ(6) * t1377;
t1937 = -pkin(5) * t2006 - qJ(6) * t1608;
t1936 = -pkin(2) * t1731 + pkin(7) * t1623;
t1850 = t1863 ^ 2;
t1935 = t1849 + t1850;
t1929 = t1860 * t1950;
t1928 = t1863 * t1950;
t1834 = -t1940 - t1996;
t1776 = -t1834 * t1860 - t1946;
t1926 = -pkin(2) * t1814 + pkin(7) * t1776 + t1717;
t1847 = t1850 * t1865;
t1836 = -t1847 - t1996;
t1774 = t1836 * t1863 - t1945;
t1816 = -0.2e1 * t1844 + t1933;
t1925 = pkin(2) * t1816 + pkin(7) * t1774 - t1718;
t1923 = -qJ(6) * t1853 - pkin(4);
t1340 = t1411 * t1853 + t1855 * t1412;
t1477 = t1560 * t1859 + t1862 * t1561;
t1919 = -t1829 * t1861 - t1864 * t1830;
t1918 = t1854 * t1838;
t1917 = t1856 * t1838;
t1320 = t1377 * t1853 - t1387 * t1855;
t1321 = t1377 * t1855 + t1387 * t1853;
t1294 = t1320 * t1862 + t1321 * t1859;
t1434 = (-pkin(5) * t1839 - 0.2e1 * qJD(6)) * t1754 + t2059;
t1296 = qJ(5) * t1321 + (t1923 - t1986) * t1434;
t1304 = -qJ(5) * t1320 + (pkin(5) * t1853 - t1981) * t1434;
t1253 = -pkin(8) * t1294 - t1296 * t1859 + t1304 * t1862;
t1905 = pkin(4) * t1320 + t1939;
t1265 = -pkin(3) * t1294 - t1905;
t1295 = -t1320 * t1859 + t1321 * t1862;
t1290 = t1295 * t1863 + t1434 * t1860;
t1916 = -pkin(2) * t1294 + pkin(7) * t1290 + t1860 * t1253 + t1863 * t1265;
t1306 = t1340 * t1859 + t1973;
t1326 = -pkin(4) * t1555 + qJ(5) * t1340;
t1276 = -pkin(8) * t1306 - qJ(5) * t1973 - t1326 * t1859;
t1287 = -pkin(3) * t1306 - t1990;
t1307 = t1340 * t1862 - t1974;
t1303 = t1307 * t1863 + t1555 * t1860;
t1915 = -pkin(2) * t1306 + pkin(7) * t1303 + t1860 * t1276 + t1863 * t1287;
t1365 = -pkin(5) * t1624 + t1377;
t1366 = -qJ(6) * t1624 + t1387;
t1313 = t1365 * t1855 + t1366 * t1853 + t2080;
t1315 = -t1365 * t1853 + t1366 * t1855 - t2087;
t1278 = -t1313 * t1859 + t1315 * t1862 - t2102;
t1903 = t1937 + t1989;
t1351 = -t1903 - t2104;
t1914 = t1860 * t1278 + t1863 * t1351 + t2115;
t1319 = t1340 + t2080;
t1324 = -t1339 - t2087;
t1285 = -t1319 * t1859 + t1324 * t1862 - t2102;
t1368 = -t1989 - t2104;
t1913 = t1860 * t1285 + t1863 * t1368 + t2115;
t1866 = 0.2e1 * qJD(6) * t1754 - t2059;
t1407 = pkin(5) * t1724 + t1866 + t2064;
t1352 = -t2116 + t1853 * t1407 + (pkin(4) + t1986) * t2012;
t1357 = -pkin(5) * t1965 + t1407 * t1855 + t2117;
t1309 = -t1352 * t1859 + t1357 * t1862 + t2138;
t1744 = 0.2e1 * t1977;
t1879 = t1893 + t2015 - t2120;
t1323 = t1744 - t1879 + 0.2e1 * t1975 + t2139;
t1912 = t1860 * t1309 + t1863 * t1323 + t2141;
t1408 = (-t2013 + t1724) * pkin(5) + t1866;
t1355 = t1855 * t1408 + t1923 * t2013 + t2066;
t1359 = -t1408 * t1853 - t1981 * t2013 - t2065;
t1312 = -t1355 * t1859 + t1359 * t1862 - t2088;
t1329 = t1387 + t2014 + t2082;
t1911 = t1860 * t1312 + t1863 * t1329 + t2097;
t1431 = -pkin(4) * t2013 - t1970 + t2066;
t1474 = t1971 - t2065;
t1333 = -t1431 * t1859 + t1474 * t1862 - t2088;
t1347 = t1411 + t2082;
t1910 = t1860 * t1333 + t1863 * t1347 + t2097;
t1435 = -pkin(4) * t2012 + t1971 + t2116;
t1478 = t1970 - t2117;
t1349 = -t1435 * t1859 + t1478 * t1862 - t2138;
t1904 = -t1938 + t2120;
t1354 = t1743 - t1904 - t2139;
t1909 = t1860 * t1349 + t1863 * t1354 - t2141;
t1748 = -t1993 - t1994;
t1641 = t1748 * t1859 + t2017;
t1521 = -pkin(3) * t1641 + t1560;
t1559 = -pkin(8) * t1641 + t1958;
t1642 = t1748 * t1862 - t2019;
t1789 = t1839 * t1809;
t1689 = t1789 - t1867;
t1570 = t1642 * t1863 - t1689 * t1860;
t1908 = -pkin(2) * t1641 + pkin(7) * t1570 + t1863 * t1521 + t1860 * t1559;
t1763 = -t1802 - t1993;
t1658 = t1763 * t1862 + t1955;
t1523 = -pkin(3) * t1658 + t1561;
t1565 = -pkin(8) * t1658 + t1957;
t1659 = -t1763 * t1859 + t1954;
t1693 = (qJD(4) - t1839) * t1808 + t1870;
t1579 = t1659 * t1863 - t1693 * t1860;
t1907 = -pkin(2) * t1658 + pkin(7) * t1579 + t1863 * t1523 + t1860 * t1565;
t1820 = t1935 * qJDD(1);
t1823 = t1847 + t1940;
t1906 = pkin(2) * t1823 + pkin(7) * t1820 + t1623;
t1902 = -pkin(3) * t1668 + pkin(8) * t1477;
t1822 = qJDD(1) * t1864 - t1861 * t1865;
t1900 = -pkin(6) * t1822 - g(3) * t1861;
t1898 = t1860 * t1508 - t1931;
t1476 = -t1560 * t1862 + t1561 * t1859;
t1622 = t1700 * t1863 - t1702 * t1860;
t1892 = t1829 * t1864 - t1830 * t1861;
t1690 = (-qJD(4) - t1839) * t1809 + t1869;
t1616 = t1690 * t1859 - t1692 * t1862;
t1433 = -pkin(8) * t1616 - t1476;
t1618 = t1690 * t1862 + t1692 * t1859;
t1725 = t1802 + t1994;
t1546 = t1618 * t1863 - t1725 * t1860;
t1888 = pkin(7) * t1546 + t1860 * t1433 + (-pkin(2) - t1991) * t1616;
t1885 = pkin(3) * t1689 + pkin(8) * t1642 - t1957;
t1884 = pkin(3) * t1693 + pkin(8) * t1659 + t1958;
t1882 = pkin(3) * t1725 + pkin(8) * t1618 + t1477;
t1444 = t1477 * t1863 + t1668 * t1860;
t1881 = pkin(7) * t1444 + (-pkin(2) + t1901) * t1476;
t1878 = -pkin(3) * t1434 + pkin(8) * t1295 + t1296 * t1862 + t1304 * t1859;
t1877 = t1313 * t1862 + t1315 * t1859 + t2096;
t1876 = t1319 * t1862 + t1324 * t1859 + t2096;
t1875 = t1352 * t1862 + t1357 * t1859 + t2135;
t1874 = t1355 * t1862 + t1359 * t1859 + t2081;
t1873 = t1431 * t1862 + t1474 * t1859 + t2081;
t1872 = t1435 * t1862 + t1478 * t1859 - t2135;
t1871 = -pkin(3) * t1555 + pkin(8) * t1307 - qJ(5) * t1974 + t1326 * t1862;
t1745 = -0.2e1 * t1976;
t1868 = -t1670 * t1754 + t1745 - t1886 - t2014;
t1835 = t1847 - t1996;
t1833 = -t1940 + t1996;
t1824 = -t1847 + t1940;
t1821 = qJDD(1) * t1861 + t1864 * t1865;
t1807 = t1863 * t1826;
t1806 = t1935 * t1980;
t1790 = -pkin(6) * t1821 + g(3) * t1864;
t1783 = -t1802 + t1993;
t1782 = -t1993 + t1994;
t1781 = -t1849 * t1980 + t1863 * t1889;
t1780 = -t1815 * t1860 - t1850 * t1980;
t1779 = qJDD(3) * t1854 + t1806 * t1856;
t1778 = -qJDD(3) * t1856 + t1806 * t1854;
t1775 = -t1833 * t1860 + t1807;
t1773 = t1835 * t1863 - t1947;
t1772 = t1834 * t1863 - t1947;
t1771 = t1833 * t1863 + t1945;
t1770 = t1836 * t1860 + t1807;
t1769 = t1835 * t1860 + t1946;
t1767 = t1899 * t1863;
t1764 = t1802 - t1994;
t1762 = t1820 * t1856 - t1823 * t1854;
t1761 = t1820 * t1854 + t1823 * t1856;
t1758 = t1816 * t1863 - t1768;
t1757 = t1814 * t1863 + t1816 * t1860;
t1741 = t1781 * t1856 - t1918;
t1740 = t1780 * t1856 + t1918;
t1739 = t1781 * t1854 + t1917;
t1738 = t1780 * t1854 - t1917;
t1737 = t1775 * t1856 + t1854 * t1934;
t1736 = t1773 * t1856 + t1854 * t1933;
t1735 = t1775 * t1854 - t1856 * t1934;
t1734 = t1773 * t1854 - t1856 * t1933;
t1711 = t1776 * t1856 + t1814 * t1854;
t1710 = t1774 * t1856 - t1816 * t1854;
t1709 = t1776 * t1854 - t1814 * t1856;
t1708 = t1774 * t1854 + t1816 * t1856;
t1707 = -pkin(1) * t1818 - t1756;
t1706 = pkin(1) * t1819 - t1755;
t1704 = (t1808 * t1862 - t1809 * t1859) * t1839;
t1703 = (t1808 * t1859 + t1809 * t1862) * t1839;
t1701 = t1758 * t1856 + t1824 * t1854;
t1699 = t1758 * t1854 - t1824 * t1856;
t1691 = t1788 + t1747;
t1688 = t1789 + t1867;
t1687 = t1747 * t1862 + t1809 * t1942;
t1686 = t1747 * t1859 - t1809 * t1941;
t1685 = -t1808 * t1941 + t1859 * t1867;
t1684 = t1808 * t1942 + t1862 * t1867;
t1683 = t1704 * t1863 - t1951;
t1682 = t1704 * t1860 + t1791;
t1678 = -pkin(7) * t1772 + t1718;
t1677 = t1782 * t1862 + t1955;
t1676 = -t1783 * t1859 + t2017;
t1675 = -pkin(7) * t1770 + t1717;
t1674 = t1782 * t1859 - t1954;
t1673 = t1783 * t1862 + t2019;
t1663 = pkin(1) * t1666;
t1661 = -pkin(2) * t1772 + t1702;
t1660 = -pkin(2) * t1770 + t1700;
t1650 = pkin(1) * t1984 + qJ(2) * t1920;
t1629 = t1687 * t1863 + t1929;
t1628 = t1685 * t1863 - t1929;
t1627 = t1687 * t1860 - t1928;
t1626 = t1685 * t1860 + t1928;
t1617 = t1689 * t1862 - t1691 * t1859;
t1615 = t1689 * t1859 + t1691 * t1862;
t1603 = t1683 * t1856 + t1703 * t1854;
t1602 = t1683 * t1854 - t1703 * t1856;
t1593 = pkin(1) * t1708 + t1925;
t1592 = pkin(1) * t1709 + t1926;
t1591 = t1677 * t1863 - t1688 * t1860;
t1590 = t1676 * t1863 + t1692 * t1860;
t1585 = t1677 * t1860 + t1688 * t1863;
t1584 = t1676 * t1860 - t1692 * t1863;
t1576 = t1659 * t1860 + t1693 * t1863;
t1572 = -qJ(2) * t1761 + t1622 * t1856;
t1571 = qJ(2) * t1762 + t1622 * t1854;
t1569 = t1642 * t1860 + t1689 * t1863;
t1567 = t1623 * t1856 + t1731 * t1854;
t1566 = t1623 * t1854 - t1731 * t1856;
t1564 = t1617 * t1863 + t1764 * t1860;
t1563 = t1617 * t1860 - t1764 * t1863;
t1558 = pkin(1) * t1761 + t1906;
t1550 = t1629 * t1856 + t1686 * t1854;
t1549 = t1628 * t1856 - t1684 * t1854;
t1548 = t1629 * t1854 - t1686 * t1856;
t1547 = t1628 * t1854 + t1684 * t1856;
t1545 = t1618 * t1860 + t1725 * t1863;
t1539 = -qJ(2) * t1709 - t1661 * t1854 + t1678 * t1856;
t1538 = -qJ(2) * t1708 - t1660 * t1854 + t1675 * t1856;
t1533 = t1591 * t1856 + t1674 * t1854;
t1532 = t1590 * t1856 + t1673 * t1854;
t1531 = t1591 * t1854 - t1674 * t1856;
t1530 = t1590 * t1854 - t1673 * t1856;
t1529 = t1579 * t1856 + t1658 * t1854;
t1528 = t1579 * t1854 - t1658 * t1856;
t1527 = -pkin(1) * t1772 + qJ(2) * t1711 + t1661 * t1856 + t1678 * t1854;
t1526 = -pkin(1) * t1770 + qJ(2) * t1710 + t1660 * t1856 + t1675 * t1854;
t1525 = t1570 * t1856 + t1641 * t1854;
t1524 = t1570 * t1854 - t1641 * t1856;
t1489 = t1564 * t1856 + t1615 * t1854;
t1488 = t1564 * t1854 - t1615 * t1856;
t1481 = t1546 * t1856 + t1616 * t1854;
t1480 = t1546 * t1854 - t1616 * t1856;
t1479 = pkin(1) * t1566 + t1936;
t1459 = -pkin(2) * t1576 - t1884;
t1458 = -pkin(2) * t1569 - t1885;
t1457 = -qJ(2) * t1566 - (pkin(2) * t1854 - pkin(7) * t1856) * t1622;
t1443 = t1477 * t1860 - t1668 * t1863;
t1424 = qJ(2) * t1567 - (-pkin(2) * t1856 - pkin(7) * t1854 - pkin(1)) * t1622;
t1413 = -pkin(7) * t1576 - t1523 * t1860 + t1565 * t1863;
t1409 = -pkin(7) * t1569 - t1521 * t1860 + t1559 * t1863;
t1388 = -pkin(2) * t1545 - t1882;
t1371 = t1444 * t1856 + t1476 * t1854;
t1370 = t1444 * t1854 - t1476 * t1856;
t1369 = -pkin(7) * t1545 + t1433 * t1863 + t1616 * t1992;
t1364 = -pkin(2) * t1443 - t1902;
t1358 = pkin(1) * t1528 + t1907;
t1356 = pkin(1) * t1524 + t1908;
t1345 = -pkin(7) * t1443 + (-pkin(8) * t1863 + t1992) * t1476;
t1338 = -qJ(2) * t1528 + t1413 * t1856 - t1459 * t1854;
t1331 = pkin(1) * t1480 + t1888;
t1330 = -qJ(2) * t1524 + t1409 * t1856 - t1458 * t1854;
t1327 = -pkin(1) * t1576 + qJ(2) * t1529 + t1413 * t1854 + t1459 * t1856;
t1325 = -pkin(1) * t1569 + qJ(2) * t1525 + t1409 * t1854 + t1458 * t1856;
t1318 = -t1872 + t2143;
t1317 = -qJ(2) * t1480 + t1369 * t1856 - t1388 * t1854;
t1316 = -t1873 - t2105;
t1314 = -pkin(1) * t1545 + qJ(2) * t1481 + t1369 * t1854 + t1388 * t1856;
t1310 = pkin(1) * t1370 + t1881;
t1302 = t1307 * t1860 - t1555 * t1863;
t1300 = t1349 * t1863 - t1354 * t1860 + t2142;
t1299 = -qJ(2) * t1370 + t1345 * t1856 - t1364 * t1854;
t1298 = -t1874 - t2105;
t1297 = t1333 * t1863 - t1347 * t1860 - t2103;
t1292 = -t1875 - t2143;
t1291 = -pkin(1) * t1443 + qJ(2) * t1371 + t1345 * t1854 + t1364 * t1856;
t1289 = t1295 * t1860 - t1434 * t1863;
t1283 = t1909 - t2153;
t1282 = t1312 * t1863 - t1329 * t1860 - t2103;
t1281 = t1910 + t2122;
t1280 = t1309 * t1863 - t1323 * t1860 - t2142;
t1279 = -t1876 - t2121;
t1274 = t1303 * t1856 + t1306 * t1854;
t1273 = t1303 * t1854 - t1306 * t1856;
t1272 = t1285 * t1863 - t1368 * t1860 - t2119;
t1271 = t1300 * t1856 - t1318 * t1854 + t2152;
t1270 = t1911 + t2122;
t1269 = -t1877 - t2121;
t1268 = t1300 * t1854 + t1318 * t1856 - t2151;
t1267 = t1912 + t2153;
t1266 = t1297 * t1856 - t1316 * t1854 - t2118;
t1263 = t1297 * t1854 + t1316 * t1856 + t2114;
t1262 = t1278 * t1863 - t1351 * t1860 - t2119;
t1261 = t1290 * t1856 + t1294 * t1854;
t1260 = t1290 * t1854 - t1294 * t1856;
t1259 = t1913 + t2134;
t1258 = t1282 * t1856 - t1298 * t1854 - t2118;
t1257 = -pkin(2) * t1302 - t1871;
t1256 = t1282 * t1854 + t1298 * t1856 + t2114;
t1255 = t1914 + t2134;
t1254 = t1280 * t1856 - t1292 * t1854 - t2152;
t1251 = t1280 * t1854 + t1292 * t1856 + t2151;
t1250 = t1272 * t1856 - t1279 * t1854 - t2133;
t1249 = -pkin(7) * t1302 + t1276 * t1863 - t1287 * t1860;
t1248 = t1272 * t1854 + t1279 * t1856 + t2128;
t1247 = -pkin(2) * t1289 - t1878;
t1246 = t1262 * t1856 - t1269 * t1854 - t2133;
t1245 = t1262 * t1854 + t1269 * t1856 + t2128;
t1244 = -pkin(7) * t1289 + t1253 * t1863 - t1265 * t1860;
t1243 = pkin(1) * t1273 + t1915;
t1242 = -qJ(2) * t1273 + t1249 * t1856 - t1257 * t1854;
t1241 = -pkin(1) * t1302 + qJ(2) * t1274 + t1249 * t1854 + t1257 * t1856;
t1240 = pkin(1) * t1260 + t1916;
t1239 = -qJ(2) * t1260 + t1244 * t1856 - t1247 * t1854;
t1238 = -pkin(1) * t1289 + qJ(2) * t1261 + t1244 * t1854 + t1247 * t1856;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t1822, 0, -t1821, 0, t1900, -t1790, -t1892, -pkin(6) * t1892, 0, 0, t1759, 0, -t2009, 0, t2052, -t2051, t2055, pkin(6) * t2055 + qJ(2) * t1959 - t1861 * t1650, -t1739 * t1861 + t1741 * t1864, -t1699 * t1861 + t1701 * t1864, -t1735 * t1861 + t1737 * t1864, -t1738 * t1861 + t1740 * t1864, -t1734 * t1861 + t1736 * t1864, -t1778 * t1861 + t1779 * t1864, t1864 * t1538 - t1861 * t1526 - pkin(6) * (t1708 * t1864 + t1710 * t1861), t1864 * t1539 - t1861 * t1527 - pkin(6) * (t1709 * t1864 + t1711 * t1861), t1864 * t1572 - t1861 * t1571 - pkin(6) * (t1761 * t1864 + t1762 * t1861), t1864 * t1457 - t1861 * t1424 - pkin(6) * (t1566 * t1864 + t1567 * t1861), -t1548 * t1861 + t1550 * t1864, -t1488 * t1861 + t1489 * t1864, -t1530 * t1861 + t1532 * t1864, -t1547 * t1861 + t1549 * t1864, -t1531 * t1861 + t1533 * t1864, -t1602 * t1861 + t1603 * t1864, t1864 * t1330 - t1861 * t1325 - pkin(6) * (t1524 * t1864 + t1525 * t1861), t1864 * t1338 - t1861 * t1327 - pkin(6) * (t1528 * t1864 + t1529 * t1861), t1864 * t1317 - t1861 * t1314 - pkin(6) * (t1480 * t1864 + t1481 * t1861), t1864 * t1299 - t1861 * t1291 - pkin(6) * (t1370 * t1864 + t1371 * t1861), t2022, t2147, t2127, t2073, t2145, t2071, -t1861 * t1263 + t1864 * t1266 - t2108, -t1861 * t1268 + t1864 * t1271 + t2154, -t1861 * t1248 + t1864 * t1250 - t2124, t1864 * t1242 - t1861 * t1241 - pkin(6) * (t1273 * t1864 + t1274 * t1861), t2022, t2127, -t2147, t2071, -t2145, t2073, -t1861 * t1256 + t1864 * t1258 - t2108, -t1861 * t1245 + t1864 * t1246 - t2124, -t1861 * t1251 + t1864 * t1254 - t2154, t1864 * t1239 - t1861 * t1238 - pkin(6) * (t1260 * t1864 + t1261 * t1861); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t1821, 0, t1822, 0, t1790, t1900, t1919, pkin(6) * t1919, 0, 0, t2009, 0, t1759, 0, t2051, t2052, t2056, pkin(6) * t2056 + qJ(2) * t1960 + t1864 * t1650, t1739 * t1864 + t1741 * t1861, t1699 * t1864 + t1701 * t1861, t1735 * t1864 + t1737 * t1861, t1738 * t1864 + t1740 * t1861, t1734 * t1864 + t1736 * t1861, t1778 * t1864 + t1779 * t1861, t1861 * t1538 + t1864 * t1526 + pkin(6) * (-t1708 * t1861 + t1710 * t1864), t1861 * t1539 + t1864 * t1527 + pkin(6) * (-t1709 * t1861 + t1711 * t1864), t1861 * t1572 + t1864 * t1571 + pkin(6) * (-t1761 * t1861 + t1762 * t1864), t1861 * t1457 + t1864 * t1424 + pkin(6) * (-t1566 * t1861 + t1567 * t1864), t1548 * t1864 + t1550 * t1861, t1488 * t1864 + t1489 * t1861, t1530 * t1864 + t1532 * t1861, t1547 * t1864 + t1549 * t1861, t1531 * t1864 + t1533 * t1861, t1602 * t1864 + t1603 * t1861, t1861 * t1330 + t1864 * t1325 + pkin(6) * (-t1524 * t1861 + t1525 * t1864), t1861 * t1338 + t1864 * t1327 + pkin(6) * (-t1528 * t1861 + t1529 * t1864), t1861 * t1317 + t1864 * t1314 + pkin(6) * (-t1480 * t1861 + t1481 * t1864), t1861 * t1299 + t1864 * t1291 + pkin(6) * (-t1370 * t1861 + t1371 * t1864), t2021, -t2148, t2126, t2072, -t2146, t2070, t1864 * t1263 + t1861 * t1266 + t2109, t1864 * t1268 + t1861 * t1271 + t2155, t1864 * t1248 + t1861 * t1250 + t2125, t1861 * t1242 + t1864 * t1241 + pkin(6) * (-t1273 * t1861 + t1274 * t1864), t2021, t2126, t2148, t2070, t2146, t2072, t1864 * t1256 + t1861 * t1258 + t2109, t1864 * t1245 + t1861 * t1246 + t2125, t1864 * t1251 + t1861 * t1254 - t2155, t1861 * t1239 + t1864 * t1238 + pkin(6) * (-t1260 * t1861 + t1261 * t1864); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1829, t1830, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1706, t1707, 0, -t1663, t1768, t1757, t1771, -t1767, t1769, 0, t1593, t1592, t1558, t1479, t1627, t1563, t1584, t1626, t1585, t1682, t1356, t1358, t1331, t1310, t1898, -t1403, t2093, t2027, -t1452, t2029, t1281, t1283, t1259, t1243, t1898, t2093, t1403, t2029, t1452, t2027, t1270, t1255, t1267, t1240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1865, 0, 0, -g(3), -t1829, 0, 0, 0, t1819, 0, -t1818, 0, t2002, -t1786, t1666, qJ(2) * t1666, t1741, t1701, t1737, t1740, t1736, t1779, t1538, t1539, t1572, t1457, t1550, t1489, t1532, t1549, t1533, t1603, t1330, t1338, t1317, t1299, t2001, -t1343, t2110, t2047, -t1385, t2048, t1266, t1271, t1250, t1242, t2001, t2110, t1343, t2048, t1385, t2047, t1258, t1246, t1254, t1239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1865, 0, qJDD(1), 0, g(3), 0, -t1830, 0, 0, 0, t1818, 0, t1819, 0, t1786, t2002, t1920, t1650, t1739, t1699, t1735, t1738, t1734, t1778, t1526, t1527, t1571, t1424, t1548, t1488, t1530, t1547, t1531, t1602, t1325, t1327, t1314, t1291, t2003, -t1341, t2112, t2049, -t1381, t2050, t1263, t1268, t1248, t1241, t2003, t2112, t1341, t2050, t1381, t2049, t1256, t1245, t1251, t1238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1829, t1830, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t1706, t1707, 0, -t1663, t1768, t1757, t1771, -t1767, t1769, 0, t1593, t1592, t1558, t1479, t1627, t1563, t1584, t1626, t1585, t1682, t1356, t1358, t1331, t1310, t1898, -t1403, t2093, t2027, -t1452, t2029, t1281, t1283, t1259, t1243, t1898, t2093, t1403, t2029, t1452, t2027, t1270, t1255, t1267, t1240; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t1865, 0, 0, -t1984, t1755, 0, t1781, t1758, t1775, t1780, t1773, t1806, t1675, t1678, t1622, pkin(7) * t1622, t1629, t1564, t1590, t1628, t1591, t1683, t1409, t1413, t1369, t1345, t1897, -t1405, t2091, t2025, -t1456, t2024, t1297, t1300, t1272, t1249, t1897, t2091, t1405, t2024, t1456, t2025, t1282, t1262, t1280, t1244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1865, 0, qJDD(1), 0, t1984, 0, t1756, 0, t1838, -t1824, -t1934, -t1838, -t1933, -qJDD(3), t1660, t1661, 0, pkin(2) * t1622, -t1686, -t1615, -t1673, t1684, -t1674, -t1703, t1458, t1459, t1388, t1364, -t1505, t1416, -t2075, -t1997, -t1492, -t2000, t1316, t1318, t1279, t1257, -t1505, -t2075, -t1416, -t2000, t1492, -t1997, t1298, t1269, t1292, t1247; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t1755, -t1756, 0, 0, t1768, t1757, t1771, -t1767, t1769, 0, t1925, t1926, t1906, t1936, t1627, t1563, t1584, t1626, t1585, t1682, t1908, t1907, t1888, t1881, t1898, -t1403, t2093, t2027, -t1452, t2029, t1910, t1909, t1913, t1915, t1898, t2093, t1403, t2029, t1452, t2027, t1911, t1914, t1912, t1916; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1889, t1816, t1826, -t1927, t1835, t1927, 0, t1731, t1700, 0, t1687, t1617, t1676, t1685, t1677, t1704, t1559, t1565, t1433, -pkin(8) * t1476, t1508, -t1420, t2074, t1998, -t1497, t1999, t1333, t1349, t1285, t1276, t1508, t2074, t1420, t1999, t1497, t1998, t1312, t1278, t1309, t1253; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1844, t1814, t1833, t1815, t1825, -t1844, -t1731, 0, t1702, 0, -t1950, -t1764, -t1692, t1950, t1688, t1803, t1521, t1523, -pkin(3) * t1616, -pkin(3) * t1476, -t1952, -t2010, -t2006, t1952, t1608, t1803, t1347, t1354, t1368, t1287, -t1952, -t2006, t2010, t1803, -t1608, t1952, t1329, t1351, t1323, t1265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1838, t1824, t1934, t1838, t1933, qJDD(3), -t1700, -t1702, 0, 0, t1686, t1615, t1673, -t1684, t1674, t1703, t1885, t1884, t1882, t1902, t1505, -t1416, t2075, t1997, t1492, t2000, t1873, t1872, t1876, t1871, t1505, t2075, t1416, t2000, -t1492, t1997, t1874, t1877, t1875, t1878; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1747, t1689, t1883, -t1788, t1782, t1788, 0, t1668, t1560, 0, t1599, -t1514, t2058, t1890, t1588, t1895, t1474, t1478, t1324, -qJ(5) * t1339, t1599, t2058, t1514, t1895, -t1588, t1890, t1359, t1315, t1357, t1304; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1789, t1691, t1783, -t1867, -t1727, t1789, -t1668, 0, t1561, 0, t1598, t1512, t2057, t1896, t1582, t1887, t1431, t1435, t1319, t1326, t1598, t2057, -t1512, t1887, -t1582, t1896, t1355, t1313, t1352, t1296; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1950, t1764, t1692, -t1950, -t1688, -t1803, -t1560, -t1561, 0, 0, t1952, t2010, t2006, -t1952, -t1608, -t1803, t1745 - t1921 + t2067, t1744 + t1904, t1989, t1990, t1952, t2006, -t2010, -t1803, t1608, -t1952, t1868 + t2067, t1903, t1879 + t2011, t1905; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1653, -t2013, t2008, -t1953, t1713, t1953, 0, t1555, t1411, 0, t1653, t2008, t2013, t1953, -t1713, -t1953, -qJ(6) * t2013, t1366, t1407, -qJ(6) * t1434; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1724, t2012, -t1714, -t1652, -t2007, t1724, -t1555, 0, t1412, 0, -t1724, -t1714, -t2012, t1724, t2007, -t1652, t1408, t1365, pkin(5) * t2012, -pkin(5) * t1434; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1952, t2010, t2006, -t1952, -t1608, -t1803, -t1411, -t1412, 0, 0, t1952, t2006, -t2010, -t1803, t1608, -t1952, t1868, t1937, t1377 + t2015, t1939; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1653, t2008, t2013, t1953, -t1713, -t1953, 0, t1387, -t1434, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1952, t2006, -t2010, -t1803, t1608, -t1952, -t1387, 0, t1377, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1724, t1714, t2012, -t1724, -t2007, t1652, t1434, -t1377, 0, 0;];
m_new_reg  = t1;