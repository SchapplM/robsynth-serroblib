% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6PRRPRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 03:50
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6PRRPRP2_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP2_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_invdynm_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 03:48:53
% EndTime: 2019-05-05 03:50:07
% DurationCPUTime: 79.07s
% Computational Cost: add. (241333->1002), mult. (525274->1397), div. (0->0), fcn. (385106->12), ass. (0->689)
t1947 = cos(qJ(5));
t1944 = sin(qJ(5));
t1936 = sin(pkin(11));
t1939 = cos(pkin(11));
t1948 = cos(qJ(3));
t2115 = qJD(2) * t1948;
t1945 = sin(qJ(3));
t2116 = qJD(2) * t1945;
t1894 = t1936 * t2115 + t1939 * t2116;
t1855 = -qJD(3) * t1947 + t1894 * t1944;
t1924 = qJD(3) * t2115;
t2064 = t1945 * qJDD(2);
t1902 = t1924 + t2064;
t2054 = qJD(3) * t2116;
t2062 = t1948 * qJDD(2);
t1983 = -t2054 + t2062;
t1845 = t1902 * t1939 + t1936 * t1983;
t1993 = -t1944 * qJDD(3) - t1947 * t1845;
t1959 = -qJD(5) * t1855 - t1993;
t1892 = t1936 * t2116 - t1939 * t2115;
t1885 = qJD(5) + t1892;
t2086 = t1855 * t1885;
t2144 = -t2086 + t1959;
t2101 = t2144 * t1944;
t1857 = qJD(3) * t1944 + t1894 * t1947;
t2047 = qJDD(3) * t1947 - t1944 * t1845;
t1779 = qJD(5) * t1857 - t2047;
t1820 = t1857 * t1885;
t2150 = t1779 + t1820;
t1618 = -t1947 * t2150 - t2101;
t1854 = t1857 ^ 2;
t2129 = t1855 ^ 2;
t2147 = t1854 - t2129;
t1579 = t1618 * t1936 - t1939 * t2147;
t1581 = t1618 * t1939 + t1936 * t2147;
t1491 = t1579 * t1948 + t1581 * t1945;
t1938 = sin(pkin(6));
t1941 = cos(pkin(6));
t1492 = t1579 * t1945 - t1581 * t1948;
t2100 = t2144 * t1947;
t1614 = -t1944 * t2150 + t2100;
t1946 = sin(qJ(2));
t1949 = cos(qJ(2));
t2022 = t1492 * t1946 + t1614 * t1949;
t1407 = t1938 * t1491 + t1941 * t2022;
t1457 = t1492 * t1949 - t1614 * t1946;
t1937 = sin(pkin(10));
t1940 = cos(pkin(10));
t2281 = t1407 * t1937 - t1457 * t1940;
t2280 = t1407 * t1940 + t1457 * t1937;
t2128 = t1885 ^ 2;
t1811 = t2129 - t2128;
t1799 = t1857 * t1855;
t2046 = t1902 * t1936 - t1939 * t1983;
t1842 = qJDD(5) + t2046;
t2148 = t1799 + t1842;
t2170 = t2148 * t1944;
t1675 = t1811 * t1947 - t2170;
t1707 = t1779 - t1820;
t1599 = t1675 * t1936 + t1707 * t1939;
t1603 = t1675 * t1939 - t1707 * t1936;
t1510 = t1599 * t1948 + t1603 * t1945;
t1513 = t1599 * t1945 - t1603 * t1948;
t2169 = t2148 * t1947;
t1671 = t1811 * t1944 + t2169;
t2014 = t1513 * t1946 + t1671 * t1949;
t1432 = t1938 * t1510 + t1941 * t2014;
t1473 = t1513 * t1949 - t1671 * t1946;
t2279 = t1432 * t1937 - t1473 * t1940;
t2278 = t1432 * t1940 + t1473 * t1937;
t1428 = -t1941 * t1510 + t1938 * t2014;
t1405 = -t1941 * t1491 + t1938 * t2022;
t1812 = t1854 - t2128;
t2149 = -t1799 + t1842;
t2097 = t2149 * t1944;
t2191 = -t1812 * t1947 + t2097;
t2143 = t2086 + t1959;
t2096 = t2149 * t1947;
t2190 = t1812 * t1944 + t2096;
t2210 = t1936 * t2143 + t1939 * t2190;
t2213 = t1936 * t2190 - t1939 * t2143;
t2225 = -t1945 * t2213 + t1948 * t2210;
t2250 = t1946 * t2191 + t1949 * t2225;
t2224 = t1945 * t2210 + t1948 * t2213;
t2251 = t1946 * t2225 - t1949 * t2191;
t2262 = -t1938 * t2224 + t1941 * t2251;
t2269 = -t1937 * t2262 + t1940 * t2250;
t2268 = t1937 * t2250 + t1940 * t2262;
t2142 = -t2128 - t2129;
t2166 = t1947 * t2142 - t2097;
t2188 = t1936 * t2150 + t1939 * t2166;
t2189 = t1936 * t2166 - t1939 * t2150;
t2206 = t1945 * t2188 + t1948 * t2189;
t2167 = t1944 * t2142 + t2096;
t2207 = -t1945 * t2189 + t1948 * t2188;
t2227 = t1946 * t2207 - t1949 * t2167;
t2252 = -t1938 * t2206 + t1941 * t2227;
t2267 = pkin(1) * t2252;
t2253 = t1938 * t2227 + t1941 * t2206;
t2266 = pkin(1) * t2253;
t2263 = t1938 * t2251 + t1941 * t2224;
t2226 = t1946 * t2167 + t1949 * t2207;
t2261 = qJ(1) * (-t1937 * t2252 + t1940 * t2226);
t2260 = qJ(1) * (t1937 * t2226 + t1940 * t2252);
t2259 = (-t1938 * t2253 - t1941 * t2252) * pkin(7);
t2258 = pkin(7) * t2226;
t1776 = t2128 + t1854;
t1654 = t1776 * t1947 + t2170;
t2249 = pkin(2) * t1654;
t2248 = pkin(2) * t2206;
t2247 = pkin(3) * t1654;
t2246 = pkin(4) * t1654;
t2245 = pkin(8) * t2206;
t2244 = pkin(9) * t1654;
t1665 = t1776 * t1944 - t2169;
t2243 = pkin(9) * t1665;
t2241 = t1654 * t1949;
t2240 = t1665 * t1936;
t2239 = t1665 * t1939;
t2234 = t1946 * t1654;
t2232 = -pkin(2) * t2167 + pkin(8) * t2207;
t2084 = t1885 * t1947;
t1808 = t1857 * t2084;
t2085 = t1885 * t1944;
t2057 = t1855 * t2085;
t2033 = -t1808 - t2057;
t1807 = t1857 * t2085;
t2056 = t1855 * t2084;
t2034 = t1807 - t2056;
t2135 = t1842 * t1936 + t1939 * t2034;
t2139 = -t1842 * t1939 + t1936 * t2034;
t2164 = -t1945 * t2139 + t1948 * t2135;
t2185 = t1946 * t2033 + t1949 * t2164;
t2163 = t1945 * t2135 + t1948 * t2139;
t2187 = t1946 * t2164 - t1949 * t2033;
t2208 = -t1938 * t2163 + t1941 * t2187;
t2231 = -t1937 * t2208 + t1940 * t2185;
t1985 = -t1779 * t1947 + t2057;
t1986 = t1779 * t1944 + t2056;
t2059 = t1936 * t1799;
t2136 = t1939 * t1986 - t2059;
t2058 = t1939 * t1799;
t2137 = t1936 * t1986 + t2058;
t2162 = -t1945 * t2137 + t1948 * t2136;
t2184 = t1946 * t1985 + t1949 * t2162;
t2161 = t1945 * t2136 + t1948 * t2137;
t2186 = t1946 * t2162 - t1949 * t1985;
t2209 = -t1938 * t2161 + t1941 * t2186;
t2230 = -t1937 * t2209 + t1940 * t2184;
t2229 = t1937 * t2185 + t1940 * t2208;
t2228 = t1937 * t2184 + t1940 * t2209;
t2223 = pkin(3) * t2189;
t2222 = qJ(4) * t2189;
t1698 = -t1944 * t1959 - t1808;
t1699 = t1947 * t1959 - t1807;
t2035 = t1699 * t1939 + t2059;
t2036 = t1699 * t1936 - t2058;
t2133 = -t1945 * t2036 + t1948 * t2035;
t2215 = t1698 * t1949 + t1946 * t2133;
t2214 = -pkin(3) * t2167 + qJ(4) * t2188;
t2212 = t1938 * t2186 + t1941 * t2161;
t2211 = t1938 * t2187 + t1941 * t2163;
t2203 = pkin(4) * t2167;
t2202 = pkin(9) * t2166;
t2201 = pkin(9) * t2167;
t2200 = -qJ(6) * t1944 - pkin(4);
t2199 = qJ(6) * t2144;
t2132 = t1945 * t2035 + t1948 * t2036;
t2159 = -t1938 * t2132 + t1941 * t2215;
t2165 = -t1698 * t1946 + t1949 * t2133;
t2183 = t1937 * t2165 + t1940 * t2159;
t2182 = -t1937 * t2159 + t1940 * t2165;
t2146 = t1854 + t2129;
t2181 = pkin(4) * t2146;
t1843 = t1894 * t1892;
t2141 = qJDD(3) - t1843;
t2180 = t1936 * t2141;
t2178 = t1936 * t2146;
t2175 = t1939 * t2141;
t2173 = t1939 * t2146;
t2050 = g(1) * t1937 - g(2) * t1940;
t2120 = g(3) - qJDD(1);
t2168 = -t1938 * t2120 + t1941 * t2050;
t2081 = t1894 * qJD(3);
t1800 = t2046 + t2081;
t2160 = t1938 * t2215 + t1941 * t2132;
t2158 = t1937 * t2120;
t2157 = t1940 * t2120;
t2113 = qJD(4) * t1892;
t1880 = -0.2e1 * t2113;
t1872 = t1938 * t2050 + t1941 * t2120;
t1909 = g(1) * t1940 + g(2) * t1937;
t1823 = -t1949 * t1909 + t1946 * t2168;
t2131 = qJD(2) ^ 2;
t1951 = -pkin(2) * t2131 + qJDD(2) * pkin(8) + t1823;
t1765 = -t1872 * t1945 + t1948 * t1951;
t1910 = qJD(3) * pkin(3) - qJ(4) * t2116;
t2127 = t1948 ^ 2;
t1931 = t2127 * t2131;
t1722 = -pkin(3) * t1931 + qJ(4) * t1983 - qJD(3) * t1910 + t1765;
t1764 = t1948 * t1872 + t1945 * t1951;
t1921 = t1945 * t2131 * t1948;
t1911 = qJDD(3) + t1921;
t1950 = -t1764 + (-t1902 + t1924) * qJ(4) + t1911 * pkin(3);
t2067 = t1722 * t1939 + t1936 * t1950;
t1607 = t1880 + t2067;
t1827 = pkin(4) * t1892 - pkin(9) * t1894;
t2130 = qJD(3) ^ 2;
t1573 = -pkin(4) * t2130 + qJDD(3) * pkin(9) - t1827 * t1892 + t1607;
t1822 = -t1946 * t1909 - t1949 * t2168;
t1809 = -qJDD(2) * pkin(2) - pkin(8) * t2131 + t1822;
t1749 = -pkin(3) * t1983 - qJ(4) * t1931 + t1910 * t2116 + qJDD(4) + t1809;
t1883 = t1892 * qJD(3);
t1803 = -t1883 + t1845;
t1637 = pkin(4) * t1800 - pkin(9) * t1803 + t1749;
t1522 = t1573 * t1944 - t1637 * t1947;
t1523 = t1573 * t1947 + t1637 * t1944;
t1440 = t1522 * t1944 + t1523 * t1947;
t2111 = qJD(6) * t1885;
t1873 = 0.2e1 * t2111;
t1794 = pkin(5) * t1855 - qJ(6) * t1857;
t2028 = -pkin(5) * t2128 + qJ(6) * t1842 - t1794 * t1855 + t1523;
t1482 = t1873 + t2028;
t1487 = -pkin(5) * t1842 - qJ(6) * t2128 + t1794 * t1857 + qJDD(6) + t1522;
t1420 = t1482 * t1947 + t1487 * t1944;
t2048 = t1936 * t1722 - t1939 * t1950;
t1989 = -qJDD(3) * pkin(4) - pkin(9) * t2130 + t2048;
t1572 = (0.2e1 * qJD(4) + t1827) * t1894 + t1989;
t2140 = t1779 * pkin(5) - t2199;
t1505 = (pkin(5) * t1885 - 0.2e1 * qJD(6)) * t1857 + t1572 + t2140;
t1390 = t1420 * t1936 - t1505 * t1939;
t1956 = pkin(9) * t1420 + (-pkin(5) * t1947 + t2200) * t1505;
t2145 = pkin(3) * t1390 + t1956;
t2138 = -t1909 * t1940 - t1937 * t2050;
t2134 = -t1937 * t1909 + t1940 * t2050;
t1886 = t1892 ^ 2;
t1887 = t1894 ^ 2;
t2112 = qJD(4) * t1894;
t1606 = t2048 + 0.2e1 * t2112;
t1525 = -t1606 * t1939 + t1607 * t1936;
t2125 = pkin(3) * t1525;
t1804 = t1883 + t1845;
t1984 = -t2046 + t2081;
t1719 = -t1804 * t1939 + t1936 * t1984;
t2124 = pkin(3) * t1719;
t1830 = -t1886 - t2130;
t1751 = t1830 * t1936 + t2175;
t2123 = pkin(3) * t1751;
t2122 = pkin(4) * t1936;
t1743 = t1822 * t1946 + t1823 * t1949;
t2121 = pkin(7) * t1743;
t2118 = qJ(6) * t1947;
t2117 = qJD(2) * qJD(3);
t1933 = t1945 ^ 2;
t2114 = t2131 * t1933;
t2110 = t1525 * t1945;
t2109 = t1525 * t1948;
t2103 = t2143 * t1944;
t2102 = t2143 * t1947;
t2095 = t1749 * t1936;
t2094 = t1749 * t1939;
t2093 = t1809 * t1945;
t2092 = t1809 * t1948;
t1833 = qJDD(3) + t1843;
t2089 = t1833 * t1936;
t2088 = t1833 * t1939;
t2083 = t1892 * t1936;
t2082 = t1892 * t1939;
t2080 = t1894 * t1936;
t2079 = t1894 * t1939;
t1903 = -0.2e1 * t2054 + t2062;
t1858 = t1903 * t1948;
t2078 = t1911 * t1945;
t1912 = qJDD(3) - t1921;
t2077 = t1912 * t1945;
t2076 = t1912 * t1948;
t1568 = t1944 * t1572;
t2071 = t1946 * t1872;
t1569 = t1947 * t1572;
t2070 = t1949 * t1872;
t2068 = -pkin(4) * t1572 + pkin(9) * t1440;
t2066 = qJDD(3) * t1949;
t2065 = t1938 * qJDD(2);
t2063 = t1946 * qJDD(2);
t2061 = t1933 + t2127;
t2060 = -pkin(4) * t1939 - pkin(3);
t2055 = t1949 * t1843;
t2053 = t1946 * t1843;
t1714 = (qJD(5) + t1885) * t1855 + t1993;
t2052 = pkin(4) * t1714 + t1568 + t2243;
t2051 = -pkin(4) * t2150 - t1569 + t2202;
t1526 = t1606 * t1936 + t1607 * t1939;
t1660 = t1764 * t1945 + t1765 * t1948;
t2045 = t1949 * t1921;
t2044 = t1946 * t1921;
t1462 = pkin(5) * t2146 + t1482;
t1466 = qJ(6) * t2146 + t1487;
t1617 = -t1707 * t1947 + t2103;
t2043 = pkin(9) * t1617 + t1462 * t1947 + t1466 * t1944 + t2181;
t1709 = (-qJD(5) + t1885) * t1857 + t2047;
t1616 = t1709 * t1947 + t2103;
t2042 = pkin(9) * t1616 + t1440 + t2181;
t1416 = t1440 * t1936 - t1572 * t1939;
t2041 = pkin(3) * t1416 + t2068;
t1877 = -t1887 - t2130;
t1784 = t1877 * t1939 - t2089;
t2040 = pkin(3) * t1784 - t2067;
t2039 = -pkin(5) * t1487 + qJ(6) * t1482;
t2038 = -pkin(5) * t2143 - qJ(6) * t1707;
t1659 = t1764 * t1948 - t1765 * t1945;
t1904 = t2061 * qJDD(2);
t1907 = t1931 + t2114;
t1848 = t1904 * t1949 - t1907 * t1946;
t2031 = pkin(7) * t1848 + t1659 * t1946;
t1905 = qJDD(2) * t1949 - t1946 * t2131;
t2030 = -pkin(7) * t1905 - t2071;
t1995 = t1949 * t2131 + t2063;
t2029 = -pkin(7) * t1995 + t2070;
t1391 = t1420 * t1939 + t1505 * t1936;
t1359 = -t1390 * t1945 + t1391 * t1948;
t1419 = t1482 * t1944 - t1487 * t1947;
t2027 = t1359 * t1946 - t1419 * t1949;
t1417 = t1440 * t1939 + t1572 * t1936;
t1369 = -t1416 * t1945 + t1417 * t1948;
t1439 = -t1522 * t1947 + t1523 * t1944;
t2026 = t1369 * t1946 - t1439 * t1949;
t1442 = t1526 * t1948 - t2110;
t2025 = t1442 * t1946 - t1749 * t1949;
t1564 = t1616 * t1936 + t2173;
t1566 = t1616 * t1939 - t2178;
t1485 = -t1564 * t1945 + t1566 * t1948;
t1612 = t1709 * t1944 - t2102;
t2024 = t1485 * t1946 - t1612 * t1949;
t1565 = t1617 * t1936 + t2173;
t1567 = t1617 * t1939 - t2178;
t1486 = -t1565 * t1945 + t1567 * t1948;
t1613 = -t1707 * t1944 - t2102;
t2023 = t1486 * t1946 - t1613 * t1949;
t1583 = t1939 * t2144 - t2240;
t1585 = -t1936 * t2144 - t2239;
t1497 = -t1583 * t1945 + t1585 * t1948;
t2020 = t1497 * t1946 - t2241;
t1591 = t1714 * t1939 + t2240;
t1593 = -t1714 * t1936 + t2239;
t1503 = -t1591 * t1945 + t1593 * t1948;
t2018 = t1503 * t1946 + t2241;
t1718 = -t1800 * t1936 + t1803 * t1939;
t1720 = -t1800 * t1939 - t1803 * t1936;
t1621 = -t1718 * t1945 + t1720 * t1948;
t1839 = t1887 - t1886;
t2008 = t1621 * t1946 - t1839 * t1949;
t1721 = t1804 * t1936 + t1939 * t1984;
t1622 = -t1719 * t1945 + t1721 * t1948;
t1797 = -t1886 - t1887;
t2007 = t1622 * t1946 - t1797 * t1949;
t1752 = t1830 * t1939 - t2180;
t1649 = -t1751 * t1945 + t1752 * t1948;
t2006 = t1649 * t1946 - t1800 * t1949;
t2005 = t1660 * t1946 - t1809 * t1949;
t1875 = t1886 - t2130;
t1782 = t1875 * t1936 + t2088;
t1785 = t1875 * t1939 - t2089;
t1683 = -t1782 * t1945 + t1785 * t1948;
t2004 = t1683 * t1946 - t1949 * t1984;
t1876 = -t1887 + t2130;
t1783 = t1876 * t1939 + t2180;
t1786 = -t1876 * t1936 + t2175;
t1684 = -t1783 * t1945 + t1786 * t1948;
t2003 = t1684 * t1946 - t1804 * t1949;
t1787 = -t1877 * t1936 - t2088;
t1685 = -t1784 * t1945 + t1787 * t1948;
t2002 = t1685 * t1946 - t1803 * t1949;
t1742 = t1822 * t1949 - t1823 * t1946;
t1901 = 0.2e1 * t1924 + t2064;
t1847 = -t1901 * t1945 + t1858;
t1908 = -t1931 + t2114;
t2001 = t1847 * t1946 - t1908 * t1949;
t1920 = -t1931 - t2130;
t1865 = t1920 * t1948 - t2078;
t2000 = t1865 * t1946 + t1903 * t1949;
t1918 = -t2114 - t2130;
t1867 = -t1918 * t1945 - t2076;
t1999 = t1867 * t1946 - t1901 * t1949;
t1890 = t1995 * t1941;
t1998 = t1890 * t1940 + t1905 * t1937;
t1997 = t1890 * t1937 - t1905 * t1940;
t1996 = t1904 * t1946 + t1907 * t1949;
t1816 = (-t2079 - t2083) * qJD(3);
t1817 = (t2080 - t2082) * qJD(3);
t1731 = -t1816 * t1945 + t1817 * t1948;
t1994 = t1731 * t1946 - t2066;
t1899 = t2061 * t2117;
t1992 = t1899 * t1946 - t2066;
t1991 = pkin(3) * t1591 + t2052;
t1990 = t2051 + t2223;
t1790 = qJD(3) * t2083 - t1939 * t2046;
t1791 = qJD(3) * t2082 + t1936 * t2046;
t1690 = -t1790 * t1945 + t1791 * t1948;
t1988 = t1690 * t1946 + t2055;
t1792 = qJD(3) * t2079 + t1845 * t1936;
t1793 = -qJD(3) * t2080 + t1845 * t1939;
t1691 = -t1792 * t1945 + t1793 * t1948;
t1987 = t1691 * t1946 - t2055;
t1919 = t1931 - t2130;
t1864 = t1919 * t1948 - t2077;
t1982 = t1864 * t1946 - t1949 * t2062;
t1900 = t1948 * t1911;
t1917 = -t2114 + t2130;
t1866 = -t1917 * t1945 + t1900;
t1981 = t1866 * t1946 - t1949 * t2064;
t1882 = -0.2e1 * t2112;
t1953 = 0.2e1 * qJD(6) * t1857 - t1894 * t1827 + t1882 - t1989 - t2140;
t1469 = -pkin(5) * t1820 + t1953 + t2199;
t1980 = pkin(4) * t2144 + pkin(5) * t2100 + t1469 * t1944 - t2243;
t1979 = pkin(3) * t1565 + t2043;
t1978 = pkin(3) * t1564 + t2042;
t1470 = (-t2150 - t1820) * pkin(5) + t1953;
t1977 = t1947 * t1470 + t2150 * t2200 + t2202;
t1869 = -t1945 * t1983 - t2117 * t2127;
t1976 = t1869 * t1946 - t2045;
t1870 = t1902 * t1948 - t1933 * t2117;
t1975 = t1870 * t1946 + t2045;
t1370 = -pkin(4) * t1419 - t2039;
t1372 = -pkin(9) * t1419 + (pkin(5) * t1944 - t2118) * t1505;
t1337 = -pkin(3) * t1419 + qJ(4) * t1391 + t1370 * t1939 + t1372 * t1936;
t1341 = -qJ(4) * t1390 - t1370 * t1936 + t1372 * t1939;
t1358 = t1390 * t1948 + t1391 * t1945;
t1313 = -pkin(8) * t1358 - t1337 * t1945 + t1341 * t1948;
t1336 = -pkin(2) * t1358 - t2145;
t1351 = t1359 * t1949 + t1419 * t1946;
t1974 = pkin(7) * t1351 + t1313 * t1946 + t1336 * t1949;
t1354 = qJ(4) * t1417 + (-pkin(9) * t1936 + t2060) * t1439;
t1364 = -qJ(4) * t1416 + (-pkin(9) * t1939 + t2122) * t1439;
t1368 = t1416 * t1948 + t1417 * t1945;
t1328 = -pkin(8) * t1368 - t1354 * t1945 + t1364 * t1948;
t1347 = -pkin(2) * t1368 - t2041;
t1361 = t1369 * t1949 + t1439 * t1946;
t1973 = pkin(7) * t1361 + t1328 * t1946 + t1347 * t1949;
t1393 = -pkin(9) * t1613 - t1462 * t1944 + t1466 * t1947;
t1532 = -pkin(4) * t1613 - t2038;
t1367 = -pkin(3) * t1613 + qJ(4) * t1567 + t1393 * t1936 + t1532 * t1939;
t1371 = -qJ(4) * t1565 + t1393 * t1939 - t1532 * t1936;
t1484 = t1565 * t1948 + t1567 * t1945;
t1345 = -pkin(8) * t1484 - t1367 * t1945 + t1371 * t1948;
t1366 = -pkin(2) * t1484 - t1979;
t1448 = t1486 * t1949 + t1613 * t1946;
t1972 = pkin(7) * t1448 + t1345 * t1946 + t1366 * t1949;
t1954 = pkin(5) * t1776 + qJ(6) * t2148 + t2028;
t1443 = -t1954 - 0.2e1 * t2111 - t2246;
t1444 = -pkin(5) * t2101 + t1469 * t1947 - t2244;
t1373 = qJ(4) * t1585 + t1443 * t1939 + t1444 * t1936 - t2247;
t1378 = -qJ(4) * t1583 - t1443 * t1936 + t1444 * t1939;
t1495 = t1583 * t1948 + t1585 * t1945;
t1349 = -pkin(8) * t1495 - t1373 * t1945 + t1378 * t1948;
t1958 = pkin(3) * t1583 + t1980;
t1387 = -pkin(2) * t1495 - t1958;
t1459 = t1497 * t1949 + t2234;
t1971 = pkin(7) * t1459 + t1349 * t1946 + t1387 * t1949;
t1445 = -t1470 * t1944 - t2118 * t2150 - t2201;
t1952 = pkin(5) * t2149 + qJ(6) * t2142 - t1487;
t1446 = -t1952 - t2203;
t1374 = t1445 * t1936 + t1446 * t1939 + t2214;
t1379 = t1445 * t1939 - t1446 * t1936 - t2222;
t1350 = -t1374 * t1945 + t1379 * t1948 - t2245;
t1955 = t1977 + t2223;
t1388 = -t1955 - t2248;
t1970 = t1350 * t1946 + t1388 * t1949 + t2258;
t1425 = -pkin(9) * t1612 - t1439;
t1386 = qJ(4) * t1566 + t1936 * t1425 + t1612 * t2060;
t1392 = -qJ(4) * t1564 + t1425 * t1939 + t1612 * t2122;
t1483 = t1564 * t1948 + t1566 * t1945;
t1353 = -pkin(8) * t1483 - t1386 * t1945 + t1392 * t1948;
t1377 = -pkin(2) * t1483 - t1978;
t1447 = t1485 * t1949 + t1612 * t1946;
t1969 = pkin(7) * t1447 + t1353 * t1946 + t1377 * t1949;
t1479 = t1522 - t2203;
t1531 = t1568 - t2201;
t1395 = t1479 * t1939 + t1531 * t1936 + t2214;
t1409 = -t1479 * t1936 + t1531 * t1939 - t2222;
t1362 = -t1395 * t1945 + t1409 * t1948 - t2245;
t1404 = -t1990 - t2248;
t1968 = t1362 * t1946 + t1404 * t1949 + t2258;
t1481 = t1523 + t2246;
t1533 = t1569 + t2244;
t1396 = qJ(4) * t1593 + t1481 * t1939 + t1533 * t1936 + t2247;
t1411 = -qJ(4) * t1591 - t1481 * t1936 + t1533 * t1939;
t1501 = t1591 * t1948 + t1593 * t1945;
t1363 = -pkin(8) * t1501 - t1396 * t1945 + t1411 * t1948;
t1410 = -pkin(2) * t1501 - t1991;
t1463 = t1503 * t1949 - t2234;
t1967 = pkin(7) * t1463 + t1363 * t1946 + t1410 * t1949;
t1441 = t1526 * t1945 + t2109;
t1494 = -pkin(3) * t1749 + qJ(4) * t1526;
t1383 = -pkin(8) * t1441 - qJ(4) * t2109 - t1494 * t1945;
t1403 = -pkin(2) * t1441 - t2125;
t1437 = t1442 * t1949 + t1749 * t1946;
t1966 = pkin(7) * t1437 + t1383 * t1946 + t1403 * t1949;
t1480 = -pkin(3) * t1797 + qJ(4) * t1721 + t1526;
t1489 = -qJ(4) * t1719 - t1525;
t1620 = t1719 * t1948 + t1721 * t1945;
t1402 = -pkin(8) * t1620 - t1480 * t1945 + t1489 * t1948;
t1557 = -pkin(2) * t1620 - t2124;
t1582 = t1622 * t1949 + t1797 * t1946;
t1965 = pkin(7) * t1582 + t1402 * t1946 + t1557 * t1949;
t1624 = -pkin(3) * t1800 + qJ(4) * t1752 - t2094;
t1647 = -qJ(4) * t1751 + t2095;
t1648 = t1751 * t1948 + t1752 * t1945;
t1499 = -pkin(8) * t1648 - t1624 * t1945 + t1647 * t1948;
t1524 = -pkin(2) * t1648 + t1606 - t2123;
t1623 = t1649 * t1949 + t1800 * t1946;
t1964 = pkin(7) * t1623 + t1499 * t1946 + t1524 * t1949;
t1631 = -pkin(3) * t1803 + qJ(4) * t1787 + t2095;
t1667 = -qJ(4) * t1784 + t2094;
t1682 = t1784 * t1948 + t1787 * t1945;
t1519 = -pkin(8) * t1682 - t1631 * t1945 + t1667 * t1948;
t1534 = -pkin(2) * t1682 + t1880 - t2040;
t1634 = t1685 * t1949 + t1803 * t1946;
t1963 = pkin(7) * t1634 + t1519 * t1946 + t1534 * t1949;
t1861 = t1920 * t1945 + t1900;
t1728 = -pkin(2) * t1861 + t1764;
t1762 = -pkin(8) * t1861 + t2093;
t1814 = t1865 * t1949 - t1903 * t1946;
t1962 = pkin(7) * t1814 + t1728 * t1949 + t1762 * t1946;
t1863 = t1918 * t1948 - t2077;
t1729 = -pkin(2) * t1863 + t1765;
t1763 = -pkin(8) * t1863 + t2092;
t1815 = t1867 * t1949 + t1901 * t1946;
t1961 = pkin(7) * t1815 + t1729 * t1949 + t1763 * t1946;
t1627 = t1660 * t1949 + t1809 * t1946;
t1957 = pkin(7) * t1627 - (-pkin(2) * t1949 - pkin(8) * t1946) * t1659;
t1927 = t1946 * qJDD(3);
t1926 = t1941 * qJDD(2);
t1891 = t1905 * t1941;
t1889 = t1905 * t1938;
t1888 = t1995 * t1938;
t1871 = t1899 * t1949 + t1927;
t1862 = t1917 * t1948 + t2078;
t1860 = t1919 * t1945 + t2076;
t1859 = (t1902 + t1924) * t1945;
t1850 = t1992 * t1941;
t1849 = t1992 * t1938;
t1846 = t1901 * t1948 + t1903 * t1945;
t1841 = t1996 * t1941;
t1840 = t1996 * t1938;
t1837 = -t1891 * t1937 - t1940 * t1995;
t1836 = t1891 * t1940 - t1937 * t1995;
t1829 = t1870 * t1949 - t2044;
t1828 = t1869 * t1949 + t2044;
t1826 = t1866 * t1949 + t1945 * t2063;
t1825 = t1864 * t1949 + t1946 * t2062;
t1806 = t1847 * t1949 + t1908 * t1946;
t1781 = -t2070 + (t1888 * t1938 + t1890 * t1941) * pkin(7);
t1780 = -t2071 + (-t1889 * t1938 - t1891 * t1941) * pkin(7);
t1774 = -t1938 * t1859 + t1941 * t1975;
t1773 = -t1938 * t1858 + t1941 * t1976;
t1772 = t1941 * t1859 + t1938 * t1975;
t1771 = t1941 * t1858 + t1938 * t1976;
t1770 = -t1938 * t1862 + t1941 * t1981;
t1769 = -t1938 * t1860 + t1941 * t1982;
t1768 = t1941 * t1862 + t1938 * t1981;
t1767 = t1941 * t1860 + t1938 * t1982;
t1758 = -t1938 * t1863 + t1941 * t1999;
t1757 = -t1938 * t1861 + t1941 * t2000;
t1756 = t1941 * t1863 + t1938 * t1999;
t1755 = t1941 * t1861 + t1938 * t2000;
t1741 = -t1938 * t1846 + t1941 * t2001;
t1740 = t1941 * t1846 + t1938 * t2001;
t1739 = pkin(2) * t1903 + pkin(8) * t1865 - t2092;
t1738 = -pkin(2) * t1901 + pkin(8) * t1867 + t2093;
t1737 = t1743 * t1941;
t1734 = t1743 * t1938;
t1730 = t1816 * t1948 + t1817 * t1945;
t1727 = t1731 * t1949 + t1927;
t1726 = -pkin(1) * t1889 + t1938 * t1822 + t1941 * t2029;
t1725 = pkin(1) * t1888 + t1938 * t1823 + t1941 * t2030;
t1724 = pkin(1) * t1891 - t1941 * t1822 + t1938 * t2029;
t1723 = -pkin(1) * t1890 - t1941 * t1823 + t1938 * t2030;
t1701 = -t1742 * t1941 + t1938 * t1872;
t1700 = -t1742 * t1938 - t1941 * t1872;
t1689 = t1792 * t1948 + t1793 * t1945;
t1688 = t1790 * t1948 + t1791 * t1945;
t1681 = t1783 * t1948 + t1786 * t1945;
t1680 = t1782 * t1948 + t1785 * t1945;
t1653 = t1691 * t1949 + t2053;
t1652 = t1690 * t1949 - t2053;
t1638 = pkin(2) * t1907 + pkin(8) * t1904 + t1660;
t1633 = t1684 * t1949 + t1804 * t1946;
t1632 = t1683 * t1949 + t1946 * t1984;
t1630 = -pkin(2) * t1809 + pkin(8) * t1660;
t1629 = -t1938 * t1730 + t1941 * t1994;
t1628 = t1730 * t1941 + t1938 * t1994;
t1626 = pkin(1) * t1701 + t1938 * t2121;
t1625 = -pkin(1) * t1700 + t1941 * t2121;
t1619 = t1718 * t1948 + t1720 * t1945;
t1608 = t1949 * t1659 + (-t1840 * t1938 - t1841 * t1941) * pkin(7);
t1604 = (-t1700 * t1938 - t1701 * t1941) * pkin(7);
t1595 = t1621 * t1949 + t1839 * t1946;
t1577 = -t1938 * t1689 + t1941 * t1987;
t1576 = -t1938 * t1688 + t1941 * t1988;
t1575 = t1941 * t1689 + t1938 * t1987;
t1574 = t1941 * t1688 + t1938 * t1988;
t1563 = -t1938 * t1682 + t1941 * t2002;
t1562 = -t1938 * t1681 + t1941 * t2003;
t1561 = -t1938 * t1680 + t1941 * t2004;
t1560 = t1941 * t1682 + t1938 * t2002;
t1559 = t1941 * t1681 + t1938 * t2003;
t1558 = t1941 * t1680 + t1938 * t2004;
t1556 = -t1946 * t1729 + t1949 * t1763 + (-t1756 * t1938 - t1758 * t1941) * pkin(7);
t1555 = -t1946 * t1728 + t1949 * t1762 + (-t1755 * t1938 - t1757 * t1941) * pkin(7);
t1546 = t1659 * t1938 + t1941 * t2005;
t1545 = -t1659 * t1941 + t1938 * t2005;
t1542 = -t1938 * t1648 + t1941 * t2006;
t1541 = t1941 * t1648 + t1938 * t2006;
t1540 = -pkin(1) * t1756 - t1938 * t1738 + t1941 * t1961;
t1539 = -pkin(1) * t1755 - t1938 * t1739 + t1941 * t1962;
t1538 = pkin(1) * t1758 + t1941 * t1738 + t1938 * t1961;
t1537 = pkin(1) * t1757 + t1941 * t1739 + t1938 * t1962;
t1536 = -pkin(1) * t1840 - t1938 * t1638 + t1941 * t2031;
t1535 = pkin(1) * t1841 + t1941 * t1638 + t1938 * t2031;
t1518 = -t1938 * t1619 + t1941 * t2008;
t1517 = t1941 * t1619 + t1938 * t2008;
t1516 = -t1938 * t1620 + t1941 * t2007;
t1515 = t1941 * t1620 + t1938 * t2007;
t1506 = -pkin(2) * t1803 + pkin(8) * t1685 + t1631 * t1948 + t1667 * t1945;
t1488 = -pkin(2) * t1800 + pkin(8) * t1649 + t1624 * t1948 + t1647 * t1945;
t1436 = -(pkin(2) * t1946 - pkin(8) * t1949) * t1659 + (-t1545 * t1938 - t1546 * t1941) * pkin(7);
t1435 = -pkin(1) * t1545 - t1938 * t1630 + t1941 * t1957;
t1434 = pkin(1) * t1546 + t1941 * t1630 + t1938 * t1957;
t1423 = -t1938 * t1501 + t1941 * t2018;
t1421 = t1941 * t1501 + t1938 * t2018;
t1414 = -t1938 * t1495 + t1941 * t2020;
t1412 = t1941 * t1495 + t1938 * t2020;
t1401 = -t1938 * t1484 + t1941 * t2023;
t1400 = -t1938 * t1483 + t1941 * t2024;
t1399 = t1941 * t1484 + t1938 * t2023;
t1398 = t1941 * t1483 + t1938 * t2024;
t1397 = -pkin(2) * t1797 + pkin(8) * t1622 + t1480 * t1948 + t1489 * t1945;
t1394 = t1949 * t1519 - t1946 * t1534 + (-t1560 * t1938 - t1563 * t1941) * pkin(7);
t1389 = t1949 * t1499 - t1946 * t1524 + (-t1541 * t1938 - t1542 * t1941) * pkin(7);
t1385 = -pkin(1) * t1560 - t1938 * t1506 + t1941 * t1963;
t1384 = pkin(1) * t1563 + t1941 * t1506 + t1938 * t1963;
t1382 = -t1938 * t1441 + t1941 * t2025;
t1381 = t1941 * t1441 + t1938 * t2025;
t1380 = -pkin(2) * t1749 + pkin(8) * t1442 - qJ(4) * t2110 + t1494 * t1948;
t1376 = -pkin(1) * t1541 - t1938 * t1488 + t1941 * t1964;
t1375 = pkin(1) * t1542 + t1941 * t1488 + t1938 * t1964;
t1365 = t1949 * t1402 - t1946 * t1557 + (-t1515 * t1938 - t1516 * t1941) * pkin(7);
t1360 = pkin(8) * t1503 + t1396 * t1948 + t1411 * t1945 + t2249;
t1357 = t1395 * t1948 + t1409 * t1945 + t2232;
t1356 = -pkin(1) * t1515 - t1938 * t1397 + t1941 * t1965;
t1355 = pkin(1) * t1516 + t1941 * t1397 + t1938 * t1965;
t1352 = -pkin(2) * t1612 + pkin(8) * t1485 + t1386 * t1948 + t1392 * t1945;
t1348 = t1374 * t1948 + t1379 * t1945 + t2232;
t1346 = pkin(8) * t1497 + t1373 * t1948 + t1378 * t1945 - t2249;
t1344 = -pkin(2) * t1613 + pkin(8) * t1486 + t1367 * t1948 + t1371 * t1945;
t1343 = -t1938 * t1368 + t1941 * t2026;
t1342 = t1941 * t1368 + t1938 * t2026;
t1340 = t1949 * t1363 - t1946 * t1410 + (-t1421 * t1938 - t1423 * t1941) * pkin(7);
t1339 = t1949 * t1362 - t1946 * t1404 + t2259;
t1338 = t1949 * t1383 - t1946 * t1403 + (-t1381 * t1938 - t1382 * t1941) * pkin(7);
t1335 = -t1938 * t1358 + t1941 * t2027;
t1334 = t1941 * t1358 + t1938 * t2027;
t1333 = -pkin(1) * t1381 - t1938 * t1380 + t1941 * t1966;
t1332 = pkin(1) * t1382 + t1941 * t1380 + t1938 * t1966;
t1331 = t1949 * t1350 - t1946 * t1388 + t2259;
t1330 = t1949 * t1349 - t1946 * t1387 + (-t1412 * t1938 - t1414 * t1941) * pkin(7);
t1329 = t1949 * t1353 - t1946 * t1377 + (-t1398 * t1938 - t1400 * t1941) * pkin(7);
t1327 = -pkin(1) * t1421 - t1938 * t1360 + t1941 * t1967;
t1326 = pkin(1) * t1423 + t1941 * t1360 + t1938 * t1967;
t1325 = -t1938 * t1357 + t1941 * t1968 - t2266;
t1324 = t1941 * t1357 + t1938 * t1968 + t2267;
t1323 = -pkin(2) * t1439 + pkin(8) * t1369 + t1354 * t1948 + t1364 * t1945;
t1322 = t1949 * t1345 - t1946 * t1366 + (-t1399 * t1938 - t1401 * t1941) * pkin(7);
t1321 = -pkin(1) * t1398 - t1938 * t1352 + t1941 * t1969;
t1320 = pkin(1) * t1400 + t1941 * t1352 + t1938 * t1969;
t1319 = -t1938 * t1348 + t1941 * t1970 - t2266;
t1318 = t1941 * t1348 + t1938 * t1970 + t2267;
t1317 = -pkin(1) * t1412 - t1938 * t1346 + t1941 * t1971;
t1316 = pkin(1) * t1414 + t1941 * t1346 + t1938 * t1971;
t1315 = -pkin(1) * t1399 - t1938 * t1344 + t1941 * t1972;
t1314 = pkin(1) * t1401 + t1941 * t1344 + t1938 * t1972;
t1312 = -pkin(2) * t1419 + pkin(8) * t1359 + t1337 * t1948 + t1341 * t1945;
t1311 = t1949 * t1328 - t1946 * t1347 + (-t1342 * t1938 - t1343 * t1941) * pkin(7);
t1310 = -pkin(1) * t1342 - t1938 * t1323 + t1941 * t1973;
t1309 = pkin(1) * t1343 + t1941 * t1323 + t1938 * t1973;
t1308 = t1949 * t1313 - t1946 * t1336 + (-t1334 * t1938 - t1335 * t1941) * pkin(7);
t1307 = -pkin(1) * t1334 - t1938 * t1312 + t1941 * t1974;
t1306 = pkin(1) * t1335 + t1941 * t1312 + t1938 * t1974;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t2158, -t2157, -t2134, -qJ(1) * t2134, 0, 0, -t1997, 0, t1837, t1937 * t2065, -qJ(1) * t1836 - t1726 * t1937 + t1780 * t1940, qJ(1) * t1998 - t1937 * t1725 + t1940 * t1781, -t1737 * t1937 + t1742 * t1940, t1940 * t1604 - t1937 * t1625 - qJ(1) * (t1701 * t1940 + t1743 * t1937), -t1774 * t1937 + t1829 * t1940, -t1741 * t1937 + t1806 * t1940, -t1770 * t1937 + t1826 * t1940, -t1773 * t1937 + t1828 * t1940, -t1769 * t1937 + t1825 * t1940, -t1850 * t1937 + t1871 * t1940, t1940 * t1555 - t1937 * t1539 - qJ(1) * (t1757 * t1940 + t1814 * t1937), t1940 * t1556 - t1937 * t1540 - qJ(1) * (t1758 * t1940 + t1815 * t1937), t1940 * t1608 - t1937 * t1536 - qJ(1) * (t1841 * t1940 + t1848 * t1937), t1940 * t1436 - t1937 * t1435 - qJ(1) * (t1546 * t1940 + t1627 * t1937), -t1577 * t1937 + t1653 * t1940, -t1518 * t1937 + t1595 * t1940, -t1562 * t1937 + t1633 * t1940, -t1576 * t1937 + t1652 * t1940, -t1561 * t1937 + t1632 * t1940, -t1629 * t1937 + t1727 * t1940, t1940 * t1389 - t1937 * t1376 - qJ(1) * (t1542 * t1940 + t1623 * t1937), t1940 * t1394 - t1937 * t1385 - qJ(1) * (t1563 * t1940 + t1634 * t1937), t1940 * t1365 - t1937 * t1356 - qJ(1) * (t1516 * t1940 + t1582 * t1937), t1940 * t1338 - t1937 * t1333 - qJ(1) * (t1382 * t1940 + t1437 * t1937), t2182, t2281, t2269, t2230, t2279, t2231, -t1937 * t1325 + t1940 * t1339 - t2260, t1940 * t1340 - t1937 * t1327 - qJ(1) * (t1423 * t1940 + t1463 * t1937), t1940 * t1329 - t1937 * t1321 - qJ(1) * (t1400 * t1940 + t1447 * t1937), t1940 * t1311 - t1937 * t1310 - qJ(1) * (t1343 * t1940 + t1361 * t1937), t2182, t2269, -t2281, t2231, -t2279, t2230, -t1937 * t1319 + t1940 * t1331 - t2260, t1940 * t1322 - t1937 * t1315 - qJ(1) * (t1401 * t1940 + t1448 * t1937), t1940 * t1330 - t1937 * t1317 - qJ(1) * (t1414 * t1940 + t1459 * t1937), t1940 * t1308 - t1937 * t1307 - qJ(1) * (t1335 * t1940 + t1351 * t1937); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t2157, -t2158, t2138, qJ(1) * t2138, 0, 0, t1998, 0, t1836, -t1940 * t2065, qJ(1) * t1837 + t1726 * t1940 + t1780 * t1937, qJ(1) * t1997 + t1940 * t1725 + t1937 * t1781, t1737 * t1940 + t1742 * t1937, t1937 * t1604 + t1940 * t1625 + qJ(1) * (-t1701 * t1937 + t1743 * t1940), t1774 * t1940 + t1829 * t1937, t1741 * t1940 + t1806 * t1937, t1770 * t1940 + t1826 * t1937, t1773 * t1940 + t1828 * t1937, t1769 * t1940 + t1825 * t1937, t1850 * t1940 + t1871 * t1937, t1937 * t1555 + t1940 * t1539 + qJ(1) * (-t1757 * t1937 + t1814 * t1940), t1937 * t1556 + t1940 * t1540 + qJ(1) * (-t1758 * t1937 + t1815 * t1940), t1937 * t1608 + t1940 * t1536 + qJ(1) * (-t1841 * t1937 + t1848 * t1940), t1937 * t1436 + t1940 * t1435 + qJ(1) * (-t1546 * t1937 + t1627 * t1940), t1577 * t1940 + t1653 * t1937, t1518 * t1940 + t1595 * t1937, t1562 * t1940 + t1633 * t1937, t1576 * t1940 + t1652 * t1937, t1561 * t1940 + t1632 * t1937, t1629 * t1940 + t1727 * t1937, t1937 * t1389 + t1940 * t1376 + qJ(1) * (-t1542 * t1937 + t1623 * t1940), t1937 * t1394 + t1940 * t1385 + qJ(1) * (-t1563 * t1937 + t1634 * t1940), t1937 * t1365 + t1940 * t1356 + qJ(1) * (-t1516 * t1937 + t1582 * t1940), t1937 * t1338 + t1940 * t1333 + qJ(1) * (-t1382 * t1937 + t1437 * t1940), t2183, -t2280, t2268, t2228, -t2278, t2229, t1940 * t1325 + t1937 * t1339 + t2261, t1937 * t1340 + t1940 * t1327 + qJ(1) * (-t1423 * t1937 + t1463 * t1940), t1937 * t1329 + t1940 * t1321 + qJ(1) * (-t1400 * t1937 + t1447 * t1940), t1937 * t1311 + t1940 * t1310 + qJ(1) * (-t1343 * t1937 + t1361 * t1940), t2183, t2268, t2280, t2229, t2278, t2228, t1940 * t1319 + t1937 * t1331 + t2261, t1937 * t1322 + t1940 * t1315 + qJ(1) * (-t1401 * t1937 + t1448 * t1940), t1937 * t1330 + t1940 * t1317 + qJ(1) * (-t1414 * t1937 + t1459 * t1940), t1937 * t1308 + t1940 * t1307 + qJ(1) * (-t1335 * t1937 + t1351 * t1940); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t2050, t1909, 0, 0, 0, 0, t1888, 0, t1889, t1926, t1724, t1723, t1734, t1626, t1772, t1740, t1768, t1771, t1767, t1849, t1537, t1538, t1535, t1434, t1575, t1517, t1559, t1574, t1558, t1628, t1375, t1384, t1355, t1332, t2160, -t1405, t2263, t2212, -t1428, t2211, t1324, t1326, t1320, t1309, t2160, t2263, t1405, t2211, t1428, t2212, t1318, t1314, t1316, t1306; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2120, -t2050, 0, 0, 0, t1905, 0, -t1995, 0, t1780, t1781, t1742, t1604, t1829, t1806, t1826, t1828, t1825, t1871, t1555, t1556, t1608, t1436, t1653, t1595, t1633, t1652, t1632, t1727, t1389, t1394, t1365, t1338, t2165, -t1457, t2250, t2184, -t1473, t2185, t1339, t1340, t1329, t1311, t2165, t2250, t1457, t2185, t1473, t2184, t1331, t1322, t1330, t1308; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2120, 0, -t1909, 0, 0, 0, t1890, 0, t1891, -t2065, t1726, t1725, t1737, t1625, t1774, t1741, t1770, t1773, t1769, t1850, t1539, t1540, t1536, t1435, t1577, t1518, t1562, t1576, t1561, t1629, t1376, t1385, t1356, t1333, t2159, -t1407, t2262, t2209, -t1432, t2208, t1325, t1327, t1321, t1310, t2159, t2262, t1407, t2208, t1432, t2209, t1319, t1315, t1317, t1307; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2050, t1909, 0, 0, 0, 0, t1888, 0, t1889, t1926, t1724, t1723, t1734, t1626, t1772, t1740, t1768, t1771, t1767, t1849, t1537, t1538, t1535, t1434, t1575, t1517, t1559, t1574, t1558, t1628, t1375, t1384, t1355, t1332, t2160, -t1405, t2263, t2212, -t1428, t2211, t1324, t1326, t1320, t1309, t2160, t2263, t1405, t2211, t1428, t2212, t1318, t1314, t1316, t1306; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), 0, -t2131, 0, 0, -t1872, t1822, 0, t1870, t1847, t1866, t1869, t1864, t1899, t1762, t1763, t1659, pkin(8) * t1659, t1691, t1621, t1684, t1690, t1683, t1731, t1499, t1519, t1402, t1383, t2133, -t1492, t2225, t2162, -t1513, t2164, t1362, t1363, t1353, t1328, t2133, t2225, t1492, t2164, t1513, t2162, t1350, t1345, t1349, t1313; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2131, 0, qJDD(2), 0, t1872, 0, t1823, 0, t1921, -t1908, -t2064, -t1921, -t2062, -qJDD(3), t1728, t1729, 0, pkin(2) * t1659, -t1843, -t1839, -t1804, t1843, -t1984, -qJDD(3), t1524, t1534, t1557, t1403, t1698, -t1614, -t2191, -t1985, -t1671, -t2033, t1404, t1410, t1377, t1347, t1698, -t2191, t1614, -t2033, t1671, -t1985, t1388, t1366, t1387, t1336; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t1822, -t1823, 0, 0, t1859, t1846, t1862, t1858, t1860, 0, t1739, t1738, t1638, t1630, t1689, t1619, t1681, t1688, t1680, t1730, t1488, t1506, t1397, t1380, t2132, t1491, t2224, t2161, t1510, t2163, t1357, t1360, t1352, t1323, t2132, t2224, -t1491, t2163, -t1510, t2161, t1348, t1344, t1346, t1312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1902, t1903, t1911, -t1924, t1919, t1924, 0, t1809, t1764, 0, t1793, t1720, t1786, t1791, t1785, t1817, t1647, t1667, t1489, -qJ(4) * t1525, t2035, t1581, t2210, t2136, t1603, t2135, t1409, t1411, t1392, t1364, t2035, t2210, -t1581, t2135, -t1603, t2136, t1379, t1371, t1378, t1341; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2054, t1901, t1917, t1983, t1912, -t2054, -t1809, 0, t1765, 0, t1792, t1718, t1783, t1790, t1782, t1816, t1624, t1631, t1480, t1494, t2036, t1579, t2213, t2137, t1599, t2139, t1395, t1396, t1386, t1354, t2036, t2213, -t1579, t2139, -t1599, t2137, t1374, t1367, t1373, t1337; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1921, t1908, t2064, t1921, t2062, qJDD(3), -t1764, -t1765, 0, 0, t1843, t1839, t1804, -t1843, t1984, qJDD(3), t1882 - t2048 + t2123, t2040 + 0.2e1 * t2113, t2124, t2125, -t1698, t1614, t2191, t1985, t1671, t2033, t1990, t1991, t1978, t2041, -t1698, t2191, -t1614, t2033, -t1671, t1985, t1955, t1979, t1958, t2145; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1845, -t1800, t2141, t1883, t1875, -t1883, 0, t1749, t1606, 0, t1699, t1618, t2190, t1986, t1675, t2034, t1531, t1533, t1425, -pkin(9) * t1439, t1699, t2190, -t1618, t2034, -t1675, t1986, t1445, t1393, t1444, t1372; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2081, t1803, t1876, -t2046, t1833, -t2081, -t1749, 0, t1607, 0, -t1799, -t2147, -t2143, t1799, t1707, -t1842, t1479, t1481, -pkin(4) * t1612, -pkin(4) * t1439, -t1799, -t2143, t2147, -t1842, -t1707, t1799, t1446, t1532, t1443, t1370; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1843, t1839, t1804, -t1843, t1984, qJDD(3), -t1606, -t1607, 0, 0, -t1698, t1614, t2191, t1985, t1671, t2033, t2051, t2052, t2042, t2068, -t1698, t2191, -t1614, t2033, -t1671, t1985, t1977, t2043, t1980, t1956; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1959, -t2150, t2149, t2086, t1811, -t2086, 0, t1572, t1522, 0, t1959, t2149, t2150, -t2086, -t1811, t2086, -qJ(6) * t2150, t1466, t1469, -qJ(6) * t1505; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1820, t2144, -t1812, -t1779, t2148, -t1820, -t1572, 0, t1523, 0, t1820, -t1812, -t2144, -t1820, -t2148, -t1779, t1470, t1462, pkin(5) * t2144, -pkin(5) * t1505; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1799, t2147, t2143, -t1799, -t1707, t1842, -t1522, -t1523, 0, 0, t1799, t2143, -t2147, t1842, t1707, -t1799, t1952, t2038, t1873 + t1954, t2039; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1959, t2149, t2150, -t2086, -t1811, t2086, 0, t1487, -t1505, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1799, t2143, -t2147, t1842, t1707, -t1799, -t1487, 0, t1482, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1820, t1812, t2144, t1820, t2148, t1779, t1505, -t1482, 0, 0;];
m_new_reg  = t1;
