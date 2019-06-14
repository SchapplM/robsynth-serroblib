% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6RRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-08 04:55
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6RRRRRP4_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP4_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP4_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP4_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP4_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP4_invdynm_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 04:53:52
% EndTime: 2019-05-08 04:54:58
% DurationCPUTime: 70.77s
% Computational Cost: add. (391137->951), mult. (780917->1223), div. (0->0), fcn. (579007->10), ass. (0->638)
t2067 = sin(qJ(3));
t2068 = sin(qJ(2));
t2072 = cos(qJ(3));
t2073 = cos(qJ(2));
t2032 = (t2067 * t2073 + t2068 * t2072) * qJD(1);
t2062 = qJD(2) + qJD(3);
t2066 = sin(qJ(4));
t2071 = cos(qJ(4));
t1997 = t2032 * t2066 - t2071 * t2062;
t1998 = t2032 * t2071 + t2062 * t2066;
t2065 = sin(qJ(5));
t2070 = cos(qJ(5));
t1952 = -t1997 * t2065 + t1998 * t2070;
t1949 = t1952 ^ 2;
t2199 = qJD(1) * t2073;
t2200 = qJD(1) * t2068;
t2030 = t2067 * t2200 - t2072 * t2199;
t2026 = qJD(4) + t2030;
t2024 = qJD(5) + t2026;
t2211 = t2024 ^ 2;
t1869 = t2211 + t1949;
t2197 = qJD(2) * t2068;
t2138 = qJD(1) * t2197;
t2151 = t2073 * qJDD(1);
t2102 = t2138 - t2151;
t2137 = qJD(2) * t2199;
t2152 = t2068 * qJDD(1);
t2103 = t2137 + t2152;
t2088 = t2067 * t2103 + t2072 * t2102;
t1971 = -t2032 * qJD(3) - t2088;
t1968 = qJDD(4) - t1971;
t1967 = qJDD(5) + t1968;
t1950 = t2070 * t1997 + t1998 * t2065;
t2175 = t1952 * t1950;
t2228 = t1967 + t2175;
t2182 = t2228 * t2065;
t1770 = t1869 * t2070 + t2182;
t2181 = t2228 * t2070;
t1789 = t1869 * t2065 - t2181;
t1698 = t1770 * t2066 + t1789 * t2071;
t2082 = t2067 * t2102 - t2072 * t2103;
t2079 = -t2030 * qJD(3) - t2082;
t2150 = qJDD(2) + qJDD(3);
t2077 = -t2066 * t2150 - t2071 * t2079;
t1910 = -t1997 * qJD(4) - t2077;
t2078 = -t2066 * t2079 + t2071 * t2150;
t2076 = t1998 * qJD(4) - t2078;
t1808 = -t1950 * qJD(5) + t2070 * t1910 - t2065 * t2076;
t2176 = t1950 * t2024;
t2231 = t1808 - t2176;
t1623 = t1698 * t2067 - t2072 * t2231;
t1625 = t1698 * t2072 + t2067 * t2231;
t1545 = t1623 * t2068 - t1625 * t2073;
t1680 = t1770 * t2071 - t1789 * t2066;
t2069 = sin(qJ(1));
t2074 = cos(qJ(1));
t2337 = pkin(6) * (t1545 * t2074 + t1680 * t2069);
t2336 = pkin(6) * (t1545 * t2069 - t1680 * t2074);
t1550 = t1623 * t2073 + t1625 * t2068;
t2335 = pkin(1) * t1550;
t2334 = pkin(7) * t1550;
t2333 = -pkin(1) * t1680 + pkin(7) * t1545;
t2133 = -t2065 * t1910 - t2070 * t2076;
t1807 = qJD(5) * t1952 - t2133;
t1922 = t2024 * t1952;
t2232 = t1807 + t1922;
t1669 = -t2232 * t2065 + t2070 * t2231;
t2185 = t2231 * t2065;
t1673 = -t2232 * t2070 - t2185;
t1599 = -t1669 * t2066 + t1673 * t2071;
t2213 = t1950 ^ 2;
t2230 = t1949 - t2213;
t1581 = t1599 * t2067 - t2072 * t2230;
t1583 = t1599 * t2072 + t2067 * t2230;
t1517 = t1581 * t2068 - t1583 * t2073;
t1595 = t1669 * t2071 + t1673 * t2066;
t2332 = t1517 * t2069 + t1595 * t2074;
t2331 = t1517 * t2074 - t1595 * t2069;
t1916 = t2213 - t2211;
t1799 = t1916 * t2065 + t2181;
t1803 = t1916 * t2070 - t2182;
t1711 = t1799 * t2066 - t1803 * t2071;
t1760 = t1807 - t1922;
t1629 = t1711 * t2067 - t1760 * t2072;
t1633 = t1711 * t2072 + t1760 * t2067;
t1561 = t1629 * t2068 - t1633 * t2073;
t1708 = t1799 * t2071 + t1803 * t2066;
t2328 = t1561 * t2069 - t1708 * t2074;
t2327 = t1561 * t2074 + t1708 * t2069;
t2326 = pkin(2) * t1623;
t2325 = pkin(8) * t1623;
t2324 = -pkin(2) * t1680 - pkin(8) * t1625;
t1516 = t1581 * t2073 + t1583 * t2068;
t1556 = t1629 * t2073 + t1633 * t2068;
t2321 = pkin(3) * t1680;
t2320 = pkin(9) * t1680;
t2313 = pkin(3) * t2231 - pkin(9) * t1698;
t1917 = t1949 - t2211;
t2225 = -t2175 + t1967;
t2180 = t2225 * t2065;
t2266 = -t2070 * t1917 + t2180;
t1834 = t2070 * t2225;
t2267 = t1917 * t2065 + t1834;
t2276 = t2066 * t2267 + t2071 * t2266;
t2224 = t2176 + t1808;
t2275 = -t2066 * t2266 + t2071 * t2267;
t2294 = t2067 * t2224 + t2072 * t2275;
t2295 = t2067 * t2275 - t2072 * t2224;
t2305 = -t2068 * t2295 + t2073 * t2294;
t2312 = t2069 * t2305 - t2074 * t2276;
t2311 = t2069 * t2276 + t2074 * t2305;
t2223 = -t2211 - t2213;
t2236 = t2070 * t2223 - t2180;
t2239 = t2065 * t2223 + t1834;
t2265 = -t2066 * t2239 + t2071 * t2236;
t2279 = t2067 * t2232 + t2072 * t2265;
t2282 = t2067 * t2265 - t2072 * t2232;
t2292 = t2068 * t2279 + t2073 * t2282;
t2310 = pkin(1) * t2292;
t2309 = pkin(7) * t2292;
t2264 = t2066 * t2236 + t2071 * t2239;
t2293 = -t2068 * t2282 + t2073 * t2279;
t2308 = -pkin(1) * t2264 + pkin(7) * t2293;
t2307 = pkin(6) * (t2069 * t2293 - t2074 * t2264);
t2306 = pkin(6) * (t2069 * t2264 + t2074 * t2293);
t2304 = t2068 * t2294 + t2073 * t2295;
t2303 = pkin(2) * t2282;
t2302 = pkin(4) * t1770;
t2301 = pkin(8) * t2282;
t2300 = pkin(10) * t1770;
t2299 = pkin(10) * t1789;
t2296 = -pkin(2) * t2264 + pkin(8) * t2279;
t2289 = pkin(3) * t2264;
t2288 = pkin(9) * t2264;
t2283 = -pkin(3) * t2232 + pkin(9) * t2265;
t2092 = (-t1950 * t2065 - t1952 * t2070) * t2024;
t2169 = t2024 * t2065;
t1913 = t1952 * t2169;
t2168 = t2024 * t2070;
t2143 = t1950 * t2168;
t2111 = t1913 - t2143;
t2220 = t2066 * t2111 + t2071 * t2092;
t2219 = -t2066 * t2092 + t2071 * t2111;
t2234 = t1967 * t2067 + t2072 * t2219;
t2238 = -t2072 * t1967 + t2067 * t2219;
t2263 = -t2068 * t2238 + t2073 * t2234;
t2281 = t2069 * t2263 - t2074 * t2220;
t2104 = t1807 * t2065 + t2143;
t2112 = -t2070 * t1807 + t1950 * t2169;
t2217 = t2066 * t2104 + t2071 * t2112;
t2145 = t2067 * t2175;
t2218 = -t2066 * t2112 + t2071 * t2104;
t2235 = t2072 * t2218 - t2145;
t2144 = t2072 * t2175;
t2237 = t2067 * t2218 + t2144;
t2261 = -t2068 * t2237 + t2073 * t2235;
t2280 = t2069 * t2261 - t2074 * t2217;
t2278 = t2069 * t2220 + t2074 * t2263;
t2277 = t2069 * t2217 + t2074 * t2261;
t1826 = -t2213 - t1949;
t2274 = pkin(3) * t1826;
t2273 = pkin(4) * t1826;
t2272 = pkin(4) * t2239;
t2271 = pkin(10) * t2236;
t2270 = pkin(10) * t2239;
t2269 = t1826 * t2067;
t2268 = t1826 * t2072;
t2262 = t2068 * t2234 + t2073 * t2238;
t2260 = t2068 * t2235 + t2073 * t2237;
t2257 = qJ(6) * t2231;
t1965 = t1998 * t1997;
t2229 = -t1965 + t1968;
t2256 = t2066 * t2229;
t1990 = t2032 * t2030;
t2227 = -t1990 + t2150;
t2254 = t2067 * t2227;
t2248 = t2071 * t2229;
t2246 = t2072 * t2227;
t2048 = t2069 * g(1) - t2074 * g(2);
t2064 = t2073 ^ 2;
t2148 = pkin(2) * t2073 + pkin(1);
t1973 = ((pkin(8) * t2064 + pkin(7)) * qJD(1) - pkin(2) * t2197) * qJD(1) + t2148 * qJDD(1) - (qJD(2) * pkin(2) - pkin(8) * t2200) * t2200 + t2048;
t1749 = t1808 * t2065 + t1952 * t2168;
t1750 = t1808 * t2070 - t1913;
t1662 = t1749 * t2071 + t1750 * t2066;
t1665 = -t1749 * t2066 + t1750 * t2071;
t2113 = t2072 * t1665 + t2145;
t2114 = t2067 * t1665 - t2144;
t2216 = -t2068 * t2114 + t2073 * t2113;
t2240 = -t2074 * t1662 + t2069 * t2216;
t2233 = t1662 * t2069 + t2074 * t2216;
t2157 = qJD(3) + t2062;
t2160 = t2062 * t2032;
t1819 = (t2030 * t2157 + t2082) * pkin(9) + (-t1971 + t2160) * pkin(3) - t1973;
t2049 = t2074 * g(1) + t2069 * g(2);
t2214 = qJD(1) ^ 2;
t2089 = -pkin(1) * t2214 + qJDD(1) * pkin(7) - t2049;
t2011 = -t2068 * g(3) + t2073 * t2089;
t2060 = t2064 * t2214;
t2075 = qJD(2) ^ 2;
t2053 = -t2060 - t2075;
t1963 = pkin(2) * t2053 + pkin(8) * t2151 + t2011;
t2203 = t2073 * g(3);
t2080 = qJDD(2) * pkin(2) - t2203 + ((-pkin(8) - pkin(7)) * qJDD(1) + t2148 * t2214 + t2049) * t2068;
t1896 = t2072 * t1963 + t2067 * t2080;
t1988 = pkin(3) * t2030 - pkin(9) * t2032;
t2210 = t2062 ^ 2;
t1833 = -pkin(3) * t2210 + pkin(9) * t2150 - t2030 * t1988 + t1896;
t1730 = -t2071 * t1819 + t2066 * t1833;
t1731 = t2066 * t1819 + t2071 * t1833;
t1649 = t2066 * t1730 + t2071 * t1731;
t1981 = t2026 * t1997;
t1885 = t1981 + t1910;
t2021 = t2062 * t2030;
t2226 = -t2021 + t2079;
t2215 = t2068 * t2113 + t2073 * t2114;
t2212 = t1997 ^ 2;
t1995 = t1998 ^ 2;
t2025 = t2026 ^ 2;
t2028 = t2030 ^ 2;
t2029 = t2032 ^ 2;
t1895 = t2067 * t1963 - t2072 * t2080;
t1809 = -t1895 * t2072 + t1896 * t2067;
t2209 = pkin(2) * t1809;
t1938 = t2021 + t2079;
t2085 = (-qJD(3) + t2062) * t2032 - t2088;
t1862 = -t1938 * t2072 + t2067 * t2085;
t2208 = pkin(2) * t1862;
t2207 = pkin(3) * t2067;
t1688 = pkin(4) * t2229 - pkin(10) * t1885 - t1730;
t1974 = pkin(4) * t2026 - pkin(10) * t1998;
t1702 = -pkin(4) * t2212 - pkin(10) * t2076 - t2026 * t1974 + t1731;
t1604 = -t2070 * t1688 + t1702 * t2065;
t1605 = t2065 * t1688 + t2070 * t1702;
t1530 = -t1604 * t2070 + t1605 * t2065;
t2206 = pkin(4) * t1530;
t1752 = t2070 * t2224;
t1762 = (-qJD(5) + t2024) * t1952 + t2133;
t1670 = t1762 * t2065 - t1752;
t2205 = pkin(4) * t1670;
t2204 = pkin(5) * t2070;
t2202 = qJ(6) * t2070;
t2201 = qJD(1) * qJD(2);
t2063 = t2068 ^ 2;
t2198 = t2214 * t2063;
t2196 = qJD(6) * t2024;
t2195 = t1530 * t2066;
t2194 = t1530 * t2071;
t1832 = -t2150 * pkin(3) - t2210 * pkin(9) + t2032 * t1988 + t1895;
t1733 = t2076 * pkin(4) - t2212 * pkin(10) + t1998 * t1974 + t1832;
t2192 = t1733 * t2065;
t2191 = t1733 * t2070;
t2186 = t2224 * t2065;
t2184 = t1809 * t2068;
t2183 = t1809 * t2073;
t1898 = t1965 + t1968;
t2179 = t1898 * t2066;
t2178 = t1898 * t2071;
t2173 = t1973 * t2067;
t2172 = t1973 * t2072;
t1985 = t1990 + t2150;
t2171 = t1985 * t2067;
t2170 = t1985 * t2072;
t2167 = t2026 * t2066;
t2166 = t2026 * t2071;
t2033 = qJDD(1) * pkin(1) + pkin(7) * t2214 + t2048;
t2165 = t2033 * t2068;
t2164 = t2033 * t2073;
t2039 = 0.2e1 * t2137 + t2152;
t2000 = t2039 * t2068;
t2040 = -0.2e1 * t2138 + t2151;
t1999 = t2040 * t2073;
t2055 = t2073 * t2214 * t2068;
t2046 = qJDD(2) + t2055;
t2163 = t2046 * t2068;
t2047 = qJDD(2) - t2055;
t2162 = t2047 * t2068;
t2161 = t2047 * t2073;
t2159 = t2062 * t2067;
t2158 = t2062 * t2072;
t1828 = t2066 * t1832;
t1829 = t2071 * t1832;
t2016 = 0.2e1 * t2196;
t1887 = pkin(5) * t1950 - qJ(6) * t1952;
t2110 = -pkin(5) * t2211 + t1967 * qJ(6) - t1950 * t1887 + t1605;
t1579 = t2016 + t2110;
t1585 = -t1967 * pkin(5) - qJ(6) * t2211 + t1887 * t1952 + qJDD(6) + t1604;
t2156 = -pkin(5) * t1585 + qJ(6) * t1579;
t2155 = -pkin(5) * t2224 - qJ(6) * t1760;
t2154 = -pkin(3) * t1832 + pkin(9) * t1649;
t2153 = t2063 + t2064;
t2147 = -pkin(3) * t2072 - pkin(2);
t2142 = t2067 * t1965;
t2141 = t2072 * t1965;
t2140 = t2069 * t1990;
t2139 = t2074 * t1990;
t1946 = -t1995 - t2025;
t1839 = -t1946 * t2066 - t2178;
t1886 = (qJD(4) + t2026) * t1997 + t2077;
t2136 = pkin(3) * t1886 + pkin(9) * t1839 + t1828;
t1932 = -t2025 - t2212;
t1825 = t1932 * t2071 - t2256;
t1982 = t2026 * t1998;
t1882 = -t1982 - t2076;
t2135 = pkin(3) * t1882 + pkin(9) * t1825 - t1829;
t2134 = -qJ(6) * t2065 - pkin(4);
t1531 = t1604 * t2065 + t2070 * t1605;
t1810 = t1895 * t2067 + t2072 * t1896;
t2010 = t2068 * t2089 + t2203;
t1962 = t2010 * t2068 + t2073 * t2011;
t2132 = -t2048 * t2069 - t2074 * t2049;
t2131 = t2069 * t2055;
t2130 = t2074 * t2055;
t1520 = t1579 * t2070 + t1585 * t2065;
t2083 = t1807 * pkin(5) + t1733 - t2257;
t1609 = (pkin(5) * t2024 - 0.2e1 * qJD(6)) * t1952 + t2083;
t1476 = pkin(10) * t1520 + (t2134 - t2204) * t1609;
t1519 = t1579 * t2065 - t1585 * t2070;
t1480 = -t1519 * t2066 + t1520 * t2071;
t1490 = -pkin(10) * t1519 + (pkin(5) * t2065 - t2202) * t1609;
t2129 = -pkin(3) * t1609 + pkin(9) * t1480 + t2071 * t1476 + t2066 * t1490;
t1563 = -pkin(5) * t1826 + t1579;
t1564 = -qJ(6) * t1826 + t1585;
t1672 = -t1760 * t2070 + t2186;
t1498 = pkin(10) * t1672 + t1563 * t2070 + t1564 * t2065 - t2273;
t1668 = -t1760 * t2065 - t1752;
t1500 = -pkin(10) * t1668 - t1563 * t2065 + t1564 * t2070;
t1598 = -t1668 * t2066 + t1672 * t2071;
t2128 = pkin(9) * t1598 + t2071 * t1498 + t2066 * t1500 - t2274;
t1674 = t1762 * t2070 + t2186;
t1508 = pkin(10) * t1674 + t1531 - t2273;
t1514 = -pkin(10) * t1670 - t1530;
t1600 = -t1670 * t2066 + t1674 * t2071;
t2127 = pkin(9) * t1600 + t2071 * t1508 + t2066 * t1514 - t2274;
t2081 = 0.2e1 * qJD(6) * t1952 - t2083;
t1588 = -pkin(5) * t1922 + t2081 + t2257;
t1528 = -t2299 + t2065 * t1588 + (pkin(4) + t2204) * t2231;
t1538 = -pkin(5) * t2185 + t1588 * t2070 - t2300;
t2126 = t2071 * t1528 + t2066 * t1538 + t2313;
t1589 = t2081 + (-t2232 - t1922) * pkin(5);
t1533 = t2070 * t1589 + t2134 * t2232 + t2271;
t1540 = -t1589 * t2065 - t2202 * t2232 - t2270;
t2125 = t2071 * t1533 + t2066 * t1540 + t2283;
t1614 = -pkin(4) * t2232 - t2191 + t2271;
t1654 = t2192 - t2270;
t2124 = t2071 * t1614 + t2066 * t1654 + t2283;
t1616 = -pkin(4) * t2231 + t2192 + t2299;
t1675 = t2191 + t2300;
t2123 = t2071 * t1616 + t2066 * t1675 - t2313;
t1883 = (-qJD(4) + t2026) * t1998 + t2078;
t1795 = t1883 * t2071 + t1885 * t2066;
t1914 = t1995 + t2212;
t2122 = pkin(3) * t1914 + pkin(9) * t1795 + t1649;
t1621 = t1649 * t2067 - t1832 * t2072;
t2121 = pkin(2) * t1621 + t2154;
t2009 = -t2029 - t2210;
t1939 = t2009 * t2072 - t2171;
t2120 = pkin(2) * t1939 - t1896;
t2119 = pkin(4) * t1519 + t2156;
t2118 = pkin(4) * t1668 + t2155;
t2117 = -t1605 - t2302;
t2043 = qJDD(1) * t2074 - t2069 * t2214;
t2115 = -pkin(6) * t2043 - g(3) * t2069;
t1648 = -t1730 * t2071 + t1731 * t2066;
t1961 = t2010 * t2073 - t2011 * t2068;
t2109 = t2048 * t2074 - t2049 * t2069;
t1768 = t1825 * t2067 + t1882 * t2072;
t2108 = pkin(2) * t1768 + t2135;
t1774 = t1839 * t2067 + t1886 * t2072;
t2107 = pkin(2) * t1774 + t2136;
t1984 = -t2210 - t2028;
t1919 = t1984 * t2067 + t2246;
t2106 = pkin(2) * t1919 - t1895;
t2105 = -t1604 + t2272;
t1494 = t1531 * t2071 - t2195;
t1523 = -pkin(4) * t1733 + pkin(10) * t1531;
t2101 = -pkin(3) * t1733 + pkin(9) * t1494 - pkin(10) * t2195 + t2071 * t1523;
t1473 = t1480 * t2067 - t1609 * t2072;
t2100 = pkin(2) * t1473 + t2129;
t1573 = t1598 * t2067 - t2268;
t2099 = pkin(2) * t1573 + t2128;
t1574 = t1600 * t2067 - t2268;
t2098 = pkin(2) * t1574 + t2127;
t2097 = t2126 - t2326;
t2096 = t2125 + t2303;
t2095 = t2124 + t2303;
t2094 = t2123 + t2326;
t1736 = t1795 * t2067 + t1914 * t2072;
t2093 = pkin(2) * t1736 + t2122;
t2091 = pkin(5) * t1869 + qJ(6) * t2228 + t2110;
t1487 = t1494 * t2067 - t1733 * t2072;
t2090 = pkin(2) * t1487 + t2101;
t2087 = t2091 + t2302;
t2086 = pkin(5) * t2225 + qJ(6) * t2223 - t1585;
t2084 = t2086 + t2272;
t2052 = t2060 - t2075;
t2051 = -t2075 - t2198;
t2050 = t2075 - t2198;
t2045 = -t2060 + t2198;
t2044 = t2060 + t2198;
t2042 = qJDD(1) * t2069 + t2074 * t2214;
t2041 = t2153 * qJDD(1);
t2037 = t2073 * t2046;
t2036 = t2153 * t2201;
t2027 = -pkin(6) * t2042 + g(3) * t2074;
t2015 = -t2029 + t2210;
t2014 = t2028 - t2210;
t2013 = -t2063 * t2201 + t2073 * t2103;
t2012 = -t2064 * t2201 + t2068 * t2102;
t2008 = -t2051 * t2068 - t2161;
t2007 = -t2050 * t2068 + t2037;
t2006 = t2053 * t2073 - t2163;
t2005 = t2052 * t2073 - t2162;
t2004 = t2051 * t2073 - t2162;
t2003 = t2050 * t2073 + t2163;
t2002 = t2053 * t2068 + t2037;
t2001 = t2052 * t2068 + t2161;
t1992 = t1999 - t2000;
t1991 = t2039 * t2073 + t2040 * t2068;
t1989 = t2029 - t2028;
t1980 = -t1995 + t2025;
t1979 = -t2025 + t2212;
t1978 = (-t2030 * t2072 + t2032 * t2067) * t2062;
t1977 = (-t2030 * t2067 - t2032 * t2072) * t2062;
t1976 = -pkin(7) * t2004 - t2164;
t1975 = -pkin(7) * t2002 - t2165;
t1972 = -t2028 - t2029;
t1970 = -pkin(1) * t2004 + t2011;
t1969 = -pkin(1) * t2002 + t2010;
t1959 = t1995 - t2212;
t1954 = pkin(1) * t2040 + pkin(7) * t2006 + t2164;
t1953 = -pkin(1) * t2039 + pkin(7) * t2008 - t2165;
t1944 = t2014 * t2072 - t2171;
t1943 = -t2015 * t2067 + t2246;
t1942 = t2014 * t2067 + t2170;
t1941 = t2015 * t2072 + t2254;
t1940 = -t2009 * t2067 - t2170;
t1933 = t2032 * t2157 + t2088;
t1931 = pkin(1) * t2033 + pkin(7) * t1962;
t1930 = -t2032 * t2159 + t2072 * t2079;
t1929 = t2032 * t2158 + t2067 * t2079;
t1928 = -t1971 * t2067 + t2030 * t2158;
t1927 = t1971 * t2072 + t2030 * t2159;
t1926 = pkin(1) * t2044 + pkin(7) * t2041 + t1962;
t1920 = t1984 * t2072 - t2254;
t1906 = (-t1997 * t2071 + t1998 * t2066) * t2026;
t1905 = (-t1997 * t2066 - t1998 * t2071) * t2026;
t1904 = -t1977 * t2068 + t1978 * t2073;
t1903 = t1977 * t2073 + t1978 * t2068;
t1890 = -pkin(8) * t1939 - t2172;
t1884 = -t1981 + t1910;
t1881 = -t1982 + t2076;
t1880 = -pkin(8) * t1919 - t2173;
t1877 = -t1942 * t2068 + t1944 * t2073;
t1876 = -t1941 * t2068 + t1943 * t2073;
t1875 = t1942 * t2073 + t1944 * t2068;
t1874 = t1941 * t2073 + t1943 * t2068;
t1873 = t1910 * t2071 - t1998 * t2167;
t1872 = t1910 * t2066 + t1998 * t2166;
t1871 = t1997 * t2166 + t2066 * t2076;
t1870 = -t1997 * t2167 + t2071 * t2076;
t1866 = -t1939 * t2068 + t1940 * t2073;
t1865 = t1939 * t2073 + t1940 * t2068;
t1864 = t1938 * t2067 + t2072 * t2085;
t1863 = -t1933 * t2072 - t2067 * t2226;
t1861 = -t1933 * t2067 + t2072 * t2226;
t1860 = t1906 * t2072 + t1968 * t2067;
t1859 = t1906 * t2067 - t1968 * t2072;
t1858 = t1979 * t2071 - t2179;
t1857 = -t1980 * t2066 + t2248;
t1856 = t1979 * t2066 + t2178;
t1855 = t1980 * t2071 + t2256;
t1854 = -t1929 * t2068 + t1930 * t2073;
t1853 = -t1927 * t2068 + t1928 * t2073;
t1852 = t1929 * t2073 + t1930 * t2068;
t1851 = t1927 * t2073 + t1928 * t2068;
t1846 = -t1919 * t2068 + t1920 * t2073;
t1845 = t1919 * t2073 + t1920 * t2068;
t1838 = t1946 * t2071 - t2179;
t1824 = t1932 * t2066 + t2248;
t1818 = -pkin(2) * t2226 + pkin(8) * t1940 - t2173;
t1815 = t1873 * t2072 + t2142;
t1814 = t1871 * t2072 - t2142;
t1813 = t1873 * t2067 - t2141;
t1812 = t1871 * t2067 + t2141;
t1811 = -pkin(2) * t1933 + pkin(8) * t1920 + t2172;
t1794 = t1882 * t2071 - t1884 * t2066;
t1793 = t1883 * t2066 - t1885 * t2071;
t1792 = t1882 * t2066 + t1884 * t2071;
t1786 = pkin(2) * t1973 + pkin(8) * t1810;
t1785 = -t1862 * t2068 + t1864 * t2073;
t1784 = -t1861 * t2068 + t1863 * t2073;
t1783 = t1862 * t2073 + t1864 * t2068;
t1782 = t1861 * t2073 + t1863 * t2068;
t1781 = t1858 * t2072 - t1881 * t2067;
t1780 = t1857 * t2072 + t1885 * t2067;
t1779 = t1858 * t2067 + t1881 * t2072;
t1778 = t1857 * t2067 - t1885 * t2072;
t1777 = -t1859 * t2068 + t1860 * t2073;
t1776 = t1859 * t2073 + t1860 * t2068;
t1775 = t1839 * t2072 - t1886 * t2067;
t1769 = t1825 * t2072 - t1882 * t2067;
t1751 = -pkin(1) * t1865 - t2120;
t1744 = t1794 * t2072 + t1959 * t2067;
t1743 = t1794 * t2067 - t1959 * t2072;
t1738 = -pkin(9) * t1838 + t1829;
t1737 = t1795 * t2072 - t1914 * t2067;
t1735 = -pkin(1) * t1845 - t2106;
t1734 = -pkin(9) * t1824 + t1828;
t1727 = -pkin(8) * t1862 - t1809;
t1722 = -t1813 * t2068 + t1815 * t2073;
t1721 = -t1812 * t2068 + t1814 * t2073;
t1720 = t1813 * t2073 + t1815 * t2068;
t1719 = t1812 * t2073 + t1814 * t2068;
t1718 = -pkin(2) * t1972 + pkin(8) * t1864 + t1810;
t1717 = -pkin(1) * t1783 - t2208;
t1716 = -pkin(7) * t1865 - t1818 * t2068 + t1890 * t2073;
t1715 = t1810 * t2073 - t2184;
t1714 = t1810 * t2068 + t2183;
t1713 = -pkin(7) * t1845 - t1811 * t2068 + t1880 * t2073;
t1704 = -pkin(1) * t2226 + pkin(7) * t1866 + t1818 * t2073 + t1890 * t2068;
t1703 = -pkin(3) * t1838 + t1731;
t1701 = -pkin(3) * t1824 + t1730;
t1693 = -pkin(1) * t1933 + pkin(7) * t1846 + t1811 * t2073 + t1880 * t2068;
t1692 = -t1779 * t2068 + t1781 * t2073;
t1691 = -t1778 * t2068 + t1780 * t2073;
t1690 = t1779 * t2073 + t1781 * t2068;
t1689 = t1778 * t2073 + t1780 * t2068;
t1685 = -t1774 * t2068 + t1775 * t2073;
t1684 = t1774 * t2073 + t1775 * t2068;
t1677 = -t1768 * t2068 + t1769 * t2073;
t1676 = t1768 * t2073 + t1769 * t2068;
t1656 = -t1743 * t2068 + t1744 * t2073;
t1655 = t1743 * t2073 + t1744 * t2068;
t1652 = -t1736 * t2068 + t1737 * t2073;
t1651 = t1736 * t2073 + t1737 * t2068;
t1650 = -pkin(1) * t1714 - t2209;
t1622 = t1649 * t2072 + t1832 * t2067;
t1612 = -pkin(9) * t1793 - t1648;
t1611 = -pkin(7) * t1714 - pkin(8) * t2183 - t1786 * t2068;
t1610 = pkin(1) * t1973 + pkin(7) * t1715 - pkin(8) * t2184 + t1786 * t2073;
t1607 = -pkin(7) * t1783 - t1718 * t2068 + t1727 * t2073;
t1606 = -pkin(1) * t1972 + pkin(7) * t1785 + t1718 * t2073 + t1727 * t2068;
t1602 = -pkin(8) * t1774 - t1703 * t2067 + t1738 * t2072;
t1601 = -pkin(8) * t1768 - t1701 * t2067 + t1734 * t2072;
t1596 = t1670 * t2071 + t1674 * t2066;
t1594 = t1668 * t2071 + t1672 * t2066;
t1590 = -pkin(2) * t1838 + pkin(8) * t1775 + t1703 * t2072 + t1738 * t2067;
t1587 = -pkin(1) * t1684 - t2107;
t1586 = -pkin(2) * t1824 + pkin(8) * t1769 + t1701 * t2072 + t1734 * t2067;
t1577 = -pkin(1) * t1676 - t2108;
t1576 = t1600 * t2072 + t2269;
t1575 = t1598 * t2072 + t2269;
t1562 = -pkin(8) * t1736 + t1612 * t2072 + t1793 * t2207;
t1549 = -t1621 * t2068 + t1622 * t2073;
t1548 = t1621 * t2073 + t1622 * t2068;
t1547 = pkin(8) * t1737 + t2067 * t1612 + t1793 * t2147;
t1542 = -pkin(3) * t1596 - t2205;
t1541 = -t2117 + t2321;
t1536 = -t2105 - t2289;
t1535 = -pkin(1) * t1651 - t2093;
t1534 = -t1616 * t2066 + t1675 * t2071 + t2320;
t1529 = -t1614 * t2066 + t1654 * t2071 - t2288;
t1526 = -pkin(8) * t1621 + (-pkin(9) * t2072 + t2207) * t1648;
t1525 = -t2084 - t2289;
t1524 = -pkin(3) * t1594 - t2118;
t1521 = -t2087 - 0.2e1 * t2196 - t2321;
t1512 = -t1574 * t2068 + t1576 * t2073;
t1511 = -t1573 * t2068 + t1575 * t2073;
t1510 = t1574 * t2073 + t1576 * t2068;
t1509 = t1573 * t2073 + t1575 * t2068;
t1506 = pkin(8) * t1622 + (-pkin(9) * t2067 + t2147) * t1648;
t1505 = -pkin(7) * t1684 - t1590 * t2068 + t1602 * t2073;
t1504 = -pkin(7) * t1676 - t1586 * t2068 + t1601 * t2073;
t1503 = -pkin(1) * t1838 + pkin(7) * t1685 + t1590 * t2073 + t1602 * t2068;
t1502 = -pkin(1) * t1824 + pkin(7) * t1677 + t1586 * t2073 + t1601 * t2068;
t1501 = -pkin(1) * t1548 - t2121;
t1496 = -pkin(7) * t1651 - t1547 * t2068 + t1562 * t2073;
t1495 = -pkin(1) * t1793 + pkin(7) * t1652 + t1547 * t2073 + t1562 * t2068;
t1493 = t1531 * t2066 + t2194;
t1491 = -t1533 * t2066 + t1540 * t2071 - t2288;
t1488 = t1494 * t2072 + t1733 * t2067;
t1486 = -t1528 * t2066 + t1538 * t2071 - t2320;
t1485 = t1534 * t2072 - t1541 * t2067 - t2325;
t1484 = -t2094 - t2335;
t1483 = t1529 * t2072 - t1536 * t2067 - t2301;
t1482 = -t2095 - t2310;
t1481 = t1534 * t2067 + t1541 * t2072 - t2324;
t1479 = t1519 * t2071 + t1520 * t2066;
t1477 = t1529 * t2067 + t1536 * t2072 + t2296;
t1474 = t1480 * t2072 + t1609 * t2067;
t1472 = -pkin(3) * t1493 - t2206;
t1471 = -pkin(9) * t1596 - t1508 * t2066 + t1514 * t2071;
t1470 = -pkin(7) * t1548 - t1506 * t2068 + t1526 * t2073;
t1469 = -t2096 - t2310;
t1468 = -pkin(1) * t1648 + pkin(7) * t1549 + t1506 * t2073 + t1526 * t2068;
t1467 = -t2097 + t2335;
t1466 = t1491 * t2072 - t1525 * t2067 - t2301;
t1465 = t1486 * t2072 - t1521 * t2067 + t2325;
t1464 = t1491 * t2067 + t1525 * t2072 + t2296;
t1463 = -pkin(9) * t1594 - t1498 * t2066 + t1500 * t2071;
t1462 = t1486 * t2067 + t1521 * t2072 + t2324;
t1461 = -pkin(9) * t1493 - pkin(10) * t2194 - t1523 * t2066;
t1460 = -t1487 * t2068 + t1488 * t2073;
t1459 = t1487 * t2073 + t1488 * t2068;
t1458 = -pkin(8) * t1574 + t1471 * t2072 - t1542 * t2067;
t1457 = -pkin(3) * t1479 - t2119;
t1456 = -pkin(1) * t1510 - t2098;
t1455 = -pkin(2) * t1596 + pkin(8) * t1576 + t1471 * t2067 + t1542 * t2072;
t1454 = -t1481 * t2068 + t1485 * t2073 - t2334;
t1453 = t1481 * t2073 + t1485 * t2068 - t2333;
t1452 = -pkin(8) * t1573 + t1463 * t2072 - t1524 * t2067;
t1451 = -t1477 * t2068 + t1483 * t2073 - t2309;
t1450 = -t1473 * t2068 + t1474 * t2073;
t1449 = t1473 * t2073 + t1474 * t2068;
t1448 = -pkin(1) * t1509 - t2099;
t1447 = t1477 * t2073 + t1483 * t2068 + t2308;
t1446 = -pkin(2) * t1594 + pkin(8) * t1575 + t1463 * t2067 + t1524 * t2072;
t1445 = -pkin(9) * t1479 - t1476 * t2066 + t1490 * t2071;
t1444 = -t1464 * t2068 + t1466 * t2073 - t2309;
t1443 = t1464 * t2073 + t1466 * t2068 + t2308;
t1442 = -t1462 * t2068 + t1465 * t2073 + t2334;
t1441 = t1462 * t2073 + t1465 * t2068 + t2333;
t1440 = -pkin(8) * t1487 + t1461 * t2072 - t1472 * t2067;
t1439 = -pkin(1) * t1459 - t2090;
t1438 = -pkin(7) * t1510 - t1455 * t2068 + t1458 * t2073;
t1437 = -pkin(2) * t1493 + pkin(8) * t1488 + t1461 * t2067 + t1472 * t2072;
t1436 = -pkin(1) * t1596 + pkin(7) * t1512 + t1455 * t2073 + t1458 * t2068;
t1435 = -pkin(7) * t1509 - t1446 * t2068 + t1452 * t2073;
t1434 = -pkin(1) * t1594 + pkin(7) * t1511 + t1446 * t2073 + t1452 * t2068;
t1433 = -pkin(8) * t1473 + t1445 * t2072 - t1457 * t2067;
t1432 = -pkin(1) * t1449 - t2100;
t1431 = -pkin(2) * t1479 + pkin(8) * t1474 + t1445 * t2067 + t1457 * t2072;
t1430 = -pkin(7) * t1459 - t1437 * t2068 + t1440 * t2073;
t1429 = -pkin(1) * t1493 + pkin(7) * t1460 + t1437 * t2073 + t1440 * t2068;
t1428 = -pkin(7) * t1449 - t1431 * t2068 + t1433 * t2073;
t1427 = -pkin(1) * t1479 + pkin(7) * t1450 + t1431 * t2073 + t1433 * t2068;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t2043, 0, -t2042, 0, t2115, -t2027, -t2109, -pkin(6) * t2109, t2013 * t2074 - t2131, t1992 * t2074 + t2045 * t2069, t2007 * t2074 + t2069 * t2152, t2012 * t2074 + t2131, t2005 * t2074 + t2069 * t2151, qJDD(2) * t2069 + t2036 * t2074, t2074 * t1975 - t2069 * t1969 - pkin(6) * (t2006 * t2069 + t2040 * t2074), t2074 * t1976 - t2069 * t1970 - pkin(6) * (t2008 * t2069 - t2039 * t2074), t2074 * t1961 - pkin(6) * (t2041 * t2069 + t2044 * t2074), -pkin(6) * (t1962 * t2069 + t2033 * t2074) - (pkin(1) * t2069 - pkin(7) * t2074) * t1961, t1854 * t2074 + t2140, t1784 * t2074 + t1989 * t2069, t1876 * t2074 + t1938 * t2069, t1853 * t2074 - t2140, t1877 * t2074 + t2069 * t2085, t2074 * t1904 + t2069 * t2150, t2074 * t1713 - t2069 * t1735 - pkin(6) * (t1846 * t2069 - t1933 * t2074), t2074 * t1716 - t2069 * t1751 - pkin(6) * (t1866 * t2069 - t2074 * t2226), t2074 * t1607 - t2069 * t1717 - pkin(6) * (t1785 * t2069 - t1972 * t2074), t2074 * t1611 - t2069 * t1650 - pkin(6) * (t1715 * t2069 + t1973 * t2074), t1722 * t2074 + t1872 * t2069, t1656 * t2074 + t1792 * t2069, t1691 * t2074 + t1855 * t2069, t1721 * t2074 - t1870 * t2069, t1692 * t2074 + t1856 * t2069, t1777 * t2074 + t1905 * t2069, t2074 * t1504 - t2069 * t1577 - pkin(6) * (t1677 * t2069 - t1824 * t2074), t2074 * t1505 - t2069 * t1587 - pkin(6) * (t1685 * t2069 - t1838 * t2074), t2074 * t1496 - t2069 * t1535 - pkin(6) * (t1652 * t2069 - t1793 * t2074), t2074 * t1470 - t2069 * t1501 - pkin(6) * (t1549 * t2069 - t1648 * t2074), t2233, -t2331, t2311, t2277, t2327, t2278, t2074 * t1451 - t2069 * t1482 - t2307, t2074 * t1454 - t2069 * t1484 + t2336, t2074 * t1438 - t2069 * t1456 - pkin(6) * (t1512 * t2069 - t1596 * t2074), t2074 * t1430 - t2069 * t1439 - pkin(6) * (t1460 * t2069 - t1493 * t2074), t2233, t2311, t2331, t2278, -t2327, t2277, t2074 * t1444 - t2069 * t1469 - t2307, t2074 * t1435 - t2069 * t1448 - pkin(6) * (t1511 * t2069 - t1594 * t2074), t2074 * t1442 - t2069 * t1467 - t2336, t2074 * t1428 - t2069 * t1432 - pkin(6) * (t1450 * t2069 - t1479 * t2074); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t2042, 0, t2043, 0, t2027, t2115, t2132, pkin(6) * t2132, t2013 * t2069 + t2130, t1992 * t2069 - t2045 * t2074, t2007 * t2069 - t2074 * t2152, t2012 * t2069 - t2130, t2005 * t2069 - t2074 * t2151, -qJDD(2) * t2074 + t2036 * t2069, t2069 * t1975 + t2074 * t1969 + pkin(6) * (t2006 * t2074 - t2040 * t2069), t2069 * t1976 + t2074 * t1970 + pkin(6) * (t2008 * t2074 + t2039 * t2069), t2069 * t1961 + pkin(6) * (t2041 * t2074 - t2044 * t2069), pkin(6) * (t1962 * t2074 - t2033 * t2069) - (-pkin(1) * t2074 - pkin(7) * t2069) * t1961, t1854 * t2069 - t2139, t1784 * t2069 - t1989 * t2074, t1876 * t2069 - t1938 * t2074, t1853 * t2069 + t2139, t1877 * t2069 - t2074 * t2085, t2069 * t1904 - t2074 * t2150, t2069 * t1713 + t2074 * t1735 + pkin(6) * (t1846 * t2074 + t1933 * t2069), t2069 * t1716 + t2074 * t1751 + pkin(6) * (t1866 * t2074 + t2069 * t2226), t2069 * t1607 + t2074 * t1717 + pkin(6) * (t1785 * t2074 + t1972 * t2069), t2069 * t1611 + t2074 * t1650 + pkin(6) * (t1715 * t2074 - t1973 * t2069), t1722 * t2069 - t1872 * t2074, t1656 * t2069 - t1792 * t2074, t1691 * t2069 - t1855 * t2074, t1721 * t2069 + t1870 * t2074, t1692 * t2069 - t1856 * t2074, t1777 * t2069 - t1905 * t2074, t2069 * t1504 + t2074 * t1577 + pkin(6) * (t1677 * t2074 + t1824 * t2069), t2069 * t1505 + t2074 * t1587 + pkin(6) * (t1685 * t2074 + t1838 * t2069), t2069 * t1496 + t2074 * t1535 + pkin(6) * (t1652 * t2074 + t1793 * t2069), t2069 * t1470 + t2074 * t1501 + pkin(6) * (t1549 * t2074 + t1648 * t2069), t2240, -t2332, t2312, t2280, t2328, t2281, t2069 * t1451 + t2074 * t1482 + t2306, t2069 * t1454 + t2074 * t1484 - t2337, t2069 * t1438 + t2074 * t1456 + pkin(6) * (t1512 * t2074 + t1596 * t2069), t2069 * t1430 + t2074 * t1439 + pkin(6) * (t1460 * t2074 + t1493 * t2069), t2240, t2312, t2332, t2281, -t2328, t2280, t2069 * t1444 + t2074 * t1469 + t2306, t2069 * t1435 + t2074 * t1448 + pkin(6) * (t1511 * t2074 + t1594 * t2069), t2069 * t1442 + t2074 * t1467 + t2337, t2069 * t1428 + t2074 * t1432 + pkin(6) * (t1450 * t2074 + t1479 * t2069); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t2048, t2049, 0, 0, t2000, t1991, t2003, t1999, t2001, 0, t1954, t1953, t1926, t1931, t1852, t1782, t1874, t1851, t1875, t1903, t1693, t1704, t1606, t1610, t1720, t1655, t1689, t1719, t1690, t1776, t1502, t1503, t1495, t1468, t2215, t1516, t2304, t2260, -t1556, t2262, t1447, t1453, t1436, t1429, t2215, t2304, -t1516, t2262, t1556, t2260, t1443, t1434, t1441, t1427; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t2214, 0, 0, -g(3), -t2048, 0, t2013, t1992, t2007, t2012, t2005, t2036, t1975, t1976, t1961, pkin(7) * t1961, t1854, t1784, t1876, t1853, t1877, t1904, t1713, t1716, t1607, t1611, t1722, t1656, t1691, t1721, t1692, t1777, t1504, t1505, t1496, t1470, t2216, -t1517, t2305, t2261, t1561, t2263, t1451, t1454, t1438, t1430, t2216, t2305, t1517, t2263, -t1561, t2261, t1444, t1435, t1442, t1428; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2214, 0, qJDD(1), 0, g(3), 0, -t2049, 0, t2055, -t2045, -t2152, -t2055, -t2151, -qJDD(2), t1969, t1970, 0, pkin(1) * t1961, -t1990, -t1989, -t1938, t1990, -t2085, -t2150, t1735, t1751, t1717, t1650, -t1872, -t1792, -t1855, t1870, -t1856, -t1905, t1577, t1587, t1535, t1501, -t1662, -t1595, -t2276, -t2217, -t1708, -t2220, t1482, t1484, t1456, t1439, -t1662, -t2276, t1595, -t2220, t1708, -t2217, t1469, t1448, t1467, t1432; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t2048, t2049, 0, 0, t2000, t1991, t2003, t1999, t2001, 0, t1954, t1953, t1926, t1931, t1852, t1782, t1874, t1851, t1875, t1903, t1693, t1704, t1606, t1610, t1720, t1655, t1689, t1719, t1690, t1776, t1502, t1503, t1495, t1468, t2215, t1516, t2304, t2260, -t1556, t2262, t1447, t1453, t1436, t1429, t2215, t2304, -t1516, t2262, t1556, t2260, t1443, t1434, t1441, t1427; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2103, t2040, t2046, -t2137, t2052, t2137, 0, -t2033, t2010, 0, t1930, t1863, t1943, t1928, t1944, t1978, t1880, t1890, t1727, -pkin(8) * t1809, t1815, t1744, t1780, t1814, t1781, t1860, t1601, t1602, t1562, t1526, t2113, t1583, t2294, t2235, -t1633, t2234, t1483, t1485, t1458, t1440, t2113, t2294, -t1583, t2234, t1633, t2235, t1466, t1452, t1465, t1433; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2138, t2039, t2050, -t2102, t2047, -t2138, t2033, 0, t2011, 0, t1929, t1861, t1941, t1927, t1942, t1977, t1811, t1818, t1718, t1786, t1813, t1743, t1778, t1812, t1779, t1859, t1586, t1590, t1547, t1506, t2114, t1581, t2295, t2237, -t1629, t2238, t1477, t1481, t1455, t1437, t2114, t2295, -t1581, t2238, t1629, t2237, t1464, t1446, t1462, t1431; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2055, t2045, t2152, t2055, t2151, qJDD(2), -t2010, -t2011, 0, 0, t1990, t1989, t1938, -t1990, t2085, t2150, t2106, t2120, t2208, t2209, t1872, t1792, t1855, -t1870, t1856, t1905, t2108, t2107, t2093, t2121, t1662, t1595, t2276, t2217, t1708, t2220, t2095, t2094, t2098, t2090, t1662, t2276, -t1595, t2220, -t1708, t2217, t2096, t2099, t2097, t2100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2079, -t1933, t2227, t2021, t2014, -t2021, 0, -t1973, t1895, 0, t1873, t1794, t1857, t1871, t1858, t1906, t1734, t1738, t1612, -pkin(9) * t1648, t1665, t1599, t2275, t2218, -t1711, t2219, t1529, t1534, t1471, t1461, t1665, t2275, -t1599, t2219, t1711, t2218, t1491, t1463, t1486, t1445; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2160, t2226, t2015, t1971, t1985, -t2160, t1973, 0, t1896, 0, -t1965, -t1959, -t1885, t1965, t1881, -t1968, t1701, t1703, -pkin(3) * t1793, -pkin(3) * t1648, -t2175, -t2230, -t2224, t2175, t1760, -t1967, t1536, t1541, t1542, t1472, -t2175, -t2224, t2230, -t1967, -t1760, t2175, t1525, t1524, t1521, t1457; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1990, t1989, t1938, -t1990, t2085, t2150, -t1895, -t1896, 0, 0, t1872, t1792, t1855, -t1870, t1856, t1905, t2135, t2136, t2122, t2154, t1662, t1595, t2276, t2217, t1708, t2220, t2124, t2123, t2127, t2101, t1662, t2276, -t1595, t2220, -t1708, t2217, t2125, t2128, t2126, t2129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1910, t1882, t2229, t1981, t1979, -t1981, 0, t1832, t1730, 0, t1750, t1673, t2267, t2104, t1803, t2111, t1654, t1675, t1514, -pkin(10) * t1530, t1750, t2267, -t1673, t2111, -t1803, t2104, t1540, t1500, t1538, t1490; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1982, t1884, t1980, -t2076, t1898, -t1982, -t1832, 0, t1731, 0, t1749, t1669, t2266, t2112, t1799, t2092, t1614, t1616, t1508, t1523, t1749, t2266, -t1669, t2092, -t1799, t2112, t1533, t1498, t1528, t1476; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1965, t1959, t1885, -t1965, -t1881, t1968, -t1730, -t1731, 0, 0, t2175, t2230, t2224, -t2175, -t1760, t1967, t2105, t2117, t2205, t2206, t2175, t2224, -t2230, t1967, t1760, -t2175, t2084, t2118, t2016 + t2087, t2119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1808, -t2232, t2225, t2176, t1916, -t2176, 0, t1733, t1604, 0, t1808, t2225, t2232, -t2176, -t1916, t2176, -qJ(6) * t2232, t1564, t1588, -qJ(6) * t1609; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1922, t2231, -t1917, -t1807, t2228, -t1922, -t1733, 0, t1605, 0, t1922, -t1917, -t2231, -t1922, -t2228, -t1807, t1589, t1563, pkin(5) * t2231, -pkin(5) * t1609; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2175, t2230, t2224, -t2175, -t1760, t1967, -t1604, -t1605, 0, 0, t2175, t2224, -t2230, t1967, t1760, -t2175, t2086, t2155, t2016 + t2091, t2156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1808, t2225, t2232, -t2176, -t1916, t2176, 0, t1585, -t1609, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2175, t2224, -t2230, t1967, t1760, -t2175, -t1585, 0, t1579, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1922, t1917, t2231, t1922, t2228, t1807, t1609, -t1579, 0, 0;];
m_new_reg  = t1;