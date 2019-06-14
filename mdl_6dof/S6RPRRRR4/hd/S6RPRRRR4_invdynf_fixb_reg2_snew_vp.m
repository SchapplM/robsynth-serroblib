% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RPRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 03:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RPRRRR4_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR4_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR4_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_invdynf_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 03:17:26
% EndTime: 2019-05-06 03:17:37
% DurationCPUTime: 10.73s
% Computational Cost: add. (108098->330), mult. (270856->470), div. (0->0), fcn. (221754->12), ass. (0->260)
t2208 = cos(pkin(11));
t2220 = qJD(1) ^ 2;
t2213 = sin(qJ(1));
t2218 = cos(qJ(1));
t2190 = t2213 * g(1) - t2218 * g(2);
t2230 = -qJDD(2) + t2190;
t2207 = sin(pkin(11));
t2204 = t2207 ^ 2;
t2205 = t2208 ^ 2;
t2242 = t2204 + t2205;
t2164 = t2220 * (pkin(7) * t2242 + qJ(2)) + (pkin(2) * t2208 + pkin(1)) * qJDD(1) + t2230;
t2201 = t2208 * qJDD(1);
t2212 = sin(qJ(3));
t2217 = cos(qJ(3));
t2241 = t2207 * qJDD(1);
t2147 = t2217 * t2201 - t2212 * t2241;
t2227 = t2207 * t2217 + t2208 * t2212;
t2182 = t2227 * qJD(1);
t2249 = t2182 * qJD(3);
t2168 = t2147 - t2249;
t2255 = qJD(1) * t2208;
t2256 = qJD(1) * t2207;
t2180 = t2212 * t2256 - t2217 * t2255;
t2178 = t2180 ^ 2;
t2232 = qJD(3) * pkin(3) - pkin(8) * t2182;
t2119 = t2168 * pkin(3) + t2178 * pkin(8) - t2182 * t2232 + t2164;
t2211 = sin(qJ(4));
t2216 = cos(qJ(4));
t2160 = -t2180 * t2211 + t2182 * t2216;
t2251 = t2180 * qJD(3);
t2274 = t2227 * qJDD(1);
t2170 = t2274 - t2251;
t2233 = -t2216 * t2168 + t2211 * t2170;
t2122 = -qJD(4) * t2160 - t2233;
t2158 = t2180 * t2216 + t2182 * t2211;
t2157 = t2158 ^ 2;
t2206 = qJD(3) + qJD(4);
t2231 = pkin(4) * t2206 - pkin(9) * t2160;
t2076 = t2122 * pkin(4) + t2157 * pkin(9) - t2160 * t2231 + t2119;
t2187 = t2242 * t2220;
t2210 = sin(qJ(5));
t2215 = cos(qJ(5));
t2133 = t2215 * t2158 + t2160 * t2210;
t2131 = qJD(6) + t2133;
t2275 = qJD(6) + t2131;
t2135 = -t2158 * t2210 + t2160 * t2215;
t2202 = qJD(5) + t2206;
t2209 = sin(qJ(6));
t2214 = cos(qJ(6));
t2124 = t2135 * t2209 - t2214 * t2202;
t2269 = t2124 ^ 2;
t2126 = t2135 * t2214 + t2202 * t2209;
t2268 = t2126 ^ 2;
t2267 = t2131 ^ 2;
t2266 = t2133 ^ 2;
t2265 = t2135 ^ 2;
t2264 = t2160 ^ 2;
t2263 = t2182 ^ 2;
t2262 = t2202 ^ 2;
t2261 = t2206 ^ 2;
t2254 = t2124 * t2126;
t2253 = t2133 * t2135;
t2252 = t2158 * t2160;
t2250 = t2180 * t2182;
t2248 = t2205 * t2220;
t2247 = t2206 * t2158;
t2246 = t2208 * t2220;
t2245 = qJD(4) - t2206;
t2244 = qJD(5) - t2202;
t2243 = qJD(6) - t2131;
t2191 = -g(1) * t2218 - g(2) * t2213;
t2183 = -pkin(1) * t2220 + qJDD(1) * qJ(2) + t2191;
t2238 = -t2208 * g(3) - 0.2e1 * qJD(2) * t2256;
t2153 = (pkin(2) * t2246 - pkin(7) * qJDD(1) - t2183) * t2207 + t2238;
t2172 = -g(3) * t2207 + 0.2e1 * qJD(2) * t2255 + t2208 * t2183;
t2155 = -pkin(2) * t2248 + pkin(7) * t2201 + t2172;
t2128 = t2217 * t2153 - t2212 * t2155;
t2165 = qJDD(3) - t2250;
t2097 = (-t2170 - t2251) * pkin(8) + t2165 * pkin(3) + t2128;
t2129 = t2212 * t2153 + t2217 * t2155;
t2106 = -t2178 * pkin(3) + t2168 * pkin(8) - qJD(3) * t2232 + t2129;
t2078 = t2211 * t2097 + t2216 * t2106;
t2063 = -t2157 * pkin(4) + t2122 * pkin(9) - t2206 * t2231 + t2078;
t2077 = t2216 * t2097 - t2211 * t2106;
t2228 = -t2211 * t2168 - t2216 * t2170;
t2123 = -qJD(4) * t2158 - t2228;
t2240 = qJDD(3) + qJDD(4);
t2137 = t2240 - t2252;
t2223 = (-t2123 - t2247) * pkin(9) + t2137 * pkin(4) + t2077;
t2039 = t2215 * t2063 + t2210 * t2223;
t2237 = qJDD(5) + t2240;
t2038 = -t2063 * t2210 + t2215 * t2223;
t2229 = -t2210 * t2122 - t2215 * t2123;
t2083 = -qJD(5) * t2133 - t2229;
t2236 = t2202 * t2133 - t2083;
t2235 = -t2209 * t2083 + t2214 * t2237;
t2234 = -t2215 * t2122 + t2210 * t2123;
t2225 = -qJD(5) * t2135 - qJDD(6) - t2234;
t2070 = (qJD(5) + t2202) * t2135 + t2234;
t2224 = -t2214 * t2083 - t2209 * t2237;
t2219 = qJD(3) ^ 2;
t2194 = t2207 * t2246;
t2189 = -qJDD(1) * t2213 - t2218 * t2220;
t2188 = qJDD(1) * t2218 - t2213 * t2220;
t2186 = t2242 * qJDD(1);
t2185 = t2208 * t2187;
t2184 = t2207 * t2187;
t2177 = qJDD(1) * pkin(1) + t2220 * qJ(2) + t2230;
t2173 = -t2219 - t2263;
t2171 = -t2207 * t2183 + t2238;
t2169 = t2274 - 0.2e1 * t2251;
t2167 = -t2147 + 0.2e1 * t2249;
t2166 = -qJDD(3) - t2250;
t2163 = -t2219 - t2178;
t2151 = -t2261 - t2264;
t2148 = -t2178 - t2263;
t2146 = t2166 * t2217 - t2173 * t2212;
t2145 = t2166 * t2212 + t2173 * t2217;
t2144 = -t2171 * t2207 + t2172 * t2208;
t2143 = t2171 * t2208 + t2172 * t2207;
t2142 = t2147 * t2217 + t2212 * t2274;
t2141 = t2147 * t2212 - t2217 * t2274;
t2140 = t2163 * t2217 - t2165 * t2212;
t2139 = t2163 * t2212 + t2165 * t2217;
t2138 = -t2240 - t2252;
t2136 = -t2261 - t2157;
t2127 = -t2262 - t2265;
t2121 = -t2157 - t2264;
t2118 = -t2145 * t2207 + t2146 * t2208;
t2117 = t2145 * t2208 + t2146 * t2207;
t2116 = t2138 * t2216 - t2151 * t2211;
t2115 = t2138 * t2211 + t2151 * t2216;
t2114 = -t2141 * t2207 + t2142 * t2208;
t2113 = t2141 * t2208 + t2142 * t2207;
t2112 = t2158 * t2245 + t2228;
t2111 = t2123 - t2247;
t2110 = -t2160 * t2245 - t2233;
t2109 = (qJD(4) + t2206) * t2160 + t2233;
t2108 = -t2139 * t2207 + t2140 * t2208;
t2107 = t2139 * t2208 + t2140 * t2207;
t2105 = t2136 * t2216 - t2137 * t2211;
t2104 = t2136 * t2211 + t2137 * t2216;
t2103 = pkin(5) * t2133 - pkin(10) * t2135;
t2101 = -t2237 - t2253;
t2100 = t2237 - t2253;
t2099 = -t2262 - t2266;
t2095 = -t2128 * t2212 + t2129 * t2217;
t2094 = t2128 * t2217 + t2129 * t2212;
t2093 = -t2265 - t2266;
t2092 = -t2267 - t2268;
t2091 = t2101 * t2215 - t2127 * t2210;
t2090 = t2101 * t2210 + t2127 * t2215;
t2089 = -t2115 * t2212 + t2116 * t2217;
t2088 = t2115 * t2217 + t2116 * t2212;
t2087 = -t2267 - t2269;
t2086 = -t2268 - t2269;
t2085 = t2110 * t2216 - t2112 * t2211;
t2084 = t2110 * t2211 + t2112 * t2216;
t2082 = -t2104 * t2212 + t2105 * t2217;
t2081 = t2104 * t2217 + t2105 * t2212;
t2080 = t2099 * t2215 - t2100 * t2210;
t2079 = t2099 * t2210 + t2100 * t2215;
t2075 = -t2094 * t2207 + t2095 * t2208;
t2074 = t2094 * t2208 + t2095 * t2207;
t2073 = t2133 * t2244 + t2229;
t2071 = -t2135 * t2244 - t2234;
t2069 = t2225 - t2254;
t2068 = -t2225 - t2254;
t2067 = -t2090 * t2211 + t2091 * t2216;
t2066 = t2090 * t2216 + t2091 * t2211;
t2065 = -t2088 * t2207 + t2089 * t2208;
t2064 = t2088 * t2208 + t2089 * t2207;
t2059 = t2124 * t2243 + t2224;
t2058 = -t2124 * t2275 - t2224;
t2057 = -t2126 * t2243 + t2235;
t2056 = t2126 * t2275 - t2235;
t2055 = -t2084 * t2212 + t2085 * t2217;
t2054 = t2084 * t2217 + t2085 * t2212;
t2053 = -t2081 * t2207 + t2082 * t2208;
t2052 = t2081 * t2208 + t2082 * t2207;
t2051 = -t2079 * t2211 + t2080 * t2216;
t2050 = t2079 * t2216 + t2080 * t2211;
t2049 = -t2077 * t2211 + t2078 * t2216;
t2048 = t2077 * t2216 + t2078 * t2211;
t2047 = t2069 * t2214 - t2092 * t2209;
t2046 = t2069 * t2209 + t2092 * t2214;
t2045 = -t2068 * t2209 + t2087 * t2214;
t2044 = t2068 * t2214 + t2087 * t2209;
t2043 = t2071 * t2215 - t2073 * t2210;
t2042 = t2071 * t2210 + t2073 * t2215;
t2041 = -t2066 * t2212 + t2067 * t2217;
t2040 = t2066 * t2217 + t2067 * t2212;
t2037 = t2070 * pkin(5) + t2236 * pkin(10) - t2076;
t2036 = t2057 * t2214 - t2059 * t2209;
t2035 = t2057 * t2209 + t2059 * t2214;
t2034 = -t2054 * t2207 + t2055 * t2208;
t2033 = t2054 * t2208 + t2055 * t2207;
t2032 = t2047 * t2215 + t2058 * t2210;
t2031 = t2047 * t2210 - t2058 * t2215;
t2030 = t2045 * t2215 + t2056 * t2210;
t2029 = t2045 * t2210 - t2056 * t2215;
t2028 = -pkin(5) * t2262 + pkin(10) * t2237 - t2133 * t2103 + t2039;
t2027 = -pkin(5) * t2237 - pkin(10) * t2262 + t2103 * t2135 - t2038;
t2026 = -t2050 * t2212 + t2051 * t2217;
t2025 = t2050 * t2217 + t2051 * t2212;
t2024 = -t2048 * t2212 + t2049 * t2217;
t2023 = t2048 * t2217 + t2049 * t2212;
t2022 = t2036 * t2215 + t2086 * t2210;
t2021 = t2036 * t2210 - t2086 * t2215;
t2020 = -t2042 * t2211 + t2043 * t2216;
t2019 = t2042 * t2216 + t2043 * t2211;
t2018 = -t2040 * t2207 + t2041 * t2208;
t2017 = t2040 * t2208 + t2041 * t2207;
t2016 = -t2038 * t2210 + t2039 * t2215;
t2015 = t2038 * t2215 + t2039 * t2210;
t2014 = t2028 * t2214 + t2037 * t2209;
t2013 = -t2028 * t2209 + t2037 * t2214;
t2012 = -t2031 * t2211 + t2032 * t2216;
t2011 = t2031 * t2216 + t2032 * t2211;
t2010 = -t2029 * t2211 + t2030 * t2216;
t2009 = t2029 * t2216 + t2030 * t2211;
t2008 = -t2025 * t2207 + t2026 * t2208;
t2007 = t2025 * t2208 + t2026 * t2207;
t2006 = -t2023 * t2207 + t2024 * t2208;
t2005 = t2023 * t2208 + t2024 * t2207;
t2004 = -t2021 * t2211 + t2022 * t2216;
t2003 = t2021 * t2216 + t2022 * t2211;
t2002 = -t2019 * t2212 + t2020 * t2217;
t2001 = t2019 * t2217 + t2020 * t2212;
t2000 = -t2015 * t2211 + t2016 * t2216;
t1999 = t2015 * t2216 + t2016 * t2211;
t1998 = -t2013 * t2209 + t2014 * t2214;
t1997 = t2013 * t2214 + t2014 * t2209;
t1996 = -t2011 * t2212 + t2012 * t2217;
t1995 = t2011 * t2217 + t2012 * t2212;
t1994 = -t2009 * t2212 + t2010 * t2217;
t1993 = t2009 * t2217 + t2010 * t2212;
t1992 = -t2003 * t2212 + t2004 * t2217;
t1991 = t2003 * t2217 + t2004 * t2212;
t1990 = -t2001 * t2207 + t2002 * t2208;
t1989 = t2001 * t2208 + t2002 * t2207;
t1988 = t1998 * t2215 + t2027 * t2210;
t1987 = t1998 * t2210 - t2027 * t2215;
t1986 = -t1999 * t2212 + t2000 * t2217;
t1985 = t1999 * t2217 + t2000 * t2212;
t1984 = -t1995 * t2207 + t1996 * t2208;
t1983 = t1995 * t2208 + t1996 * t2207;
t1982 = -t1993 * t2207 + t1994 * t2208;
t1981 = t1993 * t2208 + t1994 * t2207;
t1980 = -t1991 * t2207 + t1992 * t2208;
t1979 = t1991 * t2208 + t1992 * t2207;
t1978 = -t1987 * t2211 + t1988 * t2216;
t1977 = t1987 * t2216 + t1988 * t2211;
t1976 = -t1985 * t2207 + t1986 * t2208;
t1975 = t1985 * t2208 + t1986 * t2207;
t1974 = -t1977 * t2212 + t1978 * t2217;
t1973 = t1977 * t2217 + t1978 * t2212;
t1972 = -t1973 * t2207 + t1974 * t2208;
t1971 = t1973 * t2208 + t1974 * t2207;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2189, -t2188, 0, -t2190 * t2213 + t2191 * t2218, 0, 0, 0, 0, 0, 0, -t2185 * t2218 - t2201 * t2213, t2184 * t2218 + t2213 * t2241, t2186 * t2218 - t2187 * t2213, t2144 * t2218 - t2177 * t2213, 0, 0, 0, 0, 0, 0, t2108 * t2218 + t2167 * t2213, t2118 * t2218 + t2169 * t2213, t2114 * t2218 + t2148 * t2213, t2075 * t2218 - t2164 * t2213, 0, 0, 0, 0, 0, 0, t2053 * t2218 + t2109 * t2213, t2065 * t2218 + t2111 * t2213, t2034 * t2218 + t2121 * t2213, t2006 * t2218 - t2119 * t2213, 0, 0, 0, 0, 0, 0, t2008 * t2218 + t2070 * t2213, t2018 * t2218 - t2213 * t2236, t1990 * t2218 + t2093 * t2213, t1976 * t2218 - t2076 * t2213, 0, 0, 0, 0, 0, 0, t1982 * t2218 + t2044 * t2213, t1984 * t2218 + t2046 * t2213, t1980 * t2218 + t2035 * t2213, t1972 * t2218 + t1997 * t2213; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2188, t2189, 0, t2190 * t2218 + t2191 * t2213, 0, 0, 0, 0, 0, 0, -t2185 * t2213 + t2201 * t2218, t2184 * t2213 - t2218 * t2241, t2186 * t2213 + t2187 * t2218, t2144 * t2213 + t2177 * t2218, 0, 0, 0, 0, 0, 0, t2108 * t2213 - t2167 * t2218, t2118 * t2213 - t2169 * t2218, t2114 * t2213 - t2148 * t2218, t2075 * t2213 + t2164 * t2218, 0, 0, 0, 0, 0, 0, t2053 * t2213 - t2109 * t2218, t2065 * t2213 - t2111 * t2218, t2034 * t2213 - t2121 * t2218, t2006 * t2213 + t2119 * t2218, 0, 0, 0, 0, 0, 0, t2008 * t2213 - t2070 * t2218, t2018 * t2213 + t2218 * t2236, t1990 * t2213 - t2093 * t2218, t1976 * t2213 + t2076 * t2218, 0, 0, 0, 0, 0, 0, t1982 * t2213 - t2044 * t2218, t1984 * t2213 - t2046 * t2218, t1980 * t2213 - t2035 * t2218, t1972 * t2213 - t1997 * t2218; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t2143, 0, 0, 0, 0, 0, 0, t2107, t2117, t2113, t2074, 0, 0, 0, 0, 0, 0, t2052, t2064, t2033, t2005, 0, 0, 0, 0, 0, 0, t2007, t2017, t1989, t1975, 0, 0, 0, 0, 0, 0, t1981, t1983, t1979, t1971; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2220, -qJDD(1), 0, t2191, 0, 0, 0, 0, 0, 0, -t2185, t2184, t2186, t2144, 0, 0, 0, 0, 0, 0, t2108, t2118, t2114, t2075, 0, 0, 0, 0, 0, 0, t2053, t2065, t2034, t2006, 0, 0, 0, 0, 0, 0, t2008, t2018, t1990, t1976, 0, 0, 0, 0, 0, 0, t1982, t1984, t1980, t1972; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2220, 0, t2190, 0, 0, 0, 0, 0, 0, t2201, -t2241, t2187, t2177, 0, 0, 0, 0, 0, 0, -t2167, -t2169, -t2148, t2164, 0, 0, 0, 0, 0, 0, -t2109, -t2111, -t2121, t2119, 0, 0, 0, 0, 0, 0, -t2070, t2236, -t2093, t2076, 0, 0, 0, 0, 0, 0, -t2044, -t2046, -t2035, -t1997; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t2143, 0, 0, 0, 0, 0, 0, t2107, t2117, t2113, t2074, 0, 0, 0, 0, 0, 0, t2052, t2064, t2033, t2005, 0, 0, 0, 0, 0, 0, t2007, t2017, t1989, t1975, 0, 0, 0, 0, 0, 0, t1981, t1983, t1979, t1971; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2248, t2194, t2201, t2172, 0, 0, 0, 0, 0, 0, t2140, t2146, t2142, t2095, 0, 0, 0, 0, 0, 0, t2082, t2089, t2055, t2024, 0, 0, 0, 0, 0, 0, t2026, t2041, t2002, t1986, 0, 0, 0, 0, 0, 0, t1994, t1996, t1992, t1974; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2194, -t2204 * t2220, -t2241, t2171, 0, 0, 0, 0, 0, 0, t2139, t2145, t2141, t2094, 0, 0, 0, 0, 0, 0, t2081, t2088, t2054, t2023, 0, 0, 0, 0, 0, 0, t2025, t2040, t2001, t1985, 0, 0, 0, 0, 0, 0, t1993, t1995, t1991, t1973; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2201, t2241, -t2187, -t2177, 0, 0, 0, 0, 0, 0, t2167, t2169, t2148, -t2164, 0, 0, 0, 0, 0, 0, t2109, t2111, t2121, -t2119, 0, 0, 0, 0, 0, 0, t2070, -t2236, t2093, -t2076, 0, 0, 0, 0, 0, 0, t2044, t2046, t2035, t1997; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2163, t2166, t2147, t2129, 0, 0, 0, 0, 0, 0, t2105, t2116, t2085, t2049, 0, 0, 0, 0, 0, 0, t2051, t2067, t2020, t2000, 0, 0, 0, 0, 0, 0, t2010, t2012, t2004, t1978; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2165, t2173, -t2274, t2128, 0, 0, 0, 0, 0, 0, t2104, t2115, t2084, t2048, 0, 0, 0, 0, 0, 0, t2050, t2066, t2019, t1999, 0, 0, 0, 0, 0, 0, t2009, t2011, t2003, t1977; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2167, t2169, t2148, -t2164, 0, 0, 0, 0, 0, 0, t2109, t2111, t2121, -t2119, 0, 0, 0, 0, 0, 0, t2070, -t2236, t2093, -t2076, 0, 0, 0, 0, 0, 0, t2044, t2046, t2035, t1997; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2136, t2138, t2110, t2078, 0, 0, 0, 0, 0, 0, t2080, t2091, t2043, t2016, 0, 0, 0, 0, 0, 0, t2030, t2032, t2022, t1988; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2137, t2151, t2112, t2077, 0, 0, 0, 0, 0, 0, t2079, t2090, t2042, t2015, 0, 0, 0, 0, 0, 0, t2029, t2031, t2021, t1987; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2109, t2111, t2121, -t2119, 0, 0, 0, 0, 0, 0, t2070, -t2236, t2093, -t2076, 0, 0, 0, 0, 0, 0, t2044, t2046, t2035, t1997; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2099, t2101, t2071, t2039, 0, 0, 0, 0, 0, 0, t2045, t2047, t2036, t1998; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2100, t2127, t2073, t2038, 0, 0, 0, 0, 0, 0, -t2056, -t2058, -t2086, -t2027; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2070, -t2236, t2093, -t2076, 0, 0, 0, 0, 0, 0, t2044, t2046, t2035, t1997; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2087, t2069, t2057, t2014; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2068, t2092, t2059, t2013; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2056, t2058, t2086, t2027;];
f_new_reg  = t1;