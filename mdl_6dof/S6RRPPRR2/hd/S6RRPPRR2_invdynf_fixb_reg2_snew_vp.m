% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RRPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 09:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RRPPRR2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR2_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR2_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_invdynf_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:49:23
% EndTime: 2019-05-06 09:49:33
% DurationCPUTime: 9.82s
% Computational Cost: add. (98425->341), mult. (235040->482), div. (0->0), fcn. (173900->12), ass. (0->272)
t2271 = cos(qJ(2));
t2309 = qJD(1) * t2271;
t2291 = qJD(2) * t2309;
t2267 = sin(qJ(2));
t2293 = t2267 * qJDD(1);
t2236 = t2291 + t2293;
t2255 = t2271 * qJDD(1);
t2310 = qJD(1) * t2267;
t2292 = qJD(2) * t2310;
t2237 = t2255 - t2292;
t2262 = sin(pkin(10));
t2264 = cos(pkin(10));
t2295 = -t2262 * t2236 + t2264 * t2237;
t2230 = t2262 * t2309 + t2264 * t2310;
t2261 = sin(pkin(11));
t2263 = cos(pkin(11));
t2214 = -qJD(2) * t2263 + t2230 * t2261;
t2215 = qJD(2) * t2261 + t2230 * t2263;
t2304 = t2214 * t2215;
t2170 = -t2295 - t2304;
t2307 = qJD(2) * t2230;
t2194 = -t2295 + t2307;
t2228 = t2262 * t2310 - t2264 * t2309;
t2227 = qJD(5) + t2228;
t2223 = qJD(6) + t2227;
t2323 = qJD(6) + t2223;
t2259 = t2271 ^ 2;
t2273 = qJD(1) ^ 2;
t2281 = qJD(2) * pkin(2) - qJ(3) * t2310;
t2268 = sin(qJ(1));
t2272 = cos(qJ(1));
t2245 = t2268 * g(1) - t2272 * g(2);
t2282 = qJDD(1) * pkin(1) + t2245;
t2192 = t2237 * pkin(2) + (qJ(3) * t2259 + pkin(7)) * t2273 - t2281 * t2310 - qJDD(3) + t2282;
t2322 = qJD(2) ^ 2;
t2266 = sin(qJ(5));
t2270 = cos(qJ(5));
t2185 = t2270 * t2214 + t2215 * t2266;
t2187 = -t2214 * t2266 + t2215 * t2270;
t2265 = sin(qJ(6));
t2269 = cos(qJ(6));
t2160 = t2269 * t2185 + t2187 * t2265;
t2321 = t2160 ^ 2;
t2162 = -t2185 * t2265 + t2187 * t2269;
t2320 = t2162 ^ 2;
t2319 = t2185 ^ 2;
t2318 = t2187 ^ 2;
t2317 = t2214 ^ 2;
t2316 = t2215 ^ 2;
t2315 = t2223 ^ 2;
t2314 = t2227 ^ 2;
t2207 = t2228 ^ 2;
t2313 = t2230 ^ 2;
t2312 = -2 * qJD(3);
t2311 = -2 * qJD(4);
t2308 = qJD(2) * t2228;
t2306 = t2160 * t2162;
t2305 = t2185 * t2187;
t2303 = t2215 * t2228;
t2302 = t2227 * t2185;
t2301 = t2228 * t2214;
t2300 = t2228 * t2230;
t2299 = t2259 * t2273;
t2298 = t2267 * t2273;
t2297 = qJD(5) - t2227;
t2296 = qJD(6) - t2223;
t2246 = -g(1) * t2272 - g(2) * t2268;
t2278 = -pkin(1) * t2273 + qJDD(1) * pkin(7) + t2246;
t2221 = -t2267 * g(3) + t2271 * t2278;
t2189 = -pkin(2) * t2299 + t2237 * qJ(3) - qJD(2) * t2281 + t2221;
t2276 = t2267 * t2278;
t2274 = -t2276 - t2236 * qJ(3) + qJDD(2) * pkin(2) + (qJ(3) * qJD(1) * qJD(2) + pkin(2) * t2298 - g(3)) * t2271;
t2156 = t2264 * t2189 + t2228 * t2312 + t2262 * t2274;
t2203 = pkin(3) * t2228 - qJ(4) * t2230;
t2135 = -pkin(3) * t2322 + qJDD(2) * qJ(4) - t2203 * t2228 + t2156;
t2209 = t2236 * t2264 + t2237 * t2262;
t2288 = -t2209 + t2308;
t2146 = pkin(3) * t2194 + t2288 * qJ(4) - t2192;
t2108 = -t2261 * t2135 + t2263 * t2146 + t2215 * t2311;
t2202 = qJDD(2) * t2261 + t2209 * t2263;
t2169 = -t2202 - t2301;
t2093 = pkin(4) * t2170 + t2169 * pkin(8) + t2108;
t2109 = t2263 * t2135 + t2261 * t2146 + t2214 * t2311;
t2198 = pkin(4) * t2228 - pkin(8) * t2215;
t2284 = t2263 * qJDD(2) - t2209 * t2261;
t2097 = -pkin(4) * t2317 + pkin(8) * t2284 - t2228 * t2198 + t2109;
t2067 = t2266 * t2093 + t2270 * t2097;
t2258 = t2267 ^ 2;
t2294 = t2258 + t2259;
t2289 = -qJDD(5) + t2295;
t2066 = t2270 * t2093 - t2266 * t2097;
t2279 = -t2270 * t2202 - t2266 * t2284;
t2153 = -t2185 * qJD(5) - t2279;
t2285 = t2266 * t2202 - t2270 * t2284;
t2280 = qJD(5) * t2187 + t2285;
t2287 = -t2265 * t2153 - t2269 * t2280;
t2286 = t2262 * t2189 - t2264 * t2274;
t2283 = -qJDD(6) + t2289;
t2151 = -t2289 - t2305;
t2275 = -t2269 * t2153 + t2265 * t2280;
t2134 = qJDD(4) - t2322 * qJ(4) - qJDD(2) * pkin(3) + ((2 * qJD(3)) + t2203) * t2230 + t2286;
t2111 = -t2284 * pkin(4) - t2317 * pkin(8) + t2215 * t2198 + t2134;
t2252 = t2271 * t2298;
t2251 = -t2299 - t2322;
t2250 = -t2258 * t2273 - t2322;
t2244 = -qJDD(2) + t2252;
t2243 = qJDD(2) + t2252;
t2242 = t2294 * t2273;
t2241 = -qJDD(1) * t2268 - t2272 * t2273;
t2240 = qJDD(1) * t2272 - t2268 * t2273;
t2239 = t2294 * qJDD(1);
t2238 = t2255 - 0.2e1 * t2292;
t2235 = 0.2e1 * t2291 + t2293;
t2233 = t2273 * pkin(7) + t2282;
t2222 = -t2313 - t2322;
t2220 = -t2271 * g(3) - t2276;
t2219 = t2244 * t2271 - t2250 * t2267;
t2218 = -t2243 * t2267 + t2251 * t2271;
t2217 = t2244 * t2267 + t2250 * t2271;
t2216 = t2243 * t2271 + t2251 * t2267;
t2206 = -qJDD(2) - t2300;
t2205 = qJDD(2) - t2300;
t2204 = -t2207 - t2322;
t2197 = -t2209 - t2308;
t2195 = t2295 + t2307;
t2193 = -t2313 - t2207;
t2191 = -t2220 * t2267 + t2221 * t2271;
t2190 = t2220 * t2271 + t2221 * t2267;
t2182 = -t2207 - t2316;
t2179 = t2206 * t2264 - t2222 * t2262;
t2178 = t2206 * t2262 + t2222 * t2264;
t2177 = -t2207 - t2317;
t2175 = -t2316 - t2317;
t2174 = t2204 * t2264 - t2205 * t2262;
t2173 = t2204 * t2262 + t2205 * t2264;
t2172 = pkin(5) * t2227 - pkin(9) * t2187;
t2171 = t2295 - t2304;
t2168 = t2202 - t2301;
t2167 = t2284 + t2303;
t2166 = -t2284 + t2303;
t2165 = -t2314 - t2318;
t2164 = t2195 * t2264 - t2197 * t2262;
t2163 = t2195 * t2262 + t2197 * t2264;
t2159 = -t2178 * t2267 + t2179 * t2271;
t2158 = t2178 * t2271 + t2179 * t2267;
t2157 = -t2314 - t2319;
t2155 = t2230 * t2312 - t2286;
t2152 = t2289 - t2305;
t2150 = t2171 * t2263 - t2182 * t2261;
t2149 = t2171 * t2261 + t2182 * t2263;
t2143 = -t2170 * t2261 + t2177 * t2263;
t2142 = t2170 * t2263 + t2177 * t2261;
t2141 = -t2315 - t2320;
t2140 = -t2173 * t2267 + t2174 * t2271;
t2139 = t2173 * t2271 + t2174 * t2267;
t2138 = t2167 * t2263 - t2169 * t2261;
t2137 = t2167 * t2261 + t2169 * t2263;
t2136 = -t2318 - t2319;
t2132 = -t2163 * t2267 + t2164 * t2271;
t2131 = t2163 * t2271 + t2164 * t2267;
t2130 = t2185 * t2297 + t2279;
t2129 = t2153 - t2302;
t2128 = -t2187 * t2297 - t2285;
t2127 = (qJD(5) + t2227) * t2187 + t2285;
t2126 = t2150 * t2264 + t2168 * t2262;
t2125 = t2150 * t2262 - t2168 * t2264;
t2124 = t2143 * t2264 + t2166 * t2262;
t2123 = t2143 * t2262 - t2166 * t2264;
t2122 = t2152 * t2270 - t2165 * t2266;
t2121 = t2152 * t2266 + t2165 * t2270;
t2120 = t2138 * t2264 + t2175 * t2262;
t2119 = t2138 * t2262 - t2175 * t2264;
t2118 = -t2315 - t2321;
t2117 = -t2155 * t2262 + t2156 * t2264;
t2116 = t2155 * t2264 + t2156 * t2262;
t2115 = t2283 - t2306;
t2114 = -t2283 - t2306;
t2113 = -t2151 * t2266 + t2157 * t2270;
t2112 = t2151 * t2270 + t2157 * t2266;
t2110 = -t2320 - t2321;
t2107 = t2115 * t2269 - t2141 * t2265;
t2106 = t2115 * t2265 + t2141 * t2269;
t2105 = t2128 * t2270 - t2130 * t2266;
t2104 = t2128 * t2266 + t2130 * t2270;
t2103 = -t2125 * t2267 + t2126 * t2271;
t2102 = t2125 * t2271 + t2126 * t2267;
t2101 = -t2123 * t2267 + t2124 * t2271;
t2100 = t2123 * t2271 + t2124 * t2267;
t2099 = -t2121 * t2261 + t2122 * t2263;
t2098 = t2121 * t2263 + t2122 * t2261;
t2095 = -t2119 * t2267 + t2120 * t2271;
t2094 = t2119 * t2271 + t2120 * t2267;
t2092 = -t2114 * t2265 + t2118 * t2269;
t2091 = t2114 * t2269 + t2118 * t2265;
t2088 = -t2116 * t2267 + t2117 * t2271;
t2087 = t2116 * t2271 + t2117 * t2267;
t2086 = t2160 * t2296 + t2275;
t2085 = -t2160 * t2323 - t2275;
t2084 = -t2162 * t2296 + t2287;
t2083 = t2162 * t2323 - t2287;
t2082 = -t2112 * t2261 + t2113 * t2263;
t2081 = t2112 * t2263 + t2113 * t2261;
t2080 = pkin(5) * t2280 - pkin(9) * t2319 + t2187 * t2172 + t2111;
t2079 = t2099 * t2264 + t2129 * t2262;
t2078 = t2099 * t2262 - t2129 * t2264;
t2077 = -t2108 * t2261 + t2109 * t2263;
t2076 = t2108 * t2263 + t2109 * t2261;
t2075 = t2082 * t2264 + t2127 * t2262;
t2074 = t2082 * t2262 - t2127 * t2264;
t2073 = -t2106 * t2266 + t2107 * t2270;
t2072 = t2106 * t2270 + t2107 * t2266;
t2071 = -t2104 * t2261 + t2105 * t2263;
t2070 = t2104 * t2263 + t2105 * t2261;
t2069 = t2077 * t2264 + t2134 * t2262;
t2068 = t2077 * t2262 - t2134 * t2264;
t2065 = t2071 * t2264 + t2136 * t2262;
t2064 = t2071 * t2262 - t2136 * t2264;
t2063 = -t2091 * t2266 + t2092 * t2270;
t2062 = t2091 * t2270 + t2092 * t2266;
t2061 = t2084 * t2269 - t2086 * t2265;
t2060 = t2084 * t2265 + t2086 * t2269;
t2059 = -pkin(5) * t2319 - pkin(9) * t2280 - t2227 * t2172 + t2067;
t2058 = (-t2153 - t2302) * pkin(9) + t2151 * pkin(5) + t2066;
t2057 = -t2078 * t2267 + t2079 * t2271;
t2056 = t2078 * t2271 + t2079 * t2267;
t2055 = -t2074 * t2267 + t2075 * t2271;
t2054 = t2074 * t2271 + t2075 * t2267;
t2053 = -t2072 * t2261 + t2073 * t2263;
t2052 = t2072 * t2263 + t2073 * t2261;
t2051 = -t2068 * t2267 + t2069 * t2271;
t2050 = t2068 * t2271 + t2069 * t2267;
t2049 = -t2066 * t2266 + t2067 * t2270;
t2048 = t2066 * t2270 + t2067 * t2266;
t2047 = -t2064 * t2267 + t2065 * t2271;
t2046 = t2064 * t2271 + t2065 * t2267;
t2045 = -t2062 * t2261 + t2063 * t2263;
t2044 = t2062 * t2263 + t2063 * t2261;
t2043 = -t2060 * t2266 + t2061 * t2270;
t2042 = t2060 * t2270 + t2061 * t2266;
t2041 = t2053 * t2264 + t2085 * t2262;
t2040 = t2053 * t2262 - t2085 * t2264;
t2039 = t2045 * t2264 + t2083 * t2262;
t2038 = t2045 * t2262 - t2083 * t2264;
t2037 = t2058 * t2265 + t2059 * t2269;
t2036 = t2058 * t2269 - t2059 * t2265;
t2035 = -t2048 * t2261 + t2049 * t2263;
t2034 = t2048 * t2263 + t2049 * t2261;
t2033 = -t2042 * t2261 + t2043 * t2263;
t2032 = t2042 * t2263 + t2043 * t2261;
t2031 = t2035 * t2264 + t2111 * t2262;
t2030 = t2035 * t2262 - t2111 * t2264;
t2029 = -t2040 * t2267 + t2041 * t2271;
t2028 = t2040 * t2271 + t2041 * t2267;
t2027 = t2033 * t2264 + t2110 * t2262;
t2026 = t2033 * t2262 - t2110 * t2264;
t2025 = -t2038 * t2267 + t2039 * t2271;
t2024 = t2038 * t2271 + t2039 * t2267;
t2023 = -t2036 * t2265 + t2037 * t2269;
t2022 = t2036 * t2269 + t2037 * t2265;
t2021 = -t2030 * t2267 + t2031 * t2271;
t2020 = t2030 * t2271 + t2031 * t2267;
t2019 = -t2026 * t2267 + t2027 * t2271;
t2018 = t2026 * t2271 + t2027 * t2267;
t2017 = -t2022 * t2266 + t2023 * t2270;
t2016 = t2022 * t2270 + t2023 * t2266;
t2015 = -t2016 * t2261 + t2017 * t2263;
t2014 = t2016 * t2263 + t2017 * t2261;
t2013 = t2015 * t2264 + t2080 * t2262;
t2012 = t2015 * t2262 - t2080 * t2264;
t2011 = -t2012 * t2267 + t2013 * t2271;
t2010 = t2012 * t2271 + t2013 * t2267;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2241, -t2240, 0, -t2245 * t2268 + t2246 * t2272, 0, 0, 0, 0, 0, 0, t2218 * t2272 - t2238 * t2268, t2219 * t2272 + t2235 * t2268, t2239 * t2272 - t2242 * t2268, t2191 * t2272 - t2233 * t2268, 0, 0, 0, 0, 0, 0, t2140 * t2272 + t2194 * t2268, t2159 * t2272 - t2268 * t2288, t2132 * t2272 + t2193 * t2268, t2088 * t2272 - t2192 * t2268, 0, 0, 0, 0, 0, 0, t2101 * t2272 + t2142 * t2268, t2103 * t2272 + t2149 * t2268, t2095 * t2272 + t2137 * t2268, t2051 * t2272 + t2076 * t2268, 0, 0, 0, 0, 0, 0, t2055 * t2272 + t2081 * t2268, t2057 * t2272 + t2098 * t2268, t2047 * t2272 + t2070 * t2268, t2021 * t2272 + t2034 * t2268, 0, 0, 0, 0, 0, 0, t2025 * t2272 + t2044 * t2268, t2029 * t2272 + t2052 * t2268, t2019 * t2272 + t2032 * t2268, t2011 * t2272 + t2014 * t2268; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2240, t2241, 0, t2245 * t2272 + t2246 * t2268, 0, 0, 0, 0, 0, 0, t2218 * t2268 + t2238 * t2272, t2219 * t2268 - t2235 * t2272, t2239 * t2268 + t2242 * t2272, t2191 * t2268 + t2233 * t2272, 0, 0, 0, 0, 0, 0, t2140 * t2268 - t2194 * t2272, t2159 * t2268 + t2272 * t2288, t2132 * t2268 - t2193 * t2272, t2088 * t2268 + t2192 * t2272, 0, 0, 0, 0, 0, 0, t2101 * t2268 - t2142 * t2272, t2103 * t2268 - t2149 * t2272, t2095 * t2268 - t2137 * t2272, t2051 * t2268 - t2076 * t2272, 0, 0, 0, 0, 0, 0, t2055 * t2268 - t2081 * t2272, t2057 * t2268 - t2098 * t2272, t2047 * t2268 - t2070 * t2272, t2021 * t2268 - t2034 * t2272, 0, 0, 0, 0, 0, 0, t2025 * t2268 - t2044 * t2272, t2029 * t2268 - t2052 * t2272, t2019 * t2268 - t2032 * t2272, t2011 * t2268 - t2014 * t2272; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2216, t2217, 0, t2190, 0, 0, 0, 0, 0, 0, t2139, t2158, t2131, t2087, 0, 0, 0, 0, 0, 0, t2100, t2102, t2094, t2050, 0, 0, 0, 0, 0, 0, t2054, t2056, t2046, t2020, 0, 0, 0, 0, 0, 0, t2024, t2028, t2018, t2010; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2273, -qJDD(1), 0, t2246, 0, 0, 0, 0, 0, 0, t2218, t2219, t2239, t2191, 0, 0, 0, 0, 0, 0, t2140, t2159, t2132, t2088, 0, 0, 0, 0, 0, 0, t2101, t2103, t2095, t2051, 0, 0, 0, 0, 0, 0, t2055, t2057, t2047, t2021, 0, 0, 0, 0, 0, 0, t2025, t2029, t2019, t2011; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2273, 0, t2245, 0, 0, 0, 0, 0, 0, t2238, -t2235, t2242, t2233, 0, 0, 0, 0, 0, 0, -t2194, t2288, -t2193, t2192, 0, 0, 0, 0, 0, 0, -t2142, -t2149, -t2137, -t2076, 0, 0, 0, 0, 0, 0, -t2081, -t2098, -t2070, -t2034, 0, 0, 0, 0, 0, 0, -t2044, -t2052, -t2032, -t2014; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2216, t2217, 0, t2190, 0, 0, 0, 0, 0, 0, t2139, t2158, t2131, t2087, 0, 0, 0, 0, 0, 0, t2100, t2102, t2094, t2050, 0, 0, 0, 0, 0, 0, t2054, t2056, t2046, t2020, 0, 0, 0, 0, 0, 0, t2024, t2028, t2018, t2010; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2251, t2244, t2255, t2221, 0, 0, 0, 0, 0, 0, t2174, t2179, t2164, t2117, 0, 0, 0, 0, 0, 0, t2124, t2126, t2120, t2069, 0, 0, 0, 0, 0, 0, t2075, t2079, t2065, t2031, 0, 0, 0, 0, 0, 0, t2039, t2041, t2027, t2013; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2243, t2250, -t2293, t2220, 0, 0, 0, 0, 0, 0, t2173, t2178, t2163, t2116, 0, 0, 0, 0, 0, 0, t2123, t2125, t2119, t2068, 0, 0, 0, 0, 0, 0, t2074, t2078, t2064, t2030, 0, 0, 0, 0, 0, 0, t2038, t2040, t2026, t2012; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2238, t2235, -t2242, -t2233, 0, 0, 0, 0, 0, 0, t2194, -t2288, t2193, -t2192, 0, 0, 0, 0, 0, 0, t2142, t2149, t2137, t2076, 0, 0, 0, 0, 0, 0, t2081, t2098, t2070, t2034, 0, 0, 0, 0, 0, 0, t2044, t2052, t2032, t2014; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2204, t2206, t2195, t2156, 0, 0, 0, 0, 0, 0, t2143, t2150, t2138, t2077, 0, 0, 0, 0, 0, 0, t2082, t2099, t2071, t2035, 0, 0, 0, 0, 0, 0, t2045, t2053, t2033, t2015; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2205, t2222, t2197, t2155, 0, 0, 0, 0, 0, 0, -t2166, -t2168, -t2175, -t2134, 0, 0, 0, 0, 0, 0, -t2127, -t2129, -t2136, -t2111, 0, 0, 0, 0, 0, 0, -t2083, -t2085, -t2110, -t2080; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2194, -t2288, t2193, -t2192, 0, 0, 0, 0, 0, 0, t2142, t2149, t2137, t2076, 0, 0, 0, 0, 0, 0, t2081, t2098, t2070, t2034, 0, 0, 0, 0, 0, 0, t2044, t2052, t2032, t2014; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2177, t2171, t2167, t2109, 0, 0, 0, 0, 0, 0, t2113, t2122, t2105, t2049, 0, 0, 0, 0, 0, 0, t2063, t2073, t2043, t2017; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2170, t2182, t2169, t2108, 0, 0, 0, 0, 0, 0, t2112, t2121, t2104, t2048, 0, 0, 0, 0, 0, 0, t2062, t2072, t2042, t2016; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2166, t2168, t2175, t2134, 0, 0, 0, 0, 0, 0, t2127, t2129, t2136, t2111, 0, 0, 0, 0, 0, 0, t2083, t2085, t2110, t2080; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2157, t2152, t2128, t2067, 0, 0, 0, 0, 0, 0, t2092, t2107, t2061, t2023; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2151, t2165, t2130, t2066, 0, 0, 0, 0, 0, 0, t2091, t2106, t2060, t2022; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2127, t2129, t2136, t2111, 0, 0, 0, 0, 0, 0, t2083, t2085, t2110, t2080; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2118, t2115, t2084, t2037; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2114, t2141, t2086, t2036; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2083, t2085, t2110, t2080;];
f_new_reg  = t1;
