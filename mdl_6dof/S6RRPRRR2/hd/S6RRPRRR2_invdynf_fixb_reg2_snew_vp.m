% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 19:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RRPRRR2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRRR2_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_invdynf_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 19:59:15
% EndTime: 2019-05-06 19:59:25
% DurationCPUTime: 10.67s
% Computational Cost: add. (100862->345), mult. (234072->488), div. (0->0), fcn. (180019->12), ass. (0->275)
t2303 = cos(qJ(2));
t2288 = t2303 * qJDD(1);
t2298 = sin(qJ(2));
t2348 = qJD(1) * t2298;
t2329 = qJD(2) * t2348;
t2269 = t2288 - t2329;
t2292 = t2303 ^ 2;
t2306 = qJD(1) ^ 2;
t2318 = qJD(2) * pkin(2) - qJ(3) * t2348;
t2299 = sin(qJ(1));
t2304 = cos(qJ(1));
t2277 = t2299 * g(1) - t2304 * g(2);
t2319 = qJDD(1) * pkin(1) + t2277;
t2228 = t2269 * pkin(2) - t2318 * t2348 - qJDD(3) + t2319 + (qJ(3) * t2292 + pkin(7)) * t2306;
t2347 = qJD(1) * t2303;
t2328 = qJD(2) * t2347;
t2331 = t2298 * qJDD(1);
t2268 = t2328 + t2331;
t2293 = sin(pkin(11));
t2294 = cos(pkin(11));
t2249 = -t2268 * t2293 + t2294 * t2269;
t2261 = -t2293 * t2348 + t2294 * t2347;
t2260 = t2261 ^ 2;
t2262 = (t2293 * t2303 + t2294 * t2298) * qJD(1);
t2322 = qJD(2) * pkin(3) - pkin(8) * t2262;
t2199 = t2249 * pkin(3) + t2260 * pkin(8) - t2262 * t2322 + t2228;
t2297 = sin(qJ(4));
t2302 = cos(qJ(4));
t2241 = -t2302 * t2261 + t2262 * t2297;
t2240 = qJD(5) + t2241;
t2239 = qJD(6) + t2240;
t2360 = qJD(6) + t2239;
t2243 = t2261 * t2297 + t2262 * t2302;
t2290 = qJD(2) + qJD(4);
t2296 = sin(qJ(5));
t2301 = cos(qJ(5));
t2229 = t2243 * t2296 - t2301 * t2290;
t2231 = t2243 * t2301 + t2290 * t2296;
t2295 = sin(qJ(6));
t2300 = cos(qJ(6));
t2207 = t2300 * t2229 + t2231 * t2295;
t2359 = t2207 ^ 2;
t2209 = -t2229 * t2295 + t2231 * t2300;
t2358 = t2209 ^ 2;
t2357 = t2229 ^ 2;
t2356 = t2231 ^ 2;
t2355 = t2239 ^ 2;
t2354 = t2240 ^ 2;
t2353 = t2241 ^ 2;
t2352 = t2243 ^ 2;
t2351 = t2262 ^ 2;
t2350 = t2290 ^ 2;
t2349 = t2298 * g(3);
t2346 = qJD(2) * t2261;
t2345 = qJD(2) * t2262;
t2344 = t2207 * t2209;
t2343 = t2229 * t2231;
t2342 = t2240 * t2229;
t2341 = t2241 * t2243;
t2340 = t2261 * t2262;
t2339 = t2292 * t2306;
t2278 = -g(1) * t2304 - g(2) * t2299;
t2265 = -pkin(1) * t2306 + qJDD(1) * pkin(7) + t2278;
t2338 = t2298 * t2265;
t2337 = t2298 * t2306;
t2336 = t2303 * t2265;
t2335 = qJD(4) - t2290;
t2334 = qJD(5) - t2240;
t2333 = qJD(6) - t2239;
t2308 = -pkin(2) * t2339 + t2269 * qJ(3) - qJD(2) * t2318 - t2349;
t2314 = qJ(3) * qJD(1) * qJD(2) + pkin(2) * t2337 - g(3);
t2315 = qJDD(2) * pkin(2) - t2268 * qJ(3) - t2338;
t2201 = t2294 * (t2308 + t2336) + t2293 * (t2303 * t2314 + t2315) + 0.2e1 * qJD(3) * t2261;
t2181 = -t2260 * pkin(3) + t2249 * pkin(8) - qJD(2) * t2322 + t2201;
t2200 = -t2293 * t2308 + t2294 * t2315 - 0.2e1 * qJD(3) * t2262 + (-t2293 * t2265 + t2294 * t2314) * t2303;
t2250 = t2294 * t2268 + t2269 * t2293;
t2236 = -t2250 + t2346;
t2246 = qJDD(2) + t2340;
t2307 = pkin(3) * t2246 + pkin(8) * t2236 + t2200;
t2153 = t2302 * t2181 + t2297 * t2307;
t2217 = pkin(4) * t2241 - pkin(9) * t2243;
t2330 = qJDD(2) + qJDD(4);
t2142 = -pkin(4) * t2350 + pkin(9) * t2330 - t2241 * t2217 + t2153;
t2323 = -t2302 * t2249 + t2297 * t2250;
t2190 = (qJD(4) + t2290) * t2243 + t2323;
t2320 = -t2297 * t2249 - t2302 * t2250;
t2206 = -qJD(4) * t2241 - t2320;
t2325 = t2290 * t2241 - t2206;
t2151 = pkin(4) * t2190 + pkin(9) * t2325 - t2199;
t2118 = t2301 * t2142 + t2296 * t2151;
t2291 = t2298 ^ 2;
t2332 = t2291 + t2292;
t2117 = -t2296 * t2142 + t2301 * t2151;
t2152 = -t2181 * t2297 + t2302 * t2307;
t2316 = -t2301 * t2206 - t2296 * t2330;
t2185 = -t2229 * qJD(5) - t2316;
t2324 = t2296 * t2206 - t2301 * t2330;
t2317 = qJD(5) * t2231 + t2324;
t2326 = -t2295 * t2185 - t2300 * t2317;
t2313 = qJD(4) * t2243 + qJDD(5) + t2323;
t2141 = -t2330 * pkin(4) - t2350 * pkin(9) + t2217 * t2243 - t2152;
t2312 = -qJDD(6) - t2313;
t2178 = t2313 - t2343;
t2310 = -t2300 * t2185 + t2295 * t2317;
t2305 = qJD(2) ^ 2;
t2282 = t2303 * t2337;
t2281 = -t2305 - t2339;
t2280 = -t2291 * t2306 - t2305;
t2276 = -qJDD(2) + t2282;
t2275 = qJDD(2) + t2282;
t2274 = t2332 * t2306;
t2273 = -qJDD(1) * t2299 - t2304 * t2306;
t2272 = qJDD(1) * t2304 - t2299 * t2306;
t2271 = t2332 * qJDD(1);
t2270 = t2288 - 0.2e1 * t2329;
t2267 = 0.2e1 * t2328 + t2331;
t2264 = t2306 * pkin(7) + t2319;
t2257 = -t2305 - t2351;
t2256 = t2336 - t2349;
t2255 = -t2303 * g(3) - t2338;
t2254 = t2276 * t2303 - t2280 * t2298;
t2253 = -t2275 * t2298 + t2281 * t2303;
t2252 = t2276 * t2298 + t2280 * t2303;
t2251 = t2275 * t2303 + t2281 * t2298;
t2247 = -qJDD(2) + t2340;
t2245 = -t2305 - t2260;
t2237 = -t2350 - t2352;
t2235 = t2250 + t2346;
t2234 = t2249 + t2345;
t2233 = -t2249 + t2345;
t2232 = -t2260 - t2351;
t2225 = -t2255 * t2298 + t2256 * t2303;
t2224 = t2255 * t2303 + t2256 * t2298;
t2221 = t2247 * t2294 - t2257 * t2293;
t2220 = t2247 * t2293 + t2257 * t2294;
t2219 = t2245 * t2294 - t2246 * t2293;
t2218 = t2245 * t2293 + t2246 * t2294;
t2216 = -t2330 - t2341;
t2215 = t2330 - t2341;
t2214 = -t2350 - t2353;
t2212 = pkin(5) * t2240 - pkin(10) * t2231;
t2211 = t2234 * t2294 - t2236 * t2293;
t2210 = t2234 * t2293 + t2236 * t2294;
t2205 = -t2352 - t2353;
t2204 = -t2220 * t2298 + t2221 * t2303;
t2203 = t2220 * t2303 + t2221 * t2298;
t2202 = -t2354 - t2356;
t2198 = t2216 * t2302 - t2237 * t2297;
t2197 = t2216 * t2297 + t2237 * t2302;
t2196 = -t2354 - t2357;
t2194 = -t2356 - t2357;
t2193 = t2241 * t2335 + t2320;
t2191 = -t2243 * t2335 - t2323;
t2189 = -t2218 * t2298 + t2219 * t2303;
t2188 = t2218 * t2303 + t2219 * t2298;
t2187 = t2214 * t2302 - t2215 * t2297;
t2186 = t2214 * t2297 + t2215 * t2302;
t2182 = -t2355 - t2358;
t2179 = -t2313 - t2343;
t2177 = -t2210 * t2298 + t2211 * t2303;
t2176 = t2210 * t2303 + t2211 * t2298;
t2173 = -t2355 - t2359;
t2172 = t2229 * t2334 + t2316;
t2171 = t2185 - t2342;
t2170 = -t2231 * t2334 - t2324;
t2169 = (qJD(5) + t2240) * t2231 + t2324;
t2168 = -t2200 * t2293 + t2201 * t2294;
t2167 = t2200 * t2294 + t2201 * t2293;
t2166 = -t2197 * t2293 + t2198 * t2294;
t2165 = t2197 * t2294 + t2198 * t2293;
t2164 = t2191 * t2302 - t2193 * t2297;
t2163 = t2191 * t2297 + t2193 * t2302;
t2162 = -t2358 - t2359;
t2161 = t2312 - t2344;
t2160 = -t2312 - t2344;
t2159 = t2179 * t2301 - t2202 * t2296;
t2158 = t2179 * t2296 + t2202 * t2301;
t2157 = -t2186 * t2293 + t2187 * t2294;
t2156 = t2186 * t2294 + t2187 * t2293;
t2155 = -t2178 * t2296 + t2196 * t2301;
t2154 = t2178 * t2301 + t2196 * t2296;
t2148 = t2161 * t2300 - t2182 * t2295;
t2147 = t2161 * t2295 + t2182 * t2300;
t2146 = t2170 * t2301 - t2172 * t2296;
t2145 = t2170 * t2296 + t2172 * t2301;
t2144 = -t2167 * t2298 + t2168 * t2303;
t2143 = t2167 * t2303 + t2168 * t2298;
t2139 = -t2165 * t2298 + t2166 * t2303;
t2138 = t2165 * t2303 + t2166 * t2298;
t2137 = -t2160 * t2295 + t2173 * t2300;
t2136 = t2160 * t2300 + t2173 * t2295;
t2135 = t2207 * t2333 + t2310;
t2134 = -t2207 * t2360 - t2310;
t2133 = -t2209 * t2333 + t2326;
t2132 = t2209 * t2360 - t2326;
t2131 = t2159 * t2302 + t2171 * t2297;
t2130 = t2159 * t2297 - t2171 * t2302;
t2129 = t2155 * t2302 + t2169 * t2297;
t2128 = t2155 * t2297 - t2169 * t2302;
t2127 = -t2163 * t2293 + t2164 * t2294;
t2126 = t2163 * t2294 + t2164 * t2293;
t2125 = t2146 * t2302 + t2194 * t2297;
t2124 = t2146 * t2297 - t2194 * t2302;
t2123 = -t2156 * t2298 + t2157 * t2303;
t2122 = t2156 * t2303 + t2157 * t2298;
t2121 = -t2152 * t2297 + t2153 * t2302;
t2120 = t2152 * t2302 + t2153 * t2297;
t2119 = pkin(5) * t2317 - pkin(10) * t2357 + t2212 * t2231 + t2141;
t2116 = -t2147 * t2296 + t2148 * t2301;
t2115 = t2147 * t2301 + t2148 * t2296;
t2114 = -t2136 * t2296 + t2137 * t2301;
t2113 = t2136 * t2301 + t2137 * t2296;
t2112 = t2133 * t2300 - t2135 * t2295;
t2111 = t2133 * t2295 + t2135 * t2300;
t2110 = -t2130 * t2293 + t2131 * t2294;
t2109 = t2130 * t2294 + t2131 * t2293;
t2108 = -t2128 * t2293 + t2129 * t2294;
t2107 = t2128 * t2294 + t2129 * t2293;
t2106 = -t2126 * t2298 + t2127 * t2303;
t2105 = t2126 * t2303 + t2127 * t2298;
t2104 = -pkin(5) * t2357 - pkin(10) * t2317 - t2240 * t2212 + t2118;
t2103 = -t2124 * t2293 + t2125 * t2294;
t2102 = t2124 * t2294 + t2125 * t2293;
t2101 = (-t2185 - t2342) * pkin(10) + t2178 * pkin(5) + t2117;
t2100 = t2116 * t2302 + t2134 * t2297;
t2099 = t2116 * t2297 - t2134 * t2302;
t2098 = -t2120 * t2293 + t2121 * t2294;
t2097 = t2120 * t2294 + t2121 * t2293;
t2096 = t2114 * t2302 + t2132 * t2297;
t2095 = t2114 * t2297 - t2132 * t2302;
t2094 = -t2117 * t2296 + t2118 * t2301;
t2093 = t2117 * t2301 + t2118 * t2296;
t2092 = -t2111 * t2296 + t2112 * t2301;
t2091 = t2111 * t2301 + t2112 * t2296;
t2090 = -t2109 * t2298 + t2110 * t2303;
t2089 = t2109 * t2303 + t2110 * t2298;
t2088 = -t2107 * t2298 + t2108 * t2303;
t2087 = t2107 * t2303 + t2108 * t2298;
t2086 = t2094 * t2302 + t2141 * t2297;
t2085 = t2094 * t2297 - t2141 * t2302;
t2084 = t2092 * t2302 + t2162 * t2297;
t2083 = t2092 * t2297 - t2162 * t2302;
t2082 = -t2102 * t2298 + t2103 * t2303;
t2081 = t2102 * t2303 + t2103 * t2298;
t2080 = t2101 * t2295 + t2104 * t2300;
t2079 = t2101 * t2300 - t2104 * t2295;
t2078 = -t2099 * t2293 + t2100 * t2294;
t2077 = t2099 * t2294 + t2100 * t2293;
t2076 = -t2097 * t2298 + t2098 * t2303;
t2075 = t2097 * t2303 + t2098 * t2298;
t2074 = -t2095 * t2293 + t2096 * t2294;
t2073 = t2095 * t2294 + t2096 * t2293;
t2072 = -t2085 * t2293 + t2086 * t2294;
t2071 = t2085 * t2294 + t2086 * t2293;
t2070 = -t2083 * t2293 + t2084 * t2294;
t2069 = t2083 * t2294 + t2084 * t2293;
t2068 = -t2079 * t2295 + t2080 * t2300;
t2067 = t2079 * t2300 + t2080 * t2295;
t2066 = -t2077 * t2298 + t2078 * t2303;
t2065 = t2077 * t2303 + t2078 * t2298;
t2064 = -t2073 * t2298 + t2074 * t2303;
t2063 = t2073 * t2303 + t2074 * t2298;
t2062 = -t2071 * t2298 + t2072 * t2303;
t2061 = t2071 * t2303 + t2072 * t2298;
t2060 = -t2069 * t2298 + t2070 * t2303;
t2059 = t2069 * t2303 + t2070 * t2298;
t2058 = -t2067 * t2296 + t2068 * t2301;
t2057 = t2067 * t2301 + t2068 * t2296;
t2056 = t2058 * t2302 + t2119 * t2297;
t2055 = t2058 * t2297 - t2119 * t2302;
t2054 = -t2055 * t2293 + t2056 * t2294;
t2053 = t2055 * t2294 + t2056 * t2293;
t2052 = -t2053 * t2298 + t2054 * t2303;
t2051 = t2053 * t2303 + t2054 * t2298;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2273, -t2272, 0, -t2277 * t2299 + t2278 * t2304, 0, 0, 0, 0, 0, 0, t2253 * t2304 - t2270 * t2299, t2254 * t2304 + t2267 * t2299, t2271 * t2304 - t2274 * t2299, t2225 * t2304 - t2264 * t2299, 0, 0, 0, 0, 0, 0, t2189 * t2304 + t2233 * t2299, t2204 * t2304 + t2235 * t2299, t2177 * t2304 + t2232 * t2299, t2144 * t2304 - t2228 * t2299, 0, 0, 0, 0, 0, 0, t2123 * t2304 + t2190 * t2299, t2139 * t2304 - t2299 * t2325, t2106 * t2304 + t2205 * t2299, t2076 * t2304 - t2199 * t2299, 0, 0, 0, 0, 0, 0, t2088 * t2304 + t2154 * t2299, t2090 * t2304 + t2158 * t2299, t2082 * t2304 + t2145 * t2299, t2062 * t2304 + t2093 * t2299, 0, 0, 0, 0, 0, 0, t2064 * t2304 + t2113 * t2299, t2066 * t2304 + t2115 * t2299, t2060 * t2304 + t2091 * t2299, t2052 * t2304 + t2057 * t2299; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2272, t2273, 0, t2277 * t2304 + t2278 * t2299, 0, 0, 0, 0, 0, 0, t2253 * t2299 + t2270 * t2304, t2254 * t2299 - t2267 * t2304, t2271 * t2299 + t2274 * t2304, t2225 * t2299 + t2264 * t2304, 0, 0, 0, 0, 0, 0, t2189 * t2299 - t2233 * t2304, t2204 * t2299 - t2235 * t2304, t2177 * t2299 - t2232 * t2304, t2144 * t2299 + t2228 * t2304, 0, 0, 0, 0, 0, 0, t2123 * t2299 - t2190 * t2304, t2139 * t2299 + t2304 * t2325, t2106 * t2299 - t2205 * t2304, t2076 * t2299 + t2199 * t2304, 0, 0, 0, 0, 0, 0, t2088 * t2299 - t2154 * t2304, t2090 * t2299 - t2158 * t2304, t2082 * t2299 - t2145 * t2304, t2062 * t2299 - t2093 * t2304, 0, 0, 0, 0, 0, 0, t2064 * t2299 - t2113 * t2304, t2066 * t2299 - t2115 * t2304, t2060 * t2299 - t2091 * t2304, t2052 * t2299 - t2057 * t2304; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2251, t2252, 0, t2224, 0, 0, 0, 0, 0, 0, t2188, t2203, t2176, t2143, 0, 0, 0, 0, 0, 0, t2122, t2138, t2105, t2075, 0, 0, 0, 0, 0, 0, t2087, t2089, t2081, t2061, 0, 0, 0, 0, 0, 0, t2063, t2065, t2059, t2051; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2306, -qJDD(1), 0, t2278, 0, 0, 0, 0, 0, 0, t2253, t2254, t2271, t2225, 0, 0, 0, 0, 0, 0, t2189, t2204, t2177, t2144, 0, 0, 0, 0, 0, 0, t2123, t2139, t2106, t2076, 0, 0, 0, 0, 0, 0, t2088, t2090, t2082, t2062, 0, 0, 0, 0, 0, 0, t2064, t2066, t2060, t2052; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2306, 0, t2277, 0, 0, 0, 0, 0, 0, t2270, -t2267, t2274, t2264, 0, 0, 0, 0, 0, 0, -t2233, -t2235, -t2232, t2228, 0, 0, 0, 0, 0, 0, -t2190, t2325, -t2205, t2199, 0, 0, 0, 0, 0, 0, -t2154, -t2158, -t2145, -t2093, 0, 0, 0, 0, 0, 0, -t2113, -t2115, -t2091, -t2057; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2251, t2252, 0, t2224, 0, 0, 0, 0, 0, 0, t2188, t2203, t2176, t2143, 0, 0, 0, 0, 0, 0, t2122, t2138, t2105, t2075, 0, 0, 0, 0, 0, 0, t2087, t2089, t2081, t2061, 0, 0, 0, 0, 0, 0, t2063, t2065, t2059, t2051; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2281, t2276, t2288, t2256, 0, 0, 0, 0, 0, 0, t2219, t2221, t2211, t2168, 0, 0, 0, 0, 0, 0, t2157, t2166, t2127, t2098, 0, 0, 0, 0, 0, 0, t2108, t2110, t2103, t2072, 0, 0, 0, 0, 0, 0, t2074, t2078, t2070, t2054; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2275, t2280, -t2331, t2255, 0, 0, 0, 0, 0, 0, t2218, t2220, t2210, t2167, 0, 0, 0, 0, 0, 0, t2156, t2165, t2126, t2097, 0, 0, 0, 0, 0, 0, t2107, t2109, t2102, t2071, 0, 0, 0, 0, 0, 0, t2073, t2077, t2069, t2053; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2270, t2267, -t2274, -t2264, 0, 0, 0, 0, 0, 0, t2233, t2235, t2232, -t2228, 0, 0, 0, 0, 0, 0, t2190, -t2325, t2205, -t2199, 0, 0, 0, 0, 0, 0, t2154, t2158, t2145, t2093, 0, 0, 0, 0, 0, 0, t2113, t2115, t2091, t2057; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2245, t2247, t2234, t2201, 0, 0, 0, 0, 0, 0, t2187, t2198, t2164, t2121, 0, 0, 0, 0, 0, 0, t2129, t2131, t2125, t2086, 0, 0, 0, 0, 0, 0, t2096, t2100, t2084, t2056; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2246, t2257, t2236, t2200, 0, 0, 0, 0, 0, 0, t2186, t2197, t2163, t2120, 0, 0, 0, 0, 0, 0, t2128, t2130, t2124, t2085, 0, 0, 0, 0, 0, 0, t2095, t2099, t2083, t2055; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2233, t2235, t2232, -t2228, 0, 0, 0, 0, 0, 0, t2190, -t2325, t2205, -t2199, 0, 0, 0, 0, 0, 0, t2154, t2158, t2145, t2093, 0, 0, 0, 0, 0, 0, t2113, t2115, t2091, t2057; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2214, t2216, t2191, t2153, 0, 0, 0, 0, 0, 0, t2155, t2159, t2146, t2094, 0, 0, 0, 0, 0, 0, t2114, t2116, t2092, t2058; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2215, t2237, t2193, t2152, 0, 0, 0, 0, 0, 0, -t2169, -t2171, -t2194, -t2141, 0, 0, 0, 0, 0, 0, -t2132, -t2134, -t2162, -t2119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2190, -t2325, t2205, -t2199, 0, 0, 0, 0, 0, 0, t2154, t2158, t2145, t2093, 0, 0, 0, 0, 0, 0, t2113, t2115, t2091, t2057; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2196, t2179, t2170, t2118, 0, 0, 0, 0, 0, 0, t2137, t2148, t2112, t2068; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2178, t2202, t2172, t2117, 0, 0, 0, 0, 0, 0, t2136, t2147, t2111, t2067; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2169, t2171, t2194, t2141, 0, 0, 0, 0, 0, 0, t2132, t2134, t2162, t2119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2173, t2161, t2133, t2080; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2160, t2182, t2135, t2079; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2132, t2134, t2162, t2119;];
f_new_reg  = t1;