% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RRRRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-08 09:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RRRRRR3_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR3_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR3_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRR3_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR3_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR3_invdynf_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 09:01:20
% EndTime: 2019-05-08 09:01:31
% DurationCPUTime: 11.50s
% Computational Cost: add. (133175->346), mult. (265920->485), div. (0->0), fcn. (200662->12), ass. (0->277)
t2332 = sin(qJ(3));
t2338 = cos(qJ(3));
t2339 = cos(qJ(2));
t2381 = qJD(1) * t2339;
t2333 = sin(qJ(2));
t2382 = qJD(1) * t2333;
t2296 = t2332 * t2382 - t2338 * t2381;
t2295 = qJD(4) + t2296;
t2293 = qJD(5) + t2295;
t2290 = qJD(6) + t2293;
t2395 = qJD(6) + t2290;
t2324 = t2339 * qJDD(1);
t2364 = qJD(2) * t2382;
t2305 = t2324 - t2364;
t2328 = t2339 ^ 2;
t2342 = qJD(1) ^ 2;
t2334 = sin(qJ(1));
t2340 = cos(qJ(1));
t2313 = t2334 * g(1) - g(2) * t2340;
t2353 = qJDD(1) * pkin(1) + t2313;
t2354 = qJD(2) * pkin(2) - pkin(8) * t2382;
t2272 = t2305 * pkin(2) + (pkin(8) * t2328 + pkin(7)) * t2342 - t2354 * t2382 + t2353;
t2298 = (t2332 * t2339 + t2333 * t2338) * qJD(1);
t2326 = qJD(2) + qJD(3);
t2331 = sin(qJ(4));
t2337 = cos(qJ(4));
t2281 = t2298 * t2331 - t2326 * t2337;
t2282 = t2298 * t2337 + t2326 * t2331;
t2330 = sin(qJ(5));
t2336 = cos(qJ(5));
t2261 = t2281 * t2336 + t2282 * t2330;
t2263 = -t2281 * t2330 + t2282 * t2336;
t2329 = sin(qJ(6));
t2335 = cos(qJ(6));
t2235 = t2261 * t2335 + t2263 * t2329;
t2394 = t2235 ^ 2;
t2237 = -t2261 * t2329 + t2263 * t2335;
t2393 = t2237 ^ 2;
t2392 = t2261 ^ 2;
t2391 = t2263 ^ 2;
t2390 = t2281 ^ 2;
t2389 = t2282 ^ 2;
t2388 = t2290 ^ 2;
t2387 = t2293 ^ 2;
t2386 = t2295 ^ 2;
t2385 = t2296 ^ 2;
t2384 = t2298 ^ 2;
t2383 = t2326 ^ 2;
t2380 = t2235 * t2237;
t2379 = t2261 * t2263;
t2378 = t2281 * t2282;
t2377 = t2293 * t2261;
t2376 = t2295 * t2281;
t2375 = t2296 * t2298;
t2374 = t2298 * t2326;
t2373 = t2328 * t2342;
t2372 = t2333 * t2342;
t2363 = qJD(2) * t2381;
t2367 = t2333 * qJDD(1);
t2304 = t2363 + t2367;
t2362 = -qJD(3) * t2298 - t2304 * t2332 + t2305 * t2338;
t2252 = -t2362 + t2374;
t2356 = -t2338 * t2304 - t2332 * t2305;
t2270 = -qJD(3) * t2296 - t2356;
t2359 = t2296 * t2326 - t2270;
t2211 = pkin(3) * t2252 + pkin(9) * t2359 - t2272;
t2314 = -g(1) * t2340 - g(2) * t2334;
t2347 = -pkin(1) * t2342 + qJDD(1) * pkin(7) + t2314;
t2289 = -t2333 * g(3) + t2339 * t2347;
t2269 = -pkin(2) * t2373 + t2305 * pkin(8) - qJD(2) * t2354 + t2289;
t2346 = t2333 * t2347;
t2343 = -t2346 - t2304 * pkin(8) + qJDD(2) * pkin(2) + (pkin(8) * qJD(1) * qJD(2) + pkin(2) * t2372 - g(3)) * t2339;
t2239 = t2269 * t2338 + t2332 * t2343;
t2278 = pkin(3) * t2296 - pkin(9) * t2298;
t2366 = qJDD(2) + qJDD(3);
t2217 = -pkin(3) * t2383 + pkin(9) * t2366 - t2278 * t2296 + t2239;
t2180 = t2211 * t2337 - t2331 * t2217;
t2357 = -qJDD(4) + t2362;
t2240 = -t2357 - t2378;
t2349 = -t2337 * t2270 - t2331 * t2366;
t2245 = -t2281 * qJD(4) - t2349;
t2168 = (-t2245 - t2376) * pkin(10) + t2240 * pkin(4) + t2180;
t2181 = t2211 * t2331 + t2217 * t2337;
t2273 = pkin(4) * t2295 - pkin(10) * t2282;
t2358 = t2331 * t2270 - t2337 * t2366;
t2352 = -qJD(4) * t2282 - t2358;
t2174 = -pkin(4) * t2390 + pkin(10) * t2352 - t2273 * t2295 + t2181;
t2140 = t2168 * t2330 + t2174 * t2336;
t2371 = -t2290 + qJD(6);
t2370 = -t2293 + qJD(5);
t2369 = -t2295 + qJD(4);
t2327 = t2333 ^ 2;
t2368 = t2327 + t2328;
t2139 = t2168 * t2336 - t2330 * t2174;
t2344 = -t2336 * t2245 - t2330 * t2352;
t2206 = -t2261 * qJD(5) - t2344;
t2360 = t2330 * t2245 - t2336 * t2352;
t2351 = qJD(5) * t2263 + t2360;
t2361 = -t2329 * t2206 - t2335 * t2351;
t2238 = -t2269 * t2332 + t2338 * t2343;
t2355 = -qJDD(5) + t2357;
t2350 = -qJDD(6) + t2355;
t2220 = -t2355 - t2379;
t2216 = -pkin(3) * t2366 - pkin(9) * t2383 + t2278 * t2298 - t2238;
t2345 = -t2335 * t2206 + t2329 * t2351;
t2182 = -pkin(4) * t2352 - pkin(10) * t2390 + t2273 * t2282 + t2216;
t2341 = qJD(2) ^ 2;
t2318 = t2339 * t2372;
t2316 = -t2341 - t2373;
t2315 = -t2327 * t2342 - t2341;
t2312 = -qJDD(2) + t2318;
t2311 = qJDD(2) + t2318;
t2310 = t2368 * t2342;
t2309 = -qJDD(1) * t2334 - t2340 * t2342;
t2308 = qJDD(1) * t2340 - t2334 * t2342;
t2307 = t2368 * qJDD(1);
t2306 = t2324 - 0.2e1 * t2364;
t2303 = 0.2e1 * t2363 + t2367;
t2299 = t2342 * pkin(7) + t2353;
t2288 = -t2339 * g(3) - t2346;
t2287 = -t2383 - t2384;
t2286 = t2312 * t2339 - t2315 * t2333;
t2285 = -t2311 * t2333 + t2316 * t2339;
t2284 = t2312 * t2333 + t2315 * t2339;
t2283 = t2311 * t2339 + t2316 * t2333;
t2277 = -t2366 - t2375;
t2276 = t2366 - t2375;
t2275 = -t2383 - t2385;
t2271 = -t2384 - t2385;
t2268 = -t2288 * t2333 + t2289 * t2339;
t2267 = t2288 * t2339 + t2289 * t2333;
t2259 = -t2386 - t2389;
t2257 = t2277 * t2338 - t2287 * t2332;
t2256 = t2277 * t2332 + t2287 * t2338;
t2255 = (qJD(3) - t2326) * t2296 + t2356;
t2253 = t2362 + t2374;
t2251 = -t2386 - t2390;
t2249 = t2275 * t2338 - t2276 * t2332;
t2248 = t2275 * t2332 + t2276 * t2338;
t2247 = -t2389 - t2390;
t2246 = pkin(5) * t2293 - pkin(11) * t2263;
t2242 = -t2387 - t2391;
t2241 = t2357 - t2378;
t2234 = t2281 * t2369 + t2349;
t2233 = t2245 - t2376;
t2232 = -t2282 * t2369 - t2358;
t2231 = (qJD(4) + t2295) * t2282 + t2358;
t2230 = -t2387 - t2392;
t2229 = -t2256 * t2333 + t2257 * t2339;
t2228 = t2256 * t2339 + t2257 * t2333;
t2227 = t2253 * t2338 - t2255 * t2332;
t2226 = t2253 * t2332 + t2255 * t2338;
t2224 = -t2248 * t2333 + t2249 * t2339;
t2223 = t2248 * t2339 + t2249 * t2333;
t2222 = -t2388 - t2393;
t2221 = t2355 - t2379;
t2219 = t2241 * t2337 - t2259 * t2331;
t2218 = t2241 * t2331 + t2259 * t2337;
t2214 = -t2391 - t2392;
t2213 = -t2240 * t2331 + t2251 * t2337;
t2212 = t2240 * t2337 + t2251 * t2331;
t2208 = -t2238 * t2332 + t2239 * t2338;
t2207 = t2238 * t2338 + t2239 * t2332;
t2203 = t2232 * t2337 - t2234 * t2331;
t2202 = t2232 * t2331 + t2234 * t2337;
t2201 = -t2388 - t2394;
t2200 = t2221 * t2336 - t2242 * t2330;
t2199 = t2221 * t2330 + t2242 * t2336;
t2198 = -t2226 * t2333 + t2227 * t2339;
t2197 = t2226 * t2339 + t2227 * t2333;
t2196 = t2350 - t2380;
t2195 = -t2350 - t2380;
t2194 = t2219 * t2338 + t2233 * t2332;
t2193 = t2219 * t2332 - t2233 * t2338;
t2192 = -t2220 * t2330 + t2230 * t2336;
t2191 = t2220 * t2336 + t2230 * t2330;
t2190 = t2213 * t2338 + t2231 * t2332;
t2189 = t2213 * t2332 - t2231 * t2338;
t2188 = t2261 * t2370 + t2344;
t2187 = t2206 - t2377;
t2186 = -t2263 * t2370 - t2360;
t2185 = (qJD(5) + t2293) * t2263 + t2360;
t2184 = t2203 * t2338 + t2247 * t2332;
t2183 = t2203 * t2332 - t2247 * t2338;
t2179 = -t2393 - t2394;
t2178 = t2196 * t2335 - t2222 * t2329;
t2177 = t2196 * t2329 + t2222 * t2335;
t2176 = -t2207 * t2333 + t2208 * t2339;
t2175 = t2207 * t2339 + t2208 * t2333;
t2172 = -t2199 * t2331 + t2200 * t2337;
t2171 = t2199 * t2337 + t2200 * t2331;
t2170 = -t2195 * t2329 + t2201 * t2335;
t2169 = t2195 * t2335 + t2201 * t2329;
t2165 = -t2193 * t2333 + t2194 * t2339;
t2164 = t2193 * t2339 + t2194 * t2333;
t2163 = -t2191 * t2331 + t2192 * t2337;
t2162 = t2191 * t2337 + t2192 * t2331;
t2161 = -t2189 * t2333 + t2190 * t2339;
t2160 = t2189 * t2339 + t2190 * t2333;
t2159 = t2186 * t2336 - t2188 * t2330;
t2158 = t2186 * t2330 + t2188 * t2336;
t2157 = -t2183 * t2333 + t2184 * t2339;
t2156 = t2183 * t2339 + t2184 * t2333;
t2155 = -t2180 * t2331 + t2181 * t2337;
t2154 = t2180 * t2337 + t2181 * t2331;
t2153 = t2235 * t2371 + t2345;
t2152 = -t2235 * t2395 - t2345;
t2151 = -t2237 * t2371 + t2361;
t2150 = t2237 * t2395 - t2361;
t2149 = pkin(5) * t2351 - pkin(11) * t2392 + t2246 * t2263 + t2182;
t2148 = t2172 * t2338 + t2187 * t2332;
t2147 = t2172 * t2332 - t2187 * t2338;
t2146 = t2155 * t2338 + t2216 * t2332;
t2145 = t2155 * t2332 - t2216 * t2338;
t2144 = -t2177 * t2330 + t2178 * t2336;
t2143 = t2177 * t2336 + t2178 * t2330;
t2142 = t2163 * t2338 + t2185 * t2332;
t2141 = t2163 * t2332 - t2185 * t2338;
t2138 = -t2169 * t2330 + t2170 * t2336;
t2137 = t2169 * t2336 + t2170 * t2330;
t2136 = -t2158 * t2331 + t2159 * t2337;
t2135 = t2158 * t2337 + t2159 * t2331;
t2134 = t2136 * t2338 + t2214 * t2332;
t2133 = t2136 * t2332 - t2214 * t2338;
t2132 = t2151 * t2335 - t2153 * t2329;
t2131 = t2151 * t2329 + t2153 * t2335;
t2130 = -pkin(5) * t2392 - pkin(11) * t2351 - t2246 * t2293 + t2140;
t2129 = (-t2206 - t2377) * pkin(11) + t2220 * pkin(5) + t2139;
t2128 = -t2147 * t2333 + t2148 * t2339;
t2127 = t2147 * t2339 + t2148 * t2333;
t2126 = -t2145 * t2333 + t2146 * t2339;
t2125 = t2145 * t2339 + t2146 * t2333;
t2124 = -t2143 * t2331 + t2144 * t2337;
t2123 = t2143 * t2337 + t2144 * t2331;
t2122 = -t2141 * t2333 + t2142 * t2339;
t2121 = t2141 * t2339 + t2142 * t2333;
t2120 = -t2139 * t2330 + t2140 * t2336;
t2119 = t2139 * t2336 + t2140 * t2330;
t2118 = -t2137 * t2331 + t2138 * t2337;
t2117 = t2137 * t2337 + t2138 * t2331;
t2116 = t2124 * t2338 + t2152 * t2332;
t2115 = t2124 * t2332 - t2152 * t2338;
t2114 = -t2133 * t2333 + t2134 * t2339;
t2113 = t2133 * t2339 + t2134 * t2333;
t2112 = -t2131 * t2330 + t2132 * t2336;
t2111 = t2131 * t2336 + t2132 * t2330;
t2110 = t2118 * t2338 + t2150 * t2332;
t2109 = t2118 * t2332 - t2150 * t2338;
t2108 = t2129 * t2329 + t2130 * t2335;
t2107 = t2129 * t2335 - t2130 * t2329;
t2106 = -t2119 * t2331 + t2120 * t2337;
t2105 = t2119 * t2337 + t2120 * t2331;
t2104 = t2106 * t2338 + t2182 * t2332;
t2103 = t2106 * t2332 - t2182 * t2338;
t2102 = -t2115 * t2333 + t2116 * t2339;
t2101 = t2115 * t2339 + t2116 * t2333;
t2100 = -t2111 * t2331 + t2112 * t2337;
t2099 = t2111 * t2337 + t2112 * t2331;
t2098 = -t2109 * t2333 + t2110 * t2339;
t2097 = t2109 * t2339 + t2110 * t2333;
t2096 = -t2107 * t2329 + t2108 * t2335;
t2095 = t2107 * t2335 + t2108 * t2329;
t2094 = t2100 * t2338 + t2179 * t2332;
t2093 = t2100 * t2332 - t2179 * t2338;
t2092 = -t2103 * t2333 + t2104 * t2339;
t2091 = t2103 * t2339 + t2104 * t2333;
t2090 = -t2095 * t2330 + t2096 * t2336;
t2089 = t2095 * t2336 + t2096 * t2330;
t2088 = -t2093 * t2333 + t2094 * t2339;
t2087 = t2093 * t2339 + t2094 * t2333;
t2086 = -t2089 * t2331 + t2090 * t2337;
t2085 = t2089 * t2337 + t2090 * t2331;
t2084 = t2086 * t2338 + t2149 * t2332;
t2083 = t2086 * t2332 - t2149 * t2338;
t2082 = -t2083 * t2333 + t2084 * t2339;
t2081 = t2083 * t2339 + t2084 * t2333;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2309, -t2308, 0, -t2313 * t2334 + t2314 * t2340, 0, 0, 0, 0, 0, 0, t2285 * t2340 - t2306 * t2334, t2286 * t2340 + t2303 * t2334, t2307 * t2340 - t2310 * t2334, t2268 * t2340 - t2299 * t2334, 0, 0, 0, 0, 0, 0, t2224 * t2340 + t2252 * t2334, t2229 * t2340 - t2334 * t2359, t2198 * t2340 + t2271 * t2334, t2176 * t2340 - t2272 * t2334, 0, 0, 0, 0, 0, 0, t2161 * t2340 + t2212 * t2334, t2165 * t2340 + t2218 * t2334, t2157 * t2340 + t2202 * t2334, t2126 * t2340 + t2154 * t2334, 0, 0, 0, 0, 0, 0, t2122 * t2340 + t2162 * t2334, t2128 * t2340 + t2171 * t2334, t2114 * t2340 + t2135 * t2334, t2092 * t2340 + t2105 * t2334, 0, 0, 0, 0, 0, 0, t2098 * t2340 + t2117 * t2334, t2102 * t2340 + t2123 * t2334, t2088 * t2340 + t2099 * t2334, t2082 * t2340 + t2085 * t2334; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2308, t2309, 0, t2313 * t2340 + t2314 * t2334, 0, 0, 0, 0, 0, 0, t2285 * t2334 + t2306 * t2340, t2286 * t2334 - t2303 * t2340, t2307 * t2334 + t2310 * t2340, t2268 * t2334 + t2299 * t2340, 0, 0, 0, 0, 0, 0, t2224 * t2334 - t2252 * t2340, t2229 * t2334 + t2340 * t2359, t2198 * t2334 - t2271 * t2340, t2176 * t2334 + t2272 * t2340, 0, 0, 0, 0, 0, 0, t2161 * t2334 - t2212 * t2340, t2165 * t2334 - t2218 * t2340, t2157 * t2334 - t2202 * t2340, t2126 * t2334 - t2154 * t2340, 0, 0, 0, 0, 0, 0, t2122 * t2334 - t2162 * t2340, t2128 * t2334 - t2171 * t2340, t2114 * t2334 - t2135 * t2340, t2092 * t2334 - t2105 * t2340, 0, 0, 0, 0, 0, 0, t2098 * t2334 - t2117 * t2340, t2102 * t2334 - t2123 * t2340, t2088 * t2334 - t2099 * t2340, t2082 * t2334 - t2085 * t2340; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2283, t2284, 0, t2267, 0, 0, 0, 0, 0, 0, t2223, t2228, t2197, t2175, 0, 0, 0, 0, 0, 0, t2160, t2164, t2156, t2125, 0, 0, 0, 0, 0, 0, t2121, t2127, t2113, t2091, 0, 0, 0, 0, 0, 0, t2097, t2101, t2087, t2081; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2342, -qJDD(1), 0, t2314, 0, 0, 0, 0, 0, 0, t2285, t2286, t2307, t2268, 0, 0, 0, 0, 0, 0, t2224, t2229, t2198, t2176, 0, 0, 0, 0, 0, 0, t2161, t2165, t2157, t2126, 0, 0, 0, 0, 0, 0, t2122, t2128, t2114, t2092, 0, 0, 0, 0, 0, 0, t2098, t2102, t2088, t2082; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2342, 0, t2313, 0, 0, 0, 0, 0, 0, t2306, -t2303, t2310, t2299, 0, 0, 0, 0, 0, 0, -t2252, t2359, -t2271, t2272, 0, 0, 0, 0, 0, 0, -t2212, -t2218, -t2202, -t2154, 0, 0, 0, 0, 0, 0, -t2162, -t2171, -t2135, -t2105, 0, 0, 0, 0, 0, 0, -t2117, -t2123, -t2099, -t2085; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2283, t2284, 0, t2267, 0, 0, 0, 0, 0, 0, t2223, t2228, t2197, t2175, 0, 0, 0, 0, 0, 0, t2160, t2164, t2156, t2125, 0, 0, 0, 0, 0, 0, t2121, t2127, t2113, t2091, 0, 0, 0, 0, 0, 0, t2097, t2101, t2087, t2081; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2316, t2312, t2324, t2289, 0, 0, 0, 0, 0, 0, t2249, t2257, t2227, t2208, 0, 0, 0, 0, 0, 0, t2190, t2194, t2184, t2146, 0, 0, 0, 0, 0, 0, t2142, t2148, t2134, t2104, 0, 0, 0, 0, 0, 0, t2110, t2116, t2094, t2084; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2311, t2315, -t2367, t2288, 0, 0, 0, 0, 0, 0, t2248, t2256, t2226, t2207, 0, 0, 0, 0, 0, 0, t2189, t2193, t2183, t2145, 0, 0, 0, 0, 0, 0, t2141, t2147, t2133, t2103, 0, 0, 0, 0, 0, 0, t2109, t2115, t2093, t2083; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2306, t2303, -t2310, -t2299, 0, 0, 0, 0, 0, 0, t2252, -t2359, t2271, -t2272, 0, 0, 0, 0, 0, 0, t2212, t2218, t2202, t2154, 0, 0, 0, 0, 0, 0, t2162, t2171, t2135, t2105, 0, 0, 0, 0, 0, 0, t2117, t2123, t2099, t2085; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2275, t2277, t2253, t2239, 0, 0, 0, 0, 0, 0, t2213, t2219, t2203, t2155, 0, 0, 0, 0, 0, 0, t2163, t2172, t2136, t2106, 0, 0, 0, 0, 0, 0, t2118, t2124, t2100, t2086; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2276, t2287, t2255, t2238, 0, 0, 0, 0, 0, 0, -t2231, -t2233, -t2247, -t2216, 0, 0, 0, 0, 0, 0, -t2185, -t2187, -t2214, -t2182, 0, 0, 0, 0, 0, 0, -t2150, -t2152, -t2179, -t2149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2252, -t2359, t2271, -t2272, 0, 0, 0, 0, 0, 0, t2212, t2218, t2202, t2154, 0, 0, 0, 0, 0, 0, t2162, t2171, t2135, t2105, 0, 0, 0, 0, 0, 0, t2117, t2123, t2099, t2085; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2251, t2241, t2232, t2181, 0, 0, 0, 0, 0, 0, t2192, t2200, t2159, t2120, 0, 0, 0, 0, 0, 0, t2138, t2144, t2112, t2090; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2240, t2259, t2234, t2180, 0, 0, 0, 0, 0, 0, t2191, t2199, t2158, t2119, 0, 0, 0, 0, 0, 0, t2137, t2143, t2111, t2089; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2231, t2233, t2247, t2216, 0, 0, 0, 0, 0, 0, t2185, t2187, t2214, t2182, 0, 0, 0, 0, 0, 0, t2150, t2152, t2179, t2149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2230, t2221, t2186, t2140, 0, 0, 0, 0, 0, 0, t2170, t2178, t2132, t2096; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2220, t2242, t2188, t2139, 0, 0, 0, 0, 0, 0, t2169, t2177, t2131, t2095; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2185, t2187, t2214, t2182, 0, 0, 0, 0, 0, 0, t2150, t2152, t2179, t2149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2201, t2196, t2151, t2108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2195, t2222, t2153, t2107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2150, t2152, t2179, t2149;];
f_new_reg  = t1;