% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
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
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 21:57
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RRRRPR8_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR8_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR8_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR8_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR8_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR8_invdynf_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 21:57:12
% EndTime: 2019-05-07 21:57:22
% DurationCPUTime: 10.15s
% Computational Cost: add. (45808->330), mult. (92126->390), div. (0->0), fcn. (66620->10), ass. (0->228)
t2386 = sin(qJ(2));
t2428 = qJD(1) * t2386;
t2375 = qJD(2) * t2428;
t2391 = cos(qJ(2));
t2414 = t2391 * qJDD(1);
t2409 = t2375 - t2414;
t2407 = -qJDD(3) - t2409;
t2405 = -qJDD(4) + t2407;
t2385 = sin(qJ(3));
t2390 = cos(qJ(3));
t2350 = -t2390 * qJD(2) + t2385 * t2428;
t2351 = qJD(2) * t2385 + t2390 * t2428;
t2384 = sin(qJ(4));
t2389 = cos(qJ(4));
t2329 = t2350 * t2389 + t2384 * t2351;
t2331 = -t2350 * t2384 + t2351 * t2389;
t2423 = t2329 * t2331;
t2284 = t2405 - t2423;
t2328 = t2331 ^ 2;
t2427 = qJD(1) * t2391;
t2372 = -qJD(3) + t2427;
t2368 = qJD(4) - t2372;
t2431 = t2368 ^ 2;
t2442 = -t2328 - t2431;
t2264 = t2284 * t2384 + t2389 * t2442;
t2266 = t2284 * t2389 - t2384 * t2442;
t2222 = t2264 * t2385 - t2266 * t2390;
t2413 = qJD(2) * t2427;
t2415 = t2386 * qJDD(1);
t2354 = t2413 + t2415;
t2406 = -t2385 * qJDD(2) - t2390 * t2354;
t2325 = -qJD(3) * t2350 - t2406;
t2410 = -t2390 * qJDD(2) + t2385 * t2354;
t2404 = -qJD(3) * t2351 - t2410;
t2400 = -t2389 * t2325 - t2384 * t2404;
t2399 = -t2329 * qJD(4) - t2400;
t2422 = t2329 * t2368;
t2398 = t2399 - t2422;
t2207 = t2222 * t2391 - t2386 * t2398;
t2230 = t2264 * t2390 + t2266 * t2385;
t2387 = sin(qJ(1));
t2392 = cos(qJ(1));
t2465 = t2207 * t2387 + t2230 * t2392;
t2464 = t2207 * t2392 - t2230 * t2387;
t2205 = t2222 * t2386 + t2391 * t2398;
t2285 = t2405 + t2423;
t2435 = t2329 ^ 2;
t2441 = -t2431 - t2435;
t2448 = t2285 * t2384 + t2389 * t2441;
t2449 = -t2389 * t2285 + t2384 * t2441;
t2451 = -t2385 * t2449 + t2390 * t2448;
t2461 = t2386 * t2451;
t2450 = t2385 * t2448 + t2390 * t2449;
t2460 = t2387 * t2450;
t2459 = t2391 * t2451;
t2458 = t2392 * t2450;
t2259 = t2399 + t2422;
t2411 = t2384 * t2325 - t2389 * t2404;
t2418 = qJD(4) - t2368;
t2401 = -t2331 * t2418 - t2411;
t2439 = t2259 * t2384 + t2389 * t2401;
t2440 = -t2389 * t2259 + t2384 * t2401;
t2446 = t2385 * t2439 + t2390 * t2440;
t2280 = t2328 + t2435;
t2447 = -t2385 * t2440 + t2390 * t2439;
t2452 = -t2280 * t2386 + t2391 * t2447;
t2457 = t2387 * t2452 - t2392 * t2446;
t2456 = t2387 * t2446 + t2392 * t2452;
t2453 = t2280 * t2391 + t2386 * t2447;
t2362 = qJD(6) - t2368;
t2443 = qJD(6) + t2362;
t2438 = qJD(2) ^ 2;
t2383 = sin(qJ(6));
t2388 = cos(qJ(6));
t2294 = -t2388 * t2329 + t2331 * t2383;
t2437 = t2294 ^ 2;
t2296 = t2329 * t2383 + t2331 * t2388;
t2436 = t2296 ^ 2;
t2434 = t2350 ^ 2;
t2433 = t2351 ^ 2;
t2432 = t2362 ^ 2;
t2430 = t2372 ^ 2;
t2429 = t2391 * g(3);
t2424 = t2294 * t2296;
t2421 = t2350 * t2351;
t2420 = t2350 * t2372;
t2419 = qJD(3) + t2372;
t2417 = qJD(6) - t2362;
t2365 = t2387 * g(1) - t2392 * g(2);
t2393 = qJD(1) ^ 2;
t2346 = qJDD(1) * pkin(1) + t2393 * pkin(7) + t2365;
t2301 = (-t2354 - t2413) * pkin(8) + (t2409 + t2375) * pkin(2) - t2346;
t2366 = -g(1) * t2392 - g(2) * t2387;
t2347 = -pkin(1) * t2393 + qJDD(1) * pkin(7) + t2366;
t2338 = -g(3) * t2386 + t2391 * t2347;
t2352 = (-pkin(2) * t2391 - pkin(8) * t2386) * qJD(1);
t2312 = -pkin(2) * t2438 + qJDD(2) * pkin(8) + t2352 * t2427 + t2338;
t2277 = t2385 * t2301 + t2390 * t2312;
t2339 = -pkin(3) * t2372 - pkin(9) * t2351;
t2244 = -pkin(3) * t2434 + pkin(9) * t2404 + t2372 * t2339 + t2277;
t2276 = t2390 * t2301 - t2385 * t2312;
t2320 = -t2407 - t2421;
t2395 = (-t2325 + t2420) * pkin(9) + t2320 * pkin(3) + t2276;
t2215 = t2389 * t2244 + t2384 * t2395;
t2379 = t2386 ^ 2;
t2380 = t2391 ^ 2;
t2416 = t2379 + t2380;
t2412 = pkin(4) * t2368 - (2 * qJD(5));
t2214 = -t2384 * t2244 + t2389 * t2395;
t2408 = -pkin(5) * t2368 - pkin(10) * t2331;
t2403 = qJD(4) * t2331 + t2411;
t2402 = -qJDD(6) - t2405;
t2297 = pkin(4) * t2329 - qJ(5) * t2331;
t2201 = -pkin(4) * t2431 - qJ(5) * t2405 + 0.2e1 * qJD(5) * t2368 - t2329 * t2297 + t2215;
t2202 = pkin(4) * t2405 - qJ(5) * t2431 + t2331 * t2297 + qJDD(5) - t2214;
t2311 = t2429 - qJDD(2) * pkin(2) - t2438 * pkin(8) + (qJD(1) * t2352 + t2347) * t2386;
t2270 = -t2404 * pkin(3) - t2434 * pkin(9) + t2351 * t2339 + t2311;
t2397 = -t2383 * t2399 + t2388 * t2403;
t2396 = t2403 * pkin(4) - qJ(5) * t2398 + t2270;
t2394 = -t2383 * t2403 - t2388 * t2399;
t2371 = t2391 * t2393 * t2386;
t2370 = -t2380 * t2393 - t2438;
t2369 = -t2379 * t2393 - t2438;
t2364 = -qJDD(2) + t2371;
t2363 = qJDD(2) + t2371;
t2360 = t2416 * t2393;
t2359 = -qJDD(1) * t2387 - t2392 * t2393;
t2358 = qJDD(1) * t2392 - t2387 * t2393;
t2357 = t2416 * qJDD(1);
t2355 = -0.2e1 * t2375 + t2414;
t2353 = 0.2e1 * t2413 + t2415;
t2337 = -t2386 * t2347 - t2429;
t2336 = t2364 * t2391 - t2369 * t2386;
t2335 = -t2363 * t2386 + t2370 * t2391;
t2334 = t2364 * t2386 + t2369 * t2391;
t2333 = t2363 * t2391 + t2370 * t2386;
t2332 = -t2430 - t2433;
t2326 = -t2430 - t2434;
t2319 = t2407 - t2421;
t2318 = -t2433 - t2434;
t2307 = -t2337 * t2386 + t2338 * t2391;
t2306 = t2337 * t2391 + t2338 * t2386;
t2305 = t2350 * t2419 + t2406;
t2304 = t2325 + t2420;
t2303 = -t2351 * t2419 - t2410;
t2302 = (qJD(3) - t2372) * t2351 + t2410;
t2293 = t2319 * t2390 - t2332 * t2385;
t2292 = t2319 * t2385 + t2332 * t2390;
t2288 = -t2320 * t2385 + t2326 * t2390;
t2287 = t2320 * t2390 + t2326 * t2385;
t2281 = -t2432 - t2436;
t2272 = t2303 * t2390 - t2305 * t2385;
t2271 = t2303 * t2385 + t2305 * t2390;
t2269 = t2293 * t2391 + t2304 * t2386;
t2268 = t2293 * t2386 - t2304 * t2391;
t2263 = t2288 * t2391 + t2302 * t2386;
t2262 = t2288 * t2386 - t2302 * t2391;
t2256 = t2368 * t2331 + t2403;
t2255 = (qJD(4) + t2368) * t2331 + t2411;
t2254 = -t2432 - t2437;
t2252 = t2402 - t2424;
t2251 = -t2402 - t2424;
t2246 = t2272 * t2391 + t2318 * t2386;
t2245 = t2272 * t2386 - t2318 * t2391;
t2240 = -t2436 - t2437;
t2239 = -t2276 * t2385 + t2277 * t2390;
t2238 = t2276 * t2390 + t2277 * t2385;
t2237 = t2252 * t2388 - t2281 * t2383;
t2236 = t2252 * t2383 + t2281 * t2388;
t2235 = t2239 * t2391 + t2311 * t2386;
t2234 = t2239 * t2386 - t2311 * t2391;
t2225 = -t2251 * t2383 + t2254 * t2388;
t2224 = t2251 * t2388 + t2254 * t2383;
t2219 = t2294 * t2417 + t2394;
t2218 = -t2294 * t2443 - t2394;
t2217 = -t2296 * t2417 + t2397;
t2216 = t2296 * t2443 - t2397;
t2213 = t2256 * t2386 + t2459;
t2211 = -t2256 * t2391 + t2461;
t2209 = t2412 * t2331 + t2396;
t2208 = t2255 * t2386 + t2459;
t2206 = -t2255 * t2391 + t2461;
t2204 = t2236 * t2384 + t2237 * t2389;
t2203 = -t2236 * t2389 + t2237 * t2384;
t2196 = t2224 * t2384 + t2225 * t2389;
t2195 = -t2224 * t2389 + t2225 * t2384;
t2194 = t2411 * pkin(5) + t2435 * pkin(10) + t2396 + (pkin(5) * qJD(4) - t2408 + t2412) * t2331;
t2189 = -pkin(5) * t2435 + pkin(10) * t2403 + t2368 * t2408 + t2201;
t2188 = (t2329 * t2418 + t2400) * pkin(10) + t2285 * pkin(5) + t2202;
t2187 = t2217 * t2388 - t2219 * t2383;
t2186 = t2217 * t2383 + t2219 * t2388;
t2185 = -t2214 * t2384 + t2215 * t2389;
t2184 = t2214 * t2389 + t2215 * t2384;
t2183 = -t2203 * t2385 + t2204 * t2390;
t2182 = t2203 * t2390 + t2204 * t2385;
t2181 = t2201 * t2389 + t2202 * t2384;
t2180 = t2201 * t2384 - t2202 * t2389;
t2179 = -t2195 * t2385 + t2196 * t2390;
t2178 = t2195 * t2390 + t2196 * t2385;
t2177 = t2183 * t2391 - t2218 * t2386;
t2176 = t2183 * t2386 + t2218 * t2391;
t2175 = t2188 * t2383 + t2189 * t2388;
t2174 = t2188 * t2388 - t2189 * t2383;
t2173 = t2186 * t2384 + t2187 * t2389;
t2172 = -t2186 * t2389 + t2187 * t2384;
t2171 = -t2184 * t2385 + t2185 * t2390;
t2170 = t2184 * t2390 + t2185 * t2385;
t2169 = t2179 * t2391 - t2216 * t2386;
t2168 = t2179 * t2386 + t2216 * t2391;
t2167 = t2171 * t2391 + t2270 * t2386;
t2166 = t2171 * t2386 - t2270 * t2391;
t2165 = -t2180 * t2385 + t2181 * t2390;
t2164 = t2180 * t2390 + t2181 * t2385;
t2163 = t2165 * t2391 + t2209 * t2386;
t2162 = t2165 * t2386 - t2209 * t2391;
t2161 = -t2174 * t2383 + t2175 * t2388;
t2160 = t2174 * t2388 + t2175 * t2383;
t2159 = -t2172 * t2385 + t2173 * t2390;
t2158 = t2172 * t2390 + t2173 * t2385;
t2157 = t2159 * t2391 - t2240 * t2386;
t2156 = t2159 * t2386 + t2240 * t2391;
t2155 = t2160 * t2384 + t2161 * t2389;
t2154 = -t2160 * t2389 + t2161 * t2384;
t2153 = -t2154 * t2385 + t2155 * t2390;
t2152 = t2154 * t2390 + t2155 * t2385;
t2151 = t2153 * t2391 + t2194 * t2386;
t2150 = t2153 * t2386 - t2194 * t2391;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2359, -t2358, 0, -t2365 * t2387 + t2366 * t2392, 0, 0, 0, 0, 0, 0, t2335 * t2392 - t2355 * t2387, t2336 * t2392 + t2353 * t2387, t2357 * t2392 - t2360 * t2387, t2307 * t2392 - t2346 * t2387, 0, 0, 0, 0, 0, 0, t2263 * t2392 + t2287 * t2387, t2269 * t2392 + t2292 * t2387, t2246 * t2392 + t2271 * t2387, t2235 * t2392 + t2238 * t2387, 0, 0, 0, 0, 0, 0, t2208 * t2392 + t2460, -t2464, t2456, t2167 * t2392 + t2170 * t2387, 0, 0, 0, 0, 0, 0, t2213 * t2392 + t2460, t2456, t2464, t2163 * t2392 + t2164 * t2387, 0, 0, 0, 0, 0, 0, t2169 * t2392 + t2178 * t2387, t2177 * t2392 + t2182 * t2387, t2157 * t2392 + t2158 * t2387, t2151 * t2392 + t2152 * t2387; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2358, t2359, 0, t2365 * t2392 + t2366 * t2387, 0, 0, 0, 0, 0, 0, t2335 * t2387 + t2355 * t2392, t2336 * t2387 - t2353 * t2392, t2357 * t2387 + t2360 * t2392, t2307 * t2387 + t2346 * t2392, 0, 0, 0, 0, 0, 0, t2263 * t2387 - t2287 * t2392, t2269 * t2387 - t2292 * t2392, t2246 * t2387 - t2271 * t2392, t2235 * t2387 - t2238 * t2392, 0, 0, 0, 0, 0, 0, t2208 * t2387 - t2458, -t2465, t2457, t2167 * t2387 - t2170 * t2392, 0, 0, 0, 0, 0, 0, t2213 * t2387 - t2458, t2457, t2465, t2163 * t2387 - t2164 * t2392, 0, 0, 0, 0, 0, 0, t2169 * t2387 - t2178 * t2392, t2177 * t2387 - t2182 * t2392, t2157 * t2387 - t2158 * t2392, t2151 * t2387 - t2152 * t2392; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2333, t2334, 0, t2306, 0, 0, 0, 0, 0, 0, t2262, t2268, t2245, t2234, 0, 0, 0, 0, 0, 0, t2206, -t2205, t2453, t2166, 0, 0, 0, 0, 0, 0, t2211, t2453, t2205, t2162, 0, 0, 0, 0, 0, 0, t2168, t2176, t2156, t2150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2393, -qJDD(1), 0, t2366, 0, 0, 0, 0, 0, 0, t2335, t2336, t2357, t2307, 0, 0, 0, 0, 0, 0, t2263, t2269, t2246, t2235, 0, 0, 0, 0, 0, 0, t2208, -t2207, t2452, t2167, 0, 0, 0, 0, 0, 0, t2213, t2452, t2207, t2163, 0, 0, 0, 0, 0, 0, t2169, t2177, t2157, t2151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2393, 0, t2365, 0, 0, 0, 0, 0, 0, t2355, -t2353, t2360, t2346, 0, 0, 0, 0, 0, 0, -t2287, -t2292, -t2271, -t2238, 0, 0, 0, 0, 0, 0, -t2450, -t2230, -t2446, -t2170, 0, 0, 0, 0, 0, 0, -t2450, -t2446, t2230, -t2164, 0, 0, 0, 0, 0, 0, -t2178, -t2182, -t2158, -t2152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2333, t2334, 0, t2306, 0, 0, 0, 0, 0, 0, t2262, t2268, t2245, t2234, 0, 0, 0, 0, 0, 0, t2206, -t2205, t2453, t2166, 0, 0, 0, 0, 0, 0, t2211, t2453, t2205, t2162, 0, 0, 0, 0, 0, 0, t2168, t2176, t2156, t2150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2370, t2364, t2414, t2338, 0, 0, 0, 0, 0, 0, t2288, t2293, t2272, t2239, 0, 0, 0, 0, 0, 0, t2451, -t2222, t2447, t2171, 0, 0, 0, 0, 0, 0, t2451, t2447, t2222, t2165, 0, 0, 0, 0, 0, 0, t2179, t2183, t2159, t2153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2363, t2369, -t2415, t2337, 0, 0, 0, 0, 0, 0, -t2302, -t2304, -t2318, -t2311, 0, 0, 0, 0, 0, 0, -t2255, -t2398, t2280, -t2270, 0, 0, 0, 0, 0, 0, -t2256, t2280, t2398, -t2209, 0, 0, 0, 0, 0, 0, t2216, t2218, t2240, -t2194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2355, t2353, -t2360, -t2346, 0, 0, 0, 0, 0, 0, t2287, t2292, t2271, t2238, 0, 0, 0, 0, 0, 0, t2450, t2230, t2446, t2170, 0, 0, 0, 0, 0, 0, t2450, t2446, -t2230, t2164, 0, 0, 0, 0, 0, 0, t2178, t2182, t2158, t2152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2326, t2319, t2303, t2277, 0, 0, 0, 0, 0, 0, t2448, t2266, t2439, t2185, 0, 0, 0, 0, 0, 0, t2448, t2439, -t2266, t2181, 0, 0, 0, 0, 0, 0, t2196, t2204, t2173, t2155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2320, t2332, t2305, t2276, 0, 0, 0, 0, 0, 0, t2449, t2264, t2440, t2184, 0, 0, 0, 0, 0, 0, t2449, t2440, -t2264, t2180, 0, 0, 0, 0, 0, 0, t2195, t2203, t2172, t2154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2302, t2304, t2318, t2311, 0, 0, 0, 0, 0, 0, t2255, t2398, -t2280, t2270, 0, 0, 0, 0, 0, 0, t2256, -t2280, -t2398, t2209, 0, 0, 0, 0, 0, 0, -t2216, -t2218, -t2240, t2194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2441, t2284, t2401, t2215, 0, 0, 0, 0, 0, 0, t2441, t2401, -t2284, t2201, 0, 0, 0, 0, 0, 0, t2225, t2237, t2187, t2161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2285, t2442, -t2259, t2214, 0, 0, 0, 0, 0, 0, -t2285, -t2259, -t2442, -t2202, 0, 0, 0, 0, 0, 0, -t2224, -t2236, -t2186, -t2160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2255, t2398, -t2280, t2270, 0, 0, 0, 0, 0, 0, t2256, -t2280, -t2398, t2209, 0, 0, 0, 0, 0, 0, -t2216, -t2218, -t2240, t2194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2441, t2401, -t2284, t2201, 0, 0, 0, 0, 0, 0, t2225, t2237, t2187, t2161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2256, -t2280, -t2398, t2209, 0, 0, 0, 0, 0, 0, -t2216, -t2218, -t2240, t2194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2285, t2259, t2442, t2202, 0, 0, 0, 0, 0, 0, t2224, t2236, t2186, t2160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2254, t2252, t2217, t2175; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2251, t2281, t2219, t2174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2216, t2218, t2240, -t2194;];
f_new_reg  = t1;