% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 17:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RRPRPR14_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR14_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR14_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_invdynf_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 17:02:05
% EndTime: 2019-05-06 17:02:15
% DurationCPUTime: 9.95s
% Computational Cost: add. (23159->342), mult. (52910->394), div. (0->0), fcn. (38682->10), ass. (0->235)
t2364 = sin(pkin(6));
t2368 = sin(qJ(2));
t2372 = cos(qJ(2));
t2329 = (qJD(1) * qJD(2) * t2372 + qJDD(1) * t2368) * t2364;
t2321 = qJDD(4) + t2329;
t2365 = cos(pkin(6));
t2359 = qJD(1) * t2365 + qJD(2);
t2367 = sin(qJ(4));
t2371 = cos(qJ(4));
t2419 = t2364 * t2372;
t2409 = qJD(1) * t2419;
t2314 = t2359 * t2367 + t2371 * t2409;
t2316 = t2359 * t2371 - t2367 * t2409;
t2427 = t2316 * t2314;
t2280 = t2321 + t2427;
t2420 = t2364 * t2368;
t2410 = qJD(1) * t2420;
t2346 = qJD(4) + t2410;
t2340 = t2346 ^ 2;
t2436 = t2316 ^ 2;
t2408 = -t2340 - t2436;
t2252 = t2280 * t2371 + t2367 * t2408;
t2250 = t2280 * t2367 - t2371 * t2408;
t2411 = qJDD(1) * t2364;
t2330 = -qJD(2) * t2410 + t2372 * t2411;
t2358 = qJDD(1) * t2365 + qJDD(2);
t2384 = t2367 * t2330 - t2371 * t2358;
t2415 = qJD(4) + t2346;
t2379 = -t2314 * t2415 - t2384;
t2391 = t2250 * t2372 + t2368 * t2379;
t2209 = t2364 * t2252 + t2365 * t2391;
t2225 = t2250 * t2368 - t2372 * t2379;
t2369 = sin(qJ(1));
t2373 = cos(qJ(1));
t2468 = t2209 * t2369 + t2225 * t2373;
t2467 = t2209 * t2373 - t2225 * t2369;
t2207 = -t2365 * t2252 + t2364 * t2391;
t2313 = t2314 ^ 2;
t2287 = -t2340 - t2313;
t2405 = -t2321 + t2427;
t2245 = -t2287 * t2367 + t2371 * t2405;
t2462 = t2245 * t2368;
t2461 = t2245 * t2372;
t2416 = qJD(4) - t2346;
t2270 = t2314 * t2416 + t2384;
t2460 = t2270 * t2367;
t2459 = t2270 * t2371;
t2246 = t2287 * t2371 + t2367 * t2405;
t2458 = t2364 * t2246;
t2457 = t2365 * t2246;
t2403 = t2359 * t2409;
t2383 = -t2329 - t2403;
t2355 = t2359 ^ 2;
t2374 = qJD(1) ^ 2;
t2423 = t2364 ^ 2 * t2374;
t2356 = t2368 ^ 2 * t2423;
t2311 = -t2356 - t2355;
t2344 = t2372 * t2368 * t2423;
t2328 = -t2344 + t2358;
t2386 = t2311 * t2372 - t2328 * t2368;
t2263 = t2364 * t2383 + t2365 * t2386;
t2290 = t2311 * t2368 + t2328 * t2372;
t2456 = t2263 * t2369 + t2290 * t2373;
t2455 = t2263 * t2373 - t2290 * t2369;
t2335 = t2359 * t2410;
t2302 = t2330 - t2335;
t2357 = t2372 ^ 2 * t2423;
t2331 = -t2357 - t2355;
t2440 = -t2358 - t2344;
t2385 = t2331 * t2368 - t2372 * t2440;
t2275 = t2364 * t2302 + t2365 * t2385;
t2297 = -t2331 * t2372 - t2368 * t2440;
t2454 = t2275 * t2369 + t2297 * t2373;
t2453 = t2275 * t2373 - t2297 * t2369;
t2442 = -t2313 - t2436;
t2448 = t2368 * t2442;
t2447 = t2372 * t2442;
t2404 = t2371 * t2330 + t2367 * t2358;
t2285 = qJD(4) * t2316 + t2404;
t2306 = t2346 * t2316;
t2265 = t2285 + t2306;
t2441 = pkin(2) * t2335 - 0.2e1 * qJD(3) * t2410;
t2273 = -t2365 * t2302 + t2364 * t2385;
t2261 = t2364 * t2386 - t2365 * t2383;
t2366 = sin(qJ(6));
t2370 = cos(qJ(6));
t2298 = -t2370 * t2314 + t2346 * t2366;
t2439 = t2298 ^ 2;
t2300 = t2314 * t2366 + t2346 * t2370;
t2438 = t2300 ^ 2;
t2312 = qJD(6) + t2316;
t2437 = t2312 ^ 2;
t2435 = 0.2e1 * qJD(3);
t2434 = 2 * qJD(5);
t2433 = g(3) * t2372;
t2432 = t2365 * g(3);
t2429 = t2298 * t2300;
t2428 = t2314 * t2346;
t2348 = t2369 * g(1) - g(2) * t2373;
t2323 = pkin(8) * t2364 * t2374 + qJDD(1) * pkin(1) + t2348;
t2426 = t2323 * t2365;
t2325 = (-pkin(2) * t2372 - qJ(3) * t2368) * t2364 * qJD(1);
t2425 = t2325 * t2368;
t2424 = t2359 * t2372;
t2414 = pkin(5) * t2316 - pkin(10) * t2346 + t2434;
t2413 = qJD(6) - t2312;
t2412 = qJD(6) + t2312;
t2349 = -g(1) * t2373 - g(2) * t2369;
t2324 = -pkin(1) * t2374 + pkin(8) * t2411 + t2349;
t2289 = -g(3) * t2420 + t2372 * t2324 + t2368 * t2426;
t2326 = pkin(3) * t2410 - pkin(9) * t2359;
t2237 = -pkin(3) * t2357 - t2432 - t2329 * qJ(3) + (-pkin(2) - pkin(9)) * t2330 + (-t2323 + (-qJ(3) * t2424 - t2326 * t2368) * qJD(1)) * t2364 + t2441;
t2406 = t2368 * t2324 - t2372 * t2426;
t2381 = -t2358 * pkin(2) - t2355 * qJ(3) + qJDD(3) + t2406;
t2239 = t2329 * pkin(3) + t2440 * pkin(9) + (t2433 + (-pkin(3) * t2424 + t2425) * qJD(1)) * t2364 + t2381;
t2213 = -t2367 * t2237 + t2371 * t2239;
t2407 = t2370 * t2285 - t2366 * t2321;
t2307 = -t2364 * t2323 - t2432;
t2292 = pkin(4) * t2314 - qJ(5) * t2316;
t2197 = -t2321 * pkin(4) - t2340 * qJ(5) + t2316 * t2292 + qJDD(5) - t2213;
t2286 = -qJD(4) * t2314 - t2384;
t2185 = t2405 * pkin(10) + (t2286 + t2428) * pkin(5) + t2197;
t2377 = -t2355 * pkin(2) + t2358 * qJ(3) + t2325 * t2409 + t2289;
t2236 = t2330 * pkin(3) - pkin(9) * t2357 + (t2435 + t2326) * t2359 + t2377;
t2375 = (-t2286 + t2428) * qJ(5) + t2236;
t2187 = -t2313 * pkin(5) + (pkin(4) + pkin(10)) * t2285 + (pkin(4) * t2346 - t2414) * t2316 + t2375;
t2173 = t2185 * t2370 - t2187 * t2366;
t2174 = t2185 * t2366 + t2187 * t2370;
t2162 = t2173 * t2370 + t2174 * t2366;
t2214 = t2371 * t2237 + t2367 * t2239;
t2378 = -t2340 * pkin(4) + t2321 * qJ(5) - t2314 * t2292 + t2214;
t2186 = -t2285 * pkin(5) - t2313 * pkin(10) + t2346 * t2414 + t2378;
t2160 = -t2162 * t2371 + t2186 * t2367;
t2163 = -t2173 * t2366 + t2174 * t2370;
t2402 = -t2160 * t2372 + t2163 * t2368;
t2196 = t2346 * t2434 + t2378;
t2178 = t2196 * t2367 - t2197 * t2371;
t2202 = pkin(4) * t2265 - 0.2e1 * qJD(5) * t2316 + t2375;
t2401 = -t2178 * t2372 + t2202 * t2368;
t2183 = t2213 * t2371 + t2214 * t2367;
t2400 = -t2183 * t2372 + t2236 * t2368;
t2227 = -t2300 * t2413 + t2407;
t2389 = -t2366 * t2285 - t2370 * t2321;
t2229 = t2298 * t2413 + t2389;
t2200 = t2227 * t2366 + t2229 * t2370;
t2248 = -t2438 - t2439;
t2188 = -t2200 * t2371 + t2248 * t2367;
t2201 = t2227 * t2370 - t2229 * t2366;
t2399 = -t2188 * t2372 + t2201 * t2368;
t2380 = qJDD(6) + t2286;
t2242 = t2380 - t2429;
t2256 = -t2437 - t2439;
t2217 = t2242 * t2370 + t2256 * t2366;
t2226 = t2300 * t2412 - t2407;
t2194 = -t2217 * t2371 + t2226 * t2367;
t2218 = -t2242 * t2366 + t2256 * t2370;
t2398 = -t2194 * t2372 + t2218 * t2368;
t2243 = -t2380 - t2429;
t2272 = -t2437 - t2438;
t2219 = t2243 * t2366 + t2272 * t2370;
t2228 = -t2298 * t2412 - t2389;
t2198 = -t2219 * t2371 + t2228 * t2367;
t2220 = t2243 * t2370 - t2272 * t2366;
t2397 = -t2198 * t2372 + t2220 * t2368;
t2266 = -t2285 + t2306;
t2230 = t2266 * t2367 + t2459;
t2396 = -t2230 * t2372 + t2448;
t2271 = -t2316 * t2416 - t2404;
t2231 = t2271 * t2367 + t2459;
t2395 = -t2231 * t2372 + t2448;
t2264 = t2316 * t2415 + t2404;
t2394 = t2264 * t2368 + t2461;
t2393 = -t2265 * t2368 - t2461;
t2254 = t2359 * t2435 + t2377;
t2257 = (qJD(1) * t2425 + t2433) * t2364 + t2381;
t2390 = t2254 * t2368 - t2257 * t2372;
t2288 = -g(3) * t2419 - t2406;
t2388 = t2288 * t2372 + t2289 * t2368;
t2304 = t2329 - t2403;
t2305 = t2330 + t2335;
t2387 = -t2304 * t2372 + t2305 * t2368;
t2342 = -qJDD(1) * t2369 - t2373 * t2374;
t2341 = qJDD(1) * t2373 - t2369 * t2374;
t2332 = -t2356 - t2357;
t2277 = t2304 * t2368 + t2305 * t2372;
t2259 = -t2364 * t2332 + t2365 * t2387;
t2258 = t2365 * t2332 + t2364 * t2387;
t2255 = -t2330 * pkin(2) + qJ(3) * t2383 + t2307 + t2441;
t2249 = -t2288 * t2368 + t2289 * t2372;
t2241 = -t2364 * t2307 + t2365 * t2388;
t2240 = t2365 * t2307 + t2364 * t2388;
t2235 = -t2259 * t2369 + t2277 * t2373;
t2234 = t2259 * t2373 + t2277 * t2369;
t2233 = t2271 * t2371 - t2460;
t2232 = t2266 * t2371 - t2460;
t2223 = t2254 * t2372 + t2257 * t2368;
t2222 = -t2265 * t2372 + t2462;
t2221 = t2264 * t2372 - t2462;
t2216 = t2230 * t2368 + t2447;
t2215 = t2231 * t2368 + t2447;
t2212 = -t2364 * t2255 + t2365 * t2390;
t2211 = t2365 * t2255 + t2364 * t2390;
t2206 = t2365 * t2393 + t2458;
t2205 = t2365 * t2394 - t2458;
t2204 = t2364 * t2393 - t2457;
t2203 = t2364 * t2394 + t2457;
t2199 = t2219 * t2367 + t2228 * t2371;
t2195 = t2217 * t2367 + t2226 * t2371;
t2193 = -t2364 * t2232 + t2365 * t2396;
t2192 = -t2364 * t2233 + t2365 * t2395;
t2191 = t2365 * t2232 + t2364 * t2396;
t2190 = t2365 * t2233 + t2364 * t2395;
t2189 = t2200 * t2367 + t2248 * t2371;
t2184 = -t2213 * t2367 + t2214 * t2371;
t2182 = t2198 * t2368 + t2220 * t2372;
t2181 = t2194 * t2368 + t2218 * t2372;
t2180 = t2183 * t2368 + t2236 * t2372;
t2179 = t2196 * t2371 + t2197 * t2367;
t2177 = t2188 * t2368 + t2201 * t2372;
t2176 = -t2364 * t2199 + t2365 * t2397;
t2175 = t2365 * t2199 + t2364 * t2397;
t2172 = -t2364 * t2195 + t2365 * t2398;
t2171 = t2365 * t2195 + t2364 * t2398;
t2170 = t2178 * t2368 + t2202 * t2372;
t2169 = -t2364 * t2189 + t2365 * t2399;
t2168 = t2365 * t2189 + t2364 * t2399;
t2167 = -t2364 * t2184 + t2365 * t2400;
t2166 = t2365 * t2184 + t2364 * t2400;
t2165 = -t2364 * t2179 + t2365 * t2401;
t2164 = t2365 * t2179 + t2364 * t2401;
t2161 = t2162 * t2367 + t2186 * t2371;
t2159 = t2160 * t2368 + t2163 * t2372;
t2158 = -t2364 * t2161 + t2365 * t2402;
t2157 = t2365 * t2161 + t2364 * t2402;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2342, -t2341, 0, -t2348 * t2369 + t2349 * t2373, 0, 0, 0, 0, 0, 0, -t2454, -t2456, t2235, -t2241 * t2369 + t2249 * t2373, 0, 0, 0, 0, 0, 0, t2235, t2454, t2456, -t2212 * t2369 + t2223 * t2373, 0, 0, 0, 0, 0, 0, -t2205 * t2369 + t2221 * t2373, -t2468, -t2193 * t2369 + t2216 * t2373, -t2167 * t2369 + t2180 * t2373, 0, 0, 0, 0, 0, 0, -t2192 * t2369 + t2215 * t2373, -t2206 * t2369 + t2222 * t2373, t2468, -t2165 * t2369 + t2170 * t2373, 0, 0, 0, 0, 0, 0, -t2172 * t2369 + t2181 * t2373, -t2176 * t2369 + t2182 * t2373, -t2169 * t2369 + t2177 * t2373, -t2158 * t2369 + t2159 * t2373; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2341, t2342, 0, t2348 * t2373 + t2349 * t2369, 0, 0, 0, 0, 0, 0, t2453, t2455, t2234, t2241 * t2373 + t2249 * t2369, 0, 0, 0, 0, 0, 0, t2234, -t2453, -t2455, t2212 * t2373 + t2223 * t2369, 0, 0, 0, 0, 0, 0, t2205 * t2373 + t2221 * t2369, t2467, t2193 * t2373 + t2216 * t2369, t2167 * t2373 + t2180 * t2369, 0, 0, 0, 0, 0, 0, t2192 * t2373 + t2215 * t2369, t2206 * t2373 + t2222 * t2369, -t2467, t2165 * t2373 + t2170 * t2369, 0, 0, 0, 0, 0, 0, t2172 * t2373 + t2181 * t2369, t2176 * t2373 + t2182 * t2369, t2169 * t2373 + t2177 * t2369, t2158 * t2373 + t2159 * t2369; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2273, t2261, t2258, t2240, 0, 0, 0, 0, 0, 0, t2258, -t2273, -t2261, t2211, 0, 0, 0, 0, 0, 0, t2203, t2207, t2191, t2166, 0, 0, 0, 0, 0, 0, t2190, t2204, -t2207, t2164, 0, 0, 0, 0, 0, 0, t2171, t2175, t2168, t2157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2374, -qJDD(1), 0, t2349, 0, 0, 0, 0, 0, 0, -t2297, -t2290, t2277, t2249, 0, 0, 0, 0, 0, 0, t2277, t2297, t2290, t2223, 0, 0, 0, 0, 0, 0, t2221, -t2225, t2216, t2180, 0, 0, 0, 0, 0, 0, t2215, t2222, t2225, t2170, 0, 0, 0, 0, 0, 0, t2181, t2182, t2177, t2159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2374, 0, t2348, 0, 0, 0, 0, 0, 0, t2275, t2263, t2259, t2241, 0, 0, 0, 0, 0, 0, t2259, -t2275, -t2263, t2212, 0, 0, 0, 0, 0, 0, t2205, t2209, t2193, t2167, 0, 0, 0, 0, 0, 0, t2192, t2206, -t2209, t2165, 0, 0, 0, 0, 0, 0, t2172, t2176, t2169, t2158; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2273, t2261, t2258, t2240, 0, 0, 0, 0, 0, 0, t2258, -t2273, -t2261, t2211, 0, 0, 0, 0, 0, 0, t2203, t2207, t2191, t2166, 0, 0, 0, 0, 0, 0, t2190, t2204, -t2207, t2164, 0, 0, 0, 0, 0, 0, t2171, t2175, t2168, t2157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2331, -t2328, t2305, t2289, 0, 0, 0, 0, 0, 0, t2305, -t2331, t2328, t2254, 0, 0, 0, 0, 0, 0, t2264, t2379, t2442, t2236, 0, 0, 0, 0, 0, 0, t2442, -t2265, -t2379, t2202, 0, 0, 0, 0, 0, 0, t2218, t2220, t2201, t2163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2440, t2311, -t2304, t2288, 0, 0, 0, 0, 0, 0, -t2304, t2440, -t2311, -t2257, 0, 0, 0, 0, 0, 0, t2245, t2250, -t2230, -t2183, 0, 0, 0, 0, 0, 0, -t2231, -t2245, -t2250, -t2178, 0, 0, 0, 0, 0, 0, -t2194, -t2198, -t2188, -t2160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2302, -t2383, t2332, t2307, 0, 0, 0, 0, 0, 0, t2332, t2302, t2383, t2255, 0, 0, 0, 0, 0, 0, t2246, -t2252, t2232, t2184, 0, 0, 0, 0, 0, 0, t2233, -t2246, t2252, t2179, 0, 0, 0, 0, 0, 0, t2195, t2199, t2189, t2161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2332, t2302, t2383, t2255, 0, 0, 0, 0, 0, 0, t2246, -t2252, t2232, t2184, 0, 0, 0, 0, 0, 0, t2233, -t2246, t2252, t2179, 0, 0, 0, 0, 0, 0, t2195, t2199, t2189, t2161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2305, t2331, -t2328, -t2254, 0, 0, 0, 0, 0, 0, -t2264, -t2379, -t2442, -t2236, 0, 0, 0, 0, 0, 0, -t2442, t2265, t2379, -t2202, 0, 0, 0, 0, 0, 0, -t2218, -t2220, -t2201, -t2163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2304, -t2440, t2311, t2257, 0, 0, 0, 0, 0, 0, -t2245, -t2250, t2230, t2183, 0, 0, 0, 0, 0, 0, t2231, t2245, t2250, t2178, 0, 0, 0, 0, 0, 0, t2194, t2198, t2188, t2160; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2287, -t2280, t2266, t2214, 0, 0, 0, 0, 0, 0, t2271, -t2287, t2280, t2196, 0, 0, 0, 0, 0, 0, t2226, t2228, t2248, t2186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2405, t2408, t2270, t2213, 0, 0, 0, 0, 0, 0, t2270, t2405, -t2408, -t2197, 0, 0, 0, 0, 0, 0, -t2217, -t2219, -t2200, -t2162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2264, t2379, t2442, t2236, 0, 0, 0, 0, 0, 0, t2442, -t2265, -t2379, t2202, 0, 0, 0, 0, 0, 0, t2218, t2220, t2201, t2163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2442, -t2265, -t2379, t2202, 0, 0, 0, 0, 0, 0, t2218, t2220, t2201, t2163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2271, t2287, -t2280, -t2196, 0, 0, 0, 0, 0, 0, -t2226, -t2228, -t2248, -t2186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2270, -t2405, t2408, t2197, 0, 0, 0, 0, 0, 0, t2217, t2219, t2200, t2162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2256, t2243, t2227, t2174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2242, t2272, t2229, t2173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2226, t2228, t2248, t2186;];
f_new_reg  = t1;
