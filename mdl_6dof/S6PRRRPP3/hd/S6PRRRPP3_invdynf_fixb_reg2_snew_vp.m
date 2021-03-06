% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6PRRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 06:56
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6PRRRPP3_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPP3_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPP3_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_invdynf_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:56:00
% EndTime: 2019-05-05 06:56:09
% DurationCPUTime: 9.53s
% Computational Cost: add. (13883->297), mult. (27064->353), div. (0->0), fcn. (19143->10), ass. (0->210)
t2433 = cos(qJ(3));
t2429 = sin(qJ(4));
t2432 = cos(qJ(4));
t2430 = sin(qJ(3));
t2485 = qJD(2) * t2430;
t2391 = -qJD(3) * t2432 + t2429 * t2485;
t2484 = qJD(2) * t2433;
t2415 = -qJD(4) + t2484;
t2374 = t2391 * t2415;
t2470 = qJD(3) * t2484;
t2472 = t2430 * qJDD(2);
t2397 = t2470 + t2472;
t2444 = t2429 * qJDD(3) + t2432 * t2397;
t2439 = qJD(4) * t2391 - t2444;
t2494 = t2374 - t2439;
t2393 = qJD(3) * t2429 + t2432 * t2485;
t2389 = t2393 ^ 2;
t2489 = t2415 ^ 2;
t2355 = t2489 + t2389;
t2417 = qJD(3) * t2485;
t2471 = t2433 * qJDD(2);
t2466 = t2417 - t2471;
t2463 = -qJDD(4) - t2466;
t2475 = t2391 * t2393;
t2441 = t2463 - t2475;
t2320 = t2355 * t2429 + t2432 * t2441;
t2509 = t2320 * t2430;
t2283 = t2433 * t2494 - t2509;
t2424 = sin(pkin(6));
t2426 = cos(pkin(6));
t2507 = t2320 * t2433;
t2292 = t2430 * t2494 + t2507;
t2431 = sin(qJ(2));
t2434 = cos(qJ(2));
t2495 = t2355 * t2432 - t2441 * t2429;
t2506 = t2434 * t2495;
t2452 = t2292 * t2431 + t2506;
t2254 = t2424 * t2283 + t2426 * t2452;
t2508 = t2431 * t2495;
t2269 = t2292 * t2434 - t2508;
t2423 = sin(pkin(10));
t2425 = cos(pkin(10));
t2541 = t2254 * t2423 - t2269 * t2425;
t2540 = t2254 * t2425 + t2269 * t2423;
t2376 = t2393 * t2415;
t2467 = -t2432 * qJDD(3) + t2429 * t2397;
t2442 = qJD(4) * t2393 + t2467;
t2493 = -t2376 + t2442;
t2440 = t2463 + t2475;
t2388 = t2391 ^ 2;
t2491 = -t2388 - t2489;
t2505 = t2440 * t2429 + t2491 * t2432;
t2515 = t2430 * t2505;
t2291 = -t2433 * t2493 + t2515;
t2513 = t2433 * t2505;
t2294 = t2430 * t2493 + t2513;
t2504 = t2491 * t2429 - t2440 * t2432;
t2512 = t2434 * t2504;
t2450 = t2294 * t2431 - t2512;
t2256 = -t2424 * t2291 + t2426 * t2450;
t2514 = t2431 * t2504;
t2271 = t2294 * t2434 + t2514;
t2537 = t2256 * t2423 - t2271 * t2425;
t2536 = t2256 * t2425 + t2271 * t2423;
t2251 = -t2426 * t2283 + t2424 * t2452;
t2329 = (-qJD(4) - t2415) * t2393 - t2467;
t2334 = -t2439 - t2374;
t2481 = t2334 * t2432;
t2510 = t2329 * t2429 - t2481;
t2492 = -t2388 - t2389;
t2503 = t2430 * t2492;
t2482 = t2334 * t2429;
t2511 = t2432 * t2329 + t2482;
t2516 = t2433 * t2511 + t2503;
t2520 = t2431 * t2510 + t2434 * t2516;
t2500 = t2433 * t2492;
t2517 = t2430 * t2511 - t2500;
t2521 = t2431 * t2516 - t2434 * t2510;
t2526 = -t2424 * t2517 + t2426 * t2521;
t2531 = -t2423 * t2526 + t2425 * t2520;
t2253 = t2426 * t2291 + t2424 * t2450;
t2530 = t2423 * t2520 + t2425 * t2526;
t2527 = t2424 * t2521 + t2426 * t2517;
t2465 = g(1) * t2423 - g(2) * t2425;
t2443 = t2426 * t2465;
t2486 = -g(3) + qJDD(1);
t2468 = t2424 * t2486;
t2497 = t2443 + t2468;
t2490 = qJD(3) ^ 2;
t2488 = -2 * qJD(5);
t2487 = 2 * qJD(6);
t2474 = qJD(4) - t2415;
t2403 = -g(1) * t2425 - g(2) * t2423;
t2350 = t2434 * t2403 + t2431 * t2497;
t2435 = qJD(2) ^ 2;
t2338 = -pkin(2) * t2435 + qJDD(2) * pkin(8) + t2350;
t2438 = -t2424 * t2465 + t2426 * t2486;
t2324 = t2433 * t2338 + t2430 * t2438;
t2419 = t2430 ^ 2;
t2420 = t2433 ^ 2;
t2473 = t2419 + t2420;
t2469 = -pkin(4) * t2415 + t2488;
t2395 = (-pkin(3) * t2433 - pkin(9) * t2430) * qJD(2);
t2304 = -pkin(3) * t2490 + qJDD(3) * pkin(9) + t2395 * t2484 + t2324;
t2464 = t2431 * t2403 - t2434 * t2497;
t2337 = -qJDD(2) * pkin(2) - t2435 * pkin(8) + t2464;
t2436 = (-t2397 - t2470) * pkin(9) + (t2466 + t2417) * pkin(3) + t2337;
t2273 = -t2429 * t2304 + t2432 * t2436;
t2358 = pkin(4) * t2391 - qJ(5) * t2393;
t2264 = t2463 * pkin(4) - t2489 * qJ(5) + t2393 * t2358 + qJDD(5) - t2273;
t2257 = -t2439 * pkin(5) + (-pkin(5) * t2391 + t2487) * t2415 + t2264 + t2440 * qJ(6);
t2274 = t2432 * t2304 + t2429 * t2436;
t2263 = -t2489 * pkin(4) - qJ(5) * t2463 - t2391 * t2358 + t2415 * t2488 + t2274;
t2370 = pkin(5) * t2393 + qJ(6) * t2415;
t2258 = -pkin(5) * t2442 - t2388 * qJ(6) - t2415 * t2370 + qJDD(6) + t2263;
t2230 = t2257 * t2429 + t2258 * t2432;
t2369 = t2433 * t2438;
t2303 = -t2369 - qJDD(3) * pkin(3) - t2490 * pkin(9) + (qJD(2) * t2395 + t2338) * t2430;
t2437 = t2442 * pkin(4) - qJ(5) * t2494 + t2303;
t2259 = -t2388 * pkin(5) + t2467 * qJ(6) + t2391 * t2487 + (qJ(6) * qJD(4) - t2370 + t2469) * t2393 + t2437;
t2225 = t2230 * t2433 + t2259 * t2430;
t2229 = -t2257 * t2432 + t2258 * t2429;
t2462 = t2225 * t2431 - t2229 * t2434;
t2232 = t2263 * t2432 + t2264 * t2429;
t2268 = t2393 * t2469 + t2437;
t2228 = t2232 * t2433 + t2268 * t2430;
t2231 = t2263 * t2429 - t2264 * t2432;
t2461 = t2228 * t2431 - t2231 * t2434;
t2248 = -t2273 * t2429 + t2274 * t2432;
t2234 = t2248 * t2433 + t2303 * t2430;
t2247 = t2273 * t2432 + t2274 * t2429;
t2460 = t2234 * t2431 - t2247 * t2434;
t2330 = -t2376 - t2442;
t2300 = t2330 * t2432 + t2482;
t2280 = t2300 * t2433 + t2503;
t2297 = t2330 * t2429 - t2481;
t2457 = t2280 * t2431 - t2297 * t2434;
t2323 = -t2430 * t2338 + t2369;
t2282 = -t2323 * t2430 + t2324 * t2433;
t2456 = t2282 * t2431 - t2337 * t2434;
t2327 = t2393 * t2474 + t2467;
t2287 = t2327 * t2430 + t2513;
t2454 = t2287 * t2431 - t2512;
t2331 = -t2391 * t2474 + t2444;
t2293 = -t2331 * t2430 - t2507;
t2451 = t2293 * t2431 - t2506;
t2449 = t2431 * t2350 - t2434 * t2464;
t2414 = t2430 * t2435 * t2433;
t2404 = qJDD(3) + t2414;
t2413 = -t2420 * t2435 - t2490;
t2366 = -t2404 * t2430 + t2413 * t2433;
t2398 = -0.2e1 * t2417 + t2471;
t2448 = t2366 * t2431 + t2398 * t2434;
t2405 = -qJDD(3) + t2414;
t2412 = -t2419 * t2435 - t2490;
t2367 = t2405 * t2433 - t2412 * t2430;
t2396 = 0.2e1 * t2470 + t2472;
t2447 = t2367 * t2431 - t2396 * t2434;
t2399 = t2473 * qJDD(2);
t2402 = t2473 * t2435;
t2446 = t2399 * t2431 + t2402 * t2434;
t2445 = qJDD(2) * t2434 - t2431 * t2435;
t2401 = -qJDD(2) * t2431 - t2434 * t2435;
t2380 = t2445 * t2426;
t2379 = t2401 * t2426;
t2378 = t2445 * t2424;
t2377 = t2401 * t2424;
t2365 = t2405 * t2430 + t2412 * t2433;
t2364 = t2404 * t2433 + t2413 * t2430;
t2362 = t2399 * t2434 - t2402 * t2431;
t2357 = t2446 * t2426;
t2356 = t2446 * t2424;
t2340 = t2367 * t2434 + t2396 * t2431;
t2339 = t2366 * t2434 - t2398 * t2431;
t2316 = -t2424 * t2365 + t2426 * t2447;
t2315 = -t2424 * t2364 + t2426 * t2448;
t2314 = t2426 * t2365 + t2424 * t2447;
t2313 = t2426 * t2364 + t2424 * t2448;
t2306 = t2350 * t2434 + t2431 * t2464;
t2302 = t2424 ^ 2 * t2465 + (-t2468 + t2449) * t2426;
t2301 = t2426 ^ 2 * t2486 + (t2449 - t2443) * t2424;
t2290 = t2331 * t2433 - t2509;
t2284 = -t2327 * t2433 + t2515;
t2281 = t2323 * t2433 + t2324 * t2430;
t2277 = t2300 * t2430 - t2500;
t2272 = t2282 * t2434 + t2337 * t2431;
t2270 = t2293 * t2434 + t2508;
t2266 = t2287 * t2434 + t2514;
t2262 = t2280 * t2434 + t2297 * t2431;
t2255 = -t2424 * t2290 + t2426 * t2451;
t2252 = t2426 * t2290 + t2424 * t2451;
t2250 = -t2424 * t2281 + t2426 * t2456;
t2249 = t2426 * t2281 + t2424 * t2456;
t2245 = -t2424 * t2284 + t2426 * t2454;
t2242 = t2426 * t2284 + t2424 * t2454;
t2240 = -t2424 * t2277 + t2426 * t2457;
t2237 = t2426 * t2277 + t2424 * t2457;
t2233 = t2248 * t2430 - t2303 * t2433;
t2227 = t2232 * t2430 - t2268 * t2433;
t2226 = t2234 * t2434 + t2247 * t2431;
t2224 = t2230 * t2430 - t2259 * t2433;
t2223 = t2228 * t2434 + t2231 * t2431;
t2222 = -t2424 * t2233 + t2426 * t2460;
t2221 = t2426 * t2233 + t2424 * t2460;
t2220 = t2225 * t2434 + t2229 * t2431;
t2219 = -t2424 * t2227 + t2426 * t2461;
t2218 = t2426 * t2227 + t2424 * t2461;
t2217 = -t2424 * t2224 + t2426 * t2462;
t2216 = t2426 * t2224 + t2424 * t2462;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t2425 * t2403 - t2423 * t2465, 0, 0, 0, 0, 0, 0, -t2380 * t2423 + t2401 * t2425, -t2379 * t2423 - t2425 * t2445, 0, -t2302 * t2423 + t2306 * t2425, 0, 0, 0, 0, 0, 0, -t2315 * t2423 + t2339 * t2425, -t2316 * t2423 + t2340 * t2425, -t2357 * t2423 + t2362 * t2425, -t2250 * t2423 + t2272 * t2425, 0, 0, 0, 0, 0, 0, -t2245 * t2423 + t2266 * t2425, -t2541, -t2240 * t2423 + t2262 * t2425, -t2222 * t2423 + t2226 * t2425, 0, 0, 0, 0, 0, 0, t2531, t2537, -t2255 * t2423 + t2270 * t2425, -t2219 * t2423 + t2223 * t2425, 0, 0, 0, 0, 0, 0, t2531, t2541, -t2537, -t2217 * t2423 + t2220 * t2425; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t2423 * t2403 + t2425 * t2465, 0, 0, 0, 0, 0, 0, t2380 * t2425 + t2401 * t2423, t2379 * t2425 - t2423 * t2445, 0, t2302 * t2425 + t2306 * t2423, 0, 0, 0, 0, 0, 0, t2315 * t2425 + t2339 * t2423, t2316 * t2425 + t2340 * t2423, t2357 * t2425 + t2362 * t2423, t2250 * t2425 + t2272 * t2423, 0, 0, 0, 0, 0, 0, t2245 * t2425 + t2266 * t2423, t2540, t2240 * t2425 + t2262 * t2423, t2222 * t2425 + t2226 * t2423, 0, 0, 0, 0, 0, 0, t2530, -t2536, t2255 * t2425 + t2270 * t2423, t2219 * t2425 + t2223 * t2423, 0, 0, 0, 0, 0, 0, t2530, -t2540, t2536, t2217 * t2425 + t2220 * t2423; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t2486, 0, 0, 0, 0, 0, 0, t2378, t2377, 0, t2301, 0, 0, 0, 0, 0, 0, t2313, t2314, t2356, t2249, 0, 0, 0, 0, 0, 0, t2242, t2251, t2237, t2221, 0, 0, 0, 0, 0, 0, t2527, -t2253, t2252, t2218, 0, 0, 0, 0, 0, 0, t2527, -t2251, t2253, t2216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2403, 0, 0, 0, 0, 0, 0, t2401, -t2445, 0, t2306, 0, 0, 0, 0, 0, 0, t2339, t2340, t2362, t2272, 0, 0, 0, 0, 0, 0, t2266, t2269, t2262, t2226, 0, 0, 0, 0, 0, 0, t2520, -t2271, t2270, t2223, 0, 0, 0, 0, 0, 0, t2520, -t2269, t2271, t2220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2465, 0, 0, 0, 0, 0, 0, t2380, t2379, 0, t2302, 0, 0, 0, 0, 0, 0, t2315, t2316, t2357, t2250, 0, 0, 0, 0, 0, 0, t2245, t2254, t2240, t2222, 0, 0, 0, 0, 0, 0, t2526, -t2256, t2255, t2219, 0, 0, 0, 0, 0, 0, t2526, -t2254, t2256, t2217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2486, 0, 0, 0, 0, 0, 0, t2378, t2377, 0, t2301, 0, 0, 0, 0, 0, 0, t2313, t2314, t2356, t2249, 0, 0, 0, 0, 0, 0, t2242, t2251, t2237, t2221, 0, 0, 0, 0, 0, 0, t2527, -t2253, t2252, t2218, 0, 0, 0, 0, 0, 0, t2527, -t2251, t2253, t2216; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2435, -qJDD(2), 0, t2350, 0, 0, 0, 0, 0, 0, t2366, t2367, t2399, t2282, 0, 0, 0, 0, 0, 0, t2287, t2292, t2280, t2234, 0, 0, 0, 0, 0, 0, t2516, -t2294, t2293, t2228, 0, 0, 0, 0, 0, 0, t2516, -t2292, t2294, t2225; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t2435, 0, -t2464, 0, 0, 0, 0, 0, 0, t2398, -t2396, t2402, -t2337, 0, 0, 0, 0, 0, 0, -t2504, t2495, -t2297, -t2247, 0, 0, 0, 0, 0, 0, -t2510, t2504, -t2495, -t2231, 0, 0, 0, 0, 0, 0, -t2510, -t2495, -t2504, -t2229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2438, 0, 0, 0, 0, 0, 0, t2364, t2365, 0, t2281, 0, 0, 0, 0, 0, 0, t2284, -t2283, t2277, t2233, 0, 0, 0, 0, 0, 0, t2517, -t2291, t2290, t2227, 0, 0, 0, 0, 0, 0, t2517, t2283, t2291, t2224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2413, t2405, t2471, t2324, 0, 0, 0, 0, 0, 0, t2505, t2320, t2300, t2248, 0, 0, 0, 0, 0, 0, t2511, -t2505, -t2320, t2232, 0, 0, 0, 0, 0, 0, t2511, -t2320, t2505, t2230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2404, t2412, -t2472, t2323, 0, 0, 0, 0, 0, 0, -t2327, -t2494, -t2492, -t2303, 0, 0, 0, 0, 0, 0, -t2492, t2493, t2331, -t2268, 0, 0, 0, 0, 0, 0, -t2492, t2494, -t2493, -t2259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2398, t2396, -t2402, t2337, 0, 0, 0, 0, 0, 0, t2504, -t2495, t2297, t2247, 0, 0, 0, 0, 0, 0, t2510, -t2504, t2495, t2231, 0, 0, 0, 0, 0, 0, t2510, t2495, t2504, t2229; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2491, t2441, t2330, t2274, 0, 0, 0, 0, 0, 0, t2329, -t2491, -t2441, t2263, 0, 0, 0, 0, 0, 0, t2329, -t2441, t2491, t2258; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2440, -t2355, -t2334, t2273, 0, 0, 0, 0, 0, 0, -t2334, t2440, t2355, -t2264, 0, 0, 0, 0, 0, 0, -t2334, t2355, -t2440, -t2257; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2327, t2494, t2492, t2303, 0, 0, 0, 0, 0, 0, t2492, -t2493, -t2331, t2268, 0, 0, 0, 0, 0, 0, t2492, -t2494, t2493, t2259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2492, -t2493, -t2331, t2268, 0, 0, 0, 0, 0, 0, t2492, -t2494, t2493, t2259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2329, t2491, t2441, -t2263, 0, 0, 0, 0, 0, 0, -t2329, t2441, -t2491, -t2258; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2334, -t2440, -t2355, t2264, 0, 0, 0, 0, 0, 0, t2334, -t2355, t2440, t2257; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2492, -t2494, t2493, t2259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2334, -t2355, t2440, t2257; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2329, -t2441, t2491, t2258;];
f_new_reg  = t1;
