% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RRRPRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 09:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RRRPRP11_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP11_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP11_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_invdynf_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 09:21:38
% EndTime: 2019-05-07 09:21:45
% DurationCPUTime: 7.81s
% Computational Cost: add. (33033->313), mult. (71774->409), div. (0->0), fcn. (56182->10), ass. (0->247)
t2462 = sin(pkin(6));
t2466 = sin(qJ(2));
t2526 = t2462 * t2466;
t2517 = qJD(1) * t2526;
t2453 = qJD(2) * t2517;
t2470 = cos(qJ(2));
t2519 = qJDD(1) * t2470;
t2507 = t2462 * t2519 - t2453;
t2486 = -qJDD(3) + t2507;
t2465 = sin(qJ(3));
t2469 = cos(qJ(3));
t2463 = cos(pkin(6));
t2512 = qJD(1) * t2463 + qJD(2);
t2425 = t2465 * t2517 - t2469 * t2512;
t2427 = t2465 * t2512 + t2469 * t2517;
t2529 = t2425 * t2427;
t2388 = t2486 + t2529;
t2424 = t2425 ^ 2;
t2525 = t2462 * t2470;
t2516 = qJD(1) * t2525;
t2449 = -qJD(3) + t2516;
t2445 = t2449 ^ 2;
t2396 = -t2445 - t2424;
t2353 = t2388 * t2465 + t2396 * t2469;
t2554 = t2353 * t2466;
t2553 = t2353 * t2470;
t2352 = t2388 * t2469 - t2396 * t2465;
t2552 = t2462 * t2352;
t2551 = t2463 * t2352;
t2482 = t2486 - t2529;
t2537 = t2427 ^ 2;
t2514 = -t2445 - t2537;
t2361 = t2465 * t2514 - t2469 * t2482;
t2550 = t2361 * t2466;
t2549 = t2361 * t2470;
t2520 = qJDD(1) * t2462;
t2436 = qJD(2) * t2516 + t2466 * t2520;
t2518 = t2463 * qJDD(1) + qJDD(2);
t2483 = -t2469 * t2436 - t2465 * t2518;
t2523 = qJD(3) + t2449;
t2379 = t2425 * t2523 + t2483;
t2548 = t2379 * t2465;
t2547 = t2379 * t2469;
t2358 = t2465 * t2482 + t2469 * t2514;
t2546 = t2462 * t2358;
t2545 = t2463 * t2358;
t2472 = qJD(1) ^ 2;
t2544 = t2462 * t2472;
t2540 = -t2424 - t2537;
t2543 = t2466 * t2540;
t2542 = t2470 * t2540;
t2480 = -t2425 * qJD(3) - t2483;
t2378 = t2425 * t2449 + t2480;
t2508 = t2512 ^ 2;
t2423 = qJD(5) + t2427;
t2541 = qJD(5) + t2423;
t2464 = sin(qJ(5));
t2468 = cos(qJ(5));
t2405 = -t2468 * t2425 - t2449 * t2464;
t2404 = t2405 ^ 2;
t2407 = t2425 * t2464 - t2449 * t2468;
t2539 = t2407 ^ 2;
t2538 = t2423 ^ 2;
t2536 = -2 * qJD(4);
t2535 = -2 * qJD(6);
t2534 = t2463 * g(3);
t2533 = (-pkin(2) * t2470 - pkin(9) * t2466) * t2544;
t2530 = t2407 * t2405;
t2527 = t2462 ^ 2 * t2472;
t2524 = qJD(3) - t2449;
t2522 = qJD(5) - t2423;
t2467 = sin(qJ(1));
t2471 = cos(qJ(1));
t2452 = -g(1) * t2471 - g(2) * t2467;
t2432 = -pkin(1) * t2472 + pkin(8) * t2520 + t2452;
t2451 = t2467 * g(1) - t2471 * g(2);
t2479 = qJDD(1) * pkin(1) + pkin(8) * t2544 + t2451;
t2476 = t2463 * t2479;
t2521 = t2470 * t2432 + t2466 * t2476;
t2367 = t2518 * pkin(9) - t2508 * pkin(2) + (-g(3) * t2466 + t2470 * t2533) * t2462 + t2521;
t2506 = qJD(1) * t2512;
t2484 = t2466 * t2506;
t2485 = t2470 * t2506;
t2368 = t2453 * pkin(2) - t2436 * pkin(9) - t2534 + (-pkin(9) * t2485 + (t2484 - t2519) * pkin(2) - t2479) * t2462;
t2337 = -t2465 * t2367 + t2469 * t2368;
t2400 = pkin(3) * t2425 - qJ(4) * t2427;
t2325 = pkin(3) * t2486 - t2445 * qJ(4) + t2427 * t2400 + qJDD(4) - t2337;
t2300 = -pkin(4) * t2379 + pkin(10) * t2388 + t2325;
t2509 = t2465 * t2436 - t2469 * t2518;
t2395 = qJD(3) * t2427 + t2509;
t2408 = pkin(4) * t2427 + pkin(10) * t2449;
t2510 = t2466 * t2432 - t2470 * t2476;
t2366 = -t2518 * pkin(2) - t2508 * pkin(9) + (g(3) * t2470 + t2466 * t2533) * t2462 + t2510;
t2474 = t2395 * pkin(3) - qJ(4) * t2378 + t2366;
t2513 = -pkin(3) * t2449 + t2536;
t2307 = -t2424 * pkin(4) + t2395 * pkin(10) + (-t2408 + t2513) * t2427 + t2474;
t2289 = t2464 * t2300 + t2468 * t2307;
t2515 = -t2538 - t2539;
t2511 = -t2468 * t2395 - t2464 * t2486;
t2355 = -qJD(5) * t2407 - t2511;
t2385 = pkin(5) * t2423 - qJ(6) * t2407;
t2280 = -pkin(5) * t2404 + qJ(6) * t2355 - t2385 * t2423 + t2405 * t2535 + t2289;
t2478 = -qJDD(5) - t2480;
t2475 = -t2478 - t2530;
t2491 = -t2464 * t2395 + t2468 * t2486;
t2481 = t2405 * t2522 + t2491;
t2501 = t2468 * t2300 - t2464 * t2307;
t2473 = pkin(5) * t2475 + qJ(6) * t2481 + t2407 * t2535 + t2501;
t2263 = t2464 * t2280 + t2468 * t2473;
t2338 = t2469 * t2367 + t2465 * t2368;
t2477 = -t2445 * pkin(3) - qJ(4) * t2486 - t2425 * t2400 + t2338;
t2306 = -t2395 * pkin(4) - t2424 * pkin(10) + (t2536 - t2408) * t2449 + t2477;
t2291 = -t2355 * pkin(5) - t2404 * qJ(6) + t2407 * t2385 + qJDD(6) + t2306;
t2262 = t2263 * t2465 + t2291 * t2469;
t2264 = t2468 * t2280 - t2464 * t2473;
t2505 = t2262 * t2466 - t2264 * t2470;
t2271 = t2464 * t2289 + t2468 * t2501;
t2268 = t2271 * t2465 + t2306 * t2469;
t2272 = t2468 * t2289 - t2464 * t2501;
t2504 = t2268 * t2466 - t2272 * t2470;
t2324 = t2449 * t2536 + t2477;
t2295 = t2324 * t2469 + t2325 * t2465;
t2326 = t2427 * t2513 + t2474;
t2503 = t2295 * t2466 - t2326 * t2470;
t2341 = -t2407 * t2522 - t2511;
t2314 = t2464 * t2341 + t2468 * t2481;
t2356 = -t2404 - t2539;
t2298 = t2314 * t2465 + t2356 * t2469;
t2315 = t2468 * t2341 - t2464 * t2481;
t2502 = t2298 * t2466 - t2315 * t2470;
t2362 = -t2538 - t2404;
t2329 = t2464 * t2362 + t2468 * t2475;
t2340 = t2407 * t2541 + t2511;
t2309 = t2329 * t2465 + t2340 * t2469;
t2330 = t2468 * t2362 - t2464 * t2475;
t2500 = t2309 * t2466 - t2330 * t2470;
t2350 = t2478 - t2530;
t2331 = t2464 * t2350 + t2468 * t2515;
t2342 = -t2405 * t2541 - t2491;
t2311 = t2331 * t2465 + t2342 * t2469;
t2332 = t2468 * t2350 - t2464 * t2515;
t2499 = t2311 * t2466 - t2332 * t2470;
t2313 = -t2337 * t2465 + t2338 * t2469;
t2498 = t2313 * t2466 - t2366 * t2470;
t2373 = t2427 * t2523 + t2509;
t2345 = -t2373 * t2469 - t2548;
t2497 = t2345 * t2466 - t2542;
t2415 = t2449 * t2427;
t2375 = -t2395 - t2415;
t2346 = t2375 * t2469 - t2548;
t2496 = t2346 * t2466 - t2542;
t2372 = t2427 * t2524 + t2509;
t2495 = -t2372 * t2470 + t2554;
t2374 = t2395 - t2415;
t2494 = t2374 * t2470 - t2554;
t2493 = -t2378 * t2470 - t2550;
t2376 = -t2425 * t2524 - t2483;
t2492 = t2376 * t2470 + t2550;
t2397 = -g(3) * t2525 - t2510;
t2398 = -g(3) * t2526 + t2521;
t2490 = t2397 * t2470 + t2398 * t2466;
t2411 = t2462 * t2485 - t2436;
t2439 = t2462 * t2484;
t2412 = t2439 + t2507;
t2489 = t2411 * t2470 + t2412 * t2466;
t2460 = t2466 ^ 2;
t2421 = -t2460 * t2527 - t2508;
t2448 = t2470 * t2466 * t2527;
t2434 = t2448 - t2518;
t2488 = t2421 * t2470 + t2434 * t2466;
t2433 = t2448 + t2518;
t2461 = t2470 ^ 2;
t2437 = -t2461 * t2527 - t2508;
t2487 = t2433 * t2470 + t2437 * t2466;
t2447 = -qJDD(1) * t2467 - t2471 * t2472;
t2446 = qJDD(1) * t2471 - t2467 * t2472;
t2438 = (-t2460 - t2461) * t2527;
t2416 = -t2462 * t2479 - t2534;
t2413 = t2439 - t2507;
t2410 = t2512 * t2516 + t2436;
t2403 = -t2433 * t2466 + t2437 * t2470;
t2399 = -t2421 * t2466 + t2434 * t2470;
t2382 = -t2411 * t2466 + t2412 * t2470;
t2381 = -t2462 * t2413 + t2463 * t2487;
t2380 = t2463 * t2413 + t2462 * t2487;
t2370 = -t2462 * t2410 + t2463 * t2488;
t2369 = t2463 * t2410 + t2462 * t2488;
t2365 = -t2462 * t2438 + t2463 * t2489;
t2364 = t2463 * t2438 + t2462 * t2489;
t2357 = -t2397 * t2466 + t2398 * t2470;
t2348 = -t2462 * t2416 + t2463 * t2490;
t2347 = t2463 * t2416 + t2462 * t2490;
t2344 = t2375 * t2465 + t2547;
t2343 = -t2373 * t2465 + t2547;
t2336 = -t2376 * t2466 + t2549;
t2335 = t2378 * t2466 - t2549;
t2334 = -t2374 * t2466 - t2553;
t2333 = t2372 * t2466 + t2553;
t2328 = t2346 * t2470 + t2543;
t2327 = t2345 * t2470 + t2543;
t2323 = t2463 * t2492 + t2546;
t2322 = t2463 * t2493 - t2546;
t2321 = t2462 * t2492 - t2545;
t2320 = t2462 * t2493 + t2545;
t2319 = t2463 * t2494 - t2552;
t2318 = t2463 * t2495 + t2552;
t2317 = t2462 * t2494 + t2551;
t2316 = t2462 * t2495 - t2551;
t2312 = t2337 * t2469 + t2338 * t2465;
t2310 = -t2331 * t2469 + t2342 * t2465;
t2308 = -t2329 * t2469 + t2340 * t2465;
t2304 = -t2462 * t2344 + t2463 * t2496;
t2303 = -t2462 * t2343 + t2463 * t2497;
t2302 = t2463 * t2344 + t2462 * t2496;
t2301 = t2463 * t2343 + t2462 * t2497;
t2297 = -t2314 * t2469 + t2356 * t2465;
t2296 = t2313 * t2470 + t2366 * t2466;
t2294 = t2324 * t2465 - t2325 * t2469;
t2293 = t2311 * t2470 + t2332 * t2466;
t2292 = t2309 * t2470 + t2330 * t2466;
t2290 = t2298 * t2470 + t2315 * t2466;
t2287 = -t2462 * t2312 + t2463 * t2498;
t2286 = t2463 * t2312 + t2462 * t2498;
t2285 = t2295 * t2470 + t2326 * t2466;
t2284 = -t2462 * t2310 + t2463 * t2499;
t2283 = t2463 * t2310 + t2462 * t2499;
t2282 = -t2462 * t2308 + t2463 * t2500;
t2281 = t2463 * t2308 + t2462 * t2500;
t2278 = -t2462 * t2297 + t2463 * t2502;
t2277 = t2463 * t2297 + t2462 * t2502;
t2276 = -t2462 * t2294 + t2463 * t2503;
t2275 = t2463 * t2294 + t2462 * t2503;
t2274 = -t2284 * t2467 + t2293 * t2471;
t2273 = t2284 * t2471 + t2293 * t2467;
t2270 = -t2282 * t2467 + t2292 * t2471;
t2269 = t2282 * t2471 + t2292 * t2467;
t2267 = -t2271 * t2469 + t2306 * t2465;
t2266 = -t2278 * t2467 + t2290 * t2471;
t2265 = t2278 * t2471 + t2290 * t2467;
t2261 = -t2263 * t2469 + t2291 * t2465;
t2260 = t2268 * t2470 + t2272 * t2466;
t2259 = t2262 * t2470 + t2264 * t2466;
t2258 = -t2462 * t2267 + t2463 * t2504;
t2257 = t2463 * t2267 + t2462 * t2504;
t2256 = -t2462 * t2261 + t2463 * t2505;
t2255 = t2463 * t2261 + t2462 * t2505;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2447, -t2446, 0, -t2451 * t2467 + t2452 * t2471, 0, 0, 0, 0, 0, 0, -t2381 * t2467 + t2403 * t2471, -t2370 * t2467 + t2399 * t2471, -t2365 * t2467 + t2382 * t2471, -t2348 * t2467 + t2357 * t2471, 0, 0, 0, 0, 0, 0, -t2318 * t2467 + t2333 * t2471, -t2322 * t2467 + t2335 * t2471, -t2304 * t2467 + t2328 * t2471, -t2287 * t2467 + t2296 * t2471, 0, 0, 0, 0, 0, 0, -t2303 * t2467 + t2327 * t2471, -t2319 * t2467 + t2334 * t2471, -t2323 * t2467 + t2336 * t2471, -t2276 * t2467 + t2285 * t2471, 0, 0, 0, 0, 0, 0, t2270, t2274, t2266, -t2258 * t2467 + t2260 * t2471, 0, 0, 0, 0, 0, 0, t2270, t2274, t2266, -t2256 * t2467 + t2259 * t2471; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2446, t2447, 0, t2451 * t2471 + t2452 * t2467, 0, 0, 0, 0, 0, 0, t2381 * t2471 + t2403 * t2467, t2370 * t2471 + t2399 * t2467, t2365 * t2471 + t2382 * t2467, t2348 * t2471 + t2357 * t2467, 0, 0, 0, 0, 0, 0, t2318 * t2471 + t2333 * t2467, t2322 * t2471 + t2335 * t2467, t2304 * t2471 + t2328 * t2467, t2287 * t2471 + t2296 * t2467, 0, 0, 0, 0, 0, 0, t2303 * t2471 + t2327 * t2467, t2319 * t2471 + t2334 * t2467, t2323 * t2471 + t2336 * t2467, t2276 * t2471 + t2285 * t2467, 0, 0, 0, 0, 0, 0, t2269, t2273, t2265, t2258 * t2471 + t2260 * t2467, 0, 0, 0, 0, 0, 0, t2269, t2273, t2265, t2256 * t2471 + t2259 * t2467; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2380, t2369, t2364, t2347, 0, 0, 0, 0, 0, 0, t2316, t2320, t2302, t2286, 0, 0, 0, 0, 0, 0, t2301, t2317, t2321, t2275, 0, 0, 0, 0, 0, 0, t2281, t2283, t2277, t2257, 0, 0, 0, 0, 0, 0, t2281, t2283, t2277, t2255; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2472, -qJDD(1), 0, t2452, 0, 0, 0, 0, 0, 0, t2403, t2399, t2382, t2357, 0, 0, 0, 0, 0, 0, t2333, t2335, t2328, t2296, 0, 0, 0, 0, 0, 0, t2327, t2334, t2336, t2285, 0, 0, 0, 0, 0, 0, t2292, t2293, t2290, t2260, 0, 0, 0, 0, 0, 0, t2292, t2293, t2290, t2259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2472, 0, t2451, 0, 0, 0, 0, 0, 0, t2381, t2370, t2365, t2348, 0, 0, 0, 0, 0, 0, t2318, t2322, t2304, t2287, 0, 0, 0, 0, 0, 0, t2303, t2319, t2323, t2276, 0, 0, 0, 0, 0, 0, t2282, t2284, t2278, t2258, 0, 0, 0, 0, 0, 0, t2282, t2284, t2278, t2256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2380, t2369, t2364, t2347, 0, 0, 0, 0, 0, 0, t2316, t2320, t2302, t2286, 0, 0, 0, 0, 0, 0, t2301, t2317, t2321, t2275, 0, 0, 0, 0, 0, 0, t2281, t2283, t2277, t2257, 0, 0, 0, 0, 0, 0, t2281, t2283, t2277, t2255; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2437, t2434, t2412, t2398, 0, 0, 0, 0, 0, 0, t2353, -t2361, t2346, t2313, 0, 0, 0, 0, 0, 0, t2345, -t2353, t2361, t2295, 0, 0, 0, 0, 0, 0, t2309, t2311, t2298, t2268, 0, 0, 0, 0, 0, 0, t2309, t2311, t2298, t2262; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2433, t2421, t2411, t2397, 0, 0, 0, 0, 0, 0, -t2372, -t2378, -t2540, -t2366, 0, 0, 0, 0, 0, 0, -t2540, t2374, t2376, -t2326, 0, 0, 0, 0, 0, 0, -t2330, -t2332, -t2315, -t2272, 0, 0, 0, 0, 0, 0, -t2330, -t2332, -t2315, -t2264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2413, t2410, t2438, t2416, 0, 0, 0, 0, 0, 0, -t2352, t2358, t2344, t2312, 0, 0, 0, 0, 0, 0, t2343, t2352, -t2358, t2294, 0, 0, 0, 0, 0, 0, t2308, t2310, t2297, t2267, 0, 0, 0, 0, 0, 0, t2308, t2310, t2297, t2261; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2396, t2482, t2375, t2338, 0, 0, 0, 0, 0, 0, -t2373, -t2396, -t2482, t2324, 0, 0, 0, 0, 0, 0, t2340, t2342, t2356, t2306, 0, 0, 0, 0, 0, 0, t2340, t2342, t2356, t2291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2388, t2514, t2379, t2337, 0, 0, 0, 0, 0, 0, t2379, t2388, -t2514, -t2325, 0, 0, 0, 0, 0, 0, -t2329, -t2331, -t2314, -t2271, 0, 0, 0, 0, 0, 0, -t2329, -t2331, -t2314, -t2263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2372, t2378, t2540, t2366, 0, 0, 0, 0, 0, 0, t2540, -t2374, -t2376, t2326, 0, 0, 0, 0, 0, 0, t2330, t2332, t2315, t2272, 0, 0, 0, 0, 0, 0, t2330, t2332, t2315, t2264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2540, -t2374, -t2376, t2326, 0, 0, 0, 0, 0, 0, t2330, t2332, t2315, t2272, 0, 0, 0, 0, 0, 0, t2330, t2332, t2315, t2264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2373, t2396, t2482, -t2324, 0, 0, 0, 0, 0, 0, -t2340, -t2342, -t2356, -t2306, 0, 0, 0, 0, 0, 0, -t2340, -t2342, -t2356, -t2291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2379, -t2388, t2514, t2325, 0, 0, 0, 0, 0, 0, t2329, t2331, t2314, t2271, 0, 0, 0, 0, 0, 0, t2329, t2331, t2314, t2263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2362, t2350, t2341, t2289, 0, 0, 0, 0, 0, 0, t2362, t2350, t2341, t2280; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2475, t2515, t2481, t2501, 0, 0, 0, 0, 0, 0, t2475, t2515, t2481, t2473; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2340, t2342, t2356, t2306, 0, 0, 0, 0, 0, 0, t2340, t2342, t2356, t2291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2362, t2350, t2341, t2280; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2475, t2515, t2481, t2473; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2340, t2342, t2356, t2291;];
f_new_reg  = t1;
