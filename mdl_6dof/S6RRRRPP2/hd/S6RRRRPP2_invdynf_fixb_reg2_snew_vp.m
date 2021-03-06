% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RRRRPP2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 18:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RRRRPP2_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP2_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP2_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_invdynf_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:06:13
% EndTime: 2019-05-07 18:06:21
% DurationCPUTime: 8.56s
% Computational Cost: add. (21096->304), mult. (42917->342), div. (0->0), fcn. (30826->8), ass. (0->198)
t2457 = sin(qJ(2));
t2461 = cos(qJ(2));
t2456 = sin(qJ(3));
t2460 = cos(qJ(3));
t2499 = t2456 * t2461;
t2481 = t2457 * t2460 + t2499;
t2432 = t2481 * qJD(1);
t2455 = sin(qJ(4));
t2459 = cos(qJ(4));
t2497 = qJD(2) + qJD(3);
t2413 = t2455 * t2432 - t2459 * t2497;
t2500 = t2456 * t2457;
t2480 = -t2460 * t2461 + t2500;
t2430 = qJD(1) * t2480;
t2511 = qJD(2) * t2457;
t2487 = qJD(1) * t2511;
t2492 = t2461 * qJDD(1);
t2478 = -t2487 + t2492;
t2486 = t2461 * qJD(1) * qJD(2);
t2493 = t2457 * qJDD(1);
t2479 = t2486 + t2493;
t2471 = -t2430 * qJD(3) + t2456 * t2478 + t2460 * t2479;
t2491 = qJDD(2) + qJDD(3);
t2468 = -t2455 * t2491 - t2459 * t2471;
t2466 = -t2413 * qJD(4) - t2468;
t2428 = qJD(4) + t2430;
t2502 = t2413 * t2428;
t2465 = t2466 - t2502;
t2408 = t2428 ^ 2;
t2415 = t2459 * t2432 + t2455 * t2497;
t2412 = t2415 ^ 2;
t2383 = t2412 + t2408;
t2392 = t2415 * t2413;
t2434 = t2460 * t2478;
t2475 = t2456 * t2479 - t2434;
t2473 = t2432 * qJD(3) + qJDD(4) + t2475;
t2523 = t2392 + t2473;
t2536 = -t2383 * t2455 + t2459 * t2523;
t2549 = -t2465 * t2456 + t2460 * t2536;
t2550 = t2456 * t2536 + t2460 * t2465;
t2285 = t2457 * t2550 - t2461 * t2549;
t2458 = sin(qJ(1));
t2462 = cos(qJ(1));
t2535 = t2383 * t2459 + t2455 * t2523;
t2569 = t2285 * t2458 + t2462 * t2535;
t2568 = t2285 * t2462 - t2458 * t2535;
t2517 = t2413 ^ 2;
t2521 = -t2517 - t2408;
t2524 = -t2392 + t2473;
t2540 = t2455 * t2521 + t2524 * t2459;
t2552 = t2462 * t2540;
t2470 = t2455 * t2471 - t2459 * t2491;
t2469 = t2415 * qJD(4) + t2470;
t2349 = t2428 * t2415 + t2469;
t2539 = -t2524 * t2455 + t2459 * t2521;
t2554 = t2460 * t2539;
t2548 = t2349 * t2456 + t2554;
t2557 = t2456 * t2539;
t2551 = -t2460 * t2349 + t2557;
t2561 = -t2457 * t2551 + t2461 * t2548;
t2567 = t2458 * t2561 - t2552;
t2555 = t2458 * t2540;
t2565 = t2462 * t2561 + t2555;
t2558 = t2457 * t2549 + t2461 * t2550;
t2351 = t2466 + t2502;
t2495 = t2428 - qJD(4);
t2467 = t2415 * t2495 - t2470;
t2543 = t2455 * t2467;
t2520 = -t2351 * t2459 + t2543;
t2542 = t2459 * t2467;
t2519 = t2351 * t2455 + t2542;
t2369 = t2412 + t2517;
t2544 = t2369 * t2460;
t2537 = t2456 * t2519 + t2544;
t2545 = t2369 * t2456;
t2538 = t2460 * t2519 - t2545;
t2546 = -t2457 * t2537 + t2461 * t2538;
t2563 = t2458 * t2546 - t2462 * t2520;
t2562 = t2458 * t2520 + t2462 * t2546;
t2560 = t2457 * t2548 + t2461 * t2551;
t2547 = t2457 * t2538 + t2461 * t2537;
t2453 = t2457 ^ 2;
t2454 = t2461 ^ 2;
t2522 = -t2454 - t2453;
t2541 = (pkin(8) * t2522 - pkin(7)) * qJD(1);
t2490 = t2497 ^ 2;
t2525 = qJD(4) + t2428;
t2463 = qJD(2) ^ 2;
t2518 = qJD(1) ^ 2;
t2446 = -t2518 * t2454 - t2463;
t2516 = t2430 ^ 2;
t2515 = t2432 ^ 2;
t2514 = t2461 * g(3);
t2513 = t2415 * qJ(6);
t2501 = t2430 * t2432;
t2496 = 0.2e1 * qJD(3) + qJD(2);
t2483 = t2462 * g(1) + t2458 * g(2);
t2476 = -pkin(1) * t2518 + qJDD(1) * pkin(7) - t2483;
t2422 = -t2457 * g(3) + t2461 * t2476;
t2390 = pkin(2) * t2446 + pkin(8) * t2492 + t2422;
t2488 = pkin(2) * t2461 + pkin(1);
t2472 = qJDD(2) * pkin(2) - t2514 + ((-pkin(7) - pkin(8)) * qJDD(1) + t2488 * t2518 + t2483) * t2457;
t2357 = t2460 * t2390 + t2456 * t2472;
t2407 = pkin(3) * t2430 - pkin(9) * t2432;
t2332 = -pkin(3) * t2490 + pkin(9) * t2491 - t2430 * t2407 + t2357;
t2443 = t2458 * g(1) - t2462 * g(2);
t2484 = t2496 * t2432;
t2464 = t2496 * t2430 * pkin(9) + (-t2434 + t2484) * pkin(3) + (pkin(3) * t2500 - pkin(9) * t2481 - t2488) * qJDD(1) + (t2541 + (0.2e1 * t2457 * pkin(2) + pkin(3) * t2499 + pkin(9) * t2480) * qJD(2)) * qJD(1) - t2443;
t2293 = t2459 * t2332 + t2455 * t2464;
t2485 = pkin(4) * t2428 - (2 * qJD(5));
t2292 = -t2455 * t2332 + t2459 * t2464;
t2356 = -t2456 * t2390 + t2460 * t2472;
t2387 = pkin(4) * t2413 - qJ(5) * t2415;
t2482 = t2473 * qJ(5) + 0.2e1 * qJD(5) * t2428 - t2413 * t2387 + t2293;
t2477 = -t2473 * pkin(4) - t2408 * qJ(5) + qJDD(5) - t2292;
t2331 = -t2491 * pkin(3) - t2490 * pkin(9) + t2432 * t2407 - t2356;
t2474 = t2469 * pkin(4) - qJ(5) * t2465 + t2331;
t2448 = t2461 * t2518 * t2457;
t2445 = -t2453 * t2518 - t2463;
t2442 = -qJDD(2) + t2448;
t2441 = qJDD(2) + t2448;
t2440 = t2522 * t2518;
t2439 = -qJDD(1) * t2458 - t2462 * t2518;
t2438 = qJDD(1) * t2462 - t2458 * t2518;
t2437 = t2522 * qJDD(1);
t2436 = -0.2e1 * t2487 + t2492;
t2435 = 0.2e1 * t2486 + t2493;
t2433 = qJDD(1) * pkin(1) + pkin(7) * t2518 + t2443;
t2421 = -t2457 * t2476 - t2514;
t2420 = -t2515 - t2490;
t2419 = t2442 * t2461 - t2445 * t2457;
t2418 = -t2441 * t2457 + t2446 * t2461;
t2417 = t2442 * t2457 + t2445 * t2461;
t2416 = t2441 * t2461 + t2446 * t2457;
t2406 = -t2491 - t2501;
t2405 = t2491 - t2501;
t2404 = -t2490 - t2516;
t2397 = t2488 * qJDD(1) + (-0.2e1 * pkin(2) * t2511 - t2541) * qJD(1) + t2443;
t2396 = -t2515 - t2516;
t2389 = -t2421 * t2457 + t2422 * t2461;
t2388 = t2421 * t2461 + t2422 * t2457;
t2381 = t2406 * t2460 - t2420 * t2456;
t2380 = t2406 * t2456 + t2420 * t2460;
t2379 = t2481 * qJDD(1);
t2378 = -t2430 * t2497 + t2471;
t2377 = -t2456 * t2493 + t2434 + (-qJD(1) * t2499 + t2432) * qJD(2);
t2376 = t2484 + t2475;
t2371 = t2404 * t2460 - t2405 * t2456;
t2370 = t2404 * t2456 + t2405 * t2460;
t2354 = -t2413 * t2495 + t2468;
t2347 = t2415 * t2525 + t2470;
t2344 = -t2380 * t2457 + t2381 * t2461;
t2343 = t2380 * t2461 + t2381 * t2457;
t2342 = t2377 * t2460 + t2379 * t2456;
t2341 = t2377 * t2456 - t2379 * t2460;
t2340 = -t2370 * t2457 + t2371 * t2461;
t2339 = t2370 * t2461 + t2371 * t2457;
t2321 = -t2356 * t2456 + t2357 * t2460;
t2320 = t2356 * t2460 + t2357 * t2456;
t2317 = t2354 * t2455 - t2542;
t2314 = -t2354 * t2459 - t2543;
t2313 = -t2341 * t2457 + t2342 * t2461;
t2312 = t2341 * t2461 + t2342 * t2457;
t2304 = t2347 * t2456 + t2554;
t2301 = -t2347 * t2460 + t2557;
t2297 = t2317 * t2460 + t2545;
t2294 = t2317 * t2456 - t2544;
t2291 = t2485 * t2415 + t2474;
t2290 = -t2320 * t2457 + t2321 * t2461;
t2289 = t2320 * t2461 + t2321 * t2457;
t2288 = t2387 * t2415 + t2477;
t2287 = -pkin(4) * t2408 + t2482;
t2279 = -t2301 * t2457 + t2304 * t2461;
t2276 = t2301 * t2461 + t2304 * t2457;
t2274 = t2470 * pkin(5) + t2517 * qJ(6) - qJDD(6) + t2474 + (pkin(5) * t2525 + t2485 + t2513) * t2415;
t2271 = -t2294 * t2457 + t2297 * t2461;
t2268 = t2294 * t2461 + t2297 * t2457;
t2267 = -t2292 * t2455 + t2293 * t2459;
t2266 = t2292 * t2459 + t2293 * t2455;
t2265 = -t2517 * pkin(5) + t2469 * qJ(6) + 0.2e1 * qJD(6) * t2413 + (-t2513 + (-pkin(4) - pkin(5)) * t2428) * t2428 + t2482;
t2264 = -t2473 * pkin(5) + (pkin(5) * t2413 - 0.2e1 * qJD(6) + t2387) * t2415 + t2477 - t2351 * qJ(6);
t2263 = t2267 * t2460 + t2331 * t2456;
t2262 = t2267 * t2456 - t2331 * t2460;
t2261 = t2287 * t2459 + t2288 * t2455;
t2260 = t2287 * t2455 - t2288 * t2459;
t2259 = t2264 * t2455 + t2265 * t2459;
t2258 = -t2264 * t2459 + t2265 * t2455;
t2257 = t2261 * t2460 + t2291 * t2456;
t2256 = t2261 * t2456 - t2291 * t2460;
t2255 = -t2262 * t2457 + t2263 * t2461;
t2254 = t2262 * t2461 + t2263 * t2457;
t2253 = t2259 * t2460 + t2274 * t2456;
t2252 = t2259 * t2456 - t2274 * t2460;
t2251 = -t2256 * t2457 + t2257 * t2461;
t2250 = t2256 * t2461 + t2257 * t2457;
t2249 = -t2252 * t2457 + t2253 * t2461;
t2248 = t2252 * t2461 + t2253 * t2457;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2439, -t2438, 0, -t2443 * t2458 - t2462 * t2483, 0, 0, 0, 0, 0, 0, t2418 * t2462 - t2436 * t2458, t2419 * t2462 + t2435 * t2458, -t2437 * t2462 + t2440 * t2458, t2389 * t2462 - t2433 * t2458, 0, 0, 0, 0, 0, 0, t2340 * t2462 + t2376 * t2458, t2344 * t2462 + t2378 * t2458, t2313 * t2462 + t2396 * t2458, t2290 * t2462 - t2397 * t2458, 0, 0, 0, 0, 0, 0, t2279 * t2462 + t2555, t2568, t2562, t2255 * t2462 + t2266 * t2458, 0, 0, 0, 0, 0, 0, t2565, t2562, -t2568, t2251 * t2462 + t2260 * t2458, 0, 0, 0, 0, 0, 0, t2565, -t2568, t2271 * t2462 + t2314 * t2458, t2249 * t2462 + t2258 * t2458; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2438, t2439, 0, t2443 * t2462 - t2458 * t2483, 0, 0, 0, 0, 0, 0, t2418 * t2458 + t2436 * t2462, t2419 * t2458 - t2435 * t2462, -t2437 * t2458 - t2440 * t2462, t2389 * t2458 + t2433 * t2462, 0, 0, 0, 0, 0, 0, t2340 * t2458 - t2376 * t2462, t2344 * t2458 - t2378 * t2462, t2313 * t2458 - t2396 * t2462, t2290 * t2458 + t2397 * t2462, 0, 0, 0, 0, 0, 0, t2279 * t2458 - t2552, t2569, t2563, t2255 * t2458 - t2266 * t2462, 0, 0, 0, 0, 0, 0, t2567, t2563, -t2569, t2251 * t2458 - t2260 * t2462, 0, 0, 0, 0, 0, 0, t2567, -t2569, t2271 * t2458 - t2314 * t2462, t2249 * t2458 - t2258 * t2462; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2416, t2417, 0, t2388, 0, 0, 0, 0, 0, 0, t2339, t2343, t2312, t2289, 0, 0, 0, 0, 0, 0, t2276, -t2558, t2547, t2254, 0, 0, 0, 0, 0, 0, t2560, t2547, t2558, t2250, 0, 0, 0, 0, 0, 0, t2560, t2558, t2268, t2248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2518, -qJDD(1), 0, -t2483, 0, 0, 0, 0, 0, 0, t2418, t2419, -t2437, t2389, 0, 0, 0, 0, 0, 0, t2340, t2344, t2313, t2290, 0, 0, 0, 0, 0, 0, t2279, t2285, t2546, t2255, 0, 0, 0, 0, 0, 0, t2561, t2546, -t2285, t2251, 0, 0, 0, 0, 0, 0, t2561, -t2285, t2271, t2249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2518, 0, t2443, 0, 0, 0, 0, 0, 0, t2436, -t2435, -t2440, t2433, 0, 0, 0, 0, 0, 0, -t2376, -t2378, -t2396, t2397, 0, 0, 0, 0, 0, 0, -t2540, t2535, -t2520, -t2266, 0, 0, 0, 0, 0, 0, -t2540, -t2520, -t2535, -t2260, 0, 0, 0, 0, 0, 0, -t2540, -t2535, -t2314, -t2258; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2416, t2417, 0, t2388, 0, 0, 0, 0, 0, 0, t2339, t2343, t2312, t2289, 0, 0, 0, 0, 0, 0, t2276, -t2558, t2547, t2254, 0, 0, 0, 0, 0, 0, t2560, t2547, t2558, t2250, 0, 0, 0, 0, 0, 0, t2560, t2558, t2268, t2248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2446, t2442, t2492, t2422, 0, 0, 0, 0, 0, 0, t2371, t2381, t2342, t2321, 0, 0, 0, 0, 0, 0, t2304, -t2549, t2538, t2263, 0, 0, 0, 0, 0, 0, t2548, t2538, t2549, t2257, 0, 0, 0, 0, 0, 0, t2548, t2549, t2297, t2253; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2441, t2445, -t2493, t2421, 0, 0, 0, 0, 0, 0, t2370, t2380, t2341, t2320, 0, 0, 0, 0, 0, 0, t2301, -t2550, t2537, t2262, 0, 0, 0, 0, 0, 0, t2551, t2537, t2550, t2256, 0, 0, 0, 0, 0, 0, t2551, t2550, t2294, t2252; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2436, t2435, t2440, -t2433, 0, 0, 0, 0, 0, 0, t2376, t2378, t2396, -t2397, 0, 0, 0, 0, 0, 0, t2540, -t2535, t2520, t2266, 0, 0, 0, 0, 0, 0, t2540, t2520, t2535, t2260, 0, 0, 0, 0, 0, 0, t2540, t2535, t2314, t2258; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2404, t2406, t2377, t2357, 0, 0, 0, 0, 0, 0, t2539, -t2536, t2519, t2267, 0, 0, 0, 0, 0, 0, t2539, t2519, t2536, t2261, 0, 0, 0, 0, 0, 0, t2539, t2536, t2317, t2259; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2405, t2420, -t2379, t2356, 0, 0, 0, 0, 0, 0, -t2347, -t2465, t2369, -t2331, 0, 0, 0, 0, 0, 0, -t2349, t2369, t2465, -t2291, 0, 0, 0, 0, 0, 0, -t2349, t2465, -t2369, -t2274; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2376, t2378, t2396, -t2397, 0, 0, 0, 0, 0, 0, t2540, -t2535, t2520, t2266, 0, 0, 0, 0, 0, 0, t2540, t2520, t2535, t2260, 0, 0, 0, 0, 0, 0, t2540, t2535, t2314, t2258; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2521, -t2523, t2467, t2293, 0, 0, 0, 0, 0, 0, t2521, t2467, t2523, t2287, 0, 0, 0, 0, 0, 0, t2521, t2523, -t2467, t2265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2524, -t2383, -t2351, t2292, 0, 0, 0, 0, 0, 0, t2524, -t2351, t2383, -t2288, 0, 0, 0, 0, 0, 0, t2524, t2383, -t2354, -t2264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2347, t2465, -t2369, t2331, 0, 0, 0, 0, 0, 0, t2349, -t2369, -t2465, t2291, 0, 0, 0, 0, 0, 0, t2349, -t2465, t2369, t2274; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2521, t2467, t2523, t2287, 0, 0, 0, 0, 0, 0, t2521, t2523, -t2467, t2265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2349, -t2369, -t2465, t2291, 0, 0, 0, 0, 0, 0, t2349, -t2465, t2369, t2274; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2524, t2351, -t2383, t2288, 0, 0, 0, 0, 0, 0, -t2524, -t2383, t2354, t2264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2521, t2523, -t2467, t2265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2524, -t2383, t2354, t2264; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2349, t2465, -t2369, -t2274;];
f_new_reg  = t1;
