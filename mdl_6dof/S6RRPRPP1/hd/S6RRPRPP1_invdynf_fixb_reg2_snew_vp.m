% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RRPRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 12:21
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RRPRPP1_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP1_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP1_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_invdynf_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 12:21:32
% EndTime: 2019-05-06 12:21:42
% DurationCPUTime: 11.07s
% Computational Cost: add. (43488->311), mult. (101475->392), div. (0->0), fcn. (73812->10), ass. (0->222)
t2599 = cos(qJ(2));
t2640 = qJD(1) * t2599;
t2623 = qJD(2) * t2640;
t2596 = sin(qJ(2));
t2625 = t2596 * qJDD(1);
t2566 = t2623 + t2625;
t2584 = t2599 * qJDD(1);
t2641 = qJD(1) * t2596;
t2624 = qJD(2) * t2641;
t2567 = t2584 - t2624;
t2591 = sin(pkin(9));
t2593 = cos(pkin(9));
t2617 = t2566 * t2591 - t2593 * t2567;
t2616 = qJDD(4) + t2617;
t2561 = t2591 * t2640 + t2593 * t2641;
t2595 = sin(qJ(4));
t2598 = cos(qJ(4));
t2544 = -t2598 * qJD(2) + t2561 * t2595;
t2545 = qJD(2) * t2595 + t2561 * t2598;
t2590 = sin(pkin(10));
t2592 = cos(pkin(10));
t2521 = t2544 * t2592 + t2590 * t2545;
t2523 = -t2544 * t2590 + t2545 * t2592;
t2635 = t2521 * t2523;
t2476 = t2616 + t2635;
t2519 = t2523 ^ 2;
t2559 = t2591 * t2641 - t2593 * t2640;
t2558 = qJD(4) + t2559;
t2646 = t2558 ^ 2;
t2653 = -t2519 - t2646;
t2436 = t2476 * t2590 - t2592 * t2653;
t2438 = t2476 * t2592 + t2590 * t2653;
t2425 = t2436 * t2595 - t2438 * t2598;
t2541 = t2566 * t2593 + t2567 * t2591;
t2614 = -t2595 * qJDD(2) - t2598 * t2541;
t2510 = -qJD(4) * t2544 - t2614;
t2615 = t2598 * qJDD(2) - t2595 * t2541;
t2609 = -qJD(4) * t2545 + t2615;
t2605 = t2592 * t2510 + t2590 * t2609;
t2634 = t2521 * t2558;
t2604 = t2605 - t2634;
t2401 = t2425 * t2591 - t2593 * t2604;
t2403 = t2425 * t2593 + t2591 * t2604;
t2376 = t2401 * t2596 - t2403 * t2599;
t2409 = t2436 * t2598 + t2438 * t2595;
t2597 = sin(qJ(1));
t2600 = cos(qJ(1));
t2683 = t2376 * t2597 - t2409 * t2600;
t2682 = t2376 * t2600 + t2409 * t2597;
t2378 = t2401 * t2599 + t2403 * t2596;
t2477 = t2616 - t2635;
t2488 = t2521 ^ 2;
t2652 = -t2646 - t2488;
t2660 = -t2477 * t2590 + t2592 * t2652;
t2661 = t2592 * t2477 + t2590 * t2652;
t2662 = t2595 * t2660 + t2598 * t2661;
t2619 = -t2510 * t2590 + t2592 * t2609;
t2632 = t2558 * t2523;
t2610 = -t2619 + t2632;
t2663 = -t2595 * t2661 + t2598 * t2660;
t2670 = t2591 * t2610 + t2593 * t2663;
t2671 = t2591 * t2663 - t2593 * t2610;
t2675 = -t2596 * t2671 + t2599 * t2670;
t2679 = t2597 * t2675 - t2600 * t2662;
t2678 = t2597 * t2662 + t2600 * t2675;
t2459 = t2605 + t2634;
t2611 = t2619 + t2632;
t2650 = t2459 * t2590 + t2592 * t2611;
t2651 = -t2592 * t2459 + t2590 * t2611;
t2658 = t2595 * t2650 + t2598 * t2651;
t2466 = t2519 + t2488;
t2659 = -t2595 * t2651 + t2598 * t2650;
t2664 = -t2466 * t2591 + t2593 * t2659;
t2665 = t2466 * t2593 + t2591 * t2659;
t2669 = -t2596 * t2665 + t2599 * t2664;
t2677 = t2597 * t2669 - t2600 * t2658;
t2676 = t2597 * t2658 + t2600 * t2669;
t2674 = t2596 * t2670 + t2599 * t2671;
t2668 = t2596 * t2664 + t2599 * t2665;
t2588 = t2599 ^ 2;
t2601 = qJD(1) ^ 2;
t2612 = qJD(2) * pkin(2) - qJ(3) * t2641;
t2575 = t2597 * g(1) - t2600 * g(2);
t2613 = qJDD(1) * pkin(1) + t2575;
t2527 = t2567 * pkin(2) + (qJ(3) * t2588 + pkin(7)) * t2601 - t2612 * t2641 - qJDD(3) + t2613;
t2649 = qJD(2) ^ 2;
t2648 = t2544 ^ 2;
t2647 = t2545 ^ 2;
t2645 = t2559 ^ 2;
t2644 = t2561 ^ 2;
t2643 = -2 * qJD(3);
t2642 = -2 * qJD(5);
t2639 = qJD(2) * t2559;
t2638 = qJD(2) * t2561;
t2633 = t2544 * t2545;
t2631 = t2558 * t2544;
t2630 = t2559 * t2561;
t2629 = t2588 * t2601;
t2628 = t2596 * t2601;
t2627 = qJD(4) - t2558;
t2576 = -g(1) * t2600 - g(2) * t2597;
t2608 = -pkin(1) * t2601 + qJDD(1) * pkin(7) + t2576;
t2551 = -t2596 * g(3) + t2599 * t2608;
t2524 = -pkin(2) * t2629 + t2567 * qJ(3) - qJD(2) * t2612 + t2551;
t2606 = t2596 * t2608;
t2602 = -t2606 - t2566 * qJ(3) + qJDD(2) * pkin(2) + (qJ(3) * qJD(1) * qJD(2) + pkin(2) * t2628 - g(3)) * t2599;
t2482 = t2593 * t2524 + t2559 * t2643 + t2591 * t2602;
t2535 = pkin(3) * t2559 - pkin(8) * t2561;
t2464 = -pkin(3) * t2649 + qJDD(2) * pkin(8) - t2535 * t2559 + t2482;
t2529 = t2617 + t2638;
t2621 = -t2541 + t2639;
t2472 = pkin(3) * t2529 + pkin(8) * t2621 - t2527;
t2433 = t2598 * t2464 + t2595 * t2472;
t2587 = t2596 ^ 2;
t2626 = t2587 + t2588;
t2533 = pkin(4) * t2558 - qJ(5) * t2545;
t2422 = -pkin(4) * t2648 + qJ(5) * t2609 - t2558 * t2533 + t2433;
t2432 = -t2595 * t2464 + t2598 * t2472;
t2497 = t2616 - t2633;
t2603 = (-t2510 - t2631) * qJ(5) + t2497 * pkin(4) + t2432;
t2385 = t2592 * t2422 + t2521 * t2642 + t2590 * t2603;
t2620 = t2590 * t2422 - t2592 * t2603;
t2618 = t2591 * t2524 - t2593 * t2602;
t2463 = -qJDD(2) * pkin(3) - t2649 * pkin(8) + ((2 * qJD(3)) + t2535) * t2561 + t2618;
t2431 = -t2609 * pkin(4) - t2648 * qJ(5) + t2545 * t2533 + qJDD(5) + t2463;
t2582 = t2599 * t2628;
t2581 = -t2629 - t2649;
t2580 = -t2587 * t2601 - t2649;
t2574 = -qJDD(2) + t2582;
t2573 = qJDD(2) + t2582;
t2572 = t2626 * t2601;
t2571 = -qJDD(1) * t2597 - t2600 * t2601;
t2570 = qJDD(1) * t2600 - t2597 * t2601;
t2569 = t2626 * qJDD(1);
t2568 = t2584 - 0.2e1 * t2624;
t2565 = 0.2e1 * t2623 + t2625;
t2563 = t2601 * pkin(7) + t2613;
t2552 = -t2644 - t2649;
t2550 = -t2599 * g(3) - t2606;
t2549 = t2574 * t2599 - t2580 * t2596;
t2548 = -t2573 * t2596 + t2581 * t2599;
t2547 = t2574 * t2596 + t2580 * t2599;
t2546 = t2573 * t2599 + t2581 * t2596;
t2539 = -qJDD(2) - t2630;
t2538 = qJDD(2) - t2630;
t2536 = -t2645 - t2649;
t2532 = -t2541 - t2639;
t2530 = -t2617 + t2638;
t2528 = -t2644 - t2645;
t2526 = -t2550 * t2596 + t2551 * t2599;
t2525 = t2550 * t2599 + t2551 * t2596;
t2513 = -t2646 - t2647;
t2512 = t2539 * t2593 - t2552 * t2591;
t2511 = t2539 * t2591 + t2552 * t2593;
t2508 = -t2646 - t2648;
t2502 = -t2647 - t2648;
t2500 = t2536 * t2593 - t2538 * t2591;
t2499 = t2536 * t2591 + t2538 * t2593;
t2498 = -t2616 - t2633;
t2494 = t2530 * t2593 - t2532 * t2591;
t2493 = t2530 * t2591 + t2532 * t2593;
t2492 = t2544 * t2627 + t2614;
t2491 = t2510 - t2631;
t2490 = -t2545 * t2627 + t2615;
t2489 = (qJD(4) + t2558) * t2545 - t2615;
t2487 = pkin(5) * t2521 - qJ(6) * t2523;
t2486 = -t2511 * t2596 + t2512 * t2599;
t2485 = t2511 * t2599 + t2512 * t2596;
t2481 = t2561 * t2643 - t2618;
t2475 = t2498 * t2598 - t2513 * t2595;
t2474 = t2498 * t2595 + t2513 * t2598;
t2470 = -t2497 * t2595 + t2508 * t2598;
t2469 = t2497 * t2598 + t2508 * t2595;
t2468 = -t2499 * t2596 + t2500 * t2599;
t2467 = t2499 * t2599 + t2500 * t2596;
t2454 = -t2493 * t2596 + t2494 * t2599;
t2453 = t2493 * t2599 + t2494 * t2596;
t2451 = t2490 * t2598 - t2492 * t2595;
t2450 = t2490 * t2595 + t2492 * t2598;
t2445 = t2475 * t2593 + t2491 * t2591;
t2444 = t2475 * t2591 - t2491 * t2593;
t2443 = t2470 * t2593 + t2489 * t2591;
t2442 = t2470 * t2591 - t2489 * t2593;
t2441 = -t2481 * t2591 + t2482 * t2593;
t2440 = t2481 * t2593 + t2482 * t2591;
t2435 = t2451 * t2593 + t2502 * t2591;
t2434 = t2451 * t2591 - t2502 * t2593;
t2420 = -t2444 * t2596 + t2445 * t2599;
t2419 = t2444 * t2599 + t2445 * t2596;
t2416 = -t2442 * t2596 + t2443 * t2599;
t2415 = t2442 * t2599 + t2443 * t2596;
t2414 = -t2440 * t2596 + t2441 * t2599;
t2413 = t2440 * t2599 + t2441 * t2596;
t2408 = -t2434 * t2596 + t2435 * t2599;
t2407 = t2434 * t2599 + t2435 * t2596;
t2406 = -t2432 * t2595 + t2433 * t2598;
t2405 = t2432 * t2598 + t2433 * t2595;
t2396 = -t2619 * pkin(5) + (pkin(5) * t2558 - (2 * qJD(6))) * t2523 + t2431 - t2604 * qJ(6);
t2391 = t2406 * t2593 + t2463 * t2591;
t2390 = t2406 * t2591 - t2463 * t2593;
t2384 = t2523 * t2642 - t2620;
t2383 = qJDD(6) - t2616 * pkin(5) - t2646 * qJ(6) + ((2 * qJD(5)) + t2487) * t2523 + t2620;
t2382 = -pkin(5) * t2646 + qJ(6) * t2616 + 0.2e1 * qJD(6) * t2558 - t2521 * t2487 + t2385;
t2373 = -t2390 * t2596 + t2391 * t2599;
t2372 = t2390 * t2599 + t2391 * t2596;
t2367 = -t2384 * t2590 + t2385 * t2592;
t2366 = t2384 * t2592 + t2385 * t2590;
t2365 = t2382 * t2592 + t2383 * t2590;
t2364 = t2382 * t2590 - t2383 * t2592;
t2363 = -t2366 * t2595 + t2367 * t2598;
t2362 = t2366 * t2598 + t2367 * t2595;
t2361 = t2363 * t2593 + t2431 * t2591;
t2360 = t2363 * t2591 - t2431 * t2593;
t2359 = -t2364 * t2595 + t2365 * t2598;
t2358 = t2364 * t2598 + t2365 * t2595;
t2357 = t2359 * t2593 + t2396 * t2591;
t2356 = t2359 * t2591 - t2396 * t2593;
t2355 = -t2360 * t2596 + t2361 * t2599;
t2354 = t2360 * t2599 + t2361 * t2596;
t2353 = -t2356 * t2596 + t2357 * t2599;
t2352 = t2356 * t2599 + t2357 * t2596;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2571, -t2570, 0, -t2575 * t2597 + t2576 * t2600, 0, 0, 0, 0, 0, 0, t2548 * t2600 - t2568 * t2597, t2549 * t2600 + t2565 * t2597, t2569 * t2600 - t2572 * t2597, t2526 * t2600 - t2563 * t2597, 0, 0, 0, 0, 0, 0, t2468 * t2600 + t2529 * t2597, t2486 * t2600 - t2597 * t2621, t2454 * t2600 + t2528 * t2597, t2414 * t2600 - t2527 * t2597, 0, 0, 0, 0, 0, 0, t2416 * t2600 + t2469 * t2597, t2420 * t2600 + t2474 * t2597, t2408 * t2600 + t2450 * t2597, t2373 * t2600 + t2405 * t2597, 0, 0, 0, 0, 0, 0, t2678, -t2682, t2676, t2355 * t2600 + t2362 * t2597, 0, 0, 0, 0, 0, 0, t2678, t2676, t2682, t2353 * t2600 + t2358 * t2597; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2570, t2571, 0, t2575 * t2600 + t2576 * t2597, 0, 0, 0, 0, 0, 0, t2548 * t2597 + t2568 * t2600, t2549 * t2597 - t2565 * t2600, t2569 * t2597 + t2572 * t2600, t2526 * t2597 + t2563 * t2600, 0, 0, 0, 0, 0, 0, t2468 * t2597 - t2529 * t2600, t2486 * t2597 + t2600 * t2621, t2454 * t2597 - t2528 * t2600, t2414 * t2597 + t2527 * t2600, 0, 0, 0, 0, 0, 0, t2416 * t2597 - t2469 * t2600, t2420 * t2597 - t2474 * t2600, t2408 * t2597 - t2450 * t2600, t2373 * t2597 - t2405 * t2600, 0, 0, 0, 0, 0, 0, t2679, -t2683, t2677, t2355 * t2597 - t2362 * t2600, 0, 0, 0, 0, 0, 0, t2679, t2677, t2683, t2353 * t2597 - t2358 * t2600; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2546, t2547, 0, t2525, 0, 0, 0, 0, 0, 0, t2467, t2485, t2453, t2413, 0, 0, 0, 0, 0, 0, t2415, t2419, t2407, t2372, 0, 0, 0, 0, 0, 0, t2674, t2378, t2668, t2354, 0, 0, 0, 0, 0, 0, t2674, t2668, -t2378, t2352; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2601, -qJDD(1), 0, t2576, 0, 0, 0, 0, 0, 0, t2548, t2549, t2569, t2526, 0, 0, 0, 0, 0, 0, t2468, t2486, t2454, t2414, 0, 0, 0, 0, 0, 0, t2416, t2420, t2408, t2373, 0, 0, 0, 0, 0, 0, t2675, -t2376, t2669, t2355, 0, 0, 0, 0, 0, 0, t2675, t2669, t2376, t2353; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2601, 0, t2575, 0, 0, 0, 0, 0, 0, t2568, -t2565, t2572, t2563, 0, 0, 0, 0, 0, 0, -t2529, t2621, -t2528, t2527, 0, 0, 0, 0, 0, 0, -t2469, -t2474, -t2450, -t2405, 0, 0, 0, 0, 0, 0, -t2662, t2409, -t2658, -t2362, 0, 0, 0, 0, 0, 0, -t2662, -t2658, -t2409, -t2358; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2546, t2547, 0, t2525, 0, 0, 0, 0, 0, 0, t2467, t2485, t2453, t2413, 0, 0, 0, 0, 0, 0, t2415, t2419, t2407, t2372, 0, 0, 0, 0, 0, 0, t2674, t2378, t2668, t2354, 0, 0, 0, 0, 0, 0, t2674, t2668, -t2378, t2352; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2581, t2574, t2584, t2551, 0, 0, 0, 0, 0, 0, t2500, t2512, t2494, t2441, 0, 0, 0, 0, 0, 0, t2443, t2445, t2435, t2391, 0, 0, 0, 0, 0, 0, t2670, t2403, t2664, t2361, 0, 0, 0, 0, 0, 0, t2670, t2664, -t2403, t2357; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2573, t2580, -t2625, t2550, 0, 0, 0, 0, 0, 0, t2499, t2511, t2493, t2440, 0, 0, 0, 0, 0, 0, t2442, t2444, t2434, t2390, 0, 0, 0, 0, 0, 0, t2671, t2401, t2665, t2360, 0, 0, 0, 0, 0, 0, t2671, t2665, -t2401, t2356; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2568, t2565, -t2572, -t2563, 0, 0, 0, 0, 0, 0, t2529, -t2621, t2528, -t2527, 0, 0, 0, 0, 0, 0, t2469, t2474, t2450, t2405, 0, 0, 0, 0, 0, 0, t2662, -t2409, t2658, t2362, 0, 0, 0, 0, 0, 0, t2662, t2658, t2409, t2358; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2536, t2539, t2530, t2482, 0, 0, 0, 0, 0, 0, t2470, t2475, t2451, t2406, 0, 0, 0, 0, 0, 0, t2663, t2425, t2659, t2363, 0, 0, 0, 0, 0, 0, t2663, t2659, -t2425, t2359; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2538, t2552, t2532, t2481, 0, 0, 0, 0, 0, 0, -t2489, -t2491, -t2502, -t2463, 0, 0, 0, 0, 0, 0, -t2610, -t2604, t2466, -t2431, 0, 0, 0, 0, 0, 0, -t2610, t2466, t2604, -t2396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2529, -t2621, t2528, -t2527, 0, 0, 0, 0, 0, 0, t2469, t2474, t2450, t2405, 0, 0, 0, 0, 0, 0, t2662, -t2409, t2658, t2362, 0, 0, 0, 0, 0, 0, t2662, t2658, t2409, t2358; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2508, t2498, t2490, t2433, 0, 0, 0, 0, 0, 0, t2660, -t2438, t2650, t2367, 0, 0, 0, 0, 0, 0, t2660, t2650, t2438, t2365; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2497, t2513, t2492, t2432, 0, 0, 0, 0, 0, 0, t2661, -t2436, t2651, t2366, 0, 0, 0, 0, 0, 0, t2661, t2651, t2436, t2364; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2489, t2491, t2502, t2463, 0, 0, 0, 0, 0, 0, t2610, t2604, -t2466, t2431, 0, 0, 0, 0, 0, 0, t2610, -t2466, -t2604, t2396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2652, -t2476, t2611, t2385, 0, 0, 0, 0, 0, 0, t2652, t2611, t2476, t2382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2477, t2653, -t2459, t2384, 0, 0, 0, 0, 0, 0, t2477, -t2459, -t2653, -t2383; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2610, t2604, -t2466, t2431, 0, 0, 0, 0, 0, 0, t2610, -t2466, -t2604, t2396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2652, t2611, t2476, t2382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2610, -t2466, -t2604, t2396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2477, t2459, t2653, t2383;];
f_new_reg  = t1;
