% Calculate inertial parameters regressor of inverse dynamics cutting forces vector with Newton-Euler for
% S6RRPPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
%
% Output:
% f_new_reg [(3*7)x(7*10)]
%   inertial parameter regressor of inverse dynamics cutting forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 10:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new_reg = S6RRPPRR4_invdynf_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_invdynf_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_invdynf_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRR4_invdynf_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR4_invdynf_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_invdynf_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_f_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:27:07
% EndTime: 2019-05-06 10:27:18
% DurationCPUTime: 12.68s
% Computational Cost: add. (49068->348), mult. (133015->471), div. (0->0), fcn. (103888->12), ass. (0->266)
t2688 = sin(qJ(2));
t2683 = sin(pkin(6));
t2692 = cos(qJ(2));
t2740 = t2683 * t2692;
t2733 = qJD(1) * t2740;
t2735 = qJDD(1) * t2683;
t2657 = qJD(2) * t2733 + t2688 * t2735;
t2741 = t2683 * t2688;
t2734 = qJD(1) * t2741;
t2658 = -qJD(2) * t2734 + t2692 * t2735;
t2682 = sin(pkin(11));
t2684 = cos(pkin(11));
t2615 = t2684 * t2657 + t2682 * t2658;
t2644 = t2682 * t2734 - t2684 * t2733;
t2685 = cos(pkin(6));
t2677 = qJD(1) * t2685 + qJD(2);
t2746 = t2644 * t2677;
t2591 = -t2615 + t2746;
t2676 = qJDD(1) * t2685 + qJDD(2);
t2646 = (t2682 * t2692 + t2684 * t2688) * qJD(1) * t2683;
t2745 = t2646 * t2644;
t2603 = t2676 + t2745;
t2612 = t2646 ^ 2;
t2675 = t2677 ^ 2;
t2730 = -t2675 - t2612;
t2571 = t2603 * t2682 - t2684 * t2730;
t2573 = t2603 * t2684 + t2682 * t2730;
t2709 = t2571 * t2692 + t2573 * t2688;
t2529 = -t2683 * t2591 + t2685 * t2709;
t2543 = t2571 * t2688 - t2573 * t2692;
t2689 = sin(qJ(1));
t2693 = cos(qJ(1));
t2787 = t2529 * t2689 + t2543 * t2693;
t2786 = t2529 * t2693 - t2543 * t2689;
t2614 = t2657 * t2682 - t2684 * t2658;
t2743 = t2677 * t2646;
t2587 = t2614 - t2743;
t2589 = t2615 + t2746;
t2770 = -t2587 * t2684 + t2589 * t2682;
t2771 = -t2587 * t2682 - t2589 * t2684;
t2775 = -t2688 * t2771 + t2692 * t2770;
t2643 = t2644 ^ 2;
t2760 = -t2612 - t2643;
t2774 = t2688 * t2770 + t2692 * t2771;
t2778 = -t2683 * t2760 + t2685 * t2774;
t2785 = -t2689 * t2778 + t2693 * t2775;
t2784 = t2689 * t2775 + t2693 * t2778;
t2586 = t2614 + t2743;
t2602 = -t2675 - t2643;
t2701 = t2676 - t2745;
t2560 = t2602 * t2682 + t2684 * t2701;
t2563 = -t2602 * t2684 + t2682 * t2701;
t2711 = t2560 * t2692 - t2563 * t2688;
t2524 = t2683 * t2586 - t2685 * t2711;
t2535 = t2560 * t2688 + t2563 * t2692;
t2781 = t2524 * t2689 - t2535 * t2693;
t2780 = t2524 * t2693 + t2535 * t2689;
t2779 = t2683 * t2774 + t2685 * t2760;
t2527 = t2685 * t2591 + t2683 * t2709;
t2521 = t2685 * t2586 + t2683 * t2711;
t2687 = sin(qJ(5));
t2691 = cos(qJ(5));
t2620 = -t2691 * t2644 + t2677 * t2687;
t2619 = qJD(6) + t2620;
t2761 = qJD(6) + t2619;
t2622 = t2644 * t2687 + t2677 * t2691;
t2641 = qJD(5) + t2646;
t2686 = sin(qJ(6));
t2690 = cos(qJ(6));
t2597 = t2622 * t2686 - t2690 * t2641;
t2759 = t2597 ^ 2;
t2599 = t2622 * t2690 + t2641 * t2686;
t2758 = t2599 ^ 2;
t2757 = t2619 ^ 2;
t2756 = t2620 ^ 2;
t2755 = t2622 ^ 2;
t2754 = t2641 ^ 2;
t2753 = 2 * qJD(4);
t2752 = -0.2e1 * t2646;
t2750 = t2597 * t2599;
t2747 = t2620 * t2622;
t2667 = t2689 * g(1) - g(2) * t2693;
t2694 = qJD(1) ^ 2;
t2653 = pkin(8) * t2683 * t2694 + qJDD(1) * pkin(1) + t2667;
t2744 = t2653 * t2685;
t2742 = t2683 ^ 2 * t2694;
t2608 = pkin(3) * t2644 - qJ(4) * t2646;
t2739 = (2 * qJD(3)) + t2608;
t2738 = qJD(5) - t2641;
t2737 = qJD(5) + t2641;
t2736 = qJD(6) - t2619;
t2623 = pkin(4) * t2646 - pkin(9) * t2677;
t2630 = -t2685 * g(3) - t2683 * t2653;
t2652 = pkin(2) * t2677 - qJ(3) * t2734;
t2681 = t2692 ^ 2;
t2732 = t2681 * t2742;
t2584 = -t2658 * pkin(2) - qJ(3) * t2732 + t2652 * t2734 + qJDD(3) + t2630;
t2695 = pkin(3) * t2743 + qJ(4) * t2591 + qJD(4) * t2752 + t2584;
t2510 = -t2643 * pkin(4) - t2646 * t2623 + (pkin(3) + pkin(9)) * t2614 + t2695;
t2668 = -g(1) * t2693 - g(2) * t2689;
t2654 = -pkin(1) * t2694 + pkin(8) * t2735 + t2668;
t2725 = -t2688 * t2654 + t2692 * t2744;
t2731 = t2688 * t2742;
t2568 = t2676 * pkin(2) - t2657 * qJ(3) + (pkin(2) * t2731 + (qJ(3) * qJD(1) * t2677 - g(3)) * t2683) * t2692 + t2725;
t2601 = -g(3) * t2741 + t2692 * t2654 + t2688 * t2744;
t2569 = -pkin(2) * t2732 + qJ(3) * t2658 - t2652 * t2677 + t2601;
t2728 = -t2684 * t2568 + t2682 * t2569;
t2700 = -t2676 * pkin(3) - t2675 * qJ(4) + qJDD(4) + t2728;
t2696 = -t2676 * pkin(9) + t2589 * pkin(4) + (pkin(9) * t2644 + t2739) * t2646 + t2700;
t2480 = t2691 * t2510 + t2687 * t2696;
t2533 = -0.2e1 * qJD(3) * t2644 + t2682 * t2568 + t2684 * t2569;
t2729 = qJDD(5) + t2615;
t2479 = -t2687 * t2510 + t2691 * t2696;
t2706 = t2687 * t2614 + t2691 * t2676;
t2578 = -qJD(5) * t2620 + t2706;
t2727 = -t2686 * t2578 + t2690 * t2729;
t2726 = -t2691 * t2614 + t2687 * t2676;
t2724 = t2677 * t2733;
t2585 = pkin(5) * t2620 - pkin(10) * t2622;
t2470 = -pkin(5) * t2754 + pkin(10) * t2729 - t2620 * t2585 + t2480;
t2697 = -t2675 * pkin(3) + t2676 * qJ(4) - t2644 * t2608 + t2533;
t2496 = -t2614 * pkin(4) - t2643 * pkin(9) + (t2753 + t2623) * t2677 + t2697;
t2546 = t2622 * t2737 + t2726;
t2482 = (t2620 * t2641 - t2578) * pkin(10) + t2546 * pkin(5) + t2496;
t2452 = -t2470 * t2686 + t2482 * t2690;
t2453 = t2470 * t2690 + t2482 * t2686;
t2442 = -t2452 * t2686 + t2453 * t2690;
t2469 = -pkin(5) * t2729 - pkin(10) * t2754 + t2622 * t2585 - t2479;
t2436 = t2442 * t2687 - t2469 * t2691;
t2441 = t2452 * t2690 + t2453 * t2686;
t2432 = -t2436 * t2684 + t2441 * t2682;
t2433 = t2436 * t2682 + t2441 * t2684;
t2723 = t2432 * t2692 + t2433 * t2688;
t2454 = t2479 * t2691 + t2480 * t2687;
t2450 = -t2454 * t2684 + t2496 * t2682;
t2451 = t2454 * t2682 + t2496 * t2684;
t2722 = t2450 * t2692 + t2451 * t2688;
t2516 = -t2599 * t2736 + t2727;
t2698 = -t2690 * t2578 - t2686 * t2729;
t2518 = t2597 * t2736 + t2698;
t2491 = t2516 * t2690 - t2518 * t2686;
t2544 = -t2758 - t2759;
t2477 = t2491 * t2687 - t2544 * t2691;
t2490 = t2516 * t2686 + t2518 * t2690;
t2458 = -t2477 * t2684 + t2490 * t2682;
t2459 = t2477 * t2682 + t2490 * t2684;
t2721 = t2458 * t2692 + t2459 * t2688;
t2699 = -qJD(5) * t2622 - qJDD(6) - t2726;
t2536 = -t2699 - t2750;
t2545 = -t2757 - t2759;
t2500 = -t2536 * t2686 + t2545 * t2690;
t2515 = t2599 * t2761 - t2727;
t2484 = t2500 * t2687 - t2515 * t2691;
t2499 = t2536 * t2690 + t2545 * t2686;
t2463 = -t2484 * t2684 + t2499 * t2682;
t2464 = t2484 * t2682 + t2499 * t2684;
t2720 = t2463 * t2692 + t2464 * t2688;
t2537 = t2699 - t2750;
t2556 = -t2757 - t2758;
t2509 = t2537 * t2690 - t2556 * t2686;
t2517 = -t2597 * t2761 - t2698;
t2486 = t2509 * t2687 - t2517 * t2691;
t2508 = t2537 * t2686 + t2556 * t2690;
t2467 = -t2486 * t2684 + t2508 * t2682;
t2468 = t2486 * t2682 + t2508 * t2684;
t2719 = t2467 * t2692 + t2468 * t2688;
t2513 = t2677 * t2753 + t2697;
t2514 = t2739 * t2646 + t2700;
t2488 = t2513 * t2682 - t2514 * t2684;
t2489 = t2513 * t2684 + t2514 * t2682;
t2718 = t2488 * t2692 + t2489 * t2688;
t2532 = qJD(3) * t2752 - t2728;
t2492 = t2532 * t2684 + t2533 * t2682;
t2493 = -t2532 * t2682 + t2533 * t2684;
t2717 = t2492 * t2692 + t2493 * t2688;
t2547 = -t2622 * t2738 - t2726;
t2549 = t2620 * t2738 - t2706;
t2519 = t2547 * t2687 + t2549 * t2691;
t2567 = -t2755 - t2756;
t2497 = -t2519 * t2684 + t2567 * t2682;
t2498 = t2519 * t2682 + t2567 * t2684;
t2716 = t2497 * t2692 + t2498 * t2688;
t2557 = t2729 - t2747;
t2575 = -t2754 - t2756;
t2538 = t2557 * t2691 + t2575 * t2687;
t2505 = -t2538 * t2684 + t2546 * t2682;
t2506 = t2538 * t2682 + t2546 * t2684;
t2715 = t2505 * t2692 + t2506 * t2688;
t2558 = -t2729 - t2747;
t2581 = -t2754 - t2755;
t2540 = t2558 * t2687 + t2581 * t2691;
t2548 = -t2620 * t2737 + t2706;
t2511 = -t2540 * t2684 + t2548 * t2682;
t2512 = t2540 * t2682 + t2548 * t2684;
t2714 = t2511 * t2692 + t2512 * t2688;
t2600 = -g(3) * t2740 + t2725;
t2707 = t2600 * t2692 + t2601 * t2688;
t2661 = t2677 * t2734;
t2624 = t2658 + t2661;
t2626 = t2724 - t2657;
t2705 = t2624 * t2688 + t2626 * t2692;
t2680 = t2688 ^ 2;
t2642 = -t2680 * t2742 - t2675;
t2665 = t2692 * t2731;
t2656 = t2665 - t2676;
t2704 = t2642 * t2692 + t2656 * t2688;
t2655 = t2665 + t2676;
t2659 = -t2675 - t2732;
t2703 = t2655 * t2692 + t2659 * t2688;
t2664 = -qJDD(1) * t2689 - t2693 * t2694;
t2663 = qJDD(1) * t2693 - t2689 * t2694;
t2660 = (-t2680 - t2681) * t2742;
t2627 = -t2658 + t2661;
t2625 = t2724 + t2657;
t2616 = -t2655 * t2688 + t2659 * t2692;
t2607 = -t2642 * t2688 + t2656 * t2692;
t2594 = t2624 * t2692 - t2626 * t2688;
t2583 = -t2683 * t2627 + t2685 * t2703;
t2582 = t2685 * t2627 + t2683 * t2703;
t2580 = -t2683 * t2625 + t2685 * t2704;
t2579 = t2685 * t2625 + t2683 * t2704;
t2577 = -t2683 * t2660 + t2685 * t2705;
t2576 = t2685 * t2660 + t2683 * t2705;
t2559 = -t2600 * t2688 + t2601 * t2692;
t2555 = -t2683 * t2630 + t2685 * t2707;
t2554 = t2685 * t2630 + t2683 * t2707;
t2541 = t2558 * t2691 - t2581 * t2687;
t2539 = -t2557 * t2687 + t2575 * t2691;
t2531 = t2614 * pkin(3) + t2695;
t2520 = t2547 * t2691 - t2549 * t2687;
t2487 = t2509 * t2691 + t2517 * t2687;
t2485 = t2500 * t2691 + t2515 * t2687;
t2483 = -t2511 * t2688 + t2512 * t2692;
t2481 = -t2505 * t2688 + t2506 * t2692;
t2478 = t2491 * t2691 + t2544 * t2687;
t2476 = -t2497 * t2688 + t2498 * t2692;
t2475 = -t2683 * t2541 + t2685 * t2714;
t2474 = t2685 * t2541 + t2683 * t2714;
t2473 = -t2683 * t2539 + t2685 * t2715;
t2472 = t2685 * t2539 + t2683 * t2715;
t2471 = -t2492 * t2688 + t2493 * t2692;
t2466 = -t2683 * t2584 + t2685 * t2717;
t2465 = t2685 * t2584 + t2683 * t2717;
t2462 = -t2683 * t2520 + t2685 * t2716;
t2461 = t2685 * t2520 + t2683 * t2716;
t2460 = -t2488 * t2688 + t2489 * t2692;
t2457 = -t2683 * t2531 + t2685 * t2718;
t2456 = t2685 * t2531 + t2683 * t2718;
t2455 = -t2479 * t2687 + t2480 * t2691;
t2449 = -t2467 * t2688 + t2468 * t2692;
t2448 = -t2463 * t2688 + t2464 * t2692;
t2447 = -t2683 * t2487 + t2685 * t2719;
t2446 = t2685 * t2487 + t2683 * t2719;
t2445 = -t2458 * t2688 + t2459 * t2692;
t2444 = -t2683 * t2485 + t2685 * t2720;
t2443 = t2685 * t2485 + t2683 * t2720;
t2440 = -t2683 * t2478 + t2685 * t2721;
t2439 = t2685 * t2478 + t2683 * t2721;
t2438 = -t2450 * t2688 + t2451 * t2692;
t2437 = t2442 * t2691 + t2469 * t2687;
t2435 = -t2683 * t2455 + t2685 * t2722;
t2434 = t2685 * t2455 + t2683 * t2722;
t2431 = -t2432 * t2688 + t2433 * t2692;
t2430 = -t2683 * t2437 + t2685 * t2723;
t2429 = t2685 * t2437 + t2683 * t2723;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, t2664, -t2663, 0, -t2667 * t2689 + t2668 * t2693, 0, 0, 0, 0, 0, 0, -t2583 * t2689 + t2616 * t2693, -t2580 * t2689 + t2607 * t2693, -t2577 * t2689 + t2594 * t2693, -t2555 * t2689 + t2559 * t2693, 0, 0, 0, 0, 0, 0, t2781, t2787, t2785, -t2466 * t2689 + t2471 * t2693, 0, 0, 0, 0, 0, 0, t2785, -t2781, -t2787, -t2457 * t2689 + t2460 * t2693, 0, 0, 0, 0, 0, 0, -t2473 * t2689 + t2481 * t2693, -t2475 * t2689 + t2483 * t2693, -t2462 * t2689 + t2476 * t2693, -t2435 * t2689 + t2438 * t2693, 0, 0, 0, 0, 0, 0, -t2444 * t2689 + t2448 * t2693, -t2447 * t2689 + t2449 * t2693, -t2440 * t2689 + t2445 * t2693, -t2430 * t2689 + t2431 * t2693; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, t2663, t2664, 0, t2667 * t2693 + t2668 * t2689, 0, 0, 0, 0, 0, 0, t2583 * t2693 + t2616 * t2689, t2580 * t2693 + t2607 * t2689, t2577 * t2693 + t2594 * t2689, t2555 * t2693 + t2559 * t2689, 0, 0, 0, 0, 0, 0, -t2780, -t2786, t2784, t2466 * t2693 + t2471 * t2689, 0, 0, 0, 0, 0, 0, t2784, t2780, t2786, t2457 * t2693 + t2460 * t2689, 0, 0, 0, 0, 0, 0, t2473 * t2693 + t2481 * t2689, t2475 * t2693 + t2483 * t2689, t2462 * t2693 + t2476 * t2689, t2435 * t2693 + t2438 * t2689, 0, 0, 0, 0, 0, 0, t2444 * t2693 + t2448 * t2689, t2447 * t2693 + t2449 * t2689, t2440 * t2693 + t2445 * t2689, t2430 * t2693 + t2431 * t2689; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2582, t2579, t2576, t2554, 0, 0, 0, 0, 0, 0, t2521, -t2527, t2779, t2465, 0, 0, 0, 0, 0, 0, t2779, -t2521, t2527, t2456, 0, 0, 0, 0, 0, 0, t2472, t2474, t2461, t2434, 0, 0, 0, 0, 0, 0, t2443, t2446, t2439, t2429; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2694, -qJDD(1), 0, t2668, 0, 0, 0, 0, 0, 0, t2616, t2607, t2594, t2559, 0, 0, 0, 0, 0, 0, -t2535, t2543, t2775, t2471, 0, 0, 0, 0, 0, 0, t2775, t2535, -t2543, t2460, 0, 0, 0, 0, 0, 0, t2481, t2483, t2476, t2438, 0, 0, 0, 0, 0, 0, t2448, t2449, t2445, t2431; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t2694, 0, t2667, 0, 0, 0, 0, 0, 0, t2583, t2580, t2577, t2555, 0, 0, 0, 0, 0, 0, -t2524, -t2529, t2778, t2466, 0, 0, 0, 0, 0, 0, t2778, t2524, t2529, t2457, 0, 0, 0, 0, 0, 0, t2473, t2475, t2462, t2435, 0, 0, 0, 0, 0, 0, t2444, t2447, t2440, t2430; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, t2582, t2579, t2576, t2554, 0, 0, 0, 0, 0, 0, t2521, -t2527, t2779, t2465, 0, 0, 0, 0, 0, 0, t2779, -t2521, t2527, t2456, 0, 0, 0, 0, 0, 0, t2472, t2474, t2461, t2434, 0, 0, 0, 0, 0, 0, t2443, t2446, t2439, t2429; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2659, t2656, t2624, t2601, 0, 0, 0, 0, 0, 0, -t2563, -t2573, t2770, t2493, 0, 0, 0, 0, 0, 0, t2770, t2563, t2573, t2489, 0, 0, 0, 0, 0, 0, t2506, t2512, t2498, t2451, 0, 0, 0, 0, 0, 0, t2464, t2468, t2459, t2433; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2655, t2642, t2626, t2600, 0, 0, 0, 0, 0, 0, t2560, -t2571, t2771, t2492, 0, 0, 0, 0, 0, 0, t2771, -t2560, t2571, t2488, 0, 0, 0, 0, 0, 0, t2505, t2511, t2497, t2450, 0, 0, 0, 0, 0, 0, t2463, t2467, t2458, t2432; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2627, t2625, t2660, t2630, 0, 0, 0, 0, 0, 0, t2586, -t2591, t2760, t2584, 0, 0, 0, 0, 0, 0, t2760, -t2586, t2591, t2531, 0, 0, 0, 0, 0, 0, t2539, t2541, t2520, t2455, 0, 0, 0, 0, 0, 0, t2485, t2487, t2478, t2437; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2602, -t2603, -t2587, t2533, 0, 0, 0, 0, 0, 0, -t2587, -t2602, t2603, t2513, 0, 0, 0, 0, 0, 0, t2546, t2548, t2567, t2496, 0, 0, 0, 0, 0, 0, t2499, t2508, t2490, t2441; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2701, t2730, -t2589, t2532, 0, 0, 0, 0, 0, 0, -t2589, -t2701, -t2730, -t2514, 0, 0, 0, 0, 0, 0, -t2538, -t2540, -t2519, -t2454, 0, 0, 0, 0, 0, 0, -t2484, -t2486, -t2477, -t2436; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2586, -t2591, t2760, t2584, 0, 0, 0, 0, 0, 0, t2760, -t2586, t2591, t2531, 0, 0, 0, 0, 0, 0, t2539, t2541, t2520, t2455, 0, 0, 0, 0, 0, 0, t2485, t2487, t2478, t2437; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2760, -t2586, t2591, t2531, 0, 0, 0, 0, 0, 0, t2539, t2541, t2520, t2455, 0, 0, 0, 0, 0, 0, t2485, t2487, t2478, t2437; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2587, t2602, -t2603, -t2513, 0, 0, 0, 0, 0, 0, -t2546, -t2548, -t2567, -t2496, 0, 0, 0, 0, 0, 0, -t2499, -t2508, -t2490, -t2441; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2589, t2701, t2730, t2514, 0, 0, 0, 0, 0, 0, t2538, t2540, t2519, t2454, 0, 0, 0, 0, 0, 0, t2484, t2486, t2477, t2436; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2575, t2558, t2547, t2480, 0, 0, 0, 0, 0, 0, t2500, t2509, t2491, t2442; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2557, t2581, t2549, t2479, 0, 0, 0, 0, 0, 0, -t2515, -t2517, -t2544, -t2469; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2546, t2548, t2567, t2496, 0, 0, 0, 0, 0, 0, t2499, t2508, t2490, t2441; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2545, t2537, t2516, t2453; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2536, t2556, t2518, t2452; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2515, t2517, t2544, t2469;];
f_new_reg  = t1;