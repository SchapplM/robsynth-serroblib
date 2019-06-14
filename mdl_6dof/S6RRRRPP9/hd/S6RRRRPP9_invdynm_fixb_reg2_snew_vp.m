% Calculate inertial parameters regressor of inverse dynamics cutting torque vector with Newton-Euler for
% S6RRRRPP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
%
% Output:
% m_new_reg [(3*7)x(%Nl%*10)]
%   inertial parameter regressor of inverse dynamics cutting torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 19:24
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new_reg = S6RRRRPP9_invdynm_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP9_invdynm_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP9_invdynm_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP9_invdynm_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP9_invdynm_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP9_invdynm_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_m_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 19:22:57
% EndTime: 2019-05-07 19:24:16
% DurationCPUTime: 84.98s
% Computational Cost: add. (277301->1032), mult. (592807->1261), div. (0->0), fcn. (469010->10), ass. (0->686)
t2231 = sin(qJ(1));
t2235 = cos(qJ(1));
t2230 = sin(qJ(2));
t2234 = cos(qJ(2));
t2228 = sin(qJ(4));
t2232 = cos(qJ(4));
t2229 = sin(qJ(3));
t2233 = cos(qJ(3));
t2227 = cos(pkin(6));
t2443 = qJD(1) * t2227;
t2361 = qJD(2) + t2443;
t2226 = sin(pkin(6));
t2392 = t2226 * t2230;
t2375 = qJD(1) * t2392;
t2191 = t2229 * t2361 + t2233 * t2375;
t2440 = qJD(2) * qJD(1);
t2365 = t2234 * t2440;
t2379 = t2230 * qJDD(1);
t2197 = (t2365 + t2379) * t2226;
t2355 = qJDD(1) * t2227 + qJDD(2);
t2357 = -t2229 * t2197 + t2233 * t2355;
t2284 = t2191 * qJD(3) - t2357;
t2135 = qJDD(4) + t2284;
t2391 = t2226 * t2234;
t2219 = qJD(1) * t2391;
t2382 = t2219 - qJD(3);
t2160 = t2191 * t2228 + t2232 * t2382;
t2162 = t2232 * t2191 - t2228 * t2382;
t2402 = t2162 * t2160;
t2492 = -t2402 + t2135;
t2415 = t2492 * t2232;
t2189 = t2229 * t2375 - t2233 * t2361;
t2185 = qJD(4) + t2189;
t2183 = t2185 ^ 2;
t2467 = t2160 ^ 2;
t2494 = -t2183 - t2467;
t2537 = t2494 * t2228 + t2415;
t2261 = -t2233 * t2197 - t2229 * t2355;
t2137 = -t2189 * qJD(3) - t2261;
t2380 = qJDD(1) * t2234;
t2283 = t2230 * t2440 - t2380;
t2265 = t2283 * t2226;
t2257 = qJDD(3) + t2265;
t2054 = qJD(4) * t2162 + t2137 * t2228 - t2232 * t2257;
t2401 = t2185 * t2162;
t2496 = t2054 + t2401;
t2416 = t2492 * t2228;
t2538 = t2494 * t2232 - t2416;
t2561 = t2229 * t2496 + t2233 * t2538;
t2592 = t2230 * t2537 + t2234 * t2561;
t2566 = t2229 * t2538 - t2233 * t2496;
t2596 = t2230 * t2561 - t2234 * t2537;
t2620 = -t2226 * t2566 + t2227 * t2596;
t2647 = pkin(7) * (t2231 * t2592 + t2235 * t2620);
t2649 = pkin(7) * (t2231 * t2620 - t2235 * t2592);
t2241 = t2232 * t2137 + t2228 * t2257;
t2021 = (qJD(4) + t2185) * t2160 - t2241;
t2466 = t2162 ^ 2;
t2093 = -t2466 - t2183;
t2491 = t2402 + t2135;
t2518 = t2232 * t2491;
t2534 = t2093 * t2228 + t2518;
t2582 = t2229 * t2534;
t1871 = t2021 * t2233 - t2582;
t2576 = t2233 * t2534;
t1873 = t2021 * t2229 + t2576;
t2516 = t2491 * t2228;
t1957 = t2093 * t2232 - t2516;
t2609 = t1957 * t2234;
t2320 = t1873 * t2230 + t2609;
t1766 = t2226 * t1871 + t2227 * t2320;
t2610 = t1957 * t2230;
t1818 = t1873 * t2234 - t2610;
t2697 = pkin(7) * (t1766 * t2235 + t2231 * t1818);
t2696 = pkin(7) * (t2231 * t1766 - t1818 * t2235);
t1763 = -t2227 * t1871 + t2226 * t2320;
t2695 = pkin(8) * (t1763 * t2226 + t1766 * t2227);
t2624 = t2226 * t2596 + t2227 * t2566;
t2645 = (t2226 * t2624 + t2227 * t2620) * pkin(8);
t2694 = pkin(1) * t1763;
t2693 = pkin(1) * t1766;
t2659 = pkin(1) * t2620;
t2657 = pkin(1) * t2624;
t2119 = t2183 - t2466;
t1978 = t2119 * t2228 - t2415;
t2055 = -t2160 * qJD(4) + t2241;
t2403 = t2160 * t2185;
t2489 = t2403 + t2055;
t1880 = t1978 * t2229 + t2489 * t2233;
t1886 = t1978 * t2233 - t2489 * t2229;
t1972 = -t2119 * t2232 - t2416;
t2314 = t1886 * t2230 - t1972 * t2234;
t1779 = -t2226 * t1880 + t2227 * t2314;
t1823 = t1886 * t2234 + t1972 * t2230;
t2653 = t2231 * t1779 - t1823 * t2235;
t2690 = t1779 * t2235 + t2231 * t1823;
t2121 = t2467 - t2183;
t1981 = t2121 * t2232 - t2516;
t2497 = t2054 - t2401;
t1883 = t1981 * t2229 + t2497 * t2233;
t1889 = t1981 * t2233 - t2497 * t2229;
t1975 = t2121 * t2228 + t2518;
t2311 = t1889 * t2230 - t1975 * t2234;
t1782 = -t2226 * t1883 + t2227 * t2311;
t1826 = t1889 * t2234 + t1975 * t2230;
t2654 = t2231 * t1782 - t1826 * t2235;
t2689 = t1782 * t2235 + t2231 * t1826;
t2490 = -t2403 + t2055;
t2515 = t2496 * t2232;
t1896 = t2490 * t2228 + t2515;
t2101 = t2467 - t2466;
t1854 = t1896 * t2229 - t2101 * t2233;
t1857 = t1896 * t2233 + t2101 * t2229;
t2536 = -t2496 * t2228 + t2490 * t2232;
t2326 = t1857 * t2230 + t2234 * t2536;
t1744 = -t2226 * t1854 + t2227 * t2326;
t1793 = t1857 * t2234 - t2230 * t2536;
t2684 = t2231 * t1744 - t1793 * t2235;
t2683 = t1744 * t2235 + t2231 * t1793;
t2682 = pkin(8) * t1818;
t2643 = pkin(8) * t2592;
t1776 = t2227 * t1883 + t2226 * t2311;
t1773 = t2227 * t1880 + t2226 * t2314;
t1741 = t2227 * t1854 + t2226 * t2326;
t2671 = pkin(2) * t1871;
t2616 = pkin(2) * t2566;
t2670 = pkin(9) * t1871;
t2612 = pkin(9) * t2566;
t2618 = pkin(2) * t1957;
t2661 = pkin(9) * t1873 + t2618;
t2598 = -pkin(2) * t2537 + pkin(9) * t2561;
t2059 = t2466 + t2467;
t2543 = t2489 * t2228 - t2232 * t2497;
t2565 = t2059 * t2233 + t2229 * t2543;
t2535 = -t2497 * t2228 - t2489 * t2232;
t2560 = -t2059 * t2229 + t2233 * t2543;
t2597 = t2230 * t2560 - t2234 * t2535;
t2619 = -t2226 * t2565 + t2227 * t2597;
t2660 = pkin(1) * t2619;
t2623 = t2226 * t2597 + t2227 * t2565;
t2658 = pkin(1) * t2623;
t2593 = t2230 * t2535 + t2234 * t2560;
t2650 = pkin(7) * (-t2231 * t2619 + t2235 * t2593);
t2648 = pkin(7) * (t2231 * t2593 + t2235 * t2619);
t2646 = (-t2226 * t2623 - t2227 * t2619) * pkin(8);
t2586 = pkin(3) * t2537;
t2642 = pkin(8) * t2593;
t2583 = pkin(10) * t2537;
t2641 = pkin(10) * t2538;
t2617 = pkin(2) * t2565;
t2615 = pkin(3) * t1957;
t2614 = pkin(9) * t2560;
t2613 = pkin(9) * t2565;
t2611 = pkin(10) * t1957;
t2585 = pkin(10) * t2534;
t2438 = qJD(5) * t2185;
t2178 = -0.2e1 * t2438;
t2600 = t2178 + t2615;
t2599 = -pkin(2) * t2535 + t2614;
t2587 = pkin(3) * t2535;
t2584 = pkin(10) * t2535;
t2570 = -pkin(3) * t2496 + t2641;
t2569 = pkin(3) * t2059 + pkin(10) * t2543;
t2400 = t2185 * t2228;
t2116 = t2162 * t2400;
t2399 = t2185 * t2232;
t2371 = t2160 * t2399;
t2340 = t2116 - t2371;
t2404 = t2135 * t2229;
t2476 = t2233 * t2340 + t2404;
t2498 = (t2160 * t2228 + t2162 * t2232) * t2185;
t2508 = -t2230 * t2498 + t2234 * t2476;
t2131 = t2233 * t2135;
t2480 = t2229 * t2340 - t2131;
t2507 = t2230 * t2476 + t2234 * t2498;
t2532 = -t2226 * t2480 + t2227 * t2507;
t2564 = -t2231 * t2532 + t2235 * t2508;
t2559 = t2231 * t2508 + t2235 * t2532;
t2557 = pkin(4) * (t2496 + t2401);
t2556 = qJ(5) * t2059;
t2374 = t2229 * t2402;
t2383 = t2232 * t2055 - t2116;
t2475 = t2233 * t2383 + t2374;
t2117 = t2162 * t2399;
t2495 = t2055 * t2228 + t2117;
t2509 = t2230 * t2495 + t2234 * t2475;
t2548 = t2231 * t2509;
t2372 = t2160 * t2400;
t2002 = -t2054 * t2232 + t2372;
t2287 = t2054 * t2228 + t2371;
t2477 = t2233 * t2287 - t2374;
t2510 = t2002 * t2230 + t2234 * t2477;
t2547 = t2231 * t2510;
t2545 = t2235 * t2509;
t2544 = t2235 * t2510;
t2533 = t2226 * t2507 + t2227 * t2480;
t2531 = -2 * qJD(6);
t2529 = qJ(5) * t2490;
t2528 = qJ(5) * t2494;
t2527 = qJ(5) * t2496;
t2526 = qJ(6) * t2492;
t2373 = t2233 * t2402;
t2478 = t2229 * t2383 - t2373;
t2525 = t2226 * t2478;
t2479 = t2229 * t2287 + t2373;
t2524 = t2226 * t2479;
t2522 = t2227 * t2478;
t2521 = t2227 * t2479;
t2216 = g(1) * t2235 + t2231 * g(2);
t2451 = pkin(8) * t2226;
t2468 = qJD(1) ^ 2;
t2193 = -pkin(1) * t2468 + qJDD(1) * t2451 - t2216;
t2462 = pkin(2) * t2234;
t2353 = -pkin(9) * t2230 - t2462;
t2444 = qJD(1) * t2226;
t2196 = t2353 * t2444;
t2215 = t2231 * g(1) - t2235 * g(2);
t2259 = qJDD(1) * pkin(1) + t2451 * t2468 + t2215;
t2251 = t2227 * t2259;
t2245 = -g(3) * t2392 + t2230 * t2251;
t2354 = t2361 ^ 2;
t2077 = t2355 * pkin(9) - t2354 * pkin(2) + (t2196 * t2444 + t2193) * t2234 + t2245;
t2450 = pkin(9) * t2234;
t2463 = pkin(2) * t2230;
t2352 = -t2450 + t2463;
t2448 = t2227 * g(3);
t2237 = -t2197 * pkin(9) - t2448 + ((-pkin(1) - t2462) * qJDD(1) + ((t2227 * t2352 - t2451) * qJD(1) + (-t2450 + 0.2e1 * t2463) * qJD(2)) * qJD(1) - t2215) * t2226;
t1984 = t2233 * t2077 + t2229 * t2237;
t2146 = pkin(3) * t2189 - pkin(10) * t2191;
t2378 = t2382 ^ 2;
t1922 = -pkin(3) * t2378 + pkin(10) * t2257 - t2189 * t2146 + t1984;
t2358 = t2230 * t2193 - t2234 * t2251;
t2442 = qJD(1) * t2230;
t2076 = -t2355 * pkin(2) - t2354 * pkin(9) + (g(3) * t2234 + t2196 * t2442) * t2226 + t2358;
t2173 = t2382 * t2189;
t2088 = t2173 + t2137;
t2174 = t2382 * t2191;
t2513 = t2284 - t2174;
t1927 = pkin(3) * t2513 - t2088 * pkin(10) + t2076;
t1829 = t2228 * t1922 - t2232 * t1927;
t2100 = pkin(4) * t2160 - qJ(5) * t2162;
t1799 = -t2135 * pkin(4) - t2183 * qJ(5) + t2162 * t2100 + qJDD(5) + t1829;
t2246 = t2055 * pkin(5) + t1799 - t2526;
t2243 = 0.2e1 * qJD(6) * t2185 - t2246;
t2514 = t2243 - pkin(5) * (t2489 + t2403);
t2347 = t2233 * t2173;
t2348 = t2229 * t2174;
t2112 = t2347 - t2348;
t2512 = t2230 * t2112 - t2234 * t2257;
t2389 = t2227 * t2234;
t2390 = t2227 * t2230;
t2469 = -t2389 * t2495 + t2390 * t2475 - t2525;
t2506 = t2235 * t2469 + t2548;
t2505 = -t2231 * t2469 + t2545;
t2470 = -t2002 * t2389 + t2390 * t2477 - t2524;
t2504 = t2235 * t2470 + t2547;
t2503 = -t2231 * t2470 + t2544;
t2398 = t2191 * t2189;
t2248 = t2257 - t2398;
t2501 = t2229 * t2248;
t2500 = t2233 * t2248;
t1983 = t2229 * t2077 - t2233 * t2237;
t1921 = -t2257 * pkin(3) - t2378 * pkin(10) + t2191 * t2146 + t1983;
t2249 = t2054 * pkin(4) + t1921 - t2529;
t2493 = pkin(5) * t2467 + t2160 * t2531 - t2249;
t2488 = -(t2496 + t2054) * qJ(6) + pkin(5) * t2494;
t1756 = (pkin(5) * t2160 + t2531) * t2185 + t2246;
t2179 = 0.2e1 * t2438;
t2115 = pkin(5) * t2162 - qJ(6) * t2185;
t1830 = t2232 * t1922 + t2228 * t1927;
t2333 = t2183 * pkin(4) - t2135 * qJ(5) + t2160 * t2100 - t1830;
t2250 = -t2054 * pkin(5) - qJ(6) * t2467 + t2185 * t2115 + qJDD(6) - t2333;
t1769 = t2179 + t2250;
t2464 = pkin(4) + qJ(6);
t2487 = qJ(5) * t1769 - t1756 * t2464;
t1798 = t2179 - t2333;
t1722 = t1798 * t2232 + t1799 * t2228;
t2460 = pkin(4) * t2185;
t2465 = -0.2e1 * qJD(5);
t1813 = (t2465 + t2460) * t2162 + t2249;
t2364 = qJ(5) * t2228 + pkin(3);
t2458 = pkin(4) * t2232;
t2486 = -t1813 * (t2364 + t2458) + pkin(10) * t1722;
t2010 = qJ(5) * t2497;
t2485 = -t2464 * t2489 - t2010;
t2439 = qJD(5) * t2162;
t1789 = t2249 - 0.2e1 * t2439 + t2557;
t2484 = t2232 * t1789 + t2364 * t2496 - t2641;
t2153 = 0.2e1 * t2439;
t1788 = -pkin(4) * t2401 - qJ(5) * t2021 + t2153 - t2249;
t2483 = t2585 + t2228 * t1788 - t2021 * (pkin(3) + t2458);
t2033 = qJ(5) * t2491;
t2482 = -t2093 * t2464 + t2033;
t2481 = t2464 * t2492 + t2528;
t2472 = -t2002 * t2391 + t2392 * t2477 + t2521;
t2471 = -t2391 * t2495 + t2392 * t2475 + t2522;
t2187 = t2189 ^ 2;
t2188 = t2191 ^ 2;
t2223 = t2226 ^ 2;
t2461 = pkin(3) * t2229;
t2459 = pkin(4) * t2228;
t2457 = pkin(5) * t1756;
t2456 = pkin(5) * t1769;
t2455 = pkin(5) * t2491;
t2454 = pkin(5) * t2492;
t2453 = pkin(5) * t2093;
t2445 = t2054 * qJ(6);
t2441 = qJD(1) * t2234;
t2435 = t1921 * t2228;
t2434 = t1921 * t2232;
t2411 = t2076 * t2229;
t2410 = t2076 * t2233;
t2123 = -t2257 - t2398;
t2406 = t2123 * t2229;
t2405 = t2123 * t2233;
t2366 = t2230 * t2234 * t2468;
t2214 = t2223 * t2366;
t2194 = t2214 + t2355;
t2397 = t2194 * t2230;
t2396 = t2194 * t2234;
t2195 = -t2214 + t2355;
t2395 = t2195 * t2230;
t2394 = t2195 * t2234;
t2393 = t2223 * t2468;
t2176 = t2226 * t2259 + t2448;
t2387 = t2230 * t2176;
t2386 = t2234 * t2176;
t2224 = t2230 ^ 2;
t2225 = t2234 ^ 2;
t2381 = t2224 + t2225;
t2377 = -pkin(3) * t2233 - pkin(2);
t2370 = t2230 * t2398;
t2369 = t2234 * t2398;
t2368 = t2224 * t2393;
t2367 = t2225 * t2393;
t2363 = -pkin(4) * t2489 - t2010;
t2362 = t2115 - t2460;
t2360 = 0.2e1 * qJD(2) + t2443;
t1753 = t1829 * t2228 + t2232 * t1830;
t1877 = t1983 * t2229 + t2233 * t1984;
t2356 = -t2215 * t2231 - t2235 * t2216;
t2351 = -pkin(3) * t1921 + pkin(10) * t1753;
t2213 = qJDD(1) * t2235 - t2231 * t2468;
t2350 = -pkin(7) * t2213 - g(3) * t2231;
t2349 = t2229 * t2173;
t2346 = t2233 * t2174;
t2345 = -pkin(4) * t1799 + qJ(5) * t1798;
t2184 = -t2368 - t2354;
t2145 = -t2184 * t2230 - t2394;
t2339 = pkin(8) * t2145 - t2387;
t2201 = -t2354 - t2367;
t2151 = t2201 * t2234 - t2397;
t2338 = pkin(8) * t2151 + t2386;
t2335 = qJD(1) * t2361;
t2334 = t2223 * t2230 * t2365;
t1705 = t1756 * t2228 + t1769 * t2232;
t2238 = (t2465 - t2362) * t2162 - t2493;
t1770 = t2238 + t2445;
t1685 = t1705 * t2233 + t1770 * t2229;
t1704 = -t1756 * t2232 + t1769 * t2228;
t2332 = t1685 * t2230 - t1704 * t2234;
t1702 = t1722 * t2233 + t1813 * t2229;
t1721 = t1798 * t2228 - t1799 * t2232;
t2331 = t1702 * t2230 - t1721 * t2234;
t1727 = t1753 * t2233 + t1921 * t2229;
t1752 = -t1829 * t2232 + t1830 * t2228;
t2330 = t1727 * t2230 - t1752 * t2234;
t1865 = -t2229 * t2490 + t2576;
t2323 = t1865 * t2230 + t2609;
t2317 = t1877 * t2230 - t2076 * t2234;
t2310 = -t2002 * t2234 + t2230 * t2477;
t2309 = t2230 * t2475 - t2234 * t2495;
t1876 = -t1983 * t2233 + t1984 * t2229;
t2046 = (-t2160 * t2232 + t2162 * t2228) * t2185;
t2000 = t2046 * t2233 + t2404;
t2043 = t2117 + t2372;
t2306 = t2000 * t2230 + t2043 * t2234;
t2024 = -t2088 * t2229 - t2233 * t2513;
t2147 = t2188 - t2187;
t2305 = t2024 * t2230 - t2147 * t2234;
t2087 = -t2191 * t2219 + t2357;
t2089 = -t2173 + t2137;
t2025 = t2087 * t2233 + t2089 * t2229;
t2114 = t2187 + t2188;
t2304 = t2025 * t2230 + t2114 * t2234;
t2142 = -t2378 - t2187;
t2050 = t2142 * t2233 - t2501;
t2303 = t2050 * t2230 - t2234 * t2513;
t2148 = -t2188 - t2378;
t2062 = -t2148 * t2229 + t2405;
t2090 = (0.2e1 * qJD(3) - t2219) * t2189 + t2261;
t2302 = t2062 * t2230 + t2090 * t2234;
t2167 = -t2188 + t2378;
t2065 = -t2167 * t2229 + t2500;
t2301 = t2065 * t2230 - t2089 * t2234;
t2166 = t2187 - t2378;
t2066 = t2166 * t2233 + t2406;
t2085 = t2174 + t2284;
t2300 = t2066 * t2230 + t2085 * t2234;
t2143 = g(3) * t2391 + t2358;
t2144 = t2234 * t2193 + t2245;
t2299 = -t2234 * t2143 + t2230 * t2144;
t2060 = t2143 * t2230 + t2144 * t2234;
t2288 = t2226 * t2335;
t2205 = t2234 * t2288;
t2169 = t2205 + t2197;
t2204 = t2230 * t2288;
t2172 = -t2204 - t2265;
t2297 = t2169 * t2234 + t2172 * t2230;
t2170 = -t2205 + t2197;
t2171 = t2204 - t2265;
t2296 = -t2170 * t2234 + t2171 * t2230;
t2295 = t2184 * t2234 - t2395;
t2199 = t2354 - t2368;
t2294 = t2199 * t2234 + t2397;
t2293 = t2201 * t2230 + t2396;
t2200 = -t2354 + t2367;
t2292 = t2200 * t2230 + t2394;
t2291 = t2215 * t2235 - t2231 * t2216;
t2289 = t2226 * t2355;
t2080 = t2229 * t2284 - t2347;
t2286 = t2080 * t2230 + t2369;
t2082 = t2233 * t2137 + t2348;
t2285 = t2082 * t2230 - t2369;
t2282 = -t2434 + t2570;
t2281 = pkin(3) * t2021 + t2435 - t2585;
t1687 = -t1770 * t2464 + t2456;
t1707 = -qJ(5) * t1770 + t2457;
t1652 = -pkin(10) * t1704 - t1687 * t2228 + t1707 * t2232;
t1663 = -pkin(3) * t1704 - t2487;
t1684 = t1705 * t2229 - t1770 * t2233;
t1627 = -pkin(9) * t1684 + t1652 * t2233 - t1663 * t2229;
t2256 = -pkin(3) * t1770 + pkin(10) * t1705 + t1687 * t2232 + t1707 * t2228;
t1637 = -pkin(2) * t1684 - t2256;
t1664 = t1685 * t2234 + t1704 * t2230;
t2280 = pkin(8) * t1664 + t1627 * t2230 + t1637 * t2234;
t1686 = -pkin(3) * t1721 - t2345;
t1689 = -pkin(10) * t1721 + (-qJ(5) * t2232 + t2459) * t1813;
t1701 = t1722 * t2229 - t1813 * t2233;
t1651 = -pkin(9) * t1701 - t1686 * t2229 + t1689 * t2233;
t1662 = -pkin(2) * t1701 - t2486;
t1676 = t1702 * t2234 + t1721 * t2230;
t2279 = pkin(8) * t1676 + t1651 * t2230 + t1662 * t2234;
t2244 = -pkin(5) * t2497 + t2250;
t1724 = t2059 * t2464 + t2179 + t2244;
t1730 = -t2514 + t2556;
t1681 = -t1724 * t2228 + t1730 * t2232 - t2584;
t1783 = -t2485 - t2587;
t1669 = t1681 * t2233 - t1783 * t2229 - t2613;
t2255 = t1724 * t2232 + t1730 * t2228 + t2569;
t1677 = -t2255 - t2617;
t2278 = t1669 * t2230 + t1677 * t2234 + t2642;
t2240 = t2153 + t2493;
t1725 = t2162 * t2115 + t2240 + t2488 - t2557;
t1910 = -t2454 - t2527;
t1712 = -t1725 * t2228 + t1910 * t2232 - t2583;
t1718 = t1756 - t2481 - t2586;
t1672 = t1712 * t2233 - t1718 * t2229 - t2612;
t2254 = t1725 * t2232 + t1910 * t2228 + t2570;
t1699 = -t2254 - t2616;
t2277 = t1672 * t2230 + t1699 * t2234 + t2643;
t1747 = t2162 * t2362 + t2240 - t2445 + t2453 + t2529;
t1840 = t2464 * t2490 + t2455;
t1710 = t1747 * t2232 - t1840 * t2228 + t2611;
t1720 = -t2250 - t2482 + t2600;
t1862 = t2233 * t2490 + t2582;
t1673 = -pkin(9) * t1862 + t1710 * t2233 - t1720 * t2229;
t2253 = pkin(3) * t2490 + t1747 * t2228 + t1840 * t2232 + t2585;
t1696 = -pkin(2) * t1862 - t2253;
t1814 = t1865 * t2234 - t2610;
t2276 = pkin(8) * t1814 + t1673 * t2230 + t1696 * t2234;
t1726 = t1753 * t2229 - t1921 * t2233;
t1679 = -pkin(9) * t1726 + (-pkin(10) * t2233 + t2461) * t1752;
t1691 = -pkin(2) * t1726 - t2351;
t1698 = t1727 * t2234 + t1752 * t2230;
t2275 = pkin(8) * t1698 + t1679 * t2230 + t1691 * t2234;
t1784 = pkin(4) * t2059 + t1798;
t1787 = t1799 + t2556;
t1703 = -t1784 * t2228 + t1787 * t2232 - t2584;
t1812 = -t2363 - t2587;
t1680 = t1703 * t2233 - t1812 * t2229 - t2613;
t2252 = t1784 * t2232 + t1787 * t2228 + t2569;
t1692 = -t2252 - t2617;
t2274 = t1680 * t2230 + t1692 * t2234 + t2642;
t1739 = qJ(5) * t2515 - t1789 * t2228 + t2583;
t2242 = -pkin(4) * t2492 + t1799 - t2528;
t1749 = -t2242 + t2586;
t1688 = t1739 * t2233 - t1749 * t2229 + t2612;
t1714 = -t2484 + t2616;
t2273 = t1688 * t2230 + t1714 * t2234 - t2643;
t1740 = t1788 * t2232 + t2021 * t2459 + t2611;
t2262 = -pkin(4) * t2093 + t2033 - t2333;
t1750 = -t2262 + t2600;
t1690 = t1740 * t2233 - t1750 * t2229 + t2670;
t1715 = -t2483 + t2671;
t2272 = t1690 * t2230 + t1715 * t2234 + t2682;
t1723 = -t1752 - t2584;
t1700 = t1723 * t2233 + t2461 * t2535 - t2613;
t2263 = t1753 + t2569;
t1706 = -t2263 - t2617;
t2271 = t1700 * t2230 + t1706 * t2234 + t2642;
t1796 = t1829 - t2586;
t1831 = t2435 - t2583;
t1716 = -t1796 * t2229 + t1831 * t2233 - t2612;
t1754 = -t2282 - t2616;
t2270 = t1716 * t2230 + t1754 * t2234 + t2643;
t1797 = t1830 - t2615;
t1838 = t2434 - t2611;
t1717 = -t1797 * t2229 + t1838 * t2233 - t2670;
t1755 = -t2281 - t2671;
t2269 = t1717 * t2230 + t1755 * t2234 - t2682;
t2049 = t2142 * t2229 + t2500;
t1909 = -pkin(2) * t2049 + t1983;
t1962 = -pkin(9) * t2049 + t2411;
t1966 = t2050 * t2234 + t2230 * t2513;
t2268 = pkin(8) * t1966 + t1909 * t2234 + t1962 * t2230;
t2061 = t2148 * t2233 + t2406;
t1913 = -pkin(2) * t2061 + t1984;
t1967 = -pkin(9) * t2061 + t2410;
t1968 = t2062 * t2234 - t2090 * t2230;
t2267 = pkin(8) * t1968 + t1913 * t2234 + t1967 * t2230;
t2109 = t2170 * t2230 + t2171 * t2234;
t2266 = pkin(8) * t2109 + t2060;
t2023 = t2087 * t2229 - t2089 * t2233;
t1827 = -pkin(9) * t2023 - t1876;
t1942 = t2025 * t2234 - t2114 * t2230;
t2264 = pkin(8) * t1942 + t1827 * t2230 - t2023 * t2462;
t1839 = t1877 * t2234 + t2076 * t2230;
t2260 = pkin(8) * t1839 + t1876 * t2353;
t2239 = -pkin(5) * t2403 + t2243;
t2222 = t2226 * t2223;
t2212 = t2231 * qJDD(1) + t2235 * t2468;
t2209 = t2227 * t2355;
t2203 = t2381 * t2393;
t2202 = (t2224 - t2225) * t2393;
t2198 = -pkin(7) * t2212 + g(3) * t2235;
t2175 = t2361 * t2381 * t2444;
t2168 = (t2360 * t2441 + t2379) * t2226;
t2165 = t2234 * t2197 - t2224 * t2288;
t2164 = (-t2225 * t2335 + t2230 * t2283) * t2226;
t2150 = t2200 * t2234 - t2395;
t2149 = -t2199 * t2230 + t2396;
t2141 = (t2227 * t2197 + (qJD(2) * t2227 * t2226 + (t2226 * t2227 ^ 2 + t2222) * qJD(1)) * t2441) * t2230;
t2140 = t2197 * t2392 + t2334;
t2139 = -t2223 * t2234 * t2283 - t2334;
t2138 = -t2222 * t2366 + (-t2360 * t2442 + t2380) * t2226 * t2389;
t2111 = t2349 + t2346;
t2110 = -t2169 * t2230 + t2172 * t2234;
t2099 = t2226 * t2172 + t2227 * t2293;
t2098 = -t2226 * t2171 + t2227 * t2292;
t2097 = -t2226 * t2170 + t2227 * t2294;
t2096 = -t2227 * t2172 + t2226 * t2293;
t2095 = t2227 * t2171 + t2226 * t2292;
t2094 = t2227 * t2170 + t2226 * t2294;
t2084 = -t2226 * t2168 + t2227 * t2295;
t2083 = t2227 * t2168 + t2226 * t2295;
t2081 = t2229 * t2137 - t2346;
t2079 = -t2233 * t2284 - t2349;
t2078 = t2234 * t2112 + t2230 * t2257;
t2075 = -t2226 * t2202 + t2227 * t2297;
t2074 = t2226 * t2203 + t2227 * t2296;
t2073 = t2227 * t2202 + t2226 * t2297;
t2072 = -t2227 * t2203 + t2226 * t2296;
t2064 = t2166 * t2229 - t2405;
t2063 = t2167 * t2233 + t2501;
t2031 = t2082 * t2234 + t2370;
t2030 = t2080 * t2234 - t2370;
t2029 = t2226 * t2176 + t2227 * t2299;
t2028 = -t2227 * t2176 + t2226 * t2299;
t2027 = -t2226 * t2111 + t2227 * t2512;
t2026 = t2227 * t2111 + t2226 * t2512;
t2022 = t2088 * t2233 - t2229 * t2513;
t1997 = t2046 * t2229 - t2131;
t1986 = t2066 * t2234 - t2085 * t2230;
t1985 = t2065 * t2234 + t2089 * t2230;
t1969 = -t2387 + (-t2096 * t2226 - t2099 * t2227) * pkin(8);
t1965 = -t2386 + (-t2083 * t2226 - t2084 * t2227) * pkin(8);
t1964 = -pkin(1) * t2096 + t2226 * t2143 + t2227 * t2338;
t1963 = pkin(1) * t2099 - t2227 * t2143 + t2226 * t2338;
t1955 = t2024 * t2234 + t2147 * t2230;
t1954 = -pkin(1) * t2083 + t2226 * t2144 + t2227 * t2339;
t1953 = pkin(1) * t2084 - t2227 * t2144 + t2226 * t2339;
t1946 = -t2226 * t2081 + t2227 * t2285;
t1945 = -t2226 * t2079 + t2227 * t2286;
t1944 = t2227 * t2081 + t2226 * t2285;
t1943 = t2227 * t2079 + t2226 * t2286;
t1941 = pkin(1) * t2029 + t2060 * t2451;
t1940 = pkin(8) * t2060 * t2227 - pkin(1) * t2028;
t1924 = -pkin(1) * t2072 + t2227 * t2266;
t1923 = pkin(1) * t2074 + t2226 * t2266;
t1920 = pkin(2) * t2090 + pkin(9) * t2062 + t2411;
t1918 = -t2226 * t2064 + t2227 * t2300;
t1917 = -t2226 * t2063 + t2227 * t2301;
t1916 = t2227 * t2064 + t2226 * t2300;
t1915 = t2227 * t2063 + t2226 * t2301;
t1914 = (-t2028 * t2226 - t2029 * t2227) * pkin(8);
t1912 = -pkin(2) * t2513 + pkin(9) * t2050 - t2410;
t1911 = (-t2072 * t2226 - t2074 * t2227) * pkin(8) - t2299;
t1908 = -t2226 * t2061 + t2227 * t2302;
t1907 = t2227 * t2061 + t2226 * t2302;
t1906 = t2000 * t2234 - t2043 * t2230;
t1903 = -t2226 * t2049 + t2227 * t2303;
t1902 = t2227 * t2049 + t2226 * t2303;
t1869 = -t2226 * t2022 + t2227 * t2305;
t1868 = t2227 * t2022 + t2226 * t2305;
t1861 = -t2226 * t2023 + t2227 * t2304;
t1860 = t2227 * t2023 + t2226 * t2304;
t1841 = -pkin(2) * t2076 + pkin(9) * t1877;
t1837 = -t2226 * t1997 + t2227 * t2306;
t1834 = t2227 * t1997 + t2226 * t2306;
t1817 = pkin(2) * t2114 + pkin(9) * t2025 + t1877;
t1811 = t2227 * t2309 - t2525;
t1806 = t2227 * t2310 - t2524;
t1805 = t2226 * t2309 + t2522;
t1800 = t2226 * t2310 + t2521;
t1786 = -t2226 * t1876 + t2227 * t2317;
t1785 = t2227 * t1876 + t2226 * t2317;
t1760 = -t2226 * t1862 + t2227 * t2323;
t1757 = t2227 * t1862 + t2226 * t2323;
t1751 = -t2230 * t1913 + t2234 * t1967 + (-t1907 * t2226 - t1908 * t2227) * pkin(8);
t1748 = -t2230 * t1909 + t2234 * t1962 + (-t1902 * t2226 - t1903 * t2227) * pkin(8);
t1738 = -pkin(1) * t1907 - t2226 * t1920 + t2227 * t2267;
t1737 = pkin(1) * t1908 + t2227 * t1920 + t2226 * t2267;
t1729 = -pkin(1) * t1902 - t2226 * t1912 + t2227 * t2268;
t1728 = pkin(1) * t1903 + t2227 * t1912 + t2226 * t2268;
t1719 = t2023 * t2463 + t2234 * t1827 + (-t1860 * t2226 - t1861 * t2227) * pkin(8);
t1713 = t1797 * t2233 + t1838 * t2229 - t2661;
t1711 = t1796 * t2233 + t1831 * t2229 + t2598;
t1709 = -pkin(1) * t1860 - t2226 * t1817 + t2227 * t2264;
t1708 = pkin(1) * t1861 + t2227 * t1817 + t2226 * t2264;
t1697 = t2352 * t1876 + (-t1785 * t2226 - t1786 * t2227) * pkin(8);
t1695 = -pkin(1) * t1785 - t2226 * t1841 + t2227 * t2260;
t1694 = pkin(1) * t1786 + t2227 * t1841 + t2226 * t2260;
t1693 = t2229 * t1723 + t2377 * t2535 + t2614;
t1683 = t1740 * t2229 + t1750 * t2233 + t2661;
t1682 = t1739 * t2229 + t1749 * t2233 - t2598;
t1678 = t1703 * t2229 + t1812 * t2233 + t2599;
t1675 = -t2226 * t1726 + t2227 * t2330;
t1674 = t2227 * t1726 + t2226 * t2330;
t1671 = pkin(9) * t1865 + t1710 * t2229 + t1720 * t2233 + t2618;
t1670 = t1712 * t2229 + t1718 * t2233 + t2598;
t1668 = pkin(9) * t1727 + (-pkin(10) * t2229 + t2377) * t1752;
t1667 = t2234 * t1717 - t2230 * t1755 + t2695;
t1666 = t1681 * t2229 + t1783 * t2233 + t2599;
t1665 = t2234 * t1716 - t2230 * t1754 - t2645;
t1661 = -t2226 * t1701 + t2227 * t2331;
t1660 = t2227 * t1701 + t2226 * t2331;
t1659 = -t2226 * t1713 + t2227 * t2269 + t2694;
t1658 = t2227 * t1713 + t2226 * t2269 - t2693;
t1657 = -t2226 * t1711 + t2227 * t2270 - t2657;
t1656 = t2227 * t1711 + t2226 * t2270 + t2659;
t1655 = t2234 * t1690 - t2230 * t1715 - t2695;
t1654 = t2234 * t1688 - t2230 * t1714 + t2645;
t1653 = t2234 * t1700 - t2230 * t1706 + t2646;
t1650 = -t2226 * t1684 + t2227 * t2332;
t1649 = t2227 * t1684 + t2226 * t2332;
t1648 = t2234 * t1672 - t2230 * t1699 - t2645;
t1647 = t2234 * t1680 - t2230 * t1692 + t2646;
t1646 = t2234 * t1673 - t2230 * t1696 + (-t1757 * t2226 - t1760 * t2227) * pkin(8);
t1645 = -t2226 * t1693 + t2227 * t2271 - t2658;
t1644 = t2227 * t1693 + t2226 * t2271 + t2660;
t1643 = -t2226 * t1683 + t2227 * t2272 - t2694;
t1642 = t2227 * t1683 + t2226 * t2272 + t2693;
t1641 = -t2226 * t1682 + t2227 * t2273 + t2657;
t1640 = t2227 * t1682 + t2226 * t2273 - t2659;
t1639 = -pkin(2) * t1721 + pkin(9) * t1702 + t1686 * t2233 + t1689 * t2229;
t1638 = t2234 * t1669 - t2230 * t1677 + t2646;
t1636 = -t2226 * t1678 + t2227 * t2274 - t2658;
t1635 = t2227 * t1678 + t2226 * t2274 + t2660;
t1634 = -t2226 * t1670 + t2227 * t2277 - t2657;
t1633 = t2227 * t1670 + t2226 * t2277 + t2659;
t1632 = -pkin(1) * t1757 - t2226 * t1671 + t2227 * t2276;
t1631 = pkin(1) * t1760 + t2227 * t1671 + t2226 * t2276;
t1630 = t2234 * t1679 - t2230 * t1691 + (-t1674 * t2226 - t1675 * t2227) * pkin(8);
t1629 = -t2226 * t1666 + t2227 * t2278 - t2658;
t1628 = t2227 * t1666 + t2226 * t2278 + t2660;
t1626 = -pkin(2) * t1704 + pkin(9) * t1685 + t1652 * t2229 + t1663 * t2233;
t1625 = -pkin(1) * t1674 - t2226 * t1668 + t2227 * t2275;
t1624 = pkin(1) * t1675 + t2227 * t1668 + t2226 * t2275;
t1623 = t2234 * t1651 - t2230 * t1662 + (-t1660 * t2226 - t1661 * t2227) * pkin(8);
t1622 = -pkin(1) * t1660 - t2226 * t1639 + t2227 * t2279;
t1621 = pkin(1) * t1661 + t2227 * t1639 + t2226 * t2279;
t1620 = t2234 * t1627 - t2230 * t1637 + (-t1649 * t2226 - t1650 * t2227) * pkin(8);
t1619 = -pkin(1) * t1649 - t2226 * t1626 + t2227 * t2280;
t1618 = pkin(1) * t1650 + t2227 * t1626 + t2226 * t2280;
t1 = [0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, t2213, 0, -t2212, 0, t2350, -t2198, -t2291, -pkin(7) * t2291, -t2231 * t2141 + t2165 * t2235, -t2231 * t2075 + t2110 * t2235, -t2231 * t2097 + t2149 * t2235, -t2231 * t2138 + t2164 * t2235, -t2231 * t2098 + t2150 * t2235, t2235 * t2175 + t2231 * t2289, t2235 * t1969 - t2231 * t1964 - pkin(7) * (t2099 * t2235 + t2231 * t2151), t2235 * t1965 - t2231 * t1954 - pkin(7) * (t2084 * t2235 + t2231 * t2145), t2235 * t1911 - t2231 * t1924 - pkin(7) * (t2074 * t2235 + t2231 * t2109), t2235 * t1914 - t2231 * t1940 - pkin(7) * (t2029 * t2235 + t2231 * t2060), -t2231 * t1946 + t2031 * t2235, -t2231 * t1869 + t1955 * t2235, -t2231 * t1917 + t1985 * t2235, -t2231 * t1945 + t2030 * t2235, -t2231 * t1918 + t1986 * t2235, -t2231 * t2027 + t2078 * t2235, t2235 * t1748 - t2231 * t1729 - pkin(7) * (t1903 * t2235 + t2231 * t1966), t2235 * t1751 - t2231 * t1738 - pkin(7) * (t1908 * t2235 + t2231 * t1968), t2235 * t1719 - t2231 * t1709 - pkin(7) * (t1861 * t2235 + t2231 * t1942), t2235 * t1697 - t2231 * t1695 - pkin(7) * (t1786 * t2235 + t2231 * t1839), t2505, t2684, t2653, -t2231 * t1806 + t2544, -t2654, t2564, -t2231 * t1657 + t2235 * t1665 - t2647, -t2231 * t1659 + t2235 * t1667 + t2697, -t2231 * t1645 + t2235 * t1653 - t2648, t2235 * t1630 - t2231 * t1625 - pkin(7) * (t1675 * t2235 + t2231 * t1698), -t2231 * t1837 + t1906 * t2235, -t2653, t2654, -t2231 * t1811 + t2545, t2684, t2503, -t2231 * t1636 + t2235 * t1647 - t2648, -t2231 * t1641 + t2235 * t1654 + t2647, -t2231 * t1643 + t2235 * t1655 - t2697, t2235 * t1623 - t2231 * t1622 - pkin(7) * (t1661 * t2235 + t2231 * t1676), t2564, t2654, t2653, t2503, -t2684, t2505, -t2231 * t1629 + t2235 * t1638 - t2648, t2235 * t1646 - t2231 * t1632 - pkin(7) * (t1760 * t2235 + t2231 * t1814), -t2231 * t1634 + t2235 * t1648 - t2647, t2235 * t1620 - t2231 * t1619 - pkin(7) * (t1650 * t2235 + t2231 * t1664); 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, t2212, 0, t2213, 0, t2198, t2350, t2356, pkin(7) * t2356, t2141 * t2235 + t2231 * t2165, t2075 * t2235 + t2231 * t2110, t2097 * t2235 + t2231 * t2149, t2138 * t2235 + t2231 * t2164, t2098 * t2235 + t2231 * t2150, t2231 * t2175 - t2235 * t2289, t2231 * t1969 + t2235 * t1964 + pkin(7) * (-t2231 * t2099 + t2151 * t2235), t2231 * t1965 + t2235 * t1954 + pkin(7) * (-t2231 * t2084 + t2145 * t2235), t2231 * t1911 + t2235 * t1924 + pkin(7) * (-t2231 * t2074 + t2109 * t2235), t2231 * t1914 + t2235 * t1940 + pkin(7) * (-t2231 * t2029 + t2060 * t2235), t1946 * t2235 + t2231 * t2031, t1869 * t2235 + t2231 * t1955, t1917 * t2235 + t2231 * t1985, t1945 * t2235 + t2231 * t2030, t1918 * t2235 + t2231 * t1986, t2027 * t2235 + t2231 * t2078, t2231 * t1748 + t2235 * t1729 + pkin(7) * (-t2231 * t1903 + t1966 * t2235), t2231 * t1751 + t2235 * t1738 + pkin(7) * (-t2231 * t1908 + t1968 * t2235), t2231 * t1719 + t2235 * t1709 + pkin(7) * (-t2231 * t1861 + t1942 * t2235), t2231 * t1697 + t2235 * t1695 + pkin(7) * (-t2231 * t1786 + t1839 * t2235), t2506, -t2683, -t2690, t1806 * t2235 + t2547, t2689, t2559, t2235 * t1657 + t2231 * t1665 - t2649, t2235 * t1659 + t2231 * t1667 + t2696, t2235 * t1645 + t2231 * t1653 + t2650, t2231 * t1630 + t2235 * t1625 + pkin(7) * (-t2231 * t1675 + t1698 * t2235), t1837 * t2235 + t2231 * t1906, t2690, -t2689, t1811 * t2235 + t2548, -t2683, t2504, t2235 * t1636 + t2231 * t1647 + t2650, t2235 * t1641 + t2231 * t1654 + t2649, t2235 * t1643 + t2231 * t1655 - t2696, t2231 * t1623 + t2235 * t1622 + pkin(7) * (-t2231 * t1661 + t1676 * t2235), t2559, -t2689, -t2690, t2504, t2683, t2506, t2235 * t1629 + t2231 * t1638 + t2650, t2231 * t1646 + t2235 * t1632 + pkin(7) * (-t2231 * t1760 + t1814 * t2235), t2235 * t1634 + t2231 * t1648 - t2649, t2231 * t1620 + t2235 * t1619 + pkin(7) * (-t2231 * t1650 + t1664 * t2235); 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, qJDD(1), t2215, t2216, 0, 0, t2140, t2073, t2094, t2139, t2095, t2209, t1963, t1953, t1923, t1941, t1944, t1868, t1915, t1943, t1916, t2026, t1728, t1737, t1708, t1694, t2471, -t1741, -t1773, t1800, t1776, t2533, t1656, t1658, t1644, t1624, t1834, t1773, -t1776, t1805, -t1741, t2472, t1635, t1640, t1642, t1621, t2533, -t1776, -t1773, t2472, t1741, t2471, t1628, t1631, t1633, t1618; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0, -t2468, 0, 0, -g(3), -t2215, 0, t2165, t2110, t2149, t2164, t2150, t2175, t1969, t1965, t1911, t1914, t2031, t1955, t1985, t2030, t1986, t2078, t1748, t1751, t1719, t1697, t2509, -t1793, -t1823, t2510, t1826, t2508, t1665, t1667, t1653, t1630, t1906, t1823, -t1826, t2509, -t1793, t2510, t1647, t1654, t1655, t1623, t2508, -t1826, -t1823, t2510, t1793, t2509, t1638, t1646, t1648, t1620; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2468, 0, qJDD(1), 0, g(3), 0, -t2216, 0, t2141, t2075, t2097, t2138, t2098, -t2289, t1964, t1954, t1924, t1940, t1946, t1869, t1917, t1945, t1918, t2027, t1729, t1738, t1709, t1695, t2469, -t1744, -t1779, t1806, t1782, t2532, t1657, t1659, t1645, t1625, t1837, t1779, -t1782, t1811, -t1744, t2470, t1636, t1641, t1643, t1622, t2532, -t1782, -t1779, t2470, t1744, t2469, t1629, t1632, t1634, t1619; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t2215, t2216, 0, 0, t2140, t2073, t2094, t2139, t2095, t2209, t1963, t1953, t1923, t1941, t1944, t1868, t1915, t1943, t1916, t2026, t1728, t1737, t1708, t1694, t2471, -t1741, -t1773, t1800, t1776, t2533, t1656, t1658, t1644, t1624, t1834, t1773, -t1776, t1805, -t1741, t2472, t1635, t1640, t1642, t1621, t2533, -t1776, -t1773, t2472, t1741, t2471, t1628, t1631, t1633, t1618; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2197, t2172, t2194, -t2205, t2200, t2205, 0, -t2176, t2143, 0, t2082, t2024, t2065, t2080, t2066, t2112, t1962, t1967, t1827, -pkin(9) * t1876, t2475, -t1857, -t1886, t2477, t1889, t2476, t1716, t1717, t1700, t1679, t2000, t1886, -t1889, t2475, -t1857, t2477, t1680, t1688, t1690, t1651, t2476, -t1889, -t1886, t2477, t1857, t2475, t1669, t1673, t1672, t1627; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2204, t2169, t2199, -t2265, t2195, -t2204, t2176, 0, t2144, 0, -t2398, -t2147, -t2089, t2398, t2085, -t2257, t1909, t1913, -pkin(2) * t2023, -pkin(2) * t1876, -t2495, -t2536, t1972, -t2002, -t1975, t2498, t1754, t1755, t1706, t1691, t2043, -t1972, t1975, -t2495, -t2536, -t2002, t1692, t1714, t1715, t1662, t2498, t1975, t1972, -t2002, t2536, -t2495, t1677, t1696, t1699, t1637; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2214, t2202, t2170, t2214, t2171, t2355, -t2143, -t2144, 0, 0, t2081, t2022, t2063, t2079, t2064, t2111, t1912, t1920, t1817, t1841, t2478, -t1854, -t1880, t2479, t1883, t2480, t1711, t1713, t1693, t1668, t1997, t1880, -t1883, t2478, -t1854, t2479, t1678, t1682, t1683, t1639, t2480, -t1883, -t1880, t2479, t1854, t2478, t1666, t1671, t1670, t1626; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2137, -t2513, t2248, -t2173, t2166, t2173, 0, t2076, t1983, 0, t2383, -t1896, -t1978, t2287, t1981, t2340, t1831, t1838, t1723, -pkin(10) * t1752, t2046, t1978, -t1981, t2383, -t1896, t2287, t1703, t1739, t1740, t1689, t2340, -t1981, -t1978, t2287, t1896, t2383, t1681, t1710, t1712, t1652; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2174, t2088, t2167, -t2284, -t2123, t2174, -t2076, 0, t1984, 0, -t2402, t2101, -t2489, t2402, t2497, -t2135, t1796, t1797, -t2587, -pkin(3) * t1752, -t2135, t2489, -t2497, -t2402, t2101, t2402, t1812, t1749, t1750, t1686, -t2135, -t2497, -t2489, t2402, -t2101, -t2402, t1783, t1720, t1718, t1663; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2398, t2147, t2089, -t2398, -t2085, t2257, -t1983, -t1984, 0, 0, t2495, t2536, -t1972, t2002, t1975, -t2498, t2282, t2281, t2263, t2351, -t2043, t1972, -t1975, t2495, t2536, t2002, t2252, t2484, t2483, t2486, -t2498, -t1975, -t1972, t2002, -t2536, t2495, t2255, t2253, t2254, t2256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2055, -t2496, t2492, t2403, t2121, -t2403, 0, t1921, t1829, 0, -t2403, -t2492, -t2121, t2055, -t2496, t2403, t1787, t2527, t1788, -qJ(5) * t1813, -t2403, -t2121, t2492, t2403, t2496, t2055, t1730, t1747, t1910, t1707; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2401, t2490, t2119, -t2054, t2491, -t2401, -t1921, 0, t1830, 0, -t2401, -t2119, -t2491, t2401, t2490, -t2054, t1784, t1789, -pkin(4) * t2021, -pkin(4) * t1813, -t2401, -t2491, t2119, -t2054, -t2490, t2401, t1724, t1840, t1725, t1687; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2402, -t2101, t2489, -t2402, -t2497, t2135, -t1829, -t1830, 0, 0, t2135, -t2489, t2497, t2402, -t2101, -t2402, t2363, t2242, t2179 + t2262, t2345, t2135, t2497, t2489, -t2402, t2101, t2402, t2485, t1769 + t2482, t2239 + t2481, t2487; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2135, -t2489, t2497, t2402, -t2101, -t2402, 0, t1799, t1798, 0, t2135, t2497, t2489, -t2402, t2101, t2402, -qJ(6) * t2489, -qJ(6) * t2093 + t1769, t2239 + t2526, -qJ(6) * t1756; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2403, t2492, t2121, -t2055, t2496, -t2403, -t1799, 0, t1813, 0, t2403, t2121, -t2492, -t2403, -t2496, -t2055, t2514, t1770 - t2453, t2454, -t2457; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2401, t2119, t2491, -t2401, -t2490, t2054, -t1798, -t1813, 0, 0, t2401, t2491, -t2119, t2054, t2490, -t2401, -qJ(6) * t2059 + t2178 - t2244, -qJ(6) * t2490 - t2455, t2238 - t2488, qJ(6) * t1770 - t2456; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2135, t2497, t2489, -t2402, t2101, t2402, 0, t1769, -t1756, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2401, t2491, -t2119, t2054, t2490, -t2401, -t1769, 0, t1770, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2403, -t2121, t2492, t2403, t2496, t2055, t1756, -t1770, 0, 0;];
m_new_reg  = t1;