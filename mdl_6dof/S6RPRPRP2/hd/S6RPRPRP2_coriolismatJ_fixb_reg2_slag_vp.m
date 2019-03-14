% Calculate inertial parameters regressor of coriolis matrix for
% S6RPRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRPRP2_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP2_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP2_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP2_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:06:20
% EndTime: 2019-03-09 03:06:32
% DurationCPUTime: 8.67s
% Computational Cost: add. (9346->462), mult. (17307->562), div. (0->0), fcn. (18490->8), ass. (0->369)
t355 = cos(qJ(5));
t352 = sin(pkin(10));
t354 = sin(qJ(3));
t356 = cos(qJ(3));
t549 = cos(pkin(10));
t330 = t352 * t354 - t356 * t549;
t328 = t330 ^ 2;
t412 = t549 * t354;
t504 = t352 * t356;
t332 = t412 + t504;
t329 = t332 ^ 2;
t438 = -t329 - t328;
t589 = t438 * t355;
t595 = qJD(4) * t589;
t353 = sin(qJ(5));
t590 = t438 * t353;
t594 = qJD(4) * t590;
t593 = t589 * qJD(1);
t592 = t590 * qJD(1);
t345 = sin(pkin(9)) * pkin(1) + pkin(7);
t494 = qJ(4) + t345;
t325 = t494 * t354;
t327 = t494 * t356;
t190 = t549 * t325 + t352 * t327;
t587 = t190 * t353;
t591 = t587 / 0.2e1;
t323 = t332 * qJ(6);
t560 = t354 * pkin(3);
t229 = t332 * pkin(4) + t330 * pkin(8) + t560;
t188 = t353 * t229;
t586 = t190 * t355;
t491 = -t586 / 0.2e1 + t188 / 0.2e1;
t588 = -t323 - t491;
t584 = -t352 * t325 + t549 * t327;
t525 = t584 * t353;
t524 = t584 * t355;
t446 = t332 * qJD(1);
t423 = t355 * t446;
t251 = t330 * t423;
t222 = t353 * t330;
t462 = t222 * qJD(3);
t585 = t251 + t462;
t322 = t332 * qJD(3);
t304 = t355 * t322;
t460 = t222 * qJD(5);
t484 = -t304 + t460;
t437 = t329 - t328;
t226 = t355 * t330;
t211 = t226 * qJD(3);
t479 = qJD(5) * t353;
t433 = t332 * t479;
t583 = t433 + t211;
t478 = qJD(5) * t355;
t305 = t332 * t478;
t489 = t462 - t305;
t463 = t222 * qJD(1);
t169 = t463 + t479;
t160 = t437 * t353;
t472 = t160 * qJD(1);
t582 = t472 - t304;
t581 = t472 - t460;
t350 = t353 ^ 2;
t351 = t355 ^ 2;
t513 = t330 * t332;
t114 = (-t350 - t351 + 0.1e1) * t513;
t475 = t114 * qJD(2);
t347 = -cos(pkin(9)) * pkin(1) - pkin(2);
t334 = -t356 * pkin(3) + t347;
t360 = t330 * pkin(4) - t332 * pkin(8) + t334;
t120 = t353 * t360 + t524;
t536 = t120 * t355;
t119 = -t355 * t360 + t525;
t541 = t119 * t353;
t377 = -t541 / 0.2e1 - t536 / 0.2e1;
t126 = -t586 + t188;
t532 = t126 * t355;
t500 = t355 * t229;
t125 = t500 + t587;
t535 = t125 * t353;
t572 = t190 / 0.2e1;
t573 = t584 / 0.2e1;
t7 = (t532 / 0.2e1 - t535 / 0.2e1 + t572) * t332 + (t573 + t377) * t330;
t580 = -t7 * qJD(1) - t475;
t501 = t355 * qJ(6);
t563 = pkin(5) * t353;
t335 = t501 - t563;
t130 = t330 * t335 + t584;
t224 = t353 * t332;
t227 = t355 * t332;
t402 = pkin(5) * t224 - qJ(6) * t227;
t131 = t190 + t402;
t562 = t330 * pkin(5);
t80 = t119 - t562;
t552 = t80 * t353;
t514 = t330 * qJ(6);
t79 = t120 + t514;
t554 = t79 * t355;
t379 = -t554 / 0.2e1 - t552 / 0.2e1;
t561 = t332 * pkin(5);
t93 = -t125 - t561;
t550 = t93 * t353;
t92 = t323 + t126;
t551 = t92 * t355;
t5 = (t551 / 0.2e1 + t131 / 0.2e1 + t550 / 0.2e1) * t332 + (t130 / 0.2e1 + t379) * t330;
t579 = t5 * qJD(1) + t475;
t521 = t190 * t332;
t31 = t521 + (-t536 - t541) * t330;
t566 = -t351 / 0.2e1;
t567 = -t350 / 0.2e1;
t101 = (t566 + t567 + 0.1e1 / 0.2e1) * t513;
t476 = t101 * qJD(2);
t578 = qJD(1) * t31 + t476;
t21 = t131 * t332 + (-t552 - t554) * t330;
t577 = qJD(1) * t21 + t476;
t449 = t330 * qJD(1);
t576 = -t449 - qJD(5);
t310 = t350 * t330;
t311 = t351 * t330;
t575 = 0.2e1 * t353 * t227 * (-qJD(5) + t449) - (-t310 + t311) * qJD(3);
t574 = t120 / 0.2e1;
t503 = t353 * qJ(6);
t559 = t355 * pkin(5);
t398 = t503 + t559;
t214 = t398 * t332;
t571 = t214 / 0.2e1;
t570 = -t229 / 0.2e1;
t569 = -t335 / 0.2e1;
t568 = t335 / 0.2e1;
t565 = -t353 / 0.2e1;
t564 = t355 / 0.2e1;
t558 = qJD(3) * pkin(3);
t555 = t79 * t330;
t436 = -t119 / 0.2e1 + t80 / 0.2e1;
t361 = t436 * t355 + (-t79 / 0.2e1 + t574) * t353;
t8 = t330 * t571 + t332 * t361;
t553 = t8 * qJD(1);
t538 = t120 * t330;
t57 = t190 * t227 - t538;
t546 = qJD(1) * t57;
t73 = -t330 * t584 + t521;
t545 = qJD(1) * t73;
t10 = t80 * t226 - t93 * t227 + (t92 * t332 - t555) * t353;
t544 = t10 * qJD(1);
t11 = ((-t120 + t79) * t355 + (-t119 + t80) * t353) * t332;
t543 = t11 * qJD(1);
t542 = t119 * t330;
t540 = t119 * t355;
t381 = t562 / 0.2e1 - t436;
t403 = t514 / 0.2e1 + t79 / 0.2e1;
t537 = t120 * t353;
t12 = -t537 / 0.2e1 + t403 * t353 + t381 * t355;
t539 = t12 * qJD(1);
t534 = t125 * t355;
t533 = t126 * t353;
t13 = (t574 - t403) * t355 + t381 * t353;
t531 = t13 * qJD(1);
t530 = t130 * t353;
t529 = t130 * t355;
t528 = t131 * t353;
t527 = t131 * t355;
t17 = (t79 - t529) * t332 + (t92 + t527) * t330;
t526 = t17 * qJD(1);
t19 = (t533 + t534) * t332 + (-t537 + t540) * t330;
t523 = t19 * qJD(1);
t24 = (-t119 + t525) * t332 + (t125 - t587) * t330;
t518 = t24 * qJD(1);
t344 = t352 * pkin(3) + pkin(8);
t275 = t344 * t310;
t276 = t344 * t311;
t414 = -t275 / 0.2e1 - t276 / 0.2e1;
t346 = -pkin(3) * t549 - pkin(4);
t324 = -t398 + t346;
t511 = t332 * t324;
t421 = t511 / 0.2e1;
t375 = t421 + t414;
t378 = t564 * t93 + t565 * t92;
t29 = t375 + t378;
t517 = t29 * qJD(1);
t32 = -t538 + (t214 * t353 + t527) * t332;
t516 = t32 * qJD(1);
t33 = -t542 + (-t214 * t355 + t528) * t332;
t515 = t33 * qJD(1);
t512 = t330 * t352;
t510 = t332 * t344;
t509 = t332 * t346;
t508 = t335 * t353;
t507 = t344 * t330;
t420 = t509 / 0.2e1;
t374 = t420 + t414;
t376 = -t534 / 0.2e1 - t533 / 0.2e1;
t35 = t374 + t376;
t506 = t35 * qJD(1);
t415 = -t190 / 0.2e1 + t572;
t38 = t415 * t332 + (-t584 / 0.2e1 + t573) * t330;
t499 = t38 * qJD(1);
t41 = -t131 * t227 + t555;
t498 = t41 * qJD(1);
t56 = -t190 * t224 + t542;
t497 = t56 * qJD(1);
t416 = t227 / 0.2e1;
t493 = t528 / 0.2e1 + t324 * t416;
t492 = t591 + t500 / 0.2e1;
t210 = t226 * qJD(5);
t303 = t353 * t322;
t488 = t210 + t303;
t487 = -t251 - t305;
t486 = -t275 - t276;
t485 = (t567 + t351 / 0.2e1) * t332;
t482 = qJD(1) * t356;
t481 = qJD(3) * t353;
t480 = qJD(3) * t355;
t477 = t101 * qJD(1);
t473 = t119 * qJD(5);
t162 = t437 * t355;
t470 = t162 * qJD(1);
t413 = t549 * t332;
t366 = -t512 / 0.2e1 - t413 / 0.2e1;
t173 = (-t354 / 0.2e1 + t366) * pkin(3);
t468 = t173 * qJD(1);
t174 = t330 * t560 + t334 * t332;
t467 = t174 * qJD(1);
t175 = -t334 * t330 + t332 * t560;
t466 = t175 * qJD(1);
t465 = t437 * qJD(1);
t220 = (t350 / 0.2e1 + t566) * t332;
t464 = t220 * qJD(5);
t461 = t222 * qJD(4);
t459 = t222 * qJD(6);
t458 = t224 * qJD(1);
t212 = t226 * qJD(1);
t457 = t227 * qJD(1);
t231 = t310 + t311;
t455 = t231 * qJD(1);
t454 = t231 * qJD(3);
t451 = t438 * qJD(1);
t326 = t412 / 0.2e1 + t504 / 0.2e1;
t450 = t326 * qJD(1);
t320 = t330 * qJD(3);
t448 = t330 * qJD(4);
t447 = t330 * qJD(5);
t318 = t330 * qJD(6);
t445 = t332 * qJD(4);
t338 = t351 - t350;
t444 = t338 * qJD(5);
t339 = -t354 ^ 2 + t356 ^ 2;
t443 = t339 * qJD(1);
t442 = t353 * qJD(6);
t441 = t354 * qJD(3);
t440 = t355 * qJD(6);
t439 = t356 * qJD(3);
t313 = t561 / 0.2e1;
t417 = t507 / 0.2e1;
t435 = -t528 / 0.2e1 - t324 * t227 / 0.2e1 + t355 * t417;
t434 = t313 + t492;
t432 = t344 * t479;
t431 = t344 * t478;
t430 = t330 * t446;
t261 = t330 * t322;
t429 = t347 * t354 * qJD(1);
t428 = t347 * t482;
t341 = t353 * t478;
t427 = t353 * t446;
t426 = t332 * t442;
t340 = t353 * t480;
t425 = t353 * t440;
t424 = t354 * t439;
t422 = t332 * t440;
t419 = -t508 / 0.2e1;
t418 = -t507 / 0.2e1;
t156 = t326 + t485;
t411 = t156 * qJD(1) + t340;
t166 = t220 * qJD(1) - t340;
t279 = t355 * t329 * t353 * qJD(1);
t139 = qJD(3) * t220 + t279;
t408 = t329 * t341;
t407 = t353 * t423;
t406 = t353 * t304;
t405 = -t561 / 0.2e1 - t587 / 0.2e1;
t401 = 0.2e1 * t406;
t400 = -t214 / 0.2e1 + t418;
t6 = t130 * t131 + t79 * t92 + t80 * t93;
t397 = t6 * qJD(1) + t5 * qJD(2);
t9 = -t119 * t79 + t120 * t80 + t131 * t214;
t396 = t9 * qJD(1) + t8 * qJD(2);
t395 = t550 + t551;
t16 = -t119 * t125 + t120 * t126 + t190 * t584;
t394 = t16 * qJD(1) + t7 * qJD(2);
t18 = (-t80 + t530) * t332 + (-t93 - t528) * t330;
t393 = t18 * qJD(1);
t25 = (-t120 + t524) * t332 + (-t126 - t586) * t330;
t392 = t25 * qJD(1);
t53 = t334 * t560;
t391 = t53 * qJD(1) + t38 * qJD(2);
t390 = t532 - t535;
t389 = t324 * t330 + t510;
t388 = -t330 * t346 - t510;
t22 = (t419 - pkin(5) / 0.2e1) * t332 + (t570 + t400) * t355 + t405 + t493;
t239 = -t324 * t353 - t335 * t355;
t386 = -t22 * qJD(1) + t239 * qJD(3);
t238 = t324 * t355 - t508;
t358 = (t332 * t568 - t131 / 0.2e1) * t355 + (t421 + t400) * t353;
t27 = t358 + t588;
t385 = -t27 * qJD(1) + t238 * qJD(3);
t384 = t576 * t353;
t383 = -t335 * qJD(5) - t442;
t382 = -t92 * qJ(6) / 0.2e1 + t93 * pkin(5) / 0.2e1;
t380 = -t501 / 0.2e1 + t563 / 0.2e1;
t37 = t434 + t435;
t373 = -t37 * qJD(1) + t324 * t481;
t271 = t355 * t418;
t51 = t271 + (t420 + t570) * t355 + t415 * t353;
t372 = -t51 * qJD(1) - t346 * t481;
t359 = (t417 - t509 / 0.2e1) * t353 + t586 / 0.2e1;
t49 = t359 + t491;
t371 = -t49 * qJD(1) - t346 * t480;
t370 = t332 * t384;
t164 = t326 * qJD(5) + t430;
t312 = t351 * t329;
t232 = t350 * t329 - t312;
t142 = -t232 * qJD(1) + t401;
t252 = -t338 * qJD(3) + 0.2e1 * t407;
t368 = t160 * qJD(3) + t305 * t330;
t367 = t380 * t330;
t365 = t232 * qJD(5) + t330 * t401;
t137 = (t568 + t380) * t330;
t357 = t131 * t569 + t324 * t571 + t344 * t361;
t4 = t357 + t382;
t364 = t324 * t335 * qJD(3) - t4 * qJD(1) + t137 * qJD(2);
t363 = -qJD(5) * t398 + t440;
t240 = t312 + t328;
t362 = t240 * qJD(1) + t406 + t447;
t342 = t354 * t482;
t307 = t326 * qJD(3);
t280 = t353 * t422;
t278 = t576 * qJ(6);
t264 = t350 * qJD(3) + t407;
t217 = t231 * qJD(4);
t172 = t560 / 0.2e1 + t366 * pkin(3);
t171 = t212 + t478;
t157 = -t326 + t485;
t154 = -t261 * t351 - t408;
t153 = -t261 * t350 + t408;
t144 = t355 * t370;
t138 = t330 * t569 + t367;
t134 = t351 * t430 - t464;
t133 = t330 * t427 - t211;
t132 = t350 * t430 + t464;
t124 = t162 * qJD(3) - t330 * t433;
t117 = t120 * qJD(5);
t116 = t210 - t470;
t102 = t114 * qJD(3);
t100 = -t464 + (-t351 * t446 - t340) * t330;
t99 = t464 + (-t350 * t446 + t340) * t330;
t94 = t101 * qJD(4);
t68 = t303 + t470;
t52 = t346 * t416 + t271 + t492 + t591;
t50 = t359 - t491;
t36 = -t500 / 0.2e1 + t405 + t435;
t34 = t374 - t376;
t30 = qJD(3) * t38;
t28 = t375 - t378;
t26 = t358 - t588;
t23 = t332 * t419 + t355 * t400 + t313 + t434 + t493;
t15 = -t540 / 0.2e1 + t79 * t565 + t537 / 0.2e1 + t80 * t564 + (t503 / 0.2e1 + t559 / 0.2e1) * t330;
t14 = t367 + t377 - t379;
t3 = t357 - t382;
t2 = qJD(3) * t7 + t94;
t1 = qJD(3) * t5 + qJD(5) * t8 + t94;
t20 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t424, t339 * qJD(3), 0, -t424, 0, 0, t347 * t441, t347 * t439, 0, 0, -t261, -t437 * qJD(3), 0, t261, 0, 0, t174 * qJD(3), t175 * qJD(3), -qJD(4) * t438, qJD(3) * t53 + qJD(4) * t73, t154, t365, t124, t153, -t368, t261, qJD(3) * t24 + qJD(5) * t57 - t594, qJD(3) * t25 + qJD(5) * t56 - t595, -qJD(3) * t19, qJD(3) * t16 + qJD(4) * t31, t154, t124, -t365, t261, t368, t153, t18 * qJD(3) + t32 * qJD(5) - t329 * t425 - t594, -t10 * qJD(3) - t11 * qJD(5) - t330 * t426, qJD(3) * t17 + qJD(5) * t33 + qJD(6) * t240 + t595, qJD(3) * t6 + qJD(4) * t21 + qJD(5) * t9 + qJD(6) * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t342, t443, t439, -t342, -t441, 0, -t345 * t439 + t429, t345 * t441 + t428, 0, 0, -t430, -t465, -t320, t430, -t322, 0, -qJD(3) * t584 + t467, qJD(3) * t190 + t466 (t330 * t549 - t332 * t352) * t558 (-t190 * t352 - t549 * t584) * t558 + t172 * qJD(4) + t391, t100, t575, t68, t99, -t582, t164, t518 + (t353 * t388 - t524) * qJD(3) + t52 * qJD(5) (t355 * t388 + t525) * qJD(3) + t50 * qJD(5) + t392, qJD(3) * t390 - t523 (t344 * t390 + t346 * t584) * qJD(3) + t34 * qJD(4) + t394, t100, t68, -t575, t164, t582, t99 (-t353 * t389 - t529) * qJD(3) + t23 * qJD(5) + t157 * qJD(6) + t393, qJD(3) * t395 + t15 * qJD(5) - t544, t526 + (t355 * t389 - t530) * qJD(3) + t26 * qJD(5) + t280 (t130 * t324 + t344 * t395) * qJD(3) + t28 * qJD(4) + t3 * qJD(5) + t36 * qJD(6) + t397; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t451, qJD(3) * t172 + t545, 0, 0, 0, 0, 0, 0, -t592, -t593, 0, qJD(3) * t34 + t578, 0, 0, 0, 0, 0, 0, -t592, 0, t593, qJD(3) * t28 + qJD(5) * t14 + t577; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, -t142, t370, t139, t487, t307, qJD(3) * t52 - t117 + t546, qJD(3) * t50 + t473 + t497, 0, 0, -t139, t370, t142, t307, -t487, t139, qJD(3) * t23 - t117 + t516, t15 * qJD(3) + qJD(5) * t402 - t426 - t543, t26 * qJD(3) + t318 - t473 + t515, t3 * qJD(3) + t14 * qJD(4) + (-pkin(5) * t120 - qJ(6) * t119) * qJD(5) + t79 * qJD(6) + t396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t157 - t279, t370, t362, qJD(3) * t36 + qJD(5) * t79 + t498; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t441, -t439, 0, 0, 0, 0, 0, 0, 0, 0, -t322, t320, 0, t499 + (-t413 - t512) * t558, 0, 0, 0, 0, 0, 0, t484, t488, -t454 (t486 + t509) * qJD(3) - t580, 0, 0, 0, 0, 0, 0, t484, -t454, -t488 (t486 + t511) * qJD(3) + t138 * qJD(5) - t459 + t579; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t477, 0, 0, 0, 0, 0, 0, 0, 0, 0, t477; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t489, t583, 0, 0, 0, 0, 0, 0, 0, 0, t489, 0, -t583, t138 * qJD(3) - t214 * qJD(5) + t422 + t553; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t489; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t342, -t443, 0, t342, 0, 0, -t429, -t428, 0, 0, t430, t465, 0, -t430, 0, 0, -t445 - t467, t448 - t466, 0, qJD(4) * t173 - t391, t134, 0.2e1 * t144, t116, t132, t581, -t164, t51 * qJD(5) - t355 * t445 - t518, qJD(4) * t224 + qJD(5) * t49 - t392, -t217 + t523, qJD(4) * t35 - t394, t134, t116, -0.2e1 * t144, -t164, -t581, t132, -qJD(4) * t227 + qJD(5) * t22 + qJD(6) * t156 - t393, -qJD(5) * t12 + qJD(6) * t226 - t217 + t544, t27 * qJD(5) - t353 * t445 + t280 - t526, qJD(4) * t29 + qJD(5) * t4 + qJD(6) * t37 - t397; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t499, 0, 0, 0, 0, 0, 0, 0, 0, 0, t580, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5) * t137 - t579; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t341, t444, 0, -t341, 0, 0, t346 * t479, t346 * t478, 0, 0, t341, 0, -t444, 0, 0, -t341, -t239 * qJD(5) + t425, 0, -t238 * qJD(5) + t350 * qJD(6), t383 * t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t446, t449, 0, t468, 0, 0, 0, 0, 0, 0, -t423, t458, -t455, t506, 0, 0, 0, 0, 0, 0, -t457, -t455, -t427, t517; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, -t252, t171, t166, -t169, -t450, -t372 - t431, -t371 + t432, 0, 0, -t166, t171, t252, -t450, t169, t166, -t386 - t431, t363 - t539, -t385 - t432, t344 * t363 - t364; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t411, t171, t264, -t373 + t431; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t322, -t320, t451, -qJD(3) * t173 - t545, 0, 0, 0, 0, 0, 0, t592 - t484, -t224 * qJD(3) - t355 * t447 + t593, t454, -qJD(3) * t35 - t578, 0, 0, 0, 0, 0, 0, t227 * qJD(3) - t353 * t447 + t592, t454, -t593 + t488, -qJD(3) * t29 - qJD(5) * t13 + t459 - t577; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t477, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t477; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t446, -t449, 0, -t468, 0, 0, 0, 0, 0, 0, t423, -t458, t455, -t506, 0, 0, 0, 0, 0, 0, t457, t455, t427, -t517; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t169, t576 * t355, 0, 0, 0, 0, 0, 0, 0, 0, t384, 0, t171, -t383 - t531; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, t142, t133, -t139, t585, t307, -qJD(3) * t51 + t461 - t546, -t49 * qJD(3) + t355 * t448 - t497, 0, 0, t139, t133, -t142, t307, -t585, -t139, -t22 * qJD(3) + t353 * t448 - t516, qJD(3) * t12 + t543, -t27 * qJD(3) - t226 * qJD(4) + t318 - t515, qJ(6) * t318 - t4 * qJD(3) + t13 * qJD(4) - t396; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t137 - t553; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, t252, -t212, -t166, t463, t450, t372, t371, 0, 0, t166, -t212, -t252, t450, -t463, -t166, t386, t539, t385, t364; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t463, t355 * t449, 0, 0, 0, 0, 0, 0, 0, 0, t353 * t449, 0, -t212, t531; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(6), qJ(6) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t576, -t278; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t156 + t279, t133, -t362, -qJ(6) * t447 - t37 * qJD(3) - t461 - t498; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t411, -t212, -t264, t373; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t463; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t576, t278; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t20;