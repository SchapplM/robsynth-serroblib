% Calculate inertial parameters regressor of coriolis matrix for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRRR1_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:37:51
% EndTime: 2019-03-08 18:38:06
% DurationCPUTime: 10.40s
% Computational Cost: add. (7450->497), mult. (16249->688), div. (0->0), fcn. (18935->8), ass. (0->418)
t484 = qJD(2) + qJD(3);
t632 = qJD(4) + t484;
t357 = cos(qJ(5));
t355 = sin(qJ(3));
t356 = sin(qJ(2));
t358 = cos(qJ(2));
t582 = cos(qJ(3));
t321 = t355 * t358 + t356 * t582;
t581 = cos(qJ(4));
t316 = t581 * t321;
t320 = t355 * t356 - t358 * t582;
t354 = sin(qJ(4));
t546 = t354 * t320;
t612 = -t316 + t546;
t250 = t612 ^ 2;
t435 = t581 * t320 + t354 * t321;
t626 = -t435 ^ 2 + t250;
t630 = t626 * t357;
t634 = t630 * qJD(1);
t353 = sin(qJ(5));
t631 = t626 * t353;
t633 = t631 * qJD(1);
t629 = t626 * qJD(1);
t351 = t353 ^ 2;
t522 = qJD(1) * t612;
t621 = t435 * t522;
t628 = t351 * t621;
t352 = t357 ^ 2;
t627 = t352 * t621;
t575 = t612 * pkin(4);
t576 = pkin(6) * t435;
t128 = t575 - t576;
t624 = -t435 / 0.2e1;
t602 = t612 / 0.2e1;
t344 = t358 * pkin(2) + pkin(1);
t285 = -pkin(3) * t320 + t344;
t416 = -pkin(4) * t435 - pkin(6) * t612;
t88 = t285 + t416;
t568 = t612 * t88;
t613 = t484 * t357;
t623 = t353 * t613;
t622 = t357 * t612;
t524 = t351 + t352;
t438 = t524 * t88;
t420 = t435 * t438;
t432 = t484 * t435;
t491 = qJD(4) * t435;
t620 = t432 + t491;
t589 = -t352 / 0.2e1;
t590 = t351 / 0.2e1;
t444 = t590 + t589;
t109 = t444 * t612;
t483 = -qJD(3) - qJD(4);
t437 = qJD(2) - t483;
t547 = t353 * t357;
t614 = t437 * t547;
t79 = -t109 * qJD(1) + t614;
t421 = -t316 / 0.2e1;
t611 = t421 + t546 / 0.2e1;
t619 = qJD(5) * t611;
t618 = t611 * qJD(1);
t617 = -0.2e1 * t612;
t257 = t484 * t321;
t616 = t320 * t257;
t337 = t352 - t351;
t615 = t337 * t437;
t433 = t484 * t612;
t422 = t582 * t581;
t545 = t354 * t355;
t319 = (t422 - t545) * pkin(2);
t588 = t352 / 0.2e1;
t610 = (t588 + t590) * t319;
t504 = t109 * qJD(5);
t583 = t357 / 0.2e1;
t448 = t435 * t583;
t111 = t353 * t435;
t452 = t111 / 0.2e1;
t97 = t353 * t448 + t357 * t452;
t609 = -t97 * qJD(4) + t504;
t601 = t435 / 0.2e1;
t98 = (t601 + t624) * t547;
t608 = -qJD(4) * t98 + t504;
t349 = qJD(5) * t357;
t339 = t353 * t349;
t528 = t98 * qJD(1);
t85 = t528 - t339;
t84 = t528 + t339;
t607 = t484 * t97 - t504;
t606 = t484 * t98 + t504;
t467 = qJD(1) * t547;
t24 = t109 * t437 + t250 * t467;
t603 = pkin(4) / 0.2e1;
t598 = -t612 / 0.2e1;
t423 = pkin(2) * t582 + pkin(3);
t327 = t581 * t423;
t310 = pkin(2) * t545 - t327;
t302 = -pkin(4) + t310;
t597 = -t302 / 0.2e1;
t398 = t354 * t423;
t470 = t581 * t355;
t311 = pkin(2) * t470 + t398;
t303 = pkin(6) + t311;
t596 = -t303 / 0.2e1;
t595 = t311 / 0.2e1;
t307 = t316 / 0.2e1;
t594 = -t319 / 0.2e1;
t269 = t354 * pkin(3);
t342 = pkin(6) + t269;
t593 = t342 / 0.2e1;
t480 = t581 * pkin(3);
t343 = -t480 - pkin(4);
t592 = -t343 / 0.2e1;
t591 = -t351 / 0.2e1;
t587 = -t353 / 0.2e1;
t586 = t353 / 0.2e1;
t585 = t354 / 0.2e1;
t584 = -t357 / 0.2e1;
t580 = pkin(2) * t355;
t579 = pkin(2) * t356;
t578 = pkin(3) * t321;
t577 = pkin(4) * t353;
t474 = t582 * t354;
t317 = (t474 + t470) * pkin(2);
t574 = t317 * pkin(4);
t573 = pkin(2) * qJD(2);
t572 = pkin(3) * qJD(3);
t571 = pkin(3) * qJD(4);
t436 = t524 * t612;
t92 = -t578 + t128;
t91 = t92 - t579;
t1 = t436 * t91 + t420;
t570 = t1 * qJD(1);
t2 = t436 * t92 + t420;
t569 = t2 * qJD(1);
t3 = t128 * t436 + t420;
t567 = t3 * qJD(1);
t7 = t91 * t438;
t566 = t7 * qJD(1);
t8 = t92 * t438;
t565 = t8 * qJD(1);
t9 = t128 * t438;
t564 = t9 * qJD(1);
t390 = t317 * t602 - t435 * t594;
t361 = (t302 / 0.2e1 + t592) * t435 + (t596 + t593) * t612 + t390;
t12 = t361 * t353;
t561 = t12 * qJD(1);
t445 = t310 / 0.2e1 + t597;
t368 = (t596 + t595) * t612 - t445 * t435;
t397 = pkin(6) * t602 + t435 * t603;
t362 = t368 + t397;
t18 = t362 * t353;
t560 = t18 * qJD(1);
t550 = t311 * t612;
t551 = t310 * t435;
t365 = t551 / 0.2e1 - t550 / 0.2e1 + t390;
t473 = t581 * t435;
t557 = t612 * t354;
t374 = (-t557 / 0.2e1 - t473 / 0.2e1) * pkin(3);
t22 = t374 - t365;
t559 = t22 * qJD(1);
t558 = t612 * t303;
t556 = t435 * t302;
t410 = -t435 * t91 + t568;
t27 = t410 * t353;
t555 = t27 * qJD(1);
t28 = t410 * t357;
t554 = t28 * qJD(1);
t409 = -t435 * t92 + t568;
t29 = t409 * t353;
t553 = t29 * qJD(1);
t30 = t409 * t357;
t552 = t30 * qJD(1);
t408 = -t128 * t435 + t568;
t35 = t408 * t353;
t549 = t35 * qJD(1);
t548 = t353 * t612;
t118 = t357 * t435;
t36 = t408 * t357;
t542 = t36 * qJD(1);
t198 = t524 * t310;
t66 = -t198 * t303 + t302 * t311;
t538 = t66 * qJD(2);
t149 = t285 * t612;
t288 = -t578 - t579;
t68 = -t288 * t435 + t149;
t537 = t68 * qJD(1);
t150 = t285 * t435;
t69 = t288 * t612 + t150;
t536 = t69 * qJD(1);
t235 = t524 * t319;
t74 = t235 * t303 + t302 * t317;
t533 = t74 * qJD(2);
t75 = -t435 * t578 - t149;
t532 = t75 * qJD(1);
t76 = t578 * t612 - t150;
t531 = t76 * qJD(1);
t451 = t435 * t587;
t453 = -t548 / 0.2e1;
t527 = pkin(4) * t451 + pkin(6) * t453;
t348 = pkin(4) * t584;
t449 = -t622 / 0.2e1;
t526 = pkin(6) * t449 + t348 * t435;
t299 = t311 * qJD(4);
t312 = t317 * qJD(3);
t525 = -t312 - t299;
t523 = qJD(1) * t435;
t521 = qJD(1) * t285;
t520 = qJD(1) * t321;
t519 = qJD(1) * t344;
t518 = qJD(1) * t358;
t517 = qJD(2) * t198;
t516 = qJD(2) * t235;
t515 = qJD(2) * t310;
t514 = qJD(2) * t311;
t513 = qJD(2) * t317;
t512 = qJD(2) * t319;
t511 = qJD(2) * t353;
t510 = qJD(2) * t357;
t509 = qJD(3) * t343;
t508 = qJD(3) * t344;
t507 = qJD(4) * t357;
t506 = qJD(5) * t353;
t447 = t602 + t598;
t110 = t447 * t353;
t503 = t110 * qJD(1);
t502 = t111 * qJD(1);
t446 = 0.2e1 * t624;
t113 = t446 * t353;
t103 = t113 * qJD(1);
t116 = t447 * t357;
t501 = t116 * qJD(1);
t500 = t118 * qJD(1);
t120 = t446 * t357;
t106 = t120 * qJD(1);
t125 = t337 * t250;
t499 = t125 * qJD(1);
t132 = t310 * t317 + t311 * t319;
t498 = t132 * qJD(2);
t200 = t320 ^ 2 - t321 ^ 2;
t497 = t200 * qJD(1);
t233 = -t320 * t579 + t344 * t321;
t496 = t233 * qJD(1);
t234 = -t344 * t320 - t321 * t579;
t495 = t234 * qJD(1);
t248 = t307 + t421;
t493 = t248 * qJD(1);
t490 = t612 * qJD(4);
t489 = t269 * qJD(2);
t427 = -t480 / 0.2e1;
t399 = -t327 / 0.2e1 + t427;
t400 = -t422 / 0.2e1;
t271 = pkin(2) * t400 - t399;
t488 = t271 * qJD(2);
t338 = -t356 ^ 2 + t358 ^ 2;
t487 = t338 * qJD(1);
t486 = t356 * qJD(2);
t485 = t358 * qJD(2);
t481 = t344 * t579;
t479 = pkin(1) * t356 * qJD(1);
t478 = pkin(1) * t518;
t477 = t354 * t572;
t476 = t354 * t571;
t475 = t88 * t523;
t472 = t581 * t351;
t471 = t581 * t352;
t469 = t351 * t522;
t468 = t352 * t522;
t466 = t353 * t507;
t465 = t435 * t506;
t464 = t435 * t349;
t130 = t612 * t523;
t462 = t435 * t521;
t461 = t612 * t521;
t460 = t288 * t521;
t459 = t320 * t520;
t458 = t320 * t519;
t457 = t321 * t519;
t456 = t356 * t485;
t324 = t343 * t586;
t443 = t589 + t591;
t442 = t582 * qJD(2);
t441 = t582 * qJD(3);
t440 = t581 * qJD(3);
t439 = t581 * qJD(4);
t430 = t285 * pkin(3) * t520;
t429 = qJD(1) * t481;
t428 = -qJD(5) + t523;
t426 = t250 * t339;
t424 = t612 * t467;
t419 = t269 / 0.2e1 + t595;
t418 = t603 + t445;
t417 = t483 * t269;
t414 = t443 * t310;
t413 = -0.2e1 * t424;
t412 = 0.2e1 * t424;
t411 = t594 + t592 + t597;
t407 = t435 * t424;
t393 = -t317 / 0.2e1 + t419;
t175 = t393 * t357;
t363 = t342 * t598 - t435 * t592 + (t581 * t601 + t585 * t612) * pkin(3);
t360 = t363 + t397;
t31 = t360 * t353;
t405 = t31 * qJD(1) - t175 * qJD(2);
t378 = t472 + t471;
t244 = (t342 * t378 + t343 * t354) * pkin(3);
t375 = t471 / 0.2e1 + t472 / 0.2e1;
t359 = (t302 * t585 + t303 * t375) * pkin(3) + t342 * t414 + t343 * t595;
t41 = t574 / 0.2e1 + t443 * t319 * pkin(6) + t359;
t404 = -t41 * qJD(2) - t244 * qJD(3);
t318 = t378 * pkin(3);
t81 = t524 * (t480 / 0.2e1 - t310 / 0.2e1 + t594);
t403 = -qJD(2) * t81 - qJD(3) * t318;
t402 = t428 * t353;
t401 = t428 * t357;
t396 = -t576 / 0.2e1 + t575 / 0.2e1;
t392 = qJD(4) * t416;
t391 = t303 * t624 + t597 * t612;
t389 = -t435 * t593 + t592 * t612;
t15 = t361 * t357;
t388 = -qJD(1) * t15 - t317 * t511;
t21 = t362 * t357;
t387 = -qJD(1) * t21 - t311 * t511;
t377 = t91 / 0.2e1 + t391;
t46 = t377 * t357;
t386 = qJD(1) * t46 - t302 * t511;
t43 = t377 * t353;
t385 = -qJD(1) * t43 - t302 * t510;
t72 = -t130 + t619;
t384 = -t621 + t619;
t383 = qJD(2) * (t556 - t558);
t382 = qJD(3) * (-t342 * t612 + t343 * t435);
t381 = t427 + t603 + t592;
t379 = t128 / 0.2e1 + t396;
t376 = t92 / 0.2e1 + t389;
t373 = t620 * t612;
t372 = (t433 + t490) * t435;
t161 = t411 * t353;
t50 = t376 * t357;
t371 = qJD(1) * t50 + qJD(2) * t161 - t353 * t509;
t162 = t411 * t357;
t47 = t376 * t353;
t370 = -qJD(1) * t47 + qJD(2) * t162 - t357 * t509;
t174 = t393 * t353;
t34 = t360 * t357;
t369 = -qJD(1) * t34 - qJD(2) * t174 - t353 * t477;
t169 = t418 * t353;
t268 = t381 * t353;
t56 = t379 * t357;
t367 = qJD(1) * t56 + qJD(2) * t169 + qJD(3) * t268 + qJD(4) * t577;
t170 = t418 * t357;
t270 = t381 * t357;
t53 = t379 * t353;
t366 = pkin(4) * t507 - qJD(1) * t53 + qJD(2) * t170 + qJD(3) * t270;
t364 = -t558 / 0.2e1 + t556 / 0.2e1 + t390;
t347 = -t577 / 0.2e1;
t340 = t356 * t518;
t336 = t353 * t476;
t326 = t337 * qJD(5);
t325 = t343 * t583;
t314 = t319 * qJD(3);
t313 = t318 * qJD(4);
t298 = t310 * qJD(4);
t295 = t353 * t312;
t284 = t353 * t299;
t277 = t302 * t583;
t276 = t302 * t586;
t273 = t357 * t427 + t325 + t348;
t272 = t353 * t427 + t324 + t347;
t259 = (t400 + t545) * pkin(2) + t399;
t258 = -t269 / 0.2e1 - t398 / 0.2e1 + (-t470 - t474 / 0.2e1) * pkin(2);
t256 = t484 * t320;
t242 = t248 * qJD(4);
t215 = t235 * qJD(3);
t199 = t339 * t617;
t191 = t198 * qJD(4);
t176 = t317 * t584 - t357 * t419;
t173 = t317 * t586 + t353 * t419;
t172 = t310 * t583 + t277 + t348;
t171 = t310 * t586 + t276 + t347;
t166 = -t546 + 0.2e1 * t307;
t164 = t319 * t584 + t277 + t325;
t163 = t319 * t587 + t276 + t324;
t127 = t413 + t615;
t126 = t412 - t615;
t119 = -t118 / 0.2e1 + t448;
t115 = t612 * t586 + t548 / 0.2e1;
t114 = t452 + t451;
t94 = t349 + t106;
t93 = -t103 - t506;
t80 = pkin(3) * t375 + t414 + t610;
t67 = t632 * t611;
t64 = t166 * qJD(4) - t433;
t58 = (-t444 + t588 + t591) * t435;
t55 = t128 * t583 - t357 * t396;
t54 = t128 * t587 + t353 * t396;
t49 = -t357 * t389 + t583 * t92;
t48 = t353 * t389 + t587 * t92;
t45 = -t357 * t391 + t583 * t91;
t44 = t353 * t391 + t587 * t91;
t40 = -t574 / 0.2e1 + t359 + t610 * pkin(6);
t38 = -t608 - t627;
t37 = t608 - t628;
t33 = t357 * t363 + t526;
t32 = t353 * t363 + t527;
t26 = t199 + 0.2e1 * t407;
t23 = t374 + t365;
t20 = t357 * t368 + t526;
t19 = t353 * t368 + t527;
t17 = -qJD(4) * t110 + qJD(5) * t120 - t634;
t16 = -qJD(4) * t116 - qJD(5) * t113 + t633;
t14 = t342 * t449 + t343 * t448 + t357 * t364;
t13 = t324 * t435 + t342 * t453 + t353 * t364;
t11 = (t468 + t623) * t435 - t609;
t10 = (t469 - t623) * t435 + t609;
t6 = t115 * qJD(4) + t119 * qJD(5) + t353 * t433 + t634;
t5 = qJD(4) * t622 + t114 * qJD(5) + t612 * t613 - t633;
t4 = t58 * qJD(4) + t337 * t432 + t199 - 0.2e1 * t407;
t25 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t456, t338 * qJD(2), 0, -t456, 0, 0, -pkin(1) * t486, -pkin(1) * t485, 0, 0, -t616, t484 * t200, 0, t616, 0, 0, -qJD(2) * t233 - t321 * t508, -qJD(2) * t234 + t320 * t508, 0, -qJD(2) * t481, t373, -t632 * t626, 0, -t372, 0, 0, qJD(2) * t68 - qJD(3) * t75 + t285 * t490, qJD(2) * t69 - qJD(3) * t76 + t285 * t491, 0 (qJD(2) * t288 - t321 * t572) * t285, t352 * t373 - t426, -0.2e1 * t353 * t620 * t622 - t125 * qJD(5), t465 * t612 + t630 * t632, t351 * t373 + t426, t464 * t612 - t631 * t632, -t372, qJD(2) * t28 + qJD(3) * t30 + qJD(4) * t36 + t465 * t88, -qJD(2) * t27 - qJD(3) * t29 - qJD(4) * t35 + t464 * t88, -qJD(2) * t1 - qJD(3) * t2 - qJD(4) * t3, qJD(2) * t7 + qJD(3) * t8 + qJD(4) * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t340, t487, -t485, -t340, t486, 0, -t479, -t478, 0, 0, -t459, t497, t256, t459, t257, 0, -t496, -t495 (-t320 * t582 + t321 * t355) * t573, -t429, t621, -t629, t620, -t130, t64, 0, t537, t536 (-t550 + t551) * qJD(2) + t23 * qJD(3), t460, t11, t4, t6, t10, t5, t72, t13 * qJD(3) + t19 * qJD(4) + t45 * qJD(5) + t353 * t383 + t554, t14 * qJD(3) + t20 * qJD(4) + t44 * qJD(5) + t357 * t383 - t555, -t570, t566; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t459, t497, t256, t459, t257, 0, -t457, t458, 0, 0, t621, -t629, t620, -t130, t64, 0, -t532, -t531, t23 * qJD(2) + (-t473 - t557) * t572, -t430, t11, t4, t6, t10, t5, t72, t13 * qJD(2) + t32 * qJD(4) + t49 * qJD(5) + t353 * t382 + t552, t14 * qJD(2) + t33 * qJD(4) + t48 * qJD(5) + t357 * t382 - t553, -t569, t565; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t621, -t629, t620, -t621, t166 * t484 - t490, 0, t461, t462, 0, 0 -(-t466 - t468) * t435 + t607, t199 + t484 * t58 - (-qJD(4) * t337 + t412) * t435, t115 * t484 + t353 * t490 + t634 -(t466 - t469) * t435 - t607, t357 * t490 + t484 * t622 - t633, t384, t19 * qJD(2) + t32 * qJD(3) + t55 * qJD(5) + t353 * t392 + t542, t20 * qJD(2) + t33 * qJD(3) + t54 * qJD(5) + t357 * t392 - t549, -t567, t564; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t614 * t617 - t499, t119 * t484 + t402 * t612, t24, t114 * t484 + t401 * t612, t67, t45 * qJD(2) + t49 * qJD(3) + t55 * qJD(4) + t402 * t88, t44 * qJD(2) + t48 * qJD(3) + t54 * qJD(4) + t401 * t88, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t340, -t487, 0, t340, 0, 0, t479, t478, 0, 0, t459, -t497, 0, -t459, 0, 0, t496, t495, 0, t429, -t621, t629, 0, t130, -t242, 0, -t537, -t536, -qJD(3) * t22, -t460, t38, t26, t17, t37, t16, -t72, qJD(3) * t12 + qJD(4) * t18 - qJD(5) * t46 - t554, qJD(3) * t15 + qJD(4) * t21 + qJD(5) * t43 + t555, t570, -t566; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t580, -pkin(2) * t441, 0, 0, 0, 0, 0, 0, 0, 0, t525, -t314 + t298, 0, qJD(3) * t132, t339, t326, 0, -t339, 0, 0, t302 * t506 + t357 * t525, t302 * t349 + t284 + t295, t215 - t191, qJD(3) * t74 + qJD(4) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t484 * t580 (-t442 - t441) * pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t258 - t312 - t513, qJD(4) * t259 - t314 - t512, -t559, t498 + (-t317 * t581 + t319 * t354) * t572, t339, t326, 0, -t339, 0, 0, t176 * qJD(4) + t163 * qJD(5) - t317 * t613 + t561, qJD(4) * t173 + qJD(5) * t164 + t295 - t388, qJD(4) * t80 + t215 + t516, t533 + (t235 * t342 + t317 * t343) * qJD(3) + t40 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t493, 0, qJD(3) * t258 - t299 - t514, qJD(3) * t259 + t298 + t515, 0, 0, t84, t326, -t503, -t84, -t501, 0, t560 + t176 * qJD(3) + t171 * qJD(5) + (-qJD(2) - qJD(4)) * t357 * t311, qJD(3) * t173 + qJD(5) * t172 + t284 - t387, qJD(3) * t80 - t191 - t517, t538 + t40 * qJD(3) + (-t311 * pkin(4) - pkin(6) * t198) * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t127, t94, -t79, t93, -t618, qJD(3) * t163 + qJD(4) * t171 - t303 * t349 - t386, qJD(3) * t164 + qJD(4) * t172 + t303 * t506 - t385, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t459, -t497, 0, -t459, 0, 0, t457, -t458, 0, 0, -t621, t629, 0, t130, -t242, 0, t532, t531, qJD(2) * t22, t430, t38, t26, t17, t37, t16, -t72, -qJD(2) * t12 + qJD(4) * t31 - qJD(5) * t50 - t552, -qJD(2) * t15 + qJD(4) * t34 + qJD(5) * t47 + t553, t569, -t565; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t355 * t573, pkin(2) * t442, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t269 + t513, -qJD(4) * t271 + t512, t559, -t498, t339, t326, 0, -t339, 0, 0, -qJD(4) * t175 - qJD(5) * t161 + t317 * t510 - t561, qJD(4) * t174 - qJD(5) * t162 + t388, qJD(4) * t81 - t516, qJD(4) * t41 - t533; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t476, -pkin(3) * t439, 0, 0, t339, t326, 0, -t339, 0, 0, t343 * t506 - t357 * t476, t343 * t349 + t336, t313, t244 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t493, 0, t417 - t489, -t488 + (-t440 - t439) * pkin(3), 0, 0, t84, t326, -t503, -t84, -t501, 0, t272 * qJD(5) + t357 * t417 + t405, qJD(5) * t273 + t336 - t369, t313 - t403 (-pkin(4) * t354 + pkin(6) * t378) * t571 - t404; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t127, t94, -t79, t93, -t618, qJD(4) * t272 - t342 * t349 - t371, qJD(4) * t273 + t342 * t506 - t370, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t621, t629, 0, t621, t484 * t248, 0, -t461, -t462, 0, 0, -t606 - t627, -t413 * t435 + t199, -qJD(5) * t118 + t110 * t484 - t634, t606 - t628, qJD(5) * t111 + t116 * t484 + t633, -t384, -qJD(2) * t18 - qJD(3) * t31 - qJD(5) * t56 - t542, -qJD(2) * t21 - qJD(3) * t34 + qJD(5) * t53 + t549, t567, -t564; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t493, 0, qJD(3) * t269 + t514, qJD(3) * t271 - t515, 0, 0, -t85, t326, t503, t85, t501, 0, qJD(3) * t175 - qJD(5) * t169 + t311 * t510 - t560, -qJD(3) * t174 - qJD(5) * t170 + t387, -qJD(3) * t81 + t517, -qJD(3) * t41 - t538; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t493, 0, t477 + t489, pkin(3) * t440 + t488, 0, 0, -t85, t326, t503, t85, t501, 0, -qJD(5) * t268 + t357 * t477 - t405, -qJD(5) * t270 + t369, t403, t404; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t339, t326, 0, -t339, 0, 0, -pkin(4) * t506, -pkin(4) * t349, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t127, t349 - t500, -t79, t502 - t506, -t618, -pkin(6) * t349 - t367, pkin(6) * t506 - t366, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0.2e1 * t612 * t614 + t499, t118 * qJD(4) - t120 * t484 - t353 * t621, -t24, -t111 * qJD(4) + t113 * t484 - t357 * t621, t67, qJD(2) * t46 + qJD(3) * t50 + qJD(4) * t56 - t353 * t475, -qJD(2) * t43 - qJD(3) * t47 - qJD(4) * t53 - t357 * t475, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t126, -t106, t79, t103, t618, qJD(3) * t161 + qJD(4) * t169 + t386, qJD(3) * t162 + qJD(4) * t170 + t385, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t126, -t106, t79, t103, t618, qJD(4) * t268 + t371, qJD(4) * t270 + t370, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t126, t500, t79, -t502, t618, t367, t366, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t25;
