% Calculate minimal parameter regressor of coriolis matrix for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x24]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6PRRPRP2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:32:18
% EndTime: 2019-03-08 21:32:36
% DurationCPUTime: 8.96s
% Computational Cost: add. (8097->519), mult. (18167->720), div. (0->0), fcn. (20807->10), ass. (0->410)
t363 = cos(qJ(5));
t358 = sin(pkin(11));
t361 = sin(qJ(3));
t364 = cos(qJ(3));
t591 = cos(pkin(11));
t392 = -t358 * t361 + t364 * t591;
t327 = t392 ^ 2;
t456 = t591 * t361;
t554 = t358 * t364;
t331 = t456 + t554;
t328 = t331 ^ 2;
t490 = -t328 - t327;
t637 = t490 * t363;
t642 = qJD(2) * t637;
t360 = sin(qJ(5));
t638 = t490 * t360;
t641 = qJD(2) * t638;
t640 = qJD(4) * t637;
t639 = qJD(4) * t638;
t324 = t331 * qJ(6);
t597 = t361 * pkin(3);
t245 = pkin(4) * t331 - pkin(9) * t392 + t597;
t224 = t360 * t245;
t595 = -qJ(4) - pkin(8);
t340 = t595 * t361;
t342 = t595 * t364;
t276 = -t591 * t340 - t342 * t358;
t634 = t276 * t363;
t539 = t224 / 0.2e1 - t634 / 0.2e1;
t636 = -t324 - t539;
t359 = sin(pkin(6));
t362 = sin(qJ(2));
t553 = t359 * t362;
t592 = cos(pkin(6));
t321 = t361 * t553 - t364 * t592;
t322 = t361 * t592 + t364 * t553;
t217 = t591 * t321 + t322 * t358;
t614 = t217 / 0.2e1;
t457 = t614 - t217 / 0.2e1;
t242 = t363 * t392;
t463 = t242 / 0.2e1;
t635 = t217 * t463;
t632 = t358 * t340 - t591 * t342;
t566 = t632 * t360;
t565 = t632 * t363;
t633 = 0.2e1 * t360;
t238 = t360 * t392;
t498 = t238 * qJD(2);
t504 = qJD(5) * t360;
t631 = -t498 + t504;
t548 = t363 * qJ(6);
t600 = pkin(5) * t360;
t341 = t548 - t600;
t551 = t360 * qJ(6);
t596 = t363 * pkin(5);
t438 = t551 + t596;
t235 = t438 * t331;
t351 = pkin(3) * t358 + pkin(9);
t557 = t351 * t392;
t439 = -t235 / 0.2e1 + t557 / 0.2e1;
t352 = -pkin(3) * t591 - pkin(4);
t325 = -t438 + t352;
t608 = t331 / 0.2e1;
t467 = t325 * t608;
t240 = t360 * t331;
t243 = t363 * t331;
t444 = pkin(5) * t240 - qJ(6) * t243;
t158 = t276 + t444;
t617 = -t158 / 0.2e1;
t369 = (t341 * t608 + t617) * t363 + (t467 + t439) * t360;
t25 = t369 + t636;
t558 = t341 * t360;
t253 = t325 * t363 - t558;
t109 = t457 * t363;
t500 = t109 * qJD(1);
t630 = -qJD(2) * t25 + qJD(3) * t253 - t500;
t513 = qJD(3) * t363;
t556 = t352 * t331;
t374 = (-t557 / 0.2e1 - t556 / 0.2e1) * t360 + t634 / 0.2e1;
t64 = t374 + t539;
t629 = -qJD(2) * t64 - t352 * t513 + t500;
t107 = t457 * t360;
t501 = t107 * qJD(1);
t503 = qJD(5) * t363;
t442 = -t351 * t503 - t501;
t393 = -t358 * t321 + t322 * t591;
t365 = cos(qJ(2));
t552 = t359 * t365;
t185 = t360 * t393 + t363 * t552;
t441 = t457 * t392;
t609 = -t331 / 0.2e1;
t372 = (t393 * t608 - t441) * t360 + t185 * t609;
t280 = t331 * t552;
t561 = t280 * t363;
t368 = t561 / 0.2e1 + t372;
t628 = qJD(1) * t368;
t461 = -t243 / 0.2e1;
t186 = -t360 * t552 + t363 * t393;
t576 = t186 * t392;
t391 = -t576 / 0.2e1 + t217 * t461;
t465 = t553 / 0.2e1;
t281 = t392 * t552;
t549 = t360 * t281;
t534 = -t549 / 0.2e1 + t363 * t465;
t379 = t391 + t534;
t627 = qJD(1) * t379;
t546 = t363 * t281;
t226 = t360 * t553 + t546;
t626 = qJD(2) * (t226 * t392 + t243 * t280);
t625 = -qJD(2) * t368 - t107 * qJD(5);
t98 = t107 * qJD(3);
t624 = qJD(2) * t379 + t98;
t603 = t360 / 0.2e1;
t468 = t217 * t603;
t604 = -t360 / 0.2e1;
t108 = -t217 * t604 + t468;
t370 = -t561 / 0.2e1 + t372;
t623 = qJD(2) * t370 + t108 * qJD(5) - t393 * t513;
t460 = t243 / 0.2e1;
t446 = t217 * t460 + t576 / 0.2e1;
t407 = t446 + t534;
t593 = t108 * qJD(3) - t186 * qJD(5);
t622 = qJD(2) * t407 + t593;
t621 = qJD(3) * t368 - qJD(5) * t379;
t225 = -t363 * t553 + t549;
t620 = qJD(3) * t370 + qJD(5) * t407 + (t225 * t392 + t280 * t240) * qJD(2);
t619 = -qJ(6) / 0.2e1;
t547 = t363 * t245;
t550 = t276 * t360;
t598 = t331 * pkin(5);
t119 = -t547 - t550 - t598;
t618 = t119 / 0.2e1;
t616 = t393 / 0.2e1;
t612 = t235 / 0.2e1;
t611 = -t245 / 0.2e1;
t610 = t392 / 0.2e1;
t356 = t360 ^ 2;
t606 = -t356 / 0.2e1;
t357 = t363 ^ 2;
t605 = -t357 / 0.2e1;
t602 = -t363 / 0.2e1;
t601 = t363 / 0.2e1;
t599 = t392 * pkin(5);
t594 = qJD(3) * pkin(3);
t353 = -pkin(3) * t364 - pkin(2);
t386 = -pkin(4) * t392 - pkin(9) * t331 + t353;
t137 = t360 * t386 + t565;
t560 = t392 * qJ(6);
t116 = t137 - t560;
t459 = -t116 / 0.2e1 + t137 / 0.2e1;
t400 = -t560 / 0.2e1 - t459;
t136 = -t363 * t386 + t566;
t117 = t136 + t599;
t458 = -t136 / 0.2e1 + t117 / 0.2e1;
t408 = -t599 / 0.2e1 - t458;
t12 = t360 * t408 - t363 * t400;
t590 = qJD(2) * t12;
t586 = t116 * t363;
t10 = -t137 * t243 + (t586 + (t117 - t136) * t360) * t331;
t589 = t10 * qJD(2);
t11 = t360 * t400 + t363 * t408;
t588 = t11 * qJD(2);
t587 = t116 * t392;
t585 = t117 * t360;
t584 = t136 * t392;
t583 = t137 * t392;
t157 = -t341 * t392 + t632;
t582 = t157 * t360;
t581 = t157 * t363;
t580 = t158 * t360;
t579 = t158 * t363;
t578 = t185 * t392;
t577 = t185 * t360;
t575 = t186 * t363;
t573 = t217 * t280;
t572 = t217 * t341;
t571 = t225 * t360;
t570 = t225 * t363;
t569 = t226 * t360;
t568 = t226 * t363;
t26 = (t393 - t575 - t577) * t217;
t567 = t26 * qJD(1);
t563 = t280 * t325;
t562 = t280 * t360;
t559 = t331 * t351;
t378 = -(t605 + t606) * t557 + t467;
t536 = -t634 + t224;
t118 = t324 + t536;
t404 = t118 * t604 + t119 * t601;
t40 = t378 + t404;
t545 = t40 * qJD(2);
t41 = t185 * t225 + t186 * t226 + t573;
t544 = t41 * qJD(1);
t45 = (-t393 / 0.2e1 + t616) * t331 - t441;
t543 = t45 * qJD(2);
t71 = -t359 ^ 2 * t362 * t365 + t281 * t393 + t573;
t541 = t71 * qJD(1);
t540 = t580 / 0.2e1 + t325 * t460;
t537 = t547 / 0.2e1 + t550 / 0.2e1;
t466 = -t553 / 0.2e1;
t535 = t549 / 0.2e1 + t363 * t466;
t533 = t546 / 0.2e1 + t360 * t465;
t532 = -t546 / 0.2e1 + t360 * t466;
t530 = t331 * t606 + t357 * t608;
t346 = t357 - t356;
t489 = t328 - t327;
t192 = t489 * t360;
t528 = qJD(2) * t192;
t194 = t489 * t363;
t526 = qJD(2) * t194;
t522 = qJD(2) * t392;
t521 = qJD(2) * t331;
t520 = qJD(2) * t360;
t519 = qJD(2) * t362;
t518 = qJD(2) * t363;
t517 = qJD(2) * t364;
t516 = qJD(3) * t109;
t515 = qJD(3) * t360;
t514 = qJD(3) * t361;
t512 = qJD(3) * t364;
t511 = qJD(3) * t365;
t510 = qJD(4) * t238;
t509 = qJD(4) * t360;
t508 = qJD(4) * t363;
t507 = qJD(5) * t109;
t506 = qJD(5) * t136;
t505 = qJD(5) * t392;
t502 = qJD(6) * t360;
t485 = -t591 / 0.2e1;
t384 = t331 * t485 + t358 * t610;
t205 = (-t361 / 0.2e1 + t384) * pkin(3);
t499 = t205 * qJD(2);
t497 = t240 * qJD(2);
t234 = t242 * qJD(2);
t496 = t243 * qJD(2);
t316 = t356 * t392;
t317 = t357 * t392;
t246 = -t316 - t317;
t495 = t246 * qJD(2);
t494 = t490 * qJD(2);
t326 = t456 / 0.2e1 + t554 / 0.2e1;
t493 = t326 * qJD(2);
t323 = t392 * qJD(6);
t347 = -t361 ^ 2 + t364 ^ 2;
t492 = t347 * qJD(2);
t491 = t363 * qJD(6);
t488 = pkin(2) * t361 * qJD(2);
t487 = pkin(2) * t517;
t318 = t598 / 0.2e1;
t486 = t361 * t552;
t462 = -t242 / 0.2e1;
t484 = -t580 / 0.2e1 + t325 * t461 + t351 * t462;
t483 = t318 + t537;
t482 = t357 * t521;
t481 = t392 * t504;
t480 = t392 * t503;
t479 = t351 * t504;
t477 = t360 * t491;
t476 = t392 * t331 * qJD(3);
t475 = qJD(2) * t552;
t474 = t360 * t503;
t473 = t331 * t520;
t472 = t331 * t502;
t348 = t360 * t513;
t471 = t361 * t517;
t470 = t331 * t518;
t313 = t331 * t513;
t469 = t217 * t608;
t464 = -t240 / 0.2e1;
t455 = (t356 + t357) * t217;
t188 = t326 + t530;
t454 = qJD(2) * t188 + t348;
t236 = (t356 / 0.2e1 + t605) * t331;
t453 = qJD(2) * t236 - t348;
t292 = t328 * t360 * t518;
t452 = qJD(3) * t236 + t292;
t451 = qJD(5) - t522;
t450 = t360 * t470;
t449 = t360 * t313;
t448 = t331 * t468 + t578 / 0.2e1;
t447 = t217 * t464 - t578 / 0.2e1;
t445 = -t598 / 0.2e1 - t550 / 0.2e1;
t440 = 0.2e1 * t449;
t405 = t586 / 0.2e1 + t585 / 0.2e1;
t366 = -t405 * t217 + t185 * t618 + t186 * t118 / 0.2e1 + t158 * t616 + t157 * t614;
t401 = -t568 / 0.2e1 - t571 / 0.2e1;
t2 = -t563 / 0.2e1 + t401 * t351 + t366;
t7 = t116 * t118 + t117 * t119 + t157 * t158;
t437 = t2 * qJD(1) + t7 * qJD(2);
t371 = t185 * t459 + t186 * t458 + t217 * t612;
t410 = t225 * pkin(5) / 0.2e1 + t226 * t619;
t4 = t371 + t410;
t8 = -t116 * t136 + t117 * t137 + t158 * t235;
t436 = t4 * qJD(1) + t8 * qJD(2);
t383 = (t185 * t602 + t186 * t603) * t392;
t43 = -t383 + t401;
t9 = -t117 * t242 - t119 * t243 + (t118 * t331 + t587) * t360;
t435 = t43 * qJD(1) - t9 * qJD(2);
t15 = (t116 - t581) * t331 - (t118 + t579) * t392;
t269 = t562 / 0.2e1;
t390 = t186 * t609 - t635;
t380 = t217 * t462 + t393 * t461 - t390;
t29 = t269 + t380;
t434 = t29 * qJD(1) + t15 * qJD(2);
t16 = (-t117 + t582) * t331 - (-t119 - t580) * t392;
t433 = t16 * qJD(2) + t628;
t373 = t457 * t632;
t385 = t281 * t358 / 0.2e1 + t280 * t485;
t17 = (t486 / 0.2e1 + t385) * pkin(3) + t373;
t81 = t353 * t597;
t432 = -t17 * qJD(1) + t81 * qJD(2);
t21 = t158 * t331 - (-t585 - t586) * t392;
t375 = -(-t575 / 0.2e1 - t577 / 0.2e1) * t392 + t469;
t402 = t570 / 0.2e1 - t569 / 0.2e1;
t23 = t375 + t402;
t431 = qJD(1) * t23 + qJD(2) * t21;
t37 = (-t136 + t566) * t331 - t547 * t392;
t430 = t37 * qJD(2) + t628;
t270 = -t562 / 0.2e1;
t381 = t393 * t460 + t390 + t635;
t32 = t270 + t381;
t38 = (-t137 + t565) * t331 - (-t536 - t634) * t392;
t429 = t32 * qJD(1) + t38 * qJD(2);
t428 = t45 * qJD(1);
t46 = t583 + (t235 * t360 + t579) * t331;
t427 = qJD(2) * t46 - t627;
t47 = t584 + (-t235 * t363 + t580) * t331;
t56 = t447 + t533;
t426 = qJD(1) * t56 - qJD(2) * t47;
t51 = -t158 * t243 - t587;
t53 = t446 + t535;
t425 = -qJD(1) * t53 + qJD(2) * t51;
t78 = t243 * t276 + t583;
t424 = qJD(2) * t78 - t627;
t57 = t448 + t532;
t77 = -t240 * t276 - t584;
t423 = -qJD(1) * t57 + qJD(2) * t77;
t422 = t118 * t363 + t119 * t360;
t421 = -t325 * t392 + t559;
t420 = t352 * t392 - t559;
t120 = t276 * t331 + t392 * t632;
t403 = t393 * t610 + t469;
t82 = t465 - t403;
t418 = -qJD(1) * t82 + qJD(2) * t120;
t19 = (-t558 / 0.2e1 - pkin(5) / 0.2e1) * t331 + (t611 + t439) * t363 + t445 + t540;
t254 = -t325 * t360 - t341 * t363;
t417 = -qJD(2) * t19 + qJD(3) * t254;
t415 = t451 * t360;
t414 = t451 * t363;
t110 = (t601 - t602) * t217;
t413 = -qJD(3) * t110 - qJD(5) * t185;
t412 = -qJD(5) * t341 - t502;
t411 = pkin(5) * t618 + t118 * t619;
t409 = t548 / 0.2e1 - t600 / 0.2e1;
t50 = t483 + t484;
t399 = -qJD(2) * t50 + t325 * t515;
t287 = t351 * t463;
t66 = t287 + (t556 / 0.2e1 + t611) * t363;
t398 = -qJD(2) * t66 - t352 * t515;
t396 = t331 * t414;
t395 = qJD(5) * t326 - t392 * t521;
t394 = -qJD(5) * t110 - t393 * t515;
t247 = t346 * t328;
t389 = qJD(2) * t247 + t440;
t388 = -qJD(3) * t346 + 0.2e1 * t450;
t387 = t409 * t217;
t35 = t572 / 0.2e1 - t387;
t367 = (t360 * t459 + t363 * t458) * t351 + t341 * t617 + t325 * t612;
t6 = t367 + t411;
t382 = t325 * t341 * qJD(3) + t35 * qJD(1) - t6 * qJD(2);
t377 = -qJD(5) * t438 + t491;
t260 = t328 * t357 + t327;
t376 = qJD(2) * t260 + t449 - t505;
t315 = t326 * qJD(3);
t312 = t331 * t515;
t293 = t363 * t472;
t291 = t451 * qJ(6);
t279 = qJD(3) * t356 + t450;
t233 = t242 * qJD(5);
t229 = t238 * qJD(5);
t227 = t236 * qJD(5);
t204 = t597 / 0.2e1 + t384 * pkin(3);
t203 = -t234 + t503;
t189 = -t326 + t530;
t138 = qJD(3) * t242 - t392 * t473;
t129 = t137 * qJD(5);
t112 = t331 * t415;
t83 = t465 + t403;
t67 = t276 * t603 + t352 * t460 + t287 + t537;
t65 = t374 - t539;
t60 = t391 + t535;
t59 = t447 + t532;
t58 = t448 + t533;
t49 = -t547 / 0.2e1 + t445 + t484;
t44 = t45 * qJD(3);
t42 = -t383 - t401;
t39 = t378 - t404;
t36 = -t572 / 0.2e1 - t387;
t31 = t269 + t381;
t30 = t270 + t380;
t24 = t369 - t636;
t22 = t375 - t402;
t20 = t341 * t464 + t363 * t439 + t318 + t483 + t540;
t18 = -t373 + (-t486 / 0.2e1 + t385) * pkin(3);
t14 = t136 * t602 + t116 * t604 + t137 * t603 + t117 * t601 - (t551 / 0.2e1 + t596 / 0.2e1) * t392;
t13 = t136 * t604 + t137 * t602 + t392 * t409 + t405;
t5 = t367 - t411;
t3 = t371 - t410;
t1 = t563 / 0.2e1 + t366 + (t568 + t571) * t351 / 0.2e1;
t27 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t71, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t41 + qJD(3) * t26; 0, 0, -t359 * t519, -t475, 0, 0, 0, 0, 0 (-t361 * t511 - t362 * t517) * t359 (t361 * t519 - t364 * t511) * t359 (t280 * t331 + t281 * t392) * qJD(2) + t44, t541 + (t276 * t280 + t281 * t632 + t353 * t553) * qJD(2) + t18 * qJD(3) + t83 * qJD(4), 0, 0, 0, 0, 0, t620, t31 * qJD(3) + t59 * qJD(5) + t626, t620, t42 * qJD(3) + (-t569 + t570) * t521, t30 * qJD(3) + t58 * qJD(5) - t626, t544 + (t116 * t226 + t117 * t225 + t158 * t280) * qJD(2) + t1 * qJD(3) + t22 * qJD(4) + t3 * qJD(5) + t60 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t322 * qJD(3) - t361 * t475, t321 * qJD(3) - t364 * t475, t543, t18 * qJD(2) + (-t217 * t358 - t393 * t591) * t594, 0, 0, 0, 0, 0, t623, qJD(2) * t31 - t394, t623, t42 * qJD(2) - qJD(3) * t455, qJD(2) * t30 + t394, t567 + t1 * qJD(2) + (t325 * t393 - t351 * t455) * qJD(3) + t36 * qJD(5) - t108 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t622, qJD(2) * t59 - t413, t622, 0, qJD(2) * t58 + t413, t3 * qJD(2) + t36 * qJD(3) + (-pkin(5) * t186 - qJ(6) * t185) * qJD(5) + t186 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t60 - t593; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -qJD(3) * t17 - qJD(4) * t82 - t541, 0, 0, 0, 0, 0, t621, qJD(3) * t32 - qJD(5) * t57, t621, qJD(3) * t43, qJD(3) * t29 - qJD(5) * t56, qJD(3) * t2 + qJD(4) * t23 + qJD(5) * t4 - qJD(6) * t53 - t544; 0, 0, 0, 0, t361 * t512, t347 * qJD(3), 0, 0, 0, -pkin(2) * t514, -pkin(2) * t512, -qJD(4) * t490, qJD(3) * t81 + qJD(4) * t120, -t328 * t474 + t357 * t476, -qJD(5) * t247 - t392 * t440, qJD(3) * t194 + t331 * t481, -qJD(3) * t192 + t331 * t480, -t476, qJD(3) * t37 + qJD(5) * t78 - t639, qJD(3) * t38 + qJD(5) * t77 - t640, qJD(3) * t16 + qJD(5) * t46 - t328 * t477 - t639, -qJD(3) * t9 - qJD(5) * t10 + t392 * t472, qJD(3) * t15 + qJD(5) * t47 + qJD(6) * t260 + t640, qJD(3) * t7 + qJD(4) * t21 + qJD(5) * t8 + qJD(6) * t51; 0, 0, 0, 0, t471, t492, t512, -t514, 0, -pkin(8) * t512 - t488, pkin(8) * t514 - t487 (-t331 * t358 - t392 * t591) * t594 + t428 (-t276 * t358 - t591 * t632) * t594 + t204 * qJD(4) + t432, -t227 - (-t348 - t482) * t392 (-t316 + t317) * qJD(3) + (-qJD(5) - t522) * t243 * t633, t312 + t526, t313 - t528, t395 (t360 * t420 - t565) * qJD(3) + t67 * qJD(5) + t430 (t363 * t420 + t566) * qJD(3) + t65 * qJD(5) + t429 (-t360 * t421 - t581) * qJD(3) + t20 * qJD(5) + t189 * qJD(6) + t433, qJD(3) * t422 + t14 * qJD(5) + t435 (t363 * t421 - t582) * qJD(3) + t24 * qJD(5) + t293 + t434 (t157 * t325 + t351 * t422) * qJD(3) + t39 * qJD(4) + t5 * qJD(5) + t49 * qJD(6) + t437; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t494, qJD(3) * t204 + t418, 0, 0, 0, 0, 0, -t641, -t642, -t641, 0, t642, qJD(3) * t39 + qJD(5) * t13 + t431; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t452, -t389, -t112, -t396, t315, qJD(3) * t67 - t129 + t424, qJD(3) * t65 + t423 + t506, qJD(3) * t20 - t129 + t427, t14 * qJD(3) + qJD(5) * t444 - t472 - t589, qJD(3) * t24 - t323 - t426 - t506, t5 * qJD(3) + t13 * qJD(4) + (-pkin(5) * t137 - qJ(6) * t136) * qJD(5) + t116 * qJD(6) + t436; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t189 - t292, -t112, t376, qJD(3) * t49 + qJD(5) * t116 + t425; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t543, qJD(2) * t17, 0, 0, 0, 0, 0, t625, -qJD(2) * t32 - t507, t625, -qJD(2) * t43, -qJD(2) * t29 + t507, -qJD(2) * t2 - qJD(5) * t35 + qJD(6) * t107 - t567; 0, 0, 0, 0, -t471, -t492, 0, 0, 0, t488, t487, -t428, qJD(4) * t205 - t432, -t392 * t482 - t227, -t396 * t633, -t233 - t526, t229 + t528, -t395, qJD(5) * t66 - t331 * t508 - t430, qJD(4) * t240 + qJD(5) * t64 - t429, -qJD(4) * t243 + qJD(5) * t19 + qJD(6) * t188 - t433, -qJD(4) * t246 - qJD(5) * t11 - qJD(6) * t242 - t435, qJD(5) * t25 - t331 * t509 + t293 - t434, qJD(4) * t40 + qJD(5) * t6 + qJD(6) * t50 - t437; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t474, t346 * qJD(5), 0, 0, 0, t352 * t504, t352 * t503, -qJD(5) * t254 + t477, 0, -qJD(5) * t253 + qJD(6) * t356, t412 * t325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t499, 0, 0, 0, 0, 0, -t470, t497, -t496, -t495, -t473, t545; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t453, -t388, t203, -t631, -t493, -t398 + t442, t479 - t629, -t417 + t442, t377 - t588, -t479 - t630, t351 * t377 - t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t454, t203, t279, -t399 - t442; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t494, -qJD(3) * t205 - t418, 0, 0, 0, 0, 0, t229 + t313 + t641, -qJD(3) * t240 + t480 + t642, qJD(3) * t243 + t481 + t641, t246 * qJD(3), -t233 + t312 - t642, -qJD(3) * t40 - qJD(5) * t12 - qJD(6) * t238 - t431; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t499, 0, 0, 0, 0, 0, t470, -t497, t496, t495, t473, -t545; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t631, -t414, -t415, 0, t203, -t412 - t590; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t631; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t624, qJD(2) * t57 + t516, t624, 0, qJD(2) * t56 - t516, -qJD(2) * t4 + qJD(3) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t452, t389, t138, -qJD(3) * t238 - t392 * t470, t315, -qJD(3) * t66 - t424 - t510, -qJD(3) * t64 - t392 * t508 - t423, -qJD(3) * t19 - t392 * t509 - t427, qJD(3) * t11 + t589, -qJD(3) * t25 + qJD(4) * t242 - t323 + t426, -qJ(6) * t323 - qJD(3) * t6 + qJD(4) * t12 - t436; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t453, t388, t234, -t498, t493, t398 + t501, t629, t417 + t501, t588, t630, t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t498, -t392 * t518, -t392 * t520, 0, t234, t590; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(6), qJ(6) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t451, t291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t53 - t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t188 + t292, t138, -t376, qJ(6) * t505 - qJD(3) * t50 - t425 + t510; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t454, t234, -t279, t399 - t501; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t498; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t451, -t291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t27;