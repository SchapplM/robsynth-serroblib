% Calculate inertial parameters regressor of coriolis matrix for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6PRPRRP1_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:58:36
% EndTime: 2019-03-08 19:58:48
% DurationCPUTime: 8.63s
% Computational Cost: add. (6913->468), mult. (17465->683), div. (0->0), fcn. (18695->10), ass. (0->388)
t432 = sin(qJ(5));
t427 = t432 ^ 2;
t435 = cos(qJ(5));
t429 = t435 ^ 2;
t406 = t429 + t427;
t433 = sin(qJ(4));
t436 = cos(qJ(4));
t568 = t433 * t436;
t314 = (-0.1e1 + t406) * t568;
t431 = sin(pkin(11));
t612 = sin(pkin(6));
t613 = cos(pkin(11));
t465 = t613 * t612;
t434 = sin(qJ(2));
t480 = t434 * t612;
t631 = cos(qJ(2));
t333 = t431 * t480 - t465 * t631;
t593 = t333 * t432;
t614 = cos(pkin(6));
t479 = t614 * t433;
t470 = t612 * t631;
t334 = t431 * t470 + t434 * t465;
t588 = t334 * t436;
t277 = t479 + t588;
t599 = t277 * t435;
t185 = t593 + t599;
t606 = t185 * t435;
t591 = t333 * t435;
t600 = t277 * t432;
t184 = -t591 + t600;
t608 = t184 * t432;
t453 = t608 / 0.2e1 + t606 / 0.2e1;
t635 = -t429 / 0.2e1;
t518 = 0.1e1 / 0.2e1 + t635;
t276 = t334 * t433 - t436 * t614;
t601 = t276 * t433;
t639 = -t277 / 0.2e1;
t33 = (t639 + t453) * t436 + (-t427 / 0.2e1 + t518) * t601;
t616 = t33 * qJD(1) + t314 * qJD(3);
t416 = pkin(2) * t431 + pkin(8);
t428 = t433 ^ 2;
t430 = t436 ^ 2;
t626 = t436 * pkin(9);
t629 = t433 * pkin(4);
t393 = -t626 + t629;
t379 = t432 * t393;
t569 = t433 * t435;
t513 = t416 * t569;
t302 = t379 - t513;
t597 = t302 * t433;
t380 = t435 * t393;
t571 = t433 * t416;
t382 = t432 * t571;
t301 = t382 + t380;
t598 = t301 * t433;
t417 = -pkin(2) * t613 - pkin(3);
t469 = -t436 * pkin(4) - t433 * pkin(9);
t440 = t469 + t417;
t563 = t436 * t416;
t511 = t435 * t563;
t274 = t432 * t440 + t511;
t602 = t274 * t436;
t339 = t435 * t440;
t516 = t432 * t563;
t273 = -t339 + t516;
t603 = t273 * t436;
t70 = (t597 / 0.2e1 + t602 / 0.2e1) * t435 + (-t598 / 0.2e1 + t603 / 0.2e1) * t432 + (t428 / 0.2e1 - t430 / 0.2e1) * t416;
t660 = -t70 * qJD(2) - t616;
t564 = t436 * qJ(6);
t628 = t433 * pkin(5);
t252 = -t435 * t564 + t301 + t628;
t272 = -t432 * t564 + t302;
t573 = t432 * t433;
t238 = -qJ(6) * t573 + t274;
t630 = t432 * pkin(5);
t481 = t416 + t630;
t343 = t481 * t436;
t567 = t435 * qJ(6);
t512 = t433 * t567;
t223 = -t512 + t339 + (-t416 * t432 - pkin(5)) * t436;
t634 = -t432 / 0.2e1;
t496 = t223 * t634;
t632 = t435 / 0.2e1;
t441 = t238 * t632 + t496 - t343 / 0.2e1;
t342 = t481 * t433;
t638 = t342 / 0.2e1;
t52 = t441 * t436 + (t252 * t634 + t272 * t632 + t638) * t433;
t659 = t52 * qJD(2) + t616;
t658 = 0.2e1 * t432;
t237 = t273 + t512;
t657 = t223 + t237;
t133 = t302 * t436 + (-t274 + 0.2e1 * t511) * t433;
t644 = -t185 / 0.2e1;
t467 = t593 / 0.2e1 + t644;
t59 = (-t599 / 0.2e1 - t467) * t433;
t560 = t59 * qJD(1);
t656 = t133 * qJD(2) - t560;
t132 = t273 * t433 + (t301 - 0.2e1 * t382) * t436;
t645 = t184 / 0.2e1;
t466 = t591 / 0.2e1 + t645;
t60 = (-t600 / 0.2e1 + t466) * t433;
t559 = t60 * qJD(1);
t655 = -t132 * qJD(2) - t559;
t585 = t343 * t435;
t586 = t342 * t435;
t105 = (t272 + t586) * t436 + (-t238 + t585) * t433;
t654 = t105 * qJD(2) - t560;
t566 = t435 * t428;
t210 = -t416 * t566 - t602;
t468 = t334 / 0.2e1 - t601 / 0.2e1;
t79 = t435 * t468 + t436 * t467;
t556 = t79 * qJD(1);
t653 = qJD(2) * t210 + t556;
t143 = t429 * t428 * pkin(5) - t237 * t436 - t342 * t573;
t80 = -t432 * t468 + t436 * t466;
t555 = t80 * qJD(1);
t652 = qJD(2) * t143 - t555;
t515 = t432 * t566;
t517 = t342 * t569;
t604 = t238 * t436;
t137 = -pkin(5) * t515 - t517 - t604;
t651 = qJD(2) * t137 + t556;
t572 = t432 * t436;
t92 = -t223 * t433 + t252 * t436 - t342 * t572 - t343 * t573;
t650 = -t92 * qJD(2) - t559;
t575 = t428 * t432;
t209 = -t416 * t575 - t603;
t649 = -qJD(2) * t209 + t555;
t565 = t435 * t436;
t590 = t334 * t432;
t203 = -t333 * t565 + t590;
t605 = t185 * t436;
t607 = t184 * t436;
t589 = t334 * t435;
t202 = t333 * t572 + t589;
t643 = t202 / 0.2e1;
t45 = (t607 / 0.2e1 - t203 / 0.2e1) * t435 + (-t605 / 0.2e1 + t643) * t432;
t562 = t45 * qJD(1);
t91 = (t598 - t603) * t435 + (t597 + t602) * t432;
t648 = t91 * qJD(2) - t562;
t68 = (t223 * t436 + t252 * t433) * t435 + (t272 * t433 + t604) * t432;
t647 = -t68 * qJD(2) + t562;
t592 = t333 * t433;
t36 = -t184 * t202 + t185 * t203 - t276 * t592;
t452 = t202 * t634 + t203 * t632;
t208 = t333 * t436;
t495 = t208 / 0.2e1;
t90 = (t495 + t452) * t433;
t622 = t36 * qJD(1) + t90 * qJD(3);
t41 = (t277 - t606 - t608) * t276;
t624 = t41 * qJD(1) + t33 * qJD(3);
t646 = pkin(5) / 0.2e1;
t642 = -t223 / 0.2e1;
t641 = -t237 / 0.2e1;
t640 = t252 / 0.2e1;
t625 = pkin(9) + qJ(6);
t388 = t625 * t435;
t637 = -t388 / 0.2e1;
t636 = t427 / 0.2e1;
t633 = -t433 / 0.2e1;
t627 = t436 * pkin(5);
t471 = t642 - t627 / 0.2e1;
t64 = (t641 + t471) * t569;
t618 = t52 * qJD(4) + t64 * qJD(5);
t617 = t70 * qJD(4);
t483 = t237 / 0.2e1 + t223 / 0.2e1;
t53 = (t627 / 0.2e1 - t483) * t435;
t611 = qJD(2) * t53;
t610 = qJD(2) * t64;
t67 = t657 * t573;
t609 = qJD(2) * t67;
t152 = t276 * t432;
t154 = t276 * t435;
t594 = t333 * t416;
t587 = t342 * t432;
t387 = t625 * t432;
t583 = t387 * t433;
t582 = t387 * t436;
t581 = t388 * t436;
t424 = -pkin(5) * t435 - pkin(4);
t579 = t424 * t432;
t578 = t424 * t435;
t577 = t427 * t433;
t576 = t427 * t436;
t574 = t429 * t433;
t570 = t433 * t424;
t83 = (-t277 * t436 + t334 - t601) * t333;
t554 = t83 * qJD(1);
t553 = t90 * qJD(1);
t407 = t429 - t427;
t408 = t430 - t428;
t545 = qJD(4) * t432;
t544 = qJD(4) * t435;
t543 = qJD(5) * t185;
t542 = qJD(5) * t238;
t541 = qJD(5) * t388;
t540 = qJD(5) * t432;
t425 = qJD(5) * t435;
t539 = qJD(5) * t436;
t538 = qJD(6) * t432;
t537 = qJD(6) * t435;
t536 = qJD(6) * t436;
t531 = t333 * qJD(2);
t482 = t636 + t635;
t361 = t482 * t433;
t530 = t361 * qJD(5);
t375 = t406 * t428;
t529 = t375 * qJD(2);
t377 = t408 * t432;
t528 = t377 * qJD(2);
t378 = t430 * t435 - t566;
t527 = t378 * qJD(2);
t526 = t406 * qJD(4);
t525 = t408 * qJD(2);
t524 = t433 * qJD(2);
t523 = t433 * qJD(4);
t522 = t433 * qJD(5);
t521 = t436 * qJD(2);
t520 = t436 * qJD(4);
t365 = t380 / 0.2e1;
t519 = t365 + t382 / 0.2e1;
t514 = t432 * t569;
t510 = t417 * t521;
t509 = t435 * t524;
t508 = t432 * t544;
t507 = t432 * t522;
t506 = t432 * t539;
t505 = t435 * t522;
t504 = t435 * t539;
t503 = t433 * t537;
t411 = t432 * t425;
t502 = t432 * t521;
t501 = t432 * t536;
t500 = t433 * t520;
t499 = t433 * t521;
t498 = t435 * t523;
t497 = t435 * t536;
t309 = t592 / 0.2e1;
t494 = -t587 / 0.2e1;
t493 = -t573 / 0.2e1;
t492 = t573 / 0.2e1;
t491 = t572 / 0.2e1;
t490 = t570 / 0.2e1;
t489 = -t569 / 0.2e1;
t488 = t569 / 0.2e1;
t487 = -t565 / 0.2e1;
t486 = t564 / 0.2e1;
t485 = -t563 / 0.2e1;
t484 = t563 / 0.2e1;
t478 = t406 * t276;
t477 = pkin(5) * t505;
t476 = -qJD(5) + t521;
t475 = t432 * t498;
t474 = t428 * t411;
t473 = t276 * t488;
t472 = -pkin(5) * t514 + t582 / 0.2e1;
t464 = t184 * t488 + t185 * t493;
t463 = -t301 * t432 + t302 * t435;
t292 = t387 * t432 + t388 * t435;
t119 = (t223 * t435 + t238 * t432) * t433;
t74 = t309 + t464;
t462 = qJD(1) * t74 - qJD(2) * t119;
t194 = pkin(5) * t579;
t448 = t424 * t489 + t494;
t25 = t483 * t388 + (t640 + t448) * pkin(5);
t461 = -qJD(2) * t25 + qJD(4) * t194;
t460 = t476 * t433;
t112 = t494 + (-t567 / 0.2e1 + t637) * t436 + (-t578 / 0.2e1 + (0.1e1 - t482) * pkin(5)) * t433 + t519;
t349 = t435 * t630 - t579;
t459 = -qJD(2) * t112 - qJD(4) * t349;
t364 = -t379 / 0.2e1;
t120 = t364 + (t571 / 0.2e1 - t342 / 0.2e1) * t435 + (t486 + t490) * t432 + t472;
t362 = t427 * pkin(5) + t578;
t458 = -qJD(2) * t120 + qJD(4) * t362;
t457 = t626 / 0.2e1 - t629 / 0.2e1;
t450 = t457 * t432;
t278 = t379 / 0.2e1 - t450;
t456 = pkin(4) * t544 - qJD(2) * t278;
t449 = t457 * t435;
t279 = -t380 / 0.2e1 + t449;
t455 = pkin(4) * t545 - qJD(2) * t279;
t454 = t301 * t645 + t302 * t644;
t338 = t435 * t460;
t308 = -t592 / 0.2e1;
t205 = t309 + t308;
t451 = -t205 * qJD(1) - t417 * t524;
t305 = -qJD(2) * t361 + t508;
t373 = t509 + t545;
t371 = t432 * t524 - t544;
t284 = qJD(2) * t515 + qJD(4) * t361;
t376 = t407 * t428;
t298 = qJD(2) * t376 + 0.2e1 * t475;
t331 = -qJD(4) * t407 + t509 * t658;
t447 = t452 * pkin(9);
t437 = -t441 * t276 - t184 * t252 / 0.2e1 + t185 * t272 / 0.2e1 + t277 * t638;
t439 = t203 * t637 + t333 * t490 + t387 * t643;
t3 = t437 + t439;
t65 = t223 * t252 + t238 * t272 + t342 * t343;
t446 = t3 * qJD(1) + t65 * qJD(2) + t52 * qJD(3);
t66 = pkin(5) * t517 - t238 * t657;
t7 = t483 * t185 + (t276 * t489 + t643) * pkin(5);
t445 = -qJD(1) * t7 + qJD(2) * t66 + qJD(3) * t64;
t444 = t588 / 0.2e1 + t479 / 0.2e1;
t102 = t416 ^ 2 * t568 - t273 * t301 + t274 * t302;
t5 = pkin(4) * t309 + t571 * t639 + t447 + (t274 * t632 + t273 * t432 / 0.2e1 + t485) * t276 + t454;
t443 = -t5 * qJD(1) + t102 * qJD(2) + t70 * qJD(3);
t347 = t388 * t493;
t442 = t347 + (t238 / 0.2e1 + t583 / 0.2e1) * t435;
t340 = -t577 / 0.2e1 + t518 * t433;
t71 = t444 - t453;
t89 = t432 * t471 + t442 + t485;
t438 = -qJD(1) * t71 + qJD(2) * t89 - qJD(3) * t340 + qJD(4) * t292;
t422 = t429 * t436;
t421 = -t524 / 0.2e1;
t420 = t524 / 0.2e1;
t419 = t523 / 0.2e1;
t413 = t435 * t521;
t412 = t432 * t523;
t391 = t407 * qJD(5);
t374 = -t413 + t425;
t372 = t476 * t432;
t366 = (t521 - qJD(5) / 0.2e1) * t433;
t358 = t373 * pkin(5);
t357 = (t422 + t576) * qJD(4);
t355 = t412 - t504;
t354 = -t432 * t520 - t505;
t353 = -t498 - t506;
t352 = -t435 * t520 + t507;
t341 = t574 / 0.2e1 + (t636 + 0.1e1 / 0.2e1) * t433;
t337 = t373 * t436;
t336 = t371 * t436;
t335 = t432 * t460;
t330 = t429 * t500 - t474;
t329 = t427 * t500 + t474;
t321 = t412 - t527;
t320 = -t504 + t527;
t319 = t506 - t528;
t318 = t498 + t528;
t300 = -t378 * qJD(4) + t433 * t506;
t299 = t377 * qJD(4) + t433 * t504;
t296 = t338 * t658;
t293 = t314 * qJD(4);
t283 = -t429 * t499 - t530;
t282 = -t427 * t499 + t530;
t281 = -t376 * qJD(5) - 0.2e1 * t436 * t475;
t280 = t428 * t594;
t251 = -t530 + (t429 * t524 + t508) * t436;
t250 = t530 + (t427 * t524 - t508) * t436;
t232 = (t422 - t576) * qJD(4) + 0.2e1 * (-qJD(5) - t521) * t514;
t229 = -pkin(5) * t572 + t388 * t492 + t347;
t225 = t382 + t365 + t449;
t224 = t364 - t450 + t513;
t206 = 0.2e1 * t309;
t121 = t424 * t493 + t586 / 0.2e1 + t416 * t488 + t364 + t432 * t486 - t472;
t113 = t581 / 0.2e1 + t577 * t646 - pkin(5) * t574 / 0.2e1 + t628 + qJ(6) * t487 - t448 + t519;
t88 = pkin(5) * t491 + t442 + t484 + t496;
t82 = t605 / 0.2e1 + t473 + t333 * t491 + t589 / 0.2e1;
t81 = -t607 / 0.2e1 + t276 * t493 + t435 * t495 - t590 / 0.2e1;
t73 = t308 + t464;
t72 = t444 + t453;
t62 = t185 * t633 + t277 * t488 + t333 * t493;
t61 = t184 * t633 + t277 * t492 + t333 * t488;
t54 = pkin(5) * t487 - t435 * t483;
t51 = qJD(2) * t80;
t50 = qJD(2) * t79;
t44 = (t184 * t632 + t185 * t634) * t436 + t452;
t28 = qJD(2) * t60;
t27 = qJD(2) * t59;
t26 = t587 * t646 + t657 * t637 + (t424 * t488 + t640) * pkin(5);
t24 = qJD(2) * t82 + qJD(4) * t152 - t543;
t23 = qJD(2) * t81 + qJD(4) * t154 + qJD(5) * t184;
t22 = qJD(2) * t62 + qJD(5) * t154 + t277 * t545;
t21 = qJD(2) * t61 + qJD(5) * t152 - t277 * t544;
t20 = pkin(5) * t152;
t18 = -qJD(4) * t59 - qJD(5) * t80;
t17 = -qJD(4) * t60 - qJD(5) * t79;
t14 = qJD(2) * t45;
t13 = qJD(4) * t45;
t12 = t44 * qJD(2) - qJD(4) * t478;
t11 = (t203 * t436 - t333 * t566) * qJD(2) + t62 * qJD(4) + t81 * qJD(5);
t10 = (-t202 * t436 - t333 * t575) * qJD(2) + t61 * qJD(4) + t82 * qJD(5);
t9 = t44 * qJD(4) + (-t202 * t435 - t203 * t432) * t524;
t8 = (t641 + t642) * t185 + (t473 + t643) * pkin(5);
t6 = -t274 * t154 / 0.2e1 - t273 * t152 / 0.2e1 + t276 * t484 + (t277 * t416 / 0.2e1 + t333 * pkin(4) / 0.2e1) * t433 + t447 - t454;
t4 = qJD(2) * t90 + qJD(4) * t33;
t2 = t437 - t439;
t1 = qJD(2) * t36 + qJD(4) * t41;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t83, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t480, -qJD(2) * t470, 0, 0, 0, 0, 0, 0, 0, 0, -t334 * qJD(2), t531, 0 (-t333 * t431 - t334 * t613) * qJD(2) * pkin(2), 0, 0, 0, 0, 0, 0, t206 * qJD(4) - t334 * t521, qJD(4) * t208 + t334 * t524 (-t428 - t430) * t531, t554 + (t334 * t417 - t430 * t594 - t280) * qJD(2), 0, 0, 0, 0, 0, 0, t10, t11, t9 (-t202 * t273 + t203 * t274 - t280) * qJD(2) + t6 * qJD(4) + t622, 0, 0, 0, 0, 0, 0, t10, t11, t9 (t202 * t223 + t203 * t238 - t342 * t592) * qJD(2) + t2 * qJD(4) + t8 * qJD(5) + t73 * qJD(6) + t622; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t206 - qJD(4) * t277, qJD(2) * t208 + qJD(4) * t276, 0, 0, 0, 0, 0, 0, 0, 0, t21, t22, t12, t6 * qJD(2) + (-pkin(4) * t277 - pkin(9) * t478) * qJD(4) + t624, 0, 0, 0, 0, 0, 0, t21, t22, t12, t2 * qJD(2) + (-t276 * t292 + t277 * t424) * qJD(4) + t20 * qJD(5) + t72 * qJD(6) + t624; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t23, 0, 0, 0, 0, 0, 0, 0, 0, t24, t23, 0, -pkin(5) * t543 + qJD(2) * t8 + qJD(4) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t73 + qJD(4) * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t205 * qJD(4), 0, 0, -t554, 0, 0, 0, 0, 0, 0, t17, t18, t13, -qJD(4) * t5 - t622, 0, 0, 0, 0, 0, 0, t17, t18, t13, qJD(4) * t3 - qJD(5) * t7 + qJD(6) * t74 - t622; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t500, t408 * qJD(4), 0, -t500, 0, 0, t417 * t523, t417 * t520, 0, 0, t330, t281, t300, t329, t299, -t500, -qJD(4) * t132 - qJD(5) * t210, qJD(4) * t133 + qJD(5) * t209, -qJD(4) * t91, qJD(4) * t102, t330, t281, t300, t329, t299, -t500, -t92 * qJD(4) - t137 * qJD(5) + t433 * t497, t105 * qJD(4) + t143 * qJD(5) - t433 * t501, -qJD(4) * t68 + qJD(5) * t67 + qJD(6) * t375, qJD(4) * t65 + qJD(5) * t66 - qJD(6) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t553 + t617, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t553 + t618; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t499, t525, t520, -t499, -t523, 0, -t416 * t520 - t451, t416 * t523 + t510, 0, 0, t251, t232, t321, t250, t318, -t366 (t432 * t469 - t511) * qJD(4) + t225 * qJD(5) + t655 (t435 * t469 + t516) * qJD(4) + t224 * qJD(5) + t656, qJD(4) * t463 - t648 (-pkin(4) * t563 + pkin(9) * t463) * qJD(4) + t443, t251, t232, t321, t250, t318, -t366 (t424 * t572 - t583 - t585) * qJD(4) + t113 * qJD(5) + t501 + t650 (t343 * t432 - t388 * t433 + t424 * t565) * qJD(4) + t121 * qJD(5) + t497 + t654 ((t272 + t582) * t435 + (-t252 - t581) * t432) * qJD(4) + t54 * qJD(5) + t647 (-t252 * t387 + t272 * t388 + t343 * t424) * qJD(4) + t26 * qJD(5) + t88 * qJD(6) + t446; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t284, -t298, t335, t284, t338, t419, qJD(4) * t225 - qJD(5) * t274 - t653, qJD(4) * t224 + qJD(5) * t273 - t649, 0, 0, -t284, -t298, t335, t284, t338, t419, qJD(4) * t113 - t542 - t651, qJD(4) * t121 + qJD(5) * t237 + t652, pkin(5) * t507 + qJD(4) * t54 + t609, -pkin(5) * t542 + qJD(4) * t26 + t445; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t337, -t336, t529, qJD(4) * t88 + t462; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t553 + t617, 0, 0, 0, 0, 0, 0, 0, 0, 0, t553 + t618; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t293, 0, 0, 0, 0, 0, 0, 0, 0, 0, t293; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t523, -t520, 0, 0, 0, 0, 0, 0, 0, 0, t353, t355, t357 (t406 * t626 - t629) * qJD(4) - t660, 0, 0, 0, 0, 0, 0, t353, t355, t357 (t292 * t436 + t570) * qJD(4) + t229 * qJD(5) + t341 * qJD(6) + t659; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t354, t352, 0, 0, 0, 0, 0, 0, 0, 0, t354, t352, 0, qJD(4) * t229 - t477 + t610; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t341 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t205 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t27, -t14, qJD(2) * t5 - t624, 0, 0, 0, 0, 0, 0, t28, t27, -t14, -qJD(2) * t3 - qJD(6) * t71 - t624; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t499, -t525, 0, t499, 0, 0, t451, -t510, 0, 0, t283, t296, t320, t282, t319, t366, qJD(5) * t279 - t655, qJD(5) * t278 - t656, t648, -t443, t283, t296, t320, t282, t319, t366, -qJD(5) * t112 - t650, -qJD(5) * t120 - t654, qJD(5) * t53 - t647, -qJD(5) * t25 + qJD(6) * t89 - t446; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t660, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(6) * t340 - t659; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t411, t391, 0, -t411, 0, 0, -pkin(4) * t540, -pkin(4) * t425, 0, 0, t411, t391, 0, -t411, 0, 0, -t349 * qJD(5), t362 * qJD(5), qJD(6) * t406, qJD(5) * t194 + qJD(6) * t292; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t305, -t331, t374, -t305, t372, t421, -pkin(9) * t425 - t455, pkin(9) * t540 - t456, 0, 0, t305, -t331, t374, -t305, t372, t421, t459 - t541, qJD(5) * t387 + t458, -pkin(5) * t425 + t611, -pkin(5) * t541 + t461; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t526, t438; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, t51, 0, 0, 0, 0, 0, 0, 0, 0, t50, t51, 0, qJD(2) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t284, t298, -t336, -t284, -t337, t419, -qJD(4) * t279 + t653, -qJD(4) * t278 + t649, 0, 0, t284, t298, -t336, -t284, -t337, t419, qJD(4) * t112 - t503 + t651, qJD(4) * t120 + t433 * t538 - t652, -qJD(4) * t53 - t609, -pkin(5) * t503 + qJD(4) * t25 - t445; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t610; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t305, t331, t413, t305, -t502, t420, t455, t456, 0, 0, -t305, t331, t413, t305, -t502, t420, -t459 - t538, -t458 - t537, -t611, -pkin(5) * t538 - t461; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t373, t371, 0, -t358; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t74 + qJD(4) * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t338, t335, -t529, -qJD(4) * t89 - t462 + t477; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t340 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t540, t425, -t526, pkin(5) * t540 - t438; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t373, -t371, 0, t358; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t15;
