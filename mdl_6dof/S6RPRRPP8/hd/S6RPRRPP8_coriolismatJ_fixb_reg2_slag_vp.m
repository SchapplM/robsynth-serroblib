% Calculate inertial parameters regressor of coriolis matrix for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRRPP8_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:55:59
% EndTime: 2019-03-09 04:56:16
% DurationCPUTime: 10.98s
% Computational Cost: add. (5474->722), mult. (10442->765), div. (0->0), fcn. (8543->4), ass. (0->485)
t394 = sin(qJ(4));
t390 = t394 ^ 2;
t396 = cos(qJ(4));
t391 = t396 ^ 2;
t534 = t390 + t391;
t393 = pkin(4) + qJ(6);
t553 = t396 * t393;
t561 = t394 * qJ(5);
t667 = -t553 - t561;
t666 = 0.2e1 * t394;
t397 = cos(qJ(3));
t632 = pkin(8) * t397;
t395 = sin(qJ(3));
t633 = pkin(3) * t395;
t445 = -t632 + t633;
t425 = qJ(2) + t445;
t285 = t396 * t425;
t398 = -pkin(1) - pkin(7);
t372 = t395 * t398;
t202 = t372 * t394 - t285;
t552 = t396 * t397;
t502 = pkin(5) * t552;
t137 = -t202 - t502;
t557 = t394 * t398;
t461 = -pkin(4) + t557;
t86 = t502 - t285 + (-qJ(6) + t461) * t395;
t665 = t86 + t137;
t508 = t395 * qJD(1);
t355 = qJD(4) + t508;
t300 = t355 * qJ(5);
t495 = t86 / 0.2e1 + t137 / 0.2e1;
t664 = t396 * t495;
t381 = qJD(5) * t396;
t629 = t396 * pkin(4);
t440 = -t561 - t629;
t408 = t440 * qJD(4);
t663 = t408 + t381;
t589 = t202 * t396;
t481 = -t589 / 0.2e1;
t171 = t395 * t461 - t285;
t595 = t171 * t396;
t662 = t481 + t595 / 0.2e1;
t505 = t397 * qJD(3);
t344 = t394 * t505;
t392 = t397 ^ 2;
t371 = t392 * t396;
t656 = t395 ^ 2;
t293 = t396 * t656 - t371;
t510 = t293 * qJD(1);
t661 = t510 - t344;
t380 = t395 * qJD(4);
t352 = t396 * t380;
t206 = t510 + t352;
t348 = t396 * t505;
t370 = t392 * t394;
t559 = t394 * t656;
t291 = -t370 + t559;
t512 = t291 * qJD(1);
t660 = t512 + t348;
t351 = t394 * t380;
t204 = t512 + t351;
t551 = t396 * t398;
t317 = t395 * t551;
t203 = t394 * t425 + t317;
t503 = t656 / 0.2e1;
t636 = -t396 / 0.2e1;
t640 = t371 / 0.2e1;
t250 = t396 * t503 + t636 + t640;
t515 = t250 * qJD(2);
t659 = qJD(4) * t203 + t515;
t462 = t503 + 0.1e1 / 0.2e1;
t249 = t396 * t462 + t640;
t517 = t249 * qJD(1);
t105 = t517 + t344 + t352;
t641 = -t370 / 0.2e1;
t248 = -t394 * t462 + t641;
t518 = t248 * qJD(1);
t104 = t518 + t348 - t351;
t548 = t397 * t398;
t318 = t394 * t548;
t628 = t397 * pkin(3);
t630 = t395 * pkin(8);
t312 = t628 + t630;
t554 = t396 * t312;
t211 = -t318 + t554;
t294 = t394 * t312;
t319 = t396 * t548;
t212 = t319 + t294;
t585 = t203 * t397;
t588 = t202 * t397;
t637 = t395 / 0.2e1;
t31 = -t395 * t548 + (t212 * t637 + t585 / 0.2e1) * t396 + (-t211 * t395 / 0.2e1 + t588 / 0.2e1) * t394;
t555 = t395 * t397;
t200 = (-0.1e1 + t534) * t555;
t520 = t200 * qJD(2);
t658 = -t31 * qJD(1) - t520;
t507 = t395 * qJD(3);
t653 = pkin(5) + pkin(8);
t311 = t653 * t394;
t570 = t311 * t395;
t558 = t394 * t397;
t357 = pkin(4) * t558;
t619 = qJ(5) * t396;
t437 = qJ(6) * t394 - t619;
t155 = t357 + (-t398 + t437) * t397;
t598 = t155 * t396;
t657 = -t570 / 0.2e1 - t598 / 0.2e1;
t333 = t391 - t390;
t87 = t552 * t666 * (-qJD(4) + t508) - t333 * t507;
t428 = t355 * t397;
t241 = t396 * t428;
t187 = t241 * t666;
t655 = -pkin(4) / 0.2e1;
t139 = -pkin(5) * t558 + t203;
t382 = t395 * qJ(5);
t113 = t139 + t382;
t88 = t113 * t394;
t654 = -t88 / 0.2e1;
t652 = t139 / 0.2e1;
t651 = -t171 / 0.2e1;
t383 = t397 * qJ(5);
t264 = pkin(4) * t552 + t394 * t383;
t550 = t397 * qJ(6);
t210 = t396 * t550 + t264;
t650 = -t210 / 0.2e1;
t218 = t357 + (-t398 - t619) * t397;
t649 = -t218 / 0.2e1;
t263 = -pkin(3) + t667;
t648 = -t263 / 0.2e1;
t647 = -t264 / 0.2e1;
t384 = pkin(4) * t394;
t309 = t384 - t619;
t646 = -t309 / 0.2e1;
t645 = -t311 / 0.2e1;
t644 = -t312 / 0.2e1;
t313 = t653 * t396;
t643 = t313 / 0.2e1;
t642 = -t357 / 0.2e1;
t639 = -t393 / 0.2e1;
t638 = t394 / 0.2e1;
t635 = t396 / 0.2e1;
t634 = t397 / 0.2e1;
t631 = t395 * pkin(5);
t627 = t397 * pkin(4);
t560 = t394 * t395;
t149 = t203 * t560;
t474 = -t560 / 0.2e1;
t62 = t149 / 0.2e1 + t203 * t474;
t626 = t31 * qJD(3) - t62 * qJD(4);
t599 = t155 * t395;
t556 = t395 * t396;
t537 = qJ(5) * t556 + t372;
t154 = -t393 * t560 + t537;
t600 = t154 * t397;
t417 = t600 / 0.2e1 - t599 / 0.2e1;
t182 = -t383 - t212;
t118 = pkin(5) * t560 - t182;
t604 = t118 * t395;
t605 = t113 * t397;
t93 = t318 - t393 * t397 + (-t312 - t631) * t396;
t620 = t93 * t395;
t621 = t86 * t397;
t3 = (t645 - t604 / 0.2e1 - t605 / 0.2e1) * t396 + (t643 - t620 / 0.2e1 - t621 / 0.2e1) * t394 + t417;
t625 = t3 * qJD(1);
t7 = t113 * t118 + t154 * t155 + t86 * t93;
t624 = t7 * qJD(1);
t562 = t393 * t395;
t8 = t654 + (t652 - t382 / 0.2e1) * t394 + (-t562 / 0.2e1 + t495) * t396;
t623 = t8 * qJD(1);
t622 = t86 * t395;
t583 = t210 * t394;
t45 = -t139 * t395 + (t583 + t598) * t397;
t618 = qJD(1) * t45;
t114 = t155 * t558;
t46 = -t137 * t395 + t210 * t552 - t114;
t617 = qJD(1) * t46;
t48 = t86 * t396 - t88;
t616 = qJD(1) * t48;
t51 = t114 - t622;
t615 = qJD(1) * t51;
t586 = t203 * t395;
t52 = -t586 + (t218 * t396 + t264 * t394) * t397;
t614 = qJD(1) * t52;
t590 = t202 * t395;
t53 = -t218 * t558 + t264 * t552 + t590;
t613 = qJD(1) * t53;
t606 = t113 * t395;
t54 = t155 * t552 - t606;
t612 = qJD(1) * t54;
t167 = -t382 - t203;
t597 = t167 * t394;
t60 = t595 + t597;
t611 = qJD(1) * t60;
t596 = t167 * t395;
t65 = t218 * t552 + t596;
t610 = qJD(1) * t65;
t587 = t203 * t394;
t69 = t587 - t589;
t609 = qJD(1) * t69;
t608 = qJD(2) * t62;
t470 = t553 / 0.2e1;
t475 = t561 / 0.2e1;
t412 = t475 + t470;
t465 = -t139 / 0.2e1 + t113 / 0.2e1;
t582 = t210 * t397;
t10 = t582 / 0.2e1 + (t394 * t465 - t664) * t395 + t412;
t607 = t10 * qJD(1);
t13 = t113 * t137 + t139 * t86 + t155 * t210;
t603 = t13 * qJD(1);
t217 = -pkin(4) * t560 + t537;
t188 = -t211 - t627;
t592 = t188 * t394;
t593 = t182 * t396;
t14 = (t217 / 0.2e1 + t167 * t635 + t394 * t651) * t397 + (t649 + t593 / 0.2e1 - t592 / 0.2e1) * t395;
t602 = t14 * qJD(1);
t15 = (t397 * t93 - t622) * t396 + (-t118 * t397 + t606) * t394;
t601 = t15 * qJD(1);
t18 = ((-t113 + t139) * t396 - t665 * t394) * t397;
t594 = t18 * qJD(1);
t464 = t202 / 0.2e1 + t651;
t483 = -t597 / 0.2e1;
t404 = t396 * t464 + t483;
t421 = t475 + t629 / 0.2e1;
t444 = -t149 / 0.2e1 + t264 * t634;
t19 = t395 * t404 + t421 + t444;
t591 = t19 * qJD(1);
t21 = -t620 - t621 + (-t599 + t600) * t394;
t584 = t21 * qJD(1);
t581 = t217 * t396;
t580 = t218 * t394;
t25 = t154 * t552 - t155 * t556 - t604 - t605;
t579 = t25 * qJD(1);
t26 = t167 * t202 + t171 * t203 + t218 * t264;
t578 = t26 * qJD(1);
t577 = t263 * t396;
t473 = t560 / 0.2e1;
t314 = qJ(5) * t473;
t480 = -t587 / 0.2e1;
t414 = t483 + t480;
t500 = pkin(4) * t637;
t27 = t314 + (t500 + t464) * t396 + t414;
t576 = t27 * qJD(1);
t274 = t384 + t437;
t575 = t274 * t394;
t29 = (-t171 * t395 + t188 * t397) * t396 + (t182 * t397 - t596) * t394;
t574 = t29 * qJD(1);
t297 = -pkin(3) + t440;
t573 = t297 * t395;
t572 = t309 * t394;
t569 = t313 * t395;
t32 = ((t167 + t203) * t396 + (-t171 + t202) * t394) * t397;
t568 = t32 * qJD(1);
t33 = t182 * t395 - t218 * t556 + (t167 + t581) * t397;
t567 = t33 * qJD(1);
t34 = t202 * t556 - t149 + (t211 * t396 + t212 * t394) * t397;
t566 = t34 * qJD(1);
t35 = -t171 * t397 - t188 * t395 + t217 * t558 - t218 * t560;
t565 = t35 * qJD(1);
t564 = t390 * t392;
t563 = t393 * t394;
t549 = t397 * t274;
t58 = -t585 + (-t212 + 0.2e1 * t319) * t395;
t547 = t58 * qJD(1);
t59 = -t588 + (t211 + 0.2e1 * t318) * t395;
t546 = t59 * qJD(1);
t545 = t62 * qJD(1);
t522 = qJD(5) * t394;
t349 = t397 * t522;
t544 = t520 + t349;
t247 = -t559 / 0.2e1 + t641 + t638;
t224 = t247 * qJD(2);
t379 = t395 * qJD(5);
t543 = t224 + t379;
t276 = t294 / 0.2e1;
t542 = t276 + t319 / 0.2e1;
t280 = t554 / 0.2e1;
t541 = t280 - t318 / 0.2e1;
t369 = t391 * t392;
t302 = t369 + t656;
t531 = qJD(2) * t395;
t340 = t394 * t531;
t540 = t302 * qJD(5) + t340;
t536 = t534 * t632;
t226 = t248 * qJD(2);
t535 = t379 - t226;
t102 = t392 * t557 + t590;
t533 = qJD(1) * t102;
t103 = -t392 * t551 - t586;
t532 = qJD(1) * t103;
t530 = qJD(3) * t394;
t529 = qJD(3) * t396;
t528 = qJD(4) * t202;
t526 = qJD(4) * t313;
t525 = qJD(4) * t394;
t524 = qJD(4) * t396;
t463 = -t390 / 0.2e1 + t391 / 0.2e1;
t245 = (0.1e1 / 0.2e1 + t463) * t397;
t523 = qJD(5) * t245;
t521 = qJD(6) * t396;
t246 = (-0.1e1 / 0.2e1 + t463) * t397;
t519 = t246 * qJD(5);
t516 = t249 * qJD(2);
t275 = t463 * t397;
t514 = t275 * qJD(4);
t292 = t534 * t397;
t511 = t292 * qJD(1);
t332 = -t392 + t656;
t509 = t332 * qJD(1);
t310 = t333 * qJD(4);
t378 = t395 * qJD(6);
t506 = t397 * qJD(1);
t504 = t397 * qJD(4);
t501 = pkin(8) * t525;
t499 = -t631 / 0.2e1;
t498 = -t630 / 0.2e1;
t497 = -t627 / 0.2e1;
t496 = t655 - qJ(6) / 0.2e1;
t469 = -t552 / 0.2e1;
t471 = t556 / 0.2e1;
t494 = -t580 / 0.2e1 + t297 * t469 + pkin(8) * t471;
t493 = t383 + t542;
t492 = qJ(2) * t508;
t491 = qJ(2) * t506;
t490 = t396 * t506;
t489 = t394 * t504;
t488 = t396 * t504;
t487 = t394 * t381;
t486 = t394 * t521;
t343 = t394 * t524;
t342 = t394 * t529;
t350 = t395 * t505;
t485 = t395 * t506;
t484 = t619 / 0.2e1;
t479 = t263 * t634;
t478 = t297 * t634;
t477 = -t572 / 0.2e1;
t476 = -t569 / 0.2e1;
t472 = t558 / 0.2e1;
t468 = t552 / 0.2e1;
t467 = -t549 / 0.2e1;
t466 = t397 * t646;
t277 = -t294 / 0.2e1;
t460 = t277 - t319 / 0.2e1 - t383;
t162 = qJD(1) * t245 + t342;
t459 = qJD(1) * t246 + t342;
t198 = qJD(1) * t275 + t342;
t290 = -t369 + t564;
t305 = t394 * t348;
t190 = qJD(1) * t290 - 0.2e1 * t305;
t303 = t394 * qJD(1) * t371;
t131 = qJD(3) * t245 - t303;
t458 = qJD(3) * t246 - t303;
t175 = qJD(3) * t275 - t303;
t341 = t394 * t508;
t158 = qJD(4) * t247 + t341;
t159 = qJD(4) * t248 - t341;
t347 = t396 * t508;
t457 = qJD(4) * t249 + t347;
t456 = qJD(4) * t250 - t347;
t306 = t318 / 0.2e1;
t454 = t306 + t494;
t453 = t394 * t490;
t452 = t392 * t343;
t451 = t397 * t486;
t450 = t496 * t397;
t449 = t644 + t499;
t448 = t639 + t496;
t320 = pkin(5) * t473;
t446 = t320 + t383 / 0.2e1 + t542;
t443 = -qJ(5) * t380 - t516;
t301 = t564 + t656;
t442 = qJD(1) * t301 - t305 + t380;
t156 = qJD(1) * t302 + t305 + t380;
t441 = pkin(3) * t472 + pkin(8) * t473;
t439 = t573 + t632;
t438 = t263 * t472 + t657;
t22 = t167 * t182 + t171 * t188 + t217 * t218;
t436 = t22 * qJD(1) - t14 * qJD(2);
t44 = -t398 ^ 2 * t555 - t202 * t211 + t203 * t212;
t435 = t44 * qJD(1) + t31 * qJD(2);
t415 = t155 * t638 + t476;
t407 = t306 + t415;
t16 = (t650 + t449) * t396 + (t575 / 0.2e1 + t577 / 0.2e1 + t448) * t397 + t407;
t99 = -t263 * t394 + t274 * t396;
t434 = -qJD(1) * t16 + qJD(3) * t99;
t402 = t396 * t467 + t438;
t24 = (t650 + t499) * t394 + t402 + t460;
t98 = t575 + t577;
t433 = -qJD(1) * t24 + qJD(3) * t98;
t432 = t592 - t593;
t431 = -t211 * t394 + t212 * t396;
t168 = t297 * t396 + t572;
t399 = pkin(8) * t474 + (t478 + t647) * t394 + (t466 + t649) * t396;
t37 = t399 + t460;
t430 = -qJD(1) * t37 + qJD(3) * t168;
t169 = -t297 * t394 + t309 * t396;
t38 = t264 * t635 + (t477 + pkin(4)) * t397 + t494 + t541;
t429 = -qJD(1) * t38 - qJD(3) * t169;
t242 = t448 * t395;
t427 = -qJD(1) * t242 + qJD(4) * t393;
t426 = -t521 - t522;
t424 = -t628 / 0.2e1 + t498;
t346 = t396 * t531;
t423 = t392 * t487 - t346;
t422 = -t182 * qJ(5) / 0.2e1 + t188 * t655;
t420 = t118 * qJ(5) / 0.2e1 + t93 * t639;
t143 = t276 + t441;
t419 = pkin(3) * t529 - qJD(1) * t143;
t144 = (t644 + t424) * t396;
t418 = pkin(3) * t530 - qJD(1) * t144;
t416 = -t155 * t274 / 0.2e1 + t210 * t648;
t413 = -t563 / 0.2e1 + t484;
t40 = t450 + (t479 + t449) * t396 + t407;
t411 = qJD(1) * t40 + t263 * t530;
t42 = t558 * t648 + t446 - t657;
t410 = qJD(1) * t42 + t263 * t529;
t49 = t306 + t497 + t580 / 0.2e1 + (t644 + t478 + t498) * t396;
t409 = qJD(1) * t49 + t297 * t530;
t219 = -qJD(3) * t333 + 0.2e1 * t453;
t192 = -qJD(3) * t291 + t395 * t488;
t406 = qJD(3) * t293 + t395 * t489;
t1 = t311 * t465 - t313 * t495 + t416 + t420;
t55 = t263 * t274;
t56 = (t274 / 0.2e1 + t413) * t397;
t405 = t1 * qJD(1) + t56 * qJD(2) - t55 * qJD(3);
t153 = qJD(4) * t290 + 0.2e1 * t305 * t395;
t138 = t642 + (t484 + t309 / 0.2e1) * t397;
t400 = pkin(8) * t480 + t218 * t646 + t297 * t647;
t5 = pkin(8) * t404 + t400 + t422;
t403 = t297 * t309 * qJD(3) - t5 * qJD(1) - t138 * qJD(2);
t401 = qJD(4) * t667 - qJD(6) * t394;
t389 = qJ(2) * qJD(2);
t388 = qJ(5) * qJD(5);
t387 = qJD(1) * qJ(2);
t377 = t390 * qJD(5);
t376 = pkin(8) * t524;
t368 = qJ(5) * t379;
t363 = -t506 / 0.2e1;
t362 = t506 / 0.2e1;
t361 = t505 / 0.2e1;
t345 = t396 * t379;
t304 = t396 * t349;
t296 = t311 * qJD(4);
t288 = t347 + t524;
t287 = t341 + t525;
t286 = (t508 + qJD(4) / 0.2e1) * t397;
t269 = t292 * qJD(2);
t267 = t292 * qJD(3);
t262 = -t394 * t507 + t488;
t261 = t396 * t507 + t489;
t260 = -qJD(3) * t391 + t453;
t259 = qJD(3) * t390 + t453;
t240 = (t490 + t530) * t395;
t239 = t394 * t428;
t238 = (t394 * t506 - t529) * t395;
t231 = t250 * qJD(5);
t228 = t249 * qJD(5);
t216 = -t350 * t391 - t452;
t215 = -t350 * t390 + t452;
t184 = t200 * qJD(3);
t173 = t391 * t485 + t514;
t172 = t390 * t485 - t514;
t140 = qJ(5) * t468 + t466 + t642;
t97 = t514 + (-t391 * t506 - t342) * t395;
t96 = -t514 + (-t390 * t506 + t342) * t395;
t92 = t396 * t424 + t280 - t318;
t91 = -t319 + t277 + t441;
t71 = t562 / 0.2e1 + qJ(6) * t637 + t500 + t137;
t57 = t397 * t413 + t467;
t50 = -t554 / 0.2e1 + t497 + t454;
t43 = t438 + t446;
t41 = t263 * t469 + t396 * t449 + t306 - t415 + t450;
t39 = (t477 - pkin(4)) * t397 + (t264 / 0.2e1 + t644) * t396 + t454;
t36 = t399 + t493;
t28 = pkin(4) * t471 + t314 - t414 + t662;
t23 = -t583 / 0.2e1 + t320 + t402 + t493;
t20 = t167 * t473 + t171 * t471 + t395 * t481 + t421 - t444;
t17 = t476 + t393 * t634 + t627 / 0.2e1 + t550 / 0.2e1 + pkin(5) * t471 + (t479 + t650) * t396 + (t549 / 0.2e1 + t155 / 0.2e1) * t394 + t541;
t12 = t14 * qJD(3);
t11 = -t582 / 0.2e1 + t139 * t473 + t113 * t474 + t412 + t665 * t471;
t9 = t139 * t638 + t395 * t470 + t314 + t654 + t664;
t6 = -t400 + t422 + (t597 / 0.2e1 + t662) * pkin(8);
t4 = t113 * t468 + t118 * t471 + t311 * t636 + t313 * t638 + t472 * t86 + t473 * t93 - t417;
t2 = t113 * t645 + t311 * t652 + t643 * t665 - t416 + t420;
t30 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t389, -t350, t332 * qJD(3), 0, t350, 0, 0, qJ(2) * t505 + t531, -qJ(2) * t507 + qJD(2) * t397, 0, t389, t216, t153, -t406, t215, -t192, t350, qJD(3) * t59 + qJD(4) * t103 + t346, qJD(3) * t58 + qJD(4) * t102 - t340, -qJD(3) * t34 - t269, qJD(2) * t69 + qJD(3) * t44, t350, t406, t192, t216, t153, t215, qJD(3) * t29 + qJD(4) * t32 - t349 * t395 - t269, -qJD(3) * t35 - qJD(4) * t52 + t423, -qJD(3) * t33 - qJD(4) * t53 + t540, -qJD(2) * t60 + qJD(3) * t22 + qJD(4) * t26 - qJD(5) * t65, t350, t192, -t406, t215, -t153, t216, t15 * qJD(3) + t18 * qJD(4) + t426 * t555 - t269, -qJD(3) * t25 - qJD(4) * t46 - t392 * t486 + t540, qJD(3) * t21 + qJD(4) * t45 + qJD(6) * t301 - t423, -qJD(2) * t48 + qJD(3) * t7 + qJD(4) * t13 - qJD(5) * t54 + qJD(6) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), t387, 0, 0, 0, 0, 0, 0, t508, t506, 0, t387, 0, 0, 0, 0, 0, 0, -t456, -t158, -t511, t609 + t626, 0, 0, 0, 0, 0, 0, -t511, t456, t158, qJD(4) * t20 - t12 + t231 - t611, 0, 0, 0, 0, 0, 0, -t511, t158, -t456, qJD(3) * t4 + qJD(4) * t11 + qJD(6) * t247 + t231 - t616; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t485, t509, -t507, t485, -t505, 0, -t398 * t507 + t491, -t398 * t505 - t492, 0, 0, t97, t87, -t661, t96, t660, t286, t546 + (t394 * t445 - t317) * qJD(3) + t92 * qJD(4), t547 + (-pkin(8) * t552 + (pkin(3) * t396 + t557) * t395) * qJD(3) + t91 * qJD(4), qJD(3) * t431 - t566 (-pkin(3) * t372 + pkin(8) * t431) * qJD(3) + t435, t286, t661, -t660, t97, t87, t96, qJD(3) * t432 + t28 * qJD(4) + t574, -t565 + (t394 * t439 + t581) * qJD(3) + t39 * qJD(4) - t519, -t567 + (-t217 * t394 + t396 * t439) * qJD(3) + t36 * qJD(4) + t304 (pkin(8) * t432 + t217 * t297) * qJD(3) + t6 * qJD(4) + t50 * qJD(5) + t436, t286, -t660, -t661, t96, -t87, t97, t601 + ((t118 - t570) * t396 + (t93 + t569) * t394) * qJD(3) + t9 * qJD(4), -t579 + (-t154 * t394 + t263 * t556 + t313 * t397) * qJD(3) + t23 * qJD(4) + t304 + t245 * qJD(6), t584 + (-t154 * t396 - t263 * t560 - t311 * t397) * qJD(3) + t17 * qJD(4) + t519 - t451, t624 + t4 * qJD(2) + (t118 * t313 + t154 * t263 + t311 * t93) * qJD(3) + t2 * qJD(4) + t41 * qJD(5) + t43 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175, t190, -t239, -t175, -t241, t361, qJD(3) * t92 + t532 - t659, qJD(3) * t91 - t224 + t528 + t533, 0, -t608, t361, t239, t241, t175, t190, -t175, t568 + t28 * qJD(3) + (-t383 * t396 + t357) * qJD(4) - t349, qJD(3) * t39 - t614 + t659, qJD(3) * t36 - t528 + t543 - t613, t578 + t20 * qJD(2) + t6 * qJD(3) + (-pkin(4) * t203 - qJ(5) * t202) * qJD(4) - t167 * qJD(5), t361, t241, -t239, -t175, -t190, t175, t594 + t9 * qJD(3) - t349 + ((t563 - t619) * qJD(4) - t521) * t397, qJD(3) * t23 + qJD(4) * t137 + t543 - t617, qJD(3) * t17 - qJD(4) * t139 + t378 - t515 + t618, t603 + t11 * qJD(2) + t2 * qJD(3) + (qJ(5) * t137 - t139 * t393) * qJD(4) + t113 * qJD(5) + t71 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t239, -t458, t156, qJD(3) * t50 - qJD(4) * t167 + t515 - t610, 0, 0, 0, 0, 0, 0, -t239, t156, t458, qJD(3) * t41 + qJD(4) * t113 + t515 - t612; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t241, t131, t442, qJD(3) * t43 + qJD(4) * t71 + t224 + t615; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), -t387, 0, 0, 0, 0, 0, 0, -t508, -t506, 0, -t387, 0, 0, 0, 0, 0, 0, -t457, -t159, t511, -t609 + t626, 0, 0, 0, 0, 0, 0, t511, t457, t159, -qJD(4) * t19 - t12 + t228 + t611, 0, 0, 0, 0, 0, 0, t511, t159, -t457, -qJD(3) * t3 - qJD(4) * t10 + qJD(6) * t248 + t228 + t616; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t507, -t505, 0, 0, 0, 0, 0, 0, 0, 0, -t261, -t262, t267 (t536 - t633) * qJD(3) - t658, 0, 0, 0, 0, 0, 0, t267, t261, t262, -t602 + (t536 + t573) * qJD(3) + t140 * qJD(4) + t544, 0, 0, 0, 0, 0, 0, t267, t262, -t261, t263 * t507 - t625 + t57 * qJD(4) + ((t311 * t394 + t313 * t396) * qJD(3) + t521) * t397 + t544; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, -t104, 0, -t545, 0, 0, 0, 0, 0, 0, 0, t105, t104, t140 * qJD(3) + t395 * t408 + t345 - t591, 0, 0, 0, 0, 0, 0, 0, t104, -t105, t57 * qJD(3) + t395 * t401 + t345 - t607; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t485, -t509, 0, -t485, 0, 0, -t491, t492, 0, 0, t173, -t187, t206, t172, -t204, -t286, qJD(4) * t144 - t546, qJD(4) * t143 - t547, t566, -t435, -t286, -t206, t204, t173, -t187, t172, -qJD(4) * t27 + t345 - t574, qJD(4) * t38 - t523 + t565, qJD(4) * t37 + t304 + t567, -qJD(4) * t5 - qJD(5) * t49 - t436, -t286, t204, t206, t172, t187, t173, qJD(4) * t8 - t378 * t394 + t345 - t601, qJD(4) * t24 + qJD(6) * t246 + t304 + t579, qJD(4) * t16 - t451 + t523 - t584, qJD(2) * t3 - qJD(4) * t1 - qJD(5) * t40 - qJD(6) * t42 - t624; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t658, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t138 - t520 + t602, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t56 - t520 + t625; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t343, t310, 0, -t343, 0, 0, -pkin(3) * t525, -pkin(3) * t524, 0, 0, 0, 0, 0, t343, t310, -t343, 0, qJD(4) * t169 - t487, -qJD(4) * t168 + t377 (qJD(4) * t309 - t522) * t297, 0, 0, 0, -t343, -t310, t343, 0, -qJD(4) * t98 + t377 + t486, -qJD(4) * t99 + qJD(6) * t391 + t487, t55 * qJD(4) + t263 * t426; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, -t219, t288, -t198, -t287, t363, -t376 - t418, -t419 + t501, 0, 0, t363, -t288, t287, t198, -t219, -t198, -t576 + t663, t376 - t429, -t430 - t501, pkin(8) * t663 + t403, t363, t287, t288, -t198, t219, t198, t381 + t401 + t623, -t296 - t433, -t434 - t526 (-qJ(5) * t311 - t313 * t393) * qJD(4) + t313 * qJD(5) - t311 * qJD(6) - t405; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t288, -t162, t259, t376 - t409, 0, 0, 0, 0, 0, 0, t288, t259, t162, -t411 + t526; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t287, t459, -t260, -t296 - t410; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t175, -t190, t238, t175, t240, t361, -qJD(3) * t144 + t516 - t532, -qJD(3) * t143 + t226 - t533, 0, t608, t361, -t238, -t240, -t175, -t190, t175, qJD(3) * t27 - t568, -qJD(3) * t38 - t516 + t614, -qJD(3) * t37 + t535 + t613, qJD(2) * t19 + qJD(3) * t5 + t368 - t578, t361, -t240, t238, t175, t190, -t175, -qJD(3) * t8 - t594, -qJD(3) * t24 + t535 + t617, -qJD(3) * t16 + t378 + t516 - t618, qJD(2) * t10 + qJD(3) * t1 - qJD(6) * t242 + t368 - t603; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t517, t518, 0, t545, 0, 0, 0, 0, 0, 0, 0, -t517, -t518, qJD(3) * t138 + t591, 0, 0, 0, 0, 0, 0, 0, -t518, t517, qJD(3) * t56 + t607; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t198, t219, -t347, t198, t341, t362, t418, t419, 0, 0, t362, t347, -t341, -t198, t219, t198, t576, t429, t430, -t403, t362, -t341, -t347, t198, -t219, -t198, -t623, t433, t434, t405; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t388, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJD(6), qJD(6) * t393 + t388; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t355, t300, 0, 0, 0, 0, 0, 0, 0, t355, 0, t300; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t355, t427; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t238, t131, -t156, qJD(3) * t49 + t443 + t610, 0, 0, 0, 0, 0, 0, t238, -t156, -t131, qJD(3) * t40 - t378 + t443 + t612; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t517, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t517; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t347, t162, -t259, t409, 0, 0, 0, 0, 0, 0, -t347, -t259, -t162, t411; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t355, -t300, 0, 0, 0, 0, 0, 0, 0, -t355, 0, -qJD(6) - t300; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t355; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t240, -t458, -t442, qJD(3) * t42 + qJD(4) * t242 + t535 - t615; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t518; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t341, -t459, t260, t410; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t355, qJD(5) - t427; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t355; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t30;
