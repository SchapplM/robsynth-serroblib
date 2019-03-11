% Calculate inertial parameters regressor of coriolis matrix for
% S6RRPPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRPPRP2_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:32:05
% EndTime: 2019-03-09 08:32:20
% DurationCPUTime: 11.26s
% Computational Cost: add. (9719->519), mult. (18058->609), div. (0->0), fcn. (19880->6), ass. (0->407)
t374 = sin(qJ(2));
t595 = -qJ(3) - pkin(7);
t344 = t595 * t374;
t576 = cos(pkin(9));
t339 = t576 * t344;
t372 = sin(pkin(9));
t376 = cos(qJ(2));
t529 = t372 * t376;
t415 = t595 * t529;
t253 = t339 + t415;
t337 = t374 * t576 + t529;
t598 = t337 * pkin(4);
t189 = t253 - t598;
t373 = sin(qJ(5));
t550 = t189 * t373;
t375 = cos(qJ(5));
t636 = t189 * t375;
t358 = t576 * t376;
t530 = t372 * t374;
t335 = -t358 + t530;
t601 = t335 * pkin(4);
t531 = t372 * t344;
t624 = t595 * t358;
t630 = -t624 + t531;
t632 = t630 - t601;
t633 = t375 * t632;
t437 = t633 / 0.2e1;
t634 = t373 * t632;
t635 = -t634 / 0.2e1;
t360 = -t376 * pkin(2) - pkin(1);
t400 = -t337 * qJ(4) + t360;
t608 = pkin(3) + pkin(8);
t154 = t335 * t608 + t400;
t86 = t636 + (qJ(6) * t335 + t154) * t373;
t631 = -t86 / 0.2e1;
t389 = -t624 / 0.2e1;
t222 = t373 * t337;
t482 = t222 * qJD(1);
t493 = qJD(5) * t373;
t179 = t482 + t493;
t610 = t337 ^ 2;
t611 = t335 ^ 2;
t617 = t611 - t610;
t629 = t617 * qJD(1);
t628 = t617 * qJD(2);
t618 = t610 + t611;
t627 = t618 * qJD(1);
t626 = qJD(3) * t618;
t224 = t375 * t335;
t453 = t373 * t224;
t424 = 0.2e1 * t453;
t597 = t337 * pkin(5);
t72 = t597 - t86;
t625 = t72 + t86;
t101 = t154 * t373 + t636;
t102 = t375 * t154 - t550;
t39 = (t101 * t373 + t102 * t375) * t337;
t586 = t72 * t373;
t87 = qJ(6) * t224 + t102;
t412 = -t87 * t375 + t586;
t31 = t412 * t337;
t473 = t335 * qJD(1);
t298 = t373 * t473;
t365 = t375 * qJD(2);
t623 = t365 + t298;
t364 = t373 * qJD(2);
t622 = -t375 * t473 + t364;
t250 = t624 / 0.2e1 + t389;
t357 = pkin(2) * t372 + qJ(4);
t465 = t357 * qJD(2);
t620 = qJD(1) * t250 + t465;
t381 = t531 / 0.2e1 + t389;
t377 = -t601 / 0.2e1 + t381;
t369 = t374 * pkin(2);
t537 = t335 * qJ(4);
t423 = t369 + t537;
t159 = t337 * t608 + t423;
t151 = t375 * t159;
t106 = t151 + t634;
t561 = t106 * t373;
t527 = t373 * t159;
t105 = t633 - t527;
t562 = t105 * t375;
t396 = -t562 / 0.2e1 - t561 / 0.2e1;
t37 = t377 + t396;
t619 = qJD(1) * t37 + t465;
t407 = -t253 * t337 - t335 * t630;
t616 = qJD(1) * t407;
t614 = qJD(3) * t407;
t518 = -t151 / 0.2e1 + t635;
t325 = t530 / 0.2e1 - t358 / 0.2e1;
t612 = t373 * t631 - t586 / 0.2e1;
t600 = t335 * pkin(5);
t73 = -t600 + t633 + (-t337 * qJ(6) - t159) * t373;
t609 = t73 / 0.2e1;
t606 = -t335 / 0.2e1;
t605 = -t337 / 0.2e1;
t370 = t373 ^ 2;
t604 = -t370 / 0.2e1;
t371 = t375 ^ 2;
t603 = t371 / 0.2e1;
t602 = pkin(5) * t373;
t599 = t337 * pkin(3);
t596 = t375 * pkin(5);
t581 = t87 * t373;
t585 = t72 * t375;
t594 = -t585 / 0.2e1 - t581 / 0.2e1;
t450 = -pkin(4) - t596;
t143 = t335 * t450 + t630;
t535 = t335 * t373;
t7 = pkin(5) * t143 * t535 - t625 * t87;
t592 = qJD(1) * t7;
t8 = t625 * t224;
t591 = qJD(1) * t8;
t590 = qJD(2) * pkin(2);
t142 = t337 * t450 + t253;
t526 = t375 * qJ(6);
t88 = t337 * t526 + t106;
t3 = t142 * t143 + t72 * t73 + t87 * t88;
t589 = t3 * qJD(1);
t431 = -t222 / 0.2e1;
t421 = pkin(5) * t431;
t4 = t421 + t612;
t588 = t4 * qJD(1);
t578 = t88 * t375;
t584 = t73 * t373;
t6 = -t31 + (t578 - t584) * t335;
t587 = t6 * qJD(1);
t583 = t73 * t375;
t579 = t88 * t373;
t456 = t597 / 0.2e1;
t420 = t373 * t456;
t9 = t420 - t612;
t577 = t9 * qJD(1);
t555 = t143 * t335;
t22 = -t555 + (t581 + t585) * t337;
t575 = qJD(1) * t22;
t341 = t357 + t602;
t435 = t341 * t606;
t359 = -pkin(2) * t576 - pkin(3);
t356 = -pkin(8) + t359;
t315 = (-qJ(6) + t356) * t375;
t540 = t315 * t375;
t528 = t373 * qJ(6);
t314 = t356 * t373 - t528;
t544 = t314 * t373;
t379 = (t544 / 0.2e1 + t540 / 0.2e1) * t337 + t435;
t398 = t584 / 0.2e1 - t578 / 0.2e1;
t24 = t379 + t398;
t574 = qJD(1) * t24;
t29 = -t632 * t335 + (-t101 * t375 + t102 * t373) * t337;
t573 = qJD(1) * t29;
t30 = t412 * t335;
t572 = qJD(1) * t30;
t571 = qJD(1) * t31;
t569 = qJD(1) * t39;
t538 = t611 * t375;
t40 = -t87 * t337 + (-pkin(5) * t538 + t555) * t373;
t568 = qJD(1) * t40;
t56 = t101 * t337 + t224 * t632;
t567 = qJD(1) * t56;
t57 = -t102 * t337 + t535 * t632;
t566 = qJD(1) * t57;
t565 = qJD(5) * t87;
t454 = t86 / 0.2e1 + t72 / 0.2e1;
t10 = (t456 + t454) * t375;
t564 = t10 * qJD(1);
t563 = t105 * t373;
t560 = t106 * t375;
t13 = -t101 * t105 + t102 * t106 + t189 * t632;
t559 = t13 * qJD(1);
t14 = t39 + (t560 - t563) * t335;
t558 = t14 * qJD(1);
t557 = t142 * t373;
t556 = t142 * t375;
t554 = t143 * t373;
t553 = t143 * t375;
t15 = (t73 - t553) * t337 + (-t72 - t556) * t335;
t552 = t15 * qJD(1);
t16 = (-t88 + t554) * t337 + (t87 + t557) * t335;
t551 = t16 * qJD(1);
t25 = (t105 - t633) * t337 + (t101 - t636) * t335;
t546 = t25 * qJD(1);
t26 = (-t106 + t634) * t337 + (t102 + t550) * t335;
t545 = t26 * qJD(1);
t543 = t314 * t375;
t542 = t315 * t337;
t541 = t315 * t373;
t427 = t603 + t370 / 0.2e1;
t414 = t427 * t337;
t536 = t335 * t357;
t434 = -t536 / 0.2e1;
t382 = t356 * t414 + t434;
t397 = t563 / 0.2e1 - t560 / 0.2e1;
t33 = t382 + t397;
t539 = t33 * qJD(1);
t534 = t337 * t357;
t533 = t341 * t373;
t532 = t341 * t375;
t309 = t370 * t337;
t310 = t371 * t337;
t41 = pkin(5) * t370 * t611 + t143 * t224 + t337 * t86;
t525 = t41 * qJD(1);
t220 = t335 * pkin(3) + t400;
t226 = t423 + t599;
t55 = t220 * t226;
t524 = t55 * qJD(1);
t85 = t360 * t369;
t521 = t85 * qJD(1);
t140 = t618 * t373;
t367 = qJD(4) * t375;
t520 = -t140 * qJD(3) + t367 * t610;
t519 = -t527 / 0.2e1 + t437;
t497 = qJD(3) * t335;
t517 = t373 * t497;
t471 = t335 * qJD(4);
t516 = t373 * t471;
t320 = t337 * qJD(3);
t515 = t222 * qJD(4) + t375 * t320;
t351 = t370 - t371;
t513 = t370 + t371;
t103 = -t220 * t337 - t226 * t335;
t512 = qJD(1) * t103;
t104 = t220 * t335 - t226 * t337;
t511 = qJD(1) * t104;
t508 = qJD(1) * t140;
t164 = t617 * t373;
t507 = qJD(1) * t164;
t165 = t617 * t375;
t506 = qJD(1) * t165;
t196 = t335 * t369 + t337 * t360;
t505 = qJD(1) * t196;
t197 = -t335 * t360 + t337 * t369;
t504 = qJD(1) * t197;
t502 = qJD(1) * t373;
t501 = qJD(1) * t376;
t500 = qJD(2) * t224;
t499 = qJD(2) * t630;
t498 = qJD(2) * t253;
t496 = qJD(4) * t337;
t366 = qJD(4) * t373;
t495 = qJD(5) * t222;
t494 = qJD(5) * t314;
t492 = qJD(5) * t375;
t491 = qJD(6) * t373;
t490 = qJD(6) * t375;
t361 = t369 / 0.2e1;
t116 = t361 + (pkin(3) / 0.2e1 - t359 / 0.2e1) * t337 + (qJ(4) / 0.2e1 + t357 / 0.2e1) * t335;
t489 = t116 * qJD(1);
t125 = t414 + t309 / 0.2e1 + t310 / 0.2e1;
t488 = t125 * qJD(1);
t141 = t513 * t337 * t335;
t487 = t141 * qJD(1);
t432 = t371 * t606;
t393 = t335 * t604 + t432;
t160 = t393 - t325;
t486 = t160 * qJD(1);
t166 = t618 * t375;
t163 = t166 * qJD(1);
t383 = (t372 * t606 + t576 * t605) * pkin(2);
t183 = -t369 / 0.2e1 + t383;
t485 = t183 * qJD(1);
t428 = t603 + t604;
t223 = t428 * t335;
t481 = t223 * qJD(5);
t480 = t224 * qJD(1);
t227 = -t309 - t310;
t479 = t227 * qJD(1);
t228 = t513 * t611;
t478 = t228 * qJD(1);
t475 = t325 * qJD(1);
t474 = t610 * qJD(1);
t472 = t335 * qJD(2);
t470 = t337 * qJD(1);
t469 = t337 * qJD(2);
t343 = -0.1e1 / 0.2e1 - t427;
t468 = t343 * qJD(2);
t467 = t513 * qJD(2);
t353 = -t374 ^ 2 + t376 ^ 2;
t466 = t353 * qJD(1);
t464 = t374 * qJD(2);
t463 = t376 * qJD(2);
t461 = pkin(1) * t374 * qJD(1);
t460 = pkin(1) * t501;
t459 = pkin(5) * t493;
t458 = pkin(5) * t492;
t457 = t600 / 0.2e1;
t455 = t596 / 0.2e1;
t449 = t220 * t470;
t447 = t610 * t502;
t446 = t335 * t364;
t445 = t373 * t365;
t444 = t337 * t493;
t443 = t337 * t492;
t303 = t335 * t492;
t442 = t335 * t491;
t441 = t335 * t470;
t248 = t335 * t469;
t440 = t374 * t463;
t355 = t373 * t492;
t301 = t375 * t470;
t439 = t337 * t490;
t438 = -t553 / 0.2e1;
t436 = t314 * t605;
t433 = -t533 / 0.2e1;
t430 = -t526 / 0.2e1;
t172 = qJD(5) * t325 + t441;
t419 = qJD(5) + t470;
t418 = t335 * t301;
t417 = t611 * t355;
t416 = pkin(4) / 0.2e1 + t455;
t413 = qJD(2) * t424;
t387 = t335 * t433 + t438;
t1 = t454 * t314 + (t609 + t387) * pkin(5);
t115 = pkin(5) * t532;
t411 = -qJD(1) * t1 + qJD(2) * t115;
t409 = t561 + t562;
t406 = t335 * t356 + t534;
t395 = t543 / 0.2e1 - t541 / 0.2e1;
t380 = t335 * t395 + t594;
t388 = -t339 / 0.2e1 - t415 / 0.2e1;
t19 = t337 * t416 + t380 + t388;
t194 = t540 + t544;
t405 = qJD(1) * t19 - qJD(2) * t194;
t21 = (t542 / 0.2e1 - t88 / 0.2e1) * t373 + (t436 - t600 / 0.2e1 - t73 / 0.2e1) * t375 + t377;
t404 = qJD(1) * t21 + qJD(2) * t341;
t27 = t438 + (-t528 / 0.2e1 + t314 / 0.2e1) * t337 + (t433 + (-0.1e1 + t428) * pkin(5)) * t335 + t519;
t270 = (t341 + t602) * t375;
t403 = -qJD(1) * t27 + qJD(2) * t270;
t281 = t371 * pkin(5) - t533;
t384 = -pkin(5) * t453 + t554 / 0.2e1 + t375 * t435;
t34 = (t430 + t315 / 0.2e1) * t337 + t384 + t518;
t402 = -qJD(1) * t34 + qJD(2) * t281;
t267 = t419 * t375;
t401 = -qJD(6) * t335 + t496;
t157 = t303 + t418;
t394 = t356 * t605 + t536 / 0.2e1;
t378 = t394 * t375 + t635;
t53 = t378 - t518;
t392 = -qJD(1) * t53 + t357 * t364;
t51 = t437 - t633 / 0.2e1 + (t159 / 0.2e1 + t394) * t373;
t391 = -qJD(1) * t51 - t357 * t365;
t174 = -qJD(1) * t223 + t445;
t144 = qJD(2) * t223 + t502 * t538;
t229 = t351 * t611;
t148 = qJD(1) * t229 + t413;
t240 = qJD(1) * t424 - qJD(2) * t351;
t354 = t374 * t501;
t346 = t357 * qJD(4);
t345 = t351 * qJD(5);
t342 = 0.1e1 / 0.2e1 - t427;
t306 = t325 * qJD(2);
t299 = t373 * t470;
t266 = -t301 - t492;
t264 = -t299 - t493;
t256 = t623 * pkin(5);
t218 = t227 * qJD(2);
t216 = t227 * qJD(3);
t215 = (-qJD(5) * t337 - t474) * t375;
t214 = t375 * t474 + t446;
t209 = t224 * qJD(3);
t208 = t224 * qJD(4);
t198 = t222 * qJD(3);
t191 = t531 + 0.2e1 * t389;
t182 = t361 + t383;
t162 = t166 * qJD(3);
t161 = -t393 - t325;
t158 = t622 * t337;
t156 = t419 * t535;
t155 = t623 * t337;
t153 = t248 * t371 - t417;
t152 = t248 * t370 + t417;
t146 = -0.2e1 * t267 * t535;
t137 = t141 * qJD(4);
t132 = -t443 + t506;
t131 = t446 - t506;
t130 = -t444 + t507;
t129 = -t335 * t365 - t507;
t128 = -t371 * t441 - t481;
t127 = -t370 * t441 + t481;
t124 = -qJD(5) * t229 + t337 * t413;
t122 = t125 * qJD(4);
t118 = t125 * qJD(3);
t117 = t434 + t359 * t337 / 0.2e1 + t361 + t537 / 0.2e1 + t599 / 0.2e1;
t114 = -qJD(2) * t165 - t335 * t444;
t113 = -qJD(2) * t164 + t303 * t337;
t110 = -t443 + t446 - t163;
t109 = -t481 + (t371 * t473 - t445) * t337;
t108 = t481 + (t370 * t473 + t445) * t337;
t107 = (t310 - t309) * qJD(2) + (-qJD(5) + t470) * t424;
t92 = t447 - t500;
t91 = -t447 - t495;
t62 = t495 + t500 + t508;
t54 = t378 + t518;
t52 = t373 * t394 + t437 + t519;
t36 = t377 - t396;
t35 = -t542 / 0.2e1 + t337 * t430 - t384 + t518;
t32 = t382 - t397;
t28 = pkin(5) * t432 + qJ(6) * t431 + t370 * t457 - t387 + t436 + t519 - t600;
t23 = t379 - t398;
t20 = t579 / 0.2e1 + t583 / 0.2e1 - t395 * t337 - t416 * t335 + t381;
t18 = -t598 / 0.2e1 + t596 * t605 + t380 - t388;
t12 = t375 * t631 + t581 / 0.2e1 + t337 * t455 + t594;
t11 = t420 + t612;
t5 = t421 - t612;
t2 = pkin(5) * t609 + t143 * t455 + t457 * t533 - t625 * t314 / 0.2e1;
t17 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t440, t353 * qJD(2), 0, -t440, 0, 0, -pkin(1) * t464, -pkin(1) * t463, 0, 0, -t248, t628, 0, t248, 0, 0, t196 * qJD(2), t197 * qJD(2), t626, qJD(2) * t85 + t614, 0, 0, 0, -t248, t628, t248, t626, qJD(2) * t103 + t337 * t471, qJD(2) * t104 + qJD(4) * t610, qJD(2) * t55 - t220 * t496 + t614, t152, t124, t113, t153, t114, -t248, qJD(2) * t25 + qJD(5) * t57 + t366 * t610 + t162, qJD(2) * t26 + qJD(5) * t56 + t520, qJD(2) * t14 - t137, qJD(2) * t13 + qJD(3) * t29 - qJD(4) * t39, t152, t124, t113, t153, t114, -t248, t15 * qJD(2) + t40 * qJD(5) + t222 * t401 + t162, qJD(2) * t16 + qJD(5) * t41 - t335 * t439 + t520, qJD(2) * t6 - qJD(5) * t8 + qJD(6) * t228 - t137, qJD(2) * t3 + qJD(3) * t22 + qJD(4) * t31 + qJD(5) * t7 - qJD(6) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t354, t466, t463, -t354, -t464, 0, -pkin(7) * t463 - t461, pkin(7) * t464 - t460, 0, 0, -t441, t629, -t472, t441, -t469, 0, -t499 + t505, -t498 + t504 (t335 * t576 - t337 * t372) * t590, t521 + (t253 * t372 - t576 * t630) * t590 + t182 * qJD(3), 0, t472, t469, -t441, t629, t441 (-t335 * t359 - t534) * qJD(2) - t471, t499 + t512, t498 + t511, t524 + (t253 * t357 + t359 * t630) * qJD(2) + t117 * qJD(3) + t191 * qJD(4), t108, t107, t129, t109, t131, -t172, t546 + (-t375 * t406 + t550) * qJD(2) - t208 + t52 * qJD(5), t545 + (t373 * t406 + t636) * qJD(2) + t54 * qJD(5) + t516, -qJD(2) * t409 + t558, t559 + (t189 * t357 + t356 * t409) * qJD(2) + t32 * qJD(3) + t36 * qJD(4), t108, t107, t129, t109, t131, -t172, t552 + (-t315 * t335 - t337 * t532 + t557) * qJD(2) - t208 + t28 * qJD(5) - t439, t551 + (t222 * t341 + t314 * t335 + t556) * qJD(2) + t35 * qJD(5) + t222 * qJD(6) + t516, t587 + (-t579 - t583 + (-t541 + t543) * t337) * qJD(2) + t5 * qJD(5), t589 + (t142 * t341 + t314 * t88 + t315 * t73) * qJD(2) + t23 * qJD(3) + t20 * qJD(4) + t2 * qJD(5) + t18 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t627, qJD(2) * t182 + t616, 0, 0, 0, 0, 0, 0, t627, 0, 0, qJD(2) * t117 + t616, 0, 0, 0, 0, 0, 0, t163, -t508, 0, qJD(2) * t32 + t573, 0, 0, 0, 0, 0, 0, t163, -t508, 0, qJD(2) * t23 + qJD(5) * t12 + qJD(6) * t161 + t575; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t472, t441, t474, qJD(2) * t191 - t449, 0, 0, 0, 0, 0, 0, t92, t214, -t487, qJD(2) * t36 - t569, 0, 0, 0, 0, 0, 0, t92, t214, -t487, qJD(2) * t20 + qJD(5) * t11 + t571; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, -t148, t157, -t144, -t156, -t306, qJD(2) * t52 - qJD(5) * t102 + t566, qJD(2) * t54 + qJD(5) * t101 + t567, 0, 0, t144, -t148, t157, -t144, -t156, -t306, qJD(2) * t28 - t565 + t568, qJD(2) * t35 + qJD(5) * t86 + t525, -pkin(5) * t303 + qJD(2) * t5 - t591, -pkin(5) * t565 + qJD(2) * t2 + qJD(3) * t12 + qJD(4) * t11 + t592; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, qJD(2) * t222 - t418, t478, qJD(2) * t18 + qJD(3) * t161 - t572; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t354, -t466, 0, t354, 0, 0, t461, t460, 0, 0, t441, -t629, 0, -t441, 0, 0, -t320 - t505, t497 - t504, 0, qJD(3) * t183 - t521, 0, 0, 0, t441, -t629, -t441, 0, t320 - t512, -t497 - t511, -qJD(3) * t116 + qJD(4) * t250 - t524, t127, t146, t130, t128, t132, t172, qJD(5) * t51 - t517 - t546, qJD(5) * t53 - t209 - t545, t216 - t558, qJD(3) * t33 + qJD(4) * t37 - t559, t127, t146, t130, t128, t132, t172, -qJD(5) * t27 - t517 - t552, -qJD(5) * t34 - t209 - t551, -qJD(5) * t4 + t216 - t587, qJD(3) * t24 + qJD(4) * t21 - qJD(5) * t1 + qJD(6) * t19 - t589; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t346, -t355, t345, 0, t355, 0, 0, t357 * t492 + t366, -t357 * t493 + t367, 0, t346, -t355, t345, 0, t355, 0, 0, qJD(5) * t270 + t366, qJD(5) * t281 + t367, t513 * qJD(6), qJD(4) * t341 + qJD(5) * t115 - qJD(6) * t194; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t470, t473, 0, t485, 0, 0, 0, 0, 0, 0, 0, t470, -t473, -t489, 0, 0, 0, 0, 0, 0, -t298, -t480, t479, t539, 0, 0, 0, 0, 0, 0, -t298, -t480, t479, t574; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t620, 0, 0, 0, 0, 0, 0, t364, t365, 0, t619, 0, 0, 0, 0, 0, 0, t364, t365, 0, qJD(6) * t342 + t404; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t174, -t240, t264, t174, t266, t475, -t356 * t493 - t391, -t356 * t492 - t392, 0, 0, -t174, -t240, t264, t174, t266, t475, t403 - t494, -qJD(5) * t315 + t402, t459 - t588, -pkin(5) * t494 + t411; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t467, qJD(4) * t342 + t405; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t469, -t472, -t627, -qJD(2) * t183 - t616, 0, 0, 0, 0, 0, 0, -t627, -t469, t472, qJD(2) * t116 - t496 - t616, 0, 0, 0, 0, 0, 0, t110, t62, -t218, -qJD(2) * t33 - t122 - t573, 0, 0, 0, 0, 0, 0, t110, t62, -t218, -qJD(2) * t24 - qJD(5) * t10 - qJD(6) * t160 - t122 - t575; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t470, -t473, 0, -t485, 0, 0, 0, 0, 0, 0, 0, -t470, t473, t489, 0, 0, 0, 0, 0, 0, t298, t480, -t479, -t539, 0, 0, 0, 0, 0, 0, t298, t480, -t479, -t574; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t470, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t488, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t488; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t266, t179, 0, 0, 0, 0, 0, 0, 0, 0, t266, t179, 0, -t458 - t564; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t486; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t441, -t474, -qJD(2) * t250 + t320 + t449, 0, 0, 0, 0, 0, 0, t91, t215, t487, -qJD(2) * t37 + t118 + t569, 0, 0, 0, 0, 0, 0, t91, t215, t487, -qJD(2) * t21 - qJD(5) * t9 + t118 - t571; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t620, 0, 0, 0, 0, 0, 0, -t364, -t365, 0, -t619, 0, 0, 0, 0, 0, 0, -t364, -t365, 0, qJD(6) * t343 - t404; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t470, 0, 0, 0, 0, 0, 0, 0, 0, 0, t488, 0, 0, 0, 0, 0, 0, 0, 0, 0, t488; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179, -t267, 0, 0, 0, 0, 0, 0, 0, 0, -t179, -t267, 0, -t459 - t577; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t468; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, t148, t158, t144, t155, -t306, -qJD(2) * t51 + t515 - t566, -qJD(2) * t53 + t337 * t367 - t198 - t567, 0, 0, -t144, t148, t158, t144, t155, -t306, qJD(2) * t27 - t442 + t515 - t568, t34 * qJD(2) + t375 * t401 - t198 - t525, qJD(2) * t4 + t591, -pkin(5) * t442 + qJD(2) * t1 + qJD(3) * t10 + qJD(4) * t9 - t592; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t174, t240, t299, -t174, t301, -t475, t391, t392, 0, 0, t174, t240, t299, -t174, t301, -t475, -t403 - t490, -t402 + t491, t588, -pkin(5) * t490 - t411; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t301, -t482, 0, 0, 0, 0, 0, 0, 0, 0, t301, -t482, 0, t564; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t482, t301, 0, 0, 0, 0, 0, 0, 0, 0, t482, t301, 0, t577; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t623, t622, 0, -t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, t157, -t478, -qJD(2) * t19 + qJD(3) * t160 + t335 * t459 + t572; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t492, -t493, -t467, -qJD(4) * t343 - t405 + t458; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t486; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t468; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t623, -t622, 0, t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t17;
