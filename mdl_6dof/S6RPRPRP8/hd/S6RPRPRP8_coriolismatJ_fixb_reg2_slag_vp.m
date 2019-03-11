% Calculate inertial parameters regressor of coriolis matrix for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRPRP8_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:26:11
% EndTime: 2019-03-09 03:26:26
% DurationCPUTime: 10.13s
% Computational Cost: add. (9537->536), mult. (16753->623), div. (0->0), fcn. (18065->6), ass. (0->410)
t390 = cos(qJ(5));
t387 = cos(pkin(9));
t389 = sin(qJ(3));
t391 = cos(qJ(3));
t612 = sin(pkin(9));
t356 = t387 * t389 + t391 * t612;
t354 = t356 ^ 2;
t370 = t612 * t389;
t554 = t387 * t391;
t358 = -t370 + t554;
t355 = t358 ^ 2;
t484 = -t355 - t354;
t649 = t484 * t390;
t655 = qJD(1) * t649;
t388 = sin(qJ(5));
t650 = t484 * t388;
t654 = qJD(1) * t650;
t653 = t650 * qJD(4);
t525 = qJD(2) * t356;
t652 = -qJD(4) * t649 - t388 * t525;
t392 = -pkin(1) - pkin(7);
t546 = -qJ(4) + t392;
t360 = t546 * t389;
t361 = t546 * t391;
t277 = t360 * t612 - t387 * t361;
t647 = t277 * t388;
t651 = t647 / 0.2e1;
t349 = t358 * qJ(6);
t617 = t391 * pkin(3);
t245 = pkin(4) * t358 + pkin(8) * t356 + t617;
t201 = t388 * t245;
t646 = t277 * t390;
t544 = t201 / 0.2e1 - t646 / 0.2e1;
t648 = -t349 - t544;
t642 = t387 * t360 + t612 * t361;
t568 = t642 * t388;
t567 = t642 * t390;
t486 = -t355 / 0.2e1;
t645 = t486 - t354 / 0.2e1;
t446 = t354 / 0.2e1 + t355 / 0.2e1;
t238 = t388 * t358;
t517 = qJD(5) * t390;
t438 = -t238 * qJD(3) - t356 * t517;
t424 = 0.1e1 / 0.2e1 + t446;
t408 = t424 * t390;
t637 = qJD(1) * t408;
t644 = t637 - t438;
t616 = qJD(3) * pkin(3);
t643 = (-t356 * t387 + t358 * t612) * t616;
t483 = t355 - t354;
t239 = t390 * t358;
t518 = qJD(5) * t388;
t477 = t356 * t518;
t437 = t239 * qJD(3) - t477;
t233 = t388 * t356;
t507 = t233 * qJD(1);
t187 = t507 + t518;
t205 = t233 * qJD(5);
t173 = t483 * t388;
t511 = t173 * qJD(1);
t641 = t511 - t205;
t347 = t358 * qJD(3);
t331 = t390 * t347;
t640 = t511 - t331;
t385 = t388 ^ 2;
t386 = t390 ^ 2;
t538 = t385 + t386;
t108 = (-0.1e1 + t538) * t358 * t356;
t515 = t108 * qJD(2);
t373 = pkin(3) * t389 + qJ(2);
t404 = pkin(4) * t356 - pkin(8) * t358 + t373;
t130 = t388 * t404 + t567;
t586 = t130 * t390;
t129 = -t390 * t404 + t568;
t591 = t129 * t388;
t415 = t591 / 0.2e1 + t586 / 0.2e1;
t132 = -t646 + t201;
t582 = t132 * t390;
t550 = t390 * t245;
t131 = t550 + t647;
t585 = t131 * t388;
t629 = -t642 / 0.2e1;
t630 = t277 / 0.2e1;
t9 = (t629 + t415) * t358 + (t582 / 0.2e1 - t585 / 0.2e1 + t630) * t356;
t639 = -qJD(1) * t9 - t515;
t332 = pkin(5) * t238;
t444 = -qJ(6) * t239 + t332;
t139 = t277 + t444;
t551 = t390 * qJ(6);
t621 = pkin(5) * t388;
t362 = t551 - t621;
t140 = t356 * t362 + t642;
t620 = t356 * pkin(5);
t107 = t129 - t620;
t598 = t107 * t388;
t562 = t356 * qJ(6);
t106 = t130 + t562;
t599 = t106 * t390;
t417 = t599 / 0.2e1 + t598 / 0.2e1;
t619 = t358 * pkin(5);
t112 = -t131 - t619;
t594 = t112 * t388;
t111 = t349 + t132;
t595 = t111 * t390;
t4 = (-t140 / 0.2e1 + t417) * t358 + (t595 / 0.2e1 + t139 / 0.2e1 + t594 / 0.2e1) * t356;
t638 = qJD(1) * t4 + t515;
t636 = qJD(2) * t408;
t153 = 0.1e1 / 0.2e1 - t446;
t497 = t356 * qJD(1);
t635 = -t497 - qJD(5);
t336 = t385 * t356;
t337 = t386 * t356;
t634 = 0.2e1 * t388 * t239 * (-qJD(5) + t497) - (-t336 + t337) * qJD(3);
t330 = t390 * t497;
t633 = -qJD(5) * t408 - t330;
t632 = t130 / 0.2e1;
t631 = -t245 / 0.2e1;
t628 = -t332 / 0.2e1;
t627 = -t358 / 0.2e1;
t626 = t358 / 0.2e1;
t625 = -t362 / 0.2e1;
t624 = -t386 / 0.2e1;
t381 = -t388 / 0.2e1;
t623 = t390 / 0.2e1;
t622 = t9 * qJD(3);
t618 = t390 * pkin(5);
t553 = t388 * qJ(6);
t435 = t553 + t618;
t222 = t435 * t358;
t460 = -t129 / 0.2e1 + t107 / 0.2e1;
t397 = t460 * t390 + (-t106 / 0.2e1 + t632) * t388;
t395 = t222 * t627 + t356 * t397;
t419 = t553 / 0.2e1 + t618 / 0.2e1;
t7 = t395 - t419;
t614 = t7 * qJD(1);
t24 = t139 * t358 + (-t598 - t599) * t356;
t611 = qJD(1) * t24;
t571 = t277 * t358;
t31 = t571 + (-t586 - t591) * t356;
t610 = qJD(1) * t31;
t580 = t139 * t390;
t588 = t130 * t356;
t32 = -t588 + (t222 * t388 + t580) * t358;
t609 = qJD(1) * t32;
t581 = t139 * t388;
t592 = t129 * t356;
t33 = -t592 + (-t222 * t390 + t581) * t358;
t608 = qJD(1) * t33;
t601 = t106 * t356;
t38 = -t139 * t239 + t601;
t607 = qJD(1) * t38;
t597 = t107 * t390;
t600 = t106 * t388;
t39 = -t597 + t600;
t606 = qJD(1) * t39;
t587 = t130 * t388;
t590 = t129 * t390;
t430 = -t587 + t590;
t605 = qJD(1) * t430;
t57 = -t238 * t277 + t592;
t604 = qJD(1) * t57;
t58 = t239 * t277 - t588;
t603 = qJD(1) * t58;
t10 = -t106 * t129 + t107 * t130 + t139 * t222;
t602 = t10 * qJD(1);
t236 = t390 * t356;
t11 = t107 * t236 - t112 * t239 + (t111 * t358 - t601) * t388;
t596 = t11 * qJD(1);
t12 = ((t106 - t130) * t390 + (t107 - t129) * t388) * t358;
t593 = t12 * qJD(1);
t418 = t620 / 0.2e1 - t460;
t441 = t562 / 0.2e1 + t106 / 0.2e1;
t13 = -t587 / 0.2e1 + t441 * t388 + t418 * t390;
t589 = t13 * qJD(1);
t584 = t131 * t390;
t583 = t132 * t388;
t14 = (t632 - t441) * t390 + t418 * t388;
t579 = t14 * qJD(1);
t578 = t140 * t388;
t577 = t140 * t390;
t17 = (t106 - t577) * t358 + (t111 + t580) * t356;
t576 = t17 * qJD(1);
t18 = (-t107 + t578) * t358 + (-t112 - t581) * t356;
t575 = t18 * qJD(1);
t19 = (t583 + t584) * t358 + t430 * t356;
t574 = t19 * qJD(1);
t27 = (-t129 + t568) * t358 + (t131 - t647) * t356;
t573 = t27 * qJD(1);
t28 = (-t130 + t567) * t358 + (-t132 - t646) * t356;
t566 = t28 * qJD(1);
t458 = t624 - t385 / 0.2e1;
t371 = pkin(3) * t612 + pkin(8);
t557 = t371 * t356;
t421 = t458 * t557;
t372 = -pkin(3) * t387 - pkin(4);
t350 = t372 - t435;
t467 = t350 * t626;
t401 = t421 + t467;
t416 = t111 * t381 + t112 * t623;
t30 = t401 + t416;
t565 = t30 * qJD(1);
t559 = t358 * t372;
t465 = t559 / 0.2e1;
t400 = t421 + t465;
t414 = -t584 / 0.2e1 - t583 / 0.2e1;
t35 = t400 + t414;
t564 = t35 * qJD(1);
t563 = t350 * t356;
t561 = t356 * t372;
t560 = t358 * t371;
t558 = t362 * t388;
t556 = t385 * t358;
t555 = t386 * t358;
t459 = -t277 / 0.2e1 + t630;
t44 = (t642 / 0.2e1 + t629) * t358 + t459 * t356;
t549 = t44 * qJD(1);
t407 = t354 * t458 + t486;
t84 = t407 + t458;
t547 = t84 * qJD(1);
t461 = t239 / 0.2e1;
t545 = t581 / 0.2e1 + t350 * t461;
t542 = t550 / 0.2e1 + t651;
t541 = (t555 + t556) * t371;
t540 = -t556 / 0.2e1 + t555 / 0.2e1;
t105 = -t356 * t642 + t571;
t537 = qJD(1) * t105;
t449 = t645 * t388;
t141 = t381 + t449;
t536 = qJD(1) * t141;
t175 = t483 * t390;
t533 = qJD(1) * t175;
t199 = t356 * t617 + t358 * t373;
t531 = qJD(1) * t199;
t200 = -t356 * t373 + t358 * t617;
t530 = qJD(1) * t200;
t527 = qJD(1) * t373;
t167 = t153 * t390;
t526 = qJD(2) * t167;
t524 = qJD(3) * t388;
t523 = qJD(3) * t390;
t522 = qJD(4) * t388;
t521 = qJD(4) * t390;
t520 = qJD(5) * t129;
t519 = qJD(5) * t356;
t516 = qJD(6) * t388;
t514 = t424 * qJD(1);
t380 = t388 / 0.2e1;
t448 = t446 * t388;
t164 = t380 + t448;
t513 = t164 * qJD(1);
t405 = -t612 * t356 / 0.2e1 + t387 * t627;
t191 = (-t391 / 0.2e1 + t405) * pkin(3);
t510 = t191 * qJD(1);
t509 = t483 * qJD(1);
t231 = (t385 / 0.2e1 + t624) * t358;
t508 = t231 * qJD(5);
t213 = t236 * qJD(1);
t506 = t238 * qJD(1);
t504 = t239 * qJD(1);
t247 = t336 + t337;
t501 = t247 * qJD(1);
t249 = t538 * t358;
t500 = t249 * qJD(1);
t499 = t484 * qJD(1);
t351 = t554 / 0.2e1 - t370 / 0.2e1;
t498 = t351 * qJD(1);
t344 = t356 * qJD(3);
t342 = t356 * qJD(6);
t496 = t358 * qJD(1);
t365 = t389 ^ 2 - t391 ^ 2;
t495 = t365 * qJD(1);
t366 = t386 - t385;
t494 = t366 * qJD(5);
t493 = t389 * qJD(1);
t492 = t389 * qJD(3);
t491 = t390 * qJD(6);
t490 = t391 * qJD(1);
t489 = t391 * qJD(3);
t339 = t619 / 0.2e1;
t462 = t557 / 0.2e1;
t482 = -t581 / 0.2e1 - t350 * t239 / 0.2e1 + t390 * t462;
t481 = t339 + t542;
t480 = qJ(2) * t493;
t479 = qJ(2) * t490;
t475 = t358 * t517;
t474 = t371 * t518;
t473 = t371 * t517;
t472 = t356 * t496;
t288 = t356 * t347;
t368 = t388 * t517;
t328 = t388 * t497;
t471 = t388 * t496;
t470 = t358 * t516;
t367 = t388 * t523;
t469 = t390 * t496;
t468 = t389 * t489;
t466 = t362 * t626;
t464 = -t558 / 0.2e1;
t463 = -t557 / 0.2e1;
t169 = t351 + t540;
t456 = qJD(1) * t169 + t367;
t179 = qJD(1) * t231 - t367;
t306 = t390 * t355 * t388 * qJD(1);
t151 = qJD(3) * t231 + t306;
t273 = t356 * t469;
t455 = qJD(3) * t233 + t273;
t452 = t388 * t469;
t451 = t388 * t331;
t450 = t355 * t368;
t447 = -t619 / 0.2e1 - t647 / 0.2e1;
t443 = 0.2e1 * t451;
t442 = -t222 / 0.2e1 + t463;
t440 = -qJD(5) * t238 - t344 * t390;
t436 = t273 + t475;
t5 = t106 * t111 + t107 * t112 + t139 * t140;
t434 = t5 * qJD(1) + t4 * qJD(2);
t20 = -t129 * t131 + t130 * t132 + t277 * t642;
t433 = t20 * qJD(1) + t9 * qJD(2);
t59 = t373 * t617;
t432 = t59 * qJD(1) + t44 * qJD(2);
t431 = t594 + t595;
t429 = t582 - t585;
t428 = t560 + t563;
t427 = -t560 - t561;
t21 = (t464 - pkin(5) / 0.2e1) * t358 + (t631 + t442) * t390 + t447 + t545;
t257 = -t350 * t388 - t362 * t390;
t426 = -qJD(1) * t21 + qJD(3) * t257;
t256 = t350 * t390 - t558;
t394 = (t466 - t139 / 0.2e1) * t390 + (t467 + t442) * t388;
t26 = t394 + t648;
t425 = -qJD(1) * t26 + qJD(3) * t256;
t423 = qJD(4) * t233 + t636;
t422 = -qJD(5) * t362 - t516;
t420 = -t111 * qJ(6) / 0.2e1 + t112 * pkin(5) / 0.2e1;
t37 = t481 + t482;
t413 = -qJD(1) * t37 + t350 * t524;
t298 = t390 * t463;
t54 = t298 + (t465 + t631) * t390 + t459 * t388;
t412 = -qJD(1) * t54 - t372 * t524;
t396 = (t462 - t559 / 0.2e1) * t388 + t646 / 0.2e1;
t52 = t396 + t544;
t411 = -qJD(1) * t52 - t372 * t523;
t109 = t635 * t238;
t177 = qJD(5) * t351 + t472;
t410 = -qJD(5) * t239 + t344 * t388;
t338 = t386 * t355;
t248 = t355 * t385 - t338;
t156 = -qJD(1) * t248 + t443;
t274 = -qJD(3) * t366 + 0.2e1 * t452;
t406 = qJD(3) * t173 + t356 * t475;
t403 = qJD(5) * t248 + t356 * t443;
t149 = t628 + (t551 / 0.2e1 + t625) * t358;
t393 = t397 * t371 + t139 * t625 + t222 * t350 / 0.2e1;
t2 = t393 + t420;
t402 = t350 * t362 * qJD(3) - t2 * qJD(1) + t149 * qJD(2);
t399 = -qJD(5) * t435 + t491;
t258 = t338 + t354;
t398 = qJD(1) * t258 + t451 + t519;
t384 = qJ(2) * qJD(2);
t383 = qJD(1) * qJ(2);
t369 = t389 * t490;
t333 = t351 * qJD(3);
t329 = t388 * t347;
t307 = t390 * t470;
t305 = t635 * qJ(6);
t291 = qJD(3) * t385 + t452;
t229 = t247 * qJD(4);
t227 = t249 * qJD(2);
t225 = t249 * qJD(3);
t223 = t247 * qJD(3);
t212 = t236 * qJD(5);
t190 = t617 / 0.2e1 + t405 * pkin(3);
t189 = t213 + t517;
t170 = -t351 + t540;
t165 = t380 + t449;
t163 = -t288 * t386 - t450;
t162 = -t288 * t385 + t450;
t157 = t390 * t109;
t150 = qJ(6) * t461 + t466 + t628;
t144 = t390 * t645 + t623;
t143 = t381 + t448;
t138 = t386 * t472 - t508;
t137 = -qJD(3) * t236 + t356 * t471;
t136 = t385 * t472 + t508;
t128 = t130 * qJD(5);
t116 = qJD(3) * t175 - t358 * t477;
t110 = t212 - t533;
t102 = t108 * qJD(3);
t98 = -t508 + (-t386 * t496 - t367) * t356;
t97 = t508 + (-t385 * t496 + t367) * t356;
t83 = t407 - t458;
t79 = t84 * qJD(2);
t77 = t84 * qJD(4);
t76 = t83 * qJD(2);
t75 = t83 * qJD(4);
t70 = t329 + t533;
t55 = t372 * t461 + t298 + t542 + t651;
t53 = t396 - t544;
t42 = t44 * qJD(3);
t36 = -t550 / 0.2e1 + t447 + t482;
t34 = t400 - t414;
t29 = t401 - t416;
t25 = t394 - t648;
t22 = t358 * t464 + t390 * t442 + t339 + t481 + t545;
t16 = -t590 / 0.2e1 - t600 / 0.2e1 + t587 / 0.2e1 + t597 / 0.2e1 + t419 * t356;
t15 = (-t551 / 0.2e1 + t621 / 0.2e1) * t356 - t415 + t417;
t6 = t395 + t419;
t3 = t4 * qJD(3);
t1 = t393 - t420;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t384, -t468, t365 * qJD(3), 0, t468, 0, 0, qJ(2) * t489 + qJD(2) * t389, -qJ(2) * t492 + qJD(2) * t391, 0, t384, -t288, -t483 * qJD(3), 0, t288, 0, 0, qJD(3) * t199 + t525, qJD(2) * t358 + qJD(3) * t200, -qJD(4) * t484, qJD(2) * t373 + qJD(3) * t59 + qJD(4) * t105, t163, t403, t116, t162, -t406, t288, qJD(3) * t27 + qJD(5) * t58 + t390 * t525 - t653, qJD(3) * t28 + qJD(5) * t57 + t652, -qJD(3) * t19 - t227, -qJD(2) * t430 + qJD(3) * t20 + qJD(4) * t31, t163, t116, -t403, t288, t406, t162, t18 * qJD(3) - t653 + t32 * qJD(5) + (-t355 * t516 + t525) * t390, -qJD(3) * t11 - qJD(5) * t12 - t356 * t470 - t227, qJD(3) * t17 + qJD(5) * t33 + qJD(6) * t258 - t652, qJD(2) * t39 + qJD(3) * t5 + qJD(4) * t24 + qJD(5) * t10 + qJD(6) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), t383, 0, 0, 0, 0, 0, 0, t493, t490, 0, t383, 0, 0, 0, 0, 0, 0, t497, t496, 0, qJD(4) * t153 + t42 + t527, 0, 0, 0, 0, 0, 0, qJD(5) * t167 + t330, qJD(5) * t143 - t328, -t500, t75 - t605 + t622, 0, 0, 0, 0, 0, 0, qJD(5) * t144 + t330, -t500, qJD(5) * t165 + t328, qJD(5) * t6 - qJD(6) * t167 + t3 + t606 + t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t369, t495, -t492, t369, -t489, 0, -t392 * t492 + t479, -t392 * t489 - t480, 0, 0, -t472, -t509, -t344, t472, -t347, 0, -qJD(3) * t642 + t531, qJD(3) * t277 + t530, -t643 (-t277 * t612 - t387 * t642) * t616 + t190 * qJD(4) + t432, t98, t634, t70, t97, -t640, t177, t573 + (t388 * t427 - t567) * qJD(3) + t55 * qJD(5), t566 + (t390 * t427 + t568) * qJD(3) + t53 * qJD(5), qJD(3) * t429 - t574 (t371 * t429 + t372 * t642) * qJD(3) + t34 * qJD(4) + t433, t98, t70, -t634, t177, t640, t97, t575 + (-t388 * t428 - t577) * qJD(3) + t22 * qJD(5) + t170 * qJD(6), qJD(3) * t431 + qJD(5) * t16 - t596, t576 + (t390 * t428 - t578) * qJD(3) + t25 * qJD(5) + t307 (t140 * t350 + t371 * t431) * qJD(3) + t29 * qJD(4) + t1 * qJD(5) + t36 * qJD(6) + t434; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t499, qJD(2) * t153 + qJD(3) * t190 + t537, 0, 0, 0, 0, 0, 0, -t654, -t655, 0, qJD(3) * t34 + t610 + t76, 0, 0, 0, 0, 0, 0, -t654, 0, t655, qJD(3) * t29 + qJD(5) * t15 + t611 + t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151, -t156, t109, t151, -t436, t333, qJD(3) * t55 - t128 + t526 + t603, qJD(2) * t143 + qJD(3) * t53 + t520 + t604, 0, 0, -t151, t109, t156, t333, t436, t151, qJD(2) * t144 + qJD(3) * t22 - t128 + t609, t16 * qJD(3) + qJD(5) * t444 - t470 - t593, qJD(2) * t165 + qJD(3) * t25 + t342 - t520 + t608, t602 + t6 * qJD(2) + t1 * qJD(3) + t15 * qJD(4) + (-pkin(5) * t130 - qJ(6) * t129) * qJD(5) + t106 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t170 - t306, t109, t398, qJD(3) * t36 + qJD(5) * t106 - t526 + t607; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), -t383, 0, 0, 0, 0, 0, 0, -t493, -t490, 0, -t383, 0, 0, 0, 0, 0, 0, -t497, -t496, 0, -qJD(4) * t424 + t42 - t527, 0, 0, 0, 0, 0, 0, t633, -qJD(5) * t141 + t328, t500, t77 + t605 + t622, 0, 0, 0, 0, 0, 0, t633, t500, -qJD(5) * t164 - t328, qJD(5) * t7 + qJD(6) * t408 + t3 - t606 + t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t492, -t489, 0, 0, 0, 0, 0, 0, 0, 0, -t344, -t347, 0, t549 + t643, 0, 0, 0, 0, 0, 0, t440, t410, t225 (t541 + t561) * qJD(3) - t639, 0, 0, 0, 0, 0, 0, t440, t225, -t410 (t541 + t563) * qJD(3) + t150 * qJD(5) + t238 * qJD(6) + t638; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t514, 0, 0, 0, 0, 0, 0, 0, 0, 0, t547, 0, 0, 0, 0, 0, 0, 0, 0, 0, t547; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t644, -t536 - t437, 0, 0, 0, 0, 0, 0, 0, 0, -t644, 0, t437 - t513, t150 * qJD(3) + t356 * t399 + t614; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t644; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t369, -t495, 0, -t369, 0, 0, -t479, t480, 0, 0, t472, t509, 0, -t472, 0, 0, -qJD(4) * t358 - t531, qJD(4) * t356 - t530, 0, qJD(4) * t191 - t432, t138, 0.2e1 * t157, t110, t136, t641, -t177, qJD(5) * t54 - t358 * t521 - t573, qJD(4) * t238 + qJD(5) * t52 - t566, -t229 + t574, qJD(4) * t35 - t433, t138, t110, -0.2e1 * t157, -t177, -t641, t136, -qJD(4) * t239 + qJD(5) * t21 + qJD(6) * t169 - t575, -qJD(5) * t13 + qJD(6) * t236 - t229 + t596, qJD(5) * t26 - t358 * t522 + t307 - t576, qJD(4) * t30 + qJD(5) * t2 + qJD(6) * t37 - t434; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t549, 0, 0, 0, 0, 0, 0, 0, 0, 0, t639, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5) * t149 - t638; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t368, t494, 0, -t368, 0, 0, t372 * t518, t372 * t517, 0, 0, t368, 0, -t494, 0, 0, -t368, -qJD(5) * t257 + t388 * t491, 0, -qJD(5) * t256 + qJD(6) * t385, t422 * t350; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t496, t497, 0, t510, 0, 0, 0, 0, 0, 0, -t469, t506, -t501, t564, 0, 0, 0, 0, 0, 0, -t504, -t501, -t471, t565; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t179, -t274, t189, t179, -t187, -t498, -t412 - t473, -t411 + t474, 0, 0, -t179, t189, t274, -t498, t187, t179, -t426 - t473, t399 - t589, -t425 - t474, t371 * t399 - t402; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t456, t189, t291, -t413 + t473; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t347, -t344, t499, qJD(2) * t424 - qJD(3) * t191 - t537, 0, 0, 0, 0, 0, 0, -t205 + t331 + t654, t438 + t655, t223, -qJD(3) * t35 - t610 - t79, 0, 0, 0, 0, 0, 0, t437 + t654, t223, t212 + t329 - t655, -qJD(3) * t30 - qJD(5) * t14 + qJD(6) * t233 - t611 - t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t514, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t547, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t547; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t496, -t497, 0, -t510, 0, 0, 0, 0, 0, 0, t469, -t506, t501, -t564, 0, 0, 0, 0, 0, 0, t504, t501, t471, -t565; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t187, -t330 - t517, 0, 0, 0, 0, 0, 0, 0, 0, -t328 - t518, 0, t189, -t422 - t579; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, t156, t137, -t151, t455, t333, -qJD(3) * t54 + t423 - t603, qJD(2) * t141 - qJD(3) * t52 + t356 * t521 - t604, 0, 0, t151, t137, -t156, t333, -t455, -t151, -qJD(3) * t21 + t356 * t522 - t609 + t636, qJD(3) * t13 + t593, qJD(2) * t164 - qJD(3) * t26 - qJD(4) * t236 + t342 - t608, qJ(6) * t342 - qJD(2) * t7 - qJD(3) * t2 + qJD(4) * t14 - t602; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t637, t536, 0, 0, 0, 0, 0, 0, 0, 0, t637, 0, t513, qJD(3) * t149 - t614; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, t274, -t213, -t179, t507, t498, t412, t411, 0, 0, t179, -t213, -t274, t498, -t507, -t179, t426, t589, t425, t402; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t507, t330, 0, 0, 0, 0, 0, 0, 0, 0, t328, 0, -t213, t579; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(6), qJ(6) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t635, -t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t169 + t306, t137, -t398, -qJ(6) * t519 - qJD(3) * t37 - t423 - t607; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t637; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t456, -t213, -t291, t413; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t507; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t635, t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t8;
