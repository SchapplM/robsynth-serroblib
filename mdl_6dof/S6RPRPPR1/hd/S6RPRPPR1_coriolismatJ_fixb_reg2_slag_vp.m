% Calculate inertial parameters regressor of coriolis matrix for
% S6RPRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRPPR1_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:39:42
% EndTime: 2019-03-09 02:39:57
% DurationCPUTime: 10.67s
% Computational Cost: add. (15386->462), mult. (29280->621), div. (0->0), fcn. (33809->10), ass. (0->372)
t402 = sin(pkin(10));
t405 = sin(qJ(3));
t406 = cos(qJ(3));
t616 = cos(pkin(10));
t376 = t402 * t405 - t406 * t616;
t401 = sin(pkin(11));
t629 = cos(qJ(6));
t500 = t629 * t401;
t403 = cos(pkin(11));
t404 = sin(qJ(6));
t566 = t404 * t403;
t560 = t376 * (-t566 / 0.2e1 - t500 / 0.2e1);
t475 = t616 * t405;
t570 = t402 * t406;
t372 = t475 / 0.2e1 + t570 / 0.2e1;
t380 = t475 + t570;
t399 = t401 ^ 2;
t400 = t403 ^ 2;
t478 = t400 / 0.2e1 + t399 / 0.2e1;
t643 = t478 * t380;
t255 = t643 - t372;
t655 = qJD(5) * t255;
t254 = t643 + t372;
t654 = t254 * qJD(5);
t393 = pkin(3) * t402 + qJ(5);
t626 = pkin(8) + t393;
t364 = t626 * t403;
t477 = t626 * t401;
t280 = t364 * t404 + t477 * t629;
t281 = t364 * t629 - t404 * t477;
t382 = t500 + t566;
t632 = t382 / 0.2e1;
t200 = t382 * t380;
t637 = -t200 / 0.2e1;
t499 = t629 * t403;
t567 = t404 * t401;
t378 = -t499 + t567;
t575 = t380 * t378;
t638 = -t575 / 0.2e1;
t652 = -t378 / 0.2e1;
t466 = sin(pkin(9)) * pkin(1) + pkin(7);
t442 = qJ(4) + t466;
t427 = t402 * t442;
t420 = t405 * t427;
t426 = t406 * t442;
t409 = t426 * t616 - t420;
t397 = -cos(pkin(9)) * pkin(1) - pkin(2);
t385 = -t406 * pkin(3) + t397;
t411 = t376 * pkin(4) - t380 * qJ(5) + t385;
t145 = t401 * t411 + t403 * t409;
t288 = t401 * t380;
t126 = -pkin(8) * t288 + t145;
t144 = -t401 * t409 + t403 * t411;
t628 = pkin(8) * t403;
t407 = t376 * pkin(5) - t380 * t628 + t144;
t67 = t126 * t404 - t407 * t629;
t68 = t126 * t629 + t404 * t407;
t653 = -t280 * t638 - t281 * t637 - t67 * t632 - t652 * t68;
t483 = t200 * t632;
t342 = t380 * t567;
t277 = t380 * t499 - t342;
t587 = t277 * t378;
t118 = -t483 - t587 / 0.2e1;
t650 = t118 * qJD(6);
t649 = t200 * qJD(1);
t540 = qJD(3) * t382;
t648 = -qJD(1) * t118 + t378 * t540;
t549 = qJD(1) * t277;
t647 = qJD(3) * t118 - t200 * t549;
t646 = t200 ^ 2;
t341 = t376 * t567;
t276 = -t376 * t499 + t341;
t237 = t276 * t378;
t571 = t382 * t376;
t572 = t382 * t571;
t81 = -t237 - t572;
t645 = qJD(3) * t81;
t644 = t200 * t378;
t510 = t376 * qJD(1);
t496 = t200 * t510;
t373 = t376 ^ 2;
t375 = t380 ^ 2;
t642 = -t375 - t373;
t504 = t375 - t373;
t640 = qJD(3) * t378 + t649;
t482 = t571 / 0.2e1;
t461 = t482 - t560;
t544 = qJD(3) * t461;
t534 = qJD(6) * t200;
t179 = t461 * qJD(6);
t374 = t378 ^ 2;
t639 = t382 ^ 2;
t636 = t277 / 0.2e1;
t635 = t376 / 0.2e1;
t634 = t378 / 0.2e1;
t633 = -t382 / 0.2e1;
t396 = -pkin(3) * t616 - pkin(4);
t384 = -t403 * pkin(5) + t396;
t631 = t384 / 0.2e1;
t76 = (t636 + t638) * t382 + t644;
t630 = t76 * qJD(5);
t627 = t405 * pkin(3);
t624 = qJD(3) * pkin(3);
t623 = t67 * t571;
t621 = t68 * t276;
t417 = t442 * t616;
t359 = t405 * t417;
t284 = -t406 * t427 - t359;
t289 = pkin(4) * t380 + qJ(5) * t376 + t627;
t162 = -t284 * t401 + t403 * t289;
t110 = pkin(5) * t380 + t376 * t628 + t162;
t502 = t629 * t110;
t163 = t403 * t284 + t401 * t289;
t285 = t401 * t376;
t133 = pkin(8) * t285 + t163;
t568 = t404 * t133;
t83 = t502 - t568;
t501 = t629 * t133;
t569 = t404 * t110;
t84 = t501 + t569;
t9 = -t200 * t84 + t276 * t67 - t277 * t83 + t571 * t68;
t619 = t9 * qJD(1);
t590 = t276 * t382;
t595 = t571 * t378;
t115 = t590 / 0.2e1 - t595 / 0.2e1;
t618 = t115 * qJD(3);
t283 = t402 * t426 + t359;
t202 = pkin(5) * t288 + t283;
t31 = -t200 * t202 + t376 * t67;
t614 = qJD(1) * t31;
t32 = t202 * t277 - t376 * t68;
t613 = qJD(1) * t32;
t612 = qJD(1) * t81;
t455 = -t144 * t403 - t145 * t401;
t90 = t455 * t380;
t611 = qJD(1) * t90;
t581 = t376 * t380;
t588 = t277 * t276;
t594 = t200 * t571;
t91 = t588 / 0.2e1 - t594 / 0.2e1 + t581 / 0.2e1;
t610 = qJD(1) * t91;
t96 = -t200 * t636 + t575 * t637;
t609 = qJD(1) * t96;
t589 = t277 * t571;
t593 = t276 * t200;
t98 = -t589 - t593;
t608 = qJD(1) * t98;
t606 = t144 * t401;
t605 = t145 * t403;
t360 = t406 * t417;
t282 = t360 - t420;
t503 = pkin(5) * t285;
t201 = t282 - t503;
t16 = t200 * t201 - t202 * t571 + t376 * t83 - t380 * t67;
t604 = t16 * qJD(1);
t603 = t162 * t401;
t602 = t162 * t403;
t601 = t163 * t401;
t600 = t163 * t403;
t17 = t201 * t277 + t202 * t276 - t376 * t84 - t380 * t68;
t599 = t17 * qJD(1);
t598 = t202 * t380;
t26 = (t601 + t602) * t380 + t455 * t376;
t597 = t26 * qJD(1);
t596 = t571 * t280;
t592 = t276 * t281;
t591 = t276 * t376;
t586 = t277 * t380;
t585 = t282 * t403;
t584 = t283 * t282;
t583 = t283 * t380;
t582 = t376 * t571;
t580 = t376 * t393;
t579 = t376 * t402;
t578 = t378 * t376;
t38 = (t282 * t401 + t144) * t380 + (-t283 * t401 + t162) * t376;
t577 = t38 * qJD(1);
t576 = t380 * t200;
t574 = t380 * t384;
t573 = t380 * t396;
t416 = -t478 * t580 + t573 / 0.2e1;
t435 = -t602 / 0.2e1 - t601 / 0.2e1;
t64 = t416 + t435;
t565 = t64 * qJD(1);
t408 = t409 * t376;
t88 = t284 * t380 / 0.2e1 - t408 / 0.2e1 + t282 * t635 + t583 / 0.2e1;
t564 = t88 * qJD(1);
t99 = t589 - t593;
t563 = t99 * qJD(1);
t467 = t499 / 0.2e1;
t562 = -t341 / 0.2e1 + t376 * t467;
t468 = -t499 / 0.2e1;
t561 = t341 / 0.2e1 + t376 * t468;
t386 = t399 + t400;
t433 = -t200 * t633 + t575 * t634;
t102 = t433 + t372;
t558 = qJD(1) * t102;
t119 = -t576 + t582;
t557 = qJD(1) * t119;
t120 = t576 + t582;
t556 = qJD(1) * t120;
t121 = t586 + t591;
t555 = qJD(1) * t121;
t122 = t586 - t591;
t554 = qJD(1) * t122;
t127 = -t408 + t583;
t553 = qJD(1) * t127;
t262 = t376 * t627 + t380 * t385;
t552 = qJD(1) * t262;
t263 = -t376 * t385 + t380 * t627;
t551 = qJD(1) * t263;
t291 = t642 * t403;
t548 = qJD(1) * t291;
t547 = qJD(1) * t406;
t185 = t482 + t560;
t545 = qJD(3) * t185;
t484 = t578 / 0.2e1;
t188 = t484 + t561;
t543 = qJD(3) * t188;
t542 = qJD(3) * t282;
t539 = qJD(3) * t384;
t538 = qJD(3) * t403;
t536 = qJD(4) * t380;
t535 = qJD(5) * t376;
t533 = qJD(6) * t382;
t137 = (0.1e1 / 0.2e1 - t478) * t581;
t531 = t137 * qJD(1);
t529 = t185 * qJD(1);
t176 = t461 * qJD(1);
t528 = t188 * qJD(1);
t485 = -t578 / 0.2e1;
t189 = t485 + t561;
t527 = t189 * qJD(1);
t190 = t484 + t562;
t526 = t190 * qJD(3);
t191 = t485 + t562;
t525 = t191 * qJD(1);
t193 = t386 * t375;
t524 = t193 * qJD(1);
t523 = t575 * qJD(1);
t521 = t254 * qJD(1);
t256 = t504 * t401;
t520 = t256 * qJD(1);
t257 = t642 * t401;
t519 = t257 * qJD(1);
t258 = t504 * t403;
t518 = t258 * qJD(1);
t476 = t616 * t380;
t423 = -t579 / 0.2e1 - t476 / 0.2e1;
t261 = (-t405 / 0.2e1 + t423) * pkin(3);
t517 = t261 * qJD(1);
t516 = t504 * qJD(1);
t268 = t277 * qJD(6);
t515 = t288 * qJD(1);
t361 = t399 * t376;
t362 = t400 * t376;
t290 = t361 + t362;
t514 = t290 * qJD(1);
t513 = t290 * qJD(3);
t512 = t642 * qJD(1);
t511 = t372 * qJD(1);
t369 = t376 * qJD(3);
t370 = t378 * qJD(6);
t509 = t380 * qJD(1);
t371 = t380 * qJD(3);
t508 = t386 * qJD(3);
t391 = -t405 ^ 2 + t406 ^ 2;
t507 = t391 * qJD(1);
t506 = t405 * qJD(3);
t505 = t406 * qJD(3);
t498 = t575 * t510;
t495 = t277 * t510;
t494 = t378 * t371;
t492 = t401 * t538;
t491 = t380 * t535;
t316 = t376 * t509;
t315 = t376 * t371;
t490 = t382 * t370;
t489 = t397 * t405 * qJD(1);
t488 = t397 * t547;
t487 = t403 * t509;
t358 = t403 * t371;
t486 = t405 * t505;
t481 = t567 / 0.2e1;
t474 = qJD(6) * t372 + t316;
t473 = qJD(5) + t539;
t472 = t376 * t487;
t471 = t401 * t316;
t10 = t201 * t202 - t67 * t83 + t68 * t84;
t4 = t84 * t636 + t621 / 0.2e1 + t83 * t637 - t623 / 0.2e1 + t201 * t635 + t598 / 0.2e1;
t463 = t10 * qJD(1) + t4 * qJD(2);
t462 = t403 * t471;
t20 = t598 + t621 - t623;
t460 = t20 * qJD(1) + t91 * qJD(2);
t23 = -t200 * t68 - t575 * t67;
t459 = -qJD(1) * t23 - qJD(2) * t96;
t436 = t606 / 0.2e1 - t605 / 0.2e1;
t24 = (t600 / 0.2e1 - t603 / 0.2e1 + t283 / 0.2e1) * t380 + (t282 / 0.2e1 + t436) * t376;
t25 = t144 * t162 + t145 * t163 + t584;
t458 = t25 * qJD(1) + t24 * qJD(2);
t93 = t284 * t409 + t385 * t627 + t584;
t457 = t93 * qJD(1) + t88 * qJD(2);
t97 = -t277 * t575 + t646;
t456 = qJD(1) * t97 + qJD(3) * t76;
t454 = t600 - t603;
t453 = -t376 * t396 - t380 * t393;
t413 = -t596 / 0.2e1 + t592 / 0.2e1 + t574 / 0.2e1;
t438 = t633 * t84 + t634 * t83;
t19 = t413 + t438;
t452 = -qJD(1) * t19 + qJD(2) * t115;
t39 = (-t145 + t585) * t380 + (-t283 * t403 - t163) * t376;
t451 = t39 * qJD(1);
t53 = t583 + (-t605 + t606) * t376;
t450 = qJD(1) * t53 + qJD(2) * t137;
t138 = (0.1e1 - t386) * t581;
t449 = t24 * qJD(1) + t138 * qJD(2);
t89 = (t282 - t409) * t380 + (-t284 - t283) * t376;
t448 = t89 * qJD(1);
t186 = -t571 / 0.2e1 - t560;
t414 = t202 * t632 - t281 * t376 / 0.2e1 + t277 * t631;
t424 = -t568 / 0.2e1 + t502 / 0.2e1;
t27 = -t414 + t424;
t447 = t27 * qJD(1) + t186 * qJD(2);
t415 = -t200 * t631 + t202 * t652 + t280 * t635;
t425 = -t569 / 0.2e1 - t501 / 0.2e1;
t28 = -t415 + t425;
t446 = t28 * qJD(1) + t190 * qJD(2);
t107 = -t277 ^ 2 + t646;
t95 = -t382 * t277 + t644;
t445 = qJD(1) * t107 + qJD(3) * t95;
t265 = t374 - t639;
t444 = qJD(1) * t95 + qJD(3) * t265;
t297 = t374 + t639;
t443 = qJD(1) * t76 + qJD(3) * t297;
t148 = t342 / 0.2e1 + (t481 - t499) * t380;
t441 = qJD(1) * t148 - t540;
t437 = qJD(3) * t466;
t432 = t587 / 0.2e1 - t483;
t92 = t581 + t588 - t594;
t422 = t4 * qJD(1) + t92 * qJD(2) + t115 * qJD(4);
t103 = t432 + t372;
t123 = t280 * t382 - t281 * t378;
t412 = t360 / 0.2e1 - t420 / 0.2e1;
t410 = -t503 / 0.2e1 + t412;
t14 = t410 + t653;
t419 = qJD(1) * t14 + qJD(2) * t103 - qJD(3) * t123;
t322 = t386 * t393;
t73 = t412 + t436;
t418 = qJD(1) * t73 - qJD(2) * t255 - qJD(3) * t322;
t392 = t405 * t547;
t363 = t372 * qJD(3);
t357 = t401 * t371;
t317 = t382 * t371;
t260 = t627 / 0.2e1 + t423 * pkin(3);
t180 = t185 * qJD(6);
t178 = t188 * qJD(6);
t177 = t190 * qJD(6);
t151 = -t342 / 0.2e1 + (t467 + t481 + t468) * t380;
t141 = -t176 - t533;
t105 = -t432 + t372;
t104 = -t433 + t372;
t94 = t95 * qJD(6);
t74 = t412 - t436;
t63 = t416 - t435;
t46 = qJD(3) * t88;
t30 = t414 + t424;
t29 = t415 + t425;
t21 = qJD(3) * t24 + qJD(4) * t137;
t18 = t413 - t438;
t15 = t410 - t653;
t1 = qJD(3) * t4 + qJD(4) * t91 + qJD(5) * t96;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t486, t391 * qJD(3), 0, -t486, 0, 0, t397 * t506, t397 * t505, 0, 0, -t315, -t504 * qJD(3), 0, t315, 0, 0, t262 * qJD(3), t263 * qJD(3), qJD(3) * t89 - qJD(4) * t642, qJD(3) * t93 + qJD(4) * t127, -t400 * t315, 0.2e1 * t358 * t285, t258 * qJD(3), -t399 * t315, -t256 * qJD(3), t315, qJD(3) * t38 - qJD(4) * t257 - t403 * t491, qJD(3) * t39 - qJD(4) * t291 + t401 * t491, -qJD(3) * t26 + qJD(5) * t193, qJD(3) * t25 + qJD(4) * t53 + qJD(5) * t90 (qJD(3) * t276 - t534) * t277, qJD(3) * t99 + qJD(6) * t107, qJD(3) * t121 - t376 * t534 (-qJD(3) * t571 + t268) * t200, qJD(3) * t119 - t268 * t376, t315, qJD(3) * t16 + qJD(4) * t120 + qJD(6) * t32 + t535 * t575, qJD(3) * t17 + qJD(4) * t122 + qJD(6) * t31 + t200 * t535, qJD(3) * t9 + qJD(4) * t98 + qJD(5) * t97, qJD(3) * t10 + qJD(4) * t20 + qJD(5) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t392, t507, t505, -t392, -t506, 0, -t406 * t437 + t489, t405 * t437 + t488, 0, 0, -t316, -t516, -t369, t316, -t371, 0, -t542 + t552, -qJD(3) * t284 + t551 (t376 * t616 - t380 * t402) * t624 + t448 (-t282 * t616 + t284 * t402) * t624 + t260 * qJD(4) + t457 (-t400 * t509 - t492) * t376, 0.2e1 * t462 + (t361 - t362) * qJD(3), t357 + t518 (-t399 * t509 + t492) * t376, t358 - t520, t316, t577 + (t401 * t453 - t585) * qJD(3) - t285 * qJD(5), t401 * t542 + (qJD(3) * t453 - t535) * t403 + t451, qJD(3) * t454 - t597 (t282 * t396 + t393 * t454) * qJD(3) + t63 * qJD(4) + t74 * qJD(5) + t458, t650 + (t540 + t549) * t276, t563 + (-t237 + t572) * qJD(3) + t94, -t177 + t317 + t555, -t571 * t640 - t650, -t180 - t494 + t557, t474, t604 + (t201 * t378 - t280 * t380 - t384 * t571) * qJD(3) - t461 * qJD(5) + t30 * qJD(6), t599 + (t201 * t382 + t276 * t384 - t281 * t380) * qJD(3) - t191 * qJD(5) + t29 * qJD(6), t619 + (t276 * t280 + t281 * t571 - t378 * t84 - t382 * t83) * qJD(3) + t630 (t201 * t384 - t280 * t83 + t281 * t84) * qJD(3) + t18 * qJD(4) + t15 * qJD(5) + t463; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t512, qJD(3) * t260 + t553, 0, 0, 0, 0, 0, 0, -t519, -t548, 0, qJD(3) * t63 + t450 - t655, 0, 0, 0, 0, 0, 0, -t180 + t556, -qJD(6) * t189 + t554, t608, t18 * qJD(3) + (t590 - t595) * qJD(4) + t104 * qJD(5) + t460; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t285 - t472 (t401 * t509 - t538) * t376, t524, qJD(3) * t74 - qJD(4) * t255 + t611, 0, 0, 0, 0, 0, 0, qJD(6) * t151 + t498 - t544, -qJD(3) * t191 + t496, t456, qJD(3) * t15 + qJD(4) * t104 - t459; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t647, t445, -t526 + (-qJD(6) - t510) * t200, -t647, -t268 - t495 - t545, t363, qJD(3) * t30 - qJD(4) * t185 + qJD(5) * t151 - qJD(6) * t68 + t613, qJD(3) * t29 - qJD(4) * t189 + qJD(6) * t67 + t614, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138 * qJD(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t506, -t505, 0, 0, 0, 0, 0, 0, 0, 0, -t371, t369, 0, t564 + (-t476 - t579) * t624, 0, 0, 0, 0, 0, 0, -t358, t357, -t513 (-t386 * t580 + t573) * qJD(3) + t654 + t449, 0, 0, 0, 0, 0, 0, t494 + t179, -t178 + t317, t645 (t574 + t592 - t596) * qJD(3) + t105 * qJD(5) + t422; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t531, 0, 0, 0, 0, 0, 0, 0, 0, 0, t610 + t618; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t254 * qJD(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t105 + t609; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t268 + t544, t534 - t543, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t392, -t507, 0, t392, 0, 0, -t489, -t488, 0, 0, t316, t516, 0, -t316, 0, 0, -t536 - t552, qJD(4) * t376 - t551, -t448, qJD(4) * t261 - t457, t400 * t316, -0.2e1 * t462, -t518, t399 * t316, t520, -t316, -t403 * t536 - t577, qJD(4) * t288 - t451, -qJD(4) * t290 + t597, qJD(4) * t64 - qJD(5) * t73 - t458, -t276 * t549 + t650, t94 - t563, -t178 - t555, t571 * t649 - t650, -t179 - t557, -t474, qJD(4) * t575 - qJD(5) * t185 - qJD(6) * t27 - t604, qJD(4) * t200 - qJD(5) * t189 - qJD(6) * t28 - t599, qJD(4) * t81 - t619 + t630, qJD(4) * t19 - qJD(5) * t14 - t463; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t564, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t449 + t655, 0, 0, 0, 0, 0, 0, -qJD(6) * t186, -t177, 0, -qJD(5) * t103 - t422; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t386 * qJD(5), t322 * qJD(5), -t490, t265 * qJD(6), 0, t490, 0, 0, t384 * t533, -t384 * t370, qJD(5) * t297, qJD(5) * t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t509, t510, 0, t517, 0, 0, 0, 0, 0, 0, -t487, t515, -t514, t565, 0, 0, 0, 0, 0, 0, t523, t649, t612, -t452; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t508, -t418, 0, 0, 0, 0, 0, 0, -t529, -t527, t443, -t419; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t648, t444, -t370 - t528, t648, t141, -t511, -qJD(6) * t281 + t382 * t539 - t447, qJD(6) * t280 - t378 * t539 - t446, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t371, -t369, t512, -qJD(3) * t261 - t553, 0, 0, 0, 0, 0, 0, t358 + t519, -qJD(3) * t288 + t548, t513, -qJD(3) * t64 - t450 - t654, 0, 0, 0, 0, 0, 0, -qJD(3) * t575 - t179 - t556, -qJD(3) * t200 - qJD(6) * t191 - t554, -t608 - t645, -qJD(3) * t19 - qJD(5) * t102 - t460; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t531, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t610 + t618; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t509, -t510, 0, -t517, 0, 0, 0, 0, 0, 0, t487, -t515, t514, -t565, 0, 0, 0, 0, 0, 0, -t523, -t649, -t612, t452; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t521, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t558; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, t370 - t525, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t472, -t471, -t524, qJD(3) * t73 + qJD(4) * t254 - t611, 0, 0, 0, 0, 0, 0, -qJD(6) * t148 - t498 + t545, qJD(3) * t189 - t496 - t534, -t456, qJD(3) * t14 + qJD(4) * t102 + t459; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t255 * qJD(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t103 - t609; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t508, t418, 0, 0, 0, 0, 0, 0, t529 + t533, -t370 + t527, -t443, t419; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t521, 0, 0, 0, 0, 0, 0, 0, 0, 0, t558; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t441, -t640, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t647, -t445, t543 + t496, t647, t495 + t544, t363, qJD(3) * t27 + qJD(4) * t461 + qJD(5) * t148 - t613, qJD(3) * t28 + qJD(4) * t191 + qJD(5) * t200 - t614, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t186, t526, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t648, -t444, t528, -t648, t176, t511, -t382 * t473 + t447, t378 * t473 + t446, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, t525, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t441, t640, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t2;
