% Calculate inertial parameters regressor of coriolis matrix for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6PRPRPR2_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:33:01
% EndTime: 2019-03-08 19:33:16
% DurationCPUTime: 10.73s
% Computational Cost: add. (10999->508), mult. (26980->777), div. (0->0), fcn. (30920->12), ass. (0->389)
t433 = cos(pkin(12));
t432 = sin(pkin(11));
t421 = pkin(2) * t432 + pkin(8);
t437 = cos(qJ(4));
t587 = t437 * t421;
t395 = t433 * t587;
t431 = sin(pkin(12));
t632 = cos(pkin(11));
t422 = -pkin(2) * t632 - pkin(3);
t637 = t437 * pkin(4);
t471 = t422 - t637;
t435 = sin(qJ(4));
t589 = t435 * qJ(5);
t454 = t471 - t589;
t325 = t431 * t454 + t395;
t599 = t431 * t435;
t272 = -pkin(9) * t599 + t325;
t434 = sin(qJ(6));
t636 = pkin(9) + qJ(5);
t439 = (-t421 * t431 - pkin(5)) * t437 + (-t435 * t636 + t471) * t433;
t640 = cos(qJ(6));
t162 = t434 * t272 - t439 * t640;
t163 = t272 * t640 + t434 * t439;
t404 = t636 * t433;
t523 = t636 * t431;
t331 = t434 * t404 + t523 * t640;
t332 = t404 * t640 - t434 * t523;
t542 = t640 * t431;
t590 = t434 * t433;
t390 = t542 + t590;
t365 = t390 * t435;
t648 = t390 / 0.2e1;
t541 = t640 * t433;
t591 = t434 * t431;
t663 = t541 - t591;
t650 = t663 / 0.2e1;
t368 = t663 * t435;
t653 = -t368 / 0.2e1;
t667 = t332 * t365 / 0.2e1 + t331 * t653 - t162 * t648 - t163 * t650;
t656 = -t365 / 0.2e1;
t652 = t368 / 0.2e1;
t633 = cos(pkin(6));
t520 = t633 * t435;
t436 = sin(qJ(2));
t631 = sin(pkin(6));
t510 = t632 * t631;
t641 = cos(qJ(2));
t511 = t631 * t641;
t363 = t432 * t511 + t436 * t510;
t610 = t363 * t437;
t328 = t520 + t610;
t521 = t436 * t631;
t362 = t432 * t521 - t510 * t641;
t612 = t362 * t433;
t500 = -t328 * t431 + t612;
t666 = -t500 / 0.2e1;
t466 = t590 / 0.2e1 + t542 / 0.2e1;
t665 = t437 * t466;
t427 = t431 ^ 2;
t428 = t433 ^ 2;
t410 = t428 + t427;
t643 = -t437 / 0.2e1;
t527 = t390 * t643;
t276 = t527 + t665;
t576 = qJD(2) * t437;
t354 = t368 * t576;
t361 = t368 * qJD(6);
t664 = qJD(4) * t276 - t354 + t361;
t613 = t362 * t431;
t619 = t328 * t433;
t214 = t613 + t619;
t123 = t434 * t214 - t500 * t640;
t124 = t214 * t640 + t434 * t500;
t634 = t123 * t652 + t124 * t656;
t364 = t365 ^ 2;
t662 = t368 ^ 2;
t661 = t390 ^ 2;
t660 = t124 / 0.2e1;
t327 = t363 * t435 - t437 * t633;
t659 = -t327 / 0.2e1;
t658 = t327 / 0.2e1;
t657 = t328 / 0.2e1;
t367 = t390 * t437;
t655 = -t367 / 0.2e1;
t654 = t367 / 0.2e1;
t408 = t437 * t591;
t370 = t437 * t541 - t408;
t651 = t370 / 0.2e1;
t649 = -t390 / 0.2e1;
t423 = -pkin(5) * t433 - pkin(4);
t647 = t423 / 0.2e1;
t646 = -t431 / 0.2e1;
t645 = -t435 / 0.2e1;
t644 = t435 / 0.2e1;
t642 = t437 / 0.2e1;
t639 = t431 * pkin(5);
t638 = t435 * pkin(4);
t602 = t390 * t368;
t609 = t365 * t663;
t200 = t602 - t609;
t635 = t200 * qJD(5);
t598 = t431 * t437;
t231 = t362 * t598 + t363 * t433;
t546 = t640 * t231;
t596 = t433 * t437;
t232 = -t362 * t596 + t363 * t431;
t594 = t434 * t232;
t144 = t546 - t594;
t545 = t640 * t232;
t595 = t434 * t231;
t145 = t545 + t595;
t246 = t362 * t437;
t530 = t246 / 0.2e1;
t48 = t144 * t656 + t145 * t652 + t435 * t530;
t630 = qJD(1) * t48;
t626 = t214 * t433;
t324 = -t431 * t587 + t433 * t454;
t625 = t324 * t431;
t624 = t324 * t437;
t623 = t325 * t433;
t622 = t325 * t437;
t621 = t327 * t435;
t620 = t328 * t421;
t618 = t328 * t437;
t394 = t421 * t599;
t405 = -qJ(5) * t437 + t638;
t338 = t433 * t405 + t394;
t616 = t338 * t435;
t597 = t433 * t435;
t339 = t431 * t405 - t421 * t597;
t615 = t339 * t435;
t614 = t362 * t421;
t611 = t362 * t435;
t608 = t367 * t390;
t606 = t368 * t663;
t605 = t368 * t435;
t522 = t421 + t639;
t378 = t522 * t435;
t604 = t378 * t435;
t603 = t390 * t365;
t601 = t427 * t437;
t429 = t435 ^ 2;
t600 = t429 * t431;
t275 = t435 * pkin(5) - pkin(9) * t596 + t338;
t593 = t434 * t275;
t318 = -pkin(9) * t598 + t339;
t592 = t434 * t318;
t588 = t435 * t365;
t425 = t435 * t437;
t88 = (t363 - t618 - t621) * t362;
t585 = t88 * qJD(1);
t430 = t437 ^ 2;
t415 = t430 - t429;
t480 = t232 * t433 / 0.2e1 + t231 * t646;
t103 = (t530 + t480) * t435;
t583 = qJD(1) * t103;
t251 = -t367 * t437 + t588;
t580 = qJD(2) * t251;
t252 = t370 * t437 - t605;
t579 = qJD(2) * t252;
t578 = qJD(2) * t365;
t577 = qJD(2) * t368;
t574 = qJD(4) * t663;
t573 = qJD(4) * t390;
t572 = qJD(4) * t423;
t571 = qJD(4) * t431;
t570 = qJD(4) * t433;
t569 = qJD(5) * t437;
t568 = qJD(6) * t390;
t172 = -t324 * t435 + (t338 - 0.2e1 * t394) * t437;
t567 = t172 * qJD(2);
t184 = -t370 * t365 - t368 * t367;
t566 = t184 * qJD(2);
t565 = t276 * qJD(2);
t278 = (t649 - t466) * t437;
t564 = t278 * qJD(2);
t512 = t541 / 0.2e1;
t279 = -t663 * t643 + t437 * t512 - t408 / 0.2e1;
t563 = t279 * qJD(2);
t562 = t279 * qJD(4);
t280 = t408 / 0.2e1 + (t650 - t541 / 0.2e1) * t437;
t561 = t280 * qJD(2);
t560 = t280 * qJD(4);
t559 = t362 * qJD(2);
t558 = t363 * qJD(2);
t557 = t365 * qJD(6);
t385 = t410 * t429;
t556 = t385 * qJD(2);
t381 = t663 * qJD(6);
t388 = -t430 * t431 + t600;
t555 = t388 * qJD(2);
t389 = t415 * t433;
t554 = t389 * qJD(2);
t553 = t410 * qJD(4);
t552 = t415 * qJD(2);
t551 = t435 * qJD(2);
t550 = t435 * qJD(4);
t549 = t437 * qJD(4);
t544 = t640 * t275;
t177 = t544 - t592;
t543 = t640 * t318;
t178 = t543 + t593;
t379 = t522 * t437;
t25 = t178 * t652 + t163 * t651 + t177 * t656 + t162 * t654 + t379 * t643 + t604 / 0.2e1;
t548 = t25 * qJD(4);
t547 = t327 * t611;
t540 = t365 * t577;
t539 = t663 * t573;
t538 = t663 * t550;
t537 = t421 * t549;
t536 = t431 * t570;
t535 = t435 * t569;
t534 = t390 * t381;
t533 = t422 * t576;
t532 = t433 * t550;
t531 = t435 * t549;
t416 = t437 * t551;
t344 = -t611 / 0.2e1;
t345 = t611 / 0.2e1;
t529 = -t606 / 0.2e1;
t528 = -t603 / 0.2e1;
t526 = t597 / 0.2e1;
t524 = t587 / 0.2e1;
t519 = t410 * t327;
t518 = t410 * t437;
t516 = qJD(4) * t278 - t354;
t515 = qJD(5) + t572;
t514 = t433 * t416;
t30 = t162 * t370 - t163 * t367 - t177 * t368 - t178 * t365;
t185 = t327 * t390;
t186 = t327 * t663;
t440 = t123 * t651 + t124 * t655 + t185 * t653 - t186 * t656;
t482 = t144 * t649 + t145 * t650;
t7 = t440 - t482;
t509 = t7 * qJD(1) + t30 * qJD(2);
t508 = t431 * t514;
t452 = -t186 * t643 + t328 * t653 + t370 * t659;
t35 = (t362 * t649 + t660) * t435 + t452;
t53 = -t163 * t435 + t178 * t437 + t379 * t368 + t378 * t370;
t507 = -t35 * qJD(1) + t53 * qJD(2);
t40 = t344 - t634;
t57 = t162 * t368 - t163 * t365;
t506 = qJD(1) * t40 - qJD(2) * t57;
t444 = t214 * t646 + t433 * t666;
t441 = t444 * t437;
t56 = t441 - t480;
t89 = (t616 + t624) * t433 + (t615 + t622) * t431;
t505 = -t56 * qJD(1) + t89 * qJD(2);
t17 = -t186 * t652 + t124 * t651 + t185 * t656 + t123 * t654 - t618 / 0.2e1 + t621 / 0.2e1;
t26 = -t123 * t185 - t124 * t186 + t327 * t328;
t504 = t26 * qJD(1) + t17 * qJD(3);
t20 = -t123 * t144 + t124 * t145 - t547;
t503 = t20 * qJD(1) + t48 * qJD(3);
t372 = (0.1e1 / 0.2e1 - t428 / 0.2e1 - t427 / 0.2e1) * t435;
t472 = t500 * t431;
t443 = t626 / 0.2e1 - t472 / 0.2e1;
t50 = (-t328 / 0.2e1 + t443) * t437 + t327 * t372;
t54 = (t328 + t472 - t626) * t327;
t502 = t54 * qJD(1) + t50 * qJD(3);
t499 = -t338 * t431 + t339 * t433;
t115 = -t162 * t437 - t378 * t365;
t470 = -t595 / 0.2e1 - t545 / 0.2e1;
t484 = t123 * t643 + t327 * t656;
t45 = t470 - t484;
t498 = qJD(1) * t45 - qJD(2) * t115;
t116 = -t163 * t437 - t378 * t368;
t469 = -t594 / 0.2e1 + t546 / 0.2e1;
t483 = t124 * t642 + t327 * t652;
t44 = t469 - t483;
t497 = qJD(1) * t44 + qJD(2) * t116;
t173 = t339 * t437 + (-t325 + 0.2e1 * t395) * t435;
t64 = (-t613 / 0.2e1 + t214 / 0.2e1 - t619 / 0.2e1) * t435;
t496 = t64 * qJD(1) - t173 * qJD(2);
t182 = (t324 * t433 + t325 * t431) * t435;
t442 = t444 * t435;
t83 = t345 + t442;
t495 = -qJD(1) * t83 + qJD(2) * t182;
t51 = t214 * t232 + t231 * t500 - t547;
t494 = t51 * qJD(1) + t103 * qJD(3);
t153 = -t602 - t609;
t213 = t364 - t662;
t493 = qJD(2) * t213 + qJD(4) * t153;
t382 = t663 ^ 2;
t255 = t382 - t661;
t492 = qJD(2) * t153 + qJD(4) * t255;
t263 = t364 + t662;
t491 = qJD(2) * t263 + qJD(4) * t200;
t329 = t382 + t661;
t490 = qJD(2) * t200 + qJD(4) * t329;
t489 = t574 - t578;
t488 = t573 + t577;
t487 = (t421 / 0.2e1 + t639 / 0.2e1) * t437;
t485 = t123 * t648 + t124 * t650;
t479 = t623 / 0.2e1 - t625 / 0.2e1;
t243 = t345 + t344;
t477 = -t243 * qJD(1) - t422 * t551;
t311 = t603 / 0.2e1;
t196 = t529 + t311;
t476 = -qJD(4) * t196 - t540;
t475 = qJD(2) * t196 - t539;
t312 = t606 / 0.2e1;
t197 = t528 + t312;
t474 = qJD(4) * t197 - t540;
t473 = -qJD(2) * t197 - t539;
t468 = -t593 / 0.2e1 - t543 / 0.2e1;
t467 = -t592 / 0.2e1 + t544 / 0.2e1;
t465 = t512 - t591 / 0.2e1;
t438 = -t123 * t177 / 0.2e1 + t178 * t660 - t185 * t162 / 0.2e1 - t186 * t163 / 0.2e1 + t379 * t658 + t378 * t657;
t448 = t144 * t331 / 0.2e1 - t145 * t332 / 0.2e1 + t423 * t345;
t3 = t438 + t448;
t32 = -t162 * t177 + t163 * t178 + t378 * t379;
t464 = t3 * qJD(1) + t32 * qJD(2) + t25 * qJD(3);
t463 = t610 / 0.2e1 + t520 / 0.2e1;
t462 = t480 * qJ(5);
t114 = t421 ^ 2 * t425 + t324 * t338 + t325 * t339;
t168 = t524 - t479;
t445 = -t214 * t339 / 0.2e1 + t338 * t666;
t18 = pkin(4) * t345 - t168 * t327 + t620 * t645 + t445 + t462;
t76 = (t615 / 0.2e1 + t622 / 0.2e1) * t433 + (-t616 / 0.2e1 - t624 / 0.2e1) * t431 + (t429 / 0.2e1 - t430 / 0.2e1) * t421;
t461 = -t18 * qJD(1) + t114 * qJD(2) + t76 * qJD(3);
t453 = t185 * t642 + t327 * t655 + t328 * t656;
t34 = (t362 * t650 + t123 / 0.2e1) * t435 + t453;
t52 = t162 * t435 + t177 * t437 - t379 * t365 - t378 * t367;
t460 = -t34 * qJD(1) - t52 * qJD(2);
t170 = t365 * t367 + t368 * t370 - t425;
t459 = t17 * qJD(1) + t25 * qJD(2) + t170 * qJD(3);
t341 = t435 * t518 - t425;
t458 = t50 * qJD(1) + t76 * qJD(2) + t341 * qJD(3);
t455 = t649 + t466;
t132 = t455 * t327;
t277 = t455 * t437;
t450 = t332 * t642 + t368 * t647 + t378 * t648;
t70 = -t450 + t467;
t457 = t132 * qJD(1) + t70 * qJD(2) - t277 * qJD(3);
t133 = (-t663 / 0.2e1 + t465) * t327;
t451 = t331 * t643 - t365 * t647 + t378 * t650;
t71 = -t451 + t468;
t456 = t133 * qJD(1) + t71 * qJD(2) + t280 * qJD(3);
t449 = (-t589 - t637) * qJD(4) + t569;
t167 = t331 * t390 + t332 * t663;
t187 = t644 + t529 + t528;
t38 = t487 + t667;
t42 = t463 - t485;
t447 = qJD(1) * t42 + qJD(2) * t38 + qJD(3) * t187 - qJD(4) * t167;
t399 = t410 * qJ(5);
t80 = -t443 + t463;
t446 = qJD(1) * t80 + qJD(2) * t168 + qJD(3) * t372 - qJD(4) * t399;
t424 = t550 / 0.2e1;
t417 = t428 * t437;
t414 = t431 * t550;
t383 = qJD(6) * t645 + t416;
t377 = t390 * t550;
t373 = (0.1e1 + t410) * t644;
t330 = t429 * t614;
t321 = t370 * t663;
t281 = t527 - t665;
t262 = t279 * qJD(6);
t261 = t280 * qJD(6);
t244 = 0.2e1 * t345;
t207 = -t365 * t576 + t562;
t192 = t197 * qJD(6);
t191 = t196 * qJD(6);
t188 = t312 + t311 + t644;
t169 = t524 + t479;
t166 = -t560 + (-qJD(6) + t576) * t365;
t152 = t153 * qJD(6);
t135 = (t466 + t648) * t327;
t134 = (t465 + t650) * t327;
t82 = t344 + t442;
t81 = t443 + t463;
t75 = t76 * qJD(4);
t73 = t450 + t467;
t72 = t451 + t468;
t67 = t214 * t645 + t328 * t526 + t344 * t431;
t66 = t362 * t526 + t500 * t644 + t599 * t657;
t55 = t441 + t480;
t47 = t469 + t483;
t46 = t470 + t484;
t43 = t463 + t485;
t41 = t344 + t634;
t39 = t487 - t667;
t37 = t124 * t645 + t344 * t390 - t452;
t36 = t123 * t645 - t344 * t663 - t453;
t31 = qJD(2) * t103 + qJD(4) * t50;
t19 = t623 * t659 + t625 * t658 + t327 * t524 + (t620 / 0.2e1 + t362 * pkin(4) / 0.2e1) * t435 + t462 - t445;
t6 = t440 + t482;
t2 = t438 - t448;
t1 = qJD(2) * t48 + qJD(4) * t17;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t88, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t51 + qJD(4) * t54, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t20 + qJD(4) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t521, -qJD(2) * t511, 0, 0, 0, 0, 0, 0, 0, 0, -t558, t559, 0 (-t362 * t432 - t363 * t632) * qJD(2) * pkin(2), 0, 0, 0, 0, 0, 0, t244 * qJD(4) - t437 * t558, qJD(4) * t246 + t363 * t551 (-t429 - t430) * t559, t585 + (t363 * t422 - t430 * t614 - t330) * qJD(2), 0, 0, 0, 0, 0, 0 (-t231 * t437 - t362 * t600) * qJD(2) + t66 * qJD(4) (t232 * t437 - t429 * t612) * qJD(2) + t67 * qJD(4), t55 * qJD(4) + (-t231 * t433 - t232 * t431) * t551 (t231 * t324 + t232 * t325 - t330) * qJD(2) + t19 * qJD(4) + t82 * qJD(5) + t494, 0, 0, 0, 0, 0, 0 (-t144 * t437 - t362 * t588) * qJD(2) + t36 * qJD(4) + t47 * qJD(6) (t145 * t437 - t362 * t605) * qJD(2) + t37 * qJD(4) + t46 * qJD(6) (-t144 * t368 - t145 * t365) * qJD(2) + t6 * qJD(4) (-t144 * t162 + t145 * t163 - t362 * t604) * qJD(2) + t2 * qJD(4) + t41 * qJD(5) + t503; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t244 - qJD(4) * t328, qJD(2) * t246 + qJD(4) * t327, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t66 - t328 * t570, qJD(2) * t67 + t328 * t571, t55 * qJD(2) - qJD(4) * t519, t19 * qJD(2) + (-t328 * pkin(4) - qJ(5) * t519) * qJD(4) + t81 * qJD(5) + t502, 0, 0, 0, 0, 0, 0, qJD(2) * t36 + qJD(6) * t135 - t328 * t574, qJD(2) * t37 + qJD(6) * t134 + t328 * t573, t6 * qJD(2) + (-t185 * t390 - t186 * t663) * qJD(4), t2 * qJD(2) + (-t185 * t331 - t186 * t332 + t328 * t423) * qJD(4) + t43 * qJD(5) + t504; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t82 + qJD(4) * t81, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t41 + qJD(4) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t47 + qJD(4) * t135 - qJD(6) * t124, qJD(2) * t46 + qJD(4) * t134 + qJD(6) * t123, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t243 * qJD(4), 0, 0, -t585, 0, 0, 0, 0, 0, 0, 0, -t64 * qJD(4), t56 * qJD(4), -qJD(4) * t18 + qJD(5) * t83 - t494, 0, 0, 0, 0, 0, 0, -qJD(4) * t34 - qJD(6) * t44, -qJD(4) * t35 - qJD(6) * t45, qJD(4) * t7, qJD(4) * t3 - qJD(5) * t40 - t503; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t531, t415 * qJD(4), 0, -t531, 0, 0, t422 * t550, t422 * t549, 0, 0, t428 * t531, -0.2e1 * t433 * t431 * t531, -t389 * qJD(4), t427 * t531, -t388 * qJD(4), -t531, -t172 * qJD(4) + t433 * t535, t173 * qJD(4) - t431 * t535, -qJD(4) * t89 + qJD(5) * t385, qJD(4) * t114 - qJD(5) * t182 (qJD(4) * t370 - t557) * t368, qJD(4) * t184 + qJD(6) * t213, -t252 * qJD(4) + t437 * t557 (qJD(4) * t367 + t361) * t365, -t251 * qJD(4) + t361 * t437, -t531, -t52 * qJD(4) - t116 * qJD(6) + t368 * t569, t53 * qJD(4) + t115 * qJD(6) - t365 * t569, qJD(4) * t30 + qJD(5) * t263, qJD(4) * t32 + qJD(5) * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75 - t583, 0, 0, 0, 0, 0, 0, 0, 0, 0, t548 - t630; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t416, t552, t549, -t416, -t550, 0, -t477 - t537, t421 * t550 + t533, 0, 0 (t428 * t551 + t536) * t437, -0.2e1 * t508 + (t417 - t601) * qJD(4), t414 - t554 (t427 * t551 - t536) * t437, t532 - t555, -t416, -t395 * qJD(4) + t431 * t449 - t567, t431 * t537 + t433 * t449 - t496, qJD(4) * t499 - t505 (-pkin(4) * t587 + qJ(5) * t499) * qJD(4) + t169 * qJD(5) + t461, t370 * t488 + t192, t566 + (t321 - t608) * qJD(4) + t152, -t261 + t377 - t579, -t367 * t489 + t191, -qJD(6) * t276 + t538 - t580, -t383 (-t331 * t435 + t367 * t423 - t379 * t663) * qJD(4) - t278 * qJD(5) + t73 * qJD(6) + t460 (-t332 * t435 + t370 * t423 + t379 * t390) * qJD(4) + t279 * qJD(5) + t72 * qJD(6) + t507 (-t177 * t390 + t178 * t663 + t331 * t370 - t332 * t367) * qJD(4) + t509 + t635 (-t177 * t331 + t178 * t332 + t379 * t423) * qJD(4) + t39 * qJD(5) + t464; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t433 * t551 + t571) * t437 (-t431 * t551 + t570) * t437, t556, qJD(4) * t169 - t495, 0, 0, 0, 0, 0, 0, -t516, t207, t491, qJD(4) * t39 - t506; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t474, t493, t166, -t476, -t664, t424, qJD(4) * t73 - qJD(6) * t163 - t497, qJD(4) * t72 + qJD(6) * t162 - t498, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75 + t583, 0, 0, 0, 0, 0, 0, 0, 0, 0, t548 + t630; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t341 * qJD(4), 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t550, -t549, 0, 0, 0, 0, 0, 0, 0, 0, -t532, t414 (t417 + t601) * qJD(4) (qJ(5) * t518 - t638) * qJD(4) + t373 * qJD(5) + t458, 0, 0, 0, 0, 0, 0, qJD(6) * t281 - t538, t377 - t262 (t321 + t608) * qJD(4) (t331 * t367 + t332 * t370 + t423 * t435) * qJD(4) + t188 * qJD(5) + t459; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t373 * qJD(4), 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t281 - t361, t557 - t562, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t243 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64 * qJD(2), -t56 * qJD(2), qJD(2) * t18 - qJD(5) * t80 - t502, 0, 0, 0, 0, 0, 0, qJD(2) * t34 - qJD(6) * t132, qJD(2) * t35 - qJD(6) * t133, -qJD(2) * t7, -qJD(2) * t3 - qJD(5) * t42 - t504; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t416, -t552, 0, t416, 0, 0, t477, -t533, 0, 0, -t428 * t416, 0.2e1 * t508, t554, -t427 * t416, t555, t416, t567, t496, t505, -qJD(5) * t168 - t461, -t370 * t577 + t192, t152 - t566, -t262 + t579, -t367 * t578 + t191, -qJD(6) * t278 + t580, t383, -qJD(5) * t276 - qJD(6) * t70 - t460, qJD(5) * t280 - qJD(6) * t71 - t507, -t509 + t635, -qJD(5) * t38 - t464; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(5) * t372 - t458, 0, 0, 0, 0, 0, 0, qJD(6) * t277, -t261, 0, -qJD(5) * t187 - t459; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t410 * qJD(5), t399 * qJD(5), t534, t255 * qJD(6), 0, -t534, 0, 0, t423 * t568, t423 * t381, qJD(5) * t329, qJD(5) * t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t553, -t446, 0, 0, 0, 0, 0, 0, -t565, t561, t490, -t447; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t473, t492, t381 - t563, t475, -t564 - t568, -t551 / 0.2e1, -qJD(6) * t332 + t390 * t572 - t457, qJD(6) * t331 + t572 * t663 - t456, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2) * t83 + qJD(4) * t80, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t40 + qJD(4) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t514, t431 * t416, -t556, qJD(4) * t168 + t495, 0, 0, 0, 0, 0, 0, t664, t166, -t491, qJD(4) * t38 + t506; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t372 * qJD(4), 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t187; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t553, t446, 0, 0, 0, 0, 0, 0, t565 + t568, t381 - t561, -t490, t447; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t488, t489, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t44 + qJD(4) * t132, qJD(2) * t45 + qJD(4) * t133, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t474, -t493, t207, t476, t516, t424, qJD(4) * t70 - qJD(5) * t368 + t497, qJD(4) * t71 + qJD(5) * t365 + t498, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4) * t277, t560, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t473, -t492, t563, -t475, t564, t551 / 0.2e1, -t390 * t515 + t457, -t515 * t663 + t456, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t488, -t489, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t4;
