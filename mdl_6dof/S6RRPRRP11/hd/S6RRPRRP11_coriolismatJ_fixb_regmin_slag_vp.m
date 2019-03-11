% Calculate minimal parameter regressor of coriolis matrix for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x30]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRPRRP11_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP11_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:48:40
% EndTime: 2019-03-09 12:49:03
% DurationCPUTime: 12.37s
% Computational Cost: add. (11373->598), mult. (21003->806), div. (0->0), fcn. (21370->6), ass. (0->459)
t420 = sin(qJ(4));
t422 = cos(qJ(5));
t596 = t422 * t420;
t419 = sin(qJ(5));
t423 = cos(qJ(4));
t602 = t419 * t423;
t349 = t596 + t602;
t424 = cos(qJ(2));
t319 = t349 * t424;
t671 = t319 / 0.2e1;
t425 = -pkin(2) - pkin(8);
t698 = -pkin(9) + t425;
t359 = t698 * t420;
t342 = t422 * t359;
t360 = t698 * t423;
t604 = t419 * t360;
t274 = t342 + t604;
t397 = pkin(4) * t420 + qJ(3);
t421 = sin(qJ(2));
t712 = t397 * t671 + t274 * t421 / 0.2e1;
t661 = -t421 / 0.2e1;
t591 = t424 * qJ(3);
t372 = t421 * pkin(2) - t591;
t344 = pkin(8) * t421 + t372;
t413 = t424 * pkin(7);
t414 = t424 * pkin(3);
t374 = t413 + t414;
t358 = t374 * t423;
t655 = pkin(9) * t421;
t211 = pkin(4) * t424 + t358 + (-t344 - t655) * t420;
t192 = t422 * t211;
t326 = t423 * t344;
t357 = t374 * t420;
t583 = t326 + t357;
t231 = t423 * t655 + t583;
t605 = t419 * t231;
t493 = t192 / 0.2e1 - t605 / 0.2e1;
t711 = t493 + t712;
t595 = t422 * t423;
t603 = t419 * t420;
t352 = t595 - t603;
t592 = t423 * t424;
t475 = pkin(4) * t592 + t374;
t448 = t475 * t352;
t52 = t448 / 0.2e1 + t493 - t712;
t506 = t419 * t359 - t422 * t360;
t692 = -t352 * qJ(6) - t506;
t678 = -t692 / 0.2e1;
t699 = t692 / 0.2e1;
t704 = t699 + t678;
t35 = t704 * t349;
t198 = -qJ(6) * t349 + t274;
t378 = t422 * t592;
t598 = t420 * t424;
t317 = t419 * t598 - t378;
t708 = -t397 * t317 / 0.2e1 + t506 * t661;
t647 = t422 * pkin(4);
t534 = t647 / 0.2e1;
t403 = pkin(5) + t647;
t665 = -t403 / 0.2e1;
t496 = t534 + t665;
t657 = pkin(5) * t349;
t170 = t657 / 0.2e1 - t496 * t349;
t561 = qJD(5) * t349;
t707 = -pkin(5) * t561 - qJD(4) * t170;
t649 = t419 * pkin(4);
t706 = -t170 * qJD(5) + qJD(4) * (-t403 * t349 + t352 * t649);
t337 = t352 * qJD(5);
t518 = -t603 / 0.2e1;
t377 = t421 * t595;
t667 = t377 / 0.2e1;
t668 = t352 / 0.2e1;
t437 = t667 + (t518 + t668) * t421;
t690 = t437 * qJD(1);
t705 = -t352 * qJD(4) - t337 - t690;
t599 = t420 * t421;
t316 = t419 * t599 - t377;
t318 = t421 * t349;
t481 = t316 * t349 + t318 * t352;
t703 = qJD(3) * t481;
t684 = pkin(5) / 0.2e1;
t464 = t684 + t496;
t695 = t464 * t349;
t702 = qJD(4) * t695;
t701 = qJD(5) * t695;
t547 = qJD(4) + qJD(5);
t191 = t419 * t211;
t214 = t422 * t231;
t513 = -t191 / 0.2e1 - t214 / 0.2e1;
t700 = t513 + t708;
t449 = t475 * t349;
t51 = -t449 / 0.2e1 + t513 - t708;
t436 = t352 * t661 + t421 * t518 + t667;
t697 = t436 * t547;
t696 = t437 * t547;
t636 = qJ(3) * t421;
t485 = -t424 * t425 + t636;
t339 = -pkin(1) - t485;
t682 = pkin(3) + pkin(7);
t373 = t682 * t421;
t356 = t423 * t373;
t267 = t339 * t420 - t356;
t227 = pkin(9) * t598 - t267;
t648 = t421 * pkin(4);
t208 = t227 + t648;
t600 = t420 * t373;
t268 = t339 * t423 + t600;
t228 = -pkin(9) * t592 + t268;
t213 = t422 * t228;
t125 = t419 * t208 + t213;
t635 = qJ(6) * t317;
t108 = t125 + t635;
t693 = t704 * t108;
t190 = t422 * t208;
t606 = t419 * t228;
t124 = -t190 + t606;
t308 = t319 * qJ(6);
t107 = -t124 + t308;
t620 = t198 * t319;
t99 = pkin(5) * t421 + t107;
t638 = t99 * t349;
t672 = t317 / 0.2e1;
t625 = t108 * t352;
t84 = -t625 / 0.2e1;
t689 = t692 * t672 + t317 * t678 + t620 / 0.2e1 + t638 / 0.2e1 + t84;
t415 = t420 ^ 2;
t417 = t423 ^ 2;
t386 = t415 - t417;
t550 = t424 * qJD(1);
t527 = t423 * t550;
t443 = qJD(2) * t386 + 0.2e1 * t420 * t527;
t549 = t424 * qJD(3);
t688 = qJD(2) * t485 - t549;
t687 = t547 * t506;
t597 = t422 * t227;
t127 = t597 - t606;
t114 = t127 + t308;
t536 = t649 / 0.2e1;
t664 = t403 / 0.2e1;
t446 = t316 * t536 + t318 * t664;
t607 = t419 * t227;
t126 = -t213 - t607;
t113 = t126 - t635;
t680 = t113 / 0.2e1;
t683 = -t99 / 0.2e1;
t11 = -t349 * (t114 / 0.2e1 + t683) - t352 * (t108 / 0.2e1 + t680) + t446;
t686 = t319 ^ 2;
t346 = t349 ^ 2;
t685 = t352 ^ 2;
t508 = t192 - t605;
t101 = pkin(5) * t424 - qJ(6) * t318 + t508;
t681 = t101 / 0.2e1;
t679 = -t190 / 0.2e1;
t676 = -t198 / 0.2e1;
t675 = t198 / 0.2e1;
t673 = -t208 / 0.2e1;
t516 = -t213 / 0.2e1;
t515 = -t342 / 0.2e1;
t670 = t349 / 0.2e1;
t669 = -t352 / 0.2e1;
t666 = t378 / 0.2e1;
t663 = -t419 / 0.2e1;
t662 = t420 / 0.2e1;
t660 = t422 / 0.2e1;
t659 = t423 / 0.2e1;
t658 = pkin(5) * t319;
t656 = pkin(5) * t352;
t654 = t113 * pkin(5);
t653 = t198 * pkin(5);
t652 = t316 * pkin(5);
t651 = t317 * pkin(5);
t650 = t318 * pkin(5);
t646 = t423 * pkin(4);
t511 = t675 + t676;
t36 = t352 * t511 + t35;
t645 = t36 * qJD(4) + t35 * qJD(5);
t643 = -t638 / 0.2e1 + t625 / 0.2e1;
t642 = pkin(4) * qJD(4);
t641 = pkin(4) * qJD(5);
t640 = pkin(5) * qJD(5);
t639 = pkin(5) * qJD(6);
t637 = t107 - t99;
t16 = t637 * t317;
t633 = qJD(1) * t16;
t245 = t475 - t651;
t19 = t108 * t637 - t245 * t658;
t632 = qJD(1) * t19;
t535 = -t648 / 0.2e1;
t497 = t535 + t227 / 0.2e1;
t47 = t422 * t497 + t679;
t631 = qJD(1) * t47;
t450 = t475 * t319;
t546 = pkin(4) * t598;
t63 = t126 * t421 + t317 * t546 - t450;
t630 = qJD(1) * t63;
t249 = t475 * t317;
t64 = t127 * t421 - t319 * t546 - t249;
t629 = qJD(1) * t64;
t77 = -t124 * t421 - t249;
t628 = qJD(1) * t77;
t78 = -t125 * t421 - t450;
t627 = qJD(1) * t78;
t626 = t101 * t352;
t586 = t214 + t191;
t110 = -qJ(6) * t316 + t586;
t624 = t110 * t349;
t25 = t108 * t316 + t318 * t99;
t15 = t101 * t319 + t110 * t317 - t25;
t623 = t15 * qJD(1);
t328 = (-t646 - t682) * t421;
t244 = t328 + t652;
t17 = t101 * t99 + t108 * t110 + t244 * t245;
t622 = t17 * qJD(1);
t18 = (t108 + t113) * t319 + (t114 - t99) * t317;
t621 = t18 * qJD(1);
t465 = -t546 - t658;
t20 = t108 * t114 + t99 * t113 + t245 * t465;
t619 = t20 * qJD(1);
t618 = t25 * qJD(1);
t27 = -t124 * t424 + t316 * t475 - t328 * t317 + t421 * t508;
t617 = t27 * qJD(1);
t28 = t125 * t424 - t318 * t475 + t328 * t319 + t421 * t586;
t614 = t28 * qJD(1);
t613 = t316 * t198;
t612 = t318 * t692;
t611 = t349 * t317;
t610 = t352 * t319;
t601 = t420 * t344;
t594 = t423 * t349;
t593 = t423 * t352;
t121 = t610 - t611;
t589 = t547 * t121;
t142 = t317 * t668 + t319 * t670;
t588 = t547 * t142;
t460 = t596 / 0.2e1 + t602 / 0.2e1;
t241 = t349 * t661 + t421 * t460;
t585 = t547 * t241;
t238 = (t670 + t460) * t421;
t584 = t547 * t238;
t416 = t421 ^ 2;
t418 = t424 ^ 2;
t387 = t418 - t416;
t130 = t316 * t317 + t318 * t319;
t581 = qJD(1) * t130;
t177 = -t267 * t421 + t374 * t592;
t580 = qJD(1) * t177;
t178 = -t268 * t421 - t374 * t598;
t579 = qJD(1) * t178;
t200 = -t316 * t421 + t317 * t424;
t578 = qJD(1) * t200;
t201 = t318 * t421 - t319 * t424;
t577 = qJD(1) * t201;
t576 = qJD(1) * t319;
t351 = t387 * t420;
t575 = qJD(1) * t351;
t355 = t387 * t423;
t574 = qJD(1) * t355;
t573 = qJD(1) * t421;
t572 = qJD(1) * t423;
t571 = qJD(2) * qJ(3);
t570 = qJD(2) * t352;
t569 = qJD(2) * t397;
t568 = qJD(3) * t420;
t567 = qJD(3) * t421;
t566 = qJD(3) * t423;
t565 = qJD(4) * t420;
t564 = qJD(4) * t421;
t563 = qJD(4) * t423;
t562 = qJD(4) * t425;
t117 = (-t267 - t356) * t424 - t601 * t421;
t560 = t117 * qJD(1);
t118 = t583 * t421 - t374 * t599 + (t268 - t600) * t424;
t559 = t118 * qJD(1);
t131 = t316 * t319 + t317 * t318;
t558 = t131 * qJD(1);
t220 = t238 * qJD(1);
t492 = -pkin(2) * t424 - t636;
t361 = -pkin(1) + t492;
t282 = t361 * t424 + t372 * t421;
t557 = t282 * qJD(1);
t283 = -t361 * t421 + t372 * t424;
t556 = t283 * qJD(1);
t555 = t387 * qJD(1);
t554 = t416 * qJD(1);
t553 = t420 * qJD(2);
t552 = t421 * qJD(2);
t551 = t423 * qJD(2);
t410 = t424 * qJD(2);
t548 = t424 * qJD(4);
t545 = pkin(1) * t573;
t544 = pkin(1) * t550;
t542 = pkin(7) * t552;
t541 = t419 * t641;
t540 = -t658 / 0.2e1;
t539 = t656 / 0.2e1;
t538 = -t650 / 0.2e1;
t537 = t650 / 0.2e1;
t533 = -t646 / 0.2e1;
t532 = t664 + t684;
t531 = t107 / 0.2e1 + t683;
t529 = t316 * t573;
t528 = t318 * t573;
t526 = t420 * t551;
t525 = t420 * t410;
t524 = t421 * t548;
t523 = t361 * t372 * qJD(1);
t522 = t361 * t573;
t392 = t421 * t410;
t391 = t421 * t550;
t394 = t423 * t410;
t521 = t420 * t563;
t520 = t419 * t672;
t519 = t349 * t663;
t517 = -t598 / 0.2e1;
t509 = pkin(4) * t547;
t504 = t547 * t349;
t503 = t547 * t421;
t502 = pkin(4) * t517;
t501 = qJD(4) + t573;
t500 = t420 * t394;
t498 = -t346 / 0.2e1 - t685 / 0.2e1;
t495 = t646 + t656;
t461 = t107 * t670 + t84;
t427 = -t620 / 0.2e1 - t461 + t689;
t8 = t537 + t427;
t491 = qJD(1) * t8;
t5 = -t317 * t704 - t319 * t511 - t11;
t490 = t5 * qJD(1);
t489 = qJD(5) + t501;
t484 = -t11 * qJD(1) + t36 * qJD(2);
t438 = t461 + t643;
t14 = t538 + t438;
t483 = -qJD(1) * t14 - qJD(2) * t35;
t432 = -t108 * t349 / 0.2e1 + t692 * t671 + t198 * t672 + t99 * t669;
t434 = (-pkin(7) / 0.2e1 - pkin(3) / 0.2e1 + t533) * t421 + t652 / 0.2e1;
t21 = -t432 + t434;
t95 = -t198 * t349 - t352 * t692;
t482 = -qJD(1) * t21 + qJD(2) * t95;
t428 = -t413 / 0.2e1 - t414 / 0.2e1 - t613 / 0.2e1 + t651 / 0.2e1 - t612 / 0.2e1 + t424 * t533;
t462 = t626 / 0.2e1 + t624 / 0.2e1;
t23 = t428 + t462;
t307 = t397 + t657;
t480 = -qJD(1) * t23 + qJD(2) * t307;
t232 = pkin(4) * t594 + t352 * t397;
t40 = t374 * t669 + (t317 * t659 + (t660 + t349 * t662 - t593 / 0.2e1) * t424) * pkin(4) + t711;
t479 = -t40 * qJD(1) + t232 * qJD(2);
t233 = pkin(4) * t593 - t349 * t397;
t39 = t374 * t670 + (t319 * t659 + (t663 + t352 * t662 + t594 / 0.2e1) * t424) * pkin(4) + t700;
t478 = t39 * qJD(1) - t233 * qJD(2);
t272 = t515 + t342 / 0.2e1;
t45 = t516 + t213 / 0.2e1 + (t673 + t497) * t419;
t477 = qJD(1) * t45 + qJD(2) * t272;
t141 = t611 / 0.2e1 + t610 / 0.2e1;
t26 = t108 * t317 + t319 * t99;
t476 = qJD(1) * t26 + qJD(3) * t141;
t310 = t317 ^ 2;
t165 = t310 - t686;
t55 = qJD(1) * t165 + qJD(2) * t121;
t215 = t346 - t685;
t65 = qJD(1) * t121 + qJD(2) * t215;
t128 = t532 * t319 + (t520 + t598 / 0.2e1) * pkin(4);
t163 = -t532 * t352 + (t519 - t423 / 0.2e1) * pkin(4);
t474 = qJD(1) * t128 + qJD(2) * t163;
t133 = t464 * t317;
t473 = qJD(1) * t133 - qJD(2) * t695;
t167 = -0.1e1 / 0.2e1 + t498;
t472 = qJD(1) * t141 + qJD(2) * t167;
t153 = -t610 - t611;
t229 = t310 + t686;
t471 = qJD(1) * t229 + qJD(2) * t153;
t275 = t346 + t685;
t470 = qJD(1) * t153 + qJD(2) * t275;
t237 = t666 + (t518 + t669) * t424;
t469 = qJD(1) * t237 + qJD(2) * t349;
t239 = (t670 - t460) * t424;
t468 = qJD(1) * t239 + t570;
t467 = -t570 + t576;
t466 = -t554 - t564;
t459 = t425 * t661 - t591 / 0.2e1;
t50 = t449 / 0.2e1 + t700;
t457 = qJD(1) * t50 + t349 * t569;
t49 = -t448 / 0.2e1 + t711;
t456 = qJD(1) * t49 - t352 * t569;
t455 = t501 * t598;
t246 = (t344 / 0.2e1 + t459) * t420;
t454 = -qJ(3) * t551 - t246 * qJD(1);
t440 = t459 * t423;
t247 = t326 / 0.2e1 + t440;
t453 = qJ(3) * t553 - t247 * qJD(1);
t112 = -qJD(2) * t142 + t317 * t576;
t123 = -qJD(1) * t142 + t349 * t570;
t341 = (t417 / 0.2e1 - t415 / 0.2e1) * t424;
t452 = qJD(1) * t341 + t526;
t447 = t101 * t664 + t110 * t536;
t445 = t418 * t420 * t572 - qJD(2) * t341;
t354 = t386 * t418;
t444 = -qJD(1) * t354 + 0.2e1 * t500;
t426 = t108 * t699 + t692 * t680 + t114 * t675 + t245 * t495 / 0.2e1 + t99 * t676 + t465 * t307 / 0.2e1;
t1 = -t426 + t447;
t43 = t307 * t495;
t442 = -t1 * qJD(1) + t43 * qJD(2) + t36 * qJD(3);
t3 = -t531 * t198 - t693 + (t245 * t669 + t307 * t671 + t681) * pkin(5);
t44 = t307 * t656;
t441 = -qJD(1) * t3 + qJD(2) * t44 + qJD(3) * t35;
t435 = qJD(2) * t492 + t549;
t431 = (t108 * t660 + t419 * t531) * pkin(4) + t108 * t665;
t10 = -t654 / 0.2e1 + t431;
t429 = (t198 * t660 + t419 * t704) * pkin(4) + t198 * t665;
t30 = t653 / 0.2e1 + t429;
t320 = (-t403 + t647) * t649;
t433 = -qJD(1) * t10 - qJD(2) * t30 - qJD(3) * t695 - qJD(4) * t320;
t430 = t113 * t668 + t114 * t670 + t446;
t409 = pkin(7) * t410;
t401 = -t550 / 0.2e1;
t400 = t550 / 0.2e1;
t399 = t410 / 0.2e1;
t393 = t421 * t572;
t390 = t420 * t573;
t348 = -t393 - t563;
t347 = -t390 - t565;
t345 = t391 + t548 / 0.2e1;
t327 = t341 * qJD(4);
t309 = t391 + (qJD(4) / 0.2e1 + qJD(5) / 0.2e1) * t424;
t240 = -t319 / 0.2e1 - t460 * t424;
t236 = t419 * t517 + t424 * t668 + t666;
t230 = t467 * pkin(5);
t226 = t436 * qJD(3);
t223 = t241 * qJD(3);
t219 = t238 * qJD(3);
t216 = t437 * qJD(3);
t194 = 0.2e1 * t515 - t604;
t184 = -t357 - t326 / 0.2e1 + t440;
t183 = t358 - t601 / 0.2e1 + t459 * t420;
t166 = 0.1e1 / 0.2e1 + t498;
t162 = pkin(4) * t519 + t352 * t665 + t646 / 0.2e1 + t539;
t155 = qJD(2) * t437 - t319 * t573;
t154 = qJD(2) * t238 - t317 * t573;
t140 = t153 * qJD(6);
t137 = t141 * qJD(6);
t134 = -t504 - t220;
t132 = -t651 / 0.2e1 + t496 * t317;
t129 = pkin(4) * t520 + t319 * t664 + t502 + t540;
t116 = qJD(2) * t436 + t319 * t489;
t115 = qJD(2) * t241 + t317 * t489;
t48 = t422 * t535 + t606 + t679 - t597 / 0.2e1;
t46 = 0.2e1 * t516 - t607 / 0.2e1 + (t535 + t673) * t419;
t42 = t319 * t533 + t352 * t502 - t424 * t649 / 0.2e1 + t51;
t41 = t317 * t533 + t349 * t502 + t424 * t534 + t52;
t29 = -t653 / 0.2e1 + t429;
t24 = -t428 + t462;
t22 = t432 + t434;
t13 = t537 + t438;
t12 = t430 + t643;
t9 = t654 / 0.2e1 + t431;
t7 = t538 + t427;
t6 = -t198 * t671 - t430 + t689;
t4 = pkin(5) * t681 + t107 * t675 + t198 * t683 + t245 * t539 + t307 * t540 + t693;
t2 = t426 + t447;
t31 = [0, 0, 0, t392, t387 * qJD(2), 0, 0, 0, -pkin(1) * t552, -pkin(1) * t410, 0, qJD(2) * t283 - t421 * t549, -qJD(2) * t282 + qJD(3) * t416 (qJD(2) * t372 - t567) * t361, -t392 * t415 + t418 * t521, -qJD(4) * t354 - 0.2e1 * t421 * t500, -qJD(2) * t351 - t423 * t524, -qJD(2) * t355 + t420 * t524, t392, qJD(2) * t117 + qJD(4) * t178 + t416 * t568, -qJD(2) * t118 - qJD(4) * t177 + t416 * t566 (-qJD(2) * t318 - t317 * t547) * t319, qJD(2) * t131 + t165 * t547, qJD(2) * t201 + t317 * t503, qJD(2) * t200 + t319 * t503, t392, qJD(2) * t27 + qJD(4) * t63 + qJD(5) * t78 + t318 * t567, -qJD(2) * t28 - qJD(4) * t64 - qJD(5) * t77 - t316 * t567, qJD(2) * t15 + qJD(3) * t130 + qJD(4) * t18 + qJD(5) * t16 + qJD(6) * t229, qJD(2) * t17 + qJD(3) * t25 + qJD(4) * t20 + qJD(5) * t19 + qJD(6) * t26; 0, 0, 0, t391, t555, t410, -t552, 0, -t409 - t545, t542 - t544, t435, t409 + t556, -t542 - t557, pkin(7) * t435 + t523, -t327 + (-t415 * t550 + t526) * t421, -t421 * t443 + 0.2e1 * t424 * t521, t394 - t575, -t525 - t574, t345, t183 * qJD(4) - t373 * t553 - t423 * t688 + t560, t184 * qJD(4) - t373 * t551 + t420 * t688 - t559, -t318 * t467 + t588, t558 + (-t316 * t352 - t318 * t349) * qJD(2) + t589, t352 * t410 + t577 + t585, -t349 * t410 + t578 + t697, t309, t617 + (t316 * t397 + t328 * t349 - t424 * t506) * qJD(2) + t236 * qJD(3) + t41 * qJD(4) + t52 * qJD(5), -t614 + (-t274 * t424 + t318 * t397 + t328 * t352) * qJD(2) + t240 * qJD(3) + t42 * qJD(4) + t51 * qJD(5), t623 + (-t612 - t613 - t624 - t626) * qJD(2) - t703 + t6 * qJD(4) + t7 * qJD(5) + t140, t622 + (t101 * t692 + t110 * t198 + t244 * t307) * qJD(2) + t24 * qJD(3) + t2 * qJD(4) + t4 * qJD(5) + t22 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t410, -t391, t554, t409 - t522, 0, 0, 0, 0, 0, t420 * t554 + t394, t423 * t554 - t525, 0, 0, 0, 0, 0, qJD(2) * t236 + t528 + t585, qJD(2) * t240 - t529 + t697, -qJD(2) * t481 + t581, t24 * qJD(2) + t12 * qJD(4) + t13 * qJD(5) + t137 + t618 + t703; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t445, t444, -t501 * t592, t455, t399, qJD(2) * t183 - qJD(4) * t268 + t579, qJD(2) * t184 + qJD(4) * t267 - t580, -t112, t55, t115, t116, t399, qJD(2) * t41 + qJD(4) * t126 + qJD(5) * t46 + t223 + t630, qJD(2) * t42 - qJD(4) * t127 + qJD(5) * t48 + t226 - t629, t621 + t6 * qJD(2) + (-t403 * t317 + t319 * t649) * qJD(4) + t132 * qJD(5), t619 + t2 * qJD(2) + t12 * qJD(3) + (t113 * t403 + t114 * t649) * qJD(4) + t9 * qJD(5) + t129 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112, t55, t115, t116, t399, qJD(2) * t52 + qJD(4) * t46 - qJD(5) * t125 + t223 + t627, qJD(2) * t51 + qJD(4) * t48 + qJD(5) * t124 + t226 - t628, qJD(2) * t7 + qJD(4) * t132 - t317 * t640 + t633, qJD(2) * t4 + qJD(3) * t13 + qJD(4) * t9 - t108 * t640 + t632; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t471, qJD(2) * t22 + qJD(4) * t129 + t476; 0, 0, 0, -t391, -t555, 0, 0, 0, t545, t544, 0, -t556, t557, -t523, t391 * t415 - t327, 0.2e1 * t423 * t455, -t420 * t564 + t575, -t421 * t563 + t574, -t345, qJD(4) * t246 - t560, qJD(4) * t247 + t559, t318 * t576 + t588, -t558 + t589, -t577 - t584, -t578 - t696, -t309, qJD(3) * t237 - qJD(4) * t40 - qJD(5) * t49 - t617, qJD(3) * t239 - qJD(4) * t39 - qJD(5) * t50 + t614, -qJD(4) * t5 + qJD(5) * t8 + t140 - t623, -qJD(3) * t23 - qJD(4) * t1 - qJD(5) * t3 - qJD(6) * t21 - t622; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3), -t521, t386 * qJD(4), 0, 0, 0, qJ(3) * t563 + t568, -qJ(3) * t565 + t566, -t352 * t504, t547 * t215, 0, 0, 0, qJD(3) * t349 + qJD(4) * t232 + t337 * t397, qJD(3) * t352 + qJD(4) * t233 - t397 * t561, qJD(6) * t275, qJD(3) * t307 + qJD(4) * t43 + qJD(5) * t44 + qJD(6) * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t571, 0, 0, 0, 0, 0, t553, t551, 0, 0, 0, 0, 0, t469, t468, 0, qJD(6) * t166 + t480 + t645; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t452, t443, t347, t348, t401, -t420 * t562 - t454, -t423 * t562 - t453, -t123, t65, t134, t705, t401, -qJD(4) * t274 + t194 * qJD(5) + t479, -t478 + t687, -t490 - t706 (-t198 * t403 + t649 * t692) * qJD(4) + t29 * qJD(5) + t162 * qJD(6) + t442; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, t65, t134, t705, t401, qJD(4) * t194 - qJD(5) * t274 - t456, -t457 + t687, t491 - t707, qJD(4) * t29 - t198 * t640 + t441; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t470, qJD(3) * t166 + qJD(4) * t162 + t482; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t391, -t554, t522, 0, 0, 0, 0, 0, t466 * t420, t466 * t423, 0, 0, 0, 0, 0, -qJD(2) * t237 - t528 - t584, -qJD(2) * t239 + t529 - t696, -t581, qJD(2) * t23 - qJD(4) * t11 + qJD(5) * t14 + t137 - t618; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t571, 0, 0, 0, 0, 0, -t553, -t551, 0, 0, 0, 0, 0, -t469, -t468, 0, qJD(6) * t167 - t480 + t645; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t347, t348, 0, 0, 0, 0, 0, t134, t705, 0, t484 + t706; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, t705, 0, -t483 + t707; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t472; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t445, -t444 (t527 + t553) * t421 (-t420 * t550 + t551) * t421, t399, -qJD(2) * t246 + t420 * t567 - t579, -qJD(2) * t247 + t421 * t566 + t580, t112, -t55, t154, t155, t399, qJD(2) * t40 + qJD(5) * t45 + t219 - t630, qJD(2) * t39 + qJD(5) * t47 + t216 + t629, qJD(2) * t5 + qJD(5) * t133 - t621, qJD(2) * t1 + qJD(3) * t11 + qJD(5) * t10 + qJD(6) * t128 - t619; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t452, -t443, t390, t393, t400, t454, t453, t123, -t65, t220, t690, t400, qJD(5) * t272 - t479, t478, t490 - t701, qJD(5) * t30 + qJD(6) * t163 - t442; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t390, t393, 0, 0, 0, 0, 0, t220, t690, 0, -t484 + t701; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t541, -t422 * t641, 0, t320 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t419 * t509 + t477, -t422 * t509 + t631, t473, -pkin(5) * t541 - t433; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t474; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, -t55, t154, t155, t399, qJD(2) * t49 - qJD(4) * t45 + t219 - t627, qJD(2) * t50 - qJD(4) * t47 + t216 + t628, -qJD(2) * t8 - qJD(4) * t133 - t633, qJD(2) * t3 - qJD(3) * t14 - qJD(4) * t10 + t319 * t639 - t632; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123, -t65, t220, t690, t400, -qJD(4) * t272 + t456, t457, -t491 + t702, -qJD(4) * t30 - t352 * t639 - t441; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t220, t690, 0, t483 - t702; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t419 * t642 - t477, t422 * t642 - t631, -t473, t433; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t471, qJD(2) * t21 - qJD(4) * t128 - t319 * t640 - t476; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t470, pkin(5) * t337 - qJD(3) * t167 - qJD(4) * t163 - t482; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t472; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t474; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t230; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t31;
