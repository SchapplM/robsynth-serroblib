% Calculate minimal parameter regressor of coriolis matrix for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x26]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRRPR10_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:29:54
% EndTime: 2019-12-31 21:30:22
% DurationCPUTime: 13.90s
% Computational Cost: add. (11751->632), mult. (28992->935), div. (0->0), fcn. (32312->10), ass. (0->480)
t411 = cos(pkin(5));
t413 = sin(qJ(3));
t416 = cos(qJ(3));
t410 = sin(pkin(5));
t414 = sin(qJ(2));
t617 = t410 * t414;
t359 = -t411 * t416 + t413 * t617;
t361 = t411 * t413 + t416 * t617;
t409 = sin(pkin(10));
t655 = cos(pkin(10));
t263 = t655 * t359 + t361 * t409;
t415 = cos(qJ(5));
t693 = t415 * t263;
t695 = t693 / 0.2e1;
t696 = -t693 / 0.2e1;
t697 = t695 + t696;
t698 = qJD(5) * t697;
t683 = t263 / 0.2e1;
t659 = -qJ(4) - pkin(8);
t384 = t659 * t416;
t506 = t655 * t413;
t305 = -t384 * t409 - t659 * t506;
t694 = -t305 / 0.2e1;
t338 = t655 * t361;
t620 = t409 * t359;
t690 = t338 - t620;
t685 = t690 / 0.2e1;
t412 = sin(qJ(5));
t151 = t412 * t690;
t417 = cos(qJ(2));
t616 = t410 * t417;
t394 = pkin(7) * t616;
t395 = t411 * t414 * pkin(1);
t364 = t394 + t395;
t335 = pkin(8) * t411 + t364;
t495 = -pkin(2) * t417 - pkin(8) * t414;
t336 = (-pkin(1) + t495) * t410;
t238 = t413 * t335 - t416 * t336;
t193 = -t361 * qJ(4) - t238;
t168 = -pkin(3) * t616 + t193;
t239 = t335 * t416 + t336 * t413;
t194 = -qJ(4) * t359 + t239;
t191 = t655 * t194;
t102 = t409 * t168 + t191;
t692 = t102 / 0.2e1;
t619 = t409 * t413;
t444 = -t655 * t384 + t659 * t619;
t691 = -t444 / 0.2e1;
t630 = t444 * t412;
t629 = t444 * t415;
t505 = t655 * t416;
t368 = -t505 + t619;
t509 = 0.2e1 * t683;
t494 = t509 * t368;
t366 = t368 ^ 2;
t618 = t409 * t416;
t371 = t506 + t618;
t367 = t371 ^ 2;
t689 = -t367 - t366;
t404 = t412 ^ 2;
t406 = t415 ^ 2;
t392 = t406 - t404;
t688 = qJD(3) * t392;
t572 = qJD(2) * t415;
t537 = t412 * t572;
t498 = t371 * t537;
t687 = 0.2e1 * t498 - t688;
t621 = t409 * t194;
t108 = t655 * t193 - t621;
t686 = t108 / 0.2e1;
t684 = -t263 / 0.2e1;
t680 = t305 / 0.2e1;
t662 = pkin(1) * t417;
t362 = pkin(7) * t617 - t411 * t662;
t334 = -pkin(2) * t411 + t362;
t677 = t334 / 0.2e1;
t676 = -t361 / 0.2e1;
t675 = -t368 / 0.2e1;
t674 = t368 / 0.2e1;
t673 = -t371 / 0.2e1;
t672 = t371 / 0.2e1;
t399 = pkin(3) * t409 + pkin(9);
t671 = t399 / 0.2e1;
t400 = -t655 * pkin(3) - pkin(4);
t670 = -t400 / 0.2e1;
t669 = t400 / 0.2e1;
t668 = -t409 / 0.2e1;
t667 = t409 / 0.2e1;
t666 = -t412 / 0.2e1;
t665 = t412 / 0.2e1;
t664 = -t413 / 0.2e1;
t663 = t415 / 0.2e1;
t661 = t361 * pkin(3);
t660 = t413 * pkin(3);
t658 = pkin(3) * qJD(3);
t107 = t193 * t409 + t191;
t225 = t415 * t616 + t151;
t274 = pkin(3) * t359 + t334;
t422 = pkin(4) * t263 - pkin(9) * t690 + t274;
t96 = -pkin(9) * t616 + t102;
t36 = t412 * t96 - t415 * t422;
t139 = pkin(4) * t690 + pkin(9) * t263 + t661;
t602 = t415 * t139;
t615 = t412 * t108;
t481 = t602 - t615;
t610 = t412 * t263;
t101 = t655 * t168 - t621;
t95 = pkin(4) * t616 - t101;
t5 = t107 * t225 + t263 * t481 - t36 * t690 - t95 * t610;
t657 = t5 * qJD(1);
t227 = -t412 * t616 + t415 * t690;
t37 = t412 * t422 + t415 * t96;
t604 = t415 * t108;
t613 = t412 * t139;
t482 = t604 + t613;
t6 = t107 * t227 - t263 * t482 - t37 * t690 - t693 * t95;
t656 = t6 * qJD(1);
t21 = -t225 * t95 + t263 * t36;
t654 = qJD(1) * t21;
t22 = t227 * t95 - t263 * t37;
t653 = qJD(1) * t22;
t29 = -t101 * t690 - t102 * t263;
t652 = qJD(1) * t29;
t550 = t263 * t610;
t633 = t690 * t225;
t58 = -t550 - t633;
t651 = qJD(1) * t58;
t59 = -t550 + t633;
t650 = qJD(1) * t59;
t549 = t263 * t693;
t636 = t227 * t690;
t60 = -t549 + t636;
t649 = qJD(1) * t60;
t61 = t549 + t636;
t648 = qJD(1) * t61;
t528 = t227 * t675;
t445 = -t263 * t406 * t672 + t415 * t528;
t548 = t413 * t616;
t314 = -t409 * t548 + t505 * t616;
t594 = t415 * t314;
t270 = t412 * t617 + t594;
t631 = t270 * t412;
t67 = -t631 / 0.2e1 + t445;
t647 = qJD(1) * t67;
t521 = -t617 / 0.2e1;
t529 = t225 * t675;
t632 = t263 * t371;
t76 = t529 - t594 / 0.2e1 + (-t632 / 0.2e1 + t521) * t412;
t646 = qJD(1) * t76;
t453 = t371 * t696 + t528;
t386 = t415 * t617;
t606 = t412 * t314;
t491 = t386 / 0.2e1 - t606 / 0.2e1;
t77 = t453 - t491;
t645 = qJD(1) * t77;
t644 = t107 * t412;
t643 = t107 * t415;
t363 = (pkin(2) * t414 - pkin(8) * t417) * t410;
t348 = t416 * t363;
t605 = t413 * t362;
t504 = t348 + t605;
t593 = t416 * t417;
t224 = (pkin(3) * t414 - qJ(4) * t593) * t410 + t504;
t346 = t413 * t363;
t347 = t416 * t362;
t585 = t346 - t347;
t234 = -qJ(4) * t548 + t585;
t132 = t655 * t224 - t409 * t234;
t119 = -pkin(4) * t617 - t132;
t269 = -t386 + t606;
t313 = t371 * t616;
t383 = pkin(3) * t548;
t311 = t383 + t364;
t165 = pkin(4) * t313 - pkin(9) * t314 + t311;
t601 = t415 * t165;
t133 = t409 * t224 + t655 * t234;
t120 = pkin(9) * t617 + t133;
t614 = t412 * t120;
t54 = t601 - t614;
t11 = t119 * t225 + t263 * t54 + t269 * t95 - t313 * t36;
t642 = t11 * qJD(1);
t603 = t415 * t120;
t612 = t412 * t165;
t55 = t603 + t612;
t12 = t119 * t227 - t263 * t55 + t270 * t95 - t313 * t37;
t641 = t12 * qJD(1);
t17 = (-t102 + t107) * t690 + (t101 - t108) * t263;
t640 = t17 * qJD(1);
t18 = -t101 * t314 - t102 * t313 - t132 * t690 - t133 * t263;
t639 = t18 * qJD(1);
t19 = -t101 * t107 + t102 * t108 + t274 * t661;
t638 = t19 * qJD(1);
t20 = t101 * t132 + t102 * t133 + t274 * t311;
t637 = t20 * qJD(1);
t635 = t227 * t412;
t634 = t227 * t415;
t628 = t334 * t416;
t280 = t412 * t371;
t600 = t415 * t225;
t456 = t600 / 0.2e1 + t635 / 0.2e1;
t434 = t280 * t693 + t368 * t456;
t458 = t269 * t666 + t270 * t663;
t35 = t434 - t458;
t627 = t35 * qJD(1);
t626 = t368 * t371;
t625 = t371 * t225;
t624 = t371 * t415;
t403 = t410 ^ 2;
t408 = t417 ^ 2;
t623 = t403 * t408;
t622 = t403 * t414;
t611 = t412 * t225;
t285 = pkin(4) * t371 + pkin(9) * t368 + t660;
t609 = t412 * t285;
t608 = t305 * t412;
t607 = t412 * t313;
t597 = t415 * t285;
t596 = t305 * t415;
t595 = t415 * t313;
t65 = -t600 - t635;
t56 = t65 * t263;
t592 = t56 * qJD(1);
t57 = -t225 * t270 - t227 * t269;
t591 = t57 * qJD(1);
t85 = t238 * t617 - t334 * t548 - t364 * t359 + t504 * t616;
t590 = t85 * qJD(1);
t86 = t364 * t361 + (-t239 * t414 + (t585 + t628) * t417) * t410;
t589 = t86 * qJD(1);
t87 = -t225 * t313 - t263 * t269;
t588 = t87 * qJD(1);
t88 = t227 * t313 + t263 * t270;
t587 = t88 * qJD(1);
t393 = -t413 ^ 2 + t416 ^ 2;
t148 = -t238 * t616 - t334 * t359;
t584 = qJD(1) * t148;
t149 = -t239 * t616 - t334 * t361;
t583 = qJD(1) * t149;
t156 = 0.2e1 * t696;
t582 = qJD(1) * t156;
t581 = qJD(1) * t227;
t580 = qJD(1) * t263;
t365 = t506 / 0.2e1 + t618 / 0.2e1;
t310 = t365 * t616;
t579 = qJD(1) * t310;
t578 = qJD(1) * t361;
t577 = qJD(1) * t411;
t576 = qJD(1) * t417;
t575 = qJD(2) * t368;
t574 = qJD(2) * t371;
t573 = qJD(2) * t413;
t571 = qJD(2) * t416;
t570 = qJD(2) * t417;
t569 = qJD(3) * t412;
t568 = qJD(3) * t413;
t567 = qJD(3) * t415;
t566 = qJD(3) * t416;
t565 = qJD(4) * t415;
t564 = qJD(5) * t263;
t563 = qJD(5) * t310;
t562 = qJD(5) * t412;
t561 = qJD(5) * t415;
t199 = -t359 * t416 - t361 * t413;
t232 = t199 * t616;
t560 = t232 * qJD(1);
t288 = -t359 * t617 + t413 * t623;
t559 = t288 * qJD(1);
t289 = -t361 * t617 + t416 * t623;
t558 = t289 * qJD(1);
t308 = pkin(1) * t622 + t364 * t411;
t557 = t308 * qJD(1);
t309 = -t362 * t411 + t403 * t662;
t556 = t309 * qJD(1);
t370 = (-t414 ^ 2 + t408) * t403;
t555 = t370 * qJD(1);
t554 = t410 * qJD(3);
t553 = t367 - t366;
t552 = t661 / 0.2e1;
t551 = t660 / 0.2e1;
t547 = t410 * t593;
t546 = t95 * t665;
t545 = t95 * t663;
t544 = -t655 / 0.2e1;
t543 = t95 / 0.2e1 + t686;
t401 = -pkin(3) * t416 - pkin(2);
t542 = t410 * t576;
t541 = t368 * t574;
t540 = t406 * t574;
t539 = t371 * t572;
t538 = t410 * t570;
t536 = t417 * t554;
t535 = t412 * t567;
t534 = t368 * t561;
t533 = qJD(3) * t626;
t532 = t403 * t576;
t531 = qJD(2) * t617;
t530 = t412 * t561;
t357 = t371 * t567;
t527 = t227 * t672;
t526 = t690 * t672;
t525 = t624 / 0.2e1;
t524 = t263 * t671;
t523 = t225 * t670;
t522 = t227 * t669;
t520 = t617 / 0.2e1;
t519 = -t614 / 0.2e1;
t518 = t610 / 0.2e1;
t517 = -t607 / 0.2e1;
t516 = t607 / 0.2e1;
t515 = -t280 / 0.2e1;
t514 = t280 / 0.2e1;
t513 = -t603 / 0.2e1;
t510 = -t595 / 0.2e1;
t508 = t680 + t694;
t507 = t346 / 0.2e1 - t347 / 0.2e1;
t503 = qJD(2) + t577;
t502 = -qJD(5) - t580;
t501 = -qJD(5) - t575;
t500 = t570 * t622;
t499 = t414 * t532;
t497 = t416 * t542;
t496 = t165 / 0.2e1 + t95 * t673;
t493 = 0.2e1 * t412 * t357;
t490 = t524 + t139 / 0.2e1;
t489 = t383 / 0.2e1 + t394 / 0.2e1 + t395 / 0.2e1;
t488 = -qJD(3) + t542;
t452 = pkin(4) * t368 - pkin(9) * t371 + t401;
t178 = -t415 * t452 + t630;
t419 = t178 * t685 + t225 * t691 + t481 * t675 + (t597 + t608) * t684 + t305 * t518 + t368 * t546;
t427 = -t119 * t415 / 0.2e1 + t269 * t669 + t399 * t517;
t1 = (t36 / 0.2e1 - t644 / 0.2e1) * t371 + t419 + t427;
t48 = (-t178 + t630) * t371 + t597 * t368;
t487 = -t1 * qJD(1) + t48 * qJD(2);
t179 = t412 * t452 + t629;
t418 = t179 * t685 + t227 * t691 + t482 * t674 + (-t596 + t609) * t683 + t305 * t695 + t368 * t545;
t428 = t119 * t665 + t270 * t669 + t399 * t510;
t2 = (t37 / 0.2e1 - t643 / 0.2e1) * t371 + t418 + t428;
t49 = (-t179 + t629) * t371 - t609 * t368;
t486 = -t2 * qJD(1) + t49 * qJD(2);
t423 = t444 * t685 + t690 * t691;
t441 = (t313 * t668 + t314 * t544) * pkin(3);
t9 = (t692 - t107 / 0.2e1) * t371 + (t686 - t101 / 0.2e1) * t368 + t441 + t423;
t485 = t9 * qJD(1);
t118 = t401 * t660;
t424 = t108 * t691 + t305 * t692 + t107 * t694 + t101 * t444 / 0.2e1;
t450 = t133 * t667 + t132 * t655 / 0.2e1;
t7 = (t274 * t664 + t401 * t676 + t450) * pkin(3) + t424;
t484 = -t7 * qJD(1) + t118 * qJD(2);
t51 = 0.2e1 * t685 * t371 + t494;
t79 = t263 ^ 2 + t690 ^ 2;
t483 = qJD(1) * t79 + qJD(2) * t51;
t480 = -t263 * t400 - t399 * t690;
t479 = -t368 * t400 - t371 * t399;
t111 = t178 * t368 - t305 * t280;
t436 = t178 * t684 + t225 * t680 + t36 * t675;
t14 = -t412 * t496 + t436 + t513;
t478 = -qJD(1) * t14 + qJD(2) * t111;
t112 = -t179 * t368 + t305 * t624;
t435 = t179 * t683 + t227 * t694 + t37 * t674;
t13 = t415 * t496 + t435 + t519;
t477 = qJD(1) * t13 - qJD(2) * t112;
t158 = t305 * t371 - t368 * t444;
t425 = t101 * t672 + t102 * t674 - t263 * t691 + t690 * t694;
t27 = t425 + t489;
t476 = -qJD(1) * t27 + qJD(2) * t158;
t235 = t553 * t412;
t462 = t412 * t494;
t429 = t690 * t515 + t462 - t625 / 0.2e1;
t39 = t510 + t429;
t475 = -qJD(1) * t39 + qJD(2) * t235;
t236 = t689 * t412;
t430 = t690 * t514 + t462 + t625 / 0.2e1;
t41 = t510 + t430;
t474 = qJD(1) * t41 - qJD(2) * t236;
t237 = t553 * t415;
t421 = (t526 - t494) * t415 + t527;
t43 = t517 + t421;
t473 = -qJD(1) * t43 - qJD(2) * t237;
t287 = t689 * t415;
t420 = (t526 + t494) * t415 + t527;
t45 = t516 + t420;
t472 = qJD(1) * t45 - qJD(2) * t287;
t471 = qJD(1) * t51 - qJD(2) * t689;
t470 = t501 * t415;
t449 = -t263 * t667 + t544 * t690;
t135 = (t676 + t449) * pkin(3);
t448 = t368 * t668 + t371 * t544;
t254 = (t664 + t448) * pkin(3);
t469 = qJD(1) * t135 + qJD(2) * t254;
t468 = qJD(1) * t151 + qJD(2) * t280;
t152 = t509 * t412;
t277 = t412 * t368;
t110 = qJD(1) * t152 + qJD(2) * t277;
t155 = 0.2e1 * t695;
t282 = t415 * t368;
t467 = -qJD(1) * t155 - qJD(2) * t282;
t233 = t359 ^ 2 - t361 ^ 2;
t466 = qJD(1) * t233 + qJD(2) * t199;
t465 = qJD(1) * t199 + qJD(2) * t393;
t260 = -t620 / 0.2e1 + t338 / 0.2e1;
t464 = qJD(1) * t260 + qJD(2) * t365;
t463 = t573 + t578;
t461 = t410 * t495;
t426 = pkin(2) * t359 / 0.2e1 + t628 / 0.2e1 - pkin(8) * t548 / 0.2e1;
t142 = t426 + t507;
t460 = pkin(2) * t571 - qJD(1) * t142;
t451 = pkin(2) * t676 + pkin(8) * t547 / 0.2e1;
t144 = -t348 / 0.2e1 + (t677 - t362 / 0.2e1) * t413 + t451;
t459 = pkin(2) * t573 - qJD(1) * t144;
t457 = t368 * t671 + t371 * t670;
t261 = t359 * t664 + t361 * t416 / 0.2e1;
t455 = -qJD(2) * t261 + t359 * t578;
t454 = qJD(1) * t261 + t413 * t571;
t53 = (-t611 + t634) * t371;
t80 = t225 ^ 2 - t227 ^ 2;
t447 = qJD(1) * t80 - qJD(2) * t53 + qJD(3) * t65;
t446 = t285 / 0.2e1 + t457;
t23 = t412 * t490 + t415 * t543 + t523;
t97 = t412 * t446 + t415 * t508;
t443 = -qJD(1) * t23 - qJD(2) * t97 - t400 * t567;
t25 = t412 * t543 - t415 * t490 + t522;
t99 = t412 * t508 - t415 * t446;
t442 = -qJD(1) * t25 - qJD(2) * t99 - t400 * t569;
t106 = t456 * t371;
t124 = -t611 / 0.2e1 + t634 / 0.2e1;
t440 = qJD(2) * t106 - qJD(3) * t124 + t225 * t581;
t276 = (t404 / 0.2e1 - t406 / 0.2e1) * t371;
t439 = qJD(1) * t124 - qJD(2) * t276 + t535;
t134 = t690 * t674 + t632 / 0.2e1;
t438 = qJD(2) * t134 + qJD(5) * t260 + t580 * t690;
t437 = qJD(1) * t134 + qJD(5) * t365 + t541;
t286 = t392 * t367;
t433 = qJD(1) * t53 + qJD(2) * t286 + t493;
t432 = -qJD(1) * t65 + t687;
t431 = qJD(1) * t106 + qJD(3) * t276 + t367 * t537;
t385 = qJD(2) * t520;
t358 = t365 * qJD(3);
t339 = (t532 - t554 / 0.2e1) * t414;
t323 = -0.2e1 * t371 * t530;
t300 = t595 / 0.2e1;
t272 = t277 * qJD(5);
t271 = t276 * qJD(5);
t256 = t261 * qJD(3);
t253 = pkin(3) * t448 + t551;
t252 = t690 * t567;
t196 = t199 * qJD(3);
t195 = (qJD(1) * t690 + t574) * t415;
t167 = qJD(2) * t310 + qJD(3) * t260;
t153 = t263 * t666 + t518;
t147 = t153 * qJD(5);
t146 = t152 * qJD(5);
t145 = t413 * t677 + t605 / 0.2e1 + t348 / 0.2e1 + t451;
t143 = t426 - t507;
t136 = pkin(3) * t449 + t552;
t125 = t134 * qJD(3);
t121 = t124 * qJD(5);
t109 = -t110 - t562;
t103 = t106 * qJD(5);
t100 = t305 * t665 + t608 / 0.2e1 + t597 / 0.2e1 - t457 * t415;
t98 = t305 * t663 + t596 / 0.2e1 - t609 / 0.2e1 + t457 * t412;
t78 = t453 + t491;
t75 = t263 * t515 + t529 + t594 / 0.2e1 + t412 * t520;
t66 = t631 / 0.2e1 + t445;
t64 = t65 * qJD(5);
t52 = t53 * qJD(5);
t50 = t51 * qJD(4);
t44 = t517 + t420;
t42 = t516 + t421;
t40 = t300 + t430;
t38 = t300 + t429;
t34 = t434 + t458;
t28 = -t425 + t489;
t26 = t399 * t696 + t522 + t546 - t615 / 0.2e1 + t602 / 0.2e1;
t24 = t412 * t524 + t523 + t545 - t604 / 0.2e1 - t613 / 0.2e1;
t16 = t95 * t525 + t519 + t601 / 0.2e1 - t435;
t15 = t95 * t515 + t513 - t612 / 0.2e1 - t436;
t10 = t101 * t674 + t102 * t673 + t107 * t672 + t108 * t675 - t423 + t441;
t8 = pkin(3) * t450 + t274 * t551 + t401 * t552 - t424;
t4 = t107 * t525 + t37 * t673 - t418 + t428;
t3 = t107 * t514 + t36 * t673 - t419 + t427;
t30 = [0, 0, 0, t500, t370 * qJD(2), t411 * t538, -t411 * t531, 0, -t308 * qJD(2), -t309 * qJD(2), (-t359 * qJD(3) + t416 * t538) * t361, qJD(2) * t232 + qJD(3) * t233, -t289 * qJD(2) + t359 * t536, t288 * qJD(2) + t361 * t536, -t500, -qJD(2) * t85 - qJD(3) * t149, qJD(2) * t86 + qJD(3) * t148, qJD(2) * t18 + qJD(3) * t17 + qJD(4) * t79, qJD(2) * t20 + qJD(3) * t19 + qJD(4) * t29, (qJD(2) * t270 - qJD(5) * t225 - t263 * t567) * t227, qJD(2) * t57 - qJD(3) * t56 + qJD(5) * t80, qJD(2) * t88 + qJD(3) * t60 - t225 * t564, qJD(2) * t87 - qJD(3) * t59 - t227 * t564, (qJD(2) * t313 + qJD(3) * t690) * t263, qJD(2) * t11 + qJD(3) * t5 - qJD(4) * t58 + qJD(5) * t22, qJD(2) * t12 + qJD(3) * t6 + qJD(4) * t61 + qJD(5) * t21; 0, 0, 0, t499, t555, t503 * t616, -t503 * t617, 0, -qJD(2) * t364 - t557, qJD(2) * t362 - t556, t463 * t547 + t256, t393 * t538 + t196 + t560, t413 * t531 - t558, t416 * t531 + t559, -t339, -t590 + (-t364 * t416 + t413 * t461) * qJD(2) + t145 * qJD(3), t589 + (t364 * t413 + t416 * t461) * qJD(2) + t143 * qJD(3), t639 + (-t132 * t371 - t133 * t368 + t305 * t314 - t313 * t444) * qJD(2) + t10 * qJD(3) + t50, t637 + (-t132 * t305 + t133 * t444 + t311 * t401) * qJD(2) + t8 * qJD(3) + t28 * qJD(4), qJD(3) * t66 - t103 + (t539 + t581) * t270, t591 + t34 * qJD(3) - t52 + (-t269 * t415 - t631) * t574, t587 + (t270 * t368 + t371 * t595) * qJD(2) + t42 * qJD(3) + t75 * qJD(5), t588 + (-t269 * t368 - t280 * t313) * qJD(2) + t38 * qJD(3) + t78 * qJD(5), t563 + t125 + (t575 + t580) * t313, t642 + (t119 * t280 - t178 * t313 + t269 * t305 + t368 * t54) * qJD(2) + t3 * qJD(3) + t40 * qJD(4) + t16 * qJD(5), t641 + (t119 * t624 - t179 * t313 + t270 * t305 - t368 * t55) * qJD(2) + t4 * qJD(3) + t44 * qJD(4) + t15 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t455, t466, t488 * t359, t488 * t361, t385, qJD(2) * t145 - qJD(3) * t239 - t583, qJD(2) * t143 + qJD(3) * t238 + t584, t640 + t10 * qJD(2) + (t263 * t655 - t409 * t690) * t658, t638 + t8 * qJD(2) + (-t107 * t655 + t108 * t409) * t658 + t136 * qJD(4), qJD(2) * t66 + t121 - (t569 + t581) * t693, t34 * qJD(2) - t263 * t688 - t592 + t64, qJD(2) * t42 + t569 * t690 + t649 + t698, qJD(2) * t38 + t147 + t252 - t650, t438, t657 + t3 * qJD(2) + (t412 * t480 - t643) * qJD(3) + t26 * qJD(5), t656 + t4 * qJD(2) + (t415 * t480 + t644) * qJD(3) + t24 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t483, qJD(2) * t28 + qJD(3) * t136 + t652, 0, 0, 0, 0, 0, qJD(2) * t40 + t147 - t651, qJD(2) * t44 + t648 + t698; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t440, t447, qJD(2) * t75 + qJD(3) * t697 + t225 * t502, qJD(2) * t78 + qJD(3) * t153 + t227 * t502, t167, qJD(2) * t16 + qJD(3) * t26 + qJD(4) * t153 - qJD(5) * t37 + t653, qJD(2) * t15 + qJD(3) * t24 + qJD(4) * t697 + qJD(5) * t36 + t654; 0, 0, 0, -t499, -t555, -t411 * t542, t577 * t617, 0, t557, t556, -t361 * t497 + t256, t196 - t560, -t416 * t536 + t558, t413 * t536 - t559, t339, qJD(3) * t144 + t590, qJD(3) * t142 - t589, -qJD(3) * t9 + t50 - t639, -qJD(3) * t7 - qJD(4) * t27 - t637, qJD(3) * t67 - t270 * t581 - t103, qJD(3) * t35 - t52 - t591, qJD(3) * t43 + qJD(5) * t76 - t587, qJD(3) * t39 + qJD(5) * t77 - t588, -t313 * t580 + t125 - t563, -qJD(3) * t1 + qJD(4) * t41 - qJD(5) * t13 - t642, -qJD(3) * t2 + qJD(4) * t45 - qJD(5) * t14 - t641; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t413 * t566, t393 * qJD(3), 0, 0, 0, -pkin(2) * t568, -pkin(2) * t566, -qJD(4) * t689, qJD(3) * t118 + qJD(4) * t158, -t367 * t530 - t406 * t533, -qJD(5) * t286 + t368 * t493, qJD(3) * t237 - t562 * t626, -qJD(3) * t235 - t371 * t534, t533, qJD(3) * t48 - qJD(4) * t236 + qJD(5) * t112, qJD(3) * t49 - qJD(4) * t287 + qJD(5) * t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t454, t465, -t488 * t416, t488 * t413, qJD(1) * t521, -pkin(8) * t566 - t459, pkin(8) * t568 - t460, (t368 * t655 - t371 * t409) * t658 - t485, (-t305 * t409 - t444 * t655) * t658 + t253 * qJD(4) + t484, t647 - t271 + (-t535 - t540) * t368, t368 * t687 + t323 + t627, t371 * t569 - t473, t357 - t475, t437, (t412 * t479 - t629) * qJD(3) + t100 * qJD(5) + t487, (t415 * t479 + t630) * qJD(3) + t98 * qJD(5) + t486; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t471, qJD(3) * t253 + t476, 0, 0, 0, 0, 0, t474, t472; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t431, -t433, t280 * t501 + t646, t371 * t470 + t645, t358 - t579, qJD(3) * t100 - qJD(5) * t179 - t477, qJD(3) * t98 + qJD(5) * t178 + t478; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t455, -t466, (-qJD(1) * t359 + t571) * t616, -t463 * t616, t385, -qJD(2) * t144 + t583, -qJD(2) * t142 - t584, qJD(2) * t9 - t640, qJD(2) * t7 + qJD(4) * t135 - t638, -qJD(2) * t67 + t581 * t693 + t121, -qJD(2) * t35 + t592 + t64, -qJD(2) * t43 + qJD(5) * t155 - t649, -qJD(2) * t39 - t146 + t650, -t438, qJD(2) * t1 + qJD(5) * t25 - t565 * t690 - t657, qJD(2) * t2 + qJD(4) * t151 + qJD(5) * t23 - t656; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t454, -t465, t497, -t413 * t542, qJD(1) * t520, t459, t460, t485, qJD(4) * t254 - t484, t368 * t540 - t271 - t647, -0.2e1 * t368 * t498 + t323 - t627, qJD(5) * t282 + t473, -t272 + t475, -t437, qJD(5) * t99 - t371 * t565 - t487, qJD(4) * t280 + qJD(5) * t97 - t486; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t530, t392 * qJD(5), 0, 0, 0, t400 * t562, t400 * t561; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t469, 0, 0, 0, 0, 0, -t195, t468; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t439, -t432, -t467 + t561, t109, -t464, -t399 * t561 - t442, t399 * t562 - t443; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t483, qJD(2) * t27 - qJD(3) * t135 - t652, 0, 0, 0, 0, 0, -qJD(2) * t41 - t146 + t252 + t651, -qJD(2) * t45 - qJD(3) * t151 + qJD(5) * t156 - t648; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t471, -qJD(3) * t254 - t476, 0, 0, 0, 0, 0, -t272 + t357 - t474, -qJD(3) * t280 - t472 - t534; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t469, 0, 0, 0, 0, 0, t195, -t468; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, t470 + t582; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t440, -t447, -qJD(2) * t76 - qJD(3) * t155 + t225 * t580, -qJD(2) * t77 + qJD(3) * t152 + t227 * t580, t167, qJD(2) * t13 - qJD(3) * t25 + qJD(4) * t152 - t653, qJD(2) * t14 - qJD(3) * t23 - qJD(4) * t156 - t654; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t431, t433, -qJD(3) * t282 + t412 * t541 - t646, qJD(3) * t277 + t368 * t539 - t645, t358 + t579, -qJD(3) * t99 + qJD(4) * t277 + t477, -qJD(3) * t97 + t368 * t565 - t478; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t439, t432, t467, t110, t464, t442, t443; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, t368 * t572 - t582; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t30;