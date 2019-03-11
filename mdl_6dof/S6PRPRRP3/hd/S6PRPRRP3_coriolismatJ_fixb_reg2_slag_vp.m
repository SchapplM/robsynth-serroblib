% Calculate inertial parameters regressor of coriolis matrix for
% S6PRPRRP3
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
% Datum: 2019-03-08 20:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6PRPRRP3_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:07:52
% EndTime: 2019-03-08 20:08:08
% DurationCPUTime: 11.36s
% Computational Cost: add. (10840->545), mult. (24983->713), div. (0->0), fcn. (28667->10), ass. (0->450)
t479 = sin(pkin(11));
t713 = cos(qJ(4));
t576 = t713 * t479;
t481 = cos(pkin(11));
t483 = sin(qJ(4));
t639 = t483 * t481;
t447 = t576 + t639;
t480 = sin(pkin(6));
t485 = cos(qJ(2));
t646 = t480 * t485;
t387 = t447 * t646;
t720 = -t387 / 0.2e1;
t558 = t646 / 0.2e1;
t622 = t447 * t558;
t484 = cos(qJ(5));
t712 = sin(qJ(2));
t579 = t480 * t712;
t575 = t713 * t481;
t640 = t483 * t479;
t507 = t575 - t640;
t388 = t507 * t646;
t482 = sin(qJ(5));
t642 = t482 * t388;
t320 = t484 * t579 - t642;
t543 = t482 * t579;
t635 = t484 * t388;
t321 = t543 + t635;
t715 = -t482 / 0.2e1;
t728 = t321 * t484 / 0.2e1 + t320 * t715;
t736 = pkin(4) * t720 + pkin(9) * t728;
t477 = t482 ^ 2;
t478 = t484 ^ 2;
t716 = -t478 / 0.2e1;
t735 = t716 - t477 / 0.2e1;
t641 = t482 * t484;
t580 = t447 * t641;
t547 = 0.2e1 * t580;
t696 = cos(pkin(6));
t499 = t479 * t579 - t481 * t696;
t734 = t499 * t479;
t472 = -pkin(3) * t481 - pkin(2);
t536 = -pkin(4) * t507 - pkin(9) * t447;
t503 = t472 + t536;
t706 = pkin(8) + qJ(3);
t550 = t706 * t479;
t538 = t483 * t550;
t458 = t706 * t481;
t577 = t713 * t458;
t372 = t577 - t538;
t643 = t482 * t372;
t219 = -t484 * t503 + t643;
t638 = t484 * qJ(6);
t190 = t447 * t638 + t219;
t709 = t507 * pkin(5);
t165 = -t190 - t709;
t733 = t165 + t190;
t587 = t447 * qJD(2);
t424 = t484 * t587;
t613 = qJD(4) * t482;
t393 = t424 + t613;
t442 = t507 ^ 2;
t443 = t447 ^ 2;
t732 = -t443 - t442;
t583 = t443 - t442;
t552 = t190 / 0.2e1 + t165 / 0.2e1;
t581 = -t709 / 0.2e1;
t518 = t581 + t552;
t26 = t518 * t482;
t495 = t483 * t499;
t430 = t479 * t696 + t481 * t579;
t578 = t713 * t430;
t317 = t578 - t495;
t637 = t484 * t317;
t274 = -t482 * t646 + t637;
t675 = t274 * t484;
t260 = t675 / 0.2e1;
t565 = -t675 / 0.2e1;
t157 = t260 + t565;
t602 = t157 * qJD(1);
t731 = qJD(2) * t26 - t602;
t644 = t482 * t317;
t273 = t484 * t646 + t644;
t677 = t273 * t482;
t510 = t565 - t677 / 0.2e1;
t316 = t430 * t483 + t499 * t713;
t668 = t316 * t447;
t563 = t668 / 0.2e1;
t492 = -t507 * t510 + t563;
t664 = t321 * t482;
t665 = t320 * t484;
t509 = -t665 / 0.2e1 - t664 / 0.2e1;
t52 = t492 + t509;
t634 = t52 * qJD(1);
t645 = t482 * t219;
t371 = t458 * t483 + t550 * t713;
t660 = t371 * t447;
t636 = t484 * t372;
t220 = t482 * t503 + t636;
t681 = t220 * t484;
t87 = t660 - (-t645 - t681) * t507;
t730 = -qJD(2) * t87 - t634;
t342 = t482 * t447;
t305 = pkin(5) * t342 + t371;
t674 = t305 * t447;
t191 = -qJ(6) * t342 + t220;
t684 = t191 * t484;
t687 = t165 * t482;
t45 = t674 - (-t684 + t687) * t507;
t729 = -qJD(2) * t45 - t634;
t438 = t576 / 0.2e1 + t639 / 0.2e1;
t707 = t447 * pkin(4);
t708 = t507 * pkin(9);
t363 = t707 - t708;
t353 = t484 * t363;
t659 = t371 * t482;
t238 = t353 + t659;
t345 = t484 * t507;
t437 = t447 * pkin(5);
t171 = -qJ(6) * t345 + t238 + t437;
t727 = t171 / 0.2e1;
t726 = -t273 / 0.2e1;
t725 = -t274 / 0.2e1;
t724 = t274 / 0.2e1;
t723 = t317 / 0.2e1;
t722 = t320 / 0.2e1;
t721 = -t353 / 0.2e1;
t360 = t371 * t484;
t355 = t360 / 0.2e1;
t705 = pkin(9) + qJ(6);
t459 = t705 * t482;
t719 = t459 / 0.2e1;
t460 = t705 * t484;
t718 = -t460 / 0.2e1;
t714 = -t484 / 0.2e1;
t711 = pkin(5) * t443;
t51 = t492 - t509;
t669 = t316 * t387;
t75 = -t273 * t320 + t274 * t321 + t669;
t628 = t75 * qJD(1);
t704 = qJD(3) * t51 + t628;
t703 = qJD(3) * t52 - t628;
t514 = (t725 + t637 / 0.2e1) * t447;
t657 = t387 * t482;
t64 = -t657 / 0.2e1 + t514;
t631 = t64 * qJD(1);
t700 = qJD(3) * t342 - t631;
t431 = t477 * t507;
t432 = t478 * t507;
t348 = -t431 - t432;
t511 = t273 * t714 + t482 * t724;
t502 = t511 * t507;
t82 = -t502 - t728;
t627 = t82 * qJD(1);
t698 = -qJD(3) * t348 - t627;
t340 = t482 * t507;
t537 = t579 / 0.2e1;
t505 = t537 - t668 / 0.2e1;
t556 = -t642 / 0.2e1;
t676 = t274 * t507;
t103 = t556 - t676 / 0.2e1 + t505 * t484;
t605 = t103 * qJD(1);
t697 = -qJD(3) * t340 + t605;
t16 = t518 * t484;
t695 = qJD(2) * t16;
t21 = t733 * t342;
t694 = qJD(2) * t21;
t473 = -pkin(5) * t484 - pkin(4);
t562 = t447 * t473 / 0.2e1;
t649 = t460 * t484;
t652 = t459 * t482;
t491 = -(-t649 / 0.2e1 - t652 / 0.2e1) * t507 + t562;
t352 = t482 * t363;
t239 = -t360 + t352;
t201 = -qJ(6) * t340 + t239;
t682 = t201 * t482;
t685 = t171 * t484;
t513 = t685 / 0.2e1 + t682 / 0.2e1;
t65 = t491 - t513;
t691 = qJD(2) * t65;
t93 = t191 * t507 + (t482 * t711 + t674) * t484;
t689 = qJD(2) * t93;
t94 = -t190 * t507 - t305 * t342 + t478 * t711;
t688 = qJD(2) * t94;
t686 = t165 * t484;
t525 = t191 * t482 + t686;
t20 = (t682 + t685) * t447 + t525 * t507;
t683 = t20 * qJD(2);
t680 = t238 * t484;
t679 = t239 * t482;
t678 = t273 * t507;
t673 = t305 * t482;
t672 = t305 * t484;
t306 = pkin(5) * t340 + t372;
t671 = t306 * t482;
t670 = t306 * t484;
t188 = t316 * t482;
t667 = t317 * t507;
t33 = (t165 + t671) * t447 - (t171 - t673) * t507;
t662 = t33 * qJD(2);
t37 = (t679 + t680) * t447 - (t219 * t484 - t220 * t482) * t507;
t661 = t37 * qJD(2);
t658 = t387 * t371;
t656 = t387 * t484;
t39 = (-t191 + t670) * t447 - (-t201 - t672) * t507;
t655 = t39 * qJD(2);
t654 = t430 * t481;
t653 = t447 * t484;
t651 = t459 * t484;
t650 = t460 * t482;
t648 = t473 * t482;
t647 = t473 * t484;
t56 = (t317 - t675 - t677) * t316;
t633 = t56 * qJD(1);
t515 = (t726 + t644 / 0.2e1) * t447;
t61 = t656 / 0.2e1 + t515;
t632 = t61 * qJD(1);
t71 = (-t219 + t643) * t447 - (t238 - t659) * t507;
t630 = t71 * qJD(2);
t72 = (-t220 + t636) * t447 - (-t239 - t360) * t507;
t629 = t72 * qJD(2);
t498 = -t735 * t708 - t707 / 0.2e1;
t512 = t680 / 0.2e1 + t679 / 0.2e1;
t91 = t498 - t512;
t626 = t91 * qJD(2);
t625 = -t687 / 0.2e1 + t684 / 0.2e1;
t623 = -t352 / 0.2e1 + t355;
t464 = t479 ^ 2 + t481 ^ 2;
t468 = t477 + t478;
t469 = t478 - t477;
t133 = -t219 * t507 - t342 * t371;
t620 = qJD(2) * t133;
t134 = t220 * t507 + t371 * t653;
t619 = qJD(2) * t134;
t288 = t583 * t482;
t618 = qJD(2) * t288;
t289 = t732 * t482;
t617 = qJD(2) * t289;
t290 = t583 * t484;
t616 = qJD(2) * t290;
t351 = t732 * t484;
t210 = qJD(2) * t351;
t615 = qJD(2) * t480;
t614 = qJD(3) * t484;
t612 = qJD(4) * t484;
t611 = qJD(5) * t191;
t610 = qJD(5) * t274;
t609 = qJD(5) * t460;
t608 = qJD(5) * t482;
t474 = qJD(5) * t484;
t607 = qJD(6) * t482;
t606 = qJD(6) * t484;
t554 = -t635 / 0.2e1;
t104 = t554 + t678 / 0.2e1 - t505 * t482;
t604 = t104 * qJD(1);
t463 = t480 ^ 2 * t712 * t485;
t126 = t317 * t388 - t463 + t669;
t603 = t126 * qJD(1);
t601 = t157 * qJD(3);
t244 = -t463 + (t654 + t734) * t646;
t600 = t244 * qJD(1);
t540 = t735 * t447;
t284 = t540 - t438;
t599 = t284 * qJD(2);
t504 = -t575 / 0.2e1 + t640 / 0.2e1;
t292 = (t507 / 0.2e1 + t504) * t646;
t597 = t292 * qJD(1);
t596 = t583 * qJD(2);
t551 = t477 / 0.2e1 + t716;
t339 = t551 * t447;
t595 = t339 * qJD(5);
t594 = t340 * qJD(2);
t593 = t342 * qJD(2);
t333 = t345 * qJD(2);
t592 = t348 * qJD(2);
t349 = t468 * t443;
t591 = t349 * qJD(2);
t590 = t732 * qJD(2);
t589 = t438 * qJD(2);
t588 = t507 * qJD(2);
t434 = t507 * qJD(4);
t586 = t447 * qJD(4);
t585 = t464 * qJD(2);
t584 = t468 * qJD(4);
t582 = pkin(5) * t608;
t574 = t482 * t612;
t573 = t447 * t608;
t572 = t447 * t474;
t571 = t447 * t607;
t570 = t447 * t606;
t569 = t507 * t587;
t370 = t507 * t586;
t470 = t482 * t474;
t568 = t484 * t588;
t567 = t507 * t606;
t566 = -t684 / 0.2e1;
t564 = t673 / 0.2e1;
t354 = t659 / 0.2e1;
t561 = -t653 / 0.2e1;
t560 = t653 / 0.2e1;
t559 = t647 / 0.2e1;
t557 = qJ(6) * t715;
t555 = -t342 / 0.2e1;
t553 = -t345 / 0.2e1;
t549 = t464 * t485;
t548 = t468 * t316;
t359 = t507 * t424;
t233 = -qJD(4) * t340 - t359;
t546 = t482 * t581;
t545 = qJD(2) * t472 + qJD(3);
t542 = t443 * t470;
t541 = t316 * t560;
t539 = qJD(2) * t579;
t535 = -t447 * t614 - t632;
t534 = qJD(4) * t547;
t533 = -t507 * t614 + t604;
t532 = -t359 + t572;
t15 = t165 * t171 + t191 * t201 + t305 * t306;
t487 = (t566 + t687 / 0.2e1 + t306 / 0.2e1) * t316 + t171 * t726 + t201 * t724 + t305 * t723;
t497 = t320 * t719 + t321 * t718 + t473 * t720;
t2 = t487 + t497;
t531 = qJD(1) * t2 + qJD(2) * t15;
t23 = pkin(5) * t305 * t653 - t191 * t733;
t6 = t552 * t274 + (t316 * t561 + t722) * pkin(5);
t530 = -qJD(1) * t6 + qJD(2) * t23;
t38 = -t219 * t238 + t220 * t239 + t371 * t372;
t486 = (-t681 / 0.2e1 - t645 / 0.2e1 + t372 / 0.2e1) * t316 + t238 * t726 + t239 * t724 + t371 * t723;
t5 = t486 - t736;
t529 = t5 * qJD(1) + t38 * qJD(2);
t80 = t525 * t447;
t96 = t447 * t511 + t622;
t527 = -qJD(1) * t96 - qJD(2) * t80;
t278 = pkin(5) * t648;
t8 = t552 * t460 + (t727 + t473 * t561 - t673 / 0.2e1) * pkin(5);
t526 = -qJD(2) * t8 + qJD(4) * t278;
t524 = -t238 * t482 + t239 * t484;
t382 = t649 + t652;
t429 = pkin(5) * t641 - t648;
t496 = (pkin(5) * t551 + t559) * t447;
t77 = -t437 + t721 + (t305 / 0.2e1 - t371 / 0.2e1) * t482 - (t718 - t638 / 0.2e1) * t507 + t496;
t523 = qJD(2) * t77 - qJD(4) * t429;
t441 = t477 * pkin(5) + t647;
t500 = -pkin(5) * t580 - t672 / 0.2e1 + t482 * t562;
t89 = -(t557 - t459 / 0.2e1) * t507 + t500 + t623;
t522 = -qJD(2) * t89 + qJD(4) * t441;
t143 = -t667 / 0.2e1 + t505;
t200 = t372 * t507 + t660;
t521 = qJD(1) * t143 - qJD(2) * t200;
t489 = t654 / 0.2e1 + t734 / 0.2e1;
t281 = t537 - t489;
t455 = t464 * qJ(3);
t520 = qJD(1) * t281 - qJD(2) * t455;
t519 = -t708 / 0.2e1 + t707 / 0.2e1;
t493 = t482 * t519 + t355;
t128 = -t360 / 0.2e1 + t352 / 0.2e1 + t493;
t517 = pkin(4) * t612 - qJD(2) * t128;
t506 = t519 * t484;
t130 = t721 - t506;
t516 = pkin(4) * t613 - qJD(2) * t130;
t195 = (-qJD(5) + t588) * t342;
t300 = -qJD(2) * t339 + t574;
t391 = t482 * t587 - t612;
t295 = qJD(5) * t438 - t569;
t248 = qJD(2) * t443 * t641 + qJD(4) * t339;
t350 = t469 * t443;
t265 = qJD(2) * t350 + t534;
t361 = qJD(2) * t547 - qJD(4) * t469;
t488 = -t495 / 0.2e1 + t578 / 0.2e1;
t109 = t488 + t510;
t490 = -t577 / 0.2e1 + t546 + t538 / 0.2e1;
t494 = (-t650 / 0.2e1 + t651 / 0.2e1) * t447 + t625;
t43 = t490 + t494;
t501 = -qJD(1) * t109 + qJD(2) * t43 + qJD(4) * t382;
t461 = t469 * qJD(5);
t426 = t438 * qJD(4);
t423 = t484 * t586;
t392 = -t474 + t568;
t383 = t393 * pkin(5);
t347 = t353 / 0.2e1;
t338 = t351 * qJD(3);
t334 = t348 * qJD(4);
t324 = t340 * qJD(5);
t310 = t474 - t333;
t309 = t594 - t608;
t294 = t720 - t622;
t293 = t504 * t646 - t507 * t558;
t285 = t289 * qJD(3);
t283 = t540 + t438;
t282 = t537 + t489;
t276 = t370 * t478 - t542;
t275 = t370 * t477 + t542;
t263 = 0.2e1 * t484 * t195;
t232 = -t478 * t569 - t595;
t231 = qJD(4) * t345 - t482 * t569;
t230 = -t477 * t569 + t595;
t221 = -qJD(5) * t350 - t507 * t534;
t203 = qJD(4) * t290 + t507 * t573;
t202 = -qJD(4) * t288 + t507 * t572;
t199 = -qJD(5) * t345 - t616;
t197 = t324 + t618;
t189 = t316 * t484;
t185 = -t595 - (-t478 * t587 - t574) * t507;
t184 = t595 - (-t477 * t587 + t574) * t507;
t183 = (-t431 + t432) * qJD(4) + (-qJD(5) - t588) * t547;
t178 = -qJD(4) * t342 + t474 * t507 + t210;
t164 = t482 * t586 + t616;
t163 = t423 - t618;
t162 = t324 + t423 + t617;
t148 = t157 * qJD(5);
t144 = t667 / 0.2e1 + t563 + t537;
t131 = 0.2e1 * t354 + t347 - t506;
t129 = t493 + t623;
t110 = t260 + t677 / 0.2e1 + t488;
t106 = t676 / 0.2e1 + t541 + t556 + t484 * t537;
t105 = -t678 / 0.2e1 + t316 * t555 + t554 - t543 / 0.2e1;
t97 = t273 * t560 + t274 * t555 + t622;
t92 = t498 + t512;
t90 = -t500 + t623 - (t557 + t719) * t507;
t81 = -t502 + t728;
t76 = qJ(6) * t553 - t507 * t718 + t347 + t354 + t437 + t496 + t564;
t70 = qJD(2) * t104;
t69 = qJD(2) * t103;
t66 = t491 + t513;
t63 = t657 / 0.2e1 + t514;
t62 = -t656 / 0.2e1 + t515;
t42 = -t490 + t494;
t41 = qJD(2) * t106 + qJD(4) * t188 - t610;
t40 = qJD(2) * t105 + qJD(4) * t189 + qJD(5) * t273;
t36 = qJD(2) * t82;
t35 = pkin(5) * t188;
t32 = qJD(2) * t61;
t31 = qJD(2) * t64;
t30 = qJD(4) * t82;
t29 = -qJD(2) * t52 + t148;
t28 = qJD(2) * t51 + t148;
t27 = t190 * t715 + t546 + t566 + t625;
t25 = qJD(2) * t63 + qJD(5) * t189 + t317 * t613;
t24 = qJD(2) * t62 + qJD(5) * t188 - t317 * t612;
t22 = t81 * qJD(2) - qJD(4) * t548;
t19 = qJD(4) * t64 - qJD(5) * t104;
t18 = qJD(4) * t61 - qJD(5) * t103;
t17 = t190 * t714 - t686 / 0.2e1 + pkin(5) * t553;
t14 = t81 * qJD(4) + (-t664 - t665) * t587;
t11 = (t321 * t507 + t387 * t653) * qJD(2) + t63 * qJD(4) + t105 * qJD(5);
t10 = (-t320 * t507 + t342 * t387) * qJD(2) + t62 * qJD(4) + t106 * qJD(5);
t9 = t559 * t437 + t733 * t718 + (t564 + t727) * pkin(5);
t7 = t733 * t725 + (t541 + t722) * pkin(5);
t4 = t486 + t736;
t3 = qJD(2) * t75 + qJD(4) * t56;
t1 = t487 - t497;
t12 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t126, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t539, -t485 * t615, 0, 0, 0, 0, 0, 0, 0, 0, -t481 * t539, t479 * t539, t549 * t615, t600 + t282 * qJD(3) + (-pkin(2) * t712 + qJ(3) * t549) * t615, 0, 0, 0, 0, 0, 0, qJD(4) * t294 - t507 * t539, qJD(4) * t293 + t447 * t539 (t387 * t447 + t388 * t507) * qJD(2), t603 + (t372 * t388 + t472 * t579 + t658) * qJD(2) + t144 * qJD(3), 0, 0, 0, 0, 0, 0, t10, t11, t14 (-t219 * t320 + t220 * t321 + t658) * qJD(2) + t4 * qJD(4) + t704, 0, 0, 0, 0, 0, 0, t10, t11, t14 (t165 * t320 + t191 * t321 + t305 * t387) * qJD(2) + t1 * qJD(4) + t7 * qJD(5) + t97 * qJD(6) + t704; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t282 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t144 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t294 - qJD(4) * t317, qJD(2) * t293 + qJD(4) * t316, 0, 0, 0, 0, 0, 0, 0, 0, t24, t25, t22, t633 + t4 * qJD(2) + (-t317 * pkin(4) - pkin(9) * t548) * qJD(4), 0, 0, 0, 0, 0, 0, t24, t25, t22, t633 + t1 * qJD(2) + (-t316 * t382 + t317 * t473) * qJD(4) + t35 * qJD(5) + t110 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t40, 0, t601, 0, 0, 0, 0, 0, 0, t41, t40, 0, -pkin(5) * t610 + qJD(2) * t7 + qJD(4) * t35 + t601; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t97 + qJD(4) * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t281 - t600, 0, 0, 0, 0, 0, 0, 0, -t292 * qJD(4), 0, -qJD(3) * t143 - t603, 0, 0, 0, 0, 0, 0, t18, t19, t30, qJD(4) * t5 + t703, 0, 0, 0, 0, 0, 0, t18, t19, t30, qJD(4) * t2 - qJD(5) * t6 - qJD(6) * t96 + t703; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t464 * qJD(3), t455 * qJD(3), t370, -t583 * qJD(4), 0, -t370, 0, 0, t472 * t586, t472 * t434, -t732 * qJD(3), qJD(3) * t200, t276, t221, t203, t275, t202, -t370, qJD(4) * t71 + qJD(5) * t134 - t285, qJD(4) * t72 + qJD(5) * t133 - t338, -qJD(4) * t37, qJD(3) * t87 + qJD(4) * t38, t276, t221, t203, t275, t202, -t370, qJD(4) * t33 + qJD(5) * t93 + t447 * t567 - t285, qJD(4) * t39 + qJD(5) * t94 - t507 * t571 - t338, -qJD(4) * t20 + qJD(5) * t21 + qJD(6) * t349, qJD(3) * t45 + qJD(4) * t15 + qJD(5) * t23 - qJD(6) * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t585, -t520, 0, 0, 0, 0, 0, 0, 0, 0, -t590, -t521, 0, 0, 0, 0, 0, 0, -t617, -t210, 0, qJD(4) * t92 - t730, 0, 0, 0, 0, 0, 0, -t617, -t210, 0, qJD(4) * t66 + qJD(5) * t27 + qJD(6) * t283 - t729; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t569, -t596, t434, -t569, -t586, 0, -qJD(4) * t372 + t472 * t587, qJD(4) * t371 + t472 * t588 - t597, 0, 0, t185, t183, t164, t184, t163, t295, t632 + t630 + (t482 * t536 - t636) * qJD(4) + t131 * qJD(5), t629 + (t484 * t536 + t643) * qJD(4) + t129 * qJD(5) + t631, qJD(4) * t524 + t627 - t661, t92 * qJD(3) + (-t372 * pkin(4) + pkin(9) * t524) * qJD(4) + t529, t185, t183, t164, t184, t163, t295, t632 + t662 + (t340 * t473 - t447 * t459 - t670) * qJD(4) + t76 * qJD(5) + t340 * qJD(6), t655 + (t345 * t473 - t460 * t447 + t671) * qJD(4) + t90 * qJD(5) + t567 + t631, t627 - t683 + (-t171 * t482 + t201 * t484 - (t650 - t651) * t507) * qJD(4) + t17 * qJD(5), t66 * qJD(3) + (-t171 * t459 + t201 * t460 + t306 * t473) * qJD(4) + t9 * qJD(5) + t42 * qJD(6) + t531; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t248, -t265, t195, t248, -t532, t426, qJD(4) * t131 - qJD(5) * t220 - t605 + t619, qJD(4) * t129 + qJD(5) * t219 - t604 + t620, 0, 0, -t248, -t265, t195, t248, -t532, t426, qJD(4) * t76 - t605 - t611 + t689, qJD(4) * t90 + qJD(5) * t190 - t604 + t688, pkin(5) * t573 + qJD(4) * t17 + t694, -pkin(5) * t611 + qJD(3) * t27 + qJD(4) * t9 + t530; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t233, -t391 * t507, t591, qJD(3) * t283 + qJD(4) * t42 + t527; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t281 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t143 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t585, t520, 0, 0, 0, 0, 0, 0, t586, t434, t590, t521, 0, 0, 0, 0, 0, 0, t162, t178, t334, -qJD(4) * t91 + t730, 0, 0, 0, 0, 0, 0, t162, t178, t334, -qJD(4) * t65 - qJD(5) * t26 + qJD(6) * t284 + t729; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t587, t588, 0, 0, 0, 0, 0, 0, 0, 0, t424, -t593, t592, -t626, 0, 0, 0, 0, 0, 0, t424, -t593, t592, -t691; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t309, t392, 0, t602, 0, 0, 0, 0, 0, 0, t309, t392, 0, -t582 - t731; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t599; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t292 * qJD(2), 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t31, -t36, -qJD(2) * t5 - t633, 0, 0, 0, 0, 0, 0, -t32, -t31, -t36, -qJD(2) * t2 - qJD(6) * t109 - t633; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t569, t596, 0, t569, 0, 0, -t447 * t545, -t507 * t545 + t597, 0, 0, t232, t263, t199, t230, t197, -t295, qJD(5) * t130 + t535 - t630, qJD(5) * t128 - t629 + t700, t661 + t698, qJD(3) * t91 - t529, t232, t263, t199, t230, t197, -t295, qJD(5) * t77 + t535 - t662, -qJD(5) * t89 - t655 + t700, -qJD(5) * t16 + t683 + t698, qJD(3) * t65 - qJD(5) * t8 + qJD(6) * t43 - t531; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t587, -t588, 0, 0, 0, 0, 0, 0, 0, 0, -t424, t593, -t592, t626, 0, 0, 0, 0, 0, 0, -t424, t593, -t592, t691; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t470, t461, 0, -t470, 0, 0, -pkin(4) * t608, -pkin(4) * t474, 0, 0, t470, t461, 0, -t470, 0, 0, -t429 * qJD(5), t441 * qJD(5), qJD(6) * t468, qJD(5) * t278 + qJD(6) * t382; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t300, -t361, t310, -t300, t309, -t589, -pkin(9) * t474 - t516, pkin(9) * t608 - t517, 0, 0, t300, -t361, t310, -t300, t309, -t589, t523 - t609, qJD(5) * t459 + t522, -pkin(5) * t474 - t695, -pkin(5) * t609 + t526; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t584, t501; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, t70, 0, -t601, 0, 0, 0, 0, 0, 0, t69, t70, 0, qJD(2) * t6 - t601; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t248, t265, t231, -t248, t233, t426, -qJD(4) * t130 - t619 + t697, -qJD(4) * t128 + t533 - t620, 0, 0, t248, t265, t231, -t248, t233, t426, -qJD(4) * t77 - t570 - t689 + t697, qJD(4) * t89 + t533 + t571 - t688, qJD(4) * t16 - t694, -pkin(5) * t570 + qJD(3) * t26 + qJD(4) * t8 - t530; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t594, -t568, 0, -t602, 0, 0, 0, 0, 0, 0, -t594, -t568, 0, t731; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t300, t361, t333, t300, -t594, t589, t516, t517, 0, 0, -t300, t361, t333, t300, -t594, t589, -t523 - t607, -t522 - t606, t695, -pkin(5) * t607 - t526; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t393, t391, 0, -t383; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t96 + qJD(4) * t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t532, t195, -t591, pkin(5) * t572 - qJD(3) * t284 - qJD(4) * t43 - t527; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t599; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t608, t474, -t584, -t501 + t582; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t393, -t391, 0, t383; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t12;
