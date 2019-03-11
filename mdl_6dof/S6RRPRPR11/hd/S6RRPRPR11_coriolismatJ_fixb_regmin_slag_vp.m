% Calculate minimal parameter regressor of coriolis matrix for
% S6RRPRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x30]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRPRPR11_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:15:07
% EndTime: 2019-03-09 11:15:35
% DurationCPUTime: 14.31s
% Computational Cost: add. (10529->565), mult. (20701->752), div. (0->0), fcn. (22227->8), ass. (0->448)
t751 = qJD(4) + qJD(6);
t425 = sin(qJ(4));
t430 = -pkin(2) - pkin(8);
t713 = -qJ(5) + t430;
t371 = t713 * t425;
t428 = cos(qJ(4));
t372 = t713 * t428;
t422 = sin(pkin(10));
t423 = cos(pkin(10));
t282 = t423 * t371 + t422 * t372;
t498 = t422 * t428 + t423 * t425;
t220 = -pkin(9) * t498 + t282;
t424 = sin(qJ(6));
t427 = cos(qJ(6));
t360 = -t422 * t425 + t423 * t428;
t707 = -t422 * t371 + t423 * t372;
t725 = -pkin(9) * t360 + t707;
t111 = t427 * t220 + t424 * t725;
t426 = sin(qJ(2));
t675 = -t426 / 0.2e1;
t766 = t111 * t675;
t186 = -t424 * t360 - t427 * t498;
t610 = t751 * t186;
t331 = t498 * t426;
t319 = t427 * t331;
t619 = t426 * t428;
t621 = t425 * t426;
t329 = -t422 * t621 + t423 * t619;
t628 = t424 * t329;
t519 = t628 / 0.2e1 + t319 / 0.2e1;
t746 = t186 * t426;
t749 = -t746 / 0.2e1 + t519;
t758 = t749 * qJD(1);
t757 = -t758 + t610;
t503 = -t424 * t220 + t427 * t725;
t429 = cos(qJ(2));
t614 = t428 * t429;
t620 = t425 * t429;
t330 = t422 * t620 - t423 * t614;
t318 = t427 * t330;
t332 = t498 * t429;
t625 = t424 * t332;
t526 = t318 + t625;
t674 = t426 / 0.2e1;
t402 = pkin(4) * t425 + qJ(3);
t322 = pkin(5) * t498 + t402;
t684 = -t322 / 0.2e1;
t765 = -t503 * t674 - t526 * t684;
t764 = t751 * t111;
t763 = t751 * t503;
t500 = t427 * t360 - t424 * t498;
t483 = t751 * t500;
t317 = t427 * t329;
t626 = t424 * t331;
t705 = t317 / 0.2e1 - t626 / 0.2e1;
t715 = t500 * t426;
t741 = t715 / 0.2e1 + t705;
t752 = t741 * qJD(1);
t756 = -t752 - t483;
t748 = t746 / 0.2e1 + t519;
t762 = qJD(3) * t748;
t761 = qJD(3) * t749;
t760 = qJD(6) * t748;
t759 = qJD(6) * t749;
t651 = qJ(3) * t426;
t510 = -t429 * t430 + t651;
t355 = -pkin(1) - t510;
t695 = pkin(3) + pkin(7);
t378 = t695 * t426;
t368 = t428 * t378;
t278 = t355 * t425 - t368;
t250 = qJ(5) * t620 - t278;
t233 = pkin(4) * t426 + t250;
t622 = t425 * t378;
t279 = t355 * t428 + t622;
t251 = -qJ(5) * t614 + t279;
t634 = t423 * t251;
t118 = t422 * t233 + t634;
t121 = -t250 * t422 - t634;
t712 = t118 + t121;
t740 = -t715 / 0.2e1 + t705;
t755 = qJD(3) * t740;
t754 = qJD(6) * t740;
t753 = qJD(6) * t741;
t234 = t422 * t251;
t122 = t423 * t250 - t234;
t117 = t423 * t233 - t234;
t694 = -t117 / 0.2e1;
t750 = -t694 - t122 / 0.2e1;
t575 = qJD(5) * t526;
t747 = qJD(3) * t741 - t575;
t320 = t427 * t332;
t627 = t424 * t330;
t708 = -t320 + t627;
t719 = t526 ^ 2 - t708 ^ 2;
t738 = qJD(1) * t719;
t720 = t186 ^ 2 - t500 ^ 2;
t736 = qJD(2) * t720;
t735 = -t707 / 0.2e1;
t734 = t707 / 0.2e1;
t416 = t429 * pkin(7);
t417 = t429 * pkin(3);
t379 = t416 + t417;
t398 = pkin(4) * t614;
t348 = t398 + t379;
t255 = -pkin(5) * t330 + t348;
t733 = t255 * t526;
t732 = t255 * t708;
t591 = qJD(1) * t426;
t731 = t526 * t591;
t594 = qJD(1) * t708;
t730 = t526 * t594;
t205 = t708 * t591;
t587 = qJD(2) * t186;
t726 = qJD(1) * t526 + t587;
t588 = qJD(2) * t500;
t486 = t588 + t594;
t679 = t498 / 0.2e1;
t680 = t360 / 0.2e1;
t677 = -t423 / 0.2e1;
t678 = t422 / 0.2e1;
t716 = (t329 * t678 + t331 * t677) * pkin(4);
t20 = t122 * t679 + t498 * t694 + t680 * t712 - t716;
t19 = -t716 - (t118 / 0.2e1 + t121 / 0.2e1) * t360 + t498 * t750;
t502 = -t329 * t498 + t331 * t360;
t723 = qJD(3) * t502;
t721 = qJD(5) * t186;
t717 = -t500 / 0.2e1;
t664 = pkin(4) * qJD(4);
t714 = (t360 * t422 - t423 * t498) * t664;
t709 = -t317 + t626;
t702 = qJD(4) * t526;
t701 = qJD(5) * t500;
t418 = t425 ^ 2;
t420 = t428 ^ 2;
t387 = t418 - t420;
t562 = t429 * qJD(1);
t548 = t428 * t562;
t452 = qJD(2) * t387 + 0.2e1 * t425 * t548;
t561 = t429 * qJD(3);
t699 = qJD(2) * t510 - t561;
t697 = t360 ^ 2;
t696 = t498 ^ 2;
t693 = t708 / 0.2e1;
t692 = t186 / 0.2e1;
t689 = t500 / 0.2e1;
t688 = -t186 / 0.2e1;
t534 = -t320 / 0.2e1;
t667 = t428 * pkin(4);
t323 = pkin(5) * t360 + t667;
t683 = -t323 / 0.2e1;
t682 = -t330 / 0.2e1;
t676 = t423 / 0.2e1;
t673 = -t428 / 0.2e1;
t672 = -t429 / 0.2e1;
t671 = pkin(4) * t422;
t670 = pkin(9) * t330;
t668 = t332 * pkin(9);
t32 = 0.2e1 * t708 * t717 + (-t688 + t692) * t526;
t37 = t186 * t526 - t500 * t708;
t666 = t32 * qJD(4) + t37 * qJD(6);
t56 = t186 * t693 + t526 * t689;
t57 = -t526 * t717 + t692 * t708;
t665 = t56 * qJD(4) + t57 * qJD(6);
t613 = t429 * qJ(3);
t376 = t426 * pkin(2) - t613;
t358 = pkin(8) * t426 + t376;
t370 = t379 * t428;
t238 = pkin(4) * t429 + t370 + (-qJ(5) * t426 - t358) * t425;
t344 = t428 * t358;
t369 = t379 * t425;
t609 = t344 + t369;
t253 = qJ(5) * t619 + t609;
t119 = t423 * t238 - t253 * t422;
t80 = pkin(5) * t429 - pkin(9) * t331 + t119;
t663 = t424 * t80;
t87 = t118 + t670;
t662 = t424 * t87;
t120 = t422 * t238 + t423 * t253;
t88 = pkin(9) * t329 + t120;
t661 = t424 * t88;
t89 = t121 - t670;
t660 = t424 * t89;
t90 = t122 + t668;
t659 = t424 * t90;
t658 = t427 * t80;
t657 = t427 * t87;
t656 = t427 * t88;
t655 = t427 * t89;
t654 = t427 * t90;
t347 = (-t667 - t695) * t426;
t254 = -pkin(5) * t329 + t347;
t446 = t426 * pkin(5) + t117 + t668;
t71 = t427 * t446;
t38 = -t71 + t662;
t5 = (t658 - t661) * t426 - t38 * t429 - t254 * t526 + t255 * t709;
t653 = t5 * qJD(1);
t440 = t424 * t446;
t39 = t440 + t657;
t525 = t319 + t628;
t6 = (t656 + t663) * t426 + t39 * t429 - t254 * t708 - t255 * t525;
t652 = t6 * qJD(1);
t40 = t655 - t659;
t559 = pkin(4) * t620;
t481 = -pkin(5) * t332 - t559;
t24 = t40 * t426 - t481 * t526 + t732;
t650 = qJD(1) * t24;
t41 = t654 + t660;
t25 = t41 * t426 - t481 * t708 - t733;
t649 = qJD(1) * t25;
t33 = -t38 * t426 - t733;
t648 = qJD(1) * t33;
t34 = -t39 * t426 + t732;
t647 = qJD(1) * t34;
t646 = t119 * t360;
t645 = t120 * t498;
t43 = t117 * t331 - t118 * t329;
t21 = t119 * t332 + t120 * t330 - t43;
t644 = t21 * qJD(1);
t22 = t712 * t332 + (-t117 + t122) * t330;
t643 = t22 * qJD(1);
t23 = t117 * t119 + t118 * t120 + t347 * t348;
t642 = t23 * qJD(1);
t26 = t117 * t121 + t118 * t122 - t348 * t559;
t641 = t26 * qJD(1);
t638 = t282 * t329;
t637 = t331 * t707;
t636 = t360 * t332;
t635 = t498 * t330;
t623 = t425 * t358;
t612 = t43 * qJD(1);
t52 = t525 * t526 - t708 * t709;
t611 = t52 * qJD(1);
t419 = t426 ^ 2;
t421 = t429 ^ 2;
t388 = t421 - t419;
t112 = -t426 * t709 + t429 * t526;
t608 = qJD(1) * t112;
t113 = t426 * t525 + t429 * t708;
t607 = qJD(1) * t113;
t540 = t500 * t674;
t123 = t540 - t705;
t606 = qJD(1) * t123;
t124 = t540 + t705;
t605 = qJD(1) * t124;
t537 = t186 * t674;
t140 = t537 + t519;
t600 = qJD(1) * t140;
t141 = t537 - t519;
t599 = qJD(1) * t141;
t157 = -t329 * t330 + t331 * t332;
t598 = qJD(1) * t157;
t195 = -t278 * t426 + t379 * t614;
t597 = qJD(1) * t195;
t196 = -t279 * t426 - t379 * t620;
t596 = qJD(1) * t196;
t211 = t534 + t320 / 0.2e1;
t595 = qJD(1) * t211;
t365 = t388 * t425;
t593 = qJD(1) * t365;
t367 = t388 * t428;
t592 = qJD(1) * t367;
t590 = qJD(1) * t428;
t589 = qJD(2) * qJ(3);
t586 = qJD(2) * t322;
t585 = qJD(2) * t426;
t413 = qJD(2) * t429;
t584 = qJD(3) * t425;
t583 = qJD(3) * t426;
t582 = qJD(3) * t428;
t581 = qJD(4) * t708;
t580 = qJD(4) * t425;
t579 = qJD(4) * t426;
t578 = qJD(4) * t428;
t577 = qJD(4) * t430;
t576 = qJD(5) * t708;
t573 = qJD(6) * t322;
t108 = (-t278 - t368) * t429 - t623 * t426;
t572 = t108 * qJD(1);
t109 = t609 * t426 - t379 * t621 + (t279 - t622) * t429;
t571 = t109 * qJD(1);
t514 = -pkin(2) * t429 - t651;
t373 = -pkin(1) + t514;
t291 = t373 * t429 + t376 * t426;
t568 = t291 * qJD(1);
t292 = -t373 * t426 + t376 * t429;
t567 = t292 * qJD(1);
t566 = t388 * qJD(1);
t565 = t419 * qJD(1);
t564 = t425 * qJD(2);
t563 = t428 * qJD(2);
t560 = t429 * qJD(4);
t558 = pkin(1) * t591;
t557 = pkin(1) * t562;
t556 = pkin(7) * t585;
t555 = t667 / 0.2e1;
t554 = t89 / 0.2e1 + t87 / 0.2e1;
t553 = -t659 / 0.2e1;
t552 = -t654 / 0.2e1;
t551 = pkin(4) * t423 + pkin(5);
t550 = t525 * t591;
t549 = t709 * t591;
t547 = t425 * t563;
t546 = t429 * t564;
t545 = t426 * t560;
t544 = t373 * t376 * qJD(1);
t543 = t373 * t591;
t393 = t426 * t413;
t392 = t426 * t562;
t395 = t429 * t563;
t542 = t425 * t578;
t535 = t620 / 0.2e1;
t528 = t735 + t734;
t524 = -t559 / 0.2e1;
t523 = qJD(4) + t591;
t522 = -qJD(6) - t591;
t521 = t425 * t395;
t518 = -t697 / 0.2e1 - t696 / 0.2e1;
t342 = t424 * t671 - t427 * t551;
t517 = t71 / 0.2e1 + t342 * t675;
t515 = -t318 / 0.2e1 - t625 / 0.2e1;
t3 = t330 * t528 - t19;
t513 = t3 * qJD(1);
t432 = t255 * t717 + t481 * t692 - t526 * t683 + t684 * t708 - t766;
t479 = -t661 / 0.2e1 + t658 / 0.2e1;
t445 = t342 * t672 + t479;
t7 = t432 + t445;
t81 = -t186 * t323 + t322 * t500;
t512 = -t7 * qJD(1) + t81 * qJD(2);
t431 = t255 * t688 + t481 * t717 + t683 * t708 - t765;
t343 = t424 * t551 + t427 * t671;
t480 = -t663 / 0.2e1 - t656 / 0.2e1;
t444 = t343 * t672 + t480;
t8 = t431 + t444;
t82 = t186 * t322 + t323 * t500;
t511 = -t8 * qJD(1) + t82 * qJD(2);
t50 = t498 * t528;
t509 = t19 * qJD(1) + t50 * qJD(2);
t508 = qJD(2) * t32 + t738;
t507 = qJD(1) * t32 + t736;
t506 = qJD(2) * t37 + t738;
t505 = qJD(1) * t37 + t736;
t116 = -t282 * t498 - t360 * t707;
t435 = t117 * t680 + t118 * t679 + t282 * t682 + t332 * t735;
t455 = (-pkin(7) / 0.2e1 - pkin(3) / 0.2e1 - t667 / 0.2e1) * t426;
t27 = t455 + t435;
t497 = -qJD(1) * t27 + qJD(2) * t116;
t437 = t398 / 0.2e1 + t416 / 0.2e1 + t417 / 0.2e1 - t638 / 0.2e1 + t637 / 0.2e1;
t478 = t646 / 0.2e1 + t645 / 0.2e1;
t31 = t437 - t478;
t496 = qJD(1) * t31 + qJD(2) * t402;
t171 = t635 / 0.2e1 + t636 / 0.2e1;
t42 = t117 * t332 + t118 * t330;
t495 = -qJD(1) * t42 - qJD(3) * t171;
t458 = t500 * t429;
t127 = -t458 / 0.2e1 + t515;
t494 = qJD(1) * t127 - t587;
t456 = t186 * t429;
t473 = t627 / 0.2e1 + t534;
t131 = -t456 / 0.2e1 + t473;
t493 = qJD(1) * t131 + t588;
t166 = 0.2e1 * t534 + t627;
t491 = qJD(1) * t166 + t588;
t178 = -0.1e1 / 0.2e1 + t518;
t490 = qJD(1) * t171 + qJD(2) * t178;
t174 = -t635 - t636;
t240 = t330 ^ 2 + t332 ^ 2;
t489 = qJD(1) * t240 + qJD(2) * t174;
t284 = t696 + t697;
t488 = qJD(1) * t174 + qJD(2) * t284;
t475 = t330 * t678 + t332 * t676;
t179 = (t535 + t475) * pkin(4);
t474 = t360 * t677 - t498 * t678;
t232 = (t673 + t474) * pkin(4);
t487 = qJD(1) * t179 + qJD(2) * t232;
t484 = -t565 - t579;
t482 = qJD(6) * t526 + t702;
t477 = t119 * t676 + t120 * t678;
t471 = t430 * t675 - t613 / 0.2e1;
t441 = t255 * t689 + t322 * t693 + t766;
t15 = -t441 + t479;
t470 = qJD(1) * t15 - t500 * t586;
t442 = t255 * t692 + t765;
t16 = -t442 + t480;
t469 = qJD(1) * t16 - t186 * t586;
t468 = -qJD(2) * t57 - t730;
t467 = qJD(2) * t56 + t730;
t466 = -qJD(1) * t57 - t186 * t588;
t465 = qJD(1) * t56 + t500 * t587;
t464 = t523 * t620;
t257 = (t358 / 0.2e1 + t471) * t425;
t463 = -qJ(3) * t563 - t257 * qJD(1);
t448 = t471 * t428;
t258 = t344 / 0.2e1 + t448;
t462 = qJ(3) * t564 - t258 * qJD(1);
t357 = (t420 / 0.2e1 - t418 / 0.2e1) * t429;
t461 = qJD(1) * t357 + t547;
t460 = qJD(2) * t124 - qJD(6) * t211 + t205;
t454 = t421 * t425 * t590 - qJD(2) * t357;
t366 = t387 * t421;
t453 = -qJD(1) * t366 + 0.2e1 * t521;
t436 = t750 * t282 + t712 * t735;
t1 = (t348 * t673 + t402 * t535 + t477) * pkin(4) + t436;
t64 = t402 * t667;
t449 = -t1 * qJD(1) + t64 * qJD(2) - t50 * qJD(3);
t12 = -t424 * t554 + t517 + t552;
t447 = qJD(1) * t12 - qJD(4) * t342;
t443 = qJD(2) * t123 + qJD(6) * t166 + t205 + t581;
t439 = qJD(2) * t514 + t561;
t433 = t343 * t674 + t440 / 0.2e1;
t11 = t427 * t554 + t433 + t553;
t438 = qJD(1) * t11 + qJD(4) * t343;
t412 = pkin(7) * t413;
t406 = t413 / 0.2e1;
t405 = -t562 / 0.2e1;
t404 = t562 / 0.2e1;
t394 = t426 * t590;
t391 = t425 * t591;
t364 = -t394 - t578;
t363 = -t391 - t580;
t359 = t392 + t560 / 0.2e1;
t346 = t357 * qJD(4);
t328 = t392 + (qJD(4) / 0.2e1 + qJD(6) / 0.2e1) * t429;
t325 = t343 * qJD(6);
t324 = t342 * qJD(6);
t231 = pkin(4) * t474 + t555;
t207 = -t369 - t344 / 0.2e1 + t448;
t206 = t370 - t623 / 0.2e1 + t471 * t425;
t180 = pkin(4) * t475 + t524;
t177 = 0.1e1 / 0.2e1 + t518;
t164 = t174 * qJD(5);
t163 = t171 * qJD(5);
t132 = t456 / 0.2e1 + t473;
t126 = t458 / 0.2e1 + t515;
t49 = t50 * qJD(4);
t48 = -qJD(2) * t141 - t731;
t35 = qJD(2) * t140 + (qJD(6) + t523) * t526;
t30 = t437 + t478;
t28 = t455 - t435;
t18 = t441 + t479;
t17 = t442 + t480;
t14 = -t657 / 0.2e1 + t553 + t655 / 0.2e1 - t433;
t13 = t662 / 0.2e1 + t552 - t660 / 0.2e1 - t517;
t10 = -t431 + t444;
t9 = -t432 + t445;
t4 = t330 * t734 + t707 * t682 - t20;
t2 = pkin(4) * t477 + t348 * t555 + t402 * t524 - t436;
t29 = [0, 0, 0, t393, t388 * qJD(2), 0, 0, 0, -pkin(1) * t585, -pkin(1) * t413, 0, qJD(2) * t292 - t426 * t561, -qJD(2) * t291 + qJD(3) * t419 (qJD(2) * t376 - t583) * t373, -t393 * t418 + t421 * t542, -qJD(4) * t366 - 0.2e1 * t426 * t521, -qJD(2) * t365 - t428 * t545, -qJD(2) * t367 + t425 * t545, t393, qJD(2) * t108 + qJD(4) * t196 + t419 * t584, -qJD(2) * t109 - qJD(4) * t195 + t419 * t582, qJD(2) * t21 + qJD(3) * t157 + qJD(4) * t22 + qJD(5) * t240, qJD(2) * t23 + qJD(3) * t43 + qJD(4) * t26 + qJD(5) * t42 (qJD(2) * t525 + t482) * t708, qJD(2) * t52 + t719 * t751, qJD(2) * t113 + t426 * t482, qJD(2) * t112 + (-qJD(6) * t708 - t581) * t426, t393, qJD(2) * t5 + qJD(4) * t24 + qJD(6) * t34 + (qJD(3) * t525 - t576) * t426, -qJD(2) * t6 - qJD(4) * t25 - qJD(6) * t33 + (-qJD(3) * t709 - t575) * t426; 0, 0, 0, t392, t566, t413, -t585, 0, -t412 - t558, t556 - t557, t439, t412 + t567, -t556 - t568, pkin(7) * t439 + t544, -t346 + (-t418 * t562 + t547) * t426, -t426 * t452 + 0.2e1 * t429 * t542, t395 - t593, -t546 - t592, t359, t206 * qJD(4) - t378 * t564 - t428 * t699 + t572, t207 * qJD(4) - t378 * t563 + t425 * t699 - t571, t644 + (-t637 + t638 - t645 - t646) * qJD(2) - t723 + t4 * qJD(4) + t164, t642 + (t119 * t707 + t120 * t282 + t347 * t402) * qJD(2) + t30 * qJD(3) + t2 * qJD(4) + t28 * qJD(5), t486 * t525 + t665, t611 + (t186 * t525 - t500 * t709) * qJD(2) + t666, qJD(4) * t140 + t413 * t500 + t607 + t760, -qJD(4) * t123 + t186 * t413 + t608 + t754, t328, t653 + (-t186 * t254 + t322 * t709 + t429 * t503) * qJD(2) + t126 * qJD(3) + t9 * qJD(4) - t124 * qJD(5) + t18 * qJD(6), -t652 + (-t111 * t429 + t254 * t500 + t322 * t525) * qJD(2) + t132 * qJD(3) + t10 * qJD(4) - t141 * qJD(5) + t17 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t413, -t392, t565, t412 - t543, 0, 0, 0, 0, 0, t425 * t565 + t395, t428 * t565 - t546, -qJD(2) * t502 + t598, t30 * qJD(2) + t20 * qJD(4) + t163 + t612 + t723, 0, 0, 0, 0, 0, qJD(2) * t126 + qJD(4) * t748 + t550 + t760, qJD(2) * t132 + qJD(4) * t740 - t549 + t754; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t454, t453, -t523 * t614, t464, t406, qJD(2) * t206 - qJD(4) * t279 + t596, qJD(2) * t207 + qJD(4) * t278 - t597, t643 + t4 * qJD(2) + (-t330 * t423 + t332 * t422) * t664, t641 + t2 * qJD(2) + t20 * qJD(3) + t180 * qJD(5) + (t121 * t423 + t122 * t422) * t664, t467, t508, t35, -t443, t406, qJD(2) * t9 + qJD(4) * t40 + qJD(6) * t14 + t650 + t762, qJD(2) * t10 - qJD(4) * t41 + qJD(6) * t13 - t649 + t755; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t489, qJD(2) * t28 + qJD(4) * t180 - t495, 0, 0, 0, 0, 0, -t460, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t468, t506, qJD(2) * t748 - t522 * t526 + t702, qJD(2) * t740 - qJD(4) * t166 + t522 * t708, t406, qJD(2) * t18 + qJD(4) * t14 + qJD(5) * t211 - qJD(6) * t39 + t647 + t762, qJD(2) * t17 + qJD(4) * t13 + qJD(6) * t38 - t648 + t755; 0, 0, 0, -t392, -t566, 0, 0, 0, t558, t557, 0, -t567, t568, -t544, t392 * t418 - t346, 0.2e1 * t428 * t464, -t425 * t579 + t593, -t426 * t578 + t592, -t359, qJD(4) * t257 - t572, qJD(4) * t258 + t571, -qJD(4) * t3 + t164 - t644, qJD(3) * t31 - qJD(4) * t1 - qJD(5) * t27 - t642, -t525 * t594 + t665, -t611 + t666, qJD(4) * t141 - t607 - t759, -qJD(4) * t124 - t608 - t753, -t328, qJD(3) * t127 - qJD(4) * t7 - qJD(5) * t123 - qJD(6) * t15 - t653, qJD(3) * t131 - qJD(4) * t8 - qJD(5) * t140 - qJD(6) * t16 + t652; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), qJ(3) * qJD(3), -t542, t387 * qJD(4), 0, 0, 0, qJ(3) * t578 + t584, -qJ(3) * t580 + t582, qJD(5) * t284, qJD(3) * t402 + qJD(4) * t64 + qJD(5) * t116, t610 * t500, t751 * t720, 0, 0, 0, -qJD(3) * t186 + qJD(4) * t81 + t500 * t573, qJD(3) * t500 + qJD(4) * t82 + t186 * t573; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), t589, 0, 0, 0, 0, 0, t564, t563, 0, qJD(5) * t177 - t49 + t496, 0, 0, 0, 0, 0, t494, t493; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t461, t452, t363, t364, t405, -t425 * t577 - t463, -t428 * t577 - t462, -t513 - t714, t231 * qJD(5) + (-t282 * t423 + t422 * t707) * t664 + t449, t465, t507, t599 + t610, -t483 - t605, t405, t512 - t764, t511 - t763; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t488, qJD(3) * t177 + qJD(4) * t231 + t497, 0, 0, 0, 0, 0, -t606, -t600; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t466, t505, t757, t756, t405, -t470 - t764, -t469 - t763; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t392, -t565, t543, 0, 0, 0, 0, 0, t484 * t425, t484 * t428, -t598, -qJD(2) * t31 - qJD(4) * t19 + t163 - t612, 0, 0, 0, 0, 0, -qJD(2) * t127 - qJD(4) * t749 - t550 - t759, -qJD(2) * t131 - qJD(4) * t741 + t549 - t753; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(2), -t589, 0, 0, 0, 0, 0, -t564, -t563, 0, qJD(5) * t178 - t49 - t496, 0, 0, 0, 0, 0, -t494, -t493; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t363, t364, 0, -t509 + t714, 0, 0, 0, 0, 0, t757, t756; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t490, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t757, t756; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t454, -t453 (t548 + t564) * t426 (-t425 * t562 + t563) * t426, t406, -qJD(2) * t257 + t425 * t583 - t596, -qJD(2) * t258 + t426 * t582 + t597, qJD(2) * t3 - t643, qJD(2) * t1 + qJD(3) * t19 + qJD(5) * t179 - t641, -t467, -t508, t48, t460, t406, qJD(2) * t7 - qJD(6) * t11 - t576 - t650 + t761, qJD(2) * t8 - qJD(6) * t12 + t649 + t747; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t461, -t452, t391, t394, t404, t463, t462, t513, qJD(5) * t232 - t449, -t465, -t507, -t599, t605, t404, -t512 - t701, -t511 - t721; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t391, t394, 0, t509, 0, 0, 0, 0, 0, t758, t752; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t325, t324; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t487, 0, 0, 0, 0, 0, -t486, -t726; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t595, 0, -t325 - t438, t324 - t447; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t489, qJD(2) * t27 - qJD(4) * t179 + t495, 0, 0, 0, 0, 0, t443, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t488, -qJD(3) * t178 - qJD(4) * t232 - t497, 0, 0, 0, 0, 0, t483 + t606, t600 + t610; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t490, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t487, 0, 0, 0, 0, 0, t486, t726; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t491, t726; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t468, -t506, qJD(2) * t749 - t731, qJD(2) * t741 + qJD(4) * t211 + t205, t406, qJD(2) * t15 + qJD(4) * t11 - qJD(5) * t166 - t647 + t761, qJD(2) * t16 + qJD(4) * t12 + t648 + t747; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t466, -t505, t758, t752, t404, t470 - t701, t469 - t721; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t758, t752; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t595, 0, t438, t447; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t491, -t726; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t29;
