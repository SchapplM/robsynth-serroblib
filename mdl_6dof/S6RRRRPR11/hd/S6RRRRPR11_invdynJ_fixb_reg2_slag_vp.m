% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRPR11_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR11_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPR11_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPR11_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_invdynJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:23:43
% EndTime: 2019-03-09 23:24:49
% DurationCPUTime: 39.42s
% Computational Cost: add. (45134->1076), mult. (109087->1436), div. (0->0), fcn. (88496->18), ass. (0->432)
t442 = sin(qJ(2));
t631 = cos(pkin(6));
t547 = pkin(1) * t631;
t414 = t442 * t547;
t441 = sin(qJ(3));
t444 = cos(qJ(3));
t500 = pkin(3) * t441 - pkin(10) * t444;
t436 = sin(pkin(6));
t445 = cos(qJ(2));
t599 = t436 * t445;
t701 = -(t414 + (pkin(8) + t500) * t599) * qJD(1) + t500 * qJD(3);
t439 = sin(qJ(6));
t526 = t631 * qJD(1);
t486 = t526 + qJD(2);
t573 = qJD(1) * t436;
t543 = t442 * t573;
t313 = t441 * t486 + t444 * t543;
t440 = sin(qJ(4));
t443 = cos(qJ(4));
t572 = qJD(1) * t445;
t406 = t436 * t572;
t493 = t406 - qJD(3);
t258 = t313 * t440 + t443 * t493;
t260 = t443 * t313 - t440 * t493;
t435 = sin(pkin(12));
t437 = cos(pkin(12));
t524 = t258 * t437 + t435 * t260;
t658 = cos(qJ(6));
t478 = t658 * t524;
t490 = -t258 * t435 + t437 * t260;
t561 = qJD(1) * qJD(2);
t532 = t442 * t561;
t509 = t436 * t532;
t560 = qJDD(1) * t445;
t405 = t436 * t560;
t557 = qJDD(3) - t405;
t467 = t509 + t557;
t566 = qJD(4) * t440;
t468 = qJD(3) * t486;
t531 = t445 * t561;
t558 = t442 * qJDD(1);
t473 = t531 + t558;
t454 = t436 * t473 + t468;
t521 = t631 * qJDD(1);
t480 = t521 + qJDD(2);
t514 = t441 * t543;
t460 = qJD(3) * t514 - t441 * t480;
t665 = t454 * t444 - t460;
t680 = qJD(4) * t493 - t665;
t136 = t313 * t566 - t440 * t467 + t443 * t680;
t309 = t443 * t467;
t567 = qJD(4) * t260;
t137 = t440 * t665 - t309 + t567;
t525 = -t435 * t136 + t137 * t437;
t562 = qJD(6) * t439;
t71 = -t136 * t437 - t137 * t435;
t25 = qJD(6) * t478 + t439 * t525 + t490 * t562 - t658 * t71;
t311 = -t444 * t486 + t514;
t303 = qJD(4) + t311;
t294 = qJD(6) + t303;
t92 = t439 * t490 + t478;
t640 = t294 * t92;
t700 = -t25 + t640;
t647 = t92 ^ 2;
t685 = -t439 * t524 + t490 * t658;
t648 = t685 ^ 2;
t699 = -t647 + t648;
t646 = t92 * t685;
t512 = t444 * t406;
t298 = t440 * t512 - t443 * t543;
t568 = qJD(3) * t444;
t496 = -t440 * t568 + t298;
t593 = t444 * t445;
t597 = t440 * t442;
t299 = (t443 * t593 + t597) * t573;
t698 = -t443 * t568 + t299;
t511 = pkin(1) * t526;
t337 = -pkin(8) * t543 + t445 * t511;
t502 = pkin(2) * t442 - pkin(9) * t445;
t338 = t502 * t573;
t248 = t444 * t337 + t441 * t338;
t223 = pkin(10) * t543 + t248;
t569 = qJD(3) * t441;
t553 = pkin(9) * t569;
t697 = t701 * t443 + (t223 + t553) * t440;
t656 = pkin(3) * t444;
t501 = pkin(10) * t441 + t656;
t386 = -pkin(2) - t501;
t564 = qJD(4) * t443;
t696 = -t443 * t223 + t386 * t564 + t440 * t701;
t291 = -pkin(2) * t486 - t337;
t182 = t311 * pkin(3) - t313 * pkin(10) + t291;
t576 = pkin(8) * t599 + t414;
t329 = t631 * pkin(9) + t576;
t292 = qJD(2) * pkin(9) + qJD(1) * t329;
t484 = -pkin(2) * t445 - pkin(9) * t442 - pkin(1);
t301 = t484 * t573;
t209 = t444 * t292 + t441 * t301;
t186 = -pkin(10) * t493 + t209;
t117 = t443 * t182 - t186 * t440;
t84 = -qJ(5) * t260 + t117;
t76 = pkin(4) * t303 + t84;
t118 = t182 * t440 + t186 * t443;
t85 = -qJ(5) * t258 + t118;
t80 = t435 * t85;
t47 = t437 * t76 - t80;
t671 = pkin(11) * t490;
t34 = pkin(5) * t303 + t47 - t671;
t639 = t437 * t85;
t48 = t435 * t76 + t639;
t687 = pkin(11) * t524;
t38 = t48 - t687;
t534 = qJD(6) * t658;
t570 = qJD(2) * t445;
t539 = t441 * t570;
t210 = t436 * (qJD(1) * (t442 * t568 + t539) + t441 * t558) + t441 * t468 - t444 * t480;
t201 = qJDD(4) + t210;
t482 = qJD(2) * t511;
t508 = pkin(1) * t521;
t530 = t436 * t558;
t516 = t442 * t482 - t445 * t508 + (t436 * t531 + t530) * pkin(8);
t233 = -pkin(2) * t480 + t516;
t101 = t210 * pkin(3) - pkin(10) * t665 + t233;
t549 = pkin(8) * t405 + t442 * t508 + t445 * t482;
t261 = -pkin(8) * t509 + t549;
t232 = pkin(9) * t480 + t261;
t479 = qJD(2) * t502;
t244 = (qJD(1) * t479 + qJDD(1) * t484) * t436;
t472 = -t444 * t232 - t441 * t244 + t292 * t569 - t301 * t568;
t89 = pkin(10) * t467 - t472;
t37 = -qJD(4) * t118 + t443 * t101 - t440 * t89;
t24 = pkin(4) * t201 + qJ(5) * t136 - qJD(5) * t260 + t37;
t475 = -t440 * t101 - t182 * t564 + t186 * t566 - t443 * t89;
t30 = -qJ(5) * t137 - qJD(5) * t258 - t475;
t8 = t437 * t24 - t435 * t30;
t6 = pkin(5) * t201 - pkin(11) * t71 + t8;
t9 = t435 * t24 + t437 * t30;
t7 = -pkin(11) * t525 + t9;
t1 = t34 * t534 - t38 * t562 + t439 * t6 + t658 * t7;
t657 = sin(qJ(1));
t504 = t631 * t657;
t659 = cos(qJ(1));
t359 = -t442 * t504 + t445 * t659;
t545 = t436 * t657;
t284 = t359 * t444 + t441 * t545;
t358 = t442 * t659 + t445 * t504;
t432 = qJ(4) + pkin(12);
t427 = qJ(6) + t432;
t420 = sin(t427);
t421 = cos(t427);
t205 = t284 * t421 + t358 * t420;
t600 = t436 * t442;
t355 = t441 * t631 + t444 * t600;
t505 = t631 * t659;
t357 = t442 * t505 + t445 * t657;
t546 = t436 * t659;
t280 = t357 * t444 - t441 * t546;
t356 = t442 * t657 - t445 * t505;
t677 = t280 * t421 + t356 * t420;
t208 = -t441 * t292 + t444 * t301;
t185 = pkin(3) * t493 - t208;
t147 = t258 * pkin(4) + qJD(5) + t185;
t79 = pkin(5) * t524 + t147;
t695 = t79 * t92 + g(1) * t205 + g(2) * t677 - g(3) * (-t355 * t421 + t420 * t599) - t1;
t236 = pkin(3) * t313 + pkin(10) * t311;
t145 = -t208 * t440 + t443 * t236;
t438 = -qJ(5) - pkin(10);
t527 = qJD(4) * t438;
t694 = -pkin(4) * t313 - qJD(5) * t440 - t145 + (-qJ(5) * t311 + t527) * t443;
t146 = t443 * t208 + t440 * t236;
t563 = qJD(5) * t443;
t615 = t311 * t440;
t693 = qJ(5) * t615 - t440 * t527 + t146 - t563;
t595 = t443 * t444;
t416 = pkin(9) * t595;
t513 = t441 * t406;
t692 = -pkin(4) * t513 + qJ(5) * t299 - t441 * t563 + (pkin(4) * t441 - qJ(5) * t595) * qJD(3) + (-t416 + (qJ(5) * t441 - t386) * t440) * qJD(4) + t697;
t596 = t441 * t443;
t691 = -qJ(5) * t298 - (-pkin(9) * qJD(3) - qJ(5) * qJD(4)) * t596 - (-qJD(5) * t441 + (-pkin(9) * qJD(4) - qJ(5) * qJD(3)) * t444) * t440 - t696;
t26 = qJD(6) * t685 + t439 * t71 + t658 * t525;
t638 = t685 * t294;
t690 = -t26 + t638;
t366 = t435 * t440 - t437 * t443;
t367 = t435 * t443 + t437 * t440;
t565 = qJD(4) * t441;
t585 = t437 * t298 + t299 * t435 + t366 * t565 - t367 * t568;
t352 = t367 * qJD(4);
t584 = t352 * t441 - t496 * t435 + t698 * t437;
t581 = t367 * t311 + t352;
t580 = t303 * t366;
t689 = t441 * t564 - t496;
t662 = t513 - t569;
t12 = t439 * t34 + t38 * t658;
t2 = -qJD(6) * t12 - t439 * t7 + t658 * t6;
t204 = -t284 * t420 + t358 * t421;
t678 = t280 * t420 - t356 * t421;
t688 = -t685 * t79 - g(1) * t204 + g(2) * t678 - g(3) * (-t355 * t420 - t421 * t599) + t2;
t637 = t435 * t691 + t437 * t692;
t636 = t435 * t692 - t437 * t691;
t635 = t693 * t435 + t437 * t694;
t634 = t435 * t694 - t693 * t437;
t686 = t524 * t490;
t684 = -t662 * pkin(5) + t584 * pkin(11) + t637;
t683 = -t585 * pkin(11) - t636;
t682 = -pkin(5) * t313 + t580 * pkin(11) + t635;
t681 = -t581 * pkin(11) + t634;
t679 = -t117 * t303 - t475;
t362 = -pkin(8) * t600 + t445 * t547;
t341 = qJD(2) * t362;
t425 = sin(t432);
t426 = cos(t432);
t676 = t280 * t425 - t356 * t426;
t675 = t280 * t426 + t356 * t425;
t674 = t280 * t440 - t356 * t443;
t673 = t280 * t443 + t356 * t440;
t431 = t436 ^ 2;
t672 = 0.2e1 * t431;
t670 = -t118 * t303 - t37;
t669 = t313 * t493;
t476 = t493 * t441;
t506 = -t209 + (t566 + t615) * pkin(4);
t328 = -pkin(2) * t631 - t362;
t354 = t441 * t600 - t444 * t631;
t219 = t354 * pkin(3) - t355 * pkin(10) + t328;
t577 = pkin(2) * t599 + pkin(9) * t600;
t330 = -pkin(1) * t436 - t577;
t238 = t444 * t329 + t441 * t330;
t221 = -pkin(10) * t599 + t238;
t144 = t440 * t219 + t443 * t221;
t666 = t248 + t553;
t247 = -t441 * t337 + t338 * t444;
t222 = -pkin(3) * t543 - t247;
t424 = pkin(9) * t568;
t582 = pkin(4) * t689 - t222 + t424;
t321 = t440 * t386 + t416;
t664 = (qJDD(2) + 0.2e1 * t521) * t436;
t663 = t512 - t568;
t342 = t576 * qJD(2);
t224 = -t284 * t440 + t358 * t443;
t277 = t355 * t440 + t443 * t599;
t661 = -g(1) * t224 + g(2) * t674 + g(3) * t277;
t446 = qJD(1) ^ 2;
t655 = pkin(4) * t435;
t654 = pkin(4) * t440;
t653 = pkin(4) * t443;
t652 = pkin(9) * t444;
t649 = g(3) * t436;
t428 = t441 * pkin(9);
t378 = pkin(5) * t425 + t654;
t645 = pkin(9) + t378;
t369 = t443 * t386;
t270 = -qJ(5) * t596 + t369 + (-pkin(9) * t440 - pkin(4)) * t444;
t598 = t440 * t441;
t287 = -qJ(5) * t598 + t321;
t188 = t437 * t270 - t287 * t435;
t335 = t366 * t441;
t153 = -pkin(5) * t444 + pkin(11) * t335 + t188;
t189 = t435 * t270 + t437 * t287;
t334 = t367 * t441;
t161 = -pkin(11) * t334 + t189;
t86 = t153 * t658 - t439 * t161;
t644 = qJD(6) * t86 + t439 * t684 - t683 * t658;
t87 = t439 * t153 + t161 * t658;
t643 = -qJD(6) * t87 + t683 * t439 + t658 * t684;
t540 = t436 * t570;
t276 = -qJD(3) * t354 + t444 * t540;
t571 = qJD(2) * t442;
t541 = t436 * t571;
t181 = -qJD(4) * t277 + t276 * t443 + t440 * t541;
t275 = qJD(3) * t355 + t436 * t539;
t551 = t440 * t599;
t278 = t355 * t443 - t551;
t339 = t436 * t479;
t154 = -t329 * t569 + t330 * t568 + t441 * t339 + t444 * t341;
t149 = pkin(10) * t541 + t154;
t171 = t275 * pkin(3) - t276 * pkin(10) + t342;
t65 = -qJD(4) * t144 - t149 * t440 + t443 * t171;
t45 = pkin(4) * t275 - qJ(5) * t181 - qJD(5) * t278 + t65;
t180 = -qJD(4) * t551 + t276 * t440 + t355 * t564 - t443 * t541;
t64 = t443 * t149 + t440 * t171 + t219 * t564 - t221 * t566;
t49 = -qJ(5) * t180 - qJD(5) * t277 + t64;
t17 = t435 * t45 + t437 * t49;
t390 = t438 * t440;
t391 = t438 * t443;
t289 = t437 * t390 + t391 * t435;
t255 = -pkin(11) * t367 + t289;
t290 = t435 * t390 - t437 * t391;
t256 = -pkin(11) * t366 + t290;
t156 = t255 * t658 - t439 * t256;
t642 = qJD(6) * t156 + t439 * t682 + t658 * t681;
t157 = t439 * t255 + t256 * t658;
t641 = -qJD(6) * t157 - t439 * t681 + t658 * t682;
t53 = t437 * t84 - t80;
t143 = t443 * t219 - t221 * t440;
t104 = pkin(4) * t354 - qJ(5) * t278 + t143;
t119 = -qJ(5) * t277 + t144;
t60 = t435 * t104 + t437 * t119;
t422 = pkin(4) * t437 + pkin(5);
t343 = t422 * t658 - t439 * t655;
t52 = -t435 * t84 - t639;
t39 = t52 + t687;
t40 = t53 - t671;
t633 = t343 * qJD(6) - t439 * t39 - t40 * t658;
t344 = t439 * t422 + t655 * t658;
t632 = -t344 * qJD(6) - t39 * t658 + t439 * t40;
t628 = t136 * t440;
t627 = t137 * t443;
t626 = t490 * t303;
t625 = t524 * t303;
t624 = t490 ^ 2;
t623 = t201 * t440;
t622 = t201 * t443;
t621 = t258 * t303;
t619 = t258 * t440;
t618 = t260 * t258;
t617 = t260 * t303;
t616 = t303 * t313;
t614 = t313 * t311;
t608 = t356 * t441;
t606 = t358 * t441;
t605 = t420 * t444;
t604 = t421 * t444;
t603 = t425 * t444;
t602 = t426 * t444;
t601 = t431 * t446;
t594 = t444 * t210;
t592 = t334 * t534 - t335 * t562 - t439 * t585 + t584 * t658;
t243 = -t439 * t334 - t335 * t658;
t591 = qJD(6) * t243 - t439 * t584 - t585 * t658;
t268 = -t439 * t366 + t367 * t658;
t590 = -qJD(6) * t268 + t439 * t580 - t581 * t658;
t589 = t366 * t534 + t367 * t562 + t439 * t581 + t580 * t658;
t588 = -pkin(5) * t585 + t582;
t587 = (-t443 * t569 - t444 * t566) * pkin(9) + t696;
t586 = -qJD(4) * t321 + t697;
t583 = pkin(5) * t581 + t506;
t376 = pkin(4) * t598 + t428;
t575 = t659 * pkin(1) + pkin(8) * t545;
t433 = t442 ^ 2;
t434 = t445 ^ 2;
t574 = t433 - t434;
t556 = g(3) * t600;
t555 = g(3) * t599;
t552 = t445 * t601;
t550 = t359 * pkin(2) + t575;
t423 = pkin(3) + t653;
t548 = pkin(9) + t654;
t542 = t441 * t572;
t535 = g(3) * t577;
t533 = pkin(1) * t672;
t379 = pkin(5) * t426 + t653;
t16 = -t435 * t49 + t437 * t45;
t59 = t437 * t104 - t119 * t435;
t237 = -t441 * t329 + t330 * t444;
t523 = t445 * t493;
t522 = t303 * t443;
t520 = qJD(3) * t493;
t518 = t442 * t552;
t517 = t441 * t232 - t444 * t244 + t292 * t568 + t301 * t569;
t510 = t442 * t531;
t507 = -pkin(1) * t657 + pkin(8) * t546;
t503 = t436 * t446 * t631;
t279 = t357 * t441 + t444 * t546;
t283 = t359 * t441 - t444 * t545;
t499 = -g(1) * t279 + g(2) * t283;
t498 = -g(1) * t356 + g(2) * t358;
t497 = g(1) * t359 + g(2) * t357;
t220 = pkin(3) * t599 - t237;
t494 = t525 * pkin(5);
t492 = t524 ^ 2;
t491 = -t117 * t443 - t118 * t440;
t374 = pkin(3) + t379;
t430 = -pkin(11) + t438;
t489 = t374 * t444 - t430 * t441;
t488 = t423 * t444 - t438 * t441;
t485 = 0.2e1 * t526 + qJD(2);
t483 = pkin(9) * t358 + t550;
t481 = -t357 * pkin(2) + t507;
t155 = -t329 * t568 - t330 * t569 + t339 * t444 - t441 * t341;
t192 = -t277 * t435 + t278 * t437;
t50 = pkin(5) * t354 - pkin(11) * t192 + t59;
t191 = t437 * t277 + t278 * t435;
t54 = -pkin(11) * t191 + t60;
t18 = -t439 * t54 + t50 * t658;
t19 = t439 * t50 + t54 * t658;
t125 = -t439 * t191 + t192 * t658;
t474 = -pkin(10) * t201 + t185 * t303;
t471 = g(1) * t659 + g(2) * t657;
t470 = g(1) * t283 + g(2) * t279 + g(3) * t354;
t469 = g(1) * t284 + g(2) * t280 + g(3) * t355;
t165 = pkin(4) * t277 + t220;
t465 = -t356 * pkin(9) + t481;
t90 = -pkin(3) * t467 + t517;
t464 = t470 - t90;
t463 = -g(1) * t358 - g(2) * t356 + t555;
t462 = -t497 - t556;
t459 = t463 * t441;
t457 = t441 * t517 - t444 * t472 - t497;
t150 = -pkin(3) * t541 - t155;
t456 = pkin(10) * qJD(4) * t303 - t464;
t96 = pkin(4) * t180 + t150;
t61 = t137 * pkin(4) + qJDD(5) + t90;
t449 = -t470 + t61;
t347 = t358 * pkin(2);
t345 = t356 * pkin(2);
t340 = t576 * qJD(1);
t322 = pkin(5) * t366 - t423;
t320 = -t440 * t652 + t369;
t271 = pkin(5) * t334 + t376;
t267 = t366 * t658 + t367 * t439;
t242 = t334 * t658 - t335 * t439;
t225 = t284 * t443 + t358 * t440;
t217 = t284 * t426 + t358 * t425;
t216 = -t284 * t425 + t358 * t426;
t195 = qJDD(6) + t201;
t138 = -t201 * t444 - t303 * t476;
t134 = pkin(4) * t260 + pkin(5) * t490;
t131 = t201 * t354 + t275 * t303;
t124 = t191 * t658 + t192 * t439;
t115 = pkin(5) * t191 + t165;
t114 = -t180 * t435 + t181 * t437;
t112 = t437 * t180 + t181 * t435;
t57 = pkin(5) * t112 + t96;
t42 = qJD(6) * t125 + t112 * t658 + t439 * t114;
t41 = t439 * t112 - t114 * t658 + t191 * t534 + t192 * t562;
t33 = t494 + t61;
t13 = -pkin(11) * t112 + t17;
t11 = t34 * t658 - t439 * t38;
t10 = pkin(5) * t275 - pkin(11) * t114 + t16;
t4 = -qJD(6) * t19 + t10 * t658 - t439 * t13;
t3 = qJD(6) * t18 + t439 * t10 + t13 * t658;
t5 = [0, 0, 0, 0, 0, qJDD(1), g(1) * t657 - g(2) * t659, t471, 0, 0 (qJDD(1) * t433 + 0.2e1 * t510) * t431 (t445 * t558 - t561 * t574) * t672, t442 * t664 + t485 * t540 (qJDD(1) * t434 - 0.2e1 * t510) * t431, t445 * t664 - t485 * t541, t480 * t631, -t342 * t486 + t362 * t480 - t516 * t631 + g(1) * t357 - g(2) * t359 + (-t532 + t560) * t533, -t261 * t631 - t341 * t486 - t473 * t533 - t480 * t576 + t498 ((-t337 * qJD(2) + qJDD(1) * t576 + t261) * t445 + (-qJD(2) * t340 - qJDD(1) * t362 + t516) * t442 - t471) * t436, t431 * qJDD(1) * pkin(1) ^ 2 - g(1) * t507 - g(2) * t575 + t261 * t576 - t337 * t342 + t340 * t341 - t362 * t516, t313 * t276 + t355 * t665, -t355 * t210 - t313 * t275 - t276 * t311 - t354 * t665, -t276 * t493 + t355 * t557 + ((qJD(1) * t355 + t313) * t571 - t665 * t445) * t436, t210 * t354 + t275 * t311, t275 * t493 - t354 * t557 + (t210 * t445 + (-qJD(1) * t354 - t311) * t571) * t436 (-t557 * t445 + (-t406 - t493) * t571) * t436, -t155 * t493 + t237 * t557 + t342 * t311 + t328 * t210 + t233 * t354 + t291 * t275 + g(1) * t280 - g(2) * t284 + (t517 * t445 + (qJD(1) * t237 + t208) * t571) * t436, t154 * t493 - t209 * t541 + t233 * t355 - t238 * t467 + t291 * t276 + t342 * t313 + t328 * t665 - t472 * t599 + t499, -t154 * t311 - t155 * t313 - t208 * t276 - t209 * t275 - t238 * t210 - t237 * t665 + t354 * t472 + t355 * t517 - t498, -g(1) * t465 - g(2) * t483 + t209 * t154 + t208 * t155 + t233 * t328 - t237 * t517 - t238 * t472 + t291 * t342, -t136 * t278 + t181 * t260, t136 * t277 - t137 * t278 - t180 * t260 - t181 * t258, -t136 * t354 + t181 * t303 + t201 * t278 + t260 * t275, t137 * t277 + t180 * t258, -t137 * t354 - t180 * t303 - t201 * t277 - t258 * t275, t131, g(1) * t673 - g(2) * t225 + t117 * t275 + t220 * t137 + t143 * t201 + t150 * t258 + t185 * t180 + t90 * t277 + t65 * t303 + t37 * t354, -g(1) * t674 - g(2) * t224 - t118 * t275 - t220 * t136 - t144 * t201 + t150 * t260 + t185 * t181 + t90 * t278 - t64 * t303 + t475 * t354, -t117 * t181 - t118 * t180 + t136 * t143 - t137 * t144 - t258 * t64 - t260 * t65 + t277 * t475 - t278 * t37 - t499, -t475 * t144 + t118 * t64 + t37 * t143 + t117 * t65 + t90 * t220 + t185 * t150 - g(1) * (-pkin(3) * t280 - pkin(10) * t279 + t465) - g(2) * (pkin(3) * t284 + pkin(10) * t283 + t483) t114 * t490 + t192 * t71, -t112 * t490 - t114 * t524 - t71 * t191 - t192 * t525, t114 * t303 + t192 * t201 + t275 * t490 + t354 * t71, t112 * t524 + t191 * t525, -t112 * t303 - t191 * t201 - t275 * t524 - t354 * t525, t131, g(1) * t675 - g(2) * t217 + t147 * t112 + t16 * t303 + t165 * t525 + t61 * t191 + t59 * t201 + t47 * t275 + t8 * t354 + t96 * t524, -g(1) * t676 - g(2) * t216 + t147 * t114 + t165 * t71 - t17 * t303 + t61 * t192 - t60 * t201 - t48 * t275 - t9 * t354 + t96 * t490, -t48 * t112 - t47 * t114 - t16 * t490 - t17 * t524 - t9 * t191 - t8 * t192 - t525 * t60 - t59 * t71 - t499, t9 * t60 + t48 * t17 + t8 * t59 + t47 * t16 + t61 * t165 + t147 * t96 - g(1) * (t279 * t438 - t280 * t423 - t356 * t548 + t481) - g(2) * (-t283 * t438 + t284 * t423 + t358 * t548 + t550) -t125 * t25 - t41 * t685, t124 * t25 - t125 * t26 + t41 * t92 - t42 * t685, t125 * t195 - t25 * t354 + t275 * t685 - t294 * t41, t124 * t26 + t42 * t92, -t124 * t195 - t26 * t354 - t275 * t92 - t294 * t42, t195 * t354 + t275 * t294, g(1) * t677 - g(2) * t205 + t11 * t275 + t115 * t26 + t33 * t124 + t18 * t195 + t2 * t354 + t4 * t294 + t79 * t42 + t57 * t92, -g(1) * t678 - g(2) * t204 - t1 * t354 - t115 * t25 - t12 * t275 + t33 * t125 - t19 * t195 - t3 * t294 - t79 * t41 + t57 * t685, -t1 * t124 + t11 * t41 - t12 * t42 - t125 * t2 + t18 * t25 - t19 * t26 - t3 * t92 - t4 * t685 - t499, t1 * t19 + t12 * t3 + t2 * t18 + t11 * t4 + t33 * t115 + t79 * t57 - g(1) * (t279 * t430 - t280 * t374 - t356 * t645 + t481) - g(2) * (-t283 * t430 + t284 * t374 + t358 * t645 + t550); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t518, t574 * t601, -t445 * t503 + t530, t518, t442 * t503 + t405, t480, pkin(1) * t442 * t601 + t340 * t486 - t463 - t516, pkin(1) * t552 + t337 * t486 + (pkin(8) * t561 + g(3)) * t600 + t497 - t549, 0, 0, -t460 * t441 + (t441 * t454 - t669) * t444, -t441 * t210 + t311 * t663 + t313 * t662 + t444 * t665, -t444 * t520 + t441 * t557 + (t444 * t523 + (qJD(2) * t441 - t313) * t442) * t573, -t311 * t476 - t594, t441 * t520 + t444 * t557 + (-t441 * t523 + (qJD(2) * t444 + t311) * t442) * t573, t493 * t543, -pkin(2) * t210 + t247 * t493 - t208 * t543 - t340 * t311 + (-pkin(9) * t467 - t291 * t493) * t441 + (pkin(9) * t520 - t233 - t463) * t444, -pkin(2) * t665 - g(1) * t606 - g(2) * t608 + t209 * t543 - t340 * t313 - t467 * t652 - t666 * t493 + (t233 + t555) * t441 - t663 * t291, -pkin(9) * t594 + t428 * t665 + t457 - t556 + (t247 + t424) * t313 + t666 * t311 + t662 * t209 + t663 * t208, -t233 * pkin(2) - t209 * t248 - t208 * t247 - t291 * t340 + g(1) * t347 + g(2) * t345 - t535 + ((-t208 * t444 - t209 * t441) * qJD(3) + t457) * pkin(9), -t136 * t596 + (-t440 * t565 - t698) * t260, t258 * t299 + t260 * t298 + (-t258 * t443 - t260 * t440) * t568 + (t628 - t627 + (-t260 * t443 + t619) * qJD(4)) * t441, t136 * t444 - t698 * t303 + (-t260 * t493 - t303 * t566 + t622) * t441, t137 * t598 + t258 * t689, t137 * t444 + t496 * t303 + (t258 * t493 - t303 * t564 - t623) * t441, t138, -t185 * t298 + t320 * t201 - t222 * t258 + t586 * t303 + t462 * t440 + (-t37 + (pkin(9) * t258 + t185 * t440) * qJD(3) - t463 * t443) * t444 + (pkin(9) * t137 - t117 * t493 + t185 * t564 + t90 * t440) * t441, -t185 * t299 - t321 * t201 - t222 * t260 - t587 * t303 + t462 * t443 + (-t475 + (pkin(9) * t260 + t185 * t443) * qJD(3) + t463 * t440) * t444 + (-pkin(9) * t136 + t118 * t493 - t185 * t566 + t90 * t443) * t441, t117 * t299 + t118 * t298 + t136 * t320 - t137 * t321 - t586 * t260 - t587 * t258 + t491 * t568 + (t475 * t440 - t37 * t443 + (t117 * t440 - t118 * t443) * qJD(4) - t463) * t441, -t475 * t321 + t37 * t320 - t185 * t222 - g(1) * (-pkin(10) * t606 - t358 * t656 - t347) - g(2) * (-pkin(10) * t608 - t356 * t656 - t345) - g(3) * (t501 * t599 + t577) + t587 * t118 + t586 * t117 + (t185 * t568 + t90 * t441 - t497) * pkin(9), -t335 * t71 - t490 * t584, -t71 * t334 + t335 * t525 + t490 * t585 + t524 * t584, -t201 * t335 - t303 * t584 - t444 * t71 - t476 * t490, t334 * t525 - t524 * t585, -t334 * t201 + t585 * t303 + t525 * t444 + t476 * t524, t138, t188 * t201 - t8 * t444 + t47 * t569 + t376 * t525 + t61 * t334 - g(1) * (-t358 * t602 + t359 * t425) - g(2) * (-t356 * t602 + t357 * t425) + t637 * t303 - t585 * t147 + (-t47 * t542 - g(3) * (t425 * t442 + t426 * t593)) * t436 + t582 * t524, -t189 * t201 + t9 * t444 - t48 * t569 + t376 * t71 - t61 * t335 - g(1) * (t358 * t603 + t359 * t426) - g(2) * (t356 * t603 + t357 * t426) - t636 * t303 + t582 * t490 - t584 * t147 + (t48 * t542 - g(3) * (-t425 * t593 + t426 * t442)) * t436, -t188 * t71 - t189 * t525 - t9 * t334 + t8 * t335 + t47 * t584 + t48 * t585 - t490 * t637 - t524 * t636 - t459, t9 * t189 + t8 * t188 + t61 * t376 - g(1) * (-t358 * t488 + t359 * t548 - t347) - g(2) * (-t356 * t488 + t357 * t548 - t345) - t535 + t636 * t48 + t637 * t47 - (pkin(4) * t597 + t445 * t488) * t649 + t582 * t147, -t243 * t25 - t592 * t685, t242 * t25 - t243 * t26 - t591 * t685 + t592 * t92, t195 * t243 + t25 * t444 - t294 * t592 - t476 * t685, t242 * t26 + t591 * t92, -t195 * t242 + t26 * t444 - t294 * t591 + t476 * t92, -t195 * t444 - t294 * t476, t86 * t195 - t2 * t444 + t11 * t569 + t271 * t26 + t33 * t242 - g(1) * (-t358 * t604 + t359 * t420) - g(2) * (-t356 * t604 + t357 * t420) + t588 * t92 + t591 * t79 + t643 * t294 + (-t11 * t542 - g(3) * (t420 * t442 + t421 * t593)) * t436, -t87 * t195 + t1 * t444 - t12 * t569 - t271 * t25 + t33 * t243 - g(1) * (t358 * t605 + t359 * t421) - g(2) * (t356 * t605 + t357 * t421) + t588 * t685 - t592 * t79 - t644 * t294 + (t12 * t542 - g(3) * (-t420 * t593 + t421 * t442)) * t436, -t1 * t242 + t11 * t592 - t12 * t591 - t2 * t243 + t25 * t86 - t26 * t87 - t643 * t685 - t644 * t92 - t459, t1 * t87 + t2 * t86 + t33 * t271 - g(1) * (-t358 * t489 + t359 * t645 - t347) - g(2) * (-t356 * t489 + t357 * t645 - t345) - t535 + t588 * t79 - (t378 * t442 + t445 * t489) * t649 + t644 * t12 + t643 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t614, -t311 ^ 2 + t313 ^ 2, -t311 * t493 + t665, -t614, -t210 - t669, t467, -t209 * t493 - t291 * t313 + t470 - t517, -t208 * t493 + t291 * t311 + t469 + t472, 0, 0, t260 * t522 - t628 (-t136 - t621) * t443 + (-t137 - t617) * t440, -t260 * t313 + t303 * t522 + t623, t303 * t619 - t627, -t303 ^ 2 * t440 + t258 * t313 + t622, -t616, -pkin(3) * t137 - t117 * t313 - t145 * t303 - t209 * t258 + t440 * t474 - t443 * t456, pkin(3) * t136 + t118 * t313 + t146 * t303 - t209 * t260 + t440 * t456 + t443 * t474, t145 * t260 + t146 * t258 + ((-t137 + t567) * pkin(10) + t679) * t443 + ((qJD(4) * t258 - t136) * pkin(10) + t670) * t440 - t469, -t117 * t145 - t118 * t146 - t185 * t209 + t464 * pkin(3) + (qJD(4) * t491 - t37 * t440 - t443 * t475 - t469) * pkin(10), t367 * t71 - t490 * t580, -t71 * t366 - t367 * t525 - t490 * t581 + t524 * t580, t201 * t367 - t303 * t580 - t313 * t490, t366 * t525 + t524 * t581, -t366 * t201 - t303 * t581 + t313 * t524, -t616, t147 * t581 + t289 * t201 + t303 * t635 - t47 * t313 + t61 * t366 - t423 * t525 + t426 * t470 + t506 * t524, -t147 * t580 - t201 * t290 - t303 * t634 + t313 * t48 + t367 * t61 - t423 * t71 - t425 * t470 + t490 * t506, -t289 * t71 - t290 * t525 - t9 * t366 - t8 * t367 + t47 * t580 - t48 * t581 - t490 * t635 - t524 * t634 - t469, t9 * t290 + t8 * t289 - t61 * t423 - g(1) * (-t283 * t423 - t284 * t438) - g(2) * (-t279 * t423 - t280 * t438) - g(3) * (-t354 * t423 - t355 * t438) + t634 * t48 + t635 * t47 + t506 * t147, -t25 * t268 - t589 * t685, t25 * t267 - t26 * t268 + t589 * t92 + t590 * t685, t195 * t268 - t294 * t589 - t313 * t685, t26 * t267 - t590 * t92, -t195 * t267 + t294 * t590 + t313 * t92, -t294 * t313, -t11 * t313 + t156 * t195 + t26 * t322 + t267 * t33 + t294 * t641 + t421 * t470 + t583 * t92 - t590 * t79, t12 * t313 - t157 * t195 - t25 * t322 + t268 * t33 - t294 * t642 - t420 * t470 + t583 * t685 - t589 * t79, -t1 * t267 + t11 * t589 + t12 * t590 + t156 * t25 - t157 * t26 - t2 * t268 - t641 * t685 - t642 * t92 - t469, t1 * t157 + t2 * t156 + t33 * t322 - g(1) * (-t283 * t374 - t284 * t430) - g(2) * (-t279 * t374 - t280 * t430) - g(3) * (-t354 * t374 - t355 * t430) + t583 * t79 + t642 * t12 + t641 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t618, -t258 ^ 2 + t260 ^ 2, -t136 + t621, -t618, -t313 * t564 + t440 * t680 + t309 + t617, t201, -t185 * t260 + t661 - t670, g(1) * t225 + g(2) * t673 + g(3) * t278 + t185 * t258 - t679, 0, 0, t686, -t492 + t624, t71 + t625, -t686, -t525 + t626, t201, -t52 * t303 - t147 * t490 - g(1) * t216 + g(2) * t676 - g(3) * (-t355 * t425 - t426 * t599) + (t437 * t201 - t260 * t524) * pkin(4) + t8, t53 * t303 + t147 * t524 + g(1) * t217 + g(2) * t675 - g(3) * (-t355 * t426 + t425 * t599) + (-t201 * t435 - t260 * t490) * pkin(4) - t9 (-t435 * t525 - t437 * t71) * pkin(4) + (t48 + t52) * t490 + (t53 - t47) * t524, -t47 * t52 - t48 * t53 + (-t147 * t260 + t9 * t435 + t8 * t437 + t661) * pkin(4), t646, t699, t700, -t646, t690, t195, -t134 * t92 + t343 * t195 + t294 * t632 + t688, -t134 * t685 - t344 * t195 - t294 * t633 + t695, t25 * t343 - t26 * t344 + (t12 - t632) * t685 + (-t11 - t633) * t92, t1 * t344 + t2 * t343 - t79 * t134 - g(1) * (-t284 * t378 + t358 * t379) - g(2) * (-t280 * t378 + t356 * t379) - g(3) * (-t355 * t378 - t379 * t599) + t633 * t12 + t632 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t525 + t626, t71 - t625, -t492 - t624, t47 * t490 + t48 * t524 + t449, 0, 0, 0, 0, 0, 0, t26 + t638, -t25 - t640, -t647 - t648, t11 * t685 + t12 * t92 + t449 + t494; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t646, t699, t700, -t646, t690, t195, t12 * t294 + t688, t11 * t294 + t695, 0, 0;];
tau_reg  = t5;
