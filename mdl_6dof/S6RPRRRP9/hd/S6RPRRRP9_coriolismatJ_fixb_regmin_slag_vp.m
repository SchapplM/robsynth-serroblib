% Calculate minimal parameter regressor of coriolis matrix for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x29]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRRRP9_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP9_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:28:35
% EndTime: 2019-03-09 06:28:57
% DurationCPUTime: 12.83s
% Computational Cost: add. (10458->639), mult. (21050->859), div. (0->0), fcn. (21520->6), ass. (0->486)
t452 = cos(qJ(4));
t733 = pkin(8) + pkin(9);
t398 = t733 * t452;
t451 = cos(qJ(5));
t384 = t451 * t398;
t449 = sin(qJ(4));
t396 = t733 * t449;
t448 = sin(qJ(5));
t636 = t448 * t396;
t304 = t384 - t636;
t435 = t451 * t452;
t635 = t448 * t449;
t375 = -t435 + t635;
t207 = -qJ(6) * t375 + t304;
t629 = t451 * t449;
t634 = t448 * t452;
t379 = t629 + t634;
t453 = cos(qJ(3));
t624 = t453 * t379;
t747 = -t624 / 0.2e1;
t532 = t451 * t396 + t448 * t398;
t742 = -t379 * qJ(6) - t532;
t744 = t742 / 0.2e1;
t746 = t744 - t742 / 0.2e1;
t583 = qJD(4) + qJD(5);
t745 = 0.2e1 * t449;
t450 = sin(qJ(3));
t522 = pkin(3) * t450 - pkin(8) * t453;
t387 = qJ(2) + t522;
t371 = t452 * t387;
t454 = -pkin(1) - pkin(7);
t631 = t449 * t454;
t537 = pkin(4) - t631;
t627 = t452 * t453;
t581 = pkin(9) * t627;
t258 = t450 * t537 + t371 - t581;
t626 = t452 * t454;
t567 = t450 * t626;
t316 = t449 * t387 + t567;
t632 = t449 * t453;
t285 = -pkin(9) * t632 + t316;
t261 = t451 * t285;
t138 = t448 * t258 + t261;
t682 = qJ(6) * t624;
t123 = t138 - t682;
t743 = t746 * t123;
t237 = t451 * t258;
t638 = t448 * t285;
t137 = -t237 + t638;
t400 = t448 * t632;
t348 = t451 * t627 - t400;
t337 = t348 * qJ(6);
t122 = -t137 - t337;
t690 = t453 * pkin(3);
t692 = t450 * pkin(8);
t397 = t690 + t692;
t385 = t452 * t397;
t264 = t452 * t450 * pkin(9) + t453 * t537 + t385;
t249 = t451 * t264;
t382 = t449 * t397;
t623 = t453 * t454;
t402 = t452 * t623;
t633 = t449 * t450;
t288 = pkin(9) * t633 + t382 + t402;
t637 = t448 * t288;
t517 = t249 / 0.2e1 - t637 / 0.2e1;
t703 = pkin(4) * t452;
t438 = -pkin(3) - t703;
t640 = t438 * t348;
t468 = -t640 / 0.2e1 + t517;
t693 = t449 * pkin(4);
t536 = -t454 + t693;
t515 = t536 * t453;
t496 = -t515 / 0.2e1;
t660 = t304 * t450;
t49 = t660 / 0.2e1 + t379 * t496 + t468;
t487 = t634 / 0.2e1 + t629 / 0.2e1;
t713 = t379 / 0.2e1;
t241 = (t713 - t487) * t453;
t593 = t241 * qJD(2);
t605 = qJD(3) * t438;
t741 = qJD(1) * t49 - t379 * t605 + t593;
t238 = t375 * t693 + t379 * t438;
t572 = t693 / 0.2e1;
t707 = -t450 / 0.2e1;
t462 = -t304 * t707 - t572 * t624;
t691 = t451 * pkin(4);
t570 = t691 / 0.2e1;
t578 = -t703 / 0.2e1;
t714 = -t379 / 0.2e1;
t43 = (t375 * t578 + t536 * t714 + t570) * t453 + t462 + t468;
t740 = t43 * qJD(1) - t238 * qJD(3) + t593;
t653 = t348 * t207;
t669 = t123 * t379;
t739 = -t653 / 0.2e1 - t669 / 0.2e1;
t665 = t742 * t624;
t117 = pkin(5) * t450 + t122;
t671 = t117 * t375;
t738 = t742 * t747 + t665 / 0.2e1 + t671 / 0.2e1 + t739;
t344 = t379 * t450;
t663 = t207 * t344;
t399 = t448 * t633;
t347 = t435 * t450 - t399;
t664 = t742 * t347;
t719 = t347 / 0.2e1;
t737 = t742 * t719 - t663 / 0.2e1 - t664 / 0.2e1;
t444 = t449 ^ 2;
t446 = t452 ^ 2;
t414 = t446 - t444;
t534 = t627 * t745;
t471 = qJD(1) * t534 - qJD(3) * t414;
t548 = -t435 / 0.2e1;
t711 = t400 / 0.2e1;
t716 = -t375 / 0.2e1;
t243 = t711 + (t548 + t716) * t453;
t592 = t243 * qJD(2);
t736 = t532 * t583 - t592;
t735 = t348 ^ 2;
t734 = t379 ^ 2;
t732 = -t117 / 0.2e1;
t533 = t249 - t637;
t121 = pkin(5) * t453 + qJ(6) * t347 + t533;
t731 = t121 / 0.2e1;
t439 = t450 * t454;
t315 = t439 * t449 - t371;
t284 = -t315 - t581;
t639 = t448 * t284;
t147 = -t261 - t639;
t126 = t147 + t682;
t730 = t126 / 0.2e1;
t729 = -t207 / 0.2e1;
t725 = t207 / 0.2e1;
t724 = -t237 / 0.2e1;
t723 = -t258 / 0.2e1;
t550 = -t261 / 0.2e1;
t722 = -t344 / 0.2e1;
t721 = t624 / 0.2e1;
t718 = -t348 / 0.2e1;
t717 = t348 / 0.2e1;
t715 = t375 / 0.2e1;
t549 = -t384 / 0.2e1;
t712 = t399 / 0.2e1;
t437 = pkin(5) + t691;
t710 = -t437 / 0.2e1;
t709 = -t449 / 0.2e1;
t706 = t450 / 0.2e1;
t705 = t451 / 0.2e1;
t704 = t453 / 0.2e1;
t702 = pkin(5) * t624;
t701 = pkin(5) * t347;
t700 = pkin(5) * t348;
t699 = pkin(5) * t375;
t698 = pkin(5) * t379;
t697 = t126 * pkin(5);
t696 = t207 * pkin(5);
t695 = t344 * pkin(5);
t694 = t448 * pkin(4);
t670 = t123 * t344;
t672 = t117 * t347;
t687 = -t672 / 0.2e1 - t670 / 0.2e1;
t686 = pkin(4) * qJD(4);
t685 = pkin(4) * qJD(5);
t684 = pkin(5) * qJD(5);
t683 = pkin(5) * qJD(6);
t466 = t122 * t719 + t670 / 0.2e1 + t687;
t625 = t453 * t348;
t546 = -t625 / 0.2e1;
t14 = (t546 + t715) * pkin(5) + t466;
t680 = qJD(1) * t14;
t620 = t117 - t122;
t17 = t620 * t624;
t679 = qJD(1) * t17;
t263 = t515 + t702;
t21 = -t123 * t620 + t263 * t700;
t678 = qJD(1) * t21;
t79 = t147 * t450 + (t348 * t536 + t624 * t703) * t453;
t677 = qJD(1) * t79;
t630 = t451 * t284;
t148 = t630 - t638;
t276 = t624 * t515;
t582 = pkin(4) * t627;
t80 = t148 * t450 - t348 * t582 + t276;
t676 = qJD(1) * t80;
t92 = -t137 * t450 + t276;
t675 = qJD(1) * t92;
t93 = -t138 * t450 + t348 * t515;
t674 = qJD(1) * t93;
t499 = t582 + t700;
t476 = t499 * t453;
t573 = t694 / 0.2e1;
t642 = t437 * t375;
t477 = -t642 / 0.2e1 + t379 * t573;
t127 = t148 - t337;
t542 = t127 / 0.2e1 + t732;
t543 = t123 / 0.2e1 + t730;
t11 = t476 / 0.2e1 - t542 * t347 + t543 * t344 + t477;
t673 = t11 * qJD(1);
t248 = t448 * t264;
t279 = t451 * t288;
t616 = t279 + t248;
t125 = qJ(6) * t344 + t616;
t18 = -t121 * t348 - t125 * t624 + t670 + t672;
t668 = t18 * qJD(1);
t373 = -pkin(4) * t633 + t439;
t262 = t373 - t695;
t19 = t117 * t121 + t123 * t125 + t262 * t263;
t667 = t19 * qJD(1);
t20 = (-t123 - t126) * t348 + (t117 - t127) * t624;
t666 = t20 * qJD(1);
t22 = t117 * t126 + t123 * t127 + t263 * t499;
t662 = t22 * qJD(1);
t661 = t532 * t450;
t659 = t344 * t450;
t658 = t624 * t375;
t657 = t624 * t453;
t656 = t347 * t624;
t655 = t347 * t375;
t654 = t347 * t450;
t652 = t348 * t344;
t651 = t348 * t375;
t36 = t669 - t671;
t650 = t36 * qJD(1);
t37 = t533 * t450 + t373 * t624 + (-t344 * t536 - t137) * t453;
t649 = t37 * qJD(1);
t648 = t379 * t344;
t647 = t379 * t624;
t646 = t379 * t348;
t38 = t138 * t453 + t347 * t515 - t373 * t348 + t450 * t616;
t645 = t38 * qJD(1);
t644 = t437 * t624;
t643 = t437 * t347;
t641 = t438 * t624;
t447 = t453 ^ 2;
t436 = t447 * t452;
t445 = t450 ^ 2;
t628 = t452 * t445;
t571 = pkin(4) * t707;
t525 = t571 + t284 / 0.2e1;
t60 = t451 * t525 + t724;
t622 = t60 * qJD(1);
t115 = -t646 + t658;
t621 = t583 * t115;
t152 = -t647 / 0.2e1 - t651 / 0.2e1;
t619 = t583 * t152;
t413 = t445 - t447;
t615 = qJD(1) * qJ(2);
t134 = -t656 / 0.2e1 + t652 / 0.2e1;
t614 = qJD(1) * t134;
t154 = -t647 + t651;
t613 = qJD(1) * t154;
t250 = -t315 * t450 - t447 * t631;
t612 = qJD(1) * t250;
t251 = -t316 * t450 - t447 * t626;
t611 = qJD(1) * t251;
t610 = qJD(1) * t348;
t378 = t413 * t449;
t609 = qJD(1) * t378;
t381 = -t436 + t628;
t608 = qJD(1) * t381;
t607 = qJD(2) * t450;
t606 = qJD(3) * t379;
t604 = qJD(3) * t452;
t603 = qJD(4) * t449;
t602 = qJD(4) * t452;
t601 = qJD(5) * t438;
t136 = t652 + t656;
t600 = t136 * qJD(1);
t566 = t449 * t623;
t155 = -t315 * t453 + (t385 + t566) * t450;
t599 = t155 * qJD(1);
t156 = t316 * t453 + (-t402 + t382) * t450;
t598 = t156 * qJD(1);
t488 = -t654 / 0.2e1 + t546;
t516 = t435 / 0.2e1 - t635 / 0.2e1;
t179 = -t488 + t516;
t597 = t179 * qJD(1);
t490 = t659 / 0.2e1 + t657 / 0.2e1;
t180 = -t487 - t490;
t596 = t180 * qJD(1);
t199 = -t657 + t659;
t595 = t199 * qJD(1);
t200 = t625 - t654;
t594 = t200 * qJD(1);
t240 = (t713 + t487) * t450;
t213 = t240 * qJD(1);
t242 = t712 + (t548 + t715) * t450;
t215 = t242 * qJD(1);
t580 = 0.1e1 / 0.2e1 + t445 / 0.2e1;
t351 = (-t447 / 0.2e1 - t580) * t449;
t591 = t351 * qJD(1);
t352 = t436 / 0.2e1 + t580 * t452;
t590 = t352 * qJD(1);
t589 = t413 * qJD(1);
t588 = t450 * qJD(1);
t587 = t450 * qJD(3);
t586 = t453 * qJD(1);
t585 = t453 * qJD(3);
t584 = t453 * qJD(4);
t579 = t448 * t685;
t577 = -t701 / 0.2e1;
t576 = t700 / 0.2e1;
t575 = t698 / 0.2e1;
t574 = -t694 / 0.2e1;
t569 = pkin(4) * t704;
t568 = t710 - pkin(5) / 0.2e1;
t565 = qJ(2) * t588;
t564 = qJ(2) * t586;
t563 = t449 * t604;
t562 = t449 * t585;
t561 = t452 * t585;
t560 = t450 * t603;
t559 = t449 * t584;
t558 = t450 * t602;
t557 = t452 * t584;
t556 = t375 * t588;
t555 = t379 * t588;
t554 = t449 * t602;
t428 = t450 * t585;
t553 = t450 * t586;
t552 = t448 * t747;
t551 = t448 * t716;
t547 = -t627 / 0.2e1;
t544 = t122 / 0.2e1 + t732;
t539 = t725 + t729;
t538 = -t248 / 0.2e1 - t279 / 0.2e1;
t535 = pkin(4) * t583;
t529 = t583 * t379;
t528 = t583 * t450;
t527 = t452 * t569;
t526 = -qJD(4) - t588;
t524 = t588 + qJD(4) / 0.2e1;
t523 = t570 + t710;
t521 = t693 + t698;
t519 = qJD(3) * t534;
t457 = t122 * t716 + t738 - t739;
t8 = t577 + t457;
t514 = qJD(1) * t8;
t478 = t643 / 0.2e1 + t344 * t573;
t5 = t348 * t539 + t375 * t542 + t379 * t543 + t624 * t746 + t478;
t513 = t5 * qJD(1);
t512 = -qJD(5) + t526;
t511 = t375 * t496 - t641 / 0.2e1 + t538;
t456 = t117 * t747 + t121 * t722 + t123 * t717 + t125 * t719 - t262 * t453 / 0.2e1 + t263 * t706;
t491 = t207 * t714 + t715 * t742;
t10 = t456 + t491;
t128 = t344 * t624 + t347 * t348 - t450 * t453;
t510 = t10 * qJD(1) + t128 * qJD(2);
t33 = -t117 * t348 - t123 * t624;
t509 = qJD(1) * t33 + qJD(2) * t134;
t239 = -t375 * t438 + t379 * t693;
t463 = t532 * t707 + t693 * t718;
t486 = t641 / 0.2e1 + t538;
t42 = (t379 * t578 + t536 * t715 + t574) * t453 + t463 + t486;
t507 = -t42 * qJD(1) + t239 * qJD(3);
t295 = t549 + t384 / 0.2e1;
t58 = t550 + t261 / 0.2e1 + (t723 + t525) * t448;
t506 = qJD(1) * t58 + qJD(3) * t295;
t505 = t526 * t453;
t339 = t624 ^ 2;
t183 = t339 - t735;
t47 = qJD(1) * t183 + qJD(3) * t115;
t374 = t375 ^ 2;
t208 = t374 - t734;
t62 = qJD(1) * t115 + qJD(3) * t208;
t131 = t568 * t348 + (t552 + t547) * pkin(4);
t178 = t568 * t379 + (t551 + t709) * pkin(4);
t504 = qJD(1) * t131 + qJD(3) * t178;
t497 = pkin(5) / 0.2e1 + t523;
t140 = t497 * t624;
t188 = t497 * t375;
t503 = -qJD(1) * t140 - qJD(3) * t188;
t153 = t646 + t658;
t219 = t339 + t735;
t502 = qJD(1) * t219 + qJD(3) * t153;
t291 = t374 + t734;
t501 = qJD(1) * t153 + qJD(3) * t291;
t500 = -t606 - t610;
t498 = t692 / 0.2e1 + t690 / 0.2e1;
t495 = t515 / 0.2e1;
t494 = t453 * t521;
t474 = t498 * t449;
t286 = t382 / 0.2e1 + t474;
t493 = pkin(3) * t604 - qJD(1) * t286;
t475 = t498 * t452;
t287 = -t385 / 0.2e1 - t475;
t492 = pkin(3) * qJD(3) * t449 - qJD(1) * t287;
t489 = t655 / 0.2e1 - t648 / 0.2e1;
t50 = -t661 / 0.2e1 + t375 * t495 + t486;
t485 = qJD(1) * t50 + t375 * t605;
t483 = t452 * t505;
t100 = qJD(3) * t152 - t610 * t624;
t119 = -qJD(1) * t152 + t375 * t606;
t364 = (t444 / 0.2e1 - t446 / 0.2e1) * t453;
t482 = -qJD(1) * t364 + t563;
t481 = t379 * t495 + t640 / 0.2e1 + t517;
t480 = t125 * t573 + t437 * t731;
t479 = -t644 / 0.2e1 + t348 * t573;
t473 = qJD(1) * t436 * t449 + qJD(3) * t364;
t377 = t414 * t447;
t472 = qJD(1) * t377 + t519;
t336 = t438 + t699;
t455 = t117 * t729 + t123 * t744 + t742 * t730 + t127 * t725 + t263 * t521 / 0.2e1 + t499 * t336 / 0.2e1;
t1 = -t455 + t480;
t25 = t494 / 0.2e1 - t746 * t347 + t539 * t344 + t479;
t41 = t336 * t521;
t470 = -t1 * qJD(1) - t25 * qJD(2) + t41 * qJD(3);
t3 = -t544 * t207 - t743 + (t263 * t714 + t336 * t718 + t731) * pkin(5);
t464 = t663 / 0.2e1 + t737;
t30 = (t747 + t721) * pkin(5) + t464;
t46 = t336 * t698;
t469 = -qJD(1) * t3 + qJD(2) * t30 + qJD(3) * t46;
t145 = t706 + t489;
t460 = t117 * t713 + t123 * t715 + t207 * t721 + t717 * t742;
t465 = t439 / 0.2e1 - t695 / 0.2e1 + t449 * t571;
t23 = t460 + t465;
t91 = -t207 * t375 - t379 * t742;
t467 = -qJD(1) * t23 - qJD(2) * t145 + qJD(3) * t91;
t142 = t497 * t347;
t459 = (t123 * t705 + t448 * t544) * pkin(4) + t123 * t710;
t16 = -t697 / 0.2e1 + t459;
t458 = (t207 * t705 + t448 * t746) * pkin(4) + t207 * t710;
t32 = t696 / 0.2e1 + t458;
t350 = (-t437 + t691) * t694;
t461 = -qJD(1) * t16 - qJD(2) * t142 - qJD(3) * t32 - qJD(4) * t350;
t434 = -t586 / 0.2e1;
t433 = t586 / 0.2e1;
t432 = t585 / 0.2e1;
t427 = t452 * t588;
t426 = t449 * t587;
t425 = t449 * t588;
t372 = t524 * t453;
t359 = t364 * qJD(4);
t354 = -t628 / 0.2e1 - t436 / 0.2e1 + t452 / 0.2e1;
t353 = t709 + (t445 + t447) * t449 / 0.2e1;
t338 = (qJD(5) / 0.2e1 + t524) * t453;
t247 = t450 * t487 + t722;
t246 = t375 * t704 + t451 * t547 + t711;
t245 = t375 * t707 + t450 * t548 + t712;
t244 = -t453 * t487 + t747;
t236 = t500 * pkin(5);
t233 = -t566 + t385 / 0.2e1 - t475;
t232 = -t402 - t382 / 0.2e1 + t474;
t210 = 0.2e1 * t549 + t636;
t187 = t699 / 0.2e1 - t523 * t375;
t182 = t488 + t516;
t181 = -t487 + t490;
t177 = pkin(4) * t551 + t379 * t710 + t572 + t575;
t176 = t182 * qJD(2);
t175 = t181 * qJD(2);
t174 = t180 * qJD(2);
t173 = t179 * qJD(2);
t158 = qJD(3) * t240 + t348 * t588;
t157 = qJD(3) * t242 + t588 * t624;
t151 = t153 * qJD(6);
t146 = t706 - t489;
t144 = -t529 - t213;
t143 = -t375 * t583 - t215;
t141 = t347 * t523 + t577;
t139 = t702 / 0.2e1 - t523 * t624;
t130 = pkin(4) * t552 + t348 * t710 + t527 + t576;
t129 = t134 * qJD(6);
t109 = qJD(3) * t243 + t596;
t108 = qJD(3) * t241 + t597;
t103 = t247 * qJD(3) + t348 * t512;
t102 = t245 * qJD(3) + t512 * t624;
t65 = t244 * qJD(3) - t347 * t583 - t597;
t64 = t246 * qJD(3) + t344 * t583 - t596;
t61 = t451 * t571 + t638 + t724 - t630 / 0.2e1;
t59 = 0.2e1 * t550 - t639 / 0.2e1 + (t571 + t723) * t448;
t52 = -t660 / 0.2e1 + t481;
t51 = t661 / 0.2e1 + t511;
t45 = t379 * t527 + t453 * t574 - t463 + t511;
t44 = t375 * t527 + t451 * t569 - t462 + t481;
t31 = -t696 / 0.2e1 + t458;
t29 = pkin(5) * t747 - t702 / 0.2e1 + t464;
t26 = -t207 * t722 - t494 / 0.2e1 + t479 + t737;
t24 = -t460 + t465;
t15 = t697 / 0.2e1 + t459;
t13 = pkin(5) * t546 - t699 / 0.2e1 + t466;
t12 = t127 * t719 + t126 * t722 - t476 / 0.2e1 + t477 + t687;
t9 = t456 - t491;
t7 = t701 / 0.2e1 + t457;
t6 = t126 * t714 + t127 * t716 - t207 * t718 + t478 + t738;
t4 = pkin(5) * t731 + t122 * t725 + t207 * t732 + t263 * t575 + t336 * t576 + t743;
t2 = t455 + t480;
t27 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t428, t413 * qJD(3), 0, 0, 0, qJ(2) * t585 + t607, -qJ(2) * t587 + qJD(2) * t453, -t428 * t446 - t447 * t554, -qJD(4) * t377 + t450 * t519, -qJD(3) * t381 - t450 * t559, qJD(3) * t378 - t450 * t557, t428, qJD(3) * t155 + qJD(4) * t251 + t452 * t607, -qJD(3) * t156 - qJD(4) * t250 - t449 * t607 (-qJD(3) * t347 - t583 * t624) * t348, t136 * qJD(3) + t183 * t583, t200 * qJD(3) - t528 * t624, t199 * qJD(3) - t348 * t528, t428, qJD(3) * t37 + qJD(4) * t79 + qJD(5) * t93 - t375 * t607, -qJD(3) * t38 - qJD(4) * t80 - qJD(5) * t92 - t379 * t607, qJD(2) * t154 + qJD(3) * t18 + qJD(4) * t20 + qJD(5) * t17 + qJD(6) * t219, qJD(2) * t36 + qJD(3) * t19 + qJD(4) * t22 + qJD(5) * t21 + qJD(6) * t33; 0, 0, 0, 0, qJD(1), t615, 0, 0, 0, 0, 0, t588, t586, 0, 0, 0, 0, 0, qJD(4) * t354 + t427, qJD(4) * t353 - t425, 0, 0, 0, 0, 0, t182 * t583 - t556, t181 * t583 - t555, t613, t650 + (t344 * t375 + t379 * t347) * qJD(2) + t9 * qJD(3) + t12 * qJD(4) + t13 * qJD(5) + t129; 0, 0, 0, 0, 0, 0, -t553, t589, -t587, -t585, 0, -t454 * t587 + t564, -t454 * t585 - t565, -t359 + (-t446 * t586 - t563) * t450, t450 * t471 - 0.2e1 * t453 * t554, t562 - t608, t561 + t609, t372, t599 + (t449 * t522 - t567) * qJD(3) + t233 * qJD(4), -t598 + (-pkin(8) * t627 + (pkin(3) * t452 + t631) * t450) * qJD(3) + t232 * qJD(4), t347 * t500 + t619, t600 + (t648 + t655) * qJD(3) + t621, t245 * t583 + t379 * t585 + t594, t247 * t583 - t375 * t585 + t595, t338, t649 + (-t344 * t438 + t373 * t375 - t453 * t532) * qJD(3) + t44 * qJD(4) + t52 * qJD(5), -t645 + (-t304 * t453 - t347 * t438 + t373 * t379) * qJD(3) + t45 * qJD(4) + t51 * qJD(5), t668 + (-t121 * t379 - t125 * t375 + t663 + t664) * qJD(3) + t6 * qJD(4) + t7 * qJD(5) + t151, t667 + t9 * qJD(2) + (t121 * t742 + t125 * t207 + t262 * t336) * qJD(3) + t2 * qJD(4) + t4 * qJD(5) + t24 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t473, -t472, t449 * t505, t483, t432, qJD(2) * t354 + qJD(3) * t233 - qJD(4) * t316 + t611, qJD(2) * t353 + qJD(3) * t232 + qJD(4) * t315 - t612, t100, t47, t102, t103, t432, qJD(3) * t44 + qJD(4) * t147 + qJD(5) * t59 + t176 + t677, qJD(3) * t45 - qJD(4) * t148 + qJD(5) * t61 + t175 - t676, t666 + t6 * qJD(3) + (-t348 * t694 + t644) * qJD(4) + t139 * qJD(5), t662 + t12 * qJD(2) + t2 * qJD(3) + (t126 * t437 + t127 * t694) * qJD(4) + t15 * qJD(5) + t130 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t47, t102, t103, t432, qJD(3) * t52 + qJD(4) * t59 - qJD(5) * t138 + t176 + t674, qJD(3) * t51 + qJD(4) * t61 + qJD(5) * t137 + t175 - t675, qJD(3) * t7 + qJD(4) * t139 + t624 * t684 + t679, qJD(2) * t13 + qJD(3) * t4 + qJD(4) * t15 - t123 * t684 + t678; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t502, qJD(3) * t24 + qJD(4) * t130 + t509; 0, 0, 0, 0, -qJD(1), -t615, 0, 0, 0, 0, 0, -t588, -t586, 0, 0, 0, 0, 0, -qJD(4) * t352 - t427, -qJD(4) * t351 + t425, 0, 0, 0, 0, 0, -t179 * t583 + t556, -t180 * t583 + t555, -t613, qJD(3) * t10 - qJD(4) * t11 + qJD(5) * t14 + t129 - t650; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t587, -t585, 0, 0, 0, 0, 0, -t452 * t587 - t559, t426 - t557, 0, 0, 0, 0, 0, t244 * t583 + t375 * t587, t246 * t583 + t379 * t587, -qJD(3) * t154 (t336 * t450 + t653 - t665) * qJD(3) + t26 * qJD(4) + t29 * qJD(5) + t146 * qJD(6) + t510; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t558 - t562 - t590, t560 - t561 - t591, 0, 0, 0, 0, 0, t65, t64, 0, -t673 + t26 * qJD(3) + (-t344 * t694 - t643) * qJD(4) + t141 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t64, 0, qJD(3) * t29 + qJD(4) * t141 - t347 * t684 + t680; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t146 + t614; 0, 0, 0, 0, 0, 0, t553, -t589, 0, 0, 0, -t564, t565, t446 * t553 - t359, t483 * t745, t558 + t608, -t560 - t609, -t372, qJD(4) * t287 - t599, qJD(4) * t286 + t598, t347 * t610 + t619, -t600 + t621, -t242 * t583 - t594, -t240 * t583 - t595, -t338, -qJD(4) * t43 - qJD(5) * t49 - t649, -qJD(4) * t42 - qJD(5) * t50 + t645, -qJD(4) * t5 + qJD(5) * t8 + t151 - t668, -qJD(2) * t10 - qJD(4) * t1 - qJD(5) * t3 - qJD(6) * t23 - t667; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t583 * t241, -t583 * t243, 0, -qJD(4) * t25 + qJD(5) * t30 - qJD(6) * t145 - t510; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t554, t414 * qJD(4), 0, 0, 0, -pkin(3) * t603, -pkin(3) * t602, -t375 * t529, t583 * t208, 0, 0, 0, qJD(4) * t238 + t379 * t601, qJD(4) * t239 - t375 * t601, qJD(6) * t291, qJD(4) * t41 + qJD(5) * t46 + qJD(6) * t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t482, -t471, t427 + t602, -t425 - t603, t434, -pkin(8) * t602 - t492, pkin(8) * t603 - t493, -t119, t62, t143, t144, t434, -qJD(4) * t304 + t210 * qJD(5) - t740, t507 + t736 (-t379 * t694 + t642) * qJD(4) + t187 * qJD(5) - t513 (-t207 * t437 + t694 * t742) * qJD(4) + t31 * qJD(5) + t177 * qJD(6) + t470; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, t62, t143, t144, t434, qJD(4) * t210 - qJD(5) * t304 - t741, -t485 + t736, qJD(4) * t187 + t375 * t684 + t514, qJD(4) * t31 - t207 * t684 + t469; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t501, qJD(4) * t177 + t467; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t473, t472 (t449 * t586 - t604) * t450, t452 * t553 + t426, t432, qJD(2) * t352 - qJD(3) * t287 - t611, qJD(2) * t351 - qJD(3) * t286 + t612, -t100, -t47, t157, t158, t432, qJD(3) * t43 + qJD(5) * t58 + t173 - t677, qJD(3) * t42 + qJD(5) * t60 + t174 + t676, qJD(3) * t5 - qJD(5) * t140 - t666, qJD(2) * t11 + qJD(3) * t1 + qJD(5) * t16 + qJD(6) * t131 - t662; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t590, t591, 0, 0, 0, 0, 0, t108, t109, 0, qJD(3) * t25 + qJD(5) * t142 + t673; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t482, t471, -t427, t425, t433, t492, t493, t119, -t62, t215, t213, t433, qJD(5) * t295 + t740, t592 - t507, -qJD(5) * t188 + t513, qJD(5) * t32 + qJD(6) * t178 - t470; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t579, -t451 * t685, 0, t350 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t448 * t535 + t506, -t451 * t535 + t622, t503, -pkin(5) * t579 - t461; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t504; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t47, t157, t158, t432, qJD(3) * t49 - qJD(4) * t58 + t173 - t674, qJD(3) * t50 - qJD(4) * t60 + t174 + t675, -qJD(3) * t8 + qJD(4) * t140 - t679, -qJD(2) * t14 + qJD(3) * t3 - qJD(4) * t16 - t348 * t683 - t678; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, t109, 0, -qJD(3) * t30 - qJD(4) * t142 - t680; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, -t62, t215, t213, t433, -qJD(4) * t295 + t741, t592 + t485, qJD(4) * t188 - t514, -qJD(4) * t32 - t379 * t683 - t469; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t448 * t686 - t506, t451 * t686 - t622, -t503, t461; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t502, qJD(3) * t23 - qJD(4) * t131 + t348 * t684 - t509; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t145 - t614; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t501, -qJD(4) * t178 + t379 * t684 - t467; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t504; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t27;
