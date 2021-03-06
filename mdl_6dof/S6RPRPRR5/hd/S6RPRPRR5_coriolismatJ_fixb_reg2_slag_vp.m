% Calculate inertial parameters regressor of coriolis matrix for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRPRR5_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:49:46
% EndTime: 2019-03-09 03:50:10
% DurationCPUTime: 20.11s
% Computational Cost: add. (18067->621), mult. (34051->799), div. (0->0), fcn. (39900->8), ass. (0->496)
t576 = qJD(3) - qJD(5);
t442 = sin(qJ(3));
t439 = cos(pkin(10));
t730 = pkin(7) + qJ(2);
t533 = t730 * t439;
t438 = sin(pkin(10));
t534 = t730 * t438;
t740 = cos(qJ(3));
t363 = -t442 * t534 + t533 * t740;
t405 = -t438 * t442 + t439 * t740;
t272 = -t405 * pkin(8) + t363;
t441 = sin(qJ(5));
t362 = t442 * t533 + t534 * t740;
t563 = t740 * t438;
t674 = t442 * t439;
t406 = t563 + t674;
t473 = pkin(8) * t406 - t362;
t739 = cos(qJ(5));
t179 = -t272 * t739 + t441 * t473;
t836 = t576 * t179;
t443 = cos(qJ(6));
t433 = qJD(6) * t443;
t781 = t405 * t739 + t406 * t441;
t800 = t443 * t781;
t816 = t800 / 0.2e1;
t823 = 0.2e1 * t816;
t824 = t823 * qJD(1);
t835 = t433 + t824;
t690 = t781 ^ 2;
t382 = t739 * t406;
t676 = t441 * t405;
t346 = -t382 + t676;
t806 = t346 ^ 2;
t575 = t806 - t690;
t830 = t575 * t443;
t834 = qJD(1) * t830;
t440 = sin(qJ(6));
t831 = t575 * t440;
t833 = t831 * qJD(1);
t798 = t441 * t272 + t473 * t739;
t832 = t576 * t798;
t821 = t179 * t443;
t760 = -t821 / 0.2e1;
t811 = t798 * t440;
t164 = t811 / 0.2e1;
t165 = -t811 / 0.2e1;
t829 = t164 + t165;
t810 = t798 * t443;
t167 = t810 / 0.2e1;
t168 = -t810 / 0.2e1;
t828 = t167 + t168;
t817 = -t800 / 0.2e1;
t819 = t816 + t817;
t827 = qJD(6) * t819;
t826 = qJD(6) * t823;
t825 = t575 * qJD(1);
t736 = t179 * pkin(5);
t822 = t179 * t440;
t820 = t179 * t798;
t342 = t676 / 0.2e1 - t382 / 0.2e1;
t793 = t346 * qJD(1);
t813 = t781 * t793;
t476 = qJD(6) * t342 + t813;
t497 = -t179 * t346 + t781 * t798;
t818 = 0.2e1 * t440;
t437 = t443 ^ 2;
t801 = t437 * t781;
t317 = -t801 / 0.2e1;
t318 = t801 / 0.2e1;
t802 = t346 * pkin(5);
t815 = -t802 / 0.2e1;
t743 = t441 / 0.2e1;
t166 = t798 * t743;
t428 = -pkin(2) * t439 - pkin(1);
t721 = qJ(4) * t406;
t486 = -t428 + t721;
t761 = pkin(3) + pkin(4);
t237 = t405 * t761 + t486;
t737 = pkin(9) * t346;
t512 = pkin(5) * t781 + t737;
t448 = t237 + t512;
t80 = t440 * t448 - t821;
t726 = t80 * t443;
t79 = -t443 * t448 - t822;
t727 = t79 * t440;
t809 = t726 + t727;
t734 = t781 * pkin(9);
t231 = -t802 + t734;
t401 = t405 ^ 2;
t765 = t406 ^ 2;
t273 = t765 - t401;
t808 = t273 * qJD(1);
t807 = t273 * qJD(3);
t805 = -t346 / 0.2e1;
t803 = -t781 / 0.2e1;
t757 = t781 / 0.2e1;
t503 = t440 * t80 - t443 * t79;
t799 = t503 * t781;
t637 = qJD(5) * t346;
t489 = qJD(3) * t346 - t637;
t638 = qJD(5) * t781;
t488 = qJD(3) * t781 - t638;
t651 = qJD(2) * t346;
t780 = t401 + t765;
t797 = qJD(2) * t780;
t796 = qJD(2) * t781;
t794 = t342 * qJD(1);
t792 = t780 * qJD(1);
t791 = t781 * qJD(1);
t412 = qJ(4) * t739 - t441 * t761;
t410 = -pkin(9) + t412;
t564 = t437 * t739;
t436 = t440 ^ 2;
t565 = t436 * t739;
t461 = t564 / 0.2e1 + t565 / 0.2e1;
t566 = t412 * t739;
t790 = t566 / 0.2e1 - t461 * t410;
t788 = -0.2e1 * t346;
t411 = qJ(4) * t441 + t739 * t761;
t409 = pkin(5) + t411;
t423 = t437 - t436;
t786 = t423 * t576;
t222 = 0.2e1 * t817;
t562 = t739 * t781;
t328 = -t562 / 0.2e1;
t574 = -t739 / 0.2e1;
t514 = t781 * t574;
t744 = -t441 / 0.2e1;
t782 = t346 * t744 + t328 + t514;
t663 = t436 + t437;
t452 = t461 * pkin(9);
t748 = t436 / 0.2e1;
t536 = t437 / 0.2e1 + t748;
t755 = -t409 / 0.2e1;
t764 = -pkin(5) / 0.2e1;
t134 = t452 + (t411 * t536 + t755 + t764) * t441 + t790;
t408 = t564 + t565;
t365 = (-t739 + t408) * t441;
t364 = t365 * qJD(4);
t680 = t440 * t231;
t107 = t680 - t810;
t706 = t107 * t443;
t671 = t443 * t231;
t106 = t671 + t811;
t709 = t106 * t440;
t481 = t706 / 0.2e1 - t709 / 0.2e1;
t568 = t179 * t739;
t518 = t568 / 0.2e1;
t767 = t809 * t739 / 0.2e1;
t573 = t166 + t767;
t8 = t441 * t481 + t518 + t573;
t779 = -t8 * qJD(1) + t134 * qJD(3) - t364;
t398 = pkin(3) * t406;
t689 = t405 * qJ(4);
t353 = t398 - t689;
t732 = t406 * pkin(4);
t271 = -t353 - t732;
t138 = t271 - t231;
t683 = t440 * t138;
t86 = t810 + t683;
t722 = t86 * t443;
t673 = t443 * t138;
t85 = t673 - t811;
t725 = t85 * t440;
t483 = t722 / 0.2e1 - t725 / 0.2e1;
t537 = t411 / 0.2e1 + t755;
t750 = t412 / 0.2e1;
t753 = -t410 / 0.2e1;
t451 = -(t753 + t750) * t346 + t537 * t781;
t778 = -t737 / 0.2e1 + t451;
t581 = t408 * qJD(3);
t777 = -qJD(5) * t408 + t581;
t196 = t222 * t440;
t747 = -t437 / 0.2e1;
t206 = (t748 + t747) * t346;
t610 = t206 * qJD(6);
t776 = -t196 * qJD(3) + t610;
t775 = -t196 * qJD(5) - t610;
t677 = t440 * t443;
t195 = (t803 + t757) * t677;
t774 = -qJD(3) * t195 - t610;
t425 = t440 * t433;
t614 = t195 * qJD(1);
t773 = t614 - t425;
t772 = t614 + t425;
t771 = qJD(5) * t195 + t610;
t494 = t362 * t406 + t363 * t405;
t770 = qJD(1) * t494;
t769 = qJD(2) * t494;
t768 = 0.2e1 * t409;
t540 = t346 * t743;
t123 = -t540 + t782;
t555 = qJD(1) * t677;
t88 = -t206 * t576 - t555 * t806;
t763 = t85 / 0.2e1;
t762 = -t86 / 0.2e1;
t759 = -t798 / 0.2e1;
t572 = -t398 / 0.2e1;
t754 = t409 / 0.2e1;
t752 = t410 / 0.2e1;
t751 = -t411 / 0.2e1;
t749 = -t436 / 0.2e1;
t746 = -t440 / 0.2e1;
t745 = t440 / 0.2e1;
t742 = t443 / 0.2e1;
t741 = t8 * qJD(5);
t738 = pkin(3) * t405;
t731 = t441 * pkin(5);
t5 = -t79 * t85 + t80 * t86 + t820;
t729 = t5 * qJD(1);
t723 = t86 * t440;
t724 = t85 * t443;
t6 = -(t723 + t724) * t346 + t799;
t728 = t6 * qJD(1);
t13 = (-t107 / 0.2e1 + t762) * t443 + (t106 / 0.2e1 + t763) * t440;
t720 = qJD(1) * t13;
t699 = t798 * t346;
t18 = t781 * t809 + t699;
t719 = qJD(1) * t18;
t25 = t503 * t406;
t718 = qJD(1) * t25;
t445 = t737 / 0.2e1 + pkin(5) * t757 + t451;
t31 = t445 * t443;
t717 = qJD(1) * t31;
t678 = t346 * t440;
t51 = t678 * t798 + t781 * t79;
t716 = qJD(1) * t51;
t668 = t443 * t346;
t52 = -t668 * t798 - t781 * t80;
t715 = qJD(1) * t52;
t60 = -t179 * t781 + t699;
t714 = qJD(1) * t60;
t171 = -t568 / 0.2e1;
t90 = t171 + t518;
t713 = qJD(1) * t90;
t711 = qJD(4) * t90;
t707 = t107 * t440;
t708 = t106 * t443;
t10 = -(t707 + t708) * t346 - t799;
t710 = t10 * qJD(1);
t14 = -t346 * t79 + t440 * t497 + t781 * t85;
t705 = t14 * qJD(1);
t15 = -t80 * t346 + t443 * t497 - t781 * t86;
t704 = t15 * qJD(1);
t16 = -(-t79 - t822) * t346 + (t106 - t811) * t781;
t703 = t16 * qJD(1);
t511 = t536 * t781;
t456 = t346 * t754 + t410 * t511;
t484 = t724 / 0.2e1 + t723 / 0.2e1;
t19 = t456 + t484;
t696 = t19 * qJD(1);
t27 = t237 * t271;
t694 = t27 * qJD(1);
t28 = t760 + t821 / 0.2e1 + t445 * t440;
t693 = t28 * qJD(1);
t458 = pkin(9) * t511 + t815;
t482 = -t708 / 0.2e1 - t707 / 0.2e1;
t32 = t458 - t482;
t692 = t32 * qJD(1);
t688 = t409 * t441;
t687 = t411 * t441;
t686 = t436 * t781;
t679 = t440 * t781;
t478 = -t346 * t755 + t752 * t781;
t464 = t138 / 0.2e1 + t478;
t47 = t440 * t464 + t828;
t667 = t47 * qJD(1);
t49 = -t443 * t464 + t829;
t666 = t49 * qJD(1);
t340 = -t486 - t738;
t99 = t340 * t353;
t664 = t99 * qJD(1);
t420 = t438 ^ 2 + t439 ^ 2;
t110 = t237 * t346 + t271 * t781;
t662 = qJD(1) * t110;
t111 = t237 * t781 - t271 * t346;
t661 = qJD(1) * t111;
t306 = t686 / 0.2e1;
t543 = t781 * t749;
t477 = t543 + t317;
t120 = t318 + t306 + t477;
t660 = qJD(1) * t120;
t130 = -t690 - t806;
t127 = t130 * t443;
t658 = qJD(1) * t127;
t173 = t340 * t406 - t353 * t405;
t656 = qJD(1) * t173;
t174 = -t340 * t405 - t353 * t406;
t655 = qJD(1) * t174;
t650 = qJD(3) * qJ(4);
t646 = qJD(3) * t362;
t645 = qJD(3) * t440;
t644 = qJD(3) * t443;
t515 = t562 / 0.2e1;
t225 = t328 + t515;
t643 = qJD(4) * t225;
t642 = qJD(4) * t406;
t641 = qJD(4) * t408;
t640 = qJD(4) * t441;
t635 = qJD(5) * t440;
t634 = qJD(5) * t443;
t633 = qJD(6) * t440;
t539 = t346 / 0.2e1 + t805;
t122 = t441 * t539 + t514 + t515;
t112 = t122 * t440;
t632 = t112 * qJD(1);
t114 = t122 * t443;
t631 = t114 * qJD(1);
t567 = t346 * t739;
t450 = t441 * t511 - t567 / 0.2e1;
t116 = t406 * t536 + t450;
t630 = t116 * qJD(1);
t307 = -t686 / 0.2e1;
t119 = t317 + t307 + t477;
t629 = t119 * qJD(1);
t628 = t122 * qJD(1);
t125 = t130 * t440;
t626 = t125 * qJD(1);
t384 = t398 / 0.2e1;
t470 = t384 - t689 / 0.2e1 + t732 / 0.2e1;
t479 = t346 * t751 + t412 * t803;
t129 = t470 + t479;
t625 = t129 * qJD(1);
t624 = t130 * qJD(1);
t149 = t663 * t406 * t346;
t620 = t149 * qJD(1);
t466 = t739 * t805 + t743 * t781;
t457 = t406 / 0.2e1 + t466;
t150 = t457 * t440;
t619 = t150 * qJD(1);
t153 = t457 * t443;
t618 = t153 * qJD(1);
t465 = -t674 / 0.2e1 - t563 / 0.2e1;
t468 = t567 / 0.2e1 + t781 * t744;
t185 = t465 + t468;
t617 = t185 * qJD(1);
t282 = t443 * t540;
t194 = -t668 * t744 - t282;
t615 = t194 * qJD(1);
t611 = t206 * qJD(1);
t330 = t346 * t745;
t207 = 0.2e1 * t330;
t609 = t207 * qJD(1);
t209 = t757 * t818;
t201 = t209 * qJD(1);
t310 = -t679 / 0.2e1;
t211 = 0.2e1 * t310;
t608 = t211 * qJD(1);
t214 = t539 * t440;
t607 = t214 * qJD(1);
t538 = 0.2e1 * t805;
t216 = t538 * t440;
t606 = t216 * qJD(1);
t217 = t539 * t443;
t605 = t217 * qJD(1);
t218 = t538 * t443;
t604 = t218 * qJD(1);
t602 = t222 * qJD(1);
t600 = t225 * qJD(1);
t229 = t663 * t781;
t599 = t229 * qJD(1);
t230 = t423 * t806;
t598 = t230 * qJD(1);
t261 = t689 + 0.2e1 * t572;
t595 = t261 * qJD(1);
t355 = t363 * qJD(3);
t585 = t765 * qJD(1);
t584 = t405 * qJD(1);
t391 = t405 * qJD(3);
t583 = t405 * qJD(4);
t582 = t406 * qJD(1);
t393 = t406 * qJD(3);
t413 = t420 * qJ(2);
t580 = t413 * qJD(1);
t579 = t420 * qJD(1);
t578 = t423 * qJD(6);
t577 = t441 * qJD(3);
t561 = t237 * t791;
t560 = t237 * t793;
t559 = t237 * t582;
t558 = t436 * t793;
t557 = t437 * t793;
t556 = t428 * t582;
t554 = t440 * t644;
t553 = t781 * t642;
t552 = t440 * t634;
t551 = qJD(6) * t781 * t346;
t548 = t346 * t791;
t547 = t781 * t582;
t546 = t346 * t582;
t361 = t405 * t582;
t360 = t405 * t393;
t545 = t443 * t793;
t542 = t679 / 0.2e1;
t541 = t781 * t746;
t535 = qJD(6) * t739;
t532 = t739 * qJD(3);
t531 = t739 * qJD(4);
t356 = t663 * t411;
t530 = t576 * t440;
t529 = t576 * t443;
t525 = t440 * t547;
t524 = t443 * t547;
t523 = t806 * t425;
t522 = t346 * t425;
t521 = t346 * t555;
t520 = t318 + t543;
t519 = t179 * t574;
t513 = t764 + t537;
t510 = 0.2e1 * t521;
t509 = -0.2e1 * t521;
t507 = -qJD(5) * t218 - t346 * t644;
t506 = t443 * t530;
t505 = t440 * t529;
t9 = -t106 * t79 + t107 * t80 - t820;
t504 = t9 * qJD(1) + t8 * qJD(4);
t502 = -t722 + t725;
t480 = -t179 * t754 + t750 * t798;
t2 = t736 / 0.2e1 + (pkin(9) * t762 + t107 * t752 + t751 * t80) * t443 + (pkin(9) * t763 + t106 * t753 + t751 * t79) * t440 + t480;
t205 = -t356 * t410 + t409 * t412;
t500 = t2 * qJD(1) + t205 * qJD(3);
t236 = t408 * t410 + t688;
t3 = t519 + (t759 + t483) * t441 - t767;
t499 = qJD(1) * t3 - qJD(3) * t236;
t498 = t706 - t709;
t495 = -t346 * t410 + t409 * t781;
t357 = t566 + t687;
t53 = t519 + t518 + (t798 / 0.2e1 + t759) * t441;
t493 = -qJD(1) * t53 + qJD(3) * t357;
t17 = -(-t80 - t821) * t346 + (-t107 - t810) * t781;
t492 = t17 * qJD(1) + t194 * qJD(4);
t491 = t346 * (-qJD(6) - t791);
t490 = qJD(3) * t412 + t640;
t487 = qJD(5) * t412 + t640;
t485 = t734 / 0.2e1 + t815;
t474 = qJD(3) * t409 - t531;
t472 = t489 * t781;
t471 = t488 * t346;
t469 = t231 / 0.2e1 + t485;
t415 = -qJD(5) * t739 + t532;
t336 = t513 * t440;
t58 = -t443 * t469 + t829;
t460 = pkin(5) * t635 - qJD(1) * t58 + qJD(3) * t336;
t337 = t513 * t443;
t56 = t440 * t469 + t828;
t459 = pkin(5) * t634 - qJD(1) * t56 + qJD(3) * t337;
t424 = t440 * t577;
t454 = -t441 * t635 + t443 * t535 + t424;
t426 = t443 * t577;
t453 = -t440 * t535 - t441 * t634 + t426;
t414 = t576 * t441;
t395 = t406 * qJD(2);
t352 = t384 + t572;
t339 = t768 * t742;
t338 = t768 * t745;
t276 = 0.2e1 * t522;
t275 = -0.2e1 * t522;
t233 = t510 - t786;
t232 = t509 + t786;
t228 = t576 * t342;
t215 = -t678 / 0.2e1 + t330;
t213 = t678 / 0.2e1 + t330;
t212 = t310 + t542;
t210 = t542 + t541;
t204 = t225 * qJD(5);
t203 = t217 * qJD(5);
t200 = t210 * qJD(6);
t199 = t209 * qJD(6);
t190 = t194 * qJD(5);
t186 = t201 + t633;
t184 = t465 - t468;
t178 = -t506 + t611;
t177 = t505 - t611;
t152 = t406 * t742 - t443 * t466;
t151 = t406 * t746 + t440 * t466;
t135 = t688 / 0.2e1 - t731 / 0.2e1 + t452 - t663 * t687 / 0.2e1 - t790;
t128 = t470 - t479;
t121 = t318 + t307 + t520;
t118 = t317 + t306 + t520;
t117 = t450 + (t747 + t749) * t406;
t115 = t443 * t782 - t282;
t113 = t123 * t440;
t89 = t90 * qJD(5);
t59 = 0.2e1 * t164 + t671 / 0.2e1 - t485 * t443;
t57 = 0.2e1 * t167 - t680 / 0.2e1 + t485 * t440;
t54 = 0.2e1 * t166 + t171 + t519;
t50 = 0.2e1 * t165 + t673 / 0.2e1 - t478 * t443;
t48 = 0.2e1 * t168 - t683 / 0.2e1 + t478 * t440;
t33 = t458 + t482;
t30 = pkin(5) * t817 + t443 * t778 + t822;
t29 = pkin(5) * t541 + t440 * t778 + 0.2e1 * t760;
t20 = t456 - t484;
t12 = -t481 + t483;
t4 = t441 * t483 + t519 + t573;
t1 = -t736 / 0.2e1 + (-t726 / 0.2e1 - t727 / 0.2e1) * t411 + t481 * t410 + t480 + t483 * pkin(9);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t420 * qJD(2), t413 * qJD(2), t360, -t807, 0, -t360, 0, 0, t428 * t393, t428 * t391, t797, t769, t360, 0, t807, 0, 0, -t360, qJD(3) * t173 + t406 * t583, t797, qJD(3) * t174 + qJD(4) * t765, qJD(3) * t99 - t340 * t642 + t769, -t471, t576 * t575, 0, t472, 0, 0, qJD(3) * t110 - t237 * t637 + t553, qJD(3) * t111 - t237 * t638 - t346 * t642, qJD(2) * t130, qJD(2) * t60 + qJD(3) * t27 + t237 * t642, -t437 * t471 - t523, t488 * t668 * t818 - t230 * qJD(6), t440 * t551 - t576 * t830, -t436 * t471 + t523, t443 * t551 + t576 * t831, t472, qJD(2) * t125 + qJD(3) * t14 + qJD(5) * t16 + qJD(6) * t52 + t443 * t553, qJD(2) * t127 + qJD(3) * t15 + qJD(5) * t17 + qJD(6) * t51 - t440 * t553, -qJD(3) * t6 + qJD(4) * t149 - qJD(5) * t10, qJD(2) * t18 + qJD(3) * t5 + qJD(4) * t25 + qJD(5) * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t579, t580, 0, 0, 0, 0, 0, 0, 0, 0, t792, t770, 0, 0, 0, 0, 0, 0, 0, t792, 0, qJD(3) * t352 + t770, 0, 0, 0, 0, 0, 0, 0, 0, t624, qJD(3) * t128 + qJD(4) * t184 + t714, 0, 0, 0, 0, 0, 0, t200 + t203 + t626, qJD(5) * t215 + t658 + t827, t118 * qJD(5), qJD(3) * t20 + qJD(4) * t117 + qJD(5) * t33 + t719; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t361, -t808, t391, -t361, -t393, 0, -t355 + t556, t428 * t584 + t646, 0, 0, t361, t391, t808, 0, t393, -t361, -t355 + t656 (-t721 - t738) * qJD(3) + t583, -t646 + t655, t664 + t352 * qJD(2) + (-pkin(3) * t363 - qJ(4) * t362) * qJD(3) + t363 * qJD(4), -t548, t825, -t488, t813, t489, 0, t662 + t836, t661 + t832 (-t346 * t412 + t411 * t781) * qJD(3) + t123 * qJD(4), t694 + t128 * qJD(2) + (t179 * t411 + t412 * t798) * qJD(3) + t54 * qJD(4) (-t554 - t557) * t781 + t775, t121 * qJD(5) + t275 + (-qJD(3) * t423 + t510) * t781, qJD(5) * t213 - t346 * t645 + t827 - t834 (t554 - t558) * t781 - t775, t200 + t507 + t833, t476, t705 + (t440 * t495 + t821) * qJD(3) + t113 * qJD(4) + t29 * qJD(5) + t50 * qJD(6), t704 + (t443 * t495 - t822) * qJD(3) + t115 * qJD(4) + t30 * qJD(5) + t48 * qJD(6), qJD(3) * t502 + qJD(5) * t12 - t728, t729 + t20 * qJD(2) + (t179 * t409 - t410 * t502) * qJD(3) + t4 * qJD(4) + t1 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t361, t391, t585, -t340 * t582 + t355, 0, 0, 0, 0, 0, 0, t547, -t546, qJD(3) * t123 + t204, qJD(2) * t184 + qJD(3) * t54 + t559 + t89, 0, 0, 0, 0, 0, 0, qJD(3) * t113 + qJD(6) * t152 + t524, qJD(3) * t115 + qJD(6) * t151 + t190 - t525, t620, qJD(2) * t117 + qJD(3) * t4 + t718 + t741; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t813, -t825, t488, -t813, -t489, 0, -t560 - t836, -t561 - t832, t643, t711 (-t552 + t557) * t781 + t776, t121 * qJD(3) + t276 + (-qJD(5) * t423 + t509) * t781, qJD(3) * t213 - t346 * t635 + t827 + t834 (t552 + t558) * t781 - t776, -qJD(3) * t218 + qJD(6) * t212 - t346 * t634 - t833, -t476, t703 + t217 * qJD(2) + t29 * qJD(3) + (t440 * t512 + t821) * qJD(5) + t59 * qJD(6), t215 * qJD(2) + t30 * qJD(3) + (t443 * t512 - t822) * qJD(5) + t57 * qJD(6) + t492, qJD(2) * t118 + qJD(3) * t12 + qJD(5) * t498 - t710, t33 * qJD(2) + t1 * qJD(3) + (pkin(9) * t498 + t736) * qJD(5) + t504; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, t505 * t788 - t598, -t440 * t491 + (qJD(3) + qJD(5)) * t819, -t88, t210 * qJD(3) + t212 * qJD(5) - t443 * t491, t228, qJD(2) * t210 + qJD(3) * t50 + qJD(4) * t152 + qJD(5) * t59 - qJD(6) * t80 + t715, qJD(2) * t819 + qJD(3) * t48 + qJD(4) * t151 + qJD(5) * t57 + qJD(6) * t79 + t716, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t579, -t580, 0, 0, 0, 0, 0, 0, t393, t391, -t792, -t770, 0, 0, 0, 0, 0, 0, t393, -t792, -t391, -qJD(3) * t261 - t642 - t770, 0, 0, 0, 0, 0, 0, -t489, -t488, -t624, qJD(3) * t129 + qJD(4) * t185 - t714, 0, 0, 0, 0, 0, 0, t199 + t507 - t626, qJD(3) * t207 + qJD(5) * t216 - t658 + t826, qJD(3) * t229 + qJD(5) * t119, -qJD(3) * t19 - qJD(4) * t116 - qJD(5) * t32 - t719; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t582, t584, 0, 0, 0, 0, 0, 0, 0, 0, t582, 0, -t584, -t595, 0, 0, 0, 0, 0, 0, -t793, -t791, 0, t625, 0, 0, 0, 0, 0, 0, -t545, t609, t599, -t696; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t582, 0, 0, 0, 0, 0, 0, 0, 0, 0, t617, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t630; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t793, t791, 0, 0, 0, 0, 0, 0, 0, 0, -t604, t606, t629, -t692; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t186, t835, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t361, t808, 0, t361, 0, 0, -t395 - t556 (-qJD(1) * t428 - qJD(2)) * t405, 0, 0, -t361, 0, -t808, 0, 0, t361, -t395 - t656, 0, qJD(2) * t405 - t655, qJD(2) * t261 - t664, t548, -t825, 0, -t813, 0, 0, t651 - t662, -t661 + t796, -qJD(4) * t122, -qJD(2) * t129 - qJD(4) * t53 - t694, t437 * t548 - t771, qJD(5) * t120 + t509 * t781 + t275, qJD(5) * t214 + qJD(6) * t222 + t834, t436 * t548 + t771, t199 + t203 - t833, -t476, -qJD(4) * t112 + qJD(5) * t28 + qJD(6) * t49 + t443 * t651 - t705, -qJD(2) * t207 - qJD(4) * t114 + qJD(5) * t31 + qJD(6) * t47 - t704, -qJD(2) * t229 + qJD(5) * t13 + t728, qJD(2) * t19 - qJD(4) * t3 + qJD(5) * t2 - t729; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t582, -t584, 0, 0, 0, 0, 0, 0, 0, 0, -t582, 0, t584, t595, 0, 0, 0, 0, 0, 0, t793, t791, 0, -t625, 0, 0, 0, 0, 0, 0, t545, -t609, -t599, t696; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4), 0, 0, 0, 0, 0, 0, t487, -qJD(5) * t411 + t531, 0, qJD(4) * t357, t425, t578, 0, -t425, 0, 0, -t409 * t633 + t443 * t487, -t409 * t433 - t440 * t487, qJD(5) * t356 - t641, qJD(4) * t236 + qJD(5) * t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t650, 0, 0, 0, 0, 0, 0, t577, t532, -t628, t493, 0, 0, 0, 0, 0, 0, t426 - t632, -t424 - t631, -t581, qJD(5) * t135 + t364 - t499; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t412 * t576, -t411 * t576, 0, 0, -t772, -t578 + t660, t607, t772, t605, 0, t338 * qJD(6) + t412 * t529 + t693, qJD(6) * t339 - t412 * t530 + t717, t356 * t576 + t720, t135 * qJD(4) + (-pkin(5) * t412 - pkin(9) * t356) * qJD(5) + t500; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, t232, -t433 + t602, t178, t186, -t794, qJD(5) * t338 - t409 * t645 - t410 * t433 + t666, qJD(5) * t339 - t409 * t644 + t410 * t633 + t667, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t361, 0, -t585 (qJD(1) * t340 + qJD(2)) * t406, 0, 0, 0, 0, 0, 0, -t547, t546, qJD(3) * t122 + t204, -qJD(2) * t185 + qJD(3) * t53 - t559 + t89, 0, 0, 0, 0, 0, 0, qJD(3) * t112 - qJD(6) * t153 - t524, qJD(3) * t114 + qJD(6) * t150 + t190 + t525, -t620, qJD(2) * t116 + qJD(3) * t3 - t718 + t741; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t582, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t617, 0, 0, 0, 0, 0, 0, 0, 0, 0, t630; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t650, 0, 0, 0, 0, 0, 0, -t414, -t415, t628, -t493, 0, 0, 0, 0, 0, 0, -t453 + t632, t454 + t631, t777, -qJD(5) * t134 + t499; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t365 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t414, t415, t600, t713, 0, 0, 0, 0, 0, 0, t453, -t454 + t615, -t777 (pkin(9) * t408 - t731) * qJD(5) - t779; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t415 * t440 - t433 * t441 - t618, t415 * t443 + t441 * t633 + t619, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t813, t825, 0, t813, 0, 0, t560 - t651, t561 - t796, -t643, -t711, -t437 * t813 - t774, -qJD(3) * t120 + t510 * t781 + t276, -qJD(3) * t214 + t826 - t834, -t436 * t813 + t774, -qJD(3) * t217 + qJD(6) * t211 + t833, t476, qJD(2) * t218 - qJD(3) * t28 + qJD(6) * t58 - t703, -qJD(2) * t216 - qJD(3) * t31 + qJD(6) * t56 - t492, -qJD(2) * t119 - qJD(3) * t13 + t710, qJD(2) * t32 - qJD(3) * t2 - t504; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t793, -t791, 0, 0, 0, 0, 0, 0, 0, 0, t604, -t606, -t629, t692; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t490, qJD(3) * t411 - t531, 0, 0, t773, -t578 - t660, -t607, -t773, -t605, 0, -t336 * qJD(6) - t443 * t490 - t693, -qJD(6) * t337 + t440 * t490 - t717, -qJD(3) * t356 + t641 - t720, qJD(4) * t134 - t500; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t577, -t532, -t600, -t713, 0, 0, 0, 0, 0, 0, -t426, t424 - t615, t581, t779; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t425, t578, 0, -t425, 0, 0, -pkin(5) * t633, -pkin(5) * t433, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178, t233, t835, t177, t608 - t633, t794, -pkin(9) * t433 - t460, pkin(9) * t633 - t459, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, -t506 * t788 + t598, -qJD(3) * t222 - qJD(5) * t823 - t440 * t813, t88, -qJD(3) * t209 - qJD(5) * t211 - t443 * t813, t228, -qJD(2) * t209 - qJD(3) * t49 + qJD(4) * t153 - qJD(5) * t58 - t715, -qJD(2) * t823 - qJD(3) * t47 - qJD(4) * t150 - qJD(5) * t56 - t716, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t201, -t824, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178, t233, -t602, t177, -t201, t794, t336 * qJD(5) + t440 * t474 - t666, t337 * qJD(5) + t443 * t474 - t667, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t440 * t532 + t618, -t443 * t532 - t619, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t177, t232, -t824, t178, -t608, -t794, t460, t459, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t7;
