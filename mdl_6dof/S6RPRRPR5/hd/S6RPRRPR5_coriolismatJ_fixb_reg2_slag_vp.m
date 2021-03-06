% Calculate inertial parameters regressor of coriolis matrix for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPRRPR5_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:13:41
% EndTime: 2019-03-09 05:14:06
% DurationCPUTime: 18.79s
% Computational Cost: add. (19807->615), mult. (37808->716), div. (0->0), fcn. (45700->8), ass. (0->461)
t423 = cos(qJ(6));
t418 = t423 ^ 2;
t419 = sin(pkin(10));
t420 = cos(pkin(10));
t721 = sin(qJ(3));
t723 = cos(qJ(3));
t376 = -t419 * t723 - t420 * t721;
t722 = cos(qJ(4));
t369 = t722 * t376;
t375 = -t419 * t721 + t420 * t723;
t422 = sin(qJ(4));
t654 = t422 * t375;
t759 = -t369 + t654;
t778 = t418 * t759;
t303 = t778 / 0.2e1;
t421 = sin(qJ(6));
t417 = t421 ^ 2;
t779 = t417 * t759;
t799 = -t779 / 0.2e1;
t640 = t799 + t303;
t798 = t779 / 0.2e1;
t800 = -t778 / 0.2e1;
t829 = t640 + t798 + t800;
t834 = qJD(2) * t829;
t833 = t829 * qJD(4);
t714 = pkin(7) + qJ(2);
t515 = t714 * t420;
t516 = t714 * t419;
t354 = t515 * t721 + t516 * t723;
t261 = t376 * pkin(8) - t354;
t355 = t515 * t723 - t516 * t721;
t436 = t375 * pkin(8) + t355;
t808 = -t722 * t261 + t422 * t436;
t428 = pkin(5) * t759 + t808;
t832 = t421 * t428;
t831 = t423 * t428;
t368 = t722 * t375;
t653 = t422 * t376;
t350 = t368 + t653;
t718 = pkin(5) * t350;
t431 = t422 * t261;
t761 = t722 * t436;
t774 = t761 + t431;
t793 = t774 + t718;
t830 = t428 * t793;
t610 = qJD(6) * t421;
t777 = t421 * t759;
t802 = -t777 / 0.2e1;
t815 = 0.2e1 * t802;
t818 = t815 * qJD(1);
t828 = t818 - t610;
t415 = qJD(3) + qJD(4);
t715 = t422 * pkin(3);
t398 = qJ(5) + t715;
t729 = t398 / 0.2e1;
t811 = t423 * t793;
t735 = t811 / 0.2e1;
t787 = t759 ^ 2;
t788 = t350 ^ 2;
t792 = -t788 + t787;
t809 = t792 * t423;
t826 = qJD(1) * t809;
t810 = t792 * t421;
t825 = qJD(1) * t810;
t780 = t350 * t423;
t796 = t780 / 0.2e1;
t807 = -t780 / 0.2e1 + t796;
t824 = qJD(3) * t807;
t823 = qJD(4) * t807;
t822 = qJD(5) * t807;
t801 = t777 / 0.2e1;
t806 = t802 + t801;
t821 = qJD(6) * t806;
t820 = qJD(6) * t815;
t580 = t807 * qJD(1);
t819 = t807 * qJD(2);
t469 = t415 * t808;
t816 = 0.2e1 * t796;
t817 = qJD(4) * t816 + t821;
t767 = t759 * qJD(1);
t794 = t350 * t767;
t814 = t417 * t794;
t813 = t418 * t794;
t812 = t421 * t793;
t805 = t792 * qJD(1);
t554 = t722 * pkin(3);
t402 = -t554 - pkin(4);
t396 = -pkin(9) + t402;
t738 = pkin(4) + pkin(9);
t775 = t738 * t759;
t781 = t350 * qJ(5);
t159 = -t781 + t775;
t650 = t423 * t159;
t91 = t650 + t812;
t698 = t91 * t421;
t659 = t421 * t159;
t90 = t811 - t659;
t699 = t90 * t423;
t713 = t699 / 0.2e1 + t698 / 0.2e1;
t803 = -t713 * t396 + t428 * t729;
t795 = t774 * pkin(4);
t727 = t418 / 0.2e1;
t212 = (t727 - t417 / 0.2e1) * t350;
t763 = t415 * t423;
t776 = t421 * t763;
t173 = t212 * qJD(1) + t776;
t394 = t417 - t418;
t225 = t394 * t788;
t790 = t225 * qJD(1) - 0.2e1 * t350 * t776;
t789 = t415 * t792;
t786 = -t350 / 0.2e1;
t731 = -t759 / 0.2e1;
t730 = t759 / 0.2e1;
t784 = t759 * pkin(4);
t783 = qJ(5) * t759;
t399 = -t420 * pkin(2) - pkin(1);
t362 = -t375 * pkin(3) + t399;
t448 = t362 - t783;
t720 = pkin(4) * t350;
t176 = t448 - t720;
t782 = t176 * t350;
t673 = t176 * t759;
t528 = t350 * t729;
t648 = t738 * t350;
t146 = t448 - t648;
t70 = t146 * t421 - t831;
t71 = t423 * t146 + t832;
t24 = (t421 * t70 + t71 * t423) * t759;
t517 = t417 / 0.2e1 + t727;
t500 = t759 * t517;
t711 = qJD(3) * pkin(3);
t749 = t422 * t711;
t616 = qJD(3) * t350;
t184 = qJD(4) * t350 + t616;
t499 = (t786 + t350 / 0.2e1) * t421;
t743 = t499 * qJD(1);
t773 = t763 + t743;
t619 = qJD(2) * t350;
t341 = -t653 / 0.2e1 - t368 / 0.2e1;
t771 = qJD(6) * t341;
t770 = t341 * qJD(1);
t769 = t350 * qJD(1);
t566 = t350 * qJD(5);
t496 = -t761 / 0.2e1 - t431 / 0.2e1;
t728 = t402 / 0.2e1;
t764 = t394 * t415;
t760 = qJ(5) + t398;
t531 = -t783 / 0.2e1;
t758 = -t648 / 0.2e1 + t531;
t425 = t718 / 0.2e1 - t496;
t716 = t376 * pkin(3);
t490 = -t781 - t716;
t153 = t490 + t775;
t651 = t423 * t153;
t77 = t651 + t812;
t702 = t77 * t421;
t660 = t421 * t153;
t76 = t811 - t660;
t703 = t76 * t423;
t461 = t702 / 0.2e1 + t703 / 0.2e1;
t18 = t425 - t461;
t558 = t398 * qJD(3);
t756 = -qJD(1) * t18 - t558;
t724 = t423 / 0.2e1;
t525 = t759 * t724;
t191 = t421 * t525 + t423 * t801;
t582 = t212 * qJD(6);
t755 = -t191 * qJD(4) + t582;
t754 = -t191 * qJD(3) + t582;
t655 = t421 * t423;
t192 = (t730 + t731) * t655;
t753 = qJD(4) * t192 - t582;
t752 = -qJD(3) * t192 - t582;
t609 = qJD(6) * t423;
t395 = t421 * t609;
t590 = t192 * qJD(1);
t751 = t590 - t395;
t750 = t590 + t395;
t406 = qJD(4) * t715;
t748 = t406 + t749;
t480 = t350 * t774 + t759 * t808;
t747 = qJD(1) * t480;
t746 = qJD(2) * t480;
t479 = t788 + t787;
t744 = t479 * qJD(1);
t742 = t499 * qJD(3);
t110 = -0.2e1 * t496;
t502 = -t369 / 0.2e1;
t343 = t502 + t369 / 0.2e1;
t551 = -t343 * qJD(2) + t110 * qJD(3) + qJD(4) * t774;
t555 = -t722 / 0.2e1;
t503 = t350 * t555;
t725 = t422 / 0.2e1;
t435 = pkin(3) * (t725 * t759 - t503) - t759 * t729 - t396 * t786;
t620 = qJD(1) * t423;
t543 = t421 * t620;
t96 = -t212 * t415 + t543 * t788;
t741 = qJD(2) * t479;
t739 = t376 ^ 2;
t737 = -qJ(5) / 0.2e1;
t736 = qJ(5) / 0.2e1;
t726 = -t421 / 0.2e1;
t3 = -t70 * t76 + t71 * t77 - t830;
t710 = t3 * qJD(1);
t4 = -t70 * t90 + t71 * t91 - t830;
t709 = t4 * qJD(1);
t701 = t77 * t423;
t704 = t76 * t421;
t5 = -(t701 - t704) * t350 + t24;
t708 = t5 * qJD(1);
t697 = t91 * t423;
t700 = t90 * t421;
t6 = t24 - (t697 - t700) * t350;
t707 = t6 * qJD(1);
t706 = t70 * t423;
t705 = t71 * t421;
t14 = t793 * t350 + (t705 - t706) * t759;
t695 = qJD(1) * t14;
t693 = qJD(1) * t24;
t52 = t70 * t759 - t780 * t793;
t692 = qJD(1) * t52;
t657 = t421 * t350;
t53 = -t657 * t793 - t71 * t759;
t691 = qJD(1) * t53;
t186 = t490 + t784;
t82 = t186 * t350 - t673;
t688 = qJD(1) * t82;
t83 = -t186 * t759 - t782;
t687 = qJD(1) * t83;
t226 = -t781 + t784;
t92 = t226 * t350 - t673;
t686 = qJD(1) * t92;
t93 = -t226 * t759 - t782;
t685 = qJD(1) * t93;
t482 = t350 * t428 + t759 * t793;
t10 = -t70 * t350 - t423 * t482 + t759 * t76;
t684 = t10 * qJD(1);
t11 = -t71 * t350 + t421 * t482 - t759 * t77;
t683 = t11 * qJD(1);
t12 = (t90 - t811) * t759 - (t70 + t831) * t350;
t682 = t12 * qJD(1);
t13 = (-t91 + t812) * t759 - (t71 - t832) * t350;
t681 = t13 * qJD(1);
t442 = t396 * t500 + t528;
t462 = t704 / 0.2e1 - t701 / 0.2e1;
t15 = t442 + t462;
t675 = t15 * qJD(1);
t529 = t781 / 0.2e1;
t441 = -t500 * t738 + t529;
t460 = t700 / 0.2e1 - t697 / 0.2e1;
t21 = t441 + t460;
t672 = t21 * qJD(1);
t27 = t176 * t186;
t671 = t27 * qJD(1);
t669 = t759 * t398;
t37 = t176 * t226;
t664 = t37 * qJD(1);
t54 = t362 * t716;
t645 = t54 * qJD(1);
t552 = t715 / 0.2e1;
t501 = -t398 / 0.2e1 + t552;
t438 = -pkin(3) * t503 + t350 * t728 + t501 * t759;
t465 = pkin(4) * t786 + t531;
t94 = t438 - t465;
t644 = t94 * qJD(1);
t240 = 0.2e1 * t502 + t654;
t641 = t240 * qJD(2);
t407 = qJD(4) * t554;
t411 = qJD(5) * t421;
t639 = t421 * t407 + t411;
t412 = qJD(5) * t423;
t638 = t423 * t407 + t412;
t409 = t421 * qJD(4);
t410 = t421 * qJD(3);
t636 = t409 + t410;
t392 = t419 ^ 2 + t420 ^ 2;
t635 = t417 + t418;
t634 = qJ(5) * qJD(4);
t633 = qJD(1) * t829;
t124 = t479 * t421;
t631 = qJD(1) * t124;
t127 = t479 * t423;
t629 = qJD(1) * t127;
t160 = -t350 * t716 - t362 * t759;
t626 = qJD(1) * t160;
t161 = -t362 * t350 + t716 * t759;
t625 = qJD(1) * t161;
t622 = qJD(1) * t362;
t621 = qJD(1) * t421;
t617 = qJD(3) * t759;
t613 = qJD(4) * t759;
t612 = qJD(4) * t362;
t611 = qJD(4) * t423;
t608 = qJD(6) * t738;
t112 = t716 / 0.2e1 + (t729 + t736) * t350 + (t728 - pkin(4) / 0.2e1) * t759;
t605 = t112 * qJD(1);
t115 = 0.2e1 * t799 + 0.2e1 * t800;
t604 = t115 * qJD(1);
t116 = t798 + t303 + t500;
t603 = t116 * qJD(1);
t518 = 0.2e1 * t786;
t519 = 0.2e1 * t730;
t121 = pkin(4) * t519 + qJ(5) * t518;
t602 = t121 * qJD(1);
t125 = t635 * t759 * t350;
t601 = t125 * qJD(1);
t452 = t350 * t725 + t555 * t759;
t174 = (t376 / 0.2e1 + t452) * pkin(3);
t594 = t174 * qJD(1);
t188 = -t354 * t376 + t355 * t375;
t593 = t188 * qJD(1);
t202 = 0.2e1 * t801;
t589 = t202 * qJD(1);
t207 = t518 * t421;
t585 = t207 * qJD(1);
t213 = t519 * t423;
t198 = t213 * qJD(1);
t215 = t518 * t423;
t581 = t215 * qJD(1);
t579 = t780 * qJD(1);
t224 = -t779 - t778;
t578 = t224 * qJD(1);
t575 = t240 * qJD(1);
t374 = t375 ^ 2;
t269 = t374 - t739;
t574 = t269 * qJD(1);
t571 = t343 * qJD(1);
t570 = t343 * qJD(3);
t569 = t343 * qJD(4);
t568 = t787 * qJD(1);
t353 = t374 + t739;
t564 = t353 * qJD(1);
t563 = t375 * qJD(1);
t372 = t375 * qJD(3);
t562 = t376 * qJD(1);
t561 = t376 * qJD(3);
t379 = t392 * qJ(2);
t560 = t379 * qJD(1);
t559 = t392 * qJD(1);
t557 = t407 + qJD(5);
t553 = -t715 / 0.2e1;
t550 = t722 * t398;
t549 = t176 * t767;
t548 = t350 * t622;
t547 = t417 * t769;
t546 = t418 * t769;
t545 = t759 * t622;
t544 = t787 * t621;
t542 = t423 * t410;
t541 = t350 * t409;
t540 = t423 * t409;
t539 = t759 * t610;
t538 = t759 * t609;
t229 = t759 * t769;
t535 = t375 * t562;
t534 = t375 * t561;
t533 = t421 * t769;
t532 = t759 * t620;
t530 = -t781 / 0.2e1;
t522 = t648 / 0.2e1;
t514 = t422 * t635;
t513 = -t794 + t771;
t512 = t229 - t771;
t510 = qJD(1) * t399 + qJD(2);
t509 = qJD(6) + t767;
t405 = qJD(3) * t554;
t508 = t554 / 0.2e1;
t507 = t788 * t395;
t505 = t350 * t543;
t373 = pkin(3) * t514;
t498 = 0.2e1 * t505;
t497 = -0.2e1 * t505;
t495 = t207 * qJD(4) - t350 * t410;
t494 = -t636 + t580;
t493 = t636 + t580;
t488 = t702 + t703;
t487 = t698 + t699;
t486 = t759 * t505;
t485 = t783 + t648;
t437 = -t428 * t736 - t461 * t738;
t1 = (t793 * t555 + (-t705 / 0.2e1 + t706 / 0.2e1) * t422) * pkin(3) + t437 + t803;
t314 = (t396 * t514 + t550) * pkin(3);
t484 = -t1 * qJD(1) + t314 * qJD(3);
t8 = -t461 + t713;
t483 = qJD(1) * t8 + qJD(3) * t373;
t477 = -t350 * t396 + t669;
t439 = t808 * t501 + (t508 + t728) * t774;
t464 = t795 / 0.2e1 - t808 * t737;
t25 = t439 + t464;
t363 = -pkin(3) * t550 - t402 * t715;
t476 = t25 * qJD(1) - t363 * qJD(3);
t475 = t509 * t421;
t474 = t509 * t423;
t471 = qJD(2) * t759;
t468 = -qJD(3) * t774 - qJD(4) * t110;
t182 = qJD(4) * t240 + t617;
t466 = qJD(3) * t240 + t613;
t463 = -t784 / 0.2e1 + t529;
t459 = t737 + t501;
t457 = t396 * t731 - t528;
t456 = -t731 * t738 + t530;
t450 = t153 / 0.2e1 + t457;
t50 = t450 * t423;
t455 = -qJD(1) * t50 + t398 * t410;
t48 = t735 - t811 / 0.2e1 + t450 * t421;
t454 = -qJD(1) * t48 - t423 * t558;
t453 = (t613 + t617) * t350;
t128 = t184 * t759;
t449 = t159 / 0.2e1 + t456;
t20 = t425 - t713;
t360 = -qJ(5) + (-0.1e1 / 0.2e1 + t517) * t715;
t447 = qJD(1) * t20 - qJD(3) * t360 + t634;
t41 = (-t435 + t758) * t421;
t446 = -t41 * qJD(1) - t405 * t423;
t43 = (t522 + t783 / 0.2e1 + t435) * t423;
t445 = -t43 * qJD(1) - t405 * t421;
t100 = t449 * t421;
t357 = t459 * t423;
t444 = -qJ(5) * t611 - t100 * qJD(1) + t357 * qJD(3);
t101 = t449 * t423;
t356 = t459 * t421;
t443 = qJ(5) * t409 - t101 * qJD(1) - t356 * qJD(3);
t416 = qJ(5) * qJD(5);
t408 = qJ(5) * t554;
t393 = t415 * qJ(5);
t389 = t398 * qJD(5);
t382 = t394 * qJD(6);
t371 = -t716 / 0.2e1;
t370 = t373 * qJD(4);
t361 = t517 * t715 + qJ(5) + t552;
t359 = t423 * t552 + t724 * t760;
t358 = t421 * t553 + t726 * t760;
t270 = 0.2e1 * t350 * t395;
t231 = t498 + t764;
t230 = t497 - t764;
t223 = t415 * t341;
t222 = t423 * t731 + t525;
t204 = -t657 / 0.2e1 + t350 * t726;
t200 = t213 * qJD(6);
t199 = t222 * qJD(6);
t194 = t499 * qJD(4);
t187 = -t198 - t609;
t175 = pkin(3) * t452 + t371;
t122 = t530 + t784 / 0.2e1 + t463;
t119 = 0.2e1 * t640;
t113 = t728 * t759 + t371 - t463 + t528;
t95 = t438 + t465;
t57 = -t812 - t650 / 0.2e1 + t456 * t423;
t56 = t811 - t659 / 0.2e1 + t456 * t421;
t51 = -t812 - t651 / 0.2e1 + t457 * t423;
t49 = 0.2e1 * t735 - t660 / 0.2e1 + t457 * t421;
t42 = -t832 + (t435 + t758) * t423;
t40 = -t831 + qJ(5) * t801 + (t522 - t435) * t421;
t26 = t439 - t464;
t22 = t441 - t460;
t19 = t425 + t713;
t17 = t425 + t461;
t16 = t442 - t462;
t9 = -t461 - t713;
t2 = t793 * t508 + t552 * t705 + t553 * t706 + t437 - t803;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t392 * qJD(2), t379 * qJD(2), -t534, t269 * qJD(3), 0, t534, 0, 0, -t399 * t561, t399 * t372, qJD(2) * t353, qJD(2) * t188, t128, -t789, 0, -t453, 0, 0, -qJD(3) * t160 + t612 * t759, -qJD(3) * t161 + t350 * t612, t741, -qJD(3) * t54 + t746, 0, 0, 0, t128, -t789, -t453, t741, qJD(3) * t82 + qJD(4) * t92 - t566 * t759, qJD(3) * t83 + qJD(4) * t93 + qJD(5) * t787, qJD(3) * t27 + qJD(4) * t37 - qJD(5) * t673 + t746, -t417 * t453 + t507, -t225 * qJD(6) - 0.2e1 * t453 * t655, -t350 * t538 + t415 * t810, -t418 * t453 - t507, t350 * t539 + t415 * t809, t128, qJD(2) * t127 + qJD(3) * t10 + qJD(4) * t12 + qJD(6) * t53 + t411 * t787, -qJD(2) * t124 + qJD(3) * t11 + qJD(4) * t13 + qJD(6) * t52 + t412 * t787, qJD(3) * t5 + qJD(4) * t6 + qJD(5) * t125, qJD(2) * t14 + qJD(3) * t3 + qJD(4) * t4 - qJD(5) * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t559, t560, 0, 0, 0, 0, 0, 0, 0, 0, t564, t593, 0, 0, 0, 0, 0, 0, t569, 0, t744, qJD(3) * t175 + t747, 0, 0, 0, 0, 0, 0, t744, -t569, 0, qJD(3) * t113 + qJD(4) * t122 - qJD(5) * t343 + t747, 0, 0, 0, 0, 0, 0, t194 + t199 + t629, -t631 + t821 + t823 + t824, t833, qJD(3) * t16 + qJD(4) * t22 + qJD(5) * t829 + t695; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t535, t574, t372, t535, t561, 0, -qJD(3) * t355 - t399 * t562, qJD(3) * t354 + t399 * t563, 0, 0, t229, -t805, t184, -t794, -t182, 0, t468 - t626, t469 - t625 (-t350 * t722 - t422 * t759) * t711, -t645 + t175 * qJD(2) + (-t422 * t808 - t722 * t774) * t711, 0, -t184, t182, t229, -t805, -t794 (t350 * t402 - t669) * qJD(3) + t95 * qJD(4) + t566, -t468 + t688, -t469 + t687, t671 + t113 * qJD(2) + (-t398 * t808 + t402 * t774) * qJD(3) + t26 * qJD(4) + t110 * qJD(5) (t542 - t547) * t759 - t755, -0.2e1 * t486 + (t778 - t779) * qJD(3) + t119 * qJD(4) + t270, t423 * t616 + t817 + t825 (-t542 - t546) * t759 + t755, t199 + t495 + t826, t512, t684 + (-t423 * t477 - t832) * qJD(3) + t42 * qJD(4) - t215 * qJD(5) + t49 * qJD(6), t683 + t819 + (t421 * t477 - t831) * qJD(3) + t40 * qJD(4) + t204 * qJD(5) + t51 * qJD(6), -qJD(3) * t488 + t9 * qJD(4) + t708, t710 + t16 * qJD(2) + (t396 * t488 - t428 * t398) * qJD(3) + t2 * qJD(4) + t17 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t794, -t805, t184, -t794, -t466, 0, t545 - t551, t469 + t548, 0, 0, 0, -t184, t466, t794, -t805, -t794, t95 * qJD(3) + (-t783 - t720) * qJD(4) + t566, t551 + t686, -t469 + t685, t664 + t122 * qJD(2) + t26 * qJD(3) + (-qJ(5) * t808 - t795) * qJD(4) + t774 * qJD(5) (t540 - t547) * t759 - t754, t119 * qJD(3) + t270 + (-qJD(4) * t394 + t497) * t759, qJD(3) * t816 + t350 * t611 + t825 (-t540 - t546) * t759 + t754, qJD(3) * t207 - t541 + t826, -t513, t682 + t499 * qJD(2) + t42 * qJD(3) + (-t423 * t485 - t832) * qJD(4) + t816 * qJD(5) + t56 * qJD(6), -t428 * t611 + t681 + t819 + t40 * qJD(3) + t57 * qJD(6) + (qJD(4) * t485 - t566) * t421, t9 * qJD(3) - qJD(4) * t487 + t707 + t834, t709 + t22 * qJD(2) + t2 * qJD(3) + (-qJ(5) * t428 - t487 * t738) * qJD(4) + t19 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, -t794, t568, -t549 + t551, 0, 0, 0, 0, 0, 0, -qJD(3) * t215 + t544 + t817, qJD(3) * t204 + t423 * t568 - t541, t601, qJD(3) * t17 + qJD(4) * t19 - t693 + t834; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, -t790, qJD(3) * t806 - t509 * t780, -t96, t222 * qJD(3) + t350 * t475, -t223, qJD(2) * t222 + qJD(3) * t49 + qJD(4) * t56 + qJD(5) * t806 - qJD(6) * t71 + t691, qJD(2) * t806 + qJD(3) * t51 + qJD(4) * t57 + qJD(6) * t70 + t692, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t559, -t560, 0, 0, 0, 0, 0, 0, -t561, t372, -t564, -t593, 0, 0, 0, 0, 0, 0, t182, t184, -t744, -qJD(3) * t174 - t747, 0, 0, 0, 0, 0, 0, -t744, -t182, -t184, -qJD(3) * t112 + qJD(4) * t121 - qJD(5) * t240 - t747, 0, 0, 0, 0, 0, 0, -t200 + t495 - t629, -qJD(3) * t780 + qJD(4) * t215 + t631 - t820, -qJD(3) * t224 - qJD(4) * t115, -qJD(3) * t15 - qJD(4) * t21 - qJD(5) * t116 - t695; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t562, t563, 0, 0, 0, 0, 0, 0, 0, 0, t767, t769, 0, -t594, 0, 0, 0, 0, 0, 0, 0, -t767, -t769, -t605, 0, 0, 0, 0, 0, 0, -t533, -t579, -t578, -t675; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t575, t769, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t575, -t769, t602, 0, 0, 0, 0, 0, 0, t585, t581, -t604, -t672; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t575, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t603; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187, -t828, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t535, -t574, 0, -t535, 0, 0, t510 * t376, -t510 * t375, 0, 0, -t229, t805, 0, t794, -t569, 0, -t471 + t626, -t619 + t625, 0, qJD(2) * t174 + t645, 0, 0, t569, -t229, t805, t794, qJD(4) * t94, t471 - t688, t619 - t687, qJD(2) * t112 + qJD(4) * t25 - t671, t753 + t814, t270 + 0.2e1 * t486 - t833, -qJD(6) * t202 + t823 - t825, -t753 + t813, t194 - t200 - t826, -t512, qJD(4) * t43 + qJD(6) * t48 + t421 * t619 - t684 - t822, qJD(2) * t780 + qJD(4) * t41 + qJD(5) * t499 + qJD(6) * t50 - t683, qJD(2) * t224 - qJD(4) * t8 - t708, qJD(2) * t15 - qJD(4) * t1 + qJD(5) * t18 - t710; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t562, -t563, 0, 0, 0, 0, 0, 0, 0, 0, -t767, -t769, 0, t594, 0, 0, 0, 0, 0, 0, 0, t767, t769, t605, 0, 0, 0, 0, 0, 0, t533, t579, t578, t675; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t406, -t407, 0, 0, 0, 0, 0, 0, 0, 0, 0, t406, t557, -qJD(4) * t363 + t389, -t395, t382, 0, t395, 0, 0, t398 * t609 + t639, -t398 * t610 + t638, -t370, qJD(4) * t314 + t389; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t571, 0, -t748, -t407 - t405, 0, 0, 0, 0, t571, 0, 0, 0, t644, t748, t557 + t405 (-pkin(4) * t715 + t408) * qJD(4) + t389 + t476, t751, t382 - t633, t580, -t751, t743, 0, t359 * qJD(6) - t445 + t639, t358 * qJD(6) - t446 + t638, -t370 - t483 (-t373 * t738 + t408) * qJD(4) + t361 * qJD(5) + t484; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t415, qJD(4) * t398 + t558, 0, 0, 0, 0, 0, 0, -t494, t773, 0, qJD(4) * t361 - t756; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t173, t231, -t589 - t610, t173, t187, t770, qJD(4) * t359 - t396 * t610 - t454, qJD(4) * t358 - t396 * t609 - t455, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t794, t805, 0, t794, t570, 0, -t545 - t641, -t619 - t548, 0, 0, 0, 0, -t570, -t794, t805, t794, -qJD(3) * t94, t641 - t686, t619 - t685, -qJD(2) * t121 - qJD(3) * t25 - t664, t752 + t814, qJD(3) * t829 + t498 * t759 + t270, -t539 - t825 - t824, -t752 + t813, -t538 - t826 - t742, t513, -qJD(2) * t207 - qJD(3) * t43 + qJD(6) * t100 - t682 + t822, -qJD(2) * t215 - qJD(3) * t41 + qJD(6) * t101 - t681, qJD(2) * t115 + qJD(3) * t8 - t707, qJD(2) * t21 + qJD(3) * t1 + qJD(5) * t20 - t709; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t575, -t769, 0, 0, 0, 0, 0, 0, 0, 0, 0, t575, t769, -t602, 0, 0, 0, 0, 0, 0, -t585, -t581, t604, t672; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t571, 0, t749, t405, 0, 0, 0, 0, -t571, 0, 0, 0, -t644, -t749, qJD(5) - t405, t416 - t476, -t750, t382 + t633, -t580, t750, -t743, 0, -t357 * qJD(6) + t411 + t445, t356 * qJD(6) + t412 + t446, t483, -qJD(5) * t360 - t484; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), t416, -t395, t382, 0, t395, 0, 0, qJ(5) * t609 + t411, -qJ(5) * t610 + t412, 0, t416; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t415, t393, 0, 0, 0, 0, 0, 0, t493, t763, 0, t447; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t173, t231, -t475, t173, -t474, t770, t421 * t608 - t444, t423 * t608 - t443, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t794, -t568, t549 + t641, 0, 0, 0, 0, 0, 0, -t544 + t820 - t823 + t824, -t742 + (-qJD(6) * t759 - t568) * t423, -t601, qJD(2) * t116 - qJD(3) * t18 - qJD(4) * t20 + t693; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t575, 0, 0, 0, 0, 0, 0, 0, 0, 0, t603; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t415, -t634 - t558, 0, 0, 0, 0, 0, 0, t494, -t773, 0, qJD(4) * t360 + t756; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t415, -t393, 0, 0, 0, 0, 0, 0, -t493, -t763, 0, -t447; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t828, -t474, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, t790, t202 * qJD(3) + (t350 * t620 + t409) * t759, t96, t213 * qJD(3) + (-t350 * t621 + t611) * t759, -t223, qJD(2) * t213 - qJD(3) * t48 - qJD(4) * t100 - qJD(5) * t815 - t691, qJD(2) * t815 - qJD(3) * t50 - qJD(4) * t101 + t412 * t759 - t692, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198, t818, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t230, t589, -t173, t198, -t770, qJD(4) * t357 + t454, -qJD(4) * t356 + t455, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t230, t759 * t621, -t173, t532, -t770, t444, t443, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t818, t532, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t7;
