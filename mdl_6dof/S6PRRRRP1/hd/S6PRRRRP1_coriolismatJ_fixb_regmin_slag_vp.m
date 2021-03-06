% Calculate minimal parameter regressor of coriolis matrix for
% S6PRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x27]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6PRRRRP1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:59:21
% EndTime: 2019-03-08 23:59:44
% DurationCPUTime: 12.83s
% Computational Cost: add. (11054->593), mult. (24467->846), div. (0->0), fcn. (27167->10), ass. (0->479)
t481 = cos(qJ(5));
t478 = sin(qJ(5));
t479 = sin(qJ(3));
t482 = cos(qJ(3));
t477 = sin(pkin(6));
t480 = sin(qJ(2));
t683 = t477 * t480;
t718 = cos(pkin(6));
t522 = t479 * t683 - t482 * t718;
t729 = sin(qJ(4));
t502 = t729 * t522;
t414 = t479 * t718 + t482 * t683;
t730 = cos(qJ(4));
t621 = t730 * t414;
t493 = t621 - t502;
t681 = t478 * t493;
t483 = cos(qJ(2));
t682 = t477 * t483;
t216 = t481 * t682 + t681;
t709 = t216 * t478;
t194 = t709 / 0.2e1;
t676 = t481 * t493;
t217 = -t478 * t682 + t676;
t707 = t217 * t481;
t381 = t621 / 0.2e1;
t789 = t381 - t502 / 0.2e1;
t770 = t194 + t707 / 0.2e1 + t789;
t565 = t414 * t729 + t730 * t522;
t796 = (t493 - t707 - t709) * t565;
t799 = t796 * qJD(1);
t801 = t770 * qJD(6) + t799;
t771 = -t709 / 0.2e1 - t707 / 0.2e1 + t789;
t800 = -qJD(6) * t771 - t799;
t635 = qJD(3) + qJD(4);
t475 = t478 ^ 2;
t736 = -t475 / 0.2e1;
t759 = -pkin(9) - pkin(8);
t452 = t759 * t479;
t453 = t759 * t482;
t564 = -t730 * t452 - t453 * t729;
t791 = t564 * t481;
t608 = t791 / 0.2e1;
t618 = t729 * t479;
t619 = t730 * t482;
t436 = t618 - t619;
t302 = t481 * t436;
t464 = t729 * t482;
t465 = t730 * t479;
t438 = -t465 - t464;
t792 = t564 * t478;
t798 = -t438 * pkin(5) + qJ(6) * t302 + t792;
t512 = -t619 / 0.2e1 + t618 / 0.2e1;
t601 = t682 / 0.2e1;
t240 = t436 * t601 + t512 * t682;
t653 = qJD(4) * t565;
t659 = qJD(3) * t565;
t797 = t240 * qJD(2) + t653 + t659;
t690 = t438 * t478;
t263 = -pkin(5) * t690 + t564;
t470 = -pkin(3) * t482 - pkin(2);
t569 = t436 * pkin(4) + t438 * pkin(10);
t521 = t470 + t569;
t440 = t729 * t452;
t620 = t730 * t453;
t566 = -t620 + t440;
t776 = t566 * t478;
t173 = -t481 * t521 + t776;
t689 = t438 * t481;
t149 = qJ(6) * t689 - t173;
t726 = t436 * pkin(5);
t121 = t149 + t726;
t775 = t566 * t481;
t174 = t478 * t521 + t775;
t150 = qJ(6) * t690 + t174;
t731 = t481 / 0.2e1;
t734 = -t478 / 0.2e1;
t536 = t121 * t734 + t150 * t731;
t299 = t478 * t436;
t410 = pkin(5) * t299;
t540 = -t410 + t566;
t753 = t493 / 0.2e1;
t795 = t263 * t753 + (t540 / 0.2e1 - t536) * t565;
t780 = 0.2e1 * t438;
t774 = t635 * t438;
t562 = t436 * t774;
t793 = t438 * t601;
t647 = qJD(5) * t481;
t758 = -t121 / 0.2e1;
t594 = t149 / 0.2e1 + t758;
t523 = (-t726 / 0.2e1 + t594) * t481;
t767 = qJD(2) * t523;
t788 = -pkin(5) * t647 + t767;
t787 = qJD(1) * t771;
t784 = -t410 / 0.2e1 + t440 / 0.2e1;
t783 = (t173 - t776) * t438;
t782 = (t174 - t775) * t438;
t747 = -t436 / 0.2e1;
t239 = (t747 + t512) * t682;
t643 = t239 * qJD(1);
t781 = t564 * t635 - t643;
t779 = 0.2e1 * t478;
t361 = t436 * t682;
t678 = t478 * t361;
t286 = t481 * t683 + t678;
t752 = t286 / 0.2e1;
t631 = t730 * pkin(3);
t777 = t263 * t540;
t723 = t481 * pkin(5);
t469 = -pkin(4) - t723;
t704 = t493 * t469;
t330 = t635 * t436;
t476 = t481 ^ 2;
t458 = t476 + t475;
t773 = t635 * t458;
t772 = t121 - t149;
t431 = t438 ^ 2;
t634 = -t436 ^ 2 + t431;
t424 = t465 / 0.2e1 + t464 / 0.2e1;
t746 = -t438 / 0.2e1;
t238 = (t746 - t424) * t682;
t644 = t238 * qJD(1);
t663 = qJD(2) * t470;
t769 = t438 * t663 + t644;
t728 = pkin(3) * t479;
t288 = t436 * t728 - t438 * t470;
t768 = -qJD(2) * t288 + t644;
t430 = -t620 / 0.2e1;
t766 = t430 + t784;
t429 = t620 / 0.2e1;
t765 = t429 - t784;
t735 = t476 / 0.2e1;
t298 = (t736 + t735) * t438;
t661 = qJD(2) * t481;
t615 = t478 * t661;
t764 = t298 * t635 + t431 * t615;
t459 = t476 - t475;
t276 = t459 * t635 + t615 * t780;
t732 = -t481 / 0.2e1;
t733 = t478 / 0.2e1;
t535 = t216 * t732 + t217 * t733;
t518 = t535 * t436;
t624 = t478 * t683;
t674 = t481 * t361;
t287 = t624 - t674;
t534 = t286 * t733 + t287 * t732;
t55 = t518 + t534;
t673 = t55 * qJD(1);
t762 = qJD(5) * t523 - t673;
t595 = t302 / 0.2e1;
t537 = pkin(5) * t595 + t481 * t594;
t761 = t537 * qJD(5) + t673;
t724 = t438 * pkin(4);
t725 = t436 * pkin(10);
t332 = -t724 + t725;
t305 = t332 + t728;
t291 = t481 * t305;
t123 = t291 + t798;
t757 = t123 / 0.2e1;
t311 = t481 * t332;
t130 = t311 + t798;
t756 = t130 / 0.2e1;
t755 = -t216 / 0.2e1;
t754 = t217 / 0.2e1;
t630 = t729 * pkin(3);
t467 = t630 + pkin(10);
t672 = qJ(6) + t467;
t417 = t672 * t478;
t750 = -t417 / 0.2e1;
t418 = t672 * t481;
t749 = -t418 / 0.2e1;
t748 = t418 / 0.2e1;
t745 = t438 / 0.2e1;
t468 = -t631 - pkin(4);
t447 = t468 - t723;
t743 = -t447 / 0.2e1;
t742 = t447 / 0.2e1;
t722 = -qJ(6) - pkin(10);
t448 = t722 * t478;
t741 = -t448 / 0.2e1;
t449 = t722 * t481;
t740 = t449 / 0.2e1;
t739 = -t449 / 0.2e1;
t738 = t467 / 0.2e1;
t737 = -t469 / 0.2e1;
t727 = pkin(5) * t478;
t721 = pkin(5) * qJD(5);
t54 = t518 - t534;
t720 = t54 * qJD(2);
t719 = t55 * qJD(2);
t25 = t772 * t690;
t715 = qJD(2) * t25;
t714 = t123 * t478;
t713 = t130 * t478;
t625 = qJ(6) * t299;
t290 = t478 * t305;
t671 = -t290 + t791;
t158 = t625 - t671;
t712 = t158 * t481;
t310 = t478 * t332;
t670 = t791 - t310;
t159 = t625 - t670;
t711 = t159 * t481;
t710 = t216 * t436;
t708 = t217 * t436;
t550 = -t121 * t481 - t150 * t478;
t527 = t550 * t436;
t23 = (-t123 * t481 - t158 * t478) * t438 + t527;
t706 = t23 * qJD(2);
t24 = (-t130 * t481 - t159 * t478) * t438 + t527;
t705 = t24 * qJD(2);
t360 = t438 * t682;
t698 = t360 * t478;
t697 = t360 * t481;
t696 = t417 * t478;
t695 = t417 * t481;
t694 = t418 * t478;
t693 = t418 * t481;
t692 = t436 * t468;
t691 = t438 * t467;
t688 = t448 * t478;
t687 = t448 * t481;
t686 = t449 * t478;
t685 = t449 * t481;
t45 = -t216 * t286 + t217 * t287 - t360 * t565;
t684 = t45 * qJD(1);
t677 = t481 * t565;
t416 = t458 * t631;
t450 = t458 * qJD(6);
t669 = t416 * qJD(4) + t450;
t289 = -t436 * t470 - t438 * t728;
t666 = qJD(2) * t289;
t665 = qJD(2) * t436;
t664 = qJD(2) * t438;
t662 = qJD(2) * t480;
t660 = qJD(2) * t482;
t658 = qJD(3) * t478;
t657 = qJD(3) * t479;
t656 = qJD(3) * t481;
t655 = qJD(3) * t482;
t654 = qJD(3) * t483;
t652 = qJD(4) * t470;
t651 = qJD(4) * t478;
t650 = qJD(4) * t481;
t649 = qJD(5) * t217;
t648 = qJD(5) * t478;
t236 = t634 * t478;
t646 = t236 * qJD(2);
t237 = t634 * t481;
t645 = t237 * qJD(2);
t642 = t634 * qJD(2);
t641 = t298 * qJD(2);
t640 = t299 * qJD(2);
t296 = t302 * qJD(2);
t308 = t458 * t431;
t639 = t308 * qJD(2);
t309 = t459 * t431;
t638 = t309 * qJD(2);
t637 = t424 * qJD(2);
t460 = -t479 ^ 2 + t482 ^ 2;
t636 = t460 * qJD(2);
t633 = pkin(5) * t689;
t632 = t729 / 0.2e1;
t629 = pkin(2) * t479 * qJD(2);
t628 = pkin(2) * t660;
t627 = qJD(6) * t727;
t626 = t727 / 0.2e1;
t623 = t478 * t730;
t622 = t481 * t730;
t617 = t436 * t663;
t614 = t438 * t648;
t613 = t438 * t647;
t340 = t436 * t664;
t612 = t477 * t662;
t611 = qJD(2) * t682;
t463 = t478 * t647;
t610 = t479 * t660;
t609 = t263 * t734;
t260 = t565 * t733;
t261 = t565 * t731;
t607 = t698 / 0.2e1;
t606 = -t697 / 0.2e1;
t605 = t690 / 0.2e1;
t604 = -t689 / 0.2e1;
t603 = t689 / 0.2e1;
t602 = t683 / 0.2e1;
t600 = t565 * t734;
t599 = t678 / 0.2e1;
t598 = t299 / 0.2e1;
t597 = -t676 / 0.2e1;
t596 = t674 / 0.2e1;
t592 = t753 - t493 / 0.2e1;
t591 = -t290 / 0.2e1 + t608;
t590 = t310 / 0.2e1 - t791 / 0.2e1;
t587 = t730 * qJD(3);
t586 = t730 * qJD(4);
t585 = t729 * qJD(3);
t584 = t729 * qJD(4);
t582 = t635 * t481;
t581 = t635 * t478;
t580 = pkin(5) * t604;
t244 = pkin(5) * t260;
t579 = pkin(3) * t584;
t578 = pkin(3) * t585;
t577 = -t631 / 0.2e1;
t576 = t630 / 0.2e1;
t573 = t565 * t604;
t572 = -t623 / 0.2e1;
t571 = t622 / 0.2e1;
t382 = -t621 / 0.2e1;
t563 = t478 * t582;
t561 = t481 * t581;
t560 = pkin(3) * t572;
t17 = t121 * t123 + t150 * t158 + t777;
t486 = t123 * t755 + t158 * t754 + t795;
t504 = t287 * t749 - t360 * t743 + t417 * t752;
t2 = t486 + t504;
t559 = t2 * qJD(1) + t17 * qJD(2);
t20 = t121 * t130 + t150 * t159 + t777;
t487 = t130 * t755 + t159 * t754 + t795;
t503 = t286 * t741 + t287 * t740 - t360 * t737;
t4 = t487 + t503;
t558 = t4 * qJD(1) + t20 * qJD(2);
t26 = -t150 * t772 - t263 * t633;
t7 = -t594 * t217 + (t565 * t603 + t752) * pkin(5);
t557 = -qJD(1) * t7 + qJD(2) * t26;
t539 = (t216 / 0.2e1 - t681 / 0.2e1) * t438;
t37 = t606 + t539;
t48 = t311 * t436 + t783;
t556 = t37 * qJD(1) + t48 * qJD(2);
t513 = t493 * t746;
t495 = t216 * t745 + t478 * t513;
t39 = t606 + t495;
t46 = t291 * t436 + t783;
t555 = t39 * qJD(1) + t46 * qJD(2);
t538 = (t754 + t597) * t438;
t42 = t607 + t538;
t49 = t782 + (t670 - t791) * t436;
t554 = t42 * qJD(1) + t49 * qJD(2);
t494 = t217 * t745 + t481 * t513;
t44 = t607 + t494;
t47 = t782 + (t671 - t791) * t436;
t553 = t44 * qJD(1) + t47 * qJD(2);
t56 = t550 * t438;
t64 = -t438 * t535 - t793;
t552 = -qJD(1) * t64 - qJD(2) * t56;
t529 = t565 * t745 + t602;
t67 = t596 - t710 / 0.2e1 - t529 * t478;
t99 = t173 * t436 + t564 * t690;
t551 = qJD(1) * t67 - qJD(2) * t99;
t293 = t693 + t696;
t548 = t691 - t692;
t547 = t685 + t688;
t100 = -t174 * t436 - t564 * t689;
t66 = t599 + t708 / 0.2e1 + t529 * t481;
t546 = qJD(1) * t66 - qJD(2) * t100;
t143 = t592 * t481;
t542 = -t468 / 0.2e1 + t577;
t488 = (t738 - t630 / 0.2e1 - pkin(10) / 0.2e1) * t438 + (-pkin(4) / 0.2e1 + t542) * t436;
t58 = t478 * t488;
t545 = t143 * qJD(1) - t58 * qJD(2);
t544 = t438 * (qJD(5) + t665);
t272 = t382 + t381;
t342 = t429 + t430;
t543 = -t272 * qJD(1) - t342 * qJD(2);
t541 = t725 / 0.2e1 - t724 / 0.2e1;
t533 = t694 / 0.2e1 - t695 / 0.2e1;
t532 = -t693 / 0.2e1 - t696 / 0.2e1;
t531 = t685 / 0.2e1 + t688 / 0.2e1;
t530 = t436 * t738 + t468 * t745;
t528 = t481 * t544;
t247 = -qJD(5) * t424 + t340;
t526 = pkin(4) / 0.2e1 + t542;
t525 = t541 * t481;
t524 = t563 * t780;
t484 = t532 * t565 + (t194 * t730 + t217 * t571 + t565 * t632) * pkin(3) + t493 * t742;
t16 = -t704 / 0.2e1 - t531 * t565 + t484;
t193 = (t417 * t623 + t418 * t622 + t447 * t729) * pkin(3);
t485 = (t121 * t572 + t150 * t571 + t263 * t632) * pkin(3) + t130 * t750 + t159 * t748 + t540 * t742;
t505 = t123 * t741 + t158 * t740 + t540 * t737;
t6 = t485 + t505;
t520 = t16 * qJD(1) + t6 * qJD(2) + t193 * qJD(3);
t160 = t447 * t727;
t9 = -t594 * t418 + (t447 * t603 + t609 + t757) * pkin(5);
t519 = -qJD(2) * t9 + qJD(3) * t160;
t517 = t530 * t481;
t13 = (t159 / 0.2e1 - t158 / 0.2e1) * t481 + (-t130 / 0.2e1 + t757) * t478 + ((t750 + t741) * t481 + (t748 + t740) * t478) * t436;
t516 = -qJD(2) * t13 - qJD(3) * t416;
t492 = t438 * t533 + t536;
t32 = t492 + t765;
t515 = qJD(2) * t32 + qJD(3) * t293 - t787;
t11 = t594 * t449 + (t469 * t603 + t609 + t756) * pkin(5);
t222 = t469 * t727;
t89 = (t577 + t737 + t743) * t727;
t514 = -qJD(2) * t11 - qJD(3) * t89 + qJD(4) * t222;
t370 = t526 * t478;
t96 = -t311 / 0.2e1 - t525;
t510 = pkin(4) * t651 - qJD(2) * t96 + qJD(3) * t370;
t371 = t526 * t481;
t499 = t478 * t541 + t608;
t94 = t499 + t590;
t509 = pkin(4) * t650 - qJD(2) * t94 + qJD(3) * t371;
t508 = (-t585 - t584) * pkin(3);
t81 = -t291 / 0.2e1 - t517;
t507 = -qJD(2) * t81 - t468 * t658;
t497 = t478 * t530 + t608;
t79 = t497 - t591;
t506 = -qJD(2) * t79 - t468 * t656;
t167 = t576 + (t740 + t749) * t481 + (t448 / 0.2e1 + t750) * t478;
t491 = (-t686 / 0.2e1 + t687 / 0.2e1) * t438 + t536;
t34 = t491 + t765;
t501 = qJD(2) * t34 - qJD(3) * t167 - qJD(4) * t547 - t787;
t140 = t592 * t478;
t61 = t481 * t488;
t500 = -t140 * qJD(1) - t61 * qJD(2) - t478 * t578;
t490 = t691 / 0.2e1 - t692 / 0.2e1 + (t729 * t746 + t730 * t747) * pkin(3);
t474 = pkin(5) * t648;
t457 = t478 * t579;
t451 = t459 * qJD(5);
t373 = pkin(4) * t732 + t468 * t731 + t481 * t577;
t372 = pkin(4) * t734 + t468 * t733 + t560;
t338 = t697 / 0.2e1;
t337 = -t698 / 0.2e1;
t329 = (-t438 * t661 + t581) * pkin(5);
t306 = t635 * t424;
t294 = t298 * qJD(5);
t278 = 0.2e1 * t429 - t440;
t268 = t296 + t647;
t267 = -t640 - t648;
t241 = -t424 * t682 + t793;
t235 = t241 * qJD(2);
t230 = t239 * qJD(2);
t228 = t238 * qJD(2);
t225 = t263 * t626;
t219 = t561 - t641;
t218 = -t563 + t641;
t210 = t528 * t779;
t189 = 0.2e1 * t382 + t502;
t179 = -t340 * t476 - t294;
t168 = t576 - t531 - t532;
t152 = qJD(5) * t302 - t645;
t151 = -qJD(5) * t299 + t646;
t148 = t261 + t677 / 0.2e1;
t147 = 0.2e1 * t261;
t146 = t260 - t600;
t145 = 0.2e1 * t260;
t144 = t493 * t732 + t597;
t139 = t681 / 0.2e1 + t493 * t733;
t112 = -t294 + (t476 * t664 - t563) * t436;
t102 = -t478 * t774 + t645;
t101 = -t438 * t582 - t646;
t98 = (qJD(5) - t665) * t689 * t779 - t459 * t330;
t97 = t792 + t311 / 0.2e1 - t525;
t95 = t499 - t590;
t90 = pkin(5) * t560 + (t447 + t469) * t626;
t82 = t792 + t291 / 0.2e1 - t517;
t80 = t497 + t591;
t69 = -t708 / 0.2e1 + t573 + t599 + t481 * t602;
t68 = t710 / 0.2e1 + t565 * t605 + t596 - t624 / 0.2e1;
t65 = t216 * t604 + t217 * t605 - t793;
t62 = (-t476 / 0.2e1 + 0.2e1 * t736 - t735) * t565;
t60 = pkin(4) * t595 + pkin(10) * t603 + t490 * t481 + t776;
t59 = pkin(4) * t598 + pkin(10) * t605 + t490 * t478 - t775;
t43 = t337 + t494;
t41 = t337 + t538;
t40 = t338 + t495;
t38 = t338 + t539;
t33 = t491 + t766;
t31 = t492 + t766;
t30 = 0.2e1 * t244;
t28 = -pkin(5) * t600 + t244;
t15 = -t677 * t739 - t448 * t600 + t704 / 0.2e1 + t484;
t14 = t711 / 0.2e1 - t713 / 0.2e1 - t449 * t598 + t712 / 0.2e1 + t448 * t595 - t714 / 0.2e1 + t533 * t436;
t12 = pkin(5) * t756 + t121 * t740 + t149 * t739 + t469 * t580 + t225;
t10 = pkin(5) * t757 + t121 * t749 + t149 * t748 + t447 * t580 + t225;
t8 = t149 * t754 + t217 * t758 + (t573 + t752) * pkin(5);
t5 = t485 - t505;
t3 = t487 - t503;
t1 = t486 - t504;
t18 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t45 + t635 * t796; 0, 0, -t612, -t611, 0, 0, 0, 0, 0 (-t479 * t654 - t480 * t660) * t477 (t479 * t662 - t482 * t654) * t477, 0, 0, 0, 0, 0, t241 * t635 + t436 * t612, t240 * t635 - t438 * t612, 0, 0, 0, 0, 0 (t286 * t436 + t360 * t690) * qJD(2) + t40 * qJD(3) + t38 * qJD(4) + t69 * qJD(5) (-t287 * t436 + t360 * t689) * qJD(2) + t43 * qJD(3) + t41 * qJD(4) + t68 * qJD(5), t635 * t54 + (t286 * t481 + t287 * t478) * t664, t684 + (t121 * t286 + t150 * t287 - t263 * t360) * qJD(2) + t1 * qJD(3) + t3 * qJD(4) + t8 * qJD(5) + t65 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t414 * qJD(3) - t479 * t611, qJD(3) * t522 - t482 * t611, 0, 0, 0, 0, 0, -qJD(3) * t493 + qJD(4) * t189 + t235, t797, 0, 0, 0, 0, 0, qJD(2) * t40 + qJD(4) * t144 + qJD(5) * t146 - t493 * t656, qJD(2) * t43 + qJD(4) * t139 + qJD(5) * t148 + t493 * t658, t62 * qJD(4) - t458 * t659 + t720, t1 * qJD(2) + (-t293 * t565 + t447 * t493) * qJD(3) + t15 * qJD(4) + t28 * qJD(5) + t801; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3) * t189 - qJD(4) * t493 + t235, t797, 0, 0, 0, 0, 0, qJD(2) * t38 + qJD(3) * t144 + qJD(5) * t145 - t493 * t650, qJD(2) * t41 + qJD(3) * t139 + qJD(5) * t147 + t493 * t651, t62 * qJD(3) - t458 * t653 + t720, t3 * qJD(2) + t15 * qJD(3) + (t547 * t565 + t704) * qJD(4) + t30 * qJD(5) + t801; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t69 + qJD(3) * t146 + qJD(4) * t145 - t649, qJD(2) * t68 + qJD(3) * t148 + qJD(4) * t147 + qJD(5) * t216, 0, -pkin(5) * t649 + qJD(2) * t8 + qJD(3) * t28 + qJD(4) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t65 + t635 * t770; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t635 * t238, -t635 * t239, 0, 0, 0, 0, 0, qJD(3) * t39 + qJD(4) * t37 - qJD(5) * t66, qJD(3) * t44 + qJD(4) * t42 - qJD(5) * t67, t55 * t635, qJD(3) * t2 + qJD(4) * t4 - qJD(5) * t7 - qJD(6) * t64 - t684; 0, 0, 0, 0, t479 * t655, t460 * qJD(3), 0, 0, 0, -pkin(2) * t657, -pkin(2) * t655, t562, -t635 * t634, 0, 0, 0, qJD(3) * t288 - t438 * t652, qJD(3) * t289 - t436 * t652, -t431 * t463 + t476 * t562, -t309 * qJD(5) - t436 * t524, t237 * t635 + t436 * t614, -t236 * t635 + t436 * t613, -t562, qJD(3) * t46 + qJD(4) * t48 + qJD(5) * t100, qJD(3) * t47 + qJD(4) * t49 + qJD(5) * t99, -qJD(3) * t23 - qJD(4) * t24 - qJD(5) * t25 + qJD(6) * t308, qJD(3) * t17 + qJD(4) * t20 + qJD(5) * t26 - qJD(6) * t56; 0, 0, 0, 0, t610, t636, t655, -t657, 0, -pkin(8) * t655 - t629, pkin(8) * t657 - t628, t340, -t642, -t330, t774, 0, -qJD(3) * t566 + qJD(4) * t278 - t768, t666 + t781, t112, t98, t102, t101, -t247 (t478 * t548 - t775) * qJD(3) + t59 * qJD(4) + t82 * qJD(5) + t555 (t481 * t548 + t776) * qJD(3) + t60 * qJD(4) + t80 * qJD(5) + t553, -t706 + (-t714 + t712 + (t694 - t695) * t436) * qJD(3) + t14 * qJD(4) + t761 (-t123 * t417 + t158 * t418 + t447 * t540) * qJD(3) + t5 * qJD(4) + t10 * qJD(5) + t31 * qJD(6) + t559; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t340, -t642, -t330, t774, 0, qJD(3) * t278 - qJD(4) * t566 - t769, -t617 + t781, t112, t98, t102, t101, -t247, t59 * qJD(3) + (t478 * t569 - t775) * qJD(4) + t97 * qJD(5) + t556, t60 * qJD(3) + (t481 * t569 + t776) * qJD(4) + t95 * qJD(5) + t554, -t705 + t14 * qJD(3) + (-t713 + t711 + (-t686 + t687) * t436) * qJD(4) + t761, t5 * qJD(3) + (t130 * t448 - t159 * t449 + t469 * t540) * qJD(4) + t12 * qJD(5) + t33 * qJD(6) + t558; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t764, t561 * t780 - t638, t478 * t544, t528, t306, qJD(3) * t82 + qJD(4) * t97 - qJD(5) * t174 - t546, qJD(3) * t80 + qJD(4) * t95 + qJD(5) * t173 - t551, -pkin(5) * t614 + t537 * t635 - t715, qJD(3) * t10 + qJD(4) * t12 - t150 * t721 + t557; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t639, qJD(3) * t31 + qJD(4) * t33 + t552; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4) * t272 + t228, t230, 0, 0, 0, 0, 0, -qJD(2) * t39 - qJD(4) * t143, -qJD(2) * t44 + qJD(4) * t140, -t719, -qJD(2) * t2 + qJD(4) * t16 + t800; 0, 0, 0, 0, -t610, -t636, 0, 0, 0, t629, t628, -t340, t642, 0, 0, 0, qJD(4) * t342 + t768, t643 - t666, t179, t210, t152, t151, t247, qJD(4) * t58 + qJD(5) * t81 - t555, qJD(4) * t61 + qJD(5) * t79 - t553, qJD(4) * t13 + t706 + t762, qJD(4) * t6 - qJD(5) * t9 + qJD(6) * t32 - t559; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t579, -pkin(3) * t586, t463, t451, 0, 0, 0, t468 * t648 - t481 * t579, t468 * t647 + t457, t669, qJD(4) * t193 + qJD(5) * t160 + qJD(6) * t293; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t508 - t543 (-t587 - t586) * pkin(3), t463, t451, 0, 0, 0, t372 * qJD(5) + t481 * t508 - t545, t373 * qJD(5) + t457 - t500, -t516 + t669 (-t448 * t623 - t449 * t622 + t469 * t729) * pkin(3) * qJD(4) + t90 * qJD(5) + t168 * qJD(6) + t520; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t219, t276, t268, t267, -t637, qJD(4) * t372 - t467 * t647 - t507, qJD(4) * t373 + t467 * t648 - t506, t788, qJD(4) * t90 - t418 * t721 + t519; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t773, qJD(4) * t168 + t515; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3) * t272 + t228, t230, 0, 0, 0, 0, 0, -qJD(2) * t37 + qJD(3) * t143, -qJD(2) * t42 - qJD(3) * t140, -t719, -qJD(2) * t4 - qJD(3) * t16 + t800; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t340, t642, 0, 0, 0, -qJD(3) * t342 + t769, t643 + t617, t179, t210, t152, t151, t247, -qJD(3) * t58 + qJD(5) * t96 - t556, -qJD(3) * t61 + qJD(5) * t94 - t554, -qJD(3) * t13 + t705 + t762, -qJD(3) * t6 - qJD(5) * t11 + qJD(6) * t34 - t558; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t578 + t543, pkin(3) * t587, t463, t451, 0, 0, 0, -t370 * qJD(5) + t481 * t578 + t545, -t371 * qJD(5) + t500, t450 + t516, -qJD(5) * t89 - qJD(6) * t167 - t520; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t463, t451, 0, 0, 0, -pkin(4) * t648, -pkin(4) * t647, t450, qJD(5) * t222 - qJD(6) * t547; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t219, t276, t268, t267, -t637, -pkin(10) * t647 - t510, pkin(10) * t648 - t509, t788, t449 * t721 + t514; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t773, t501; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t66, qJD(2) * t67, 0, qJD(2) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t764, -t524 + t638, -t302 * t635 - t340 * t478, t299 * t635 - t340 * t481, t306, -qJD(3) * t81 - qJD(4) * t96 + t546, -qJD(3) * t79 - qJD(4) * t94 + t551, -t523 * t635 + t715, qJD(3) * t9 + qJD(4) * t11 + qJD(6) * t633 - t557; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t218, -t276, -t296, t640, t637, qJD(4) * t370 + t507, qJD(4) * t371 + t506, -t767, qJD(4) * t89 - t519 - t627; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t218, -t276, -t296, t640, t637, t510, t509, -t767, -t514 - t627; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t64 + t635 * t771; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t639, -pkin(5) * t613 - qJD(3) * t32 - qJD(4) * t34 - t552; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t773, qJD(4) * t167 + t474 - t515; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t773, t474 - t501; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t18;
