% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRRP11_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP11_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP11_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_invdynJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:53:49
% EndTime: 2019-03-10 02:55:11
% DurationCPUTime: 46.30s
% Computational Cost: add. (55142->1168), mult. (158005->1541), div. (0->0), fcn. (132833->14), ass. (0->500)
t481 = cos(qJ(2));
t727 = cos(pkin(6));
t660 = pkin(1) * t727;
t465 = t481 * t660;
t452 = qJD(1) * t465;
t478 = sin(qJ(2));
t474 = sin(pkin(6));
t726 = cos(pkin(7));
t575 = t474 * (-pkin(10) * t726 - pkin(9));
t555 = t478 * t575;
t342 = qJD(1) * t555 + t452;
t464 = t478 * t660;
t513 = t481 * t575 - t464;
t343 = t513 * qJD(1);
t725 = sin(pkin(7));
t638 = t481 * t725;
t552 = pkin(2) * t478 - pkin(10) * t638;
t682 = qJD(1) * t474;
t381 = t552 * t682;
t761 = sin(qJ(3));
t614 = t725 * t761;
t763 = cos(qJ(3));
t617 = t726 * t763;
t408 = pkin(2) * t617 - pkin(10) * t614;
t616 = t726 * t761;
t692 = t408 * qJD(3) - t763 * t342 - t343 * t616 - t381 * t614;
t534 = -t478 * t616 + t481 * t763;
t377 = t534 * t474;
t361 = qJD(1) * t377;
t615 = t725 * t763;
t579 = qJD(3) * t615;
t550 = t579 - t361;
t260 = -t343 * t725 + t726 * t381;
t535 = t478 * t617 + t481 * t761;
t376 = t535 * t474;
t360 = qJD(1) * t376;
t804 = -t360 * pkin(3) + t361 * pkin(11) - t260 + (pkin(3) * t614 - pkin(11) * t615) * qJD(3);
t613 = t725 * t682;
t584 = t478 * t613;
t803 = pkin(11) * t584 - t692;
t517 = t535 * qJD(2);
t587 = t481 * t616;
t533 = t478 * t763 + t587;
t490 = qJD(3) * t533 + t517;
t658 = t761 * t478;
t621 = qJDD(1) * t658;
t488 = (qJD(1) * t490 + t621) * t474;
t636 = t727 * qJD(1);
t596 = t636 + qJD(2);
t542 = t596 * t725;
t527 = t761 * t542;
t631 = t727 * qJDD(1);
t582 = t631 + qJDD(2);
t541 = t582 * t725;
t588 = t481 * t617;
t570 = t474 * t588;
t688 = -qJDD(1) * t570 - t763 * t541;
t509 = qJD(3) * t527 + t688;
t200 = t488 + t509;
t802 = -qJDD(4) - t200;
t477 = sin(qJ(4));
t480 = cos(qJ(4));
t402 = t477 * t614 - t480 * t726;
t776 = -qJD(4) * t402 - t477 * t584 + t480 * t550;
t403 = t477 * t726 + t480 * t614;
t693 = qJD(4) * t403 + t477 * t550 + t480 * t584;
t410 = pkin(2) * t616 + pkin(10) * t615;
t691 = t410 * qJD(3) - t761 * t342 + t343 * t617 + t381 * t615;
t578 = qJD(3) * t614;
t549 = t578 - t360;
t764 = cos(qJ(1));
t619 = t727 * t764;
t762 = sin(qJ(1));
t404 = t478 * t619 + t481 * t762;
t591 = t474 * t614;
t551 = t762 * t478 - t481 * t619;
t780 = t551 * t726;
t280 = -t404 * t763 + t764 * t591 + t761 * t780;
t643 = t474 * t726;
t788 = t551 * t725 - t764 * t643;
t217 = t280 * t480 - t477 * t788;
t592 = t474 * t615;
t277 = t404 * t761 + t592 * t764 + t763 * t780;
t476 = sin(qJ(5));
t479 = cos(qJ(5));
t801 = t217 * t476 + t277 * t479;
t800 = t217 * t479 - t277 * t476;
t387 = pkin(11) * t726 + t410;
t659 = t725 * pkin(2);
t388 = -pkin(3) * t615 - pkin(11) * t614 - t659;
t678 = qJD(4) * t480;
t680 = qJD(4) * t477;
t778 = -t387 * t680 + t388 * t678 + t477 * t804 - t803 * t480;
t775 = pkin(3) * t584 + t691;
t627 = t474 * t658;
t687 = -qJD(1) * t570 - t763 * t542;
t309 = qJD(1) * t627 + t687;
t537 = qJD(4) + t309;
t797 = -pkin(12) * t549 - t778;
t796 = t693 * pkin(4) - t776 * pkin(12) + t775;
t216 = t280 * t477 + t480 * t788;
t702 = t474 * t481;
t685 = pkin(9) * t702 + t464;
t391 = t685 * qJD(1);
t622 = t481 * t643;
t303 = t391 + (qJD(1) * t622 + t542) * pkin(10);
t520 = pkin(2) * t727 + t555;
t308 = qJD(2) * pkin(2) + qJD(1) * t520 + t452;
t641 = t478 * t725;
t547 = -pkin(2) * t481 - pkin(10) * t641 - pkin(1);
t356 = t547 * t682;
t173 = t303 * t763 + t308 * t616 + t356 * t614;
t793 = -t173 + t537 * (pkin(4) * t477 - pkin(12) * t480);
t777 = -t387 * t678 - t388 * t680 + t803 * t477 + t480 * t804;
t519 = t533 * t474;
t311 = qJD(1) * t519 + t527;
t700 = t476 * t480;
t220 = -t309 * t700 - t479 * t311;
t792 = -t476 * t678 + t220;
t348 = t476 * t403 + t479 * t615;
t695 = qJD(5) * t348 - t476 * t549 - t479 * t776;
t590 = t476 * t615;
t676 = qJD(5) * t479;
t694 = -qJD(5) * t590 + t403 * t676 + t476 * t776 - t479 * t549;
t386 = -pkin(3) * t726 - t408;
t273 = t402 * pkin(4) - t403 * pkin(12) + t386;
t287 = t480 * t387 + t477 * t388;
t276 = -pkin(12) * t615 + t287;
t677 = qJD(5) * t476;
t746 = t273 * t676 - t276 * t677 + t476 * t796 - t479 * t797;
t186 = t476 * t273 + t479 * t276;
t745 = -qJD(5) * t186 + t476 * t797 + t479 * t796;
t696 = -pkin(4) * t549 - t777;
t790 = t477 * t676 - t792;
t618 = t727 * t762;
t536 = t478 * t764 + t481 * t618;
t789 = t536 * t725 + t762 * t643;
t623 = t474 * t638;
t518 = -qJDD(1) * t623 + t582 * t726 + qJDD(3);
t624 = t474 * t641;
t583 = qJD(2) * t624;
t500 = qJD(1) * t583 + t518;
t492 = t311 * t680 - t477 * t500;
t670 = qJDD(1) * t478;
t649 = t474 * t670;
t671 = qJD(1) * qJD(2);
t650 = t481 * t671;
t652 = t761 * qJD(3);
t657 = t478 * t682;
t199 = (qJD(2) * t616 + t652) * t657 - t761 * t541 - t763 * t649 + (-qJDD(1) * t587 - t650 * t763) * t474 + t687 * qJD(3);
t672 = t481 * t613 - qJD(3);
t521 = -t596 * t726 + t672;
t508 = qJD(4) * t521 + t199;
t767 = t480 * t508 + t492;
t787 = -qJD(5) * t537 + t767;
t405 = -t478 * t618 + t481 * t764;
t512 = t536 * t726;
t282 = t405 * t763 - t512 * t761 + t591 * t762;
t218 = t282 * t477 - t480 * t789;
t601 = t727 * t725;
t553 = t761 * t601;
t337 = t553 + t519;
t602 = t727 * t726;
t526 = t623 - t602;
t271 = t337 * t477 + t480 * t526;
t561 = g(1) * t218 - g(2) * t216 + g(3) * t271;
t471 = t474 ^ 2;
t786 = 0.2e1 * t471;
t349 = t403 * t479 - t590;
t747 = pkin(5) * t693 + qJ(6) * t695 - qJD(6) * t349 + t745;
t785 = -qJ(6) * t694 - qJD(6) * t348 + t746;
t735 = pkin(5) * t694 + t696;
t172 = -t761 * t303 + t308 * t617 + t356 * t615;
t226 = pkin(3) * t311 + pkin(11) * t309;
t132 = t480 * t172 + t477 * t226;
t115 = pkin(12) * t311 + t132;
t667 = pkin(11) * t680;
t784 = t793 * t479 + (t115 + t667) * t476;
t442 = -pkin(4) * t480 - pkin(12) * t477 - pkin(3);
t783 = -t479 * t115 + t442 * t676 + t476 * t793;
t246 = t480 * t311 - t477 * t521;
t70 = t246 * t677 + t476 * t802 + t479 * t787;
t187 = t246 * t476 - t479 * t537;
t244 = t311 * t477 + t480 * t521;
t243 = qJD(5) + t244;
t722 = t187 * t243;
t782 = -t70 - t722;
t71 = t246 * t676 - t476 * t787 + t479 * t802;
t189 = t479 * t246 + t476 * t537;
t719 = t189 * t243;
t781 = t71 + t719;
t341 = t465 + t520;
t686 = pkin(2) * t702 + pkin(10) * t624;
t370 = -pkin(1) * t474 - t686;
t254 = -t341 * t725 + t726 * t370;
t554 = t763 * t601;
t336 = t627 - t570 - t554;
t334 = t336 * pkin(3);
t168 = -t337 * pkin(11) + t254 + t334;
t325 = (t622 + t601) * pkin(10) + t685;
t203 = t763 * t325 + t341 * t616 + t370 * t614;
t175 = -pkin(11) * t526 + t203;
t110 = t477 * t168 + t480 * t175;
t107 = pkin(12) * t336 + t110;
t202 = -t761 * t325 + t341 * t617 + t370 * t615;
t174 = pkin(3) * t526 - t202;
t272 = t337 * t480 - t477 * t526;
t119 = t271 * pkin(4) - t272 * pkin(12) + t174;
t60 = t479 * t107 + t476 * t119;
t779 = t132 + t667;
t698 = t479 * t480;
t466 = pkin(11) * t698;
t384 = t476 * t442 + t466;
t506 = -qJDD(4) - t509;
t774 = t488 - t506;
t773 = (qJDD(2) + 0.2e1 * t631) * t474;
t708 = t309 * t480;
t772 = -t678 - t708;
t395 = t685 * qJD(2);
t770 = t71 * pkin(5) + qJDD(6);
t769 = qJD(3) * t311 + t474 * (qJD(1) * t517 + t621) + qJDD(4) + t688;
t219 = t282 * t480 + t477 * t789;
t281 = t405 * t761 + t512 * t763 - t592 * t762;
t152 = -t219 * t476 + t281 * t479;
t210 = t272 * t476 - t336 * t479;
t768 = -g(1) * t152 - g(2) * t801 + g(3) * t210;
t766 = t189 ^ 2;
t482 = qJD(1) ^ 2;
t635 = -t477 * t199 - t480 * t500;
t673 = t246 * qJD(4);
t124 = t635 + t673;
t123 = qJDD(5) + t124;
t760 = pkin(5) * t123;
t759 = pkin(9) * t474;
t758 = pkin(11) * t477;
t757 = pkin(11) * t480;
t750 = t187 * pkin(5);
t475 = -qJ(6) - pkin(12);
t154 = pkin(3) * t521 - t172;
t104 = t244 * pkin(4) - t246 * pkin(12) + t154;
t238 = -t308 * t725 + t726 * t356;
t146 = t309 * pkin(3) - t311 * pkin(11) + t238;
t155 = -pkin(11) * t521 + t173;
t97 = t477 * t146 + t480 * t155;
t89 = pkin(12) * t537 + t97;
t46 = t479 * t104 - t476 * t89;
t38 = -qJ(6) * t189 + t46;
t32 = pkin(5) * t243 + t38;
t749 = -t38 + t32;
t744 = qJ(6) * t70;
t743 = qJ(6) * t71;
t629 = pkin(9) * t657;
t605 = qJD(2) * t636;
t593 = pkin(1) * t605;
t625 = pkin(1) * t631;
t669 = qJDD(1) * t481;
t648 = t474 * t669;
t662 = pkin(9) * t648 + t478 * t625 + t481 * t593;
t306 = -qJD(2) * t629 + t662;
t606 = t726 * t671;
t634 = qJDD(1) * t726;
t233 = (t541 + (-t478 * t606 + t481 * t634) * t474) * pkin(10) + t306;
t451 = t481 * t625;
t544 = -t478 * t593 + t451;
t567 = -t650 - t670;
t545 = t567 * pkin(9);
t239 = t582 * pkin(2) + ((-t478 * t634 - t481 * t606) * pkin(10) + t545) * t474 + t544;
t540 = qJD(2) * t552;
t283 = (qJD(1) * t540 + qJDD(1) * t547) * t474;
t581 = qJD(3) * t617;
t91 = t763 * t233 + t239 * t616 + t283 * t614 - t303 * t652 + t308 * t581 + t356 * t579;
t82 = pkin(11) * t500 + t91;
t162 = -t239 * t725 + t726 * t283;
t90 = t200 * pkin(3) + t199 * pkin(11) + t162;
t630 = t146 * t680 + t155 * t678 + t477 * t82 - t480 * t90;
t22 = -pkin(4) * t774 + t630;
t12 = t22 + t770;
t742 = t12 * t476;
t741 = t22 * t476;
t47 = t104 * t476 + t479 * t89;
t39 = -qJ(6) * t187 + t47;
t740 = t243 * t39;
t739 = t243 * t46;
t738 = t243 * t47;
t737 = t476 * t70;
t736 = t479 * t71;
t221 = -t309 * t698 + t311 * t476;
t675 = qJD(6) * t479;
t709 = t309 * t477;
t734 = pkin(5) * t709 + qJ(6) * t221 - t477 * t675 + (pkin(5) * t477 - qJ(6) * t698) * qJD(4) + (-t466 + (qJ(6) * t477 - t442) * t476) * qJD(5) + t784;
t699 = t477 * t479;
t733 = qJ(6) * t220 + (-qJD(4) * pkin(11) - qJ(6) * qJD(5)) * t699 + (-qJD(6) * t477 + (-pkin(11) * qJD(5) - qJ(6) * qJD(4)) * t480) * t476 + t783;
t679 = qJD(4) * t479;
t732 = (-t477 * t679 - t480 * t677) * pkin(11) + t783;
t731 = -qJD(5) * t384 + t784;
t131 = -t477 * t172 + t226 * t480;
t114 = -pkin(4) * t311 - t131;
t666 = pkin(11) * t678;
t730 = pkin(5) * t790 - t114 + t666;
t646 = qJD(5) * t475;
t159 = pkin(4) * t246 + pkin(12) * t244;
t96 = t480 * t146 - t477 * t155;
t69 = t476 * t159 + t479 * t96;
t716 = t244 * t476;
t729 = -qJ(6) * t716 + t476 * t646 + t675 - t69;
t68 = t479 * t159 - t476 * t96;
t715 = t244 * t479;
t728 = -pkin(5) * t246 - qJ(6) * t715 - qJD(6) * t476 + t479 * t646 - t68;
t724 = t123 * t476;
t723 = t123 * t479;
t721 = t187 * t476;
t720 = t189 * t187;
t718 = t189 * t476;
t717 = t243 * t246;
t714 = t246 * t244;
t713 = t277 * t477;
t712 = t277 * t480;
t711 = t281 * t477;
t710 = t281 * t480;
t707 = t311 * t309;
t706 = t336 * t477;
t705 = t336 * t480;
t704 = t471 * t482;
t703 = t474 * t478;
t701 = t476 * t477;
t697 = t480 * t124;
t684 = t764 * pkin(1) + t762 * t759;
t472 = t478 ^ 2;
t473 = t481 ^ 2;
t683 = t472 - t473;
t681 = qJD(2) * t474;
t674 = t154 * qJD(4);
t665 = t481 * t704;
t664 = t561 * t476;
t663 = t377 * pkin(3) + t686;
t661 = pkin(5) * t476 + pkin(11);
t656 = t478 * t681;
t654 = qJD(4) + t687;
t653 = t763 * qJD(3);
t651 = pkin(1) * t786;
t645 = t404 * t725;
t644 = t405 * t725;
t642 = t477 * t725;
t639 = t480 * t725;
t637 = qJD(6) + t750;
t571 = -t146 * t678 + t155 * t680 - t477 * t90 - t480 * t82;
t21 = pkin(12) * t774 - t571;
t580 = qJD(3) * t616;
t574 = t761 * t233 - t239 * t617 - t283 * t615 + t303 * t653 + t308 * t580 + t356 * t578;
t83 = -pkin(3) * t500 + t574;
t35 = t124 * pkin(4) + pkin(12) * t767 + t83;
t3 = t104 * t676 + t479 * t21 + t476 * t35 - t89 * t677;
t59 = -t107 * t476 + t479 * t119;
t109 = t168 * t480 - t477 * t175;
t185 = t479 * t273 - t276 * t476;
t286 = -t477 * t387 + t388 * t480;
t633 = t477 * t537;
t632 = t243 * t479;
t628 = t478 * t665;
t626 = t478 * t650;
t620 = -pkin(1) * t762 + t764 * t759;
t612 = t474 * t482 * t727;
t610 = g(1) * t801 - g(2) * t152;
t153 = t219 * t479 + t281 * t476;
t609 = -g(1) * t800 - g(2) * t153;
t608 = g(1) * t216 + g(2) * t218;
t607 = -g(1) * t277 + g(2) * t281;
t275 = pkin(4) * t615 - t286;
t603 = t479 * t678 - t221;
t600 = -t46 * t479 - t47 * t476;
t345 = t513 * qJD(2);
t382 = t474 * t540;
t261 = -t345 * t725 + t726 * t382;
t211 = t272 * t479 + t336 * t476;
t469 = pkin(5) * t479 + pkin(4);
t599 = -t469 * t480 + t475 * t477;
t595 = 0.2e1 * t636 + qJD(2);
t594 = pkin(11) * t376 + t663;
t453 = qJD(2) * t465;
t344 = qJD(2) * t555 + t453;
t134 = -t325 * t652 + t341 * t581 + t763 * t344 + t345 * t616 + t370 * t579 + t382 * t614;
t127 = pkin(11) * t583 + t134;
t258 = qJD(3) * t553 + t474 * t490;
t259 = qJD(3) * t554 + ((t588 - t658) * qJD(3) + t534 * qJD(2)) * t474;
t140 = t258 * pkin(3) - t259 * pkin(11) + t261;
t45 = -t477 * t127 + t140 * t480 - t168 * t680 - t175 * t678;
t577 = -t551 * pkin(2) + pkin(10) * t645;
t576 = -t536 * pkin(2) + pkin(10) * t644;
t106 = -pkin(4) * t336 - t109;
t88 = -pkin(4) * t537 - t96;
t572 = -pkin(12) * t123 + t243 * t88;
t44 = t480 * t127 + t477 * t140 + t168 * t678 - t175 * t680;
t41 = pkin(12) * t258 + t44;
t135 = -t325 * t653 - t341 * t580 - t761 * t344 + t345 * t617 - t370 * t578 + t382 * t615;
t128 = -pkin(3) * t583 - t135;
t160 = qJD(4) * t272 + t259 * t477 - t480 * t583;
t161 = -qJD(4) * t271 + t259 * t480 + t477 * t583;
t64 = t160 * pkin(4) - t161 * pkin(12) + t128;
t9 = -t107 * t677 + t119 * t676 + t479 * t41 + t476 * t64;
t300 = -t404 * t616 - t551 * t763;
t569 = t300 * pkin(3) + t577;
t302 = -t405 * t616 - t536 * t763;
t568 = t302 * pkin(3) + t576;
t566 = g(1) * t764 + g(2) * t762;
t565 = -g(1) * (t281 * t700 + t282 * t479) - g(2) * (t277 * t700 - t280 * t479) - g(3) * (t336 * t700 + t337 * t479);
t564 = -g(1) * (-t281 * t698 + t282 * t476) - g(2) * (-t277 * t698 - t280 * t476) - g(3) * (-t336 * t698 + t337 * t476);
t248 = t300 * t480 + t404 * t642;
t250 = t302 * t480 + t405 * t642;
t299 = t404 * t617 - t551 * t761;
t301 = t405 * t617 - t536 * t761;
t313 = t377 * t480 + t477 * t624;
t563 = -g(1) * (-t250 * t476 + t301 * t479) - g(2) * (-t248 * t476 + t299 * t479) - g(3) * (-t313 * t476 + t376 * t479);
t562 = -g(1) * (t250 * t479 + t301 * t476) - g(2) * (t248 * t479 + t299 * t476) - g(3) * (t313 * t479 + t376 * t476);
t560 = -g(1) * t219 + g(2) * t217 - g(3) * t272;
t247 = t300 * t477 - t404 * t639;
t249 = t302 * t477 - t405 * t639;
t312 = t377 * t477 - t480 * t624;
t559 = -g(1) * t249 - g(2) * t247 - g(3) * t312;
t558 = g(1) * t281 + g(2) * t277 + g(3) * t336;
t557 = g(1) * t282 - g(2) * t280 + g(3) * t337;
t556 = g(1) * t301 + g(2) * t299 + g(3) * t376;
t546 = -t22 + t561;
t42 = -pkin(4) * t258 - t45;
t539 = t299 * pkin(11) + t569;
t538 = t301 * pkin(11) + t568;
t529 = qJD(2) * t617 + t653;
t525 = g(1) * t153 - g(2) * t800 + g(3) * t211 - t3;
t4 = -qJD(5) * t47 - t21 * t476 + t479 * t35;
t524 = t561 - t630;
t523 = qJD(4) * t537;
t10 = -qJD(5) * t60 - t41 * t476 + t479 * t64;
t516 = t481 * (qJD(2) * t761 + t580);
t515 = t477 * t630 - t480 * t571 - t557;
t510 = t521 * t725;
t507 = qJD(3) * t510;
t503 = t4 + t768;
t502 = -t404 * pkin(2) - pkin(10) * t788 + t620;
t501 = t405 * pkin(2) + pkin(10) * t789 + t684;
t499 = t280 * pkin(3) + t502;
t498 = t282 * pkin(3) + t501;
t495 = -pkin(11) * t277 + t499;
t494 = t281 * pkin(11) + t498;
t493 = t500 * t725;
t491 = t309 * t537 + t523;
t445 = t475 * t479;
t444 = t475 * t476;
t430 = t661 * t477;
t424 = t479 * t442;
t409 = -pkin(9) * t703 + t465;
t394 = -pkin(9) * t656 + t453;
t389 = t452 - t629;
t383 = -pkin(11) * t700 + t424;
t351 = -qJ(6) * t701 + t384;
t335 = -qJ(6) * t699 + t424 + (-pkin(11) * t476 - pkin(5)) * t480;
t307 = t474 * t545 + t544;
t269 = t281 * pkin(3);
t267 = t277 * pkin(3);
t222 = pkin(5) * t348 + t275;
t184 = t187 ^ 2;
t156 = -qJ(6) * t348 + t186;
t142 = pkin(5) * t402 - qJ(6) * t349 + t185;
t101 = t258 * t476 - t272 * t677 + (qJD(5) * t336 + t161) * t479;
t100 = qJD(5) * t211 + t161 * t476 - t258 * t479;
t87 = -t184 + t766;
t81 = pkin(5) * t210 + t106;
t80 = -pkin(5) * t716 + t97;
t76 = -t123 * t480 + t243 * t633;
t74 = t123 * t402 + t243 * t693;
t65 = t637 + t88;
t63 = t123 * t271 + t160 * t243;
t52 = -qJ(6) * t210 + t60;
t51 = -t189 * t246 + t243 * t632 + t724;
t50 = -t243 ^ 2 * t476 + t187 * t246 + t723;
t49 = -t71 + t719;
t48 = -t70 + t722;
t43 = pkin(5) * t271 - qJ(6) * t211 + t59;
t37 = t243 * t721 - t736;
t36 = t189 * t632 - t737;
t31 = t187 * t694 + t348 * t71;
t30 = -t189 * t695 - t349 * t70;
t29 = t187 * t790 + t71 * t701;
t28 = -t70 * t699 + (-t477 * t677 + t603) * t189;
t27 = t100 * t187 + t210 * t71;
t26 = t101 * t189 - t211 * t70;
t25 = pkin(5) * t100 + t42;
t19 = t480 * t71 + t792 * t243 + (-t187 * t537 - t243 * t676 - t724) * t477;
t18 = t480 * t70 + t603 * t243 + (t189 * t537 - t243 * t677 + t723) * t477;
t17 = -t123 * t348 - t187 * t693 - t243 * t694 - t402 * t71;
t16 = t123 * t349 + t189 * t693 - t243 * t695 - t402 * t70;
t15 = t101 * t243 + t123 * t211 + t160 * t189 - t271 * t70;
t14 = -t100 * t243 - t123 * t210 - t160 * t187 - t271 * t71;
t13 = -t476 * t781 + t479 * t782;
t11 = t187 * t695 - t189 * t694 + t348 * t70 - t349 * t71;
t8 = t187 * t221 + t189 * t220 + (-t187 * t479 - t718) * t678 + (t737 - t736 + (-t189 * t479 + t721) * qJD(5)) * t477;
t7 = -t100 * t189 - t101 * t187 + t210 * t70 - t211 * t71;
t6 = -qJ(6) * t100 - qJD(6) * t210 + t9;
t5 = pkin(5) * t160 - qJ(6) * t101 - qJD(6) * t211 + t10;
t2 = -qJD(6) * t187 + t3 - t743;
t1 = -qJD(6) * t189 + t4 + t744 + t760;
t20 = [0, 0, 0, 0, 0, qJDD(1), g(1) * t762 - g(2) * t764, t566, 0, 0 (qJDD(1) * t472 + 0.2e1 * t626) * t471 (t478 * t669 - t671 * t683) * t786, t481 * t595 * t681 + t478 * t773 (qJDD(1) * t473 - 0.2e1 * t626) * t471, t481 * t773 - t595 * t656, t582 * t727, -t395 * t596 + t409 * t582 + t307 * t727 + g(1) * t404 - g(2) * t405 + (-t478 * t671 + t669) * t651, -g(1) * t551 + g(2) * t536 - t306 * t727 - t394 * t596 + t567 * t651 - t582 * t685 ((-t389 * qJD(2) + qJDD(1) * t685 + t306 + (-qJD(2) * t409 + t394) * qJD(1)) * t481 + (-t391 * qJD(2) - qJDD(1) * t409 - t307) * t478 - t566) * t474, t471 * qJDD(1) * pkin(1) ^ 2 - g(1) * t620 - g(2) * t684 + t306 * t685 + t307 * t409 - t389 * t395 + t391 * t394, -t199 * t337 + t259 * t311, t199 * t336 - t200 * t337 - t258 * t311 - t259 * t309, t199 * t526 - t259 * t521 + t311 * t583 + t337 * t500, t200 * t336 + t258 * t309, t200 * t526 + t258 * t521 - t309 * t583 - t336 * t500, -t500 * t526 - t510 * t656, -g(1) * t280 - g(2) * t282 - t135 * t521 + t162 * t336 + t172 * t583 + t254 * t200 + t202 * t500 + t238 * t258 + t261 * t309 + t526 * t574, t134 * t521 + t162 * t337 - t173 * t583 - t254 * t199 - t203 * t500 + t238 * t259 + t261 * t311 + t526 * t91 + t607, g(1) * t788 - g(2) * t789 - t134 * t309 - t135 * t311 - t172 * t259 - t173 * t258 + t202 * t199 - t203 * t200 - t91 * t336 + t337 * t574, -g(1) * t502 - g(2) * t501 + t173 * t134 + t172 * t135 + t162 * t254 - t202 * t574 + t91 * t203 + t238 * t261, t246 * t161 - t272 * t767, -t272 * t124 - t246 * t160 - t161 * t244 + t271 * t767, t161 * t537 + t246 * t258 - t272 * t802 - t336 * t767, t124 * t271 + t160 * t244, -t160 * t654 + t271 * t506 - t124 * t336 - t244 * t258 + (-t271 * t621 + (-t271 * t516 + (-t160 * t761 - t271 * t529) * t478) * qJD(1)) * t474, -t506 * t336 + t654 * t258 + (t336 * t621 + (t336 * t516 + (t258 * t761 + t336 * t529) * t478) * qJD(1)) * t474, -g(1) * t217 - g(2) * t219 - t109 * t802 + t174 * t124 + t128 * t244 + t154 * t160 + t96 * t258 + t83 * t271 - t336 * t630 + t45 * t537, t110 * t802 + t128 * t246 + t154 * t161 - t174 * t767 - t97 * t258 + t83 * t272 + t336 * t571 - t44 * t537 + t608, t109 * t767 - t110 * t124 - t97 * t160 - t96 * t161 - t44 * t244 - t45 * t246 + t271 * t571 + t272 * t630 - t607, -g(1) * t495 - g(2) * t494 - t109 * t630 - t110 * t571 + t154 * t128 + t83 * t174 + t97 * t44 + t96 * t45, t26, t7, t15, t27, t14, t63, t10 * t243 + t100 * t88 + t106 * t71 + t123 * t59 + t160 * t46 + t187 * t42 + t210 * t22 + t271 * t4 + t609, t101 * t88 - t106 * t70 - t123 * t60 - t160 * t47 + t189 * t42 + t211 * t22 - t243 * t9 - t271 * t3 + t610, -t10 * t189 - t100 * t47 - t101 * t46 - t187 * t9 - t210 * t3 - t211 * t4 + t59 * t70 - t60 * t71 - t608, t3 * t60 + t47 * t9 + t4 * t59 + t46 * t10 + t22 * t106 + t88 * t42 - g(1) * (t217 * pkin(4) + t216 * pkin(12) + t495) - g(2) * (t219 * pkin(4) + t218 * pkin(12) + t494) t26, t7, t15, t27, t14, t63, t1 * t271 + t100 * t65 + t12 * t210 + t123 * t43 + t160 * t32 + t187 * t25 + t243 * t5 + t71 * t81 + t609, t101 * t65 + t12 * t211 - t123 * t52 - t160 * t39 + t189 * t25 - t2 * t271 - t243 * t6 - t70 * t81 + t610, -t1 * t211 - t100 * t39 - t101 * t32 - t187 * t6 - t189 * t5 - t2 * t210 + t43 * t70 - t52 * t71 - t608, t2 * t52 + t39 * t6 + t1 * t43 + t32 * t5 + t12 * t81 + t65 * t25 - g(1) * (-t216 * t475 + t217 * t469 - t277 * t661 + t499) - g(2) * (-t218 * t475 + t219 * t469 + t281 * t661 + t498); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t628, t683 * t704, -t481 * t612 + t649, t628, t478 * t612 + t648, t582, t451 + t391 * t596 + g(1) * t536 + g(2) * t551 + (-t605 + t704) * t478 * pkin(1) + (-g(3) * t481 + t545) * t474, pkin(1) * t665 + t389 * t596 + g(1) * t405 + g(2) * t404 + (pkin(9) * t671 + g(3)) * t703 - t662, 0, 0, -t199 * t614 + t311 * t550, -t199 * t615 - t200 * t614 + t361 * t309 + t311 * t360 + (-t309 * t615 - t311 * t614) * qJD(3), -t199 * t726 - t311 * t584 + t361 * t521 + t493 * t761 - t507 * t763, -t200 * t615 + t309 * t549, -t200 * t726 + t309 * t584 - t360 * t521 + t493 * t763 + t507 * t761, t518 * t726 - (qJD(1) * t602 - t672) * t584, -g(1) * t302 - g(2) * t300 - g(3) * t377 - t162 * t615 - t172 * t584 - t200 * t659 + t549 * t238 - t260 * t309 + t408 * t500 + t521 * t691 - t574 * t726, t162 * t614 + t173 * t584 + t199 * t659 + t550 * t238 - t260 * t311 - t410 * t500 + t521 * t692 - t726 * t91 + t556, -g(3) * t624 + t91 * t615 + t574 * t614 - g(1) * t644 - g(2) * t645 + t172 * t361 + t173 * t360 + t408 * t199 - t410 * t200 + t691 * t311 - t692 * t309 + (-t172 * t615 - t173 * t614) * qJD(3), -g(1) * t576 - g(2) * t577 - g(3) * t686 - t162 * t659 - t172 * t691 + t173 * t692 - t238 * t260 - t408 * t574 + t91 * t410, t246 * t776 - t403 * t767, -t403 * t124 - t244 * t776 - t246 * t693 + t402 * t767, t246 * t549 - t403 * t802 + t537 * t776 + t615 * t767, t124 * t402 + t244 * t693, t402 * t506 + t124 * t615 - t549 * t244 + (-t402 * t621 + (-t402 * t516 + (-t402 * t529 - t693 * t761) * t478) * qJD(1)) * t474 - t693 * t654, t537 * t549 + t615 * t802, -g(1) * t250 - g(2) * t248 - g(3) * t313 + t386 * t124 + t154 * t693 + t244 * t775 - t286 * t802 + t83 * t402 + t537 * t777 + t549 * t96 + t615 * t630, t154 * t776 + t246 * t775 + t287 * t802 - t386 * t767 + t83 * t403 - t537 * t778 - t549 * t97 - t571 * t615 - t559, -t287 * t124 - t244 * t778 - t246 * t777 + t286 * t767 + t402 * t571 + t403 * t630 - t693 * t97 - t776 * t96 - t556, -g(1) * t538 - g(2) * t539 - g(3) * t594 + t154 * t775 - t286 * t630 - t287 * t571 + t83 * t386 + t777 * t96 + t778 * t97, t30, t11, t16, t31, t17, t74, t123 * t185 + t187 * t696 + t22 * t348 + t243 * t745 + t275 * t71 + t4 * t402 + t46 * t693 + t694 * t88 + t562, -t123 * t186 + t189 * t696 + t22 * t349 - t243 * t746 - t275 * t70 - t3 * t402 - t47 * t693 - t695 * t88 + t563, t185 * t70 - t186 * t71 - t187 * t746 - t189 * t745 - t3 * t348 - t349 * t4 + t46 * t695 - t47 * t694 + t559, t3 * t186 + t4 * t185 + t22 * t275 - g(1) * (t250 * pkin(4) + t249 * pkin(12) + t538) - g(2) * (t248 * pkin(4) + t247 * pkin(12) + t539) - g(3) * (pkin(4) * t313 + pkin(12) * t312 + t594) + t696 * t88 + t746 * t47 + t745 * t46, t30, t11, t16, t31, t17, t74, t1 * t402 + t12 * t348 + t123 * t142 + t187 * t735 + t222 * t71 + t243 * t747 + t32 * t693 + t65 * t694 + t562, t12 * t349 - t123 * t156 + t189 * t735 - t2 * t402 - t222 * t70 - t243 * t785 - t39 * t693 - t65 * t695 + t563, -t1 * t349 + t142 * t70 - t156 * t71 - t187 * t785 - t189 * t747 - t2 * t348 + t32 * t695 - t39 * t694 + t559, t2 * t156 + t1 * t142 + t12 * t222 - g(1) * (-t249 * t475 + t250 * t469 + t301 * t661 + t568) - g(2) * (-t247 * t475 + t248 * t469 + t299 * t661 + t569) - g(3) * (-t312 * t475 + t313 * t469 + t376 * t661 + t663) + t735 * t65 + t785 * t39 + t747 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t707, -t309 ^ 2 + t311 ^ 2, -t309 * t521 - t199, -t707, -t311 * t521 - t200, t500, -t173 * t521 - t238 * t311 + t558 - t574, -t172 * t521 + t238 * t309 + t557 - t91, 0, 0, -t492 * t477 + (t246 * t537 - t477 * t508) * t480, -t246 * t709 - t480 * t767 + (-t124 - t673) * t477 + t772 * t244, -t246 * t311 + t477 * t769 + t480 * t491, t244 * t633 - t697, t244 * t311 - t477 * t491 + t480 * t769, -t537 * t311, -pkin(3) * t124 + g(1) * t710 + g(2) * t712 + g(3) * t705 - t131 * t537 + t154 * t709 - t173 * t244 - t96 * t311 + t477 * t674 - t83 * t480 - t523 * t757 + t758 * t802, pkin(3) * t767 - g(1) * t711 - g(2) * t713 - g(3) * t706 + t154 * t708 - t173 * t246 + t97 * t311 + t83 * t477 + t480 * t674 + t537 * t779 + t757 * t802, -pkin(11) * t697 - t767 * t758 + t515 + (-t680 - t709) * t97 + t772 * t96 + (t131 + t666) * t246 + t779 * t244, -t83 * pkin(3) + g(1) * t269 + g(2) * t267 + g(3) * t334 - t96 * t131 - t97 * t132 - t154 * t173 + ((-t477 * t97 - t480 * t96) * qJD(4) + t515) * pkin(11), t28, t8, t18, t29, t19, t76, -t114 * t187 + t123 * t383 - t220 * t88 + t731 * t243 + (-t4 + (pkin(11) * t187 + t476 * t88) * qJD(4)) * t480 + (pkin(11) * t71 + t46 * t537 + t676 * t88 + t741) * t477 + t564, -t114 * t189 - t123 * t384 - t221 * t88 - t732 * t243 + (t3 + (pkin(11) * t189 + t479 * t88) * qJD(4)) * t480 + (-pkin(11) * t70 + t22 * t479 - t47 * t537 - t677 * t88) * t477 + t565, t220 * t47 + t221 * t46 + t383 * t70 - t384 * t71 - t731 * t189 - t732 * t187 + t600 * t678 + (-t3 * t476 - t4 * t479 + (t46 * t476 - t47 * t479) * qJD(5) + t558) * t477, t3 * t384 + t4 * t383 - t88 * t114 - g(1) * (-pkin(4) * t710 - pkin(12) * t711 - t269) - g(2) * (-pkin(4) * t712 - pkin(12) * t713 - t267) - g(3) * (-pkin(4) * t705 - pkin(12) * t706 - t334) + t732 * t47 + t731 * t46 + (t22 * t477 + t678 * t88 - t557) * pkin(11), t28, t8, t18, t29, t19, t76, t123 * t335 - t220 * t65 + t430 * t71 + (qJD(4) * t476 * t65 - t1) * t480 + t734 * t243 + t730 * t187 + (t32 * t537 + t65 * t676 + t742) * t477 + t564, -t123 * t351 - t221 * t65 - t430 * t70 + (t65 * t679 + t2) * t480 - t733 * t243 + t730 * t189 + (t12 * t479 - t39 * t537 - t65 * t677) * t477 + t565, t220 * t39 + t221 * t32 + t335 * t70 - t351 * t71 - t734 * t189 - t733 * t187 + (-t32 * t479 - t39 * t476) * t678 + (-t1 * t479 - t2 * t476 + (t32 * t476 - t39 * t479) * qJD(5) + t558) * t477, t2 * t351 + t1 * t335 + t12 * t430 - g(1) * (t281 * t599 + t282 * t661 - t269) - g(2) * (t277 * t599 - t280 * t661 - t267) - g(3) * (t336 * t599 + t337 * t661 - t334) + t730 * t65 + t733 * t39 + t734 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t714, -t244 ^ 2 + t246 ^ 2, t244 * t537 - t767, -t714, t246 * t309 - t635, -t802, -t154 * t246 + t537 * t97 + t524, t154 * t244 + t537 * t96 - t560 + t571, 0, 0, t36, t13, t51, t37, t50, -t717, -pkin(4) * t71 - t187 * t97 - t243 * t68 - t246 * t46 + t572 * t476 + (-pkin(12) * qJD(5) * t243 + t546) * t479, pkin(4) * t70 - t189 * t97 + t741 + t246 * t47 + (pkin(12) * t677 + t69) * t243 + t572 * t479 - t664, t187 * t69 + t189 * t68 + (t3 - t739 + (qJD(5) * t189 - t71) * pkin(12)) * t479 + (-t4 - t738 + (qJD(5) * t187 - t70) * pkin(12)) * t476 + t560, -t46 * t68 - t47 * t69 - t88 * t97 + t546 * pkin(4) + (qJD(5) * t600 + t3 * t479 - t4 * t476 + t560) * pkin(12), t36, t13, t51, t37, t50, -t717, t123 * t444 - t187 * t80 - t246 * t32 - t469 * t71 + t728 * t243 + (t244 * t65 + (t65 + t750) * qJD(5)) * t476 + (-t12 + t561) * t479, t65 * t715 + t742 + t123 * t445 - t189 * t80 + t246 * t39 + t469 * t70 - t729 * t243 + (pkin(5) * t718 + t479 * t65) * qJD(5) - t664, t444 * t70 + t445 * t71 - t728 * t189 - t729 * t187 + (-t243 * t32 + t2) * t479 + (-t1 - t740) * t476 + t560, -t2 * t445 + t1 * t444 - t12 * t469 - g(1) * (-t218 * t469 - t219 * t475) - g(2) * (t216 * t469 + t217 * t475) - g(3) * (-t271 * t469 - t272 * t475) + (pkin(5) * t677 - t80) * t65 + t729 * t39 + t728 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t720, t87, t48, -t720, t49, t123, -t189 * t88 + t503 + t738, t187 * t88 + t525 + t739, 0, 0, t720, t87, t48, -t720, t49, t123, 0.2e1 * t760 + t744 + t740 + (-t637 - t65) * t189 + t503, -pkin(5) * t766 + t743 + t243 * t38 + (qJD(6) + t65) * t187 + t525, pkin(5) * t70 - t187 * t749, t749 * t39 + (-t65 * t189 + t1 + t768) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t781, t782, -t184 - t766, pkin(4) * t802 + t39 * t187 + t32 * t189 - t524 + t770;];
tau_reg  = t20;
