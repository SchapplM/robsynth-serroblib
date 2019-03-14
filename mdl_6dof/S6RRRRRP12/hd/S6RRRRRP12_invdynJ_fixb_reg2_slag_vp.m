% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RRRRRP12
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
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRRRRP12_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP12_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP12_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_invdynJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:20:48
% EndTime: 2019-03-10 03:22:06
% DurationCPUTime: 42.30s
% Computational Cost: add. (55289->1166), mult. (158418->1524), div. (0->0), fcn. (133157->14), ass. (0->500)
t473 = cos(qJ(2));
t734 = cos(pkin(6));
t674 = pkin(1) * t734;
t459 = t473 * t674;
t448 = qJD(1) * t459;
t470 = sin(qJ(2));
t467 = sin(pkin(6));
t733 = cos(pkin(7));
t572 = t467 * (-pkin(10) * t733 - pkin(9));
t552 = t470 * t572;
t340 = qJD(1) * t552 + t448;
t458 = t470 * t674;
t506 = t473 * t572 - t458;
t341 = t506 * qJD(1);
t732 = sin(pkin(7));
t652 = t473 * t732;
t549 = pkin(2) * t470 - pkin(10) * t652;
t693 = qJD(1) * t467;
t381 = t549 * t693;
t759 = sin(qJ(3));
t627 = t732 * t759;
t761 = cos(qJ(3));
t630 = t733 * t761;
t408 = pkin(2) * t630 - pkin(10) * t627;
t629 = t733 * t759;
t701 = t408 * qJD(3) - t761 * t340 - t341 * t629 - t381 * t627;
t529 = -t470 * t629 + t473 * t761;
t377 = t529 * t467;
t362 = qJD(1) * t377;
t628 = t732 * t761;
t577 = qJD(3) * t628;
t548 = t577 - t362;
t530 = t470 * t630 + t473 * t759;
t511 = t530 * qJD(2);
t585 = t473 * t629;
t528 = t470 * t761 + t585;
t484 = qJD(3) * t528 + t511;
t672 = t759 * t470;
t633 = qJDD(1) * t672;
t482 = (qJD(1) * t484 + t633) * t467;
t650 = t734 * qJD(1);
t593 = t650 + qJD(2);
t539 = t593 * t732;
t522 = t759 * t539;
t644 = t734 * qJDD(1);
t581 = t644 + qJDD(2);
t538 = t581 * t732;
t586 = t473 * t630;
t565 = t467 * t586;
t699 = -qJDD(1) * t565 - t761 * t538;
t499 = qJD(3) * t522 + t699;
t479 = -t499 - t482;
t805 = -qJDD(4) + t479;
t242 = -t341 * t732 + t733 * t381;
t376 = t530 * t467;
t361 = qJD(1) * t376;
t804 = -t361 * pkin(3) + t362 * pkin(11) - t242 + (pkin(3) * t627 - pkin(11) * t628) * qJD(3);
t626 = t732 * t693;
t583 = t470 * t626;
t803 = pkin(11) * t583 - t701;
t469 = sin(qJ(4));
t472 = cos(qJ(4));
t402 = t469 * t627 - t472 * t733;
t778 = qJD(4) * t402 + t469 * t583 - t472 * t548;
t403 = t469 * t733 + t472 * t627;
t702 = qJD(4) * t403 + t469 * t548 + t472 * t583;
t410 = pkin(2) * t629 + pkin(10) * t628;
t700 = t410 * qJD(3) - t759 * t340 + t341 * t630 + t381 * t628;
t576 = qJD(3) * t627;
t547 = t576 - t361;
t640 = t467 * t672;
t698 = -qJD(1) * t565 - t761 * t539;
t306 = qJD(1) * t640 + t698;
t532 = qJD(4) + t306;
t754 = pkin(11) * t472;
t709 = t467 * t473;
t696 = pkin(9) * t709 + t458;
t392 = t696 * qJD(1);
t656 = t467 * t733;
t634 = t473 * t656;
t299 = t392 + (qJD(1) * t634 + t539) * pkin(10);
t514 = pkin(2) * t734 + t552;
t305 = qJD(2) * pkin(2) + qJD(1) * t514 + t448;
t654 = t470 * t732;
t544 = -pkin(2) * t473 - pkin(10) * t654 - pkin(1);
t354 = t544 * t693;
t155 = -t759 * t299 + t305 * t630 + t354 * t628;
t513 = t528 * t467;
t308 = qJD(1) * t513 + t522;
t209 = pkin(3) * t308 + pkin(11) * t306;
t118 = t472 * t155 + t469 * t209;
t691 = qJD(4) * t469;
t783 = pkin(11) * t691 + t118;
t387 = pkin(11) * t733 + t410;
t673 = t732 * pkin(2);
t388 = -pkin(3) * t628 - pkin(11) * t627 - t673;
t690 = qJD(4) * t472;
t782 = -t387 * t691 + t388 * t690 + t469 * t804 - t803 * t472;
t777 = pkin(3) * t583 + t700;
t474 = cos(qJ(1));
t651 = t474 * t734;
t760 = sin(qJ(1));
t404 = t470 * t760 - t473 * t651;
t405 = t470 * t651 + t473 * t760;
t588 = t467 * t627;
t270 = -t404 * t629 + t405 * t761 - t474 * t588;
t344 = t404 * t732 - t474 * t656;
t200 = t270 * t472 + t344 * t469;
t589 = t467 * t628;
t269 = t404 * t630 + t405 * t759 + t474 * t589;
t468 = sin(qJ(5));
t471 = cos(qJ(5));
t131 = t200 * t468 - t269 * t471;
t132 = t200 * t471 + t269 * t468;
t801 = -t547 * pkin(12) - t782;
t800 = t702 * pkin(4) + t778 * pkin(12) + t777;
t799 = pkin(12) * t308 + t783;
t156 = t299 * t761 + t305 * t629 + t354 * t627;
t797 = -qJD(5) * t754 - t156 + t532 * (pkin(4) * t469 - pkin(12) * t472);
t781 = -t387 * t690 - t388 * t691 + t803 * t469 + t472 * t804;
t347 = t468 * t403 + t471 * t628;
t704 = qJD(5) * t347 - t468 * t547 + t471 * t778;
t587 = t468 * t628;
t688 = qJD(5) * t471;
t703 = -qJD(5) * t587 + t403 * t688 - t468 * t778 - t471 * t547;
t631 = t734 * t760;
t535 = t474 * t470 + t473 * t631;
t794 = t535 * t732 + t760 * t656;
t719 = t306 * t469;
t793 = t691 + t719;
t636 = t467 * t652;
t512 = -qJDD(1) * t636 + t581 * t733 + qJDD(3);
t637 = t467 * t654;
t582 = qJD(2) * t637;
t491 = qJD(1) * t582 + t512;
t487 = t308 * t691 - t469 * t491;
t682 = qJDD(1) * t470;
t663 = t467 * t682;
t683 = qJD(1) * qJD(2);
t664 = t473 * t683;
t666 = t759 * qJD(3);
t671 = t470 * t693;
t182 = (qJD(2) * t629 + t666) * t671 - t759 * t538 - t761 * t663 + (-qJDD(1) * t585 - t664 * t761) * t467 + t698 * qJD(3);
t684 = t473 * t626 - qJD(3);
t518 = -t593 * t733 + t684;
t498 = qJD(4) * t518 + t182;
t768 = t472 * t498 + t487;
t792 = -qJD(5) * t532 + t768;
t227 = t308 * t469 + t472 * t518;
t224 = qJD(5) + t227;
t649 = t469 * t182 + t472 * t491;
t229 = t472 * t308 - t469 * t518;
t685 = t229 * qJD(4);
t111 = -t649 + t685;
t110 = qJDD(5) + t111;
t791 = -t270 * t469 + t344 * t472;
t169 = t229 * t468 - t471 * t532;
t171 = t471 * t229 + t468 * t532;
t708 = t468 * t472;
t205 = -t306 * t708 - t471 * t308;
t707 = t471 * t472;
t206 = -t306 * t707 + t308 * t468;
t59 = t229 * t688 - t468 * t792 + t471 * t805;
t743 = t471 * t59;
t689 = qJD(5) * t468;
t58 = t229 * t689 + t468 * t805 + t471 * t792;
t744 = t468 * t58;
t790 = t469 * ((t169 * t468 - t171 * t471) * qJD(5) - t743 + t744) - (t169 * t471 + t171 * t468) * t690 + t169 * t206 + t171 * t205;
t617 = t468 * t690 - t205;
t729 = t110 * t468;
t789 = t617 * t224 + t469 * (t169 * t532 + t224 * t688 + t729) - t472 * t59;
t646 = t171 * t224;
t727 = t169 * t224;
t788 = t468 * (t59 + t646) + t471 * (t58 + t727);
t464 = t467 ^ 2;
t787 = 0.2e1 * t464;
t386 = -pkin(3) * t733 - t408;
t265 = t402 * pkin(4) - t403 * pkin(12) + t386;
t281 = t472 * t387 + t469 * t388;
t268 = -pkin(12) * t628 + t281;
t749 = t265 * t688 - t268 * t689 + t468 * t800 - t471 * t801;
t741 = -pkin(4) * t547 - t781;
t440 = -pkin(4) * t472 - pkin(12) * t469 - pkin(3);
t737 = t440 * t688 + t468 * t797 - t471 * t799;
t613 = t734 * t732;
t325 = (t634 + t613) * pkin(10) + t696;
t339 = t459 + t514;
t697 = pkin(2) * t709 + pkin(10) * t637;
t370 = -pkin(1) * t467 - t697;
t185 = -t759 * t325 + t339 * t630 + t370 * t628;
t614 = t734 * t733;
t521 = t636 - t614;
t157 = pkin(3) * t521 - t185;
t550 = t759 * t613;
t335 = t550 + t513;
t263 = t335 * t469 + t472 * t521;
t258 = t263 * pkin(4);
t264 = t335 * t472 - t469 * t521;
t106 = -t264 * pkin(12) + t157 + t258;
t551 = t761 * t613;
t334 = t640 - t565 - t551;
t237 = -t339 * t732 + t733 * t370;
t333 = t334 * pkin(3);
t660 = -t335 * pkin(11) + t333;
t151 = t237 + t660;
t186 = t761 * t325 + t339 * t629 + t370 * t627;
t158 = -pkin(11) * t521 + t186;
t97 = t469 * t151 + t472 * t158;
t93 = pkin(12) * t334 + t97;
t784 = t468 * t106 + t471 * t93;
t780 = t468 * t265 + t471 * t268;
t715 = t334 * t472;
t716 = t334 * t469;
t779 = -pkin(4) * t715 - pkin(12) * t716;
t495 = -qJDD(4) - t499;
t776 = t482 - t495;
t775 = (qJDD(2) + 0.2e1 * t644) * t467;
t718 = t306 * t472;
t774 = -t690 - t718;
t611 = pkin(5) * t471 + qJ(6) * t468;
t396 = t696 * qJD(2);
t772 = qJD(3) * t308 + t467 * (qJD(1) * t511 + t633) + qJDD(4) + t699;
t220 = -t305 * t732 + t733 * t354;
t129 = t306 * pkin(3) - t308 * pkin(11) + t220;
t138 = -pkin(11) * t518 + t156;
t82 = t472 * t129 - t469 * t138;
t75 = -pkin(4) * t532 - t82;
t38 = t169 * pkin(5) - t171 * qJ(6) + t75;
t753 = pkin(12) * t110;
t771 = t224 * t38 - t753;
t240 = qJD(3) * t550 + t467 * t484;
t449 = qJD(2) * t459;
t342 = qJD(2) * t552 + t449;
t343 = t506 * qJD(2);
t536 = qJD(2) * t549;
t382 = t467 * t536;
t579 = qJD(3) * t630;
t120 = -t325 * t666 + t339 * t579 + t761 * t342 + t343 * t629 + t370 * t577 + t382 * t627;
t114 = pkin(11) * t582 + t120;
t241 = qJD(3) * t551 + ((t586 - t672) * qJD(3) + t529 * qJD(2)) * t467;
t243 = -t343 * t732 + t733 * t382;
t125 = t240 * pkin(3) - t241 * pkin(11) + t243;
t32 = t472 * t114 + t469 * t125 + t151 * t690 - t158 * t691;
t30 = pkin(12) * t240 + t32;
t578 = qJD(3) * t629;
t667 = t761 * qJD(3);
t121 = -t325 * t667 - t339 * t578 - t759 * t342 + t343 * t630 - t370 * t576 + t382 * t628;
t115 = -pkin(3) * t582 - t121;
t143 = qJD(4) * t264 + t241 * t469 - t472 * t582;
t144 = -qJD(4) * t263 + t241 * t472 + t469 * t582;
t51 = t143 * pkin(4) - t144 * pkin(12) + t115;
t9 = -qJD(5) * t784 - t30 * t468 + t471 * t51;
t736 = -t440 * t689 + t468 * t799 + t471 * t797;
t748 = -qJD(5) * t780 + t468 * t801 + t471 * t800;
t348 = t403 * t471 - t587;
t765 = t169 * t704 - t171 * t703 + t347 * t58 - t348 * t59;
t764 = t110 * t347 + t169 * t702 + t224 * t703 + t402 * t59;
t763 = t171 ^ 2;
t762 = t224 ^ 2;
t475 = qJD(1) ^ 2;
t758 = pkin(5) * t110;
t756 = pkin(9) * t467;
t755 = pkin(11) * t469;
t752 = pkin(12) * t171;
t751 = qJ(6) * t702 + qJD(6) * t402 + t749;
t750 = -pkin(5) * t702 - t748;
t747 = pkin(5) * t703 + qJ(6) * t704 - t348 * qJD(6) + t741;
t83 = t469 * t129 + t472 * t138;
t76 = pkin(12) * t532 + t83;
t137 = pkin(3) * t518 - t155;
t90 = t227 * pkin(4) - t229 * pkin(12) + t137;
t35 = t468 * t90 + t471 * t76;
t29 = qJ(6) * t224 + t35;
t746 = t224 * t29;
t745 = t224 * t35;
t740 = qJ(6) * t793 - qJD(6) * t472 + t737;
t739 = -pkin(5) * t793 - t736;
t117 = -t469 * t155 + t209 * t472;
t101 = -pkin(4) * t308 - t117;
t610 = pkin(5) * t468 - qJ(6) * t471;
t580 = pkin(11) + t610;
t738 = -pkin(5) * t205 + qJ(6) * t206 - t101 + (qJD(5) * t611 - qJD(6) * t471) * t469 + t580 * t690;
t735 = -qJD(6) * t468 + t224 * t610 - t83;
t140 = pkin(4) * t229 + pkin(12) * t227;
t57 = t468 * t140 + t471 * t82;
t731 = qJ(6) * t110;
t728 = t110 * t471;
t726 = t171 * t169;
t725 = t224 * t229;
t724 = t229 * t227;
t723 = t269 * t469;
t722 = t269 * t472;
t406 = -t470 * t631 + t473 * t474;
t509 = t535 * t733;
t273 = t406 * t759 + t509 * t761 - t589 * t760;
t721 = t273 * t469;
t720 = t273 * t472;
t717 = t308 * t306;
t712 = t440 * t471;
t711 = t464 * t475;
t710 = t467 * t470;
t706 = t472 * t111;
t34 = -t468 * t76 + t471 * t90;
t705 = qJD(6) - t34;
t384 = pkin(11) * t707 + t468 * t440;
t695 = t474 * pkin(1) + t760 * t756;
t465 = t470 ^ 2;
t466 = t473 ^ 2;
t694 = t465 - t466;
t692 = qJD(2) * t467;
t686 = t137 * qJD(4);
t681 = qJDD(1) * t473;
t678 = t473 * t711;
t259 = t269 * pkin(3);
t677 = -pkin(4) * t722 - pkin(12) * t723 - t259;
t261 = t273 * pkin(3);
t676 = -pkin(4) * t720 - pkin(12) * t721 - t261;
t618 = qJD(2) * t650;
t590 = pkin(1) * t618;
t638 = pkin(1) * t644;
t662 = t467 * t681;
t675 = pkin(9) * t662 + t470 * t638 + t473 * t590;
t670 = t470 * t692;
t669 = qJD(4) + t698;
t668 = t169 ^ 2 - t763;
t665 = pkin(1) * t787;
t658 = t405 * t732;
t657 = t406 * t732;
t655 = t469 * t732;
t653 = t472 * t732;
t642 = pkin(9) * t671;
t303 = -qJD(2) * t642 + t675;
t619 = t733 * t683;
t648 = qJDD(1) * t733;
t215 = (t538 + (-t470 * t619 + t473 * t648) * t467) * pkin(10) + t303;
t447 = t473 * t638;
t541 = -t470 * t590 + t447;
t563 = -t664 - t682;
t542 = t563 * pkin(9);
t221 = t581 * pkin(2) + ((-t470 * t648 - t473 * t619) * pkin(10) + t542) * t467 + t541;
t276 = (qJD(1) * t536 + qJDD(1) * t544) * t467;
t78 = t761 * t215 + t221 * t629 + t276 * t627 - t299 * t666 + t305 * t579 + t354 * t577;
t69 = pkin(11) * t491 + t78;
t145 = -t221 * t732 + t733 * t276;
t77 = -pkin(3) * t479 + t182 * pkin(11) + t145;
t566 = -t129 * t690 + t138 * t691 - t469 * t77 - t472 * t69;
t16 = pkin(12) * t776 - t566;
t571 = t759 * t215 - t221 * t630 - t276 * t628 + t299 * t667 + t305 * t578 + t354 * t576;
t70 = -pkin(3) * t491 + t571;
t25 = t111 * pkin(4) + pkin(12) * t768 + t70;
t4 = -t468 * t16 + t471 * t25 - t76 * t688 - t90 * t689;
t96 = t151 * t472 - t469 * t158;
t280 = -t469 * t387 + t388 * t472;
t647 = t532 * t469;
t645 = t224 * t468;
t643 = t129 * t691 + t138 * t690 + t469 * t69 - t472 * t77;
t641 = t470 * t678;
t639 = t470 * t664;
t632 = -pkin(1) * t760 + t474 * t756;
t625 = t467 * t475 * t734;
t274 = t406 * t761 - t509 * t759 + t588 * t760;
t204 = t274 * t472 + t469 * t794;
t135 = t204 * t468 - t273 * t471;
t623 = g(1) * t131 - g(2) * t135;
t136 = t204 * t471 + t273 * t468;
t622 = g(1) * t132 - g(2) * t136;
t203 = t274 * t469 - t472 * t794;
t621 = g(1) * t791 + g(2) * t203;
t620 = -g(1) * t269 + g(2) * t273;
t267 = pkin(4) * t628 - t280;
t616 = t471 * t690 - t206;
t615 = (qJD(5) * t169 - t58) * pkin(12);
t191 = t264 * t468 - t334 * t471;
t192 = t264 * t471 + t334 * t468;
t612 = -pkin(5) * t191 + qJ(6) * t192;
t86 = qJD(5) * t192 + t144 * t468 - t240 * t471;
t609 = t169 * t86 + t191 * t59;
t28 = -pkin(5) * t224 + t705;
t608 = t28 * t471 - t29 * t468;
t606 = -t34 * t471 - t35 * t468;
t43 = t106 * t471 - t468 * t93;
t56 = t140 * t471 - t468 * t82;
t601 = -t169 * t229 - t728;
t167 = t265 * t471 - t268 * t468;
t592 = 0.2e1 * t650 + qJD(2);
t591 = t377 * pkin(3) + pkin(11) * t376 + t697;
t33 = -t469 * t114 + t125 * t472 - t151 * t691 - t158 * t690;
t575 = -t404 * pkin(2) + pkin(10) * t658;
t574 = -t535 * pkin(2) + pkin(10) * t657;
t573 = g(1) * t474 + g(2) * t760;
t92 = -pkin(4) * t334 - t96;
t3 = t471 * t16 + t468 * t25 + t90 * t688 - t689 * t76;
t8 = t106 * t688 + t471 * t30 + t468 * t51 - t689 * t93;
t567 = t224 * t75 - t753;
t564 = t169 * t703 + t347 * t59;
t159 = -t269 * t708 - t270 * t471;
t161 = -t273 * t708 - t274 * t471;
t222 = -t334 * t708 - t335 * t471;
t562 = g(1) * t161 + g(2) * t159 + g(3) * t222;
t160 = -t269 * t707 + t270 * t468;
t162 = -t273 * t707 + t274 * t468;
t223 = -t334 * t707 + t335 * t468;
t561 = -g(1) * t162 - g(2) * t160 - g(3) * t223;
t296 = -t404 * t761 - t405 * t629;
t231 = t296 * t472 + t405 * t655;
t295 = -t404 * t759 + t405 * t630;
t163 = t231 * t468 - t295 * t471;
t298 = -t406 * t629 - t535 * t761;
t233 = t298 * t472 + t406 * t655;
t297 = t406 * t630 - t535 * t759;
t165 = t233 * t468 - t297 * t471;
t310 = t377 * t472 + t469 * t637;
t234 = t310 * t468 - t376 * t471;
t560 = -g(1) * t165 - g(2) * t163 - g(3) * t234;
t164 = t231 * t471 + t295 * t468;
t166 = t233 * t471 + t297 * t468;
t235 = t310 * t471 + t376 * t468;
t559 = -g(1) * t166 - g(2) * t164 - g(3) * t235;
t558 = g(1) * t203 - g(2) * t791 + g(3) * t263;
t557 = -g(1) * t204 - g(2) * t200 - g(3) * t264;
t230 = t296 * t469 - t405 * t653;
t232 = t298 * t469 - t406 * t653;
t309 = t377 * t469 - t472 * t637;
t556 = -g(1) * t232 - g(2) * t230 - g(3) * t309;
t555 = g(1) * t273 + g(2) * t269 + g(3) * t334;
t554 = g(1) * t274 + g(2) * t270 + g(3) * t335;
t553 = g(1) * t297 + g(2) * t295 + g(3) * t376;
t545 = t169 * t645 - t743;
t543 = -pkin(12) * t743 + t557;
t31 = -pkin(4) * t240 - t33;
t537 = t310 * pkin(4) + pkin(12) * t309 + t591;
t534 = t296 * pkin(3) + t295 * pkin(11) + t575;
t533 = t298 * pkin(3) + t297 * pkin(11) + t574;
t87 = t240 * t468 - t264 * t689 + (qJD(5) * t334 + t144) * t471;
t531 = t169 * t87 + t171 * t86 - t191 * t58 + t192 * t59;
t527 = t110 * t191 + t143 * t169 + t224 * t86 + t263 * t59;
t526 = qJD(2) * t630 + t667;
t524 = pkin(12) * qJD(5) * t224 - t558;
t523 = g(1) * t135 + g(2) * t131 + g(3) * t191 + t4;
t520 = qJD(4) * t532;
t17 = -pkin(4) * t776 + t643;
t5 = t59 * pkin(5) + t58 * qJ(6) - t171 * qJD(6) + t17;
t517 = -t5 - t524;
t516 = -t405 * pkin(2) - pkin(10) * t344 + t632;
t515 = t17 + t524;
t510 = t473 * (qJD(2) * t759 + t578);
t507 = t469 * t643 - t472 * t566 - t554;
t505 = t231 * pkin(4) + t230 * pkin(12) + t534;
t504 = t233 * pkin(4) + t232 * pkin(12) + t533;
t503 = -g(1) * t136 - g(2) * t132 - g(3) * t192 + t3;
t502 = t518 * t732;
t501 = t468 * t469 * t59 + (t469 * t688 + t617) * t169;
t500 = t171 * t38 + qJDD(6) - t523;
t497 = -pkin(3) * t270 - pkin(11) * t269 + t516;
t496 = qJD(3) * t502;
t494 = t406 * pkin(2) + pkin(10) * t794 + t695;
t493 = -pkin(4) * t200 + pkin(12) * t791 + t497;
t489 = t274 * pkin(3) + t273 * pkin(11) + t494;
t488 = t491 * t732;
t486 = t306 * t532 + t520;
t485 = t204 * pkin(4) + t203 * pkin(12) + t489;
t476 = t59 - t646;
t433 = -pkin(4) - t611;
t409 = -pkin(9) * t710 + t459;
t395 = -pkin(9) * t670 + t449;
t390 = t448 - t642;
t389 = t580 * t469;
t383 = -pkin(11) * t708 + t712;
t357 = -t712 + (pkin(11) * t468 + pkin(5)) * t472;
t356 = -qJ(6) * t472 + t384;
t304 = t467 * t542 + t541;
t197 = t203 * pkin(4);
t195 = t791 * pkin(4);
t173 = pkin(5) * t347 - qJ(6) * t348 + t267;
t142 = -pkin(5) * t402 - t167;
t141 = qJ(6) * t402 + t780;
t105 = pkin(5) * t171 + qJ(6) * t169;
t65 = -t110 * t472 + t224 * t647;
t63 = t110 * t402 + t224 * t702;
t52 = -t612 + t92;
t50 = t110 * t263 + t143 * t224;
t42 = -pkin(5) * t229 - t56;
t41 = qJ(6) * t229 + t57;
t40 = -pkin(5) * t263 - t43;
t39 = qJ(6) * t263 + t784;
t37 = -t171 * t229 + t471 * t762 + t729;
t36 = -t58 + t727;
t26 = t471 * t646 - t744;
t22 = -t171 * t704 - t348 * t58;
t21 = -t469 * t471 * t58 + (-t469 * t689 + t616) * t171;
t20 = t171 * t87 - t192 * t58;
t13 = t472 * t58 + t616 * t224 + (t171 * t532 - t224 * t689 + t728) * t469;
t12 = t110 * t348 + t171 * t702 - t224 * t704 - t402 * t58;
t11 = pkin(5) * t86 - qJ(6) * t87 - qJD(6) * t192 + t31;
t10 = t110 * t192 + t143 * t171 + t224 * t87 - t263 * t58;
t7 = -pkin(5) * t143 - t9;
t6 = qJ(6) * t143 + qJD(6) * t263 + t8;
t2 = qJDD(6) - t4 - t758;
t1 = qJD(6) * t224 + t3 + t731;
t14 = [0, 0, 0, 0, 0, qJDD(1), g(1) * t760 - g(2) * t474, t573, 0, 0 (qJDD(1) * t465 + 0.2e1 * t639) * t464 (t470 * t681 - t683 * t694) * t787, t473 * t592 * t692 + t470 * t775 (qJDD(1) * t466 - 0.2e1 * t639) * t464, t473 * t775 - t592 * t670, t581 * t734, -t396 * t593 + t409 * t581 + t304 * t734 + g(1) * t405 - g(2) * t406 + (-t470 * t683 + t681) * t665, -g(1) * t404 + g(2) * t535 - t303 * t734 - t395 * t593 + t563 * t665 - t581 * t696 ((-t390 * qJD(2) + qJDD(1) * t696 + t303 + (-qJD(2) * t409 + t395) * qJD(1)) * t473 + (-t392 * qJD(2) - qJDD(1) * t409 - t304) * t470 - t573) * t467, t464 * qJDD(1) * pkin(1) ^ 2 - g(1) * t632 - g(2) * t695 + t303 * t696 + t304 * t409 - t390 * t396 + t392 * t395, -t182 * t335 + t241 * t308, t182 * t334 - t240 * t308 - t241 * t306 + t335 * t479, t182 * t521 - t241 * t518 + t308 * t582 + t335 * t491, t240 * t306 - t334 * t479, t240 * t518 - t306 * t582 - t334 * t491 - t479 * t521, -t491 * t521 - t502 * t670, g(1) * t270 - g(2) * t274 - t121 * t518 + t145 * t334 + t155 * t582 + t185 * t491 + t220 * t240 - t237 * t479 + t243 * t306 + t521 * t571, t120 * t518 + t145 * t335 - t156 * t582 - t237 * t182 - t186 * t491 + t220 * t241 + t243 * t308 + t521 * t78 + t620, g(1) * t344 - g(2) * t794 - t120 * t306 - t121 * t308 - t155 * t241 - t156 * t240 + t185 * t182 + t186 * t479 - t78 * t334 + t335 * t571, -g(1) * t516 - g(2) * t494 + t156 * t120 + t155 * t121 + t145 * t237 - t185 * t571 + t78 * t186 + t220 * t243, t229 * t144 - t264 * t768, -t264 * t111 - t229 * t143 - t144 * t227 + t263 * t768, t144 * t532 + t229 * t240 - t264 * t805 - t334 * t768, t111 * t263 + t143 * t227, -t143 * t669 + t263 * t495 - t111 * t334 - t227 * t240 + (-t263 * t633 + (-t263 * t510 + (-t143 * t759 - t263 * t526) * t470) * qJD(1)) * t467, -t495 * t334 + t669 * t240 + (t334 * t633 + (t334 * t510 + (t240 * t759 + t334 * t526) * t470) * qJD(1)) * t467, g(1) * t200 - g(2) * t204 + t157 * t111 + t115 * t227 + t137 * t143 + t82 * t240 + t70 * t263 + t33 * t532 - t334 * t643 - t805 * t96, t115 * t229 + t137 * t144 - t157 * t768 - t83 * t240 + t70 * t264 - t32 * t532 + t334 * t566 + t805 * t97 + t621, -t97 * t111 - t83 * t143 - t82 * t144 - t32 * t227 - t33 * t229 + t263 * t566 + t264 * t643 + t768 * t96 - t620, -g(1) * t497 - g(2) * t489 + t137 * t115 + t70 * t157 + t83 * t32 + t82 * t33 - t566 * t97 - t643 * t96, t20, -t531, t10, t609, -t527, t50, t110 * t43 + t143 * t34 + t169 * t31 + t17 * t191 + t224 * t9 + t263 * t4 + t59 * t92 + t75 * t86 + t622, -t110 * t784 - t143 * t35 + t17 * t192 + t171 * t31 - t224 * t8 - t263 * t3 - t58 * t92 + t75 * t87 - t623, -t169 * t8 - t171 * t9 - t191 * t3 - t192 * t4 - t34 * t87 - t35 * t86 + t43 * t58 - t59 * t784 - t621, -g(1) * t493 - g(2) * t485 + t17 * t92 + t3 * t784 + t75 * t31 + t34 * t9 + t35 * t8 + t4 * t43, t20, t10, t531, t50, t527, t609, t11 * t169 - t110 * t40 - t143 * t28 + t191 * t5 - t2 * t263 - t224 * t7 + t38 * t86 + t52 * t59 + t622, -t1 * t191 - t169 * t6 + t171 * t7 + t192 * t2 + t28 * t87 - t29 * t86 - t39 * t59 - t40 * t58 - t621, t1 * t263 - t11 * t171 + t110 * t39 + t143 * t29 - t192 * t5 + t224 * t6 - t38 * t87 + t52 * t58 + t623, t1 * t39 + t29 * t6 + t5 * t52 + t38 * t11 + t2 * t40 + t28 * t7 - g(1) * (-pkin(5) * t132 - qJ(6) * t131 + t493) - g(2) * (t136 * pkin(5) + t135 * qJ(6) + t485); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t641, t694 * t711, -t473 * t625 + t663, t641, t470 * t625 + t662, t581, t447 + t392 * t593 + g(1) * t535 + g(2) * t404 + (-t618 + t711) * t470 * pkin(1) + (-g(3) * t473 + t542) * t467, pkin(1) * t678 + t390 * t593 + g(1) * t406 + g(2) * t405 + (pkin(9) * t683 + g(3)) * t710 - t675, 0, 0, -t182 * t627 + t308 * t548, -t182 * t628 + t479 * t627 + t362 * t306 + t308 * t361 + (-t306 * t628 - t308 * t627) * qJD(3), -t182 * t733 - t308 * t583 + t362 * t518 + t488 * t759 - t496 * t761, t306 * t547 + t479 * t628, t306 * t583 - t361 * t518 + t479 * t733 + t488 * t761 + t496 * t759, t512 * t733 - (qJD(1) * t614 - t684) * t583, -g(1) * t298 - g(2) * t296 - g(3) * t377 - t145 * t628 - t155 * t583 + t220 * t547 - t242 * t306 + t408 * t491 + t479 * t673 + t518 * t700 - t571 * t733, t145 * t627 + t156 * t583 + t182 * t673 + t220 * t548 - t242 * t308 - t410 * t491 + t518 * t701 - t733 * t78 + t553, -g(3) * t637 + t78 * t628 + t571 * t627 - g(1) * t657 - g(2) * t658 + t155 * t362 + t156 * t361 + t408 * t182 + t410 * t479 + t700 * t308 - t701 * t306 + (-t155 * t628 - t156 * t627) * qJD(3), -g(1) * t574 - g(2) * t575 - g(3) * t697 - t145 * t673 - t155 * t700 + t156 * t701 - t220 * t242 - t408 * t571 + t78 * t410, -t229 * t778 - t403 * t768, -t403 * t111 + t227 * t778 - t229 * t702 + t402 * t768, t229 * t547 - t403 * t805 - t532 * t778 + t628 * t768, t111 * t402 + t227 * t702, t402 * t495 + t111 * t628 - t547 * t227 + (-t402 * t633 + (-t402 * t510 + (-t402 * t526 - t702 * t759) * t470) * qJD(1)) * t467 - t702 * t669, t532 * t547 + t628 * t805, -g(1) * t233 - g(2) * t231 - g(3) * t310 + t386 * t111 + t137 * t702 + t227 * t777 - t280 * t805 + t70 * t402 + t532 * t781 + t547 * t82 + t628 * t643, -t137 * t778 + t229 * t777 + t281 * t805 - t386 * t768 + t70 * t403 - t532 * t782 - t547 * t83 - t566 * t628 - t556, -t281 * t111 - t227 * t782 - t229 * t781 + t280 * t768 + t402 * t566 + t403 * t643 - t702 * t83 + t778 * t82 - t553, -g(1) * t533 - g(2) * t534 - g(3) * t591 + t137 * t777 - t280 * t643 - t281 * t566 + t70 * t386 + t781 * t82 + t782 * t83, t22, t765, t12, t564, -t764, t63, t110 * t167 + t169 * t741 + t17 * t347 + t224 * t748 + t267 * t59 + t34 * t702 + t4 * t402 + t703 * t75 + t559, -t110 * t780 + t17 * t348 + t171 * t741 - t224 * t749 - t267 * t58 - t3 * t402 - t35 * t702 - t704 * t75 - t560, t167 * t58 - t169 * t749 - t171 * t748 - t3 * t347 + t34 * t704 - t348 * t4 - t35 * t703 - t59 * t780 + t556, -g(1) * t504 - g(2) * t505 - g(3) * t537 + t4 * t167 + t17 * t267 + t3 * t780 + t34 * t748 + t35 * t749 + t741 * t75, t22, t12, -t765, t63, t764, t564, -t110 * t142 + t169 * t747 + t173 * t59 - t2 * t402 - t224 * t750 - t28 * t702 + t347 * t5 + t38 * t703 + t559, -t1 * t347 - t141 * t59 - t142 * t58 - t169 * t751 + t171 * t750 + t2 * t348 - t28 * t704 - t29 * t703 + t556, t1 * t402 + t110 * t141 - t171 * t747 + t173 * t58 + t224 * t751 + t29 * t702 - t348 * t5 + t38 * t704 + t560, t1 * t141 + t5 * t173 + t2 * t142 - g(1) * (t166 * pkin(5) + t165 * qJ(6) + t504) - g(2) * (t164 * pkin(5) + t163 * qJ(6) + t505) - g(3) * (pkin(5) * t235 + qJ(6) * t234 + t537) + t747 * t38 + t751 * t29 + t750 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t717, -t306 ^ 2 + t308 ^ 2, -t306 * t518 - t182, -t717, -t308 * t518 + t479, t491, -t156 * t518 - t220 * t308 + t555 - t571, -t155 * t518 + t220 * t306 + t554 - t78, 0, 0, -t487 * t469 + (t229 * t532 - t469 * t498) * t472, -t229 * t719 - t472 * t768 + (-t111 - t685) * t469 + t774 * t227, -t229 * t308 + t469 * t772 + t472 * t486, t227 * t647 - t706, t227 * t308 - t469 * t486 + t472 * t772, -t532 * t308, -pkin(3) * t111 + g(1) * t720 + g(2) * t722 + g(3) * t715 - t117 * t532 + t137 * t719 - t156 * t227 - t82 * t308 + t469 * t686 - t70 * t472 - t520 * t754 + t755 * t805, pkin(3) * t768 - g(1) * t721 - g(2) * t723 - g(3) * t716 + t137 * t718 - t156 * t229 + t83 * t308 + t70 * t469 + t472 * t686 + t532 * t783 + t754 * t805, -pkin(11) * t706 + t117 * t229 + t783 * t227 + t685 * t754 - t768 * t755 + t774 * t82 - t793 * t83 + t507, -t70 * pkin(3) + g(1) * t261 + g(2) * t259 + g(3) * t333 - t82 * t117 - t83 * t118 - t137 * t156 + ((-t83 * t469 - t82 * t472) * qJD(4) + t507) * pkin(11), t21, t790, t13, t501, -t789, t65, -t101 * t169 + t110 * t383 - t205 * t75 + t736 * t224 + (-t4 + (pkin(11) * t169 + t468 * t75) * qJD(4)) * t472 + (pkin(11) * t59 + t17 * t468 + t34 * t532 + t688 * t75) * t469 + t561, -t101 * t171 - t110 * t384 - t206 * t75 - t737 * t224 + (t3 + (pkin(11) * t171 + t471 * t75) * qJD(4)) * t472 + (-pkin(11) * t58 + t17 * t471 - t35 * t532 - t689 * t75) * t469 + t562, t205 * t35 + t206 * t34 + t383 * t58 - t384 * t59 - t736 * t171 - t737 * t169 + t606 * t690 + (-t3 * t468 - t4 * t471 + (t34 * t468 - t35 * t471) * qJD(5) + t555) * t469, t3 * t384 + t4 * t383 - t75 * t101 - g(1) * t676 - g(2) * t677 - g(3) * (-t333 + t779) + t737 * t35 + t736 * t34 + (t17 * t469 + t690 * t75 - t554) * pkin(11), t21, t13, -t790, t65, t789, t501, -t110 * t357 + t2 * t472 + t389 * t59 + t617 * t38 - t739 * t224 + t738 * t169 + (-t28 * t532 + t38 * t688 + t468 * t5) * t469 + t561, t205 * t29 - t206 * t28 - t356 * t59 - t357 * t58 + t739 * t171 - t740 * t169 + t608 * t690 + (-t1 * t468 + t2 * t471 + (-t28 * t468 - t29 * t471) * qJD(5) + t555) * t469, -t1 * t472 + t110 * t356 + t389 * t58 - t616 * t38 + t740 * t224 - t738 * t171 + (t29 * t532 + t38 * t689 - t471 * t5) * t469 - t562, t1 * t356 + t5 * t389 + t2 * t357 - g(1) * (pkin(5) * t162 + pkin(11) * t274 + qJ(6) * t161 + t676) - g(2) * (pkin(5) * t160 + pkin(11) * t270 + qJ(6) * t159 + t677) - g(3) * (pkin(5) * t223 + qJ(6) * t222 - t660 + t779) + t738 * t38 + t740 * t29 + t739 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t724, -t227 ^ 2 + t229 ^ 2, t227 * t532 - t768, -t724, t229 * t306 + t649, -t805, -t137 * t229 + t532 * t83 + t558 - t643, t137 * t227 + t532 * t82 - t557 + t566, 0, 0, t26, -t788, t37, t545, -t468 * t762 - t601, -t725, -pkin(4) * t59 - t169 * t83 - t224 * t56 - t229 * t34 + t468 * t567 - t471 * t515, pkin(4) * t58 - t171 * t83 + t224 * t57 + t229 * t35 + t468 * t515 + t471 * t567, t169 * t57 + t171 * t56 + (-t227 * t34 + t3 + (-t34 + t752) * qJD(5)) * t471 + (-t4 + t615 - t745) * t468 + t543, -t17 * pkin(4) + g(1) * t197 - g(2) * t195 + g(3) * t258 - t34 * t56 - t35 * t57 - t75 * t83 + (qJD(5) * t606 + t3 * t471 - t4 * t468 + t557) * pkin(12), t26, t37, t788, -t725, t224 * t645 + t601, t545, t735 * t169 + t224 * t42 + t229 * t28 + t433 * t59 + t468 * t771 + t517 * t471, t169 * t41 - t171 * t42 + (t227 * t28 + t1 + (t28 + t752) * qJD(5)) * t471 + (t2 + t615 - t746) * t468 + t543, -t735 * t171 - t224 * t41 - t229 * t29 + t433 * t58 + t517 * t468 - t471 * t771, t5 * t433 - t29 * t41 - t28 * t42 - g(1) * (-t203 * t611 - t197) - g(2) * (t611 * t791 + t195) - g(3) * (-t263 * t611 - t258) + t735 * t38 + (qJD(5) * t608 + t1 * t471 + t2 * t468 + t557) * pkin(12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t726, -t668, t36, -t726, -t476, t110, -t171 * t75 + t523 + t745, t169 * t75 + t224 * t34 - t503, 0, 0, t726, t36, t668, t110, t476, -t726, -t105 * t169 - t500 + t745 + 0.2e1 * t758, pkin(5) * t58 - qJ(6) * t59 + (t29 - t35) * t171 + (t28 - t705) * t169, 0.2e1 * t731 + t105 * t171 - t169 * t38 + (0.2e1 * qJD(6) - t34) * t224 + t503, t1 * qJ(6) - t2 * pkin(5) - t38 * t105 - t28 * t35 - g(1) * (-pkin(5) * t135 + qJ(6) * t136) - g(2) * (-pkin(5) * t131 + qJ(6) * t132) - g(3) * t612 + t705 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t726 - t110, t36, -t762 - t763, t500 - t746 - t758;];
tau_reg  = t14;