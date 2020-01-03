% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP6_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP6_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP6_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:09
% EndTime: 2019-12-31 18:43:07
% DurationCPUTime: 48.98s
% Computational Cost: add. (26254->898), mult. (34997->1193), div. (0->0), fcn. (33148->8), ass. (0->449)
t818 = -Icges(5,4) - Icges(6,4);
t775 = Icges(5,1) + Icges(6,1);
t815 = Icges(5,5) + Icges(6,5);
t795 = -Icges(5,2) - Icges(6,2);
t794 = Icges(5,6) + Icges(6,6);
t375 = cos(qJ(4));
t817 = t818 * t375;
t372 = sin(qJ(4));
t816 = t818 * t372;
t814 = Icges(5,3) + Icges(6,3);
t813 = -t794 * t372 + t815 * t375;
t812 = t795 * t372 - t817;
t811 = t775 * t375 + t816;
t370 = qJ(1) + pkin(8);
t366 = sin(t370);
t376 = cos(qJ(3));
t604 = t372 * t376;
t367 = cos(t370);
t608 = t367 * t375;
t261 = t366 * t604 + t608;
t603 = t375 * t376;
t610 = t367 * t372;
t262 = t366 * t603 - t610;
t373 = sin(qJ(3));
t613 = t366 * t373;
t132 = Icges(6,5) * t262 - Icges(6,6) * t261 + Icges(6,3) * t613;
t135 = Icges(5,5) * t262 - Icges(5,6) * t261 + Icges(5,3) * t613;
t793 = t132 + t135;
t612 = t366 * t375;
t263 = -t367 * t604 + t612;
t614 = t366 * t372;
t264 = t367 * t603 + t614;
t609 = t367 * t373;
t134 = Icges(6,5) * t264 + Icges(6,6) * t263 + Icges(6,3) * t609;
t137 = Icges(5,5) * t264 + Icges(5,6) * t263 + Icges(5,3) * t609;
t773 = t134 + t137;
t248 = Icges(6,4) * t262;
t138 = -Icges(6,2) * t261 + Icges(6,6) * t613 + t248;
t251 = Icges(5,4) * t262;
t141 = -Icges(5,2) * t261 + Icges(5,6) * t613 + t251;
t792 = t138 + t141;
t637 = Icges(6,4) * t264;
t140 = Icges(6,2) * t263 + Icges(6,6) * t609 + t637;
t640 = Icges(5,4) * t264;
t143 = Icges(5,2) * t263 + Icges(5,6) * t609 + t640;
t791 = t140 + t143;
t247 = Icges(6,4) * t261;
t145 = -Icges(6,1) * t262 - Icges(6,5) * t613 + t247;
t250 = Icges(5,4) * t261;
t148 = -Icges(5,1) * t262 - Icges(5,5) * t613 + t250;
t790 = -t145 - t148;
t249 = Icges(6,4) * t263;
t146 = Icges(6,1) * t264 + Icges(6,5) * t609 + t249;
t252 = Icges(5,4) * t263;
t149 = Icges(5,1) * t264 + Icges(5,5) * t609 + t252;
t789 = t146 + t149;
t772 = t813 * t373 - t814 * t376;
t810 = t814 * t373 + t813 * t376;
t803 = t812 * t373 - t794 * t376;
t719 = t794 * t373 + t812 * t376;
t770 = t811 * t373 - t815 * t376;
t718 = t815 * t373 + t811 * t376;
t809 = (-t815 * t372 - t794 * t375) * t373;
t808 = (t795 * t375 + t816) * t373;
t807 = (-t775 * t372 + t817) * t373;
t779 = t792 * t263 + t790 * t264 + t609 * t793;
t778 = t263 * t791 + t264 * t789 + t609 * t773;
t730 = t263 * t803 + t264 * t770 + t609 * t772;
t806 = qJD(3) * t810 + qJD(4) * t809;
t805 = qJD(3) * t719 + qJD(4) * t808;
t804 = qJD(3) * t718 + qJD(4) * t807;
t781 = -t792 * t261 + t790 * t262 + t613 * t793;
t780 = -t261 * t791 + t262 * t789 + t613 * t773;
t558 = qJD(4) * t376;
t355 = qJD(1) - t558;
t561 = qJD(3) * t373;
t402 = t355 * t375 + t372 * t561;
t564 = qJD(1) * t376;
t497 = -qJD(4) + t564;
t127 = t367 * t402 + t497 * t614;
t440 = t355 * t372;
t401 = -t375 * t561 + t440;
t128 = t367 * t401 - t497 * t612;
t560 = qJD(3) * t376;
t532 = t367 * t560;
t565 = qJD(1) * t373;
t538 = t366 * t565;
t418 = t532 - t538;
t68 = Icges(6,5) * t128 + Icges(6,6) * t127 + Icges(6,3) * t418;
t70 = Icges(5,5) * t128 + Icges(5,6) * t127 + Icges(5,3) * t418;
t802 = t68 + t70;
t129 = t366 * t402 - t497 * t610;
t130 = t366 * t401 + t497 * t608;
t420 = t366 * t560 + t367 * t565;
t69 = Icges(6,5) * t130 + Icges(6,6) * t129 + Icges(6,3) * t420;
t71 = Icges(5,5) * t130 + Icges(5,6) * t129 + Icges(5,3) * t420;
t801 = t69 + t71;
t72 = Icges(6,4) * t128 + Icges(6,2) * t127 + Icges(6,6) * t418;
t74 = Icges(5,4) * t128 + Icges(5,2) * t127 + Icges(5,6) * t418;
t800 = t72 + t74;
t73 = Icges(6,4) * t130 + Icges(6,2) * t129 + Icges(6,6) * t420;
t75 = Icges(5,4) * t130 + Icges(5,2) * t129 + Icges(5,6) * t420;
t799 = t73 + t75;
t76 = Icges(6,1) * t128 + Icges(6,4) * t127 + Icges(6,5) * t418;
t78 = Icges(5,1) * t128 + Icges(5,4) * t127 + Icges(5,5) * t418;
t798 = t76 + t78;
t77 = Icges(6,1) * t130 + Icges(6,4) * t129 + Icges(6,5) * t420;
t79 = Icges(5,1) * t130 + Icges(5,4) * t129 + Icges(5,5) * t420;
t797 = t77 + t79;
t731 = -t261 * t803 + t262 * t770 + t613 * t772;
t788 = -t372 * t803 + t375 * t770;
t759 = t141 * t372 + t148 * t375;
t760 = t138 * t372 + t145 * t375;
t787 = t759 + t760;
t559 = qJD(4) * t373;
t563 = qJD(3) * t366;
t295 = t367 * t559 + t563;
t562 = qJD(3) * t367;
t296 = -t366 * t559 + t562;
t733 = t778 * t295 - t296 * t779 + t730 * t355;
t784 = t127 * t803 + t128 * t770 + t263 * t805 + t264 * t804 + t418 * t772 + t609 * t806;
t783 = t129 * t803 + t130 * t770 - t261 * t805 + t262 * t804 + t420 * t772 + t613 * t806;
t757 = t132 * t376;
t56 = -t373 * t760 - t757;
t754 = t135 * t376;
t58 = -t373 * t759 - t754;
t777 = t56 + t58;
t786 = t295 * (t264 * t795 + t249 + t252 + t789) - t296 * (t262 * t795 - t247 - t250 + t790) + t355 * (t770 + t808);
t449 = -t143 * t372 + t149 * t375;
t451 = -t140 * t372 + t146 * t375;
t785 = t295 * (-t367 * t772 - t449 - t451) - t296 * (-t366 * t772 + t787) + t355 * (-t788 + t810);
t741 = t127 * t792 + t128 * t790 + t263 * t799 + t264 * t797 + t418 * t793 + t609 * t801;
t740 = t127 * t791 + t128 * t789 + t263 * t800 + t264 * t798 + t418 * t773 + t609 * t802;
t738 = t129 * t792 + t130 * t790 - t261 * t799 + t262 * t797 + t420 * t793 + t613 * t801;
t737 = t129 * t791 + t130 * t789 - t261 * t800 + t262 * t798 + t420 * t773 + t613 * t802;
t782 = (qJD(3) * t788 - t806) * t376 + (t804 * t375 - t805 * t372 + (-t372 * t770 - t375 * t803) * qJD(4) + t772 * qJD(3)) * t373;
t57 = -t134 * t376 + t373 * t451;
t59 = -t137 * t376 + t373 * t449;
t776 = t57 + t59;
t729 = t373 * t788 - t376 * t772;
t365 = pkin(4) * t375 + pkin(3);
t611 = t366 * t376;
t769 = -rSges(6,1) * t262 + rSges(6,2) * t261 - t365 * t611;
t557 = qJD(5) * t373;
t334 = t367 * t557;
t552 = pkin(4) * t610;
t675 = pkin(3) * t376;
t371 = -qJ(5) - pkin(7);
t671 = pkin(7) + t371;
t726 = t373 * t671;
t598 = t552 + (t675 + t726) * t366 - rSges(6,3) * t613 + t769;
t768 = t355 * t598 + t334;
t767 = -t809 * t355 + (-t261 * t815 - t262 * t794) * t296 + (-t263 * t815 + t264 * t794) * t295;
t378 = qJD(1) ^ 2;
t758 = rSges(6,1) + pkin(4);
t734 = t780 * t295 - t296 * t781 + t731 * t355;
t750 = t366 * t367;
t359 = t366 * rSges(4,3);
t607 = t367 * t376;
t233 = rSges(4,1) * t607 - rSges(4,2) * t609 + t359;
t313 = t367 * pkin(2) + t366 * pkin(6);
t377 = cos(qJ(1));
t369 = t377 * pkin(1);
t717 = t369 + t313;
t749 = t233 + t717;
t672 = pkin(3) - t365;
t725 = t376 * t672;
t244 = -t725 - t726;
t474 = rSges(6,1) * t375 - rSges(6,2) * t372;
t662 = rSges(6,3) * t373;
t748 = t376 * t474 + t244 + t662;
t524 = t671 * t376;
t525 = t672 * t373;
t747 = rSges(6,3) * t376 - t373 * t474 - t524 + t525;
t351 = pkin(3) * t607;
t292 = pkin(7) * t609 + t351;
t347 = pkin(3) * t373 - pkin(7) * t376;
t416 = (t292 + t717) * qJD(1) - t347 * t563;
t531 = t366 * t557;
t606 = t371 * t373;
t724 = t264 * rSges(6,1) + t263 * rSges(6,2) + rSges(6,3) * t609 + pkin(4) * t614 + t365 * t607;
t597 = -t367 * t606 - t292 + t724;
t45 = t295 * t747 + t355 * t597 + t416 + t531;
t504 = t45 * t747;
t673 = pkin(7) * t373;
t348 = t673 + t675;
t290 = t348 * t366;
t527 = t290 * t563 + t292 * t562 + qJD(2);
t556 = qJD(5) * t376;
t37 = -t295 * t598 + t296 * t597 + t527 - t556;
t506 = t37 * t598;
t746 = -t504 + t506;
t478 = rSges(5,1) * t262 - rSges(5,2) * t261;
t153 = rSges(5,3) * t613 + t478;
t477 = rSges(5,1) * t375 - rSges(5,2) * t372;
t286 = -rSges(5,3) * t376 + t373 * t477;
t745 = -t153 * t355 - t286 * t296;
t744 = 0.2e1 * qJD(3);
t553 = qJD(3) * qJD(4);
t521 = t376 * t553;
t223 = qJD(1) * t295 + t366 * t521;
t224 = qJD(1) * t296 + t367 * t521;
t522 = t373 * t553;
t743 = t779 * t223 + t778 * t224 + t740 * t295 - t296 * t741 + t784 * t355 + t730 * t522;
t742 = t223 * t781 + t224 * t780 + t295 * t737 - t296 * t738 + t355 * t783 + t522 * t731;
t326 = pkin(7) * t532;
t533 = t367 * t561;
t419 = -t366 * t564 - t533;
t192 = pkin(3) * t419 - pkin(7) * t538 + t326;
t325 = t366 * pkin(3) * t561;
t193 = pkin(7) * t420 + qJD(1) * t351 - t325;
t696 = qJD(1) * t290;
t544 = t193 * t563 + (t192 + t696) * t562;
t567 = qJD(1) * t366;
t476 = rSges(6,1) * t130 + rSges(6,2) * t129;
t659 = pkin(4) * qJD(4);
t547 = t375 * t659;
t615 = t365 * t373;
t667 = t325 + (qJD(1) * t244 - t547) * t367 + (t557 + pkin(4) * t440 + (-t524 - t615) * qJD(3)) * t366 + rSges(6,3) * t420 + t476;
t548 = t372 * t659;
t498 = t376 * t548;
t605 = t371 * t376;
t702 = t128 * rSges(6,1) + t127 * rSges(6,2) + rSges(6,3) * t532 + qJD(1) * t552 + t366 * t547 + t371 * t538 + t334;
t668 = -t326 + (t673 + t725) * t567 + (-t498 + (t525 - t605) * qJD(3)) * t367 - rSges(6,3) * t538 + t702;
t5 = t668 * t296 + t667 * t295 - t598 * t224 - t597 * t223 + (-t292 * t567 + t557) * qJD(3) + t544;
t739 = t5 * t597;
t16 = (-qJD(3) * t760 - t69) * t376 + (qJD(3) * t132 - t372 * t73 + t375 * t77 + (-t138 * t375 + t145 * t372) * qJD(4)) * t373;
t18 = (-qJD(3) * t759 - t71) * t376 + (qJD(3) * t135 - t372 * t75 + t375 * t79 + (-t141 * t375 + t148 * t372) * qJD(4)) * t373;
t736 = t16 + t18;
t17 = (qJD(3) * t451 - t68) * t376 + (qJD(3) * t134 - t372 * t72 + t375 * t76 + (-t140 * t375 - t146 * t372) * qJD(4)) * t373;
t19 = (qJD(3) * t449 - t70) * t376 + (qJD(3) * t137 - t372 * t74 + t375 * t78 + (-t143 * t375 - t149 * t372) * qJD(4)) * t373;
t735 = t17 + t19;
t732 = t295 * t776 - t296 * t777 + t355 * t729;
t723 = t803 * t366;
t722 = t803 * t367;
t721 = t770 * t366;
t720 = t770 * t367;
t502 = t367 * rSges(3,1) - rSges(3,2) * t366;
t368 = Icges(4,4) * t376;
t458 = -Icges(4,2) * t373 + t368;
t332 = Icges(4,1) * t373 + t368;
t716 = t369 + t502;
t715 = t785 * t373;
t713 = (-t803 + t807) * t355 + (t261 * t775 + t248 + t251 + t792) * t296 + (t263 * t775 - t637 - t640 - t791) * t295;
t712 = t767 * t373;
t711 = t295 * t773 - t296 * t793 + t355 * t772;
t710 = t366 * t777 + t367 * t776;
t709 = t366 * t776 - t367 * t777;
t708 = t366 * t779 + t367 * t778;
t707 = t366 * t778 - t367 * t779;
t706 = t366 * t781 + t367 * t780;
t705 = t366 * t780 - t367 * t781;
t320 = qJD(3) * t348;
t307 = (-rSges(6,1) * t372 - rSges(6,2) * t375) * t373;
t591 = qJD(3) * t748 + qJD(4) * t307 - t373 * t548 - t556;
t701 = -t320 - t591 - t556;
t289 = t347 * t366;
t291 = t347 * t367;
t566 = qJD(1) * t367;
t700 = t367 * t192 + t366 * t193 + t289 * t563 + t290 * t566 + t291 * t562;
t697 = t355 * t782 + t522 * t729;
t627 = Icges(4,3) * t367;
t226 = Icges(4,5) * t611 - Icges(4,6) * t613 - t627;
t342 = Icges(4,4) * t613;
t634 = Icges(4,5) * t367;
t230 = Icges(4,1) * t611 - t342 - t634;
t630 = Icges(4,6) * t367;
t228 = Icges(4,4) * t611 - Icges(4,2) * t613 - t630;
t620 = t228 * t373;
t447 = -t230 * t376 + t620;
t92 = -t226 * t367 - t366 * t447;
t329 = Icges(4,5) * t376 - Icges(4,6) * t373;
t328 = Icges(4,5) * t373 + Icges(4,6) * t376;
t422 = qJD(3) * t328;
t641 = Icges(4,4) * t373;
t333 = Icges(4,1) * t376 - t641;
t231 = Icges(4,5) * t366 + t333 * t367;
t229 = Icges(4,6) * t366 + t367 * t458;
t619 = t229 * t373;
t446 = -t231 * t376 + t619;
t694 = -t367 * t422 + (-t329 * t366 + t446 + t627) * qJD(1);
t227 = Icges(4,3) * t366 + t329 * t367;
t693 = -t366 * t422 + (t227 + t447) * qJD(1);
t330 = Icges(4,2) * t376 + t641;
t442 = t330 * t373 - t332 * t376;
t692 = qJD(1) * t442 + t329 * qJD(3);
t584 = -Icges(4,2) * t611 + t230 - t342;
t586 = t332 * t366 + t228;
t691 = -t373 * t584 - t376 * t586;
t690 = -m(6) / 0.2e1;
t689 = m(6) / 0.2e1;
t688 = t223 / 0.2e1;
t687 = t224 / 0.2e1;
t686 = -t295 / 0.2e1;
t685 = t295 / 0.2e1;
t684 = -t296 / 0.2e1;
t683 = t296 / 0.2e1;
t682 = -t355 / 0.2e1;
t681 = t355 / 0.2e1;
t677 = -rSges(5,3) - pkin(7);
t374 = sin(qJ(1));
t676 = pkin(1) * t374;
t674 = pkin(4) * t372;
t666 = rSges(4,1) * t376;
t664 = rSges(5,3) * t373;
t660 = pkin(1) * qJD(1);
t658 = t16 * t296;
t657 = t17 * t295;
t656 = t18 * t296;
t655 = t19 * t295;
t363 = t367 * pkin(6);
t312 = pkin(2) * t366 - t363;
t510 = -t312 - t676;
t536 = t347 * t562;
t393 = (-t290 + t510) * qJD(1) - t536;
t62 = t393 + t745;
t648 = t367 * t62;
t647 = t56 * t223;
t646 = t57 * t224;
t645 = t58 * t223;
t644 = t59 * t224;
t643 = -rSges(6,3) + t371;
t570 = rSges(4,2) * t613 + t367 * rSges(4,3);
t232 = rSges(4,1) * t611 - t570;
t336 = rSges(4,1) * t373 + rSges(4,2) * t376;
t534 = t336 * t562;
t113 = -t534 + (-t232 + t510) * qJD(1);
t624 = t113 * t366;
t623 = t113 * t367;
t537 = t336 * t563;
t114 = qJD(1) * t749 - t537;
t284 = t336 * t367;
t622 = t114 * t284;
t617 = t328 * t366;
t616 = t328 * t367;
t157 = t264 * rSges(5,1) + t263 * rSges(5,2) + rSges(5,3) * t609;
t596 = -t157 - t292;
t595 = t747 * t366;
t594 = t747 * t367;
t593 = -rSges(6,2) * t262 - t261 * t758;
t592 = -rSges(6,2) * t264 + t263 * t758;
t297 = t313 * qJD(1);
t590 = -t193 - t297;
t288 = t376 * t477 + t664;
t308 = (-rSges(5,1) * t372 - rSges(5,2) * t375) * t373;
t201 = qJD(3) * t288 + qJD(4) * t308;
t589 = -t201 - t320;
t588 = -t366 * t226 - t230 * t607;
t587 = t366 * t227 + t231 * t607;
t585 = -t332 * t367 - t229;
t583 = -t330 * t367 + t231;
t580 = t366 * t290 + t367 * t292;
t575 = -t286 - t347;
t573 = rSges(4,2) * t538 + rSges(4,3) * t566;
t572 = -t330 + t333;
t571 = t332 + t458;
t568 = qJD(1) * t329;
t164 = -t366 * t442 - t616;
t555 = t164 * qJD(1);
t554 = qJD(1) * qJD(3);
t551 = t378 * t676;
t550 = t378 * t369;
t549 = t374 * t660;
t545 = t128 * rSges(5,1) + t127 * rSges(5,2) + rSges(5,3) * t532;
t541 = -t347 + t747;
t539 = t347 * t566;
t526 = -pkin(2) - t666;
t523 = t366 * t554;
t519 = t566 / 0.2e1;
t518 = -t563 / 0.2e1;
t515 = t562 / 0.2e1;
t514 = t561 / 0.2e1;
t513 = t560 / 0.2e1;
t508 = t363 - t676;
t507 = -t365 * t376 - pkin(2);
t44 = t296 * t747 + t393 + t768;
t505 = t44 * t747;
t202 = t231 * t611;
t501 = t227 * t367 - t202;
t500 = -t226 + t619;
t499 = -t320 + t556;
t354 = pkin(6) * t566;
t489 = qJD(1) * (-pkin(2) * t567 + t354) - t551;
t488 = t347 * t523 - t550;
t487 = t373 * t674 - t307;
t484 = qJD(4) * t514;
t483 = -qJD(1) * t291 - t348 * t563;
t309 = rSges(3,1) * t366 + rSges(3,2) * t367;
t480 = -rSges(4,2) * t373 + t666;
t479 = rSges(5,1) * t130 + rSges(5,2) * t129;
t453 = -t114 * t366 - t623;
t448 = t153 * t367 - t157 * t366;
t121 = t228 * t376 + t230 * t373;
t122 = t229 * t376 + t231 * t373;
t445 = t232 * t366 + t233 * t367;
t441 = qJD(1) * t192 + t489;
t283 = t336 * t366;
t429 = qJD(1) * t289 - t348 * t562;
t93 = -t229 * t613 - t501;
t427 = (t366 * t93 - t367 * t92) * qJD(3);
t94 = -t228 * t609 - t588;
t95 = -t229 * t609 + t587;
t426 = (t366 * t95 - t367 * t94) * qJD(3);
t425 = t373 * t677 - pkin(2) - t675;
t424 = qJD(3) * t332;
t423 = qJD(3) * t330;
t417 = t37 * t667 - t5 * t598;
t410 = t44 * t598 + t45 * t597;
t409 = -t37 * t597 - t505;
t408 = -t373 * t583 + t376 * t585;
t300 = qJD(1) * t312;
t407 = -t300 - t536 - t549 - t696;
t394 = (-t373 * t571 + t376 * t572) * qJD(1);
t166 = rSges(4,1) * t419 - rSges(4,2) * t532 + t573;
t167 = -qJD(3) * t283 + (t367 * t480 + t359) * qJD(1);
t388 = t166 * t367 + t167 * t366 + (t232 * t367 - t233 * t366) * qJD(1);
t54 = t153 * t295 + t157 * t296 + t527;
t63 = t157 * t355 - t286 * t295 + t416;
t385 = t54 * t448 + (t366 * t62 - t367 * t63) * t286;
t315 = t458 * qJD(3);
t316 = t333 * qJD(3);
t382 = qJD(1) * t328 - t315 * t373 + t316 * t376 + (-t330 * t376 - t332 * t373) * qJD(3);
t381 = t409 * t366 - t367 * t746;
t318 = t480 * qJD(3);
t299 = t347 * t567;
t220 = t286 * t367;
t218 = t286 * t366;
t190 = rSges(5,1) * t263 - rSges(5,2) * t264;
t188 = -rSges(5,1) * t261 - rSges(5,2) * t262;
t165 = -t367 * t442 + t617;
t131 = t165 * qJD(1);
t112 = qJD(3) * t445 + qJD(2);
t87 = -t550 - t318 * t562 + (-t167 - t297 + t537) * qJD(1);
t86 = -t318 * t563 + (t166 - t534) * qJD(1) + t489;
t83 = rSges(5,3) * t420 + t479;
t81 = -rSges(5,3) * t538 + t545;
t65 = t382 * t366 - t367 * t692;
t64 = t366 * t692 + t382 * t367;
t61 = -qJD(3) * t446 + (-t367 * t423 + (-t366 * t458 + t630) * qJD(1)) * t376 + (-t367 * t424 + (-t333 * t366 + t634) * qJD(1)) * t373;
t60 = -qJD(3) * t447 + (qJD(1) * t229 - t366 * t423) * t376 + (qJD(1) * t231 - t366 * t424) * t373;
t55 = t388 * qJD(3);
t41 = t131 + t426;
t40 = t427 + t555;
t28 = -t201 * t296 + t223 * t286 - t355 * t83 + (-t153 * t559 - t320 * t367) * qJD(3) + t590 * qJD(1) + t488;
t27 = -t201 * t295 - t224 * t286 + t355 * t81 + (t157 * t559 - t320 * t366 - t539) * qJD(3) + t441;
t20 = t153 * t224 - t157 * t223 - t292 * t523 + t295 * t83 + t296 * t81 + t544;
t15 = -t667 * t355 - t591 * t296 - t747 * t223 + (-t531 + t590) * qJD(1) + (t367 * t499 + t559 * t598) * qJD(3) + t488;
t14 = qJD(1) * t334 + t668 * t355 - t591 * t295 + t747 * t224 + (t366 * t499 + t559 * t597 - t539) * qJD(3) + t441;
t1 = [(t373 * t787 + t754 + t757 + t777) * t295 * t682 + ((t122 + t165) * t367 + (t121 + t164) * t366) * t554 / 0.2e1 - (t65 + t41 + t60) * t562 / 0.2e1 + (t61 + t64) * t563 / 0.2e1 + (-(-t44 + t407 + t768) * t45 - t504 * t296 + t15 * (t508 + t769) + t44 * (-t377 * t660 - t476) + t14 * (t717 + t724) + t45 * (t354 - t549 + t702) + (-t14 * t606 + t45 * (-t605 - t615) * qJD(3) + (t15 * t372 + (t375 * t44 - t45 * t604) * qJD(4)) * pkin(4) + t44 * (t373 * t643 + t507) * qJD(1)) * t367 + (t15 * (-pkin(2) + t606 - t662) + t45 * (t507 - t662) * qJD(1) + (t498 - t557 + (t376 * t643 + t615) * qJD(3) + (-pkin(6) - t674) * qJD(1)) * t44) * t366) * m(6) + m(3) * ((-t309 * t378 - t551) * t716 + (-t550 + (-0.2e1 * t502 - t369 + t716) * t378) * (-t309 - t676)) + (-qJD(3) * t442 + t315 * t376 + t316 * t373) * qJD(1) + t647 / 0.2e1 + t644 / 0.2e1 + t645 / 0.2e1 + t646 / 0.2e1 + t784 * t685 + (t783 + t733) * t684 + t733 * t683 + (t28 * (-t478 + t508) + t62 * (t325 - t479) + t27 * (t717 - t596) + t63 * (-pkin(3) * t533 + t326 + t354 + t545) + (t560 * t62 * t677 + t28 * t425) * t366 + ((-t374 * t63 - t377 * t62) * pkin(1) + t425 * t648 + (-t62 * pkin(6) + t63 * (-pkin(2) - t348 - t664)) * t366) * qJD(1) - (t407 - t62 + t745) * t63) * m(5) + (t87 * (t366 * t526 + t508 + t570) + t86 * t749 + t114 * (t354 + t573) + (t336 * t624 - t622) * qJD(3) + ((-t113 * t377 - t114 * t374) * pkin(1) + (-pkin(2) - t480) * t623 + (t113 * (-rSges(4,3) - pkin(6)) + t114 * t526) * t366) * qJD(1) - (-t534 - t113 - t300 + (-t232 - t676) * qJD(1)) * t114) * m(4) + t657 / 0.2e1 + t697 - t658 / 0.2e1 + t655 / 0.2e1 - t656 / 0.2e1 + (t131 + ((t93 - t202 + (t227 + t620) * t367 + t588) * t367 + t587 * t366) * qJD(3)) * t515 + t731 * t688 + t730 * t687 + (t40 - t555 + ((t367 * t500 - t587 + t95) * t367 + (t366 * t500 + t501 + t94) * t366) * qJD(3)) * t518; m(4) * t55 + m(5) * t20 + m(6) * t5; qJD(1) * (t366 * t61 - t367 * t60 + (t121 * t366 + t122 * t367) * qJD(1)) / 0.2e1 + ((-t562 * t617 - t568) * t367 + (t394 + (t408 * t366 + (t616 - t691) * t367) * qJD(3)) * t366) * t515 + ((-t563 * t616 + t568) * t366 + (t394 + (-t691 * t367 + (t617 + t408) * t366) * qJD(3)) * t367) * t518 - qJD(1) * ((t373 * t572 + t376 * t571) * qJD(1) + ((t366 * t583 - t367 * t584) * t376 + (t366 * t585 + t367 * t586) * t373) * qJD(3)) / 0.2e1 + t705 * t688 + t707 * t687 + ((t730 * t373 + t611 * t779) * qJD(4) + ((qJD(4) * t778 + t711) * t376 + t715) * t367 + (t263 * t719 + t264 * t718) * t355 + (t263 * t723 + t264 * t721) * t296 + (-t263 * t722 - t264 * t720) * t295) * t686 + (qJD(1) * t708 + t366 * t740 - t367 * t741) * t685 + (qJD(1) * t706 + t366 * t737 - t367 * t738) * t684 + ((t731 * t373 + t607 * t780) * qJD(4) + ((qJD(4) * t781 + t711) * t376 + t715) * t366 + (-t261 * t719 + t262 * t718) * t355 + (-t261 * t723 + t262 * t721) * t296 + (t261 * t722 - t262 * t720) * t295) * t683 + ((t710 * qJD(4) - t785) * t376 + ((-t372 * t719 + t375 * t718 + t772) * t355 + (-t372 * t723 + t375 * t721 - t793) * t296 + (t372 * t722 - t375 * t720 + t773) * t295 + t729 * qJD(4)) * t373) * t682 + (qJD(1) * t710 + t366 * t735 - t367 * t736) * t681 - t732 * t559 / 0.2e1 + t709 * t484 + (t5 * t580 + (-t506 * qJD(1) + t15 * t541 + t739) * t367 + (-t505 * qJD(1) + t14 * t541 + t417) * t366 - (t373 * t410 + t376 * t381) * qJD(4) + (t295 * t748 - t594 * t355 + t366 * t701 + t541 * t566 - t483) * t45 + (t668 * t367 + (-t292 - t597) * t567 - t557 - t594 * t296 - t595 * t295 + t700) * t37 + (t296 * t748 + t595 * t355 + t367 * t701 + t299 - t429) * t44) * m(6) + (-t62 * (t218 * t355 - t288 * t296 + t429) - t63 * (-t220 * t355 - t288 * t295 + t483) - ((-t153 * t62 + t157 * t63) * t373 + t385 * t376) * qJD(4) + t62 * t299 + t20 * t580 + (t20 * t157 + t62 * t589 + (qJD(1) * t63 + t28) * t575) * t367 + (t62 * t286 * qJD(1) + t20 * t153 + t27 * t575 + t63 * t589) * t366 + (t218 * t295 + t220 * t296 + (t153 * qJD(1) + t81) * t367 + (t596 * qJD(1) + t83) * t366 + t700) * t54) * m(5) + (-(t113 * t283 - t622) * qJD(1) - (t112 * (-t283 * t366 - t284 * t367) + t453 * t480) * qJD(3) + t55 * t445 + t112 * t388 + t453 * t318 + (-t86 * t366 - t87 * t367 + (-t114 * t367 + t624) * qJD(1)) * t336) * m(4) + (qJD(1) * t64 + (-t693 * t750 + t694 * t366 ^ 2 + (t94 * t366 + t95 * t367) * qJD(1)) * t744 + t743) * t366 / 0.2e1 - (qJD(1) * t65 + (t693 * t367 ^ 2 - t694 * t750 + (t92 * t366 + t93 * t367) * qJD(1)) * t744 + t742) * t367 / 0.2e1 + (t427 + t40 + t734) * t567 / 0.2e1 + (t426 + t41 + t733) * t519 - (t366 * t734 + t367 * t733) * t558 / 0.2e1; (t373 * t706 - t376 * t731) * t688 + (t373 * t708 - t376 * t730) * t687 + (t263 * t786 + t713 * t264 - t712 * t367) * t686 + ((qJD(3) * t708 - t784) * t376 + (-qJD(1) * t707 + qJD(3) * t730 + t366 * t741 + t367 * t740) * t373) * t685 + ((qJD(3) * t706 - t783) * t376 + (-qJD(1) * t705 + qJD(3) * t731 + t366 * t738 + t367 * t737) * t373) * t684 + (-t261 * t786 + t262 * t713 - t366 * t712) * t683 + (t767 * t376 + (-t372 * t786 + t375 * t713) * t373) * t682 + ((qJD(3) * t710 - t782) * t376 + (-qJD(1) * t709 + qJD(3) * t729 + t366 * t736 + t367 * t735) * t373) * t681 - (t646 + t647 + t657 - t658 + t644 + t645 + t655 - t656 + t697) * t376 / 0.2e1 + t742 * t613 / 0.2e1 + t743 * t609 / 0.2e1 + t732 * t514 + (t373 * t710 - t376 * t729) * t484 + ((qJD(3) * t381 - t14 * t597 - t15 * t598 + t44 * t667 - t45 * t668) * t376 + (t410 * qJD(3) + (qJD(1) * t409 + t14 * t747 - t45 * t591 + t417) * t367 + (qJD(1) * t746 - t15 * t747 - t37 * t668 + t44 * t591 - t739) * t366) * t373 - (-t44 * t593 + t45 * t592) * t355 - (t37 * t592 + t44 * t487) * t296 - (t37 * t593 + t45 * t487) * t295) * m(6) + (-t62 * (-t188 * t355 - t296 * t308) - t63 * (t190 * t355 - t295 * t308) - t54 * (t188 * t295 + t190 * t296) + (qJD(3) * t385 + t28 * t153 - t27 * t157 + t62 * t83 - t63 * t81) * t376 + (t62 * (-qJD(3) * t153 + t201 * t366) + t63 * (qJD(3) * t157 - t201 * t367) + t20 * t448 + t54 * (-t153 * t567 - t157 * t566 - t366 * t81 + t367 * t83) + (-t27 * t367 + t28 * t366 + (t366 * t63 + t648) * qJD(1)) * t286) * t373) * m(5) + t733 * (t367 * t513 - t538 / 0.2e1) + t734 * (t366 * t513 + t373 * t519); 0.2e1 * ((t44 * t562 + t45 * t563 - t5) * t689 + (t295 * t45 + t296 * t44) * t690) * t376 + 0.2e1 * ((qJD(3) * t37 + t14 * t366 + t15 * t367 - t44 * t567 + t45 * t566) * t689 + (t37 * (t295 * t366 + t296 * t367) + (-t366 * t44 + t367 * t45) * t355) * t690) * t373;];
tauc = t1(:);
