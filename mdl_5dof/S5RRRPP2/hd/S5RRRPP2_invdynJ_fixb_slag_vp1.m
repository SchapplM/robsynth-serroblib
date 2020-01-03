% Calculate vector of inverse dynamics joint torques for
% S5RRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:30
% EndTime: 2019-12-31 20:51:57
% DurationCPUTime: 23.68s
% Computational Cost: add. (14488->664), mult. (15923->761), div. (0->0), fcn. (12369->6), ass. (0->370)
t411 = sin(qJ(3));
t400 = Icges(6,4) * t411;
t413 = cos(qJ(3));
t498 = Icges(6,1) * t413 + t400;
t399 = Icges(5,5) * t411;
t499 = Icges(5,1) * t413 + t399;
t794 = t498 + t499;
t793 = Icges(5,1) + Icges(6,1);
t792 = Icges(5,4) - Icges(6,5);
t791 = Icges(6,2) + Icges(5,3);
t790 = Icges(5,6) - Icges(6,6);
t410 = qJ(1) + qJ(2);
t396 = sin(t410);
t640 = Icges(4,4) * t411;
t335 = Icges(4,1) * t413 - t640;
t397 = cos(t410);
t469 = t335 * t397;
t213 = Icges(4,5) * t396 + t469;
t788 = t794 * t397;
t782 = t792 * t396 + t788;
t772 = t213 + t782;
t609 = t397 * t413;
t789 = (Icges(6,4) + Icges(5,5)) * t609;
t787 = Icges(4,5) + t792;
t786 = -Icges(4,6) + t790;
t638 = Icges(5,5) * t413;
t321 = Icges(5,3) * t411 + t638;
t198 = -Icges(5,6) * t397 + t321 * t396;
t639 = Icges(6,4) * t413;
t325 = Icges(6,2) * t411 + t639;
t202 = Icges(6,6) * t397 + t325 * t396;
t785 = -t198 - t202;
t610 = t397 * t411;
t779 = t790 * t396 + t791 * t610 + t789;
t323 = Icges(4,5) * t413 - Icges(4,6) * t411;
t464 = t323 * t397;
t201 = Icges(4,3) * t396 + t464;
t327 = Icges(5,4) * t413 + Icges(5,6) * t411;
t465 = t327 * t397;
t205 = Icges(5,2) * t396 + t465;
t784 = t201 + t205;
t752 = -t791 * t413 + t399 + t400;
t766 = -Icges(4,2) * t413 - t640 + t752;
t401 = Icges(4,4) * t413;
t497 = -Icges(4,2) * t411 + t401;
t778 = -t321 - t325;
t770 = t497 + t778;
t208 = Icges(6,5) * t397 + t396 * t498;
t210 = -Icges(5,4) * t397 + t396 * t499;
t613 = t396 * t411;
t363 = Icges(4,4) * t613;
t612 = t396 * t413;
t212 = Icges(4,1) * t612 - Icges(4,5) * t397 - t363;
t773 = -t208 - t210 - t212;
t334 = Icges(4,1) * t411 + t401;
t751 = t793 * t411 + t334 - t638 - t639;
t781 = -t787 * t411 + t786 * t413;
t780 = t335 + t794;
t200 = Icges(4,5) * t612 - Icges(4,6) * t613 - Icges(4,3) * t397;
t206 = Icges(4,4) * t612 - Icges(4,2) * t613 - Icges(4,6) * t397;
t626 = t206 * t411;
t485 = -t212 * t413 + t626;
t488 = t202 * t411 + t208 * t413;
t204 = -Icges(5,2) * t397 + t327 * t396;
t627 = t204 * t397;
t319 = Icges(6,5) * t413 + Icges(6,6) * t411;
t196 = Icges(6,3) * t397 + t319 * t396;
t629 = t196 * t397;
t492 = t198 * t411 + t210 * t413;
t708 = t396 * t492;
t726 = t396 * t488 - t627 + t629 + t708;
t692 = -t200 * t397 - t396 * t485 + t726;
t463 = t319 * t397;
t197 = -Icges(6,3) * t396 + t463;
t190 = t397 * t197;
t725 = -t205 * t397 + t782 * t612 + t779 * t613 + t190;
t466 = t497 * t397;
t207 = Icges(4,6) * t396 + t466;
t172 = t213 * t612;
t526 = t201 * t397 - t172;
t80 = -t207 * t613 - t526;
t691 = t80 + t725;
t777 = t784 * t396 + t772 * t609 + t779 * t610;
t188 = t396 * t204;
t776 = -t396 * t200 + t773 * t609 + t785 * t610 - t188;
t775 = -t206 - t785;
t774 = t207 - t779;
t409 = qJD(1) + qJD(2);
t768 = -qJD(3) * t766 + t786 * t409;
t767 = -t751 * qJD(3) + t787 * t409;
t763 = t411 * t766 + t751 * t413;
t630 = t196 * t396;
t690 = -t206 * t610 - t630 - t776;
t689 = -t197 * t396 - t207 * t610 + t777;
t765 = t770 * qJD(3);
t764 = t780 * qJD(3);
t757 = t319 - t323 - t327;
t762 = -t751 * t411 + t766 * t413;
t684 = t781 * t397;
t685 = t781 * t396;
t614 = t396 * t409;
t761 = t768 * t397 + t770 * t614;
t611 = t397 * t409;
t760 = -t768 * t396 + t409 * t466 + t778 * t611;
t759 = t767 * t397 - t614 * t780;
t758 = (t469 + t788) * t409 + t767 * t396;
t713 = t763 * t396 + t684;
t712 = t763 * t397 - t685;
t756 = (-Icges(5,2) - Icges(4,3) - Icges(6,3)) * t409 - t781 * qJD(3);
t625 = t207 * t411;
t755 = t779 * t411 + t772 * t413 - t625;
t754 = -t485 + t488 + t492;
t687 = t772 * t411 + t774 * t413;
t688 = t773 * t411 + t775 * t413;
t389 = t397 * pkin(7);
t286 = pkin(2) * t396 - t389;
t340 = pkin(7) * t611;
t753 = t409 * t286 + t340;
t750 = t762 * qJD(3) - t409 * t781 - t765 * t411 + t764 * t413;
t749 = t689 * t396 - t690 * t397;
t748 = t691 * t396 - t692 * t397;
t747 = t757 * qJD(3) + t763 * t409;
t746 = rSges(6,1) + pkin(4);
t745 = rSges(6,3) + qJ(5);
t398 = t411 * qJ(4);
t699 = t413 * pkin(3) + t398;
t264 = t699 * t396;
t232 = t409 * t264;
t744 = t232 + t753;
t743 = t712 * t409;
t742 = -t758 * t413 + t760 * t411 + (t196 - t200 - t204) * t409 - t688 * qJD(3);
t741 = t759 * t413 + t761 * t411 + (-t197 + t784) * t409 - t687 * qJD(3);
t740 = t751 + t770;
t739 = t766 + t780;
t738 = t766 * t397 + t772;
t737 = -t334 * t397 - t793 * t610 - t774 + t789;
t736 = t713 * t409;
t735 = (-t465 - t464 + t463 + t754) * t409 + t756 * t396;
t734 = t756 * t397 + t755 * t409 - t757 * t614;
t585 = rSges(6,2) * t610 - t396 * t745 + t746 * t609;
t562 = qJD(3) * t409;
t272 = -qJDD(3) * t397 + t396 * t562;
t408 = qJDD(1) + qJDD(2);
t341 = pkin(3) * t411 - qJ(4) * t413;
t412 = sin(qJ(1));
t414 = cos(qJ(1));
t416 = qJD(1) ^ 2;
t462 = (-qJDD(1) * t412 - t414 * t416) * pkin(1);
t557 = qJD(3) * qJD(4);
t721 = qJDD(4) * t411 + t413 * t557;
t433 = t272 * t341 + t721 * t397 + t462;
t559 = qJD(4) * t413;
t277 = qJD(3) * t699 - t559;
t403 = t411 * rSges(6,2);
t700 = t413 * rSges(6,1) + t403;
t580 = -qJD(3) * t700 - t277;
t651 = pkin(4) * t413;
t454 = (-t651 * qJD(3) + t580) * qJD(3);
t395 = qJD(4) * t411;
t536 = t396 * t395;
t471 = qJD(5) * t397 + t536;
t581 = -t264 - t286;
t220 = rSges(6,3) * t397 + t396 * t700;
t282 = pkin(4) * t612 + qJ(5) * t397;
t587 = t220 + t282;
t519 = t581 - t587;
t342 = rSges(6,1) * t411 - rSges(6,2) * t413;
t652 = pkin(4) * t411;
t527 = t342 + t652;
t560 = qJD(3) * t413;
t540 = t396 * t560;
t608 = t409 * t411;
t235 = t397 * t608 + t540;
t561 = qJD(3) * t411;
t541 = t396 * t561;
t315 = pkin(4) * t541;
t607 = t409 * t413;
t550 = t397 * t607;
t703 = rSges(6,1) * t541 + t745 * t614;
t605 = rSges(6,1) * t550 + rSges(6,2) * t235 - t315 + (pkin(4) * t607 + qJD(5)) * t397 - t703;
t316 = pkin(3) * t541;
t512 = -t316 + t536;
t147 = pkin(3) * t550 + qJ(4) * t235 + t512;
t287 = t397 * pkin(2) + t396 * pkin(7);
t240 = t287 * t409;
t706 = -t147 - t240;
t8 = -qJDD(5) * t396 + t527 * t272 + t454 * t397 + t519 * t408 + (-t471 - t605 + t706) * t409 + t433;
t710 = t8 - g(1);
t271 = qJDD(3) * t396 + t397 * t562;
t538 = t397 * t561;
t456 = -t396 * t607 - t538;
t552 = t396 * t608;
t537 = t397 * t560;
t308 = qJ(4) * t537;
t354 = t397 * t395;
t578 = t308 + t354;
t146 = pkin(3) * t456 - qJ(4) * t552 + t578;
t269 = pkin(3) * t609 + t397 * t398;
t407 = t414 * pkin(1);
t653 = pkin(1) * t412;
t516 = qJDD(1) * t407 - t416 * t653;
t472 = t409 * (-pkin(2) * t614 + t340) + t408 * t287 + t516;
t432 = t408 * t269 + t472 + (t146 + t354) * t409 + t721 * t396;
t515 = -t341 - t527;
t314 = rSges(6,2) * t537;
t558 = qJD(5) * t396;
t707 = t409 * t220;
t606 = -rSges(6,1) * t538 + pkin(4) * t456 - qJ(5) * t611 + t314 - t558 - t707;
t9 = qJDD(5) * t397 + t606 * t409 + t585 * t408 + t515 * t271 + (-qJD(5) * t409 + t454) * t396 + t432;
t709 = t9 - g(2);
t310 = rSges(5,1) * t541;
t383 = t396 * rSges(5,2);
t510 = t413 * rSges(5,1) + t411 * rSges(5,3);
t144 = rSges(5,3) * t540 - t310 + (t397 * t510 + t383) * t409;
t343 = rSges(5,1) * t411 - rSges(5,3) * t413;
t579 = -t510 * qJD(3) - t277;
t524 = qJD(3) * t579;
t385 = t397 * rSges(5,2);
t221 = t396 * t510 - t385;
t547 = -t221 + t581;
t23 = t272 * t343 + t397 * t524 + t547 * t408 + (-t144 - t536 + t706) * t409 + t433;
t733 = t23 - g(1);
t576 = rSges(5,2) * t611 + rSges(5,3) * t537;
t141 = rSges(5,1) * t456 - rSges(5,3) * t552 + t576;
t224 = rSges(5,1) * t609 + rSges(5,3) * t610 + t383;
t569 = -t341 - t343;
t24 = t141 * t409 + t224 * t408 + t271 * t569 + t396 * t524 + t432;
t732 = t24 - g(2);
t544 = rSges(4,1) * t541 + t235 * rSges(4,2);
t702 = rSges(4,1) * t609 + t396 * rSges(4,3);
t145 = t702 * t409 - t544;
t649 = rSges(4,1) * t413;
t349 = -rSges(4,2) * t411 + t649;
t302 = t349 * qJD(3);
t344 = rSges(4,1) * t411 + rSges(4,2) * t413;
t563 = qJD(3) * t397;
t567 = rSges(4,2) * t613 + t397 * rSges(4,3);
t222 = rSges(4,1) * t612 - t567;
t586 = -t222 - t286;
t39 = -t302 * t563 + t272 * t344 + (-t145 - t240) * t409 + t586 * t408 + t462;
t731 = t39 - g(1);
t236 = t537 - t552;
t470 = -t236 * rSges(4,2) + rSges(4,3) * t611;
t142 = rSges(4,1) * t456 + t470;
t225 = -rSges(4,2) * t610 + t702;
t564 = qJD(3) * t396;
t40 = t142 * t409 + t225 * t408 - t271 * t344 - t302 * t564 + t472;
t730 = t40 - g(2);
t238 = rSges(3,1) * t611 - rSges(3,2) * t614;
t284 = rSges(3,1) * t396 + rSges(3,2) * t397;
t729 = -t238 * t409 - t284 * t408 - g(1) + t462;
t285 = t397 * rSges(3,1) - rSges(3,2) * t396;
t624 = t284 * t409;
t728 = t285 * t408 - t409 * t624 - g(2) + t516;
t555 = -pkin(3) - t746;
t727 = t555 * t413 - pkin(2);
t460 = t515 * t563;
t522 = t397 * t555;
t646 = -rSges(6,2) - qJ(4);
t523 = t354 - t558;
t647 = pkin(1) * qJD(1);
t553 = t412 * t647;
t457 = t523 - t553;
t66 = t409 * t519 + t457 + t460;
t554 = t414 * t647;
t278 = t341 * t564;
t455 = -t315 + t471;
t546 = -t269 - t585;
t666 = t342 * t564 - t409 * (t287 - t546) + t278 - t455;
t67 = t554 - t666;
t698 = -t398 - t403 + t727;
t724 = ((-t66 * pkin(7) + t698 * t67) * t396 + (t698 * t66 - t67 * t745) * t397) * t409 + (t411 * t522 * t67 + t612 * t646 * t66) * qJD(3) - t460 * t67;
t193 = t409 * t222;
t723 = -rSges(4,1) * t538 + t193 + t470 + t753;
t722 = t409 * t221 + t576 + t578 + t744;
t227 = -t553 - t624;
t720 = t409 * t282 + t308 + t314 - t523 + t707 + t744;
t568 = -t699 - t510;
t719 = t748 * qJD(3) + t736;
t718 = t749 * qJD(3) + t743;
t717 = t754 * qJD(3) + t758 * t411 + t760 * t413;
t716 = t755 * qJD(3) + t759 * t411 - t761 * t413;
t715 = -t747 * t396 + t750 * t397;
t714 = t750 * t396 + t747 * t397;
t711 = t627 + t777;
t518 = t269 + t287;
t697 = t518 + t224;
t696 = t735 * t396 + t742 * t397;
t695 = -t734 * t396 + t741 * t397;
t694 = t742 * t396 - t735 * t397;
t693 = t741 * t396 + t734 * t397;
t681 = (-t740 * t411 + t739 * t413) * t409;
t680 = -t738 * t411 + t737 * t413;
t679 = -Icges(4,2) * t612 + t752 * t396 - t363 - t773;
t678 = t751 * t396 - t775;
t677 = t757 * t409;
t166 = t225 + t287;
t670 = -t166 * t409 + t344 * t564;
t665 = -t343 * t564 + t409 * t697 - t278;
t664 = t679 * t411 + t678 * t413;
t663 = m(5) / 0.2e1;
t662 = m(6) / 0.2e1;
t661 = t271 / 0.2e1;
t660 = t272 / 0.2e1;
t654 = -rSges(5,1) - pkin(3);
t650 = g(2) * t396;
t645 = -rSges(5,3) - qJ(4);
t542 = t344 * t563;
t459 = -t542 - t553;
t103 = t409 * t586 + t459;
t634 = t103 * t397;
t584 = -t224 - t269;
t583 = t396 * t264 + t397 * t269;
t358 = qJ(4) * t609;
t265 = -pkin(3) * t610 + t358;
t582 = t409 * t265 + t396 * t559;
t565 = t396 ^ 2 + t397 ^ 2;
t549 = t396 * t147 + (t146 + t232) * t397;
t356 = qJ(4) * t612;
t260 = -pkin(3) * t613 + t356;
t548 = t260 * t564 + t265 * t563 + t395;
t543 = t397 * t654;
t533 = -pkin(2) - t649;
t531 = -t564 / 0.2e1;
t530 = t564 / 0.2e1;
t529 = -t563 / 0.2e1;
t528 = t563 / 0.2e1;
t525 = -t200 + t625;
t511 = t264 * t564 + t269 * t563 - t559;
t350 = rSges(2,1) * t414 - rSges(2,2) * t412;
t345 = rSges(2,1) * t412 + rSges(2,2) * t414;
t104 = t554 - t670;
t494 = -t104 * t396 - t634;
t483 = t222 * t396 + t225 * t397;
t476 = -pkin(4) * t560 + t580;
t473 = t563 * t569 + t354;
t461 = -qJDD(4) * t413 + t146 * t563 + t147 * t564 + t271 * t264 + t411 * t557;
t458 = t536 + t554;
t165 = t396 * t533 + t389 + t567;
t449 = t411 * t645 + t413 * t654 - pkin(2);
t446 = t473 - t553;
t102 = t518 + t585;
t431 = t316 - t455 + t703;
t151 = t396 * t449 + t385 + t389;
t101 = t389 - t745 * t397 + (t411 * t646 + t727) * t396;
t420 = (t533 * t634 + (t103 * (-rSges(4,3) - pkin(7)) + t104 * t533) * t396) * t409;
t419 = (((t80 - t172 + (t201 + t626) * t397 + t776) * t397 + (-t708 + (-t197 - t488) * t396 + t711 + t726) * t396) * qJD(3) + t743) * t528 + (t763 * qJD(3) + t764 * t411 + t765 * t413) * t409 + (Icges(3,3) - t762) * t408 + (t687 + t712) * t661 + (-t688 + t713) * t660 + (((t397 * t525 + t629 + t689 - t711) * t397 + (t396 * t525 - t188 + t190 + t526 + t630 + t690 - t725) * t396) * qJD(3) + t719 - t736) * t531 + (t715 + t716) * t530 + (t714 + t717 + t718) * t529;
t73 = t409 * t547 + t446;
t74 = t458 + t665;
t418 = (t73 * t449 * t397 + (t73 * (-rSges(5,2) - pkin(7)) + t74 * (-pkin(2) + t568)) * t396) * t409 + (t411 * t543 * t74 + t612 * t645 * t73) * qJD(3);
t370 = rSges(6,2) * t609;
t369 = rSges(5,3) * t609;
t366 = rSges(6,2) * t612;
t365 = rSges(5,3) * t612;
t355 = t397 * t559;
t270 = t341 * t614;
t268 = t344 * t397;
t267 = -rSges(5,1) * t610 + t369;
t266 = -rSges(6,1) * t610 + t370;
t263 = t344 * t396;
t262 = -rSges(5,1) * t613 + t365;
t261 = -rSges(6,1) * t613 + t366;
t241 = t565 * t561;
t228 = t285 * t409 + t554;
t108 = t483 * qJD(3);
t72 = (t221 * t396 + t224 * t397) * qJD(3) + t511;
t59 = (t587 * t396 + t585 * t397) * qJD(3) + t511;
t10 = t221 * t271 + t584 * t272 + (t141 * t397 + t144 * t396) * qJD(3) + t461;
t7 = t587 * t271 + t546 * t272 + (t605 * t396 + t606 * t397) * qJD(3) + t461;
t1 = [Icges(2,3) * qJDD(1) + t419 + (t728 * (t285 + t407) + t729 * (-t284 - t653) + (-t238 - t554 + t228) * t227) * m(3) + ((t345 ^ 2 + t350 ^ 2) * qJDD(1) + g(1) * t345 - g(2) * t350) * m(2) + (t66 * (t431 - t554) + t709 * (t407 + t102) + t710 * (t101 - t653) + (t66 + t553 + t457 + t720) * t67 + t724) * m(6) + (t73 * (t310 + t316 - t458) + t418 + (-t553 + t73 - t446 + t722) * t74 + t732 * (t407 + t697) + t733 * (t151 - t653)) * m(5) + (t103 * (t544 - t554) + t420 + t730 * (t166 + t407) + t731 * (t165 - t653) + (-t553 + t103 - t459 + t723) * t104) * m(4); t419 + ((t523 + t720) * t67 + (-t666 + t431) * t66 + t709 * t102 + t710 * t101 + t724) * m(6) + (t418 + (-t473 + t722) * t74 + (t310 - t512 + t536 + t665) * t73 + t732 * t697 + t733 * t151) * m(5) + (t420 + t730 * t166 + t731 * t165 + (t542 + t723) * t104 + (-t670 + t544) * t103) * m(4) + (-t227 * t238 - t228 * t624 + (t227 * t409 + t728) * t285 + (t228 * t409 - t729) * t284) * m(3); t749 * t661 + t748 * t660 + (t715 * t409 + t712 * t408 + t690 * t272 + t689 * t271 + (t695 * t396 + t696 * t397) * qJD(3)) * t396 / 0.2e1 - (t714 * t409 + t713 * t408 + t692 * t272 + t691 * t271 + (t693 * t396 + t694 * t397) * qJD(3)) * t397 / 0.2e1 + (t687 * t396 + t688 * t397) * t408 / 0.2e1 - ((t739 * t411 + t740 * t413) * t409 + ((t738 * t396 - t679 * t397) * t413 + (t737 * t396 + t678 * t397) * t411) * qJD(3)) * t409 / 0.2e1 + ((t687 * t409 - t717) * t397 + (-t688 * t409 + t716) * t396) * t409 / 0.2e1 + t719 * t614 / 0.2e1 + t718 * t611 / 0.2e1 + ((t684 * t564 - t677) * t396 + ((t664 * t397 + (t680 - t685) * t396) * qJD(3) + t681) * t397) * t531 + ((t689 * t409 + t696) * t397 + (t690 * t409 + t695) * t396) * t530 + ((t691 * t409 + t694) * t397 + (t692 * t409 + t693) * t396) * t529 + ((t685 * t563 + t677) * t397 + ((t680 * t396 + (t664 - t684) * t397) * qJD(3) + t681) * t396) * t528 + (-g(1) * (t358 + t370) - g(2) * (t356 + t366) - (g(1) * t522 + t555 * t650) * t411 + t7 * t583 + (t9 * t515 + t7 * t587) * t396 + (t8 * t515 + t7 * t585) * t397 + (t476 * t396 - t582 + (pkin(4) * t610 + t515 * t397 - t266) * t409) * t67 + (t476 * t397 + t270 - t355 + (t342 * t396 + t260 + t261) * t409) * t66 + (-t548 - (t261 * t396 + t266 * t397 - t565 * t652) * qJD(3) + t549 + (t546 * t409 + t605) * t396 + (t587 * t409 + t606) * t397) * t59 + (-g(3) + (t67 * t396 + t66 * t397) * qJD(3)) * (t699 + t700 + t651)) * m(6) + (-g(1) * (t358 + t369) - g(2) * (t356 + t365) + g(3) * t568 - (g(1) * t543 + t650 * t654) * t411 - t73 * (t355 + (-t260 - t262) * t409) - t74 * (t267 * t409 + t582) - t72 * t548 - ((t72 * t267 + t568 * t73) * t397 + (t72 * t262 + t568 * t74) * t396) * qJD(3) + t73 * t270 + t10 * t583 + t72 * t549 + (t23 * t569 + t73 * t579 + t10 * t224 + t72 * t141 + (t72 * t221 + t569 * t74) * t409) * t397 + (t24 * t569 + t74 * t579 + t10 * t221 + t72 * t144 + (t73 * t343 + t584 * t72) * t409) * t396) * m(5) + ((t222 * t271 - t225 * t272 + (t142 * t397 + t145 * t396) * qJD(3)) * t483 + t108 * ((t142 + t193) * t397 + (-t225 * t409 + t145) * t396) + t494 * t302 + ((-t104 * t409 - t39) * t397 + (t103 * t409 - t40) * t396) * t344 - (t103 * t263 - t104 * t268) * t409 - (t108 * (-t263 * t396 - t268 * t397) + t494 * t349) * qJD(3) + g(1) * t268 + g(2) * t263 - g(3) * t349) * m(4); (-m(5) - m(6)) * (-g(3) * t413 + (g(1) * t397 + t650) * t411) - m(5) * (t235 * t74 + t236 * t73 + t241 * t72) - m(6) * (t235 * t67 + t236 * t66 + t241 * t59) + 0.2e1 * ((t563 * t73 + t564 * t74 - t10) * t663 + (t563 * t66 + t564 * t67 - t7) * t662) * t413 + 0.2e1 * ((qJD(3) * t72 + t23 * t397 + t24 * t396 + t611 * t74 - t614 * t73) * t663 + (qJD(3) * t59 + t396 * t9 + t397 * t8 + t611 * t67 - t614 * t66) * t662) * t411; (-t710 * t396 + t709 * t397) * m(6);];
tau = t1;
