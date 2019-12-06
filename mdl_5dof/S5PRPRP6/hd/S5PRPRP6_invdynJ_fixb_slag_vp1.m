% Calculate vector of inverse dynamics joint torques for
% S5PRPRP6
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP6_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP6_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP6_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP6_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP6_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:21
% EndTime: 2019-12-05 15:41:17
% DurationCPUTime: 47.88s
% Computational Cost: add. (9008->676), mult. (24768->978), div. (0->0), fcn. (23575->6), ass. (0->317)
t725 = Icges(6,4) + Icges(5,5);
t723 = Icges(5,6) - Icges(6,6);
t753 = Icges(5,4) - Icges(6,5);
t754 = Icges(5,1) + Icges(6,1);
t724 = Icges(5,2) + Icges(6,3);
t765 = Icges(6,2) + Icges(5,3);
t346 = sin(qJ(4));
t348 = cos(qJ(4));
t768 = t346 * t725 + t348 * t723;
t767 = t753 * t348;
t766 = t753 * t346;
t347 = sin(qJ(2));
t349 = cos(qJ(2));
t764 = t347 * t768 + t349 * t765;
t763 = (t346 * t723 - t348 * t725) * t349;
t762 = -t348 * t724 - t766;
t761 = t346 * t754 + t767;
t749 = -Icges(3,5) + Icges(4,4);
t747 = -Icges(3,6) + Icges(4,5);
t345 = cos(pkin(7));
t344 = sin(pkin(7));
t565 = t346 * t347;
t498 = t344 * t565;
t260 = t345 * t348 - t498;
t760 = t753 * t260;
t563 = t347 * t348;
t411 = -t344 * t346 + t345 * t563;
t759 = t753 * t411;
t259 = t344 * t563 + t345 * t346;
t758 = t753 * t259;
t258 = t344 * t348 + t345 * t565;
t757 = t753 * t258;
t515 = qJD(2) * t349;
t490 = t346 * t515;
t156 = qJD(4) * t411 + t345 * t490;
t562 = t348 * t349;
t495 = t345 * t562;
t157 = -qJD(2) * t495 + qJD(4) * t258;
t516 = qJD(2) * t347;
t491 = t345 * t516;
t756 = t156 * t725 - t157 * t723 - t491 * t765;
t158 = qJD(4) * t259 + t344 * t490;
t159 = -qJD(4) * t498 + (qJD(4) * t345 + t344 * t515) * t348;
t492 = t344 * t516;
t755 = t158 * t725 + t159 * t723 - t492 * t765;
t566 = t345 * t349;
t741 = -t411 * t724 - t566 * t723 - t757;
t568 = t344 * t349;
t740 = -t259 * t724 - t568 * t723 + t760;
t709 = t258 * t725 + t411 * t723 + t566 * t765;
t708 = t259 * t723 - t260 * t725 + t568 * t765;
t739 = t258 * t754 + t566 * t725 + t759;
t738 = -t260 * t754 + t568 * t725 + t758;
t752 = qJD(2) * t764 + qJD(4) * t763;
t661 = t347 * t762 - t349 * t723;
t660 = t347 * t761 + t349 * t725;
t751 = (-t346 * t724 + t767) * t349;
t750 = (-t348 * t754 + t766) * t349;
t704 = t347 * t723 + t349 * t762;
t735 = -t765 * t347 + t349 * t768;
t703 = -t347 * t725 + t349 * t761;
t748 = -Icges(3,2) - Icges(4,3);
t701 = t347 * t749 + t349 * t747;
t746 = (Icges(3,4) + Icges(4,6)) * t349;
t745 = -t156 * t753 + t157 * t724 + t491 * t723;
t744 = -t158 * t753 - t159 * t724 + t492 * t723;
t743 = t156 * t754 - t157 * t753 - t491 * t725;
t742 = t158 * t754 + t159 * t753 - t492 * t725;
t737 = qJD(2) * t661 + qJD(4) * t751;
t736 = qJD(2) * t660 + qJD(4) * t750;
t734 = t349 * t752 + t516 * t735;
t733 = t349 * t755 - t516 * t708;
t732 = t349 * t756 - t516 * t709;
t731 = -t346 * t703 + t348 * t704;
t730 = t346 * t738 - t348 * t740;
t729 = t346 * t739 - t348 * t741;
t728 = t735 * t349;
t727 = Icges(4,1) + Icges(3,3);
t722 = t747 * t345 + (t347 * t748 + t746) * t344;
t569 = t344 * t347;
t588 = Icges(3,4) * t347;
t721 = t344 * (Icges(3,1) * t349 - t588) + Icges(4,2) * t568 - Icges(4,6) * t569 + t749 * t345;
t720 = t701 * qJD(2);
t719 = t347 * t747 - t349 * t749;
t653 = t735 * t347;
t675 = t156 * t739 + t157 * t741 + t258 * t743 + t345 * t732 - t411 * t745;
t674 = t156 * t738 + t157 * t740 + t258 * t742 + t345 * t733 - t411 * t744;
t673 = t158 * t739 - t159 * t741 - t259 * t745 - t260 * t743 + t344 * t732;
t672 = t158 * t738 - t159 * t740 - t259 * t744 - t260 * t742 + t344 * t733;
t718 = -t156 * t703 - t157 * t704 + t258 * t736 + t345 * t734 - t411 * t737;
t717 = -t158 * t703 + t159 * t704 - t259 * t737 - t260 * t736 + t344 * t734;
t670 = t258 * t739 - t411 * t741 + t566 * t709;
t716 = t258 * t738 - t411 * t740 + t566 * t708;
t715 = -t259 * t741 - t260 * t739 + t568 * t709;
t698 = -t259 * t740 - t260 * t738 + t568 * t708;
t668 = t347 * t709 - t349 * t729;
t667 = t347 * t708 - t349 * t730;
t714 = -t258 * t703 - t345 * t728 + t411 * t704;
t713 = t259 * t704 + t260 * t703 - t344 * t728;
t640 = -t349 * t731 - t653;
t712 = t345 * t344;
t512 = qJD(4) * t349;
t520 = qJD(2) * t344;
t288 = t345 * t512 + t520;
t518 = qJD(2) * t345;
t289 = t344 * t512 - t518;
t705 = (t344 * t735 + t730) * t289 + (t345 * t735 + t729) * t288;
t702 = (t731 + t764) * t347;
t513 = qJD(4) * t347;
t700 = (t288 * t715 + t289 * t698 + t513 * t713) * t344 + (t288 * t670 + t289 * t716 + t513 * t714) * t345;
t627 = t288 * (-t411 * t754 - t741 + t757) + t289 * (-t259 * t754 - t740 - t760) + t513 * (t704 - t750);
t699 = t288 * (-t258 * t724 + t739 + t759) + t289 * (t260 * t724 + t738 + t758) + t513 * (-t703 - t751);
t343 = t345 ^ 2;
t676 = rSges(6,1) + pkin(4);
t307 = rSges(3,1) * t347 + rSges(3,2) * t349;
t342 = t344 ^ 2;
t629 = t342 + t343;
t697 = t307 * t629;
t696 = t344 * t719 - t345 * t727;
t695 = t344 * t727 + t345 * t719;
t694 = t720 * t344;
t693 = t720 * t345;
t691 = ((-Icges(4,6) * t347 + t349 * t748 - t588) * t516 + ((Icges(3,1) + Icges(4,2)) * t347 + t746) * t515) * t344 + (t347 * t721 + t349 * t722) * qJD(2);
t689 = t763 * t513 + (t259 * t725 + t260 * t723) * t289 + (-t258 * t723 + t411 * t725) * t288;
t687 = t347 * t722 - t349 * t721;
t680 = 0.2e1 * qJD(2);
t679 = 2 * qJDD(2);
t506 = qJD(2) * qJD(4);
t394 = qJDD(4) * t349 - t347 * t506;
t505 = qJDD(2) * t344;
t178 = t345 * t394 + t505;
t504 = qJDD(2) * t345;
t179 = t344 * t394 - t504;
t294 = qJDD(4) * t347 + t349 * t506;
t678 = t178 * t670 + t179 * t716 + t288 * t675 + t289 * t674 + t294 * t714 + t513 * t718;
t677 = t178 * t715 + t179 * t698 + t288 * t673 + t289 * t672 + t294 * t713 + t513 * t717;
t643 = (t745 * t348 - t743 * t346 + (-t346 * t741 - t348 * t739) * qJD(4) + t709 * qJD(2)) * t349 + (qJD(2) * t729 + t756) * t347;
t642 = (t744 * t348 - t742 * t346 + (-t346 * t740 - t348 * t738) * qJD(4) + t708 * qJD(2)) * t349 + (qJD(2) * t730 + t755) * t347;
t671 = t288 * t668 + t289 * t667 + t513 * t640;
t666 = rSges(6,3) + qJ(5);
t350 = qJD(2) ^ 2;
t635 = t350 * t629;
t665 = t704 * t344;
t664 = t704 * t345;
t663 = t703 * t344;
t662 = t703 * t345;
t659 = t701 * t344;
t658 = t701 * t345;
t657 = (qJD(4) * t702 + t705) * t349;
t655 = -t288 * t709 - t289 * t708;
t625 = t344 * t667 + t345 * t668;
t507 = qJD(2) * qJD(3);
t654 = qJDD(3) * t347 + t349 * t507;
t652 = t713 * t349;
t651 = t714 * t349;
t650 = t716 * t344;
t649 = t715 * t345;
t483 = -t516 / 0.2e1;
t641 = (t737 * t348 - t736 * t346 + (t346 * t704 + t348 * t703) * qJD(4) - t735 * qJD(2)) * t349 + (qJD(2) * t731 + t752) * t347;
t290 = pkin(3) * t344 + pkin(6) * t566;
t291 = -pkin(3) * t345 + pkin(6) * t568;
t634 = t290 * t504 + t291 * t505;
t564 = t346 * t349;
t499 = t344 * t564;
t633 = t676 * t499;
t496 = t345 * t564;
t632 = t676 * t496;
t339 = t349 * rSges(6,2);
t631 = t565 * t676 + t339;
t463 = -rSges(4,2) * t349 + rSges(4,3) * t347;
t630 = pkin(2) * t349 + qJ(3) * t347;
t626 = t689 * t349;
t624 = t344 * t698 + t649;
t623 = t345 * t670 + t650;
t622 = g(1) * t345 + g(2) * t344;
t621 = t622 * t347;
t509 = qJD(5) * t348;
t329 = t349 * t509;
t274 = t630 * t344;
t276 = t630 * t345;
t514 = qJD(3) * t349;
t416 = t274 * t520 + t276 * t518 + qJD(1) - t514;
t383 = t290 * t518 + t291 * t520 + t416;
t550 = rSges(6,2) * t568 - t259 * t666 - t260 * t676;
t551 = rSges(6,2) * t566 + t258 * t676 - t411 * t666;
t27 = t288 * t550 - t289 * t551 + t329 + t383;
t335 = qJD(3) * t347;
t317 = t344 * t335;
t305 = pkin(2) * t347 - qJ(3) * t349;
t404 = qJD(2) * t305;
t176 = -t344 * t404 + t317;
t319 = t345 * t335;
t177 = -t345 * t404 + t319;
t410 = t176 * t520 + t177 * t518 + t274 * t505 + t276 * t504 + t347 * t507 + qJDD(1);
t474 = pkin(6) * t635;
t510 = qJD(5) * t259;
t589 = -rSges(6,2) * t492 + t158 * t676 - t159 * t666 - t510;
t511 = qJD(5) * t411;
t601 = -rSges(6,2) * t491 + t156 * t676 + t157 * t666 - t511;
t5 = -t601 * t289 + t589 * t288 - t551 * t179 + t550 * t178 + (-qJD(4) * qJD(5) * t346 + qJDD(5) * t348 - qJDD(3)) * t349 + (-qJD(2) * t509 - t474) * t347 + t410 + t634;
t620 = t27 * t601 + t5 * t551;
t617 = -pkin(2) - pkin(6);
t616 = t178 / 0.2e1;
t615 = t179 / 0.2e1;
t614 = t629 * t483;
t613 = -t288 / 0.2e1;
t612 = t288 / 0.2e1;
t611 = -t289 / 0.2e1;
t610 = t289 / 0.2e1;
t609 = t294 / 0.2e1;
t605 = -t349 / 0.2e1;
t602 = g(3) * t349;
t599 = rSges(6,2) * t347;
t338 = t349 * rSges(5,3);
t567 = t345 * t347;
t286 = (-rSges(6,1) * t348 - rSges(6,3) * t346) * t349;
t387 = -t346 * t512 - t348 * t516;
t464 = rSges(6,1) * t346 - rSges(6,3) * t348;
t549 = t329 + t387 * qJ(5) + (t346 * t516 - t348 * t512) * pkin(4) + qJD(4) * t286 + (t347 * t464 + t339) * qJD(2);
t548 = t258 * t666 + t411 * t676;
t547 = t259 * t676 - t260 * t666;
t546 = t176 * t344 + t177 * t345;
t417 = -rSges(6,3) * t562 - t599;
t497 = t344 * t562;
t545 = -qJ(5) * t497 + t344 * t417 + t633;
t544 = -qJ(5) * t495 + t345 * t417 + t632;
t535 = -t563 * t666 + t631;
t534 = t599 + (-pkin(4) * t346 + qJ(5) * t348 - t464) * t349;
t533 = t274 * t344 + t276 * t345;
t253 = qJD(2) * t630 - t514;
t532 = -qJD(2) * t463 - t253;
t531 = (-pkin(4) * t348 - qJ(5) * t346) * t349 + t286;
t530 = rSges(5,1) * t499 + rSges(5,2) * t497;
t529 = rSges(5,1) * t496 + rSges(5,2) * t495;
t462 = rSges(4,2) * t347 + rSges(4,3) * t349;
t528 = -t305 + t462;
t527 = -t630 - t463;
t526 = t654 * t344;
t525 = t654 * t345;
t524 = rSges(4,2) * t569 + rSges(4,3) * t568;
t523 = rSges(4,2) * t567 + rSges(4,3) * t566;
t508 = -m(4) - m(5) - m(6);
t322 = qJ(3) * t566;
t468 = -pkin(2) * t567 + t322;
t321 = qJ(3) * t568;
t469 = -pkin(2) * t569 + t321;
t494 = t468 * t518 + t469 * t520 + t335;
t218 = rSges(5,1) * t565 + rSges(5,2) * t563 + t338;
t493 = pkin(6) * t349 + t630;
t481 = -t513 / 0.2e1;
t480 = t513 / 0.2e1;
t478 = -pkin(6) * t347 - t305;
t477 = t629 * t347;
t476 = qJD(2) * t528;
t475 = t290 * t345 + t291 * t344 + t533;
t465 = rSges(5,1) * t346 + rSges(5,2) * t348;
t220 = rSges(5,3) * t347 - t349 * t465;
t467 = -t220 + t478;
t466 = -pkin(6) * t515 - t253;
t310 = rSges(3,1) * t349 - rSges(3,2) * t347;
t39 = -t305 * t520 - t510 + t317 - t534 * t288 + (-pkin(6) * t520 + qJD(4) * t551) * t347;
t40 = -t305 * t518 - t511 + t319 + t534 * t289 + (-pkin(6) * t518 - qJD(4) * t550) * t347;
t460 = -t344 * t39 - t345 * t40;
t118 = rSges(5,1) * t258 + rSges(5,2) * t411 + rSges(5,3) * t566;
t452 = qJD(2) * t478;
t64 = t118 * t513 - t220 * t288 + t344 * t452 + t317;
t120 = -rSges(5,1) * t260 + rSges(5,2) * t259 + rSges(5,3) * t568;
t65 = -t120 * t513 + t220 * t289 + t345 * t452 + t319;
t453 = -t344 * t64 - t345 * t65;
t429 = t118 * t344 - t120 * t345;
t424 = t629 * t310;
t202 = rSges(4,1) * t344 + t345 * t463;
t203 = -rSges(4,1) * t345 + t344 * t463;
t423 = t202 * t345 + t203 * t344;
t420 = qJD(2) * t697;
t419 = t478 - t534;
t287 = (-rSges(5,1) * t348 + rSges(5,2) * t346) * t349;
t129 = qJD(4) * t287 + (t347 * t465 + t338) * qJD(2);
t418 = -t129 + t466;
t407 = t466 - t549;
t405 = qJD(2) * t462;
t41 = -t118 * t289 + t120 * t288 + t383;
t403 = t41 * t429;
t384 = t27 * t589 + t5 * t550;
t382 = t39 * t551 - t40 * t550;
t379 = qJD(2) * t532 + qJDD(2) * t528;
t376 = -qJDD(3) * t349 + t410;
t367 = -pkin(6) * t516 * t629 + t546;
t362 = -qJD(2) * t253 - qJDD(2) * t305 + (-qJDD(2) * t347 - t349 * t350) * pkin(6);
t355 = t344 * t362 + t526;
t354 = t345 * t362 + t525;
t353 = (-t27 * t550 + t39 * t534) * t345 + (t27 * t551 - t40 * t534) * t344;
t320 = t345 * t514;
t318 = t344 * t514;
t275 = t307 * t345;
t273 = t307 * t344;
t237 = t345 * t405;
t235 = t344 * t405;
t175 = -rSges(5,3) * t567 + t529;
t173 = -rSges(5,3) * t569 + t530;
t150 = rSges(5,1) * t259 + rSges(5,2) * t260;
t146 = rSges(5,1) * t411 - rSges(5,2) * t258;
t131 = t345 * t476 + t319;
t130 = t344 * t476 + t317;
t101 = rSges(5,1) * t158 + rSges(5,2) * t159 - rSges(5,3) * t492;
t99 = rSges(5,1) * t156 - rSges(5,2) * t157 - rSges(5,3) * t491;
t73 = t345 * t379 + t525;
t72 = t344 * t379 + t526;
t71 = qJD(2) * t423 + t416;
t66 = -qJD(2) * t420 + qJDD(2) * t424 + qJDD(1);
t36 = t423 * qJDD(2) + (t235 * t344 + t237 * t345) * qJD(2) + t376;
t29 = -t101 * t513 - t120 * t294 + t129 * t289 + t179 * t220 + t354;
t28 = t118 * t294 - t129 * t288 - t178 * t220 + t513 * t99 + t355;
t16 = t101 * t288 - t118 * t179 + t120 * t178 - t289 * t99 - t347 * t474 + t376 + t634;
t7 = qJD(5) * t157 - qJDD(5) * t411 + t179 * t534 + t289 * t549 - t294 * t550 - t513 * t589 + t354;
t6 = -qJD(5) * t159 - qJDD(5) * t259 - t178 * t534 - t288 * t549 + t294 * t551 + t513 * t601 + t355;
t1 = [m(2) * qJDD(1) + (-m(2) - m(3) + t508) * g(3) + m(3) * t66 + m(4) * t36 + m(5) * t16 + m(6) * t5; (t670 * t344 - t345 * t716) * t616 + (t344 * t715 - t345 * t698) * t615 + (((t258 * t660 - t411 * t661 - t650) * t347 + t651) * qJD(4) + (((t653 - t670) * qJD(4) + t655) * t347 + t657) * t345 + (t258 * t663 - t411 * t665) * t289 + (t258 * t662 - t411 * t664) * t288) * t613 + (t344 * t675 - t345 * t674) * t612 + (((-t259 * t661 - t260 * t660 - t649) * t347 + t652) * qJD(4) + (((t653 - t698) * qJD(4) + t655) * t347 + t657) * t344 + (-t259 * t665 - t260 * t663) * t289 + (-t259 * t664 - t260 * t662) * t288) * t611 + (t344 * t673 - t345 * t672) * t610 + (t344 * t668 - t345 * t667) * t609 - (t658 * qJD(2) * t342 - t659 * t344 * t518) * t520 / 0.2e1 + (t659 * qJD(2) * t343 - t345 * t658 * t520) * t518 / 0.2e1 + (((-t346 * t663 + t348 * t665 + t708) * t289 + (-t346 * t662 + t348 * t664 + t709) * t288 + t640 * qJD(4)) * t349 + (((-t346 * t660 + t348 * t661 - t735) * t349 - t625 + t702) * qJD(4) + t705) * t347) * t481 - t671 * t512 / 0.2e1 + (-t40 * (-t329 * t345 + t320) - t39 * (-t329 * t344 + t318) - t27 * (-t347 * t509 + t494) - (-t27 * t544 + t40 * t535) * t289 - (t27 * t545 - t39 * t535) * t288 - (t460 * t630 + (-t27 * t477 + t349 * t460) * pkin(6)) * qJD(2) - (t382 * t349 + (t39 * t544 - t40 * t545 + t353) * t347) * qJD(4) + t5 * t475 + t27 * t367 + (t40 * t407 + t419 * t7 + t620) * t345 + (t39 * t407 + t419 * t6 + t384) * t344 - g(1) * (t322 + t632) - g(2) * (t321 + t633) - g(3) * (t493 + t631) - (-rSges(6,2) + t617) * t621 + (g(3) * t563 + t562 * t622) * t666) * m(6) + (-t65 * (t218 * t289 + t320) - t64 * (-t218 * t288 + t318) - t41 * (t173 * t288 - t175 * t289 + t494) - (t453 * t630 + (t349 * t453 - t41 * t477) * pkin(6)) * qJD(2) - ((t118 * t64 - t120 * t65) * t349 + (t65 * (-t220 * t344 - t173) + t64 * (t220 * t345 + t175) + t403) * t347) * qJD(4) - g(1) * (t322 + t529) - g(2) * (t321 + t530) - g(3) * (t493 + t218) - (-rSges(5,3) + t617) * t621 + t16 * t475 + t41 * t367 + (t118 * t16 + t29 * t467 + t41 * t99 + t418 * t65) * t345 + (t101 * t41 + t120 * t16 + t28 * t467 + t418 * t64) * t344) * m(5) + (t36 * t533 + t71 * t546 + (t131 * t532 + t36 * t202 + t71 * t237 + t528 * t73) * t345 + (t130 * t532 + t36 * t203 + t71 * t235 + t528 * t72) * t344 - g(1) * (t468 + t523) - g(2) * (t469 + t524) + g(3) * t527 - t131 * t320 - t130 * t318 - t71 * t494 - ((t131 * t527 + t523 * t71) * t345 + (t130 * t527 + t524 * t71) * t344) * qJD(2)) * m(4) + (g(1) * t275 + g(2) * t273 - g(3) * t310 + t66 * t424 + (qJDD(2) * t697 + t310 * t635) * t307 + (-t420 - (-t273 * t344 - t275 * t345) * qJD(2)) * (qJD(2) * t424 + qJD(1))) * m(3) + (t678 + (t691 * t343 + (t693 * t344 - t345 * t694) * t344) * t680 + (t687 * t343 + (t695 * t344 - t345 * t696) * t344) * t679) * t344 / 0.2e1 - (t677 + (t694 * t343 + (t691 - t693) * t712) * t680 + (t696 * t343 + (t687 - t695) * t712) * t679) * t345 / 0.2e1 + (t344 * t643 - t345 * t642 + t700) * t480; -t508 * t602 + 0.2e1 * (t27 * t614 + t5 * t605) * m(6) + 0.2e1 * (t16 * t605 + t41 * t614) * m(5) + 0.2e1 * (t36 * t605 + t614 * t71) * m(4) + (t508 * t622 + m(4) * (qJD(2) * t71 + t344 * t72 + t345 * t73) + m(5) * (qJD(2) * t41 + t28 * t344 + t29 * t345) + m(6) * (qJD(2) * t27 + t344 * t6 + t345 * t7)) * t347; (t347 * t714 + t349 * t623) * t616 + (t347 * t713 + t349 * t624) * t615 + (-t258 * t627 + t345 * t626 + t411 * t699) * t613 + ((t344 * t674 + t345 * t675) * t349 + t718 * t347 + (-t347 * t623 + t651) * qJD(2)) * t612 + (t259 * t699 + t260 * t627 + t344 * t626) * t611 + ((t344 * t672 + t345 * t673) * t349 + t717 * t347 + (-t347 * t624 + t652) * qJD(2)) * t610 + (t347 * t640 + t349 * t625) * t609 + (t178 * t668 + t179 * t667 + t288 * t643 + t289 * t642 + t294 * t640 + t513 * t641) * t347 / 0.2e1 + t677 * t568 / 0.2e1 + t678 * t566 / 0.2e1 + t671 * t515 / 0.2e1 + ((t346 * t627 - t348 * t699) * t349 + t689 * t347) * t481 + ((t344 * t642 + t345 * t643) * t349 + t641 * t347 + (-t347 * t625 + t349 * t640) * qJD(2)) * t480 + (-(t258 * t40 - t260 * t39 - t27 * t564) * qJD(5) - (-t27 * t548 + t40 * t531) * t289 - (t27 * t547 - t39 * t531) * t288 - (t39 * t548 - t40 * t547) * t513 - g(1) * t548 - g(2) * t547 - (-t346 * t666 - t348 * t676) * t602 + (qJD(2) * t353 + t39 * t601 - t40 * t589 - t550 * t7 + t551 * t6) * t347 + (t382 * qJD(2) + (-t39 * t549 - t534 * t6 + t384) * t345 + (t40 * t549 + t534 * t7 - t620) * t344) * t349) * m(6) + (-g(1) * t146 - g(2) * t150 - g(3) * t287 - t65 * (-t150 * t513 + t287 * t289) - t64 * (t146 * t513 - t287 * t288) - t41 * (-t146 * t289 + t150 * t288) + (-t65 * t101 + t28 * t118 - t29 * t120 + t64 * t99 + (t403 + (-t344 * t65 + t345 * t64) * t220) * qJD(2)) * t347 + (t65 * (-qJD(2) * t120 + t129 * t344) + t64 * (qJD(2) * t118 - t129 * t345) - t16 * t429 + t41 * (t101 * t345 - t344 * t99) + (-t28 * t345 + t29 * t344) * t220) * t349) * m(5) + t700 * t483; (t157 * t40 - t159 * t39 + t387 * t27 + (t288 * t39 - t289 * t40 - g(3) + t5) * t562 + (t27 * t288 - t40 * t513 + g(2) - t6) * t259 - (t27 * t289 - t39 * t513 - g(1) + t7) * t411) * m(6);];
tau = t1;
