% Calculate vector of inverse dynamics joint torques for
% S5PRRRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:31
% EndTime: 2019-12-05 16:44:12
% DurationCPUTime: 28.94s
% Computational Cost: add. (19246->684), mult. (16285->864), div. (0->0), fcn. (12840->6), ass. (0->381)
t746 = Icges(5,4) + Icges(6,4);
t747 = Icges(5,1) + Icges(6,1);
t745 = Icges(5,2) + Icges(6,2);
t739 = Icges(5,6) + Icges(6,6);
t740 = Icges(5,5) + Icges(6,5);
t375 = qJ(3) + qJ(4);
t369 = cos(t375);
t744 = t746 * t369;
t368 = sin(t375);
t600 = Icges(6,4) * t368;
t601 = Icges(5,4) * t368;
t743 = t747 * t369 - t600 - t601;
t742 = -t745 * t368 + t744;
t373 = pkin(8) + qJ(2);
t367 = cos(t373);
t741 = t739 * t367;
t738 = Icges(5,3) + Icges(6,3);
t366 = sin(t373);
t566 = t366 * t369;
t567 = t366 * t368;
t718 = t746 * t566 - t745 * t567 - t741;
t734 = t747 * t368 + t744;
t737 = t746 * t567;
t736 = t740 * t367;
t717 = t739 * t366 + t742 * t367;
t716 = t747 * t566 - t736 - t737;
t715 = t740 * t366 + t743 * t367;
t714 = -t739 * t368 + t740 * t369;
t280 = Icges(6,2) * t369 + t600;
t282 = Icges(5,2) * t369 + t601;
t735 = t280 + t282;
t704 = t717 * t368;
t726 = t738 * t367;
t719 = -t740 * t566 + t739 * t567 + t726;
t646 = t704 + t719;
t561 = t367 * t369;
t650 = t738 * t366 + t714 * t367;
t699 = t650 * t366 + t715 * t561;
t727 = t718 * t368;
t647 = -t369 * t716 + t727;
t722 = t366 * t647;
t731 = t367 * t646 - t699 - t722;
t730 = -t280 * t367 + t715;
t729 = t740 * t368 + t739 * t369;
t728 = -t280 + t743;
t725 = t742 + t734;
t724 = -t734 * t367 - t717;
t723 = t734 * t366 + t718;
t721 = t368 * t735 - t369 * t734;
t720 = t715 * t566;
t374 = qJD(3) + qJD(4);
t713 = t729 * t374;
t570 = t282 * t374;
t712 = t728 * t374 - t570;
t711 = t725 * t374;
t710 = -t715 * qJD(2) + t723 * t374;
t709 = t724 * t374 + (-t366 * t743 + t736) * qJD(2);
t270 = t366 * t374;
t708 = t717 * qJD(2) - t270 * t280 - t366 * t570 + t716 * t374;
t707 = -t367 * t570 + t730 * t374 + (-t366 * t742 + t741) * qJD(2);
t706 = t729 * t367;
t705 = t729 * t366;
t378 = -pkin(7) - pkin(6);
t372 = -qJ(5) + t378;
t703 = t372 - rSges(6,3);
t668 = -t721 * t366 - t706;
t667 = -t721 * t367 + t705;
t701 = t717 * t567 - t720;
t700 = t650 * qJD(2);
t698 = t719 * t366 - t716 * t561;
t674 = t719 * t367 - t722;
t673 = -t650 * t367 - t701;
t562 = t367 * t368;
t672 = -t718 * t562 - t698;
t671 = -t717 * t562 + t699;
t697 = t729 * qJD(2) - t711 * t368 + t712 * t369;
t696 = -t707 * t368 + t709 * t369 + t700;
t695 = t719 * qJD(2) + t708 * t368 + t710 * t369;
t271 = t367 * t374;
t694 = t723 * t271 + t724 * t270 + (-t282 + t728) * qJD(2);
t693 = (t745 * t566 - t716 + t737) * t271 + (-t282 * t367 + t730) * t270 + t725 * qJD(2);
t692 = t721 * qJD(2) + t714 * t374;
t691 = t647 * qJD(2) - t713 * t366 + t700;
t690 = -t713 * t367 + (-t714 * t366 - t715 * t369 + t704 + t726) * qJD(2);
t689 = t667 * qJD(2);
t362 = pkin(4) * t369;
t377 = cos(qJ(3));
t371 = t377 * pkin(3);
t363 = t371 + pkin(2);
t501 = t362 + t363;
t688 = -rSges(6,1) * t566 + rSges(6,2) * t567 - t366 * t501 - t703 * t367;
t349 = t366 * rSges(6,3);
t687 = rSges(6,1) * t561 - rSges(6,2) * t562 + t367 * t501 + t349;
t376 = sin(qJ(3));
t509 = qJD(3) * t376;
t500 = pkin(3) * t509;
t558 = t368 * t374;
t267 = -pkin(4) * t558 - t500;
t338 = qJD(5) * t366;
t514 = qJD(2) * t366;
t493 = t368 * t514;
t557 = t369 * t374;
t497 = t367 * t557;
t513 = qJD(2) * t367;
t686 = rSges(6,3) * t513 + (t493 - t497) * rSges(6,2) + t367 * t267 + t338;
t685 = t668 * qJD(2);
t684 = -t691 * t366 + t695 * t367;
t683 = t690 * t366 + t696 * t367;
t682 = t695 * t366 + t691 * t367;
t681 = t696 * t366 - t690 * t367;
t680 = t673 * t270 - t674 * t271 + t685;
t679 = t671 * t270 - t672 * t271 + t689;
t678 = t710 * t368 - t708 * t369;
t677 = t709 * t368 + t707 * t369;
t676 = t692 * t366 + t697 * t367;
t675 = t697 * t366 - t692 * t367;
t670 = t716 * t368 + t718 * t369;
t669 = t715 * t368 + t717 * t369;
t498 = t367 * t558;
t403 = -t369 * t514 - t498;
t488 = t367 * t509;
t456 = pkin(3) * t488;
t458 = t363 - t501;
t519 = t378 - t372;
t607 = rSges(6,1) * t403 + t456 + (t366 * t458 + t367 * t519) * qJD(2) + t686;
t609 = t369 * rSges(6,2);
t288 = rSges(6,1) * t368 + t609;
t233 = t288 * t366;
t339 = qJD(5) * t367;
t313 = t366 * t500;
t523 = t378 * t514 + t313;
t565 = t366 * t372;
t360 = t369 * rSges(6,1);
t611 = rSges(6,2) * t368;
t649 = t360 - t611;
t606 = -t374 * t233 + t366 * t267 - t339 + t523 + (t349 - t565 + (-t458 + t649) * t367) * qJD(2);
t335 = t367 * t378;
t525 = t366 * t363 + t335;
t554 = t525 + t688;
t306 = t367 * t363;
t553 = t366 * t519 - t306 + t687;
t658 = -t693 * t368 + t694 * t369;
t657 = t714 * qJD(2) - t706 * t270 + t705 * t271;
t451 = t649 + t362;
t289 = rSges(5,1) * t368 + rSges(5,2) * t369;
t234 = t289 * t366;
t236 = t289 * t367;
t361 = t369 * rSges(5,1);
t648 = -rSges(5,2) * t368 + t361;
t204 = rSges(5,1) * t566 - rSges(5,2) * t567 - t367 * rSges(5,3);
t350 = t366 * rSges(5,3);
t206 = rSges(5,1) * t561 - rSges(5,2) * t562 + t350;
t356 = t367 * pkin(6);
t274 = pkin(2) * t366 - t356;
t185 = t274 - t525;
t355 = t366 * pkin(6);
t275 = t367 * pkin(2) + t355;
t459 = -t366 * t378 + t306;
t186 = t459 - t275;
t510 = qJD(3) * t367;
t511 = qJD(3) * t366;
t487 = -t185 * t511 + t186 * t510 + qJD(1);
t65 = t204 * t270 + t206 * t271 + t487;
t405 = -t271 * t289 - t456;
t543 = t185 - t274;
t494 = -t204 + t543;
t72 = qJD(2) * t494 + t405;
t542 = -t186 - t206;
t73 = -t270 * t289 - t313 + (t275 - t542) * qJD(2);
t656 = -(qJD(2) * t234 - t271 * t648) * t72 - t65 * (-t270 * t234 - t236 * t271) - t73 * (-qJD(2) * t236 - t270 * t648);
t504 = qJD(2) * qJD(3);
t330 = t366 * t504;
t503 = qJD(2) * qJD(4);
t190 = t366 * t503 + t330 + (-qJDD(3) - qJDD(4)) * t367;
t246 = t649 * t374;
t266 = -qJDD(3) * t367 + t330;
t556 = t377 * qJD(3) ^ 2;
t620 = pkin(3) * t376;
t420 = -pkin(3) * t367 * t556 + t266 * t620;
t454 = t543 + t554;
t614 = pkin(2) - t363;
t148 = (-t367 * t614 - t355) * qJD(2) - t523;
t264 = t275 * qJD(2);
t555 = -t148 - t264;
t12 = qJDD(5) * t366 + t190 * t288 - t246 * t271 + (t190 * t368 - t271 * t557) * pkin(4) + t454 * qJDD(2) + (t339 + t555 - t606) * qJD(2) + t420;
t654 = t12 - g(1);
t265 = qJDD(3) * t366 + t367 * t504;
t189 = qJDD(4) * t366 + t367 * t503 + t265;
t336 = pkin(6) * t513;
t147 = -t456 - t336 + (t366 * t614 - t335) * qJD(2);
t530 = qJD(2) * (-pkin(2) * t514 + t336) + qJDD(2) * t275;
t395 = qJD(2) * t147 + qJDD(2) * t186 + (-t265 * t376 - t366 * t556) * pkin(3) + t530;
t13 = -qJDD(5) * t367 - t189 * t288 - t246 * t270 + t553 * qJDD(2) + (-t189 * t368 - t270 * t557) * pkin(4) + (t338 + t607) * qJD(2) + t395;
t653 = t13 - g(2);
t619 = pkin(4) * t368;
t479 = -t288 - t619;
t495 = -t186 - t553;
t62 = -t313 - t339 + t479 * t270 + (t275 - t495) * qJD(2);
t652 = qJD(2) * t62 + t12;
t269 = qJD(2) * t274;
t651 = qJD(2) * t185 - t269;
t370 = Icges(4,4) * t377;
t437 = -Icges(4,2) * t376 + t370;
t319 = Icges(4,1) * t376 + t370;
t644 = g(1) * t367 + g(2) * t366;
t563 = t366 * t377;
t564 = t366 * t376;
t590 = Icges(4,3) * t367;
t209 = Icges(4,5) * t563 - Icges(4,6) * t564 - t590;
t326 = Icges(4,4) * t564;
t599 = Icges(4,5) * t367;
t213 = Icges(4,1) * t563 - t326 - t599;
t593 = Icges(4,6) * t367;
t211 = Icges(4,4) * t563 - Icges(4,2) * t564 - t593;
t577 = t211 * t376;
t428 = -t213 * t377 + t577;
t82 = -t209 * t367 - t366 * t428;
t316 = Icges(4,5) * t377 - Icges(4,6) * t376;
t315 = Icges(4,5) * t376 + Icges(4,6) * t377;
t407 = qJD(3) * t315;
t602 = Icges(4,4) * t376;
t320 = Icges(4,1) * t377 - t602;
t214 = Icges(4,5) * t366 + t320 * t367;
t212 = Icges(4,6) * t366 + t367 * t437;
t576 = t212 * t376;
t427 = -t214 * t377 + t576;
t639 = -t367 * t407 + (-t316 * t366 + t427 + t590) * qJD(2);
t210 = Icges(4,3) * t366 + t316 * t367;
t516 = qJD(2) * t210;
t638 = qJD(2) * t428 - t366 * t407 + t516;
t317 = Icges(4,2) * t377 + t602;
t422 = t317 * t376 - t319 * t377;
t635 = t422 * qJD(2) + t316 * qJD(3);
t535 = -Icges(4,2) * t563 + t213 - t326;
t537 = t319 * t366 + t211;
t634 = -t376 * t535 - t377 * t537;
t631 = m(2) + m(3);
t630 = t189 / 0.2e1;
t629 = t190 / 0.2e1;
t628 = t265 / 0.2e1;
t627 = t266 / 0.2e1;
t626 = -t270 / 0.2e1;
t625 = t270 / 0.2e1;
t624 = -t271 / 0.2e1;
t623 = t271 / 0.2e1;
t622 = t366 / 0.2e1;
t621 = -t367 / 0.2e1;
t616 = -qJD(2) / 0.2e1;
t615 = qJD(2) / 0.2e1;
t613 = rSges(4,1) * t377;
t351 = t366 * rSges(4,3);
t610 = t367 * t73;
t608 = qJDD(2) / 0.2e1;
t43 = -t270 * t554 + t271 * t553 + t487;
t587 = qJD(2) * t43;
t321 = rSges(4,1) * t376 + rSges(4,2) * t377;
t490 = t321 * t510;
t520 = rSges(4,2) * t564 + t367 * rSges(4,3);
t217 = rSges(4,1) * t563 - t520;
t533 = -t217 - t274;
t121 = qJD(2) * t533 - t490;
t585 = t121 * t366;
t584 = t121 * t367;
t559 = t367 * t377;
t560 = t367 * t376;
t218 = rSges(4,1) * t559 - rSges(4,2) * t560 + t351;
t165 = t218 + t275;
t122 = qJD(2) * t165 - t321 * t511;
t256 = t321 * t367;
t583 = t122 * t256;
t575 = t270 * t369;
t569 = t315 * t366;
t568 = t315 * t367;
t547 = -t366 * t185 + t367 * t186;
t546 = t366 * t204 + t367 * t206;
t545 = -t366 * t209 - t213 * t559;
t544 = t366 * t210 + t214 * t559;
t536 = -t319 * t367 - t212;
t534 = -t317 * t367 + t214;
t297 = pkin(4) * t493;
t531 = t288 * t514 + t297;
t526 = rSges(5,2) * t493 + rSges(5,3) * t513;
t512 = qJD(2) * t376;
t524 = rSges(4,2) * t366 * t512 + rSges(4,3) * t513;
t522 = -t317 + t320;
t521 = t319 + t437;
t515 = qJD(2) * t316;
t508 = qJD(3) * t377;
t141 = -t366 * t422 - t568;
t505 = t141 * qJD(2);
t118 = rSges(5,1) * t403 - rSges(5,2) * t497 + t526;
t416 = t289 * t374;
t120 = -t366 * t416 + (t367 * t648 + t350) * qJD(2);
t502 = t367 * t118 + t366 * t120 + t204 * t513;
t499 = pkin(3) * t508;
t496 = t367 * t147 + t366 * t148 - t185 * t513;
t491 = t367 * t512;
t486 = -pkin(2) - t613;
t485 = t514 / 0.2e1;
t484 = t513 / 0.2e1;
t483 = -t511 / 0.2e1;
t482 = t511 / 0.2e1;
t481 = -t510 / 0.2e1;
t480 = t510 / 0.2e1;
t404 = -t289 - t620;
t478 = t376 * (-t366 ^ 2 - t367 ^ 2);
t235 = t288 * t367;
t469 = -t270 * t233 - t235 * t271;
t182 = t214 * t563;
t468 = t210 * t367 - t182;
t465 = -t209 + t576;
t464 = -qJD(2) * t235 - t270 * t649;
t455 = -t366 * t554 + t367 * t553;
t453 = -pkin(4) * t557 - t246;
t247 = t648 * t374;
t452 = -t247 - t499;
t311 = -t619 - t620;
t448 = qJD(2) * t233 - t271 * t451;
t273 = rSges(3,1) * t367 - rSges(3,2) * t366;
t272 = rSges(3,1) * t366 + rSges(3,2) * t367;
t322 = -rSges(4,2) * t376 + t613;
t124 = t212 * t377 + t214 * t376;
t408 = qJD(3) * t317;
t137 = -t367 * t408 + (-t366 * t437 + t593) * qJD(2);
t409 = qJD(3) * t319;
t139 = -t367 * t409 + (-t320 * t366 + t599) * qJD(2);
t385 = -qJD(3) * t124 - t137 * t376 + t139 * t377 + t516;
t123 = t211 * t377 + t213 * t376;
t138 = qJD(2) * t212 - t366 * t408;
t140 = qJD(2) * t214 - t366 * t409;
t386 = qJD(2) * t209 - qJD(3) * t123 - t138 * t376 + t140 * t377;
t445 = -(t366 * t638 + t386 * t367) * t367 + (t366 * t639 + t385 * t367) * t366;
t444 = -(t386 * t366 - t367 * t638) * t367 + (t385 * t366 - t367 * t639) * t366;
t443 = -t366 * t73 - t367 * t72;
t83 = -t212 * t564 - t468;
t442 = t366 * t83 - t367 * t82;
t84 = -t211 * t560 - t545;
t85 = -t212 * t560 + t544;
t441 = t366 * t85 - t367 * t84;
t434 = -t122 * t366 - t584;
t143 = -rSges(4,2) * t367 * t508 + (-t377 * t514 - t488) * rSges(4,1) + t524;
t255 = t321 * t366;
t144 = -qJD(3) * t255 + (t322 * t367 + t351) * qJD(2);
t433 = t143 * t367 + t144 * t366;
t426 = t217 * t366 + t218 * t367;
t423 = t317 * t377 + t319 * t376;
t421 = -t501 - t360;
t419 = t606 * t366 + t607 * t367 - t554 * t513;
t418 = -t288 + t311;
t406 = t147 * t510 + t148 * t511 - t265 * t185 - t186 * t266 + qJDD(1);
t402 = t453 - t499;
t399 = -t376 * t534 + t377 * t536;
t398 = (-t376 * t521 + t377 * t522) * qJD(2);
t397 = t271 * t479 + t338 - t456;
t293 = t437 * qJD(3);
t294 = t320 * qJD(3);
t384 = qJD(2) * t315 - qJD(3) * t423 - t293 * t376 + t294 * t377;
t383 = (t671 * t366 - t672 * t367) * t630 + (t673 * t366 - t674 * t367) * t629 + (t366 * t657 + t367 * t658) * t626 + (t684 * t367 + t683 * t366 + (t672 * t366 + t671 * t367) * qJD(2)) * t625 + (t682 * t367 + t681 * t366 + (t674 * t366 + t673 * t367) * qJD(2)) * t624 + (t366 * t658 - t367 * t657) * t623 + (t676 * qJD(2) + t667 * qJDD(2) + t671 * t189 + t672 * t190 + t683 * t270 + t684 * t271) * t622 + (t675 * qJD(2) + t668 * qJDD(2) + t673 * t189 + t674 * t190 + t681 * t270 + t682 * t271) * t621 + (t694 * t368 + t693 * t369) * t616 + (t678 * t367 + t677 * t366 + (t670 * t366 + t669 * t367) * qJD(2)) * t615 + (t669 * t366 - t670 * t367) * t608 + t680 * t485 + t679 * t484;
t298 = t322 * qJD(3);
t263 = t367 * t311;
t262 = t366 * t311;
t216 = pkin(3) * t560 + t263;
t215 = pkin(3) * t564 + t262;
t142 = -t367 * t422 + t569;
t131 = t142 * qJD(2);
t100 = qJD(3) * t426 + qJD(1);
t71 = t384 * t366 - t367 * t635;
t70 = t366 * t635 + t384 * t367;
t67 = qJD(2) * t143 + qJDD(2) * t218 - t265 * t321 - t298 * t511 + t530;
t66 = -t298 * t510 + t266 * t321 + t533 * qJDD(2) + (-t144 - t264) * qJD(2);
t64 = -qJD(3) * t427 + t137 * t377 + t139 * t376;
t63 = -qJD(3) * t428 + t138 * t377 + t140 * t376;
t61 = qJD(2) * t454 + t397;
t60 = qJD(3) * t433 + t217 * t265 - t218 * t266 + qJDD(1);
t49 = qJD(3) * t441 + t131;
t48 = qJD(3) * t442 + t505;
t42 = qJD(2) * t118 + qJDD(2) * t206 - t189 * t289 - t247 * t270 + t395;
t41 = t190 * t289 - t247 * t271 + t494 * qJDD(2) + (-t120 + t555) * qJD(2) + t420;
t18 = t118 * t271 + t120 * t270 + t189 * t204 - t190 * t206 + t406;
t9 = -t189 * t554 - t190 * t553 + t270 * t606 + t271 * t607 + t406;
t1 = [t631 * qJDD(1) + m(4) * t60 + m(5) * t18 + m(6) * t9 + (-m(4) - m(5) - m(6) - t631) * g(3); (t131 + ((t83 - t182 + (t210 + t577) * t367 + t545) * t367 + t544 * t366) * qJD(3)) * t480 - m(3) * (-g(1) * t272 + g(2) * t273) + (t124 + t142) * t628 + (t141 + t123) * t627 + (((t650 + t727) * t367 + t673 + t698 + t701) * t271 + (t674 - t731) * t270 + t689) * t623 + (-t505 + ((t367 * t465 - t544 + t85) * t367 + (t366 * t465 + t468 + t84) * t366) * qJD(3) + t48) * t483 + (t64 + t70) * t482 + (t63 + t71 + t49) * t481 + (-qJD(3) * t422 + t293 * t377 + t294 * t376 + t712 * t368 + t711 * t369) * qJD(2) + (-(qJD(2) * t554 + t397 - t61 + t651) * t62 + t61 * t339 + t62 * (-rSges(6,1) * t498 + t686) + t61 * (t288 * t374 - t267) * t366 + ((t61 * (t421 + t611) - t62 * t372) * t367 + (t62 * t421 + t61 * t703) * t366) * qJD(2) + t653 * (-t565 + t687) + t654 * t688) * m(6) + (-(-qJD(2) * t204 + t405 + t651 - t72) * t73 + t72 * (t523 + (rSges(5,1) * t558 + rSges(5,2) * t557) * t366) + t73 * t526 + (-t416 - t500) * t610 + ((-t72 * rSges(5,3) + t73 * (-t363 - t361)) * t366 + (t72 * (-t363 - t648) - t73 * t378) * t367) * qJD(2) + (t42 - g(2)) * (t206 + t459) + (t41 - g(1)) * (-t204 - t525)) * m(5) + (-(-qJD(2) * t217 - t121 - t269 - t490) * t122 + t122 * (t336 + t524) + (t321 * t585 - t583) * qJD(3) + ((-pkin(2) - t322) * t584 + (t121 * (-rSges(4,3) - pkin(6)) + t122 * t486) * t366) * qJD(2) + (t67 - g(2)) * t165 + (t66 - g(1)) * (t486 * t366 + t356 + t520)) * m(4) + (t667 + t669) * t630 + (t668 + t670) * t629 + ((t671 + t731) * t271 + ((t647 + t650) * t367 + t646 * t366 + t672 - t720) * t270 + t680 - t685) * t626 + (t676 + t677) * t625 + (t423 + m(3) * (t272 ^ 2 + t273 ^ 2) + Icges(3,3) + t735 * t369 + t734 * t368) * qJDD(2) + (t675 - t678 + t679) * t624; t383 + ((t84 * t366 + t85 * t367) * qJD(2) + t445) * t482 + ((-t511 * t568 + t515) * t366 + (t398 + (-t634 * t367 + (t569 + t399) * t366) * qJD(3)) * t367) * t483 + (-t123 * t367 + t124 * t366) * t608 + t441 * t628 + t442 * t627 + t48 * t485 + ((t376 * t522 + t377 * t521) * qJD(2) + ((t366 * t534 - t367 * t535) * t377 + (t366 * t536 + t367 * t537) * t376) * qJD(3)) * t616 + ((t82 * t366 + t83 * t367) * qJD(2) + t444) * t481 + (t366 * t64 - t367 * t63 + (t123 * t366 + t124 * t367) * qJD(2)) * t615 + (qJD(2) * t70 + qJD(3) * t445 + qJDD(2) * t142 + t265 * t85 + t266 * t84) * t622 + (qJD(2) * t71 + qJD(3) * t444 + qJDD(2) * t141 + t265 * t83 + t266 * t82) * t621 + t49 * t484 + ((-t510 * t569 - t515) * t367 + (t398 + (t399 * t366 + (t568 - t634) * t367) * qJD(3)) * t366) * t480 + (-g(1) * (t263 - t235) - g(2) * (t262 - t233) - g(3) * (t371 + t451) - t61 * (-qJD(2) * t215 + t448) - t62 * (-pkin(4) * t575 + qJD(2) * t216 + t464) - t43 * (t215 * t270 + t216 * t271 + t469) - (-t62 * t491 + ((-t366 * t62 - t367 * t61) * t377 + t43 * t478) * qJD(3)) * pkin(3) + t61 * t531 + t9 * (t455 + t547) + t43 * (t419 + t496) + (t402 * t61 + t418 * t652) * t367 + (t13 * t418 + t402 * t62 + t495 * t587) * t366) * m(6) + (-g(3) * (t648 + t371) - t644 * t404 - (-t73 * t491 + (t377 * t443 + t478 * t65) * qJD(3)) * pkin(3) + t18 * (t546 + t547) + t65 * (t496 + t502) + (t452 * t72 + (qJD(2) * t73 + t41) * t404) * t367 + (t42 * t404 + t73 * t452 + (t72 * t289 + t542 * t65) * qJD(2)) * t366 + t656) * m(5) + (t60 * t426 + t100 * ((t217 * t367 - t218 * t366) * qJD(2) + t433) + t434 * t298 + (-t67 * t366 - t66 * t367 + (-t122 * t367 + t585) * qJD(2)) * t321 - (t121 * t255 - t583) * qJD(2) - (t100 * (-t255 * t366 - t256 * t367) + t434 * t322) * qJD(3) + g(1) * t256 + g(2) * t255 - g(3) * t322) * m(4); t383 + (-g(3) * t451 - t62 * t464 - t43 * t469 - (-t62 * t575 + (-t62 * t513 + t43 * (-t270 * t366 - t271 * t367)) * t368) * pkin(4) + t9 * t455 + t43 * t419 + (t453 * t62 - t553 * t587) * t366 + (t13 * t366 + t367 * t652) * t479 + (t367 * t453 - t297 - t448 + t531) * t61 - t644 * (-t609 + (-rSges(6,1) - pkin(4)) * t368)) * m(6) + (t18 * t546 + t65 * (-t206 * t514 + t502) + t443 * t247 + (-t42 * t366 - t41 * t367 + (t366 * t72 - t610) * qJD(2)) * t289 + g(1) * t236 + g(2) * t234 - g(3) * t648 + t656) * m(5); ((t270 * t43 - t653) * t367 + (-t271 * t43 + t654) * t366) * m(6);];
tau = t1;
