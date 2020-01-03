% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:24
% EndTime: 2019-12-31 20:49:51
% DurationCPUTime: 21.20s
% Computational Cost: add. (16919->594), mult. (14620->714), div. (0->0), fcn. (11299->8), ass. (0->350)
t375 = qJ(3) + pkin(8);
t368 = sin(t375);
t369 = cos(t375);
t378 = sin(qJ(3));
t380 = cos(qJ(3));
t755 = Icges(4,5) * t380 + Icges(5,5) * t369 - Icges(4,6) * t378 - Icges(5,6) * t368;
t754 = Icges(6,4) + Icges(5,5);
t753 = Icges(5,6) - Icges(6,6);
t752 = Icges(4,3) + Icges(5,3);
t376 = qJ(1) + qJ(2);
t370 = sin(t376);
t270 = Icges(6,4) * t369 + Icges(6,6) * t368;
t371 = cos(t376);
t439 = t270 * t371;
t194 = Icges(6,2) * t370 + t439;
t750 = t755 * t371;
t714 = t752 * t370 + t750;
t751 = t194 + t714;
t342 = Icges(6,5) * t368;
t475 = Icges(6,3) * t369 - t342;
t621 = Icges(5,4) * t368;
t748 = Icges(5,2) * t369 + t475 + t621;
t620 = Icges(6,5) * t369;
t273 = Icges(6,1) * t368 - t620;
t343 = Icges(5,4) * t369;
t736 = -Icges(5,1) * t368 - t273 - t343;
t587 = t370 * t380;
t588 = t370 * t378;
t592 = t369 * t370;
t594 = t368 * t370;
t715 = -Icges(4,5) * t587 - Icges(5,5) * t592 + Icges(4,6) * t588 + Icges(5,6) * t594 + t752 * t371;
t478 = Icges(6,1) * t369 + t342;
t197 = -Icges(6,4) * t371 + t370 * t478;
t308 = Icges(5,4) * t594;
t199 = Icges(5,1) * t592 - Icges(5,5) * t371 - t308;
t741 = t197 + t199;
t442 = t478 * t371;
t198 = Icges(6,4) * t370 + t442;
t276 = Icges(5,1) * t369 - t621;
t443 = t276 * t371;
t200 = Icges(5,5) * t370 + t443;
t740 = t198 + t200;
t266 = Icges(6,3) * t368 + t620;
t476 = -Icges(5,2) * t368 + t343;
t749 = t266 - t476;
t747 = t276 + t478;
t195 = Icges(5,4) * t592 - Icges(5,2) * t594 - Icges(5,6) * t371;
t209 = Icges(4,4) * t587 - Icges(4,2) * t588 - Icges(4,6) * t371;
t746 = t195 * t368 + t209 * t378;
t727 = Icges(4,5) * t378 + Icges(4,6) * t380 + t754 * t368 + t753 * t369;
t622 = Icges(4,4) * t378;
t326 = Icges(4,1) * t380 - t622;
t444 = t326 * t371;
t212 = Icges(4,5) * t370 + t444;
t745 = -t200 * t592 - t212 * t587;
t338 = Icges(4,4) * t588;
t211 = Icges(4,1) * t587 - Icges(4,5) * t371 - t338;
t744 = t199 * t369 + t211 * t380 - t746;
t323 = Icges(4,2) * t380 + t622;
t372 = Icges(4,4) * t380;
t325 = Icges(4,1) * t378 + t372;
t726 = t323 * t378 - t325 * t380 + t748 * t368 + t736 * t369;
t189 = -Icges(6,6) * t371 + t266 * t370;
t743 = -t189 + t195;
t591 = t369 * t371;
t307 = Icges(6,5) * t591;
t593 = t368 * t371;
t190 = Icges(6,6) * t370 + Icges(6,3) * t593 + t307;
t440 = t476 * t371;
t196 = Icges(5,6) * t370 + t440;
t742 = -t190 + t196;
t739 = t749 * qJD(3);
t738 = t747 * qJD(3);
t374 = qJD(1) + qJD(2);
t734 = t748 * qJD(3) - t753 * t374;
t733 = t736 * qJD(3) + t754 * t374;
t486 = -t190 * t594 + t194 * t371 - t198 * t592;
t477 = -Icges(4,2) * t378 + t372;
t441 = t477 * t371;
t210 = Icges(4,6) * t370 + t441;
t731 = -t371 * t714 - t745;
t695 = -t196 * t594 - t210 * t588 + t731;
t732 = -t486 + t695;
t584 = t371 * t380;
t730 = t190 * t593 + t212 * t584 + t751 * t370 + t740 * t591;
t193 = -Icges(6,2) * t371 + t270 * t370;
t176 = t370 * t193;
t729 = -t189 * t593 - t211 * t584 + t370 * t715 - t741 * t591 - t176;
t728 = t270 + t755;
t656 = t727 * t371;
t657 = t727 * t370;
t696 = rSges(6,1) + pkin(4);
t725 = t368 * t696;
t723 = t196 * t368 + t210 * t378;
t611 = t193 * t371;
t472 = t189 * t368 + t197 * t369;
t667 = t370 * t472;
t78 = -t611 + t667;
t722 = t744 * t370 + t715 * t371 + t78;
t585 = t371 * t378;
t678 = -t195 * t593 - t209 * t585 - t729;
t677 = -t196 * t593 - t210 * t585 + t730;
t721 = t726 * t370 + t656;
t720 = -t726 * t371 + t657;
t694 = rSges(6,3) + qJ(5);
t590 = t370 * t374;
t719 = t734 * t371 - t590 * t749;
t586 = t371 * t374;
t718 = t266 * t586 + t734 * t370 - t374 * t440;
t717 = t733 * t371 - t590 * t747;
t716 = (-t442 - t443) * t374 - t733 * t370;
t281 = rSges(6,1) * t369 + rSges(6,3) * t368;
t713 = pkin(4) * t369 + qJ(5) * t368 + t281;
t290 = t477 * qJD(3);
t291 = t326 * qJD(3);
t711 = -t290 * t378 + t291 * t380 + t738 * t369 + t739 * t368 + t727 * t374 + (-t323 * t380 - t325 * t378 + t736 * t368 - t369 * t748) * qJD(3);
t710 = (Icges(6,2) + t752) * t374 - t727 * qJD(3);
t709 = -t190 * t368 - t212 * t380 - t740 * t369 + t723;
t708 = -t472 - t744;
t673 = t210 * t380 + t212 * t378 + t740 * t368 + t742 * t369;
t672 = t209 * t380 + t211 * t378 + t741 * t368 + t743 * t369;
t707 = t728 * qJD(3) + t726 * t374;
t706 = t370 * (-t371 * t748 + t740) - t371 * (-Icges(5,2) * t592 - t475 * t370 - t308 + t741);
t705 = t743 * t371 + (-Icges(6,1) * t593 + t273 * t371 + t307 - t742) * t370;
t704 = -t736 - t749;
t703 = -t748 + t747;
t700 = t720 * t374;
t699 = (t370 * t677 - t371 * t678) * qJD(3);
t698 = (t732 * t370 - t722 * t371) * qJD(3);
t697 = t721 * t374;
t360 = t371 * rSges(6,2);
t567 = t713 * t370 - t360;
t693 = t374 * t567;
t331 = pkin(7) * t586;
t377 = -qJ(4) - pkin(7);
t341 = t371 * t377;
t344 = qJD(4) * t370;
t542 = qJD(3) * t378;
t511 = t371 * t542;
t633 = pkin(3) * t380;
t365 = pkin(2) + t633;
t632 = pkin(2) - t365;
t139 = -pkin(3) * t511 - t331 + t344 + (t370 * t632 - t341) * t374;
t363 = t371 * pkin(7);
t286 = pkin(2) * t370 - t363;
t549 = -t370 * t365 - t341;
t183 = t286 + t549;
t173 = t374 * t183;
t692 = t139 - t173;
t553 = -t694 * t369 + t725;
t691 = t370 * rSges(6,2) + pkin(4) * t591;
t539 = qJD(5) * t368;
t303 = t371 * t539;
t383 = qJD(1) ^ 2;
t379 = sin(qJ(1));
t635 = pkin(1) * t379;
t537 = t383 * t635;
t489 = t374 * (-pkin(2) * t590 + t331) - t537;
t534 = qJD(3) ^ 2 * t633;
t540 = qJD(4) * t374;
t410 = t374 * t139 + t489 + (-t534 + t540) * t370;
t634 = pkin(3) * t378;
t488 = -t553 - t634;
t455 = t488 * t371;
t538 = qJD(5) * t369;
t561 = -qJD(3) * t713 + t538;
t484 = t538 + t561;
t528 = t369 * t590;
t543 = qJD(3) * t371;
t433 = -t368 * t543 - t528;
t530 = t368 * t590;
t515 = t369 * t543;
t663 = rSges(6,2) * t586 + t515 * t694;
t583 = t433 * t696 - t694 * t530 + t303 + t663;
t17 = (t303 + t583) * t374 + (t370 * t484 + t374 * t455) * qJD(3) + t410;
t514 = t370 * t542;
t318 = pkin(3) * t514;
t381 = cos(qJ(1));
t373 = t381 * pkin(1);
t536 = t383 * t373;
t461 = t374 * t318 + t371 * t540 - t536;
t362 = t370 * pkin(7);
t589 = t370 * t377;
t302 = t374 * t589;
t548 = qJD(4) * t371 + t318;
t520 = t302 + t548;
t140 = (-t371 * t632 - t362) * t374 - t520;
t287 = t371 * pkin(2) + t362;
t244 = t287 * t374;
t581 = -t140 - t244;
t509 = t370 * t539;
t544 = qJD(3) * t370;
t516 = t369 * t544;
t517 = t368 * t544;
t662 = t696 * t517;
t686 = t368 * t586 + t516;
t582 = t509 + t686 * qJ(5) + rSges(6,3) * t516 - t662 + (t281 * t371 + t691) * t374;
t18 = (qJD(3) * t484 - t534) * t371 + ((qJD(3) * t553 - t539) * t370 + t581 - t582) * t374 + t461;
t409 = -t365 - t713;
t436 = qJD(3) * t455;
t628 = pkin(1) * qJD(1);
t531 = t379 * t628;
t551 = t303 + t344;
t460 = -t531 + t551;
t573 = t183 - t286;
t57 = t436 + (-t567 + t573) * t374 + t460;
t532 = t381 * t628;
t315 = t371 * t365;
t494 = t315 - t589;
t184 = t494 - t287;
t566 = rSges(6,1) * t591 + t593 * t694 + t691;
t524 = -t184 - t566;
t642 = -t374 * (t287 - t524) - t509 + t548 + t553 * t544;
t58 = t532 - t642;
t690 = (t58 * (-t634 - t725) * qJD(3) + (-t58 * t377 + t409 * t57) * t374) * t371 + (-t17 * t377 + (-t57 * qJD(5) - t18 * t694) * t368 + (-qJD(3) * t57 * t694 - t18 * t696) * t369 + (-t57 * rSges(6,2) + t409 * t58) * t374) * t370 - t436 * t58;
t545 = rSges(4,2) * t588 + t371 * rSges(4,3);
t213 = rSges(4,1) * t587 - t545;
t205 = t374 * t213;
t259 = t374 * t286;
t541 = qJD(3) * t380;
t510 = t371 * t541;
t445 = rSges(4,3) * t586 + (t374 * t588 - t510) * rSges(4,2);
t689 = -rSges(4,1) * t511 + t205 + t259 + t331 + t445;
t202 = rSges(5,1) * t592 - rSges(5,2) * t594 - t371 * rSges(5,3);
t558 = rSges(5,2) * t530 + rSges(5,3) * t586;
t665 = t259 - t173;
t688 = -rSges(5,1) * t528 + t374 * t202 - t365 * t590 + t344 + t558 + t665;
t687 = t715 + t723;
t284 = rSges(3,1) * t370 + rSges(3,2) * t371;
t598 = t284 * t374;
t217 = -t531 - t598;
t685 = t663 - t551 + t665 + t693;
t684 = -t697 + t698;
t683 = t699 + t700;
t415 = Icges(4,6) * t374 - qJD(3) * t323;
t148 = t370 * t415 + t374 * t441;
t418 = Icges(4,5) * t374 - qJD(3) * t325;
t150 = t370 * t418 + t374 * t444;
t682 = t708 * qJD(3) - t148 * t380 - t150 * t378 + t716 * t368 + t718 * t369;
t147 = t371 * t415 - t477 * t590;
t149 = -t326 * t590 + t371 * t418;
t681 = -t709 * qJD(3) + t147 * t380 + t149 * t378 + t717 * t368 - t719 * t369;
t680 = t707 * t370 + t711 * t371;
t679 = t711 * t370 - t707 * t371;
t676 = t148 * t378 - t150 * t380 + t716 * t369 - t718 * t368 + (-t193 + t715) * t374 + t672 * qJD(3);
t675 = -t673 * qJD(3) - t147 * t378 + t149 * t380 + t719 * t368 + t717 * t369 + t751 * t374;
t674 = t611 + t730;
t671 = (t439 + t708 + t750) * t374 + t710 * t370;
t670 = t710 * t371 + t709 * t374 - t590 * t728;
t669 = 0.2e1 * qJD(3);
t640 = t370 / 0.2e1;
t550 = rSges(5,1) * t591 + t370 * rSges(5,3);
t661 = rSges(4,1) * t584 + t370 * rSges(4,3);
t660 = t370 ^ 2 + t371 ^ 2;
t659 = -t706 * t368 + t705 * t369;
t546 = t325 + t477;
t547 = -t323 + t326;
t655 = (-t368 * t704 + t369 * t703 - t378 * t546 + t380 * t547) * t374;
t654 = t728 * t374;
t279 = rSges(5,1) * t368 + rSges(5,2) * t369;
t241 = t279 * t544;
t204 = -rSges(5,2) * t593 + t550;
t572 = -t184 - t204;
t647 = t374 * (t287 - t572) - t241;
t332 = rSges(4,1) * t378 + rSges(4,2) * t380;
t262 = t332 * t544;
t533 = rSges(4,2) * t585;
t214 = -t533 + t661;
t431 = t214 + t287;
t644 = -t374 * t431 + t262;
t563 = -Icges(4,2) * t587 + t211 - t338;
t565 = t325 * t370 + t209;
t641 = t378 * t563 + t380 * t565;
t639 = -t371 / 0.2e1;
t637 = t374 / 0.2e1;
t631 = rSges(4,1) * t380;
t630 = rSges(5,1) * t369;
t629 = rSges(5,2) * t368;
t518 = t332 * t543;
t435 = -t518 - t531;
t114 = (-t213 - t286) * t374 + t435;
t614 = t114 * t371;
t282 = -t629 + t630;
t261 = t282 * qJD(3);
t605 = t261 * t370;
t577 = -t183 * t544 + t184 * t543;
t576 = -t370 * t183 + t371 * t184;
t564 = -t325 * t371 - t210;
t562 = -t323 * t371 + t212;
t560 = t553 * t370;
t559 = t553 * t371;
t535 = pkin(3) * t585;
t526 = t140 * t544 + t543 * t692;
t525 = t370 * t140 + t371 * t692;
t522 = -rSges(5,1) * t517 - rSges(5,2) * t686;
t513 = t370 * t541;
t521 = rSges(4,1) * t514 + rSges(4,2) * t513 + t374 * t533;
t519 = t360 + t549;
t506 = -pkin(2) - t631;
t505 = -t544 / 0.2e1;
t502 = t543 / 0.2e1;
t500 = -t279 - t634;
t499 = -t282 - t633;
t285 = t371 * rSges(3,1) - rSges(3,2) * t370;
t492 = t660 * t634;
t487 = -t633 - t713;
t243 = rSges(3,1) * t586 - rSges(3,2) * t590;
t459 = t532 - t548;
t74 = t459 + t647;
t485 = t74 * t500;
t483 = t371 * t500;
t481 = -rSges(4,2) * t378 + t631;
t115 = t532 - t644;
t474 = -t115 * t370 - t614;
t467 = t202 * t370 + t204 * t371;
t458 = -pkin(3) * t541 + t561;
t456 = t494 + t550;
t453 = t315 + t566;
t138 = (t213 * t370 + t214 * t371) * qJD(3);
t434 = t467 + t576;
t432 = -t202 + t549;
t430 = qJD(3) * t483 + t344;
t426 = -t378 * t562 + t380 * t564;
t425 = t370 * t506 + t363 + t545;
t423 = -rSges(5,3) * t590 + t520 - t522;
t408 = t430 - t531;
t389 = (t506 * t614 + (t114 * (-rSges(4,3) - pkin(7)) + t115 * t506) * t370) * t374;
t135 = rSges(5,1) * t433 - rSges(5,2) * t515 + t558;
t36 = t135 * t374 + (t374 * t483 - t605) * qJD(3) + t410;
t73 = (-t202 + t573) * t374 + t408;
t386 = (-t36 * t629 + (t73 * (-t365 - t630) - t74 * t377) * t374 + qJD(3) * t485) * t371;
t385 = (((t78 - t667 + t674) * t370 + ((t714 + t746) * t371 + t695 + t729 + t745) * t371) * qJD(3) + t700) * t502 + (-t726 * qJD(3) + t290 * t380 + t291 * t378 + t738 * t368 - t739 * t369) * t374 + (t680 + t681) * t544 / 0.2e1 + (((t371 * t687 - t674 + t677) * t371 + (t370 * t687 - t176 + t486 + t678 - t731) * t370) * qJD(3) + t684 + t697) * t505 - (t679 - t682 + t683) * t543 / 0.2e1 + ((t672 - t721) * t370 + (t673 + t720) * t371) * qJD(3) * t637;
t299 = t481 * qJD(3);
t252 = t332 * t371;
t251 = t332 * t370;
t237 = t279 * t371;
t233 = t279 * t370;
t218 = t285 * t374 + t532;
t186 = -t243 * t374 - t536;
t185 = -t374 * t598 - t537;
t152 = t374 * t661 - t521;
t151 = (-t374 * t587 - t511) * rSges(4,1) + t445;
t137 = t374 * t550 + t522;
t76 = -t536 - t299 * t543 + (-t152 - t244 + t262) * t374;
t75 = t151 * t374 + (-t299 * t370 - t332 * t586) * qJD(3) + t489;
t69 = qJD(3) * t467 + t577;
t46 = -t538 + (t370 * t567 + t371 * t566) * qJD(3) + t577;
t37 = (-qJD(3) * t261 - t534) * t371 + (-t137 + t241 + t581) * t374 + t461;
t7 = (t539 + (t583 + t693) * t371 + (t374 * t524 + t582) * t370) * qJD(3) + t526;
t1 = [m(3) * (t186 * (-t284 - t635) + t185 * (t285 + t373) + (-t243 - t532 + t218) * t217) + t385 + (t18 * (t519 - t635) + t57 * (t302 - t459 + t662) + t17 * (t373 + t453) + (t57 + t531 + t460 + t685) * t58 + t690) * m(6) + (t37 * (t432 - t635) + t73 * (t423 - t532) + t36 * (t373 + t456) + t386 + (-t531 + t73 - t408 + t688) * t74) * m(5) + (t76 * (t425 - t635) + t114 * (t521 - t532) + t75 * (t373 + t431) + t389 + (-t531 + t114 - t435 + t689) * t115) * m(4); t385 + (t17 * t453 + t18 * t519 + (t551 + t685) * t58 + (-t642 + t520 + t662) * t57 + t690) * m(6) + (t36 * t456 + t37 * t432 + t386 + (-t430 + t688) * t74 + (-t548 + t647 + t423) * t73) * m(5) + (t425 * t76 + t431 * t75 + t389 + (t518 + t689) * t115 + (-t644 + t521) * t114) * m(4) + (-(-t217 * t285 - t218 * t284) * t374 + t185 * t285 - t186 * t284 - t217 * t243 - t218 * t598) * m(3); (-(-t69 * t492 + (-t69 * t237 + t499 * t73) * t371 + (-t69 * t233 + t499 * t74) * t370) * qJD(3) + t37 * t483 + t73 * (-pkin(3) * t510 - t261 * t371) + t36 * t500 * t370 + t74 * (-pkin(3) * t513 - t605) + (t135 * t543 + t137 * t544 + t526) * t434 + t69 * (t135 * t371 + t137 * t370 + t525) + (-t73 * t233 - t74 * (-t237 - t535) + (t69 * t202 + t485) * t371 + (t73 * t279 + t572 * t69) * t370 + (t202 * t371 + t370 * t572) * t434 * qJD(3)) * t374) * m(5) + (-(t114 * t251 - t115 * t252) * t374 - (t138 * (-t251 * t370 - t252 * t371) + t474 * t481) * qJD(3) + 0.2e1 * t138 * ((t151 + t205) * t371 + (-t214 * t374 + t152) * t370) + t474 * t299 + ((-t115 * t374 - t76) * t371 + (t114 * t374 - t75) * t370) * t332) * m(4) - ((t703 * t368 + t704 * t369 + t378 * t547 + t380 * t546) * t374 + ((t370 * t562 - t371 * t563) * t380 + (t370 * t564 + t371 * t565) * t378 + t706 * t369 + t705 * t368) * qJD(3)) * t374 / 0.2e1 + ((t374 * t673 + t682) * t371 + (t374 * t672 + t681) * t370) * t637 + ((-t656 * t544 + t654) * t370 + ((t641 * t371 + (t426 + t657) * t370 + t659) * qJD(3) + t655) * t371) * t505 + ((-t657 * t543 - t654) * t371 + ((t426 * t370 + (t641 + t656) * t371 + t659) * qJD(3) + t655) * t370) * t502 + (-(t368 * t46 + (t370 * t58 + t371 * t57) * t369) * qJD(5) - (t57 * t560 + t58 * (-t535 - t559)) * t374 - (-t46 * t492 + (-t46 * t559 + t487 * t57) * t371 + (-t46 * t560 + t487 * t58) * t370) * qJD(3) + t7 * t576 + t46 * t525 + (t18 * t488 + t57 * t458 + t7 * t566 + t46 * t583 + (t46 * t567 + t488 * t58) * t374) * t371 + (t17 * t488 + t58 * t458 + t7 * t567 + t46 * t582 + (t46 * t524 + t553 * t57) * t374) * t370) * m(6) + (t680 * t374 + ((t676 * t371 + t677 * t374) * t371 + (t670 * t370 + t678 * t374 + (-t671 + t675) * t371) * t370) * t669) * t640 + (t679 * t374 + ((t671 * t371 + t732 * t374) * t371 + (t675 * t370 + t722 * t374 + (-t670 + t676) * t371) * t370) * t669) * t639 + (t684 + t698) * t590 / 0.2e1 + (t683 + t699) * t586 / 0.2e1; 0.2e1 * (t17 * t639 + t18 * t640) * m(6) + 0.2e1 * (t36 * t639 + t37 * t640) * m(5); (-t369 * t7 + 0.2e1 * (t17 * t640 + t18 * t371 / 0.2e1 + (0.1e1 / 0.2e1 - t660 / 0.2e1) * qJD(3) * t46) * t368) * m(6);];
tauc = t1(:);
