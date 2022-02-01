% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRP3
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
% m [6x1]
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
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:29:35
% EndTime: 2022-01-23 09:30:11
% DurationCPUTime: 30.02s
% Computational Cost: add. (17917->644), mult. (15561->822), div. (0->0), fcn. (12116->8), ass. (0->376)
t726 = Icges(5,6) + Icges(6,6);
t735 = Icges(5,1) + Icges(6,1);
t734 = Icges(5,4) + Icges(6,4);
t731 = Icges(5,5) + Icges(6,5);
t733 = Icges(5,2) + Icges(6,2);
t359 = qJ(1) + pkin(8);
t351 = cos(t359);
t732 = t726 * t351;
t730 = Icges(5,3) + Icges(6,3);
t350 = sin(t359);
t360 = qJ(3) + qJ(4);
t353 = cos(t360);
t562 = t350 * t353;
t352 = sin(t360);
t563 = t350 * t352;
t715 = t734 * t562 - t733 * t563 - t732;
t596 = Icges(6,4) * t352;
t597 = Icges(5,4) * t352;
t729 = t735 * t353 - t596 - t597;
t728 = t734 * t563;
t727 = t731 * t351;
t701 = t735 * t562 - t727 - t728;
t699 = -t726 * t352 + t731 * t353;
t345 = Icges(6,4) * t353;
t268 = Icges(6,1) * t352 + t345;
t346 = Icges(5,4) * t353;
t270 = Icges(5,1) * t352 + t346;
t725 = -t268 - t270;
t723 = -t733 * t352 + t345 + t346;
t700 = t731 * t350 + t729 * t351;
t722 = t731 * t352 + t726 * t353;
t266 = Icges(5,2) * t353 + t597;
t721 = -t266 + t729;
t720 = t715 * t352;
t719 = t730 * t351;
t718 = t723 - t725;
t702 = -t731 * t562 + t726 * t563 + t719;
t674 = t730 * t350 + t699 * t351;
t712 = t726 * t350 + t723 * t351;
t264 = Icges(6,2) * t353 + t596;
t717 = t725 * t353 + (t264 + t266) * t352;
t696 = -t701 * t353 + t720;
t716 = t700 * t562;
t714 = t270 * t350 + t715;
t713 = -t270 * t351 - t712;
t358 = qJD(3) + qJD(4);
t566 = t264 * t358;
t711 = t358 * t721 - t566;
t710 = t718 * t358;
t691 = t722 * t351;
t690 = t722 * t350;
t709 = -t266 * t351 + t700;
t708 = t717 * t350 + t691;
t707 = -t717 * t351 + t690;
t706 = t696 * t350;
t705 = t712 * t563 - t716;
t558 = t351 * t353;
t704 = t674 * t350 + t700 * t558;
t703 = t702 * t350 - t701 * t558;
t658 = t702 * t351 - t706;
t657 = -t674 * t351 - t705;
t559 = t351 * t352;
t656 = -t715 * t559 - t703;
t655 = -t712 * t559 + t704;
t698 = qJD(1) * t722 - t710 * t352 + t711 * t353;
t498 = rSges(5,1) * t562;
t365 = -pkin(7) - pkin(6);
t322 = t351 * t365;
t363 = cos(qJ(3));
t347 = t363 * pkin(3) + pkin(2);
t524 = t350 * t347 + t322;
t362 = sin(qJ(1));
t614 = pkin(1) * t362;
t697 = -t498 - t614 - t524;
t403 = t268 * t358;
t695 = -t700 * qJD(1) + t350 * t403 + t714 * t358;
t694 = -t351 * t403 + t713 * t358 + (-t350 * t729 + t727) * qJD(1);
t255 = t350 * t358;
t693 = t712 * qJD(1) - t255 * t266 - t350 * t566 + t701 * t358;
t692 = -t351 * t566 + t709 * t358 + (-t723 * t350 + t732) * qJD(1);
t689 = t712 * t352;
t688 = t717 * qJD(1) + t699 * t358;
t687 = t707 * qJD(1);
t357 = -qJ(5) + t365;
t686 = rSges(6,3) - t357;
t685 = t708 * qJD(1);
t684 = t674 * qJD(1);
t670 = t689 + t702;
t683 = t670 * t351 - t704 - t706;
t256 = t351 * t358;
t682 = t657 * t255 - t658 * t256 - t685;
t681 = t655 * t255 - t656 * t256 + t687;
t680 = t695 * t352 - t693 * t353;
t679 = t694 * t352 + t692 * t353;
t678 = t688 * t350 + t698 * t351;
t677 = t698 * t350 - t688 * t351;
t676 = t701 * t352 + t715 * t353;
t675 = t700 * t352 + t712 * t353;
t673 = -t692 * t352 + t694 * t353 + t684;
t672 = t702 * qJD(1) + t693 * t352 + t695 * t353;
t671 = (t268 * t350 + t714) * t256 + (-t268 * t351 + t713) * t255 + (-t264 + t721) * qJD(1);
t669 = (t733 * t562 - t701 + t728) * t256 + (-t264 * t351 + t709) * t255 + t718 * qJD(1);
t668 = t696 * qJD(1) - t690 * t358 + t684;
t667 = -t691 * t358 + (-t699 * t350 - t700 * t353 + t689 + t719) * qJD(1);
t611 = pkin(4) * t353;
t284 = t347 + t611;
t335 = t350 * rSges(6,3);
t665 = rSges(6,1) * t558 - rSges(6,2) * t559 + t351 * t284 + t335;
t361 = sin(qJ(3));
t506 = qJD(3) * t361;
t495 = pkin(3) * t506;
t555 = t352 * t358;
t253 = -pkin(4) * t555 - t495;
t324 = qJD(5) * t350;
t511 = qJD(1) * t350;
t489 = t352 * t511;
t554 = t353 * t358;
t491 = t351 * t554;
t510 = qJD(1) * t351;
t664 = rSges(6,3) * t510 + (t489 - t491) * rSges(6,2) + t351 * t253 + t324;
t663 = rSges(6,2) * t563 - t350 * t284 + t351 * t686;
t367 = qJD(1) ^ 2;
t662 = -t668 * t350 + t672 * t351;
t661 = t667 * t350 + t673 * t351;
t660 = t672 * t350 + t668 * t351;
t659 = t673 * t350 - t667 * t351;
t492 = t351 * t555;
t394 = -t353 * t511 - t492;
t483 = t351 * t506;
t445 = pkin(3) * t483;
t516 = -t357 + t365;
t525 = t284 - t347;
t603 = rSges(6,1) * t394 + t445 + (-t350 * t525 + t351 * t516) * qJD(1) + t664;
t272 = rSges(6,1) * t352 + rSges(6,2) * t353;
t222 = t272 * t350;
t606 = rSges(6,1) * t353;
t274 = -rSges(6,2) * t352 + t606;
t325 = qJD(5) * t351;
t300 = t350 * t495;
t520 = t365 * t511 + t300;
t602 = -t358 * t222 + t253 * t350 - t325 + t520 + (-t350 * t357 + t335 + (t274 + t525) * t351) * qJD(1);
t550 = rSges(6,1) * t562 - t524 - t663;
t293 = t351 * t347;
t549 = t350 * t516 - t293 + t665;
t648 = -t669 * t352 + t671 * t353;
t647 = t699 * qJD(1) - t691 * t255 + t690 * t256;
t644 = 0.2e1 * qJD(3);
t236 = t274 * t358;
t444 = qJD(1) * t358;
t248 = t350 * t444;
t364 = cos(qJ(1));
t356 = t364 * pkin(1);
t500 = t367 * t356;
t553 = t363 * qJD(3) ^ 2;
t388 = -pkin(3) * t351 * t553 + qJD(1) * t300 - t500;
t342 = t350 * pkin(6);
t608 = pkin(2) - t347;
t143 = (-t351 * t608 - t342) * qJD(1) - t520;
t259 = t351 * pkin(2) + t342;
t252 = t259 * qJD(1);
t551 = -t143 - t252;
t26 = -t236 * t256 + t248 * t272 + (t248 * t352 - t256 * t554) * pkin(4) + (t325 + t551 - t602) * qJD(1) + t388;
t447 = -t350 * t365 + t293;
t177 = t447 - t259;
t472 = t259 + t356;
t439 = t177 + t472;
t612 = pkin(4) * t352;
t469 = -t272 - t612;
t61 = -t300 - t325 + t469 * t255 + (t439 + t549) * qJD(1);
t582 = qJD(1) * t61;
t643 = t26 + t582;
t343 = t351 * pkin(6);
t258 = pkin(2) * t350 - t343;
t176 = t258 - t524;
t254 = qJD(1) * t258;
t642 = qJD(1) * t176 - t254;
t337 = t350 * rSges(4,3);
t556 = t351 * t363;
t557 = t351 * t361;
t207 = rSges(4,1) * t556 - rSges(4,2) * t557 + t337;
t641 = t207 + t472;
t640 = -rSges(5,2) * t563 - t351 * rSges(5,3);
t639 = t324 - t445;
t467 = t351 * rSges(3,1) - rSges(3,2) * t350;
t354 = Icges(4,4) * t363;
t425 = -Icges(4,2) * t361 + t354;
t306 = Icges(4,1) * t361 + t354;
t638 = t356 + t467;
t637 = rSges(5,1) * t555 + rSges(5,2) * t554;
t273 = rSges(5,1) * t352 + rSges(5,2) * t353;
t223 = t273 * t350;
t225 = t273 * t351;
t604 = rSges(5,2) * t352;
t275 = rSges(5,1) * t353 - t604;
t193 = t498 + t640;
t336 = t350 * rSges(5,3);
t522 = rSges(5,1) * t558 + t336;
t195 = -rSges(5,2) * t559 + t522;
t507 = qJD(3) * t351;
t508 = qJD(3) * t350;
t482 = -t176 * t508 + t177 * t507 + qJD(2);
t65 = t193 * t255 + t195 * t256 + t482;
t396 = -t256 * t273 - t445;
t473 = -t258 - t614;
t440 = t176 + t473;
t68 = (-t193 + t440) * qJD(1) + t396;
t69 = -t255 * t273 - t300 + (t195 + t439) * qJD(1);
t635 = -(qJD(1) * t223 - t256 * t275) * t68 - t65 * (-t255 * t223 - t225 * t256) - t69 * (-qJD(1) * t225 - t255 * t275);
t560 = t350 * t363;
t561 = t350 * t361;
t586 = Icges(4,3) * t351;
t198 = Icges(4,5) * t560 - Icges(4,6) * t561 - t586;
t315 = Icges(4,4) * t561;
t595 = Icges(4,5) * t351;
t202 = Icges(4,1) * t560 - t315 - t595;
t589 = Icges(4,6) * t351;
t200 = Icges(4,4) * t560 - Icges(4,2) * t561 - t589;
t573 = t200 * t361;
t417 = -t202 * t363 + t573;
t80 = -t198 * t351 - t350 * t417;
t303 = Icges(4,5) * t363 - Icges(4,6) * t361;
t302 = Icges(4,5) * t361 + Icges(4,6) * t363;
t397 = qJD(3) * t302;
t598 = Icges(4,4) * t361;
t307 = Icges(4,1) * t363 - t598;
t203 = Icges(4,5) * t350 + t307 * t351;
t201 = Icges(4,6) * t350 + t351 * t425;
t572 = t201 * t361;
t416 = -t203 * t363 + t572;
t630 = -t351 * t397 + (-t303 * t350 + t416 + t586) * qJD(1);
t199 = Icges(4,3) * t350 + t303 * t351;
t513 = qJD(1) * t199;
t629 = qJD(1) * t417 - t350 * t397 + t513;
t304 = Icges(4,2) * t363 + t598;
t412 = t304 * t361 - t306 * t363;
t626 = t412 * qJD(1) + t303 * qJD(3);
t533 = -Icges(4,2) * t560 + t202 - t315;
t535 = t306 * t350 + t200;
t625 = -t361 * t533 - t363 * t535;
t622 = t248 / 0.2e1;
t249 = t351 * t444;
t621 = t249 / 0.2e1;
t620 = -t255 / 0.2e1;
t619 = t255 / 0.2e1;
t618 = -t256 / 0.2e1;
t617 = t256 / 0.2e1;
t616 = t350 / 0.2e1;
t615 = -t351 / 0.2e1;
t613 = pkin(3) * t361;
t610 = -qJD(1) / 0.2e1;
t609 = qJD(1) / 0.2e1;
t607 = rSges(4,1) * t363;
t41 = t255 * t550 + t256 * t549 + t482;
t583 = qJD(1) * t41;
t517 = rSges(4,2) * t561 + t351 * rSges(4,3);
t206 = rSges(4,1) * t560 - t517;
t308 = rSges(4,1) * t361 + rSges(4,2) * t363;
t484 = t308 * t507;
t103 = -t484 + (-t206 + t473) * qJD(1);
t581 = t103 * t350;
t580 = t103 * t351;
t486 = t308 * t508;
t104 = qJD(1) * t641 - t486;
t246 = t308 * t351;
t579 = t104 * t246;
t571 = t255 * t353;
t565 = t302 * t350;
t564 = t302 * t351;
t543 = -t350 * t176 + t351 * t177;
t542 = t350 * t193 + t351 * t195;
t541 = -t350 * t198 - t202 * t556;
t540 = t350 * t199 + t203 * t556;
t534 = -t306 * t351 - t201;
t532 = -t304 * t351 + t203;
t281 = pkin(4) * t489;
t530 = t272 * t511 + t281;
t526 = rSges(5,2) * t489 + rSges(5,3) * t510;
t509 = qJD(1) * t361;
t521 = rSges(4,2) * t350 * t509 + rSges(4,3) * t510;
t519 = -t304 + t307;
t518 = t306 + t425;
t512 = qJD(1) * t303;
t505 = qJD(3) * t363;
t137 = -t350 * t412 - t564;
t503 = t137 * qJD(1);
t502 = qJD(1) * qJD(3);
t501 = t367 * t614;
t118 = rSges(5,1) * t394 - rSges(5,2) * t491 + t526;
t120 = -t358 * t223 + (t275 * t351 + t336) * qJD(1);
t499 = t351 * t118 + t350 * t120 + t193 * t510;
t494 = pkin(3) * t505;
t493 = t350 * t553;
t323 = pkin(6) * t510;
t142 = -t445 - t323 + (t350 * t608 - t322) * qJD(1);
t490 = t351 * t142 + t350 * t143 - t176 * t510;
t487 = t351 * t509;
t481 = -pkin(2) - t607;
t480 = t351 * t502;
t479 = t511 / 0.2e1;
t478 = t510 / 0.2e1;
t477 = -t508 / 0.2e1;
t474 = t507 / 0.2e1;
t471 = -t273 - t613;
t298 = -t612 - t613;
t470 = t298 + t613;
t466 = t361 * (-t350 ^ 2 - t351 ^ 2);
t224 = t272 * t351;
t457 = -t255 * t222 - t224 * t256;
t173 = t203 * t560;
t456 = t199 * t351 - t173;
t453 = -t198 + t572;
t452 = -qJD(1) * t224 - t255 * t274;
t443 = t350 * t550 + t351 * t549;
t442 = qJD(1) * (-pkin(2) * t511 + t323) - t501;
t441 = -pkin(4) * t554 - t236;
t237 = t275 * t358;
t436 = -t237 - t494;
t433 = qJD(1) * t222 + (-t274 - t611) * t256;
t257 = rSges(3,1) * t350 + rSges(3,2) * t351;
t430 = -rSges(4,2) * t361 + t607;
t429 = -t350 * t69 - t351 * t68;
t422 = -t104 * t350 - t580;
t121 = t200 * t363 + t202 * t361;
t122 = t201 * t363 + t203 * t361;
t415 = t206 * t350 + t207 * t351;
t411 = qJD(1) * t142 + t442;
t410 = t350 * t602 + t351 * t603 + t510 * t550;
t409 = -t272 + t298;
t245 = t308 * t350;
t81 = -t201 * t561 - t456;
t405 = (t350 * t81 - t351 * t80) * qJD(3);
t82 = -t200 * t557 - t541;
t83 = -t201 * t557 + t540;
t404 = (t350 * t83 - t351 * t82) * qJD(3);
t399 = qJD(3) * t306;
t398 = qJD(3) * t304;
t393 = t441 - t494;
t390 = -t177 * t350 * t502 + t142 * t507 + t143 * t508 - t176 * t480;
t389 = -t361 * t532 + t363 * t534;
t387 = (-t361 * t518 + t363 * t519) * qJD(1);
t386 = t256 * t469 + t639;
t139 = -rSges(4,2) * t351 * t505 + (-t363 * t511 - t483) * rSges(4,1) + t521;
t140 = -qJD(3) * t245 + (t351 * t430 + t337) * qJD(1);
t384 = t139 * t351 + t140 * t350 + (t206 * t351 - t207 * t350) * qJD(1);
t134 = qJD(1) * t201 - t350 * t398;
t136 = qJD(1) * t203 - t350 * t399;
t373 = qJD(1) * t198 - qJD(3) * t121 - t134 * t361 + t136 * t363;
t133 = -t351 * t398 + (-t350 * t425 + t589) * qJD(1);
t135 = -t351 * t399 + (-t307 * t350 + t595) * qJD(1);
t372 = -qJD(3) * t122 - t133 * t361 + t135 * t363 + t513;
t277 = t425 * qJD(3);
t278 = t307 * qJD(3);
t371 = qJD(1) * t302 - t277 * t361 + t278 * t363 + (-t304 * t363 - t306 * t361) * qJD(3);
t370 = (t350 * t657 - t351 * t658) * t622 + (t350 * t655 - t351 * t656) * t621 + (t350 * t647 + t351 * t648) * t620 + (t662 * t351 + t661 * t350 + (t350 * t656 + t351 * t655) * qJD(1)) * t619 + (t660 * t351 + t659 * t350 + (t350 * t658 + t351 * t657) * qJD(1)) * t618 + (t350 * t648 - t351 * t647) * t617 + (t678 * qJD(1) + t656 * t248 + t655 * t249 + t661 * t255 + t662 * t256) * t616 + (t677 * qJD(1) + t658 * t248 + t657 * t249 + t659 * t255 + t660 * t256) * t615 + (t671 * t352 + t669 * t353) * t610 + (t680 * t351 + t679 * t350 + (t676 * t350 + t675 * t351) * qJD(1)) * t609 + t682 * t479 + t681 * t478;
t282 = t430 * qJD(3);
t205 = t470 * t351;
t204 = t470 * t350;
t138 = -t351 * t412 + t565;
t127 = t138 * qJD(1);
t98 = qJD(3) * t415 + qJD(2);
t79 = -t500 - t282 * t507 + (-t140 - t252 + t486) * qJD(1);
t78 = -t282 * t508 + (t139 - t484) * qJD(1) + t442;
t67 = t371 * t350 - t351 * t626;
t66 = t350 * t626 + t371 * t351;
t64 = -qJD(3) * t416 + t133 * t363 + t135 * t361;
t63 = -t417 * qJD(3) + t134 * t363 + t136 * t361;
t62 = t384 * qJD(3);
t60 = (t440 - t550) * qJD(1) + t386;
t49 = -t237 * t256 + t248 * t273 + (-t120 + t551) * qJD(1) + t388;
t48 = qJD(1) * t118 - t237 * t255 - t249 * t273 + (-t361 * t480 - t493) * pkin(3) + t411;
t43 = t127 + t404;
t42 = t405 + t503;
t25 = -pkin(3) * t493 - t236 * t255 - t249 * t272 + (-t249 * t352 - t255 * t554) * pkin(4) + (t603 + t639) * qJD(1) + t411;
t16 = t118 * t256 + t120 * t255 + t193 * t249 - t195 * t248 + t390;
t9 = -t248 * t549 + t249 * t550 + t255 * t602 + t256 * t603 + t390;
t1 = [m(3) * ((-t257 * t367 - t501) * t638 + (-t500 + (-0.2e1 * t467 - t356 + t638) * t367) * (-t257 - t614)) + (t127 + ((t81 - t173 + (t199 + t573) * t351 + t541) * t351 + t540 * t350) * qJD(3)) * t474 + (((t674 + t720) * t351 + t657 + t703 + t705) * t256 + (t658 - t683) * t255 + t687) * t617 + (t42 - t503 + ((t351 * t453 - t540 + t83) * t351 + (t350 * t453 + t456 + t82) * t350) * qJD(3)) * t477 + (t66 + t64) * t508 / 0.2e1 + (t26 * (-t614 + t663) + t60 * t325 + t25 * (t356 + t665) + (-t26 * t606 + t60 * (rSges(6,1) * t555 + rSges(6,2) * t554 - t253) - t25 * t357) * t350 - (-t550 - t614) * t582 + (-rSges(6,1) * t492 - t386 + t60 - t642 + t664) * t61) * m(6) + (-t412 * qJD(3) + t277 * t363 + t278 * t361 + m(6) * ((-t362 * t61 - t364 * t60) * pkin(1) + (t60 * (-t274 - t284) - t61 * t357) * t351 + (-t60 * t686 + t61 * (-t284 - t606)) * t350) + t710 * t353 + t711 * t352) * qJD(1) + (t49 * (-t640 + t697) + (-t604 * t351 + t356 + t447 + t522) * t48 + (t637 * t350 + t520 + (-t336 - t356 + (-t275 - t347) * t351) * qJD(1)) * t68 + (t526 + (-t495 - t637) * t351 + t68 - t396 - t642 + (t193 + t614 + t697) * qJD(1)) * t69) * m(5) + (-(-t484 - t103 - t254 + (-t206 - t614) * qJD(1)) * t104 + t79 * (t350 * t481 + t343 + t517 - t614) + t78 * t641 + t104 * (t323 + t521) + (t308 * t581 - t579) * qJD(3) + ((-t103 * t364 - t104 * t362) * pkin(1) + (-pkin(2) - t430) * t580 + (t103 * (-rSges(4,3) - pkin(6)) + t104 * t481) * t350) * qJD(1)) * m(4) - (t63 + t67 + t43) * t507 / 0.2e1 + (t676 - t708) * t622 + (t675 + t707) * t621 + ((t655 + t683) * t256 + ((t674 + t696) * t351 + t670 * t350 + t656 - t716) * t255 + t682 + t685) * t620 + (t678 + t679) * t619 + (t677 - t680 + t681) * t618 + ((t121 + t137) * t350 + (t122 + t138) * t351) * t502 / 0.2e1; m(4) * t62 + m(5) * t16 + m(6) * t9; ((-t507 * t565 - t512) * t351 + (t387 + (t389 * t350 + (t564 - t625) * t351) * qJD(3)) * t350) * t474 + ((t361 * t519 + t363 * t518) * qJD(1) + ((t350 * t532 - t351 * t533) * t363 + (t350 * t534 + t351 * t535) * t361) * qJD(3)) * t610 + ((-t508 * t564 + t512) * t350 + (t387 + (-t625 * t351 + (t565 + t389) * t350) * qJD(3)) * t351) * t477 + (t350 * t64 - t351 * t63 + (t121 * t350 + t122 * t351) * qJD(1)) * t609 + t370 + (qJD(1) * t66 + (-(t350 * t629 + t373 * t351) * t351 + (t350 * t630 + t372 * t351) * t350 + (t82 * t350 + t83 * t351) * qJD(1)) * t644) * t616 + (qJD(1) * t67 + (-(t373 * t350 - t351 * t629) * t351 + (t372 * t350 - t351 * t630) * t350 + (t80 * t350 + t81 * t351) * qJD(1)) * t644) * t615 + (t405 + t42) * t479 + (t404 + t43) * t478 + (-t60 * (-qJD(1) * t204 + t433) - t61 * (-pkin(4) * t571 + qJD(1) * t205 + t452) - t41 * (t204 * t255 + t205 * t256 + t457) - (-t61 * t487 + ((-t350 * t61 - t351 * t60) * t363 + t41 * t466) * qJD(3)) * pkin(3) + t60 * t530 + t9 * (t443 + t543) + t41 * (t410 + t490) + (t393 * t60 + t409 * t643) * t351 + (t25 * t409 + t61 * t393 + (-t177 - t549) * t583) * t350) * m(6) + (-(-t69 * t487 + (t363 * t429 + t466 * t65) * qJD(3)) * pkin(3) + t16 * (t542 + t543) + t65 * (t490 + t499) + (t436 * t68 + (qJD(1) * t69 + t49) * t471) * t351 + (t48 * t471 + t69 * t436 + (t68 * t273 + t65 * (-t177 - t195)) * qJD(1)) * t350 + t635) * m(5) + (t62 * t415 + t98 * t384 + t422 * t282 + (-t78 * t350 - t79 * t351 + (-t104 * t351 + t581) * qJD(1)) * t308 - (t103 * t245 - t579) * qJD(1) - (t98 * (-t245 * t350 - t246 * t351) + t422 * t430) * qJD(3)) * m(4); t370 + (t9 * t443 + t41 * t410 + (t441 * t61 - t549 * t583) * t350 - t61 * t452 - t41 * t457 - (-t61 * t571 + (-t61 * t510 + t41 * (-t255 * t350 - t256 * t351)) * t352) * pkin(4) + (t25 * t350 + t351 * t643) * t469 + (t351 * t441 - t281 - t433 + t530) * t60) * m(6) + (t16 * t542 + t65 * (-t195 * t511 + t499) + t429 * t237 + (-t48 * t350 - t49 * t351 + (t350 * t68 - t351 * t69) * qJD(1)) * t273 + t635) * m(5); 0.2e1 * (t25 * t615 + t26 * t616 - t41 * (-t255 * t351 + t256 * t350) / 0.2e1) * m(6);];
tauc = t1(:);
