% Calculate vector of inverse dynamics joint torques for
% S5RPRPR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:23
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:22:25
% EndTime: 2022-01-23 09:23:05
% DurationCPUTime: 35.51s
% Computational Cost: add. (18057->786), mult. (14312->974), div. (0->0), fcn. (11121->10), ass. (0->427)
t707 = Icges(4,3) + Icges(5,3);
t366 = qJ(3) + pkin(9);
t356 = sin(t366);
t358 = cos(t366);
t255 = Icges(5,5) * t358 - Icges(5,6) * t356;
t369 = sin(qJ(3));
t371 = cos(qJ(3));
t295 = Icges(4,5) * t371 - Icges(4,6) * t369;
t700 = t255 + t295;
t367 = qJ(1) + pkin(8);
t359 = cos(t367);
t706 = t707 * t359;
t357 = sin(t367);
t557 = t357 * t371;
t558 = t357 * t369;
t559 = t357 * t358;
t561 = t356 * t357;
t687 = -Icges(4,5) * t557 - Icges(5,5) * t559 + Icges(4,6) * t558 + Icges(5,6) * t561 + t706;
t696 = t707 * t357 + t700 * t359;
t584 = Icges(5,6) * t359;
t180 = Icges(5,4) * t559 - Icges(5,2) * t561 - t584;
t585 = Icges(4,6) * t359;
t197 = Icges(4,4) * t557 - Icges(4,2) * t558 - t585;
t705 = t180 * t356 + t197 * t369;
t593 = Icges(5,4) * t356;
t259 = Icges(5,1) * t358 - t593;
t183 = Icges(5,5) * t357 + t259 * t359;
t594 = Icges(4,4) * t369;
t299 = Icges(4,1) * t371 - t594;
t200 = Icges(4,5) * t357 + t299 * t359;
t704 = -t183 * t559 - t200 * t557;
t284 = Icges(5,4) * t561;
t590 = Icges(5,5) * t359;
t182 = Icges(5,1) * t559 - t284 - t590;
t307 = Icges(4,4) * t558;
t591 = Icges(4,5) * t359;
t199 = Icges(4,1) * t557 - t307 - t591;
t692 = -t182 * t358 - t199 * t371 + t705;
t256 = Icges(5,2) * t358 + t593;
t337 = Icges(5,4) * t358;
t258 = Icges(5,1) * t356 + t337;
t296 = Icges(4,2) * t371 + t594;
t361 = Icges(4,4) * t371;
t298 = Icges(4,1) * t369 + t361;
t699 = t256 * t356 - t258 * t358 + t296 * t369 - t298 * t371;
t703 = t696 * t359 + t704;
t553 = t359 * t371;
t556 = t358 * t359;
t672 = t183 * t556 + t200 * t553 + t696 * t357;
t702 = -t182 * t556 - t199 * t553 + t687 * t357;
t657 = -t692 * t357 + t687 * t359;
t430 = -Icges(5,2) * t356 + t337;
t181 = Icges(5,6) * t357 + t359 * t430;
t431 = -Icges(4,2) * t369 + t361;
t198 = Icges(4,6) * t357 + t359 * t431;
t656 = -t181 * t561 - t198 * t558 - t703;
t554 = t359 * t369;
t560 = t356 * t359;
t655 = -t180 * t560 - t197 * t554 - t702;
t654 = -t181 * t560 - t198 * t554 + t672;
t254 = Icges(5,5) * t356 + Icges(5,6) * t358;
t294 = Icges(4,5) * t369 + Icges(4,6) * t371;
t701 = t254 + t294;
t698 = -t256 * t358 - t258 * t356 - t296 * t371 - t298 * t369;
t566 = t294 * t359;
t568 = t254 * t359;
t677 = -t699 * t357 - t566 - t568;
t567 = t294 * t357;
t569 = t254 * t357;
t676 = -t699 * t359 + t567 + t569;
t697 = t181 * t356 + t198 * t369;
t234 = t430 * qJD(3);
t235 = t259 * qJD(3);
t276 = t431 * qJD(3);
t277 = t299 * qJD(3);
t695 = t701 * qJD(1) + t698 * qJD(3) - t234 * t356 + t235 * t358 - t276 * t369 + t277 * t371;
t494 = rSges(5,1) * t559;
t368 = -qJ(4) - pkin(6);
t317 = t359 * t368;
t362 = t371 * pkin(3);
t352 = t362 + pkin(2);
t521 = t357 * t352 + t317;
t370 = sin(qJ(1));
t617 = pkin(1) * t370;
t694 = -t494 - t617 - t521;
t693 = t183 * t358 + t200 * t371 - t697;
t691 = t654 * t357 - t655 * t359;
t690 = t656 * t357 - t657 * t359;
t689 = t699 * qJD(1) + t700 * qJD(3);
t372 = cos(qJ(1));
t363 = t372 * pkin(1);
t653 = t180 * t358 + t182 * t356 + t197 * t371 + t199 * t369;
t652 = t181 * t358 + t183 * t356 + t198 * t371 + t200 * t369;
t688 = t676 * qJD(1);
t263 = rSges(3,1) * t357 + rSges(3,2) * t359;
t230 = -t263 - t617;
t686 = t701 * qJD(3);
t685 = t677 * qJD(1);
t505 = qJD(3) * t359;
t508 = qJD(1) * t357;
t285 = rSges(5,2) * t561;
t507 = qJD(1) * t359;
t524 = rSges(5,3) * t507 + qJD(1) * t285;
t607 = rSges(5,2) * t358;
t110 = -t505 * t607 + (-t356 * t505 - t358 * t508) * rSges(5,1) + t524;
t339 = t357 * rSges(5,3);
t185 = rSges(5,1) * t556 - rSges(5,2) * t560 + t339;
t645 = t358 * rSges(5,1) - rSges(5,2) * t356;
t238 = t645 * qJD(3);
t502 = qJD(1) * qJD(3);
t248 = qJDD(3) * t357 + t359 * t502;
t262 = rSges(5,1) * t356 + t607;
t322 = pkin(6) * t507;
t326 = qJD(4) * t357;
t484 = t369 * t505;
t458 = pkin(3) * t484;
t612 = pkin(2) - t352;
t115 = -t458 - t322 + t326 + (t357 * t612 - t317) * qJD(1);
t346 = t357 * pkin(6);
t267 = t359 * pkin(2) + t346;
t515 = -t359 * t352 + t357 * t368;
t164 = -t267 - t515;
t374 = qJD(1) ^ 2;
t455 = qJDD(1) * t363 - t374 * t617;
t410 = qJD(1) * (-pkin(2) * t508 + t322) + qJDD(1) * t267 + t455;
t501 = qJD(1) * qJD(4);
t373 = qJD(3) ^ 2;
t552 = t371 * t373;
t377 = qJD(1) * t115 + qJDD(1) * t164 + t357 * t501 + (-t248 * t369 - t357 * t552) * pkin(3) - qJDD(4) * t359 + t410;
t506 = qJD(3) * t357;
t30 = qJD(1) * t110 + qJDD(1) * t185 - t238 * t506 - t248 * t262 + t377;
t684 = t30 - g(2);
t496 = t374 * t363;
t683 = t690 * qJD(3) + t685;
t682 = t691 * qJD(3) + t688;
t402 = qJD(3) * t256;
t107 = qJD(1) * t181 - t357 * t402;
t404 = qJD(3) * t258;
t109 = qJD(1) * t183 - t357 * t404;
t403 = qJD(3) * t296;
t123 = qJD(1) * t198 - t357 * t403;
t405 = qJD(3) * t298;
t125 = qJD(1) * t200 - t357 * t405;
t681 = t692 * qJD(3) - t107 * t358 - t109 * t356 - t123 * t371 - t125 * t369;
t106 = -t359 * t402 + (-t357 * t430 + t584) * qJD(1);
t108 = -t359 * t404 + (-t259 * t357 + t590) * qJD(1);
t122 = -t359 * t403 + (-t357 * t431 + t585) * qJD(1);
t124 = -t359 * t405 + (-t299 * t357 + t591) * qJD(1);
t680 = t693 * qJD(3) + t106 * t358 + t108 * t356 + t122 * t371 + t124 * t369;
t679 = t689 * t357 + t695 * t359;
t678 = t695 * t357 - t689 * t359;
t347 = t359 * pkin(6);
t266 = pkin(2) * t357 - t347;
t163 = t266 - t521;
t251 = qJD(1) * t266;
t675 = -qJD(1) * t163 + t251 + t326;
t674 = t687 + t697;
t673 = t696 * qJD(1);
t670 = t359 ^ 2;
t615 = pkin(4) * t356;
t616 = pkin(3) * t369;
t287 = -t615 - t616;
t237 = t359 * t287;
t349 = pkin(4) * t358;
t449 = -t362 - t349;
t278 = pkin(2) - t449;
t364 = pkin(7) - t368;
t311 = t364 * t359;
t132 = t278 * t357 - t311 - t521;
t360 = qJ(5) + t366;
t350 = sin(t360);
t565 = t350 * t357;
t273 = rSges(6,2) * t565;
t523 = t359 * rSges(6,3) + t273;
t351 = cos(t360);
t563 = t351 * t357;
t171 = rSges(6,1) * t563 - t523;
t669 = -t132 - t171;
t668 = -t652 * qJD(3) - t106 * t356 + t108 * t358 - t122 * t369 + t124 * t371 + t673;
t667 = t687 * qJD(1) + t653 * qJD(3) + t107 * t356 - t109 * t358 + t123 * t369 - t125 * t371;
t666 = t692 * qJD(1) - t686 * t357 + t673;
t665 = -t686 * t359 + (-t700 * t357 - t693 + t706) * qJD(1);
t338 = t357 * rSges(6,3);
t562 = t351 * t359;
t564 = t350 * t359;
t172 = rSges(6,1) * t562 - rSges(6,2) * t564 + t338;
t310 = t364 * t357;
t661 = t359 * t278 + t310;
t647 = t515 + t661;
t547 = t647 + t172;
t664 = rSges(4,2) * t369;
t260 = t287 * qJD(3);
t220 = t359 * t260;
t475 = -t266 - t617;
t454 = t163 + t475;
t409 = t454 + t669;
t605 = rSges(6,2) * t351;
t610 = rSges(6,1) * t350;
t246 = t605 + t610;
t365 = qJD(3) + qJD(5);
t253 = t359 * t365;
t459 = -t246 * t253 + t326;
t49 = qJD(1) * t409 + t220 + t459;
t663 = t357 * t49;
t592 = Icges(6,4) * t350;
t244 = Icges(6,1) * t351 - t592;
t170 = Icges(6,5) * t357 + t244 * t359;
t241 = Icges(6,2) * t351 + t592;
t662 = -t241 * t359 + t170;
t660 = -t241 + t244;
t321 = Icges(6,4) * t351;
t243 = Icges(6,1) * t350 + t321;
t429 = -Icges(6,2) * t350 + t321;
t659 = t243 + t429;
t474 = t267 + t363;
t453 = t164 + t474;
t658 = t453 + t547;
t340 = t357 * rSges(4,3);
t204 = rSges(4,1) * t553 - rSges(4,2) * t554 + t340;
t139 = t204 + t474;
t650 = t645 + t362;
t265 = t359 * rSges(3,1) - rSges(3,2) * t357;
t231 = t265 + t363;
t649 = -t359 * rSges(5,3) - t285;
t323 = t351 * rSges(6,1);
t646 = -rSges(6,2) * t350 + t323;
t644 = t666 * t670 + (t668 * t357 + (-t665 + t667) * t359) * t357;
t643 = t667 * t670 + (t665 * t357 + (-t666 + t668) * t359) * t357;
t240 = Icges(6,5) * t351 - Icges(6,6) * t350;
t239 = Icges(6,5) * t350 + Icges(6,6) * t351;
t570 = t239 * t359;
t168 = Icges(6,6) * t357 + t359 * t429;
t578 = t168 * t350;
t580 = Icges(6,3) * t359;
t640 = -t365 * t570 + (-t170 * t351 - t240 * t357 + t578 + t580) * qJD(1);
t583 = Icges(6,6) * t359;
t167 = Icges(6,4) * t563 - Icges(6,2) * t565 - t583;
t272 = Icges(6,4) * t565;
t589 = Icges(6,5) * t359;
t169 = Icges(6,1) * t563 - t272 - t589;
t427 = t167 * t350 - t169 * t351;
t166 = Icges(6,3) * t357 + t240 * t359;
t513 = qJD(1) * t166;
t571 = t239 * t357;
t639 = qJD(1) * t427 - t365 * t571 + t513;
t420 = t241 * t350 - t243 * t351;
t634 = qJD(1) * t420 + t240 * t365;
t631 = t357 * (-t256 * t359 + t183) - t359 * (-Icges(5,2) * t559 + t182 - t284);
t530 = -Icges(4,2) * t557 + t199 - t307;
t532 = t298 * t357 + t197;
t630 = -t369 * t530 - t371 * t532;
t252 = t357 * t365;
t629 = qJD(1) * t659 + t252 * t662 - t253 * (-Icges(6,2) * t563 + t169 - t272);
t628 = -m(5) - m(6);
t500 = qJD(1) * qJD(5);
t175 = qJDD(5) * t357 + t359 * t500 + t248;
t627 = t175 / 0.2e1;
t314 = t357 * t502;
t176 = t357 * t500 + t314 + (-qJDD(3) - qJDD(5)) * t359;
t626 = t176 / 0.2e1;
t625 = t248 / 0.2e1;
t249 = -qJDD(3) * t359 + t314;
t624 = t249 / 0.2e1;
t623 = -t252 / 0.2e1;
t622 = t252 / 0.2e1;
t621 = -t253 / 0.2e1;
t620 = t253 / 0.2e1;
t619 = t357 / 0.2e1;
t618 = -t359 / 0.2e1;
t614 = -qJD(1) / 0.2e1;
t613 = qJD(1) / 0.2e1;
t611 = rSges(4,1) * t371;
t301 = rSges(4,1) * t369 + rSges(4,2) * t371;
t228 = t301 * t359;
t99 = qJD(1) * t139 - t301 * t506;
t604 = t228 * t99;
t516 = rSges(4,2) * t558 + t359 * rSges(4,3);
t203 = rSges(4,1) * t557 - t516;
t452 = -t203 + t475;
t485 = t301 * t505;
t98 = qJD(1) * t452 - t485;
t603 = t357 * t98;
t602 = t359 * t49;
t601 = t359 * t98;
t399 = -t262 - t616;
t397 = t399 * t505 + t326;
t184 = t494 + t649;
t415 = -t184 + t454;
t65 = qJD(1) * t415 + t397;
t600 = t65 * t262;
t598 = qJDD(1) / 0.2e1;
t165 = Icges(6,5) * t563 - Icges(6,6) * t565 - t580;
t579 = t165 * t359;
t555 = t358 * t373;
t80 = -t357 * t420 - t570;
t550 = t80 * qJD(1);
t300 = t368 * t508;
t499 = pkin(3) * t558;
t292 = qJD(3) * t499;
t327 = qJD(4) * t359;
t519 = t292 + t327;
t488 = t300 + t519;
t116 = (-t359 * t612 - t346) * qJD(1) - t488;
t245 = t267 * qJD(1);
t548 = -t116 - t245;
t546 = -t357 * t165 - t169 * t562;
t545 = t357 * t166 + t170 * t562;
t541 = -t357 * t163 + t359 * t164;
t540 = t357 * t171 + t359 * t172;
t537 = -t164 - t185;
t531 = -t298 * t359 - t198;
t529 = -t296 * t359 + t200;
t528 = t364 * t507 + t220;
t526 = -t256 + t259;
t525 = t258 + t430;
t522 = t278 - t352;
t520 = rSges(4,3) * t507 + t508 * t664;
t518 = -t296 + t299;
t517 = t298 + t431;
t510 = qJD(1) * t255;
t509 = qJD(1) * t295;
t504 = qJD(3) * t371;
t498 = pkin(3) * t554;
t497 = pkin(3) * t552;
t493 = t365 * t605;
t408 = rSges(6,3) * t507 + qJD(1) * t273 - t359 * t493;
t491 = t350 * t253;
t95 = (-t351 * t508 - t491) * rSges(6,1) + t408;
t201 = t246 * t357;
t96 = -t365 * t201 + (t359 * t646 + t338) * qJD(1);
t495 = t171 * t507 + t357 * t96 + t359 * t95;
t492 = pkin(3) * t504;
t490 = t359 * t115 + t357 * t116 - t163 * t507;
t483 = -t163 * t506 + t164 * t505 + qJD(2);
t482 = -pkin(2) - t611;
t481 = t508 / 0.2e1;
t480 = t507 / 0.2e1;
t479 = -t506 / 0.2e1;
t478 = t506 / 0.2e1;
t477 = -t505 / 0.2e1;
t476 = t505 / 0.2e1;
t472 = -t278 - t323;
t407 = t243 * t365;
t471 = qJD(1) * t170 - t167 * t365 - t357 * t407;
t470 = -t168 * t365 - t359 * t407 + (-t244 * t357 + t589) * qJD(1);
t469 = qJD(1) * t168 + t169 * t365 - t241 * t252;
t468 = (-t357 * t429 + t583) * qJD(1) + t662 * t365;
t135 = t170 * t563;
t467 = t166 * t359 - t135;
t464 = -t165 + t578;
t461 = t659 * t365;
t460 = t660 * t365;
t457 = (-t357 ^ 2 - t670) * t616;
t456 = t115 * t505 + t116 * t506 - t248 * t163 + qJDD(2);
t450 = -t238 - t492;
t448 = t246 * t252 + t519;
t304 = rSges(2,1) * t372 - rSges(2,2) * t370;
t302 = rSges(2,1) * t370 + rSges(2,2) * t372;
t303 = t611 - t664;
t435 = -t357 * t99 - t601;
t128 = -rSges(4,2) * t359 * t504 + (-t371 * t508 - t484) * rSges(4,1) + t520;
t227 = t301 * t357;
t129 = -qJD(3) * t227 + (t303 * t359 + t340) * qJD(1);
t428 = t128 * t359 + t129 * t357;
t78 = t167 * t351 + t169 * t350;
t421 = t203 * t357 + t204 * t359;
t413 = -t246 + t287;
t411 = qJDD(4) * t357 + t249 * t616 + t359 * t501 - t496;
t211 = t262 * t357;
t202 = t246 * t359;
t406 = t427 * t357;
t216 = t646 * t365;
t396 = -qJD(3) * t349 - t216 - t492;
t395 = qJD(1) * t240 - t252 * t570 + t253 * t571;
t394 = t180 * t359 - t181 * t357;
t393 = -t369 * t529 + t371 * t531;
t392 = t359 * t399;
t391 = (-t356 * t525 + t358 * t526) * qJD(1);
t390 = (-t369 * t517 + t371 * t518) * qJD(1);
t386 = qJD(1) * t165 - t350 * t469 + t351 * t471;
t11 = t357 * t639 + t386 * t359;
t385 = -t350 * t468 + t351 * t470 + t513;
t12 = t357 * t640 + t385 * t359;
t13 = t386 * t357 - t359 * t639;
t14 = t385 * t357 - t359 * t640;
t61 = -t406 - t579;
t62 = -t168 * t565 - t467;
t23 = t252 * t62 - t253 * t61 + t550;
t63 = -t167 * t564 - t546;
t64 = -t168 * t564 + t545;
t81 = -t359 * t420 + t571;
t77 = t81 * qJD(1);
t24 = t252 * t64 - t253 * t63 + t77;
t387 = (-t243 * t359 - t168) * t252 - (-t243 * t357 - t167) * t253 + t660 * qJD(1);
t375 = -t350 * t629 + t387 * t351;
t384 = qJD(1) * t239 - t350 * t461 + t351 * t460;
t38 = t357 * t634 + t384 * t359;
t39 = t384 * t357 - t359 * t634;
t40 = t350 * t471 + t351 * t469;
t41 = t350 * t470 + t351 * t468;
t79 = t168 * t351 + t170 * t350;
t389 = (qJD(1) * t38 + qJDD(1) * t81 - t11 * t253 + t12 * t252 + t175 * t64 + t176 * t63) * t619 + (t387 * t350 + t351 * t629) * t614 + (qJD(1) * t39 + qJDD(1) * t80 - t13 * t253 + t14 * t252 + t175 * t62 + t176 * t61) * t618 + t23 * t481 + t24 * t480 + (-t11 * t359 + t12 * t357 + (t357 * t63 + t359 * t64) * qJD(1)) * t622 + (t357 * t64 - t359 * t63) * t627 + (t357 * t62 - t359 * t61) * t626 + (-t13 * t359 + t14 * t357 + (t357 * t61 + t359 * t62) * qJD(1)) * t621 + (t357 * t79 - t359 * t78) * t598 + (t357 * t41 - t359 * t40 + (t357 * t78 + t359 * t79) * qJD(1)) * t613 + (t357 * t395 + t359 * t375) * t623 + (t357 * t375 - t359 * t395) * t620;
t35 = t171 * t252 + t172 * t253 + (t132 * t357 + t359 * t647) * qJD(3) + t483;
t50 = t658 * qJD(1) - t506 * t615 - t448;
t388 = t35 * (-t252 * t201 - t202 * t253) + t50 * (-qJD(1) * t202 - t252 * t646);
t376 = -t356 * t631 + t394 * t358;
t279 = t303 * qJD(3);
t236 = t357 * t287;
t212 = t262 * t359;
t188 = t237 + t498;
t187 = t236 + t499;
t177 = t253 * t646;
t111 = -qJD(3) * t211 + (t359 * t645 + t339) * qJD(1);
t97 = qJD(3) * t421 + qJD(2);
t76 = t260 * t357 + t292 + t300 + (t359 * t522 + t310) * qJD(1);
t75 = t458 + (-t357 * t522 + t317) * qJD(1) + t528;
t66 = -t262 * t506 + (t185 + t453) * qJD(1) - t519;
t56 = (t184 * t357 + t185 * t359) * qJD(3) + t483;
t55 = qJD(1) * t128 + qJDD(1) * t204 - t248 * t301 - t279 * t506 + t410;
t54 = -t496 - t279 * t505 + t249 * t301 + (-t129 - t245) * qJD(1) + t452 * qJDD(1);
t48 = qJD(3) * t428 + t203 * t248 - t204 * t249 + qJDD(2);
t29 = t249 * t262 + (-qJD(3) * t238 - t497) * t359 + (-t111 + t548) * qJD(1) + t415 * qJDD(1) + t411;
t15 = t184 * t248 + t537 * t249 + (t110 * t359 + t111 * t357) * qJD(3) + t456;
t10 = -t175 * t246 - t216 * t252 + t547 * qJDD(1) + (t75 + t95) * qJD(1) + (-t248 * t356 - t357 * t555) * pkin(4) + t377;
t9 = -t359 * t497 + t176 * t246 - t216 * t253 + (t249 * t356 - t359 * t555) * pkin(4) + (-t76 - t96 + t548) * qJD(1) + t409 * qJDD(1) + t411;
t5 = t132 * t248 + t171 * t175 - t172 * t176 + t252 * t96 + t253 * t95 + (-t647 - t164) * t249 + (t357 * t76 + t359 * t75) * qJD(3) + t456;
t1 = [(t77 + (t62 + (t167 * t359 + t168 * t357) * t350 + t467 + t546) * t253 + (-t169 * t563 + t579 + t61 + (t167 * t357 - t168 * t359) * t350 + t545) * t252) * t620 - m(2) * (-g(1) * t302 + g(2) * t304) + (t79 + t81) * t627 + (t78 + t80) * t626 + (-t550 + (t64 - t406 - t545) * t253 + (t357 * t464 - t135 + t63) * t252 + ((t166 + t427) * t252 + t464 * t253) * t359 + t23) * t623 + (t41 + t38) * t622 + ((t672 * t357 + ((t696 + t705) * t359 + t656 + t702 + t704) * t359) * qJD(3) + t688) * t476 + ((-t263 * t374 - g(2) + t455) * t231 + (-t496 + (-0.2e1 * t265 - t363 + t231) * t374 - g(1)) * t230) * m(3) + (t39 + t24 + t40) * t621 + (-t699 * qJD(3) + t234 * t358 + t235 * t356 + t276 * t371 + t277 * t369 + t350 * t460 + t351 * t461) * qJD(1) + (-t615 * t663 * qJD(3) + (-g(2) + t10) * (t172 + t363 + t661) + (-rSges(6,1) * t491 - t237 * qJD(3) + t408 - t459 + t528 + t675) * t50 + (t327 + (t365 * t610 - t260 + t493) * t357 - t448) * t49 + ((-t370 * t50 - t372 * t49) * pkin(1) + (-t278 - t646) * t602 + (t49 * (-rSges(6,3) - t364) + t50 * t472) * t357 + t49 * t658 - t50 * (-t617 + t669)) * qJD(1) + (t9 - g(1)) * (t357 * t472 + t311 + t523 - t617)) * m(6) + (t600 * t506 + (t488 + (-t339 - t363 + (-t352 - t645) * t359) * qJD(1)) * t65 + (t392 * qJD(3) - t397 + t524 + t65 + (t184 + t617 + t694) * qJD(1) + t675) * t66 + t684 * (t185 + t363 - t515) + (t29 - g(1)) * (-t649 + t694)) * m(5) + (t99 * (t322 + t520) + (t301 * t603 - t604) * qJD(3) + ((-t370 * t99 - t372 * t98) * pkin(1) + (-pkin(2) - t303) * t601 + (t98 * (-rSges(4,3) - pkin(6)) + t99 * t482) * t357) * qJD(1) - (-t485 - t251 - t98 + (-t203 - t617) * qJD(1)) * t99 + (t55 - g(2)) * t139 + (t54 - g(1)) * (t482 * t357 + t347 + t516 - t617)) * m(4) + (t676 + t652) * t625 + (t677 + t653) * t624 + (((t359 * t674 + t654 - t672) * t359 + (t357 * t674 + t655 + t703) * t357) * qJD(3) + t683 - t685) * t479 + (t679 + t680) * t478 + (m(3) * (t230 ^ 2 + t265 * t231) + m(2) * (t302 ^ 2 + t304 ^ 2) + t241 * t351 + t243 * t350 + Icges(2,3) + Icges(3,3) - t698) * qJDD(1) + (t678 - t681 + t682) * t477; m(3) * qJDD(2) + m(4) * t48 + m(5) * t15 + m(6) * t5 + (-m(3) - m(4) + t628) * g(3); t389 + t691 * t625 + t690 * t624 + (t679 * qJD(1) + t643 * qJD(3) + t676 * qJDD(1) + t654 * t248 + t655 * t249) * t619 + (t678 * qJD(1) + t644 * qJD(3) + t677 * qJDD(1) + t656 * t248 + t657 * t249) * t618 + (((t357 * t529 - t359 * t530) * t371 + (t357 * t531 + t359 * t532) * t369 + t631 * t358 + t394 * t356) * qJD(3) + (t356 * t526 + t358 * t525 + t369 * t518 + t371 * t517) * qJD(1)) * t614 + (t681 * t359 + t680 * t357 + (t357 * t653 + t359 * t652) * qJD(1)) * t613 + (t357 * t652 - t359 * t653) * t598 + t683 * t481 + t682 * t480 + ((-t506 * t566 + t509) * t357 + (t390 + (-t630 * t359 + (t567 + t393) * t357) * qJD(3)) * t359 + (-t506 * t568 + t510) * t357 + (t391 + (t357 * t569 + t376) * qJD(3)) * t359) * t479 + ((t357 * t655 + t359 * t654) * qJD(1) + t643) * t478 + ((t357 * t657 + t656 * t359) * qJD(1) + t644) * t477 + ((-t505 * t567 - t509) * t359 + (t390 + (t393 * t357 + (t566 - t630) * t359) * qJD(3)) * t357 + (-t505 * t569 - t510) * t359 + (t391 + (t359 * t568 + t376) * qJD(3)) * t357) * t476 + (t5 * (t540 + t541) + t35 * (t490 + t495) + (t9 * t413 + t49 * t396 + t5 * t647 + t35 * t75 + (t35 * t132 + t413 * t50) * qJD(1)) * t359 + (t10 * t413 + t50 * t396 + t5 * t132 + t35 * t76 + (t49 * (t246 + t615) + t35 * (-t164 - t547)) * qJD(1)) * t357 - g(1) * (t237 - t202) - g(2) * (t236 - t201) - g(3) * (t646 - t449) + t49 * t177 - (t49 * (-t187 + t201) + t50 * (t188 - t498)) * qJD(1) - (t35 * t457 + (t35 * t188 + t449 * t49) * t359 + (t35 * t187 + t449 * t50) * t357) * qJD(3) - t388) * m(6) + (t48 * t421 + t97 * ((t203 * t359 - t204 * t357) * qJD(1) + t428) + t435 * t279 + (-t55 * t357 - t54 * t359 + (-t359 * t99 + t603) * qJD(1)) * t301 - (t227 * t98 - t604) * qJD(1) - (t97 * (-t227 * t357 - t228 * t359) + t435 * t303) * qJD(3) + g(1) * t228 + g(2) * t227 - g(3) * t303) * m(4) + (-(t65 * t211 + t66 * (-t212 - t498)) * qJD(1) - (t56 * t457 + (-t56 * t212 - t65 * t650) * t359) * qJD(3) + t15 * t541 + t56 * t490 + (t29 * t399 + t65 * t450 + t15 * t185 + t56 * t110 + (t56 * t184 + t399 * t66) * qJD(1)) * t359 - g(3) * t650 - g(1) * t392 + (-(-t56 * t211 - t650 * t66) * qJD(3) + t66 * t450 + t15 * t184 + t56 * t111 + (t537 * t56 + t600) * qJD(1) + t684 * t399) * t357) * m(5); t628 * (g(1) * t357 - g(2) * t359) + 0.2e1 * (t10 * t618 + t619 * t9) * m(6) + 0.2e1 * (t29 * t619 + t30 * t618) * m(5); t389 + (-t49 * (qJD(1) * t201 - t177) - t388 + t5 * t540 + t35 * (-t172 * t508 + t495) + (-t357 * t50 - t602) * t216 + (-t10 * t357 - t9 * t359 + (-t359 * t50 + t663) * qJD(1)) * t246 + g(1) * t202 + g(2) * t201 - g(3) * t646) * m(6);];
tau = t1;
