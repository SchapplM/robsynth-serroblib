% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:55:55
% EndTime: 2019-12-05 17:56:39
% DurationCPUTime: 21.00s
% Computational Cost: add. (46095->591), mult. (50999->825), div. (0->0), fcn. (54935->10), ass. (0->360)
t733 = Icges(4,3) + Icges(5,3);
t444 = cos(qJ(1));
t441 = cos(pkin(8));
t634 = cos(qJ(3));
t525 = t634 * pkin(3);
t478 = t441 * (t525 + pkin(2));
t442 = sin(qJ(3));
t443 = sin(qJ(1));
t584 = t443 * t442;
t535 = -pkin(3) * t584 - t444 * t478;
t531 = qJ(3) + pkin(9);
t498 = cos(t531);
t493 = pkin(4) * t498;
t418 = t493 + t525;
t452 = t441 * (pkin(2) + t418);
t436 = sin(t531);
t615 = pkin(4) * t436;
t616 = pkin(3) * t442;
t494 = t615 + t616;
t539 = -t443 * t494 - t444 * t452;
t440 = sin(pkin(8));
t591 = t440 * t444;
t233 = -pkin(7) * t591 - t535 + t539;
t499 = qJ(5) + t531;
t496 = cos(t499);
t472 = t443 * t496;
t433 = sin(t499);
t583 = t444 * t433;
t367 = t441 * t583 - t472;
t471 = t444 * t496;
t586 = t443 * t433;
t368 = t441 * t471 + t586;
t677 = -t368 * rSges(6,1) + t367 * rSges(6,2);
t264 = rSges(6,3) * t591 - t677;
t246 = t441 * t264;
t614 = t441 * pkin(2);
t311 = (-qJ(4) * t440 + t614) * t444 + t535;
t290 = t441 * t311;
t446 = t440 * (-t441 * rSges(6,3) + (rSges(6,1) * t496 - rSges(6,2) * t433) * t440);
t332 = t444 * t446;
t457 = t440 * (-t441 * qJ(4) + t440 * t525);
t334 = t444 * t457;
t454 = t440 * (-pkin(7) * t441 + t440 * t493);
t515 = t444 * t454 + t332 + t334;
t103 = -t233 * t441 + t246 - t290 + t515;
t724 = t264 - t233 - t311;
t580 = t724 * t441 - t103 + t515;
t731 = m(6) * t580;
t437 = t444 * qJ(2);
t613 = -qJ(4) - pkin(6);
t530 = -pkin(7) + t613;
t592 = t440 * t443;
t536 = t444 * t494 + t530 * t592;
t365 = t441 * t586 + t471;
t366 = t441 * t472 - t583;
t545 = t366 * rSges(6,1) - t365 * rSges(6,2);
t193 = -t437 + (t440 * rSges(6,3) + pkin(1) + t452) * t443 - t536 + t545;
t587 = t443 * qJ(2);
t194 = (-pkin(1) + (-rSges(6,3) + t530) * t440) * t444 - t587 + t539 + t677;
t288 = t365 * rSges(6,1) + t366 * rSges(6,2);
t289 = -t367 * rSges(6,1) - t368 * rSges(6,2);
t732 = m(6) * (-t193 * t288 - t194 * t289);
t604 = Icges(5,4) * t436;
t351 = -Icges(5,5) * t441 + (Icges(5,1) * t498 - t604) * t440;
t606 = Icges(4,4) * t442;
t382 = -Icges(4,5) * t441 + (t634 * Icges(4,1) - t606) * t440;
t384 = (-Icges(5,2) * t498 - t604) * t440;
t406 = (-t634 * Icges(4,2) - t606) * t440;
t487 = t498 * Icges(5,4);
t385 = (-Icges(5,1) * t436 - t487) * t440;
t485 = t498 * t385;
t350 = -Icges(5,6) * t441 + (-Icges(5,2) * t436 + t487) * t440;
t486 = t498 * t350;
t511 = t634 * Icges(4,4);
t407 = (-Icges(4,1) * t442 - t511) * t440;
t509 = t634 * t407;
t381 = -Icges(4,6) * t441 + (-Icges(4,2) * t442 + t511) * t440;
t510 = t634 * t381;
t730 = ((t351 / 0.2e1 + t384 / 0.2e1) * t436 + (t382 / 0.2e1 + t406 / 0.2e1) * t442 + t510 / 0.2e1 - t509 / 0.2e1 + t486 / 0.2e1 - t485 / 0.2e1) * t440;
t691 = t440 ^ 2;
t489 = t443 * t498;
t582 = t444 * t436;
t388 = t441 * t582 - t489;
t488 = t444 * t498;
t585 = t443 * t436;
t389 = t441 * t488 + t585;
t676 = -t389 * rSges(5,1) + t388 * rSges(5,2);
t277 = rSges(5,3) * t591 - t676;
t447 = t440 * (-t441 * rSges(5,3) + (rSges(5,1) * t498 - rSges(5,2) * t436) * t440);
t550 = t444 * t447 + t334;
t162 = t277 * t441 - t290 + t550;
t721 = t277 - t311;
t577 = t721 * t441 - t162 + t550;
t661 = m(5) / 0.2e1;
t717 = t577 * t661;
t204 = t246 + t332;
t668 = t443 ^ 2;
t727 = t440 * t668;
t659 = m(6) / 0.2e1;
t726 = t580 * t659;
t255 = Icges(6,5) * t368 - Icges(6,6) * t367 + Icges(6,3) * t591;
t603 = Icges(6,4) * t366;
t257 = Icges(6,2) * t365 - Icges(6,6) * t592 - t603;
t352 = Icges(6,4) * t365;
t260 = -Icges(6,1) * t366 - Icges(6,5) * t592 + t352;
t572 = -t365 * t257 + t366 * t260;
t725 = t255 * t591 - t572;
t386 = t441 * t585 + t488;
t469 = t441 * t489;
t387 = t469 - t582;
t512 = t444 * t634;
t408 = t441 * t584 + t512;
t513 = t443 * t634;
t581 = t444 * t442;
t409 = t441 * t513 - t581;
t723 = -Icges(4,5) * t409 - Icges(5,5) * t387 + Icges(4,6) * t408 + Icges(5,6) * t386 - t592 * t733;
t520 = t441 * t581;
t410 = -t513 + t520;
t411 = t441 * t512 + t584;
t722 = Icges(4,5) * t411 + Icges(5,5) * t389 - Icges(4,6) * t410 - Icges(5,6) * t388 + t733 * t591;
t605 = Icges(5,4) * t387;
t270 = Icges(5,2) * t386 - Icges(5,6) * t592 - t605;
t370 = Icges(5,4) * t386;
t273 = -Icges(5,1) * t387 - Icges(5,5) * t592 + t370;
t607 = Icges(4,4) * t409;
t304 = Icges(4,2) * t408 - Icges(4,6) * t592 - t607;
t396 = Icges(4,4) * t408;
t307 = -Icges(4,1) * t409 - Icges(4,5) * t592 + t396;
t714 = -t386 * t270 + t387 * t273 - t408 * t304 + t409 * t307;
t720 = t722 * t591 - t714;
t675 = -t411 * rSges(4,1) + t410 * rSges(4,2);
t313 = rSges(4,3) * t591 - t675;
t450 = t440 * (-t441 * rSges(4,3) + (t634 * rSges(4,1) - rSges(4,2) * t442) * t440);
t224 = t313 * t441 + t444 * t450;
t690 = -t440 / 0.2e1;
t468 = t441 * t494;
t213 = t444 * t418 + t443 * t468 + t288;
t719 = t213 * t444;
t602 = Icges(6,4) * t433;
t346 = -Icges(6,5) * t441 + (Icges(6,1) * t496 - t602) * t440;
t548 = t346 + (-Icges(6,2) * t496 - t602) * t440;
t718 = t548 * t433;
t372 = Icges(5,4) * t389;
t271 = -Icges(5,2) * t388 + Icges(5,6) * t591 + t372;
t371 = Icges(5,4) * t388;
t275 = -Icges(5,1) * t389 - Icges(5,5) * t591 + t371;
t398 = Icges(4,4) * t411;
t305 = -Icges(4,2) * t410 + Icges(4,6) * t591 + t398;
t397 = Icges(4,4) * t410;
t309 = -Icges(4,1) * t411 - Icges(4,5) * t591 + t397;
t713 = t386 * t271 + t275 * t387 + t408 * t305 + t309 * t409;
t354 = Icges(6,4) * t368;
t258 = -Icges(6,2) * t367 + Icges(6,6) * t591 + t354;
t353 = Icges(6,4) * t367;
t262 = -Icges(6,1) * t368 - Icges(6,5) * t591 + t353;
t573 = t365 * t258 + t262 * t366;
t662 = m(4) / 0.2e1;
t686 = m(6) * t440;
t709 = t686 / 0.2e1;
t254 = -Icges(6,5) * t366 + Icges(6,6) * t365 - Icges(6,3) * t592;
t112 = -t254 * t592 - t572;
t113 = -t255 * t592 + t573;
t596 = t254 * t443;
t708 = (-t255 * t727 + (-t113 + t573) * t443 + (-t112 + (-t255 * t444 - t596) * t440 + t725) * t444) * t440;
t703 = t723 * t592 + t714;
t702 = -t722 * t592 + t713;
t209 = (-pkin(1) + (-rSges(5,3) + t613) * t440) * t444 - t587 + t535 + t676;
t671 = (rSges(4,3) + pkin(6)) * t440 + pkin(1) + t614;
t239 = -t444 * t671 - t587 + t675;
t331 = t443 * t446;
t333 = t443 * t457;
t263 = -rSges(6,3) * t592 - t545;
t534 = pkin(3) * t581 + t613 * t592;
t310 = (t440 * pkin(6) - t441 * t525) * t443 + t534;
t517 = -pkin(4) * t469 + t263 + t310 - t534 + t536;
t101 = t441 * t517 - t443 * t454 - t331 - t333;
t543 = t387 * rSges(5,1) - t386 * rSges(5,2);
t556 = -rSges(5,3) * t592 + t310 - t543;
t160 = t556 * t441 - t443 * t447 - t333;
t538 = t409 * rSges(4,1) - t408 * rSges(4,2);
t312 = -rSges(4,3) * t592 - t538;
t223 = t312 * t441 - t443 * t450;
t526 = (-t223 * t444 - t224 * t443) * t662 + (-t101 * t444 - t103 * t443) * t659 + (-t160 * t444 - t162 * t443) * t661;
t529 = (t731 / 0.2e1 + t717) * t443;
t698 = t526 - t529;
t697 = t723 * t443;
t687 = m(5) * t440;
t611 = (t101 * t443 - t103 * t444) * t709 + (t160 * t443 - t162 * t444) * t687 / 0.2e1;
t612 = (t717 + t726) * t591;
t696 = t611 - t612;
t300 = -t388 * rSges(5,1) - t389 * rSges(5,2);
t421 = pkin(3) * t520;
t400 = pkin(3) * t513 - t421;
t245 = -t300 - t400;
t555 = -t386 * rSges(5,1) - t387 * rSges(5,2) - t408 * pkin(3);
t695 = t245 * t443 + t555 * t444;
t293 = Icges(5,5) * t386 + Icges(5,6) * t387;
t294 = -Icges(5,5) * t388 - Icges(5,6) * t389;
t319 = Icges(4,5) * t408 + Icges(4,6) * t409;
t320 = -Icges(4,5) * t410 - Icges(4,6) * t411;
t693 = -(t320 + t294) * t591 + (t293 + t319) * t592;
t362 = (-Icges(6,5) * t433 - Icges(6,6) * t496) * t440;
t590 = t441 * t362;
t692 = -t590 / 0.2e1 + t718 * t690;
t636 = -t441 / 0.2e1;
t506 = -t592 / 0.2e1;
t188 = (t288 * t444 + t289 * t443) * t440;
t369 = (-rSges(6,1) * t433 - rSges(6,2) * t496) * t440;
t338 = t369 * t592;
t211 = -t288 * t441 + t338;
t212 = t441 * t289 + t369 * t591;
t470 = t496 * Icges(6,4);
t345 = -Icges(6,6) * t441 + (-Icges(6,2) * t433 + t470) * t440;
t364 = (-Icges(6,1) * t433 - t470) * t440;
t549 = t345 - t364;
t117 = -t362 * t592 + t548 * t365 + t549 * t366;
t118 = t362 * t591 - t548 * t367 - t549 * t368;
t281 = Icges(6,5) * t365 + Icges(6,6) * t366;
t282 = -Icges(6,5) * t367 - Icges(6,6) * t368;
t503 = t591 / 0.2e1;
t561 = -Icges(6,2) * t368 - t262 - t353;
t562 = Icges(6,2) * t366 + t260 + t352;
t563 = Icges(6,1) * t367 + t258 + t354;
t564 = -Icges(6,1) * t365 + t257 - t603;
t159 = -t590 + (-t496 * t549 - t718) * t440;
t601 = t159 * t441;
t95 = -t441 * t281 + (-t562 * t433 - t564 * t496) * t440;
t96 = -t441 * t282 + (-t561 * t433 - t563 * t496) * t440;
t528 = (-t117 * t441 - (-t281 * t592 + t562 * t365 + t564 * t366) * t592 + (-t282 * t592 + t561 * t365 + t563 * t366) * t591) * t506 + (-t118 * t441 - (t281 * t591 - t562 * t367 - t564 * t368) * t592 + (t282 * t591 - t561 * t367 - t563 * t368) * t591) * t503 + (-t601 + (-t443 * t95 + t444 * t96) * t440) * t636;
t81 = (-t724 * t443 - t517 * t444) * t440;
t6 = t528 + m(6) * (-t101 * t211 + t103 * t212 - t81 * t188);
t683 = t6 * qJD(5);
t540 = t382 + t406;
t546 = t351 + t384;
t405 = (-Icges(4,5) * t442 - t634 * Icges(4,6)) * t440;
t588 = t441 * t405;
t383 = (-Icges(5,5) * t436 - Icges(5,6) * t498) * t440;
t589 = t441 * t383;
t680 = t588 + t589 + (t546 * t436 + t540 * t442 - t485 + t486 - t509 + t510) * t440;
t560 = -Icges(5,1) * t386 + t270 - t605;
t559 = Icges(5,1) * t388 + t271 + t372;
t554 = -Icges(4,1) * t408 + t304 - t607;
t553 = Icges(4,1) * t410 + t305 + t398;
t461 = t444 * t468;
t214 = -t443 * t418 - t289 + t461;
t578 = (-t214 * t443 + t719) * t709 - t695 * t687 / 0.2e1;
t325 = rSges(4,1) * t408 + rSges(4,2) * t409;
t326 = -rSges(4,1) * t410 - rSges(4,2) * t411;
t518 = (t443 * t213 + t214 * t444) * t659 + (t245 * t444 - t443 * t555) * t661 + (t443 * t325 - t326 * t444) * t662;
t416 = (t444 ^ 2 + t668) * t440;
t667 = 0.2e1 * t416;
t665 = 4 * qJD(1);
t664 = 2 * qJD(3);
t663 = 4 * qJD(3);
t660 = m(5) / 0.4e1;
t658 = m(6) / 0.4e1;
t655 = m(5) * t577 * t160;
t208 = -t437 + (t440 * rSges(5,3) + pkin(1) + t478) * t443 - t534 + t543;
t650 = m(5) * (t208 * t555 + t209 * t245);
t203 = t263 * t441 - t331;
t649 = t203 * t731;
t648 = t101 * t731;
t579 = -t211 * t193 + t212 * t194;
t644 = m(6) * (-t101 * t288 - t103 * t289 + t579);
t643 = m(6) * (-t203 * t213 + t204 * t214 + t579);
t638 = m(6) * (-t193 * t213 + t194 * t214);
t633 = m(3) * (-((-rSges(3,3) - qJ(2)) * t443 + rSges(3,2) * t591) * t443 - t444 * (-rSges(3,2) * t592 - t444 * rSges(3,3) - t437));
t238 = t443 * t671 - t437 + t538;
t632 = m(4) * (-t238 * t325 - t239 * t326);
t630 = m(4) * (-t444 * t238 - t239 * t443);
t628 = (t208 * t443 - t209 * t444) * t687;
t627 = m(5) * (-t444 * t208 - t209 * t443);
t624 = (t193 * t443 - t194 * t444) * t686;
t623 = m(6) * (-t444 * t193 - t194 * t443);
t622 = (t203 * t443 - t204 * t444) * t686;
t621 = m(6) * (-t203 * t444 - t204 * t443);
t618 = m(6) * t188;
t617 = m(6) * (t443 * t288 - t289 * t444);
t610 = m(6) * qJD(3);
t609 = m(6) * qJD(5);
t600 = (-(-Icges(6,3) * t441 + (Icges(6,5) * t496 - Icges(6,6) * t433) * t440) * t592 + t345 * t365 - t346 * t366) * t441;
t558 = Icges(5,2) * t387 + t273 + t370;
t557 = -Icges(5,2) * t389 - t275 - t371;
t552 = Icges(4,2) * t409 + t307 + t396;
t551 = -Icges(4,2) * t411 - t309 - t397;
t547 = t350 - t385;
t541 = t381 - t407;
t234 = (t661 + t659) * t667;
t532 = t234 * qJD(1);
t527 = t691 * t616;
t504 = -t591 / 0.2e1;
t24 = -t600 + (t573 * t444 + (t596 * t440 - t725) * t443) * t440;
t59 = -t600 + (-t112 * t443 + t113 * t444) * t440;
t497 = t24 * t503 + t59 * t504 + t708 * t506;
t473 = t440 * t496;
t459 = t473 / 0.2e1;
t460 = -t473 / 0.2e1;
t484 = t345 * t460 + t364 * t459 + t692;
t483 = -t649 / 0.2e1 + t497;
t477 = (t703 * t443 + t702 * t444) * t690 + (t713 * t444 + (t697 * t440 - t720) * t443) * t440 / 0.2e1;
t476 = (-t722 * t727 + ((-t444 * t722 - t697) * t440 + t703 + t720) * t444 + (-t702 + t713) * t443) * t690 + (t441 / 0.2e1 + t636) * (-t388 * t350 + t389 * t351 - t410 * t381 + t411 * t382 + (-t733 * t441 + (t634 * Icges(4,5) + Icges(5,5) * t498 - Icges(4,6) * t442 - Icges(5,6) * t436) * t440) * t591);
t475 = -t691 * t615 - t527;
t474 = (-rSges(5,1) * t436 - rSges(5,2) * t498) * t691 - t527;
t458 = t24 * t504 + (t117 + t95) * t506 + t708 * t592 / 0.2e1 + (t59 + t118 + t96) * t503;
t455 = t345 * t459 + t364 * t460 - t692;
t453 = t458 - t601;
t412 = (-rSges(4,1) * t442 - t634 * rSges(4,2)) * t440;
t379 = t441 * t400;
t280 = t421 - t461 + (-t525 + t418) * t443;
t243 = t441 * t326 + t412 * t591;
t242 = -t325 * t441 + t412 * t592;
t235 = (t660 + t658) * t667 - (m(5) + m(6)) * t416 / 0.2e1;
t197 = t441 * t300 + t444 * t474 + t379;
t196 = t555 * t441 + t474 * t443;
t191 = t617 / 0.2e1;
t187 = t618 / 0.2e1;
t181 = (-t263 * t444 - t264 * t443) * t440;
t168 = t695 * t440;
t164 = t405 * t591 - t540 * t410 - t541 * t411;
t163 = -t405 * t592 + t540 * t408 + t541 * t409;
t152 = t211 * t443 + t212 * t444;
t147 = t152 * t609;
t137 = t441 * t280 + t444 * t475 + t212 + t379;
t136 = -t213 * t441 + t443 * t475 + t338;
t134 = t621 / 0.2e1;
t133 = t383 * t591 - t546 * t388 - t547 * t389;
t132 = -t383 * t592 + t546 * t386 + t547 * t387;
t123 = t622 / 0.2e1;
t120 = -t441 * t320 + (-t551 * t442 - t553 * t634) * t440;
t119 = -t441 * t319 + (-t552 * t442 - t554 * t634) * t440;
t107 = (-t719 + (-t280 - t289 - t400) * t443) * t440;
t106 = -t441 * t294 + (-t557 * t436 - t498 * t559) * t440;
t105 = -t441 * t293 + (-t558 * t436 - t498 * t560) * t440;
t104 = t441 * t188 + (t211 * t444 - t212 * t443) * t440;
t99 = t104 * t609;
t74 = t624 + t628;
t71 = t484 + t732;
t53 = t623 + t627 + t630 + t633;
t49 = t134 - t617 / 0.2e1;
t48 = t191 + t134;
t47 = t191 - t621 / 0.2e1;
t45 = t643 / 0.2e1;
t43 = t123 - t618 / 0.2e1;
t42 = t187 + t123;
t41 = t187 - t622 / 0.2e1;
t35 = t644 / 0.2e1;
t28 = (-t405 / 0.2e1 - t383 / 0.2e1) * t441 + t632 + t650 + t638 - t730 + t484;
t14 = m(6) * (-t181 * t188 - t203 * t211 + t204 * t212) + t528;
t13 = t14 * qJD(5);
t12 = t526 + t529 - t518;
t11 = t518 + t698;
t10 = t518 - t698;
t9 = t578 + t696;
t8 = t611 + t612 - t578;
t7 = t578 - t696;
t4 = t45 - t644 / 0.2e1 + t483;
t3 = t35 - t643 / 0.2e1 + t483;
t2 = t35 + t45 + t649 / 0.2e1 + t453;
t1 = -t655 - t648 + (t443 * t476 + t444 * t477) * t440 + t497;
t5 = [t53 * qJD(2) + t28 * qJD(3) + t74 * qJD(4) + t71 * qJD(5), qJD(1) * t53 + qJD(3) * t11 + qJD(4) * t235 + qJD(5) * t48, t28 * qJD(1) + t11 * qJD(2) + t9 * qJD(4) + t2 * qJD(5) + (t655 / 0.4e1 + t648 / 0.4e1) * t663 + ((-t223 * t325 - t224 * t326 - t238 * t242 + t239 * t243) * t662 + (t160 * t555 + t162 * t245 - t196 * t208 + t197 * t209) * t661 + (-t101 * t213 + t103 * t214 - t136 * t193 + t137 * t194) * t659) * t664 + (t458 + (-t159 + t680) * t441 + ((t120 / 0.2e1 + t106 / 0.2e1 + t164 / 0.2e1 + t133 / 0.2e1 - t477) * t444 + (-t119 / 0.2e1 - t105 / 0.2e1 - t163 / 0.2e1 - t132 / 0.2e1 - t476) * t443) * t440) * qJD(3), qJD(1) * t74 + qJD(2) * t235 + qJD(3) * t9 + qJD(5) * t42, t71 * qJD(1) + t48 * qJD(2) + t2 * qJD(3) + t42 * qJD(4) + ((-t203 * t288 - t204 * t289 + t579) * m(6) + t453) * qJD(5); t10 * qJD(3) - t234 * qJD(4) + t47 * qJD(5) + (-t633 / 0.4e1 - t630 / 0.4e1 - t627 / 0.4e1 - t623 / 0.4e1) * t665, 0, t10 * qJD(1) + t147 + ((t242 * t443 + t243 * t444) * t662 + (t196 * t443 + t197 * t444) * t661 + (t136 * t443 + t137 * t444) * t659) * t664, -t532, t47 * qJD(1) + t152 * t610 + t147; t12 * qJD(2) + t1 * qJD(3) + t8 * qJD(4) + t3 * qJD(5) + (-t632 / 0.4e1 - t650 / 0.4e1 - t638 / 0.4e1) * t665 + (t588 / 0.2e1 + t589 / 0.2e1 + t455 - 0.2e1 * t193 * t726 - 0.2e1 * t208 * t717 + t730) * qJD(1), t12 * qJD(1), t1 * qJD(1) + (t528 + (t680 * t441 + ((t106 + t120) * t444 + (-t105 - t119) * t443) * t440) * t636 + ((t386 * t557 + t387 * t559 + t408 * t551 + t409 * t553) * t591 + (-t132 - t163) * t441 + (-t386 * t558 - t387 * t560 - t408 * t552 - t409 * t554 + t693) * t592) * t506 + ((t388 * t558 + t389 * t560 + t410 * t552 + t411 * t554) * t592 + (-t164 - t133) * t441 + (-t388 * t557 - t389 * t559 - t410 * t551 - t411 * t553 - t693) * t591) * t503) * qJD(3) + t683 + ((-t101 * t136 + t103 * t137 + t107 * t81) * t658 + ((-t721 * t443 - t556 * t444) * t440 * t168 - t160 * t196 + t162 * t197) * t660 + m(4) * (-t223 * t242 + t224 * t243 + (-t312 * t444 - t313 * t443) * t691 * (-t325 * t444 - t326 * t443)) / 0.4e1) * t663, t8 * qJD(1), t3 * qJD(1) + t6 * qJD(3) + t683; t234 * qJD(2) + t7 * qJD(3) + t41 * qJD(5) + (-t628 / 0.4e1 - t624 / 0.4e1) * t665, t532, t7 * qJD(1) + ((-t441 * t107 + (t136 * t444 - t137 * t443) * t440) * t659 + (-t441 * t168 + (t196 * t444 - t197 * t443) * t440) * t661) * t664 + t99, 0, t41 * qJD(1) + t104 * t610 + t99; (t455 - t732) * qJD(1) + t49 * qJD(2) + t4 * qJD(3) + t43 * qJD(4) + t497 * qJD(5), t49 * qJD(1), t4 * qJD(1) + ((t107 * t181 - t136 * t203 + t137 * t204) * m(6) + t528) * qJD(3) + t13, t43 * qJD(1), qJD(1) * t497 + qJD(3) * t14 + t13;];
Cq = t5;
