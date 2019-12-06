% Calculate vector of inverse dynamics joint torques for
% S5RPRPR1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-05 17:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR1_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:47:03
% EndTime: 2019-12-05 17:47:41
% DurationCPUTime: 26.85s
% Computational Cost: add. (11862->774), mult. (14954->944), div. (0->0), fcn. (11551->8), ass. (0->414)
t709 = Icges(4,3) + Icges(5,3);
t374 = qJ(3) + pkin(8);
t342 = sin(t374);
t343 = cos(t374);
t378 = sin(qJ(3));
t380 = cos(qJ(3));
t708 = Icges(4,5) * t378 + Icges(5,5) * t342 + Icges(4,6) * t380 + Icges(5,6) * t343;
t379 = sin(qJ(1));
t381 = cos(qJ(1));
t706 = t708 * t379 + t709 * t381;
t580 = Icges(5,4) * t342;
t437 = Icges(5,2) * t343 + t580;
t174 = -Icges(5,6) * t379 + t381 * t437;
t579 = Icges(5,4) * t343;
t440 = Icges(5,1) * t342 + t579;
t176 = -Icges(5,5) * t379 + t381 * t440;
t582 = Icges(4,4) * t378;
t438 = Icges(4,2) * t380 + t582;
t207 = -Icges(4,6) * t379 + t381 * t438;
t581 = Icges(4,4) * t380;
t441 = Icges(4,1) * t378 + t581;
t209 = -Icges(4,5) * t379 + t381 * t441;
t695 = t174 * t343 + t176 * t342 + t207 * t380 + t209 * t378;
t707 = t695 * t381;
t651 = -t709 * t379 + t708 * t381;
t705 = Icges(4,5) * t380 + Icges(5,5) * t343 - Icges(4,6) * t378 - Icges(5,6) * t342;
t173 = Icges(5,6) * t381 + t379 * t437;
t563 = t343 * t379;
t306 = Icges(5,4) * t563;
t565 = t342 * t379;
t575 = Icges(5,5) * t381;
t175 = Icges(5,1) * t565 + t306 + t575;
t206 = Icges(4,6) * t381 + t379 * t438;
t559 = t379 * t380;
t326 = Icges(4,4) * t559;
t561 = t378 * t379;
t576 = Icges(4,5) * t381;
t208 = Icges(4,1) * t561 + t326 + t576;
t647 = t173 * t343 + t175 * t342 + t206 * t380 + t208 * t378;
t652 = t706 * t379;
t257 = -Icges(5,2) * t342 + t579;
t259 = Icges(5,1) * t343 - t580;
t295 = -Icges(4,2) * t378 + t581;
t297 = Icges(4,1) * t380 - t582;
t703 = t257 * t343 + t259 * t342 + t295 * t380 + t297 * t378;
t668 = t173 * t563 + t175 * t565 + t206 * t559 + t208 * t561 + t706 * t381;
t667 = -t174 * t563 - t176 * t565 - t207 * t559 - t209 * t561 - t651 * t381;
t666 = -t647 * t381 + t652;
t665 = -t651 * t379 + t707;
t661 = t174 * t342 - t176 * t343 + t207 * t378 - t209 * t380;
t682 = t705 * t379;
t704 = t257 * t342 - t259 * t343 + t295 * t378 - t297 * t380;
t662 = t173 * t342 - t175 * t343 + t206 * t378 - t208 * t380;
t677 = t705 * t381;
t660 = t703 * t379 + t677;
t659 = t703 * t381 - t682;
t701 = t651 * qJD(1);
t700 = t706 * qJD(1);
t699 = t381 ^ 2;
t698 = t703 * qJD(1) - qJD(3) * t708;
t697 = t647 * qJD(1) + t682 * qJD(3) + t701;
t696 = -t695 * qJD(1) - t677 * qJD(3) + t700;
t694 = t665 * t379 + t666 * t381;
t693 = t667 * t379 + t668 * t381;
t510 = qJD(3) * t379;
t107 = qJD(1) * t174 + t257 * t510;
t214 = t259 * t379;
t109 = qJD(1) * t176 + qJD(3) * t214;
t123 = qJD(1) * t207 + t295 * t510;
t247 = t297 * t379;
t125 = qJD(1) * t209 + qJD(3) * t247;
t692 = t662 * qJD(3) - t107 * t343 - t109 * t342 - t123 * t380 - t125 * t378 + t700;
t229 = t437 * qJD(3);
t230 = t440 * qJD(3);
t267 = t438 * qJD(3);
t268 = t441 * qJD(3);
t691 = t705 * qJD(1) + t704 * qJD(3) + t229 * t343 + t230 * t342 + t267 * t380 + t268 * t378;
t213 = t257 * t381;
t106 = qJD(1) * t173 - qJD(3) * t213;
t215 = t259 * t381;
t108 = -qJD(3) * t215 + (t379 * t440 + t575) * qJD(1);
t246 = t295 * t381;
t122 = qJD(1) * t206 - qJD(3) * t246;
t248 = t297 * t381;
t124 = -qJD(3) * t248 + (t379 * t441 + t576) * qJD(1);
t690 = t661 * qJD(3) + t106 * t343 + t108 * t342 + t122 * t380 + t124 * t378 + t701;
t533 = t295 + t441;
t534 = -t438 + t297;
t535 = t257 + t440;
t536 = -t437 + t259;
t689 = (t342 * t535 - t343 * t536 + t378 * t533 - t380 * t534) * qJD(1);
t349 = qJ(5) + t374;
t332 = sin(t349);
t333 = cos(t349);
t566 = t333 * t379;
t283 = Icges(6,4) * t566;
t567 = t332 * t379;
t574 = Icges(6,5) * t381;
t163 = Icges(6,1) * t567 + t283 + t574;
t373 = qJD(3) + qJD(5);
t290 = t373 * t379;
t291 = t373 * t381;
t578 = Icges(6,4) * t332;
t238 = Icges(6,1) * t333 - t578;
t436 = Icges(6,2) * t333 + t578;
t680 = -t436 + t238;
t577 = Icges(6,4) * t333;
t439 = Icges(6,1) * t332 + t577;
t164 = -Icges(6,5) * t379 + t381 * t439;
t236 = -Icges(6,2) * t332 + t577;
t683 = t236 * t381 + t164;
t628 = qJD(1) * t680 - t290 * t683 + t291 * (-Icges(6,2) * t567 + t163 + t283);
t681 = t236 + t439;
t162 = -Icges(6,6) * t379 + t381 * t436;
t684 = -t238 * t381 + t162;
t161 = Icges(6,6) * t381 + t379 * t436;
t685 = -t238 * t379 + t161;
t629 = qJD(1) * t681 - t290 * t684 + t291 * t685;
t688 = t629 * t332 - t628 * t333;
t687 = rSges(4,2) * t378;
t596 = rSges(5,2) * t342;
t686 = t660 * qJD(1);
t452 = rSges(5,1) * t342 + rSges(5,2) * t343;
t611 = pkin(3) * t378;
t410 = t452 + t611;
t679 = t659 * qJD(1);
t678 = t708 * qJD(1);
t451 = rSges(6,1) * t332 + rSges(6,2) * t333;
t395 = t379 * (t209 + t246) - t381 * (-Icges(4,2) * t561 + t208 + t326);
t396 = t379 * (t207 - t248) - t381 * (t206 - t247);
t397 = t379 * (t176 + t213) - t381 * (-Icges(5,2) * t565 + t175 + t306);
t398 = t379 * (t174 - t215) - t381 * (t173 - t214);
t675 = -t398 * t342 + t397 * t343 - t396 * t378 + t395 * t380;
t674 = qJD(3) * t693 + t686;
t673 = qJD(3) * t694 - t679;
t672 = t379 * t698 + t381 * t691;
t671 = -t379 * t691 + t381 * t698;
t670 = -qJD(3) * t647 - t107 * t342 + t109 * t343 - t123 * t378 + t125 * t380;
t669 = t695 * qJD(3) - t106 * t342 + t108 * t343 - t122 * t378 + t124 * t380;
t664 = rSges(4,2) * t380;
t369 = t381 * rSges(5,3);
t181 = rSges(5,1) * t565 + rSges(5,2) * t563 + t369;
t509 = qJD(3) * t381;
t480 = t380 * t509;
t319 = pkin(3) * t480;
t352 = qJD(2) * t381;
t508 = qJD(4) * t379;
t461 = t352 - t508;
t454 = -t319 - t461;
t330 = pkin(3) * t561;
t377 = -qJ(4) - pkin(6);
t602 = pkin(6) + t377;
t225 = -t381 * t602 + t330;
t303 = t381 * pkin(1) + t379 * qJ(2);
t607 = pkin(6) * t381;
t470 = t303 + t607;
t456 = t225 + t470;
t600 = rSges(5,1) * t343;
t263 = -t596 + t600;
t487 = t263 * t509;
t66 = -t487 + (t181 + t456) * qJD(1) + t454;
t663 = t381 * t66;
t433 = Icges(6,5) * t332 + Icges(6,6) * t333;
t160 = -Icges(6,3) * t379 + t381 * t433;
t61 = -t381 * t160 - t162 * t566 - t164 * t567;
t421 = t236 * t333 + t238 * t332;
t234 = Icges(6,5) * t333 - Icges(6,6) * t332;
t571 = t234 * t381;
t76 = t379 * t421 + t571;
t658 = t76 * qJD(1) + t290 * t61;
t429 = t162 * t333 + t164 * t332;
t657 = t381 * t429;
t502 = t381 * t611;
t224 = t379 * t602 + t502;
t354 = t381 * qJ(2);
t299 = pkin(1) * t379 - t354;
t277 = qJD(1) * t299;
t650 = qJD(1) * t224 - t277;
t370 = t381 * rSges(4,3);
t218 = rSges(4,1) * t561 + rSges(4,2) * t559 + t370;
t649 = t218 + t470;
t278 = pkin(4) * t342 + t611;
t648 = -t278 - t451;
t304 = -rSges(3,2) * t381 + t379 * rSges(3,3);
t646 = t697 * t699 + (t690 * t379 + (-t692 + t696) * t381) * t379;
t645 = t692 * t699 + (t696 * t379 + (-t690 + t697) * t381) * t379;
t302 = rSges(4,1) * t380 - t687;
t250 = t302 * t381;
t453 = rSges(4,1) * t378 + t664;
t126 = -qJD(3) * t250 + (t379 * t453 + t370) * qJD(1);
t231 = t303 * qJD(1) - t352;
t269 = t453 * qJD(3);
t505 = qJD(1) * qJD(3);
t274 = qJDD(3) * t379 + t381 * t505;
t506 = qJD(1) * qJD(2);
t529 = qJDD(2) * t379 + t381 * t506;
t606 = pkin(6) * qJD(1) ^ 2;
t416 = -t381 * t606 + t529;
t219 = -t379 * rSges(4,3) + t381 * t453;
t608 = pkin(6) * t379;
t471 = -t299 - t608;
t458 = t219 + t471;
t50 = -t269 * t510 + t274 * t302 + (-t126 - t231) * qJD(1) + t458 * qJDD(1) + t416;
t481 = t380 * t510;
t514 = qJD(1) * t381;
t489 = t378 * t514;
t492 = t514 * t664 + (t481 + t489) * rSges(4,1);
t511 = qJD(3) * t378;
t127 = (-rSges(4,2) * t511 - rSges(4,3) * qJD(1)) * t379 + t492;
t347 = qJDD(3) * t381;
t275 = -t379 * t505 + t347;
t515 = qJD(1) * t379;
t338 = qJ(2) * t514;
t351 = qJD(2) * t379;
t528 = t338 + t351;
t497 = qJD(1) * (-pkin(1) * t515 + t528) + qJDD(1) * t303 + t379 * t506;
t409 = qJDD(1) * t607 - t379 * t606 + t497;
t51 = qJD(1) * t127 + qJDD(1) * t218 - t275 * t302 + (qJD(3) * t269 - qJDD(2)) * t381 + t409;
t644 = t379 * t50 - t381 * t51;
t375 = t379 ^ 2;
t610 = pkin(3) * t380;
t643 = (-t375 - t699) * t610;
t636 = g(1) * t379 - g(2) * t381;
t110 = -t487 + (t379 * t452 + t369) * qJD(1);
t485 = t343 * t510;
t494 = rSges(5,1) * t485 + t452 * t514;
t512 = qJD(3) * t342;
t499 = rSges(5,2) * t512;
t111 = (-rSges(5,3) * qJD(1) - t499) * t379 + t494;
t350 = qJD(4) * t381;
t317 = pkin(3) * t481;
t491 = -pkin(3) * t489 - t377 * t515 - t317;
t459 = t350 - t491;
t137 = pkin(6) * t515 + t459;
t635 = t110 * t381 + (-t111 - t137) * t379;
t464 = t681 * t373;
t465 = t680 * t373;
t632 = qJD(1) * t234 + t332 * t464 - t333 * t465;
t466 = -qJD(1) * t161 + t373 * t683;
t468 = (t379 * t439 + t574) * qJD(1) + t684 * t373;
t523 = qJD(1) * t160;
t631 = t332 * t468 - t333 * t466 + t523;
t467 = qJD(1) * t162 + t163 * t373 + t236 * t290;
t469 = -qJD(1) * t164 + t373 * t685;
t159 = Icges(6,3) * t381 + t379 * t433;
t524 = qJD(1) * t159;
t630 = t332 * t469 - t333 * t467 + t524;
t627 = -m(5) - m(6);
t626 = -pkin(1) - pkin(6);
t187 = qJD(5) * t514 + qJDD(5) * t379 + t274;
t625 = t187 / 0.2e1;
t188 = qJDD(5) * t381 - t373 * t515 + t347;
t624 = t188 / 0.2e1;
t622 = t274 / 0.2e1;
t621 = t275 / 0.2e1;
t620 = -t290 / 0.2e1;
t619 = t290 / 0.2e1;
t618 = -t291 / 0.2e1;
t617 = t291 / 0.2e1;
t616 = t379 / 0.2e1;
t615 = t381 / 0.2e1;
t614 = rSges(3,2) - pkin(1);
t613 = -rSges(5,3) - pkin(1);
t612 = -rSges(6,3) - pkin(1);
t609 = pkin(4) * t343;
t183 = t451 * t373;
t594 = rSges(6,2) * t332;
t598 = rSges(6,1) * t333;
t240 = -t594 + t598;
t382 = qJD(3) ^ 2;
t408 = qJDD(4) * t381 + t274 * t610 + t416;
t457 = t224 + t471;
t372 = -pkin(7) + t377;
t537 = t381 * t278 + t379 * t372;
t135 = t377 * t379 + t502 - t537;
t166 = -t379 * rSges(6,3) + t381 * t451;
t555 = -t135 + t166;
t413 = t457 + t555;
t325 = t377 * t514;
t530 = t319 + t325;
t138 = t508 + (t330 - t607) * qJD(1) - t530;
t455 = -t138 - t231 - t508;
t564 = t342 * t382;
t279 = t609 + t610;
t260 = t279 * qJD(3);
t562 = t372 * t381;
t83 = -t260 * t381 + (-t562 + (t278 - t611) * t379) * qJD(1) + t530;
t368 = t381 * rSges(6,3);
t570 = t240 * t291;
t98 = -t570 + (t379 * t451 + t368) * qJD(1);
t9 = -t382 * t330 - t183 * t290 + t187 * t240 + (t274 * t343 - t379 * t564) * pkin(4) + t413 * qJDD(1) + (t455 - t83 - t98) * qJD(1) + t408;
t605 = t9 * t379;
t604 = -qJD(1) / 0.2e1;
t603 = qJD(1) / 0.2e1;
t592 = rSges(3,3) * t381;
t501 = t382 * t611;
t386 = qJDD(1) * t225 + qJDD(4) * t379 + t381 * t501 + t409 + (t137 + t350) * qJD(1);
t504 = qJDD(2) * t381;
t251 = t379 * t278;
t136 = t251 - t330 + (-t372 + t377) * t381;
t165 = rSges(6,1) * t567 + rSges(6,2) * t566 + t368;
t553 = t136 + t165;
t496 = t379 * t260 + t278 * t514 + t372 * t515;
t82 = t491 + t496;
t288 = rSges(6,1) * t566;
t495 = t373 * t288 + t451 * t514;
t99 = (-rSges(6,3) * qJD(1) - t373 * t594) * t379 + t495;
t10 = t386 + (t82 + t99) * qJD(1) + t553 * qJDD(1) + (-t275 * t343 + t381 * t564) * pkin(4) - t188 * t240 - t275 * t610 + t291 * t183 - t504;
t591 = t10 * t381;
t473 = -t263 - t610;
t232 = t452 * qJD(3);
t513 = qJD(3) * t232;
t25 = t386 + qJD(1) * t111 + qJDD(1) * t181 + t473 * t275 + (-qJDD(2) + t513) * t381;
t590 = t25 * t381;
t526 = t350 + t351;
t490 = t317 + t526;
t412 = pkin(4) * t485 + t290 * t240 + t490;
t48 = qJD(1) * t413 + t412;
t587 = t381 * t48;
t365 = t379 * rSges(5,3);
t182 = t381 * t452 - t365;
t415 = t182 + t457;
t460 = t263 * t510 + t490;
t65 = qJD(1) * t415 + t460;
t585 = t381 * t65;
t584 = qJDD(1) / 0.2e1;
t583 = -t137 - t82;
t261 = t302 * t510;
t100 = qJD(1) * t458 + t261 + t351;
t573 = t100 * t381;
t560 = t379 * t160;
t192 = t379 * t234;
t77 = t381 * t421 - t192;
t558 = t77 * qJD(1);
t556 = t138 * t509 - t275 * t224;
t554 = -t136 - t225;
t544 = -t181 - t225;
t300 = rSges(3,2) * t379 + t592;
t532 = -t299 + t300;
t223 = t303 + t304;
t527 = rSges(3,2) * t515 + rSges(3,3) * t514;
t525 = t351 - t277;
t503 = -rSges(4,3) + t626;
t331 = pkin(3) * t559;
t500 = rSges(6,2) * t567;
t60 = t381 * t159 + t161 * t566 + t163 * t567;
t498 = -t225 - t553;
t484 = t343 * t509;
t483 = t378 * t510;
t479 = -t515 / 0.2e1;
t478 = t514 / 0.2e1;
t477 = -t510 / 0.2e1;
t476 = t510 / 0.2e1;
t475 = -t509 / 0.2e1;
t474 = t509 / 0.2e1;
t472 = t240 + t609;
t198 = t288 - t500;
t463 = qJD(1) * t198 + t291 * t451;
t289 = t381 * t594;
t199 = -t381 * t598 + t289;
t462 = -qJD(1) * t199 - t290 * t451;
t216 = rSges(5,1) * t563 - rSges(5,2) * t565;
t305 = rSges(2,1) * t381 - rSges(2,2) * t379;
t301 = rSges(2,1) * t379 + rSges(2,2) * t381;
t189 = t224 * t509;
t40 = -t165 * t290 - t166 * t291 - t189 + (t135 * t381 + t379 * t554) * qJD(3);
t442 = t40 * (-t198 * t290 + t291 * t199);
t101 = t649 * qJD(1) - t302 * t509 - t352;
t432 = t100 * t379 - t101 * t381;
t431 = t126 * t381 - t127 * t379;
t430 = t161 * t333 + t163 * t332;
t81 = t162 * t332 - t164 * t333;
t422 = -t218 * t379 - t219 * t381;
t414 = t60 + t560;
t411 = -t182 * t381 + t379 * t544;
t403 = -qJD(1) * t433 + t192 * t291 - t290 * t571;
t149 = t379 * t159;
t62 = -t381 * t430 + t149;
t400 = -qJD(1) * t429 - t373 * t571 + t524;
t399 = qJD(1) * t430 + t192 * t373 + t523;
t390 = qJD(1) * t421 - t433 * t373;
t11 = t399 * t379 + t630 * t381;
t12 = t400 * t379 - t631 * t381;
t13 = -t630 * t379 + t399 * t381;
t14 = t631 * t379 + t400 * t381;
t22 = t291 * t60 + t658;
t63 = -t560 + t657;
t23 = t290 * t63 + t291 * t62 - t558;
t34 = t390 * t379 + t632 * t381;
t35 = -t632 * t379 + t390 * t381;
t38 = -t332 * t467 - t333 * t469;
t39 = t332 * t466 + t333 * t468;
t80 = -t161 * t332 + t163 * t333;
t387 = (qJD(1) * t34 - qJDD(1) * t77 + t11 * t291 + t12 * t290 + t187 * t63 + t188 * t62) * t616 + (-t332 * t628 - t333 * t629) * t604 + (qJD(1) * t35 + qJDD(1) * t76 + t13 * t291 + t14 * t290 + t187 * t61 + t188 * t60) * t615 + t22 * t479 + t23 * t478 + (t11 * t381 + t12 * t379 + (-t379 * t62 + t381 * t63) * qJD(1)) * t619 + (t379 * t63 + t381 * t62) * t625 + (t379 * t61 + t381 * t60) * t624 + (t13 * t381 + t14 * t379 + (-t379 * t60 + t381 * t61) * qJD(1)) * t617 + (t379 * t81 + t381 * t80) * t584 + (t379 * t39 + t38 * t381 + (-t379 * t80 + t381 * t81) * qJD(1)) * t603 + (t403 * t379 + t381 * t688) * t620 + (-t379 * t688 + t403 * t381) * t618;
t310 = t381 * t596;
t252 = t379 * t279;
t249 = t302 * t379;
t217 = -t381 * t600 + t310;
t200 = t381 * t224;
t191 = (-t279 + t610) * t381;
t190 = -t331 + t252;
t154 = qJD(1) * t223 - t352;
t153 = qJD(1) * t532 + t351;
t152 = t381 * t166;
t129 = t381 * t138;
t112 = t422 * qJD(3);
t89 = t381 * t98;
t79 = qJD(1) * t527 + qJDD(1) * t304 + t497 - t504;
t78 = t532 * qJDD(1) + (-t304 * qJD(1) - t231) * qJD(1) + t529;
t64 = qJD(3) * t411 - t189;
t49 = -pkin(4) * t484 - t570 + (t456 + t553) * qJD(1) + t454;
t24 = t263 * t274 + (-t501 - t513) * t379 + t415 * qJDD(1) + (-t110 + t455) * qJD(1) + t408;
t5 = t135 * t275 - t165 * t187 - t166 * t188 - t290 * t99 + t291 * t98 + t554 * t274 + (t379 * t583 + t381 * t83) * qJD(3) + t556;
t1 = [-t187 * t77 / 0.2e1 + ((t414 + t63 - t657) * t291 + t658) * t620 - m(2) * (-g(1) * t301 + g(2) * t305) + t81 * t625 + (t80 + t76) * t624 - t659 * t274 / 0.2e1 + t661 * t622 + (t558 + (t379 * t429 - t149 + t61) * t291 + (t414 - t60) * t290 + ((t160 + t430) * t291 - t429 * t290) * t381 + t23) * t618 + (t38 + t35) * t617 + (((t652 - t666 + t667) * t379 + (-t707 + (-t647 + t651) * t379 + t665 + t668) * t381) * qJD(3) + t686) * t477 + (t39 + t34 + t22) * t619 + (-t703 * qJD(3) + t229 * t342 - t230 * t343 + t267 * t378 - t268 * t380 - t332 * t465 - t333 * t464) * qJD(1) + (t48 * t461 + t49 * (-t373 * t500 + t338 + t495 + t496 + t526) + (t240 * t373 + t260) * t587 + ((t372 + t612) * t587 + (t48 * (-qJ(2) + t648) + t49 * t612) * t379) * qJD(1) - (-t48 + (t555 - t608) * qJD(1) + t412 + t650) * t49 + (t10 - g(2)) * (t165 + t251 + t303 - t562) + (t9 - g(1)) * (t166 - t299 + t537)) * m(6) + (t65 * (rSges(5,1) * t484 - t509 * t596 + t325 - t454) + t66 * (-t379 * t499 + t459 + t494 + t528) + (t613 * t585 + (t65 * (-qJ(2) - t410) + t66 * t613) * t379) * qJD(1) - (-t65 + (t182 - t608) * qJD(1) + t460 + t650) * t66 + (t25 - g(2)) * (-t377 * t381 + t181 + t303 + t330) + (-g(1) + t24) * (t354 - t365 + (-pkin(1) + t377) * t379 + t410 * t381)) * m(5) + (t100 * (rSges(4,1) * t480 - t509 * t687 + t352) + t101 * (-rSges(4,2) * t483 + t492 + t528) + (t503 * t573 + (t100 * (-qJ(2) - t453) + t101 * t503) * t379) * qJD(1) - (-t100 + t261 + (t219 - t608) * qJD(1) + t525) * t101 + (-g(2) + t51) * t649 + (-g(1) + t50) * (t379 * t626 + t219 + t354)) * m(4) + (t153 * t352 + t154 * (t527 + t528) + (t153 * t614 * t381 + (t153 * (-rSges(3,3) - qJ(2)) - t154 * pkin(1)) * t379) * qJD(1) - (qJD(1) * t300 - t153 + t525) * t154 + (t79 - g(2)) * t223 + (t78 - g(1)) * (t379 * t614 + t354 + t592)) * m(3) + (t660 - t662) * t621 + ((t651 * t375 + ((t647 + t651) * t381 - t652 + t667) * t381) * qJD(3) + t673 + t679) * t475 + (t670 + t671) * t474 + (t669 + t672 + t674) * t476 + (m(2) * (t301 ^ 2 + t305 ^ 2) + Icges(2,3) + Icges(3,1) - t236 * t332 + t238 * t333 - t704) * qJDD(1); (-m(3) + t627) * t636 + 0.2e1 * (-t591 / 0.2e1 + t605 / 0.2e1) * m(6) - t590 * m(5) - t381 * t79 * m(3) + (-t636 + t644) * m(4) + 0.2e1 * (m(3) * t78 + m(5) * t24) * t616; t387 + t694 * t622 + t693 * t621 + (t672 * qJD(1) + t645 * qJD(3) - t659 * qJDD(1) + t665 * t274 + t666 * t275) * t616 + (t671 * qJD(1) + t646 * qJD(3) + t660 * qJDD(1) + t667 * t274 + t668 * t275) * t615 + ((t342 * t397 + t343 * t398 + t378 * t395 + t380 * t396) * qJD(3) + (-t342 * t536 - t343 * t535 - t378 * t534 - t380 * t533) * qJD(1)) * t604 + (t670 * t381 + t669 * t379 + (t662 * t379 + t661 * t381) * qJD(1)) * t603 + (t661 * t379 - t662 * t381) * t584 + t674 * t479 + t673 * t478 + ((-t510 * t677 - t678) * t379 + ((t379 * t682 + t675) * qJD(3) + t689) * t381) * t477 + ((-t666 * t379 + t665 * t381) * qJD(1) + t645) * t476 + ((t509 * t682 - t678) * t381 + ((-t381 * t677 - t675) * qJD(3) - t689) * t379) * t475 + ((-t668 * t379 + t667 * t381) * qJD(1) + t646) * t474 + (-t48 * (-qJD(1) * t191 + t462) - t49 * (qJD(1) * t190 + t463) - t442 - (t40 * (t191 * t381 + t643) + (-t40 * t190 - t278 * t48) * t379) * qJD(3) + t9 * t331 + t5 * (-t152 - t200) + t40 * (t129 + t89) + (t10 * (-t240 - t279) + t49 * t183 + t5 * t135 + t40 * t83 + (t40 * t498 + t472 * t48) * qJD(1)) * t381 + (t9 * t472 + t48 * (-pkin(3) * t511 - pkin(4) * t512 - t183) + t5 * t498 + t40 * (-t99 + t583) + (t49 * t472 + t40 * (t224 + t555)) * qJD(1)) * t379 - g(1) * (t198 + t252) - g(2) * (t289 + (-t279 - t598) * t381) - g(3) * t648) * m(6) + (-g(1) * (t216 + t331) - g(2) * (t310 + (-t600 - t610) * t381) + g(3) * t410 - (t452 * t663 + t64 * (t217 * t381 + t643) + (-t64 * t216 - t65 * t410) * t379) * qJD(3) + t24 * (t263 * t379 + t331) - t65 * pkin(3) * t483 + t473 * t590 + (t635 * qJD(3) - t182 * t275 + t544 * t274 + t556) * (-t200 + t411) + t64 * (t129 + t635) - (t65 * t379 - t663) * t232 + (t65 * t217 - t66 * t216 + t64 * (t544 * t381 + (t182 + t224) * t379) + (t379 * t66 + t585) * t263) * qJD(1)) * m(5) + (-(t100 * t250 + t101 * t249) * qJD(1) - (t112 * (-t249 * t379 - t250 * t381) - t432 * t453) * qJD(3) - g(1) * t249 + g(2) * t250 + g(3) * t453 + (qJD(3) * t431 - t218 * t274 - t219 * t275) * t422 + t112 * ((-t218 * t381 + t219 * t379) * qJD(1) + t431) - t432 * t269 + ((t101 * t379 + t573) * qJD(1) + t644) * t302) * m(4); t627 * (g(1) * t381 + g(2) * t379) + m(5) * (t24 * t381 + t25 * t379) + m(6) * (t10 * t379 + t381 * t9); t387 + (t5 * (-t165 * t379 - t152) + t40 * (-t379 * t99 + t89 + (-t165 * t381 + t166 * t379) * qJD(1)) - (t379 * t48 - t381 * t49) * t183 + (-t591 + t605 + (t379 * t49 + t587) * qJD(1)) * t240 - t462 * t48 - t463 * t49 - t442 - g(1) * t198 - g(2) * t199 + g(3) * t451) * m(6);];
tau = t1;
