% Calculate vector of inverse dynamics joint torques for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR6_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR6_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR6_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:23
% EndTime: 2019-12-31 17:04:53
% DurationCPUTime: 26.91s
% Computational Cost: add. (11298->720), mult. (13632->921), div. (0->0), fcn. (10741->8), ass. (0->396)
t673 = Icges(3,3) + Icges(4,3);
t348 = qJ(2) + pkin(7);
t319 = sin(t348);
t320 = cos(t348);
t352 = sin(qJ(2));
t354 = cos(qJ(2));
t661 = Icges(3,5) * t354 + Icges(4,5) * t320 - Icges(3,6) * t352 - Icges(4,6) * t319;
t355 = cos(qJ(1));
t672 = t673 * t355;
t353 = sin(qJ(1));
t526 = t353 * t354;
t530 = t352 * t353;
t533 = t320 * t353;
t535 = t319 * t353;
t662 = -Icges(3,5) * t526 - Icges(4,5) * t533 + Icges(3,6) * t530 + Icges(4,6) * t535 + t672;
t667 = t673 * t353 + t661 * t355;
t558 = Icges(4,6) * t355;
t174 = Icges(4,4) * t533 - Icges(4,2) * t535 - t558;
t559 = Icges(3,6) * t355;
t201 = Icges(3,4) * t526 - Icges(3,2) * t530 - t559;
t671 = t174 * t319 + t201 * t352;
t284 = Icges(4,4) * t535;
t564 = Icges(4,5) * t355;
t176 = Icges(4,1) * t533 - t284 - t564;
t300 = Icges(3,4) * t530;
t565 = Icges(3,5) * t355;
t203 = Icges(3,1) * t526 - t300 - t565;
t651 = -t176 * t320 - t203 * t354 + t671;
t635 = -t353 * t651 + t662 * t355;
t567 = Icges(4,4) * t319;
t246 = Icges(4,1) * t320 - t567;
t177 = Icges(4,5) * t353 + t246 * t355;
t568 = Icges(3,4) * t352;
t277 = Icges(3,1) * t354 - t568;
t204 = Icges(3,5) * t353 + t277 * t355;
t669 = -t177 * t533 - t204 * t526;
t668 = Icges(3,5) * t352 + Icges(4,5) * t319 + Icges(3,6) * t354 + Icges(4,6) * t320;
t243 = Icges(4,2) * t320 + t567;
t304 = Icges(4,4) * t320;
t245 = Icges(4,1) * t319 + t304;
t274 = Icges(3,2) * t354 + t568;
t335 = Icges(3,4) * t354;
t276 = Icges(3,1) * t352 + t335;
t659 = t243 * t319 - t245 * t320 + t274 * t352 - t276 * t354;
t666 = t355 * t667 + t669;
t525 = t354 * t355;
t532 = t320 * t355;
t607 = -t177 * t532 - t204 * t525 - t667 * t353;
t665 = -t176 * t532 - t203 * t525 + t662 * t353;
t410 = -Icges(4,2) * t319 + t304;
t175 = Icges(4,6) * t353 + t355 * t410;
t411 = -Icges(3,2) * t352 + t335;
t202 = Icges(3,6) * t353 + t355 * t411;
t663 = t175 * t319 + t202 * t352;
t634 = -t175 * t535 - t202 * t530 - t666;
t529 = t352 * t355;
t534 = t319 * t355;
t633 = -t174 * t534 - t201 * t529 - t665;
t632 = -t175 * t534 - t202 * t529 - t607;
t630 = t174 * t320 + t176 * t319 + t201 * t354 + t203 * t352;
t629 = t175 * t320 + t177 * t319 + t202 * t354 + t204 * t352;
t660 = t668 * qJD(2);
t658 = -t243 * t320 - t245 * t319 - t274 * t354 - t276 * t352;
t657 = t177 * t320 + t204 * t354 - t663;
t616 = t668 * t355;
t615 = t668 * t353;
t628 = -t353 * t659 - t616;
t627 = -t355 * t659 + t615;
t656 = t667 * qJD(1);
t655 = t355 ^ 2;
t217 = t410 * qJD(2);
t218 = t246 * qJD(2);
t252 = t411 * qJD(2);
t253 = t277 * qJD(2);
t654 = t668 * qJD(1) + t658 * qJD(2) - t217 * t319 + t218 * t320 - t252 * t352 + t253 * t354;
t387 = qJD(2) * t274;
t122 = -t355 * t387 + (-t353 * t411 + t559) * qJD(1);
t389 = qJD(2) * t276;
t124 = -t355 * t389 + (-t277 * t353 + t565) * qJD(1);
t386 = qJD(2) * t243;
t97 = -t355 * t386 + (-t353 * t410 + t558) * qJD(1);
t388 = qJD(2) * t245;
t99 = -t355 * t388 + (-t246 * t353 + t564) * qJD(1);
t653 = -qJD(2) * t629 - t122 * t352 + t124 * t354 - t319 * t97 + t320 * t99 + t656;
t100 = qJD(1) * t177 - t353 * t388;
t123 = qJD(1) * t202 - t353 * t387;
t125 = qJD(1) * t204 - t353 * t389;
t98 = qJD(1) * t175 - t353 * t386;
t652 = qJD(1) * t662 + qJD(2) * t630 - t100 * t320 + t123 * t352 - t125 * t354 + t319 * t98;
t650 = t353 * t632 - t355 * t633;
t649 = t634 * t353 - t355 * t635;
t648 = t659 * qJD(1) + qJD(2) * t661;
t647 = qJD(1) * t651 - t353 * t660 + t656;
t646 = -t660 * t355 + (-t353 * t661 - t657 + t672) * qJD(1);
t577 = rSges(4,2) * t320;
t248 = rSges(4,1) * t319 + t577;
t211 = t248 * t353;
t337 = t353 * rSges(4,3);
t305 = t320 * rSges(4,1);
t609 = -rSges(4,2) * t319 + t305;
t102 = -qJD(2) * t211 + (t355 * t609 + t337) * qJD(1);
t219 = t609 * qJD(2);
t471 = qJD(1) * qJD(2);
t311 = t353 * t471;
t260 = -qJDD(2) * t355 + t311;
t470 = qJD(1) * qJD(3);
t584 = pkin(2) * t352;
t458 = qJDD(3) * t353 + t260 * t584 + t355 * t470;
t285 = rSges(4,2) * t535;
t182 = rSges(4,1) * t533 - t355 * rSges(4,3) - t285;
t344 = t355 * pkin(5);
t287 = pkin(1) * t353 - t344;
t351 = -qJ(3) - pkin(5);
t313 = t355 * t351;
t343 = t354 * pkin(2);
t314 = t343 + pkin(1);
t487 = t353 * t314 + t313;
t170 = t287 - t487;
t503 = t170 - t287;
t459 = -t182 + t503;
t356 = qJD(2) ^ 2;
t524 = t354 * t356;
t466 = pkin(2) * t524;
t342 = t353 * pkin(5);
t476 = qJD(1) * t353;
t297 = t351 * t476;
t468 = pkin(2) * t530;
t293 = qJD(2) * t468;
t325 = qJD(3) * t355;
t485 = t293 + t325;
t457 = t297 + t485;
t580 = pkin(1) - t314;
t119 = (-t355 * t580 - t342) * qJD(1) - t457;
t288 = t355 * pkin(1) + t342;
t257 = t288 * qJD(1);
t518 = -t119 - t257;
t30 = t248 * t260 + (-qJD(2) * t219 - t466) * t355 + t459 * qJDD(1) + (-t102 + t518) * qJD(1) + t458;
t645 = t30 - g(1);
t644 = t627 * qJD(1);
t643 = t628 * qJD(1);
t473 = qJD(2) * t355;
t475 = qJD(1) * t355;
t490 = rSges(4,3) * t475 + qJD(1) * t285;
t101 = -t473 * t577 + (-t319 * t473 - t320 * t476) * rSges(4,1) + t490;
t183 = rSges(4,1) * t532 - rSges(4,2) * t534 + t337;
t259 = qJDD(2) * t353 + t355 * t471;
t318 = pkin(5) * t475;
t324 = qJD(3) * t353;
t454 = t352 * t473;
t429 = pkin(2) * t454;
t118 = -t429 - t318 + t324 + (t353 * t580 - t313) * qJD(1);
t291 = t355 * t314;
t431 = -t351 * t353 + t291;
t171 = t431 - t288;
t492 = qJD(1) * (-pkin(1) * t476 + t318) + qJDD(1) * t288;
t361 = qJD(1) * t118 + qJDD(1) * t171 + t353 * t470 + (-t259 * t352 - t353 * t524) * pkin(2) - qJDD(3) * t355 + t492;
t474 = qJD(2) * t353;
t31 = qJD(1) * t101 + qJDD(1) * t183 - t219 * t474 - t248 * t259 + t361;
t642 = t31 - g(2);
t641 = qJD(2) * t649 + t643;
t640 = qJD(2) * t650 + t644;
t639 = t353 * t648 + t355 * t654;
t638 = t353 * t654 - t355 * t648;
t637 = t651 * qJD(2) - t100 * t319 - t123 * t354 - t125 * t352 - t320 * t98;
t636 = qJD(2) * t657 + t122 * t354 + t124 * t352 + t319 * t99 + t320 * t97;
t631 = rSges(3,2) * t352;
t323 = qJ(4) + t348;
t308 = cos(t323);
t537 = t308 * t353;
t307 = sin(t323);
t539 = t307 * t353;
t557 = Icges(5,6) * t355;
t159 = Icges(5,4) * t537 - Icges(5,2) * t539 - t557;
t295 = Icges(5,4) * t308;
t224 = Icges(5,1) * t307 + t295;
t626 = -t224 * t353 - t159;
t409 = -Icges(5,2) * t307 + t295;
t160 = Icges(5,6) * t353 + t355 * t409;
t625 = -t224 * t355 - t160;
t566 = Icges(5,4) * t307;
t225 = Icges(5,1) * t308 - t566;
t162 = Icges(5,5) * t353 + t225 * t355;
t222 = Icges(5,2) * t308 + t566;
t624 = -t222 * t355 + t162;
t623 = -t222 + t225;
t622 = t224 + t409;
t378 = t201 * t355 - t202 * t353;
t379 = t174 * t355 - t175 * t353;
t596 = t353 * (-t243 * t355 + t177) - t355 * (-Icges(4,2) * t533 + t176 - t284);
t597 = t353 * (-t274 * t355 + t204) - t355 * (-Icges(3,2) * t526 + t203 - t300);
t621 = -t319 * t596 + t379 * t320 - t352 * t597 + t378 * t354;
t488 = t276 + t411;
t489 = -t274 + t277;
t493 = t245 + t410;
t494 = -t243 + t246;
t620 = (-t319 * t493 + t320 * t494 - t352 * t488 + t354 * t489) * qJD(1);
t619 = t647 * t655 + (t653 * t353 + (-t646 + t652) * t355) * t353;
t618 = t652 * t655 + (t646 * t353 + (-t647 + t653) * t355) * t353;
t617 = t661 * qJD(1);
t383 = -t248 - t584;
t614 = t355 * t383;
t263 = qJD(1) * t287;
t613 = qJD(1) * t170 - t263;
t583 = pkin(3) * t319;
t261 = -t583 - t584;
t247 = t261 * qJD(2);
t215 = t355 * t247;
t612 = t215 + t324;
t611 = t609 + t343;
t296 = t308 * rSges(5,1);
t576 = rSges(5,2) * t307;
t610 = t296 - t576;
t306 = pkin(3) * t320;
t483 = t306 + t343;
t608 = t662 + t663;
t221 = Icges(5,5) * t308 - Icges(5,6) * t307;
t347 = qJD(2) + qJD(4);
t220 = Icges(5,5) * t307 + Icges(5,6) * t308;
t544 = t220 * t355;
t550 = t160 * t307;
t554 = Icges(5,3) * t355;
t604 = -t347 * t544 + (-t162 * t308 - t221 * t353 + t550 + t554) * qJD(1);
t267 = Icges(5,4) * t539;
t563 = Icges(5,5) * t355;
t161 = Icges(5,1) * t537 - t267 - t563;
t405 = t159 * t307 - t161 * t308;
t158 = Icges(5,3) * t353 + t221 * t355;
t481 = qJD(1) * t158;
t545 = t220 * t353;
t603 = qJD(1) * t405 - t347 * t545 + t481;
t397 = t222 * t307 - t224 * t308;
t600 = qJD(1) * t397 + t221 * t347;
t270 = t347 * t353;
t271 = t347 * t355;
t595 = qJD(1) * t622 + t270 * t624 - t271 * (-Icges(5,2) * t537 + t161 - t267);
t469 = qJD(1) * qJD(4);
t187 = qJDD(4) * t353 + t355 * t469 + t259;
t594 = t187 / 0.2e1;
t188 = t353 * t469 + t311 + (-qJDD(2) - qJDD(4)) * t355;
t593 = t188 / 0.2e1;
t592 = t259 / 0.2e1;
t591 = t260 / 0.2e1;
t590 = -t270 / 0.2e1;
t589 = t270 / 0.2e1;
t588 = -t271 / 0.2e1;
t587 = t271 / 0.2e1;
t586 = t353 / 0.2e1;
t585 = -t355 / 0.2e1;
t582 = -qJD(1) / 0.2e1;
t581 = qJD(1) / 0.2e1;
t579 = rSges(3,1) * t354;
t575 = rSges(5,2) * t308;
t338 = t353 * rSges(3,3);
t336 = t353 * rSges(5,3);
t226 = rSges(5,1) * t307 + t575;
t372 = -t271 * t226 + t612;
t346 = -pkin(6) + t351;
t464 = pkin(1) + t483;
t495 = -t355 * t346 - t353 * t464;
t130 = t487 + t495;
t268 = rSges(5,2) * t539;
t165 = rSges(5,1) * t537 - t355 * rSges(5,3) - t268;
t517 = t130 - t165;
t428 = t503 + t517;
t47 = qJD(1) * t428 + t372;
t574 = t353 * t47;
t382 = qJD(2) * t614 + t324;
t63 = qJD(1) * t459 + t382;
t573 = t63 * t248;
t572 = qJDD(1) / 0.2e1;
t278 = rSges(3,1) * t352 + rSges(3,2) * t354;
t455 = t278 * t473;
t484 = rSges(3,2) * t530 + t355 * rSges(3,3);
t213 = rSges(3,1) * t526 - t484;
t497 = -t213 - t287;
t109 = qJD(1) * t497 - t455;
t553 = t109 * t353;
t552 = t109 * t355;
t214 = rSges(3,1) * t525 - rSges(3,2) * t529 + t338;
t153 = t214 + t288;
t110 = qJD(1) * t153 - t278 * t474;
t238 = t278 * t355;
t551 = t110 * t238;
t538 = t307 * t355;
t536 = t308 * t355;
t531 = t320 * t356;
t528 = t353 * t247;
t527 = t353 * t346;
t157 = Icges(5,5) * t537 - Icges(5,6) * t539 - t554;
t523 = t355 * t157;
t74 = -t353 * t397 - t544;
t520 = t74 * qJD(1);
t230 = t355 * t464;
t482 = t351 - t346;
t131 = t353 * t482 + t230 - t291;
t166 = rSges(5,1) * t536 - rSges(5,2) * t538 + t336;
t516 = t131 + t166;
t515 = -t353 * t157 - t161 * t536;
t514 = t353 * t158 + t162 * t536;
t511 = t353 * t165 + t355 * t166;
t509 = -t170 * t474 + t171 * t473;
t508 = -t353 * t170 + t355 * t171;
t502 = -t171 - t183;
t491 = rSges(5,3) * t475 + qJD(1) * t268;
t486 = rSges(3,3) * t475 + t476 * t631;
t467 = pkin(2) * t529;
t93 = -t271 * t575 + (-t271 * t307 - t308 * t476) * rSges(5,1) + t491;
t197 = t226 * t353;
t94 = -t347 * t197 + (t355 * t610 + t336) * qJD(1);
t465 = t165 * t475 + t353 * t94 + t355 * t93;
t463 = qJD(2) * t343;
t462 = t118 * t473 + t119 * t474 - t259 * t170;
t461 = t355 * t118 + t353 * t119 - t170 * t475;
t460 = -t171 - t516;
t453 = t354 * t473;
t451 = -pkin(1) - t579;
t450 = t476 / 0.2e1;
t449 = t475 / 0.2e1;
t448 = -t474 / 0.2e1;
t447 = t474 / 0.2e1;
t446 = -t473 / 0.2e1;
t445 = t473 / 0.2e1;
t443 = qJD(1) * t162 + t347 * t626;
t442 = (-t225 * t353 + t563) * qJD(1) + t625 * t347;
t441 = qJD(1) * t160 + t161 * t347 - t222 * t270;
t440 = (-t353 * t409 + t557) * qJD(1) + t624 * t347;
t133 = t162 * t537;
t439 = t355 * t158 - t133;
t437 = -t157 + t550;
t434 = t622 * t347;
t433 = t623 * t347;
t430 = t314 - t464;
t427 = (-t353 ^ 2 - t655) * t584;
t281 = rSges(2,1) * t355 - rSges(2,2) * t353;
t279 = rSges(2,1) * t353 + rSges(2,2) * t355;
t280 = t579 - t631;
t408 = t101 * t355 + t102 * t353;
t407 = -t110 * t353 - t552;
t126 = -rSges(3,2) * t453 + (-t354 * t476 - t454) * rSges(3,1) + t486;
t237 = t278 * t353;
t127 = -qJD(2) * t237 + (t280 * t355 + t338) * qJD(1);
t406 = t126 * t355 + t127 * t353;
t78 = t159 * t308 + t161 * t307;
t401 = t182 * t353 + t183 * t355;
t398 = t213 * t353 + t214 * t355;
t392 = -t464 - t296;
t391 = -t226 + t261;
t198 = t226 * t355;
t390 = t405 * t353;
t184 = t610 * t347;
t381 = -qJD(2) * t306 - t184 - t463;
t380 = qJD(1) * t221 - t270 * t544 + t271 * t545;
t370 = qJD(1) * t157 - t307 * t441 + t308 * t443;
t11 = t353 * t603 + t370 * t355;
t369 = -t307 * t440 + t308 * t442 + t481;
t12 = t353 * t604 + t369 * t355;
t13 = t370 * t353 - t355 * t603;
t14 = t369 * t353 - t355 * t604;
t57 = -t390 - t523;
t58 = -t160 * t539 - t439;
t22 = t270 * t58 - t271 * t57 + t520;
t59 = -t159 * t538 - t515;
t60 = -t160 * t538 + t514;
t75 = -t355 * t397 + t545;
t69 = t75 * qJD(1);
t23 = t270 * t60 - t271 * t59 + t69;
t368 = qJD(1) * t220 - t307 * t434 + t308 * t433;
t34 = t353 * t600 + t368 * t355;
t35 = t368 * t353 - t355 * t600;
t371 = qJD(1) * t623 + t270 * t625 - t271 * t626;
t358 = -t307 * t595 + t371 * t308;
t39 = t307 * t443 + t308 * t441;
t40 = t307 * t442 + t308 * t440;
t79 = t160 * t308 + t162 * t307;
t374 = (qJD(1) * t34 + qJDD(1) * t75 - t11 * t271 + t12 * t270 + t187 * t60 + t188 * t59) * t586 + (t371 * t307 + t308 * t595) * t582 + (qJD(1) * t35 + qJDD(1) * t74 - t13 * t271 + t14 * t270 + t187 * t58 + t188 * t57) * t585 + t22 * t450 + t23 * t449 + (-t11 * t355 + t12 * t353 + (t353 * t59 + t355 * t60) * qJD(1)) * t589 + (t353 * t60 - t355 * t59) * t594 + (t353 * t58 - t355 * t57) * t593 + (-t13 * t355 + t14 * t353 + (t353 * t57 + t355 * t58) * qJD(1)) * t588 + (t353 * t79 - t355 * t78) * t572 + (t353 * t40 - t355 * t39 + (t353 * t78 + t355 * t79) * qJD(1)) * t581 + (t353 * t380 + t355 * t358) * t590 + (t353 * t358 - t355 * t380) * t587;
t36 = t165 * t270 + t166 * t271 + (-t130 * t353 + t131 * t355) * qJD(2) + t509;
t48 = -t474 * t583 - t226 * t270 + (t288 - t460) * qJD(1) - t485;
t373 = t36 * (-t270 * t197 - t198 * t271) + t48 * (-qJD(1) * t198 - t270 * t610);
t254 = t280 * qJD(2);
t240 = t355 * t261;
t239 = t353 * t261;
t212 = t248 * t355;
t190 = t240 + t467;
t189 = t239 + t468;
t156 = t271 * t610;
t105 = t398 * qJD(2);
t77 = t528 + t293 + t297 + (-t355 * t430 - t527) * qJD(1);
t76 = t429 + t215 + (t353 * t430 + t355 * t482) * qJD(1);
t64 = -t248 * t474 + (t288 - t502) * qJD(1) - t485;
t56 = qJD(2) * t401 + t509;
t55 = qJD(1) * t126 + qJDD(1) * t214 - t254 * t474 - t259 * t278 + t492;
t54 = -t254 * t473 + t260 * t278 + t497 * qJDD(1) + (-t127 - t257) * qJD(1);
t10 = -t184 * t270 - t187 * t226 + t516 * qJDD(1) + (t76 + t93) * qJD(1) + (-t259 * t319 - t353 * t531) * pkin(3) + t361;
t9 = -t355 * t466 - t184 * t271 + t188 * t226 + (t260 * t319 - t355 * t531) * pkin(3) + t428 * qJDD(1) + (-t77 - t94 + t518) * qJD(1) + t458;
t5 = -t130 * t259 + t165 * t187 - t166 * t188 + t270 * t94 + t271 * t93 + (-t131 - t171) * t260 + (t353 * t77 + t355 * t76) * qJD(2) + t462;
t1 = [(t69 + (t58 + (t159 * t355 + t160 * t353) * t307 + t439 + t515) * t271 + (-t161 * t537 + t523 + t57 + (t159 * t353 - t160 * t355) * t307 + t514) * t270) * t587 - m(2) * (-g(1) * t279 + g(2) * t281) + (t79 + t75) * t594 + (t78 + t74) * t593 + (-t520 + (t60 - t390 - t514) * t271 + (t353 * t437 - t133 + t59) * t270 + ((t158 + t405) * t270 + t437 * t271) * t355 + t22) * t590 + (t40 + t34) * t589 + ((((t667 + t671) * t355 + t634 + t665 + t669) * t355 - t607 * t353) * qJD(2) + t644) * t445 + (t39 + t35 + t23) * t588 + (-t659 * qJD(2) + t217 * t320 + t218 * t319 + t252 * t354 + t253 * t352 + t307 * t433 + t308 * t434) * qJD(1) + (-(qJD(1) * t517 + t372 - t47 + t613) * t48 + t47 * (t325 - t528) + t48 * (t491 + t612) + (-t198 * t48 + t226 * t574) * t347 + ((t47 * (t392 + t576) - t48 * t346) * t355 + (t47 * (t346 - rSges(5,3)) + t48 * t392) * t353) * qJD(1) + (t10 - g(2)) * (t166 + t230 - t527) + (t9 - g(1)) * (-t165 + t495)) * m(5) + (-(-qJD(1) * t182 + t382 + t613 - t63) * t64 + t63 * t457 + t64 * (t324 + t490) + (t353 * t573 + t614 * t64) * qJD(2) + ((-t63 * rSges(4,3) + t64 * (-t314 - t305)) * t353 + (t63 * (-t314 - t609) - t64 * t351) * t355) * qJD(1) + t642 * (t183 + t431) + t645 * (-t182 - t487)) * m(4) + (-(-qJD(1) * t213 - t109 - t263 - t455) * t110 + t110 * (t318 + t486) + (t278 * t553 - t551) * qJD(2) + ((-pkin(1) - t280) * t552 + (t109 * (-rSges(3,3) - pkin(5)) + t110 * t451) * t353) * qJD(1) + (t55 - g(2)) * t153 + (t54 - g(1)) * (t353 * t451 + t344 + t484)) * m(3) + (t627 + t629) * t592 + (t628 + t630) * t591 + (((t355 * t608 + t607 + t632) * t355 + (t353 * t608 + t633 + t666) * t353) * qJD(2) + t641 - t643) * t448 + (t636 + t639) * t447 + (m(2) * (t279 ^ 2 + t281 ^ 2) + t222 * t308 + t224 * t307 + Icges(2,3) - t658) * qJDD(1) + (-t637 + t638 + t640) * t446; t374 + t650 * t592 + t649 * t591 + (qJD(1) * t639 + qJD(2) * t618 + qJDD(1) * t627 + t259 * t632 + t260 * t633) * t586 + (qJD(1) * t638 + qJD(2) * t619 + qJDD(1) * t628 + t259 * t634 + t260 * t635) * t585 + ((t379 * t319 + t320 * t596 + t378 * t352 + t354 * t597) * qJD(2) + (t319 * t494 + t320 * t493 + t352 * t489 + t354 * t488) * qJD(1)) * t582 + (t637 * t355 + t636 * t353 + (t353 * t630 + t355 * t629) * qJD(1)) * t581 + (t353 * t629 - t355 * t630) * t572 + t641 * t450 + t640 * t449 + ((-t474 * t616 + t617) * t353 + ((t353 * t615 + t621) * qJD(2) + t620) * t355) * t448 + ((t353 * t633 + t355 * t632) * qJD(1) + t618) * t447 + ((t353 * t635 + t355 * t634) * qJD(1) + t619) * t446 + ((-t473 * t615 - t617) * t355 + ((t355 * t616 + t621) * qJD(2) + t620) * t353) * t445 + (-g(1) * (t240 - t198) - g(2) * (t239 - t197) - g(3) * (t610 + t483) + t5 * (t508 + t511) + t36 * (t461 + t465) + (t9 * t391 + t47 * t381 + t5 * t131 + t36 * t76 + (-t36 * t130 + t391 * t48) * qJD(1)) * t355 + (t10 * t391 + t48 * t381 - t5 * t130 + t36 * t77 + (t47 * (t226 + t583) + t36 * t460) * qJD(1)) * t353 + t47 * t156 - (t47 * (-t189 + t197) + t48 * (t190 - t467)) * qJD(1) - (t36 * t427 + (t36 * t190 - t47 * t483) * t355 + (t36 * t189 - t48 * t483) * t353) * qJD(2) - t373) * m(5) + (-g(3) * t611 + t63 * (-pkin(2) * t453 - t219 * t355) + (qJD(2) * t408 + t182 * t259 + t260 * t502 + t462) * (t401 + t508) + t56 * (t408 + t461) + (t56 * t182 + t383 * t64) * t475 - (t63 * t211 + t64 * (-t212 - t467)) * qJD(1) - (t56 * t427 + (-t56 * t212 - t611 * t63) * t355) * qJD(2) + (t64 * (-t219 - t463) + (t502 * t56 + t573) * qJD(1) - (-t56 * t211 - t611 * t64) * qJD(2) + t642 * t383) * t353 + t645 * t614) * m(4) + (g(1) * t238 + g(2) * t237 - g(3) * t280 - (t109 * t237 - t551) * qJD(1) - (t105 * (-t237 * t353 - t238 * t355) + t407 * t280) * qJD(2) + (qJD(2) * t406 + t213 * t259 - t214 * t260) * t398 + t105 * ((t213 * t355 - t214 * t353) * qJD(1) + t406) + t407 * t254 + (-t55 * t353 - t54 * t355 + (-t110 * t355 + t553) * qJD(1)) * t278) * m(3); (-m(4) - m(5)) * (g(1) * t353 - g(2) * t355) + 0.2e1 * (t10 * t585 + t586 * t9) * m(5) + 0.2e1 * (t30 * t586 + t31 * t585) * m(4); t374 + (g(1) * t198 + g(2) * t197 - g(3) * t610 - t47 * (qJD(1) * t197 - t156) - t373 + t5 * t511 + t36 * (-t166 * t476 + t465) + (-t353 * t48 - t355 * t47) * t184 + (-t10 * t353 - t9 * t355 + (-t355 * t48 + t574) * qJD(1)) * t226) * m(5);];
tau = t1;
