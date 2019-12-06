% Calculate vector of inverse dynamics joint torques for
% S5PRRPR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:03
% EndTime: 2019-12-05 16:19:39
% DurationCPUTime: 27.77s
% Computational Cost: add. (17884->756), mult. (13990->962), div. (0->0), fcn. (10957->8), ass. (0->412)
t672 = Icges(4,3) + Icges(5,3);
t355 = qJ(3) + pkin(9);
t346 = sin(t355);
t348 = cos(t355);
t254 = Icges(5,5) * t348 - Icges(5,6) * t346;
t357 = sin(qJ(3));
t358 = cos(qJ(3));
t290 = Icges(4,5) * t358 - Icges(4,6) * t357;
t660 = t254 + t290;
t353 = pkin(8) + qJ(2);
t347 = cos(t353);
t671 = t672 * t347;
t345 = sin(t353);
t535 = t345 * t358;
t536 = t345 * t357;
t538 = t345 * t348;
t539 = t345 * t346;
t661 = -Icges(4,5) * t535 - Icges(5,5) * t538 + Icges(4,6) * t536 + Icges(5,6) * t539 + t671;
t666 = t672 * t345 + t660 * t347;
t564 = Icges(5,6) * t347;
t180 = Icges(5,4) * t538 - Icges(5,2) * t539 - t564;
t565 = Icges(4,6) * t347;
t197 = Icges(4,4) * t535 - Icges(4,2) * t536 - t565;
t670 = t180 * t346 + t197 * t357;
t280 = Icges(5,4) * t539;
t570 = Icges(5,5) * t347;
t182 = Icges(5,1) * t538 - t280 - t570;
t300 = Icges(4,4) * t536;
t571 = Icges(4,5) * t347;
t199 = Icges(4,1) * t535 - t300 - t571;
t650 = -t182 * t348 - t199 * t358 + t670;
t631 = -t650 * t345 + t661 * t347;
t573 = Icges(5,4) * t346;
t258 = Icges(5,1) * t348 - t573;
t183 = Icges(5,5) * t345 + t258 * t347;
t574 = Icges(4,4) * t357;
t294 = Icges(4,1) * t358 - t574;
t200 = Icges(4,5) * t345 + t294 * t347;
t668 = -t183 * t538 - t200 * t535;
t253 = Icges(5,5) * t346 + Icges(5,6) * t348;
t289 = Icges(4,5) * t357 + Icges(4,6) * t358;
t667 = t253 + t289;
t255 = Icges(5,2) * t348 + t573;
t328 = Icges(5,4) * t348;
t257 = Icges(5,1) * t346 + t328;
t291 = Icges(4,2) * t358 + t574;
t350 = Icges(4,4) * t358;
t293 = Icges(4,1) * t357 + t350;
t658 = t255 * t346 - t257 * t348 + t291 * t357 - t293 * t358;
t665 = t666 * t347 + t668;
t529 = t347 * t358;
t531 = t347 * t348;
t615 = -t183 * t531 - t200 * t529 - t666 * t345;
t664 = -t182 * t531 - t199 * t529 + t661 * t345;
t411 = -Icges(5,2) * t346 + t328;
t181 = Icges(5,6) * t345 + t347 * t411;
t412 = -Icges(4,2) * t357 + t350;
t198 = Icges(4,6) * t345 + t347 * t412;
t662 = t181 * t346 + t198 * t357;
t630 = -t181 * t539 - t198 * t536 - t665;
t530 = t347 * t357;
t534 = t346 * t347;
t629 = -t180 * t534 - t197 * t530 - t664;
t628 = -t181 * t534 - t198 * t530 - t615;
t627 = t180 * t348 + t182 * t346 + t197 * t358 + t199 * t357;
t626 = t181 * t348 + t183 * t346 + t198 * t358 + t200 * t357;
t659 = t667 * qJD(3);
t657 = -t255 * t348 - t257 * t346 - t291 * t358 - t293 * t357;
t656 = t183 * t348 + t200 * t358 - t662;
t544 = t289 * t347;
t546 = t253 * t347;
t625 = -t658 * t345 - t544 - t546;
t545 = t289 * t345;
t547 = t253 * t345;
t624 = -t658 * t347 + t545 + t547;
t655 = t666 * qJD(2);
t654 = t347 ^ 2;
t231 = t411 * qJD(3);
t232 = t258 * qJD(3);
t275 = t412 * qJD(3);
t276 = t294 * qJD(3);
t653 = qJD(2) * t667 + t657 * qJD(3) - t231 * t346 + t232 * t348 - t275 * t357 + t276 * t358;
t387 = qJD(3) * t255;
t104 = -t347 * t387 + (-t345 * t411 + t564) * qJD(2);
t389 = qJD(3) * t257;
t106 = -t347 * t389 + (-t258 * t345 + t570) * qJD(2);
t388 = qJD(3) * t291;
t122 = -t347 * t388 + (-t345 * t412 + t565) * qJD(2);
t390 = qJD(3) * t293;
t124 = -t347 * t390 + (-t294 * t345 + t571) * qJD(2);
t652 = -t626 * qJD(3) - t104 * t346 + t106 * t348 - t122 * t357 + t124 * t358 + t655;
t105 = qJD(2) * t181 - t345 * t387;
t107 = qJD(2) * t183 - t345 * t389;
t123 = qJD(2) * t198 - t345 * t388;
t125 = qJD(2) * t200 - t345 * t390;
t651 = t661 * qJD(2) + t627 * qJD(3) + t105 * t346 - t107 * t348 + t123 * t357 - t125 * t358;
t649 = t628 * t345 - t629 * t347;
t648 = t630 * t345 - t631 * t347;
t647 = t658 * qJD(2) + t660 * qJD(3);
t646 = t650 * qJD(2) - t659 * t345 + t655;
t645 = -t659 * t347 + (-t660 * t345 - t656 + t671) * qJD(2);
t644 = rSges(4,2) * t357;
t643 = t624 * qJD(2);
t349 = qJ(5) + t355;
t341 = cos(t349);
t541 = t341 * t345;
t340 = sin(t349);
t543 = t340 * t345;
t563 = Icges(6,6) * t347;
t167 = Icges(6,4) * t541 - Icges(6,2) * t543 - t563;
t312 = Icges(6,4) * t341;
t242 = Icges(6,1) * t340 + t312;
t642 = -t242 * t345 - t167;
t410 = -Icges(6,2) * t340 + t312;
t168 = Icges(6,6) * t345 + t347 * t410;
t641 = -t242 * t347 - t168;
t640 = t242 + t410;
t639 = t625 * qJD(2);
t476 = qJD(3) * t348;
t477 = qJD(3) * t347;
t480 = qJD(2) * t345;
t281 = rSges(5,2) * t539;
t479 = qJD(2) * t347;
t494 = rSges(5,3) * t479 + qJD(2) * t281;
t110 = -rSges(5,2) * t347 * t476 + (-t346 * t477 - t348 * t480) * rSges(5,1) + t494;
t330 = t345 * rSges(5,3);
t185 = rSges(5,1) * t531 - rSges(5,2) * t534 + t330;
t335 = t348 * rSges(5,1);
t619 = -rSges(5,2) * t346 + t335;
t235 = t619 * qJD(3);
t473 = qJD(2) * qJD(3);
t247 = qJDD(3) * t345 + t347 * t473;
t262 = rSges(5,1) * t346 + rSges(5,2) * t348;
t356 = -qJ(4) - pkin(6);
t308 = t347 * t356;
t313 = pkin(6) * t479;
t317 = qJD(4) * t345;
t455 = t357 * t477;
t431 = pkin(3) * t455;
t351 = t358 * pkin(3);
t342 = t351 + pkin(2);
t586 = pkin(2) - t342;
t117 = -t431 - t313 + t317 + (t345 * t586 - t308) * qJD(2);
t336 = t345 * pkin(6);
t266 = t347 * pkin(2) + t336;
t285 = t347 * t342;
t433 = -t345 * t356 + t285;
t164 = t433 - t266;
t472 = qJD(2) * qJD(4);
t499 = qJD(2) * (-pkin(2) * t480 + t313) + qJDD(2) * t266;
t359 = qJD(3) ^ 2;
t527 = t358 * t359;
t362 = qJD(2) * t117 + qJDD(2) * t164 + t345 * t472 + (-t247 * t357 - t345 * t527) * pkin(3) - qJDD(4) * t347 + t499;
t478 = qJD(3) * t345;
t32 = qJD(2) * t110 + qJDD(2) * t185 - t235 * t478 - t247 * t262 + t362;
t638 = -g(2) + t32;
t637 = t649 * qJD(3) + t643;
t636 = t648 * qJD(3) + t639;
t635 = t650 * qJD(3) - t105 * t348 - t107 * t346 - t123 * t358 - t125 * t357;
t634 = t656 * qJD(3) + t104 * t348 + t106 * t346 + t122 * t358 + t124 * t357;
t633 = t647 * t345 + t653 * t347;
t632 = t653 * t345 - t647 * t347;
t337 = t347 * pkin(6);
t265 = pkin(2) * t345 - t337;
t493 = t345 * t342 + t308;
t163 = t265 - t493;
t250 = qJD(2) * t265;
t623 = qJD(2) * t163 - t250;
t589 = pkin(4) * t346;
t590 = pkin(3) * t357;
t283 = -t589 - t590;
t259 = t283 * qJD(3);
t218 = t347 * t259;
t622 = t218 + t317;
t621 = t619 + t351;
t314 = t341 * rSges(6,1);
t583 = rSges(6,2) * t340;
t620 = t314 - t583;
t487 = pkin(4) * t348 + t351;
t618 = t646 * t654 + (t652 * t345 + (-t645 + t651) * t347) * t345;
t617 = t651 * t654 + (t645 * t345 + (-t646 + t652) * t347) * t345;
t616 = t661 + t662;
t572 = Icges(6,4) * t340;
t243 = Icges(6,1) * t341 - t572;
t170 = Icges(6,5) * t345 + t243 * t347;
t239 = Icges(6,5) * t341 - Icges(6,6) * t340;
t354 = qJD(3) + qJD(5);
t238 = Icges(6,5) * t340 + Icges(6,6) * t341;
t549 = t238 * t347;
t556 = t168 * t340;
t560 = Icges(6,3) * t347;
t614 = -t354 * t549 + (-t170 * t341 - t239 * t345 + t556 + t560) * qJD(2);
t271 = Icges(6,4) * t543;
t569 = Icges(6,5) * t347;
t169 = Icges(6,1) * t541 - t271 - t569;
t407 = t167 * t340 - t169 * t341;
t166 = Icges(6,3) * t345 + t239 * t347;
t485 = qJD(2) * t166;
t550 = t238 * t345;
t613 = qJD(2) * t407 - t354 * t550 + t485;
t240 = Icges(6,2) * t341 + t572;
t400 = t240 * t340 - t242 * t341;
t608 = qJD(2) * t400 + t239 * t354;
t505 = -t255 * t347 + t183;
t506 = -Icges(5,2) * t538 + t182 - t280;
t605 = t345 * t505 - t347 * t506;
t502 = -Icges(4,2) * t535 + t199 - t300;
t504 = t293 * t345 + t197;
t604 = -t357 * t502 - t358 * t504;
t251 = t345 * t354;
t252 = t347 * t354;
t603 = qJD(2) * t640 + t251 * (-t240 * t347 + t170) - t252 * (-Icges(6,2) * t541 + t169 - t271);
t602 = m(2) + m(3);
t601 = -m(5) - m(6);
t471 = qJD(2) * qJD(5);
t175 = qJDD(5) * t345 + t347 * t471 + t247;
t600 = t175 / 0.2e1;
t306 = t345 * t473;
t176 = t345 * t471 + t306 + (-qJDD(3) - qJDD(5)) * t347;
t599 = t176 / 0.2e1;
t598 = t247 / 0.2e1;
t248 = -qJDD(3) * t347 + t306;
t597 = t248 / 0.2e1;
t596 = -t251 / 0.2e1;
t595 = t251 / 0.2e1;
t594 = -t252 / 0.2e1;
t593 = t252 / 0.2e1;
t592 = t345 / 0.2e1;
t591 = -t347 / 0.2e1;
t588 = -qJD(2) / 0.2e1;
t587 = qJD(2) / 0.2e1;
t585 = rSges(4,1) * t358;
t582 = rSges(6,2) * t341;
t331 = t345 * rSges(4,3);
t329 = t345 * rSges(6,3);
t245 = rSges(6,1) * t340 + t582;
t373 = -t252 * t245 + t622;
t512 = t163 - t265;
t352 = -pkin(7) + t356;
t466 = pkin(2) + t487;
t498 = -t345 * t466 - t347 * t352;
t133 = t493 + t498;
t272 = rSges(6,2) * t543;
t171 = rSges(6,1) * t541 - t347 * rSges(6,3) - t272;
t523 = t133 - t171;
t429 = t512 + t523;
t49 = qJD(2) * t429 + t373;
t581 = t345 * t49;
t384 = -t262 - t590;
t383 = t384 * t477 + t317;
t184 = rSges(5,1) * t538 - t347 * rSges(5,3) - t281;
t461 = -t184 + t512;
t65 = qJD(2) * t461 + t383;
t580 = t65 * t262;
t578 = qJDD(2) / 0.2e1;
t296 = rSges(4,1) * t357 + rSges(4,2) * t358;
t457 = t296 * t477;
t488 = rSges(4,2) * t536 + t347 * rSges(4,3);
t203 = rSges(4,1) * t535 - t488;
t500 = -t203 - t265;
t108 = qJD(2) * t500 - t457;
t559 = t108 * t345;
t558 = t108 * t347;
t204 = rSges(4,1) * t529 - rSges(4,2) * t530 + t331;
t153 = t204 + t266;
t109 = qJD(2) * t153 - t296 * t478;
t226 = t296 * t347;
t557 = t109 * t226;
t548 = t240 * t354;
t542 = t340 * t347;
t540 = t341 * t347;
t537 = t345 * t352;
t165 = Icges(6,5) * t541 - Icges(6,6) * t543 - t560;
t533 = t347 * t165;
t528 = t348 * t359;
t80 = -t345 * t400 - t549;
t526 = t80 * qJD(2);
t295 = t356 * t480;
t470 = pkin(3) * t536;
t287 = qJD(3) * t470;
t318 = qJD(4) * t347;
t491 = t287 + t318;
t459 = t295 + t491;
t118 = (-t347 * t586 - t336) * qJD(2) - t459;
t244 = t266 * qJD(2);
t524 = -t118 - t244;
t229 = t347 * t466;
t486 = t356 - t352;
t134 = t345 * t486 + t229 - t285;
t172 = rSges(6,1) * t540 - rSges(6,2) * t542 + t329;
t522 = t134 + t172;
t521 = -t345 * t165 - t169 * t540;
t520 = t345 * t166 + t170 * t540;
t516 = -t345 * t163 + t347 * t164;
t515 = t345 * t171 + t347 * t172;
t511 = -t164 - t185;
t508 = -t257 * t345 - t180;
t507 = -t257 * t347 - t181;
t503 = -t293 * t347 - t198;
t501 = -t291 * t347 + t200;
t496 = -t255 + t258;
t495 = t257 + t411;
t492 = rSges(4,3) * t479 + t480 * t644;
t490 = -t291 + t294;
t489 = t293 + t412;
t482 = qJD(2) * t254;
t481 = qJD(2) * t290;
t475 = qJD(3) * t358;
t469 = pkin(3) * t530;
t468 = pkin(3) * t527;
t393 = rSges(6,3) * t479 + qJD(2) * t272 - t252 * t582;
t464 = t340 * t252;
t95 = (-t341 * t480 - t464) * rSges(6,1) + t393;
t201 = t245 * t345;
t96 = -t354 * t201 + (t347 * t620 + t329) * qJD(2);
t467 = t171 * t479 + t345 * t96 + t347 * t95;
t465 = pkin(3) * t475;
t463 = t347 * t117 + t345 * t118 - t163 * t479;
t462 = -t164 - t522;
t460 = qJDD(4) * t345 + t248 * t590 + t347 * t472;
t454 = -t163 * t478 + t164 * t477 + qJD(1);
t453 = -pkin(2) - t585;
t452 = t480 / 0.2e1;
t451 = t479 / 0.2e1;
t450 = -t478 / 0.2e1;
t449 = t478 / 0.2e1;
t448 = -t477 / 0.2e1;
t447 = t477 / 0.2e1;
t445 = qJD(2) * t170 + t642 * t354;
t444 = (-t243 * t345 + t569) * qJD(2) + t641 * t354;
t443 = qJD(2) * t168 + t169 * t354 - t345 * t548;
t442 = t170 * t354 - t347 * t548 + (-t345 * t410 + t563) * qJD(2);
t135 = t170 * t541;
t441 = t347 * t166 - t135;
t438 = -t165 + t556;
t435 = t640 * t354;
t434 = t243 * t354 - t548;
t432 = t342 - t466;
t430 = t117 * t477 + t118 * t478 - t247 * t163 + qJDD(1);
t428 = (-t345 ^ 2 - t654) * t590;
t427 = -t235 - t465;
t263 = rSges(3,1) * t347 - rSges(3,2) * t345;
t261 = rSges(3,1) * t345 + rSges(3,2) * t347;
t297 = t585 - t644;
t409 = -t109 * t345 - t558;
t128 = -rSges(4,2) * t347 * t475 + (-t358 * t480 - t455) * rSges(4,1) + t492;
t225 = t296 * t345;
t129 = -qJD(3) * t225 + (t297 * t347 + t331) * qJD(2);
t408 = t128 * t347 + t129 * t345;
t78 = t167 * t341 + t169 * t340;
t401 = t203 * t345 + t204 * t347;
t395 = -t466 - t314;
t394 = -t245 + t283;
t211 = t262 * t345;
t202 = t245 * t347;
t391 = t407 * t345;
t216 = t620 * t354;
t382 = -pkin(4) * t476 - t216 - t465;
t381 = qJD(2) * t239 - t251 * t549 + t252 * t550;
t380 = t345 * t507 - t347 * t508;
t379 = -t357 * t501 + t358 * t503;
t378 = t384 * t347;
t377 = (-t346 * t495 + t348 * t496) * qJD(2);
t376 = (-t357 * t489 + t358 * t490) * qJD(2);
t371 = qJD(2) * t165 - t340 * t443 + t341 * t445;
t11 = t345 * t613 + t371 * t347;
t370 = -t340 * t442 + t341 * t444 + t485;
t12 = t345 * t614 + t370 * t347;
t13 = t371 * t345 - t347 * t613;
t14 = t370 * t345 - t347 * t614;
t59 = -t391 - t533;
t60 = -t168 * t543 - t441;
t23 = t251 * t60 - t252 * t59 + t526;
t61 = -t167 * t542 - t521;
t62 = -t168 * t542 + t520;
t81 = -t347 * t400 + t550;
t77 = t81 * qJD(2);
t24 = t251 * t62 - t252 * t61 + t77;
t372 = t641 * t251 - t642 * t252 + (-t240 + t243) * qJD(2);
t361 = -t340 * t603 + t372 * t341;
t369 = qJD(2) * t238 - t340 * t435 + t341 * t434;
t38 = t345 * t608 + t369 * t347;
t39 = t369 * t345 - t347 * t608;
t40 = t340 * t445 + t341 * t443;
t41 = t340 * t444 + t341 * t442;
t79 = t168 * t341 + t170 * t340;
t375 = (qJD(2) * t38 + qJDD(2) * t81 - t11 * t252 + t12 * t251 + t175 * t62 + t176 * t61) * t592 + (t372 * t340 + t603 * t341) * t588 + (qJD(2) * t39 + qJDD(2) * t80 - t13 * t252 + t14 * t251 + t175 * t60 + t176 * t59) * t591 + t23 * t452 + t24 * t451 + (-t11 * t347 + t12 * t345 + (t345 * t61 + t347 * t62) * qJD(2)) * t595 + (t345 * t62 - t347 * t61) * t600 + (t345 * t60 - t347 * t59) * t599 + (-t13 * t347 + t14 * t345 + (t345 * t59 + t347 * t60) * qJD(2)) * t594 + (t345 * t79 - t347 * t78) * t578 + (t345 * t41 - t347 * t40 + (t345 * t78 + t347 * t79) * qJD(2)) * t587 + (t345 * t381 + t347 * t361) * t596 + (t345 * t361 - t347 * t381) * t593;
t35 = t171 * t251 + t172 * t252 + (-t133 * t345 + t134 * t347) * qJD(3) + t454;
t50 = -t478 * t589 - t245 * t251 + (t266 - t462) * qJD(2) - t491;
t374 = t50 * (-qJD(2) * t202 - t251 * t620) + t35 * (-t251 * t201 - t202 * t252);
t277 = t297 * qJD(3);
t234 = t347 * t283;
t233 = t345 * t283;
t212 = t262 * t347;
t188 = t234 + t469;
t187 = t233 + t470;
t177 = t252 * t620;
t111 = -qJD(3) * t211 + (t347 * t619 + t330) * qJD(2);
t97 = qJD(3) * t401 + qJD(1);
t76 = t345 * t259 + t287 + t295 + (-t347 * t432 - t537) * qJD(2);
t75 = t431 + t218 + (t345 * t432 + t347 * t486) * qJD(2);
t66 = -t262 * t478 + (t266 - t511) * qJD(2) - t491;
t56 = (t184 * t345 + t185 * t347) * qJD(3) + t454;
t55 = qJD(2) * t128 + qJDD(2) * t204 - t247 * t296 - t277 * t478 + t499;
t54 = -t277 * t477 + t248 * t296 + t500 * qJDD(2) + (-t129 - t244) * qJD(2);
t48 = qJD(3) * t408 + t203 * t247 - t204 * t248 + qJDD(1);
t31 = t248 * t262 + (-qJD(3) * t235 - t468) * t347 + t461 * qJDD(2) + (-t111 + t524) * qJD(2) + t460;
t15 = t184 * t247 + t511 * t248 + (t110 * t347 + t111 * t345) * qJD(3) + t430;
t10 = -t175 * t245 - t216 * t251 + t522 * qJDD(2) + (t75 + t95) * qJD(2) + (-t247 * t346 - t345 * t528) * pkin(4) + t362;
t9 = -t347 * t468 + t176 * t245 - t216 * t252 + (t248 * t346 - t347 * t528) * pkin(4) + t429 * qJDD(2) + (-t76 - t96 + t524) * qJD(2) + t460;
t5 = -t133 * t247 + t171 * t175 - t172 * t176 + t251 * t96 + t252 * t95 + (-t134 - t164) * t248 + (t345 * t76 + t347 * t75) * qJD(3) + t430;
t1 = [t602 * qJDD(1) + m(4) * t48 + m(5) * t15 + m(6) * t5 + (-m(4) + t601 - t602) * g(3); (t77 + (t60 + (t167 * t347 + t168 * t345) * t340 + t441 + t521) * t252 + (-t169 * t541 + t533 + t59 + (t167 * t345 - t168 * t347) * t340 + t520) * t251) * t593 - m(3) * (-g(1) * t261 + g(2) * t263) + (t79 + t81) * t600 + (t78 + t80) * t599 + (-t526 + (t62 - t391 - t520) * t252 + (t345 * t438 - t135 + t61) * t251 + ((t166 + t407) * t251 + t438 * t252) * t347 + t23) * t596 + (t41 + t38) * t595 + ((((t666 + t670) * t347 + t630 + t664 + t668) * t347 - t615 * t345) * qJD(3) + t643) * t447 + (t40 + t39 + t24) * t594 + (-t658 * qJD(3) + t231 * t348 + t232 * t346 + t275 * t358 + t276 * t357 + t340 * t434 + t341 * t435) * qJD(2) + (t49 * t318 + t50 * (-rSges(6,1) * t464 + t393 + t622) + (t245 * t354 - t259) * t581 + ((t49 * (t395 + t583) - t50 * t352) * t347 + (t49 * (t352 - rSges(6,3)) + t50 * t395) * t345) * qJD(2) - (qJD(2) * t523 + t373 - t49 + t623) * t50 + (t10 - g(2)) * (t172 + t229 - t537) + (t9 - g(1)) * (-t171 + t498)) * m(6) + (t65 * t459 + t66 * (t317 + t494) + (t345 * t580 + t378 * t66) * qJD(3) + ((-t65 * rSges(5,3) + t66 * (-t342 - t335)) * t345 + (t65 * (-t342 - t619) - t66 * t356) * t347) * qJD(2) - (-qJD(2) * t184 + t383 + t623 - t65) * t66 + t638 * (t185 + t433) + (t31 - g(1)) * (-t184 - t493)) * m(5) + (-(-qJD(2) * t203 - t108 - t250 - t457) * t109 + t109 * (t313 + t492) + (t296 * t559 - t557) * qJD(3) + ((-pkin(2) - t297) * t558 + (t108 * (-rSges(4,3) - pkin(6)) + t109 * t453) * t345) * qJD(2) + (t55 - g(2)) * t153 + (t54 - g(1)) * (t453 * t345 + t337 + t488)) * m(4) + (t624 + t626) * t598 + (t625 + t627) * t597 + (((t347 * t616 + t615 + t628) * t347 + (t345 * t616 + t629 + t665) * t345) * qJD(3) + t636 - t639) * t450 + (t633 + t634) * t449 + (m(3) * (t261 ^ 2 + t263 ^ 2) + t240 * t341 + t242 * t340 + Icges(3,3) - t657) * qJDD(2) + (t632 - t635 + t637) * t448; t375 + t649 * t598 + t648 * t597 + (qJD(2) * t633 + qJD(3) * t617 + qJDD(2) * t624 + t247 * t628 + t248 * t629) * t592 + (qJD(2) * t632 + qJD(3) * t618 + qJDD(2) * t625 + t247 * t630 + t248 * t631) * t591 + (((t345 * t501 - t347 * t502) * t358 + (t345 * t503 + t347 * t504) * t357 + t605 * t348 + t380 * t346) * qJD(3) + (t346 * t496 + t348 * t495 + t357 * t490 + t358 * t489) * qJD(2)) * t588 + (t635 * t347 + t634 * t345 + (t345 * t627 + t347 * t626) * qJD(2)) * t587 + (t345 * t626 - t347 * t627) * t578 + t636 * t452 + t637 * t451 + ((-t478 * t546 + t482) * t345 + (t377 + (-(-t346 * t506 + t348 * t508) * t347 + (-t346 * t505 + t348 * t507 + t547) * t345) * qJD(3)) * t347 + (-t478 * t544 + t481) * t345 + (t376 + (-t604 * t347 + (t545 + t379) * t345) * qJD(3)) * t347) * t450 + ((t345 * t629 + t347 * t628) * qJD(2) + t617) * t449 + ((t345 * t631 + t347 * t630) * qJD(2) + t618) * t448 + ((-t477 * t547 - t482) * t347 + (t377 + (-t346 * t605 + t347 * t546 + t380 * t348) * qJD(3)) * t345 + (-t477 * t545 - t481) * t347 + (t376 + (t379 * t345 + (t544 - t604) * t347) * qJD(3)) * t345) * t447 + (t49 * t177 - (t49 * (-t187 + t201) + t50 * (t188 - t469)) * qJD(2) - (t35 * t428 + (t35 * t188 - t487 * t49) * t347 + (t35 * t187 - t487 * t50) * t345) * qJD(3) - t374 - g(1) * (t234 - t202) - g(2) * (t233 - t201) - g(3) * (t620 + t487) + t5 * (t515 + t516) + t35 * (t463 + t467) + (t9 * t394 + t49 * t382 + t5 * t134 + t35 * t75 + (-t35 * t133 + t394 * t50) * qJD(2)) * t347 + (t10 * t394 + t50 * t382 - t5 * t133 + t35 * t76 + (t49 * (t245 + t589) + t35 * t462) * qJD(2)) * t345) * m(6) + (-(t108 * t225 - t557) * qJD(2) - (t97 * (-t225 * t345 - t226 * t347) + t409 * t297) * qJD(3) + t48 * t401 + t97 * ((t203 * t347 - t204 * t345) * qJD(2) + t408) + t409 * t277 + (-t55 * t345 - t54 * t347 + (-t109 * t347 + t559) * qJD(2)) * t296 + g(1) * t226 + g(2) * t225 - g(3) * t297) * m(4) + (-g(3) * t621 - g(1) * t378 + t15 * t516 + t56 * t463 + (t31 * t384 + t65 * t427 + t15 * t185 + t56 * t110 + (t56 * t184 + t384 * t66) * qJD(2)) * t347 - (t65 * t211 + t66 * (-t212 - t469)) * qJD(2) - (t56 * t428 + (-t56 * t212 - t621 * t65) * t347) * qJD(3) + (t66 * t427 + t15 * t184 + t56 * t111 + (t511 * t56 + t580) * qJD(2) - (-t56 * t211 - t621 * t66) * qJD(3) + t638 * t384) * t345) * m(5); t601 * (g(1) * t345 - g(2) * t347) + 0.2e1 * (t10 * t591 + t592 * t9) * m(6) + 0.2e1 * (t31 * t592 + t32 * t591) * m(5); t375 + (-t49 * (qJD(2) * t201 - t177) - t374 + g(1) * t202 + g(2) * t201 - g(3) * t620 + t5 * t515 + t35 * (-t172 * t480 + t467) + (-t345 * t50 - t347 * t49) * t216 + (-t10 * t345 - t9 * t347 + (-t347 * t50 + t581) * qJD(2)) * t245) * m(6);];
tau = t1;
