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
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:38:15
% EndTime: 2020-01-03 11:39:11
% DurationCPUTime: 32.52s
% Computational Cost: add. (18025->811), mult. (14294->985), div. (0->0), fcn. (11105->10), ass. (0->423)
t688 = Icges(4,3) + Icges(5,3);
t349 = qJ(3) + pkin(9);
t337 = sin(t349);
t339 = cos(t349);
t247 = Icges(5,5) * t339 - Icges(5,6) * t337;
t352 = sin(qJ(3));
t354 = cos(qJ(3));
t292 = Icges(4,5) * t354 - Icges(4,6) * t352;
t683 = t247 + t292;
t350 = qJ(1) + pkin(8);
t340 = cos(t350);
t528 = t339 * t340;
t532 = t337 * t340;
t338 = sin(t350);
t556 = Icges(5,6) * t338;
t169 = Icges(5,4) * t528 - Icges(5,2) * t532 + t556;
t524 = t340 * t354;
t525 = t340 * t352;
t557 = Icges(4,6) * t338;
t187 = Icges(4,4) * t524 - Icges(4,2) * t525 + t557;
t687 = -t169 * t337 - t187 * t352;
t686 = t688 * t338;
t634 = t683 * t338 - t340 * t688;
t633 = Icges(4,5) * t524 + Icges(5,5) * t528 - Icges(4,6) * t525 - Icges(5,6) * t532 + t686;
t280 = Icges(5,4) * t532;
t562 = Icges(5,5) * t338;
t171 = Icges(5,1) * t528 - t280 + t562;
t303 = Icges(4,4) * t525;
t563 = Icges(4,5) * t338;
t189 = Icges(4,1) * t524 - t303 + t563;
t673 = t171 * t339 + t189 * t354 + t687;
t565 = Icges(5,4) * t337;
t251 = Icges(5,1) * t339 - t565;
t170 = -Icges(5,5) * t340 + t251 * t338;
t566 = Icges(4,4) * t352;
t296 = Icges(4,1) * t354 - t566;
t188 = -Icges(4,5) * t340 + t296 * t338;
t529 = t338 * t354;
t531 = t338 * t339;
t636 = -t170 * t531 - t188 * t529;
t685 = t171 * t531 + t189 * t529;
t321 = Icges(5,4) * t339;
t405 = -Icges(5,2) * t337 + t321;
t168 = -Icges(5,6) * t340 + t338 * t405;
t342 = Icges(4,4) * t354;
t406 = -Icges(4,2) * t352 + t342;
t186 = -Icges(4,6) * t340 + t338 * t406;
t635 = t168 * t532 + t186 * t525;
t248 = Icges(5,2) * t339 + t565;
t293 = Icges(4,2) * t354 + t566;
t625 = Icges(4,1) * t352 + t342;
t627 = Icges(5,1) * t337 + t321;
t681 = -t248 * t337 - t293 * t352 + t339 * t627 + t354 * t625;
t621 = -t170 * t339 - t188 * t354;
t620 = t168 * t337 + t186 * t352;
t530 = t338 * t352;
t533 = t337 * t338;
t647 = -t168 * t533 - t186 * t530 - t634 * t340 - t636;
t646 = t169 * t533 + t187 * t530 + t633 * t340 - t685;
t645 = -t170 * t528 - t188 * t524 - t634 * t338 + t635;
t644 = t633 * t338 + t673 * t340;
t643 = t168 * t339 + t170 * t337 + t186 * t354 + t188 * t352;
t246 = Icges(5,5) * t337 + Icges(5,6) * t339;
t196 = t246 * t340;
t291 = Icges(4,5) * t352 + Icges(4,6) * t354;
t212 = t291 * t340;
t684 = t196 + t212;
t682 = t248 * t339 + t293 * t354 + t337 * t627 + t352 * t625;
t642 = t169 * t339 + t171 * t337 + t187 * t354 + t189 * t352;
t680 = t620 + t621;
t538 = t291 * t338;
t540 = t246 * t338;
t679 = -t538 - t540;
t641 = t681 * t338 - t684;
t640 = -t248 * t532 - t293 * t525 + t524 * t625 + t528 * t627 - t679;
t678 = t634 * qJD(1);
t677 = t340 ^ 2;
t676 = t681 * qJD(1) - t683 * qJD(3);
t675 = t679 * qJD(3) + (t683 * t340 + t680 + t686) * qJD(1);
t674 = t673 * qJD(1) + t684 * qJD(3) + t678;
t672 = -t644 * t338 - t645 * t340;
t671 = -t646 * t338 - t647 * t340;
t353 = sin(qJ(1));
t344 = t353 * pkin(1);
t422 = -t338 * rSges(3,1) - t340 * rSges(3,2);
t670 = t422 - t344;
t476 = qJD(3) * t340;
t103 = qJD(1) * t168 + t248 * t476;
t200 = t627 * t340;
t105 = qJD(1) * t170 + qJD(3) * t200;
t119 = qJD(1) * t186 + t293 * t476;
t216 = t625 * t340;
t121 = qJD(1) * t188 + qJD(3) * t216;
t669 = t633 * qJD(1) - t642 * qJD(3) + t103 * t337 - t105 * t339 + t119 * t352 - t121 * t354;
t227 = t405 * qJD(3);
t228 = t251 * qJD(3);
t272 = t406 * qJD(3);
t273 = t296 * qJD(3);
t668 = t227 * t337 - t228 * t339 + t272 * t352 - t273 * t354 + t682 * qJD(3) + (-t246 - t291) * qJD(1);
t197 = t248 * t338;
t104 = -qJD(3) * t197 + (t340 * t405 + t556) * qJD(1);
t199 = t627 * t338;
t106 = -qJD(3) * t199 + (t251 * t340 + t562) * qJD(1);
t213 = t293 * t338;
t120 = -qJD(3) * t213 + (t340 * t406 + t557) * qJD(1);
t215 = t625 * t338;
t122 = -qJD(3) * t215 + (t296 * t340 + t563) * qJD(1);
t667 = -t643 * qJD(3) - t104 * t337 + t106 * t339 - t120 * t352 + t122 * t354 + t678;
t341 = qJ(5) + t349;
t328 = sin(t341);
t329 = cos(t341);
t536 = t328 * t340;
t264 = Icges(6,4) * t536;
t534 = t329 * t340;
t561 = Icges(6,5) * t338;
t160 = Icges(6,1) * t534 - t264 + t561;
t348 = qJD(3) + qJD(5);
t244 = t338 * t348;
t526 = t340 * t348;
t315 = Icges(6,4) * t329;
t404 = -Icges(6,2) * t328 + t315;
t629 = Icges(6,1) * t328 + t315;
t656 = t629 + t404;
t564 = Icges(6,4) * t328;
t236 = Icges(6,1) * t329 - t564;
t159 = -Icges(6,5) * t340 + t236 * t338;
t233 = Icges(6,2) * t329 + t564;
t659 = -t233 * t338 + t159;
t359 = qJD(1) * t656 + t244 * (-Icges(6,2) * t534 + t160 - t264) - t526 * t659;
t658 = t233 - t236;
t555 = Icges(6,6) * t338;
t158 = Icges(6,4) * t534 - Icges(6,2) * t536 + t555;
t660 = t629 * t340 + t158;
t157 = -Icges(6,6) * t340 + t338 * t404;
t661 = t629 * t338 + t157;
t605 = qJD(1) * t658 + t244 * t660 - t526 * t661;
t666 = t359 * t328 + t329 * t605;
t375 = t338 * (-Icges(5,2) * t528 + t171 - t280) - t340 * (t170 - t197);
t610 = t338 * (t169 + t200) - t340 * (t168 + t199);
t665 = -t375 * t337 - t339 * t610;
t585 = rSges(5,1) * t337;
t254 = rSges(5,2) * t339 + t585;
t590 = pkin(3) * t352;
t444 = t254 + t590;
t355 = cos(qJ(1));
t346 = t355 * pkin(1);
t334 = qJD(1) * t346;
t473 = qJD(4) * t340;
t434 = t334 - t473;
t478 = qJD(3) * t338;
t361 = -t444 * t478 + t434;
t591 = pkin(2) * t340;
t259 = pkin(6) * t338 + t591;
t345 = t354 * pkin(3);
t330 = t345 + pkin(2);
t287 = t340 * t330;
t351 = -qJ(4) - pkin(6);
t154 = t591 - t287 + (pkin(6) + t351) * t338;
t464 = rSges(5,1) * t528;
t428 = -rSges(5,2) * t532 + t464;
t173 = rSges(5,3) * t338 + t428;
t517 = -t154 + t173;
t454 = t259 + t517;
t64 = qJD(1) * t454 + t361;
t664 = t444 * t64;
t663 = t641 * qJD(1);
t581 = rSges(6,2) * t329;
t584 = rSges(6,1) * t328;
t238 = t581 + t584;
t662 = t238 * t526;
t657 = t640 * qJD(1);
t492 = t625 + t406;
t493 = t293 - t296;
t655 = (t352 * t492 + t354 * t493) * qJD(1);
t498 = t627 + t405;
t499 = t248 - t251;
t654 = (t337 * t498 + t339 * t499) * qJD(1);
t653 = qJD(3) * t671 + t663;
t652 = qJD(3) * t672 - t657;
t651 = t680 * qJD(3) - t104 * t339 - t106 * t337 - t120 * t354 - t122 * t352;
t650 = qJD(3) * t673 - t103 * t339 - t105 * t337 - t119 * t354 - t121 * t352;
t649 = t338 * t676 + t340 * t668;
t648 = -t338 * t668 + t340 * t676;
t232 = Icges(6,5) * t329 - Icges(6,6) * t328;
t155 = -Icges(6,3) * t340 + t232 * t338;
t549 = t158 * t328;
t402 = -t160 * t329 + t549;
t388 = -t155 + t402;
t639 = t526 * t388;
t243 = qJD(1) * t259;
t632 = -qJD(1) * t154 + t243;
t583 = rSges(5,2) * t337;
t626 = t339 * rSges(5,1) - t583;
t631 = t626 + t345;
t202 = rSges(5,1) * t532 + rSges(5,2) * t528;
t307 = pkin(3) * t525;
t630 = t307 + t202;
t582 = rSges(6,2) * t328;
t628 = t329 * rSges(6,1) - t582;
t624 = -pkin(4) * t339 - t345;
t623 = t675 * t677 + (t669 * t338 + (-t667 + t674) * t340) * t338;
t622 = t667 * t677 + (t674 * t338 + (-t669 + t675) * t340) * t338;
t190 = t238 * t338;
t191 = rSges(6,1) * t536 + rSges(6,2) * t534;
t310 = t340 * t351;
t495 = t338 * t330 + t310;
t274 = pkin(2) - t624;
t347 = -pkin(7) + t351;
t502 = t338 * t274 + t340 * t347;
t128 = -t495 + t502;
t225 = t340 * t274;
t488 = t347 - t351;
t129 = t338 * t488 - t225 + t287;
t535 = t329 * t338;
t537 = t328 * t338;
t578 = rSges(6,3) * t340;
t161 = rSges(6,1) * t535 - rSges(6,2) * t537 - t578;
t462 = rSges(6,1) * t534;
t427 = -rSges(6,2) * t536 + t462;
t162 = rSges(6,3) * t338 + t427;
t326 = t338 * pkin(2);
t258 = -pkin(6) * t340 + t326;
t153 = -t258 + t495;
t471 = t153 * t478 + qJD(2);
t35 = t161 * t244 + t162 * t526 + (t128 * t338 + (-t129 - t154) * t340) * qJD(3) + t471;
t320 = qJD(4) * t338;
t589 = pkin(4) * t337;
t424 = t589 + t590;
t445 = t258 + t344;
t430 = t153 + t445;
t520 = t128 + t161;
t48 = t662 - t320 + t424 * t476 + (t430 + t520) * qJD(1);
t619 = -t48 * (-qJD(1) * t190 + t526 * t628) - t35 * (-t190 * t244 - t191 * t526);
t618 = (-t338 ^ 2 - t677) * t590;
t231 = Icges(6,5) * t328 + Icges(6,6) * t329;
t435 = t658 * t348;
t436 = t656 * t348;
t608 = -qJD(1) * t231 + t328 * t436 + t329 * t435;
t552 = Icges(6,3) * t338;
t156 = Icges(6,5) * t534 - Icges(6,6) * t536 + t552;
t440 = -qJD(1) * t157 + t160 * t348 - t233 * t526;
t442 = qJD(1) * t159 + t348 * t660;
t607 = -qJD(1) * t156 + t328 * t440 + t329 * t442;
t441 = (t340 * t404 + t555) * qJD(1) + t659 * t348;
t443 = -(t236 * t340 + t561) * qJD(1) + t661 * t348;
t487 = qJD(1) * t155;
t606 = t328 * t441 + t329 * t443 - t487;
t357 = qJD(1) ^ 2;
t604 = -m(5) - m(6);
t469 = -qJDD(3) - qJDD(5);
t163 = -qJD(1) * t526 + t338 * t469;
t603 = t163 / 0.2e1;
t470 = qJD(1) * qJD(3);
t309 = t338 * t470;
t480 = qJD(1) * t338;
t164 = qJD(5) * t480 + t340 * t469 + t309;
t602 = t164 / 0.2e1;
t240 = -qJDD(3) * t338 - t340 * t470;
t600 = t240 / 0.2e1;
t241 = -qJDD(3) * t340 + t309;
t599 = t241 / 0.2e1;
t598 = -t244 / 0.2e1;
t597 = t244 / 0.2e1;
t596 = t526 / 0.2e1;
t595 = -t526 / 0.2e1;
t594 = -t338 / 0.2e1;
t593 = -t340 / 0.2e1;
t592 = rSges(4,3) + pkin(6);
t588 = -qJD(1) / 0.2e1;
t587 = qJD(1) / 0.2e1;
t586 = rSges(4,1) * t354;
t580 = rSges(4,3) * t340;
t579 = rSges(5,3) * t340;
t297 = rSges(4,1) * t352 + rSges(4,2) * t354;
t217 = t297 * t338;
t426 = rSges(4,1) * t529 - rSges(4,2) * t530;
t192 = t426 - t580;
t95 = t297 * t476 + (t192 + t445) * qJD(1);
t577 = t217 * t95;
t423 = -t297 * t478 + t334;
t305 = rSges(4,2) * t525;
t465 = rSges(4,1) * t524;
t193 = rSges(4,3) * t338 - t305 + t465;
t504 = t193 + t259;
t96 = qJD(1) * t504 + t423;
t576 = t338 * t96;
t575 = t340 * t96;
t573 = qJDD(1) / 0.2e1;
t572 = rSges(5,3) - t351;
t571 = rSges(6,3) - t347;
t479 = qJD(1) * t340;
t317 = pkin(6) * t479;
t475 = qJD(3) * t352;
t461 = pkin(3) * t475;
t112 = t340 * t461 + t317 - t320 + (t310 + (-pkin(2) + t330) * t338) * qJD(1);
t252 = t424 * qJD(3);
t570 = -t112 - (t252 - t461) * t340 - (t488 * t340 + (t274 - t330) * t338) * qJD(1);
t551 = pkin(1) * qJDD(1);
t550 = t157 * t328;
t548 = t159 * t329;
t541 = t231 * t338;
t179 = t231 * t340;
t539 = t424 * t340;
t210 = t338 * t252;
t356 = qJD(3) ^ 2;
t527 = t339 * t356;
t79 = -t233 * t536 + t534 * t629 + t541;
t523 = t79 * qJD(1);
t521 = -t254 * t476 - (t338 * t626 - t579) * qJD(1) - t112;
t519 = -t129 + t162;
t508 = t186 + t215;
t507 = t187 + t216;
t506 = t188 - t213;
t505 = -Icges(4,2) * t524 + t189 - t303;
t503 = t274 * t479 - t210;
t497 = rSges(5,3) * t480 + qJD(1) * t464;
t494 = rSges(4,3) * t480 + qJD(1) * t465;
t491 = pkin(2) * t479 + pkin(6) * t480;
t489 = t357 * t346 + t353 * t551;
t482 = qJD(1) * t247;
t481 = qJD(1) * t292;
t477 = qJD(3) * t339;
t474 = qJD(3) * t354;
t468 = pkin(3) * t530;
t467 = t356 * t345;
t92 = t662 + (t338 * t628 - t578) * qJD(1);
t466 = -t92 + t570;
t463 = t348 * t584;
t460 = pkin(3) * t474;
t387 = rSges(6,3) * t480 + qJD(1) * t462 - t338 * t463;
t457 = t329 * t244;
t93 = (-t328 * t479 - t457) * rSges(6,2) + t387;
t458 = t161 * t479 - t162 * t480 + t338 * t93;
t270 = t330 * t479;
t113 = -t473 + t270 + (-qJD(1) * t351 - t461) * t338 - t491;
t456 = t338 * t113 + t153 * t479 + t154 * t480;
t455 = -t154 + t519;
t452 = t113 * t478 + t241 * t154 + qJDD(2);
t451 = t480 / 0.2e1;
t450 = -t479 / 0.2e1;
t449 = -t478 / 0.2e1;
t448 = t478 / 0.2e1;
t447 = -t476 / 0.2e1;
t446 = t476 / 0.2e1;
t257 = rSges(3,1) * t340 - t338 * rSges(3,2);
t439 = -t156 - t550;
t438 = -t156 + t548;
t433 = t259 + t455;
t432 = qJD(1) * t491 + qJDD(1) * t258 + t489;
t431 = -t344 * t357 + t355 * t551;
t425 = rSges(3,1) * t479 - rSges(3,2) * t480;
t300 = rSges(2,1) * t355 - t353 * rSges(2,2);
t298 = rSges(2,1) * t353 + rSges(2,2) * t355;
t299 = -rSges(4,2) * t352 + t586;
t411 = t340 * t95 - t576;
t218 = t297 * t340;
t125 = qJD(3) * t218 + (t299 * t338 - t580) * qJD(1);
t126 = -rSges(4,1) * t338 * t475 + (-t338 * t474 - t352 * t479) * rSges(4,2) + t494;
t403 = -t125 * t340 + t126 * t338;
t77 = -t158 * t329 - t160 * t328;
t395 = t192 * t338 + t193 * t340;
t394 = -t233 * t328 + t329 * t629;
t389 = -t238 - t424;
t172 = rSges(5,1) * t531 - rSges(5,2) * t533 - t579;
t381 = -qJD(1) * t232 + t179 * t244 - t526 * t541;
t378 = qJD(1) * t402 - t179 * t348 - t487;
t377 = t348 * t541 + (-t232 * t340 + t548 - t550 - t552) * qJD(1);
t374 = t352 * t506 + t354 * t508;
t373 = t352 * t505 + t354 * t507;
t368 = qJD(1) * t394 - t232 * t348;
t365 = t338 * t444;
t364 = qJD(1) * t320 - qJDD(4) * t340 + t240 * t590 + t431;
t363 = qJD(1) * t113 + qJDD(1) * t153 - qJDD(4) * t338 + t340 * t467 + t432;
t11 = t377 * t338 + t340 * t606;
t12 = t378 * t338 - t340 * t607;
t13 = -t338 * t606 + t377 * t340;
t14 = t338 * t607 + t378 * t340;
t130 = t159 * t535;
t59 = -t155 * t340 - t157 * t537 + t130;
t131 = t160 * t535;
t60 = t156 * t340 + t158 * t537 - t131;
t78 = t338 * t394 - t179;
t75 = t78 * qJD(1);
t23 = -t244 * t60 - t526 * t59 + t75;
t132 = t157 * t536;
t61 = -t155 * t338 - t159 * t534 + t132;
t62 = t156 * t338 - t340 * t402;
t24 = -t244 * t62 - t526 * t61 - t523;
t38 = t368 * t338 + t340 * t608;
t39 = -t338 * t608 + t368 * t340;
t40 = -t328 * t443 + t329 * t441;
t41 = t328 * t442 - t329 * t440;
t76 = t157 * t329 + t159 * t328;
t362 = (qJD(1) * t38 - qJDD(1) * t79 - t11 * t526 - t12 * t244 + t163 * t62 + t164 * t61) * t594 + (-t328 * t605 + t329 * t359) * t588 + (qJD(1) * t39 + qJDD(1) * t78 - t13 * t526 - t14 * t244 + t163 * t60 + t164 * t59) * t593 + t23 * t451 + t24 * t450 + (-t11 * t340 - t12 * t338 + (t338 * t61 - t340 * t62) * qJD(1)) * t598 + (-t338 * t62 - t340 * t61) * t603 + (-t338 * t60 - t340 * t59) * t602 + (-t13 * t340 - t14 * t338 + (t338 * t59 - t340 * t60) * qJD(1)) * t595 + (-t338 * t77 - t340 * t76) * t573 + (-t338 * t41 - t340 * t40 + (t338 * t76 - t340 * t77) * qJD(1)) * t587 + (t381 * t338 + t340 * t666) * t597 + (-t338 * t666 + t381 * t340) * t596;
t358 = -t244 * t238 - t210 + t434;
t275 = t299 * qJD(3);
t237 = pkin(2) * t480 - t317;
t230 = t626 * qJD(3);
t229 = t338 * t424;
t222 = t257 + t346;
t206 = t628 * t348;
t201 = t254 * t338;
t177 = -t307 + t539;
t176 = -t229 + t468;
t165 = t244 * t628;
t147 = t338 * t161;
t146 = t338 * t153;
t108 = -t478 * t585 + (-t337 * t479 - t338 * t477) * rSges(5,2) + t497;
t94 = qJD(3) * t395 + qJD(2);
t74 = -t270 + (-qJD(1) * t488 + t461) * t338 + t503;
t63 = -t320 + t444 * t476 + (t172 + t430) * qJD(1);
t55 = (t172 * t338 + t340 * t517) * qJD(3) + t471;
t54 = -t275 * t478 + t240 * t297 + t504 * qJDD(1) + (-t125 - t237) * qJD(1) + t431;
t53 = qJD(1) * t126 + qJDD(1) * t192 - t241 * t297 + t275 * t476 + t432;
t49 = qJD(1) * t433 + t358;
t47 = qJD(3) * t403 - t192 * t240 - t193 * t241 + qJDD(2);
t30 = t240 * t254 + (-qJD(3) * t230 - t467) * t338 + t454 * qJDD(1) + (-t237 + t521) * qJD(1) + t364;
t29 = t230 * t476 + qJDD(1) * t172 - t444 * t241 + (t108 - t473) * qJD(1) + t363;
t15 = -t173 * t241 + (-t153 - t172) * t240 + (t108 * t338 + t340 * t521) * qJD(3) + t452;
t10 = -t338 * t467 + t163 * t238 - t206 * t244 + (t240 * t337 - t338 * t527) * pkin(4) + t433 * qJDD(1) + (-t237 + t466) * qJD(1) + t364;
t9 = -t241 * t590 - t164 * t238 + t206 * t526 + t520 * qJDD(1) + (-t241 * t337 + t340 * t527) * pkin(4) + (t74 + t93 - t473) * qJD(1) + t363;
t5 = t129 * t241 - t161 * t163 - t162 * t164 + t244 * t93 - t526 * t92 + (-t128 - t153) * t240 + (t338 * t74 + t340 * t570) * qJD(3) + t452;
t1 = [-t163 * t79 / 0.2e1 + (t75 - (t338 * t439 + t130 + t62) * t526 + (t131 - t132 + t61 + (t155 - t549) * t338) * t244 + (t244 * t438 - t639) * t340) * t597 - m(2) * (g(2) * t300 + g(3) * t298) + t77 * t603 + (t76 + t78) * t602 - t640 * t240 / 0.2e1 - t642 * t600 + (t523 - (-t132 + t60) * t526 + (-t130 + t59) * t244 + (-t244 * t388 - t438 * t526) * t340 + (-t244 * t439 + t639) * t338 + t24) * t596 + (t40 + t39) * t595 + ((((t634 + t673) * t340 + t636 - t644) * t340 + ((t634 + t687) * t338 + (t620 - t621) * t340 - t635 + t645 + t685) * t338) * qJD(3) + t663) * t448 + ((t357 * t422 - g(2) + t431) * t222 - (qJD(1) * (t334 + t425) + t489 - g(3)) * t670) * m(3) + (t41 + t38 + t23) * t598 + (m(3) * (qJD(1) * t257 + t334 - t425) * t670 + t227 * t339 + t228 * t337 + t272 * t354 + t273 * t352 - t328 * t435 + t329 * t436 + t681 * qJD(3)) * qJD(1) + (-(qJD(1) * t519 + t358 - t49 + t632) * t48 + t49 * t320 + t48 * (-rSges(6,2) * t457 + t334 + t387 + t503) + (t49 * (-t348 * t581 - t252 - t463) - t48 * qJD(4)) * t340 + (-t49 * t344 + (-t48 * t582 + t49 * t571) * t340 + (t49 * (-t274 - t628) - t48 * t347) * t338) * qJD(1) + (-g(2) + t10) * (t338 * t571 + t225 + t346 + t427) + (-g(3) + t9) * (t161 + t344 + t502)) * m(6) + (-(qJD(1) * t173 + t361 + t632 - t64) * t63 + t64 * t320 + t63 * (t270 + t434 + t497) + (-t340 * t664 - t365 * t63) * qJD(3) + (-t64 * t344 + (t572 * t64 - t583 * t63) * t340 + (t64 * (-t330 - t626) - t63 * t351) * t338) * qJD(1) + (-g(2) + t30) * (t338 * t572 + t287 + t346 + t428) + (-g(3) + t29) * (t172 + t344 + t495)) * m(5) + (t96 * t317 + (-t297 * t575 - t577) * qJD(3) + (t96 * (t580 - t344) + (-pkin(2) - t299) * t576) * qJD(1) + (-t243 + t334 - t423 + t491 + t494 + t96 + (-t193 - t305) * qJD(1)) * t95 + (-g(2) + t54) * (-t305 + t346 + (pkin(2) + t586) * t340 + t592 * t338) + (-g(3) + t53) * (-t340 * t592 + t326 + t344 + t426)) * m(4) + (t641 + t643) * t599 + (t648 - t651) * t447 + ((((t621 + t633) * t340 + t635 - t646) * t340 + ((t620 + t633) * t338 + t636 + t647) * t338) * qJD(3) + t652 + t657) * t446 + (m(2) * (t298 ^ 2 + t300 ^ 2) + Icges(2,3) + Icges(3,3) + t233 * t329 + t629 * t328 + m(3) * (t222 * t257 + t422 * t670) + t682) * qJDD(1) + (t649 - t650 + t653) * t449; m(3) * qJDD(2) + m(4) * t47 + m(5) * t15 + m(6) * t5 + (-m(3) - m(4) + t604) * g(1); t362 + t672 * t600 + t671 * t599 + (qJD(1) * t649 + qJD(3) * t622 - qJDD(1) * t640 + t240 * t644 + t241 * t645) * t594 + (qJD(1) * t648 + qJD(3) * t623 + qJDD(1) * t641 + t240 * t646 + t241 * t647) * t593 + (((t338 * t505 - t340 * t506) * t354 + (-t338 * t507 + t340 * t508) * t352 - t337 * t610 + t339 * t375) * qJD(3) + (-t337 * t499 + t339 * t498 - t352 * t493 + t354 * t492) * qJD(1)) * t588 + (t651 * t340 + t650 * t338 + (t338 * t643 + t340 * t642) * qJD(1)) * t587 + (t338 * t642 - t340 * t643) * t573 + t653 * t451 + t652 * t450 + ((t338 * t645 - t340 * t644) * qJD(1) + t622) * t449 + ((t212 * t478 - t481) * t338 + (t655 + (-t374 * t340 + (-t538 + t373) * t338) * qJD(3)) * t340 + (t196 * t478 - t482) * t338 + (t654 + (-t338 * t540 - t665) * qJD(3)) * t340) * t448 + ((t338 * t647 - t340 * t646) * qJD(1) + t623) * t447 + ((-t476 * t538 - t481) * t340 + (-t655 + (-t373 * t338 + (t212 + t374) * t340) * qJD(3)) * t338 + (-t476 * t540 - t482) * t340 + (-t654 + (t340 * t196 + t665) * qJD(3)) * t338) * t446 + (-g(1) * (t628 - t624) - g(2) * (-t229 - t190) - g(3) * (t191 + t539) + t5 * (t146 + t147) + t35 * (t456 + t458) + t9 * t307 + (t5 * t455 + t35 * t466 + t9 * (t238 + t589) + t48 * t206 + (t35 * t128 + t389 * t49) * qJD(1)) * t340 + (t5 * t128 + t35 * t74 + t10 * t389 + t49 * (-pkin(4) * t477 - t206 - t460) + (t35 * t129 + t389 * t48) * qJD(1)) * t338 + t49 * t165 - (t49 * (-t177 - t191 - t307) + t48 * (t176 - t468)) * qJD(1) - (t35 * (-t177 * t340 + t618) + (t35 * t176 + t49 * t624) * t338) * qJD(3) + t619) * m(6) + (-g(1) * t299 + g(2) * t217 - g(3) * t218 + t47 * t395 + t94 * ((t192 * t340 - t193 * t338) * qJD(1) + t403) + t411 * t275 + (-t54 * t338 + t53 * t340 + (-t338 * t95 - t575) * qJD(1)) * t297 - (-t218 * t96 - t577) * qJD(1) - (t94 * (-t217 * t338 - t218 * t340) + t411 * t299) * qJD(3)) * m(4) + (-g(1) * t631 - g(3) * t630 + g(2) * t365 + t15 * t146 + t55 * t456 + t29 * t307 + (t15 * t517 + t55 * t521 + t29 * t254 + t63 * t230 + (t55 * t172 - t664) * qJD(1)) * t340 + (t15 * t172 + t55 * t108 - t30 * t444 + t64 * (-t230 - t460) + (-t55 * t173 - t444 * t63) * qJD(1)) * t338 - (-t64 * t630 + t63 * (-t201 - t468)) * qJD(1) - (t55 * (-t202 * t340 + t618) + t63 * t626 * t340 + (-t55 * t201 - t631 * t64) * t338) * qJD(3)) * m(5); t604 * (-g(2) * t340 - g(3) * t338) + m(5) * (-t29 * t338 - t30 * t340) + m(6) * (-t10 * t340 - t338 * t9); t362 + (t5 * (t162 * t340 + t147) + t35 * (-t340 * t92 + t458) + (-t338 * t49 + t340 * t48) * t206 + (-t10 * t338 + t9 * t340 + (-t338 * t48 - t340 * t49) * qJD(1)) * t238 - t49 * (-qJD(1) * t191 - t165) - g(1) * t628 + g(2) * t190 - g(3) * t191 + t619) * m(6);];
tau = t1;
