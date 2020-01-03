% Calculate vector of inverse dynamics joint torques for
% S5RPRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:39
% EndTime: 2019-12-31 18:09:12
% DurationCPUTime: 29.49s
% Computational Cost: add. (13249->639), mult. (12362->763), div. (0->0), fcn. (9513->8), ass. (0->338)
t707 = Icges(4,3) + Icges(5,3);
t336 = qJ(3) + pkin(8);
t329 = sin(t336);
t331 = cos(t336);
t339 = sin(qJ(3));
t341 = cos(qJ(3));
t706 = Icges(4,5) * t341 + Icges(5,5) * t331 - Icges(4,6) * t339 - Icges(5,6) * t329;
t226 = Icges(6,4) * t331 + Icges(6,6) * t329;
t337 = qJ(1) + pkin(7);
t330 = sin(t337);
t332 = cos(t337);
t154 = Icges(6,2) * t330 + t226 * t332;
t680 = t707 * t330 + t706 * t332;
t705 = t154 + t680;
t704 = t707 * t332;
t533 = t330 * t341;
t534 = t330 * t339;
t535 = t330 * t331;
t537 = t329 * t330;
t673 = -Icges(4,5) * t533 - Icges(5,5) * t535 + Icges(4,6) * t534 + Icges(5,6) * t537 + t704;
t306 = Icges(6,5) * t329;
t407 = Icges(6,1) * t331 + t306;
t158 = Icges(6,4) * t330 + t332 * t407;
t562 = Icges(5,4) * t329;
t232 = Icges(5,1) * t331 - t562;
t160 = Icges(5,5) * t330 + t232 * t332;
t688 = t158 + t160;
t227 = Icges(5,2) * t331 + t562;
t700 = Icges(6,3) * t331 + t227 - t306;
t558 = Icges(6,5) * t331;
t229 = Icges(6,1) * t329 - t558;
t310 = Icges(5,4) * t331;
t699 = Icges(5,1) * t329 + t229 + t310;
t553 = Icges(5,6) * t332;
t155 = Icges(5,4) * t535 - Icges(5,2) * t537 - t553;
t554 = Icges(4,6) * t332;
t167 = Icges(4,4) * t533 - Icges(4,2) * t534 - t554;
t703 = t155 * t329 + t167 * t339;
t563 = Icges(4,4) * t339;
t278 = Icges(4,1) * t341 - t563;
t170 = Icges(4,5) * t330 + t278 * t332;
t702 = -t160 * t535 - t170 * t533;
t222 = Icges(6,3) * t329 + t558;
t405 = -Icges(5,2) * t329 + t310;
t701 = t222 - t405;
t698 = t232 + t407;
t260 = Icges(5,4) * t537;
t559 = Icges(5,5) * t332;
t159 = Icges(5,1) * t535 - t260 - t559;
t286 = Icges(4,4) * t534;
t560 = Icges(4,5) * t332;
t169 = Icges(4,1) * t533 - t286 - t560;
t697 = t159 * t331 + t169 * t341 - t703;
t683 = Icges(4,5) * t339 + Icges(4,6) * t341 + (Icges(5,6) - Icges(6,6)) * t331 + (Icges(6,4) + Icges(5,5)) * t329;
t696 = t226 + t706;
t532 = t331 * t332;
t259 = Icges(6,5) * t532;
t536 = t329 * t332;
t552 = Icges(6,6) * t330;
t150 = Icges(6,3) * t536 + t259 + t552;
t427 = -t150 * t537 + t154 * t332 - t158 * t535;
t156 = Icges(5,6) * t330 + t332 * t405;
t333 = Icges(4,4) * t341;
t406 = -Icges(4,2) * t339 + t333;
t168 = Icges(4,6) * t330 + t332 * t406;
t695 = -t680 * t332 - t702;
t659 = -t156 * t537 - t168 * t534 + t695;
t631 = -t427 + t659;
t527 = t332 * t341;
t694 = -t159 * t532 - t169 * t527 + t673 * t330;
t693 = t150 * t536 + t170 * t527 + t705 * t330 + t688 * t532;
t275 = Icges(4,2) * t341 + t563;
t277 = Icges(4,1) * t339 + t333;
t682 = t275 * t339 - t277 * t341 + t700 * t329 - t699 * t331;
t528 = t332 * t339;
t692 = t155 * t536 + t167 * t528 + t694;
t149 = -Icges(6,6) * t332 + t222 * t330;
t691 = t149 - t155;
t690 = -t150 + t156;
t157 = -Icges(6,4) * t332 + t330 * t407;
t689 = -t157 - t159;
t687 = t701 * qJD(3);
t686 = t698 * qJD(3);
t685 = t700 * qJD(3);
t684 = t699 * qJD(3);
t153 = -Icges(6,2) * t332 + t226 * t330;
t530 = t332 * t153;
t401 = t149 * t329 + t157 * t331;
t622 = t330 * t401;
t50 = -t530 + t622;
t635 = t330 * t697 + t332 * t673 + t50;
t633 = -t156 * t536 - t168 * t528 + t693;
t681 = -t275 * t341 - t277 * t339 - t329 * t699 - t331 * t700;
t625 = t683 * t332;
t624 = t683 * t330;
t338 = -qJ(4) - pkin(6);
t291 = t332 * t338;
t334 = t341 * pkin(3);
t325 = t334 + pkin(2);
t493 = -t330 * t325 - t291;
t340 = sin(qJ(1));
t589 = pkin(1) * t340;
t679 = t493 - t589;
t678 = t156 * t329 + t168 * t339;
t647 = -t682 * t330 - t625;
t646 = -t682 * t332 + t624;
t677 = t685 * t332 + (t330 * t405 - t149 - t553) * qJD(1);
t676 = t685 * t330 + (t222 * t332 - t156 + t552) * qJD(1);
t675 = -t684 * t332 + (-t232 * t330 - t157 + t559) * qJD(1);
t674 = -t688 * qJD(1) + t684 * t330;
t465 = rSges(5,1) * t535;
t671 = -t465 + t679;
t632 = -t167 * t341 - t169 * t339 + t689 * t329 + t691 * t331;
t630 = t168 * t341 + t170 * t339 + t688 * t329 + t690 * t331;
t248 = t406 * qJD(3);
t249 = t278 * qJD(3);
t670 = t683 * qJD(1) + t681 * qJD(3) - t248 * t339 + t249 * t341 + t687 * t329 + t686 * t331;
t669 = t683 * qJD(3);
t668 = t150 * t329 + t170 * t341 + t688 * t331 - t678;
t667 = -t401 - t697;
t138 = t330 * t153;
t54 = t149 * t536 + t157 * t532 + t138;
t639 = t332 * t54;
t666 = t633 * t330 + t692 * t332 - t639;
t665 = t631 * t330 - t635 * t332;
t664 = t682 * qJD(1) + t696 * qJD(3);
t342 = cos(qJ(1));
t335 = t342 * pkin(1);
t236 = rSges(3,1) * t330 + rSges(3,2) * t332;
t206 = -t236 - t589;
t663 = t646 * qJD(1);
t662 = t647 * qJD(1);
t661 = t705 * qJD(1);
t660 = t332 ^ 2;
t344 = qJD(1) ^ 2;
t467 = t344 * t335;
t658 = -t691 * t332 + (-Icges(6,1) * t536 + t229 * t332 + t259 - t690) * t330;
t657 = t673 + t678;
t656 = t699 - t701;
t655 = t698 - t700;
t654 = (Icges(5,2) * t535 + t260 + t689) * t332 + (-t227 * t332 + t688) * t330;
t478 = qJD(3) * t330;
t653 = t665 * qJD(3) + t662;
t652 = t666 * qJD(3) + t663;
t375 = qJD(3) * t275;
t110 = qJD(1) * t168 - t330 * t375;
t378 = qJD(3) * t277;
t112 = qJD(1) * t170 - t330 * t378;
t651 = t667 * qJD(3) - t110 * t341 - t112 * t339 + t674 * t329 + t676 * t331;
t109 = -t332 * t375 + (-t330 * t406 + t554) * qJD(1);
t111 = -t332 * t378 + (-t278 * t330 + t560) * qJD(1);
t650 = t668 * qJD(3) + t109 * t341 + t111 * t339 + t675 * t329 - t677 * t331;
t649 = t664 * t330 + t670 * t332;
t648 = t670 * t330 - t664 * t332;
t634 = t54 - t692;
t645 = -t630 * qJD(3) - t109 * t339 + t111 * t341 + t677 * t329 + t675 * t331 + t661;
t610 = qJD(1) * t153;
t644 = t673 * qJD(1) - t632 * qJD(3) + t110 * t339 - t112 * t341 - t676 * t329 + t674 * t331 - t610;
t643 = t530 + t693;
t642 = t667 * qJD(1) - t669 * t330 + t661;
t641 = -t610 - t669 * t332 + (-t330 * t706 - t668 + t704) * qJD(1);
t640 = rSges(4,2) * t339;
t638 = rSges(6,3) + qJ(5);
t590 = rSges(6,1) + pkin(4);
t637 = -t654 * t329 + t331 * t658;
t615 = t331 * rSges(6,1) + t329 * rSges(6,3);
t636 = t331 * pkin(4) + t329 * qJ(5) + t615;
t489 = t277 + t406;
t490 = -t275 + t278;
t629 = (-t329 * t656 + t331 * t655 - t339 * t489 + t341 * t490) * qJD(1);
t628 = t642 * t660 + (t645 * t330 + (-t641 + t644) * t332) * t330;
t627 = t644 * t660 + (t641 * t330 + (-t642 + t645) * t332) * t330;
t626 = t696 * qJD(1);
t623 = t329 * t590;
t323 = t332 * pkin(6);
t241 = pkin(2) * t330 - t323;
t145 = t241 + t493;
t220 = qJD(1) * t241;
t621 = qJD(1) * t145 - t220;
t313 = t330 * rSges(4,3);
t172 = rSges(4,1) * t527 - rSges(4,2) * t528 + t313;
t321 = t330 * pkin(6);
t242 = t332 * pkin(2) + t321;
t445 = t242 + t335;
t119 = t172 + t445;
t614 = t331 * rSges(5,1) - rSges(5,2) * t329;
t620 = t614 + t334;
t240 = t332 * rSges(3,1) - rSges(3,2) * t330;
t207 = t240 + t335;
t477 = qJD(3) * t332;
t459 = t331 * t477;
t479 = qJD(1) * t332;
t619 = rSges(6,2) * t479 + t459 * t638;
t618 = t638 * t535;
t617 = t638 * t532;
t616 = -rSges(5,2) * t537 - t332 * rSges(5,3);
t234 = rSges(6,1) * t329 - rSges(6,3) * t331;
t497 = pkin(4) * t329 - qJ(5) * t331 + t234;
t588 = pkin(3) * t339;
t432 = -t497 - t588;
t475 = qJD(5) * t329;
t253 = t332 * t475;
t297 = qJD(4) * t330;
t494 = t253 + t297;
t612 = t432 * t477 + t494;
t586 = g(2) * t330;
t611 = -g(1) * t332 - t586;
t608 = t334 + t636;
t506 = -Icges(4,2) * t533 + t169 - t286;
t508 = t277 * t330 + t167;
t596 = -t339 * t506 - t341 * t508;
t595 = -m(5) - m(6);
t472 = qJD(1) * qJD(3);
t217 = qJDD(3) * t330 + t332 * t472;
t594 = t217 / 0.2e1;
t218 = -qJDD(3) * t332 + t330 * t472;
t593 = t218 / 0.2e1;
t592 = t330 / 0.2e1;
t591 = -t332 / 0.2e1;
t583 = pkin(2) - t325;
t480 = qJD(1) * t330;
t367 = -t329 * t477 - t331 * t480;
t461 = t329 * t480;
t582 = t367 * t590 - t461 * t638 + t253 + t619;
t314 = t330 * rSges(6,2);
t581 = (pkin(4) * t479 + qJ(5) * t478) * t331 + (qJ(5) * t479 + (-pkin(4) * qJD(3) + qJD(5)) * t330) * t329 - t234 * t478 + (t332 * t615 + t314) * qJD(1);
t580 = rSges(4,1) * t341;
t280 = rSges(4,1) * t339 + rSges(4,2) * t341;
t204 = t280 * t332;
t76 = qJD(1) * t119 - t280 * t478;
t577 = t204 * t76;
t439 = t332 * t325 - t330 * t338;
t146 = t439 - t242;
t454 = -t145 * t478 + t146 * t477 + qJD(2);
t474 = qJD(5) * t331;
t509 = t532 * t590 + t536 * t638 + t314;
t319 = t332 * rSges(6,2);
t510 = t330 * t636 - t319;
t25 = -t474 + (t330 * t510 + t332 * t509) * qJD(3) + t454;
t576 = t25 * t329;
t312 = t330 * rSges(5,3);
t488 = rSges(4,2) * t534 + t332 * rSges(4,3);
t171 = rSges(4,1) * t533 - t488;
t446 = -t241 - t589;
t433 = -t171 + t446;
t457 = t280 * t477;
t75 = qJD(1) * t433 - t457;
t574 = t330 * t75;
t572 = t332 * t75;
t235 = rSges(5,1) * t329 + rSges(5,2) * t331;
t368 = -t235 - t588;
t366 = t368 * t477 + t297;
t162 = t465 + t616;
t387 = t145 - t162 + t446;
t48 = qJD(1) * t387 + t366;
t570 = t48 * t235;
t526 = t341 * qJD(3) ^ 2;
t470 = pkin(3) * t534;
t491 = qJD(3) * t470 + qJD(4) * t332;
t462 = t338 * t480 + t491;
t103 = (-t332 * t583 - t321) * qJD(1) - t462;
t216 = t242 * qJD(1);
t522 = -t103 - t216;
t518 = -t330 * t145 + t332 * t146;
t164 = rSges(5,1) * t532 - rSges(5,2) * t536 + t312;
t515 = -t146 - t164;
t507 = -t277 * t332 - t168;
t505 = -t275 * t332 + t170;
t504 = -qJD(3) * t636 + t474;
t503 = -t537 * t590 + t618;
t502 = -t536 * t590 + t617;
t495 = rSges(5,2) * t461 + rSges(5,3) * t479;
t492 = rSges(4,3) * t479 + t480 * t640;
t487 = t330 ^ 2 + t660;
t476 = qJD(3) * t341;
t471 = qJD(1) * qJD(4);
t469 = pkin(3) * t528;
t468 = pkin(3) * t526;
t295 = pkin(6) * t479;
t456 = t339 * t477;
t102 = -pkin(3) * t456 - t295 + t297 + (t330 * t583 - t291) * qJD(1);
t466 = t332 * t102 + t330 * t103 - t145 * t479;
t464 = pkin(3) * t476;
t463 = -t146 - t509;
t455 = t330 * t475;
t453 = -pkin(2) - t580;
t450 = -t478 / 0.2e1;
t449 = t478 / 0.2e1;
t448 = -t477 / 0.2e1;
t447 = t477 / 0.2e1;
t438 = t102 * t477 + t103 * t478 - t217 * t145 + qJDD(2);
t437 = t487 * t588;
t436 = qJDD(1) * t335 - t344 * t589;
t435 = t146 + t445;
t434 = -t510 - t589;
t215 = t614 * qJD(3);
t430 = -t215 - t464;
t426 = t335 + t439;
t283 = rSges(2,1) * t342 - rSges(2,2) * t340;
t281 = rSges(2,1) * t340 + rSges(2,2) * t342;
t282 = t580 - t640;
t381 = t145 - t241 + t434;
t37 = qJD(1) * t381 + t612;
t38 = (-qJD(3) * t497 + t475) * t330 + (t435 + t509) * qJD(1) - t491;
t417 = t330 * t38 + t332 * t37;
t410 = -t330 * t76 - t572;
t115 = -rSges(4,2) * t332 * t476 + (-t341 * t480 - t456) * rSges(4,1) + t492;
t203 = t280 * t330;
t116 = -qJD(3) * t203 + (t282 * t332 + t313) * qJD(1);
t403 = t115 * t332 + t116 * t330;
t394 = t171 * t330 + t172 * t332;
t386 = -t464 + t504;
t383 = qJDD(4) * t330 + t218 * t588 + t332 * t471 - t467;
t382 = qJD(1) * (-pkin(2) * t480 + t295) + qJDD(1) * t242 + t436;
t188 = t235 * t330;
t380 = t417 * t331;
t363 = -t339 * t505 + t341 * t507;
t362 = t368 * t332;
t358 = -t325 - t636;
t357 = qJD(1) * t102 + qJDD(1) * t146 - qJDD(4) * t332 + t330 * t471 + t382;
t356 = -t468 + qJDD(5) * t329 + (t474 + t504) * qJD(3);
t250 = t282 * qJD(3);
t192 = t235 * t332;
t98 = -qJD(3) * t188 + (t332 * t614 + t312) * qJD(1);
t96 = rSges(5,1) * t367 - rSges(5,2) * t459 + t495;
t74 = qJD(3) * t394 + qJD(2);
t49 = -t235 * t478 + (t164 + t435) * qJD(1) - t491;
t43 = (t162 * t330 + t164 * t332) * qJD(3) + t454;
t42 = qJD(1) * t115 + qJDD(1) * t172 - t217 * t280 - t250 * t478 + t382;
t41 = -t467 - t250 * t477 + t218 * t280 + (-t116 - t216) * qJD(1) + t433 * qJDD(1);
t34 = qJD(3) * t403 + t171 * t217 - t172 * t218 + qJDD(2);
t18 = -t215 * t478 + qJD(1) * t96 + qJDD(1) * t164 - t217 * t235 + (-t217 * t339 - t330 * t526) * pkin(3) + t357;
t17 = t218 * t235 + (-qJD(3) * t215 - t468) * t332 + (-t98 + t522) * qJD(1) + t387 * qJDD(1) + t383;
t4 = t162 * t217 + t515 * t218 + (t330 * t98 + t332 * t96) * qJD(3) + t438;
t3 = t497 * t218 + t356 * t332 + t381 * qJDD(1) + (-t455 + t522 - t581) * qJD(1) + t383;
t2 = t509 * qJDD(1) + t432 * t217 + (t253 + t582) * qJD(1) + t356 * t330 + t357;
t1 = -qJDD(5) * t331 + t510 * t217 + t463 * t218 + (t330 * t581 + t332 * t582 + t475) * qJD(3) + t438;
t5 = [-m(2) * (-g(1) * t281 + g(2) * t283) + ((-t236 * t344 - g(2) + t436) * t207 + (-t467 + (-0.2e1 * t240 - t335 + t207) * t344 - g(1)) * t206) * m(3) + ((-t639 + ((t680 + t703) * t332 + t659 + t694 + t702) * t332 + (t50 - t622 + t643) * t330) * qJD(3) + t663) * t447 + (-t682 * qJD(3) + t248 * t341 + t249 * t339 + t686 * t329 - t687 * t331) * qJD(1) + (t37 * (-t455 + t462) + t38 * (t494 + t619) + ((-t588 - t623) * t332 * t38 + (-t331 * t638 + t623) * t330 * t37) * qJD(3) + ((-t340 * t38 - t342 * t37) * pkin(1) + (-t38 * t338 + t358 * t37) * t332 + (-t37 * rSges(6,2) + t358 * t38) * t330) * qJD(1) - (t434 * qJD(1) - t37 + t612 + t621) * t38 + (t2 - g(2)) * (t426 + t509) + (t3 - g(1)) * (t319 + (-t329 * t638 - t331 * t590) * t330 + t679)) * m(6) + (t570 * t478 + (t462 + (-t312 - t335 + (-t325 - t614) * t332) * qJD(1)) * t48 + (t362 * qJD(3) + t297 - t366 + t48 + t495 - t621 + (t162 + t589 + t671) * qJD(1)) * t49 + (t18 - g(2)) * (t164 + t426) + (t17 - g(1)) * (-t616 + t671)) * m(5) + (t76 * (t295 + t492) + (t280 * t574 - t577) * qJD(3) + ((-t340 * t76 - t342 * t75) * pkin(1) + (-pkin(2) - t282) * t572 + (t75 * (-rSges(4,3) - pkin(6)) + t76 * t453) * t330) * qJD(1) - (-t457 - t220 - t75 + (-t171 - t589) * qJD(1)) * t76 + (t42 - g(2)) * t119 + (t41 - g(1)) * (t453 * t330 + t323 + t488 - t589)) * m(4) + (m(3) * (t206 ^ 2 + t240 * t207) + m(2) * (t281 ^ 2 + t283 ^ 2) + Icges(2,3) + Icges(3,3) - t681) * qJDD(1) + (t630 + t646) * t594 + (-t632 + t647) * t593 + (((t657 * t332 + t633 - t643) * t332 + (t330 * t657 - t138 + t427 + t634 - t695) * t330) * qJD(3) + t653 - t662) * t450 + (t649 + t650) * t449 + (t648 - t651 + t652) * t448; m(3) * qJDD(2) + m(4) * t34 + m(5) * t4 + m(6) * t1 + (-m(3) - m(4) + t595) * g(3); t666 * t594 + t665 * t593 + (qJD(1) * t649 + t627 * qJD(3) + qJDD(1) * t646 + t633 * t217 + t634 * t218) * t592 + (qJD(1) * t648 + t628 * qJD(3) + qJDD(1) * t647 + t631 * t217 + t635 * t218) * t591 - (((t330 * t505 - t332 * t506) * t341 + (t330 * t507 + t332 * t508) * t339 + t654 * t331 + t658 * t329) * qJD(3) + (t329 * t655 + t331 * t656 + t339 * t490 + t341 * t489) * qJD(1)) * qJD(1) / 0.2e1 + (t651 * t332 + t650 * t330 + (-t330 * t632 + t332 * t630) * qJD(1)) * qJD(1) / 0.2e1 + (t630 * t330 + t632 * t332) * qJDD(1) / 0.2e1 + t653 * t480 / 0.2e1 + t652 * t479 / 0.2e1 + ((-t625 * t478 + t626) * t330 + ((-t596 * t332 + (t363 + t624) * t330 + t637) * qJD(3) + t629) * t332) * t450 + ((t330 * t634 + t332 * t633) * qJD(1) + t627) * t449 + ((t330 * t635 + t332 * t631) * qJD(1) + t628) * t448 + ((-t624 * t477 - t626) * t332 + ((t363 * t330 + (-t596 + t625) * t332 + t637) * qJD(3) + t629) * t330) * t447 + (g(1) * t204 + g(2) * t203 - g(3) * t282 - (t203 * t75 - t577) * qJD(1) - (t74 * (-t203 * t330 - t204 * t332) + t410 * t282) * qJD(3) + t34 * t394 + t74 * ((t171 * t332 - t172 * t330) * qJD(1) + t403) + t410 * t250 + (-t42 * t330 - t41 * t332 + (-t332 * t76 + t574) * qJD(1)) * t280) * m(4) + (-g(1) * (-t469 + t617) - g(2) * (-t470 + t618) - g(3) * t608 - t611 * t623 - (t380 + t576) * qJD(5) - (-t37 * t503 + t38 * (-t469 + t502)) * qJD(1) - (-t25 * t437 + (t25 * t502 - t37 * t608) * t332 + (t25 * t503 - t38 * t608) * t330) * qJD(3) + t1 * t518 + t25 * t466 + (t3 * t432 + t37 * t386 + t1 * t509 + t25 * t582 + (t25 * t510 + t38 * t432) * qJD(1)) * t332 + (t2 * t432 + t38 * t386 + t1 * t510 + t25 * t581 + (t25 * t463 + t37 * t497) * qJD(1)) * t330) * m(6) + (-g(1) * t362 - g(3) * t620 - t368 * t586 - (t48 * t188 + t49 * (-t192 - t469)) * qJD(1) - (-t43 * t437 + (-t43 * t192 - t48 * t620) * t332 + (-t43 * t188 - t49 * t620) * t330) * qJD(3) + t4 * t518 + t43 * t466 + (t17 * t368 + t48 * t430 + t4 * t164 + t43 * t96 + (t43 * t162 + t368 * t49) * qJD(1)) * t332 + (t18 * t368 + t49 * t430 + t4 * t162 + t43 * t98 + (t43 * t515 + t570) * qJD(1)) * t330) * m(5); t595 * (g(1) * t330 - g(2) * t332) + 0.2e1 * (t2 * t591 + t3 * t592) * m(6) + 0.2e1 * (t17 * t592 + t18 * t591) * m(5); (-(t487 * t576 + t380) * qJD(3) + (qJD(3) * t417 + g(3) - t1) * t331 + (qJD(3) * t25 + t2 * t330 + t3 * t332 + t611) * t329) * m(6);];
tau = t5;
