% Calculate vector of inverse dynamics joint torques for
% S4RRPP5
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPP5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP5_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP5_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP5_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP5_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:10
% EndTime: 2019-12-31 17:00:36
% DurationCPUTime: 23.05s
% Computational Cost: add. (4804->611), mult. (12366->711), div. (0->0), fcn. (9666->4), ass. (0->332)
t714 = Icges(4,4) - Icges(3,5);
t713 = Icges(4,5) - Icges(3,6);
t712 = Icges(4,1) + Icges(3,3);
t339 = sin(qJ(2));
t341 = cos(qJ(2));
t711 = t713 * t339 - t714 * t341;
t565 = Icges(3,4) * t339;
t244 = Icges(3,2) * t341 + t565;
t321 = Icges(5,6) * t339;
t235 = -Icges(5,2) * t341 + t321;
t554 = Icges(4,6) * t339;
t402 = Icges(4,3) * t341 + t554;
t705 = t235 - t402;
t710 = -t244 + t705;
t240 = Icges(5,4) * t339 + Icges(5,5) * t341;
t340 = sin(qJ(1));
t342 = cos(qJ(1));
t148 = Icges(5,1) * t340 + t240 * t342;
t681 = t712 * t340 + t711 * t342;
t709 = t148 + t681;
t708 = t712 * t342;
t538 = t339 * t340;
t296 = Icges(3,4) * t538;
t536 = t340 * t341;
t560 = Icges(3,5) * t342;
t138 = Icges(3,1) * t536 - t296 - t560;
t285 = Icges(5,6) * t538;
t558 = Icges(5,5) * t342;
t141 = -Icges(5,3) * t536 - t285 + t558;
t707 = -t138 + t141;
t247 = Icges(3,1) * t341 - t565;
t139 = Icges(3,5) * t340 + t247 * t342;
t400 = Icges(5,3) * t341 + t321;
t140 = Icges(5,5) * t340 + t342 * t400;
t697 = t139 + t140;
t553 = Icges(4,6) * t341;
t403 = -Icges(4,3) * t339 + t553;
t142 = Icges(4,5) * t340 - t342 * t403;
t535 = t341 * t342;
t287 = Icges(5,6) * t535;
t537 = t339 * t342;
t562 = Icges(5,4) * t340;
t144 = Icges(5,2) * t537 + t287 + t562;
t696 = t142 + t144;
t559 = Icges(4,5) * t342;
t143 = Icges(4,6) * t536 - Icges(4,3) * t538 + t559;
t286 = Icges(5,6) * t536;
t561 = Icges(5,4) * t342;
t145 = -Icges(5,2) * t538 - t286 + t561;
t706 = t143 + t145;
t674 = t714 * t536 - t713 * t538 + t708;
t324 = Icges(3,4) * t341;
t246 = Icges(3,1) * t339 + t324;
t552 = Icges(5,6) * t341;
t695 = t552 - t553 + (-Icges(4,2) - Icges(5,3)) * t339;
t704 = -t246 + t695;
t555 = Icges(3,6) * t342;
t136 = Icges(3,4) * t536 - Icges(3,2) * t538 - t555;
t288 = Icges(4,6) * t538;
t563 = Icges(4,4) * t342;
t147 = Icges(4,2) * t536 - t288 + t563;
t703 = t136 * t339 - t147 * t341;
t702 = -t139 * t536 - t142 * t538;
t700 = t138 * t341 - t143 * t339 - t703;
t699 = (-Icges(5,4) - Icges(4,5)) * t341 + (-Icges(4,4) + Icges(5,5)) * t339;
t234 = Icges(5,2) * t339 + t552;
t410 = -Icges(3,2) * t339 + t324;
t687 = t234 - t403 - t410;
t238 = Icges(3,5) * t339 + Icges(3,6) * t341;
t686 = t238 + t699;
t406 = Icges(4,2) * t341 - t554;
t698 = t247 + t400 + t406;
t432 = -t140 * t536 - t144 * t538 + t148 * t342;
t137 = Icges(3,6) * t340 + t342 * t410;
t289 = Icges(4,6) * t537;
t564 = Icges(4,4) * t340;
t146 = -Icges(4,2) * t535 + t289 + t564;
t694 = -t342 * t681 - t702;
t646 = -t137 * t538 - t146 * t536 + t694;
t614 = -t432 + t646;
t693 = t709 * t340 + t697 * t535 + t696 * t537;
t566 = Icges(5,1) * t342;
t149 = -Icges(5,4) * t538 - Icges(5,5) * t536 + t566;
t130 = t340 * t149;
t692 = t340 * t674 + t535 * t707 + t537 * t706 + t130;
t691 = t136 + t706;
t690 = t137 - t696;
t689 = t147 - t707;
t688 = -t146 + t697;
t684 = t710 * qJD(2);
t683 = t704 * qJD(2);
t384 = t244 * t339 - t246 * t341;
t676 = t339 * t705 - t341 * t695 - t384;
t679 = t699 * t342;
t645 = t695 * t536 - t538 * t705 + t679;
t539 = t238 * t342;
t73 = -t340 * t384 - t539;
t682 = t73 - t645;
t680 = t137 * t339 + t146 * t341;
t542 = t149 * t342;
t395 = t141 * t341 + t145 * t339;
t626 = t340 * t395;
t54 = t542 - t626;
t615 = t340 * t700 + t342 * t674 + t54;
t613 = -t136 * t537 + t147 * t535 - t692;
t612 = -t137 * t537 - t146 * t535 + t693;
t678 = t698 * qJD(2);
t677 = t687 * qJD(2);
t668 = t240 + t711;
t675 = t704 * t339 + t710 * t341;
t609 = t686 * t340;
t673 = -t235 + t244;
t629 = t342 * t676 + t609;
t672 = t684 * t342 + (t340 * t687 + t555 - t559 - t561) * qJD(1);
t671 = t684 * t340 + (-t234 * t342 + t137 - t142 - t562) * qJD(1);
t670 = t683 * t342 + (-t340 * t698 + t558 + t560 - t563) * qJD(1);
t669 = -t683 * t340 + (-t342 * t406 + t564 - t697) * qJD(1);
t611 = t339 * t689 + t341 * t691;
t610 = t339 * t688 + t341 * t690;
t667 = t686 * qJD(2);
t666 = t395 - t700;
t665 = t696 * t339 + t341 * t697 - t680;
t584 = rSges(5,1) + pkin(3);
t664 = qJD(1) * t686 + t675 * qJD(2) + t677 * t339 + t678 * t341;
t663 = t612 * t340 - t613 * t342;
t662 = t340 * t614 - t615 * t342;
t661 = t709 * qJD(1);
t660 = t676 * qJD(1) - t668 * qJD(2);
t659 = (Icges(4,3) * t535 + t673 * t342 + t289 - t688) * t340 + (t285 - t288 - t296 + (-Icges(3,2) - Icges(5,2) - Icges(4,3)) * t536 + t689) * t342;
t658 = t342 ^ 2;
t644 = rSges(5,3) + qJ(4);
t657 = t629 * qJD(1);
t656 = t669 * t341 + t671 * t339 + t611 * qJD(2) + (t149 + t674) * qJD(1);
t655 = -t610 * qJD(2) - t672 * t339 + t670 * t341 + t661;
t654 = (Icges(5,3) * t538 - t286 + t691) * t342 + (-Icges(5,3) * t537 + t287 - t690) * t340;
t653 = -t687 - t704;
t652 = -t402 - t673 + t698;
t649 = t682 * qJD(1);
t648 = -t667 * t342 + (-t668 * t340 + t566 - t665 + t708) * qJD(1);
t647 = t666 * qJD(1) - t667 * t340 + t661;
t515 = -rSges(5,2) * t538 + t342 * t584 - t644 * t536;
t639 = qJ(4) * t535 + t584 * t340;
t638 = t674 + t680;
t475 = qJD(2) * qJD(3);
t637 = qJDD(3) * t339 + t341 * t475;
t482 = qJD(2) * t342;
t462 = t339 * t482;
t486 = qJD(1) * t340;
t636 = t341 * t486 + t462;
t635 = t663 * qJD(2) + t657;
t634 = t662 * qJD(2) + t649;
t633 = t665 * qJD(2) + t670 * t339 + t672 * t341;
t632 = t666 * qJD(2) + t669 * t339 - t671 * t341;
t631 = -t660 * t340 + t664 * t342;
t630 = t664 * t340 + t660 * t342;
t628 = -t542 + t693;
t473 = -pkin(2) - t644;
t627 = t473 * t341 - pkin(1);
t329 = t340 * rSges(4,1);
t162 = -rSges(4,2) * t535 + rSges(4,3) * t537 + t329;
t311 = pkin(2) * t535;
t204 = qJ(3) * t537 + t311;
t259 = t342 * pkin(1) + t340 * pkin(5);
t441 = t204 + t259;
t85 = t162 + t441;
t319 = t339 * qJ(3);
t621 = t341 * pkin(2) + t319;
t199 = t621 * t340;
t335 = t342 * pkin(5);
t258 = pkin(1) * t340 - t335;
t229 = qJD(1) * t258;
t625 = -qJD(1) * t199 - t229;
t279 = qJ(3) * t536;
t196 = -pkin(2) * t538 + t279;
t197 = rSges(4,2) * t538 + rSges(4,3) * t536;
t624 = t196 + t197;
t325 = t339 * rSges(4,3);
t427 = -rSges(4,2) * t341 + t325;
t326 = t339 * rSges(5,2);
t623 = t341 * rSges(5,3) + t326;
t477 = qJD(4) * t342;
t274 = t341 * t477;
t480 = qJD(3) * t342;
t276 = t339 * t480;
t248 = pkin(2) * t339 - qJ(3) * t341;
t425 = rSges(5,2) * t341 - rSges(5,3) * t339;
t549 = qJ(4) * t339;
t449 = -t425 + t549;
t434 = -t248 - t449;
t620 = t434 * t482 + t274 + t276;
t461 = t341 * t482;
t485 = qJD(1) * t342;
t618 = rSges(5,2) * t461 + t584 * t485 + t274;
t617 = qJD(1) * t515;
t493 = t340 ^ 2 + t658;
t616 = qJD(2) * t493;
t608 = t539 + t679;
t607 = qJ(4) * t341 + t621 + t623;
t606 = t659 * t339 + t654 * t341;
t605 = (-t653 * t339 + t652 * t341) * qJD(1);
t604 = t656 * t658 + (t648 * t340 + (-t647 + t655) * t342) * t340;
t603 = t647 * t658 + (t655 * t340 + (-t648 + t656) * t342) * t340;
t602 = t668 * qJD(1);
t479 = qJD(4) * t339;
t483 = qJD(2) * t341;
t481 = qJD(3) * t341;
t168 = qJD(2) * t621 - t481;
t513 = -qJD(2) * t623 - t168;
t592 = -qJ(4) * t483 + qJD(2) * t607 - t479 + t513;
t591 = m(4) / 0.2e1;
t590 = m(5) / 0.2e1;
t476 = qJD(1) * qJD(2);
t226 = qJDD(2) * t340 + t342 * t476;
t589 = t226 / 0.2e1;
t227 = -qJDD(2) * t342 + t340 * t476;
t587 = t227 / 0.2e1;
t583 = g(2) * t340;
t580 = rSges(3,1) * t341;
t579 = rSges(4,2) * t339;
t251 = rSges(3,1) * t339 + rSges(3,2) * t341;
t203 = t251 * t342;
t327 = t340 * rSges(3,3);
t160 = rSges(3,1) * t535 - rSges(3,2) * t537 + t327;
t113 = t160 + t259;
t484 = qJD(2) * t340;
t64 = qJD(1) * t113 - t251 * t484;
t577 = t203 * t64;
t464 = t251 * t482;
t494 = rSges(3,2) * t538 + t342 * rSges(3,3);
t159 = rSges(3,1) * t536 - t494;
t518 = -t159 - t258;
t63 = qJD(1) * t518 - t464;
t575 = t340 * t63;
t573 = t342 * t63;
t571 = -rSges(5,2) - qJ(3);
t570 = -rSges(4,3) - qJ(3);
t225 = t259 * qJD(1);
t176 = t339 * t485 + t340 * t483;
t463 = t339 * t484;
t271 = pkin(2) * t463;
t318 = qJD(3) * t339;
t459 = t340 * t318;
t431 = -t271 + t459;
t83 = qJ(3) * t176 + qJD(1) * t311 + t431;
t569 = -t225 - t83;
t263 = qJ(4) * t463;
t478 = qJD(4) * t341;
t430 = t340 * t478 - t263;
t531 = t425 * t484 + t430 + (t342 * t623 + t639) * qJD(1);
t466 = t339 * t486;
t530 = -rSges(5,2) * t466 - t636 * t644 + t618;
t517 = rSges(5,2) * t537 + rSges(5,3) * t535 + t639;
t516 = -t162 - t204;
t514 = t340 * t199 + t342 * t204;
t512 = -t427 * qJD(2) - t168;
t282 = qJ(3) * t535;
t201 = -pkin(2) * t537 + t282;
t511 = qJD(1) * t201 + t340 * t481;
t510 = -t199 - t258;
t316 = pkin(5) * t485;
t509 = qJD(1) * (-pkin(1) * t486 + t316) + qJDD(1) * t259;
t426 = rSges(4,3) * t341 + t579;
t502 = -t248 + t426;
t501 = -t621 - t427;
t500 = qJ(3) * t461 + t276;
t498 = rSges(3,2) * t466 + rSges(3,3) * t485;
t202 = rSges(4,2) * t537 + rSges(4,3) * t535;
t82 = -pkin(2) * t636 - qJ(3) * t466 + t500;
t472 = t199 * t485 + t340 * t83 + t342 * t82;
t471 = t227 * t248 + t342 * t637;
t470 = t196 * t484 + t201 * t482 + t318;
t469 = -t204 - t517;
t450 = rSges(4,1) * t342 - rSges(4,3) * t538;
t164 = rSges(4,2) * t536 + t450;
t468 = t164 + t510;
t467 = t316 + t500;
t458 = -pkin(1) - t580;
t454 = -t484 / 0.2e1;
t453 = t484 / 0.2e1;
t452 = -t482 / 0.2e1;
t451 = t482 / 0.2e1;
t444 = qJD(2) * t512;
t443 = t510 + t515;
t442 = rSges(4,1) * t485 + rSges(4,2) * t636 + rSges(4,3) * t461;
t440 = t342 * t473;
t438 = (rSges(4,2) - pkin(2)) * t341 - pkin(1);
t437 = g(1) * t342 + t583;
t428 = t199 * t484 + t204 * t482 - t481;
t257 = rSges(2,1) * t342 - rSges(2,2) * t340;
t252 = rSges(2,1) * t340 + rSges(2,2) * t342;
t256 = -rSges(3,2) * t339 + t580;
t412 = -t340 * t64 - t573;
t104 = -rSges(3,1) * t636 - rSges(3,2) * t461 + t498;
t198 = t251 * t340;
t105 = -qJD(2) * t198 + (t256 * t342 + t327) * qJD(1);
t399 = t104 * t342 + t105 * t340;
t390 = t159 * t340 + t160 * t342;
t383 = -t478 - t318;
t382 = -pkin(1) - t621;
t380 = t482 * t502 + t276;
t378 = -qJDD(3) * t341 + t226 * t199 + t339 * t475 + t82 * t482 + t83 * t484;
t368 = qJDD(1) * t204 + t509 + t637 * t340 + (t276 + t82) * qJD(1);
t357 = (-qJ(4) * qJD(2) ^ 2 + qJDD(4)) * t341 + (-0.2e1 * t479 + t513) * qJD(2);
t2 = t449 * t227 + t443 * qJDD(1) + t357 * t342 + (t340 * t383 - t531 + t569) * qJD(1) + t471;
t3 = t517 * qJDD(1) + t434 * t226 + (t274 + t530) * qJD(1) + t357 * t340 + t368;
t37 = t479 + (-t340 * t515 + t342 * t517) * qJD(2) + t428;
t367 = qJD(2) * t37 + t2 * t342 + t3 * t340;
t309 = rSges(5,2) * t535;
t303 = rSges(5,2) * t536;
t277 = t341 * t480;
t220 = t256 * qJD(2);
t207 = t248 * t486;
t206 = t248 * t484;
t205 = -rSges(5,3) * t537 + t309;
t200 = -rSges(5,3) * t538 + t303;
t177 = t461 - t466;
t175 = t339 * t616;
t109 = -rSges(4,3) * t466 + t442;
t107 = t426 * t484 + (t342 * t427 + t329) * qJD(1);
t59 = t390 * qJD(2);
t44 = -t206 + (qJD(2) * t426 + t318) * t340 + t85 * qJD(1);
t43 = qJD(1) * t468 + t380;
t42 = (t162 * t342 - t164 * t340) * qJD(2) + t428;
t41 = -t206 - t263 + (qJD(2) * t425 - t383) * t340 + (t259 - t469) * qJD(1);
t40 = qJD(1) * t443 + t620;
t39 = qJD(1) * t104 + qJDD(1) * t160 - t220 * t484 - t226 * t251 + t509;
t38 = -t220 * t482 + t227 * t251 + t518 * qJDD(1) + (-t105 - t225) * qJD(1);
t18 = qJD(1) * t109 + qJDD(1) * t162 + t226 * t502 + t340 * t444 + t368;
t17 = -t227 * t426 + t342 * t444 + t468 * qJDD(1) + (-t107 - t459 + t569) * qJD(1) + t471;
t4 = -t164 * t226 + t516 * t227 + (t107 * t340 + t109 * t342) * qJD(2) + t378;
t1 = qJDD(4) * t339 - t515 * t226 + t469 * t227 + (t531 * t340 + t530 * t342 + t478) * qJD(2) + t378;
t5 = [-m(2) * (-g(1) * t252 + g(2) * t257) - t645 * t227 / 0.2e1 + (((t54 + t626 + t628) * t340 + ((t681 + t703) * t342 + t646 + t692 + t702) * t342) * qJD(2) + t657) * t451 + (t676 * qJD(2) + t678 * t339 - t677 * t341) * qJD(1) + (-(-t40 + t617 + t620 + t625) * t41 + t40 * (-t430 - t431) + t41 * (t467 + t618) + (t40 * t571 * t536 + (t40 * rSges(5,3) * t340 + t41 * t440) * t339) * qJD(2) + ((t339 * t571 + t627) * t342 * t40 + (t40 * (-pkin(5) - t584) + (-t319 - t326 + t627) * t41) * t340) * qJD(1) + (t3 - g(2)) * (t441 + t517) + (t2 - g(1)) * (t340 * t382 + t335 + t515)) * m(5) + ((t271 + (-t318 + (t341 * t570 - t579) * qJD(2)) * t340 + ((t339 * t570 + t438) * t342 + (-rSges(4,1) - pkin(5)) * t340) * qJD(1)) * t43 + (t18 - g(2)) * t85 + (t17 - g(1)) * (t335 + (t438 - t319) * t340 + t450) + (-qJD(1) * t164 - t380 + t43 - t625 - pkin(2) * t462 + t442 + t467 + (t382 - t325) * t486) * t44) * m(4) + (-(-qJD(1) * t159 - t229 - t464 - t63) * t64 + t64 * (t316 + t498) + (t251 * t575 - t577) * qJD(2) + ((-pkin(1) - t256) * t573 + (t63 * (-rSges(3,3) - pkin(5)) + t64 * t458) * t340) * qJD(1) + (t39 - g(2)) * t113 + (t38 - g(1)) * (t458 * t340 + t335 + t494)) * m(3) + (t73 + t611) * t587 + (m(2) * (t252 ^ 2 + t257 ^ 2) + Icges(2,3) - t675) * qJDD(1) + (t610 + t629) * t589 + (t631 + t633) * t453 + (((t342 * t638 + t612 - t628) * t342 + (t340 * t638 + t130 + t432 + t613 - t694) * t340) * qJD(2) + t634 - t649) * t454 + (t630 - t632 + t635) * t452; t663 * t589 + t662 * t587 + (qJD(1) * t631 + t604 * qJD(2) + qJDD(1) * t629 + t612 * t226 + t613 * t227) * t340 / 0.2e1 - (t630 * qJD(1) + t603 * qJD(2) + t682 * qJDD(1) + t614 * t226 + t615 * t227) * t342 / 0.2e1 - ((t654 * t339 - t659 * t341) * qJD(2) + (t652 * t339 + t653 * t341) * qJD(1)) * qJD(1) / 0.2e1 + (t632 * t342 + t633 * t340 + (t611 * t340 + t610 * t342) * qJD(1)) * qJD(1) / 0.2e1 + (t610 * t340 - t611 * t342) * qJDD(1) / 0.2e1 + t634 * t486 / 0.2e1 + t635 * t485 / 0.2e1 + ((-t484 * t608 + t602) * t340 + ((t340 * t609 + t606) * qJD(2) + t605) * t342) * t454 + ((t613 * t340 + t342 * t612) * qJD(1) + t604) * t453 + ((t340 * t615 + t614 * t342) * qJD(1) + t603) * t452 + ((-t482 * t609 - t602) * t342 + ((t342 * t608 + t606) * qJD(2) + t605) * t340) * t451 + (-g(1) * (t282 + t309) - g(2) * (t279 + t303) - g(3) * t607 - (g(1) * t440 + t473 * t583) * t339 + t1 * t514 + (-t1 * t515 + t3 * t434) * t340 + (t1 * t517 + t2 * t434) * t342 + (-t511 + (qJ(4) * t537 + t342 * t434 - t205) * qJD(1) + (t479 + t592) * t340) * t41 + (t339 * t477 + t207 - t277 + t592 * t342 + (-t340 * t425 + t196 + t200) * qJD(1)) * t40 + (-t470 - t478 - (t200 * t340 + t205 * t342 - t493 * t549) * qJD(2) + t472 + (qJD(1) * t469 + t531) * t340 + (t530 - t617) * t342) * t37) * m(5) + (-g(1) * (t201 + t202) - g(2) * t624 + g(3) * t501 + t43 * t207 + t4 * t514 + t42 * t472 + (t17 * t502 + t43 * t512 + t4 * t162 + t42 * t109 + (-t42 * t164 + t44 * t502) * qJD(1)) * t342 + (t18 * t502 + t44 * t512 - t4 * t164 + t42 * t107 + (t42 * t516 - t426 * t43) * qJD(1)) * t340 - t43 * (-qJD(1) * t624 + t277) - t44 * (qJD(1) * t202 + t511) - t42 * t470 - ((t42 * t202 + t43 * t501) * t342 + (t42 * t197 + t44 * t501) * t340) * qJD(2)) * m(4) + (-(t198 * t63 - t577) * qJD(1) - (t59 * (-t198 * t340 - t203 * t342) + t412 * t256) * qJD(2) + (qJD(2) * t399 + t159 * t226 - t160 * t227) * t390 + t59 * ((t159 * t342 - t160 * t340) * qJD(1) + t399) + t412 * t220 + (-t39 * t340 - t38 * t342 + (-t342 * t64 + t575) * qJD(1)) * t251 + g(1) * t203 + g(2) * t198 - g(3) * t256) * m(3); (-m(4) - m(5)) * (-g(3) * t341 + t339 * t437) - m(4) * (t175 * t42 + t176 * t44 + t177 * t43) - m(5) * (t175 * t37 + t176 * t41 + t177 * t40) + 0.2e1 * ((t43 * t482 + t44 * t484 - t4) * t591 + (t40 * t482 + t41 * t484 - t1) * t590) * t341 + 0.2e1 * ((qJD(2) * t42 + t17 * t342 + t18 * t340 - t43 * t486 + t44 * t485) * t591 + (-t40 * t486 + t41 * t485 + t367) * t590) * t339; ((t1 - g(3)) * t339 + (-t37 * t616 + t367 - t437) * t341) * m(5);];
tau = t5;
