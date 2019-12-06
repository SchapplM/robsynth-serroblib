% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:31
% EndTime: 2019-12-05 16:44:08
% DurationCPUTime: 25.38s
% Computational Cost: add. (17807->608), mult. (15317->803), div. (0->0), fcn. (12004->6), ass. (0->359)
t719 = Icges(5,4) + Icges(6,4);
t720 = Icges(5,1) + Icges(6,1);
t718 = Icges(5,2) + Icges(6,2);
t714 = Icges(5,6) + Icges(6,6);
t715 = Icges(5,5) + Icges(6,5);
t355 = qJ(3) + qJ(4);
t349 = cos(t355);
t717 = t719 * t349;
t348 = sin(t355);
t716 = t719 * t348;
t713 = t349 * t720 - t716;
t712 = -t348 * t718 + t717;
t353 = pkin(8) + qJ(2);
t347 = cos(t353);
t711 = t714 * t347;
t710 = Icges(5,3) + Icges(6,3);
t346 = sin(t353);
t546 = t346 * t349;
t547 = t346 * t348;
t690 = t546 * t719 - t547 * t718 - t711;
t706 = t349 * t718 + t716;
t709 = t348 * t720 + t717;
t708 = t719 * t547;
t707 = t715 * t347;
t689 = t346 * t714 + t347 * t712;
t688 = t546 * t720 - t707 - t708;
t687 = t346 * t715 + t347 * t713;
t686 = -t348 * t714 + t349 * t715;
t677 = t689 * t348;
t700 = t710 * t347;
t691 = -t546 * t715 + t547 * t714 + t700;
t620 = t677 + t691;
t542 = t347 * t349;
t623 = t346 * t710 + t347 * t686;
t670 = t346 * t623 + t542 * t687;
t701 = t690 * t348;
t621 = -t349 * t688 + t701;
t694 = t346 * t621;
t702 = t347 * t620 - t670 - t694;
t699 = t712 + t709;
t698 = t713 - t706;
t697 = -t347 * t706 + t687;
t696 = -t347 * t709 - t689;
t695 = t346 * t709 + t690;
t693 = t348 * t706 - t349 * t709;
t692 = t687 * t546;
t354 = qJD(3) + qJD(4);
t685 = t698 * t354;
t684 = t699 * t354;
t683 = -qJD(2) * t687 + t354 * t695;
t682 = t696 * t354 + (-t346 * t713 + t707) * qJD(2);
t254 = t346 * t354;
t681 = qJD(2) * t689 - t254 * t706 + t354 * t688;
t680 = t697 * t354 + (-t346 * t712 + t711) * qJD(2);
t262 = Icges(5,5) * t348 + Icges(5,6) * t349;
t260 = Icges(6,5) * t348 + Icges(6,6) * t349;
t552 = t260 * t347;
t679 = -t262 * t347 - t552;
t553 = t260 * t346;
t678 = t262 * t346 + t553;
t358 = -pkin(7) - pkin(6);
t352 = -qJ(5) + t358;
t676 = rSges(6,3) - t352;
t675 = t346 * t693 - t679;
t674 = -t347 * t693 + t678;
t672 = t547 * t689 - t692;
t671 = t623 * qJD(2);
t669 = t346 * t691 - t542 * t688;
t645 = t347 * t691 - t694;
t644 = -t347 * t623 - t672;
t543 = t347 * t348;
t643 = -t543 * t690 - t669;
t642 = -t543 * t689 + t670;
t668 = t685 * t349 - t684 * t348 + (t260 + t262) * qJD(2);
t667 = -t348 * t680 + t349 * t682 + t671;
t666 = qJD(2) * t691 + t348 * t681 + t349 * t683;
t255 = t347 * t354;
t665 = qJD(2) * t698 + t254 * t696 + t255 * t695;
t664 = (t546 * t718 - t688 + t708) * t255 + t697 * t254 + t699 * qJD(2);
t663 = qJD(2) * t693 + t354 * t686;
t393 = t262 * t354;
t662 = qJD(2) * t621 - t346 * t393 - t354 * t553 + t671;
t661 = -t347 * t393 - t354 * t552 + (-t346 * t686 - t349 * t687 + t677 + t700) * qJD(2);
t660 = t674 * qJD(2);
t357 = cos(qJ(3));
t343 = pkin(3) * t357 + pkin(2);
t594 = pkin(4) * t349;
t283 = t343 + t594;
t332 = t346 * rSges(6,3);
t659 = rSges(6,1) * t542 - rSges(6,2) * t543 + t283 * t347 + t332;
t356 = sin(qJ(3));
t486 = qJD(3) * t356;
t478 = pkin(3) * t486;
t539 = t348 * t354;
t252 = -pkin(4) * t539 - t478;
t321 = qJD(5) * t346;
t491 = qJD(2) * t346;
t471 = t348 * t491;
t538 = t349 * t354;
t474 = t347 * t538;
t490 = qJD(2) * t347;
t658 = rSges(6,3) * t490 + (t471 - t474) * rSges(6,2) + t347 * t252 + t321;
t657 = t675 * qJD(2);
t656 = rSges(6,2) * t547 - t346 * t283 + t347 * t676;
t655 = -t346 * t662 + t347 * t666;
t654 = t346 * t661 + t347 * t667;
t653 = t346 * t666 + t347 * t662;
t652 = t346 * t667 - t347 * t661;
t651 = t254 * t644 - t255 * t645 - t657;
t650 = t254 * t642 - t255 * t643 + t660;
t649 = t348 * t683 - t349 * t681;
t648 = t348 * t682 + t349 * t680;
t647 = t346 * t663 + t347 * t668;
t646 = t346 * t668 - t347 * t663;
t641 = t348 * t688 + t349 * t690;
t640 = t348 * t687 + t349 * t689;
t475 = t347 * t539;
t386 = -t349 * t491 - t475;
t465 = t347 * t486;
t430 = pkin(3) * t465;
t496 = -t352 + t358;
t505 = t283 - t343;
t586 = rSges(6,1) * t386 + t430 + (-t346 * t505 + t347 * t496) * qJD(2) + t658;
t272 = rSges(6,1) * t348 + rSges(6,2) * t349;
t222 = t272 * t346;
t588 = rSges(6,1) * t349;
t274 = -rSges(6,2) * t348 + t588;
t322 = qJD(5) * t347;
t299 = t346 * t478;
t500 = t358 * t491 + t299;
t585 = -t354 * t222 + t252 * t346 - t322 + t500 + (-t346 * t352 + t332 + (t274 + t505) * t347) * qJD(2);
t319 = t347 * t358;
t504 = t343 * t346 + t319;
t533 = -rSges(6,1) * t546 + t504 + t656;
t292 = t347 * t343;
t532 = t346 * t496 - t292 + t659;
t629 = -t348 * t664 + t349 * t665;
t628 = qJD(2) * t686 + t254 * t679 + t255 * t678;
t626 = 0.2e1 * qJD(3);
t234 = t274 * t354;
t429 = qJD(2) * t354;
t246 = t346 * t429;
t537 = t357 * qJD(3) ^ 2;
t400 = -pkin(3) * t347 * t537 + qJD(2) * t299;
t338 = t346 * pkin(6);
t591 = pkin(2) - t343;
t143 = (-t347 * t591 - t338) * qJD(2) - t500;
t259 = pkin(2) * t347 + t338;
t251 = t259 * qJD(2);
t534 = -t143 - t251;
t26 = -t234 * t255 + t246 * t272 + (t246 * t348 - t255 * t538) * pkin(4) + (t322 + t534 - t585) * qJD(2) + t400;
t595 = pkin(4) * t348;
t453 = -t272 - t595;
t432 = -t346 * t358 + t292;
t177 = t432 - t259;
t472 = -t177 - t532;
t61 = -t299 - t322 + t453 * t254 + (t259 - t472) * qJD(2);
t625 = qJD(2) * t61 + t26;
t339 = t347 * pkin(6);
t258 = pkin(2) * t346 - t339;
t176 = t258 - t504;
t253 = qJD(2) * t258;
t624 = qJD(2) * t176 - t253;
t622 = t321 - t430;
t350 = Icges(4,4) * t357;
t414 = -Icges(4,2) * t356 + t350;
t305 = Icges(4,1) * t356 + t350;
t619 = rSges(5,1) * t539 + rSges(5,2) * t538;
t273 = rSges(5,1) * t348 + rSges(5,2) * t349;
t223 = t273 * t346;
t225 = t273 * t347;
t587 = rSges(5,2) * t348;
t589 = rSges(5,1) * t349;
t275 = -t587 + t589;
t193 = rSges(5,1) * t546 - rSges(5,2) * t547 - rSges(5,3) * t347;
t333 = t346 * rSges(5,3);
t502 = rSges(5,1) * t542 + t333;
t195 = -rSges(5,2) * t543 + t502;
t487 = qJD(3) * t347;
t488 = qJD(3) * t346;
t464 = -t176 * t488 + t177 * t487 + qJD(1);
t65 = t193 * t254 + t195 * t255 + t464;
t387 = -t255 * t273 - t430;
t522 = t176 - t258;
t68 = (-t193 + t522) * qJD(2) + t387;
t521 = -t177 - t195;
t69 = -t254 * t273 - t299 + (t259 - t521) * qJD(2);
t617 = -(qJD(2) * t223 - t255 * t275) * t68 - t65 * (-t223 * t254 - t225 * t255) - t69 * (-qJD(2) * t225 - t254 * t275);
t544 = t346 * t357;
t545 = t346 * t356;
t569 = Icges(4,3) * t347;
t198 = Icges(4,5) * t544 - Icges(4,6) * t545 - t569;
t312 = Icges(4,4) * t545;
t578 = Icges(4,5) * t347;
t202 = Icges(4,1) * t544 - t312 - t578;
t572 = Icges(4,6) * t347;
t200 = Icges(4,4) * t544 - Icges(4,2) * t545 - t572;
t556 = t200 * t356;
t406 = -t202 * t357 + t556;
t80 = -t198 * t347 - t346 * t406;
t302 = Icges(4,5) * t357 - Icges(4,6) * t356;
t301 = Icges(4,5) * t356 + Icges(4,6) * t357;
t388 = qJD(3) * t301;
t581 = Icges(4,4) * t356;
t306 = Icges(4,1) * t357 - t581;
t203 = Icges(4,5) * t346 + t306 * t347;
t201 = Icges(4,6) * t346 + t347 * t414;
t555 = t201 * t356;
t405 = -t203 * t357 + t555;
t612 = -t347 * t388 + (-t302 * t346 + t405 + t569) * qJD(2);
t199 = Icges(4,3) * t346 + t302 * t347;
t493 = qJD(2) * t199;
t611 = qJD(2) * t406 - t346 * t388 + t493;
t303 = Icges(4,2) * t357 + t581;
t401 = t303 * t356 - t305 * t357;
t608 = t401 * qJD(2) + qJD(3) * t302;
t514 = -Icges(4,2) * t544 + t202 - t312;
t516 = t305 * t346 + t200;
t607 = -t356 * t514 - t357 * t516;
t604 = t246 / 0.2e1;
t247 = t347 * t429;
t603 = t247 / 0.2e1;
t602 = -t254 / 0.2e1;
t601 = t254 / 0.2e1;
t600 = -t255 / 0.2e1;
t599 = t255 / 0.2e1;
t598 = t346 / 0.2e1;
t597 = -t347 / 0.2e1;
t596 = pkin(3) * t356;
t593 = -qJD(2) / 0.2e1;
t592 = qJD(2) / 0.2e1;
t590 = rSges(4,1) * t357;
t334 = t346 * rSges(4,3);
t41 = -t254 * t533 + t255 * t532 + t464;
t566 = qJD(2) * t41;
t497 = rSges(4,2) * t545 + rSges(4,3) * t347;
t206 = rSges(4,1) * t544 - t497;
t307 = rSges(4,1) * t356 + rSges(4,2) * t357;
t467 = t307 * t487;
t119 = -t467 + (-t206 - t258) * qJD(2);
t564 = t119 * t346;
t563 = t119 * t347;
t468 = t307 * t488;
t540 = t347 * t357;
t541 = t347 * t356;
t207 = rSges(4,1) * t540 - rSges(4,2) * t541 + t334;
t512 = t207 + t259;
t120 = qJD(2) * t512 - t468;
t244 = t307 * t347;
t562 = t120 * t244;
t554 = t254 * t349;
t549 = t301 * t346;
t548 = t301 * t347;
t320 = pkin(6) * t490;
t142 = -t430 - t320 + (t346 * t591 - t319) * qJD(2);
t245 = qJD(2) * (-pkin(2) * t491 + t320);
t535 = qJD(2) * t142 + t245;
t526 = -t176 * t346 + t177 * t347;
t525 = t193 * t346 + t195 * t347;
t524 = -t198 * t346 - t202 * t540;
t523 = t199 * t346 + t203 * t540;
t515 = -t305 * t347 - t201;
t513 = -t303 * t347 + t203;
t281 = pkin(4) * t471;
t510 = t272 * t491 + t281;
t506 = rSges(5,2) * t471 + rSges(5,3) * t490;
t489 = qJD(2) * t356;
t501 = rSges(4,2) * t346 * t489 + rSges(4,3) * t490;
t499 = -t303 + t306;
t498 = t305 + t414;
t492 = qJD(2) * t302;
t485 = qJD(3) * t357;
t137 = -t346 * t401 - t548;
t483 = t137 * qJD(2);
t482 = qJD(2) * qJD(3);
t116 = rSges(5,1) * t386 - rSges(5,2) * t474 + t506;
t118 = -t354 * t223 + (t275 * t347 + t333) * qJD(2);
t481 = t116 * t347 + t118 * t346 + t193 * t490;
t477 = pkin(3) * t485;
t476 = t346 * t537;
t473 = t142 * t347 + t143 * t346 - t176 * t490;
t469 = t347 * t489;
t463 = -pkin(2) - t590;
t462 = t347 * t482;
t461 = t491 / 0.2e1;
t460 = t490 / 0.2e1;
t459 = -t488 / 0.2e1;
t456 = t487 / 0.2e1;
t455 = -t273 - t596;
t297 = -t595 - t596;
t454 = t297 + t596;
t451 = t356 * (-t346 ^ 2 - t347 ^ 2);
t224 = t272 * t347;
t442 = -t222 * t254 - t224 * t255;
t173 = t203 * t544;
t441 = t199 * t347 - t173;
t438 = -t198 + t555;
t437 = -qJD(2) * t224 - t254 * t274;
t428 = -t346 * t533 + t347 * t532;
t427 = -pkin(4) * t538 - t234;
t235 = t275 * t354;
t424 = -t235 - t477;
t421 = qJD(2) * t222 + (-t274 - t594) * t255;
t419 = -rSges(4,2) * t356 + t590;
t418 = -t346 * t69 - t347 * t68;
t411 = -t120 * t346 - t563;
t121 = t200 * t357 + t202 * t356;
t122 = t201 * t357 + t203 * t356;
t404 = t206 * t346 + t207 * t347;
t399 = t346 * t585 + t347 * t586 - t490 * t533;
t398 = -t272 + t297;
t243 = t307 * t346;
t81 = -t201 * t545 - t441;
t395 = (t346 * t81 - t347 * t80) * qJD(3);
t82 = -t200 * t541 - t524;
t83 = -t201 * t541 + t523;
t394 = (t346 * t83 - t347 * t82) * qJD(3);
t390 = qJD(3) * t305;
t389 = qJD(3) * t303;
t385 = t427 - t477;
t382 = -t177 * t346 * t482 + t142 * t487 + t143 * t488 - t176 * t462;
t381 = -t356 * t513 + t357 * t515;
t380 = (-t356 * t498 + t357 * t499) * qJD(2);
t379 = t255 * t453 + t622;
t139 = -rSges(4,2) * t347 * t485 + (-t357 * t491 - t465) * rSges(4,1) + t501;
t140 = -qJD(3) * t243 + (t347 * t419 + t334) * qJD(2);
t377 = t139 * t347 + t140 * t346 + (t206 * t347 - t207 * t346) * qJD(2);
t134 = qJD(2) * t201 - t346 * t389;
t136 = qJD(2) * t203 - t346 * t390;
t366 = qJD(2) * t198 - qJD(3) * t121 - t134 * t356 + t136 * t357;
t133 = -t347 * t389 + (-t346 * t414 + t572) * qJD(2);
t135 = -t347 * t390 + (-t306 * t346 + t578) * qJD(2);
t365 = -qJD(3) * t122 - t133 * t356 + t135 * t357 + t493;
t277 = t414 * qJD(3);
t278 = t306 * qJD(3);
t364 = qJD(2) * t301 - t277 * t356 + t278 * t357 + (-t303 * t357 - t305 * t356) * qJD(3);
t363 = (t346 * t644 - t347 * t645) * t604 + (t346 * t642 - t347 * t643) * t603 + (t346 * t628 + t347 * t629) * t602 + (t655 * t347 + t654 * t346 + (t346 * t643 + t347 * t642) * qJD(2)) * t601 + (t653 * t347 + t652 * t346 + (t346 * t645 + t347 * t644) * qJD(2)) * t600 + (t346 * t629 - t347 * t628) * t599 + (qJD(2) * t647 + t246 * t643 + t247 * t642 + t254 * t654 + t255 * t655) * t598 + (qJD(2) * t646 + t246 * t645 + t247 * t644 + t254 * t652 + t255 * t653) * t597 + (t348 * t665 + t349 * t664) * t593 + (t649 * t347 + t648 * t346 + (t346 * t641 + t347 * t640) * qJD(2)) * t592 + t651 * t461 + t650 * t460;
t282 = t419 * qJD(3);
t205 = t454 * t347;
t204 = t454 * t346;
t138 = -t347 * t401 + t549;
t127 = t138 * qJD(2);
t98 = qJD(3) * t404 + qJD(1);
t79 = -t282 * t487 + (-t140 - t251 + t468) * qJD(2);
t78 = -t282 * t488 + t245 + (t139 - t467) * qJD(2);
t67 = t346 * t364 - t347 * t608;
t66 = t346 * t608 + t347 * t364;
t64 = -qJD(3) * t405 + t133 * t357 + t135 * t356;
t63 = -qJD(3) * t406 + t134 * t357 + t136 * t356;
t62 = t377 * qJD(3);
t60 = (t522 + t533) * qJD(2) + t379;
t49 = -t235 * t255 + t246 * t273 + (-t118 + t534) * qJD(2) + t400;
t48 = qJD(2) * t116 - t235 * t254 - t247 * t273 + (-t356 * t462 - t476) * pkin(3) + t535;
t43 = t127 + t394;
t42 = t395 + t483;
t25 = -pkin(3) * t476 - t234 * t254 - t247 * t272 + (-t247 * t348 - t254 * t538) * pkin(4) + (t586 + t622) * qJD(2) + t535;
t16 = t116 * t255 + t118 * t254 + t193 * t247 - t195 * t246 + t382;
t9 = -t246 * t532 - t247 * t533 + t254 * t585 + t255 * t586 + t382;
t1 = [m(4) * t62 + m(5) * t16 + m(6) * t9; (t127 + ((t81 - t173 + (t199 + t556) * t347 + t524) * t347 + t523 * t346) * qJD(3)) * t456 + (((t623 + t701) * t347 + t644 + t669 + t672) * t255 + (t645 - t702) * t254 + t660) * t599 + (t42 - t483 + ((t347 * t438 - t523 + t83) * t347 + (t346 * t438 + t441 + t82) * t346) * qJD(3)) * t459 + (t66 + t64) * t488 / 0.2e1 + (-qJD(3) * t401 + t277 * t357 + t278 * t356 + t685 * t348 + t684 * t349) * qJD(2) + (-(qJD(2) * t533 + t379 - t60 + t624) * t61 + t26 * t656 + t60 * t322 + t25 * t659 + t61 * (-rSges(6,1) * t475 + t658) + (-t26 * t588 + t60 * (rSges(6,1) * t539 + rSges(6,2) * t538 - t252) - t25 * t352) * t346 + ((t60 * (-t274 - t283) - t61 * t352) * t347 + (-t60 * t676 + t61 * (-t283 - t588)) * t346) * qJD(2)) * m(6) + (-(-qJD(2) * t193 + t387 + t624 - t68) * t69 + t49 * (-t193 - t504) + t68 * (t619 * t346 + t500) + t48 * (t432 + t502) + t69 * t506 + (-t48 * t587 + t69 * (-t478 - t619)) * t347 + ((-t68 * rSges(5,3) + t69 * (-t343 - t589)) * t346 + (t68 * (-t275 - t343) - t69 * t358) * t347) * qJD(2)) * m(5) + (-(-qJD(2) * t206 - t119 - t253 - t467) * t120 + t79 * (t346 * t463 + t339 + t497) + t78 * t512 + t120 * (t320 + t501) + (t307 * t564 - t562) * qJD(3) + ((-pkin(2) - t419) * t563 + (t119 * (-rSges(4,3) - pkin(6)) + t120 * t463) * t346) * qJD(2)) * m(4) - (t63 + t67 + t43) * t487 / 0.2e1 + (t641 - t675) * t604 + (t640 + t674) * t603 + ((t642 + t702) * t255 + ((t621 + t623) * t347 + t620 * t346 + t643 - t692) * t254 + t651 + t657) * t602 + (t647 + t648) * t601 + (t646 - t649 + t650) * t600 + ((t121 + t137) * t346 + (t122 + t138) * t347) * t482 / 0.2e1; ((-t487 * t549 - t492) * t347 + (t380 + (t381 * t346 + (t548 - t607) * t347) * qJD(3)) * t346) * t456 + ((t356 * t499 + t357 * t498) * qJD(2) + ((t346 * t513 - t347 * t514) * t357 + (t346 * t515 + t347 * t516) * t356) * qJD(3)) * t593 + (t346 * t64 - t347 * t63 + (t121 * t346 + t122 * t347) * qJD(2)) * t592 + t363 + ((-t488 * t548 + t492) * t346 + (t380 + (-t607 * t347 + (t549 + t381) * t346) * qJD(3)) * t347) * t459 + (qJD(2) * t66 + (-(t346 * t611 + t347 * t366) * t347 + (t346 * t612 + t347 * t365) * t346 + (t346 * t82 + t347 * t83) * qJD(2)) * t626) * t598 + (qJD(2) * t67 + (-(t346 * t366 - t347 * t611) * t347 + (t346 * t365 - t347 * t612) * t346 + (t346 * t80 + t347 * t81) * qJD(2)) * t626) * t597 + (t395 + t42) * t461 + (t394 + t43) * t460 + (-t60 * (-qJD(2) * t204 + t421) - t61 * (-pkin(4) * t554 + qJD(2) * t205 + t437) - t41 * (t204 * t254 + t205 * t255 + t442) - (-t61 * t469 + ((-t346 * t61 - t347 * t60) * t357 + t41 * t451) * qJD(3)) * pkin(3) + t60 * t510 + t9 * (t428 + t526) + t41 * (t399 + t473) + (t385 * t60 + t398 * t625) * t347 + (t25 * t398 + t385 * t61 + t472 * t566) * t346) * m(6) + (-(-t69 * t469 + (t357 * t418 + t451 * t65) * qJD(3)) * pkin(3) + t16 * (t525 + t526) + t65 * (t473 + t481) + (t424 * t68 + (qJD(2) * t69 + t49) * t455) * t347 + (t48 * t455 + t69 * t424 + (t273 * t68 + t521 * t65) * qJD(2)) * t346 + t617) * m(5) + (t62 * t404 + t98 * t377 + t411 * t282 + (-t78 * t346 - t79 * t347 + (-t120 * t347 + t564) * qJD(2)) * t307 - (t119 * t243 - t562) * qJD(2) - (t98 * (-t243 * t346 - t244 * t347) + t411 * t419) * qJD(3)) * m(4); t363 + (t9 * t428 + t41 * t399 + (t427 * t61 - t532 * t566) * t346 - t61 * t437 - t41 * t442 - (-t61 * t554 + (-t61 * t490 + t41 * (-t254 * t346 - t255 * t347)) * t348) * pkin(4) + (t25 * t346 + t347 * t625) * t453 + (t347 * t427 - t281 - t421 + t510) * t60) * m(6) + (t16 * t525 + t65 * (-t195 * t491 + t481) + t418 * t235 + (-t48 * t346 - t49 * t347 + (t346 * t68 - t347 * t69) * qJD(2)) * t273 + t617) * m(5); 0.2e1 * (t25 * t597 + t26 * t598 - t41 * (-t254 * t347 + t255 * t346) / 0.2e1) * m(6);];
tauc = t1(:);
