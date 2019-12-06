% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1]';
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
% Datum: 2019-12-05 15:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRP6_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP6_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP6_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRPRP6_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP6_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP6_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP6_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:40:21
% EndTime: 2019-12-05 15:41:07
% DurationCPUTime: 37.32s
% Computational Cost: add. (8315->571), mult. (23329->869), div. (0->0), fcn. (22286->6), ass. (0->274)
t649 = Icges(5,4) - Icges(6,5);
t650 = Icges(5,1) + Icges(6,1);
t618 = Icges(6,4) + Icges(5,5);
t616 = Icges(5,2) + Icges(6,3);
t615 = Icges(5,6) - Icges(6,6);
t663 = Icges(6,2) + Icges(5,3);
t294 = cos(qJ(4));
t665 = t649 * t294;
t292 = sin(qJ(4));
t664 = t649 * t292;
t290 = sin(pkin(7));
t293 = sin(qJ(2));
t480 = t292 * t293;
t431 = t290 * t480;
t291 = cos(pkin(7));
t482 = t291 * t294;
t235 = -t431 + t482;
t662 = t649 * t235;
t479 = t293 * t294;
t352 = -t290 * t292 + t291 * t479;
t661 = t649 * t352;
t234 = t290 * t479 + t291 * t292;
t660 = t649 * t234;
t659 = -t616 * t294 - t664;
t658 = t618 * t292 + t615 * t294;
t657 = t650 * t292 + t665;
t233 = t290 * t294 + t291 * t480;
t656 = t649 * t233;
t295 = cos(qJ(2));
t481 = t291 * t295;
t606 = t616 * t352 + t615 * t481 + t656;
t484 = t290 * t295;
t605 = t616 * t234 + t615 * t484 - t662;
t604 = t618 * t233 + t615 * t352 + t663 * t481;
t601 = t615 * t234 - t618 * t235 + t663 * t484;
t637 = t650 * t233 + t618 * t481 + t661;
t636 = -t650 * t235 + t618 * t484 + t660;
t654 = t658 * t293 + t663 * t295;
t653 = (t615 * t292 - t618 * t294) * t295;
t633 = -t663 * t293 + t658 * t295;
t595 = t615 * t293 + t659 * t295;
t594 = -t618 * t293 + t657 * t295;
t613 = t636 * t233 + t605 * t352 + t601 * t481;
t612 = t606 * t234 - t637 * t235 + t604 * t484;
t645 = t293 * (Icges(4,4) - Icges(3,5)) + t295 * (Icges(4,5) - Icges(3,6));
t558 = t645 * t290;
t557 = t645 * t291;
t570 = t637 * t233 + t606 * t352 + t604 * t481;
t611 = t605 * t234 - t636 * t235 + t601 * t484;
t608 = -t594 * t233 + t595 * t352 - t633 * t481;
t607 = t595 * t234 + t594 * t235 - t633 * t484;
t442 = qJD(2) * t295;
t426 = t292 * t442;
t148 = qJD(4) * t352 + t291 * t426;
t149 = qJD(4) * t233 - t442 * t482;
t443 = qJD(2) * t293;
t427 = t291 * t443;
t652 = t618 * t148 - t615 * t149 - t427 * t663;
t150 = qJD(4) * t234 + t290 * t426;
t151 = -qJD(4) * t431 + (qJD(4) * t291 + t290 * t442) * t294;
t428 = t290 * t443;
t651 = t618 * t150 + t615 * t151 - t428 * t663;
t648 = t654 * qJD(2) + t653 * qJD(4);
t562 = t659 * t293 - t615 * t295;
t561 = t657 * t293 + t618 * t295;
t647 = (-t616 * t292 + t665) * t295;
t646 = (-t650 * t294 + t664) * t295;
t644 = -Icges(3,2) - Icges(4,3);
t628 = t636 * t292 + t605 * t294;
t627 = t637 * t292 + t606 * t294;
t643 = t612 * t291;
t642 = t613 * t290;
t641 = -t148 * t649 + t149 * t616 + t427 * t615;
t640 = -t150 * t649 - t151 * t616 + t428 * t615;
t639 = t148 * t650 - t149 * t649 - t427 * t618;
t638 = t150 * t650 + t151 * t649 - t428 * t618;
t635 = qJD(2) * t562 + qJD(4) * t647;
t634 = qJD(2) * t561 + qJD(4) * t646;
t632 = t295 * t648 + t443 * t633;
t631 = t295 * t651 - t443 * t601;
t630 = t295 * t652 - t443 * t604;
t629 = -t292 * t594 + t294 * t595;
t626 = t290 * t611 + t643;
t625 = t291 * t570 + t642;
t624 = (Icges(3,4) + Icges(4,6)) * t295;
t623 = t607 * t295;
t622 = t608 * t295;
t610 = t293 * t604 - t627 * t295;
t609 = t293 * t601 - t295 * t628;
t549 = t633 * t293;
t504 = Icges(3,4) * t293;
t614 = t624 * t295 + (-t504 + (Icges(3,1) + t644) * t295) * t293;
t439 = qJD(4) * t295;
t447 = qJD(2) * t290;
t261 = t291 * t439 + t447;
t526 = t261 / 0.2e1;
t445 = qJD(2) * t291;
t262 = t290 * t439 - t445;
t524 = t262 / 0.2e1;
t575 = t637 * t148 - t606 * t149 + t639 * t233 + t630 * t291 - t352 * t641;
t574 = t148 * t636 - t149 * t605 + t233 * t638 + t291 * t631 - t352 * t640;
t573 = t637 * t150 + t606 * t151 - t234 * t641 - t639 * t235 + t630 * t290;
t572 = t150 * t636 + t151 * t605 - t234 * t640 - t235 * t638 + t290 * t631;
t568 = -t295 * t629 - t549;
t400 = rSges(5,1) * t292 + rSges(5,2) * t294;
t325 = -rSges(5,3) * t293 + t295 * t400;
t603 = t261 * t325;
t602 = t262 * t325;
t288 = t290 ^ 2;
t289 = t291 ^ 2;
t448 = t288 + t289;
t598 = (t290 * t633 + t628) * t262 + (t291 * t633 + t627) * t261;
t597 = (-t150 * t594 + t151 * t595 - t234 * t635 - t235 * t634 + t290 * t632) * t293 + (-t293 * t626 + t623) * qJD(2);
t596 = (-t148 * t594 - t149 * t595 + t233 * t634 + t291 * t632 - t352 * t635) * t293 + (-t293 * t625 + t622) * qJD(2);
t593 = (t629 + t654) * t293;
t592 = t645 * qJD(2);
t590 = (-Icges(4,6) * t293 + t295 * t644 - t504) * t443 + ((Icges(3,1) + Icges(4,2)) * t293 + t624) * t442;
t440 = qJD(4) * t293;
t534 = t261 * (-t352 * t650 + t606 + t656) + t262 * (-t234 * t650 + t605 - t662) + t440 * (t595 - t646);
t589 = t261 * (-t233 * t616 + t637 + t661) + t262 * (t235 * t616 + t636 + t660) + t440 * (-t594 - t647);
t588 = t574 * t524 + t575 * t526 + t596 * qJD(4) / 0.2e1;
t586 = t653 * t440 + (t234 * t618 + t235 * t615) * t262 + (-t233 * t615 + t352 * t618) * t261;
t550 = t609 * t290 + t291 * t610;
t579 = (t261 * t612 + t262 * t611 + t440 * t607) * t290 + (t261 * t570 + t262 * t613 + t440 * t608) * t291;
t577 = qJD(4) * t597 + t261 * t573 + t262 * t572;
t576 = rSges(6,1) + pkin(4);
t542 = (t641 * t294 - t639 * t292 + (t292 * t606 - t294 * t637) * qJD(4) + t604 * qJD(2)) * t295 + (qJD(2) * t627 + t652) * t293;
t541 = (t640 * t294 - t638 * t292 + (t292 * t605 - t294 * t636) * qJD(4) + t601 * qJD(2)) * t295 + (qJD(2) * t628 + t651) * t293;
t571 = t261 * t610 + t262 * t609 + t440 * t568;
t567 = rSges(6,3) + qJ(5);
t566 = t595 * t290;
t565 = t595 * t291;
t564 = t594 * t290;
t563 = t594 * t291;
t560 = t592 * t290;
t559 = t592 * t291;
t556 = (qJD(4) * t593 + t598) * t295;
t483 = t291 * t293;
t555 = t590 * t291 + ((Icges(4,2) * t481 - Icges(4,6) * t483) * t293 + t614 * t291 - t558) * qJD(2);
t485 = t290 * t293;
t554 = t590 * t290 + ((Icges(4,2) * t484 - Icges(4,6) * t485) * t293 + t614 * t290 + t557) * qJD(2);
t397 = pkin(4) * t292 - qJ(5) * t294;
t399 = rSges(6,1) * t292 - rSges(6,3) * t294;
t553 = -rSges(6,2) * t293 + (t397 + t399) * t295;
t551 = -t261 * t604 - t262 * t601;
t398 = rSges(4,2) * t293 + rSges(4,3) * t295;
t543 = t448 * qJD(2) * t398;
t419 = -t443 / 0.2e1;
t536 = t568 * t442 + (t648 * t293 + (t635 * t294 - t634 * t292 + (t292 * t595 + t294 * t594) * qJD(4)) * t295 + (t293 * t629 - t295 * t633 - t550) * qJD(2)) * t293;
t533 = t586 * t295;
t271 = rSges(3,1) * t293 + rSges(3,2) * t295;
t346 = qJD(2) * t271;
t436 = qJD(5) * t294;
t285 = t295 * t436;
t263 = pkin(3) * t290 + pkin(6) * t481;
t264 = -pkin(3) * t291 + pkin(6) * t484;
t272 = pkin(2) * t295 + qJ(3) * t293;
t248 = t272 * t290;
t249 = t272 * t291;
t441 = qJD(3) * t295;
t358 = t248 * t447 + t249 * t445 + qJD(1) - t441;
t332 = t263 * t445 + t264 * t447 + t358;
t470 = rSges(6,2) * t484 - t234 * t567 - t235 * t576;
t471 = rSges(6,2) * t481 + t233 * t576 - t352 * t567;
t27 = t261 * t470 - t262 * t471 + t285 + t332;
t296 = qJD(2) ^ 2;
t519 = pkin(6) * t296;
t409 = t448 * t519;
t425 = t292 * t439;
t287 = qJD(3) * t293;
t279 = t290 * t287;
t269 = pkin(2) * t293 - qJ(3) * t295;
t344 = qJD(2) * t269;
t168 = -t290 * t344 + t279;
t281 = t291 * t287;
t169 = -t291 * t344 + t281;
t435 = qJD(2) * qJD(3);
t430 = t168 * t447 + t169 * t445 + t293 * t435;
t437 = qJD(5) * t234;
t517 = -rSges(6,2) * t428 + t150 * t576 - t151 * t567 - t437;
t438 = qJD(5) * t352;
t518 = -rSges(6,2) * t427 + t148 * t576 + t149 * t567 - t438;
t5 = -qJD(5) * t425 - t518 * t262 + t517 * t261 + (-t409 + (-t436 + (t290 * t471 - t291 * t470) * qJD(4)) * qJD(2)) * t293 + t430;
t531 = t27 * t518 + t471 * t5;
t528 = t448 * t419;
t527 = -t261 / 0.2e1;
t525 = -t262 / 0.2e1;
t520 = -t295 / 0.2e1;
t478 = t294 * t295;
t200 = rSges(6,2) * t295 + t293 * t399;
t259 = (-rSges(6,1) * t294 - rSges(6,3) * t292) * t295;
t336 = -t294 * t443 - t425;
t469 = t285 + t336 * qJ(5) + (t292 * t443 - t294 * t439) * pkin(4) + qJD(2) * t200 + qJD(4) * t259;
t468 = t233 * t567 + t352 * t576;
t467 = t234 * t576 - t235 * t567;
t466 = t290 * t168 + t291 * t169;
t465 = t553 * t290;
t464 = t553 * t291;
t455 = t397 * t293 + t200;
t453 = t290 * t248 + t291 * t249;
t228 = qJD(2) * t272 - t441;
t273 = -rSges(4,2) * t295 + rSges(4,3) * t293;
t452 = -t273 * qJD(2) - t228;
t451 = (-pkin(4) * t294 - qJ(5) * t292) * t295 + t259;
t450 = -t269 + t398;
t449 = -t272 - t273;
t434 = t295 * t519;
t429 = -t344 * t448 + t287;
t422 = t295 * t435;
t418 = t442 / 0.2e1;
t417 = -t440 / 0.2e1;
t416 = t440 / 0.2e1;
t414 = -pkin(6) * t293 - t269;
t413 = t293 * t448;
t412 = qJD(2) * t452;
t411 = qJD(2) * t450;
t410 = t291 * t263 + t290 * t264 + t453;
t404 = t325 + t414;
t403 = qJD(2) * t417;
t402 = qJD(4) * t418;
t401 = -pkin(6) * t442 - t228;
t274 = rSges(3,1) * t295 - rSges(3,2) * t293;
t38 = -t269 * t447 - t437 + t279 + t553 * t261 + (-pkin(6) * t447 + qJD(4) * t471) * t293;
t39 = -t269 * t445 - t438 + t281 - t553 * t262 + (-pkin(6) * t445 - qJD(4) * t470) * t293;
t396 = -t290 * t38 - t291 * t39;
t110 = rSges(5,1) * t233 + rSges(5,2) * t352 + rSges(5,3) * t481;
t388 = qJD(2) * t414;
t64 = t110 * t440 + t290 * t388 + t279 + t603;
t112 = -rSges(5,1) * t235 + rSges(5,2) * t234 + rSges(5,3) * t484;
t65 = -t112 * t440 + t291 * t388 + t281 - t602;
t389 = -t290 * t64 - t291 * t65;
t367 = t110 * t290 - t112 * t291;
t366 = t448 * t274;
t363 = t448 * t346;
t362 = t414 + t553;
t361 = t290 * t403;
t360 = t291 * t403;
t201 = rSges(5,3) * t295 + t293 * t400;
t260 = (-rSges(5,1) * t294 + rSges(5,2) * t292) * t295;
t121 = qJD(2) * t201 + qJD(4) * t260;
t359 = -t121 + t401;
t357 = -qJD(2) * t228 - t434;
t186 = rSges(4,1) * t290 + t273 * t291;
t187 = -rSges(4,1) * t291 + t273 * t290;
t70 = (t186 * t291 + t187 * t290) * qJD(2) + t358;
t351 = t70 * t398;
t347 = t401 - t469;
t40 = -t110 * t262 + t112 * t261 + t332;
t343 = t40 * t367;
t335 = t27 * t517 + t470 * t5;
t331 = t38 * t471 - t39 * t470;
t316 = -pkin(6) * t443 * t448 + t466;
t299 = (-t27 * t470 - t38 * t553) * t291 + (t27 * t471 + t39 * t553) * t290;
t282 = t291 * t441;
t280 = t290 * t441;
t276 = t291 * t422;
t275 = t290 * t422;
t142 = rSges(5,1) * t234 + rSges(5,2) * t235;
t138 = rSges(5,1) * t352 - rSges(5,2) * t233;
t123 = t291 * t411 + t281;
t122 = t290 * t411 + t279;
t94 = t363 * qJD(2);
t93 = t291 * t412 + t276;
t92 = t290 * t412 + t275;
t91 = qJD(2) * t366 + qJD(1);
t90 = rSges(5,1) * t150 + rSges(5,2) * t151 - rSges(5,3) * t428;
t88 = rSges(5,1) * t148 - rSges(5,2) * t149 - rSges(5,3) * t427;
t63 = qJD(2) * t543 + t430;
t35 = -t291 * t434 - t90 * t440 + t121 * t262 + t276 + (-t228 * t291 + (-t112 * t295 + t325 * t485) * qJD(4)) * qJD(2);
t34 = -t290 * t434 + t88 * t440 - t121 * t261 + t275 + (-t228 * t290 + (t110 * t295 - t325 * t483) * qJD(4)) * qJD(2);
t26 = t261 * t90 - t262 * t88 + (qJD(2) * qJD(4) * t367 - t409) * t293 + t430;
t7 = qJD(5) * t149 + t276 + t357 * t291 + t469 * t262 + (-t517 * t293 + (-t295 * t470 + t485 * t553) * qJD(2)) * qJD(4);
t6 = -qJD(5) * t151 + t275 + t357 * t290 - t469 * t261 + (t518 * t293 + (t295 * t471 - t483 * t553) * qJD(2)) * qJD(4);
t1 = [-m(3) * t94 + m(4) * t63 + m(5) * t26 + m(6) * t5; (((t233 * t561 - t352 * t562 - t642) * t293 + t622) * qJD(4) + (((t549 - t570) * qJD(4) + t551) * t293 + t556) * t291 + (t233 * t564 - t352 * t566) * t262 + (t233 * t563 - t352 * t565) * t261) * t527 + (t290 * t575 - t291 * t574) * t526 + (((-t234 * t562 - t235 * t561 - t643) * t293 + t623) * qJD(4) + (((t549 - t611) * qJD(4) + t551) * t293 + t556) * t290 + (-t234 * t566 - t235 * t564) * t262 + (-t234 * t565 - t235 * t563) * t261) * t525 + (t290 * t573 - t291 * t572) * t524 + t290 * t588 - t577 * t291 / 0.2e1 + (t554 * t289 + (t559 * t290 + (-t555 - t560) * t291) * t290) * t447 + (-t560 * t289 + (t555 * t290 + (-t554 + t559) * t291) * t290) * t445 - (t557 * qJD(2) * t288 - t558 * t290 * t445) * t447 / 0.2e1 + (t558 * qJD(2) * t289 - t291 * t557 * t447) * t445 / 0.2e1 + (((-t292 * t564 + t294 * t566 + t601) * t262 + (-t292 * t563 + t294 * t565 + t604) * t261 + t568 * qJD(4)) * t295 + (((-t292 * t561 + t294 * t562 - t633) * t295 - t550 + t593) * qJD(4) + t598) * t293) * t417 - t571 * t439 / 0.2e1 + (t290 * t610 - t291 * t609) * t402 + (t290 * t612 - t291 * t611) * t361 + (t570 * t290 - t291 * t613) * t360 + (-t39 * (-t285 * t291 + t282) - t38 * (-t285 * t290 + t280) - t27 * (-t293 * t436 + t429) - (-t27 * t464 + t39 * t455) * t262 - (t27 * t465 - t38 * t455) * t261 - (t396 * t272 + (-t27 * t413 + t295 * t396) * pkin(6)) * qJD(2) - (t331 * t295 + (t38 * t464 - t39 * t465 + t299) * t293) * qJD(4) + t5 * t410 + t27 * t316 + (t347 * t39 + t362 * t7 + t531) * t291 + (t347 * t38 + t362 * t6 + t335) * t290) * m(6) + (-t65 * (t201 * t262 + t282) - t64 * (-t201 * t261 + t280) - t40 * (t290 * t603 - t291 * t602 + t429) - (t389 * t272 + (t295 * t389 - t40 * t413) * pkin(6)) * qJD(2) - ((t110 * t64 - t112 * t65) * t295 + t343 * t293) * qJD(4) + t26 * t410 + t40 * t316 + (t26 * t110 + t35 * t404 + t359 * t65 + t40 * t88) * t291 + (t26 * t112 + t34 * t404 + t359 * t64 + t40 * t90) * t290) * m(5) + (t63 * t453 + (t123 * t452 + t63 * t186 + t450 * t93) * t291 + (t122 * t452 + t63 * t187 + t450 * t92) * t290 - t123 * t282 - t122 * t280 - ((t123 * t449 + t291 * t351) * t291 + (t122 * t449 + t290 * t351) * t290) * qJD(2) + (-t429 + t466 + t543) * t70) * m(4) + (-t363 * t91 - t366 * t94 + (t271 * t274 * t296 + t346 * t91) * t448) * m(3) + (t290 * t542 - t291 * t541 + t579) * t416; 0.2e1 * (t27 * t528 + t5 * t520) * m(6) + 0.2e1 * (t26 * t520 + t40 * t528) * m(5) + 0.2e1 * (t520 * t63 + t528 * t70) * m(4) + 0.2e1 * (m(4) * (qJD(2) * t70 + t290 * t92 + t291 * t93) / 0.2e1 + m(5) * (qJD(2) * t40 + t290 * t34 + t291 * t35) / 0.2e1 + m(6) * (qJD(2) * t27 + t290 * t6 + t291 * t7) / 0.2e1) * t293; (-t233 * t534 + t533 * t291 + t352 * t589) * t527 + ((t290 * t574 + t291 * t575) * t295 + t596) * t526 + (t234 * t589 + t534 * t235 + t533 * t290) * t525 + ((t290 * t572 + t291 * t573) * t295 + t597) * t524 + (qJD(4) * t536 + t261 * t542 + t262 * t541) * t293 / 0.2e1 + t577 * t484 / 0.2e1 + t481 * t588 + t571 * t418 + ((t292 * t534 - t294 * t589) * t295 + t586 * t293) * t417 + ((t290 * t541 + t291 * t542) * t295 + t536) * t416 + (t293 * t568 + t295 * t550) * t402 + (t607 * t293 + t295 * t626) * t361 + (t608 * t293 + t295 * t625) * t360 + (-(-t27 * t292 * t295 + t233 * t39 - t235 * t38) * qJD(5) - (-t27 * t468 + t39 * t451) * t262 - (t27 * t467 - t38 * t451) * t261 - (t38 * t468 - t39 * t467) * t440 + (qJD(2) * t299 + t38 * t518 - t39 * t517 - t470 * t7 + t471 * t6) * t293 + (t331 * qJD(2) + (-t38 * t469 + t553 * t6 + t335) * t291 + (t39 * t469 - t553 * t7 - t531) * t290) * t295) * m(6) + ((t34 * t110 - t35 * t112 + t64 * t88 - t65 * t90 + (t343 - (-t290 * t65 + t291 * t64) * t325) * qJD(2)) * t293 + (t65 * (-qJD(2) * t112 + t121 * t290) + t64 * (qJD(2) * t110 - t121 * t291) - t26 * t367 + t40 * (-t290 * t88 + t291 * t90) - (t290 * t35 - t291 * t34) * t325) * t295 - t65 * (-t142 * t440 + t260 * t262) - t64 * (t138 * t440 - t260 * t261) - t40 * (-t138 * t262 + t142 * t261)) * m(5) + t579 * t419; (-t352 * t7 - t234 * t6 + t5 * t478 + (-t234 * t440 - t262 * t478 + t149) * t39 + (t261 * t478 + t352 * t440 - t151) * t38 + (t234 * t261 - t262 * t352 + t336) * t27) * m(6);];
tauc = t1(:);
