% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPP1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP1_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP1_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:25
% EndTime: 2019-12-31 20:49:41
% DurationCPUTime: 9.96s
% Computational Cost: add. (30372->457), mult. (24012->570), div. (0->0), fcn. (21447->8), ass. (0->297)
t445 = qJ(1) + qJ(2);
t441 = sin(t445);
t446 = sin(qJ(3));
t597 = t446 * pkin(3);
t444 = qJ(3) + pkin(8);
t440 = cos(t444);
t439 = sin(t444);
t623 = rSges(6,1) + pkin(4);
t494 = t623 * t439;
t587 = rSges(6,3) + qJ(5);
t679 = t587 * t440 - t494;
t481 = t597 - t679;
t699 = t481 * t441;
t248 = t441 * t699;
t442 = cos(t445);
t554 = t439 * t442;
t701 = t554 * t699;
t669 = t587 * t439 + t623 * t440;
t700 = -t442 * rSges(6,2) + t669 * t441;
t421 = Icges(5,4) * t440;
t383 = Icges(5,1) * t439 + t421;
t579 = Icges(6,5) * t440;
t698 = Icges(6,1) * t439 + t383 - t579;
t420 = Icges(6,5) * t439;
t471 = Icges(6,3) * t440 - t420;
t580 = Icges(5,4) * t439;
t694 = Icges(5,2) * t440 + t471 + t580;
t262 = t481 * t442;
t651 = m(6) / 0.2e1;
t387 = rSges(5,1) * t439 + rSges(5,2) * t440;
t459 = t387 + t597;
t673 = t459 * t442;
t674 = t459 * t441;
t676 = t441 * t674 + t442 * t673;
t531 = (-t262 * t442 - t248) * t651 - m(5) * t676 / 0.2e1;
t548 = t440 * t442;
t504 = t587 * t548;
t244 = (-t494 - t597) * t442 + t504;
t652 = m(5) / 0.2e1;
t535 = (-t244 * t442 + t248) * t651 + t676 * t652;
t41 = t535 - t531;
t499 = qJD(1) + qJD(2);
t697 = t499 * t41;
t448 = cos(qJ(3));
t598 = pkin(3) * t448;
t436 = pkin(2) + t598;
t596 = -qJ(4) - pkin(7);
t503 = -t441 * t436 - t442 * t596;
t216 = t503 - t700;
t600 = sin(qJ(1)) * pkin(1);
t209 = t216 - t600;
t418 = t441 * t596;
t588 = t441 * rSges(6,2);
t217 = t588 - t418 + (t436 + t669) * t442;
t599 = cos(qJ(1)) * pkin(1);
t210 = t217 + t599;
t114 = t209 * t442 + t210 * t441;
t118 = t216 * t442 + t217 * t441;
t549 = t440 * t441;
t555 = t439 * t441;
t463 = rSges(5,1) * t549 - rSges(5,2) * t555 - t442 * rSges(5,3);
t258 = -t463 + t503;
t253 = t258 - t600;
t491 = -rSges(5,2) * t554 + t441 * rSges(5,3);
t589 = rSges(5,1) * t440;
t259 = -t418 + (t436 + t589) * t442 + t491;
t254 = t259 + t599;
t145 = t253 * t442 + t254 * t441;
t162 = t258 * t442 + t259 * t441;
t594 = (t118 + t114) * t651 + (t162 + t145) * t652;
t595 = (t114 - t118) * t651 + (t145 - t162) * t652;
t6 = t595 - t594;
t696 = t6 * qJD(1);
t380 = -Icges(5,2) * t439 + t421;
t695 = t380 + t698;
t384 = Icges(5,1) * t440 - t580;
t668 = Icges(6,1) * t440 + t420;
t693 = t668 + t384;
t691 = t209 - t216;
t443 = Icges(4,4) * t448;
t406 = -Icges(4,2) * t446 + t443;
t407 = Icges(4,1) * t446 + t443;
t689 = t406 + t407;
t306 = Icges(6,4) * t441 + t442 * t668;
t308 = Icges(5,5) * t441 + t384 * t442;
t688 = -t442 * t694 + t306 + t308;
t396 = Icges(6,5) * t548;
t298 = Icges(6,6) * t441 + Icges(6,3) * t554 + t396;
t304 = Icges(5,6) * t441 + t380 * t442;
t687 = -Icges(6,1) * t554 - t383 * t442 + t298 - t304 + t396;
t686 = -Icges(4,5) * t446 - Icges(4,6) * t448 + (-Icges(5,6) + Icges(6,6)) * t440 + (-Icges(6,4) - Icges(5,5)) * t439;
t653 = m(4) / 0.2e1;
t621 = m(3) * (t599 * (-rSges(3,1) * t441 - rSges(3,2) * t442) + (t442 * rSges(3,1) - t441 * rSges(3,2)) * t600);
t495 = t217 * t554;
t112 = -t216 * t555 + t495;
t681 = t112 * m(6) * qJD(2);
t378 = Icges(6,4) * t440 + Icges(6,6) * t439;
t564 = t378 * t441;
t301 = -Icges(6,2) * t442 + t564;
t284 = t441 * t301;
t376 = Icges(6,3) * t439 + t579;
t297 = -Icges(6,6) * t442 + t376 * t441;
t305 = -Icges(6,4) * t442 + t441 * t668;
t150 = t297 * t554 + t305 * t548 + t284;
t680 = t150 * t442;
t437 = t441 ^ 2;
t438 = t442 ^ 2;
t500 = t437 + t438;
t397 = Icges(5,4) * t555;
t307 = Icges(5,1) * t549 - Icges(5,5) * t442 - t397;
t665 = -Icges(5,2) * t549 - t471 * t441 + t305 + t307 - t397;
t303 = Icges(5,4) * t549 - Icges(5,2) * t555 - Icges(5,6) * t442;
t664 = t698 * t441 - t297 + t303;
t581 = Icges(4,4) * t446;
t405 = Icges(4,2) * t448 + t581;
t408 = Icges(4,1) * t448 - t581;
t678 = (-t405 + t408) * t448 - t689 * t446 + (t693 - t694) * t440 + (t376 - t695) * t439;
t109 = t253 * t674 - t254 * t673;
t113 = t258 * t674 - t259 * t673;
t672 = (t297 * t439 + t305 * t440) * t441;
t402 = t437 * t439;
t403 = t438 * t439;
t256 = 0.2e1 * (t402 / 0.2e1 + t403 / 0.2e1) * m(6);
t671 = t499 * t256;
t73 = -t217 * t209 + t210 * t216;
t104 = -t259 * t253 + t254 * t258;
t435 = t442 * pkin(7);
t590 = rSges(4,1) * t448;
t493 = pkin(2) + t590;
t544 = t441 * t446;
t501 = rSges(4,2) * t544 + t442 * rSges(4,3);
t280 = -t493 * t441 + t435 + t501;
t264 = t280 - t600;
t541 = t442 * t446;
t417 = rSges(4,2) * t541;
t281 = -t417 + t493 * t442 + (rSges(4,3) + pkin(7)) * t441;
t265 = t281 + t599;
t115 = -t281 * t264 + t265 * t280;
t667 = t688 * t441;
t666 = t687 * t441;
t663 = t686 * t442;
t662 = t686 * t441;
t409 = rSges(4,1) * t446 + rSges(4,2) * t448;
t498 = ((-t210 + t217) * t262 + t691 * t699) * t651 + (t109 - t113) * t652 + ((-t265 + t281) * t442 + (t264 - t280) * t441) * t409 * t653;
t361 = t409 * t441;
t362 = t409 * t442;
t141 = t264 * t361 - t265 * t362;
t154 = t280 * t361 - t281 * t362;
t82 = t209 * t699 + t210 * t244;
t83 = t216 * t699 + t217 * t244;
t661 = (t83 + t82) * t651 + (t113 + t109) * t652 + (t154 + t141) * t653;
t453 = t689 * t448 / 0.2e1 + (-t405 / 0.2e1 + t408 / 0.2e1) * t446 + (-t376 / 0.2e1 + t695 / 0.2e1) * t440 + (-t694 / 0.2e1 + t693 / 0.2e1) * t439;
t658 = 0.4e1 * qJD(1);
t656 = 0.4e1 * qJD(2);
t655 = 2 * qJD(3);
t180 = t210 * t554;
t640 = m(6) * (-t691 * t555 + t180 - t495);
t639 = m(6) * (t180 + t495 + (-t209 - t216) * t555);
t525 = t244 * t555 + t701;
t529 = t209 * t548 + t210 * t549;
t635 = m(6) * (t525 + t529);
t634 = m(6) * t73;
t527 = t216 * t548 + t217 * t549;
t633 = m(6) * (t525 + t527);
t479 = t262 * t555 - t701;
t632 = m(6) * (t479 + t529);
t630 = m(6) * (t479 + t527);
t629 = m(6) * t82;
t628 = m(6) * t83;
t627 = -t441 / 0.2e1;
t626 = t441 / 0.2e1;
t625 = -t442 / 0.2e1;
t322 = Icges(4,6) * t441 + t406 * t442;
t324 = Icges(4,5) * t441 + t408 * t442;
t543 = t441 * t448;
t287 = t324 * t543;
t404 = Icges(4,5) * t448 - Icges(4,6) * t446;
t561 = t404 * t442;
t320 = Icges(4,3) * t441 + t561;
t485 = t320 * t442 - t287;
t171 = -t322 * t544 - t485;
t319 = Icges(4,5) * t543 - Icges(4,6) * t544 - Icges(4,3) * t442;
t415 = Icges(4,4) * t544;
t323 = Icges(4,1) * t543 - Icges(4,5) * t442 - t415;
t321 = Icges(4,4) * t543 - Icges(4,2) * t544 - Icges(4,6) * t442;
t566 = t321 * t446;
t102 = -(-(-t323 * t448 + t566) * t441 - t319 * t442) * t442 + t171 * t441;
t540 = t442 * t448;
t519 = -t441 * t319 - t323 * t540;
t172 = -t321 * t541 - t519;
t518 = t441 * t320 + t324 * t540;
t173 = -t322 * t541 + t518;
t103 = -t172 * t442 + t173 * t441;
t563 = t378 * t442;
t302 = Icges(6,2) * t441 + t563;
t151 = t298 * t554 + t441 * t302 + t306 * t548;
t569 = t301 * t442;
t462 = t151 + t569;
t478 = -t298 * t555 + t302 * t442 - t306 * t549;
t25 = (t151 - t462) * t442 + (t150 + t478 - t284) * t441;
t299 = Icges(5,5) * t549 - Icges(5,6) * t555 - Icges(5,3) * t442;
t522 = -t441 * t299 - t307 * t548;
t152 = -t303 * t554 - t522;
t377 = Icges(5,5) * t440 - Icges(5,6) * t439;
t565 = t377 * t442;
t300 = Icges(5,3) * t441 + t565;
t521 = t441 * t300 + t308 * t548;
t153 = -t304 * t554 + t521;
t484 = t304 * t439 - t299;
t269 = t308 * t549;
t486 = t300 * t442 - t269;
t26 = (t442 * t484 + t153 - t521) * t442 + (t441 * t484 + t152 + t486) * t441;
t146 = -t569 + t672;
t27 = -t680 + (t146 + t462 - t672) * t441;
t149 = -t304 * t555 - t486;
t568 = t303 * t439;
t28 = (t149 - t269 + (t300 + t568) * t442 + t522) * t442 + t521 * t441;
t483 = t322 * t446 - t319;
t33 = (t442 * t483 + t173 - t518) * t442 + (t441 * t483 + t172 + t485) * t441;
t34 = (t171 - t287 + (t320 + t566) * t442 + t519) * t442 + t518 * t441;
t88 = -t146 * t442 - t441 * t478;
t89 = -(-(-t307 * t440 + t568) * t441 - t299 * t442) * t442 + t149 * t441;
t90 = t151 * t441 - t680;
t91 = -t152 * t442 + t153 * t441;
t2 = (-t34 / 0.2e1 - t28 / 0.2e1 - t27 / 0.2e1 + t91 / 0.2e1 + t90 / 0.2e1 + t103 / 0.2e1) * t442 + (t102 / 0.2e1 + t89 / 0.2e1 + t88 / 0.2e1 + t26 / 0.2e1 + t25 / 0.2e1 + t33 / 0.2e1) * t441;
t622 = t2 * qJD(3) - qJD(4) * t41;
t619 = m(4) * t115;
t617 = m(4) * t141;
t616 = m(4) * t154;
t615 = m(5) * t104;
t613 = m(5) * t109;
t612 = m(5) * t113;
t611 = m(5) * t145;
t610 = m(5) * t162;
t605 = m(6) * t114;
t604 = m(6) * t118;
t592 = m(6) * qJD(3);
t591 = m(6) * qJD(5);
t42 = t531 + t535;
t586 = t42 * qJD(3);
t585 = t41 * qJD(3) - t256 * qJD(5);
t556 = t439 * t440;
t530 = (-t244 - t262) * t699;
t523 = -t440 * t248 - t262 * t548;
t520 = -t441 * (pkin(2) * t441 - t435 + t503) + t442 * (-t441 * pkin(7) - t418 + (-pkin(2) + t436) * t442);
t509 = t407 * t441 + t321;
t508 = -t407 * t442 - t322;
t507 = -Icges(4,2) * t543 + t323 - t415;
t506 = -t405 * t442 + t324;
t505 = t500 * t556;
t502 = t402 + t403;
t108 = -t209 * t555 + t180;
t496 = m(6) * t108 * qJD(1);
t492 = rSges(5,2) * t439 - t589 - t598;
t482 = t500 * t597;
t480 = -t598 - t669;
t460 = (-t361 * t442 + t362 * t441) * t409;
t458 = t507 * t446 + t509 * t448;
t457 = -t506 * t446 + t508 * t448;
t452 = t453 + t661;
t451 = -t453 + (t626 + t627) * (t321 * t448 + t323 * t446 + (t297 + t303) * t440 + (-t305 + t307) * t439);
t450 = t42 * qJD(4) + ((t34 + t28 + t27) * t442 / 0.2e1 + (t102 + t89 + t88 + t33 + t26 + t25) * t627 + (t508 * t446 + t506 * t448 + t564 + (t377 + t404) * t441 + t688 * t440 + t687 * t439 + t678 * t442) * t626 + (-t439 * t664 + t440 * t665 + t678 * t441 - t509 * t446 + t507 * t448 + t103 - t561 - t563 - t565 + t90 + t91) * t625) * qJD(3);
t411 = -rSges(4,2) * t446 + t590;
t314 = t492 * t442;
t312 = t492 * t441;
t266 = t505 - t556;
t263 = t480 * t442;
t261 = t480 * t441;
t249 = t256 * qJD(4);
t129 = -t482 + (-t442 * t494 + t504) * t442 + t679 * t437;
t105 = (t669 * t442 + t588) * t442 + t520 + t700 * t441;
t81 = t604 + t610;
t79 = t630 / 0.2e1;
t77 = t605 + t611;
t74 = t632 / 0.2e1;
t72 = t633 / 0.2e1;
t70 = t635 / 0.2e1;
t54 = t639 / 0.2e1;
t53 = t640 / 0.2e1;
t52 = t105 * t502 + t523;
t20 = t615 + t619 + t621 + t634;
t19 = t453 + t612 + t616 + t628;
t18 = t453 + t613 + t617 + t629;
t17 = t54 - t640 / 0.2e1;
t16 = t54 + t53;
t15 = t53 - t639 / 0.2e1;
t14 = t79 - t633 / 0.2e1;
t13 = t79 + t72;
t12 = t72 - t630 / 0.2e1;
t11 = t74 - t635 / 0.2e1;
t10 = t74 + t70;
t9 = t70 - t632 / 0.2e1;
t7 = t594 + t595;
t5 = t452 + t498;
t4 = t452 - t498;
t3 = t451 + t498 - t661;
t1 = [t20 * qJD(2) + t18 * qJD(3) + t77 * qJD(4) + t108 * t591, t20 * qJD(1) + t5 * qJD(3) + t7 * qJD(4) + t16 * qJD(5) + 0.2e1 * (t621 / 0.2e1 + t104 * t652 + t115 * t653 + t73 * t651) * qJD(2), t18 * qJD(1) + t5 * qJD(2) + t10 * qJD(5) + (((-t264 * t442 - t265 * t441) * t411 + t460) * t653 + (t253 * t314 + t254 * t312) * t652 + (t209 * t263 + t210 * t261 + t530) * t651) * t655 + t450, qJD(1) * t77 + qJD(2) * t7 + t586, t16 * qJD(2) + t10 * qJD(3) + t496; t4 * qJD(3) - t6 * qJD(4) + t17 * qJD(5) + (-t634 / 0.4e1 - t615 / 0.4e1 - t619 / 0.4e1 - t621 / 0.4e1) * t658, t19 * qJD(3) + t81 * qJD(4) + t112 * t591, t4 * qJD(1) + t19 * qJD(2) + t13 * qJD(5) + (((-t280 * t442 - t281 * t441) * t411 + t460) * t653 + (t258 * t314 + t259 * t312) * t652 + (t216 * t263 + t217 * t261 + t530) * t651) * t655 + t450, qJD(2) * t81 + t586 - t696, t17 * qJD(1) + t13 * qJD(3) + t681; t451 * qJD(1) + t3 * qJD(2) + t11 * qJD(5) + (-t617 / 0.4e1 - t613 / 0.4e1 - t629 / 0.4e1) * t658 + t622, t3 * qJD(1) + t451 * qJD(2) + t14 * qJD(5) + (-t616 / 0.4e1 - t612 / 0.4e1 - t628 / 0.4e1) * t656 + t622, t499 * t2 + t52 * t591 + (m(5) * (-t673 * t314 - t674 * t312 + (t441 * t463 + t442 * (rSges(5,1) * t548 + t491) + t520) * (-t387 * t500 - t482)) + m(6) * (t105 * t129 - t261 * t699 - t262 * t263) + m(4) * ((t441 * (rSges(4,1) * t543 - t501) + t442 * (rSges(4,1) * t540 + t441 * rSges(4,3) - t417)) * (-t361 * t441 - t362 * t442) + t500 * t411 * t409) + ((t458 * t442 + (t457 - t662) * t441 + (t442 * t664 + t666) * t440 + (t442 * t665 - t667) * t439) * t442 + t663 * t437) * t626 + ((t457 * t441 + t666 * t440 - t667 * t439 + (t439 * t665 + t440 * t664 + t458 - t663) * t442) * t441 + t662 * t438) * t625) * qJD(3), -t697, t11 * qJD(1) + t14 * qJD(2) + t52 * t592 + (-t440 * t502 - t266 + t505) * t591; t6 * qJD(2) + (-t605 / 0.4e1 - t611 / 0.4e1) * t658 + t585, t696 + (-t604 / 0.4e1 - t610 / 0.4e1) * t656 + t585, ((-t312 * t442 + t314 * t441) * t652 + (-t261 * t442 + t263 * t441) * t651) * t655 + t697, 0, -t671; t15 * qJD(2) + t9 * qJD(3) + t249 - t496, t15 * qJD(1) + t12 * qJD(3) + t249 - t681, t9 * qJD(1) + t12 * qJD(2) + (-t129 * t440 + (t261 * t441 + t263 * t442 + t105) * t439 - t52 + t523) * t592 + t266 * t591, t671, t266 * t592;];
Cq = t1;
