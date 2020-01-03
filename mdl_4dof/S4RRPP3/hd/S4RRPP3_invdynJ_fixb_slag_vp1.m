% Calculate vector of inverse dynamics joint torques for
% S4RRPP3
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPP3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP3_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP3_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:23
% EndTime: 2019-12-31 16:57:50
% DurationCPUTime: 23.82s
% Computational Cost: add. (7566->586), mult. (11719->727), div. (0->0), fcn. (9163->6), ass. (0->311)
t679 = Icges(3,3) + Icges(4,3);
t321 = qJ(2) + pkin(6);
t296 = sin(t321);
t297 = cos(t321);
t325 = sin(qJ(2));
t327 = cos(qJ(2));
t678 = Icges(3,5) * t327 + Icges(4,5) * t297 - Icges(3,6) * t325 - Icges(4,6) * t296;
t214 = Icges(5,4) * t297 + Icges(5,6) * t296;
t326 = sin(qJ(1));
t328 = cos(qJ(1));
t150 = Icges(5,2) * t326 + t214 * t328;
t659 = t679 * t326 + t678 * t328;
t677 = t150 + t659;
t539 = Icges(4,4) * t296;
t215 = Icges(4,2) * t297 + t539;
t282 = Icges(5,5) * t296;
t674 = Icges(5,3) * t297 + t215 - t282;
t535 = Icges(5,5) * t297;
t217 = Icges(5,1) * t296 - t535;
t283 = Icges(4,4) * t297;
t676 = Icges(4,1) * t296 + t217 + t283;
t675 = t679 * t328;
t505 = t326 * t327;
t507 = t325 * t326;
t509 = t297 * t326;
t511 = t296 * t326;
t650 = -Icges(3,5) * t505 - Icges(4,5) * t509 + Icges(3,6) * t507 + Icges(4,6) * t511 + t675;
t390 = Icges(5,1) * t297 + t282;
t153 = -Icges(5,4) * t328 + t326 * t390;
t259 = Icges(4,4) * t511;
t536 = Icges(4,5) * t328;
t155 = Icges(4,1) * t509 - t259 - t536;
t669 = -t153 - t155;
t154 = Icges(5,4) * t326 + t328 * t390;
t220 = Icges(4,1) * t297 - t539;
t156 = Icges(4,5) * t326 + t220 * t328;
t668 = t154 + t156;
t530 = Icges(4,6) * t328;
t151 = Icges(4,4) * t509 - Icges(4,2) * t511 - t530;
t531 = Icges(3,6) * t328;
t166 = Icges(3,4) * t505 - Icges(3,2) * t507 - t531;
t673 = t151 * t296 + t166 * t325;
t540 = Icges(3,4) * t325;
t247 = Icges(3,1) * t327 - t540;
t169 = Icges(3,5) * t326 + t247 * t328;
t672 = -t156 * t509 - t169 * t505;
t210 = Icges(5,3) * t296 + t535;
t145 = -Icges(5,6) * t328 + t210 * t326;
t671 = t145 - t151;
t508 = t297 * t328;
t258 = Icges(5,5) * t508;
t510 = t296 * t328;
t529 = Icges(5,6) * t326;
t146 = Icges(5,3) * t510 + t258 + t529;
t388 = -Icges(4,2) * t296 + t283;
t152 = Icges(4,6) * t326 + t328 * t388;
t670 = t146 - t152;
t667 = t210 - t388;
t664 = t220 + t390;
t663 = t674 * qJD(2);
t662 = t676 * qJD(2);
t278 = Icges(3,4) * t507;
t537 = Icges(3,5) * t328;
t168 = Icges(3,1) * t505 - t278 - t537;
t661 = t155 * t297 + t168 * t327 - t673;
t660 = Icges(3,5) * t325 + Icges(3,6) * t327 + (Icges(4,6) - Icges(5,6)) * t297 + (Icges(5,4) + Icges(4,5)) * t296;
t409 = -t146 * t511 + t150 * t328 - t154 * t509;
t310 = Icges(3,4) * t327;
t389 = -Icges(3,2) * t325 + t310;
t167 = Icges(3,6) * t326 + t328 * t389;
t658 = -t659 * t328 - t672;
t621 = -t152 * t511 - t167 * t507 + t658;
t604 = -t409 + t621;
t657 = t152 * t296 + t167 * t325;
t504 = t327 * t328;
t656 = t146 * t510 + t169 * t504 + t677 * t326 + t668 * t508;
t149 = -Icges(5,2) * t328 + t214 * t326;
t136 = t326 * t149;
t655 = -t145 * t510 - t168 * t504 + t650 * t326 + t669 * t508 - t136;
t244 = Icges(3,2) * t327 + t540;
t246 = Icges(3,1) * t325 + t310;
t643 = t244 * t325 - t246 * t327 + t674 * t296 - t297 * t676;
t654 = t663 * t328 + (t326 * t388 - t145 - t530) * qJD(1);
t653 = t663 * t326 + (t210 * t328 - t152 + t529) * qJD(1);
t652 = -t662 * t328 + (-t220 * t326 - t153 + t536) * qJD(1);
t651 = -t668 * qJD(1) + t662 * t326;
t649 = t667 * qJD(2);
t648 = t664 * qJD(2);
t523 = t149 * t328;
t384 = t145 * t296 + t153 * t297;
t594 = t326 * t384;
t48 = -t523 + t594;
t611 = t661 * t326 + t650 * t328 + t48;
t506 = t325 * t328;
t610 = -t151 * t510 - t166 * t506 - t655;
t609 = -t152 * t510 - t167 * t506 + t656;
t646 = -t167 * t327 - t169 * t325 - t668 * t296 + t670 * t297;
t605 = -t166 * t327 - t168 * t325 + t669 * t296 + t671 * t297;
t645 = t214 + t678;
t644 = t660 * qJD(2);
t642 = -t244 * t327 - t246 * t325 - t296 * t676 - t297 * t674;
t641 = t146 * t296 + t169 * t327 + t668 * t297 - t657;
t640 = -t384 - t661;
t597 = t660 * t328;
t596 = t660 * t326;
t608 = -t643 * t326 - t597;
t607 = -t643 * t328 + t596;
t639 = t677 * qJD(1);
t638 = t328 ^ 2;
t637 = -t671 * t328 + (-Icges(5,1) * t510 + t217 * t328 + t258 + t670) * t326;
t636 = t676 - t667;
t635 = t664 - t674;
t634 = (Icges(4,2) * t509 + t259 + t669) * t328 + (-t215 * t328 + t668) * t326;
t229 = t389 * qJD(2);
t230 = t247 * qJD(2);
t633 = t660 * qJD(1) + t642 * qJD(2) - t229 * t325 + t230 * t327 + t649 * t296 + t648 * t297;
t362 = qJD(2) * t244;
t107 = -t328 * t362 + (-t326 * t389 + t531) * qJD(1);
t365 = qJD(2) * t246;
t109 = -t328 * t365 + (-t247 * t326 + t537) * qJD(1);
t632 = t646 * qJD(2) - t107 * t325 + t109 * t327 + t654 * t296 + t652 * t297 + t639;
t108 = qJD(1) * t167 - t326 * t362;
t110 = qJD(1) * t169 - t326 * t365;
t582 = qJD(1) * t149;
t631 = t650 * qJD(1) - t605 * qJD(2) + t108 * t325 - t110 * t327 - t653 * t296 + t651 * t297 - t582;
t630 = t609 * t326 - t610 * t328;
t629 = t604 * t326 - t611 * t328;
t628 = t643 * qJD(1) + t645 * qJD(2);
t627 = t640 * qJD(1) - t644 * t326 + t639;
t626 = -t582 - t644 * t328 + (-t326 * t678 - t641 + t675) * qJD(1);
t286 = t297 * rSges(4,1);
t586 = -rSges(4,2) * t296 + t286;
t199 = t586 * qJD(2);
t223 = rSges(4,1) * t296 + rSges(4,2) * t297;
t449 = qJD(1) * qJD(2);
t236 = -qJDD(2) * t328 + t326 * t449;
t448 = qJD(1) * qJD(3);
t563 = pkin(2) * t325;
t439 = qJDD(3) * t326 + t236 * t563 + t328 * t448;
t159 = rSges(4,1) * t509 - rSges(4,2) * t511 - t328 * rSges(4,3);
t319 = t328 * pkin(5);
t266 = pkin(1) * t326 - t319;
t324 = -qJ(3) - pkin(5);
t290 = t328 * t324;
t318 = t327 * pkin(2);
t291 = t318 + pkin(1);
t467 = -t326 * t291 - t290;
t143 = t266 + t467;
t492 = t143 - t266;
t441 = -t159 + t492;
t503 = t327 * qJD(2) ^ 2;
t445 = pkin(2) * t503;
t317 = t326 * pkin(5);
t456 = qJD(1) * t326;
t447 = pkin(2) * t507;
t465 = qJD(2) * t447 + qJD(3) * t328;
t438 = t324 * t456 + t465;
t558 = pkin(1) - t291;
t104 = (-t328 * t558 - t317) * qJD(1) - t438;
t267 = t328 * pkin(1) + t317;
t234 = t267 * qJD(1);
t500 = -t104 - t234;
t184 = t223 * t326;
t311 = t326 * rSges(4,3);
t89 = -qJD(2) * t184 + (t328 * t586 + t311) * qJD(1);
t20 = t223 * t236 + (-qJD(2) * t199 - t445) * t328 + t441 * qJDD(1) + (-t89 + t500) * qJD(1) + t439;
t625 = t20 - g(1);
t624 = rSges(3,2) * t325;
t620 = rSges(5,3) + qJ(4);
t623 = t607 * qJD(1);
t622 = t608 * qJD(1);
t564 = rSges(5,1) + pkin(3);
t619 = t650 + t657;
t587 = t297 * rSges(5,1) + t296 * rSges(5,3);
t618 = t297 * pkin(3) + t296 * qJ(4) + t587;
t617 = t629 * qJD(2) + t622;
t616 = t630 * qJD(2) + t623;
t615 = t628 * t326 + t633 * t328;
t614 = t633 * t326 - t628 * t328;
t613 = t640 * qJD(2) - t108 * t327 - t110 * t325 + t651 * t296 + t653 * t297;
t612 = t641 * qJD(2) + t107 * t327 + t109 * t325 + t652 * t296 - t654 * t297;
t349 = t166 * t328 - t167 * t326;
t571 = t326 * (-t244 * t328 + t169) - t328 * (-Icges(3,2) * t505 + t168 - t278);
t603 = -t634 * t296 + t637 * t297 - t325 * t571 + t349 * t327;
t469 = t246 + t389;
t470 = -t244 + t247;
t602 = (-t636 * t296 + t635 * t297 - t325 * t469 + t327 * t470) * qJD(1);
t601 = t523 + t656;
t600 = t627 * t638 + (t632 * t326 + (-t626 + t631) * t328) * t326;
t599 = t631 * t638 + (t626 * t326 + (-t627 + t632) * t328) * t326;
t598 = t645 * qJD(1);
t595 = t296 * t564;
t355 = -t223 - t563;
t593 = t328 * t355;
t238 = qJD(1) * t266;
t592 = qJD(1) * t143 - t238;
t591 = t586 + t318;
t453 = qJD(2) * t328;
t434 = t297 * t453;
t455 = qJD(1) * t328;
t590 = rSges(5,2) * t455 + t434 * t620;
t589 = t620 * t509;
t588 = t620 * t508;
t222 = rSges(5,1) * t296 - rSges(5,3) * t297;
t474 = pkin(3) * t296 - qJ(4) * t297 + t222;
t414 = -t474 - t563;
t452 = qJD(4) * t296;
t248 = t328 * t452;
t299 = qJD(3) * t326;
t468 = t248 + t299;
t584 = t414 * t453 + t468;
t561 = g(2) * t326;
t583 = -g(1) * t328 - t561;
t581 = t318 + t618;
t235 = qJDD(2) * t326 + t328 * t449;
t568 = t235 / 0.2e1;
t567 = t236 / 0.2e1;
t566 = t326 / 0.2e1;
t565 = -t328 / 0.2e1;
t354 = -t296 * t453 - t297 * t456;
t437 = t296 * t456;
t557 = t354 * t564 - t620 * t437 + t248 + t590;
t313 = t326 * rSges(5,2);
t454 = qJD(2) * t326;
t556 = (pkin(3) * t455 + qJ(4) * t454) * t297 + (qJ(4) * t455 + (-pkin(3) * qJD(2) + qJD(4)) * t326) * t296 - t222 * t454 + (t328 * t587 + t313) * qJD(1);
t555 = rSges(3,1) * t327;
t249 = rSges(3,1) * t325 + rSges(3,2) * t327;
t208 = t249 * t328;
t312 = t326 * rSges(3,3);
t191 = rSges(3,1) * t504 - rSges(3,2) * t506 + t312;
t133 = t191 + t267;
t95 = qJD(1) * t133 - t249 * t454;
t553 = t208 * t95;
t451 = qJD(4) * t297;
t484 = t508 * t564 + t620 * t510 + t313;
t316 = t328 * rSges(5,2);
t485 = t326 * t618 - t316;
t417 = t328 * t291 - t324 * t326;
t144 = t417 - t267;
t496 = -t143 * t454 + t144 * t453;
t32 = -t451 + (t326 * t485 + t328 * t484) * qJD(2) + t496;
t552 = t32 * t296;
t435 = t249 * t453;
t464 = rSges(3,2) * t507 + t328 * rSges(3,3);
t190 = rSges(3,1) * t505 - t464;
t479 = -t190 - t266;
t94 = qJD(1) * t479 - t435;
t550 = t326 * t94;
t548 = t328 * t94;
t353 = qJD(2) * t593 + t299;
t46 = qJD(1) * t441 + t353;
t546 = t46 * t223;
t495 = -t326 * t143 + t328 * t144;
t161 = rSges(4,1) * t508 - rSges(4,2) * t510 + t311;
t491 = -t144 - t161;
t486 = -qJD(2) * t618 + t451;
t481 = -t511 * t564 + t589;
t480 = -t510 * t564 + t588;
t295 = pkin(5) * t455;
t473 = qJD(1) * (-pkin(1) * t456 + t295) + qJDD(1) * t267;
t471 = rSges(4,2) * t437 + rSges(4,3) * t455;
t466 = rSges(3,3) * t455 + t456 * t624;
t463 = t326 ^ 2 + t638;
t446 = pkin(2) * t506;
t433 = t325 * t453;
t103 = -pkin(2) * t433 - t295 + t299 + (t326 * t558 - t290) * qJD(1);
t444 = t328 * t103 + t326 * t104 - t143 * t455;
t443 = t103 * t453 + t104 * t454 - t235 * t143;
t442 = qJD(2) * t318;
t440 = -t144 - t484;
t432 = t327 * t453;
t430 = t326 * t452;
t429 = -pkin(1) - t555;
t426 = -t454 / 0.2e1;
t425 = t454 / 0.2e1;
t424 = -t453 / 0.2e1;
t423 = t453 / 0.2e1;
t416 = -t485 + t492;
t415 = t463 * t563;
t252 = rSges(2,1) * t328 - rSges(2,2) * t326;
t250 = rSges(2,1) * t326 + rSges(2,2) * t328;
t251 = t555 - t624;
t35 = qJD(1) * t416 + t584;
t36 = (-qJD(2) * t474 + t452) * t326 + (t267 - t440) * qJD(1) - t465;
t401 = t326 * t36 + t328 * t35;
t87 = rSges(4,1) * t354 - rSges(4,2) * t434 + t471;
t394 = t326 * t89 + t328 * t87;
t393 = -t326 * t95 - t548;
t111 = -rSges(3,2) * t432 + (-t327 * t456 - t433) * rSges(3,1) + t466;
t207 = t249 * t326;
t112 = -qJD(2) * t207 + (t251 * t328 + t312) * qJD(1);
t386 = t111 * t328 + t112 * t326;
t379 = t159 * t326 + t161 * t328;
t376 = t190 * t326 + t191 * t328;
t369 = -t442 + t486;
t367 = t401 * t297;
t352 = qJD(1) * t103 + qJDD(1) * t144 - qJDD(3) * t328 + t326 * t448 + t473;
t344 = -t291 - t618;
t343 = -t445 + qJDD(4) * t296 + (t451 + t486) * qJD(2);
t231 = t251 * qJD(2);
t188 = t223 * t328;
t90 = t376 * qJD(2);
t47 = -t223 * t454 + (t267 - t491) * qJD(1) - t465;
t43 = qJD(2) * t379 + t496;
t42 = qJD(1) * t111 + qJDD(1) * t191 - t231 * t454 - t235 * t249 + t473;
t41 = -t231 * t453 + t236 * t249 + t479 * qJDD(1) + (-t112 - t234) * qJD(1);
t21 = -t199 * t454 + qJD(1) * t87 + qJDD(1) * t161 - t223 * t235 + (-t235 * t325 - t326 * t503) * pkin(2) + t352;
t3 = t484 * qJDD(1) + t414 * t235 + (t248 + t557) * qJD(1) + t343 * t326 + t352;
t2 = t474 * t236 + t416 * qJDD(1) + t343 * t328 + (-t430 + t500 - t556) * qJD(1) + t439;
t1 = -qJDD(4) * t297 + t485 * t235 + t440 * t236 + (t326 * t556 + t328 * t557 + t452) * qJD(2) + t443;
t4 = [-m(2) * (-g(1) * t250 + g(2) * t252) + (((t48 - t594 + t601) * t326 + ((t659 + t673) * t328 + t621 + t655 + t672) * t328) * qJD(2) + t623) * t423 + (-t643 * qJD(2) + t229 * t327 + t230 * t325 + t648 * t296 - t649 * t297) * qJD(1) + (t35 * (-t430 + t438) + t36 * (t468 + t590) + ((-t563 - t595) * t328 * t36 + (-t297 * t620 + t595) * t326 * t35) * qJD(2) + ((-t36 * t324 + t344 * t35) * t328 + (-t35 * rSges(5,2) + t344 * t36) * t326) * qJD(1) - (-qJD(1) * t485 - t35 + t584 + t592) * t36 + (t3 - g(2)) * (t417 + t484) + (t2 - g(1)) * (t316 + (-t296 * t620 - t297 * t564) * t326 + t467)) * m(5) + (t46 * t438 + t47 * (t299 + t471) + (t326 * t546 + t47 * t593) * qJD(2) + ((-t46 * rSges(4,3) + t47 * (-t291 - t286)) * t326 + (t46 * (-t291 - t586) - t47 * t324) * t328) * qJD(1) - (-qJD(1) * t159 + t353 - t46 + t592) * t47 + (t21 - g(2)) * (t161 + t417) + t625 * (-t159 + t467)) * m(4) + (-(-qJD(1) * t190 - t238 - t435 - t94) * t95 + t95 * (t295 + t466) + (t249 * t550 - t553) * qJD(2) + ((-pkin(1) - t251) * t548 + (t94 * (-rSges(3,3) - pkin(5)) + t95 * t429) * t326) * qJD(1) + (t42 - g(2)) * t133 + (t41 - g(1)) * (t429 * t326 + t319 + t464)) * m(3) + (m(2) * (t250 ^ 2 + t252 ^ 2) + Icges(2,3) - t642) * qJDD(1) + (-t646 + t607) * t568 + (-t605 + t608) * t567 + (t612 + t615) * t425 + (((t328 * t619 - t601 + t609) * t328 + (t326 * t619 - t136 + t409 + t610 - t658) * t326) * qJD(2) + t617 - t622) * t426 + (-t613 + t614 + t616) * t424; t630 * t568 + t629 * t567 + (qJD(1) * t615 + qJD(2) * t599 + qJDD(1) * t607 + t235 * t609 + t236 * t610) * t566 + (qJD(1) * t614 + qJD(2) * t600 + qJDD(1) * t608 + t235 * t604 + t236 * t611) * t565 - ((t637 * t296 + t634 * t297 + t349 * t325 + t327 * t571) * qJD(2) + (t635 * t296 + t636 * t297 + t325 * t470 + t327 * t469) * qJD(1)) * qJD(1) / 0.2e1 + (t613 * t328 + t612 * t326 + (-t326 * t605 - t328 * t646) * qJD(1)) * qJD(1) / 0.2e1 + (-t326 * t646 + t328 * t605) * qJDD(1) / 0.2e1 + t617 * t456 / 0.2e1 + t616 * t455 / 0.2e1 + ((-t454 * t597 + t598) * t326 + ((t326 * t596 + t603) * qJD(2) + t602) * t328) * t426 + ((t326 * t610 + t328 * t609) * qJD(1) + t599) * t425 + ((t326 * t611 + t328 * t604) * qJD(1) + t600) * t424 + ((-t453 * t596 - t598) * t328 + ((t328 * t597 + t603) * qJD(2) + t602) * t326) * t423 + (-g(3) * t591 - t355 * t561 + t46 * (-pkin(2) * t432 - t199 * t328) + (qJD(2) * t394 + t159 * t235 + t236 * t491 + t443) * (t379 + t495) + t43 * (t394 + t444) + (t43 * t159 + t355 * t47) * t455 + (t21 * t355 + t47 * (-t199 - t442) + (t43 * t491 + t546) * qJD(1)) * t326 - (t46 * t184 + t47 * (-t188 - t446)) * qJD(1) - (-t43 * t415 + (-t43 * t188 - t46 * t591) * t328 + (-t43 * t184 - t47 * t591) * t326) * qJD(2) + t625 * t593) * m(4) + (g(1) * t208 + g(2) * t207 - g(3) * t251 - (t207 * t94 - t553) * qJD(1) - (t90 * (-t207 * t326 - t208 * t328) + t393 * t251) * qJD(2) + (qJD(2) * t386 + t190 * t235 - t191 * t236) * t376 + t90 * ((t190 * t328 - t191 * t326) * qJD(1) + t386) + t393 * t231 + (-t42 * t326 - t41 * t328 + (-t328 * t95 + t550) * qJD(1)) * t249) * m(3) + (-g(1) * (-t446 + t588) - g(2) * (-t447 + t589) - g(3) * t581 - t583 * t595 - (t367 + t552) * qJD(4) - (-t35 * t481 + t36 * (-t446 + t480)) * qJD(1) - (-t32 * t415 + (t32 * t480 - t35 * t581) * t328 + (t32 * t481 - t36 * t581) * t326) * qJD(2) + t1 * t495 + t32 * t444 + (t2 * t414 + t35 * t369 + t1 * t484 + t32 * t557 + (t32 * t485 + t36 * t414) * qJD(1)) * t328 + (t3 * t414 + t36 * t369 + t1 * t485 + t32 * t556 + (t32 * t440 + t35 * t474) * qJD(1)) * t326) * m(5); (-m(4) - m(5)) * (g(1) * t326 - g(2) * t328) + 0.2e1 * (t2 * t566 + t3 * t565) * m(5) + 0.2e1 * (t20 * t566 + t21 * t565) * m(4); (-(t463 * t552 + t367) * qJD(2) + (qJD(2) * t401 + g(3) - t1) * t297 + (qJD(2) * t32 + t2 * t328 + t3 * t326 + t583) * t296) * m(5);];
tau = t4;
