% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRPR16_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR16_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR16_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:46
% EndTime: 2019-12-31 18:39:06
% DurationCPUTime: 14.94s
% Computational Cost: add. (18281->649), mult. (45261->926), div. (0->0), fcn. (48368->6), ass. (0->384)
t410 = cos(qJ(3));
t411 = cos(qJ(1));
t521 = t410 * t411;
t385 = qJ(4) * t521;
t397 = t411 * qJ(2);
t407 = sin(qJ(3));
t408 = sin(qJ(1));
t617 = pkin(1) + pkin(6);
t490 = rSges(5,1) + t617;
t563 = rSges(5,3) * t410;
t229 = -t385 + t397 + (-t563 + (-rSges(5,2) + pkin(3)) * t407) * t411 - t490 * t408;
t525 = t408 * t410;
t383 = qJ(4) * t525;
t528 = t407 * t408;
t498 = -rSges(5,2) * t528 - rSges(5,3) * t525;
t567 = pkin(3) * t407;
t230 = -t383 + (qJ(2) + t567) * t408 + t490 * t411 + t498;
t527 = t407 * t411;
t469 = t229 * t528 - t230 * t527;
t406 = sin(qJ(5));
t409 = cos(qJ(5));
t332 = -t406 * t411 - t409 * t525;
t333 = -t406 * t525 + t409 * t411;
t463 = t333 * rSges(6,1) + t332 * rSges(6,2);
t491 = pkin(4) + t617;
t594 = rSges(6,3) + pkin(7);
t186 = -t383 + t491 * t411 + (qJ(2) + (pkin(3) + t594) * t407) * t408 + t463;
t473 = pkin(3) * t527 - t385;
t330 = t406 * t408 - t409 * t521;
t526 = t408 * t409;
t331 = t406 * t521 + t526;
t392 = rSges(6,3) * t527;
t251 = t331 * rSges(6,1) - t330 * rSges(6,2) - t392;
t395 = pkin(7) * t527;
t658 = t251 - t395;
t634 = -t491 * t408 + t397 + t473 - t658;
t470 = -t186 * t527 + t528 * t634;
t619 = m(6) / 0.2e1;
t620 = m(5) / 0.2e1;
t376 = pkin(3) * t410 + qJ(4) * t407;
t356 = t408 * t376;
t461 = rSges(5,2) * t410 - rSges(5,3) * t407;
t282 = -t408 * t461 + t356;
t284 = (t376 - t461) * t411;
t636 = -t282 * t411 + t284 * t408;
t393 = pkin(7) * t525;
t564 = rSges(6,1) * t406;
t462 = rSges(6,2) * t409 + t564;
t421 = t462 * t407;
t317 = rSges(6,3) * t410 + t421;
t532 = t317 * t408;
t239 = t356 + t393 + t532;
t241 = (pkin(7) * t410 + t317 + t376) * t411;
t637 = -t239 * t411 + t241 * t408;
t556 = (t636 * t410 + t469) * t620 + (t637 * t410 + t470) * t619;
t351 = pkin(3) * t525 + qJ(4) * t528;
t271 = -rSges(5,2) * t525 + rSges(5,3) * t528 + t351;
t497 = -pkin(3) * t521 - qJ(4) * t527;
t272 = -t461 * t411 - t497;
t438 = t271 * t411 - t408 * t272;
t279 = t407 * rSges(6,2) * t526 + rSges(6,3) * t525 + t528 * t564;
t509 = -t279 - t351;
t227 = t393 - t509;
t228 = (t594 * t410 + t421) * t411 - t497;
t443 = t227 * t411 - t408 * t228;
t557 = (t438 * t410 + t469) * t620 + (t443 * t410 + t470) * t619;
t655 = t556 - t557;
t678 = t655 * qJD(1);
t452 = Icges(6,5) * t406 + Icges(6,6) * t409;
t299 = Icges(6,3) * t410 + t452 * t407;
t547 = Icges(6,4) * t406;
t455 = Icges(6,2) * t409 + t547;
t303 = Icges(6,6) * t410 + t455 * t407;
t546 = Icges(6,4) * t409;
t386 = t407 * t546;
t529 = t406 * t407;
t545 = Icges(6,5) * t410;
t307 = Icges(6,1) * t529 + t386 + t545;
t169 = t299 * t528 + t303 * t332 + t307 * t333;
t243 = -Icges(6,5) * t331 + Icges(6,6) * t330 + Icges(6,3) * t527;
t326 = Icges(6,4) * t331;
t245 = -Icges(6,2) * t330 - Icges(6,6) * t527 + t326;
t325 = Icges(6,4) * t330;
t249 = -Icges(6,1) * t331 + Icges(6,5) * t527 + t325;
t670 = -t332 * t245 + t333 * t249;
t124 = -t243 * t528 - t670;
t244 = Icges(6,5) * t333 + Icges(6,6) * t332 + Icges(6,3) * t528;
t548 = Icges(6,4) * t333;
t247 = Icges(6,2) * t332 + Icges(6,6) * t528 + t548;
t327 = Icges(6,4) * t332;
t250 = Icges(6,1) * t333 + Icges(6,5) * t528 + t327;
t125 = t244 * t528 + t332 * t247 + t333 * t250;
t447 = -t124 * t411 + t125 * t408;
t677 = t169 * t410 + t447 * t407;
t378 = rSges(4,1) * t410 - rSges(4,2) * t407;
t352 = t378 * t408;
t354 = t378 * t411;
t667 = -m(6) / 0.2e1;
t668 = -m(5) / 0.2e1;
t483 = t443 * t667 + t438 * t668 + m(4) * (-t352 * t411 + t408 * t354) / 0.2e1;
t518 = t636 * t668 + t637 * t667;
t42 = t518 - t483;
t676 = qJD(1) * t42;
t674 = (-Icges(5,4) + Icges(4,5)) * t410 + (Icges(5,5) - Icges(4,6)) * t407;
t673 = t124 * t408 + t125 * t411;
t419 = rSges(6,3) * t528 + t463;
t204 = t317 * t528 - t410 * t419;
t671 = -t410 * t251 - t317 * t527;
t148 = -t204 * t408 + t411 * t671;
t550 = Icges(4,4) * t407;
t457 = Icges(4,2) * t410 + t550;
t306 = -Icges(4,6) * t408 + t457 * t411;
t549 = Icges(4,4) * t410;
t459 = Icges(4,1) * t407 + t549;
t310 = -Icges(4,5) * t408 + t459 * t411;
t399 = Icges(5,6) * t407;
t541 = Icges(5,3) * t410;
t449 = t399 + t541;
t311 = Icges(5,5) * t408 + t449 * t411;
t543 = Icges(5,6) * t410;
t362 = Icges(5,2) * t407 + t543;
t313 = Icges(5,4) * t408 + t362 * t411;
t433 = t311 * t410 + t313 * t407;
t672 = (-t306 * t410 - t310 * t407 - t433) * t411;
t441 = t245 * t409 - t249 * t406;
t536 = t243 * t410;
t149 = t441 * t407 - t536;
t305 = Icges(4,6) * t411 + t457 * t408;
t387 = Icges(4,4) * t525;
t309 = Icges(4,1) * t528 + Icges(4,5) * t411 + t387;
t312 = Icges(5,5) * t411 - t449 * t408;
t314 = Icges(5,4) * t411 - t362 * t408;
t450 = -Icges(5,3) * t407 + t543;
t334 = t450 * t408;
t335 = t450 * t411;
t544 = Icges(5,2) * t410;
t451 = -t399 + t544;
t336 = t451 * t408;
t337 = t451 * t411;
t344 = -Icges(4,2) * t528 + t387;
t367 = -Icges(4,2) * t407 + t549;
t345 = t367 * t411;
t369 = Icges(4,1) * t410 - t550;
t347 = t369 * t408;
t348 = t369 * t411;
t669 = -t407 * ((-t305 + t347 + t312 + t336) * t411 + (t306 - t348 + t311 - t337) * t408) + t410 * ((-t309 - t344 + t314 - t334) * t411 + (t310 + t345 + t313 + t335) * t408);
t122 = t243 * t527 - t330 * t245 - t331 * t249;
t404 = t408 ^ 2;
t583 = m(5) * (t229 * t411 + t230 * t408);
t168 = t299 * t527 + t330 * t303 - t331 * t307;
t663 = t168 * t410;
t661 = t317 * t411;
t435 = -t305 * t410 - t309 * t407;
t660 = t435 * t411;
t365 = Icges(5,4) * t407 + Icges(5,5) * t410;
t315 = Icges(5,1) * t408 + t365 * t411;
t296 = t411 * t315;
t453 = Icges(4,5) * t407 + Icges(4,6) * t410;
t302 = -Icges(4,3) * t408 + t453 * t411;
t659 = -t411 * t302 - t306 * t525 - t310 * t528 - t433 * t408 + t296;
t657 = t674 * t408;
t656 = t674 * t411;
t297 = t411 * (Icges(5,1) * t411 - t408 * t365);
t652 = t297 - t672 + (-t302 + t315) * t408;
t640 = -t186 * t408 - t411 * t634;
t405 = t411 ^ 2;
t600 = t408 / 0.2e1;
t597 = -t411 / 0.2e1;
t596 = t411 / 0.2e1;
t648 = m(6) * t410;
t263 = -rSges(6,1) * t330 - rSges(6,2) * t331;
t264 = rSges(6,1) * t332 - rSges(6,2) * t333;
t439 = t408 * t263 + t264 * t411;
t568 = m(6) * t439;
t647 = t122 * t408;
t646 = t122 * t411;
t300 = -Icges(6,3) * t407 + t452 * t410;
t524 = t409 * t303;
t531 = t406 * t307;
t645 = (t531 / 0.2e1 + t524 / 0.2e1 + t300 / 0.2e1 - t459 / 0.2e1 - t367 / 0.2e1 - t362 / 0.2e1 - t450 / 0.2e1) * t410;
t496 = t404 + t405;
t642 = (m(5) / 0.4e1 + m(6) / 0.4e1) * (0.1e1 - t496) * t410 * t407;
t641 = t496 * t410;
t639 = t330 * t247 - t331 * t250;
t436 = t524 + t531;
t425 = t300 + t436;
t533 = t299 * t410;
t632 = t425 * t407 + t533;
t440 = t247 * t409 + t250 * t406;
t426 = t299 * t408 + t440;
t535 = t244 * t410;
t631 = t426 * t407 + t535;
t427 = -t299 * t411 + t441;
t630 = t427 * t407 - t536;
t458 = Icges(6,1) * t406 + t546;
t629 = t458 * t407 + t545;
t341 = -Icges(6,2) * t529 + t386;
t346 = (Icges(6,1) * t409 - t547) * t407;
t628 = (t346 / 0.2e1 - t303 / 0.2e1) * t406 + (t307 / 0.2e1 + t341 / 0.2e1) * t409;
t623 = 0.4e1 * qJD(1);
t622 = 2 * qJD(3);
t123 = -t244 * t527 - t639;
t448 = t408 * t123 - t646;
t56 = t448 * t407 - t663;
t618 = t56 / 0.2e1;
t126 = (t279 * t411 - t408 * t661) * t407 + (t408 * t251 + t419 * t411) * t410;
t170 = (t411 * t463 + (t251 + t392) * t408) * t407;
t190 = t671 * t528;
t562 = rSges(6,3) * t407;
t318 = t462 * t410 - t562;
t153 = (t279 - t532) * t410 + ((-t318 - t562) * t408 - t463) * t407;
t154 = (-t318 * t411 + t251) * t407;
t445 = -t153 * t411 + t154 * t408;
t539 = t204 * t411;
t614 = m(6) * (t190 + (t126 + t539) * t407 + (t170 - t445) * t410);
t613 = m(6) * (t126 * t170 - t153 * t204 + t154 * t671);
t612 = m(6) * (t153 * t186 + t154 * t634 - t204 * t227 + t228 * t671);
t349 = (rSges(6,1) * t409 - rSges(6,2) * t406) * t407;
t609 = m(6) * (-t239 * t263 - t241 * t264 + (-t186 * t411 + t408 * t634) * t349);
t323 = t411 * t473;
t350 = pkin(3) * t528 - t383;
t132 = -t323 + t658 * t411 + (-t594 * t528 - t350 - t463) * t408;
t468 = t239 * t528 + t241 * t527;
t605 = m(6) * (t132 * t641 + t468);
t604 = m(6) * (t170 * t641 + t204 * t527 + t190);
t276 = t303 * t411;
t278 = t629 * t411;
t105 = t427 * t410 + (-t276 * t409 - t278 * t406 + t243) * t407;
t603 = t105 / 0.2e1;
t275 = t303 * t408;
t277 = t629 * t408;
t106 = t426 * t410 + (t275 * t409 + t277 * t406 - t244) * t407;
t602 = t106 / 0.2e1;
t258 = Icges(6,5) * t332 - Icges(6,6) * t333;
t511 = -Icges(6,2) * t333 + t250 + t327;
t513 = Icges(6,1) * t332 - t247 - t548;
t111 = t258 * t410 + (t513 * t406 + t511 * t409) * t407;
t601 = t111 / 0.2e1;
t599 = t408 / 0.4e1;
t598 = t410 / 0.2e1;
t595 = t411 / 0.4e1;
t593 = m(3) * ((rSges(3,3) * t411 + t397) * t411 + (rSges(3,3) + qJ(2)) * t404);
t464 = rSges(4,1) * t407 + rSges(4,2) * t410;
t420 = -t408 * rSges(4,3) + t464 * t411;
t268 = -t617 * t408 + t397 + t420;
t269 = (rSges(4,3) + t617) * t411 + (qJ(2) + t464) * t408;
t592 = m(4) * (t268 * t354 + t269 * t352);
t591 = m(4) * (t268 * t411 + t269 * t408);
t373 = rSges(5,2) * t407 + t563;
t182 = t373 * t405 - t323 + (-t350 - t498) * t408;
t467 = t282 * t528 + t284 * t527;
t587 = m(5) * (t182 * t641 + t467);
t586 = m(5) * (t229 * t272 + t230 * t271);
t584 = t410 * t583;
t578 = m(6) * (t186 * t227 + t228 * t634);
t577 = t640 * t648;
t576 = m(6) * t640;
t575 = t148 * t648;
t574 = m(6) * t148;
t189 = t263 * t411 - t408 * t264;
t573 = m(6) * (t189 * t407 - t349 * t641);
t569 = t410 * t568;
t565 = m(6) * qJD(5);
t560 = t408 * t677;
t558 = t411 * t56;
t553 = t124 + (t243 * t408 + t244 * t411) * t407 + t639 + t670;
t308 = -Icges(6,5) * t407 + t458 * t410;
t530 = t406 * t308;
t304 = -Icges(6,6) * t407 + t455 * t410;
t523 = t409 * t304;
t338 = (Icges(6,5) * t409 - Icges(6,6) * t406) * t407;
t522 = t410 * t338;
t422 = m(6) * t445;
t431 = t496 * t349 * t619;
t76 = -t422 / 0.2e1 + t431;
t520 = t76 * qJD(2);
t514 = -Icges(6,1) * t330 - t245 - t326;
t512 = -Icges(6,2) * t331 - t249 - t325;
t508 = -t303 + t346;
t505 = t307 + t341;
t495 = qJD(1) * t407;
t494 = qJD(3) * t407;
t493 = qJD(3) * t410;
t492 = qJD(5) * t407;
t15 = t663 + (t553 * t408 + t646) * t407;
t488 = t15 / 0.2e1 + t618;
t115 = -t330 * t304 + t331 * t308 - t632 * t411;
t89 = t330 * t276 - t331 * t278 - t630 * t411;
t90 = -t330 * t275 + t331 * t277 - t631 * t411;
t12 = (t115 + t448) * t410 + (t408 * t90 - t411 * t89 + t168) * t407;
t257 = -Icges(6,5) * t330 - Icges(6,6) * t331;
t97 = -t257 * t527 - t512 * t330 + t514 * t331;
t98 = -t258 * t527 - t511 * t330 + t513 * t331;
t46 = t97 * t408 + t411 * t98;
t487 = -t46 / 0.2e1 + t12 / 0.2e1;
t114 = t304 * t332 + t308 * t333 + t632 * t408;
t87 = -t276 * t332 - t278 * t333 + t630 * t408;
t88 = t275 * t332 + t277 * t333 + t631 * t408;
t11 = (t114 + t447) * t410 + (t408 * t88 - t411 * t87 - t169) * t407;
t100 = t258 * t528 + t511 * t332 + t513 * t333;
t99 = t257 * t528 + t512 * t332 + t514 * t333;
t47 = t100 * t411 + t99 * t408;
t486 = t47 / 0.2e1 - t11 / 0.2e1;
t192 = t411 * (Icges(4,3) * t411 + t453 * t408) + t305 * t525 + t309 * t528;
t482 = t528 / 0.4e1;
t481 = -t527 / 0.4e1;
t476 = t365 / 0.2e1 - t453 / 0.2e1;
t475 = t464 * t496;
t110 = t257 * t410 + (t514 * t406 + t512 * t409) * t407;
t138 = -t505 * t330 + t508 * t331 - t338 * t527;
t139 = t505 * t332 + t508 * t333 + t338 * t528;
t465 = t609 / 0.2e1 + (t110 + t138) * t599 + (t111 + t139) * t595;
t151 = t440 * t407 + t535;
t446 = -t149 * t411 + t151 * t408;
t213 = t264 * t410 - t349 * t528;
t214 = -t410 * t263 - t349 * t527;
t444 = t213 * t411 - t214 * t408;
t372 = qJ(4) * t410 - t567;
t355 = t408 * t372;
t238 = t355 + (-pkin(7) * t407 + t318) * t408;
t240 = t395 + (-t318 - t372) * t411;
t442 = t238 * t408 - t240 * t411;
t281 = t373 * t408 + t355;
t283 = (-t372 - t373) * t411;
t437 = t281 * t408 - t283 * t411;
t432 = -t312 * t410 - t314 * t407;
t430 = m(6) * (t186 * t264 - t263 * t634) + t522 / 0.2e1;
t424 = -t315 + t432;
t23 = t553 * t411 - t647;
t69 = t123 * t411 + t647;
t418 = t15 * t595 + t677 * t599 - t560 / 0.4e1 + t558 / 0.4e1 + (t23 + t69) * t482 + (t481 + t527 / 0.4e1) * t673;
t199 = t432 * t408 + t297;
t417 = -t659 * t408 / 0.2e1 + (t424 * t411 + t659 - t660) * t600 + (t199 + t192) * t597 + ((t302 + t435) * t408 + t192 + t652 + t672) * t596;
t416 = t69 / 0.2e1 + t23 / 0.2e1 + t404 * t302 / 0.2e1 + (t424 * t408 - t199 + t652) * t600 + (-t296 + (t302 - t435) * t411 + t659 + t660) * t596;
t147 = t425 * t410 + (-t299 + t523 + t530) * t407;
t191 = t436 * t407 + t533;
t415 = t147 * t598 - t191 * t407 / 0.2e1 + t612 / 0.2e1 + (t106 + t114) * t482 + (t105 + t115) * t481 + (t151 + t169) * t525 / 0.4e1 - (t149 - t168) * t521 / 0.4e1;
t412 = t530 / 0.2e1 + t523 / 0.2e1 - t299 / 0.2e1 - t369 / 0.2e1 + t457 / 0.2e1 + t399 - t544 / 0.2e1 + t541 / 0.2e1;
t324 = t411 * t497;
t265 = 0.2e1 * t496 * t407 * (t620 + t619);
t221 = 0.4e1 * t642;
t202 = -t271 * t408 + t461 * t405 + t324;
t183 = -t568 / 0.2e1;
t179 = t439 * t407;
t175 = t569 / 0.2e1;
t157 = -pkin(7) * t641 + t509 * t408 - t411 * t661 + t324;
t156 = (t522 + (t508 * t406 + t505 * t409) * t407) * t410;
t152 = t573 / 0.2e1;
t141 = t574 / 0.2e1;
t130 = -t575 / 0.2e1;
t93 = t604 / 0.2e1;
t82 = t577 - t584;
t75 = t422 / 0.2e1 + t431;
t66 = t628 * t407 + t430;
t64 = t132 * t189 + (t239 * t408 + t241 * t411) * t349;
t62 = t191 * t410 + t446 * t407;
t61 = -t576 + t583 + t591 + t593;
t58 = t587 + t605;
t45 = t141 + t568 / 0.2e1;
t44 = t183 + t141;
t43 = t183 - t574 / 0.2e1;
t40 = t483 + t518;
t39 = t89 * t408 + t411 * t90;
t38 = t87 * t408 + t411 * t88;
t35 = t130 - t569 / 0.2e1;
t34 = t175 + t130;
t33 = t175 + t575 / 0.2e1;
t30 = t614 / 0.2e1;
t29 = t139 * t410 + (t100 * t408 - t411 * t99) * t407;
t28 = t138 * t410 + (t408 * t98 - t411 * t97) * t407;
t27 = t412 * t407 + t578 + t586 + t592 + t645;
t20 = (t147 + t446) * t410 + (-t105 * t411 + t106 * t408 - t191) * t407;
t19 = t93 + t30 - t573 / 0.2e1;
t18 = t152 + t93 - t614 / 0.2e1;
t17 = t152 + t30 - t604 / 0.2e1;
t9 = t556 + t557;
t7 = m(6) * t64 + t46 * t600 + t47 * t596;
t6 = t488 * t528;
t5 = t613 + (t560 / 0.2e1 - t558 / 0.2e1 + t20 / 0.2e1) * t410 + (t11 * t600 + t12 * t597 - t62 / 0.2e1) * t407;
t4 = t417 * t408 + t416 * t411;
t3 = t415 + (-t15 / 0.4e1 - t56 / 0.4e1) * t411 + (-t69 / 0.4e1 - t23 / 0.4e1) * t528 + t465;
t2 = t415 - t609 / 0.2e1 + t418 + (-t138 / 0.4e1 - t110 / 0.4e1) * t408 + (-t139 / 0.4e1 - t111 / 0.4e1) * t411;
t1 = (t191 / 0.2e1 + (t115 / 0.4e1 + t105 / 0.4e1) * t411 + (-t114 / 0.4e1 - t106 / 0.4e1) * t408) * t407 + (-t147 / 0.2e1 + (-t168 / 0.4e1 + t149 / 0.4e1) * t411 + (-t169 / 0.4e1 - t151 / 0.4e1) * t408) * t410 - t612 / 0.2e1 + t418 + t465;
t8 = [t61 * qJD(2) + t27 * qJD(3) + t82 * qJD(4) + t66 * qJD(5), qJD(1) * t61 + qJD(3) * t40 + qJD(5) * t44, t27 * qJD(1) + t40 * qJD(2) + t9 * qJD(4) + t3 * qJD(5) + ((t229 * t281 + t230 * t283 - t271 * t284 + t272 * t282) * t620 + (t186 * t240 - t227 * t241 + t228 * t239 + t238 * t634) * t619) * t622 + ((t602 + t114 / 0.2e1 + m(4) * (t269 * t464 - t352 * t378) + t476 * t411 - t416) * qJD(3) + (-t305 / 0.2e1 + t347 / 0.2e1 + t312 / 0.2e1 + t336 / 0.2e1) * t493 + (-t309 / 0.2e1 - t344 / 0.2e1 + t314 / 0.2e1 - t334 / 0.2e1) * t494) * t411 + ((t603 + t115 / 0.2e1 + m(4) * (-t268 * t464 + t354 * t378) + t476 * t408 - t417) * qJD(3) + (t306 / 0.2e1 - t348 / 0.2e1 + t311 / 0.2e1 - t337 / 0.2e1) * t493 + (t310 / 0.2e1 + t345 / 0.2e1 + t313 / 0.2e1 + t335 / 0.2e1) * t494) * t408, qJD(1) * t82 + qJD(3) * t9 + qJD(5) * t34, t66 * qJD(1) + t44 * qJD(2) + t3 * qJD(3) + t34 * qJD(4) + (m(6) * (t186 * t213 - t204 * t264 + t214 * t634 - t263 * t671) + t156) * qJD(5) + ((-t138 / 0.2e1 - t110 / 0.2e1) * t411 + (t139 / 0.2e1 + t601 - t488) * t408) * t492; -t42 * qJD(3) + t43 * qJD(5) + (t576 / 0.4e1 - t583 / 0.4e1 - t591 / 0.4e1 - t593 / 0.4e1) * t623, 0, -t676 + (t437 * t620 + t442 * t619 - m(4) * t475 / 0.2e1) * t622 + t265 * qJD(4) + t75 * qJD(5), t265 * qJD(3), t43 * qJD(1) + t75 * qJD(3) - t444 * t565; t42 * qJD(2) + t4 * qJD(3) + t655 * qJD(4) + t1 * qJD(5) + (-t592 / 0.4e1 - t586 / 0.4e1 - t578 / 0.4e1) * t623 - t412 * t495 - t645 * qJD(1), qJD(5) * t76 + t676, t4 * qJD(1) + t58 * qJD(4) + t7 * qJD(5) + (m(4) * ((-t411 * t420 + (-t411 * rSges(4,3) - t464 * t408) * t408) * (-t408 * t352 - t354 * t411) - t378 * t475) + m(6) * (t132 * t157 + t238 * t239 - t240 * t241) + m(5) * (t182 * t202 + t281 * t282 - t283 * t284) + (t39 - t656 * t404 + (t408 * t657 + t669) * t411) * t600 + (t38 + t657 * t405 + (-t411 * t656 - t669) * t408) * t596) * qJD(3), t58 * qJD(3) - 0.4e1 * t642 * qJD(4) + t18 * qJD(5) + t678, t1 * qJD(1) + t520 + t7 * qJD(3) + t18 * qJD(4) + (t62 / 0.2e1 + t487 * t411 + t486 * t408) * t492 + (m(6) * (t179 * t132 + t170 * t189 - t213 * t241 + t214 * t239 + (t408 * t671 + t539) * t349) + t29 * t596 + t28 * t600 - t613 + (-t20 / 0.2e1 + (t601 + t618) * t411 + (t110 / 0.2e1 - t677 / 0.2e1) * t408) * t410) * qJD(5); -t655 * qJD(3) + t33 * qJD(5) + (-t577 / 0.4e1 + t584 / 0.4e1) * t623, 0, -t678 + t221 * qJD(4) + t17 * qJD(5) + 0.4e1 * (-t605 / 0.4e1 - t587 / 0.4e1) * qJD(3) + ((t407 * t157 + t468) * t619 + (t407 * t202 + t467) * t620 + ((t132 - t442) * t619 + (t182 - t437) * t620) * t410) * t622, t221 * qJD(3), t33 * qJD(1) + t17 * qJD(3) + (t179 * t407 + t410 * t444) * t565; -t430 * qJD(1) + t45 * qJD(2) + t2 * qJD(3) + t35 * qJD(4) + t6 * qJD(5) - t628 * t495, qJD(1) * t45 - qJD(3) * t76, t2 * qJD(1) - t520 + (((-t69 / 0.2e1 + t602) * t410 - t486) * t411 + ((t673 / 0.2e1 + t603) * t410 + t487) * t408 + ((-t39 / 0.2e1 - t151 / 0.2e1) * t411 + (t38 / 0.2e1 - t149 / 0.2e1) * t408) * t407 + (t126 * t132 - t153 * t241 + t154 * t239 + t157 * t170 - t204 * t240 + t238 * t671 - t64) * m(6)) * qJD(3) + t19 * qJD(4) + t5 * qJD(5), qJD(1) * t35 + qJD(3) * t19, t6 * qJD(1) + t5 * qJD(3) + (m(6) * (t170 * t179 - t204 * t213 + t214 * t671) + t156 * t598 + (t29 * t600 + t28 * t597 + (-t110 * t411 + t111 * t408) * t598) * t407) * qJD(5);];
Cq = t8;
