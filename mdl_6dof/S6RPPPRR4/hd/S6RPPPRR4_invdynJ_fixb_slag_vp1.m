% Calculate vector of inverse dynamics joint torques for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPPRR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_slag_vp1: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR4_invdynJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR4_invdynJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR4_invdynJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:18
% EndTime: 2019-03-09 01:35:47
% DurationCPUTime: 27.00s
% Computational Cost: add. (15665->919), mult. (36608->1141), div. (0->0), fcn. (40940->8), ass. (0->421)
t349 = sin(qJ(5));
t351 = cos(qJ(5));
t550 = sin(pkin(9));
t551 = cos(pkin(9));
t569 = sin(qJ(1));
t570 = cos(qJ(1));
t294 = -t550 * t569 - t551 * t570;
t295 = t570 * t550 - t569 * t551;
t350 = cos(qJ(6));
t348 = sin(qJ(6));
t524 = t348 * t349;
t181 = t294 * t524 - t295 * t350;
t528 = t294 * t351;
t522 = t349 * t350;
t178 = t294 * t522 + t295 * t348;
t546 = Icges(7,4) * t178;
t100 = -Icges(7,2) * t181 - Icges(7,6) * t528 + t546;
t172 = Icges(7,4) * t181;
t103 = Icges(7,1) * t178 - Icges(7,5) * t528 - t172;
t632 = t100 * t348 - t103 * t350;
t97 = Icges(7,5) * t178 - Icges(7,6) * t181 - Icges(7,3) * t528;
t38 = t349 * t97 - t351 * t632;
t179 = -t294 * t350 - t295 * t524;
t525 = t295 * t351;
t180 = -t294 * t348 + t295 * t522;
t547 = Icges(7,4) * t180;
t101 = Icges(7,2) * t179 - Icges(7,6) * t525 + t547;
t171 = Icges(7,4) * t179;
t104 = Icges(7,1) * t180 - Icges(7,5) * t525 + t171;
t561 = t181 * t101 - t104 * t178;
t562 = -t100 * t179 - t103 * t180;
t98 = Icges(7,5) * t180 + Icges(7,6) * t179 - Icges(7,3) * t525;
t636 = -t562 - t351 * (t294 * t98 + t295 * t97) - t561;
t563 = t179 * t101 + t180 * t104;
t607 = t100 * t181 - t103 * t178;
t635 = t607 + t563 + (t294 * t97 - t295 * t98) * t351;
t608 = t178 * rSges(7,1) - t181 * rSges(7,2);
t108 = rSges(7,3) * t528 - t608;
t483 = qJD(6) * t351;
t492 = qJD(5) * t295;
t199 = t294 * t483 - t492;
t431 = rSges(7,1) * t350 - rSges(7,2) * t348;
t227 = -rSges(7,3) * t349 - t351 * t431;
t484 = qJD(6) * t349;
t326 = qJD(1) - t484;
t568 = pkin(5) * t351;
t321 = -pkin(8) * t349 - t568;
t487 = qJD(5) * t321;
t273 = qJD(4) * t295;
t337 = qJD(2) * t569;
t503 = t273 + t337;
t634 = t326 * t108 - t199 * t227 + t295 * t487 - t503;
t523 = t348 * t351;
t324 = Icges(7,4) * t523;
t521 = t350 * t351;
t542 = Icges(7,5) * t349;
t225 = -Icges(7,1) * t521 + t324 - t542;
t544 = Icges(7,4) * t350;
t413 = -Icges(7,2) * t348 + t544;
t593 = Icges(7,6) * t349 + t351 * t413;
t411 = Icges(7,5) * t350 - Icges(7,6) * t348;
t594 = Icges(7,3) * t349 + t351 * t411;
t606 = t178 * t225 + t181 * t593 + t528 * t594;
t630 = t606 * t326;
t306 = Icges(6,5) * t349 + Icges(6,6) * t351;
t154 = Icges(6,3) * t295 + t294 * t306;
t535 = t154 * t295;
t341 = Icges(6,4) * t349;
t414 = Icges(6,2) * t351 + t341;
t157 = Icges(6,6) * t295 + t294 * t414;
t627 = t157 * t349;
t626 = t157 * t351;
t493 = qJD(5) * t294;
t200 = -t295 * t483 - t493;
t220 = -Icges(7,3) * t351 + t349 * t411;
t404 = -t225 * t350 - t348 * t593;
t410 = t101 * t348 - t104 * t350;
t589 = t199 * (-t294 * t594 - t632) + t200 * (t295 * t594 + t410) + t326 * (t220 + t404);
t625 = t589 * t351;
t291 = t294 * pkin(3);
t536 = qJ(4) * t295;
t210 = -t291 + t536;
t495 = t570 * pkin(1) + t569 * qJ(2);
t305 = qJD(1) * t495;
t338 = qJD(2) * t570;
t455 = qJD(1) * t570;
t497 = -pkin(2) * t455 + t338;
t462 = t305 - t497;
t482 = t294 * qJD(4);
t623 = -qJD(1) * t210 - t462 + t482;
t412 = Icges(6,5) * t351 - Icges(6,6) * t349;
t183 = t412 * t295;
t548 = Icges(6,4) * t351;
t309 = Icges(6,2) * t349 - t548;
t311 = -Icges(6,1) * t351 + t341;
t402 = t309 * t351 + t311 * t349;
t595 = t294 * t402 - t183;
t622 = t595 * qJD(1);
t494 = qJD(1) * t294;
t254 = pkin(7) * t494;
t433 = rSges(6,1) * t351 - rSges(6,2) * t349;
t489 = qJD(5) * t433;
t621 = t294 * t489;
t312 = rSges(6,1) * t349 + rSges(6,2) * t351;
t620 = t312 * t294;
t184 = t412 * t294;
t267 = t295 * qJD(1);
t153 = -t494 * pkin(3) + t267 * qJ(4) - t482;
t387 = -t153 + t497;
t262 = -t338 + t305;
t474 = t569 * pkin(2);
t475 = t569 * pkin(1);
t380 = -t475 - t474;
t619 = t380 + t474;
t618 = t254 + t623;
t485 = qJD(5) * t351;
t392 = t295 * t485 + t494 * t349;
t617 = -rSges(7,1) * t521 + rSges(7,2) * t523;
t340 = t570 * qJ(2);
t313 = t475 - t340;
t304 = qJD(1) * t313;
t498 = t337 - t304;
t499 = qJ(2) * t455 + t337;
t616 = t499 - t498;
t345 = t349 * pkin(5);
t566 = pkin(8) * t351;
t439 = t345 - t566;
t615 = t439 * qJD(5);
t310 = Icges(6,1) * t349 + t548;
t501 = t309 - t310;
t502 = t311 + t414;
t614 = (t349 * t501 - t351 * t502) * qJD(1);
t613 = t267 * pkin(7) + t499;
t609 = t387 - t305;
t605 = t295 * rSges(4,1) - t294 * rSges(4,2);
t279 = t294 * rSges(5,3);
t604 = -t295 * rSges(5,2) + t279;
t275 = t294 * qJ(4);
t290 = t295 * pkin(3);
t603 = t290 + t275;
t343 = t570 * rSges(3,3);
t469 = t569 * rSges(3,1);
t314 = t469 - t343;
t602 = -t313 - t314;
t571 = rSges(7,3) + pkin(8);
t460 = t571 * t351;
t601 = -t460 + t345;
t509 = t180 * rSges(7,1) + t179 * rSges(7,2);
t107 = -rSges(7,3) * t525 + t509;
t600 = -t107 * t326 + t200 * t227 - t294 * t487 + t618;
t526 = t295 * t349;
t164 = rSges(6,1) * t526 + rSges(6,2) * t525 - t294 * rSges(6,3);
t599 = -qJD(1) * t164 + t618 + t621;
t206 = t294 * rSges(5,2) + t295 * rSges(5,3);
t598 = -qJD(1) * t206 + t623;
t242 = Icges(6,4) * t528;
t529 = t294 * t349;
t543 = Icges(6,5) * t295;
t162 = -Icges(6,1) * t529 - t242 - t543;
t77 = -t162 * t351 - t627;
t486 = qJD(5) * t349;
t530 = t267 * t351;
t389 = t294 * t486 + t530;
t531 = t267 * t349;
t390 = -t294 * t485 + t531;
t92 = Icges(6,4) * t390 + Icges(6,2) * t389 - Icges(6,6) * t494;
t94 = Icges(6,1) * t390 + Icges(6,4) * t389 - Icges(6,5) * t494;
t597 = qJD(5) * t77 - t349 * t94 - t351 * t92;
t158 = -Icges(6,6) * t294 + t295 * t414;
t241 = Icges(6,4) * t525;
t161 = Icges(6,1) * t526 - Icges(6,5) * t294 + t241;
t76 = t158 * t349 - t161 * t351;
t459 = t295 * t486;
t532 = t494 * t351;
t391 = t459 - t532;
t91 = Icges(6,4) * t392 - Icges(6,2) * t391 + Icges(6,6) * t267;
t93 = Icges(6,1) * t392 - Icges(6,4) * t391 + Icges(6,5) * t267;
t596 = qJD(5) * t76 - t349 * t93 - t351 * t91;
t545 = Icges(7,4) * t348;
t415 = Icges(7,1) * t350 - t545;
t592 = t351 * t415 + t542;
t297 = t414 * qJD(5);
t298 = t310 * qJD(5);
t401 = t309 * t349 - t311 * t351;
t591 = qJD(5) * t401 - t297 * t351 - t298 * t349;
t278 = (Icges(7,1) * t348 + t544) * t351;
t588 = t199 * (-Icges(7,1) * t181 - t100 - t546) + t200 * (-Icges(7,1) * t179 + t101 + t547) + t326 * (-t593 - t278);
t352 = qJD(1) ^ 2;
t26 = -t525 * t98 + t563;
t27 = t525 * t97 + t562;
t65 = -t179 * t593 + t180 * t225 + t525 * t594;
t64 = t65 * t326;
t10 = t199 * t27 + t200 * t26 + t64;
t587 = -t10 / 0.2e1;
t196 = -qJD(5) * t494 - qJDD(5) * t295;
t479 = qJDD(6) * t351;
t110 = -qJD(6) * t389 + t294 * t479 + t196;
t586 = t110 / 0.2e1;
t197 = qJD(5) * t267 - qJDD(5) * t294;
t111 = qJD(6) * t391 - t295 * t479 + t197;
t585 = t111 / 0.2e1;
t584 = t196 / 0.2e1;
t583 = t197 / 0.2e1;
t582 = -t199 / 0.2e1;
t581 = t199 / 0.2e1;
t580 = -t200 / 0.2e1;
t579 = t200 / 0.2e1;
t578 = -t494 / 0.2e1;
t293 = -qJD(5) * t483 - qJDD(6) * t349 + qJDD(1);
t576 = t293 / 0.2e1;
t573 = -t326 / 0.2e1;
t572 = t326 / 0.2e1;
t567 = pkin(7) * t294;
t289 = t295 * pkin(7);
t367 = -qJD(6) * t294 + t392;
t456 = t295 * t484;
t435 = t267 - t456;
t85 = -t348 * t367 + t350 * t435;
t86 = t348 * t435 + t350 * t367;
t42 = Icges(7,5) * t86 + Icges(7,6) * t85 + Icges(7,3) * t391;
t44 = Icges(7,4) * t86 + Icges(7,2) * t85 + Icges(7,6) * t391;
t46 = Icges(7,1) * t86 + Icges(7,4) * t85 + Icges(7,5) * t391;
t8 = (-qJD(5) * t410 - t42) * t349 + (-qJD(5) * t98 + t348 * t44 - t350 * t46 + (t101 * t350 + t104 * t348) * qJD(6)) * t351;
t565 = t8 * t200;
t366 = qJD(6) * t295 - t390;
t457 = t294 * t484;
t436 = -t494 + t457;
t87 = t348 * t366 + t350 * t436;
t88 = t348 * t436 - t350 * t366;
t43 = Icges(7,5) * t88 + Icges(7,6) * t87 - Icges(7,3) * t389;
t45 = Icges(7,4) * t88 + Icges(7,2) * t87 - Icges(7,6) * t389;
t47 = Icges(7,1) * t88 + Icges(7,4) * t87 - Icges(7,5) * t389;
t9 = (qJD(5) * t632 - t43) * t349 + (qJD(5) * t97 + t348 * t45 - t350 * t47 + (-t100 * t350 - t103 * t348) * qJD(6)) * t351;
t564 = t9 * t199;
t560 = rSges(7,3) * t351;
t559 = pkin(7) * qJD(1);
t558 = t494 * rSges(6,3);
t281 = t295 * rSges(6,3);
t55 = t154 * t294 - t157 * t525 + t162 * t526;
t556 = t295 * t55;
t37 = -t349 * t98 + t351 * t410;
t554 = t37 * t111;
t553 = t38 * t110;
t126 = t349 * t594 + t351 * t404;
t276 = (Icges(7,5) * t348 + Icges(7,6) * t350) * t351;
t167 = qJD(5) * t220 + qJD(6) * t276;
t222 = -Icges(7,6) * t351 + t349 * t413;
t168 = (Icges(7,2) * t350 + t545) * t483 + t222 * qJD(5);
t224 = -Icges(7,5) * t351 + t349 * t415;
t169 = qJD(5) * t224 + qJD(6) * t278;
t41 = (-qJD(5) * t404 - t167) * t349 + (qJD(5) * t594 + t168 * t348 - t169 * t350 + (t225 * t348 - t350 * t593) * qJD(6)) * t351;
t552 = t126 * t293 + t41 * t326;
t534 = t494 * qJ(4);
t258 = pkin(5) * t526;
t192 = -pkin(8) * t525 + t258;
t518 = t107 + t192;
t194 = t439 * t294;
t517 = t108 - t194;
t516 = -t153 - t262;
t515 = t295 * t311 + t158;
t514 = -t294 * t311 - t157;
t513 = -Icges(6,2) * t526 + t161 + t241;
t512 = Icges(6,2) * t529 + t162 - t242;
t292 = (rSges(7,1) * t348 + rSges(7,2) * t350) * t351;
t170 = qJD(6) * t292 + (t349 * t431 - t560) * qJD(5);
t511 = t170 + t615;
t507 = t227 + t321;
t505 = t267 * rSges(4,1) - rSges(4,2) * t494;
t504 = -t267 * rSges(5,2) + rSges(5,3) * t494;
t193 = pkin(5) * t525 + pkin(8) * t526;
t195 = -pkin(5) * t528 - pkin(8) * t529;
t211 = -t294 * rSges(4,1) - t295 * rSges(4,2);
t500 = qJD(1) * t338 + qJDD(2) * t569;
t496 = t570 * rSges(3,1) + t569 * rSges(3,3);
t299 = t312 * qJD(5);
t491 = qJD(5) * t299;
t490 = qJD(5) * t615;
t481 = t306 * qJD(1);
t480 = -m(5) - m(6) - m(7);
t478 = -t570 / 0.2e1;
t477 = t569 / 0.2e1;
t476 = t86 * rSges(7,1) + t85 * rSges(7,2) + rSges(7,3) * t459;
t346 = t570 * pkin(2);
t470 = m(4) - t480;
t155 = -Icges(6,3) * t294 + t295 * t306;
t56 = -t295 * t155 - t158 * t528 - t161 * t529;
t464 = t294 * t155 - t158 * t525 - t161 * t526;
t463 = pkin(5) * t392 + pkin(8) * t459;
t163 = rSges(6,1) * t529 + rSges(6,2) * t528 + t281;
t461 = t346 + t495;
t454 = qJD(1) * t569;
t453 = -t493 / 0.2e1;
t452 = t493 / 0.2e1;
t451 = -t492 / 0.2e1;
t450 = t492 / 0.2e1;
t255 = t267 * pkin(3);
t449 = t255 + t273;
t204 = qJD(1) * t603;
t448 = t204 + t273 + t498;
t447 = rSges(6,1) * t392 + rSges(6,2) * t532 + t267 * rSges(6,3);
t445 = -t291 + t461;
t444 = -t313 - t474;
t443 = rSges(7,1) * t522 - rSges(7,2) * t524;
t440 = t88 * rSges(7,1) + t87 * rSges(7,2);
t438 = qJ(4) - t460;
t437 = t560 + t566;
t434 = rSges(4,1) * t494 + t267 * rSges(4,2);
t430 = -rSges(5,2) * t494 - t267 * rSges(5,3);
t406 = t162 * t349 - t626;
t407 = t158 * t351 + t161 * t349;
t89 = Icges(6,5) * t392 - Icges(6,6) * t391 + Icges(6,3) * t267;
t90 = Icges(6,5) * t390 + Icges(6,6) * t389 - Icges(6,3) * t494;
t429 = -(t155 * t267 - t294 * t89 - t295 * t596 + t407 * t494) * t294 - (-t154 * t267 - t294 * t90 - t295 * t597 + t406 * t494) * t295;
t428 = -(-t155 * t494 + t407 * t267 + t294 * t596 - t295 * t89) * t294 - (t154 * t494 + t406 * t267 + t294 * t597 - t295 * t90) * t295;
t427 = t26 * t295 - t27 * t294;
t28 = t528 * t98 + t561;
t29 = -t528 * t97 - t607;
t426 = t28 * t295 - t29 * t294;
t265 = t295 * t559;
t400 = t444 + t603;
t384 = t194 + t400;
t35 = qJD(1) * t384 + t265 - t634;
t36 = qJD(1) * t192 - t600;
t425 = t294 * t35 + t295 * t36;
t424 = t294 * t38 - t295 * t37;
t423 = t294 * t464 - t556;
t57 = -t294 * t406 + t535;
t422 = -t294 * t56 - t295 * t57;
t165 = -t281 - t620;
t385 = -t165 + t400;
t62 = qJD(1) * t385 + t295 * t489 + t265 + t503;
t421 = -t294 * t599 - t295 * t62;
t95 = -rSges(6,2) * t459 + t447;
t96 = rSges(6,1) * t390 + rSges(6,2) * t389 - t558;
t420 = t294 * t96 - t295 * t95;
t408 = t107 * t294 + t108 * t295;
t405 = t164 * t295 - t165 * t294;
t399 = t605 + t444;
t397 = -t351 * t42 + t486 * t98;
t396 = -t351 * t43 - t486 * t97;
t388 = t449 + t534;
t386 = -t346 * t352 + t500;
t383 = t400 + t604;
t148 = rSges(7,3) * t526 - t295 * t617;
t149 = -rSges(7,3) * t529 + t294 * t617;
t382 = t445 + t536;
t379 = -t469 - t475;
t378 = qJDD(1) * t495 - qJDD(2) * t570 + (-pkin(1) * t454 + t337 + t499) * qJD(1);
t319 = rSges(2,1) * t570 - rSges(2,2) * t569;
t315 = rSges(2,1) * t569 + rSges(2,2) * t570;
t375 = -t199 * t97 + t200 * t98 - t326 * t594;
t374 = (Icges(7,5) * t179 - Icges(7,6) * t180) * t200 + (Icges(7,5) * t181 + Icges(7,6) * t178) * t199 + t276 * t326;
t373 = t340 + t380;
t372 = t349 * t515 - t351 * t513;
t371 = -t349 * t514 + t351 * t512;
t370 = qJD(4) * t494 + qJDD(4) * t295 + t386;
t369 = t275 + t373;
t368 = t290 + t373;
t363 = qJDD(1) * t289 + t494 * t559 + t370;
t23 = -t295 * t491 - t196 * t433 + (-t96 + t516) * qJD(1) + t385 * qJDD(1) + t363;
t358 = qJDD(1) * t346 - t352 * t474 + t378;
t355 = qJD(1) * t388 + qJD(4) * t267 + qJDD(1) * t210 - qJDD(4) * t294 + t358;
t354 = (qJD(1) * t267 - qJDD(1) * t294) * pkin(7) + t355;
t24 = qJD(1) * t95 + qJDD(1) * t164 + t197 * t433 + t294 * t491 + t354;
t365 = t23 * t295 - t24 * t294 - t267 * t599 + t494 * t62;
t364 = t289 + t369;
t362 = t368 + t275;
t359 = (-Icges(7,2) * t180 + t104 + t171) * t200 + (Icges(7,2) * t178 - t103 + t172) * t199 + (Icges(7,2) * t521 + t225 + t324) * t326;
t30 = -t107 * t199 + t108 * t200 - qJD(3) + (t192 * t295 + t194 * t294) * qJD(5);
t356 = -t227 * t425 + t30 * t408;
t334 = rSges(3,3) * t455;
t296 = t306 * qJD(5);
t256 = pkin(5) * t529;
t226 = t443 - t560;
t213 = qJD(1) * t602 + t337;
t191 = -pkin(8) * t528 + t256;
t190 = t433 * t294;
t189 = t433 * t295;
t160 = t294 * t310 + t543;
t150 = qJD(1) * t399 + t337;
t147 = t592 * t294;
t146 = t592 * t295;
t145 = t593 * t294;
t144 = t593 * t295;
t130 = qJDD(1) * t496 + qJD(1) * (-rSges(3,1) * t454 + t334) + t378;
t129 = -qJD(1) * t262 + qJDD(1) * t602 - t352 * t496 + t500;
t125 = rSges(7,1) * t181 + rSges(7,2) * t178;
t124 = rSges(7,1) * t179 - rSges(7,2) * t180;
t122 = t295 * t402 + t184;
t115 = pkin(5) * t390 - pkin(8) * t389;
t114 = -pkin(8) * t532 + t463;
t113 = t122 * qJD(1);
t81 = qJD(1) * t383 + t503;
t75 = qJD(1) * t505 + qJDD(1) * t211 + t358;
t74 = (-t262 + t434) * qJD(1) + t399 * qJDD(1) + t386;
t67 = qJD(5) * t405 - qJD(3);
t53 = t402 * t267 + t294 * t591 - t295 * t296 + t412 * t494;
t52 = -t267 * t412 - t294 * t296 - t295 * t591 + t402 * t494;
t51 = qJD(1) * t504 + qJDD(1) * t206 + t355;
t50 = (t430 + t516) * qJD(1) + t383 * qJDD(1) + t370;
t49 = -rSges(7,3) * t389 + t440;
t48 = -rSges(7,3) * t532 + t476;
t32 = qJD(5) * t406 + t349 * t92 - t351 * t94;
t31 = qJD(5) * t407 + t349 * t91 - t351 * t93;
t25 = qJD(5) * t420 + t164 * t196 - t165 * t197 + qJDD(3);
t22 = qJD(5) * t422 - t622;
t21 = qJD(5) * t423 + t113;
t20 = t167 * t528 + t168 * t181 - t169 * t178 + t225 * t88 + t389 * t594 - t593 * t87;
t19 = -t167 * t525 + t168 * t179 + t169 * t180 + t225 * t86 - t391 * t594 - t593 * t85;
t14 = t126 * t326 + t199 * t38 + t200 * t37;
t13 = -t295 * t490 - t293 * t108 + t110 * t227 + t199 * t170 + t196 * t321 - t326 * t49 + (-t115 + t516) * qJD(1) + t384 * qJDD(1) + t363;
t12 = qJD(1) * t114 + qJDD(1) * t192 + t293 * t107 - t111 * t227 - t200 * t170 - t197 * t321 + t294 * t490 + t326 * t48 + t354;
t11 = t199 * t29 + t200 * t28 - t630;
t7 = t107 * t110 - t111 * t108 - t114 * t492 + t115 * t493 + t192 * t196 + t197 * t194 + t199 * t48 - t200 * t49 + qJDD(3);
t6 = -t100 * t87 - t103 * t88 - t178 * t47 + t181 * t45 - t294 * t396 + t530 * t97;
t5 = t101 * t87 + t104 * t88 - t178 * t46 + t181 * t44 - t294 * t397 - t530 * t98;
t4 = -t100 * t85 - t103 * t86 + t179 * t45 + t180 * t47 + t295 * t396 + t532 * t97;
t3 = t101 * t85 + t104 * t86 + t179 * t44 + t180 * t46 + t295 * t397 - t532 * t98;
t2 = t110 * t29 + t111 * t28 + t199 * t6 + t20 * t326 + t200 * t5 - t293 * t606;
t1 = t110 * t27 + t111 * t26 + t19 * t326 + t199 * t4 + t200 * t3 + t293 * t65;
t15 = [(t32 + t53 + t21) * t451 + (t122 + t76) * t583 + (-g(1) * (t163 + t362 + t289) + (-g(2) + t24) * (t164 + t382 - t567) + (t364 + t620) * t23 - (pkin(2) * t454 + t255 + t447 - t448 + t534 + (-t163 + t380) * qJD(1) + t613) * t599 + (-rSges(6,1) * t531 - rSges(6,2) * t530 + t254 + t558 - t599 + t609 + t621) * t62 + (t23 * (rSges(6,3) + pkin(3)) - t67 * (t163 + t165) * qJD(5) - (-rSges(6,2) * t486 + qJD(4) - t489 - t559) * t599) * t295) * m(6) + (t50 * (t279 + (-rSges(5,2) + pkin(3)) * t295 + t369) - g(1) * (t362 + t604) - (t388 + t499 + t504 - t448 + (-t604 + t619) * qJD(1)) * t598 + (t430 - t598 + t609) * t81 + (t51 - g(2)) * (t382 + t206)) * m(5) + (-t595 + t77) * t584 + t65 * t585 + ((-g(2) + t75) * (t461 + t211) + (t434 - t462) * t150 + (t505 + t150 + (-t605 + t619) * qJD(1) + t616) * (qJD(1) * t211 + t462) + (-g(1) + t74) * (t373 + t605)) * m(4) - t606 * t586 + (t213 * t338 + (-g(1) + t129) * (t340 + t343 + t379) + (-qJD(1) * t213 - g(2) + t130) * (t495 + t496) + (t213 + t334 + (t314 + t379) * qJD(1) + t616) * (qJD(1) * t496 + t262)) * m(3) + (-(t294 * t438 + t256 + t289 + t368 + t608) * g(1) - t30 * (t191 - t194) * t492 + (t12 - g(2)) * (t295 * t438 + t258 + t445 + t509 - t567) + (t294 * t601 + t290 + t364 + t608) * t13 + (-t204 + t304 + t449 + t463 + t476 - (-qJ(4) + t437) * t494 + t613 + t634) * t36 + (-t600 + t254 - t440 + (t437 - t345) * t267 + (t349 * t571 + t568) * t493 + t387) * t35 + ((-t191 - t289 + t619) * t36 + (t192 - t495) * t35) * qJD(1)) * m(7) + t565 / 0.2e1 + t564 / 0.2e1 - m(2) * (-g(1) * t315 + g(2) * t319) + t553 / 0.2e1 + t554 / 0.2e1 + t19 * t579 + t20 * t581 + (qJD(5) * t402 + t297 * t349 - t298 * t351) * qJD(1) + (t64 + (t29 + t635) * t200 + (-t28 - t636) * t199) * t582 + (t630 + (t27 + t636) * t200 + t11 + (-t26 + t635) * t199) * t580 + (t622 + (t535 * t295 + (-t55 + (t154 - t407) * t294 + (-t155 + (-t160 - t162) * t349) * t295) * t294) * qJD(5) + t22) * t452 - t199 * t587 + (t31 + t52 + qJD(1) * (t160 * t351 - t627 - t77)) * t453 + (t113 + (-t556 + (t535 - t57 + (t160 * t349 + t626) * t294 + t464) * t294) * qJD(5)) * t450 + t552 + (t401 + m(2) * (t315 ^ 2 + t319 ^ 2) + Icges(2,3) + Icges(3,2) + Icges(4,3) + Icges(5,1)) * qJDD(1); (-m(3) - t470) * (g(1) * t569 - g(2) * t570) + 0.2e1 * (t12 * t478 + t13 * t477) * m(7) + 0.2e1 * (t23 * t477 + t24 * t478) * m(6) + 0.2e1 * (t477 * t50 + t478 * t51) * m(5) + 0.2e1 * (t477 * t74 + t478 * t75) * m(4) + 0.2e1 * (t129 * t477 + t130 * t478) * m(3); (m(4) + m(5)) * qJDD(3) + m(6) * t25 + m(7) * t7 + t470 * g(3); t480 * (g(1) * t295 - g(2) * t294) + m(5) * (-t267 * t598 - t294 * t51 + t295 * t50 + t494 * t81) + m(6) * t365 + (-t12 * t294 + t13 * t295 + t267 * t36 + t35 * t494) * m(7) + (-m(5) * (t294 * t81 - t295 * t598) - m(6) * (t294 * t62 - t295 * t599) - m(7) * t425) * qJD(1); (t21 + t10) * t267 / 0.2e1 + (-g(1) * t189 + g(2) * t190 - g(3) * t312 - (-t189 * t599 + t190 * t62) * qJD(1) - (t67 * (t189 * t295 + t190 * t294) + t421 * t312) * qJD(5) - t25 * t405 + t67 * (t164 * t494 + t165 * t267 - t420) + t421 * t299 + t365 * t433) * m(6) - (qJD(1) * t53 + qJD(5) * t428 - qJDD(1) * t595 + t196 * t57 + t197 * t56 + t2) * t295 / 0.2e1 + (-g(1) * (t148 + t193) - g(2) * (t149 + t195) - g(3) * (t443 + t601) - t35 * (-qJD(1) * t195 - t149 * t326 + t199 * t226 - t295 * t615) - t36 * (qJD(1) * t193 + t148 * t326 - t200 * t226 + t294 * t615) - t30 * (-t148 * t199 + t149 * t200 + t193 * t492 - t195 * t493) - ((-t107 * t36 + t108 * t35) * t351 + t356 * t349) * qJD(6) + (t30 * t517 - t36 * t507) * t267 - (-t30 * t518 + t35 * t507) * t494 + (-t13 * t507 - t35 * t511 - t7 * t518 + t30 * (t114 + t48)) * t295 + (t12 * t507 + t36 * t511 + t7 * t517 + t30 * (-t115 - t49)) * t294) * m(7) + ((t183 * t493 - t481) * t294 + (-t614 + (-t371 * t295 + (-t184 + t372) * t294) * qJD(5)) * t295) * t452 + ((-t184 * t492 - t481) * t295 + (t614 + (-t372 * t294 + (t183 + t371) * t295) * qJD(5)) * t294) * t450 + t423 * t583 + t422 * t584 + (-t26 * t294 - t27 * t295) * t585 + (-t28 * t294 - t29 * t295) * t586 + t456 * t587 + ((t144 * t181 - t146 * t178) * t200 + (-t145 * t181 + t147 * t178) * t199 + (-t178 * t224 + t181 * t222) * t326 + (t28 * t526 + t351 * t606) * qJD(6) + ((-qJD(6) * t29 - t375) * t349 + t625) * t294) * t582 + (((t144 * t348 - t146 * t350 - t98) * t200 + (-t145 * t348 + t147 * t350 + t97) * t199 + (t222 * t348 - t224 * t350 + t594) * t326 - t126 * qJD(6)) * t351 + (-qJD(6) * t424 - t589) * t349) * t573 + ((t144 * t179 + t146 * t180) * t200 + (-t145 * t179 - t147 * t180) * t199 + (t179 * t222 + t180 * t224) * t326 + (-t27 * t529 - t351 * t65) * qJD(6) + ((qJD(6) * t26 + t375) * t349 - t625) * t295) * t580 + (t457 / 0.2e1 + t578) * t11 + qJDD(1) * (-t294 * t76 - t295 * t77) / 0.2e1 + (-t294 * t37 - t295 * t38) * t576 + t22 * t578 - (qJD(1) * t52 + qJD(5) * t429 + qJDD(1) * t122 + t196 * t55 - t197 * t464 + t1) * t294 / 0.2e1 + qJD(1) * (t267 * t76 - t294 * t31 - t295 * t32 - t494 * t77) / 0.2e1 + (t267 * t37 - t294 * t8 - t295 * t9 - t38 * t494) * t572 + (t26 * t267 - t27 * t494 - t294 * t3 - t295 * t4) * t579 + (t267 * t28 - t29 * t494 - t294 * t5 - t295 * t6) * t581 + (t267 * t56 - t494 * t57 + t428) * t451 + (-t267 * t464 - t494 * t55 + t429) * t453 + t14 * t483 / 0.2e1 - qJD(1) * ((t349 * t502 + t351 * t501) * qJD(1) + ((-t294 * t515 - t295 * t514) * t351 + (-t294 * t513 - t295 * t512) * t349) * qJD(5)) / 0.2e1; -t1 * t525 / 0.2e1 + (-t349 * t65 - t351 * t427) * t585 + ((qJD(5) * t427 - t19) * t349 + (-qJD(5) * t65 - t26 * t494 - t267 * t27 + t294 * t4 - t295 * t3) * t351) * t579 + t2 * t528 / 0.2e1 + (t349 * t606 - t351 * t426) * t586 + ((qJD(5) * t426 - t20) * t349 + (qJD(5) * t606 - t267 * t29 - t28 * t494 + t294 * t6 - t295 * t5) * t351) * t581 - t14 * t485 / 0.2e1 - t349 * (t552 + t553 + t554 + t564 + t565) / 0.2e1 + (-t126 * t349 + t351 * t424) * t576 + ((-qJD(5) * t424 - t41) * t349 + (-qJD(5) * t126 - t267 * t38 + t294 * t9 - t295 * t8 - t37 * t494) * t351) * t572 + (t179 * t359 - t180 * t588 - t374 * t525) * t580 + (t178 * t588 + t181 * t359 + t374 * t528) * t582 + (-t374 * t349 + (t359 * t348 + t350 * t588) * t351) * t573 + (-t530 / 0.2e1 + t349 * t453) * t11 + (-t532 / 0.2e1 + t349 * t450) * t10 + ((qJD(5) * t356 - t12 * t107 + t13 * t108 + t35 * t49 - t36 * t48) * t349 + (t35 * (qJD(5) * t108 + t170 * t294) + t36 * (-qJD(5) * t107 + t170 * t295) + t7 * t408 + t30 * (t107 * t267 - t108 * t494 - t294 * t48 - t295 * t49) + (t12 * t295 + t13 * t294 - t267 * t35 + t36 * t494) * t227) * t351 - t35 * (-t125 * t326 + t199 * t292) - t36 * (t124 * t326 - t200 * t292) - t30 * (-t124 * t199 + t125 * t200) - g(1) * t124 - g(2) * t125 - g(3) * t292) * m(7);];
tau  = t15;
