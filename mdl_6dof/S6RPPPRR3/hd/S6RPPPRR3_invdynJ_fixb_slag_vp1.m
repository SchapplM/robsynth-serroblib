% Calculate vector of inverse dynamics joint torques for
% S6RPPPRR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPPRR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR3_invdynJ_fixb_slag_vp1: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_invdynJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_invdynJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR3_invdynJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR3_invdynJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:28
% EndTime: 2019-03-09 01:33:57
% DurationCPUTime: 24.99s
% Computational Cost: add. (22254->899), mult. (37646->1133), div. (0->0), fcn. (42116->10), ass. (0->406)
t554 = sin(pkin(9));
t555 = cos(pkin(9));
t573 = sin(qJ(1));
t574 = cos(qJ(1));
t313 = -t554 * t573 - t555 * t574;
t314 = t574 * t554 - t573 * t555;
t357 = sin(qJ(6));
t353 = pkin(10) + qJ(5);
t342 = cos(t353);
t358 = cos(qJ(6));
t528 = t342 * t358;
t200 = t313 * t357 + t314 * t528;
t529 = t342 * t357;
t201 = -t313 * t358 + t314 * t529;
t341 = sin(t353);
t533 = t314 * t341;
t102 = Icges(7,5) * t200 - Icges(7,6) * t201 + Icges(7,3) * t533;
t551 = Icges(7,4) * t200;
t105 = -Icges(7,2) * t201 + Icges(7,6) * t533 + t551;
t186 = Icges(7,4) * t201;
t108 = Icges(7,1) * t200 + Icges(7,5) * t533 - t186;
t623 = t105 * t357 - t108 * t358;
t39 = -t102 * t342 - t341 * t623;
t493 = qJD(6) * t341;
t497 = qJD(5) * t313;
t205 = t314 * t493 + t497;
t496 = qJD(5) * t314;
t206 = -t313 * t493 + t496;
t535 = t313 * t341;
t203 = t313 * t529 + t314 * t358;
t204 = -t313 * t528 + t314 * t357;
t624 = t105 * t203 + t108 * t204;
t28 = t102 * t535 - t624;
t104 = Icges(7,5) * t204 + Icges(7,6) * t203 - Icges(7,3) * t535;
t550 = Icges(7,4) * t204;
t107 = Icges(7,2) * t203 - Icges(7,6) * t535 + t550;
t187 = Icges(7,4) * t203;
t110 = Icges(7,1) * t204 - Icges(7,5) * t535 + t187;
t29 = -t104 * t535 + t203 * t107 + t204 * t110;
t492 = qJD(6) * t342;
t329 = qJD(1) + t492;
t531 = t341 * t357;
t326 = Icges(7,4) * t531;
t530 = t341 * t358;
t546 = Icges(7,5) * t342;
t226 = -Icges(7,1) * t530 + t326 + t546;
t423 = Icges(7,5) * t358 - Icges(7,6) * t357;
t373 = -Icges(7,3) * t342 + t341 * t423;
t548 = Icges(7,4) * t358;
t426 = -Icges(7,2) * t357 + t548;
t374 = -Icges(7,6) * t342 + t341 * t426;
t66 = -t203 * t374 + t204 * t226 + t373 * t535;
t13 = -t205 * t28 + t206 * t29 + t66 * t329;
t425 = Icges(6,5) * t342 - Icges(6,6) * t341;
t157 = Icges(6,3) * t313 + t314 * t425;
t552 = Icges(6,4) * t342;
t428 = -Icges(6,2) * t341 + t552;
t160 = Icges(6,6) * t313 + t314 * t428;
t553 = Icges(6,4) * t341;
t431 = Icges(6,1) * t342 - t553;
t163 = Icges(6,5) * t313 + t314 * t431;
t611 = t160 * t341 - t163 * t342;
t58 = t157 * t313 - t611 * t314;
t26 = t102 * t533 - t105 * t201 + t108 * t200;
t345 = qJD(2) * t574;
t498 = t574 * pkin(1) + t573 * qJ(2);
t278 = qJD(1) * t498 - t345;
t444 = t200 * rSges(7,1) - t201 * rSges(7,2);
t113 = -rSges(7,3) * t533 - t444;
t565 = rSges(7,1) * t358;
t443 = rSges(7,2) * t357 - t565;
t228 = rSges(7,3) * t342 + t341 * t443;
t572 = t341 * pkin(5);
t311 = pkin(8) * t342 - t572;
t290 = qJD(4) * t314;
t344 = qJD(2) * t573;
t506 = t290 + t344;
t620 = t329 * t113 + t205 * t228 + t311 * t497 - t506;
t604 = t200 * t226 + t201 * t374 - t373 * t533;
t619 = t205 * t26 + t604 * t329;
t159 = Icges(6,3) * t314 - t313 * t425;
t616 = t314 * t159;
t294 = t313 * rSges(5,3);
t354 = sin(pkin(10));
t355 = cos(pkin(10));
t446 = -rSges(5,1) * t355 + rSges(5,2) * t354;
t615 = -t314 * t446 + t294;
t467 = qJD(1) * t574;
t609 = -pkin(2) * t467 - t278;
t424 = Icges(6,5) * t341 + Icges(6,6) * t342;
t181 = t424 * t313;
t427 = Icges(6,2) * t342 + t553;
t430 = Icges(6,1) * t341 + t552;
t412 = -t341 * t427 + t342 * t430;
t118 = t314 * t412 + t181;
t613 = qJD(1) * t118;
t419 = -t160 * t342 - t163 * t341;
t585 = -t205 / 0.2e1;
t180 = t424 * t314;
t494 = qJD(5) * t342;
t284 = t314 * qJD(1);
t539 = t284 * t341;
t399 = t313 * t494 - t539;
t495 = qJD(5) * t341;
t470 = t313 * t495;
t538 = t284 * t342;
t398 = t470 + t538;
t605 = -t201 * t107 + t110 * t200;
t603 = t314 * rSges(4,1) - t313 * rSges(4,2);
t292 = t313 * qJ(4);
t458 = t314 * pkin(3) + t292;
t347 = t574 * qJ(2);
t485 = t573 * pkin(1);
t320 = t485 - t347;
t349 = t574 * rSges(3,3);
t480 = t573 * rSges(3,1);
t321 = t480 - t349;
t602 = -t320 - t321;
t484 = t573 * pkin(2);
t393 = -t485 - t484;
t601 = t393 + t484;
t283 = t313 * qJD(1);
t459 = t284 * pkin(3) + t283 * qJ(4);
t599 = t459 + t290;
t499 = t574 * rSges(3,1) + t573 * rSges(3,3);
t598 = t498 + t499;
t318 = qJD(1) * t320;
t501 = t344 - t318;
t502 = qJ(2) * t467 + t344;
t597 = t502 - t501;
t356 = -pkin(7) - qJ(4);
t537 = t284 * t356;
t595 = t537 + t609;
t281 = t313 * t356;
t334 = pkin(4) * t355 + pkin(3);
t509 = t314 * t334 - t281;
t513 = -t283 * t356 + t284 * t334;
t594 = -qJD(1) * (-t458 + t509) + t502 + t513 + t290;
t519 = t427 * t314 - t163;
t521 = -t430 * t314 - t160;
t593 = t341 * t519 + t342 * t521;
t223 = -Icges(7,3) * t341 - t342 * t423;
t414 = -t226 * t358 - t357 * t374;
t421 = t107 * t357 - t110 * t358;
t592 = -t206 * (t373 * t313 + t421) + t205 * (t373 * t314 - t623) - t329 * (t223 + t414);
t256 = (Icges(7,1) * t357 + t548) * t341;
t591 = t206 * (-Icges(7,1) * t203 + t107 + t550) - t205 * (-Icges(7,1) * t201 - t105 - t551) + t329 * (-t374 - t256);
t359 = qJD(1) ^ 2;
t208 = qJD(5) * t283 + qJDD(5) * t314;
t488 = qJDD(6) * t341;
t115 = -qJD(6) * t399 - t313 * t488 + t208;
t590 = t115 / 0.2e1;
t209 = qJD(5) * t284 - qJDD(5) * t313;
t540 = t283 * t341;
t401 = t314 * t494 + t540;
t116 = -qJD(6) * t401 - t314 * t488 + t209;
t589 = t116 / 0.2e1;
t588 = -t206 / 0.2e1;
t587 = t206 / 0.2e1;
t586 = t205 / 0.2e1;
t584 = t208 / 0.2e1;
t583 = t209 / 0.2e1;
t269 = -qJD(5) * t493 + qJDD(6) * t342 + qJDD(1);
t582 = t269 / 0.2e1;
t581 = t283 / 0.2e1;
t580 = t284 / 0.2e1;
t577 = -t329 / 0.2e1;
t576 = t329 / 0.2e1;
t575 = rSges(7,3) + pkin(8);
t571 = t342 * pkin(5);
t400 = -t283 * t342 + t314 * t495;
t379 = -qJD(6) * t313 + t400;
t448 = t314 * t492 + t284;
t88 = -t357 * t379 + t358 * t448;
t89 = t357 * t448 + t358 * t379;
t46 = Icges(7,5) * t89 + Icges(7,6) * t88 - Icges(7,3) * t401;
t48 = Icges(7,4) * t89 + Icges(7,2) * t88 - Icges(7,6) * t401;
t50 = Icges(7,1) * t89 + Icges(7,4) * t88 - Icges(7,5) * t401;
t7 = (-qJD(5) * t623 + t46) * t342 + (qJD(5) * t102 + t357 * t48 - t358 * t50 + (-t105 * t358 - t108 * t357) * qJD(6)) * t341;
t570 = t7 * t205;
t378 = qJD(6) * t314 + t398;
t449 = t313 * t492 + t283;
t90 = -t357 * t378 + t358 * t449;
t91 = t357 * t449 + t358 * t378;
t47 = Icges(7,5) * t91 + Icges(7,6) * t90 - Icges(7,3) * t399;
t49 = Icges(7,4) * t91 + Icges(7,2) * t90 - Icges(7,6) * t399;
t51 = Icges(7,1) * t91 + Icges(7,4) * t90 - Icges(7,5) * t399;
t8 = (qJD(5) * t421 + t47) * t342 + (-qJD(5) * t104 + t357 * t49 - t358 * t51 + (t107 * t358 + t110 * t357) * qJD(6)) * t341;
t569 = t8 * t206;
t568 = pkin(3) - t334;
t111 = t341 * t414 - t342 * t373;
t254 = (Icges(7,5) * t357 + Icges(7,6) * t358) * t341;
t172 = qJD(5) * t223 + qJD(6) * t254;
t225 = -Icges(7,6) * t341 - t342 * t426;
t549 = Icges(7,4) * t357;
t173 = (Icges(7,2) * t358 + t549) * t493 + t225 * qJD(5);
t429 = Icges(7,1) * t358 - t549;
t227 = -Icges(7,5) * t341 - t342 * t429;
t174 = qJD(5) * t227 + qJD(6) * t256;
t36 = (qJD(5) * t414 + t172) * t342 + (qJD(5) * t373 + t173 * t357 - t174 * t358 + (t226 * t357 - t358 * t374) * qJD(6)) * t341;
t567 = t111 * t269 + t36 * t329;
t563 = rSges(7,3) * t341;
t312 = -pkin(8) * t341 - t571;
t196 = t312 * t314;
t454 = -t320 - t484;
t411 = t454 + t458;
t396 = -t314 * t568 - t281 - t292 + t411;
t387 = -t196 + t396;
t33 = qJD(1) * t387 - t620;
t562 = t283 * t33;
t293 = t313 * rSges(6,3);
t310 = -rSges(6,1) * t342 + rSges(6,2) * t341;
t167 = t310 * t314 - t293;
t388 = -t167 + t396;
t445 = rSges(6,1) * t341 + rSges(6,2) * t342;
t56 = qJD(1) * t388 + t445 * t497 + t506;
t561 = t283 * t56;
t560 = t284 * rSges(5,3);
t559 = t284 * rSges(6,3);
t558 = t39 * t116;
t40 = t104 * t342 + t341 * t421;
t557 = t40 * t115;
t122 = pkin(5) * t398 - pkin(8) * t399;
t53 = t91 * rSges(7,1) + t90 * rSges(7,2) - rSges(7,3) * t399;
t556 = t122 + t53;
t162 = Icges(6,6) * t314 - t313 * t428;
t541 = t162 * t341;
t534 = t313 * t342;
t532 = t314 * t342;
t525 = t113 + t196;
t114 = t204 * rSges(7,1) + t203 * rSges(7,2) - rSges(7,3) * t535;
t198 = -pkin(5) * t534 - pkin(8) * t535;
t524 = t114 + t198;
t165 = Icges(6,5) * t314 - t313 * t431;
t523 = -t313 * t159 - t165 * t532;
t522 = -t165 * t534 + t616;
t520 = -t430 * t313 + t162;
t518 = t427 * t313 + t165;
t268 = t284 * qJ(4);
t490 = t313 * qJD(4);
t517 = t283 * pkin(3) - t268 - t278 + t490;
t258 = (rSges(7,1) * t357 + rSges(7,2) * t358) * t341;
t178 = qJD(6) * t258 + (t342 * t443 - t563) * qJD(5);
t282 = qJD(5) * t312;
t516 = -t178 - t282;
t514 = t228 + t311;
t512 = rSges(6,1) * t538 + t283 * rSges(6,3);
t482 = rSges(7,2) * t531;
t511 = -rSges(7,3) * t532 - t314 * t482;
t510 = -rSges(7,3) * t534 - t313 * t482;
t508 = -t313 * t334 - t314 * t356;
t194 = pkin(5) * t532 + pkin(8) * t533;
t507 = t284 * rSges(4,1) - t283 * rSges(4,2);
t218 = -t313 * rSges(4,1) - t314 * rSges(4,2);
t505 = -t427 + t431;
t504 = -t428 - t430;
t503 = qJD(1) * t345 + qJDD(2) * t573;
t491 = t425 * qJD(1);
t489 = -m(5) - m(6) - m(7);
t27 = -t104 * t533 - t605;
t487 = -t574 / 0.2e1;
t486 = t573 / 0.2e1;
t351 = t574 * pkin(2);
t483 = rSges(7,1) * t530;
t481 = m(4) - t489;
t474 = -t283 * t568 + t268 + t517 + t537;
t471 = t351 + t498;
t466 = qJD(1) * t573;
t465 = -t497 / 0.2e1;
t464 = t497 / 0.2e1;
t463 = -t496 / 0.2e1;
t462 = t496 / 0.2e1;
t461 = -t494 / 0.2e1;
t460 = t492 / 0.2e1;
t217 = -t313 * pkin(3) + qJ(4) * t314;
t213 = qJD(1) * t458;
t457 = t213 + t290 + t501;
t456 = qJD(1) * t217 - t609;
t453 = -rSges(7,1) * t528 + rSges(7,2) * t529;
t450 = t89 * rSges(7,1) + t88 * rSges(7,2);
t447 = t283 * rSges(4,1) + t284 * rSges(4,2);
t417 = t162 * t342 + t165 * t341;
t95 = Icges(6,4) * t398 + Icges(6,2) * t399 + Icges(6,6) * t283;
t97 = Icges(6,1) * t398 + Icges(6,4) * t399 + Icges(6,5) * t283;
t369 = qJD(5) * t417 + t341 * t95 - t342 * t97;
t94 = Icges(6,4) * t400 + Icges(6,2) * t401 + Icges(6,6) * t284;
t96 = Icges(6,1) * t400 + Icges(6,4) * t401 + Icges(6,5) * t284;
t370 = qJD(5) * t419 + t341 * t94 - t342 * t96;
t416 = -t165 * t342 + t541;
t92 = Icges(6,5) * t400 + Icges(6,6) * t401 + Icges(6,3) * t284;
t93 = Icges(6,5) * t398 + Icges(6,6) * t399 + Icges(6,3) * t283;
t442 = -(-t157 * t284 - t283 * t611 - t313 * t92 + t314 * t370) * t313 + (t159 * t284 + t283 * t416 - t313 * t93 + t314 * t369) * t314;
t441 = -(-t157 * t283 + t284 * t611 + t313 * t370 + t314 * t92) * t313 + (t159 * t283 - t284 * t416 + t313 * t369 + t314 * t93) * t314;
t440 = -t26 * t314 - t27 * t313;
t439 = -t28 * t314 - t29 * t313;
t438 = -t313 * t40 - t314 * t39;
t168 = -rSges(6,1) * t534 + rSges(6,2) * t535 + t314 * rSges(6,3);
t146 = -t217 + t508;
t433 = -qJD(1) * t146 - t456;
t384 = -qJD(1) * t168 + t433;
t57 = t445 * t496 - t384 - t490;
t437 = t313 * t56 + t314 * t57;
t59 = t162 * t533 + t523;
t436 = -t313 * t58 + t314 * t59;
t60 = -t157 * t314 - t160 * t535 + t163 * t534;
t61 = t162 * t535 + t522;
t435 = -t313 * t60 + t314 * t61;
t98 = rSges(6,1) * t400 + rSges(6,2) * t401 + t559;
t99 = rSges(6,1) * t470 + rSges(6,2) * t399 + t512;
t434 = -t313 * t99 - t314 * t98;
t432 = t471 + t508;
t420 = -t113 * t313 + t114 * t314;
t415 = t167 * t314 + t168 * t313;
t413 = -t341 * t430 - t342 * t427;
t410 = t603 + t454;
t409 = t283 * rSges(5,3) - t284 * t446;
t176 = t314 * rSges(5,3) + t313 * t446;
t166 = rSges(6,1) * t532 - rSges(6,2) * t533 + t293;
t408 = pkin(3) - t446;
t405 = -t310 + t334;
t397 = -t351 * t359 + t503;
t395 = t411 + t615;
t394 = -qJD(1) * t176 - t456;
t392 = -t480 - t485;
t391 = qJDD(1) * t498 - qJDD(2) * t574 + (-pkin(1) * t466 + t344 + t502) * qJD(1);
t325 = rSges(2,1) * t574 - rSges(2,2) * t573;
t322 = rSges(2,1) * t573 + rSges(2,2) * t574;
t386 = -t102 * t205 - t104 * t206 + t329 * t373;
t385 = -(Icges(7,5) * t201 + Icges(7,6) * t200) * t205 + (Icges(7,5) * t203 - Icges(7,6) * t204) * t206 + t254 * t329;
t383 = t347 + t393;
t382 = t341 * t518 + t342 * t520;
t381 = qJD(4) * t283 + qJDD(4) * t314 + t397;
t380 = -t484 + t615;
t377 = t341 * t385;
t366 = qJDD(1) * t351 - t359 * t484 + t391;
t363 = qJD(1) * t599 + qJD(4) * t284 + qJDD(1) * t217 - qJDD(4) * t313 + t366;
t361 = qJD(1) * (-t459 + t513) + qJDD(1) * t146 + t363;
t10 = qJD(1) * t122 + qJDD(1) * t198 + t269 * t114 - t115 * t228 - t206 * t178 - t208 * t311 - t282 * t496 + t329 * t53 + t361;
t121 = pkin(5) * t400 - pkin(8) * t401;
t52 = -rSges(7,3) * t401 + t450;
t11 = -t282 * t497 - t269 * t113 + t116 * t228 - t205 * t178 + t209 * t311 - t329 * t52 + (-t121 + t474) * qJD(1) + t387 * qJDD(1) + t381;
t362 = -qJD(1) * t198 - t114 * t329 + t206 * t228 + t311 * t496 + t433;
t34 = -t490 - t362;
t376 = t10 * t313 - t11 * t314 - t284 * t34 - t562;
t375 = t341 * t429 - t546;
t372 = t383 + t509;
t371 = (t341 * t504 + t342 * t505) * qJD(1);
t367 = (-Icges(7,2) * t204 + t110 + t187) * t206 - (Icges(7,2) * t200 - t108 + t186) * t205 + (Icges(7,2) * t530 + t226 + t326) * t329;
t276 = t428 * qJD(5);
t277 = t431 * qJD(5);
t365 = qJD(5) * t413 - t276 * t341 + t277 * t342;
t32 = t113 * t206 + t114 * t205 - qJD(3) + (t196 * t314 + t198 * t313) * qJD(5);
t364 = t32 * t420 + (t313 * t34 - t314 * t33) * t228;
t360 = t592 * t341;
t339 = rSges(3,3) * t467;
t279 = t310 * qJD(5);
t275 = t425 * qJD(5);
t252 = pkin(8) * t534;
t250 = pkin(8) * t532;
t229 = t453 - t563;
t220 = qJD(1) * t602 + t344;
t197 = pkin(5) * t535 - t252;
t195 = pkin(5) * t533 - t250;
t189 = t445 * t313;
t188 = t445 * t314;
t169 = qJD(1) * t410 + t344;
t154 = t313 * t483 + t510;
t153 = t314 * t483 + t511;
t152 = t375 * t313;
t151 = t375 * t314;
t150 = t374 * t313;
t149 = t374 * t314;
t144 = qJDD(1) * t499 + qJD(1) * (-rSges(3,1) * t466 + t339) + t391;
t143 = -qJD(1) * t278 + qJDD(1) * t602 - t359 * t499 + t503;
t130 = rSges(7,1) * t203 - rSges(7,2) * t204;
t129 = rSges(7,1) * t201 + rSges(7,2) * t200;
t119 = t313 * t412 - t180;
t117 = t119 * qJD(1);
t79 = -t490 - t394;
t78 = qJD(1) * t395 + t506;
t77 = qJD(1) * t507 + qJDD(1) * t218 + t366;
t76 = (-t278 + t447) * qJD(1) + t410 * qJDD(1) + t397;
t67 = qJD(5) * t415 - qJD(3);
t45 = -t275 * t314 - t283 * t424 - t284 * t412 + t313 * t365;
t44 = t275 * t313 + t283 * t412 - t284 * t424 + t314 * t365;
t43 = qJD(1) * t409 + qJDD(1) * t176 + t363;
t42 = (-t283 * t446 + t517 - t560) * qJD(1) + t395 * qJDD(1) + t381;
t38 = qJD(5) * t416 - t341 * t97 - t342 * t95;
t37 = -qJD(5) * t611 - t341 * t96 - t342 * t94;
t25 = qJD(5) * t434 - t167 * t208 + t168 * t209 + qJDD(3);
t24 = qJD(5) * t435 + t117;
t23 = qJD(5) * t436 + t613;
t22 = qJD(1) * t99 + qJDD(1) * t168 + t208 * t445 - t279 * t496 + t361;
t21 = -t279 * t497 - t209 * t445 + (-t98 + t474) * qJD(1) + t388 * qJDD(1) + t381;
t20 = -t172 * t535 + t173 * t203 + t174 * t204 + t226 * t91 + t373 * t399 - t374 * t90;
t19 = -t172 * t533 + t173 * t201 - t174 * t200 + t226 * t89 + t373 * t401 - t374 * t88;
t14 = t111 * t329 - t205 * t39 + t206 * t40;
t12 = t206 * t27 - t619;
t9 = -t115 * t113 + t114 * t116 - t121 * t496 - t122 * t497 - t208 * t196 + t198 * t209 - t205 * t53 - t206 * t52 + qJDD(3);
t6 = -t104 * t399 + t107 * t90 + t110 * t91 + t203 * t49 + t204 * t51 - t47 * t535;
t5 = t102 * t399 - t105 * t90 - t108 * t91 + t203 * t48 + t204 * t50 - t46 * t535;
t4 = -t104 * t401 + t107 * t88 + t110 * t89 - t200 * t51 + t201 * t49 - t47 * t533;
t3 = t102 * t401 - t105 * t88 - t108 * t89 - t200 * t50 + t201 * t48 - t46 * t533;
t2 = t115 * t29 + t116 * t28 + t20 * t329 - t205 * t5 + t206 * t6 + t269 * t66;
t1 = t115 * t27 + t116 * t26 + t19 * t329 - t205 * t3 + t206 * t4 - t269 * t604;
t15 = [(t405 * t561 - g(1) * (t166 + t372) + (t405 * t314 + t383) * t21 + (-rSges(6,2) * t539 + pkin(2) * t466 - t457 + t512 + (-t166 + t393) * qJD(1) + t594) * t57 + (t21 * (rSges(6,3) - t356) - t67 * (t166 + t167) * qJD(5)) * t313 + (-t384 - t559 + t595) * t56 + (t22 - g(2)) * (t168 + t432)) * m(6) + (t118 - t419) * t583 + (t119 - t417) * t584 + (t355 ^ 2 * Icges(5,2) + (Icges(5,1) * t354 + 0.2e1 * Icges(5,4) * t355) * t354 + m(2) * (t322 ^ 2 + t325 ^ 2) - t413 + Icges(2,3) + Icges(3,2) + Icges(4,3)) * qJDD(1) + t66 * t590 + (-g(1) * (-t320 + t380 + t458) + (t314 * t408 + t292 + t294 + t383) * t42 + (t409 - t457 + t502 + (-t380 + t393) * qJD(1) + t599) * t79 + (t283 * t408 - t268 - t394 - t560 + t609) * t78 + (-g(2) + t43) * (t176 + t217 + t471)) * m(5) + ((t144 - g(2)) * t598 + (t143 - g(1)) * (t347 + t349 + t392) + (-qJD(1) * t598 + t345) * t220 + (t220 + t339 + (t321 + t392) * qJD(1) + t597) * (qJD(1) * t499 + t278)) * m(3) + ((t28 + (-t102 * t313 + t104 * t314) * t341 + t605 + t624) * t206 + t12 + t619) * t588 + (t37 + t44 + t24) * t465 + t13 * t586 + t20 * t587 + (t19 + t13) * t585 + (t23 + ((t60 + (t157 - t416) * t314) * t314 + (t61 + (t157 - t541) * t313 + t616 - t522) * t313) * qJD(5) - t613) * t463 - t570 / 0.2e1 - t604 * t589 + (qJD(5) * t412 + t276 * t342 + t277 * t341) * qJD(1) + t569 / 0.2e1 + t567 + t557 / 0.2e1 + t558 / 0.2e1 + (t45 + t38) * t462 + ((t77 - g(2)) * (t471 + t218) + (t447 + t609) * t169 + (t169 + t507 + (t601 - t603) * qJD(1) + t597) * (qJD(1) * t218 - t609) + (t76 - g(1)) * (t383 + t603)) * m(4) + (t117 + ((t59 - t60 - t523) * t313 + t522 * t314) * qJD(5)) * t464 + (-g(1) * (t372 - t113 + t194) - (-t341 * t575 - t334 - t571) * t562 - t32 * (t194 + t196) * t497 + (-g(2) + t10) * (t432 + t524) + (t383 + t444 - t281 + (-t312 + t334 + t563) * t314) * t11 + (-t213 + t318 + t556 + (-t194 + t601) * qJD(1) + t594 + t620) * t34 + (-t450 + (t342 * t575 - t572) * t496 - t362 + t595) * t33) * m(7) - m(2) * (-g(1) * t322 + g(2) * t325); (-m(3) - t481) * (g(1) * t573 - g(2) * t574) + 0.2e1 * (t10 * t487 + t11 * t486) * m(7) + 0.2e1 * (t21 * t486 + t22 * t487) * m(6) + 0.2e1 * (t42 * t486 + t43 * t487) * m(5) + 0.2e1 * (t486 * t76 + t487 * t77) * m(4) + 0.2e1 * (t143 * t486 + t144 * t487) * m(3); (m(4) + m(5)) * qJDD(3) + m(6) * t25 + m(7) * t9 + t481 * g(3); t489 * (g(1) * t314 - g(2) * t313) + m(5) * (t283 * t78 + t284 * t79 - t313 * t43 + t314 * t42) - m(7) * t376 + (t21 * t314 - t22 * t313 + t284 * t57 + t561) * m(6) + (-m(5) * (t313 * t78 + t314 * t79) - m(7) * (t313 * t33 + t314 * t34) - m(6) * t437) * qJD(1); ((-t32 * t524 + t33 * t514) * t284 - (-t32 * t525 + t34 * t514) * t283 + (-t10 * t514 + t34 * t516 - t9 * t525 + t32 * (t121 + t52)) * t314 + (-t11 * t514 + t32 * t556 + t33 * t516 - t524 * t9) * t313 - g(1) * (-t252 + t510) - g(2) * (-t250 + t511) - g(3) * (t453 - t571) - (-g(3) * t575 + (g(1) * t313 + g(2) * t314) * (pkin(5) + t565)) * t341 - t33 * (-qJD(1) * t195 - t153 * t329 - t205 * t229 - t312 * t497) - t34 * (qJD(1) * t197 + t154 * t329 - t206 * t229 - t312 * t496) - t32 * (t153 * t206 + t154 * t205 + t195 * t496 + t197 * t497) - ((t113 * t33 - t114 * t34) * t341 + t364 * t342) * qJD(6)) * m(7) + ((t150 * t203 + t152 * t204) * t206 - (t149 * t203 + t151 * t204) * t205 + (t203 * t225 + t204 * t227) * t329 + (-t28 * t532 - t341 * t66) * qJD(6) + ((-qJD(6) * t29 + t386) * t342 + t360) * t313) * t588 + (-t25 * t415 + t67 * (t167 * t283 - t168 * t284 - t434) - t437 * t279 - (-t21 * t313 - t22 * t314 - t283 * t57 + t284 * t56) * t445 - (-t188 * t56 + t189 * t57) * qJD(1) - (t67 * (t188 * t314 + t189 * t313) - t437 * t310) * qJD(5) - g(1) * t189 - g(2) * t188 - g(3) * t310) * m(6) + (t314 * t460 + t580) * t12 + qJD(1) * (-t283 * t417 - t284 * t419 - t313 * t37 + t314 * t38) / 0.2e1 + qJDD(1) * (t313 * t419 - t314 * t417) / 0.2e1 + (-t28 * t313 + t29 * t314) * t590 + ((t180 * t497 + t491) * t313 + (t371 + (t382 * t314 + (-t181 - t593) * t313) * qJD(5)) * t314) * t464 + ((t181 * t496 - t491) * t314 + (t371 + (-t593 * t313 + (-t180 + t382) * t314) * qJD(5)) * t313) * t463 + (t283 * t59 + t284 * t58 + t442) * t465 + (t283 * t40 + t284 * t39 - t313 * t7 + t314 * t8) * t576 + (t26 * t284 + t27 * t283 - t3 * t313 + t314 * t4) * t585 + (t28 * t284 + t283 * t29 - t313 * t5 + t314 * t6) * t587 + (t283 * t61 + t284 * t60 + t441) * t462 + t23 * t580 + t24 * t581 + (-t313 * t39 + t314 * t40) * t582 + t436 * t583 + t435 * t584 + (-t26 * t313 + t27 * t314) * t589 + ((t150 * t201 - t152 * t200) * t206 - (t149 * t201 - t151 * t200) * t205 + (-t200 * t227 + t201 * t225) * t329 + (-t27 * t534 + t341 * t604) * qJD(6) + ((-qJD(6) * t26 + t386) * t342 + t360) * t314) * t586 + (((t150 * t357 - t152 * t358 - t104) * t206 - (t149 * t357 - t151 * t358 + t102) * t205 + (t225 * t357 - t227 * t358 + t373) * t329 - t111 * qJD(6)) * t341 + (qJD(6) * t438 - t592) * t342) * t577 + (qJD(1) * t45 + qJD(5) * t441 + qJDD(1) * t119 + t208 * t61 + t209 * t60 + t2) * t314 / 0.2e1 - (qJD(1) * t44 + qJD(5) * t442 + qJDD(1) * t118 + t208 * t59 + t209 * t58 + t1) * t313 / 0.2e1 - qJD(1) * ((t341 * t505 - t342 * t504) * qJD(1) + ((t313 * t519 - t314 * t518) * t342 + (-t313 * t521 + t314 * t520) * t341) * qJD(5)) / 0.2e1 + (t313 * t460 + t581) * t13 + t14 * t493 / 0.2e1; -t2 * t535 / 0.2e1 + (t341 * t439 + t342 * t66) * t590 + ((qJD(5) * t439 + t20) * t342 + (-qJD(5) * t66 - t28 * t283 + t284 * t29 - t313 * t6 - t314 * t5) * t341) * t587 - t1 * t533 / 0.2e1 + (t341 * t440 - t342 * t604) * t589 + ((qJD(5) * t440 + t19) * t342 + (qJD(5) * t604 - t26 * t283 + t27 * t284 - t3 * t314 - t313 * t4) * t341) * t585 - t14 * t495 / 0.2e1 + t342 * (t557 + t558 + t567 + t569 - t570) / 0.2e1 + (t111 * t342 + t341 * t438) * t582 + ((qJD(5) * t438 + t36) * t342 + (-qJD(5) * t111 - t283 * t39 + t284 * t40 - t313 * t8 - t314 * t7) * t341) * t576 + (t203 * t367 - t204 * t591 - t313 * t377) * t588 + (t200 * t591 + t201 * t367 - t314 * t377) * t586 + (t385 * t342 + (t367 * t357 + t358 * t591) * t341) * t577 + (t539 / 0.2e1 + t313 * t461) * t13 + (-t540 / 0.2e1 + t314 * t461) * t12 + ((qJD(5) * t364 + t10 * t114 - t11 * t113 - t33 * t52 + t34 * t53) * t342 + (t33 * (qJD(5) * t113 - t178 * t314) + t34 * (-qJD(5) * t114 + t178 * t313) - t9 * t420 + t32 * (t113 * t284 + t114 * t283 - t313 * t52 + t314 * t53) + t376 * t228) * t341 - t33 * (-t129 * t329 - t205 * t258) - t34 * (t130 * t329 - t206 * t258) - t32 * (t129 * t206 + t130 * t205) - g(1) * t130 - g(2) * t129 - g(3) * t258) * m(7);];
tau  = t15;
