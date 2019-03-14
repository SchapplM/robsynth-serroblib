% Calculate time derivative of joint inertia matrix for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR6_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR6_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR6_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:51:25
% EndTime: 2019-03-09 03:51:54
% DurationCPUTime: 19.40s
% Computational Cost: add. (52738->1045), mult. (44977->1457), div. (0->0), fcn. (42896->12), ass. (0->492)
t379 = sin(qJ(1));
t596 = t379 / 0.2e1;
t380 = cos(qJ(1));
t595 = -t380 / 0.2e1;
t622 = -qJD(1) / 0.2e1;
t369 = pkin(10) + qJ(3);
t361 = cos(t369);
t519 = qJD(3) * t361;
t480 = t519 / 0.2e1;
t359 = sin(t369);
t522 = qJD(1) * t379;
t494 = t359 * t522;
t621 = t380 * t480 - t494 / 0.2e1;
t521 = qJD(1) * t380;
t481 = t521 / 0.2e1;
t620 = t359 * t481 + t379 * t480;
t375 = cos(pkin(11));
t355 = t375 * pkin(4) + pkin(3);
t588 = pkin(3) - t355;
t482 = t588 * t361;
t572 = qJ(4) * t359;
t619 = t482 + t572;
t518 = qJD(3) * t379;
t487 = t361 * t518;
t393 = t359 * t521 + t487;
t618 = t359 * t588;
t580 = Icges(4,4) * t361;
t437 = -Icges(4,2) * t359 + t580;
t289 = Icges(4,6) * t379 + t380 * t437;
t581 = Icges(4,4) * t359;
t442 = Icges(4,1) * t361 - t581;
t291 = Icges(4,5) * t379 + t380 * t442;
t415 = t289 * t359 - t291 * t361;
t400 = t415 * t379;
t288 = -Icges(4,6) * t380 + t379 * t437;
t290 = -Icges(4,5) * t380 + t379 * t442;
t416 = t288 * t359 - t290 * t361;
t401 = t416 * t380;
t366 = t379 * rSges(4,3);
t563 = t359 * t380;
t617 = -rSges(4,2) * t563 + t366;
t378 = -pkin(7) - qJ(2);
t376 = cos(pkin(10));
t356 = pkin(2) * t376 + pkin(1);
t587 = rSges(4,1) * t361;
t460 = -rSges(4,2) * t359 + t587;
t408 = -t356 - t460;
t251 = (rSges(4,3) - t378) * t380 + t408 * t379;
t559 = t361 * t380;
t297 = rSges(4,1) * t559 + t617;
t475 = t380 * t356 - t379 * t378;
t252 = t297 + t475;
t616 = t251 * t380 + t252 * t379;
t432 = Icges(4,5) * t361 - Icges(4,6) * t359;
t286 = -Icges(4,3) * t380 + t379 * t432;
t523 = qJD(1) * t361;
t470 = -qJD(5) + t523;
t517 = qJD(3) * t380;
t488 = t359 * t517;
t615 = t379 * t470 + t488;
t370 = qJD(5) + qJD(6);
t473 = -t370 + t523;
t614 = t379 * t473 + t488;
t489 = t359 * t518;
t613 = t380 * t470 - t489;
t612 = t380 * t473 - t489;
t410 = rSges(3,1) * t376 - rSges(3,2) * sin(pkin(10)) + pkin(1);
t584 = rSges(3,3) + qJ(2);
t275 = t379 * t584 + t380 * t410;
t377 = -pkin(8) - qJ(4);
t367 = -pkin(9) + t377;
t526 = t367 - t377;
t477 = t526 * t359;
t373 = sin(pkin(11));
t551 = t379 * t373;
t351 = pkin(4) * t551;
t532 = -t355 * t559 - t351;
t368 = pkin(11) + qJ(5);
t360 = cos(t368);
t328 = pkin(5) * t360 + t355;
t358 = sin(t368);
t590 = pkin(5) * t358;
t591 = pkin(4) * t373;
t332 = t590 + t591;
t533 = t328 * t559 + t379 * t332;
t170 = -t380 * t477 + t532 + t533;
t362 = qJ(6) + t368;
t353 = sin(t362);
t354 = cos(t362);
t554 = t379 * t354;
t284 = -t353 * t559 + t554;
t555 = t379 * t353;
t285 = t354 * t559 + t555;
t203 = t285 * rSges(7,1) + t284 * rSges(7,2) + rSges(7,3) * t563;
t545 = t170 + t203;
t531 = t328 - t355;
t478 = t531 * t361;
t398 = -t359 * t367 + t478;
t558 = t373 * t380;
t352 = pkin(4) * t558;
t564 = t359 * t379;
t529 = t377 * t564 + t352;
t169 = -t332 * t380 + t379 * t398 + t529;
t282 = -t354 * t380 - t361 * t555;
t283 = -t353 * t380 + t361 * t554;
t452 = -t283 * rSges(7,1) - t282 * rSges(7,2);
t202 = rSges(7,3) * t564 - t452;
t546 = t169 + t202;
t611 = -t379 * t546 - t380 * t545;
t610 = 2 * m(4);
t609 = 2 * m(5);
t608 = 2 * m(6);
t607 = 2 * m(7);
t606 = m(5) / 0.2e1;
t605 = m(6) / 0.2e1;
t604 = m(7) / 0.2e1;
t435 = Icges(5,4) * t375 - Icges(5,2) * t373;
t277 = -Icges(5,6) * t361 + t359 * t435;
t603 = t277 / 0.2e1;
t440 = Icges(5,1) * t375 - Icges(5,4) * t373;
t278 = -Icges(5,5) * t361 + t359 * t440;
t602 = t278 / 0.2e1;
t550 = t379 * t375;
t312 = -t361 * t558 + t550;
t601 = t312 / 0.2e1;
t557 = t375 * t380;
t313 = t361 * t557 + t551;
t600 = t313 / 0.2e1;
t599 = -t361 / 0.2e1;
t598 = -t373 / 0.2e1;
t597 = t375 / 0.2e1;
t594 = t380 / 0.2e1;
t326 = rSges(4,1) * t359 + rSges(4,2) * t361;
t593 = m(4) * t326;
t592 = pkin(3) * t361;
t589 = qJD(1) / 0.2e1;
t586 = rSges(6,3) * t359;
t585 = rSges(7,3) * t359;
t583 = -rSges(5,3) - qJ(4);
t582 = -rSges(7,3) + t367;
t579 = Icges(6,4) * t358;
t578 = Icges(6,4) * t360;
t577 = Icges(7,4) * t353;
t576 = Icges(7,4) * t354;
t553 = t379 * t358;
t300 = -t360 * t380 - t361 * t553;
t552 = t379 * t360;
t301 = -t358 * t380 + t361 * t552;
t455 = -t301 * rSges(6,1) - t300 * rSges(6,2);
t214 = rSges(6,3) * t564 - t455;
t571 = t214 * t380;
t438 = Icges(7,1) * t354 - t577;
t264 = -Icges(7,5) * t361 + t359 * t438;
t568 = t264 * t354;
t439 = Icges(6,1) * t360 - t579;
t271 = -Icges(6,5) * t361 + t359 * t439;
t567 = t271 * t360;
t566 = t328 * t359;
t565 = t359 * t370;
t562 = t361 * t373;
t561 = t361 * t375;
t560 = t361 * t377;
t556 = t378 * t380;
t549 = qJ(4) + t377;
t474 = -t361 * t370 + qJD(1);
t414 = t379 * t474;
t176 = -t353 * t612 + t354 * t414;
t177 = t353 * t414 + t354 * t612;
t453 = t177 * rSges(7,1) + t176 * rSges(7,2);
t114 = rSges(7,3) * t393 + t453;
t486 = t361 * t517;
t548 = t114 * t563 + t202 * t486;
t413 = t380 * t474;
t174 = t353 * t614 + t354 * t413;
t175 = t353 * t413 - t354 * t614;
t502 = t175 * rSges(7,1) + t174 * rSges(7,2) + rSges(7,3) * t486;
t113 = -rSges(7,3) * t494 + t502;
t514 = qJD(5) * t360;
t511 = pkin(5) * t514;
t496 = t332 * t521 + t367 * t494 + t379 * t511;
t512 = qJD(5) * t590;
t530 = qJD(1) * t352 + t377 * t494;
t547 = -t531 * t488 + (-t531 * t522 + (-qJD(3) * t526 - t512) * t380) * t361 + t496 - t530 + t113;
t451 = rSges(7,1) * t354 - rSges(7,2) * t353;
t186 = (-rSges(7,1) * t353 - rSges(7,2) * t354) * t565 + (t361 * t451 + t585) * qJD(3);
t515 = qJD(5) * t359;
t485 = t358 * t515;
t189 = -pkin(5) * t485 + (t478 - t477) * qJD(3);
t544 = -t186 - t189;
t265 = -rSges(7,3) * t361 + t359 * t451;
t520 = qJD(3) * t359;
t543 = t203 * t520 + t265 * t494;
t348 = pkin(3) * t559;
t308 = qJ(4) * t563 + t348;
t409 = -t377 * t563 - t532;
t222 = t409 - t308;
t542 = -t222 - t308;
t236 = t313 * rSges(5,1) + t312 * rSges(5,2) + rSges(5,3) * t563;
t541 = -t236 - t308;
t237 = t359 * t531 + t361 * t526;
t540 = t237 + t265;
t148 = t361 * t202 + t265 * t564;
t450 = t572 + t592;
t293 = qJD(3) * t450 - qJD(4) * t361;
t539 = -(-t359 * t549 - t482) * qJD(3) - t293;
t257 = t361 * t549 - t618;
t325 = pkin(3) * t359 - qJ(4) * t361;
t309 = t325 * t522;
t538 = t257 * t522 + t309;
t537 = -t257 - t325;
t457 = rSges(5,1) * t375 - rSges(5,2) * t373;
t536 = -(rSges(5,3) * t359 + t361 * t457) * qJD(3) - t293;
t279 = -rSges(5,3) * t361 + t359 * t457;
t535 = -t279 - t325;
t307 = t450 * t379;
t534 = t379 * t307 + t380 * t308;
t516 = qJD(4) * t359;
t341 = t380 * t516;
t363 = qJD(2) * t379;
t528 = t341 + t363;
t364 = qJD(2) * t380;
t527 = t378 * t522 + t364;
t525 = t379 ^ 2 + t380 ^ 2;
t287 = Icges(4,3) * t379 + t380 * t432;
t524 = qJD(1) * t287;
t471 = -qJD(5) * t361 + qJD(1);
t411 = t380 * t471;
t192 = t358 * t615 + t360 * t411;
t193 = t358 * t411 - t360 * t615;
t392 = t486 - t494;
t123 = Icges(6,5) * t193 + Icges(6,6) * t192 + Icges(6,3) * t392;
t125 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t392;
t127 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t392;
t302 = -t358 * t559 + t552;
t303 = t360 * t559 + t553;
t209 = Icges(6,5) * t303 + Icges(6,6) * t302 + Icges(6,3) * t563;
t211 = Icges(6,4) * t303 + Icges(6,2) * t302 + Icges(6,6) * t563;
t213 = Icges(6,1) * t303 + Icges(6,4) * t302 + Icges(6,5) * t563;
t421 = -t211 * t358 + t213 * t360;
t35 = (qJD(3) * t421 - t123) * t361 + (qJD(3) * t209 - t125 * t358 + t127 * t360 + (-t211 * t360 - t213 * t358) * qJD(5)) * t359;
t430 = Icges(6,5) * t360 - Icges(6,6) * t358;
t204 = (-Icges(6,5) * t358 - Icges(6,6) * t360) * t515 + (Icges(6,3) * t359 + t361 * t430) * qJD(3);
t434 = -Icges(6,2) * t358 + t578;
t205 = (-Icges(6,2) * t360 - t579) * t515 + (Icges(6,6) * t359 + t361 * t434) * qJD(3);
t206 = (-Icges(6,1) * t358 - t578) * t515 + (Icges(6,5) * t359 + t361 * t439) * qJD(3);
t269 = -Icges(6,3) * t361 + t359 * t430;
t270 = -Icges(6,6) * t361 + t359 * t434;
t56 = t192 * t270 + t193 * t271 + t204 * t563 + t302 * t205 + t303 * t206 + t269 * t392;
t510 = t35 / 0.2e1 + t56 / 0.2e1;
t412 = t379 * t471;
t194 = -t358 * t613 + t360 * t412;
t195 = t358 * t412 + t360 * t613;
t124 = Icges(6,5) * t195 + Icges(6,6) * t194 + Icges(6,3) * t393;
t126 = Icges(6,4) * t195 + Icges(6,2) * t194 + Icges(6,6) * t393;
t128 = Icges(6,1) * t195 + Icges(6,4) * t194 + Icges(6,5) * t393;
t208 = Icges(6,5) * t301 + Icges(6,6) * t300 + Icges(6,3) * t564;
t210 = Icges(6,4) * t301 + Icges(6,2) * t300 + Icges(6,6) * t564;
t212 = Icges(6,1) * t301 + Icges(6,4) * t300 + Icges(6,5) * t564;
t422 = -t210 * t358 + t212 * t360;
t34 = (qJD(3) * t422 - t124) * t361 + (qJD(3) * t208 - t126 * t358 + t128 * t360 + (-t210 * t360 - t212 * t358) * qJD(5)) * t359;
t57 = t194 * t270 + t195 * t271 + t204 * t564 + t300 * t205 + t301 * t206 + t269 * t393;
t509 = t57 / 0.2e1 + t34 / 0.2e1;
t311 = t361 * t550 - t558;
t405 = t361 * t551 + t557;
t223 = Icges(5,5) * t311 - Icges(5,6) * t405 + Icges(5,3) * t564;
t508 = t223 * t564;
t507 = t223 * t563;
t224 = Icges(5,5) * t313 + Icges(5,6) * t312 + Icges(5,3) * t563;
t506 = t224 * t564;
t505 = t224 * t563;
t132 = t269 * t564 + t270 * t300 + t271 * t301;
t97 = -t208 * t361 + t359 * t422;
t504 = t97 / 0.2e1 + t132 / 0.2e1;
t133 = t269 * t563 + t302 * t270 + t303 * t271;
t98 = -t209 * t361 + t359 * t421;
t503 = t98 / 0.2e1 + t133 / 0.2e1;
t501 = t193 * rSges(6,1) + t192 * rSges(6,2) + rSges(6,3) * t486;
t454 = rSges(6,1) * t360 - rSges(6,2) * t358;
t207 = (-rSges(6,1) * t358 - rSges(6,2) * t360) * t515 + (t361 * t454 + t586) * qJD(3);
t500 = -t207 + t539;
t334 = qJ(4) * t486;
t339 = pkin(3) * t489;
t499 = t379 * (qJ(4) * t393 + qJD(1) * t348 + t379 * t516 - t339) + t380 * (-qJ(4) * t494 + t334 + t341 + (-t361 * t522 - t488) * pkin(3)) + t307 * t521;
t258 = qJD(1) * t405 + t373 * t488;
t259 = -qJD(1) * t311 - t375 * t488;
t498 = t259 * rSges(5,1) + t258 * rSges(5,2) + rSges(5,3) * t486;
t273 = -rSges(6,3) * t361 + t359 * t454;
t497 = -t273 + t537;
t215 = t303 * rSges(6,1) + t302 * rSges(6,2) + rSges(6,3) * t563;
t495 = t355 * t489 + t377 * t393;
t492 = t273 * t522;
t433 = -Icges(7,2) * t353 + t576;
t263 = -Icges(7,6) * t361 + t359 * t433;
t491 = t263 * t519;
t490 = t270 * t519;
t484 = t564 / 0.2e1;
t483 = t563 / 0.2e1;
t479 = t540 * t380;
t239 = t535 * t380;
t476 = -t328 * t361 - t356;
t472 = t361 * t512;
t469 = t361 * t114 + t186 * t564 + t265 * t393;
t468 = t539 + t544;
t221 = -t379 * t619 - t529;
t467 = t379 * t221 + t380 * t222 + t534;
t466 = t537 - t540;
t165 = t497 * t380;
t429 = Icges(7,5) * t354 - Icges(7,6) * t353;
t262 = -Icges(7,3) * t361 + t359 * t429;
t139 = -t262 * t361 + (-t263 * t353 + t568) * t359;
t121 = t262 * t564 + t263 * t282 + t264 * t283;
t196 = Icges(7,5) * t283 + Icges(7,6) * t282 + Icges(7,3) * t564;
t198 = Icges(7,4) * t283 + Icges(7,2) * t282 + Icges(7,6) * t564;
t200 = Icges(7,1) * t283 + Icges(7,4) * t282 + Icges(7,5) * t564;
t81 = t196 * t564 + t198 * t282 + t200 * t283;
t197 = Icges(7,5) * t285 + Icges(7,6) * t284 + Icges(7,3) * t563;
t199 = Icges(7,4) * t285 + Icges(7,2) * t284 + Icges(7,6) * t563;
t201 = Icges(7,1) * t285 + Icges(7,4) * t284 + Icges(7,5) * t563;
t82 = t197 * t564 + t199 * t282 + t201 * t283;
t449 = t379 * t81 + t380 * t82;
t41 = -t121 * t361 + t359 * t449;
t122 = t262 * t563 + t284 * t263 + t285 * t264;
t83 = t196 * t563 + t284 * t198 + t285 * t200;
t84 = t197 * t563 + t284 * t199 + t285 * t201;
t448 = t379 * t83 + t380 * t84;
t42 = -t122 * t361 + t359 * t448;
t424 = -t198 * t353 + t200 * t354;
t95 = -t196 * t361 + t359 * t424;
t423 = -t199 * t353 + t201 * t354;
t96 = -t197 * t361 + t359 * t423;
t444 = t379 * t95 + t380 * t96;
t108 = Icges(7,5) * t177 + Icges(7,6) * t176 + Icges(7,3) * t393;
t110 = Icges(7,4) * t177 + Icges(7,2) * t176 + Icges(7,6) * t393;
t112 = Icges(7,1) * t177 + Icges(7,4) * t176 + Icges(7,5) * t393;
t19 = t108 * t563 + t284 * t110 + t285 * t112 + t174 * t198 + t175 * t200 + t196 * t392;
t107 = Icges(7,5) * t175 + Icges(7,6) * t174 + Icges(7,3) * t392;
t109 = Icges(7,4) * t175 + Icges(7,2) * t174 + Icges(7,6) * t392;
t111 = Icges(7,1) * t175 + Icges(7,4) * t174 + Icges(7,5) * t392;
t20 = t107 * t563 + t284 * t109 + t285 * t111 + t174 * t199 + t175 * t201 + t197 * t392;
t183 = (-Icges(7,5) * t353 - Icges(7,6) * t354) * t565 + (Icges(7,3) * t359 + t361 * t429) * qJD(3);
t184 = (-Icges(7,2) * t354 - t577) * t565 + (Icges(7,6) * t359 + t361 * t433) * qJD(3);
t185 = (-Icges(7,1) * t353 - t576) * t565 + (Icges(7,5) * t359 + t361 * t438) * qJD(3);
t50 = t174 * t263 + t175 * t264 + t183 * t563 + t284 * t184 + t285 * t185 + t262 * t392;
t59 = t84 * t379 - t380 * t83;
t5 = (qJD(3) * t448 - t50) * t361 + (-qJD(1) * t59 + qJD(3) * t122 + t19 * t379 + t20 * t380) * t359;
t21 = t108 * t564 + t282 * t110 + t283 * t112 + t176 * t198 + t177 * t200 + t196 * t393;
t22 = t107 * t564 + t282 * t109 + t283 * t111 + t176 * t199 + t177 * t201 + t197 * t393;
t51 = t176 * t263 + t177 * t264 + t183 * t564 + t282 * t184 + t283 * t185 + t262 * t393;
t58 = t82 * t379 - t380 * t81;
t6 = (qJD(3) * t449 - t51) * t361 + (-qJD(1) * t58 + qJD(3) * t121 + t21 * t379 + t22 * t380) * t359;
t461 = t42 * t486 + t5 * t563 + t6 * t564 + (-t139 * t361 + t359 * t444) * t520 + t393 * t41;
t260 = qJD(1) * t312 + t373 * t489;
t261 = qJD(1) * t313 - t375 * t489;
t459 = -t261 * rSges(5,1) - t260 * rSges(5,2);
t458 = -t311 * rSges(5,1) + rSges(5,2) * t405;
t456 = t195 * rSges(6,1) + t194 * rSges(6,2);
t89 = t208 * t564 + t210 * t300 + t212 * t301;
t90 = t209 * t564 + t211 * t300 + t213 * t301;
t61 = t90 * t379 - t380 * t89;
t447 = t379 * t89 + t380 * t90;
t91 = t208 * t563 + t302 * t210 + t303 * t212;
t92 = t209 * t563 + t302 * t211 + t303 * t213;
t62 = t92 * t379 - t380 * t91;
t446 = t379 * t91 + t380 * t92;
t445 = t96 * t379 - t380 * t95;
t443 = t379 * t97 + t380 * t98;
t441 = Icges(4,1) * t359 + t580;
t436 = Icges(4,2) * t361 + t581;
t431 = Icges(5,5) * t375 - Icges(5,6) * t373;
t391 = t359 * t582 + t476;
t141 = (t332 - t378) * t380 + t391 * t379 + t452;
t142 = -t367 * t563 + t203 + t475 + t533;
t428 = t141 * t380 + t142 * t379;
t149 = -t361 * t203 - t265 * t563;
t427 = t148 * t380 + t149 * t379;
t406 = -t355 * t361 - t356 - t586;
t385 = t379 * t406 - t556;
t150 = t385 + t455 + t529;
t151 = t409 + t475 + t215;
t426 = t150 * t380 + t151 * t379;
t394 = t359 * t583 - t356 - t592;
t381 = t379 * t394 - t556;
t166 = t381 + t458;
t167 = t475 - t541;
t425 = t166 * t380 + t167 * t379;
t420 = -t215 * t379 + t571;
t419 = -t379 * t214 - t215 * t380;
t135 = t466 * t380;
t407 = t379 * (-qJ(4) * t487 + t339 + (-t380 * t619 + t351) * qJD(1) - t495) + t380 * (-t334 + (-t560 + t618) * t517 + t619 * t522 + t530) + t221 * t521 + t499;
t402 = qJD(3) * t326;
t399 = t169 * t380 - t379 * t545;
t397 = qJD(3) * t441;
t396 = qJD(3) * t436;
t395 = qJD(3) * (-Icges(4,5) * t359 - Icges(4,6) * t361);
t12 = qJD(1) * t448 - t19 * t380 + t20 * t379;
t13 = qJD(1) * t449 - t21 * t380 + t22 * t379;
t26 = (qJD(3) * t424 - t108) * t361 + (qJD(3) * t196 + (-t198 * t370 + t112) * t354 + (-t200 * t370 - t110) * t353) * t359;
t27 = (qJD(3) * t423 - t107) * t361 + (qJD(3) * t197 + (-t199 * t370 + t111) * t354 + (-t201 * t370 - t109) * t353) * t359;
t389 = t12 * t483 + t13 * t484 + t5 * t596 + t6 * t595 + (qJD(1) * t444 - t26 * t380 + t27 * t379) * t599 + t41 * t522 / 0.2e1 + t42 * t481 + t445 * t520 / 0.2e1 + t621 * t59 + t620 * t58;
t388 = rSges(4,2) * t494 + rSges(4,3) * t521 - t380 * t402;
t387 = -t361 * t183 + t262 * t520 + t519 * t568 + (t185 * t359 - t263 * t565) * t354;
t386 = -t361 * t204 + t269 * t520 + t519 * t567 + (t206 * t360 - t270 * t514) * t359;
t274 = -t379 * t410 + t380 * t584;
t137 = t139 * t520;
t60 = (-t491 + (-t264 * t370 - t184) * t359) * t353 + t387;
t9 = t137 + (qJD(3) * t444 - t60) * t361 + (-qJD(1) * t445 + t26 * t379 + t27 * t380) * t359;
t384 = -t361 * t9 - t42 * t494 + t461;
t383 = -t472 + (-t361 * t367 - t566) * qJD(3);
t382 = t137 + (t26 + t51) * t484 + (t27 + t50) * t483 + (t122 + t96) * t621 + (t121 + t95) * t620;
t317 = t460 * qJD(3);
t314 = t432 * qJD(3);
t296 = -rSges(4,3) * t380 + t379 * t460;
t268 = (Icges(5,5) * t359 + t361 * t440) * qJD(3);
t267 = (Icges(5,6) * t359 + t361 * t435) * qJD(3);
t249 = -qJD(1) * t275 + t364;
t248 = qJD(1) * t274 + t363;
t238 = t535 * t379;
t235 = rSges(5,3) * t564 - t458;
t230 = t379 * t395 + t524;
t229 = -qJD(1) * t286 + t380 * t395;
t228 = Icges(5,1) * t313 + Icges(5,4) * t312 + Icges(5,5) * t563;
t227 = Icges(5,1) * t311 - Icges(5,4) * t405 + Icges(5,5) * t564;
t226 = Icges(5,4) * t313 + Icges(5,2) * t312 + Icges(5,6) * t563;
t225 = Icges(5,4) * t311 - Icges(5,2) * t405 + Icges(5,6) * t564;
t182 = t202 * t563;
t179 = t326 * t518 + (t380 * t408 - t366) * qJD(1) + t527;
t178 = t363 + (-t556 + (-t356 - t587) * t379) * qJD(1) + t388;
t164 = t497 * t379;
t163 = t379 * t287 - t380 * t415;
t162 = t379 * t286 - t401;
t161 = -t287 * t380 - t400;
t160 = -t286 * t380 - t379 * t416;
t159 = Icges(5,1) * t261 + Icges(5,4) * t260 + Icges(5,5) * t393;
t158 = Icges(5,1) * t259 + Icges(5,4) * t258 + Icges(5,5) * t392;
t157 = Icges(5,4) * t261 + Icges(5,2) * t260 + Icges(5,6) * t393;
t156 = Icges(5,4) * t259 + Icges(5,2) * t258 + Icges(5,6) * t392;
t153 = -t361 * t215 - t273 * t563;
t152 = t214 * t361 + t273 * t564;
t145 = qJD(1) * t239 + t379 * t536;
t144 = t279 * t522 + t380 * t536 + t309;
t143 = -t269 * t361 + (-t270 * t358 + t567) * t359;
t140 = t143 * t520;
t138 = t420 * t359;
t136 = -t203 * t564 + t182;
t134 = t466 * t379;
t131 = t379 * t235 + t236 * t380 + t534;
t130 = rSges(6,3) * t393 + t456;
t129 = -rSges(6,3) * t494 + t501;
t120 = t339 + (t519 * t583 - t516) * t379 + t394 * t521 + t459 + t527;
t119 = -pkin(3) * t488 + qJD(1) * t381 + t334 + t498 + t528;
t105 = -t380 * t511 + t383 * t379 + ((t332 - t591) * t379 + t398 * t380) * qJD(1) + t495;
t102 = t312 * t226 + t313 * t228 + t505;
t101 = t312 * t225 + t313 * t227 + t507;
t100 = -t226 * t405 + t228 * t311 + t506;
t99 = -t225 * t405 + t227 * t311 + t508;
t94 = qJD(1) * t165 + t379 * t500;
t93 = t380 * t500 + t492 + t538;
t80 = (-rSges(6,3) * t519 - t516) * t379 + (t380 * t406 - t351) * qJD(1) - t456 + t495 + t527;
t79 = (-t355 * t359 - t560) * t517 + t385 * qJD(1) + t501 + t528 + t530;
t78 = -t359 * t479 - t361 * t545;
t77 = t169 * t361 + t237 * t564 + t148;
t76 = -t419 + t467;
t75 = (qJD(1) * t391 + t511) * t380 + (t472 - qJD(1) * t332 - t516 + (t361 * t582 + t566) * qJD(3)) * t379 - t453 + t527;
t74 = t383 * t380 + (-t556 + (t476 - t585) * t379) * qJD(1) + t496 + t502 + t528;
t73 = t359 * t399 + t182;
t72 = qJD(1) * t135 + t379 * t468;
t71 = t380 * t468 + t522 * t540 + t538;
t70 = (t273 * t518 + t130) * t361 + (-qJD(3) * t214 + t379 * t207 + t273 * t521) * t359;
t69 = (-t273 * t517 - t129) * t361 + (qJD(3) * t215 - t207 * t380 + t492) * t359;
t68 = -t202 * t520 + t469;
t67 = -t186 * t563 + (-t265 * t517 - t113) * t361 + t543;
t66 = (-t490 + (-qJD(5) * t271 - t205) * t359) * t358 + t386;
t64 = t467 - t611;
t63 = t379 * (rSges(5,3) * t487 - t459) + t380 * t498 + (t380 * t235 + t379 * t541) * qJD(1) + t499;
t47 = t420 * t519 + (qJD(1) * t419 - t129 * t379 + t130 * t380) * t359;
t45 = -t133 * t361 + t359 * t446;
t44 = -t132 * t361 + t359 * t447;
t43 = -t203 * t487 + (-t113 * t379 + (-t379 * t202 - t203 * t380) * qJD(1)) * t359 + t548;
t33 = (t237 * t518 + t105) * t361 + (-qJD(3) * t546 + t379 * t189 + t237 * t521) * t359 + t469;
t32 = (-qJD(3) * t479 - t547) * t361 + (qJD(3) * t170 + t237 * t522 + t380 * t544) * t359 + t543;
t31 = t123 * t564 + t300 * t125 + t301 * t127 + t194 * t211 + t195 * t213 + t209 * t393;
t30 = t124 * t564 + t300 * t126 + t301 * t128 + t194 * t210 + t195 * t212 + t208 * t393;
t29 = t123 * t563 + t302 * t125 + t303 * t127 + t192 * t211 + t193 * t213 + t209 * t392;
t28 = t124 * t563 + t302 * t126 + t303 * t128 + t192 * t210 + t193 * t212 + t208 * t392;
t25 = t129 * t380 + t379 * t130 + (t571 + (-t215 + t542) * t379) * qJD(1) + t407;
t18 = t399 * t519 + (qJD(1) * t611 + t105 * t380 - t547 * t379) * t359 + t548;
t17 = t547 * t380 + (t105 + t114) * t379 + (t546 * t380 + (t542 - t545) * t379) * qJD(1) + t407;
t16 = qJD(1) * t447 - t30 * t380 + t31 * t379;
t15 = qJD(1) * t446 - t28 * t380 + t29 * t379;
t8 = (qJD(3) * t447 - t57) * t361 + (-qJD(1) * t61 + qJD(3) * t132 + t30 * t379 + t31 * t380) * t359;
t7 = (qJD(3) * t446 - t56) * t361 + (-qJD(1) * t62 + qJD(3) * t133 + t28 * t379 + t29 * t380) * t359;
t1 = [t387 + t386 - t358 * t490 - t271 * t485 + 0.2e1 * m(3) * (t248 * t275 + t249 * t274) + (t178 * t252 + t179 * t251) * t610 + (t141 * t75 + t142 * t74) * t607 + (t150 * t80 + t151 * t79) * t608 + (t119 * t167 + t120 * t166) * t609 + (-t264 * t565 - t491) * t353 + (-Icges(5,3) * t361 - t436 + t442) * t520 + (-t184 * t353 - t205 * t358 - t267 * t373 + t268 * t375 + t431 * t520) * t359 + (-Icges(5,3) * t359 - t277 * t373 + t278 * t375 - t361 * t431 + t437 + t441) * t519; m(7) * (qJD(1) * t428 + t379 * t75 - t380 * t74) + m(6) * (qJD(1) * t426 + t379 * t80 - t380 * t79) + m(5) * (qJD(1) * t425 - t119 * t380 + t379 * t120) + m(4) * (qJD(1) * t616 - t178 * t380 + t379 * t179) + m(3) * (-t248 * t380 + t379 * t249 + (t274 * t380 + t275 * t379) * qJD(1)); 0; (-t260 * t277 / 0.2e1 - t261 * t278 / 0.2e1 + t405 * t267 / 0.2e1 - t311 * t268 / 0.2e1 - t51 / 0.2e1 - t26 / 0.2e1 + t314 * t594 - t509) * t380 + (t258 * t603 + t259 * t602 + t267 * t601 + t268 * t600 + t27 / 0.2e1 + t50 / 0.2e1 + t314 * t596 + t510) * t379 + m(4) * ((-t178 * t379 - t179 * t380) * t326 - t616 * t317) + m(5) * (t119 * t238 + t120 * t239 + t144 * t166 + t145 * t167) + m(7) * (t134 * t74 + t135 * t75 + t141 * t71 + t142 * t72) + m(6) * (t150 * t93 + t151 * t94 + t164 * t79 + t165 * t80) + ((Icges(5,5) * t261 / 0.2e1 + Icges(5,6) * t260 / 0.2e1 + Icges(5,3) * t393 / 0.2e1 + t289 * t622 + t396 * t596) * t380 + (-Icges(5,5) * t259 / 0.2e1 - Icges(5,6) * t258 / 0.2e1 - Icges(5,3) * t392 / 0.2e1 + t288 * t622 + t396 * t595) * t379) * t361 + ((-qJD(1) * t290 - t156 * t373 + t158 * t375 - t380 * t397) * t596 + (qJD(1) * t291 - t157 * t373 + t159 * t375 - t379 * t397) * t595) * t359 + ((t223 * t359 - t225 * t562 + t227 * t561) * t595 + t401 / 0.2e1 + (t224 * t359 - t226 * t562 + t228 * t561) * t596 - t400 / 0.2e1) * qJD(3) + ((-t252 * t593 + t277 * t601 + t278 * t600 + t122 / 0.2e1 + t96 / 0.2e1 + (-t224 / 0.2e1 + t289 / 0.2e1) * t361 + (t226 * t598 + t228 * t597 + t291 / 0.2e1) * t359 + t503) * t380 + (t251 * t593 + t121 / 0.2e1 - t405 * t603 + t311 * t602 + t95 / 0.2e1 + (-t223 / 0.2e1 + t288 / 0.2e1) * t361 + (t225 * t598 + t227 * t597 + t290 / 0.2e1) * t359 + t504) * t379) * qJD(1); m(5) * (t144 * t379 - t145 * t380 + (t238 * t379 + t239 * t380) * qJD(1)) + m(6) * (t93 * t379 - t380 * t94 + (t164 * t379 + t165 * t380) * qJD(1)) + m(7) * (t71 * t379 - t380 * t72 + (t134 * t379 + t135 * t380) * qJD(1)); (t134 * t72 + t135 * t71 + t17 * t64) * t607 + t379 * t12 - t380 * t13 + (t164 * t94 + t165 * t93 + t25 * t76) * t608 + t379 * t15 - t380 * t16 + (t131 * t63 + t144 * t239 + t145 * t238) * t609 + ((t379 * t296 + t297 * t380) * ((qJD(1) * t296 + t388) * t380 + (-t379 * t402 + (-t297 + t617) * qJD(1)) * t379) + t525 * t326 * t317) * t610 + t379 * ((t379 * t229 + (t162 + t400) * qJD(1)) * t379 + (t163 * qJD(1) + (t288 * t519 + t290 * t520) * t380 + (-t230 + (-t289 * t361 - t291 * t359) * qJD(3) + (t287 - t416) * qJD(1)) * t379) * t380) - t380 * ((t230 * t380 + (t161 + t401) * qJD(1)) * t380 + (t160 * qJD(1) + (-t289 * t519 - t291 * t520 + t524) * t379 + (-t229 + (t288 * t361 + t290 * t359) * qJD(3) - t415 * qJD(1)) * t380) * t379) + t379 * ((t312 * t156 + t313 * t158 + t258 * t226 + t259 * t228 + (t101 - t506) * qJD(1)) * t379 + (-t312 * t157 - t313 * t159 - t258 * t225 - t259 * t227 + (t102 + t508) * qJD(1)) * t380) - t380 * ((t405 * t157 - t311 * t159 - t260 * t225 - t261 * t227 + (t100 - t507) * qJD(1)) * t380 + (-t405 * t156 + t311 * t158 + t260 * t226 + t261 * t228 + (t99 + t505) * qJD(1)) * t379) + (t58 + t61 + (-t160 - t99) * t380 + (t100 + t161) * t379) * t522 + (t59 + t62 + (-t101 - t162) * t380 + (t102 + t163) * t379) * t521; 0.2e1 * (t425 * t606 + t426 * t605 + t428 * t604) * t519 + 0.2e1 * ((-t141 * t522 + t142 * t521 + t379 * t74 + t380 * t75) * t604 + (-t150 * t522 + t151 * t521 + t379 * t79 + t380 * t80) * t605 + (t119 * t379 + t120 * t380 - t166 * t522 + t167 * t521) * t606) * t359; 0; 0.2e1 * ((t134 * t518 + t135 * t517 - t17) * t604 + (t164 * t518 + t165 * t517 - t25) * t605 + (t238 * t518 + t239 * t517 - t63) * t606) * t361 + 0.2e1 * ((qJD(3) * t64 + t134 * t521 - t135 * t522 + t379 * t72 + t380 * t71) * t604 + (qJD(3) * t76 + t164 * t521 - t165 * t522 + t379 * t94 + t380 * t93) * t605 + (qJD(3) * t131 + t144 * t380 + t145 * t379 + t238 * t521 - t239 * t522) * t606) * t359; 0.4e1 * (t606 + t605 + t604) * (-0.1e1 + t525) * t359 * t519; (t510 * t380 + t509 * t379 + (-t379 * t503 + t380 * t504) * qJD(1)) * t359 + (-t60 - t66 + (t379 * t504 + t380 * t503) * qJD(3)) * t361 + t382 + m(7) * (t141 * t33 + t142 * t32 + t74 * t78 + t75 * t77) + m(6) * (t150 * t70 + t151 * t69 + t152 * t80 + t153 * t79) + t140; m(6) * (t70 * t379 - t380 * t69 + (t152 * t380 + t153 * t379) * qJD(1)) + m(7) * (-t32 * t380 + t33 * t379 + (t379 * t78 + t380 * t77) * qJD(1)); t389 + (qJD(3) * (t98 * t379 - t380 * t97) / 0.2e1 + t15 * t594 + t16 * t596 + (t61 * t594 - t379 * t62 / 0.2e1) * qJD(1)) * t359 + m(7) * (t134 * t32 + t135 * t33 + t17 * t73 + t18 * t64 + t71 * t77 + t72 * t78) + m(6) * (t138 * t25 + t152 * t93 + t153 * t94 + t164 * t69 + t165 * t70 + t47 * t76) + ((qJD(1) * t98 - t34) * t599 - t8 / 0.2e1 + t62 * t480 + t45 * t589) * t380 + ((qJD(1) * t97 + t35) * t599 + t7 / 0.2e1 + t61 * t480 + t44 * t589) * t379; 0.2e1 * ((t152 * t517 + t153 * t518 - t47) * t605 + (t517 * t77 + t518 * t78 - t18) * t604) * t361 + 0.2e1 * ((qJD(3) * t138 - t152 * t522 + t153 * t521 + t379 * t69 + t380 * t70) * t605 + (qJD(3) * t73 + t32 * t379 + t33 * t380 + t521 * t78 - t522 * t77) * t604) * t359; (t18 * t73 + t32 * t78 + t33 * t77) * t607 + (t138 * t47 + t152 * t70 + t153 * t69) * t608 + (t66 * t361 - t140 - t9 + (-t361 * t443 + t379 * t44 + t380 * t45) * qJD(3)) * t361 + (t380 * t7 + t379 * t8 - t361 * (t34 * t379 + t35 * t380) + (-t143 * t361 + t359 * t443) * qJD(3) + ((-t361 * t97 + t44) * t380 + (t361 * t98 - t42 - t45) * t379) * qJD(1)) * t359 + t461; -t60 * t361 + t382 + m(7) * (t141 * t68 + t142 * t67 + t148 * t75 + t149 * t74); m(7) * (qJD(1) * t427 + t68 * t379 - t380 * t67); t389 + m(7) * (t134 * t67 + t135 * t68 + t136 * t17 + t148 * t71 + t149 * t72 + t43 * t64); m(7) * ((qJD(3) * t427 - t43) * t361 + (qJD(3) * t136 + t379 * t67 + t380 * t68 + (-t148 * t379 + t149 * t380) * qJD(1)) * t359); m(7) * (t136 * t18 + t148 * t33 + t149 * t32 + t43 * t73 + t67 * t78 + t68 * t77) + t384; (t136 * t43 + t148 * t68 + t149 * t67) * t607 + t384;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;