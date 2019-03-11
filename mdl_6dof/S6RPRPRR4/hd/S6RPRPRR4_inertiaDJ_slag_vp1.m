% Calculate time derivative of joint inertia matrix for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR4_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR4_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:44:28
% EndTime: 2019-03-09 03:45:00
% DurationCPUTime: 20.99s
% Computational Cost: add. (39391->901), mult. (41605->1238), div. (0->0), fcn. (39599->10), ass. (0->440)
t624 = Icges(5,1) + Icges(4,3);
t355 = cos(qJ(3));
t352 = sin(qJ(3));
t541 = Icges(5,6) * t352;
t551 = Icges(4,4) * t352;
t623 = -t541 - t551 + (-Icges(4,2) - Icges(5,3)) * t355;
t540 = Icges(5,6) * t355;
t550 = Icges(4,4) * t355;
t622 = -t540 - t550 + (-Icges(4,1) - Icges(5,2)) * t352;
t412 = Icges(4,5) * t355 - Icges(4,6) * t352;
t415 = Icges(5,4) * t355 - Icges(5,5) * t352;
t621 = t412 - t415;
t349 = qJ(1) + pkin(10);
t344 = cos(t349);
t343 = sin(t349);
t421 = Icges(4,1) * t355 - t551;
t244 = Icges(4,5) * t343 + t344 * t421;
t409 = Icges(5,2) * t355 - t541;
t247 = Icges(5,4) * t343 - t344 * t409;
t589 = t244 - t247;
t417 = -Icges(4,2) * t352 + t550;
t242 = Icges(4,6) * t343 + t344 * t417;
t407 = -Icges(5,3) * t352 + t540;
t245 = Icges(5,5) * t343 - t344 * t407;
t591 = t242 - t245;
t613 = -t352 * t591 + t355 * t589;
t620 = t344 * t613;
t241 = -Icges(4,6) * t344 + t343 * t417;
t584 = Icges(5,5) * t344 + t343 * t407;
t592 = t241 + t584;
t243 = -Icges(4,5) * t344 + t343 * t421;
t583 = Icges(5,4) * t344 + t343 * t409;
t590 = t243 + t583;
t615 = t343 * t624 + t621 * t344;
t619 = t623 * qJD(3);
t618 = t622 * qJD(3);
t586 = -t352 * t592 + t355 * t590;
t617 = t344 * t586;
t616 = t621 * t343 - t344 * t624;
t614 = ((Icges(5,5) - Icges(4,6)) * t355 + (Icges(5,4) - Icges(4,5)) * t352) * qJD(3);
t494 = qJD(3) * t352;
t463 = -t494 / 0.2e1;
t497 = qJD(1) * t355;
t475 = t343 * t497;
t612 = t344 * t463 - t475 / 0.2e1;
t499 = qJD(1) * t344;
t464 = t499 / 0.2e1;
t611 = t343 * t463 + t355 * t464;
t472 = t344 * t494;
t610 = t472 + t475;
t609 = t615 * qJD(1);
t354 = cos(qJ(5));
t351 = sin(qJ(5));
t548 = Icges(6,4) * t351;
t414 = Icges(6,2) * t354 + t548;
t279 = Icges(6,6) * t352 - t355 * t414;
t547 = Icges(6,4) * t354;
t419 = Icges(6,1) * t351 + t547;
t280 = Icges(6,5) * t352 - t355 * t419;
t608 = t279 * t354 + t280 * t351;
t350 = qJ(5) + qJ(6);
t346 = cos(t350);
t345 = sin(t350);
t546 = Icges(7,4) * t345;
t413 = Icges(7,2) * t346 + t546;
t266 = Icges(7,6) * t352 - t355 * t413;
t545 = Icges(7,4) * t346;
t418 = Icges(7,1) * t345 + t545;
t267 = Icges(7,5) * t352 - t355 * t418;
t607 = t266 * t346 + t267 * t345;
t606 = -t615 * t344 + t617 + (t613 + t616) * t343;
t605 = t343 * t589 + t344 * t590;
t604 = t343 * t591 + t344 * t592;
t603 = t343 ^ 2;
t602 = t344 ^ 2;
t601 = qJD(3) / 0.2e1;
t580 = 2 * m(4);
t563 = rSges(4,1) * t355;
t441 = -rSges(4,2) * t352 + t563;
t560 = rSges(4,3) * t344;
t252 = t343 * t441 - t560;
t526 = t344 * t355;
t527 = t344 * t352;
t588 = -rSges(4,2) * t527 + t343 * rSges(4,3);
t253 = rSges(4,1) * t526 + t588;
t321 = rSges(4,1) * t352 + rSges(4,2) * t355;
t380 = qJD(3) * t321;
t498 = qJD(1) * t352;
t476 = t343 * t498;
t362 = rSges(4,2) * t476 + rSges(4,3) * t499 - t344 * t380;
t91 = (qJD(1) * t252 + t362) * t344 + (-t343 * t380 + (-t253 + t588) * qJD(1)) * t343;
t600 = t580 * t91;
t523 = t346 * t352;
t263 = t343 * t523 + t344 * t345;
t524 = t345 * t352;
t264 = t343 * t524 - t344 * t346;
t436 = -t264 * rSges(7,1) - t263 * rSges(7,2);
t529 = t343 * t355;
t176 = rSges(7,3) * t529 - t436;
t338 = t344 * pkin(4);
t521 = t351 * t352;
t357 = -pkin(9) - pkin(8);
t564 = -pkin(8) - t357;
t373 = pkin(5) * t521 + t355 * t564;
t340 = pkin(5) * t354 + pkin(4);
t528 = t344 * t340;
t212 = t343 * t373 + t338 - t528;
t514 = t176 + t212;
t599 = t514 * t344;
t598 = -t343 * t586 + t344 * t616;
t595 = t343 * t615 + t620;
t594 = -qJD(1) * t616 + t344 * t614;
t593 = -t343 * t614 - t609;
t587 = t343 * rSges(5,1) - rSges(5,2) * t526;
t287 = t343 * pkin(4) + pkin(8) * t526;
t500 = qJD(1) * t343;
t261 = -t343 * t345 + t344 * t523;
t262 = t343 * t346 + t344 * t524;
t175 = t262 * rSges(7,1) + t261 * rSges(7,2) + rSges(7,3) * t526;
t484 = t344 * t521;
t384 = pkin(5) * t484 + t343 * t340 - t357 * t526;
t211 = t384 - t287;
t515 = t175 + t211;
t581 = -t343 * t514 - t344 * t515;
t579 = 2 * m(5);
t578 = 2 * m(6);
t577 = 2 * m(7);
t576 = m(5) / 0.2e1;
t575 = m(6) / 0.2e1;
t574 = m(7) / 0.2e1;
t573 = t343 / 0.2e1;
t572 = -t344 / 0.2e1;
t571 = t344 / 0.2e1;
t570 = t352 / 0.2e1;
t569 = -rSges(7,3) - pkin(3);
t568 = m(4) * t321;
t567 = pkin(3) * t355;
t566 = sin(qJ(1)) * pkin(1);
t347 = cos(qJ(1)) * pkin(1);
t565 = qJD(1) / 0.2e1;
t562 = rSges(5,1) * t344;
t561 = rSges(5,2) * t352;
t559 = rSges(6,3) * t352;
t558 = pkin(5) * qJD(5);
t557 = -rSges(5,3) - qJ(4);
t556 = rSges(7,3) - t357;
t410 = Icges(7,5) * t345 + Icges(7,6) * t346;
t265 = Icges(7,3) * t352 - t355 * t410;
t146 = t265 * t352 - t355 * t607;
t348 = qJD(5) + qJD(6);
t522 = t348 * t355;
t196 = (Icges(7,2) * t345 - t545) * t522 + (Icges(7,6) * t355 + t352 * t413) * qJD(3);
t195 = (-Icges(7,5) * t346 + Icges(7,6) * t345) * t522 + (Icges(7,3) * t355 + t352 * t410) * qJD(3);
t493 = qJD(3) * t355;
t425 = t345 * t266 * t522 + t352 * t195 + t265 * t493 + t494 * t607;
t197 = (-Icges(7,1) * t346 + t546) * t522 + (Icges(7,5) * t355 + t352 * t418) * qJD(3);
t525 = t345 * t197;
t555 = t146 * t493 + ((-t525 + (-t267 * t348 - t196) * t346) * t355 + t425) * t352;
t411 = Icges(6,5) * t351 + Icges(6,6) * t354;
t278 = Icges(6,3) * t352 - t355 * t411;
t152 = t278 * t352 - t355 * t608;
t491 = qJD(5) * t355;
t219 = (Icges(6,2) * t351 - t547) * t491 + (Icges(6,6) * t355 + t352 * t414) * qJD(3);
t220 = (-Icges(6,1) * t354 + t548) * t491 + (Icges(6,5) * t355 + t352 * t419) * qJD(3);
t218 = (-Icges(6,5) * t354 + Icges(6,6) * t351) * t491 + (Icges(6,3) * t355 + t352 * t411) * qJD(3);
t424 = t351 * t279 * t491 + t352 * t218 + t278 * t493 + t494 * t608;
t554 = t152 * t493 + ((-t351 * t220 + (-qJD(5) * t280 - t219) * t354) * t355 + t424) * t352;
t473 = t343 * t494;
t474 = t344 * t497;
t371 = -t473 + t474;
t458 = t348 * t352 + qJD(1);
t366 = -t345 * t458 + t346 * t493;
t457 = t348 + t498;
t390 = t344 * t457;
t159 = t343 * t366 + t346 * t390;
t367 = t345 * t493 + t346 * t458;
t160 = t343 * t367 + t345 * t390;
t437 = t160 * rSges(7,1) + t159 * rSges(7,2);
t104 = rSges(7,3) * t371 + t437;
t553 = t104 * t526 + t175 * t473;
t538 = qJ(4) * t352;
t537 = qJ(4) * t355;
t519 = t352 * t354;
t276 = t343 * t519 + t344 * t351;
t277 = t343 * t521 - t344 * t354;
t439 = -rSges(6,1) * t277 - rSges(6,2) * t276;
t194 = rSges(6,3) * t529 - t439;
t536 = t194 * t344;
t531 = t343 * t352;
t530 = t343 * t354;
t520 = t351 * t355;
t309 = pkin(8) * t473;
t469 = t357 * t494;
t470 = t351 * t493;
t128 = t343 * t469 + t309 + (qJD(5) * t276 + t343 * t470) * pkin(5) + ((-pkin(4) + t340) * t343 + t373 * t344) * qJD(1);
t518 = t104 + t128;
t391 = t343 * t457;
t161 = t344 * t366 - t346 * t391;
t162 = t344 * t367 - t345 * t391;
t516 = t162 * rSges(7,1) + t161 * rSges(7,2);
t105 = -rSges(7,3) * t610 + t516;
t332 = pkin(4) * t499;
t453 = qJD(5) + t498;
t387 = t351 * t453;
t483 = t344 * t519;
t423 = t340 * t499 + t357 * t475 + t483 * t558 + (t470 * pkin(5) + t469) * t344;
t129 = pkin(8) * t472 - t332 + (-pkin(5) * t387 + pkin(8) * t497) * t343 + t423;
t517 = t105 + t129;
t435 = rSges(7,1) * t345 + rSges(7,2) * t346;
t210 = (-rSges(7,1) * t346 + rSges(7,2) * t345) * t522 + (rSges(7,3) * t355 + t352 * t435) * qJD(3);
t269 = rSges(7,3) * t352 - t355 * t435;
t513 = t210 * t529 + t269 * t474;
t454 = qJD(5) * t352 + qJD(1);
t364 = -t351 * t454 + t354 * t493;
t184 = t344 * t364 - t453 * t530;
t365 = t354 * t454 + t470;
t185 = -t343 * t387 + t344 * t365;
t512 = t185 * rSges(6,1) + t184 * rSges(6,2);
t468 = t354 * t491;
t234 = -pkin(5) * t468 + qJD(3) * t373;
t511 = -t210 - t234;
t432 = t538 + t567;
t281 = t432 * t343;
t282 = pkin(3) * t526 + qJ(4) * t527;
t510 = t343 * t281 + t344 * t282;
t284 = -pkin(5) * t520 + t352 * t564;
t509 = t269 + t284;
t508 = -t282 - t287;
t285 = qJD(3) * t432 - qJD(4) * t355;
t434 = -rSges(5,2) * t355 + rSges(5,3) * t352;
t507 = -t434 * qJD(3) - t285;
t319 = pkin(3) * t352 - t537;
t286 = t319 * t500;
t506 = pkin(8) * t476 + t286;
t471 = t344 * t493;
t492 = qJD(4) * t352;
t505 = qJ(4) * t471 + t344 * t492;
t433 = rSges(5,3) * t355 + t561;
t504 = -t319 + t433;
t503 = t602 + t603;
t496 = qJD(3) * t343;
t495 = qJD(3) * t344;
t490 = -rSges(6,3) - pkin(3) - pkin(8);
t487 = t351 * t558;
t388 = t344 * t453;
t182 = t343 * t364 + t354 * t388;
t183 = t343 * t365 + t351 * t388;
t106 = Icges(6,5) * t183 + Icges(6,6) * t182 + Icges(6,3) * t371;
t108 = Icges(6,4) * t183 + Icges(6,2) * t182 + Icges(6,6) * t371;
t110 = Icges(6,1) * t183 + Icges(6,4) * t182 + Icges(6,5) * t371;
t188 = Icges(6,5) * t277 + Icges(6,6) * t276 + Icges(6,3) * t529;
t190 = Icges(6,4) * t277 + Icges(6,2) * t276 + Icges(6,6) * t529;
t192 = Icges(6,1) * t277 + Icges(6,4) * t276 + Icges(6,5) * t529;
t402 = t190 * t354 + t192 * t351;
t32 = (qJD(3) * t402 + t106) * t352 + (qJD(3) * t188 - t108 * t354 - t110 * t351 + (t190 * t351 - t192 * t354) * qJD(5)) * t355;
t59 = t182 * t279 + t183 * t280 + t218 * t529 + t219 * t276 + t220 * t277 + t278 * t371;
t486 = t32 / 0.2e1 + t59 / 0.2e1;
t107 = Icges(6,5) * t185 + Icges(6,6) * t184 - Icges(6,3) * t610;
t109 = Icges(6,4) * t185 + Icges(6,2) * t184 - Icges(6,6) * t610;
t111 = Icges(6,1) * t185 + Icges(6,4) * t184 - Icges(6,5) * t610;
t274 = -t343 * t351 + t483;
t275 = t484 + t530;
t187 = Icges(6,5) * t275 + Icges(6,6) * t274 + Icges(6,3) * t526;
t189 = Icges(6,4) * t275 + Icges(6,2) * t274 + Icges(6,6) * t526;
t191 = Icges(6,1) * t275 + Icges(6,4) * t274 + Icges(6,5) * t526;
t403 = t189 * t354 + t191 * t351;
t31 = (qJD(3) * t403 + t107) * t352 + (qJD(3) * t187 - t109 * t354 - t111 * t351 + (t189 * t351 - t191 * t354) * qJD(5)) * t355;
t60 = t184 * t279 + t185 * t280 + t218 * t526 + t219 * t274 + t220 * t275 - t278 * t610;
t485 = t60 / 0.2e1 + t31 / 0.2e1;
t131 = t274 * t279 + t275 * t280 + t278 * t526;
t92 = t187 * t352 - t355 * t403;
t482 = -t92 / 0.2e1 - t131 / 0.2e1;
t132 = t276 * t279 + t277 * t280 + t278 * t529;
t93 = t188 * t352 - t355 * t402;
t481 = -t93 / 0.2e1 - t132 / 0.2e1;
t310 = pkin(3) * t473;
t480 = t343 * (pkin(3) * t474 + t343 * t492 - t310 + (t343 * t493 + t344 * t498) * qJ(4)) + t344 * (-pkin(3) * t610 - qJ(4) * t476 + t505) + t281 * t499;
t193 = t275 * rSges(6,1) + t274 * rSges(6,2) + rSges(6,3) * t526;
t331 = pkin(7) * t499;
t479 = t331 + t505;
t478 = t344 * pkin(2) + t343 * pkin(7) + t347;
t438 = rSges(6,1) * t351 + rSges(6,2) * t354;
t283 = -t355 * t438 + t559;
t477 = t283 * t500;
t467 = t529 / 0.2e1;
t466 = t526 / 0.2e1;
t465 = t412 * t601 - t415 * qJD(3) / 0.2e1;
t337 = t344 * pkin(7);
t462 = t337 - t566;
t461 = -pkin(5) * t351 - qJ(4);
t460 = -t352 * pkin(8) - t319;
t85 = t187 * t529 + t189 * t276 + t191 * t277;
t86 = t188 * t529 + t190 * t276 + t192 * t277;
t428 = t343 * t86 + t344 * t85;
t43 = t352 * t132 + t355 * t428;
t459 = t352 * t93 + t43;
t233 = t504 * t344;
t456 = t352 * t105 + t175 * t493 + t269 * t610;
t101 = Icges(7,4) * t162 + Icges(7,2) * t161 - Icges(7,6) * t610;
t103 = Icges(7,1) * t162 + Icges(7,4) * t161 - Icges(7,5) * t610;
t169 = Icges(7,5) * t262 + Icges(7,6) * t261 + Icges(7,3) * t526;
t171 = Icges(7,4) * t262 + Icges(7,2) * t261 + Icges(7,6) * t526;
t173 = Icges(7,1) * t262 + Icges(7,4) * t261 + Icges(7,5) * t526;
t405 = t171 * t346 + t173 * t345;
t99 = Icges(7,5) * t162 + Icges(7,6) * t161 - Icges(7,3) * t610;
t25 = (qJD(3) * t405 + t99) * t352 + (qJD(3) * t169 + (-t173 * t348 - t101) * t346 + (t171 * t348 - t103) * t345) * t355;
t100 = Icges(7,4) * t160 + Icges(7,2) * t159 + Icges(7,6) * t371;
t102 = Icges(7,1) * t160 + Icges(7,4) * t159 + Icges(7,5) * t371;
t170 = Icges(7,5) * t264 + Icges(7,6) * t263 + Icges(7,3) * t529;
t172 = Icges(7,4) * t264 + Icges(7,2) * t263 + Icges(7,6) * t529;
t174 = Icges(7,1) * t264 + Icges(7,4) * t263 + Icges(7,5) * t529;
t404 = t172 * t346 + t174 * t345;
t98 = Icges(7,5) * t160 + Icges(7,6) * t159 + Icges(7,3) * t371;
t26 = (qJD(3) * t404 + t98) * t352 + (qJD(3) * t170 + (-t174 * t348 - t100) * t346 + (t172 * t348 - t102) * t345) * t355;
t124 = t263 * t266 + t264 * t267 + t265 * t529;
t76 = t169 * t529 + t171 * t263 + t173 * t264;
t77 = t170 * t529 + t172 * t263 + t174 * t264;
t430 = t343 * t77 + t344 * t76;
t39 = t124 * t352 + t355 * t430;
t89 = t169 * t352 - t355 * t405;
t90 = t170 * t352 - t355 * t404;
t426 = t343 * t89 - t344 * t90;
t427 = t343 * t90 + t89 * t344;
t19 = t101 * t263 + t103 * t264 + t159 * t171 + t160 * t173 + t169 * t371 + t529 * t99;
t20 = t100 * t263 + t102 * t264 + t159 * t172 + t160 * t174 + t170 * t371 + t529 * t98;
t51 = t343 * t76 - t344 * t77;
t54 = t159 * t266 + t160 * t267 + t195 * t529 + t196 * t263 + t197 * t264 + t265 * t371;
t5 = (-qJD(3) * t430 + t54) * t352 + (-qJD(1) * t51 + qJD(3) * t124 + t19 * t344 + t20 * t343) * t355;
t123 = t261 * t266 + t262 * t267 + t265 * t526;
t21 = t101 * t261 + t103 * t262 + t161 * t171 + t162 * t173 - t169 * t610 + t526 * t99;
t22 = t100 * t261 + t102 * t262 + t161 * t172 + t162 * t174 - t170 * t610 + t526 * t98;
t74 = t169 * t526 + t171 * t261 + t173 * t262;
t75 = t170 * t526 + t172 * t261 + t174 * t262;
t431 = t343 * t75 + t344 * t74;
t50 = t343 * t74 - t344 * t75;
t55 = t161 * t266 + t162 * t267 + t195 * t526 + t196 * t261 + t197 * t262 - t265 * t610;
t6 = (-qJD(3) * t431 + t55) * t352 + (-qJD(1) * t50 + qJD(3) * t123 + t21 * t344 + t22 * t343) * t355;
t455 = t5 * t529 + t6 * t526 + t39 * t474 + (t146 * t352 + t355 * t427) * t493 + t352 * (-t427 * t494 + (-qJD(1) * t426 + t25 * t344 + t26 * t343) * t355 + t555);
t288 = pkin(8) * t529 - t338;
t452 = t344 * t287 + t343 * t288 + t510;
t451 = rSges(5,1) * t499 + rSges(5,2) * t610 + rSges(5,3) * t471;
t38 = t123 * t352 + t355 * t431;
t83 = t187 * t526 + t189 * t274 + t191 * t275;
t84 = t188 * t526 + t190 * t274 + t192 * t275;
t429 = t343 * t84 + t344 * t83;
t42 = t352 * t131 + t355 * t429;
t450 = -t352 * t92 - t38 - t42;
t445 = -t283 + t460;
t444 = -pkin(8) * t493 - t285;
t443 = t352 * t557 - pkin(2);
t440 = rSges(6,1) * t183 + rSges(6,2) * t182;
t56 = t343 * t83 - t344 * t84;
t57 = t343 * t85 - t344 * t86;
t422 = t478 + t282;
t401 = -t193 * t344 - t194 * t343;
t400 = t193 * t343 - t536;
t389 = t460 - t509;
t254 = rSges(5,3) * t527 + t587;
t221 = (-rSges(6,1) * t354 + rSges(6,2) * t351) * t491 + (rSges(6,3) * t355 + t352 * t438) * qJD(3);
t386 = -t221 + t444;
t385 = -pkin(2) - t441;
t214 = t445 * t344;
t383 = t343 * (qJD(1) * t287 - t309) + t344 * (-pkin(8) * t610 + t332) + t288 * t499 + t480;
t382 = t352 * t461 - pkin(2);
t381 = t444 + t511;
t151 = t389 * t344;
t369 = t355 * t490 - pkin(2) - t538;
t368 = (rSges(5,2) - pkin(3)) * t355 + t443;
t12 = qJD(1) * t430 + t19 * t343 - t20 * t344;
t13 = qJD(1) * t431 + t21 * t343 - t22 * t344;
t363 = t12 * t467 + t13 * t466 + t5 * t572 + t6 * t573 + (qJD(1) * t427 + t25 * t343 - t26 * t344) * t570 + t39 * t500 / 0.2e1 + t38 * t464 + t426 * t493 / 0.2e1 + t611 * t51 + t612 * t50;
t361 = (-pkin(3) - t556) * t355 + t382;
t360 = t343 * t369 - t566;
t359 = t555 + (t26 + t54) * t467 + (t25 + t55) * t466 + (t123 + t89) * t612 + (t124 + t90) * t611;
t358 = (-t343 * t39 - t344 * t38) * t494 - t38 * t475 + t455;
t299 = t441 * qJD(3);
t255 = t343 * t434 - t562;
t232 = t504 * t343;
t231 = t269 * t529;
t217 = t253 + t478;
t216 = t343 * t385 + t462 + t560;
t213 = t445 * t343;
t168 = t254 + t422;
t167 = t343 * t368 + t462 + t562;
t166 = t352 * t175;
t164 = t176 * t526;
t156 = t321 * t496 + (-t347 + (-rSges(4,3) - pkin(7)) * t343 + t385 * t344) * qJD(1);
t155 = t331 + (-t566 + (-pkin(2) - t563) * t343) * qJD(1) + t362;
t154 = qJD(1) * t233 + t343 * t507;
t153 = t344 * t507 - t433 * t500 + t286;
t150 = t389 * t343;
t148 = t193 * t352 - t283 * t526;
t147 = -t194 * t352 + t283 * t529;
t144 = -t269 * t526 + t166;
t143 = -t176 * t352 + t231;
t134 = t422 + t193 + t287;
t133 = t337 + t338 + t360 + t439;
t130 = t254 * t344 + t255 * t343 + t510;
t127 = t400 * t355;
t126 = t310 + (-t492 + (t355 * t557 - t561) * qJD(3)) * t343 + (-t347 + (-rSges(5,1) - pkin(7)) * t343 + t368 * t344) * qJD(1);
t125 = -pkin(3) * t472 + (-t566 + (t443 - t567) * t343) * qJD(1) + t451 + t479;
t122 = t384 + t422 + t175;
t121 = t343 * t361 + t436 + t462 + t528;
t120 = qJD(1) * t214 + t343 * t386;
t119 = t344 * t386 + t477 + t506;
t118 = -t175 * t529 + t164;
t113 = -rSges(6,3) * t610 + t512;
t112 = rSges(6,3) * t371 + t440;
t95 = t211 * t352 - t509 * t526 + t166;
t94 = t284 * t529 - t352 * t514 + t231;
t88 = qJD(1) * t151 + t343 * t381;
t87 = t344 * t381 + t500 * t509 + t506;
t82 = -t401 + t452;
t73 = t309 + t310 + (-t492 + (-t537 + t559) * qJD(3)) * t343 + (-t347 + (-pkin(4) - pkin(7)) * t343 + t369 * t344) * qJD(1) - t440;
t72 = qJD(1) * t360 + t472 * t490 + t332 + t479 + t512;
t71 = t164 + (t212 * t344 - t343 * t515) * t355;
t69 = (-t283 * t496 - t112) * t352 + (-qJD(3) * t194 + t221 * t343 + t283 * t499) * t355;
t68 = (t283 * t495 + t113) * t352 + (qJD(3) * t193 - t221 * t344 + t477) * t355;
t67 = -t344 * t487 + t310 + ((-t354 * t558 - qJD(4)) * t352 + (t352 * t556 + t355 * t461) * qJD(3)) * t343 + (-t347 + (-pkin(7) - t340) * t343 + t361 * t344) * qJD(1) - t437;
t66 = -t343 * t487 + t569 * t472 + (-t566 + (t355 * t569 + t382) * t343) * qJD(1) + t423 + t479 + t516;
t65 = t452 - t581;
t63 = (qJD(1) * t255 + t451) * t344 + (t433 * t496 + (-t254 - t282 + t587) * qJD(1)) * t343 + t480;
t62 = -t104 * t352 + (-t176 * t355 - t269 * t531) * qJD(3) + t513;
t61 = -t210 * t526 + t456;
t44 = t400 * t494 + (qJD(1) * t401 + t112 * t344 - t113 * t343) * t355;
t41 = (t234 * t343 + t284 * t499) * t355 - t518 * t352 + (-t355 * t514 - t509 * t531) * qJD(3) + t513;
t40 = (t284 * t495 + t129) * t352 + (qJD(3) * t211 + t284 * t500 + t344 * t511) * t355 + t456;
t37 = -t176 * t472 + (-t105 * t343 + (-t175 * t344 - t176 * t343) * qJD(1)) * t355 + t553;
t33 = t112 * t343 + t113 * t344 + (t536 + (-t193 + t508) * t343) * qJD(1) + t383;
t30 = t106 * t526 + t108 * t274 + t110 * t275 + t184 * t190 + t185 * t192 - t188 * t610;
t29 = t107 * t526 + t109 * t274 + t111 * t275 + t184 * t189 + t185 * t191 - t187 * t610;
t28 = t106 * t529 + t108 * t276 + t110 * t277 + t182 * t190 + t183 * t192 + t188 * t371;
t27 = t107 * t529 + t109 * t276 + t111 * t277 + t182 * t189 + t183 * t191 + t187 * t371;
t18 = (t211 * t343 - t599) * t494 + (qJD(1) * t581 + t128 * t344 - t517 * t343) * t355 + t553;
t17 = t517 * t344 + t518 * t343 + (t599 + (t508 - t515) * t343) * qJD(1) + t383;
t16 = qJD(1) * t429 + t29 * t343 - t30 * t344;
t15 = qJD(1) * t428 + t27 * t343 - t28 * t344;
t9 = (-qJD(3) * t429 + t60) * t352 + (-qJD(1) * t56 + qJD(3) * t131 + t29 * t344 + t30 * t343) * t355;
t8 = (-qJD(3) * t428 + t59) * t352 + (-qJD(1) * t57 + qJD(3) * t132 + t27 * t344 + t28 * t343) * t355;
t1 = [t424 + t425 - t220 * t520 - t280 * t468 - t346 * t267 * t522 + (t155 * t217 + t156 * t216) * t580 + (t125 * t168 + t126 * t167) * t579 + (t133 * t73 + t134 * t72) * t578 + (t121 * t67 + t122 * t66) * t577 + (-t196 * t346 - t219 * t354 - t525) * t355 + (t421 + t409 + t623) * t494 + (t417 + t407 - t622) * t493; 0; 0; m(5) * (t125 * t232 + t126 * t233 + t153 * t167 + t154 * t168) + m(6) * (t119 * t133 + t120 * t134 + t213 * t72 + t214 * t73) + m(7) * (t121 * t87 + t122 * t88 + t150 * t66 + t151 * t67) + (m(4) * (-t156 * t321 - t216 * t299) - t26 / 0.2e1 - t54 / 0.2e1 + t465 * t344 - t486) * t344 + (m(4) * (-t155 * t321 - t217 * t299) + t25 / 0.2e1 + t55 / 0.2e1 + t465 * t343 + t485) * t343 + ((-qJD(3) * t591 + t344 * t618) * t573 + (-qJD(3) * t592 + t343 * t618) * t572 + (t572 * t589 - t573 * t590) * qJD(1)) * t352 + ((qJD(3) * t589 + t344 * t619) * t573 + (qJD(3) * t590 + t343 * t619) * t572 + (t572 * t591 - t573 * t592) * qJD(1)) * t355 + ((-t217 * t568 + t123 / 0.2e1 + t89 / 0.2e1 + (-t245 / 0.2e1 + t242 / 0.2e1) * t355 + (-t247 / 0.2e1 + t244 / 0.2e1) * t352 - t482) * t344 + (t216 * t568 + t124 / 0.2e1 + t90 / 0.2e1 + (t584 / 0.2e1 + t241 / 0.2e1) * t355 + (t583 / 0.2e1 + t243 / 0.2e1) * t352 - t481) * t343) * qJD(1); m(4) * t91 + m(5) * t63 + m(6) * t33 + m(7) * t17; (t150 * t88 + t151 * t87 + t65 * t17) * t577 + (t119 * t214 + t120 * t213 + t33 * t82) * t578 + (t130 * t63 + t153 * t233 + t154 * t232) * t579 + t503 * t321 * t299 * t580 + (t51 + t57) * t500 + (t50 + t56) * t499 + (t253 * t600 + t598 * t500 + t593 * t602 - t12 - t15 + (-t606 + t617) * t499) * t344 + (t13 + t16 + t252 * t600 + t594 * t603 + t595 * t499 + ((t593 - t609) * t343 + t594 * t344 + t605 * t494 + t604 * t493 + (-t352 * t605 - t355 * t604) * qJD(3) + ((t586 + t615) * t343 - t620 + t595 + t598) * qJD(1)) * t344 + (-t343 * t613 + t606) * t500) * t343; 0.2e1 * ((t121 * t344 + t122 * t343) * t574 + (t133 * t344 + t134 * t343) * t575 + (t167 * t344 + t168 * t343) * t576) * t493 + 0.2e1 * ((-t121 * t500 + t122 * t499 + t343 * t66 + t344 * t67) * t574 + (-t133 * t500 + t134 * t499 + t343 * t72 + t344 * t73) * t575 + (t125 * t343 + t126 * t344 - t167 * t500 + t168 * t499) * t576) * t352; (m(5) + m(6) + m(7)) * t494; 0.2e1 * ((t150 * t496 + t151 * t495 - t17) * t574 + (t213 * t496 + t214 * t495 - t33) * t575 + (t232 * t496 + t233 * t495 - t63) * t576) * t355 + 0.2e1 * ((qJD(3) * t65 + t150 * t499 - t151 * t500 + t343 * t88 + t344 * t87) * t574 + (qJD(3) * t82 + t119 * t344 + t120 * t343 + t213 * t499 - t214 * t500) * t575 + (qJD(3) * t130 + t153 * t344 + t154 * t343 + t232 * t499 - t233 * t500) * t576) * t352; 0.4e1 * (t576 + t575 + t574) * (-0.1e1 + t503) * t352 * t493; (t343 * t481 + t344 * t482) * t494 + (t485 * t344 + t486 * t343 + (t343 * t482 - t344 * t481) * qJD(1)) * t355 + t359 + m(6) * (t133 * t69 + t134 * t68 + t147 * t73 + t148 * t72) + m(7) * (t121 * t41 + t122 * t40 + t66 * t95 + t67 * t94) + t554; m(6) * t44 + m(7) * t18; ((t343 * t92 - t344 * t93) * t601 + t16 * t571 + t15 * t573 + (t57 * t571 - t343 * t56 / 0.2e1) * qJD(1)) * t355 + m(6) * (t119 * t147 + t120 * t148 - t127 * t33 + t213 * t68 + t214 * t69 + t44 * t82) + m(7) * (t150 * t40 + t151 * t41 + t17 * t71 + t18 * t65 + t87 * t94 + t88 * t95) + t363 + (t42 * t565 + (qJD(1) * t92 - t32) * t570 + t56 * t463 - t8 / 0.2e1) * t344 + (t43 * t565 + (qJD(1) * t93 + t31) * t570 + t57 * t463 + t9 / 0.2e1) * t343; 0.2e1 * ((t147 * t495 + t148 * t496 - t44) * t575 + (t495 * t94 + t496 * t95 - t18) * t574) * t355 + 0.2e1 * ((-qJD(3) * t127 - t147 * t500 + t148 * t499 + t343 * t68 + t344 * t69) * t575 + (qJD(3) * t71 + t343 * t40 + t344 * t41 + t499 * t95 - t500 * t94) * t574) * t352; (t18 * t71 + t40 * t95 + t41 * t94) * t577 + (-t127 * t44 + t147 * t69 + t148 * t68) * t578 + ((t450 * t344 + (-t39 - t459) * t343) * qJD(3) + t554) * t352 + (t344 * t9 + t343 * t8 + t352 * (t31 * t344 + t32 * t343) + (t152 * t352 + (t343 * t93 + t344 * t92) * t355) * qJD(3) + (t343 * t450 + t344 * t459) * qJD(1)) * t355 + t455; m(7) * (t121 * t62 + t122 * t61 + t143 * t67 + t144 * t66) + t359; m(7) * t37; m(7) * (t118 * t17 + t143 * t87 + t144 * t88 + t150 * t61 + t151 * t62 + t37 * t65) + t363; m(7) * ((-t37 + (t143 * t344 + t144 * t343) * qJD(3)) * t355 + (qJD(3) * t118 + t343 * t61 + t344 * t62 + (-t143 * t343 + t144 * t344) * qJD(1)) * t352); m(7) * (t118 * t18 + t143 * t41 + t144 * t40 + t37 * t71 + t61 * t95 + t62 * t94) + t358; (t118 * t37 + t143 * t62 + t144 * t61) * t577 + t358;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
