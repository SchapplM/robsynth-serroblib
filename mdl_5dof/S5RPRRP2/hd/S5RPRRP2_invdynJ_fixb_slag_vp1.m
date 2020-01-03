% Calculate vector of inverse dynamics joint torques for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:45:07
% EndTime: 2020-01-03 11:45:27
% DurationCPUTime: 14.90s
% Computational Cost: add. (13128->497), mult. (9630->598), div. (0->0), fcn. (7393->8), ass. (0->281)
t272 = qJ(1) + pkin(8);
t262 = qJ(3) + t272;
t252 = sin(t262);
t253 = cos(t262);
t276 = cos(qJ(4));
t407 = t253 * t276;
t274 = sin(qJ(4));
t408 = t253 * t274;
t130 = Icges(6,4) * t407 - Icges(6,2) * t408 + Icges(6,6) * t252;
t132 = Icges(5,4) * t407 - Icges(5,2) * t408 + Icges(5,6) * t252;
t555 = t130 + t132;
t219 = Icges(6,5) * t276 - Icges(6,6) * t274;
t221 = Icges(5,5) * t276 - Icges(5,6) * t274;
t545 = t219 + t221;
t562 = Icges(5,5) + Icges(6,5);
t561 = -Icges(5,6) - Icges(6,6);
t560 = Icges(5,3) + Icges(6,3);
t209 = Icges(6,4) * t408;
t134 = Icges(6,1) * t407 + Icges(6,5) * t252 - t209;
t210 = Icges(5,4) * t408;
t136 = Icges(5,1) * t407 + Icges(5,5) * t252 - t210;
t559 = t134 + t136;
t263 = Icges(6,4) * t276;
t329 = -Icges(6,2) * t274 + t263;
t309 = t329 * t252;
t129 = -Icges(6,6) * t253 + t309;
t264 = Icges(5,4) * t276;
t330 = -Icges(5,2) * t274 + t264;
t131 = -Icges(5,6) * t253 + t252 * t330;
t556 = t129 + t131;
t435 = Icges(6,4) * t274;
t227 = Icges(6,1) * t276 - t435;
t133 = -Icges(6,5) * t253 + t227 * t252;
t436 = Icges(5,4) * t274;
t229 = Icges(5,1) * t276 - t436;
t312 = t229 * t252;
t135 = -Icges(5,5) * t253 + t312;
t554 = t133 + t135;
t222 = Icges(6,2) * t276 + t435;
t224 = Icges(5,2) * t276 + t436;
t553 = t222 + t224;
t558 = t545 * t252;
t557 = t555 * t274;
t477 = Icges(5,1) * t274 + t264;
t478 = Icges(6,1) * t274 + t263;
t552 = t477 + t478;
t506 = -t253 * t560 + t558;
t505 = t252 * t560 + t407 * t562 + t408 * t561;
t532 = -t276 * t559 + t557;
t410 = t252 * t276;
t508 = t554 * t410;
t551 = t559 * t410;
t507 = t556 * t408;
t218 = Icges(6,5) * t274 + Icges(6,6) * t276;
t220 = Icges(5,5) * t274 + Icges(5,6) * t276;
t550 = t218 + t220;
t544 = -t274 * t553 + t276 * t552;
t271 = qJD(1) + qJD(3);
t549 = t552 * qJD(4) - t271 * t562;
t548 = t553 * qJD(4) + t271 * t561;
t411 = t252 * t274;
t490 = -t253 * t506 - t411 * t556 + t508;
t489 = t253 * t505 + t411 * t555 - t551;
t488 = -t252 * t506 - t407 * t554 + t507;
t487 = t252 * t505 - t253 * t532;
t547 = (t329 + t330) * qJD(4);
t546 = (t227 + t229) * qJD(4);
t543 = t274 * t552 + t276 * t553;
t500 = t554 * t276;
t499 = t556 * t274;
t486 = t274 * t554 + t276 * t556;
t310 = t330 * t271;
t542 = t252 * t310 + t253 * t548 + t271 * t309;
t409 = t253 * t271;
t541 = -t252 * t548 + t253 * t310 + t329 * t409;
t311 = t227 * t271;
t540 = t252 * t311 + t253 * t549 + t271 * t312;
t539 = t229 * t409 - t252 * t549 + t253 * t311;
t273 = -qJ(5) - pkin(7);
t538 = rSges(6,3) - t273;
t150 = t218 * t253;
t152 = t220 * t253;
t510 = t252 * t544 - t150 - t152;
t414 = t220 * t252;
t416 = t218 * t252;
t509 = t407 * t552 - t408 * t553 + t414 + t416;
t260 = sin(t272);
t261 = cos(t272);
t179 = rSges(3,1) * t261 - rSges(3,2) * t260;
t277 = cos(qJ(1));
t269 = t277 * pkin(1);
t537 = t179 + t269;
t536 = -t499 + t500;
t485 = t274 * t559 + t276 * t555;
t535 = t550 * qJD(4) - t271 * t560;
t213 = rSges(6,2) * t408;
t268 = t276 * pkin(4);
t254 = t268 + pkin(3);
t534 = rSges(6,1) * t407 + t253 * t254 - t213;
t533 = -qJD(4) * t545 + t544 * t271;
t531 = -t487 * t252 - t253 * t488;
t530 = -t252 * t489 - t490 * t253;
t529 = qJD(4) * t543 - t271 * t550 + t274 * t547 - t276 * t546;
t528 = t510 * t271;
t447 = pkin(3) * t253;
t402 = t447 - t534 + (pkin(7) - t538) * t252;
t527 = -t252 * t535 - t271 * t536 + t409 * t545;
t526 = (-t532 + t558) * t271 + t535 * t253;
t480 = rSges(4,1) * t252 + rSges(4,2) * t253;
t144 = t480 * t271;
t275 = sin(qJ(1));
t267 = t275 * pkin(1);
t479 = pkin(2) * t260 + t267;
t313 = t479 * qJD(1);
t119 = t313 + t144;
t525 = t509 * t271;
t524 = qJD(4) * t485 - t271 * t505 - t274 * t542 + t276 * t540;
t523 = qJD(4) * t486 - t271 * t506 + t274 * t541 - t276 * t539;
t375 = qJD(4) * t252;
t203 = pkin(7) * t409;
t412 = t252 * t271;
t146 = pkin(3) * t412 - t203;
t373 = qJD(4) * t271;
t165 = -qJDD(4) * t252 - t253 * t373;
t444 = rSges(5,1) * t276;
t235 = -rSges(5,2) * t274 + t444;
t193 = t235 * qJD(4);
t232 = rSges(5,1) * t274 + rSges(5,2) * t276;
t270 = qJDD(1) + qJDD(3);
t426 = pkin(1) * qJDD(1);
t256 = t277 * t426;
t279 = qJD(1) ^ 2;
t425 = pkin(2) * qJDD(1);
t285 = t261 * t425 - t279 * t479 + t256;
t214 = rSges(5,2) * t408;
t140 = rSges(5,1) * t407 + rSges(5,3) * t252 - t214;
t176 = pkin(7) * t252 + t447;
t392 = t140 + t176;
t372 = qJD(4) * t274;
t359 = t253 * t372;
t405 = t271 * t276;
t368 = t252 * t405;
t305 = t359 + t368;
t371 = qJD(4) * t276;
t358 = t253 * t371;
t406 = t271 * t274;
t369 = t252 * t406;
t390 = -rSges(5,2) * t369 - rSges(5,3) * t409;
t93 = rSges(5,1) * t305 + rSges(5,2) * t358 + t390;
t27 = -t193 * t375 + t165 * t232 + (-t146 - t93) * t271 + t392 * t270 + t285;
t522 = -g(2) + t27;
t174 = rSges(4,1) * t253 - rSges(4,2) * t252;
t521 = -t144 * t271 + t174 * t270 - g(2) + t285;
t234 = rSges(6,1) * t276 - rSges(6,2) * t274;
t192 = t234 * qJD(4);
t441 = t276 * rSges(6,2);
t231 = rSges(6,1) * t274 + t441;
t237 = qJD(5) * t252;
t364 = t176 - t402;
t404 = t276 * qJD(4) ^ 2;
t230 = t253 * t273;
t363 = rSges(6,2) * t369 + rSges(6,3) * t409 + t237;
t446 = -pkin(4) * t359 - t203 - (t230 + (-pkin(3) + t254) * t252) * t271 - rSges(6,1) * t305 - rSges(6,2) * t358 + t363;
t15 = -t192 * t375 - qJDD(5) * t253 + t165 * t231 + (t165 * t274 - t252 * t404) * pkin(4) + t364 * t270 + (-t146 + t237 + t446) * t271 + t285;
t494 = g(2) - t15;
t344 = rSges(5,1) * t410 - rSges(5,2) * t411;
t138 = -rSges(5,3) * t253 + t344;
t166 = -qJDD(4) * t253 + t252 * t373;
t246 = t252 * pkin(3);
t175 = -pkin(7) * t253 + t246;
t251 = pkin(2) * t261;
t378 = t279 * t269 + t275 * t426;
t349 = t251 * t279 + t260 * t425 + t378;
t386 = pkin(3) * t409 + pkin(7) * t412;
t315 = t175 * t270 + t271 * t386 + t349;
t374 = qJD(4) * t253;
t367 = t253 * t406;
t304 = -t252 * t371 - t367;
t360 = t252 * t372;
t366 = t253 * t405;
t388 = rSges(5,1) * t366 + rSges(5,3) * t412;
t95 = -rSges(5,1) * t360 + rSges(5,2) * t304 + t388;
t26 = t138 * t270 - t166 * t232 + t193 * t374 + t271 * t95 + t315;
t520 = -g(3) + t26;
t145 = rSges(4,1) * t409 - rSges(4,2) * t412;
t519 = t145 * t271 + t270 * t480 - g(3) + t349;
t351 = pkin(4) * t274 + t231;
t504 = -rSges(6,1) * t410 - t252 * t254 - t230;
t106 = -rSges(6,2) * t411 - rSges(6,3) * t253 - t504;
t403 = -t175 + t106;
t370 = qJD(5) * t253;
t316 = rSges(6,1) * t366 + rSges(6,3) * t412 + t254 * t409 - t370;
t445 = (-pkin(4) * t372 - t271 * t273) * t252 - t386 - rSges(6,1) * t360 + rSges(6,2) * t304 + t316;
t14 = -qJDD(5) * t252 + t445 * t271 + t403 * t270 - t351 * t166 + (pkin(4) * t404 + qJD(4) * t192 - qJD(5) * t271) * t253 + t315;
t493 = g(3) - t14;
t518 = qJD(4) * t530 + t528;
t517 = qJD(4) * t531 - t525;
t516 = qJD(4) * t536 + t274 * t539 + t276 * t541;
t515 = qJD(4) * t532 + t274 * t540 + t276 * t542;
t514 = t252 * t533 + t253 * t529;
t513 = -t252 * t529 + t253 * t533;
t384 = t478 + t329;
t385 = t222 - t227;
t512 = (t274 * t384 + t276 * t385) * t271;
t382 = t477 + t330;
t383 = t224 - t229;
t511 = (t274 * t382 + t276 * t383) * t271;
t122 = t271 * t140;
t168 = t271 * t176;
t502 = -rSges(5,2) * t367 - t122 - t168 + t386 + t388;
t501 = t271 * t402 - t168 + t316;
t162 = t232 * t252;
t164 = t232 * t253;
t470 = -t232 * t374 - t271 * (t138 + t175);
t63 = t313 - t470;
t362 = t232 * t375;
t376 = qJD(1) * t261;
t442 = pkin(1) * qJD(1);
t381 = pkin(2) * t376 + t277 * t442;
t317 = -t362 + t381;
t64 = t271 * t392 + t317;
t498 = -t162 * t63 - t164 * t64;
t497 = -t351 * t375 - t370;
t496 = t252 * t527 - t253 * t523;
t495 = t252 * t526 + t253 * t524;
t492 = t252 * t523 + t253 * t527;
t491 = -t252 * t524 + t253 * t526;
t217 = pkin(4) * t408;
t482 = rSges(6,1) * t408 + rSges(6,2) * t407 + t217;
t481 = t234 + t268;
t178 = rSges(3,1) * t260 + rSges(3,2) * t261;
t471 = -t271 * (t175 + t403) + t237;
t457 = m(3) + m(4);
t455 = t165 / 0.2e1;
t454 = t166 / 0.2e1;
t448 = rSges(5,3) + pkin(7);
t415 = t219 * t271;
t413 = t221 * t271;
t401 = t252 * t478 + t129;
t400 = t253 * t478 + t130;
t399 = t252 * t477 + t131;
t398 = t253 * t477 + t132;
t397 = -t222 * t252 + t133;
t396 = -Icges(6,2) * t407 + t134 - t209;
t395 = -t224 * t252 + t135;
t394 = -Icges(5,2) * t407 + t136 - t210;
t379 = t251 + t269;
t377 = qJD(1) * t260;
t355 = -t375 / 0.2e1;
t354 = t375 / 0.2e1;
t353 = -t374 / 0.2e1;
t352 = t374 / 0.2e1;
t347 = -pkin(4) * t411 - t231 * t252;
t120 = t174 * t271 + t381;
t283 = t381 + t497;
t50 = t271 * t364 + t283;
t341 = t50 * t351;
t236 = rSges(2,1) * t277 - rSges(2,2) * t275;
t233 = rSges(2,1) * t275 + rSges(2,2) * t277;
t333 = -t252 * t64 + t253 * t63;
t322 = t138 * t252 + t140 * t253;
t314 = -t441 + (-rSges(6,1) - pkin(4)) * t274;
t306 = -pkin(2) * t377 - t275 * t442;
t299 = t252 * t314;
t290 = t274 * t397 + t276 * t401;
t289 = t274 * t396 + t276 * t400;
t288 = t274 * t395 + t276 * t399;
t287 = t274 * t394 + t276 * t398;
t108 = -t253 * t448 + t246 + t344;
t107 = t252 * t538 + t534;
t109 = -t214 + (pkin(3) + t444) * t253 + t448 * t252;
t284 = -rSges(5,1) * t368 - t146 - t390;
t282 = t498 * qJD(4);
t281 = -t509 * t165 / 0.2e1 - t485 * t455 + ((((t506 - t532) * t253 - t508 - t487) * t253 + ((t506 - t557) * t252 + (t499 + t500) * t253 - t507 + t488 + t551) * t252) * qJD(4) + t528) * t354 + (qJD(4) * t544 + t274 * t546 + t276 * t547) * t271 + (Icges(4,3) + t543) * t270 + (t510 + t486) * t454 + (t513 + t516) * t353 + ((((-t500 + t505) * t253 + t507 - t489) * t253 + ((t499 + t505) * t252 - t508 + t490) * t252) * qJD(4) + t517 + t525) * t352 + (t514 + t515 + t518) * t355;
t49 = t351 * t374 + t313 - t471;
t280 = (t50 * t504 + t49 * (-t252 * t273 - t213)) * t271 + (t253 * t314 * t50 + t299 * t49) * qJD(4);
t69 = qJD(4) * t322 + qJD(2);
t44 = qJD(2) + (t252 * t403 - t253 * t402) * qJD(4);
t25 = -t138 * t165 - t140 * t166 + qJDD(2) + (t252 * t95 - t253 * t93) * qJD(4);
t5 = qJDD(2) + t402 * t166 - t403 * t165 + (t252 * t445 + t253 * t446) * qJD(4);
t1 = [t281 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + (t521 * (t174 + t379) + t519 * (t479 + t480) + (-t120 + t145 + t381) * t119) * m(4) + ((qJDD(1) * t179 - g(2) + t256) * t537 + (-t279 * t537 + qJDD(1) * t178 - g(3) + t378 + (0.2e1 * rSges(3,1) * t376 - 0.2e1 * rSges(3,2) * t377 - qJD(1) * t179) * qJD(1)) * (t267 + t178)) * m(3) + (-g(2) * t236 - g(3) * t233 + (t233 ^ 2 + t236 ^ 2) * qJDD(1)) * m(2) + (t50 * (t306 + t363) + t280 - t494 * (t107 + t379) - t493 * (t106 + t479) + (t381 + t50 - t283 + t501) * t49) * m(6) + (t64 * (t284 + t306) + t282 + (-t317 + t64 + t381 + t502) * t63 + t522 * (t109 + t379) + t520 * (t108 + t479)) * m(5); t457 * qJDD(2) + m(5) * t25 + m(6) * t5 + (-m(5) - m(6) - t457) * g(1); t281 + (t341 * t374 + t280 + (t363 - t471) * t50 + (-t497 + t501) * t49 - t494 * t107 - t493 * t106) * m(6) + (t282 + (-t470 + t284) * t64 + (t362 + t502) * t63 + t522 * t109 + t520 * t108) * m(5) + (t119 * t145 - t120 * t144 + (-t119 * t271 + t521) * t174 + (t120 * t271 + t519) * t480) * m(4); t531 * t455 + t530 * t454 - (t514 * t271 - t509 * t270 + t488 * t166 + t487 * t165 + (t495 * t252 + t496 * t253) * qJD(4)) * t252 / 0.2e1 - (t513 * t271 + t510 * t270 + t490 * t166 + t489 * t165 + (t491 * t252 + t492 * t253) * qJD(4)) * t253 / 0.2e1 + (t485 * t252 - t486 * t253) * t270 / 0.2e1 - (((t382 + t384) * t276 + (-t383 - t385) * t274) * t271 + (((-t395 - t397) * t253 + (t394 + t396) * t252) * t276 + ((t399 + t401) * t253 + (-t398 - t400) * t252) * t274) * qJD(4)) * t271 / 0.2e1 + ((t271 * t485 - t516) * t253 + (t271 * t486 - t515) * t252) * t271 / 0.2e1 + t518 * t412 / 0.2e1 - t517 * t409 / 0.2e1 + ((-t271 * t487 + t496) * t253 + (t271 * t488 + t495) * t252) * t355 + ((t152 * t375 - t413) * t252 + (t511 + (-t288 * t253 + (-t414 + t287) * t252) * qJD(4)) * t253 + (t150 * t375 - t415) * t252 + (t512 + (-t290 * t253 + (-t416 + t289) * t252) * qJD(4)) * t253) * t354 + ((-t271 * t489 + t492) * t253 + (t271 * t490 + t491) * t252) * t353 + ((-t374 * t414 - t413) * t253 + (-t511 + (-t287 * t252 + (t152 + t288) * t253) * qJD(4)) * t252 + (-t374 * t416 - t415) * t253 + (-t512 + (-t289 * t252 + (t150 + t290) * t253) * qJD(4)) * t252) * t352 + (-t498 * t271 - (t69 * (-t162 * t252 - t164 * t253) + t333 * t235) * qJD(4) + t25 * t322 + t69 * ((t138 * t271 - t93) * t253 + (t95 - t122) * t252) + t333 * t193 + ((-t271 * t64 + t26) * t253 + (-t271 * t63 - t27) * t252) * t232 - g(1) * t235 + g(2) * t162 - g(3) * t164) * m(5) + (t14 * t217 + (-t5 * t402 + t44 * t446 + t14 * t231 + t49 * t192 + (t403 * t44 - t341) * t271) * t253 + (t5 * t403 + t44 * t445 - t15 * t351 + t50 * (-pkin(4) * t371 - t192) + (-t351 * t49 + t402 * t44) * t271) * t252 - (t347 * t49 - t482 * t50) * t271 - ((t234 * t49 - t44 * t482) * t253 + (t347 * t44 - t481 * t50) * t252) * qJD(4) - g(1) * t481 - g(3) * t482 - g(2) * t299) * m(6); (t252 * t493 + t253 * t494) * m(6);];
tau = t1;
