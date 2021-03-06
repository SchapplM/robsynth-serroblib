% Calculate time derivative of joint inertia matrix for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_inertiaDJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR4_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR4_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:04:48
% EndTime: 2019-03-09 07:05:17
% DurationCPUTime: 15.16s
% Computational Cost: add. (54160->870), mult. (40241->1210), div. (0->0), fcn. (37187->12), ass. (0->458)
t352 = sin(qJ(1));
t341 = t352 * rSges(5,3);
t354 = cos(qJ(1));
t344 = pkin(11) + qJ(3);
t335 = qJ(4) + t344;
t325 = sin(t335);
t326 = cos(t335);
t564 = rSges(5,1) * t326;
t434 = -rSges(5,2) * t325 + t564;
t259 = t354 * t434 + t341;
t328 = qJ(5) + t335;
t322 = sin(t328);
t345 = qJD(3) + qJD(4);
t334 = qJD(5) + t345;
t353 = cos(qJ(6));
t351 = sin(qJ(6));
t549 = Icges(7,4) * t353;
t415 = -Icges(7,2) * t351 + t549;
t323 = cos(t328);
t521 = t323 * t334;
t550 = Icges(7,4) * t351;
t139 = t415 * t521 + (Icges(7,6) * t334 + (-Icges(7,2) * t353 - t550) * qJD(6)) * t322;
t221 = -Icges(7,6) * t323 + t322 * t415;
t420 = Icges(7,1) * t353 - t550;
t222 = -Icges(7,5) * t323 + t322 * t420;
t596 = -t351 * t139 + (-t221 * t353 - t222 * t351) * qJD(6);
t289 = rSges(6,1) * t322 + rSges(6,2) * t323;
t385 = t289 * t334;
t332 = sin(t344);
t333 = cos(t344);
t555 = Icges(4,4) * t333;
t419 = -Icges(4,2) * t332 + t555;
t263 = Icges(4,6) * t352 + t354 * t419;
t556 = Icges(4,4) * t332;
t424 = Icges(4,1) * t333 - t556;
t265 = Icges(4,5) * t352 + t354 * t424;
t398 = t263 * t332 - t265 * t333;
t595 = t352 * t398;
t553 = Icges(5,4) * t326;
t417 = -Icges(5,2) * t325 + t553;
t255 = Icges(5,6) * t352 + t354 * t417;
t554 = Icges(5,4) * t325;
t422 = Icges(5,1) * t326 - t554;
t257 = Icges(5,5) * t352 + t354 * t422;
t400 = t255 * t325 - t257 * t326;
t594 = t352 * t400;
t551 = Icges(6,4) * t323;
t416 = -Icges(6,2) * t322 + t551;
t237 = Icges(6,6) * t352 + t354 * t416;
t552 = Icges(6,4) * t322;
t421 = Icges(6,1) * t323 - t552;
t239 = Icges(6,5) * t352 + t354 * t421;
t402 = t237 * t322 - t239 * t323;
t593 = t352 * t402;
t262 = -Icges(4,6) * t354 + t352 * t419;
t264 = -Icges(4,5) * t354 + t352 * t424;
t399 = t262 * t332 - t264 * t333;
t592 = t354 * t399;
t254 = -Icges(5,6) * t354 + t352 * t417;
t256 = -Icges(5,5) * t354 + t352 * t422;
t401 = t254 * t325 - t256 * t326;
t591 = t354 * t401;
t236 = -Icges(6,6) * t354 + t352 * t416;
t238 = -Icges(6,5) * t354 + t352 * t421;
t403 = t236 * t322 - t238 * t323;
t590 = t354 * t403;
t411 = Icges(7,5) * t353 - Icges(7,6) * t351;
t138 = t411 * t521 + (Icges(7,3) * t334 + (-Icges(7,5) * t351 - Icges(7,6) * t353) * qJD(6)) * t322;
t537 = t221 * t351;
t589 = -t334 * t537 - t138;
t505 = t353 * t354;
t507 = t352 * t351;
t281 = -t323 * t507 - t505;
t506 = t352 * t353;
t510 = t351 * t354;
t282 = t323 * t506 - t510;
t431 = -t282 * rSges(7,1) - t281 * rSges(7,2);
t523 = t322 * t352;
t194 = rSges(7,3) * t523 - t431;
t283 = -t323 * t510 + t506;
t284 = t323 * t505 + t507;
t522 = t322 * t354;
t195 = t284 * rSges(7,1) + t283 * rSges(7,2) + rSges(7,3) * t522;
t588 = -t352 * t194 - t354 * t195;
t350 = -pkin(7) - qJ(2);
t349 = cos(pkin(11));
t327 = t349 * pkin(2) + pkin(1);
t562 = rSges(4,2) * t332;
t565 = rSges(4,1) * t333;
t435 = -t562 + t565;
t389 = -t327 - t435;
t215 = (rSges(4,3) - t350) * t354 + t389 * t352;
t342 = t352 * rSges(4,3);
t490 = t354 * t565 + t342;
t216 = -t352 * t350 + (t327 - t562) * t354 + t490;
t587 = t215 * t354 + t216 * t352;
t412 = Icges(6,5) * t323 - Icges(6,6) * t322;
t234 = -Icges(6,3) * t354 + t352 * t412;
t586 = qJD(1) * t234;
t413 = Icges(5,5) * t326 - Icges(5,6) * t325;
t252 = -Icges(5,3) * t354 + t352 * t413;
t585 = qJD(1) * t252;
t414 = Icges(4,5) * t333 - Icges(4,6) * t332;
t260 = -Icges(4,3) * t354 + t352 * t414;
t439 = qJD(1) * t323 - qJD(6);
t515 = t334 * t354;
t468 = t322 * t515;
t584 = t352 * t439 + t468;
t287 = Icges(6,2) * t323 + t552;
t288 = Icges(6,1) * t322 + t551;
t397 = t287 * t322 - t288 * t323;
t582 = qJD(1) * t397 + t412 * t334;
t296 = Icges(5,2) * t326 + t554;
t297 = Icges(5,1) * t325 + t553;
t396 = t296 * t325 - t297 * t326;
t581 = qJD(1) * t396 + t413 * t345;
t392 = rSges(3,1) * t349 - rSges(3,2) * sin(pkin(11)) + pkin(1);
t557 = rSges(3,3) + qJ(2);
t245 = t352 * t557 + t354 * t392;
t343 = -pkin(8) + t350;
t487 = t343 - t350;
t307 = pkin(3) * t333 + t327;
t493 = t307 - t327;
t217 = t352 * t493 + t354 * t487;
t338 = -pkin(9) + t343;
t488 = t338 - t343;
t280 = pkin(4) * t326 + t307;
t495 = t280 - t307;
t196 = t352 * t495 + t354 * t488;
t580 = 2 * m(4);
t579 = 2 * m(5);
t578 = 2 * m(6);
t577 = 2 * m(7);
t346 = t352 ^ 2;
t347 = t354 ^ 2;
t576 = t352 / 0.2e1;
t575 = -t354 / 0.2e1;
t574 = -rSges(7,3) - pkin(10);
t105 = -t234 * t354 - t352 * t403;
t235 = Icges(6,3) * t352 + t354 * t412;
t106 = -t235 * t354 - t593;
t286 = Icges(6,5) * t322 + Icges(6,6) * t323;
t379 = t334 * t286;
t152 = -t354 * t379 - t586;
t485 = qJD(1) * t235;
t153 = -t352 * t379 + t485;
t528 = t287 * t334;
t154 = -qJD(1) * t236 - t354 * t528;
t527 = t288 * t334;
t156 = -qJD(1) * t238 - t354 * t527;
t157 = qJD(1) * t239 - t352 * t527;
t447 = t236 * t334 - t157;
t155 = qJD(1) * t237 - t352 * t528;
t449 = t238 * t334 + t155;
t524 = t322 * t334;
t13 = (t153 * t354 + (t106 + t590) * qJD(1)) * t354 + (t105 * qJD(1) + (-t154 * t322 + t156 * t323 - t237 * t521 - t239 * t524 + t485) * t352 + (-t152 + t447 * t323 + t449 * t322 + (-t234 - t402) * qJD(1)) * t354) * t352;
t440 = -qJD(6) * t323 + qJD(1);
t393 = t440 * t353;
t517 = t334 * t352;
t469 = t322 * t517;
t163 = t352 * t393 + (-t354 * t439 + t469) * t351;
t394 = t440 * t351;
t516 = t334 * t353;
t164 = t439 * t505 + (-t322 * t516 + t394) * t352;
t188 = Icges(7,5) * t282 + Icges(7,6) * t281 + Icges(7,3) * t523;
t190 = Icges(7,4) * t282 + Icges(7,2) * t281 + Icges(7,6) * t523;
t192 = Icges(7,1) * t282 + Icges(7,4) * t281 + Icges(7,5) * t523;
t481 = qJD(1) * t354;
t369 = t322 * t481 + t323 * t517;
t88 = Icges(7,5) * t164 + Icges(7,6) * t163 + Icges(7,3) * t369;
t90 = Icges(7,4) * t164 + Icges(7,2) * t163 + Icges(7,6) * t369;
t92 = Icges(7,1) * t164 + Icges(7,4) * t163 + Icges(7,5) * t369;
t19 = t163 * t190 + t164 * t192 + t188 * t369 + t281 * t90 + t282 * t92 + t523 * t88;
t189 = Icges(7,5) * t284 + Icges(7,6) * t283 + Icges(7,3) * t522;
t191 = Icges(7,4) * t284 + Icges(7,2) * t283 + Icges(7,6) * t522;
t193 = Icges(7,1) * t284 + Icges(7,4) * t283 + Icges(7,5) * t522;
t161 = t351 * t584 + t354 * t393;
t162 = -t353 * t584 + t354 * t394;
t482 = qJD(1) * t352;
t458 = t322 * t482;
t467 = t323 * t515;
t368 = -t458 + t467;
t87 = Icges(7,5) * t162 + Icges(7,6) * t161 + Icges(7,3) * t368;
t89 = Icges(7,4) * t162 + Icges(7,2) * t161 + Icges(7,6) * t368;
t91 = Icges(7,1) * t162 + Icges(7,4) * t161 + Icges(7,5) * t368;
t20 = t163 * t191 + t164 * t193 + t189 * t369 + t281 * t89 + t282 * t91 + t523 * t87;
t63 = t188 * t523 + t190 * t281 + t192 * t282;
t64 = t189 * t523 + t191 * t281 + t193 * t282;
t429 = t352 * t63 + t354 * t64;
t9 = qJD(1) * t429 - t19 * t354 + t20 * t352;
t573 = -t13 - t9;
t305 = rSges(4,1) * t332 + rSges(4,2) * t333;
t572 = m(4) * t305;
t560 = rSges(5,2) * t326;
t299 = rSges(5,1) * t325 + t560;
t571 = m(5) * t299;
t570 = m(6) * t289;
t569 = pkin(3) * t332;
t568 = pkin(4) * t325;
t567 = pkin(5) * t322;
t566 = pkin(5) * t323;
t563 = rSges(6,1) * t323;
t409 = -t190 * t351 + t192 * t353;
t23 = (t334 * t409 - t88) * t323 + (t188 * t334 - t351 * t90 + t353 * t92 + (-t190 * t353 - t192 * t351) * qJD(6)) * t322;
t559 = t23 * t354;
t408 = -t191 * t351 + t193 * t353;
t24 = (t334 * t408 - t87) * t323 + (t189 * t334 - t351 * t89 + t353 * t91 + (-t191 * t353 - t193 * t351) * qJD(6)) * t322;
t558 = t24 * t352;
t340 = t352 * rSges(6,3);
t452 = -t299 - t569;
t231 = t452 * t354;
t534 = t231 * t354;
t433 = -rSges(6,2) * t322 + t563;
t243 = t433 * t334;
t533 = t243 * t354;
t532 = t262 * t333;
t531 = t263 * t333;
t530 = t264 * t332;
t529 = t265 * t332;
t526 = t296 * t345;
t525 = t297 * t345;
t520 = t323 * t354;
t519 = t325 * t345;
t518 = t326 * t345;
t514 = t338 * t354;
t513 = t345 * t352;
t512 = t345 * t354;
t509 = t352 * t243;
t275 = t434 * t345;
t508 = t352 * t275;
t430 = rSges(7,1) * t353 - rSges(7,2) * t351;
t141 = t430 * t521 + (rSges(7,3) * t334 + (-rSges(7,1) * t351 - rSges(7,2) * t353) * qJD(6)) * t322;
t436 = pkin(10) * t322 + t566;
t504 = -t436 * t334 - t141;
t274 = t354 * t280;
t301 = t354 * t307;
t197 = -t352 * t488 + t274 - t301;
t503 = t352 * t196 + t354 * t197;
t267 = pkin(5) * t520 + pkin(10) * t522;
t502 = -t195 - t267;
t247 = rSges(6,1) * t520 - rSges(6,2) * t522 + t340;
t501 = -t197 - t247;
t218 = -t327 * t354 - t352 * t487 + t301;
t500 = t352 * t217 + t354 * t218;
t223 = -rSges(7,3) * t323 + t322 * t430;
t219 = t223 * t482;
t290 = -pkin(10) * t323 + t567;
t499 = t290 * t482 + t219;
t498 = -t223 - t290;
t246 = -rSges(6,3) * t354 + t352 * t433;
t160 = t352 * t246 + t354 * t247;
t258 = -rSges(5,3) * t354 + t352 * t434;
t173 = t352 * t258 + t354 * t259;
t457 = t325 * t482;
t309 = pkin(4) * t457;
t497 = t289 * t482 + t309;
t480 = qJD(3) * t332;
t474 = pkin(3) * t480;
t285 = -pkin(4) * t519 - t474;
t277 = t354 * t285;
t336 = qJD(2) * t352;
t496 = t277 + t336;
t494 = rSges(6,2) * t458 + rSges(6,3) * t481;
t492 = rSges(5,2) * t457 + rSges(5,3) * t481;
t491 = t343 * t482 + t352 * t474;
t319 = t338 * t482;
t337 = qJD(2) * t354;
t489 = t319 + t337;
t486 = t346 + t347;
t253 = Icges(5,3) * t352 + t354 * t413;
t484 = qJD(1) * t253;
t261 = Icges(4,3) * t352 + t354 * t414;
t483 = qJD(1) * t261;
t479 = qJD(3) * t333;
t478 = qJD(3) * t352;
t476 = pkin(4) * t518;
t475 = t354 * t562;
t473 = pkin(3) * t479;
t76 = -t188 * t323 + t322 * t409;
t220 = -Icges(7,3) * t323 + t322 * t411;
t97 = t220 * t523 + t221 * t281 + t222 * t282;
t472 = t76 / 0.2e1 + t97 / 0.2e1;
t77 = -t189 * t323 + t322 * t408;
t98 = t220 * t522 + t283 * t221 + t284 * t222;
t471 = t77 / 0.2e1 + t98 / 0.2e1;
t441 = t354 * t474;
t466 = t352 * (t352 * t285 + t481 * t495 - t319 + t491) + t354 * (-qJD(1) * t196 + t277 + t441) + t196 * t481;
t140 = t420 * t521 + (Icges(7,5) * t334 + (-Icges(7,1) * t351 - t549) * qJD(6)) * t322;
t465 = t322 * t353 * t140 + t323 * t222 * t516 + t220 * t524;
t370 = -t323 * t482 - t468;
t464 = t352 * (-t352 * t385 + (t354 * t433 + t340) * qJD(1)) + t354 * (rSges(6,1) * t370 - rSges(6,2) * t467 + t494) + t246 * t481;
t463 = t162 * rSges(7,1) + t161 * rSges(7,2) + rSges(7,3) * t467;
t384 = t299 * t345;
t462 = t352 * (qJD(1) * t259 - t352 * t384) + t354 * (-t512 * t560 + (-t325 * t512 - t326 * t482) * rSges(5,1) + t492) + t258 * t481;
t321 = t350 * t482;
t461 = t352 * (t481 * t493 + t321 - t491) + t354 * (-qJD(1) * t217 - t441) + t217 * t481;
t460 = -t197 + t502;
t459 = t309 + t499;
t456 = t332 * t482;
t455 = t521 / 0.2e1;
t454 = t482 / 0.2e1;
t453 = t481 / 0.2e1;
t451 = -t289 - t568;
t181 = t498 * t354;
t450 = t239 * t334 + t154;
t448 = -t237 * t334 + t156;
t176 = -qJD(1) * t254 - t354 * t526;
t446 = t257 * t345 + t176;
t177 = qJD(1) * t255 - t352 * t526;
t445 = t256 * t345 + t177;
t178 = -qJD(1) * t256 - t354 * t525;
t444 = -t255 * t345 + t178;
t179 = qJD(1) * t257 - t352 * t525;
t443 = t254 * t345 - t179;
t442 = -t352 * t338 + t274;
t81 = t160 + t503;
t266 = t436 * t352;
t82 = t352 * t266 + t354 * t267 - t588;
t438 = t498 - t568;
t437 = -t568 - t569;
t432 = t164 * rSges(7,1) + t163 * rSges(7,2);
t44 = t64 * t352 - t354 * t63;
t65 = t188 * t522 + t283 * t190 + t284 * t192;
t66 = t189 * t522 + t283 * t191 + t284 * t193;
t45 = t66 * t352 - t354 * t65;
t428 = t352 * t65 + t354 * t66;
t427 = t77 * t352 - t354 * t76;
t426 = t352 * t76 + t354 * t77;
t107 = t352 * t234 - t590;
t108 = t352 * t235 - t354 * t402;
t17 = t161 * t190 + t162 * t192 + t188 * t368 + t283 * t90 + t284 * t92 + t522 * t88;
t18 = t161 * t191 + t162 * t193 + t189 * t368 + t283 * t89 + t284 * t91 + t522 * t87;
t8 = qJD(1) * t428 - t17 * t354 + t18 * t352;
t425 = t44 * t482 + t45 * t481 + (-t105 * t482 - t107 * t481) * t354 + (t8 + (t108 * qJD(1) + (t155 * t322 - t157 * t323 + t236 * t521 + t238 * t524 - t586) * t354) * t354 + t106 * t482 + t108 * t481 + ((t107 + t593) * qJD(1) + (-t153 + t448 * t323 - t450 * t322 + (t235 - t403) * qJD(1)) * t354 + t352 * t152) * t352) * t352;
t423 = Icges(4,1) * t332 + t555;
t418 = Icges(4,2) * t333 + t556;
t295 = Icges(5,5) * t325 + Icges(5,6) * t326;
t387 = -t280 - t433;
t170 = (rSges(6,3) - t338) * t354 + t387 * t352;
t171 = t247 + t442;
t410 = t170 * t354 + t171 * t352;
t407 = t194 * t354 - t195 * t352;
t388 = -t307 - t434;
t205 = (rSges(5,3) - t343) * t354 + t388 * t352;
t206 = -t352 * t343 + t259 + t301;
t404 = t205 * t354 + t206 * t352;
t395 = -t476 + t504;
t391 = -t289 + t437;
t300 = pkin(10) * t467;
t95 = -rSges(7,3) * t458 + t463;
t96 = rSges(7,3) * t369 + t432;
t390 = (t194 + t266) * t481 + (pkin(5) * t370 - pkin(10) * t458 + t300 + t95) * t354 + (t96 + t369 * pkin(10) + (t323 * t481 - t469) * pkin(5)) * t352;
t149 = t438 * t354;
t386 = t464 + t466;
t51 = t82 + t503;
t383 = qJD(3) * t305;
t382 = t437 + t498;
t376 = t345 * t295;
t375 = -t473 - t476;
t374 = qJD(3) * t423;
t373 = qJD(3) * t418;
t372 = qJD(3) * (-Icges(4,5) * t332 - Icges(4,6) * t333);
t371 = t322 * t574 - t280 - t566;
t212 = t391 * t354;
t115 = -t252 * t354 - t352 * t401;
t116 = -t253 * t354 - t594;
t117 = t352 * t252 - t591;
t118 = t352 * t253 - t354 * t400;
t174 = -t354 * t376 - t585;
t175 = -t352 * t376 + t484;
t367 = t352 * ((t352 * t174 + (t117 + t594) * qJD(1)) * t352 + (t118 * qJD(1) + (t177 * t325 - t179 * t326 + t254 * t518 + t256 * t519 - t585) * t354 + (-t175 + t444 * t326 - t446 * t325 + (t253 - t401) * qJD(1)) * t352) * t354) + (-t115 * t354 + t116 * t352) * t482 + (-t117 * t354 + t118 * t352) * t481 + t425;
t366 = -t243 + t375;
t137 = t382 * t354;
t365 = t354 * t573 + t425;
t364 = t390 + t466;
t363 = t375 + t504;
t29 = t322 * t429 - t97 * t323;
t34 = t138 * t522 + t283 * t139 + t284 * t140 + t161 * t221 + t162 * t222 + t220 * t368;
t3 = (t334 * t428 - t34) * t323 + (-qJD(1) * t45 + t17 * t352 + t18 * t354 + t334 * t98) * t322;
t30 = t322 * t428 - t98 * t323;
t35 = t138 * t523 + t281 * t139 + t282 * t140 + t163 * t221 + t164 * t222 + t220 * t369;
t4 = (t334 * t429 - t35) * t323 + (-qJD(1) * t44 + t19 * t352 + t20 * t354 + t334 * t97) * t322;
t362 = t3 * t576 + t4 * t575 - t323 * (qJD(1) * t426 + t558 - t559) / 0.2e1 + t29 * t454 + t30 * t453 + t427 * t524 / 0.2e1 + t9 * t523 / 0.2e1 + t8 * t522 / 0.2e1 + (t354 * t455 - t458 / 0.2e1) * t45 + (t322 * t453 + t352 * t455) * t44;
t361 = rSges(4,2) * t456 + rSges(4,3) * t481 - t354 * t383;
t360 = t352 * t371 - t514;
t244 = -t352 * t392 + t354 * t557;
t15 = (t175 * t354 + (t116 + t591) * qJD(1)) * t354 + (t115 * qJD(1) + (-t176 * t325 + t178 * t326 - t255 * t518 - t257 * t519 + t484) * t352 + (-t174 + t443 * t326 + t445 * t325 + (-t252 - t400) * qJD(1)) * t354) * t352;
t359 = (-t15 + t573) * t354 + t367;
t241 = t416 * t334;
t242 = t421 * t334;
t358 = qJD(1) * t286 + (t242 - t528) * t323 + (-t241 - t527) * t322;
t269 = t417 * t345;
t270 = t422 * t345;
t357 = qJD(1) * t295 + (t270 - t526) * t326 + (-t269 - t525) * t325;
t356 = -t559 / 0.2e1 + t558 / 0.2e1 + (t322 * t448 + t323 * t450 + t352 * t582 + t358 * t354 + t34) * t576 + (-t322 * t447 + t323 * t449 + t358 * t352 - t354 * t582 + t35) * t575 + (t236 * t323 + t238 * t322 - t286 * t354 - t352 * t397 + t76 + t97) * t454 + (t237 * t323 + t239 * t322 + t352 * t286 - t354 * t397 + t77 + t98) * t453;
t355 = t356 + (t325 * t444 + t326 * t446 + t352 * t581 + t357 * t354) * t576 + (-t325 * t443 + t326 * t445 + t357 * t352 - t354 * t581) * t575 + (t254 * t326 + t256 * t325 - t295 * t354 - t352 * t396) * t454 + (t255 * t326 + t257 * t325 + t352 * t295 - t354 * t396) * t453;
t316 = pkin(3) * t456;
t294 = t435 * qJD(3);
t272 = -t475 + t490;
t271 = -rSges(4,3) * t354 + t352 * t435;
t230 = t452 * t352;
t225 = t451 * t354;
t224 = t451 * t352;
t214 = -qJD(1) * t245 + t337;
t213 = qJD(1) * t244 + t336;
t211 = t391 * t352;
t199 = t352 * t372 + t483;
t198 = -qJD(1) * t260 + t354 * t372;
t180 = t498 * t352;
t148 = t438 * t352;
t147 = -t299 * t481 - t508 + (-t332 * t481 - t333 * t478) * pkin(3);
t146 = t299 * t482 + t316 + (-t275 - t473) * t354;
t145 = t321 + t337 + t305 * t478 + (t354 * t389 - t342) * qJD(1);
t144 = t336 + (-t350 * t354 + (-t327 - t565) * t352) * qJD(1) + t361;
t136 = t382 * t352;
t131 = -t289 * t481 - t509 + (-t325 * t481 - t326 * t513) * pkin(4);
t130 = (-t243 - t476) * t354 + t497;
t126 = t352 * t261 - t354 * t398;
t125 = t352 * t260 - t592;
t124 = -t261 * t354 - t595;
t123 = -t260 * t354 - t352 * t399;
t120 = t337 + t299 * t513 + (t354 * t388 - t341) * qJD(1) + t491;
t119 = t336 + (-t307 - t564) * t482 + (-qJD(1) * t343 - t384 - t474) * t354 + t492;
t114 = qJD(1) * t212 + t352 * t366;
t113 = t354 * t366 + t316 + t497;
t112 = -t323 * t195 - t223 * t522;
t111 = t194 * t323 + t223 * t523;
t110 = t442 - t502;
t109 = t360 + t431;
t104 = (-t285 + t385) * t352 + (t354 * t387 - t340) * qJD(1) + t489;
t103 = -t354 * t385 + (-t514 + (-t280 - t563) * t352) * qJD(1) + t494 + t496;
t102 = -t220 * t323 + (t222 * t353 - t537) * t322;
t101 = t407 * t322;
t100 = t102 * t524;
t99 = t173 + t500;
t84 = qJD(1) * t181 + t352 * t504;
t83 = t354 * t504 + t499;
t80 = -t259 * t482 + t462;
t79 = qJD(1) * t149 + t352 * t395;
t78 = t354 * t395 + t459;
t73 = -t247 * t482 + t464;
t68 = qJD(1) * t137 + t352 * t363;
t67 = t354 * t363 + t316 + t459;
t56 = t81 + t500;
t53 = (-t285 + (t323 * t574 + t567) * t334) * t352 + t371 * t481 - t432 + t489;
t52 = -pkin(5) * t468 + qJD(1) * t360 + t300 + t463 + t496;
t50 = t51 + t500;
t49 = (-t218 - t259) * t482 + t461 + t462;
t48 = (t223 * t517 + t96) * t323 + (t352 * t141 - t194 * t334 + t223 * t481) * t322;
t47 = (-t223 * t515 - t95) * t323 + (-t141 * t354 + t195 * t334 + t219) * t322;
t37 = t322 * t596 + t589 * t323 + t465;
t36 = t482 * t501 + t386;
t31 = t407 * t521 + (qJD(1) * t588 - t352 * t95 + t354 * t96) * t322;
t28 = t482 * t502 + t390;
t25 = (-t218 + t501) * t482 + t386 + t461;
t16 = t460 * t482 + t364;
t11 = (-t218 + t460) * t482 + t364 + t461;
t1 = [t465 + t297 * t518 - t296 * t519 + t288 * t521 - t287 * t524 + t326 * t269 + t325 * t270 + 0.2e1 * m(3) * (t213 * t245 + t214 * t244) + (t144 * t216 + t145 * t215) * t580 + (t119 * t206 + t120 * t205) * t579 + (t103 * t171 + t104 * t170) * t578 + (t109 * t53 + t110 * t52) * t577 + (-t418 + t424) * t480 + (t419 + t423) * t479 + (t241 + t589) * t323 + (t242 + t596) * t322; m(7) * (t352 * t53 - t354 * t52 + (t109 * t354 + t110 * t352) * qJD(1)) + m(6) * (qJD(1) * t410 - t103 * t354 + t352 * t104) + m(5) * (qJD(1) * t404 - t119 * t354 + t352 * t120) + m(4) * (qJD(1) * t587 - t144 * t354 + t352 * t145) + m(3) * (-t213 * t354 + t352 * t214 + (t244 * t354 + t245 * t352) * qJD(1)); 0; (t346 / 0.2e1 + t347 / 0.2e1) * t414 * qJD(3) + t355 + m(4) * ((-t144 * t352 - t145 * t354) * t305 - t587 * t294) + m(5) * (t119 * t230 + t120 * t231 + t146 * t205 + t147 * t206) + m(6) * (t103 * t211 + t104 * t212 + t113 * t170 + t114 * t171) + m(7) * (t109 * t67 + t110 * t68 + t136 * t52 + t137 * t53) + ((-t216 * t572 + t531 / 0.2e1 + t529 / 0.2e1) * t354 + (t215 * t572 + t532 / 0.2e1 + t530 / 0.2e1) * t352) * qJD(1) + (-qJD(3) * t398 + (-qJD(1) * t262 - t354 * t373) * t333 + (-qJD(1) * t264 - t354 * t374) * t332) * t576 + (-qJD(3) * t399 + (qJD(1) * t263 - t352 * t373) * t333 + (qJD(1) * t265 - t352 * t374) * t332) * t575; m(5) * (t146 * t352 - t147 * t354 + (t230 * t352 + t534) * qJD(1)) + m(6) * (t113 * t352 - t114 * t354 + (t211 * t352 + t212 * t354) * qJD(1)) + m(7) * (t67 * t352 - t354 * t68 + (t136 * t352 + t137 * t354) * qJD(1)); t367 + (t11 * t50 + t136 * t68 + t137 * t67) * t577 + (t113 * t212 + t114 * t211 + t25 * t56) * t578 + (t146 * t231 + t147 * t230 + t49 * t99) * t579 + t352 * ((t352 * t198 + (t125 + t595) * qJD(1)) * t352 + (t126 * qJD(1) + (t262 * t479 + t264 * t480) * t354 + (-t199 + (-t529 - t531) * qJD(3) + (t261 - t399) * qJD(1)) * t352) * t354) + (-t125 * t354 + t126 * t352) * t481 + (-t123 * t354 + t124 * t352) * t482 - t354 * ((t354 * t199 + (t124 + t592) * qJD(1)) * t354 + (t123 * qJD(1) + (-t263 * t479 - t265 * t480 + t483) * t352 + (-t198 + (t530 + t532) * qJD(3) - t398 * qJD(1)) * t354) * t352) - t354 * t13 - t354 * t9 - t354 * t15 + ((t352 * t271 + t272 * t354) * ((qJD(1) * t271 + t361) * t354 + (-t352 * t383 + (-t272 - t475 + t342) * qJD(1)) * t352) + t486 * t305 * t294) * t580; t355 + (-t119 * t352 - t120 * t354 + (t205 * t352 - t206 * t354) * qJD(1)) * t571 - m(5) * t404 * t275 + m(6) * (t103 * t224 + t104 * t225 + t130 * t170 + t131 * t171) + m(7) * (t109 * t78 + t110 * t79 + t148 * t52 + t149 * t53); m(6) * (t130 * t352 - t131 * t354 + (t224 * t352 + t225 * t354) * qJD(1)) + m(7) * (t78 * t352 - t354 * t79 + (t148 * t352 + t149 * t354) * qJD(1)); m(5) * (t173 * t49 - t230 * t508 - t275 * t534 + t80 * t99) + m(7) * (t11 * t51 + t136 * t79 + t137 * t78 + t148 * t68 + t149 * t67 + t16 * t50) + m(6) * (t113 * t225 + t114 * t224 + t130 * t212 + t131 * t211 + t25 * t81 + t36 * t56) + t359 + (-t146 * t354 - t147 * t352 + (-t230 * t354 + t231 * t352) * qJD(1)) * t571; (t148 * t79 + t149 * t78 + t16 * t51) * t577 + (t130 * t225 + t131 * t224 + t36 * t81) * t578 + (t275 * t299 * t486 + t173 * t80) * t579 + t359; m(7) * (t109 * t83 + t110 * t84 + t180 * t52 + t181 * t53) + t356 + (-t103 * t352 - t104 * t354 + (t170 * t352 - t171 * t354) * qJD(1)) * t570 - m(6) * t410 * t243; m(7) * (t83 * t352 - t354 * t84 + (t180 * t352 + t181 * t354) * qJD(1)); m(7) * (t11 * t82 + t136 * t84 + t137 * t83 + t180 * t68 + t181 * t67 + t28 * t50) + m(6) * (t160 * t25 - t211 * t509 - t212 * t533 + t73 * t56) + (-t113 * t354 - t114 * t352 + (-t211 * t354 + t212 * t352) * qJD(1)) * t570 + t365; m(7) * (t148 * t84 + t149 * t83 + t16 * t82 + t180 * t79 + t181 * t78 + t28 * t51) + m(6) * (t160 * t36 - t224 * t509 - t225 * t533 + t73 * t81) + (-t130 * t354 - t131 * t352 + (-t224 * t354 + t225 * t352) * qJD(1)) * t570 + t365; (t243 * t289 * t486 + t160 * t73) * t578 + (t180 * t84 + t181 * t83 + t28 * t82) * t577 + t365; m(7) * (t109 * t48 + t110 * t47 + t111 * t53 + t112 * t52) + t100 + (-t37 + (t352 * t472 + t354 * t471) * t334) * t323 + ((t24 / 0.2e1 + t34 / 0.2e1) * t354 + (t23 / 0.2e1 + t35 / 0.2e1) * t352 + (-t352 * t471 + t354 * t472) * qJD(1)) * t322; m(7) * (t48 * t352 - t354 * t47 + (t111 * t354 + t112 * t352) * qJD(1)); t362 + m(7) * (t101 * t11 + t111 * t67 + t112 * t68 + t136 * t47 + t137 * t48 + t31 * t50); t362 + m(7) * (t101 * t16 + t111 * t78 + t112 * t79 + t148 * t47 + t149 * t48 + t31 * t51); t362 + m(7) * (t101 * t28 + t111 * t83 + t112 * t84 + t180 * t47 + t181 * t48 + t31 * t82); (t101 * t31 + t111 * t48 + t112 * t47) * t577 + (t37 * t323 - t100 + (t352 * t29 + t354 * t30 - t323 * t426) * t334) * t323 + (t354 * t3 + t352 * t4 + t426 * t524 + (-t102 * t334 - t23 * t352 - t24 * t354) * t323 + (t354 * t29 - t352 * t30 + t323 * t427) * qJD(1)) * t322;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
