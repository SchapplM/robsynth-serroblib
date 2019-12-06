% Calculate vector of inverse dynamics joint torques for
% S5RPRRR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:06
% EndTime: 2019-12-05 18:16:24
% DurationCPUTime: 10.36s
% Computational Cost: add. (18888->605), mult. (11804->758), div. (0->0), fcn. (9175->10), ass. (0->356)
t306 = qJD(1) ^ 2;
t298 = qJ(1) + pkin(9);
t289 = qJ(3) + t298;
t283 = cos(t289);
t273 = qJDD(4) * t283;
t296 = qJD(4) + qJD(5);
t282 = sin(t289);
t297 = qJD(1) + qJD(3);
t449 = t282 * t297;
t131 = qJDD(5) * t283 - t296 * t449 + t273;
t410 = qJD(4) * t297;
t191 = -t282 * t410 + t273;
t299 = qJ(4) + qJ(5);
t291 = cos(t299);
t281 = t291 * rSges(6,1);
t290 = sin(t299);
t228 = -rSges(6,2) * t290 + t281;
t197 = t228 * t296;
t204 = t283 * t296;
t227 = rSges(6,1) * t290 + rSges(6,2) * t291;
t295 = qJDD(1) + qJDD(3);
t300 = sin(qJ(4));
t496 = pkin(7) * t282;
t497 = pkin(3) * t283;
t207 = t496 + t497;
t287 = sin(t298);
t288 = cos(t298);
t301 = sin(qJ(1));
t303 = cos(qJ(1));
t436 = t303 * t306;
t312 = (-qJDD(1) * t287 - t288 * t306) * pkin(2) + (-qJDD(1) * t301 - t436) * pkin(1);
t495 = t283 * pkin(7);
t364 = -pkin(3) * t282 + t495;
t310 = -t207 * t297 ^ 2 + t295 * t364 + t312;
t302 = cos(qJ(4));
t437 = t302 * qJD(4) ^ 2;
t409 = qJD(4) * t300;
t395 = t282 * t409;
t250 = pkin(4) * t395;
t304 = -pkin(8) - pkin(7);
t438 = t297 * t304;
t417 = t282 * t438 + t250;
t294 = t302 * pkin(4);
t284 = t294 + pkin(3);
t494 = pkin(3) - t284;
t113 = (t283 * t494 + t496) * t297 + t417;
t443 = t283 * t291;
t404 = rSges(6,1) * t443;
t342 = rSges(6,3) * t282 + t404;
t450 = t282 * t291;
t451 = t282 * t290;
t175 = rSges(6,1) * t451 + rSges(6,2) * t450;
t444 = t283 * t290;
t242 = rSges(6,2) * t444;
t399 = t175 * t296 + t297 * t242;
t95 = -t297 * t342 + t399;
t481 = t113 + t95;
t538 = rSges(6,2) * t451 + rSges(6,3) * t283;
t493 = pkin(7) + t304;
t540 = t493 * t283;
t344 = t540 - t538;
t380 = -t284 - t281;
t367 = -pkin(3) - t380;
t552 = t282 * t367 + t344;
t563 = g(3) + t131 * t227 + t204 * t197 - t481 * t297 + t552 * t295 - (-t191 * t300 - t283 * t437) * pkin(4) - t310;
t238 = t283 * t284;
t139 = t282 * t493 - t238 + t497;
t155 = -t242 + t342;
t431 = t139 - t155;
t360 = rSges(4,1) * t282 + rSges(4,2) * t283;
t177 = t360 * t297;
t498 = pkin(2) * t287;
t500 = pkin(1) * t301;
t366 = t498 + t500;
t340 = t366 * qJD(1);
t145 = t340 + t177;
t236 = Icges(6,4) * t451;
t151 = -Icges(6,1) * t450 + Icges(6,5) * t283 + t236;
t237 = Icges(6,4) * t444;
t152 = Icges(6,1) * t443 + Icges(6,5) * t282 - t237;
t203 = t282 * t296;
t280 = Icges(6,4) * t291;
t354 = -Icges(6,2) * t290 + t280;
t537 = Icges(6,1) * t290 + t280;
t550 = t537 + t354;
t313 = t203 * (-Icges(6,2) * t443 + t152 - t237) + t204 * (Icges(6,2) * t450 + t151 + t236) + t297 * t550;
t334 = t354 * t282;
t149 = Icges(6,6) * t283 - t334;
t150 = Icges(6,4) * t443 - Icges(6,2) * t444 + Icges(6,6) * t282;
t477 = Icges(6,4) * t290;
t223 = Icges(6,2) * t291 + t477;
t226 = Icges(6,1) * t291 - t477;
t516 = t203 * (t283 * t537 + t150) + t204 * (-t282 * t537 + t149) + t297 * (t223 - t226);
t562 = t290 * t313 + t291 * t516;
t254 = pkin(3) * t449;
t405 = pkin(4) * t409;
t421 = -t283 * t438 - t284 * t449;
t112 = t254 + (-pkin(7) * t297 - t405) * t283 + t421;
t190 = qJDD(4) * t282 + t283 * t410;
t442 = t283 * t297;
t130 = qJD(5) * t442 + qJDD(5) * t282 + t190;
t179 = pkin(7) * t442 - t254;
t286 = t306 * t500;
t499 = pkin(1) * t303;
t365 = pkin(2) * t288 + t499;
t318 = -qJDD(1) * t365 + t306 * t498 + t286;
t400 = -t207 + t431;
t406 = rSges(6,1) * t450;
t554 = t204 * t227;
t398 = -t297 * t406 - t554;
t94 = t297 * t538 + t398;
t18 = t130 * t227 + t197 * t203 + (t190 * t300 + t282 * t437) * pkin(4) + (-t112 - t179 - t94) * t297 + t400 * t295 + t318;
t561 = -g(2) + t18;
t408 = qJD(4) * t302;
t392 = t283 * t408;
t393 = t283 * t409;
t447 = t282 * t302;
t407 = rSges(5,1) * t447;
t397 = -rSges(5,1) * t393 - rSges(5,2) * t392 - t297 * t407;
t448 = t282 * t300;
t415 = rSges(5,2) * t448 + rSges(5,3) * t283;
t110 = t297 * t415 + t397;
t488 = rSges(5,2) * t300;
t490 = rSges(5,1) * t302;
t270 = -t488 + t490;
t243 = t270 * qJD(4);
t268 = rSges(5,1) * t300 + rSges(5,2) * t302;
t412 = qJD(4) * t282;
t441 = t283 * t300;
t260 = rSges(5,2) * t441;
t440 = t283 * t302;
t343 = -rSges(5,1) * t440 - rSges(5,3) * t282;
t165 = -t260 - t343;
t422 = -t165 - t207;
t48 = t243 * t412 + t190 * t268 + (-t110 - t179) * t297 + t422 * t295 + t318;
t560 = -g(2) + t48;
t491 = rSges(4,1) * t283;
t206 = -rSges(4,2) * t282 + t491;
t559 = t177 * t297 - t206 * t295 - g(2) + t318;
t396 = rSges(5,1) * t395 + (t282 * t408 + t297 * t441) * rSges(5,2);
t111 = t297 * t343 + t396;
t347 = t407 - t415;
t411 = qJD(4) * t283;
t47 = t111 * t297 - t191 * t268 - t243 * t411 - t295 * t347 + t310;
t557 = -g(3) + t47;
t253 = rSges(4,2) * t449;
t178 = -rSges(4,1) * t442 + t253;
t556 = t178 * t297 - t295 * t360 - g(3) + t312;
t455 = t227 * t283;
t555 = t175 * t203 + t204 * t455 + t283 * t94;
t292 = Icges(5,4) * t302;
t355 = -Icges(5,2) * t300 + t292;
t536 = Icges(5,1) * t300 + t292;
t413 = t536 + t355;
t478 = Icges(5,4) * t300;
t264 = Icges(5,2) * t302 + t478;
t267 = Icges(5,1) * t302 - t478;
t414 = t264 - t267;
t553 = (t300 * t413 + t302 * t414) * t297;
t551 = t175 * t297 - t204 * t228 - t227 * t449;
t549 = -t282 * t494 + t540;
t532 = t365 * qJD(1);
t369 = t203 * t227 + t250;
t61 = t297 * t400 + t369 - t532;
t548 = (t197 * t282 - t203 * t228 + t227 * t442 - t297 * t455) * t61;
t148 = Icges(6,5) * t443 - Icges(6,6) * t444 + Icges(6,3) * t282;
t64 = t148 * t283 + t150 * t451 - t152 * t450;
t350 = t223 * t290 - t291 * t537;
t221 = Icges(6,5) * t290 + Icges(6,6) * t291;
t460 = t221 * t283;
t97 = t282 * t350 + t460;
t547 = t203 * t64 + t297 * t97;
t161 = Icges(5,4) * t440 - Icges(5,2) * t441 + Icges(5,6) * t282;
t258 = Icges(5,4) * t441;
t163 = Icges(5,1) * t440 + Icges(5,5) * t282 - t258;
t351 = t161 * t300 - t163 * t302;
t542 = t351 * t283;
t353 = t150 * t290 - t152 * t291;
t541 = t353 * t283;
t154 = t297 * t165;
t202 = t268 * t412;
t539 = -t154 + t202;
t279 = t287 * rSges(3,2);
t492 = rSges(3,1) * t288;
t363 = -t492 - t499;
t201 = t279 + t363;
t361 = rSges(3,1) * t287 + rSges(3,2) * t288;
t331 = t361 + t500;
t457 = t223 * t296;
t534 = -Icges(6,6) * t297 + t457;
t153 = t297 * t347;
t198 = t297 * t364;
t345 = t268 * t411 + t153 - t198;
t78 = t340 + t345;
t79 = t297 * t422 + t202 - t532;
t533 = t282 * t79 + t283 * t78;
t531 = t297 * t431 + t369;
t335 = t355 * t297;
t525 = -Icges(5,6) * t297 + qJD(4) * t264;
t107 = t282 * t525 - t283 * t335;
t336 = t267 * t297;
t523 = -Icges(5,5) * t297 + qJD(4) * t536;
t109 = t282 * t523 - t283 * t336;
t263 = Icges(5,5) * t302 - Icges(5,6) * t300;
t333 = t263 * t282;
t158 = Icges(5,3) * t283 - t333;
t160 = Icges(5,6) * t283 - t282 * t355;
t257 = Icges(5,4) * t448;
t162 = -Icges(5,1) * t447 + Icges(5,5) * t283 + t257;
t99 = t160 * t302 + t162 * t300;
t530 = qJD(4) * t99 + t107 * t300 - t109 * t302 - t158 * t297;
t529 = -Icges(6,3) * t297 + t221 * t296;
t528 = -Icges(6,5) * t297 + t296 * t537;
t100 = t161 * t302 + t163 * t300;
t106 = -t282 * t335 - t283 * t525;
t108 = -t282 * t336 - t283 * t523;
t159 = Icges(5,5) * t440 - Icges(5,6) * t441 + Icges(5,3) * t282;
t527 = qJD(4) * t100 + t106 * t300 - t108 * t302 - t159 * t297;
t262 = Icges(5,5) * t300 + Icges(5,6) * t302;
t526 = -Icges(5,3) * t297 + qJD(4) * t262;
t232 = t355 * qJD(4);
t233 = t267 * qJD(4);
t349 = t264 * t302 + t300 * t536;
t524 = qJD(4) * t349 + t232 * t300 - t233 * t302 - t262 * t297;
t423 = -Icges(5,2) * t440 + t163 - t258;
t425 = t283 * t536 + t161;
t521 = t300 * t423 + t302 * t425;
t424 = Icges(5,2) * t447 + t162 + t257;
t426 = -t282 * t536 + t160;
t520 = -t300 * t424 - t302 * t426;
t370 = -t226 * t296 + t457;
t371 = t550 * t296;
t519 = -t221 * t297 + t290 * t371 + t291 * t370;
t376 = t152 * t296 - t283 * t534 - t297 * t334;
t337 = t297 * t226;
t378 = t150 * t296 + t282 * t337 + t283 * t528;
t518 = -t148 * t297 + t290 * t376 + t291 * t378;
t222 = Icges(6,5) * t291 - Icges(6,6) * t290;
t147 = Icges(6,3) * t283 - t222 * t282;
t377 = t151 * t296 + t282 * t534 - t354 * t442;
t379 = t149 * t296 - t282 * t528 + t283 * t337;
t517 = -t147 * t297 + t290 * t377 + t291 * t379;
t515 = m(3) + m(4);
t514 = t130 / 0.2e1;
t513 = t131 / 0.2e1;
t512 = t190 / 0.2e1;
t511 = t191 / 0.2e1;
t510 = -t203 / 0.2e1;
t509 = t203 / 0.2e1;
t508 = -t204 / 0.2e1;
t507 = t204 / 0.2e1;
t506 = t282 / 0.2e1;
t505 = t283 / 0.2e1;
t504 = t295 / 0.2e1;
t503 = -t297 / 0.2e1;
t502 = t297 / 0.2e1;
t501 = -rSges(5,3) - pkin(7);
t169 = t221 * t282;
t98 = -t283 * t350 + t169;
t482 = t98 * t297;
t182 = t262 * t282;
t348 = t264 * t300 - t302 * t536;
t116 = -t283 * t348 + t182;
t468 = t116 * t297;
t465 = t149 * t290;
t464 = t151 * t291;
t463 = t160 * t300;
t462 = t162 * t302;
t458 = t222 * t297;
t456 = t227 * t282;
t454 = t262 * t283;
t453 = t263 * t297;
t188 = t268 * t282;
t452 = t268 * t283;
t446 = t283 * t112;
t445 = t283 * t139;
t439 = t297 * t206;
t435 = t147 * t283 + t149 * t451;
t434 = -t147 * t282 - t151 * t443;
t433 = t158 * t283 + t160 * t448;
t432 = t158 * t282 + t162 * t440;
t261 = pkin(4) * t448;
t391 = -t449 / 0.2e1;
t390 = t442 / 0.2e1;
t389 = -pkin(3) - t490;
t386 = -t412 / 0.2e1;
t385 = t412 / 0.2e1;
t384 = -t411 / 0.2e1;
t383 = t411 / 0.2e1;
t373 = -t159 - t462;
t368 = pkin(4) * t393;
t271 = rSges(2,1) * t303 - rSges(2,2) * t301;
t362 = rSges(2,1) * t301 + rSges(2,2) * t303;
t69 = -t162 * t447 + t433;
t70 = t159 * t283 + t161 * t448 - t163 * t447;
t359 = t282 * t70 + t283 * t69;
t71 = -t160 * t441 + t432;
t72 = t159 * t282 - t542;
t358 = t282 * t72 + t283 * t71;
t82 = t150 * t291 + t152 * t290;
t352 = -t462 + t463;
t346 = t406 - t538;
t65 = -t149 * t444 - t434;
t327 = t169 * t204 - t203 * t460 + t458;
t199 = t297 * t207;
t325 = t199 + t532;
t324 = -t282 * t458 - t283 * t529 + t297 * t353;
t323 = -t283 * t458 + t529 * t282 + (-t464 + t465) * t297;
t322 = -t283 * t526 + (-t333 + t351) * t297;
t321 = -t263 * t442 + t282 * t526 + t297 * t352;
t320 = t222 * t296 + t297 * t350;
t319 = qJD(4) * t263 + t297 * t348;
t120 = -t404 - t238 + t242 + (-rSges(6,3) + t304) * t282;
t124 = t282 * t389 + t415 + t495;
t125 = t282 * t501 + t283 * t389 + t260;
t317 = t554 - t198 + t368 + (t549 + t346) * t297;
t316 = t283 * t165 + t282 * t347;
t119 = t282 * t380 - t283 * t304 + t538;
t13 = t282 * t323 - t283 * t517;
t14 = t282 * t324 - t283 * t518;
t15 = t282 * t517 + t283 * t323;
t16 = t282 * t518 + t283 * t324;
t63 = -t151 * t450 + t435;
t30 = t204 * t63 + t547;
t66 = t148 * t282 - t541;
t31 = t203 * t66 + t204 * t65 + t482;
t40 = -t290 * t379 + t291 * t377;
t41 = -t290 * t378 + t291 * t376;
t45 = t282 * t320 - t283 * t519;
t46 = t282 * t519 + t283 * t320;
t81 = t149 * t291 + t151 * t290;
t315 = (t13 * t204 + t130 * t66 + t131 * t65 + t14 * t203 + t295 * t98 + t297 * t45) * t506 + (t327 * t282 - t283 * t562) * t510 + (t282 * t562 + t327 * t283) * t508 + (t130 * t64 + t131 * t63 + t15 * t204 + t16 * t203 + t295 * t97 + t297 * t46) * t505 + (-t290 * t516 + t291 * t313) * t503 + t30 * t391 + t31 * t390 + ((t297 * t66 + t13) * t283 + (-t297 * t65 + t14) * t282) * t509 + (t282 * t66 + t283 * t65) * t514 + (t282 * t64 + t283 * t63) * t513 + ((t297 * t64 + t15) * t283 + (-t297 * t63 + t16) * t282) * t507 + (t282 * t82 + t283 * t81) * t504 + ((t297 * t82 + t40) * t283 + (-t297 * t81 + t41) * t282) * t502;
t311 = t541 + (-t148 - t464) * t282 + t435;
t115 = t282 * t348 + t454;
t114 = t115 * t297;
t36 = qJD(4) * t359 + t114;
t37 = qJD(4) * t358 + t468;
t51 = -qJD(4) * t352 + t107 * t302 + t109 * t300;
t52 = -qJD(4) * t351 + t106 * t302 + t108 * t300;
t57 = t282 * t319 - t283 * t524;
t58 = t282 * t524 + t283 * t319;
t309 = ((t66 + t311) * t204 + t547) * t510 + (t114 + ((t433 + t72 + t542) * t283 + (-t71 + (t373 - t463) * t283 + t70 + t432) * t282) * qJD(4)) * t386 + (t82 + t98) * t514 + (t81 + t97) * t513 + (t116 + t100) * t512 + (t115 + t99) * t511 + (-t482 + (t64 + (-t148 + t465) * t283 - t353 * t282 + t434) * t204 + (-t63 + t311) * t203 + t31) * t508 + (t40 + t46) * t507 + (-t468 + ((t70 + (-t159 + t463) * t283 - t432) * t283 + (t282 * t373 + t433 - t69) * t282) * qJD(4) + t37) * t384 + (t51 + t58) * t383 + (-qJD(4) * t348 + t232 * t302 + t233 * t300 - t290 * t370 + t291 * t371) * t297 + (t41 + t45 + t30) * t509 + (t52 + t57 + t36) * t385 + (t223 * t291 + t290 * t537 + Icges(4,3) + t349) * t295;
t308 = t79 * (t254 - t397) + ((-t488 * t79 - t501 * t78) * t282 + (-t389 * t78 + t501 * t79) * t283) * t297 - t78 * t396;
t60 = t340 + t317;
t307 = t61 * (t368 - t398 - t421) + (-t61 * t538 - t60 * (-t342 - t238)) * t297 - t60 * (t399 + t417);
t146 = -t532 - t439;
t135 = t283 * t155;
t85 = qJD(4) * t316 + qJD(2);
t53 = qJD(2) + t204 * t155 + t203 * t346 + (t282 * t549 - t445) * qJD(4);
t44 = qJDD(2) + t191 * t165 + t190 * t347 + (t283 * t110 - t282 * t111) * qJD(4);
t22 = t282 * t527 + t283 * t322;
t21 = t282 * t530 + t283 * t321;
t20 = t282 * t322 - t283 * t527;
t19 = t282 * t321 - t283 * t530;
t10 = qJDD(2) - t191 * t139 + t190 * t549 + t131 * t155 + t204 * t94 + t130 * t346 - t203 * t95 + (-t282 * t113 + t446) * qJD(4);
t1 = [t309 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + ((qJDD(1) * t201 + t306 * t361 - g(2) + t286) * t201 + (t436 * pkin(1) + t331 * qJDD(1) + g(3) + (-0.2e1 * t279 + t492 - t363 + t201) * t306) * t331) * m(3) + ((qJDD(1) * t362 + g(3)) * t362 + (qJDD(1) * t271 + g(2)) * t271) * m(2) + (-(t325 + t61 - t531) * t60 + (t365 * t60 + t366 * t61) * qJD(1) + t307 + t561 * (t120 - t365) - t563 * (t119 - t366)) * m(6) + (-(t79 + t325 - t539) * t78 + (t365 * t78 + t366 * t79) * qJD(1) + t308 + t560 * (t125 - t365) + t557 * (t124 - t366)) * m(5) + (t559 * (-t206 - t365) + t556 * (-t360 - t366) + (t297 * t491 - t253 - t439) * t145) * m(4); t515 * qJDD(2) + m(5) * t44 + m(6) * t10 + (-m(5) - m(6) - t515) * g(1); t309 + (-t61 * t317 + t60 * (-t199 + t531) + t307 + t561 * t120 - t563 * t119) * m(6) + (-t79 * t345 + t78 * (-t199 + t539) + t308 + t560 * t125 + t557 * t124) * m(5) + (-t145 * t178 + t146 * t177 + (-t145 * t297 - t559) * t206 - (t146 * t297 + t556) * t360) * m(4); ((-t412 * t454 + t453) * t282 + (-t553 + (t520 * t283 + (t182 - t521) * t282) * qJD(4)) * t283) * t386 + ((t297 * t70 + t21) * t283 + (-t297 * t69 + t22) * t282) * t383 + ((t297 * t72 + t19) * t283 + (-t297 * t71 + t20) * t282) * t385 + ((t100 * t297 + t51) * t283 + (-t297 * t99 + t52) * t282) * t502 + ((-t300 * t414 + t302 * t413) * t297 + ((t282 * t423 + t283 * t424) * t302 + (-t282 * t425 - t283 * t426) * t300) * qJD(4)) * t503 + ((t182 * t411 + t453) * t283 + (t553 + (t521 * t282 + (-t454 - t520) * t283) * qJD(4)) * t282) * t384 + (t116 * t295 + t190 * t72 + t191 * t71 + t297 * t57 + (t19 * t283 + t20 * t282) * qJD(4)) * t506 + (t115 * t295 + t190 * t70 + t191 * t69 + t297 * t58 + (t21 * t283 + t22 * t282) * qJD(4)) * t505 + t315 + t358 * t512 + t359 * t511 + (t100 * t282 + t283 * t99) * t504 + t36 * t391 + t37 * t390 + (-g(1) * (t228 + t294) - g(2) * (t261 + t175) + t10 * (t282 * t552 + t135 - t445) + t18 * (t261 + t456) - t563 * t283 * (-pkin(4) * t300 - t227) + (-(-t282 ^ 2 - t283 ^ 2) * t405 + t446 - t481 * t282 + (t344 * t283 + (t283 * t367 + t431) * t282) * t297 + t555) * t53 + t548 + (-pkin(4) * t392 - (-pkin(4) * t408 - t197) * t283 + t551) * t60) * m(6) + (t44 * t316 + t85 * ((-t111 - t154) * t282 + (t110 + t153) * t283) + t48 * t188 - t47 * t452 + (t442 * t79 - t449 * t78) * t268 + t533 * t243 - (-t188 * t78 + t452 * t79) * t297 - (t85 * (-t188 * t282 - t283 * t452) + t533 * t270) * qJD(4) - g(1) * t270 - g(2) * t188 + g(3) * t452) * m(5); t315 + (t10 * (t282 * t346 + t135) + t18 * t456 - g(1) * t228 - g(2) * t175 + (-t282 * t95 + (-t155 * t282 + t283 * t346) * t297 + t555) * t53 + t548 + (t197 * t283 + t551) * t60 + t563 * t455) * m(6);];
tau = t1;
