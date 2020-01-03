% Calculate vector of inverse dynamics joint torques for
% S5RRRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 21:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRP3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP3_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRP3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRP3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:49:07
% EndTime: 2019-12-31 21:49:24
% DurationCPUTime: 12.89s
% Computational Cost: add. (17323->540), mult. (12359->626), div. (0->0), fcn. (9523->8), ass. (0->312)
t294 = qJ(1) + qJ(2);
t286 = qJ(3) + t294;
t275 = sin(t286);
t295 = sin(qJ(4));
t297 = cos(qJ(4));
t240 = Icges(5,5) * t297 - Icges(5,6) * t295;
t276 = cos(t286);
t340 = t240 * t276;
t137 = Icges(5,3) * t275 + t340;
t242 = Icges(6,4) * t297 + Icges(6,6) * t295;
t341 = t242 * t276;
t139 = Icges(6,2) * t275 + t341;
t576 = t137 + t139;
t284 = Icges(6,5) * t295;
t363 = Icges(6,1) * t297 + t284;
t142 = -Icges(6,4) * t276 + t275 * t363;
t444 = t275 * t295;
t228 = Icges(5,4) * t444;
t443 = t275 * t297;
t144 = Icges(5,1) * t443 - Icges(5,5) * t276 - t228;
t575 = -t142 - t144;
t343 = t363 * t276;
t143 = Icges(6,4) * t275 + t343;
t461 = Icges(5,4) * t295;
t248 = Icges(5,1) * t297 - t461;
t344 = t248 * t276;
t145 = Icges(5,5) * t275 + t344;
t568 = t143 + t145;
t243 = Icges(5,2) * t297 + t461;
t361 = Icges(6,3) * t297 - t284;
t574 = t243 + t361;
t460 = Icges(6,5) * t297;
t245 = Icges(6,1) * t295 - t460;
t285 = Icges(5,4) * t297;
t247 = Icges(5,1) * t295 + t285;
t573 = t245 + t247;
t238 = Icges(6,3) * t295 + t460;
t362 = -Icges(5,2) * t295 + t285;
t572 = t238 - t362;
t571 = t248 + t363;
t438 = t276 * t297;
t227 = Icges(6,5) * t438;
t439 = t276 * t295;
t135 = Icges(6,6) * t275 + Icges(6,3) * t439 + t227;
t570 = t135 * t439 + t576 * t275 + t568 * t438;
t138 = -Icges(6,2) * t276 + t242 * t275;
t127 = t275 * t138;
t134 = -Icges(6,6) * t276 + t238 * t275;
t136 = Icges(5,5) * t443 - Icges(5,6) * t444 - Icges(5,3) * t276;
t569 = -t134 * t439 - t275 * t136 + t575 * t438 - t127;
t239 = Icges(5,5) * t295 + Icges(5,6) * t297;
t241 = Icges(6,4) * t295 - Icges(6,6) * t297;
t567 = t239 + t241;
t293 = qJD(1) + qJD(2);
t280 = qJD(3) + t293;
t566 = (-Icges(5,6) + Icges(6,6)) * t280 + t574 * qJD(4);
t565 = (Icges(6,4) + Icges(5,5)) * t280 - t573 * qJD(4);
t561 = -t574 * t295 + t573 * t297;
t140 = Icges(5,4) * t443 - Icges(5,2) * t444 - Icges(5,6) * t276;
t455 = t140 * t295;
t356 = -t144 * t297 + t455;
t441 = t276 * t138;
t359 = t134 * t295 + t142 * t297;
t505 = t275 * t359;
t49 = -t441 + t505;
t511 = -t276 * t136 - t275 * t356 + t49;
t509 = -t140 * t439 - t569;
t342 = t362 * t276;
t141 = Icges(5,6) * t275 + t342;
t508 = -t141 * t439 + t570;
t376 = -t135 * t444 + t139 * t276 - t143 * t443;
t116 = t145 * t443;
t379 = t276 * t137 - t116;
t52 = -t141 * t444 - t379;
t510 = t52 - t376;
t564 = t572 * qJD(4);
t563 = t571 * qJD(4);
t562 = t240 + t242;
t560 = -t573 * t295 - t574 * t297;
t538 = rSges(6,1) + pkin(4);
t445 = t275 * t280;
t559 = t566 * t276 - t572 * t445;
t440 = t276 * t280;
t558 = -t238 * t440 - t275 * t566 + t280 * t342;
t557 = t565 * t276 - t571 * t445;
t556 = (t343 + t344) * t280 + t565 * t275;
t546 = rSges(6,3) + qJ(5);
t555 = (t135 - t141) * t297 - t568 * t295;
t507 = (t134 - t140) * t297 + t575 * t295;
t501 = t546 * t443;
t447 = t241 * t276;
t450 = t239 * t276;
t531 = t275 * t561 - t447 - t450;
t448 = t241 * t275;
t451 = t239 * t275;
t530 = t276 * t561 + t448 + t451;
t554 = (-Icges(6,2) - Icges(5,3)) * t280 + t567 * qJD(4);
t454 = t141 * t295;
t553 = t135 * t295 + t297 * t568 - t454;
t552 = -t356 + t359;
t283 = t295 * qJ(5);
t374 = t297 * rSges(6,1) + t295 * rSges(6,3);
t551 = t297 * pkin(4) + t283 + t374;
t550 = t560 * qJD(4) + t280 * t567 + t564 * t295 + t563 * t297;
t549 = t508 * t275 - t509 * t276;
t548 = t510 * t275 - t511 * t276;
t193 = rSges(4,1) * t275 + rSges(4,2) * t276;
t175 = t280 * t193;
t296 = sin(qJ(1));
t469 = pkin(1) * qJD(1);
t399 = t296 * t469;
t281 = sin(t294);
t436 = t281 * t293;
t401 = pkin(2) * t436;
t123 = -t399 - t175 - t401;
t547 = -t562 * qJD(4) + t561 * t280;
t545 = t530 * t280;
t267 = t276 * rSges(6,2);
t424 = t551 * t275 - t267;
t270 = t276 * pkin(8);
t195 = pkin(3) * t275 - t270;
t214 = pkin(8) * t440;
t544 = t280 * t195 + t214;
t411 = t538 * t295 - t546 * t297;
t543 = t555 * qJD(4) + t576 * t280 + t559 * t295 + t557 * t297;
t542 = -t556 * t297 + t558 * t295 + (-t136 - t138) * t280 - t507 * qJD(4);
t541 = t531 * t280;
t540 = (-t341 - t340 + t552) * t280 + t554 * t275;
t539 = t554 * t276 + t553 * t280 + t562 * t445;
t537 = t548 * qJD(4) + t541;
t536 = t549 * qJD(4) + t545;
t535 = t552 * qJD(4) + t556 * t295 + t558 * t297;
t534 = t553 * qJD(4) + t557 * t295 - t559 * t297;
t533 = -t547 * t275 + t550 * t276;
t532 = t550 * t275 + t547 * t276;
t528 = t275 * rSges(6,2) + pkin(4) * t438;
t416 = rSges(5,2) * t444 + t276 * rSges(5,3);
t147 = rSges(5,1) * t443 - t416;
t130 = t280 * t147;
t405 = qJD(4) * t297;
t388 = t276 * t405;
t437 = t280 * t295;
t398 = t275 * t437;
t345 = rSges(5,3) * t440 + (-t388 + t398) * rSges(5,2);
t406 = qJD(4) * t295;
t389 = t276 * t406;
t527 = -rSges(5,1) * t389 + t130 + t345 + t544;
t526 = t441 + t570;
t404 = qJD(5) * t295;
t221 = t276 * t404;
t497 = rSges(6,2) * t440 + t388 * t546 + t221;
t525 = t424 * t280 + t497 + t544;
t390 = t275 * t405;
t524 = t276 * t437 + t390;
t282 = cos(t294);
t199 = rSges(3,1) * t281 + rSges(3,2) * t282;
t452 = t199 * t293;
t156 = -t399 - t452;
t196 = t276 * pkin(3) + t275 * pkin(8);
t155 = t196 * t280;
t407 = qJD(4) * t280;
t161 = -qJDD(4) * t276 + t275 * t407;
t292 = qJDD(1) + qJDD(2);
t279 = qJDD(3) + t292;
t291 = t293 ^ 2;
t298 = cos(qJ(1));
t299 = qJD(1) ^ 2;
t339 = (-qJDD(1) * t296 - t298 * t299) * pkin(1);
t311 = t339 + (-t281 * t292 - t282 * t291) * pkin(2);
t403 = qJD(5) * t297;
t419 = -qJD(4) * t551 + t403;
t324 = qJDD(5) * t295 + (t403 + t419) * qJD(4);
t387 = t275 * t404;
t396 = -t195 - t424;
t391 = t275 * t406;
t216 = rSges(6,1) * t391;
t220 = pkin(4) * t391;
t375 = -t220 + t387;
t471 = rSges(6,3) * t390 - t216 + t524 * qJ(5) + t375 + (t276 * t374 + t528) * t280;
t6 = t411 * t161 + t396 * t279 + (-t155 - t387 - t471) * t280 + t324 * t276 + t311;
t523 = -g(1) + t6;
t160 = qJDD(4) * t275 + t276 * t407;
t274 = pkin(2) * t282;
t290 = t298 * pkin(1);
t475 = pkin(1) * t296;
t377 = qJDD(1) * t290 - t299 * t475;
t474 = pkin(2) * t281;
t328 = t292 * t274 - t291 * t474 + t377;
t314 = t280 * (-pkin(3) * t445 + t214) + t279 * t196 + t328;
t422 = rSges(6,1) * t438 + rSges(6,3) * t439 + t276 * t283 + t528;
t334 = -t280 * t443 - t389;
t472 = t538 * t334 - t398 * t546 + t497;
t7 = t422 * t279 - t411 * t160 + (t221 + t472) * t280 + t324 * t275 + t314;
t522 = -g(2) + t7;
t470 = rSges(5,1) * t297;
t256 = -rSges(5,2) * t295 + t470;
t211 = t256 * qJD(4);
t252 = rSges(5,1) * t295 + rSges(5,2) * t297;
t408 = qJD(4) * t276;
t423 = -t147 - t195;
t394 = rSges(5,1) * t391 + t524 * rSges(5,2);
t499 = rSges(5,1) * t438 + t275 * rSges(5,3);
t91 = t280 * t499 - t394;
t25 = -t211 * t408 + t161 * t252 + (-t155 - t91) * t280 + t423 * t279 + t311;
t521 = -g(1) + t25;
t154 = rSges(4,1) * t440 - rSges(4,2) * t445;
t520 = -t154 * t280 - t193 * t279 - g(1) + t311;
t149 = -rSges(5,2) * t439 + t499;
t409 = qJD(4) * t275;
t89 = rSges(5,1) * t334 + t345;
t26 = t149 * t279 - t160 * t252 - t211 * t409 + t280 * t89 + t314;
t519 = -g(2) + t26;
t194 = t276 * rSges(4,1) - rSges(4,2) * t275;
t518 = -t175 * t280 + t194 * t279 - g(2) + t328;
t517 = t275 * t540 + t276 * t542;
t516 = -t275 * t539 + t276 * t543;
t435 = t282 * t293;
t186 = rSges(3,1) * t435 - rSges(3,2) * t436;
t515 = -t186 * t293 - t199 * t292 - g(1) + t339;
t200 = t282 * rSges(3,1) - rSges(3,2) * t281;
t514 = t200 * t292 - t293 * t452 - g(2) + t377;
t513 = t275 * t542 - t276 * t540;
t512 = t275 * t543 + t276 * t539;
t503 = t411 * t409;
t502 = t196 + t422;
t500 = t546 * t438;
t329 = -t295 * t546 - t297 * t538 - pkin(3);
t393 = t276 * t538;
t347 = -t408 * t411 + t221;
t327 = t347 - t401;
t312 = t327 - t399;
t46 = t280 * t396 + t312;
t468 = t276 * t46;
t400 = t298 * t469;
t402 = pkin(2) * t435;
t337 = t400 + t402;
t315 = t387 + t337;
t47 = t280 * t502 + t315 - t503;
t300 = (t329 * t468 + (t46 * (-rSges(6,2) - pkin(8)) + t47 * (-pkin(3) - t551)) * t275) * t280 + (-t295 * t393 * t47 - t46 * t501) * qJD(4);
t466 = t280 * t46;
t498 = t466 * t502 + t300;
t495 = -t401 + t527;
t494 = t216 - t503;
t113 = t149 + t196;
t365 = -t113 * t280 + t252 * t409;
t487 = -t401 + t525;
t384 = -pkin(3) - t470;
t392 = t252 * t408;
t336 = -t392 - t401;
t316 = t336 - t399;
t61 = t280 * t423 + t316;
t467 = t276 * t61;
t62 = t337 - t365;
t303 = (t384 * t467 + (t61 * (-rSges(5,3) - pkin(8)) + t62 * t384) * t275) * t280;
t486 = t303 + (-t365 + t394) * t61;
t426 = -Icges(5,2) * t443 + t144 - t228;
t430 = t247 * t275 + t140;
t485 = -t295 * t426 - t297 * t430;
t428 = -t361 * t275 + t142;
t432 = -t245 * t275 + t134;
t484 = -t295 * t428 + t297 * t432;
t483 = t160 / 0.2e1;
t482 = t161 / 0.2e1;
t473 = g(2) * t275;
t48 = -t403 + (t424 * t275 + t422 * t276) * qJD(4);
t465 = t48 * t295;
t453 = t194 * t280;
t449 = t240 * t280;
t446 = t242 * t280;
t431 = -Icges(6,1) * t439 + t135 + t227;
t429 = -t247 * t276 - t141;
t427 = -t361 * t276 + t143;
t425 = -t243 * t276 + t145;
t421 = -t538 * t444 + t501;
t420 = -t538 * t439 + t500;
t415 = -t361 + t363;
t414 = t238 - t245;
t413 = -t243 + t248;
t412 = t247 + t362;
t383 = -t409 / 0.2e1;
t382 = t409 / 0.2e1;
t381 = -t408 / 0.2e1;
t380 = t408 / 0.2e1;
t378 = -t136 + t454;
t159 = t194 + t274;
t257 = rSges(2,1) * t298 - rSges(2,2) * t296;
t253 = rSges(2,1) * t296 + rSges(2,2) * t298;
t372 = t275 * t47 + t468;
t367 = -t275 * t62 - t467;
t354 = t147 * t275 + t149 * t276;
t346 = t372 * t297;
t158 = -t193 - t474;
t333 = -t154 - t402;
t73 = t274 + t502;
t332 = -t295 * t427 + t297 * t431;
t331 = -t295 * t425 + t297 * t429;
t112 = t275 * t384 + t270 + t416;
t111 = t113 + t274;
t326 = (t295 * t414 + t297 * t415) * t280;
t325 = (-t295 * t412 + t297 * t413) * t280;
t110 = t112 - t474;
t98 = t275 * t329 + t267 + t270;
t72 = t98 - t474;
t302 = (((t52 - t116 + (t137 + t455) * t276 + t569) * t276 + (t49 - t505 + t526) * t275) * qJD(4) + t545) * t380 + (t561 * qJD(4) + t563 * t295 - t564 * t297) * t280 + (Icges(4,3) - t560) * t279 + (t530 - t555) * t483 + (t531 - t507) * t482 + (t533 + t534) * t382 + (((t276 * t378 + t508 - t526) * t276 + (t275 * t378 - t127 + t376 + t379 + t509) * t275) * qJD(4) + t537 - t541) * t383 + (t532 + t535 + t536) * t381;
t301 = Icges(3,3) * t292 + t302;
t182 = t252 * t276;
t178 = t252 * t275;
t157 = t200 * t293 + t400;
t124 = t337 + t453;
t71 = t354 * qJD(4);
t5 = -qJDD(5) * t297 - t422 * t161 + t424 * t160 + (t275 * t471 + t276 * t472 + t404) * qJD(4);
t1 = [Icges(2,3) * qJDD(1) + t301 + (t514 * (t200 + t290) + t515 * (-t199 - t475) + (-t186 - t400 + t157) * t156) * m(3) + (g(1) * t253 - g(2) * t257 + (t253 ^ 2 + t257 ^ 2) * qJDD(1)) * m(2) + (t46 * (t216 + t220 - t315) + t300 + t522 * (t290 + t73) + t523 * (t72 - t475) + (-t399 + t46 - t312 + t487) * t47) * m(6) + (t61 * (-t337 + t394) + t303 + (-t399 - t316 + t61 + t495) * t62 + t519 * (t111 + t290) + t521 * (t110 - t475)) * m(5) + (t518 * (t159 + t290) + t520 * (t158 - t475) + (t333 - t400 + t124) * t123) * m(4); t301 + (t522 * t73 + t523 * t72 + (-t327 + t487) * t47 + (t220 + t494) * t46 + t498) * m(6) + ((-t336 + t495) * t62 + t519 * t111 + t521 * t110 + t486) * m(5) + (t518 * t159 + t520 * t158 + (t333 + t402 + t453) * t123) * m(4) + (-t156 * t186 - t157 * t452 + (t156 * t293 + t514) * t200 + (t157 * t293 - t515) * t199) * m(3); t302 + (t522 * t502 + t523 * t98 + (-t347 + t525) * t47 + (t387 - t375 + t494) * t46 + t498) * m(6) + ((t392 + t527) * t62 + t519 * t113 + t521 * t112 + t486) * m(5) + (-t123 * t154 - t124 * t175 + (t123 * t280 + t518) * t194 + (t124 * t280 - t520) * t193) * m(4); t549 * t483 + t548 * t482 + (t533 * t280 + t530 * t279 + t509 * t161 + t508 * t160 + (t516 * t275 + t517 * t276) * qJD(4)) * t275 / 0.2e1 - (t532 * t280 + t531 * t279 + t511 * t161 + t510 * t160 + (t512 * t275 + t513 * t276) * qJD(4)) * t276 / 0.2e1 + (-t275 * t555 + t507 * t276) * t279 / 0.2e1 - (((t412 - t414) * t297 + (t413 + t415) * t295) * t280 + (((-t426 - t428) * t276 + (t425 + t427) * t275) * t297 + ((t430 - t432) * t276 + (t429 + t431) * t275) * t295) * qJD(4)) * t280 / 0.2e1 + ((-t280 * t555 - t535) * t276 + (-t280 * t507 + t534) * t275) * t280 / 0.2e1 + t537 * t445 / 0.2e1 + t536 * t440 / 0.2e1 + ((-t409 * t447 + t446) * t275 + (t326 + (-t484 * t276 + (t448 + t332) * t275) * qJD(4)) * t276 + (-t409 * t450 + t449) * t275 + (t325 + (-t485 * t276 + (t451 + t331) * t275) * qJD(4)) * t276) * t383 + ((t280 * t508 + t517) * t276 + (t280 * t509 + t516) * t275) * t382 + ((t280 * t510 + t513) * t276 + (t280 * t511 + t512) * t275) * t381 + ((-t408 * t448 - t446) * t276 + (t326 + (t332 * t275 + (t447 - t484) * t276) * qJD(4)) * t275 + (-t408 * t451 - t449) * t276 + (t325 + (t331 * t275 + (t450 - t485) * t276) * qJD(4)) * t275) * t380 + (-(t346 + t465) * qJD(5) - (t420 * t47 - t421 * t46) * t280 - ((t420 * t48 - t46 * t551) * t276 + (t421 * t48 - t47 * t551) * t275) * qJD(4) - g(1) * t500 - g(2) * t501 - g(3) * t551 - (-g(1) * t393 - t473 * t538) * t295 + (-t6 * t411 + t46 * t419 + t5 * t422 + t48 * t472 + (-t411 * t47 + t424 * t48) * t280) * t276 + (-t7 * t411 + t47 * t419 + t5 * t424 + t48 * t471 + (t411 * t46 - t422 * t48) * t280) * t275) * m(6) + (g(1) * t182 + g(2) * t178 - g(3) * t256 + (t147 * t160 - t149 * t161 + (t275 * t91 + t276 * t89) * qJD(4)) * t354 + t71 * ((t89 + t130) * t276 + (-t149 * t280 + t91) * t275) + t367 * t211 + ((-t280 * t62 - t25) * t276 + (t280 * t61 - t26) * t275) * t252 - (t178 * t61 - t182 * t62) * t280 - (t71 * (-t178 * t275 - t182 * t276) + t367 * t256) * qJD(4)) * m(5); (-(-t275 * t46 + t276 * t47) * t437 - (t346 + (t275 ^ 2 + t276 ^ 2) * t465) * qJD(4) + (qJD(4) * t372 + g(3) - t5) * t297 + (qJD(4) * t48 + (t7 - t466) * t275 - t473 + (t280 * t47 + t523) * t276) * t295) * m(6);];
tau = t1;
