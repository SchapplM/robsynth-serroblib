% Calculate vector of inverse dynamics joint torques for
% S5RRRRR13
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR13_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR13_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_invdynJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR13_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR13_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:23
% EndTime: 2024-09-27 17:30:44
% DurationCPUTime: 13.53s
% Computational Cost: add. (52562->775), mult. (31407->1006), div. (0->0), fcn. (27457->16), ass. (0->371)
t376 = qJ(1) + qJ(2);
t368 = qJ(3) + t376;
t356 = cos(t368);
t355 = sin(t368);
t374 = qJD(1) + qJD(2);
t363 = qJD(3) + t374;
t377 = sin(pkin(5));
t514 = t363 * t377;
t484 = t355 * t514;
t383 = pkin(9) + pkin(10);
t378 = cos(pkin(5));
t379 = sin(qJ(4));
t508 = t378 * t379;
t318 = pkin(4) * t508 - t377 * t383;
t537 = pkin(4) * qJD(4);
t486 = t379 * t537;
t520 = t355 * t363;
t498 = -t318 * t520 - t355 * t486;
t381 = cos(qJ(4));
t541 = pkin(4) * t381;
t357 = pkin(3) + t541;
t507 = t378 * t381;
t461 = t537 * t507;
t556 = t363 * (pkin(3) - t357) - t461;
t127 = -pkin(9) * t484 - t356 * t556 + t498;
t375 = qJ(4) + qJ(5);
t364 = sin(t375);
t472 = pkin(5) + t375;
t439 = cos(t472) / 0.2e1;
t473 = pkin(5) - t375;
t451 = cos(t473);
t403 = t451 / 0.2e1 + t439;
t240 = t355 * t364 - t356 * t403;
t450 = sin(t473);
t545 = sin(t472) / 0.2e1;
t316 = t545 - t450 / 0.2e1;
t366 = cos(t375);
t241 = t316 * t356 + t355 * t366;
t516 = t356 * t377;
t163 = t241 * rSges(6,1) - t240 * rSges(6,2) - rSges(6,3) * t516;
t494 = qJD(4) * t377;
t479 = t363 * t494;
t308 = t355 * t479;
t493 = qJD(5) * t363;
t203 = t308 + (t355 * t493 + (-qJDD(4) - qJDD(5)) * t356) * t377;
t294 = pkin(3) * t355 - pkin(9) * t516;
t495 = -t356 * t318 - t355 * t357;
t205 = t294 + t495;
t373 = qJD(4) + qJD(5);
t329 = t373 * t439;
t436 = t373 * t451;
t269 = t329 - t436 / 0.2e1;
t438 = t450 / 0.2e1;
t330 = t373 * t438;
t572 = t373 * t545;
t271 = t330 + t572;
t210 = rSges(6,1) * t271 + rSges(6,2) * t269;
t315 = t545 + t438;
t317 = t439 - t451 / 0.2e1;
t222 = -rSges(6,1) * t317 + rSges(6,2) * t315 + rSges(6,3) * t378;
t492 = qJDD(4) * t377;
t250 = -t356 * t492 + t308;
t463 = t373 * t377;
t274 = t356 * t463;
t510 = t377 * t379;
t290 = pkin(4) * t510 + (-pkin(9) + t383) * t378;
t371 = qJDD(1) + qJDD(2);
t360 = qJDD(3) + t371;
t333 = qJDD(4) * t378 + t360;
t302 = qJDD(5) * t378 + t333;
t340 = qJD(4) * t378 + t363;
t314 = qJD(5) * t378 + t340;
t519 = t355 * t377;
t295 = t356 * pkin(3) + pkin(9) * t519;
t365 = sin(t376);
t367 = cos(t376);
t370 = t374 ^ 2;
t380 = sin(qJ(1));
t382 = cos(qJ(1));
t385 = qJD(1) ^ 2;
t425 = (-qJDD(1) * t380 - t382 * t385) * pkin(1);
t393 = t425 + (-t365 * t371 - t367 * t370) * pkin(2);
t388 = -t295 * t363 ^ 2 - t360 * t294 + t393;
t464 = t377 ^ 2 * qJD(4) ^ 2 * t541;
t390 = -t363 * t403 - t373 * t366;
t466 = t363 * t364 - t330 + t572;
t141 = t355 * t390 - t356 * t466;
t440 = -t316 * t363 - t364 * t373;
t467 = t363 * t366 + t329 + t436 / 0.2e1;
t142 = t355 * t440 + t356 * t467;
t94 = rSges(6,1) * t142 + rSges(6,2) * t141 + rSges(6,3) * t484;
t32 = -t340 * t127 - t302 * t163 + t203 * t222 + t333 * t205 - t274 * t210 + t250 * t290 - t314 * t94 - t356 * t464 + t388;
t603 = -g(1) + t32;
t278 = t355 * t381 + t356 * t508;
t279 = -t355 * t507 - t356 * t379;
t198 = -qJD(4) * t278 + t279 * t363;
t430 = t355 * t508 - t356 * t381;
t431 = -t355 * t379 + t356 * t507;
t199 = qJD(4) * t431 - t363 * t430;
t447 = rSges(5,1) * t199 + rSges(5,2) * t198;
t112 = rSges(5,3) * t484 + t447;
t190 = t278 * rSges(5,1) + rSges(5,2) * t431 - rSges(5,3) * t516;
t289 = rSges(5,3) * t378 + (rSges(5,1) * t379 + rSges(5,2) * t381) * t377;
t310 = (rSges(5,1) * t381 - rSges(5,2) * t379) * t377;
t296 = qJD(4) * t310;
t480 = t356 * t494;
t53 = -t112 * t340 - t190 * t333 + t250 * t289 - t296 * t480 + t388;
t602 = -g(1) + t53;
t306 = rSges(4,1) * t355 + rSges(4,2) * t356;
t281 = t363 * t306;
t538 = pkin(1) * qJD(1);
t488 = t380 * t538;
t513 = t365 * t374;
t491 = pkin(2) * t513;
t223 = -t488 - t281 - t491;
t542 = pkin(2) * t365;
t593 = -t190 - t294;
t599 = t593 - t542;
t594 = -t163 + t495;
t598 = t594 - t542;
t182 = -Icges(5,5) * t278 - Icges(5,6) * t431 + Icges(5,3) * t516;
t261 = Icges(5,4) * t278;
t184 = Icges(5,2) * t431 - Icges(5,6) * t516 + t261;
t260 = Icges(5,4) * t431;
t188 = -Icges(5,1) * t278 + Icges(5,5) * t516 - t260;
t85 = -t378 * t182 + (t184 * t381 - t188 * t379) * t377;
t584 = Icges(6,4) * t241;
t158 = Icges(6,2) * t240 + Icges(6,6) * t516 - t584;
t243 = -t355 * t403 - t356 * t364;
t245 = t316 * t355 - t356 * t366;
t583 = Icges(6,4) * t245;
t159 = Icges(6,2) * t243 + Icges(6,6) * t519 - t583;
t526 = Icges(6,4) * t317;
t220 = Icges(6,2) * t315 + Icges(6,6) * t378 - t526;
t273 = t355 * t463;
t585 = (-Icges(6,1) * t243 + t159 - t583) * t273 + t274 * (-Icges(6,1) * t240 + t158 - t584) + (-Icges(6,1) * t315 + t220 - t526) * t314;
t443 = t184 * t431 - t188 * t278;
t75 = t182 * t516 + t443;
t597 = t356 * t75;
t155 = -Icges(6,5) * t241 + Icges(6,6) * t240 + Icges(6,3) * t516;
t225 = Icges(6,4) * t240;
t161 = -Icges(6,1) * t241 + Icges(6,5) * t516 + t225;
t591 = t158 * t243 - t161 * t245;
t61 = -t155 * t519 - t591;
t156 = -Icges(6,5) * t245 + Icges(6,6) * t243 + Icges(6,3) * t519;
t226 = Icges(6,4) * t243;
t162 = -Icges(6,1) * t245 + Icges(6,5) * t519 + t226;
t62 = t156 * t519 + t243 * t159 - t162 * t245;
t219 = -Icges(6,5) * t317 + Icges(6,6) * t315 + Icges(6,3) * t378;
t299 = Icges(6,4) * t315;
t221 = -Icges(6,1) * t317 + Icges(6,5) * t378 + t299;
t96 = t219 * t519 + t220 * t243 - t221 * t245;
t27 = t273 * t62 - t274 * t61 + t96 * t314;
t59 = t155 * t516 + t158 * t240 - t161 * t241;
t71 = -t155 * t378 - t158 * t315 + t161 * t317;
t257 = t363 * t294;
t590 = -t163 * t314 + t340 * t205 - t257;
t589 = -t190 * t340 - t257;
t588 = -t279 * t184 - t188 * t430;
t586 = (-Icges(6,5) * t240 - Icges(6,6) * t241) * t274 - (Icges(6,5) * t243 + Icges(6,6) * t245) * t273 - (Icges(6,5) * t315 + Icges(6,6) * t317) * t314;
t283 = Icges(5,3) * t378 + (Icges(5,5) * t379 + Icges(5,6) * t381) * t377;
t527 = Icges(5,4) * t379;
t284 = Icges(5,6) * t378 + (Icges(5,2) * t381 + t527) * t377;
t509 = t377 * t381;
t348 = Icges(5,4) * t509;
t285 = Icges(5,1) * t510 + Icges(5,5) * t378 + t348;
t411 = -t278 * t285 + t283 * t516 - t284 * t431;
t579 = t411 * t340;
t559 = -t357 * t363 - t461;
t139 = t355 * t466 + t356 * t390;
t140 = -t355 * t467 + t356 * t440;
t483 = t356 * t514;
t93 = t140 * rSges(6,1) + t139 * rSges(6,2) + rSges(6,3) * t483;
t577 = t355 * t559 - t590 + t93;
t196 = qJD(4) * t430 - t363 * t431;
t197 = qJD(4) * t279 - t278 * t363;
t111 = t197 * rSges(5,1) + t196 * rSges(5,2) + rSges(5,3) * t483;
t313 = pkin(9) * t483;
t576 = t111 + t313 - t589;
t319 = rSges(3,1) * t365 + rSges(3,2) * t367;
t521 = t319 * t374;
t258 = -t488 - t521;
t412 = t219 * t516 + t220 * t240 - t221 * t241;
t575 = t274 * t59 + t412 * t314;
t435 = -t318 * t363 - t486;
t126 = t355 * t556 + t435 * t356 - t313;
t165 = -rSges(6,1) * t245 + t243 * rSges(6,2) + rSges(6,3) * t519;
t249 = t355 * t492 + t356 * t479;
t202 = (qJDD(5) * t355 + t356 * t493) * t377 + t249;
t465 = -t318 * t355 + t356 * t357;
t206 = t465 - t295;
t354 = pkin(2) * t367;
t369 = t382 * pkin(1);
t543 = pkin(1) * t380;
t452 = qJDD(1) * t369 - t385 * t543;
t414 = t371 * t354 - t370 * t542 + t452;
t404 = t363 * (-pkin(3) * t520 + t313) + t360 * t295 + t414;
t33 = t126 * t340 + t165 * t302 - t202 * t222 + t206 * t333 - t210 * t273 - t249 * t290 + t314 * t93 - t355 * t464 + t404;
t571 = -g(2) + t33;
t192 = -rSges(5,1) * t430 + t279 * rSges(5,2) + rSges(5,3) * t519;
t481 = t355 * t494;
t54 = t111 * t340 + t192 * t333 - t249 * t289 - t296 * t481 + t404;
t570 = -g(2) + t54;
t517 = t356 * t363;
t254 = rSges(4,1) * t517 - rSges(4,2) * t520;
t569 = -t254 * t363 - t306 * t360 - g(1) + t393;
t307 = t356 * rSges(4,1) - rSges(4,2) * t355;
t568 = -t281 * t363 + t307 * t360 - g(2) + t414;
t512 = t367 * t374;
t288 = rSges(3,1) * t512 - rSges(3,2) * t513;
t567 = -t288 * t374 - t319 * t371 - g(1) + t425;
t320 = t367 * rSges(3,1) - rSges(3,2) * t365;
t566 = t320 * t371 - t374 * t521 - g(2) + t452;
t423 = -t222 * t274 - t290 * t480;
t405 = t423 - t491;
t394 = t405 - t488;
t68 = t394 + t590;
t522 = t295 * t363;
t408 = -t165 * t314 - t206 * t340 + t273 * t222 + t290 * t481 - t522;
t490 = pkin(2) * t512;
t396 = t408 - t490;
t487 = t382 * t538;
t69 = -t396 + t487;
t565 = (t69 * t435 + t559 * t68) * t356;
t459 = t289 * t480;
t418 = -t459 - t491;
t400 = t418 - t488;
t101 = t400 + t589;
t429 = -t192 * t340 + t289 * t481 - t522;
t409 = t429 - t490;
t102 = -t409 + t487;
t525 = t101 * t356;
t564 = (-pkin(3) * t525 + (-t102 * pkin(3) + t101 * (-rSges(5,3) - pkin(9)) * t377) * t355) * t363;
t562 = t240 * t159 - t241 * t162;
t528 = Icges(5,4) * t430;
t186 = Icges(5,2) * t279 + Icges(5,6) * t519 - t528;
t262 = Icges(5,4) * t279;
t189 = -Icges(5,1) * t430 + Icges(5,5) * t519 + t262;
t505 = -t186 * t431 - t278 * t189;
t557 = -t491 + t576;
t555 = -t491 + t577;
t415 = t355 * (Icges(5,2) * t430 + t189 + t262) - t356 * (-Icges(5,2) * t278 - t188 + t260);
t397 = t273 * (Icges(6,2) * t245 + t162 + t226) - t274 * (-Icges(6,2) * t241 - t161 - t225) + t314 * (Icges(6,2) * t317 + t221 + t299);
t554 = t202 / 0.2e1;
t553 = t203 / 0.2e1;
t552 = t249 / 0.2e1;
t551 = t250 / 0.2e1;
t550 = -t273 / 0.2e1;
t549 = t273 / 0.2e1;
t548 = -t274 / 0.2e1;
t547 = t274 / 0.2e1;
t544 = t378 / 0.2e1;
t207 = Icges(6,5) * t271 + Icges(6,6) * t269;
t208 = Icges(6,4) * t271 + Icges(6,2) * t269;
t209 = Icges(6,1) * t271 + Icges(6,4) * t269;
t65 = t207 * t378 + t208 * t315 - t209 * t317 + t220 * t269 + t221 * t271;
t98 = t219 * t378 + t220 * t315 - t221 * t317;
t539 = t98 * t302 + t65 * t314;
t88 = Icges(6,5) * t142 + Icges(6,6) * t141 + Icges(6,3) * t484;
t90 = Icges(6,4) * t142 + Icges(6,2) * t141 + Icges(6,6) * t484;
t92 = Icges(6,1) * t142 + Icges(6,4) * t141 + Icges(6,5) * t484;
t30 = -t158 * t269 - t161 * t271 + t315 * t90 - t317 * t92 + t378 * t88;
t534 = t30 * t274;
t87 = Icges(6,5) * t140 + Icges(6,6) * t139 + Icges(6,3) * t483;
t89 = Icges(6,4) * t140 + Icges(6,2) * t139 + Icges(6,6) * t483;
t91 = Icges(6,1) * t140 + Icges(6,4) * t139 + Icges(6,5) * t483;
t31 = t159 * t269 + t162 * t271 + t315 * t89 - t317 * t91 + t378 * t87;
t533 = t31 * t273;
t532 = t71 * t203;
t72 = t156 * t378 + t159 * t315 - t162 * t317;
t531 = t72 * t202;
t530 = t85 * t250;
t183 = -Icges(5,5) * t430 + Icges(5,6) * t279 + Icges(5,3) * t519;
t86 = t183 * t378 + (t186 * t381 + t189 * t379) * t377;
t529 = t86 * t249;
t524 = t182 * t355;
t523 = t183 * t356;
t303 = (Icges(5,5) * t381 - Icges(5,6) * t379) * t377;
t291 = qJD(4) * t303;
t292 = (Icges(5,4) * t381 - Icges(5,2) * t379) * t494;
t305 = (Icges(5,1) * t381 - t527) * t377;
t293 = qJD(4) * t305;
t113 = t291 * t378 + (t292 * t381 + t293 * t379 + (-t284 * t379 + t285 * t381) * qJD(4)) * t377;
t166 = t283 * t378 + (t284 * t381 + t285 * t379) * t377;
t506 = t113 * t340 + t166 * t333;
t499 = -t222 - t290;
t174 = -t240 * rSges(6,1) - rSges(6,2) * t241;
t175 = t243 * rSges(6,1) + t245 * rSges(6,2);
t497 = -t284 + t305;
t496 = -Icges(5,2) * t510 + t285 + t348;
t239 = t315 * rSges(6,1) + t317 * rSges(6,2);
t489 = t163 * t483 + t93 * t516 + t94 * t519;
t77 = -t182 * t519 - t588;
t78 = t183 * t519 + t279 * t186 - t189 * t430;
t478 = t519 / 0.2e1;
t477 = -t516 / 0.2e1;
t476 = t514 / 0.2e1;
t475 = -t494 / 0.2e1;
t474 = t494 / 0.2e1;
t19 = t163 * t202 - t165 * t203 - t205 * t249 - t206 * t250 + t273 * t94 + t274 * t93 + (t126 * t356 + t127 * t355) * t494;
t471 = t19 * (t163 * t519 + t165 * t516);
t470 = t273 * t174 + t175 * t274;
t469 = t314 * t175 - t239 * t273;
t468 = -t174 * t314 - t274 * t239;
t458 = t355 * t476;
t457 = t356 * t476;
t456 = t355 * t475;
t455 = t355 * t474;
t454 = t356 * t475;
t453 = t356 * t474;
t263 = t431 * pkin(4);
t276 = t307 + t354;
t344 = rSges(2,1) * t382 - rSges(2,2) * t380;
t343 = rSges(2,1) * t380 + rSges(2,2) * t382;
t167 = t192 + t295;
t445 = -t102 * t355 - t525;
t442 = t190 * t355 + t192 * t356;
t441 = (Icges(5,5) * t431 - Icges(5,6) * t278) * t356 - (Icges(5,5) * t279 + Icges(5,6) * t430) * t355;
t60 = -t156 * t516 - t562;
t437 = -t307 * t363 - t490;
t147 = t354 + t167;
t76 = -t183 * t516 - t505;
t427 = (t355 * t76 - t597) * t377;
t426 = (t355 * t78 - t356 * t77) * t377;
t125 = t465 + t165;
t275 = -t306 - t542;
t264 = t279 * pkin(4);
t116 = t125 + t354;
t420 = -t254 - t490;
t104 = t442 * t494;
t417 = -t447 - t490;
t416 = (Icges(5,1) * t279 - t186 + t528) * t355 - (Icges(5,1) * t431 - t184 - t261) * t356;
t15 = -t139 * t158 - t140 * t161 + t243 * t90 - t245 * t92 + (-t155 * t517 + t355 * t88) * t377;
t16 = t139 * t159 + t140 * t162 + t243 * t89 - t245 * t91 + (t156 * t517 + t355 * t87) * t377;
t17 = -t141 * t158 - t142 * t161 - t240 * t90 + t241 * t92 + (-t155 * t520 - t356 * t88) * t377;
t18 = t141 * t159 + t142 * t162 - t240 * t89 + t241 * t91 + (t156 * t520 - t356 * t87) * t377;
t26 = t273 * t60 - t575;
t50 = t139 * t220 + t140 * t221 + t208 * t243 - t209 * t245 + (t207 * t355 + t219 * t517) * t377;
t51 = t141 * t220 + t142 * t221 - t208 * t240 + t209 * t241 + (-t207 * t356 + t219 * t520) * t377;
t402 = (-t15 * t274 + t16 * t273 + t202 * t62 + t203 * t61 + t302 * t96 + t314 * t50) * t478 + (t397 * t243 + t245 * t585 - t519 * t586) * t550 + (-t397 * t240 - t241 * t585 + t586 * t516) * t547 - (t397 * t315 + t317 * t585 - t378 * t586) * t314 / 0.2e1 + (-t17 * t274 + t18 * t273 + t202 * t60 + t203 * t59 - t302 * t412 + t314 * t51) * t477 + (t378 * t96 + (t355 * t62 - t356 * t61) * t377) * t554 + (-t378 * t412 + (t355 * t60 - t356 * t59) * t377) * t553 + t26 * t458 + t27 * t457 + (t378 * t50 + ((t363 * t62 - t15) * t356 + (t363 * t61 + t16) * t355) * t377) * t549 + t302 * (t378 * t98 + (t355 * t72 - t356 * t71) * t377) / 0.2e1 + (t378 * t51 + ((t363 * t60 - t17) * t356 + (t363 * t59 + t18) * t355) * t377) * t548 + (t531 + t532 + t533 - t534 + t539) * t544 + t314 * (t378 * t65 + ((t363 * t72 - t30) * t356 + (t363 * t71 + t31) * t355) * t377) / 0.2e1;
t401 = -t94 - t498;
t392 = t401 - t490;
t118 = t279 * t284 + t283 * t519 - t285 * t430;
t114 = t118 * t340;
t42 = qJD(4) * t427 - t579;
t43 = qJD(4) * t426 + t114;
t106 = Icges(5,5) * t199 + Icges(5,6) * t198 + Icges(5,3) * t484;
t108 = Icges(5,4) * t199 + Icges(5,2) * t198 + Icges(5,6) * t484;
t110 = Icges(5,1) * t199 + Icges(5,4) * t198 + Icges(5,5) * t484;
t48 = t106 * t378 + (t108 * t381 + t110 * t379 + (-t184 * t379 - t188 * t381) * qJD(4)) * t377;
t105 = Icges(5,5) * t197 + Icges(5,6) * t196 + Icges(5,3) * t483;
t107 = Icges(5,4) * t197 + Icges(5,2) * t196 + Icges(5,6) * t483;
t109 = Icges(5,1) * t197 + Icges(5,4) * t196 + Icges(5,5) * t483;
t49 = t105 * t378 + (t107 * t381 + t109 * t379 + (-t186 * t379 + t189 * t381) * qJD(4)) * t377;
t66 = t196 * t284 + t197 * t285 + t279 * t292 - t430 * t293 + (t283 * t517 + t291 * t355) * t377;
t67 = t198 * t284 + t199 * t285 + t431 * t292 + t278 * t293 + (t283 * t520 - t291 * t356) * t377;
t387 = t539 + t532 / 0.2e1 + t531 / 0.2e1 + t529 / 0.2e1 + t530 / 0.2e1 + Icges(4,3) * t360 + (t114 + ((-t443 + t75 + t78) * t355 + (t76 + (t523 - t524) * t377 - t77 + t505) * t356) * t494) * t453 - t412 * t553 - t411 * t551 + t50 * t549 + t118 * t552 + t96 * t554 - t534 / 0.2e1 + t533 / 0.2e1 + t506 + t27 * t547 + (t26 + (t61 + (t155 * t355 + t156 * t356) * t377 + t562 + t591) * t273 + t575) * t550 + (t51 + t27) * t548 + (t49 + t66) * t455 + (t42 + t579 + (t597 + (t505 + t77 + (t523 + t524) * t377 + t588) * t355) * t494) * t456 + (t48 + t67 + t43) * t454;
t386 = Icges(3,3) * t371 + t387;
t259 = t320 * t374 + t487;
t224 = -t437 + t487;
t218 = rSges(5,1) * t279 + rSges(5,2) * t430;
t217 = rSges(5,1) * t431 - rSges(5,2) * t278;
t204 = t222 * t484;
t146 = t378 * t165;
t84 = t378 * t93;
t57 = t163 * t273 + t165 * t274 + (-t205 * t355 + t206 * t356) * t494;
t37 = t107 * t431 + t109 * t278 + t186 * t198 + t189 * t199 + (-t105 * t356 + t183 * t520) * t377;
t36 = t108 * t431 + t110 * t278 + t184 * t198 - t188 * t199 + (-t106 * t356 - t182 * t520) * t377;
t35 = t107 * t279 - t109 * t430 + t186 * t196 + t189 * t197 + (t105 * t355 + t183 * t517) * t377;
t34 = t108 * t279 - t110 * t430 + t184 * t196 - t188 * t197 + (t106 * t355 - t182 * t517) * t377;
t1 = [Icges(2,3) * qJDD(1) + t386 + (t566 * (t320 + t369) + t567 * (-t319 - t543) + (-t288 - t487 + t259) * t258) * m(3) + (g(1) * t343 - g(2) * t344 + (t343 ^ 2 + t344 ^ 2) * qJDD(1)) * m(2) + (t68 * (t392 - t487) + (-t394 + t68 - t488 + t555) * t69 + t571 * (t116 + t369) + t565 + t603 * (t598 - t543)) * m(6) + (t101 * (t417 - t487) + t570 * (t369 + t147) + (t101 - t400 - t488 + t557) * t102 + t564 + t602 * (t599 - t543)) * m(5) + (t568 * (t276 + t369) + t569 * (t275 - t543) + (t420 - t487 + t224) * t223) * m(4); t386 + ((-t405 + t555) * t69 + (-t396 + t392) * t68 + t571 * t116 + t565 + t603 * t598) * m(6) + (t570 * t147 + (-t418 + t557) * t102 + (-t409 + t417) * t101 + t564 + t602 * t599) * m(5) + (t568 * t276 + t569 * t275 + (-t437 + t420) * t223) * m(4) + (-t258 * t288 - t259 * t521 + (t258 * t374 + t566) * t320 + (t259 * t374 - t567) * t319) * m(3); t387 + ((-t423 + t577) * t69 + (-t408 + t401) * t68 + t571 * t125 + t565 + t603 * t594) * m(6) + (t570 * t167 + (t459 + t576) * t102 + (-t429 - t447) * t101 + t564 + t602 * t593) * m(5) + (-t223 * t254 - t224 * t281 + (t223 * t363 + t568) * t307 + (t224 * t363 - t569) * t306) * m(4); t340 * (t113 * t378 + ((t363 * t86 - t48) * t356 + (t363 * t85 + t49) * t355) * t377) / 0.2e1 + ((t279 * t496 + t303 * t519 - t430 * t497) * t340 + (t279 * t415 - t416 * t430 - t441 * t519) * t494) * t456 + ((t278 * t497 - t303 * t516 + t431 * t496) * t340 + (t416 * t278 + t415 * t431 + t441 * t516) * t494) * t453 + t42 * t458 + t43 * t457 + (t378 * t66 + ((t363 * t78 - t34) * t356 + (t363 * t77 + t35) * t355) * t377) * t455 + (t378 * t67 + ((t363 * t76 - t36) * t356 + (t363 * t75 + t37) * t355) * t377) * t454 + (t529 + t530 + (t355 * t49 - t356 * t48) * t494 + t506) * t544 + t402 + (t118 * t378 + t426) * t552 + t333 * (t166 * t378 + (t355 * t86 - t356 * t85) * t377) / 0.2e1 + (t118 * t333 + t249 * t78 + t250 * t77 + t340 * t66 + (-t34 * t356 + t35 * t355) * t494) * t478 + (-t411 * t333 + t249 * t76 + t250 * t75 + t340 * t67 + (t355 * t37 - t356 * t36) * t494) * t477 + (-t378 * t411 + t427) * t551 - t340 * (t303 * t378 * t340 + ((t379 * t497 + t381 * t496) * t340 + ((t379 * t416 + t381 * t415) * t377 - t441 * t378) * qJD(4)) * t377) / 0.2e1 + (-g(1) * (t264 + t175) - g(2) * (t263 + t174) - g(3) * (pkin(4) * t509 + t239) - t68 * (-t263 * t340 + t468) - t69 * (t264 * t340 + t469) - t57 * ((t263 * t355 + t264 * t356) * t494 + t470) + t68 * t204 + t33 * t146 + t69 * t84 + t471 + t57 * t489 + (t32 * (-t163 + t205) + t68 * (-t127 - t94) + t33 * t206 + t69 * t126) * t378 + ((t32 * t499 - t68 * t210 + t19 * t206 + t57 * t126 + (-t57 * t205 + t499 * t69) * t363) * t356 + (t33 * t499 - t69 * t210 - t19 * t205 + t57 * t127 + (t68 * t290 + t57 * (-t165 - t206)) * t363) * t355) * t377) * m(6) + ((-t101 * t112 + t102 * t111 - t190 * t53 + t192 * t54) * t378 + ((t190 * t249 - t192 * t250) * t442 + t445 * t296 + ((-t102 * t363 - t53) * t356 + (t101 * t363 - t54) * t355) * t289 + (0.2e1 * t111 * t356 + 0.2e1 * t112 * t355 + t190 * t517 - t192 * t520) * t104) * t377 - (-t101 * t217 + t102 * t218) * t340 - (t104 * (t217 * t355 + t218 * t356) + t445 * t310) * t494 - g(1) * t218 - g(2) * t217 - g(3) * t310) * m(5); t402 + (t32 * (-t163 * t378 - t222 * t516) + t33 * (-t222 * t519 + t146) + t471 - g(1) * t175 - g(2) * t174 - g(3) * t239 + (t84 + (-t210 * t355 - t222 * t517) * t377 - t469) * t69 + (-t210 * t516 - t378 * t94 + t204 - t468) * t68 + (-t165 * t484 - t470 + t489) * t57) * m(6);];
tau = t1;
