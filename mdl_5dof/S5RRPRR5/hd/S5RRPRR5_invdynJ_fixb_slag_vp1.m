% Calculate vector of inverse dynamics joint torques for
% S5RRPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR5_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:12
% EndTime: 2022-01-20 11:02:36
% DurationCPUTime: 11.95s
% Computational Cost: add. (19945->665), mult. (13261->822), div. (0->0), fcn. (10315->10), ass. (0->382)
t333 = pkin(9) + qJ(4);
t326 = qJ(5) + t333;
t318 = sin(t326);
t336 = qJ(1) + qJ(2);
t328 = cos(t336);
t502 = t318 * t328;
t451 = rSges(6,2) * t502;
t327 = sin(t336);
t319 = cos(t326);
t500 = t319 * t328;
t468 = rSges(6,1) * t500 + t327 * rSges(6,3);
t172 = -t451 + t468;
t338 = cos(pkin(9));
t320 = t338 * pkin(3) + pkin(2);
t275 = t328 * t320;
t339 = -pkin(7) - qJ(3);
t293 = t327 * t339;
t463 = t293 - t275;
t325 = cos(t333);
t304 = pkin(4) * t325;
t257 = t304 + t320;
t332 = pkin(8) - t339;
t560 = t328 * t257 + t332 * t327;
t570 = t560 + t463;
t489 = t570 + t172;
t302 = t327 * qJ(3);
t252 = t328 * pkin(2) + t302;
t164 = -t252 - t463;
t482 = t164 + t252;
t574 = t482 + t489;
t324 = sin(t333);
t573 = rSges(5,2) * t324;
t295 = Icges(6,4) * t319;
t229 = Icges(6,1) * t318 + t295;
t399 = -Icges(6,2) * t318 + t295;
t572 = t229 + t399;
t340 = sin(qJ(1));
t524 = pkin(1) * qJD(1);
t449 = t340 * t524;
t251 = rSges(3,1) * t327 + rSges(3,2) * t328;
t335 = qJD(1) + qJD(2);
t504 = t251 * t335;
t197 = -t449 - t504;
t334 = qJD(4) + qJD(5);
t525 = rSges(6,2) * t319;
t450 = t334 * t525;
t528 = rSges(6,1) * t318;
t571 = -t334 * t528 - t450;
t249 = t328 * t334;
t503 = t318 * t327;
t264 = rSges(6,2) * t503;
t491 = t328 * t335;
t475 = rSges(6,3) * t491 + t335 * t264;
t494 = t327 * t335;
t101 = -t328 * t450 + (-t249 * t318 - t319 * t494) * rSges(6,1) + t475;
t441 = t571 * t327 - t335 * t451;
t102 = t335 * t468 + t441;
t490 = t332 * t335;
t266 = t328 * t490;
t294 = t328 * t339;
t459 = qJD(4) * t328;
t438 = t324 * t459;
t416 = pkin(4) * t438;
t470 = t257 - t320;
t104 = -t416 + t266 + (-t327 * t470 + t294) * t335;
t270 = t335 * t293;
t460 = qJD(4) * t327;
t439 = t324 * t460;
t263 = pkin(4) * t439;
t467 = t327 * t490 - t263;
t105 = t470 * t491 + t270 + t467;
t292 = t332 * t328;
t466 = t327 * t320 + t294;
t133 = t257 * t327 - t292 - t466;
t458 = qJD(4) * t335;
t216 = qJDD(4) * t327 + t328 * t458;
t457 = qJD(5) * t335;
t159 = qJDD(5) * t327 + t328 * t457 + t216;
t278 = t327 * t458;
t160 = t327 * t457 + t278 + (-qJDD(4) - qJDD(5)) * t328;
t469 = t328 * rSges(6,3) + t264;
t501 = t319 * t327;
t171 = rSges(6,1) * t501 - t469;
t217 = -qJDD(4) * t328 + t278;
t248 = t327 * t334;
t10 = t101 * t249 + t102 * t248 + t133 * t216 - t570 * t217 + t159 * t171 - t160 * t172 + (t104 * t328 + t105 * t327) * qJD(4);
t231 = t525 + t528;
t195 = t231 * t327;
t196 = t231 * t328;
t296 = t319 * rSges(6,1);
t526 = rSges(6,2) * t318;
t232 = t296 - t526;
t58 = t171 * t248 + t172 * t249 + (t133 * t327 + t328 * t570) * qJD(4);
t300 = qJD(3) * t327;
t375 = -t231 * t249 + t300 - t416;
t303 = t328 * qJ(3);
t250 = pkin(2) * t327 - t303;
t163 = t250 - t466;
t483 = t163 - t250;
t415 = -t133 - t171 + t483;
t59 = t335 * t415 + t375 - t449;
t301 = qJD(3) * t328;
t341 = cos(qJ(1));
t448 = t341 * t524;
t411 = -t301 + t448;
t561 = -t248 * t231 - t263;
t388 = t411 + t561;
t60 = t574 * t335 + t388;
t569 = -t59 * (t195 * t335 - t249 * t232) - t58 * (-t248 * t195 - t196 * t249) - t60 * (-t335 * t196 - t232 * t248) + t10 * (t327 * t171 + t328 * t172);
t496 = t325 * t328;
t390 = rSges(5,1) * t496 + t327 * rSges(5,3);
t445 = t324 * t491;
t461 = qJD(4) * t325;
t447 = rSges(5,2) * t461;
t440 = -rSges(5,1) * t439 - rSges(5,2) * t445 - t327 * t447;
t120 = t335 * t390 + t440;
t247 = rSges(5,1) * t325 - t573;
t224 = t247 * qJD(4);
t246 = rSges(5,1) * t324 + rSges(5,2) * t325;
t331 = qJDD(1) + qJDD(2);
t343 = qJD(1) ^ 2;
t380 = (-qJDD(1) * t340 - t341 * t343) * pkin(1);
t462 = qJD(3) * t335;
t370 = qJDD(3) * t327 + t328 * t462 + t380;
t497 = t325 * t327;
t454 = rSges(5,1) * t497;
t499 = t324 * t327;
t184 = -rSges(5,2) * t499 - t328 * rSges(5,3) + t454;
t444 = -t184 + t483;
t186 = t252 * t335 - t301;
t530 = pkin(2) - t320;
t488 = t270 - (-t328 * t530 - t302) * t335 - t186;
t34 = -t224 * t459 + t217 * t246 + (-t120 + t488) * t335 + t444 * t331 + t370;
t568 = t34 - g(1);
t387 = rSges(5,3) * t491 - t328 * t447 + t494 * t573;
t119 = (-t325 * t494 - t438) * rSges(5,1) + t387;
t498 = t324 * t328;
t185 = -rSges(5,2) * t498 + t390;
t282 = qJ(3) * t491;
t330 = t341 * pkin(1);
t532 = pkin(1) * t340;
t413 = qJDD(1) * t330 - t343 * t532;
t465 = t282 + t300;
t362 = -qJDD(3) * t328 + t335 * (-pkin(2) * t494 + t465) + t331 * t252 + t327 * t462 + t413;
t357 = t335 * (-t282 + (t327 * t530 - t294) * t335) + t331 * t164 + t362;
t35 = t119 * t335 + t185 * t331 - t216 * t246 - t224 * t460 + t357;
t567 = t35 - g(2);
t337 = sin(pkin(9));
t527 = rSges(4,2) * t337;
t452 = t328 * t527;
t256 = t335 * t452;
t529 = rSges(4,1) * t338;
t455 = t327 * t529;
t289 = t327 * t527;
t464 = t328 * rSges(4,3) + t289;
t193 = t455 - t464;
t472 = -t250 - t193;
t558 = -t327 * rSges(4,3) - t328 * t529;
t63 = t472 * t331 + (t335 * t558 - t186 + t256) * t335 + t370;
t566 = t63 - g(1);
t194 = -t452 - t558;
t471 = rSges(4,3) * t491 + t335 * t289;
t64 = t331 * t194 + t335 * (-t335 * t455 + t471) + t362;
t565 = t64 - g(2);
t213 = rSges(3,1) * t491 - rSges(3,2) * t494;
t564 = -t213 * t335 - t251 * t331 - g(1) + t380;
t253 = t328 * rSges(3,1) - rSges(3,2) * t327;
t563 = t253 * t331 - t335 * t504 - g(2) + t413;
t149 = t194 + t252;
t562 = t149 * t335;
t220 = t335 * t250;
t559 = -t335 * t163 + t220;
t299 = Icges(5,4) * t325;
t400 = -Icges(5,2) * t324 + t299;
t244 = Icges(5,1) * t324 + t299;
t202 = t232 * t334;
t495 = t325 * qJD(4) ^ 2;
t13 = t160 * t231 - t202 * t249 + (t217 * t324 - t328 * t495) * pkin(4) + (-t102 - t105 + t488) * t335 + t415 * t331 + t370;
t14 = -t159 * t231 - t202 * t248 + (t101 + t104) * t335 + t489 * t331 + (-t216 * t324 - t327 * t495) * pkin(4) + t357;
t429 = -t257 - t296;
t523 = t335 * t59;
t531 = pkin(4) * t324;
t557 = (t13 * t429 + (-t59 * rSges(6,3) + t429 * t60) * t335) * t327 + (t13 * t332 - t14 * t526 + t60 * (-qJD(4) * t531 + t571) + t429 * t523) * t328 + t574 * t523;
t556 = t335 * t193 + t220 + t471;
t156 = t335 * t171;
t555 = t133 * t335 + t156 + t266 - t375 + t475 + t559;
t517 = Icges(6,4) * t318;
t230 = Icges(6,1) * t319 - t517;
t385 = t230 * t328;
t170 = Icges(6,5) * t327 + t385;
t226 = Icges(6,5) * t319 - Icges(6,6) * t318;
t225 = Icges(6,5) * t318 + Icges(6,6) * t319;
t366 = Icges(6,3) * t335 - t225 * t334;
t383 = t399 * t328;
t168 = Icges(6,6) * t327 + t383;
t513 = t168 * t318;
t554 = -t226 * t494 + t328 * t366 + t335 * (-t170 * t319 + t513);
t381 = t226 * t328;
t167 = Icges(6,4) * t501 - Icges(6,2) * t503 - Icges(6,6) * t328;
t262 = Icges(6,4) * t503;
t169 = Icges(6,1) * t501 - Icges(6,5) * t328 - t262;
t398 = t167 * t318 - t169 * t319;
t553 = t327 * t366 + (t381 + t398) * t335;
t241 = Icges(5,5) * t325 - Icges(5,6) * t324;
t240 = Icges(5,5) * t324 + Icges(5,6) * t325;
t363 = Icges(5,3) * t335 - qJD(4) * t240;
t518 = Icges(5,4) * t324;
t245 = Icges(5,1) * t325 - t518;
t386 = t245 * t328;
t181 = Icges(5,5) * t327 + t386;
t384 = t400 * t328;
t179 = Icges(5,6) * t327 + t384;
t511 = t179 * t324;
t395 = -t181 * t325 + t511;
t552 = -t241 * t494 + t328 * t363 + t335 * t395;
t382 = t241 * t328;
t273 = Icges(5,4) * t499;
t180 = Icges(5,1) * t497 - Icges(5,5) * t328 - t273;
t178 = Icges(5,4) * t497 - Icges(5,2) * t499 - Icges(5,6) * t328;
t512 = t178 * t324;
t396 = -t180 * t325 + t512;
t551 = t327 * t363 + (t382 + t396) * t335;
t227 = Icges(6,2) * t319 + t517;
t393 = t227 * t318 - t229 * t319;
t550 = t226 * t334 + t335 * t393;
t242 = Icges(5,2) * t325 + t518;
t391 = t242 * t324 - t244 * t325;
t549 = t241 * qJD(4) + t335 * t391;
t176 = Icges(5,5) * t497 - Icges(5,6) * t499 - Icges(5,3) * t328;
t75 = -t328 * t176 - t327 * t396;
t162 = t335 * t184;
t548 = -rSges(5,1) * t438 + t162 + t300 + t387 + t559;
t547 = t327 * (-t242 * t328 + t181) - t328 * (-Icges(5,2) * t497 + t180 - t273);
t546 = t248 * (-t227 * t328 + t170) - t249 * (-Icges(6,2) * t501 + t169 - t262) + t335 * t572;
t545 = t159 / 0.2e1;
t544 = t160 / 0.2e1;
t543 = t216 / 0.2e1;
t542 = t217 / 0.2e1;
t541 = -t248 / 0.2e1;
t540 = t248 / 0.2e1;
t539 = -t249 / 0.2e1;
t538 = t249 / 0.2e1;
t537 = t327 / 0.2e1;
t536 = -t328 / 0.2e1;
t535 = t331 / 0.2e1;
t534 = -t335 / 0.2e1;
t533 = t335 / 0.2e1;
t408 = -t246 * t459 + t300;
t374 = t408 - t449;
t72 = t335 * t444 + t374;
t522 = t335 * t72;
t509 = t225 * t328;
t91 = -t327 * t393 - t509;
t521 = t91 * t335;
t506 = t240 * t328;
t109 = -t327 * t391 - t506;
t514 = t109 * t335;
t510 = t225 * t327;
t508 = t227 * t334;
t507 = t240 * t327;
t505 = t241 * t335;
t165 = Icges(6,5) * t501 - Icges(6,6) * t503 - Icges(6,3) * t328;
t493 = t328 * t165;
t487 = -t327 * t165 - t169 * t500;
t166 = Icges(6,3) * t327 + t381;
t486 = t327 * t166 + t170 * t500;
t485 = -t327 * t176 - t180 * t496;
t177 = Icges(5,3) * t327 + t382;
t484 = t327 * t177 + t181 * t496;
t474 = -t242 + t245;
t473 = t244 + t400;
t456 = t327 * t102 + (t101 + t156) * t328;
t442 = t185 + t482;
t437 = t494 / 0.2e1;
t436 = t491 / 0.2e1;
t435 = -pkin(2) - t529;
t434 = -t460 / 0.2e1;
t433 = t460 / 0.2e1;
t432 = -t459 / 0.2e1;
t431 = t459 / 0.2e1;
t378 = -t231 - t531;
t368 = Icges(6,5) * t335 - t229 * t334;
t428 = -t168 * t334 - t230 * t494 + t328 * t368;
t367 = Icges(6,6) * t335 - t508;
t427 = t169 * t334 + t327 * t367 + t335 * t383;
t426 = t170 * t334 + t328 * t367 - t399 * t494;
t425 = -t167 * t334 + t327 * t368 + t335 * t385;
t138 = t170 * t501;
t424 = t328 * t166 - t138;
t145 = t181 * t497;
t423 = t328 * t177 - t145;
t422 = -t165 + t513;
t420 = -t176 + t511;
t419 = t572 * t334;
t418 = t230 * t334 - t508;
t414 = t468 + t560;
t412 = t300 - t449;
t410 = -pkin(4) * t461 - t202;
t288 = rSges(2,1) * t341 - rSges(2,2) * t340;
t287 = rSges(2,1) * t340 + rSges(2,2) * t341;
t407 = -t327 * t60 - t328 * t59;
t211 = t246 * t460;
t73 = t335 * t442 - t211 + t411;
t406 = -t327 * t73 - t328 * t72;
t76 = -t179 * t499 - t423;
t405 = t327 * t76 - t328 * t75;
t77 = -t178 * t498 - t485;
t78 = -t179 * t498 + t484;
t404 = t327 * t78 - t328 * t77;
t403 = t270 + t301 - t440;
t93 = t167 * t319 + t169 * t318;
t107 = t178 * t325 + t180 * t324;
t108 = t179 * t325 + t181 * t324;
t394 = t184 * t327 + t185 * t328;
t392 = t242 * t325 + t244 * t324;
t389 = t301 - t441 - t467;
t379 = t398 * t327;
t130 = -t184 - t466;
t131 = t185 - t463;
t377 = -t226 * t335 + t248 * t509 - t249 * t510;
t373 = t178 * t328 - t179 * t327;
t148 = t327 * t435 + t303 + t464;
t372 = t327 * t429 + t292 + t469;
t371 = t172 + t560;
t369 = (-t324 * t473 + t325 * t474) * t335;
t365 = Icges(5,5) * t335 - qJD(4) * t244;
t364 = Icges(5,6) * t335 - qJD(4) * t242;
t354 = t165 * t335 - t318 * t427 + t319 * t425;
t16 = t327 * t553 + t354 * t328;
t355 = t166 * t335 - t318 * t426 + t319 * t428;
t17 = t327 * t554 + t355 * t328;
t18 = t354 * t327 - t328 * t553;
t19 = t355 * t327 - t328 * t554;
t68 = -t379 - t493;
t69 = -t168 * t503 - t424;
t30 = t248 * t69 - t249 * t68 + t521;
t70 = -t167 * t502 - t487;
t71 = -t168 * t502 + t486;
t92 = -t328 * t393 + t510;
t88 = t92 * t335;
t31 = t248 * t71 - t249 * t70 + t88;
t356 = (-t229 * t328 - t168) * t248 - (-t229 * t327 - t167) * t249 + (-t227 + t230) * t335;
t345 = -t318 * t546 + t356 * t319;
t353 = t225 * t335 - t318 * t419 + t319 * t418;
t44 = t327 * t550 + t353 * t328;
t45 = t353 * t327 - t328 * t550;
t46 = t318 * t425 + t319 * t427;
t47 = t318 * t428 + t319 * t426;
t94 = t168 * t319 + t170 * t318;
t361 = (t159 * t71 - t16 * t249 + t160 * t70 + t17 * t248 + t331 * t92 + t335 * t44) * t537 + (-t327 * t377 + t328 * t345) * t541 + (t327 * t345 + t328 * t377) * t538 + (t159 * t69 + t160 * t68 - t18 * t249 + t19 * t248 + t331 * t91 + t335 * t45) * t536 + (t356 * t318 + t319 * t546) * t534 + t30 * t437 + t31 * t436 + ((t335 * t71 - t16) * t328 + (t335 * t70 + t17) * t327) * t540 + (t327 * t71 - t328 * t70) * t545 + (t327 * t69 - t328 * t68) * t544 + ((t335 * t69 - t18) * t328 + (t335 * t68 + t19) * t327) * t539 + (t327 * t94 - t328 * t93) * t535 + ((t335 * t94 - t46) * t328 + (t335 * t93 + t47) * t327) * t533;
t115 = t328 * t364 - t400 * t494;
t117 = -t245 * t494 + t328 * t365;
t352 = -qJD(4) * t108 - t115 * t324 + t117 * t325 + t177 * t335;
t116 = t327 * t364 + t335 * t384;
t118 = t327 * t365 + t335 * t386;
t351 = -qJD(4) * t107 - t116 * t324 + t118 * t325 + t176 * t335;
t222 = t400 * qJD(4);
t223 = t245 * qJD(4);
t350 = -qJD(4) * t392 - t222 * t324 + t223 * t325 + t240 * t335;
t121 = t335 * t472 + t412;
t122 = t411 + t562;
t349 = (t121 * t435 * t328 + (t121 * (-rSges(4,3) - qJ(3)) + t122 * t435) * t327) * t335;
t348 = -t324 * t547 + t373 * t325;
t347 = (t72 * (-t390 - t275) + t73 * (-t454 - t466)) * t335;
t110 = -t328 * t391 + t507;
t103 = t110 * t335;
t36 = qJD(4) * t405 + t514;
t37 = qJD(4) * t404 + t103;
t50 = -qJD(4) * t396 + t116 * t325 + t118 * t324;
t51 = -qJD(4) * t395 + t115 * t325 + t117 * t324;
t54 = t327 * t549 + t350 * t328;
t55 = t350 * t327 - t328 * t549;
t346 = (t88 + (t69 + (t167 * t328 + t168 * t327) * t318 + t424 + t487) * t249 + (-t169 * t501 + t493 + t68 + (t167 * t327 - t168 * t328) * t318 + t486) * t248) * t538 + (t103 + ((t76 - t145 + (t177 + t512) * t328 + t485) * t328 + t484 * t327) * qJD(4)) * t431 + (t94 + t92) * t545 + (t93 + t91) * t544 + (t110 + t108) * t543 + (t109 + t107) * t542 + (t30 - t521 + (t71 - t379 - t486) * t249 + (t422 * t327 - t138 + t70) * t248 + ((t166 + t398) * t248 + t422 * t249) * t328) * t541 + (t47 + t44) * t540 + (t36 - t514 + ((t328 * t420 - t484 + t78) * t328 + (t327 * t420 + t423 + t77) * t327) * qJD(4)) * t434 + (t51 + t54) * t433 + (-qJD(4) * t391 + t222 * t325 + t223 * t324 + t318 * t418 + t319 * t419) * t335 + (t46 + t45 + t31) * t539 + (t50 + t55 + t37) * t432 + (t392 + Icges(4,2) * t338 ^ 2 + (Icges(4,1) * t337 + 0.2e1 * Icges(4,4) * t338) * t337 + Icges(3,3) + t227 * t319 + t229 * t318) * t331;
t210 = t246 * t328;
t209 = t246 * t327;
t198 = t253 * t335 + t448;
t106 = t394 * qJD(4);
t23 = t352 * t327 - t328 * t552;
t22 = t351 * t327 - t328 * t551;
t21 = t327 * t552 + t352 * t328;
t20 = t327 * t551 + t351 * t328;
t1 = [Icges(2,3) * qJDD(1) + t346 + (t563 * (t253 + t330) + t564 * (-t251 - t532) + (-t213 - t448 + t198) * t197) * m(3) + ((t287 ^ 2 + t288 ^ 2) * qJDD(1) + g(1) * t287 - g(2) * t288) * m(2) + (t13 * (t469 - t532) + t14 * (t330 + t414) - g(1) * (t372 - t532) - g(2) * (t330 + t371) + (t412 + t449 + t555) * t60 + (t389 - t448 + t388) * t59 + t557) * m(6) + (t72 * (t403 - t448) + t347 + (-t374 + t72 - t449 + t548) * t73 + t567 * (t131 + t330) + t568 * (t130 - t532)) * m(5) + (t121 * (t256 - t411) + t349 + t565 * (t149 + t330) + t566 * (t148 - t532) + (t121 + t282 + t556) * t122) * m(4); t346 + (-g(1) * t372 - g(2) * t371 + t13 * t469 + t14 * t414 + (t300 + t555) * t60 + (t389 - t301 + t561) * t59 + t557) * m(6) + (t442 * t522 + t347 + (-t408 + t548) * t73 + (-t211 - t301 + t403) * t72 + t567 * t131 + t568 * t130) * m(5) + (t349 + t565 * t149 + t566 * t148 + (t465 - t300 + t556) * t122 + (t256 + t562) * t121) * m(4) + (-t197 * t213 - t198 * t504 + (t197 * t335 + t563) * t253 + (t198 * t335 - t564) * t251) * m(3); (-m(4) - m(5) - m(6)) * (g(1) * t327 - g(2) * t328) + 0.2e1 * (t13 * t537 + t14 * t536) * m(6) + 0.2e1 * (t34 * t537 + t35 * t536) * m(5) + 0.2e1 * (t536 * t64 + t537 * t63) * m(4); ((t324 * t474 + t325 * t473) * t335 + (t373 * t324 + t325 * t547) * qJD(4)) * t534 + (t110 * t331 + t216 * t78 + t217 * t77 + t335 * t54 + (-t20 * t328 + t21 * t327) * qJD(4)) * t537 + (t109 * t331 + t216 * t76 + t217 * t75 + t335 * t55 + (-t22 * t328 + t23 * t327) * qJD(4)) * t536 + ((-t459 * t507 - t505) * t328 + (t369 + (t328 * t506 + t348) * qJD(4)) * t327) * t431 + ((-t460 * t506 + t505) * t327 + (t369 + (t327 * t507 + t348) * qJD(4)) * t328) * t434 + ((t108 * t335 - t50) * t328 + (t107 * t335 + t51) * t327) * t533 + t361 + (-t107 * t328 + t108 * t327) * t535 + ((t335 * t76 - t22) * t328 + (t335 * t75 + t23) * t327) * t432 + ((t335 * t78 - t20) * t328 + (t335 * t77 + t21) * t327) * t433 + t36 * t437 + t37 * t436 + t404 * t543 + t405 * t542 + (-(-t60 * t445 + (t407 * t325 + t58 * (-t327 ^ 2 - t328 ^ 2) * t324) * qJD(4)) * pkin(4) - g(3) * (t232 + t304) - (g(1) * t328 + g(2) * t327) * t378 + t58 * t456 + (t13 * t378 + t59 * t410 + t10 * t570 + t58 * t104 + (t58 * t133 + t378 * t60) * t335) * t328 + (t14 * t378 + t60 * t410 + t10 * t133 + t58 * t105 + (t59 * t231 - t489 * t58) * t335) * t327 + t569) * m(6) + (-(t209 * t72 - t210 * t73) * t335 - (t106 * (-t209 * t327 - t210 * t328) + t406 * t247) * qJD(4) + (t184 * t216 - t185 * t217 + (t119 * t328 + t120 * t327) * qJD(4)) * t394 + t106 * ((t119 + t162) * t328 + (-t185 * t335 + t120) * t327) + t406 * t224 + ((-t335 * t73 - t34) * t328 + (-t35 + t522) * t327) * t246 + g(1) * t210 + g(2) * t209 - g(3) * t247) * m(5); t361 + (t58 * (-t172 * t494 + t456) + t407 * t202 + ((-t335 * t60 - t13) * t328 + (-t14 + t523) * t327) * t231 + g(1) * t196 + g(2) * t195 - g(3) * t232 + t569) * m(6);];
tau = t1;
