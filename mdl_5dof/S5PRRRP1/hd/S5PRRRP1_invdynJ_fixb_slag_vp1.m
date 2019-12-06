% Calculate vector of inverse dynamics joint torques for
% S5PRRRP1
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP1_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:39:47
% EndTime: 2019-12-05 16:40:10
% DurationCPUTime: 13.70s
% Computational Cost: add. (12989->458), mult. (9330->552), div. (0->0), fcn. (7247->6), ass. (0->264)
t267 = sin(qJ(4));
t268 = cos(qJ(4));
t218 = Icges(6,5) * t268 - Icges(6,6) * t267;
t220 = Icges(5,5) * t268 - Icges(5,6) * t267;
t496 = t218 + t220;
t524 = Icges(5,5) + Icges(6,5);
t523 = Icges(5,6) + Icges(6,6);
t522 = Icges(5,3) + Icges(6,3);
t264 = pkin(8) + qJ(2);
t258 = qJ(3) + t264;
t253 = cos(t258);
t252 = sin(t258);
t392 = t252 * t268;
t393 = t252 * t267;
t132 = Icges(6,4) * t392 - Icges(6,2) * t393 - Icges(6,6) * t253;
t134 = Icges(5,4) * t392 - Icges(5,2) * t393 - Icges(5,6) * t253;
t515 = t132 + t134;
t211 = Icges(6,4) * t393;
t136 = Icges(6,1) * t392 - Icges(6,5) * t253 - t211;
t212 = Icges(5,4) * t393;
t138 = Icges(5,1) * t392 - Icges(5,5) * t253 - t212;
t520 = t136 + t138;
t413 = Icges(6,4) * t267;
t226 = Icges(6,1) * t268 - t413;
t302 = t226 * t253;
t137 = Icges(6,5) * t252 + t302;
t414 = Icges(5,4) * t267;
t228 = Icges(5,1) * t268 - t414;
t303 = t228 * t253;
t139 = Icges(5,5) * t252 + t303;
t513 = t137 + t139;
t221 = Icges(6,2) * t268 + t413;
t223 = Icges(5,2) * t268 + t414;
t519 = t221 + t223;
t259 = Icges(6,4) * t268;
t225 = Icges(6,1) * t267 + t259;
t260 = Icges(5,4) * t268;
t227 = Icges(5,1) * t267 + t260;
t518 = -t225 - t227;
t521 = t496 * t253;
t500 = t522 * t253 - t524 * t392 + t523 * t393;
t499 = t522 * t252 + t521;
t316 = -Icges(6,2) * t267 + t259;
t300 = t316 * t253;
t133 = Icges(6,6) * t252 + t300;
t317 = -Icges(5,2) * t267 + t260;
t301 = t317 * t253;
t135 = Icges(5,6) * t252 + t301;
t514 = t133 + t135;
t517 = t515 * t267;
t516 = t513 * t392;
t217 = Icges(6,5) * t267 + Icges(6,6) * t268;
t219 = Icges(5,5) * t267 + Icges(5,6) * t268;
t512 = t217 + t219;
t511 = t226 + t228;
t265 = qJD(2) + qJD(3);
t510 = -t519 * qJD(4) + t523 * t265;
t509 = t518 * qJD(4) + t524 * t265;
t491 = -t520 * t268 + t517;
t508 = t316 + t317;
t494 = t519 * t267 + t518 * t268;
t507 = -t499 * t253 + t516;
t389 = t253 * t268;
t506 = t500 * t252 - t520 * t389;
t453 = t499 * t252 + t513 * t389;
t505 = t514 * t267;
t465 = -t491 * t252 + t500 * t253;
t464 = -t514 * t393 + t507;
t390 = t253 * t267;
t463 = -t515 * t390 - t506;
t462 = -t514 * t390 + t453;
t461 = t520 * t267 + t515 * t268;
t460 = t513 * t267 + t514 * t268;
t395 = t252 * t265;
t504 = t510 * t253 - t508 * t395;
t503 = (t300 + t301) * t265 + t510 * t252;
t502 = t509 * t253 - t511 * t395;
t501 = (-t302 - t303) * t265 - t509 * t252;
t498 = t508 * qJD(4);
t497 = t511 * qJD(4);
t495 = -t512 * qJD(4) + t522 * t265;
t493 = t518 * t267 - t519 * t268;
t492 = -t513 * t268 + t505;
t397 = t219 * t253;
t400 = t217 * t253;
t459 = -t494 * t252 - t397 - t400;
t398 = t219 * t252;
t401 = t217 * t252;
t458 = -t494 * t253 + t398 + t401;
t266 = -qJ(5) - pkin(7);
t229 = t253 * t266;
t262 = t268 * pkin(4);
t254 = t262 + pkin(3);
t105 = -rSges(6,1) * t392 + rSges(6,2) * t393 + t253 * rSges(6,3) - t252 * t254 - t229;
t490 = t493 * qJD(4) + t512 * t265 - t498 * t267 + t497 * t268;
t489 = -t460 * qJD(4) + t499 * t265 - t504 * t267 + t502 * t268;
t488 = t461 * qJD(4) + t500 * t265 + t503 * t267 + t501 * t268;
t487 = t462 * t252 - t463 * t253;
t486 = t464 * t252 - t465 * t253;
t485 = t496 * qJD(4) + t494 * t265;
t484 = (t491 + t521) * t265 + t495 * t252;
t483 = t495 * t253 + t492 * t265 - t496 * t395;
t482 = t458 * t265;
t249 = t253 * pkin(7);
t178 = pkin(3) * t252 - t249;
t481 = t178 + t105;
t480 = t459 * t265;
t479 = -t252 * t484 + t253 * t488;
t478 = t252 * t483 + t253 * t489;
t248 = t252 * pkin(7);
t179 = t253 * pkin(3) + t248;
t149 = t179 * t265;
t361 = qJD(4) * t265;
t169 = -qJDD(4) * t253 + t252 * t361;
t261 = t268 * rSges(6,1);
t232 = -rSges(6,2) * t267 + t261;
t196 = t232 * qJD(4);
t418 = t268 * rSges(6,2);
t230 = rSges(6,1) * t267 + t418;
t235 = qJD(5) * t253;
t263 = qJDD(2) + qJDD(3);
t256 = sin(t264);
t257 = cos(t264);
t270 = qJD(2) ^ 2;
t297 = (-qJDD(2) * t256 - t257 * t270) * pkin(2);
t352 = -t178 + t481;
t362 = qJD(4) * t253;
t387 = t268 * qJD(4) ^ 2;
t425 = pkin(3) - t254;
t360 = qJD(4) * t267;
t345 = t252 * t360;
t369 = pkin(4) * t345 + t235;
t394 = t252 * t266;
t359 = qJD(4) * t268;
t388 = t265 * t267;
t454 = t252 * t359 + t253 * t388;
t445 = -rSges(6,1) * t345 - t454 * rSges(6,2) - t265 * t394 - t369;
t449 = rSges(6,1) * t389 + t252 * rSges(6,3);
t423 = t445 + (-t253 * t425 - t248 + t449) * t265;
t14 = -t196 * t362 + qJDD(5) * t252 + t169 * t230 + (t169 * t267 - t253 * t387) * pkin(4) + t297 + t352 * t263 + (-t149 + t235 - t423) * t265;
t452 = t14 - g(1);
t168 = qJDD(4) * t252 + t253 * t361;
t391 = t253 * t265;
t206 = pkin(7) * t391;
t251 = pkin(2) * t257;
t426 = pkin(2) * t256;
t327 = qJDD(2) * t251 - t270 * t426;
t305 = t265 * (-pkin(3) * t395 + t206) + t263 * t179 + t327;
t334 = -pkin(4) * t267 - t230;
t106 = -rSges(6,2) * t390 + t253 * t254 - t394 + t449;
t381 = -t179 + t106;
t343 = t253 * t360;
t354 = t265 * t392;
t295 = -t343 - t354;
t342 = t253 * t359;
t234 = qJD(5) * t252;
t355 = t252 * t388;
t446 = rSges(6,2) * t355 + rSges(6,3) * t391 + t234;
t424 = -pkin(4) * t343 - t206 + (t252 * t425 - t229) * t265 + rSges(6,1) * t295 - rSges(6,2) * t342 + t446;
t15 = -qJDD(5) * t253 + t424 * t265 + t381 * t263 + t334 * t168 + (-pkin(4) * t387 - qJD(4) * t196 + qJD(5) * t265) * t252 + t305;
t451 = t15 - g(2);
t422 = rSges(5,1) * t268;
t233 = -rSges(5,2) * t267 + t422;
t197 = t233 * qJD(4);
t231 = rSges(5,1) * t267 + rSges(5,2) * t268;
t368 = rSges(5,2) * t393 + t253 * rSges(5,3);
t141 = rSges(5,1) * t392 - t368;
t372 = -t141 - t178;
t348 = rSges(5,1) * t345 + t454 * rSges(5,2);
t448 = rSges(5,1) * t389 + t252 * rSges(5,3);
t96 = t265 * t448 - t348;
t26 = -t197 * t362 + t169 * t231 + (-t149 - t96) * t265 + t372 * t263 + t297;
t477 = t26 - g(1);
t143 = -rSges(5,2) * t390 + t448;
t363 = qJD(4) * t252;
t304 = rSges(5,3) * t391 + (-t342 + t355) * rSges(5,2);
t94 = rSges(5,1) * t295 + t304;
t27 = t143 * t263 - t168 * t231 - t197 * t363 + t265 * t94 + t305;
t476 = t27 - g(2);
t148 = rSges(4,1) * t391 - rSges(4,2) * t395;
t176 = rSges(4,1) * t252 + rSges(4,2) * t253;
t475 = -t148 * t265 - t176 * t263 - g(1) + t297;
t177 = t253 * rSges(4,1) - rSges(4,2) * t252;
t402 = t176 * t265;
t474 = t177 * t263 - t265 * t402 - g(2) + t327;
t473 = t252 * t488 + t253 * t484;
t472 = t252 * t489 - t253 * t483;
t471 = qJD(4) * t486 + t480;
t470 = qJD(4) * t487 + t482;
t469 = t491 * qJD(4) + t501 * t267 - t503 * t268;
t468 = -t492 * qJD(4) + t502 * t267 + t504 * t268;
t467 = t252 * t485 + t253 * t490;
t466 = t252 * t490 - t253 * t485;
t125 = t265 * t141;
t171 = t265 * t178;
t457 = -rSges(5,1) * t343 + t125 + t171 + t206 + t304;
t456 = -rSges(6,1) * t354 - t254 * t395 - t265 * t481 + t171 + t446;
t455 = t500 + t505;
t420 = pkin(2) * qJD(2);
t357 = t256 * t420;
t144 = -t357 - t402;
t447 = t232 + t262;
t110 = t143 + t179;
t439 = -t110 * t265 + t231 * t363;
t436 = t230 * t363 - t265 * (t179 + t381) + t369;
t374 = -Icges(5,2) * t392 + t138 - t212;
t378 = t227 * t252 + t134;
t435 = -t267 * t374 - t268 * t378;
t376 = -Icges(6,2) * t392 + t136 - t211;
t380 = t225 * t252 + t132;
t434 = -t267 * t376 - t268 * t380;
t433 = t168 / 0.2e1;
t432 = t169 / 0.2e1;
t346 = t231 * t362;
t296 = -t346 - t357;
t64 = t265 * t372 + t296;
t419 = t253 * t64;
t399 = t218 * t265;
t396 = t220 * t265;
t379 = -t225 * t253 - t133;
t377 = -t227 * t253 - t135;
t375 = -t221 * t253 + t137;
t373 = -t223 * t253 + t139;
t367 = -t221 + t226;
t366 = t225 + t316;
t365 = -t223 + t228;
t364 = t227 + t317;
t358 = m(2) + m(3) + m(4);
t356 = t257 * t420;
t339 = -pkin(3) - t422;
t338 = -t363 / 0.2e1;
t337 = t363 / 0.2e1;
t336 = -t362 / 0.2e1;
t335 = t362 / 0.2e1;
t326 = -pkin(4) * t390 - t230 * t253;
t325 = -pkin(4) * t359 - t196;
t181 = rSges(3,1) * t257 - rSges(3,2) * t256;
t180 = rSges(3,1) * t256 + rSges(3,2) * t257;
t65 = t356 - t439;
t320 = -t252 * t65 - t419;
t311 = t141 * t252 + t143 * t253;
t306 = -t418 + (-rSges(6,1) - pkin(4)) * t267;
t294 = t334 * t362 + t234;
t293 = -t267 * t375 + t268 * t379;
t292 = -t267 * t373 + t268 * t377;
t109 = t252 * t339 + t249 + t368;
t291 = -rSges(6,3) * t395 - t445;
t290 = (-t267 * t366 + t268 * t367) * t265;
t289 = (-t267 * t364 + t268 * t365) * t265;
t282 = t294 - t357;
t273 = (t339 * t419 + (t64 * (-rSges(5,3) - pkin(7)) + t65 * t339) * t252) * t265;
t50 = t265 * t352 + t282;
t51 = t356 - t436;
t272 = ((t50 * (-t254 - t261) - t51 * t266) * t265 + t51 * t306 * qJD(4)) * t253;
t271 = ((t453 * t252 + ((t499 + t517) * t253 + t464 + t506 - t516) * t253) * qJD(4) + t482) * t335 + (-t494 * qJD(4) + t497 * t267 + t498 * t268) * t265 + (Icges(4,3) - t493) * t263 + (t458 + t460) * t433 + (t459 + t461) * t432 + (((t455 * t253 - t453 + t462) * t253 + (t455 * t252 + t463 - t507) * t252) * qJD(4) + t471 - t480) * t338 + (t467 + t468) * t337 + (t466 - t469 + t470) * t336;
t167 = t231 * t253;
t165 = t231 * t252;
t164 = t230 * t252;
t145 = t177 * t265 + t356;
t70 = qJD(4) * t311 + qJD(1);
t44 = qJD(1) + (-t252 * t481 + t381 * t253) * qJD(4);
t25 = t141 * t168 - t143 * t169 + qJDD(1) + (t252 * t96 + t253 * t94) * qJD(4);
t5 = qJDD(1) - t381 * t169 - t481 * t168 + (t252 * t423 + t253 * t424) * qJD(4);
t1 = [m(5) * t25 + m(6) * t5 + t358 * qJDD(1) + (-m(5) - m(6) - t358) * g(3); Icges(3,3) * qJDD(2) + t271 + (t474 * (t177 + t251) + t475 * (-t176 - t426) + (-t148 - t356 + t145) * t144) * m(4) + ((t180 ^ 2 + t181 ^ 2) * qJDD(2) + g(1) * t180 - g(2) * t181) * m(3) + (t50 * (t291 - t356) + t272 + t452 * (t105 - t426) + (-t357 + t50 - t282 + t456) * t51 + t451 * (t106 + t251)) * m(6) + (t64 * (t348 - t356) + t273 + (-t357 - t296 + t64 + t457) * t65 + t476 * (t110 + t251) + t477 * (t109 - t426)) * m(5); t271 + (t272 + (-t294 + t456) * t51 + (t291 - t436) * t50 + t451 * t106 + t452 * t105) * m(6) + (t273 + (t346 + t457) * t65 + (t348 - t439) * t64 + t476 * t110 + t477 * t109) * m(5) + (-t144 * t148 - t145 * t402 + (t144 * t265 + t474) * t177 + (t145 * t265 - t475) * t176) * m(4); t487 * t433 + t486 * t432 + (t467 * t265 + t458 * t263 + t463 * t169 + t462 * t168 + (t478 * t252 + t479 * t253) * qJD(4)) * t252 / 0.2e1 - (t466 * t265 + t459 * t263 + t465 * t169 + t464 * t168 + (t472 * t252 + t473 * t253) * qJD(4)) * t253 / 0.2e1 + (t460 * t252 - t461 * t253) * t263 / 0.2e1 - (((t364 + t366) * t268 + (t365 + t367) * t267) * t265 + (((-t374 - t376) * t253 + (t373 + t375) * t252) * t268 + ((t378 + t380) * t253 + (t377 + t379) * t252) * t267) * qJD(4)) * t265 / 0.2e1 + ((t460 * t265 + t469) * t253 + (t461 * t265 + t468) * t252) * t265 / 0.2e1 + t471 * t395 / 0.2e1 + t470 * t391 / 0.2e1 + ((-t363 * t400 + t399) * t252 + (t290 + (-t434 * t253 + (t401 + t293) * t252) * qJD(4)) * t253 + (-t363 * t397 + t396) * t252 + (t289 + (-t435 * t253 + (t398 + t292) * t252) * qJD(4)) * t253) * t338 + ((t462 * t265 + t479) * t253 + (t463 * t265 + t478) * t252) * t337 + ((t464 * t265 + t473) * t253 + (t465 * t265 + t472) * t252) * t336 + ((-t362 * t401 - t399) * t253 + (t290 + (t293 * t252 + (t400 - t434) * t253) * qJD(4)) * t252 + (-t362 * t398 - t396) * t253 + (t289 + (t292 * t252 + (t397 - t435) * t253) * qJD(4)) * t252) * t335 + (t25 * t311 + t70 * ((t94 + t125) * t253 + (-t143 * t265 + t96) * t252) + t320 * t197 + ((-t265 * t65 - t26) * t253 + (t265 * t64 - t27) * t252) * t231 - (t165 * t64 - t167 * t65) * t265 - (t70 * (-t165 * t252 - t167 * t253) + t320 * t233) * qJD(4) + g(1) * t167 + g(2) * t165 - g(3) * t233) * m(5) + ((t15 * t334 + t51 * t325 - t5 * t481 + t44 * t423 + (t50 * t230 - t381 * t44) * t265) * t252 + (t14 * t334 + t50 * t325 + t5 * t381 + t44 * t424 + (t334 * t51 - t44 * t481) * t265) * t253 - (t50 * t164 + t326 * t51) * t265 - ((t326 * t44 - t447 * t50) * t253 + (-t51 * t447 + (-pkin(4) * t393 - t164) * t44) * t252) * qJD(4) - g(3) * t447 - (g(1) * t253 + g(2) * t252) * t306) * m(6); (t252 * t452 - t451 * t253) * m(6);];
tau = t1;
