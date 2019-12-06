% Calculate vector of inverse dynamics joint torques for
% S5PRRRP2
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
% Datum: 2019-12-05 16:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRP2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP2_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRP2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:41:37
% EndTime: 2019-12-05 16:41:56
% DurationCPUTime: 11.63s
% Computational Cost: add. (12996->488), mult. (9751->579), div. (0->0), fcn. (7588->6), ass. (0->276)
t269 = pkin(8) + qJ(2);
t261 = qJ(3) + t269;
t256 = sin(t261);
t271 = sin(qJ(4));
t272 = cos(qJ(4));
t225 = Icges(5,5) * t272 - Icges(5,6) * t271;
t257 = cos(t261);
t301 = t225 * t257;
t129 = Icges(5,3) * t256 + t301;
t227 = Icges(6,4) * t272 + Icges(6,6) * t271;
t302 = t227 * t257;
t131 = Icges(6,2) * t256 + t302;
t524 = t129 + t131;
t263 = Icges(6,5) * t271;
t324 = Icges(6,1) * t272 + t263;
t134 = -Icges(6,4) * t257 + t256 * t324;
t399 = t256 * t271;
t214 = Icges(5,4) * t399;
t398 = t256 * t272;
t136 = Icges(5,1) * t398 - Icges(5,5) * t257 - t214;
t523 = -t134 - t136;
t304 = t324 * t257;
t135 = Icges(6,4) * t256 + t304;
t417 = Icges(5,4) * t271;
t233 = Icges(5,1) * t272 - t417;
t305 = t233 * t257;
t137 = Icges(5,5) * t256 + t305;
t516 = t135 + t137;
t228 = Icges(5,2) * t272 + t417;
t322 = Icges(6,3) * t272 - t263;
t522 = t228 + t322;
t416 = Icges(6,5) * t272;
t230 = Icges(6,1) * t271 - t416;
t264 = Icges(5,4) * t272;
t232 = Icges(5,1) * t271 + t264;
t521 = t230 + t232;
t223 = Icges(6,3) * t271 + t416;
t323 = -Icges(5,2) * t271 + t264;
t520 = t223 - t323;
t519 = t233 + t324;
t395 = t257 * t272;
t213 = Icges(6,5) * t395;
t396 = t257 * t271;
t127 = Icges(6,6) * t256 + Icges(6,3) * t396 + t213;
t518 = t127 * t396 + t524 * t256 + t516 * t395;
t130 = -Icges(6,2) * t257 + t227 * t256;
t120 = t256 * t130;
t126 = -Icges(6,6) * t257 + t223 * t256;
t128 = Icges(5,5) * t398 - Icges(5,6) * t399 - Icges(5,3) * t257;
t517 = -t126 * t396 - t256 * t128 + t523 * t395 - t120;
t224 = Icges(5,5) * t271 + Icges(5,6) * t272;
t226 = Icges(6,4) * t271 - Icges(6,6) * t272;
t515 = t224 + t226;
t270 = qJD(2) + qJD(3);
t514 = (-Icges(5,6) + Icges(6,6)) * t270 + t522 * qJD(4);
t513 = (Icges(6,4) + Icges(5,5)) * t270 - t521 * qJD(4);
t509 = -t522 * t271 + t521 * t272;
t132 = Icges(5,4) * t398 - Icges(5,2) * t399 - Icges(5,6) * t257;
t409 = t132 * t271;
t317 = -t136 * t272 + t409;
t410 = t130 * t257;
t320 = t126 * t271 + t134 * t272;
t456 = t256 * t320;
t50 = -t410 + t456;
t463 = -t128 * t257 - t256 * t317 + t50;
t461 = -t132 * t396 - t517;
t303 = t323 * t257;
t133 = Icges(5,6) * t256 + t303;
t460 = -t133 * t396 + t518;
t336 = -t127 * t399 + t131 * t257 - t135 * t398;
t111 = t137 * t398;
t339 = t129 * t257 - t111;
t53 = -t133 * t399 - t339;
t462 = t53 - t336;
t512 = t520 * qJD(4);
t511 = t519 * qJD(4);
t510 = t225 + t227;
t508 = -t521 * t271 - t522 * t272;
t487 = rSges(6,1) + pkin(4);
t400 = t256 * t270;
t507 = t514 * t257 - t400 * t520;
t397 = t257 * t270;
t506 = -t223 * t397 - t514 * t256 + t270 * t303;
t505 = t513 * t257 - t400 * t519;
t504 = (t304 + t305) * t270 + t513 * t256;
t494 = rSges(6,3) + qJ(5);
t503 = (t127 - t133) * t272 - t516 * t271;
t459 = (t126 - t132) * t272 + t523 * t271;
t452 = t494 * t398;
t402 = t226 * t257;
t405 = t224 * t257;
t480 = t509 * t256 - t402 - t405;
t403 = t226 * t256;
t406 = t224 * t256;
t479 = t509 * t257 + t403 + t406;
t502 = (-Icges(6,2) - Icges(5,3)) * t270 + t515 * qJD(4);
t408 = t133 * t271;
t501 = t127 * t271 + t516 * t272 - t408;
t500 = -t317 + t320;
t262 = t271 * qJ(5);
t334 = t272 * rSges(6,1) + t271 * rSges(6,3);
t499 = t272 * pkin(4) + t262 + t334;
t498 = t508 * qJD(4) + t515 * t270 + t512 * t271 + t511 * t272;
t497 = t460 * t256 - t461 * t257;
t496 = t462 * t256 - t463 * t257;
t495 = -t510 * qJD(4) + t509 * t270;
t493 = t479 * t270;
t250 = t257 * rSges(6,2);
t383 = t499 * t256 - t250;
t370 = t487 * t271 - t494 * t272;
t492 = t503 * qJD(4) + t524 * t270 + t507 * t271 + t505 * t272;
t491 = -t504 * t272 + t506 * t271 + (-t128 - t130) * t270 - t459 * qJD(4);
t490 = t480 * t270;
t489 = (-t302 - t301 + t500) * t270 + t502 * t256;
t488 = t502 * t257 + t501 * t270 + t510 * t400;
t486 = qJD(4) * t496 + t490;
t485 = qJD(4) * t497 + t493;
t484 = t500 * qJD(4) + t504 * t271 + t506 * t272;
t483 = t501 * qJD(4) + t505 * t271 - t507 * t272;
t482 = -t256 * t495 + t257 * t498;
t481 = t256 * t498 + t257 * t495;
t253 = t257 * pkin(7);
t183 = pkin(3) * t256 - t253;
t206 = pkin(7) * t397;
t478 = t270 * t183 + t206;
t365 = qJD(4) * t271;
t351 = t256 * t365;
t477 = t487 * t351;
t476 = t256 * rSges(6,2) + pkin(4) * t395;
t475 = t410 + t518;
t364 = qJD(4) * t272;
t350 = t256 * t364;
t394 = t270 * t271;
t474 = t257 * t394 + t350;
t259 = sin(t269);
t425 = pkin(2) * qJD(2);
t360 = t259 * t425;
t181 = rSges(4,1) * t256 + rSges(4,2) * t257;
t407 = t181 * t270;
t142 = -t360 - t407;
t473 = t256 * t489 + t257 * t491;
t472 = -t256 * t488 + t257 * t492;
t184 = t257 * pkin(3) + t256 * pkin(7);
t148 = t184 * t270;
t366 = qJD(4) * t270;
t172 = -qJDD(4) * t257 + t256 * t366;
t268 = qJDD(2) + qJDD(3);
t362 = qJD(5) * t272;
t378 = -qJD(4) * t499 + t362;
t290 = qJDD(5) * t271 + (t362 + t378) * qJD(4);
t260 = cos(t269);
t273 = qJD(2) ^ 2;
t300 = (-qJDD(2) * t259 - t260 * t273) * pkin(2);
t363 = qJD(5) * t271;
t347 = t256 * t363;
t356 = -t183 - t383;
t453 = -t347 + t477;
t427 = rSges(6,3) * t350 + t474 * qJ(5) - t453 + (t257 * t334 + t476) * t270;
t14 = t370 * t172 + t300 + t356 * t268 + (-t148 - t347 - t427) * t270 + t290 * t257;
t471 = t14 - g(1);
t171 = qJDD(4) * t256 + t257 * t366;
t207 = t257 * t363;
t255 = pkin(2) * t260;
t430 = pkin(2) * t259;
t337 = qJDD(2) * t255 - t273 * t430;
t308 = t270 * (-pkin(3) * t400 + t206) + t268 * t184 + t337;
t381 = rSges(6,1) * t395 + rSges(6,3) * t396 + t257 * t262 + t476;
t349 = t257 * t365;
t297 = -t270 * t398 - t349;
t358 = t256 * t394;
t348 = t257 * t364;
t449 = rSges(6,2) * t397 + t348 * t494 + t207;
t428 = t297 * t487 - t358 * t494 + t449;
t15 = t381 * t268 - t370 * t171 + (t207 + t428) * t270 + t290 * t256 + t308;
t470 = t15 - g(2);
t426 = rSges(5,1) * t272;
t239 = -rSges(5,2) * t271 + t426;
t196 = t239 * qJD(4);
t236 = rSges(5,1) * t271 + rSges(5,2) * t272;
t367 = qJD(4) * t257;
t375 = rSges(5,2) * t399 + t257 * rSges(5,3);
t139 = rSges(5,1) * t398 - t375;
t382 = -t139 - t183;
t354 = rSges(5,1) * t351 + rSges(5,2) * t474;
t450 = rSges(5,1) * t395 + t256 * rSges(5,3);
t96 = t270 * t450 - t354;
t26 = -t196 * t367 + t172 * t236 + (-t148 - t96) * t270 + t382 * t268 + t300;
t469 = t26 - g(1);
t141 = -rSges(5,2) * t396 + t450;
t368 = qJD(4) * t256;
t306 = rSges(5,3) * t397 + (-t348 + t358) * rSges(5,2);
t94 = rSges(5,1) * t297 + t306;
t27 = t141 * t268 - t171 * t236 - t196 * t368 + t270 * t94 + t308;
t468 = t27 - g(2);
t147 = rSges(4,1) * t397 - rSges(4,2) * t400;
t467 = -t147 * t270 - t181 * t268 - g(1) + t300;
t182 = t257 * rSges(4,1) - rSges(4,2) * t256;
t466 = t182 * t268 - t270 * t407 - g(2) + t337;
t465 = t256 * t491 - t257 * t489;
t464 = t256 * t492 + t257 * t488;
t352 = t236 * t367;
t299 = -t352 - t360;
t62 = t270 * t382 + t299;
t457 = t270 * t62;
t455 = t370 * t368;
t454 = t184 + t381;
t451 = t494 * t395;
t123 = t270 * t139;
t448 = -rSges(5,1) * t349 + t123 + t306 + t478;
t441 = t383 * t270 + t449 + t478;
t385 = -Icges(5,2) * t398 + t136 - t214;
t389 = t232 * t256 + t132;
t440 = -t271 * t385 - t272 * t389;
t387 = -t322 * t256 + t134;
t391 = -t230 * t256 + t126;
t439 = -t271 * t387 + t272 * t391;
t438 = t171 / 0.2e1;
t437 = t172 / 0.2e1;
t429 = g(2) * t256;
t309 = -t367 * t370 + t207;
t293 = t309 - t360;
t48 = t270 * t356 + t293;
t424 = t257 * t48;
t423 = t257 * t62;
t422 = t270 * t48;
t47 = -t362 + qJD(1) + (t383 * t256 + t381 * t257) * qJD(4);
t421 = t47 * t271;
t404 = t225 * t270;
t401 = t227 * t270;
t390 = -Icges(6,1) * t396 + t127 + t213;
t388 = -t232 * t257 - t133;
t386 = -t322 * t257 + t135;
t384 = -t228 * t257 + t137;
t108 = t141 + t184;
t380 = -t399 * t487 + t452;
t379 = -t396 * t487 + t451;
t374 = -t322 + t324;
t373 = t223 - t230;
t372 = -t228 + t233;
t371 = t232 + t323;
t361 = m(2) + m(3) + m(4);
t359 = t260 * t425;
t353 = t257 * t487;
t344 = -pkin(3) - t426;
t343 = -t368 / 0.2e1;
t342 = t368 / 0.2e1;
t341 = -t367 / 0.2e1;
t340 = t367 / 0.2e1;
t338 = -t128 + t408;
t186 = rSges(3,1) * t260 - rSges(3,2) * t259;
t185 = rSges(3,1) * t259 + rSges(3,2) * t260;
t298 = t347 + t359;
t49 = t270 * t454 + t298 - t455;
t332 = t256 * t49 + t424;
t177 = t236 * t368;
t63 = t108 * t270 - t177 + t359;
t327 = -t256 * t63 - t423;
t315 = t139 * t256 + t141 * t257;
t307 = t332 * t272;
t296 = -t271 * t386 + t272 * t390;
t295 = -t271 * t384 + t272 * t388;
t107 = t256 * t344 + t253 + t375;
t294 = -t271 * t494 - t272 * t487 - pkin(3);
t292 = (t271 * t373 + t272 * t374) * t270;
t291 = (-t271 * t371 + t272 * t372) * t270;
t75 = t256 * t294 + t250 + t253;
t276 = (t344 * t423 + (t62 * (-rSges(5,3) - pkin(7)) + t63 * t344) * t256) * t270;
t275 = (((t53 - t111 + (t129 + t409) * t257 + t517) * t257 + (t50 - t456 + t475) * t256) * qJD(4) + t493) * t340 + (t509 * qJD(4) + t511 * t271 - t512 * t272) * t270 + (Icges(4,3) - t508) * t268 + (t479 - t503) * t438 + (t480 - t459) * t437 + (t482 + t483) * t342 + (((t257 * t338 + t460 - t475) * t257 + (t256 * t338 - t120 + t336 + t339 + t461) * t256) * qJD(4) + t486 - t490) * t343 + (t481 + t484 + t485) * t341;
t274 = (t294 * t424 + (t48 * (-rSges(6,2) - pkin(7)) + t49 * (-pkin(3) - t499)) * t256) * t270 + (-t271 * t353 * t49 - t48 * t452) * qJD(4);
t169 = t236 * t257;
t165 = t236 * t256;
t143 = t182 * t270 + t359;
t68 = qJD(4) * t315 + qJD(1);
t25 = t139 * t171 - t141 * t172 + qJDD(1) + (t256 * t96 + t257 * t94) * qJD(4);
t5 = -qJDD(5) * t272 + qJDD(1) - t381 * t172 + t383 * t171 + (t256 * t427 + t257 * t428 + t363) * qJD(4);
t1 = [m(5) * t25 + m(6) * t5 + t361 * qJDD(1) + (-m(5) - m(6) - t361) * g(3); Icges(3,3) * qJDD(2) + t275 + (t466 * (t182 + t255) + t467 * (-t181 - t430) + (-t147 - t359 + t143) * t142) * m(4) + ((t185 ^ 2 + t186 ^ 2) * qJDD(2) + g(1) * t185 - g(2) * t186) * m(3) + (t48 * (-t298 + t477) + t274 + t470 * (t255 + t454) + t471 * (t75 - t430) + (-t360 + t48 - t293 + t441) * t49) * m(6) + (t62 * (t354 - t359) + t276 + (-t360 - t299 + t62 + t448) * t63 + t468 * (t108 + t255) + t469 * (t107 - t430)) * m(5); t275 + (t274 + t471 * t75 + (-t309 + t441) * t49 + (t347 + t453 - t455) * t48 + (t422 + t470) * t454) * m(6) + (t276 + (t352 + t448) * t63 + (t354 - t177) * t62 + (t457 + t468) * t108 + t469 * t107) * m(5) + (-t142 * t147 - t143 * t407 + (t142 * t270 + t466) * t182 + (t143 * t270 - t467) * t181) * m(4); t497 * t438 + t496 * t437 + (t482 * t270 + t479 * t268 + t461 * t172 + t460 * t171 + (t472 * t256 + t473 * t257) * qJD(4)) * t256 / 0.2e1 - (t481 * t270 + t480 * t268 + t463 * t172 + t462 * t171 + (t464 * t256 + t465 * t257) * qJD(4)) * t257 / 0.2e1 + (-t256 * t503 + t257 * t459) * t268 / 0.2e1 - (((t371 - t373) * t272 + (t372 + t374) * t271) * t270 + (((-t385 - t387) * t257 + (t384 + t386) * t256) * t272 + ((t389 - t391) * t257 + (t388 + t390) * t256) * t271) * qJD(4)) * t270 / 0.2e1 + ((-t270 * t503 - t484) * t257 + (-t270 * t459 + t483) * t256) * t270 / 0.2e1 + t486 * t400 / 0.2e1 + t485 * t397 / 0.2e1 + ((-t368 * t402 + t401) * t256 + (t292 + (-t439 * t257 + (t403 + t296) * t256) * qJD(4)) * t257 + (-t368 * t405 + t404) * t256 + (t291 + (-t440 * t257 + (t406 + t295) * t256) * qJD(4)) * t257) * t343 + ((t270 * t460 + t473) * t257 + (t270 * t461 + t472) * t256) * t342 + ((t270 * t462 + t465) * t257 + (t270 * t463 + t464) * t256) * t341 + ((-t367 * t403 - t401) * t257 + (t292 + (t296 * t256 + (t402 - t439) * t257) * qJD(4)) * t256 + (-t367 * t406 - t404) * t257 + (t291 + (t295 * t256 + (t405 - t440) * t257) * qJD(4)) * t256) * t340 + ((-t14 * t370 + t48 * t378 + t5 * t381 + t47 * t428 + (-t370 * t49 + t383 * t47) * t270) * t257 + (-t15 * t370 + t49 * t378 + t5 * t383 + t47 * t427 + (t370 * t48 - t381 * t47) * t270) * t256 - g(1) * t451 - g(2) * t452 - g(3) * t499 - (-g(1) * t353 - t429 * t487) * t271 - (t307 + t421) * qJD(5) - (t379 * t49 - t380 * t48) * t270 - ((t379 * t47 - t48 * t499) * t257 + (t380 * t47 - t49 * t499) * t256) * qJD(4)) * m(6) + (-(t165 * t62 - t169 * t63) * t270 - (t68 * (-t165 * t256 - t169 * t257) + t327 * t239) * qJD(4) + t25 * t315 + t68 * ((t94 + t123) * t257 + (-t141 * t270 + t96) * t256) + t327 * t196 + ((-t270 * t63 - t26) * t257 + (-t27 + t457) * t256) * t236 + g(1) * t169 + g(2) * t165 - g(3) * t239) * m(5); (-(-t256 * t48 + t257 * t49) * t394 - (t307 + (t256 ^ 2 + t257 ^ 2) * t421) * qJD(4) + (qJD(4) * t332 + g(3) - t5) * t272 + (qJD(4) * t47 + (t15 - t422) * t256 - t429 + (t270 * t49 + t471) * t257) * t271) * m(6);];
tau = t1;
