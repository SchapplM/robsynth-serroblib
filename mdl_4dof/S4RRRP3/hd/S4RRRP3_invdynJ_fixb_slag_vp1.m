% Calculate vector of inverse dynamics joint torques for
% S4RRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP3_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP3_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP3_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP3_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP3_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:14:01
% EndTime: 2019-12-31 17:14:14
% DurationCPUTime: 10.70s
% Computational Cost: add. (8624->477), mult. (9572->575), div. (0->0), fcn. (7484->6), ass. (0->273)
t266 = qJ(1) + qJ(2);
t255 = sin(t266);
t267 = sin(qJ(3));
t269 = cos(qJ(3));
t203 = Icges(4,5) * t269 - Icges(4,6) * t267;
t256 = cos(t266);
t299 = t203 * t256;
t128 = Icges(4,3) * t255 + t299;
t205 = Icges(5,4) * t269 + Icges(5,6) * t267;
t300 = t205 * t256;
t130 = Icges(5,2) * t255 + t300;
t521 = t128 + t130;
t258 = Icges(5,5) * t267;
t322 = Icges(5,1) * t269 + t258;
t133 = -Icges(5,4) * t256 + t255 * t322;
t398 = t255 * t267;
t231 = Icges(4,4) * t398;
t397 = t255 * t269;
t135 = Icges(4,1) * t397 - Icges(4,5) * t256 - t231;
t520 = -t133 - t135;
t302 = t322 * t256;
t134 = Icges(5,4) * t255 + t302;
t414 = Icges(4,4) * t267;
t211 = Icges(4,1) * t269 - t414;
t303 = t211 * t256;
t136 = Icges(4,5) * t255 + t303;
t513 = t134 + t136;
t206 = Icges(4,2) * t269 + t414;
t320 = Icges(5,3) * t269 - t258;
t519 = t206 + t320;
t413 = Icges(5,5) * t269;
t208 = Icges(5,1) * t267 - t413;
t259 = Icges(4,4) * t269;
t210 = Icges(4,1) * t267 + t259;
t518 = t208 + t210;
t201 = Icges(5,3) * t267 + t413;
t321 = -Icges(4,2) * t267 + t259;
t517 = t201 - t321;
t516 = t211 + t322;
t392 = t256 * t269;
t230 = Icges(5,5) * t392;
t393 = t256 * t267;
t126 = Icges(5,6) * t255 + Icges(5,3) * t393 + t230;
t515 = t126 * t393 + t521 * t255 + t513 * t392;
t129 = -Icges(5,2) * t256 + t205 * t255;
t121 = t255 * t129;
t125 = -Icges(5,6) * t256 + t201 * t255;
t127 = Icges(4,5) * t397 - Icges(4,6) * t398 - Icges(4,3) * t256;
t514 = -t125 * t393 - t255 * t127 + t520 * t392 - t121;
t202 = Icges(4,5) * t267 + Icges(4,6) * t269;
t204 = Icges(5,4) * t267 - Icges(5,6) * t269;
t512 = t202 + t204;
t265 = qJD(1) + qJD(2);
t511 = (-Icges(4,6) + Icges(5,6)) * t265 + t519 * qJD(3);
t510 = (Icges(5,4) + Icges(4,5)) * t265 - t518 * qJD(3);
t506 = -t519 * t267 + t518 * t269;
t131 = Icges(4,4) * t397 - Icges(4,2) * t398 - Icges(4,6) * t256;
t408 = t131 * t267;
t315 = -t135 * t269 + t408;
t395 = t256 * t129;
t318 = t125 * t267 + t133 * t269;
t453 = t255 * t318;
t49 = -t395 + t453;
t460 = -t256 * t127 - t255 * t315 + t49;
t458 = -t131 * t393 - t514;
t301 = t321 * t256;
t132 = Icges(4,6) * t255 + t301;
t457 = -t132 * t393 + t515;
t334 = -t126 * t398 + t130 * t256 - t134 * t397;
t110 = t136 * t397;
t337 = t256 * t128 - t110;
t52 = -t132 * t398 - t337;
t459 = t52 - t334;
t509 = t517 * qJD(3);
t508 = t516 * qJD(3);
t507 = t203 + t205;
t505 = -t518 * t267 - t519 * t269;
t484 = rSges(5,1) + pkin(3);
t399 = t255 * t265;
t504 = t511 * t256 - t399 * t517;
t394 = t256 * t265;
t503 = -t201 * t394 - t255 * t511 + t265 * t301;
t502 = t510 * t256 - t399 * t516;
t501 = (t302 + t303) * t265 + t510 * t255;
t490 = rSges(5,3) + qJ(4);
t500 = (t126 - t132) * t269 - t513 * t267;
t456 = (t125 - t131) * t269 + t520 * t267;
t449 = t490 * t397;
t401 = t204 * t256;
t404 = t202 * t256;
t477 = t255 * t506 - t401 - t404;
t402 = t204 * t255;
t405 = t202 * t255;
t476 = t256 * t506 + t402 + t405;
t499 = (-Icges(5,2) - Icges(4,3)) * t265 + t512 * qJD(3);
t407 = t132 * t267;
t498 = t126 * t267 + t269 * t513 - t407;
t497 = -t315 + t318;
t257 = t267 * qJ(4);
t332 = t269 * rSges(5,1) + t267 * rSges(5,3);
t496 = t269 * pkin(3) + t257 + t332;
t495 = t505 * qJD(3) + t265 * t512 + t509 * t267 + t508 * t269;
t494 = t457 * t255 - t256 * t458;
t493 = t459 * t255 - t460 * t256;
t492 = -qJD(3) * t507 + t506 * t265;
t491 = t476 * t265;
t248 = t256 * rSges(5,2);
t380 = t496 * t255 - t248;
t368 = t484 * t267 - t490 * t269;
t489 = t500 * qJD(3) + t521 * t265 + t504 * t267 + t502 * t269;
t488 = -t501 * t269 + t503 * t267 + (-t127 - t129) * t265 - t456 * qJD(3);
t487 = t477 * t265;
t486 = (-t300 - t299 + t497) * t265 + t499 * t255;
t485 = t499 * t256 + t498 * t265 + t399 * t507;
t483 = t493 * qJD(3) + t487;
t482 = t494 * qJD(3) + t491;
t481 = t497 * qJD(3) + t501 * t267 + t503 * t269;
t480 = t498 * qJD(3) + t502 * t267 - t504 * t269;
t479 = -t492 * t255 + t495 * t256;
t478 = t495 * t255 + t492 * t256;
t251 = t256 * pkin(6);
t180 = pkin(2) * t255 - t251;
t215 = pkin(6) * t394;
t475 = t265 * t180 + t215;
t362 = qJD(3) * t267;
t349 = t255 * t362;
t474 = t484 * t349;
t473 = t255 * rSges(5,2) + pkin(3) * t392;
t472 = t395 + t515;
t361 = qJD(3) * t269;
t348 = t255 * t361;
t391 = t265 * t267;
t471 = t256 * t391 + t348;
t268 = sin(qJ(1));
t422 = pkin(1) * qJD(1);
t357 = t268 * t422;
t178 = rSges(3,1) * t255 + rSges(3,2) * t256;
t406 = t178 * t265;
t142 = -t357 - t406;
t470 = t255 * t486 + t256 * t488;
t469 = -t255 * t485 + t256 * t489;
t181 = t256 * pkin(2) + t255 * pkin(6);
t147 = t181 * t265;
t363 = qJD(3) * t265;
t169 = -qJDD(3) * t256 + t255 * t363;
t264 = qJDD(1) + qJDD(2);
t359 = qJD(4) * t269;
t375 = -qJD(3) * t496 + t359;
t288 = qJDD(4) * t267 + (t359 + t375) * qJD(3);
t270 = cos(qJ(1));
t271 = qJD(1) ^ 2;
t298 = (-qJDD(1) * t268 - t270 * t271) * pkin(1);
t360 = qJD(4) * t267;
t345 = t255 * t360;
t354 = -t180 - t380;
t450 = -t345 + t474;
t424 = rSges(5,3) * t348 + t471 * qJ(4) - t450 + (t256 * t332 + t473) * t265;
t14 = t368 * t169 + t298 + t354 * t264 + (-t147 - t345 - t424) * t265 + t288 * t256;
t468 = t14 - g(1);
t168 = qJDD(3) * t255 + t256 * t363;
t224 = t256 * t360;
t263 = t270 * pkin(1);
t427 = pkin(1) * t268;
t335 = qJDD(1) * t263 - t271 * t427;
t306 = t265 * (-pkin(2) * t399 + t215) + t264 * t181 + t335;
t378 = rSges(5,1) * t392 + rSges(5,3) * t393 + t256 * t257 + t473;
t347 = t256 * t362;
t295 = -t265 * t397 - t347;
t356 = t255 * t391;
t346 = t256 * t361;
t446 = rSges(5,2) * t394 + t346 * t490 + t224;
t425 = t295 * t484 - t356 * t490 + t446;
t15 = t378 * t264 - t368 * t168 + (t224 + t425) * t265 + t288 * t255 + t306;
t467 = t15 - g(2);
t423 = rSges(4,1) * t269;
t222 = -rSges(4,2) * t267 + t423;
t191 = t222 * qJD(3);
t218 = rSges(4,1) * t267 + rSges(4,2) * t269;
t364 = qJD(3) * t256;
t366 = rSges(4,2) * t398 + t256 * rSges(4,3);
t138 = rSges(4,1) * t397 - t366;
t379 = -t138 - t180;
t352 = rSges(4,1) * t349 + rSges(4,2) * t471;
t447 = rSges(4,1) * t392 + t255 * rSges(4,3);
t91 = t265 * t447 - t352;
t25 = -t191 * t364 + t169 * t218 + (-t147 - t91) * t265 + t379 * t264 + t298;
t466 = t25 - g(1);
t140 = -rSges(4,2) * t393 + t447;
t365 = qJD(3) * t255;
t304 = rSges(4,3) * t394 + (-t346 + t356) * rSges(4,2);
t89 = rSges(4,1) * t295 + t304;
t26 = t140 * t264 - t168 * t218 - t191 * t365 + t265 * t89 + t306;
t465 = t26 - g(2);
t146 = rSges(3,1) * t394 - rSges(3,2) * t399;
t464 = -t146 * t265 - t178 * t264 - g(1) + t298;
t179 = t256 * rSges(3,1) - rSges(3,2) * t255;
t463 = t179 * t264 - t265 * t406 - g(2) + t335;
t462 = t255 * t488 - t256 * t486;
t461 = t255 * t489 + t256 * t485;
t350 = t218 * t364;
t297 = -t350 - t357;
t65 = t265 * t379 + t297;
t454 = t265 * t65;
t452 = t368 * t365;
t451 = t181 + t378;
t448 = t490 * t392;
t124 = t265 * t138;
t445 = -rSges(4,1) * t347 + t124 + t304 + t475;
t438 = t380 * t265 + t446 + t475;
t382 = -Icges(4,2) * t397 + t135 - t231;
t386 = t210 * t255 + t131;
t437 = -t267 * t382 - t269 * t386;
t384 = -t320 * t255 + t133;
t388 = -t208 * t255 + t125;
t436 = -t267 * t384 + t269 * t388;
t435 = t168 / 0.2e1;
t434 = t169 / 0.2e1;
t426 = g(2) * t255;
t307 = -t364 * t368 + t224;
t291 = t307 - t357;
t47 = t265 * t354 + t291;
t421 = t256 * t47;
t420 = t256 * t65;
t419 = t265 * t47;
t46 = -t359 + (t380 * t255 + t378 * t256) * qJD(3);
t418 = t46 * t267;
t403 = t203 * t265;
t400 = t205 * t265;
t387 = -Icges(5,1) * t393 + t126 + t230;
t385 = -t210 * t256 - t132;
t383 = -t320 * t256 + t134;
t381 = -t206 * t256 + t136;
t107 = t140 + t181;
t377 = -t398 * t484 + t449;
t376 = -t393 * t484 + t448;
t372 = -t320 + t322;
t371 = t201 - t208;
t370 = -t206 + t211;
t369 = t210 + t321;
t358 = t270 * t422;
t351 = t256 * t484;
t342 = -pkin(2) - t423;
t341 = -t365 / 0.2e1;
t340 = t365 / 0.2e1;
t339 = -t364 / 0.2e1;
t338 = t364 / 0.2e1;
t336 = -t127 + t407;
t223 = rSges(2,1) * t270 - rSges(2,2) * t268;
t219 = rSges(2,1) * t268 + rSges(2,2) * t270;
t296 = t345 + t358;
t48 = t265 * t451 + t296 - t452;
t330 = t255 * t48 + t421;
t177 = t218 * t365;
t66 = t107 * t265 - t177 + t358;
t325 = -t255 * t66 - t420;
t313 = t138 * t255 + t140 * t256;
t305 = t330 * t269;
t294 = -t267 * t383 + t269 * t387;
t293 = -t267 * t381 + t269 * t385;
t106 = t255 * t342 + t251 + t366;
t292 = -t267 * t490 - t269 * t484 - pkin(2);
t290 = (t267 * t371 + t269 * t372) * t265;
t289 = (-t267 * t369 + t269 * t370) * t265;
t96 = t255 * t292 + t248 + t251;
t274 = (t342 * t420 + (t65 * (-rSges(4,3) - pkin(6)) + t66 * t342) * t255) * t265;
t273 = (((t52 - t110 + (t128 + t408) * t256 + t514) * t256 + (t49 - t453 + t472) * t255) * qJD(3) + t491) * t338 + (qJD(3) * t506 + t267 * t508 - t269 * t509) * t265 + (Icges(3,3) - t505) * t264 + (t476 - t500) * t435 + (t477 - t456) * t434 + (t479 + t480) * t340 + (((t256 * t336 + t457 - t472) * t256 + (t255 * t336 - t121 + t334 + t337 + t458) * t255) * qJD(3) + t483 - t487) * t341 + (t478 + t481 + t482) * t339;
t272 = (t292 * t421 + (t47 * (-rSges(5,2) - pkin(6)) + t48 * (-pkin(2) - t496)) * t255) * t265 + (-t267 * t351 * t48 - t47 * t449) * qJD(3);
t166 = t218 * t256;
t162 = t218 * t255;
t143 = t179 * t265 + t358;
t69 = t313 * qJD(3);
t5 = -qJDD(4) * t269 - t378 * t169 + t380 * t168 + (t255 * t424 + t256 * t425 + t360) * qJD(3);
t1 = [Icges(2,3) * qJDD(1) + t273 + (t463 * (t179 + t263) + t464 * (-t178 - t427) + (-t146 - t358 + t143) * t142) * m(3) + ((t219 ^ 2 + t223 ^ 2) * qJDD(1) + g(1) * t219 - g(2) * t223) * m(2) + (t47 * (-t296 + t474) + t272 + t467 * (t263 + t451) + t468 * (t96 - t427) + (-t357 + t47 - t291 + t438) * t48) * m(5) + (t65 * (t352 - t358) + t274 + (-t357 - t297 + t65 + t445) * t66 + t465 * (t107 + t263) + t466 * (t106 - t427)) * m(4); t273 + (t272 + t468 * t96 + (-t307 + t438) * t48 + (t345 + t450 - t452) * t47 + (t419 + t467) * t451) * m(5) + (t274 + (t350 + t445) * t66 + (t352 - t177) * t65 + (t454 + t465) * t107 + t466 * t106) * m(4) + (-t142 * t146 - t143 * t406 + (t142 * t265 + t463) * t179 + (t143 * t265 - t464) * t178) * m(3); t494 * t435 + t493 * t434 + (t479 * t265 + t476 * t264 + t458 * t169 + t457 * t168 + (t469 * t255 + t470 * t256) * qJD(3)) * t255 / 0.2e1 - (t478 * t265 + t477 * t264 + t460 * t169 + t459 * t168 + (t461 * t255 + t462 * t256) * qJD(3)) * t256 / 0.2e1 + (-t255 * t500 + t256 * t456) * t264 / 0.2e1 - (((t369 - t371) * t269 + (t370 + t372) * t267) * t265 + (((-t382 - t384) * t256 + (t381 + t383) * t255) * t269 + ((t386 - t388) * t256 + (t385 + t387) * t255) * t267) * qJD(3)) * t265 / 0.2e1 + ((-t265 * t500 - t481) * t256 + (-t265 * t456 + t480) * t255) * t265 / 0.2e1 + t483 * t399 / 0.2e1 + t482 * t394 / 0.2e1 + ((-t365 * t401 + t400) * t255 + (t290 + (-t436 * t256 + (t402 + t294) * t255) * qJD(3)) * t256 + (-t365 * t404 + t403) * t255 + (t289 + (-t437 * t256 + (t405 + t293) * t255) * qJD(3)) * t256) * t341 + ((t265 * t457 + t470) * t256 + (t265 * t458 + t469) * t255) * t340 + ((t265 * t459 + t462) * t256 + (t265 * t460 + t461) * t255) * t339 + ((-t364 * t402 - t400) * t256 + (t290 + (t294 * t255 + (t401 - t436) * t256) * qJD(3)) * t255 + (-t364 * t405 - t403) * t256 + (t289 + (t293 * t255 + (t404 - t437) * t256) * qJD(3)) * t255) * t338 + (-(t305 + t418) * qJD(4) - (t376 * t48 - t377 * t47) * t265 - ((t376 * t46 - t47 * t496) * t256 + (t377 * t46 - t48 * t496) * t255) * qJD(3) - g(1) * t448 - g(2) * t449 - g(3) * t496 - (-g(1) * t351 - t426 * t484) * t267 + (-t14 * t368 + t47 * t375 + t5 * t378 + t46 * t425 + (-t368 * t48 + t380 * t46) * t265) * t256 + (-t15 * t368 + t48 * t375 + t5 * t380 + t46 * t424 + (t368 * t47 - t378 * t46) * t265) * t255) * m(5) + ((t138 * t168 - t140 * t169 + (t255 * t91 + t256 * t89) * qJD(3)) * t313 + t69 * ((t89 + t124) * t256 + (-t140 * t265 + t91) * t255) + t325 * t191 + ((-t265 * t66 - t25) * t256 + (-t26 + t454) * t255) * t218 - (t162 * t65 - t166 * t66) * t265 - (t69 * (-t162 * t255 - t166 * t256) + t325 * t222) * qJD(3) + g(1) * t166 + g(2) * t162 - g(3) * t222) * m(4); (-(-t255 * t47 + t256 * t48) * t391 - (t305 + (t255 ^ 2 + t256 ^ 2) * t418) * qJD(3) + (qJD(3) * t330 + g(3) - t5) * t269 + (qJD(3) * t46 + (t15 - t419) * t255 - t426 + (t265 * t48 + t468) * t256) * t267) * m(5);];
tau = t1;
