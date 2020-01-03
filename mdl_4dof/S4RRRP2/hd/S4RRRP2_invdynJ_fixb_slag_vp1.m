% Calculate vector of inverse dynamics joint torques for
% S4RRRP2
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
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRP2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP2_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP2_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP2_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:53
% EndTime: 2019-12-31 17:13:08
% DurationCPUTime: 12.79s
% Computational Cost: add. (8632->446), mult. (9156->549), div. (0->0), fcn. (7147->6), ass. (0->263)
t262 = sin(qJ(3));
t264 = cos(qJ(3));
t200 = Icges(5,5) * t264 - Icges(5,6) * t262;
t202 = Icges(4,5) * t264 - Icges(4,6) * t262;
t493 = t200 + t202;
t521 = Icges(4,5) + Icges(5,5);
t520 = Icges(4,6) + Icges(5,6);
t519 = Icges(4,3) + Icges(5,3);
t260 = qJ(1) + qJ(2);
t252 = cos(t260);
t251 = sin(t260);
t391 = t251 * t264;
t392 = t251 * t262;
t130 = Icges(5,4) * t391 - Icges(5,2) * t392 - Icges(5,6) * t252;
t132 = Icges(4,4) * t391 - Icges(4,2) * t392 - Icges(4,6) * t252;
t512 = t130 + t132;
t225 = Icges(5,4) * t392;
t134 = Icges(5,1) * t391 - Icges(5,5) * t252 - t225;
t226 = Icges(4,4) * t392;
t136 = Icges(4,1) * t391 - Icges(4,5) * t252 - t226;
t517 = t134 + t136;
t408 = Icges(5,4) * t262;
t208 = Icges(5,1) * t264 - t408;
t300 = t208 * t252;
t135 = Icges(5,5) * t251 + t300;
t409 = Icges(4,4) * t262;
t210 = Icges(4,1) * t264 - t409;
t301 = t210 * t252;
t137 = Icges(4,5) * t251 + t301;
t510 = t135 + t137;
t203 = Icges(5,2) * t264 + t408;
t205 = Icges(4,2) * t264 + t409;
t516 = t203 + t205;
t253 = Icges(5,4) * t264;
t207 = Icges(5,1) * t262 + t253;
t254 = Icges(4,4) * t264;
t209 = Icges(4,1) * t262 + t254;
t515 = -t207 - t209;
t518 = t493 * t252;
t497 = t519 * t252 - t391 * t521 + t520 * t392;
t496 = t519 * t251 + t518;
t314 = -Icges(5,2) * t262 + t253;
t298 = t314 * t252;
t131 = Icges(5,6) * t251 + t298;
t315 = -Icges(4,2) * t262 + t254;
t299 = t315 * t252;
t133 = Icges(4,6) * t251 + t299;
t511 = t131 + t133;
t514 = t512 * t262;
t513 = t510 * t391;
t199 = Icges(5,5) * t262 + Icges(5,6) * t264;
t201 = Icges(4,5) * t262 + Icges(4,6) * t264;
t509 = t199 + t201;
t508 = t208 + t210;
t259 = qJD(1) + qJD(2);
t507 = -t516 * qJD(3) + t259 * t520;
t506 = t515 * qJD(3) + t259 * t521;
t488 = -t517 * t264 + t514;
t505 = t314 + t315;
t491 = t516 * t262 + t515 * t264;
t504 = -t496 * t252 + t513;
t386 = t252 * t264;
t503 = t497 * t251 - t517 * t386;
t448 = t496 * t251 + t510 * t386;
t502 = t511 * t262;
t462 = -t488 * t251 + t497 * t252;
t461 = -t511 * t392 + t504;
t387 = t252 * t262;
t460 = -t512 * t387 - t503;
t459 = -t511 * t387 + t448;
t458 = t517 * t262 + t512 * t264;
t457 = t510 * t262 + t511 * t264;
t394 = t251 * t259;
t501 = t507 * t252 - t505 * t394;
t500 = (t298 + t299) * t259 + t507 * t251;
t499 = t506 * t252 - t508 * t394;
t498 = (-t300 - t301) * t259 - t506 * t251;
t495 = t505 * qJD(3);
t494 = t508 * qJD(3);
t492 = -t509 * qJD(3) + t259 * t519;
t490 = t515 * t262 - t516 * t264;
t489 = -t510 * t264 + t502;
t396 = t201 * t252;
t399 = t199 * t252;
t456 = -t491 * t251 - t396 - t399;
t397 = t201 * t251;
t400 = t199 * t251;
t455 = -t491 * t252 + t397 + t400;
t261 = -qJ(4) - pkin(6);
t231 = t252 * t261;
t256 = t264 * pkin(3);
t248 = t256 + pkin(2);
t103 = -rSges(5,1) * t391 + rSges(5,2) * t392 + t252 * rSges(5,3) - t251 * t248 - t231;
t487 = t490 * qJD(3) + t509 * t259 - t495 * t262 + t494 * t264;
t486 = -t457 * qJD(3) + t496 * t259 - t501 * t262 + t499 * t264;
t485 = t458 * qJD(3) + t497 * t259 + t500 * t262 + t498 * t264;
t484 = t459 * t251 - t460 * t252;
t483 = t461 * t251 - t462 * t252;
t482 = t493 * qJD(3) + t491 * t259;
t481 = (t488 + t518) * t259 + t492 * t251;
t480 = t492 * t252 + t489 * t259 - t493 * t394;
t479 = t455 * t259;
t246 = t252 * pkin(6);
t174 = pkin(2) * t251 - t246;
t478 = t174 + t103;
t477 = t456 * t259;
t476 = -t251 * t481 + t252 * t485;
t475 = t251 * t480 + t252 * t486;
t474 = t251 * t485 + t252 * t481;
t245 = t251 * pkin(6);
t175 = t252 * pkin(2) + t245;
t147 = t175 * t259;
t358 = qJD(3) * t259;
t165 = -qJDD(3) * t252 + t251 * t358;
t255 = t264 * rSges(5,1);
t218 = -rSges(5,2) * t262 + t255;
t186 = t218 * qJD(3);
t414 = t264 * rSges(5,2);
t215 = rSges(5,1) * t262 + t414;
t233 = qJD(4) * t252;
t258 = qJDD(1) + qJDD(2);
t263 = sin(qJ(1));
t265 = cos(qJ(1));
t267 = qJD(1) ^ 2;
t295 = (-qJDD(1) * t263 - t265 * t267) * pkin(1);
t350 = -t174 + t478;
t359 = qJD(3) * t252;
t384 = t264 * qJD(3) ^ 2;
t421 = pkin(2) - t248;
t442 = rSges(5,1) * t386 + t251 * rSges(5,3);
t357 = qJD(3) * t262;
t343 = t251 * t357;
t366 = pkin(3) * t343 + t233;
t393 = t251 * t261;
t356 = qJD(3) * t264;
t342 = t251 * t356;
t385 = t259 * t262;
t450 = t252 * t385 + t342;
t449 = -rSges(5,1) * t343 - t450 * rSges(5,2) - t259 * t393 - t366;
t419 = t449 + (-t252 * t421 - t245 + t442) * t259;
t13 = -t186 * t359 + qJDD(4) * t251 + t165 * t215 + (t165 * t262 - t252 * t384) * pkin(3) + t295 + t350 * t258 + (-t147 + t233 - t419) * t259;
t446 = t13 - g(1);
t164 = qJDD(3) * t251 + t252 * t358;
t388 = t252 * t259;
t214 = pkin(6) * t388;
t257 = t265 * pkin(1);
t422 = pkin(1) * t263;
t325 = qJDD(1) * t257 - t267 * t422;
t303 = t259 * (-pkin(2) * t394 + t214) + t258 * t175 + t325;
t332 = -pkin(3) * t262 - t215;
t104 = -rSges(5,2) * t387 + t252 * t248 - t393 + t442;
t378 = -t175 + t104;
t341 = t252 * t357;
t352 = t259 * t391;
t293 = -t341 - t352;
t340 = t252 * t356;
t232 = qJD(4) * t251;
t353 = t251 * t385;
t454 = rSges(5,2) * t353 + rSges(5,3) * t388 + t232;
t420 = -pkin(3) * t341 - t214 + (t251 * t421 - t231) * t259 + rSges(5,1) * t293 - rSges(5,2) * t340 + t454;
t14 = -qJDD(4) * t252 + t420 * t259 + t378 * t258 + t332 * t164 + (-pkin(3) * t384 - qJD(3) * t186 + qJD(4) * t259) * t251 + t303;
t445 = t14 - g(2);
t418 = rSges(4,1) * t264;
t219 = -rSges(4,2) * t262 + t418;
t187 = t219 * qJD(3);
t216 = rSges(4,1) * t262 + rSges(4,2) * t264;
t361 = rSges(4,2) * t392 + t252 * rSges(4,3);
t139 = rSges(4,1) * t391 - t361;
t369 = -t139 - t174;
t346 = rSges(4,1) * t343 + t450 * rSges(4,2);
t441 = rSges(4,1) * t386 + t251 * rSges(4,3);
t92 = t259 * t441 - t346;
t24 = -t187 * t359 + t165 * t216 + (-t147 - t92) * t259 + t369 * t258 + t295;
t473 = t24 - g(1);
t141 = -rSges(4,2) * t387 + t441;
t360 = qJD(3) * t251;
t302 = rSges(4,3) * t388 + (-t340 + t353) * rSges(4,2);
t90 = rSges(4,1) * t293 + t302;
t25 = t141 * t258 - t164 * t216 - t187 * t360 + t259 * t90 + t303;
t472 = t25 - g(2);
t146 = rSges(3,1) * t388 - rSges(3,2) * t394;
t172 = rSges(3,1) * t251 + rSges(3,2) * t252;
t471 = -t146 * t259 - t172 * t258 - g(1) + t295;
t173 = t252 * rSges(3,1) - rSges(3,2) * t251;
t401 = t172 * t259;
t470 = t173 * t258 - t259 * t401 - g(2) + t325;
t469 = t251 * t486 - t252 * t480;
t468 = qJD(3) * t483 + t477;
t467 = qJD(3) * t484 + t479;
t466 = t488 * qJD(3) + t498 * t262 - t500 * t264;
t465 = -t489 * qJD(3) + t499 * t262 + t501 * t264;
t464 = t251 * t482 + t252 * t487;
t463 = t251 * t487 - t252 * t482;
t125 = t259 * t139;
t169 = t259 * t174;
t453 = -rSges(4,1) * t341 + t125 + t169 + t214 + t302;
t452 = -rSges(5,1) * t352 - t248 * t394 - t478 * t259 + t169 + t454;
t451 = t497 + t502;
t416 = pkin(1) * qJD(1);
t354 = t263 * t416;
t143 = -t354 - t401;
t447 = t251 * t419 + t252 * t420;
t443 = t218 + t256;
t108 = t141 + t175;
t435 = -t108 * t259 + t216 * t360;
t432 = t215 * t360 - t259 * (t175 + t378) + t366;
t371 = -Icges(4,2) * t391 + t136 - t226;
t375 = t209 * t251 + t132;
t431 = -t262 * t371 - t264 * t375;
t373 = -Icges(5,2) * t391 + t134 - t225;
t377 = t207 * t251 + t130;
t430 = -t262 * t373 - t264 * t377;
t429 = t164 / 0.2e1;
t428 = t165 / 0.2e1;
t344 = t216 * t359;
t294 = -t344 - t354;
t66 = t259 * t369 + t294;
t415 = t252 * t66;
t398 = t200 * t259;
t395 = t202 * t259;
t376 = -t207 * t252 - t131;
t374 = -t209 * t252 - t133;
t372 = -t203 * t252 + t135;
t370 = -t205 * t252 + t137;
t365 = -t203 + t208;
t364 = t207 + t314;
t363 = -t205 + t210;
t362 = t209 + t315;
t355 = t265 * t416;
t337 = -pkin(2) - t418;
t336 = -t360 / 0.2e1;
t335 = t360 / 0.2e1;
t334 = -t359 / 0.2e1;
t333 = t359 / 0.2e1;
t324 = -pkin(3) * t387 - t215 * t252;
t323 = t332 * t252;
t220 = rSges(2,1) * t265 - rSges(2,2) * t263;
t217 = rSges(2,1) * t263 + rSges(2,2) * t265;
t67 = t355 - t435;
t318 = -t251 * t67 - t415;
t309 = t139 * t251 + t141 * t252;
t304 = -t414 + (-rSges(5,1) - pkin(3)) * t262;
t292 = qJD(3) * t323 + t232;
t291 = -t251 * t478 + t252 * t378;
t290 = -t262 * t372 + t264 * t376;
t289 = -t262 * t370 + t264 * t374;
t107 = t251 * t337 + t246 + t361;
t288 = -rSges(5,3) * t394 - t449;
t287 = (-t262 * t364 + t264 * t365) * t259;
t286 = (-t262 * t362 + t264 * t363) * t259;
t279 = t292 - t354;
t270 = (t337 * t415 + (t66 * (-rSges(4,3) - pkin(6)) + t67 * t337) * t251) * t259;
t48 = t259 * t350 + t279;
t49 = t355 - t432;
t269 = ((t48 * (-t248 - t255) - t49 * t261) * t259 + t49 * t304 * qJD(3)) * t252;
t268 = ((t448 * t251 + ((t496 + t514) * t252 + t461 + t503 - t513) * t252) * qJD(3) + t479) * t333 + (-t491 * qJD(3) + t494 * t262 + t495 * t264) * t259 + (Icges(3,3) - t490) * t258 + (t455 + t457) * t429 + (t456 + t458) * t428 + (((t451 * t252 - t448 + t459) * t252 + (t451 * t251 + t460 - t504) * t251) * qJD(3) + t468 - t477) * t336 + (t464 + t465) * t335 + (t463 - t466 + t467) * t334;
t163 = t216 * t252;
t161 = t216 * t251;
t160 = t215 * t251;
t144 = t173 * t259 + t355;
t70 = t309 * qJD(3);
t42 = t291 * qJD(3);
t1 = [Icges(2,3) * qJDD(1) + t268 + (t470 * (t173 + t257) + t471 * (-t172 - t422) + (-t146 - t355 + t144) * t143) * m(3) + ((t217 ^ 2 + t220 ^ 2) * qJDD(1) + g(1) * t217 - g(2) * t220) * m(2) + (t48 * (t288 - t355) + t269 + (t48 - t279 - t354 + t452) * t49 + t445 * (t104 + t257) + t446 * (t103 - t422)) * m(5) + (t66 * (t346 - t355) + t270 + (-t294 + t66 - t354 + t453) * t67 + t472 * (t108 + t257) + t473 * (t107 - t422)) * m(4); t268 + (t269 + (-t292 + t452) * t49 + (t288 - t432) * t48 + t445 * t104 + t446 * t103) * m(5) + (t270 + (t344 + t453) * t67 + (-t435 + t346) * t66 + t472 * t108 + t473 * t107) * m(4) + (-t143 * t146 - t144 * t401 + (t143 * t259 + t470) * t173 + (t144 * t259 - t471) * t172) * m(3); t484 * t429 + t483 * t428 + (t464 * t259 + t455 * t258 + t460 * t165 + t459 * t164 + (t475 * t251 + t476 * t252) * qJD(3)) * t251 / 0.2e1 - (t463 * t259 + t456 * t258 + t462 * t165 + t461 * t164 + (t469 * t251 + t474 * t252) * qJD(3)) * t252 / 0.2e1 + (t457 * t251 - t458 * t252) * t258 / 0.2e1 - (((t362 + t364) * t264 + (t363 + t365) * t262) * t259 + (((-t371 - t373) * t252 + (t370 + t372) * t251) * t264 + ((t375 + t377) * t252 + (t374 + t376) * t251) * t262) * qJD(3)) * t259 / 0.2e1 + ((t457 * t259 + t466) * t252 + (t458 * t259 + t465) * t251) * t259 / 0.2e1 + t468 * t394 / 0.2e1 + t467 * t388 / 0.2e1 + ((-t360 * t396 + t395) * t251 + (t286 + (-t431 * t252 + (t397 + t289) * t251) * qJD(3)) * t252 + (-t360 * t399 + t398) * t251 + (t287 + (-t430 * t252 + (t400 + t290) * t251) * qJD(3)) * t252) * t336 + ((t459 * t259 + t476) * t252 + (t460 * t259 + t475) * t251) * t335 + ((t461 * t259 + t474) * t252 + (t462 * t259 + t469) * t251) * t334 + ((-t359 * t397 - t395) * t252 + (t286 + (t289 * t251 + (t396 - t431) * t252) * qJD(3)) * t251 + (-t359 * t400 - t398) * t252 + (t287 + (t290 * t251 + (t399 - t430) * t252) * qJD(3)) * t251) * t333 + (-g(3) * t443 - (g(1) * t252 + g(2) * t251) * t304 - ((t324 * t42 - t443 * t48) * t252 + (-t49 * t443 + (-pkin(3) * t392 - t160) * t42) * t251) * qJD(3) + t13 * t323 + t48 * (-pkin(3) * t340 - t186 * t252) + t14 * t332 * t251 + t49 * (-pkin(3) * t342 - t186 * t251) + (t447 * qJD(3) - t164 * t478 - t378 * t165) * t291 + t42 * t447 + (-t48 * t160 - t324 * t49 + (t48 * t215 - t378 * t42) * t251 + (t332 * t49 - t42 * t478) * t252) * t259) * m(5) + ((t139 * t164 - t141 * t165 + (t251 * t92 + t252 * t90) * qJD(3)) * t309 + t70 * ((t90 + t125) * t252 + (-t141 * t259 + t92) * t251) + t318 * t187 + ((-t259 * t67 - t24) * t252 + (t259 * t66 - t25) * t251) * t216 - (t161 * t66 - t163 * t67) * t259 - (t70 * (-t161 * t251 - t163 * t252) + t318 * t219) * qJD(3) + g(1) * t163 + g(2) * t161 - g(3) * t219) * m(4); (t251 * t446 - t445 * t252) * m(5);];
tau = t1;
