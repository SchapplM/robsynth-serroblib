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
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:27:52
% EndTime: 2022-01-23 09:28:08
% DurationCPUTime: 13.94s
% Computational Cost: add. (13128->477), mult. (9630->583), div. (0->0), fcn. (7393->8), ass. (0->275)
t276 = sin(qJ(4));
t278 = cos(qJ(4));
t222 = Icges(6,5) * t278 - Icges(6,6) * t276;
t224 = Icges(5,5) * t278 - Icges(5,6) * t276;
t510 = t222 + t224;
t538 = Icges(5,5) + Icges(6,5);
t537 = Icges(5,6) + Icges(6,6);
t536 = Icges(5,3) + Icges(6,3);
t274 = qJ(1) + pkin(8);
t266 = qJ(3) + t274;
t260 = cos(t266);
t259 = sin(t266);
t405 = t259 * t278;
t406 = t259 * t276;
t134 = Icges(6,4) * t405 - Icges(6,2) * t406 - Icges(6,6) * t260;
t136 = Icges(5,4) * t405 - Icges(5,2) * t406 - Icges(5,6) * t260;
t529 = t134 + t136;
t215 = Icges(6,4) * t406;
t138 = Icges(6,1) * t405 - Icges(6,5) * t260 - t215;
t216 = Icges(5,4) * t406;
t140 = Icges(5,1) * t405 - Icges(5,5) * t260 - t216;
t534 = t138 + t140;
t426 = Icges(6,4) * t276;
t230 = Icges(6,1) * t278 - t426;
t313 = t230 * t260;
t139 = Icges(6,5) * t259 + t313;
t427 = Icges(5,4) * t276;
t232 = Icges(5,1) * t278 - t427;
t314 = t232 * t260;
t141 = Icges(5,5) * t259 + t314;
t527 = t139 + t141;
t225 = Icges(6,2) * t278 + t426;
t227 = Icges(5,2) * t278 + t427;
t533 = t225 + t227;
t267 = Icges(6,4) * t278;
t229 = Icges(6,1) * t276 + t267;
t268 = Icges(5,4) * t278;
t231 = Icges(5,1) * t276 + t268;
t532 = -t229 - t231;
t535 = t510 * t260;
t514 = t536 * t260 - t538 * t405 + t537 * t406;
t513 = t536 * t259 + t535;
t328 = -Icges(6,2) * t276 + t267;
t311 = t328 * t260;
t135 = Icges(6,6) * t259 + t311;
t329 = -Icges(5,2) * t276 + t268;
t312 = t329 * t260;
t137 = Icges(5,6) * t259 + t312;
t528 = t135 + t137;
t531 = t529 * t276;
t530 = t527 * t405;
t221 = Icges(6,5) * t276 + Icges(6,6) * t278;
t223 = Icges(5,5) * t276 + Icges(5,6) * t278;
t526 = t221 + t223;
t525 = t230 + t232;
t273 = qJD(1) + qJD(3);
t524 = -t533 * qJD(4) + t537 * t273;
t523 = t532 * qJD(4) + t538 * t273;
t505 = -t534 * t278 + t531;
t522 = t328 + t329;
t508 = t533 * t276 + t532 * t278;
t521 = -t513 * t260 + t530;
t402 = t260 * t278;
t520 = t514 * t259 - t534 * t402;
t467 = t513 * t259 + t527 * t402;
t519 = t528 * t276;
t477 = -t505 * t259 + t514 * t260;
t476 = -t528 * t406 + t521;
t403 = t260 * t276;
t475 = -t529 * t403 - t520;
t474 = -t528 * t403 + t467;
t473 = t534 * t276 + t529 * t278;
t472 = t527 * t276 + t528 * t278;
t408 = t259 * t273;
t518 = t524 * t260 - t522 * t408;
t517 = (t311 + t312) * t273 + t524 * t259;
t516 = t523 * t260 - t525 * t408;
t515 = (-t313 - t314) * t273 - t523 * t259;
t512 = t522 * qJD(4);
t511 = t525 * qJD(4);
t509 = -t526 * qJD(4) + t536 * t273;
t507 = t532 * t276 - t533 * t278;
t506 = -t527 * t278 + t519;
t410 = t223 * t260;
t413 = t221 * t260;
t471 = -t259 * t508 - t410 - t413;
t411 = t223 * t259;
t414 = t221 * t259;
t470 = -t260 * t508 + t411 + t414;
t253 = t260 * rSges(4,1);
t179 = -rSges(4,2) * t259 + t253;
t272 = qJDD(1) + qJDD(3);
t265 = cos(t274);
t258 = pkin(2) * t265;
t279 = cos(qJ(1));
t271 = t279 * pkin(1);
t262 = qJDD(1) * t271;
t281 = qJD(1) ^ 2;
t264 = sin(t274);
t277 = sin(qJ(1));
t439 = pkin(1) * t277;
t339 = -pkin(2) * t264 - t439;
t303 = qJDD(1) * t258 + t281 * t339 + t262;
t178 = rSges(4,1) * t259 + rSges(4,2) * t260;
t401 = t273 * t178;
t504 = t179 * t272 - t273 * t401 - g(2) + t303;
t254 = t259 * pkin(7);
t181 = t260 * pkin(3) + t254;
t151 = t181 * t273;
t372 = qJD(4) * t273;
t171 = -qJDD(4) * t260 + t259 * t372;
t269 = t278 * rSges(6,1);
t237 = -rSges(6,2) * t276 + t269;
t198 = t237 * qJD(4);
t431 = t278 * rSges(6,2);
t234 = rSges(6,1) * t276 + t431;
t241 = qJD(5) * t260;
t308 = (-qJDD(1) * t277 - t279 * t281) * pkin(1);
t291 = (-qJDD(1) * t264 - t265 * t281) * pkin(2) + t308;
t255 = t260 * pkin(7);
t180 = pkin(3) * t259 - t255;
t275 = -qJ(5) - pkin(7);
t233 = t260 * t275;
t270 = t278 * pkin(4);
t261 = t270 + pkin(3);
t107 = -rSges(6,1) * t405 + rSges(6,2) * t406 + t260 * rSges(6,3) - t259 * t261 - t233;
t491 = t180 + t107;
t366 = -t180 + t491;
t373 = qJD(4) * t260;
t399 = t278 * qJD(4) ^ 2;
t438 = pkin(3) - t261;
t371 = qJD(4) * t276;
t360 = t259 * t371;
t381 = pkin(4) * t360 + t241;
t407 = t259 * t275;
t370 = qJD(4) * t278;
t400 = t273 * t276;
t468 = t259 * t370 + t260 * t400;
t460 = -rSges(6,1) * t360 - rSges(6,2) * t468 - t273 * t407 - t381;
t465 = rSges(6,1) * t402 + t259 * rSges(6,3);
t436 = t460 + (-t260 * t438 - t254 + t465) * t273;
t14 = -t198 * t373 + qJDD(5) * t259 + t171 * t234 + (t171 * t276 - t260 * t399) * pkin(4) + t366 * t272 + (-t151 + t241 - t436) * t273 + t291;
t487 = t14 - g(1);
t170 = qJDD(4) * t259 + t260 * t372;
t404 = t260 * t273;
t210 = pkin(7) * t404;
t293 = t273 * (-pkin(3) * t408 + t210) + t272 * t181 + t303;
t348 = -pkin(4) * t276 - t234;
t108 = -rSges(6,2) * t403 + t260 * t261 - t407 + t465;
t393 = -t181 + t108;
t357 = t260 * t371;
t368 = t273 * t405;
t307 = -t357 - t368;
t356 = t260 * t370;
t240 = qJD(5) * t259;
t369 = t259 * t400;
t462 = rSges(6,2) * t369 + rSges(6,3) * t404 + t240;
t437 = -pkin(4) * t357 - t210 + (t259 * t438 - t233) * t273 + rSges(6,1) * t307 - rSges(6,2) * t356 + t462;
t15 = -qJDD(5) * t260 + t437 * t273 + t393 * t272 + t348 * t170 + (-pkin(4) * t399 - qJD(4) * t198 + qJD(5) * t273) * t259 + t293;
t486 = t15 - g(2);
t435 = rSges(5,1) * t278;
t238 = -rSges(5,2) * t276 + t435;
t199 = t238 * qJD(4);
t235 = rSges(5,1) * t276 + rSges(5,2) * t278;
t380 = rSges(5,2) * t406 + t260 * rSges(5,3);
t143 = rSges(5,1) * t405 - t380;
t384 = -t143 - t180;
t362 = rSges(5,1) * t360 + rSges(5,2) * t468;
t464 = rSges(5,1) * t402 + t259 * rSges(5,3);
t96 = t273 * t464 - t362;
t26 = -t199 * t373 + t171 * t235 + (-t151 - t96) * t273 + t384 * t272 + t291;
t503 = t26 - g(1);
t145 = -rSges(5,2) * t403 + t464;
t374 = qJD(4) * t259;
t315 = rSges(5,3) * t404 + (-t356 + t369) * rSges(5,2);
t94 = rSges(5,1) * t307 + t315;
t27 = t145 * t272 - t170 * t235 - t199 * t374 + t273 * t94 + t293;
t502 = t27 - g(2);
t207 = rSges(4,2) * t408;
t150 = rSges(4,1) * t404 - t207;
t501 = t150 * t273 + t178 * t272 + g(1) - t291;
t500 = t507 * qJD(4) + t526 * t273 - t512 * t276 + t511 * t278;
t499 = -t472 * qJD(4) + t513 * t273 - t518 * t276 + t516 * t278;
t498 = t473 * qJD(4) + t514 * t273 + t517 * t276 + t515 * t278;
t497 = t474 * t259 - t475 * t260;
t496 = t476 * t259 - t477 * t260;
t495 = t510 * qJD(4) + t508 * t273;
t494 = (t505 + t535) * t273 + t509 * t259;
t493 = t509 * t260 + t506 * t273 - t510 * t408;
t492 = t470 * t273;
t490 = t471 * t273;
t489 = -t494 * t259 + t498 * t260;
t488 = t493 * t259 + t499 * t260;
t485 = t498 * t259 + t494 * t260;
t484 = t499 * t259 - t493 * t260;
t483 = t496 * qJD(4) + t490;
t482 = t497 * qJD(4) + t492;
t481 = t505 * qJD(4) + t515 * t276 - t517 * t278;
t480 = -t506 * qJD(4) + t516 * t276 + t518 * t278;
t479 = t495 * t259 + t500 * t260;
t478 = t500 * t259 - t495 * t260;
t469 = t514 + t519;
t127 = t273 * t143;
t173 = t273 * t180;
t466 = -t127 - t173;
t183 = t265 * rSges(3,1) - rSges(3,2) * t264;
t175 = t183 + t271;
t463 = t237 + t270;
t461 = t258 + t271;
t458 = t273 * t491 - t173;
t110 = t145 + t181;
t453 = -t110 * t273 + t235 * t374;
t450 = t234 * t374 - t273 * (t181 + t393) + t381;
t386 = -Icges(5,2) * t405 + t140 - t216;
t390 = t231 * t259 + t136;
t449 = -t276 * t386 - t278 * t390;
t388 = -Icges(6,2) * t405 + t138 - t215;
t392 = t229 * t259 + t134;
t448 = -t276 * t388 - t278 * t392;
t447 = m(3) + m(4);
t446 = t170 / 0.2e1;
t445 = t171 / 0.2e1;
t317 = t339 * qJD(1);
t358 = t235 * t373;
t294 = t317 - t358;
t64 = t273 * t384 + t294;
t432 = t260 * t64;
t316 = t461 * qJD(1);
t125 = t179 * t273 + t316;
t421 = t125 * t178;
t412 = t222 * t273;
t409 = t224 * t273;
t391 = -t229 * t260 - t135;
t389 = -t231 * t260 - t137;
t387 = -t225 * t260 + t139;
t385 = -t227 * t260 + t141;
t379 = -t225 + t230;
t378 = t229 + t328;
t377 = -t227 + t232;
t376 = t231 + t329;
t353 = -pkin(3) - t435;
t352 = -t374 / 0.2e1;
t351 = t374 / 0.2e1;
t350 = -t373 / 0.2e1;
t349 = t373 / 0.2e1;
t341 = -pkin(4) * t403 - t234 * t260;
t340 = -pkin(4) * t370 - t198;
t239 = rSges(2,1) * t279 - rSges(2,2) * t277;
t236 = rSges(2,1) * t277 + rSges(2,2) * t279;
t182 = rSges(3,1) * t264 + rSges(3,2) * t265;
t65 = t316 - t453;
t332 = -t259 * t65 - t432;
t323 = t143 * t259 + t145 * t260;
t318 = -t431 + (-rSges(6,1) - pkin(4)) * t276;
t306 = t348 * t373 + t240;
t305 = -t276 * t387 + t278 * t391;
t304 = -t276 * t385 + t278 * t389;
t109 = t259 * t353 + t255 + t380;
t124 = t317 - t401;
t302 = (-t276 * t378 + t278 * t379) * t273;
t301 = (-t276 * t376 + t278 * t377) * t273;
t292 = t317 + t306;
t284 = ((t467 * t259 + ((t513 + t531) * t260 + t476 + t520 - t530) * t260) * qJD(4) + t492) * t349 + (-t508 * qJD(4) + t511 * t276 + t512 * t278) * t273 + (Icges(4,3) - t507) * t272 + (t470 + t472) * t446 + (t471 + t473) * t445 + (((t260 * t469 - t467 + t474) * t260 + (t259 * t469 + t475 - t521) * t259) * qJD(4) + t483 - t490) * t352 + (t479 + t480) * t351 + (t478 - t481 + t482) * t350;
t283 = t64 * t362 + t65 * (-rSges(5,1) * t357 + t210 + t315) + (t353 * t432 + (t64 * (-rSges(5,3) - pkin(7)) + t65 * t353) * t259) * t273;
t50 = t273 * t366 + t292;
t51 = t316 - t450;
t282 = t50 * (-rSges(6,3) * t408 - t460) + t51 * (-rSges(6,1) * t368 - t261 * t408 + t462) + ((t50 * (-t261 - t269) - t51 * t275) * t273 + t51 * t318 * qJD(4)) * t260;
t169 = t235 * t260;
t167 = t235 * t259;
t166 = t234 * t259;
t70 = qJD(4) * t323 + qJD(2);
t44 = qJD(2) + (-t259 * t491 + t393 * t260) * qJD(4);
t25 = t143 * t170 - t145 * t171 + qJDD(2) + (t259 * t96 + t260 * t94) * qJD(4);
t5 = qJDD(2) - t393 * t171 - t491 * t170 + (t259 * t436 + t260 * t437) * qJD(4);
t1 = [t284 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + (t124 * t207 + (-t124 * t253 - t421) * t273 + (-t124 * t461 + t125 * t339) * qJD(1) + t504 * (t179 + t461) - t501 * (-t178 + t339)) * m(4) + ((qJDD(1) * t183 - g(2) + t262) * t175 + (-g(1) - qJDD(1) * t182 + t308 + (-0.2e1 * t183 + 0.2e1 * t175 - t271) * t281) * (-t182 - t439)) * m(3) + (g(1) * t236 - g(2) * t239 + (t236 ^ 2 + t239 ^ 2) * qJDD(1)) * m(2) + ((t339 * t51 - t461 * t50) * qJD(1) + t282 - (-t50 + t292 + t458) * t51 + t486 * (t108 + t461) + t487 * (t107 + t339)) * m(6) + ((t339 * t65 - t461 * t64) * qJD(1) + t283 - (-t64 + t294 + t466) * t65 + t502 * (t110 + t461) + t503 * (t109 + t339)) * m(5); t447 * qJDD(2) + m(5) * t25 + m(6) * t5 + (-m(5) - m(6) - t447) * g(3); t284 + (t282 - t51 * (t306 + t458) - t450 * t50 + t486 * t108 + t487 * t107) * m(6) + (t283 - t64 * t453 - t65 * (-t358 + t466) + t502 * t110 + t503 * t109) * m(5) + (-t124 * t150 - t125 * t401 + t421 * t273 + (t124 * t273 + t504) * t179 + t501 * t178) * m(4); t497 * t446 + t496 * t445 + (t479 * t273 + t470 * t272 + t475 * t171 + t474 * t170 + (t488 * t259 + t489 * t260) * qJD(4)) * t259 / 0.2e1 - (t478 * t273 + t471 * t272 + t477 * t171 + t476 * t170 + (t484 * t259 + t485 * t260) * qJD(4)) * t260 / 0.2e1 + (t259 * t472 - t260 * t473) * t272 / 0.2e1 - (((t376 + t378) * t278 + (t377 + t379) * t276) * t273 + (((-t386 - t388) * t260 + (t385 + t387) * t259) * t278 + ((t390 + t392) * t260 + (t389 + t391) * t259) * t276) * qJD(4)) * t273 / 0.2e1 + ((t273 * t472 + t481) * t260 + (t273 * t473 + t480) * t259) * t273 / 0.2e1 + t483 * t408 / 0.2e1 + t482 * t404 / 0.2e1 + ((-t374 * t413 + t412) * t259 + (t302 + (-t448 * t260 + (t414 + t305) * t259) * qJD(4)) * t260 + (-t374 * t410 + t409) * t259 + (t301 + (-t449 * t260 + (t411 + t304) * t259) * qJD(4)) * t260) * t352 + ((t273 * t474 + t489) * t260 + (t273 * t475 + t488) * t259) * t351 + ((t273 * t476 + t485) * t260 + (t273 * t477 + t484) * t259) * t350 + ((-t373 * t414 - t412) * t260 + (t302 + (t305 * t259 + (t413 - t448) * t260) * qJD(4)) * t259 + (-t373 * t411 - t409) * t260 + (t301 + (t304 * t259 + (t410 - t449) * t260) * qJD(4)) * t259) * t349 + (t25 * t323 + t70 * ((t94 + t127) * t260 + (-t145 * t273 + t96) * t259) + t332 * t199 + ((-t273 * t65 - t26) * t260 + (t273 * t64 - t27) * t259) * t235 + g(1) * t169 + g(2) * t167 - g(3) * t238 - (t167 * t64 - t169 * t65) * t273 - (t70 * (-t167 * t259 - t169 * t260) + t332 * t238) * qJD(4)) * m(5) + (-(t50 * t166 + t341 * t51) * t273 - ((t341 * t44 - t463 * t50) * t260 + (-t51 * t463 + (-pkin(4) * t406 - t166) * t44) * t259) * qJD(4) + (t15 * t348 + t51 * t340 - t5 * t491 + t44 * t436 + (t50 * t234 - t393 * t44) * t273) * t259 + (t14 * t348 + t50 * t340 + t5 * t393 + t44 * t437 + (t348 * t51 - t44 * t491) * t273) * t260 - g(3) * t463 - (g(1) * t260 + g(2) * t259) * t318) * m(6); (t259 * t487 - t260 * t486) * m(6);];
tau = t1;
