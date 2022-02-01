% Calculate vector of inverse dynamics joint torques for
% S5RRPRR4
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
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:04
% EndTime: 2022-01-20 10:48:21
% DurationCPUTime: 9.50s
% Computational Cost: add. (19139->625), mult. (11926->787), div. (0->0), fcn. (9239->10), ass. (0->365)
t312 = sin(qJ(1));
t314 = cos(qJ(1));
t317 = qJD(1) ^ 2;
t350 = (-qJDD(1) * t312 - t314 * t317) * pkin(1);
t533 = t350 - g(1);
t309 = qJ(4) + qJ(5);
t300 = cos(t309);
t287 = Icges(6,4) * t300;
t298 = sin(t309);
t225 = Icges(6,1) * t298 + t287;
t370 = -Icges(6,2) * t298 + t287;
t532 = t225 + t370;
t310 = qJ(1) + qJ(2);
t297 = pkin(9) + t310;
t290 = sin(t297);
t464 = t290 * t298;
t237 = Icges(6,4) * t464;
t291 = cos(t297);
t463 = t290 * t300;
t153 = Icges(6,1) * t463 - Icges(6,5) * t291 - t237;
t151 = Icges(6,4) * t463 - Icges(6,2) * t464 - Icges(6,6) * t291;
t475 = t151 * t298;
t369 = -t153 * t300 + t475;
t349 = t369 * t290;
t222 = Icges(6,5) * t300 - Icges(6,6) * t298;
t351 = t222 * t291;
t150 = Icges(6,3) * t290 + t351;
t479 = Icges(6,4) * t298;
t226 = Icges(6,1) * t300 - t479;
t355 = t226 * t291;
t154 = Icges(6,5) * t290 + t355;
t456 = t291 * t300;
t446 = t150 * t290 + t154 * t456;
t531 = -t349 - t446;
t419 = rSges(6,1) * t463;
t315 = -pkin(8) - pkin(7);
t271 = t291 * t315;
t313 = cos(qJ(4));
t303 = t313 * pkin(4);
t293 = t303 + pkin(3);
t432 = -t290 * t293 - t271;
t299 = sin(t310);
t489 = pkin(2) * t299;
t530 = -t419 - t489 + t432;
t484 = pkin(1) * qJD(1);
t416 = t312 * t484;
t301 = cos(t310);
t228 = rSges(3,1) * t299 + rSges(3,2) * t301;
t308 = qJD(1) + qJD(2);
t468 = t228 * t308;
t181 = -t416 - t468;
t451 = t299 * t308;
t421 = pkin(2) * t451;
t529 = t421 - t416;
t455 = t291 * t308;
t253 = pkin(7) * t455;
t311 = sin(qJ(4));
t425 = qJD(4) * t311;
t407 = t291 * t425;
t384 = pkin(4) * t407;
t488 = pkin(3) - t293;
t113 = -t384 - t253 + (t290 * t488 - t271) * t308;
t284 = t290 * pkin(7);
t408 = t290 * t425;
t247 = pkin(4) * t408;
t459 = t290 * t315;
t433 = -t308 * t459 - t247;
t114 = (-t291 * t488 - t284) * t308 + t433;
t426 = qJD(4) * t308;
t194 = qJDD(4) * t290 + t291 * t426;
t423 = qJD(5) * t308;
t133 = qJDD(5) * t290 + t291 * t423 + t194;
t248 = t290 * t426;
t134 = t290 * t423 + t248 + (-qJDD(4) - qJDD(5)) * t291;
t285 = t291 * pkin(7);
t212 = pkin(3) * t290 - t285;
t142 = t212 + t432;
t286 = t291 * pkin(3);
t213 = t286 + t284;
t239 = t291 * t293;
t385 = t239 - t459;
t143 = t385 - t213;
t240 = rSges(6,2) * t464;
t517 = -rSges(6,3) * t291 - t240;
t156 = t419 + t517;
t358 = rSges(6,1) * t456 + rSges(6,3) * t290;
t457 = t291 * t298;
t418 = rSges(6,2) * t457;
t157 = -t418 + t358;
t195 = -qJDD(4) * t291 + t248;
t307 = qJD(4) + qJD(5);
t208 = t290 * t307;
t209 = t291 * t307;
t485 = rSges(6,2) * t300;
t417 = t307 * t485;
t435 = rSges(6,3) * t455 + t240 * t308;
t452 = t298 * t307;
t462 = t290 * t308;
t95 = -t291 * t417 + (-t291 * t452 - t300 * t462) * rSges(6,1) + t435;
t411 = -t308 * t418 + (-rSges(6,1) * t452 - t417) * t290;
t96 = t308 * t358 + t411;
t10 = t133 * t156 - t134 * t157 - t142 * t194 - t143 * t195 + t208 * t96 + t209 * t95 + qJDD(3) + (t113 * t291 + t114 * t290) * qJD(4);
t227 = rSges(6,1) * t298 + t485;
t179 = t227 * t290;
t180 = t227 * t291;
t229 = rSges(6,1) * t300 - rSges(6,2) * t298;
t53 = t156 * t208 + t157 * t209 + qJD(3) + (-t142 * t290 + t143 * t291) * qJD(4);
t348 = -t209 * t227 - t384;
t334 = t348 - t416;
t398 = -t212 - t489;
t361 = t142 - t156 + t398;
t59 = t308 * t361 + t334;
t415 = t314 * t484;
t292 = pkin(2) * t301;
t397 = t213 + t292;
t443 = t143 + t157;
t507 = t308 * (t397 + t443) - t208 * t227 - t247;
t60 = t415 + t507;
t528 = -(t179 * t308 - t209 * t229) * t59 - t53 * (-t179 * t208 - t180 * t209) - t60 * (-t180 * t308 - t208 * t229) + t10 * (t156 * t290 + t157 * t291);
t183 = t213 * t308;
t201 = t229 * t307;
t306 = qJDD(1) + qJDD(2);
t305 = t308 ^ 2;
t450 = t301 * t305;
t333 = -pkin(2) * t450 + t350;
t448 = t313 * qJD(4) ^ 2;
t17 = t134 * t227 - t201 * t209 + (t195 * t311 - t291 * t448) * pkin(4) + (-t114 - t183 - t96) * t308 + t361 * t306 + t333;
t527 = t17 - g(1);
t304 = t314 * pkin(1);
t490 = pkin(1) * t312;
t383 = qJDD(1) * t304 - t317 * t490;
t343 = t292 * t306 - t305 * t489 + t383;
t332 = t308 * (-pkin(3) * t462 + t253) + t306 * t213 + t343;
t18 = -t133 * t227 - t201 * t208 + (t113 + t95) * t308 + t443 * t306 + (-t194 * t311 - t290 * t448) * pkin(4) + t332;
t526 = t18 - g(2);
t449 = t308 * t311;
t412 = t291 * t449;
t424 = qJD(4) * t313;
t414 = rSges(5,2) * t424;
t410 = rSges(5,1) * t408 + rSges(5,2) * t412 + t290 * t414;
t453 = t291 * t313;
t258 = rSges(5,1) * t453;
t516 = rSges(5,3) * t290 + t258;
t112 = t308 * t516 - t410;
t487 = rSges(5,1) * t313;
t269 = -rSges(5,2) * t311 + t487;
t242 = t269 * qJD(4);
t267 = rSges(5,1) * t311 + rSges(5,2) * t313;
t461 = t290 * t311;
t431 = rSges(5,2) * t461 + rSges(5,3) * t291;
t460 = t290 * t313;
t168 = rSges(5,1) * t460 - t431;
t382 = -t168 + t398;
t427 = qJD(4) * t291;
t45 = -t242 * t427 + t195 * t267 + (-t112 - t183) * t308 + t382 * t306 + t333;
t525 = t45 - g(1);
t357 = rSges(5,2) * t290 * t449 + rSges(5,3) * t455 - t291 * t414;
t111 = (-t308 * t460 - t407) * rSges(5,1) + t357;
t454 = t291 * t311;
t169 = -rSges(5,2) * t454 + t516;
t428 = qJD(4) * t290;
t46 = t111 * t308 + t169 * t306 - t194 * t267 - t242 * t428 + t332;
t524 = t46 - g(2);
t210 = rSges(4,1) * t290 + rSges(4,2) * t291;
t250 = rSges(4,2) * t462;
t523 = -t306 * t210 - t308 * (rSges(4,1) * t455 - t250) + (-t299 * t306 - t450) * pkin(2) + t533;
t283 = t291 * rSges(4,1);
t211 = -rSges(4,2) * t290 + t283;
t522 = -t210 * t305 + t211 * t306 - g(2) + t343;
t289 = t301 * rSges(3,1);
t203 = -rSges(3,2) * t451 + t289 * t308;
t521 = -t203 * t308 - t228 * t306 + t533;
t230 = -rSges(3,2) * t299 + t289;
t520 = t230 * t306 - t308 * t468 - g(2) + t383;
t519 = t169 + t397;
t518 = t211 + t292;
t302 = Icges(5,4) * t313;
t371 = -Icges(5,2) * t311 + t302;
t264 = Icges(5,1) * t311 + t302;
t155 = t308 * t168;
t204 = t308 * t212;
t515 = -rSges(5,1) * t407 + t155 + t204 + t253 + t357;
t221 = Icges(6,5) * t298 + Icges(6,6) * t300;
t339 = Icges(6,3) * t308 - t221 * t307;
t353 = t370 * t291;
t152 = Icges(6,6) * t290 + t353;
t474 = t152 * t298;
t514 = -t222 * t462 + t291 * t339 + t308 * (-t154 * t300 + t474);
t513 = t290 * t339 + (t351 + t369) * t308;
t261 = Icges(5,5) * t313 - Icges(5,6) * t311;
t260 = Icges(5,5) * t311 + Icges(5,6) * t313;
t336 = Icges(5,3) * t308 - qJD(4) * t260;
t480 = Icges(5,4) * t311;
t265 = Icges(5,1) * t313 - t480;
t356 = t265 * t291;
t166 = Icges(5,5) * t290 + t356;
t354 = t371 * t291;
t164 = Icges(5,6) * t290 + t354;
t472 = t164 * t311;
t366 = -t166 * t313 + t472;
t512 = -t261 * t462 + t291 * t336 + t308 * t366;
t352 = t261 * t291;
t256 = Icges(5,4) * t461;
t165 = Icges(5,1) * t460 - Icges(5,5) * t291 - t256;
t163 = Icges(5,4) * t460 - Icges(5,2) * t461 - Icges(5,6) * t291;
t473 = t163 * t311;
t367 = -t165 * t313 + t473;
t511 = t290 * t336 + (t352 + t367) * t308;
t223 = Icges(6,2) * t300 + t479;
t364 = t223 * t298 - t225 * t300;
t510 = t222 * t307 + t308 * t364;
t509 = -t267 * t428 + t308 * t519;
t262 = Icges(5,2) * t313 + t480;
t362 = t262 * t311 - t264 * t313;
t508 = qJD(4) * t261 + t308 * t362;
t161 = Icges(5,5) * t460 - Icges(5,6) * t461 - Icges(5,3) * t291;
t69 = -t161 * t291 - t290 * t367;
t144 = t308 * t156;
t506 = -t142 * t308 + t144 + t204 + t435;
t437 = -Icges(5,2) * t460 + t165 - t256;
t439 = t264 * t290 + t163;
t505 = -t311 * t437 - t313 * t439;
t504 = t208 * (-t223 * t291 + t154) - t209 * (-Icges(6,2) * t463 + t153 - t237) + t308 * t532;
t503 = t133 / 0.2e1;
t502 = t134 / 0.2e1;
t501 = t194 / 0.2e1;
t500 = t195 / 0.2e1;
t499 = -t208 / 0.2e1;
t498 = t208 / 0.2e1;
t497 = -t209 / 0.2e1;
t496 = t209 / 0.2e1;
t495 = t290 / 0.2e1;
t494 = -t291 / 0.2e1;
t493 = t306 / 0.2e1;
t492 = -t308 / 0.2e1;
t491 = t308 / 0.2e1;
t470 = t221 * t291;
t98 = -t290 * t364 - t470;
t483 = t98 * t308;
t466 = t260 * t291;
t118 = -t290 * t362 - t466;
t476 = t118 * t308;
t471 = t221 * t290;
t469 = t223 * t307;
t467 = t260 * t290;
t465 = t261 * t308;
t149 = Icges(6,5) * t463 - Icges(6,6) * t464 - Icges(6,3) * t291;
t447 = -t149 * t290 - t153 * t456;
t445 = -t290 * t161 - t165 * t453;
t162 = Icges(5,3) * t290 + t352;
t444 = t162 * t290 + t166 * t453;
t438 = -t264 * t291 - t164;
t436 = -t262 * t291 + t166;
t430 = -t262 + t265;
t429 = t264 + t371;
t79 = t415 + t509;
t422 = t79 * t489;
t420 = t290 * t96 + (t144 + t95) * t291;
t409 = t267 * t427;
t406 = t462 / 0.2e1;
t405 = t455 / 0.2e1;
t404 = -pkin(3) - t487;
t403 = -t428 / 0.2e1;
t402 = t428 / 0.2e1;
t401 = -t427 / 0.2e1;
t400 = t427 / 0.2e1;
t184 = -t210 - t489;
t346 = -pkin(4) * t311 - t227;
t341 = Icges(6,5) * t308 - t225 * t307;
t395 = -t151 * t307 + t290 * t341 + t308 * t355;
t394 = -t152 * t307 - t226 * t462 + t291 * t341;
t340 = Icges(6,6) * t308 - t469;
t393 = t153 * t307 + t290 * t340 + t308 * t353;
t392 = t154 * t307 + t291 * t340 - t370 * t462;
t139 = t166 * t460;
t391 = t162 * t291 - t139;
t390 = -t149 + t474;
t388 = -t161 + t472;
t387 = t532 * t307;
t386 = t226 * t307 - t469;
t380 = -pkin(4) * t424 - t201;
t126 = t154 * t463;
t379 = t152 * t464 - t126;
t270 = rSges(2,1) * t314 - rSges(2,2) * t312;
t268 = rSges(2,1) * t312 + rSges(2,2) * t314;
t378 = -t290 * t60 - t291 * t59;
t70 = -t164 * t461 - t391;
t377 = t290 * t70 - t291 * t69;
t71 = -t163 * t454 - t445;
t72 = -t164 * t454 + t444;
t376 = t290 * t72 - t291 * t71;
t347 = -t409 - t416;
t78 = t308 * t382 + t347;
t375 = -t290 * t79 - t291 * t78;
t374 = -t411 - t433;
t84 = t151 * t300 + t153 * t298;
t100 = t163 * t313 + t165 * t311;
t101 = t164 * t313 + t166 * t311;
t365 = t168 * t290 + t169 * t291;
t363 = t262 * t313 + t264 * t311;
t345 = t208 * t470 - t209 * t471 - t222 * t308;
t344 = -t311 * t436 + t313 * t438;
t115 = -t517 + t530;
t342 = (-t311 * t429 + t313 * t430) * t308;
t338 = Icges(5,5) * t308 - qJD(4) * t264;
t337 = Icges(5,6) * t308 - qJD(4) * t262;
t116 = t157 + t292 + t385;
t328 = t149 * t308 - t298 * t393 + t300 * t395;
t13 = t290 * t513 + t291 * t328;
t327 = t150 * t308 - t298 * t392 + t300 * t394;
t14 = t290 * t514 + t291 * t327;
t15 = t290 * t328 - t291 * t513;
t16 = t290 * t327 - t291 * t514;
t63 = -t149 * t291 - t349;
t64 = -t150 * t291 - t379;
t30 = t208 * t64 - t209 * t63 + t483;
t65 = -t151 * t457 - t447;
t66 = -t152 * t457 + t446;
t99 = -t291 * t364 + t471;
t97 = t99 * t308;
t31 = t208 * t66 - t209 * t65 + t97;
t329 = (-t225 * t291 - t152) * t208 - (-t225 * t290 - t151) * t209 + (-t223 + t226) * t308;
t319 = -t298 * t504 + t300 * t329;
t40 = t298 * t395 + t300 * t393;
t41 = t298 * t394 + t300 * t392;
t326 = t221 * t308 - t298 * t387 + t300 * t386;
t47 = t290 * t510 + t291 * t326;
t48 = t290 * t326 - t291 * t510;
t85 = t152 * t300 + t154 * t298;
t335 = (-t13 * t209 + t133 * t66 + t134 * t65 + t14 * t208 + t306 * t99 + t308 * t47) * t495 + (-t290 * t345 + t291 * t319) * t499 + (t290 * t319 + t291 * t345) * t496 + (t133 * t64 + t134 * t63 - t15 * t209 + t16 * t208 + t306 * t98 + t308 * t48) * t494 + (t298 * t329 + t300 * t504) * t492 + t30 * t406 + t31 * t405 + ((t308 * t66 - t13) * t291 + (t308 * t65 + t14) * t290) * t498 + (t290 * t66 - t291 * t65) * t503 + (t290 * t64 - t291 * t63) * t502 + ((t308 * t64 - t15) * t291 + (t308 * t63 + t16) * t290) * t497 + (t290 * t85 - t291 * t84) * t493 + ((t308 * t85 - t40) * t291 + (t308 * t84 + t41) * t290) * t491;
t124 = t290 * t404 + t285 + t431 - t489;
t107 = t291 * t337 - t371 * t462;
t109 = -t265 * t462 + t291 * t338;
t325 = -qJD(4) * t101 - t107 * t311 + t109 * t313 + t162 * t308;
t108 = t290 * t337 + t308 * t354;
t110 = t290 * t338 + t308 * t356;
t324 = -qJD(4) * t100 - t108 * t311 + t110 * t313 + t161 * t308;
t233 = t371 * qJD(4);
t234 = t265 * qJD(4);
t323 = -qJD(4) * t363 - t233 * t311 + t234 * t313 + t260 * t308;
t147 = t184 * t308 - t416;
t148 = t308 * t518 + t415;
t322 = (t147 * (-t283 - t292) + t148 * t184) * t308;
t321 = (t78 * (-t258 - t286 - t292) - t422 + (t78 * (-rSges(5,3) - pkin(7)) + t79 * t404) * t290) * t308;
t119 = -t291 * t362 + t467;
t117 = t119 * t308;
t36 = qJD(4) * t377 + t476;
t37 = qJD(4) * t376 + t117;
t51 = -qJD(4) * t367 + t108 * t313 + t110 * t311;
t52 = -qJD(4) * t366 + t107 * t313 + t109 * t311;
t57 = t290 * t508 + t291 * t323;
t58 = t290 * t323 - t291 * t508;
t320 = (t117 + ((t70 - t139 + (t162 + t473) * t291 + t445) * t291 + t444 * t290) * qJD(4)) * t400 + (t97 + (t64 + (t150 + t475) * t291 + t379 + t447) * t209 + (-t291 * t390 - t531 + t63) * t208) * t496 + (t85 + t99) * t503 + (t84 + t98) * t502 + (t119 + t101) * t501 + (t118 + t100) * t500 + (t30 - t483 + (t66 + t531) * t209 + (t290 * t390 - t126 + t65) * t208 + ((t150 + t369) * t208 + t390 * t209) * t291) * t499 + (t41 + t47) * t498 + (t36 - t476 + ((t291 * t388 - t444 + t72) * t291 + (t290 * t388 + t391 + t71) * t290) * qJD(4)) * t403 + (t52 + t57) * t402 + (-qJD(4) * t362 + t233 * t313 + t234 * t311 + t298 * t386 + t300 * t387) * t308 + (t40 + t48 + t31) * t497 + (t51 + t58 + t37) * t401 + (t223 * t300 + t225 * t298 + Icges(3,3) + Icges(4,3) + t363) * t306;
t318 = (t59 * (-t239 - t358 - t292) + t60 * t530) * t308 + t60 * (-pkin(4) * t425 - t227 * t307) * t291;
t200 = t308 * t210;
t193 = t267 * t291;
t192 = t267 * t290;
t182 = t230 * t308 + t415;
t88 = qJD(4) * t365 + qJD(3);
t44 = t168 * t194 - t169 * t195 + qJDD(3) + (t111 * t291 + t112 * t290) * qJD(4);
t22 = t290 * t325 - t291 * t512;
t21 = t290 * t324 - t291 * t511;
t20 = t290 * t512 + t291 * t325;
t19 = t290 * t511 + t291 * t324;
t1 = [Icges(2,3) * qJDD(1) + t320 + (t520 * (t230 + t304) + t521 * (-t228 - t490) + (-t203 - t415 + t182) * t181) * m(3) + ((t268 ^ 2 + t270 ^ 2) * qJDD(1) + g(1) * t268 - g(2) * t270) * m(2) + (t59 * (t374 - t415) + t318 + (-t334 + t59 + t506 + t529) * t60 + t526 * (t116 + t304) + t527 * (t115 - t490)) * m(6) + (t78 * (t410 - t415) + t321 + (-t347 + t78 + t515 + t529) * t79 + t524 * (t519 + t304) + t525 * (t124 - t490)) * m(5) + (t147 * (t250 - t415) + t322 + t522 * (t518 + t304) + t523 * (t184 - t490) + (t147 + t200 + t421) * t148) * m(4); t320 + (t318 + (t308 * t489 - t348 + t506) * t60 + (t374 + t507) * t59 + t526 * t116 + t527 * t115) * m(6) + (t422 * t308 + t321 + (t409 + t515) * t79 + (t410 + t509) * t78 + t524 * t519 + t525 * t124) * m(5) + (t147 * t250 + t322 + t148 * t200 - (-t147 * t518 - t148 * t489) * t308 + t522 * t518 + t523 * t184) * m(4) + (-t181 * t203 - t182 * t468 + (t181 * t308 + t520) * t230 + (t182 * t308 - t521) * t228) * m(3); m(4) * qJDD(3) + m(5) * t44 + m(6) * t10 + (-m(4) - m(5) - m(6)) * g(3); ((t308 * t70 - t21) * t291 + (t308 * t69 + t22) * t290) * t401 + ((-t428 * t466 + t465) * t290 + (t342 + (-t505 * t291 + (t467 + t344) * t290) * qJD(4)) * t291) * t403 + (t119 * t306 + t194 * t72 + t195 * t71 + t308 * t57 + (-t19 * t291 + t20 * t290) * qJD(4)) * t495 + ((-t427 * t467 - t465) * t291 + (t342 + (t344 * t290 + (t466 - t505) * t291) * qJD(4)) * t290) * t400 + ((t308 * t72 - t19) * t291 + (t308 * t71 + t20) * t290) * t402 + (t118 * t306 + t194 * t70 + t195 * t69 + t308 * t58 + (-t21 * t291 + t22 * t290) * qJD(4)) * t494 + t335 + t36 * t406 + t37 * t405 + ((t311 * t430 + t313 * t429) * t308 + ((t290 * t436 - t291 * t437) * t313 + (t290 * t438 + t291 * t439) * t311) * qJD(4)) * t492 + ((t101 * t308 - t51) * t291 + (t100 * t308 + t52) * t290) * t491 + t376 * t501 + t377 * t500 + (-t100 * t291 + t101 * t290) * t493 + (-g(3) * (t229 + t303) - (g(1) * t291 + g(2) * t290) * t346 - (-t60 * t412 + (t378 * t313 + t53 * (-t290 ^ 2 - t291 ^ 2) * t311) * qJD(4)) * pkin(4) + t53 * t420 + (t17 * t346 + t59 * t380 + t10 * t143 + t53 * t113 + (-t142 * t53 + t346 * t60) * t308) * t291 + (t18 * t346 + t60 * t380 - t10 * t142 + t53 * t114 + (t227 * t59 - t443 * t53) * t308) * t290 + t528) * m(6) + (-(t192 * t78 - t193 * t79) * t308 - (t88 * (-t192 * t290 - t193 * t291) + t375 * t269) * qJD(4) + t44 * t365 + t88 * ((t111 + t155) * t291 + (-t169 * t308 + t112) * t290) + t375 * t242 + ((-t308 * t79 - t45) * t291 + (t308 * t78 - t46) * t290) * t267 + g(1) * t193 + g(2) * t192 - g(3) * t269) * m(5); t335 + (t53 * (-t157 * t462 + t420) + t378 * t201 + ((-t308 * t60 - t17) * t291 + (t308 * t59 - t18) * t290) * t227 + g(1) * t180 + g(2) * t179 - g(3) * t229 + t528) * m(6);];
tau = t1;
