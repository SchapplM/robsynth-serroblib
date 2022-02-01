% Calculate vector of inverse dynamics joint torques for
% S5RRRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR5_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:01:41
% EndTime: 2022-01-20 12:02:01
% DurationCPUTime: 11.18s
% Computational Cost: add. (23835->634), mult. (14122->793), div. (0->0), fcn. (10979->10), ass. (0->373)
t317 = qJ(1) + qJ(2);
t309 = qJ(3) + t317;
t297 = sin(t309);
t298 = cos(t309);
t217 = rSges(4,1) * t297 + rSges(4,2) * t298;
t315 = qJD(1) + qJD(2);
t303 = qJD(3) + t315;
t200 = t303 * t217;
t319 = sin(qJ(1));
t500 = pkin(1) * qJD(1);
t429 = t319 * t500;
t305 = sin(t317);
t463 = t305 * t315;
t434 = pkin(2) * t463;
t146 = -t429 - t434 - t200;
t316 = qJ(4) + qJ(5);
t306 = cos(t316);
t293 = Icges(6,4) * t306;
t304 = sin(t316);
t231 = Icges(6,1) * t304 + t293;
t385 = -Icges(6,2) * t304 + t293;
t548 = t231 + t385;
t318 = sin(qJ(4));
t473 = t297 * t318;
t443 = rSges(5,2) * t473 + t298 * rSges(5,3);
t320 = cos(qJ(4));
t472 = t297 * t320;
t172 = rSges(5,1) * t472 - t443;
t152 = t303 * t172;
t290 = t298 * pkin(8);
t221 = pkin(3) * t297 - t290;
t203 = t303 * t221;
t470 = t298 * t303;
t248 = pkin(8) * t470;
t436 = qJD(4) * t320;
t427 = rSges(5,2) * t436;
t465 = t303 * t318;
t375 = rSges(5,2) * t297 * t465 + rSges(5,3) * t470 - t298 * t427;
t437 = qJD(4) * t318;
t417 = t298 * t437;
t547 = -rSges(5,1) * t417 + t152 + t203 + t248 + t375;
t475 = t297 * t304;
t254 = rSges(6,2) * t475;
t474 = t297 * t306;
t160 = rSges(6,1) * t474 - t298 * rSges(6,3) - t254;
t145 = t303 * t160;
t322 = -pkin(9) - pkin(8);
t276 = t298 * t322;
t310 = t320 * pkin(4);
t299 = t310 + pkin(3);
t444 = -t297 * t299 - t276;
t148 = t221 + t444;
t426 = t303 * t474;
t447 = rSges(6,3) * t470 + t303 * t254;
t476 = t297 * t303;
t546 = -rSges(6,1) * t426 - t303 * t148 - t299 * t476 + t145 + t203 + t447;
t251 = Icges(6,4) * t475;
t158 = Icges(6,1) * t474 - Icges(6,5) * t298 - t251;
t156 = Icges(6,4) * t474 - Icges(6,2) * t475 - Icges(6,6) * t298;
t489 = t156 * t304;
t384 = -t158 * t306 + t489;
t367 = t384 * t297;
t228 = Icges(6,5) * t306 - Icges(6,6) * t304;
t369 = t228 * t298;
t155 = Icges(6,3) * t297 + t369;
t493 = Icges(6,4) * t304;
t232 = Icges(6,1) * t306 - t493;
t373 = t232 * t298;
t159 = Icges(6,5) * t297 + t373;
t468 = t298 * t306;
t459 = t297 * t155 + t159 * t468;
t545 = -t367 - t459;
t307 = cos(t317);
t234 = rSges(3,1) * t305 + rSges(3,2) * t307;
t480 = t234 * t315;
t187 = -t429 - t480;
t289 = t297 * pkin(8);
t418 = t297 * t437;
t258 = pkin(4) * t418;
t471 = t297 * t322;
t446 = -t303 * t471 - t258;
t504 = pkin(3) - t299;
t100 = (-t504 * t298 - t289) * t303 + t446;
t438 = qJD(4) * t303;
t191 = qJDD(4) * t297 + t298 * t438;
t435 = qJD(5) * t303;
t134 = qJDD(5) * t297 + t298 * t435 + t191;
t237 = t297 * t438;
t135 = t297 * t435 + t237 + (-qJDD(4) - qJDD(5)) * t298;
t222 = t298 * pkin(3) + t289;
t398 = t298 * t299 - t471;
t149 = t398 - t222;
t469 = t298 * t304;
t431 = rSges(6,2) * t469;
t534 = rSges(6,1) * t468 + t297 * rSges(6,3);
t161 = -t431 + t534;
t192 = -qJDD(4) * t298 + t237;
t314 = qJD(4) + qJD(5);
t215 = t297 * t314;
t216 = t298 * t314;
t501 = rSges(6,2) * t306;
t430 = t314 * t501;
t464 = t304 * t314;
t91 = -t298 * t430 + (-t298 * t464 - t426) * rSges(6,1) + t447;
t421 = -t303 * t431 + (-rSges(6,1) * t464 - t430) * t297;
t92 = t534 * t303 + t421;
t397 = pkin(4) * t417;
t99 = -t397 - t248 + (t504 * t297 - t276) * t303;
t10 = t134 * t160 - t135 * t161 - t148 * t191 - t149 * t192 + t215 * t92 + t216 * t91 + (t100 * t297 + t298 * t99) * qJD(4);
t233 = rSges(6,1) * t304 + t501;
t185 = t233 * t297;
t186 = t233 * t298;
t294 = t306 * rSges(6,1);
t235 = -rSges(6,2) * t304 + t294;
t52 = t160 * t215 + t161 * t216 + (-t148 * t297 + t149 * t298) * qJD(4);
t363 = -t216 * t233 - t397;
t347 = t363 - t434;
t336 = t347 - t429;
t423 = t148 - t160 - t221;
t58 = t423 * t303 + t336;
t321 = cos(qJ(1));
t428 = t321 * t500;
t462 = t307 * t315;
t433 = pkin(2) * t462;
t365 = t428 + t433;
t456 = t149 + t161;
t422 = t222 + t456;
t453 = t215 * t233 + t258;
t59 = t422 * t303 + t365 - t453;
t544 = -t58 * (t185 * t303 - t216 * t235) - t52 * (-t215 * t185 - t186 * t216) - t59 * (-t303 * t186 - t215 * t235) + t10 * (t297 * t160 + t298 * t161);
t178 = t222 * t303;
t207 = t235 * t314;
t313 = qJDD(1) + qJDD(2);
t302 = qJDD(3) + t313;
t312 = t315 ^ 2;
t324 = qJD(1) ^ 2;
t368 = (-qJDD(1) * t319 - t321 * t324) * pkin(1);
t334 = t368 + (-t305 * t313 - t307 * t312) * pkin(2);
t461 = t320 * qJD(4) ^ 2;
t17 = t135 * t233 - t207 * t216 + (t192 * t318 - t298 * t461) * pkin(4) + (-t100 - t178 - t92) * t303 + t423 * t302 + t334;
t543 = t17 - g(1);
t296 = pkin(2) * t307;
t311 = t321 * pkin(1);
t506 = pkin(1) * t319;
t396 = qJDD(1) * t311 - t324 * t506;
t505 = pkin(2) * t305;
t357 = t313 * t296 - t312 * t505 + t396;
t344 = t303 * (-pkin(3) * t476 + t248) + t302 * t222 + t357;
t18 = -t134 * t233 - t207 * t215 + (t91 + t99) * t303 + t456 * t302 + (-t191 * t318 - t297 * t461) * pkin(4) + t344;
t542 = t18 - g(2);
t424 = t298 * t465;
t420 = rSges(5,1) * t418 + rSges(5,2) * t424 + t297 * t427;
t466 = t298 * t320;
t533 = rSges(5,1) * t466 + t297 * rSges(5,3);
t109 = t533 * t303 - t420;
t503 = rSges(5,1) * t320;
t274 = -rSges(5,2) * t318 + t503;
t245 = t274 * qJD(4);
t272 = rSges(5,1) * t318 + rSges(5,2) * t320;
t439 = qJD(4) * t298;
t448 = -t172 - t221;
t44 = -t245 * t439 + t192 * t272 + (-t109 - t178) * t303 + t448 * t302 + t334;
t541 = t44 - g(1);
t108 = (-t303 * t472 - t417) * rSges(5,1) + t375;
t467 = t298 * t318;
t173 = -rSges(5,2) * t467 + t533;
t440 = qJD(4) * t297;
t45 = t108 * t303 + t173 * t302 - t191 * t272 - t245 * t440 + t344;
t540 = t45 - g(2);
t177 = rSges(4,1) * t470 - rSges(4,2) * t476;
t539 = -t177 * t303 - t217 * t302 - g(1) + t334;
t218 = t298 * rSges(4,1) - rSges(4,2) * t297;
t538 = -t200 * t303 + t218 * t302 - g(2) + t357;
t209 = rSges(3,1) * t462 - rSges(3,2) * t463;
t537 = -t209 * t315 - t234 * t313 - g(1) + t368;
t236 = t307 * rSges(3,1) - rSges(3,2) * t305;
t536 = t236 * t313 - t315 * t480 - g(2) + t396;
t308 = Icges(5,4) * t320;
t386 = -Icges(5,2) * t318 + t308;
t269 = Icges(5,1) * t318 + t308;
t328 = (t59 * (-pkin(4) * t437 - t233 * t314) + (t58 * (-t299 - t294) - t59 * t322) * t303) * t298;
t498 = t303 * t58;
t532 = t422 * t498 + t328;
t530 = -t434 + t547;
t227 = Icges(6,5) * t304 + Icges(6,6) * t306;
t353 = Icges(6,3) * t303 - t227 * t314;
t371 = t385 * t298;
t157 = Icges(6,6) * t297 + t371;
t488 = t157 * t304;
t529 = -t228 * t476 + t353 * t298 + t303 * (-t159 * t306 + t488);
t528 = t353 * t297 + (t369 + t384) * t303;
t266 = Icges(5,5) * t320 - Icges(5,6) * t318;
t265 = Icges(5,5) * t318 + Icges(5,6) * t320;
t350 = Icges(5,3) * t303 - t265 * qJD(4);
t494 = Icges(5,4) * t318;
t270 = Icges(5,1) * t320 - t494;
t374 = t270 * t298;
t171 = Icges(5,5) * t297 + t374;
t372 = t386 * t298;
t169 = Icges(5,6) * t297 + t372;
t485 = t169 * t318;
t381 = -t171 * t320 + t485;
t527 = -t266 * t476 + t350 * t298 + t303 * t381;
t370 = t266 * t298;
t261 = Icges(5,4) * t473;
t170 = Icges(5,1) * t472 - Icges(5,5) * t298 - t261;
t168 = Icges(5,4) * t472 - Icges(5,2) * t473 - Icges(5,6) * t298;
t486 = t168 * t318;
t382 = -t170 * t320 + t486;
t526 = t350 * t297 + (t370 + t382) * t303;
t133 = t173 + t222;
t389 = -t133 * t303 + t272 * t440;
t229 = Icges(6,2) * t306 + t493;
t379 = t229 * t304 - t231 * t306;
t525 = t228 * t314 + t303 * t379;
t267 = Icges(5,2) * t320 + t494;
t377 = t318 * t267 - t269 * t320;
t524 = t266 * qJD(4) + t377 * t303;
t166 = Icges(5,5) * t472 - Icges(5,6) * t473 - Icges(5,3) * t298;
t70 = -t166 * t298 - t382 * t297;
t523 = -t434 + t546;
t414 = -pkin(3) - t503;
t419 = t272 * t439;
t364 = -t419 - t434;
t345 = t364 - t429;
t77 = t448 * t303 + t345;
t499 = t298 * t77;
t78 = t365 - t389;
t329 = (t414 * t499 + (t77 * (-rSges(5,3) - pkin(8)) + t78 * t414) * t297) * t303;
t522 = t329 + (-t389 + t420) * t77;
t450 = -Icges(5,2) * t472 + t170 - t261;
t452 = t269 * t297 + t168;
t521 = -t318 * t450 - t320 * t452;
t520 = (-t229 * t298 + t159) * t215 - (-Icges(6,2) * t474 + t158 - t251) * t216 + t548 * t303;
t519 = t134 / 0.2e1;
t518 = t135 / 0.2e1;
t517 = t191 / 0.2e1;
t516 = t192 / 0.2e1;
t515 = -t215 / 0.2e1;
t514 = t215 / 0.2e1;
t513 = -t216 / 0.2e1;
t512 = t216 / 0.2e1;
t511 = t297 / 0.2e1;
t510 = -t298 / 0.2e1;
t509 = t302 / 0.2e1;
t508 = -t303 / 0.2e1;
t507 = t303 / 0.2e1;
t482 = t227 * t298;
t97 = -t379 * t297 - t482;
t497 = t97 * t303;
t478 = t265 * t298;
t117 = -t377 * t297 - t478;
t490 = t117 * t303;
t484 = t218 * t303;
t483 = t227 * t297;
t481 = t229 * t314;
t479 = t265 * t297;
t477 = t266 * t303;
t154 = Icges(6,5) * t474 - Icges(6,6) * t475 - Icges(6,3) * t298;
t460 = -t297 * t154 - t158 * t468;
t458 = -t297 * t166 - t170 * t466;
t167 = Icges(5,3) * t297 + t370;
t457 = t297 * t167 + t171 * t466;
t451 = -t269 * t298 - t169;
t449 = -t267 * t298 + t171;
t442 = -t267 + t270;
t441 = t269 + t386;
t432 = t297 * t92 + (t145 + t91) * t298;
t416 = t476 / 0.2e1;
t415 = t470 / 0.2e1;
t413 = -t440 / 0.2e1;
t412 = t440 / 0.2e1;
t411 = -t439 / 0.2e1;
t410 = t439 / 0.2e1;
t362 = -pkin(4) * t318 - t233;
t355 = Icges(6,5) * t303 - t231 * t314;
t408 = -t156 * t314 + t355 * t297 + t303 * t373;
t407 = -t157 * t314 - t232 * t476 + t355 * t298;
t354 = Icges(6,6) * t303 - t481;
t406 = t158 * t314 + t354 * t297 + t303 * t371;
t405 = t159 * t314 + t354 * t298 - t385 * t476;
t142 = t171 * t472;
t404 = t167 * t298 - t142;
t403 = -t154 + t488;
t401 = -t166 + t485;
t400 = t548 * t314;
t399 = t232 * t314 - t481;
t395 = -pkin(4) * t436 - t207;
t190 = t218 + t296;
t128 = t159 * t474;
t394 = t157 * t475 - t128;
t275 = rSges(2,1) * t321 - rSges(2,2) * t319;
t273 = rSges(2,1) * t319 + rSges(2,2) * t321;
t393 = -t297 * t59 - t298 * t58;
t71 = -t169 * t473 - t404;
t392 = t297 * t71 - t298 * t70;
t72 = -t168 * t467 - t458;
t73 = -t169 * t467 + t457;
t391 = t297 * t73 - t298 * t72;
t390 = -t297 * t78 - t499;
t93 = t156 * t306 + t158 * t304;
t110 = t168 * t320 + t170 * t318;
t111 = t169 * t320 + t171 * t318;
t380 = t172 * t297 + t173 * t298;
t378 = t267 * t320 + t269 * t318;
t189 = -t217 - t505;
t123 = -t160 + t444;
t361 = t215 * t482 - t216 * t483 - t228 * t303;
t360 = -t177 - t433;
t359 = -t449 * t318 + t451 * t320;
t132 = t414 * t297 + t290 + t443;
t126 = t133 + t296;
t358 = -rSges(6,3) * t476 - t421 - t446;
t124 = t161 + t398;
t115 = t123 - t505;
t356 = (-t441 * t318 + t442 * t320) * t303;
t352 = Icges(5,5) * t303 - qJD(4) * t269;
t351 = Icges(5,6) * t303 - t267 * qJD(4);
t116 = t124 + t296;
t338 = t154 * t303 - t406 * t304 + t408 * t306;
t13 = t528 * t297 + t338 * t298;
t337 = t155 * t303 - t405 * t304 + t407 * t306;
t14 = t529 * t297 + t337 * t298;
t15 = t338 * t297 - t528 * t298;
t16 = t337 * t297 - t529 * t298;
t62 = -t154 * t298 - t367;
t63 = -t155 * t298 - t394;
t30 = t215 * t63 - t216 * t62 + t497;
t64 = -t156 * t469 - t460;
t65 = -t157 * t469 + t459;
t98 = -t379 * t298 + t483;
t95 = t98 * t303;
t31 = t215 * t65 - t216 * t64 + t95;
t340 = (-t231 * t298 - t157) * t215 - (-t231 * t297 - t156) * t216 + (-t229 + t232) * t303;
t325 = -t520 * t304 + t340 * t306;
t40 = t304 * t408 + t306 * t406;
t41 = t304 * t407 + t306 * t405;
t335 = t227 * t303 - t400 * t304 + t399 * t306;
t46 = t525 * t297 + t335 * t298;
t47 = t335 * t297 - t525 * t298;
t94 = t157 * t306 + t159 * t304;
t346 = (-t13 * t216 + t134 * t65 + t135 * t64 + t14 * t215 + t302 * t98 + t303 * t46) * t511 + (-t361 * t297 + t325 * t298) * t515 + (t325 * t297 + t361 * t298) * t512 + (t134 * t63 + t135 * t62 - t15 * t216 + t16 * t215 + t302 * t97 + t303 * t47) * t510 + (t340 * t304 + t520 * t306) * t508 + t30 * t416 + t31 * t415 + ((t303 * t65 - t13) * t298 + (t303 * t64 + t14) * t297) * t514 + (t297 * t65 - t298 * t64) * t519 + (t297 * t63 - t298 * t62) * t518 + ((t303 * t63 - t15) * t298 + (t303 * t62 + t16) * t297) * t513 + (t297 * t94 - t298 * t93) * t509 + ((t303 * t94 - t40) * t298 + (t303 * t93 + t41) * t297) * t507;
t125 = t132 - t505;
t342 = t358 - t433;
t103 = t351 * t298 - t386 * t476;
t105 = -t270 * t476 + t352 * t298;
t333 = -t111 * qJD(4) - t103 * t318 + t105 * t320 + t167 * t303;
t104 = t351 * t297 + t303 * t372;
t106 = t352 * t297 + t303 * t374;
t332 = -t110 * qJD(4) - t104 * t318 + t106 * t320 + t166 * t303;
t240 = t386 * qJD(4);
t241 = t270 * qJD(4);
t331 = -t378 * qJD(4) - t240 * t318 + t241 * t320 + t265 * t303;
t118 = -t377 * t298 + t479;
t114 = t118 * t303;
t36 = qJD(4) * t392 + t490;
t37 = qJD(4) * t391 + t114;
t50 = -qJD(4) * t382 + t104 * t320 + t106 * t318;
t51 = -t381 * qJD(4) + t103 * t320 + t105 * t318;
t56 = t524 * t297 + t331 * t298;
t57 = t331 * t297 - t524 * t298;
t327 = (t114 + ((t71 - t142 + (t167 + t486) * t298 + t458) * t298 + t457 * t297) * qJD(4)) * t410 + (t95 + (t63 + (t155 + t489) * t298 + t394 + t460) * t216 + (-t403 * t298 - t545 + t62) * t215) * t512 + (t94 + t98) * t519 + (t93 + t97) * t518 + (t118 + t111) * t517 + (t117 + t110) * t516 + (-t497 + (t65 + t545) * t216 + (t403 * t297 - t128 + t64) * t215 + ((t155 + t384) * t215 + t403 * t216) * t298 + t30) * t515 + (t41 + t46) * t514 + (-t490 + ((t401 * t298 - t457 + t73) * t298 + (t401 * t297 + t404 + t72) * t297) * qJD(4) + t36) * t413 + (t51 + t56) * t412 + (-t377 * qJD(4) + t240 * t320 + t241 * t318 + t399 * t304 + t400 * t306) * t303 + (t40 + t47 + t31) * t513 + (t50 + t57 + t37) * t411 + (t229 * t306 + t231 * t304 + Icges(4,3) + t378) * t302;
t326 = Icges(3,3) * t313 + t327;
t202 = t272 * t298;
t201 = t272 * t297;
t188 = t236 * t315 + t428;
t147 = t365 + t484;
t96 = t380 * qJD(4);
t22 = t333 * t297 - t527 * t298;
t21 = t332 * t297 - t526 * t298;
t20 = t527 * t297 + t333 * t298;
t19 = t526 * t297 + t332 * t298;
t1 = [Icges(2,3) * qJDD(1) + t326 + (t536 * (t236 + t311) + t537 * (-t234 - t506) + (-t209 - t428 + t188) * t187) * m(3) + ((t273 ^ 2 + t275 ^ 2) * qJDD(1) + g(1) * t273 - g(2) * t275) * m(2) + (t58 * (t342 - t428) + t328 + (-t336 + t58 - t429 + t523) * t59 + t542 * (t116 + t311) + t543 * (t115 - t506)) * m(6) + (t77 * (-t365 + t420) + t329 + (-t345 + t77 - t429 + t530) * t78 + t540 * (t126 + t311) + t541 * (t125 - t506)) * m(5) + (t538 * (t190 + t311) + t539 * (t189 - t506) + (t360 - t428 + t147) * t146) * m(4); t326 + ((-t347 + t523) * t59 + (t433 - t453 + t342) * t58 + t542 * t116 + t543 * t115 + t532) * m(6) + ((-t364 + t530) * t78 + t540 * t126 + t541 * t125 + t522) * m(5) + (t538 * t190 + t539 * t189 + (t360 + t433 + t484) * t146) * m(4) + (-t187 * t209 - t188 * t480 + (t187 * t315 + t536) * t236 + (t188 * t315 - t537) * t234) * m(3); t327 + ((-t363 + t546) * t59 + (-t453 + t358) * t58 + t542 * t124 + t543 * t123 + t532) * m(6) + ((t419 + t547) * t78 + t540 * t133 + t541 * t132 + t522) * m(5) + (-t146 * t177 - t147 * t200 + (t146 * t303 + t538) * t218 + (t147 * t303 - t539) * t217) * m(4); ((-t439 * t479 - t477) * t298 + (t356 + (t359 * t297 + (t478 - t521) * t298) * qJD(4)) * t297) * t410 + t346 + ((-t440 * t478 + t477) * t297 + (t356 + (-t521 * t298 + (t479 + t359) * t297) * qJD(4)) * t298) * t413 + ((t303 * t71 - t21) * t298 + (t303 * t70 + t22) * t297) * t411 + ((t303 * t73 - t19) * t298 + (t303 * t72 + t20) * t297) * t412 + ((t318 * t442 + t320 * t441) * t303 + ((t297 * t449 - t298 * t450) * t320 + (t297 * t451 + t298 * t452) * t318) * qJD(4)) * t508 + (t117 * t302 + t191 * t71 + t192 * t70 + t303 * t57 + (-t21 * t298 + t22 * t297) * qJD(4)) * t510 + (-t110 * t298 + t111 * t297) * t509 + ((t111 * t303 - t50) * t298 + (t110 * t303 + t51) * t297) * t507 + t36 * t416 + t37 * t415 + (t118 * t302 + t191 * t73 + t192 * t72 + t303 * t56 + (-t19 * t298 + t20 * t297) * qJD(4)) * t511 + t391 * t517 + t392 * t516 + (-g(3) * (t235 + t310) - (g(1) * t298 + g(2) * t297) * t362 + t52 * t432 + (t17 * t362 + t58 * t395 + t10 * t149 + t52 * t99 + (-t52 * t148 + t362 * t59) * t303) * t298 + (t18 * t362 + t59 * t395 - t10 * t148 + t52 * t100 + (t58 * t233 - t456 * t52) * t303) * t297 - (-t59 * t424 + (t393 * t320 + t52 * (-t297 ^ 2 - t298 ^ 2) * t318) * qJD(4)) * pkin(4) + t544) * m(6) + (-(t201 * t77 - t202 * t78) * t303 - (t96 * (-t201 * t297 - t202 * t298) + t390 * t274) * qJD(4) + (t172 * t191 - t173 * t192 + (t108 * t298 + t109 * t297) * qJD(4)) * t380 + t96 * ((t108 + t152) * t298 + (-t173 * t303 + t109) * t297) + t390 * t245 + ((-t303 * t78 - t44) * t298 + (t303 * t77 - t45) * t297) * t272 + g(1) * t202 + g(2) * t201 - g(3) * t274) * m(5); t346 + (t52 * (-t161 * t476 + t432) + t393 * t207 + ((-t303 * t59 - t17) * t298 + (-t18 + t498) * t297) * t233 + g(1) * t186 + g(2) * t185 - g(3) * t235 + t544) * m(6);];
tau = t1;
