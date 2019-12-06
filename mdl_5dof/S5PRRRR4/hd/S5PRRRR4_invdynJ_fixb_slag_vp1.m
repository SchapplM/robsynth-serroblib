% Calculate vector of inverse dynamics joint torques for
% S5PRRRR4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:08
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRRR4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR4_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:07:40
% EndTime: 2019-12-05 17:07:52
% DurationCPUTime: 7.50s
% Computational Cost: add. (18747->583), mult. (11500->748), div. (0->0), fcn. (9027->8), ass. (0->338)
t293 = qJ(4) + qJ(5);
t286 = cos(t293);
t276 = Icges(6,4) * t286;
t285 = sin(t293);
t220 = Icges(6,1) * t285 + t276;
t346 = -Icges(6,2) * t285 + t276;
t496 = t220 + t346;
t290 = pkin(9) + qJ(2);
t284 = qJ(3) + t290;
t278 = sin(t284);
t434 = t278 * t285;
t230 = Icges(6,4) * t434;
t279 = cos(t284);
t433 = t278 * t286;
t149 = Icges(6,1) * t433 - Icges(6,5) * t279 - t230;
t147 = Icges(6,4) * t433 - Icges(6,2) * t434 - Icges(6,6) * t279;
t445 = t147 * t285;
t345 = -t149 * t286 + t445;
t328 = t345 * t278;
t217 = Icges(6,5) * t286 - Icges(6,6) * t285;
t330 = t217 * t279;
t146 = Icges(6,3) * t278 + t330;
t449 = Icges(6,4) * t285;
t221 = Icges(6,1) * t286 - t449;
t334 = t221 * t279;
t150 = Icges(6,5) * t278 + t334;
t426 = t279 * t286;
t418 = t278 * t146 + t150 * t426;
t495 = -t328 - t418;
t282 = sin(t290);
t456 = pkin(2) * qJD(2);
t389 = t282 * t456;
t203 = rSges(4,1) * t278 + rSges(4,2) * t279;
t292 = qJD(2) + qJD(3);
t441 = t203 * t292;
t166 = -t389 - t441;
t425 = t279 * t292;
t244 = pkin(7) * t425;
t296 = -pkin(8) - pkin(7);
t258 = t279 * t296;
t294 = sin(qJ(4));
t396 = qJD(4) * t294;
t377 = t279 * t396;
t357 = pkin(4) * t377;
t295 = cos(qJ(4));
t288 = t295 * pkin(4);
t280 = t288 + pkin(3);
t460 = pkin(3) - t280;
t111 = -t357 - t244 + (t278 * t460 - t258) * t292;
t272 = t278 * pkin(7);
t378 = t278 * t396;
t238 = pkin(4) * t378;
t429 = t278 * t296;
t404 = -t292 * t429 - t238;
t112 = (-t279 * t460 - t272) * t292 + t404;
t397 = qJD(4) * t292;
t190 = qJDD(4) * t278 + t279 * t397;
t394 = qJD(5) * t292;
t131 = qJDD(5) * t278 + t279 * t394 + t190;
t239 = t278 * t397;
t132 = t278 * t394 + t239 + (-qJDD(4) - qJDD(5)) * t279;
t273 = t279 * pkin(7);
t205 = pkin(3) * t278 - t273;
t403 = -t278 * t280 - t258;
t140 = t205 + t403;
t206 = t279 * pkin(3) + t272;
t358 = t279 * t280 - t429;
t141 = t358 - t206;
t233 = rSges(6,2) * t434;
t152 = rSges(6,1) * t433 - t279 * rSges(6,3) - t233;
t427 = t279 * t285;
t391 = rSges(6,2) * t427;
t486 = rSges(6,1) * t426 + t278 * rSges(6,3);
t153 = -t391 + t486;
t191 = -qJDD(4) * t279 + t239;
t291 = qJD(4) + qJD(5);
t199 = t278 * t291;
t200 = t279 * t291;
t432 = t278 * t292;
t386 = t286 * t432;
t457 = rSges(6,2) * t286;
t390 = t291 * t457;
t406 = rSges(6,3) * t425 + t292 * t233;
t422 = t285 * t291;
t93 = -t279 * t390 + (-t279 * t422 - t386) * rSges(6,1) + t406;
t381 = -t292 * t391 + (-rSges(6,1) * t422 - t390) * t278;
t94 = t292 * t486 + t381;
t10 = t131 * t152 - t132 * t153 - t140 * t190 - t141 * t191 + t199 * t94 + t200 * t93 + qJDD(1) + (t111 * t279 + t112 * t278) * qJD(4);
t222 = rSges(6,1) * t285 + t457;
t175 = t222 * t278;
t176 = t222 * t279;
t277 = t286 * rSges(6,1);
t223 = -rSges(6,2) * t285 + t277;
t53 = t152 * t199 + t153 * t200 + qJD(1) + (-t140 * t278 + t141 * t279) * qJD(4);
t327 = -t200 * t222 - t357;
t312 = t327 - t389;
t383 = t140 - t152 - t205;
t61 = t292 * t383 + t312;
t415 = t141 + t153;
t382 = t206 + t415;
t283 = cos(t290);
t388 = t283 * t456;
t412 = t199 * t222 + t238;
t62 = t292 * t382 + t388 - t412;
t494 = -(t175 * t292 - t200 * t223) * t61 - t53 * (-t199 * t175 - t176 * t200) - t62 * (-t292 * t176 - t199 * t223) + t10 * (t278 * t152 + t279 * t153);
t421 = t292 * t294;
t384 = t279 * t421;
t395 = qJD(4) * t295;
t387 = rSges(5,2) * t395;
t380 = rSges(5,1) * t378 + rSges(5,2) * t384 + t278 * t387;
t423 = t279 * t295;
t485 = rSges(5,1) * t423 + t278 * rSges(5,3);
t110 = t292 * t485 - t380;
t179 = t206 * t292;
t459 = rSges(5,1) * t295;
t257 = -rSges(5,2) * t294 + t459;
t235 = t257 * qJD(4);
t256 = rSges(5,1) * t294 + rSges(5,2) * t295;
t289 = qJDD(2) + qJDD(3);
t298 = qJD(2) ^ 2;
t329 = (-qJDD(2) * t282 - t283 * t298) * pkin(2);
t398 = qJD(4) * t279;
t431 = t278 * t294;
t402 = rSges(5,2) * t431 + t279 * rSges(5,3);
t430 = t278 * t295;
t164 = rSges(5,1) * t430 - t402;
t407 = -t164 - t205;
t47 = -t235 * t398 + t191 * t256 + (-t110 - t179) * t292 + t407 * t289 + t329;
t493 = -g(1) + t47;
t196 = t223 * t291;
t420 = t295 * qJD(4) ^ 2;
t21 = t132 * t222 - t196 * t200 + (t191 * t294 - t279 * t420) * pkin(4) + t329 + (-t112 - t179 - t94) * t292 + t383 * t289;
t492 = t21 - g(1);
t275 = pkin(2) * t283;
t461 = pkin(2) * t282;
t356 = qJDD(2) * t275 - t298 * t461;
t337 = t292 * (-pkin(3) * t432 + t244) + t289 * t206 + t356;
t22 = -t131 * t222 - t196 * t199 + (t111 + t93) * t292 + t415 * t289 + (-t190 * t294 - t278 * t420) * pkin(4) + t337;
t491 = t22 - g(2);
t336 = rSges(5,2) * t278 * t421 + rSges(5,3) * t425 - t279 * t387;
t109 = (-t292 * t430 - t377) * rSges(5,1) + t336;
t424 = t279 * t294;
t165 = -rSges(5,2) * t424 + t485;
t399 = qJD(4) * t278;
t48 = t109 * t292 + t165 * t289 - t190 * t256 - t235 * t399 + t337;
t490 = t48 - g(2);
t178 = rSges(4,1) * t425 - rSges(4,2) * t432;
t489 = -t178 * t292 - t203 * t289 - g(1) + t329;
t204 = t279 * rSges(4,1) - rSges(4,2) * t278;
t488 = t204 * t289 - t292 * t441 - g(2) + t356;
t379 = t256 * t398;
t326 = -t379 - t389;
t80 = t292 * t407 + t326;
t487 = t292 * t80;
t287 = Icges(5,4) * t295;
t347 = -Icges(5,2) * t294 + t287;
t254 = Icges(5,1) * t294 + t287;
t151 = t292 * t164;
t197 = t292 * t205;
t484 = -rSges(5,1) * t377 + t151 + t197 + t244 + t336;
t216 = Icges(6,5) * t285 + Icges(6,6) * t286;
t318 = Icges(6,3) * t292 - t216 * t291;
t332 = t346 * t279;
t148 = Icges(6,6) * t278 + t332;
t444 = t148 * t285;
t483 = -t217 * t432 + t279 * t318 + t292 * (-t150 * t286 + t444);
t482 = t278 * t318 + (t330 + t345) * t292;
t251 = Icges(5,5) * t295 - Icges(5,6) * t294;
t250 = Icges(5,5) * t294 + Icges(5,6) * t295;
t315 = Icges(5,3) * t292 - qJD(4) * t250;
t450 = Icges(5,4) * t294;
t255 = Icges(5,1) * t295 - t450;
t335 = t255 * t279;
t162 = Icges(5,5) * t278 + t335;
t333 = t347 * t279;
t160 = Icges(5,6) * t278 + t333;
t442 = t160 * t294;
t342 = -t162 * t295 + t442;
t481 = -t251 * t432 + t279 * t315 + t292 * t342;
t331 = t251 * t279;
t247 = Icges(5,4) * t431;
t161 = Icges(5,1) * t430 - Icges(5,5) * t279 - t247;
t159 = Icges(5,4) * t430 - Icges(5,2) * t431 - Icges(5,6) * t279;
t443 = t159 * t294;
t343 = -t161 * t295 + t443;
t480 = t278 * t315 + (t331 + t343) * t292;
t218 = Icges(6,2) * t286 + t449;
t340 = t218 * t285 - t220 * t286;
t479 = t217 * t291 + t292 * t340;
t252 = Icges(5,2) * t295 + t450;
t338 = t252 * t294 - t254 * t295;
t478 = t251 * qJD(4) + t292 * t338;
t157 = Icges(5,5) * t430 - Icges(5,6) * t431 - Icges(5,3) * t279;
t69 = -t279 * t157 - t278 * t343;
t142 = t292 * t152;
t477 = -rSges(6,1) * t386 - t292 * t140 - t280 * t432 + t142 + t197 + t406;
t409 = -Icges(5,2) * t430 + t161 - t247;
t411 = t254 * t278 + t159;
t476 = -t294 * t409 - t295 * t411;
t475 = t199 * (-t218 * t279 + t150) - t200 * (-Icges(6,2) * t433 + t149 - t230) + t292 * t496;
t474 = t131 / 0.2e1;
t473 = t132 / 0.2e1;
t472 = t190 / 0.2e1;
t471 = t191 / 0.2e1;
t470 = -t199 / 0.2e1;
t469 = t199 / 0.2e1;
t468 = -t200 / 0.2e1;
t467 = t200 / 0.2e1;
t466 = t278 / 0.2e1;
t465 = -t279 / 0.2e1;
t464 = t289 / 0.2e1;
t463 = -t292 / 0.2e1;
t462 = t292 / 0.2e1;
t455 = t279 * t80;
t454 = t292 * t61;
t439 = t216 * t279;
t96 = -t278 * t340 - t439;
t453 = t96 * t292;
t436 = t250 * t279;
t116 = -t278 * t338 - t436;
t446 = t116 * t292;
t440 = t216 * t278;
t438 = t218 * t291;
t437 = t250 * t278;
t435 = t251 * t292;
t145 = Icges(6,5) * t433 - Icges(6,6) * t434 - Icges(6,3) * t279;
t419 = -t278 * t145 - t149 * t426;
t417 = -t278 * t157 - t161 * t423;
t158 = Icges(5,3) * t278 + t331;
t416 = t278 * t158 + t162 * t423;
t410 = -t254 * t279 - t160;
t408 = -t252 * t279 + t162;
t126 = t165 + t206;
t401 = -t252 + t255;
t400 = t254 + t347;
t393 = m(2) + m(3) + m(4);
t392 = t278 * t94 + (t142 + t93) * t279;
t376 = t432 / 0.2e1;
t375 = t425 / 0.2e1;
t374 = -pkin(3) - t459;
t373 = -t399 / 0.2e1;
t372 = t399 / 0.2e1;
t371 = -t398 / 0.2e1;
t370 = t398 / 0.2e1;
t325 = -pkin(4) * t294 - t222;
t320 = Icges(6,5) * t292 - t220 * t291;
t368 = -t147 * t291 + t278 * t320 + t292 * t334;
t367 = -t148 * t291 - t221 * t432 + t279 * t320;
t319 = Icges(6,6) * t292 - t438;
t366 = t149 * t291 + t278 * t319 + t292 * t332;
t365 = t150 * t291 + t279 * t319 - t346 * t432;
t137 = t162 * t430;
t364 = t279 * t158 - t137;
t363 = -t145 + t444;
t361 = -t157 + t442;
t360 = t496 * t291;
t359 = t221 * t291 - t438;
t355 = -pkin(4) * t395 - t196;
t122 = t150 * t433;
t354 = t148 * t434 - t122;
t212 = rSges(3,1) * t283 - rSges(3,2) * t282;
t211 = rSges(3,1) * t282 + rSges(3,2) * t283;
t353 = -t278 * t62 - t279 * t61;
t70 = -t160 * t431 - t364;
t352 = t70 * t278 - t69 * t279;
t71 = -t159 * t424 - t417;
t72 = -t160 * t424 + t416;
t351 = t278 * t72 - t279 * t71;
t198 = t256 * t399;
t81 = t126 * t292 - t198 + t388;
t350 = -t278 * t81 - t455;
t82 = t147 * t286 + t149 * t285;
t98 = t159 * t295 + t161 * t294;
t99 = t160 * t295 + t162 * t294;
t341 = t164 * t278 + t165 * t279;
t339 = t252 * t295 + t254 * t294;
t118 = -t152 + t403;
t324 = t199 * t439 - t200 * t440 - t217 * t292;
t323 = -t294 * t408 + t295 * t410;
t125 = t278 * t374 + t273 + t402;
t322 = -rSges(6,3) * t432 - t381 - t404;
t119 = t153 + t358;
t321 = (-t294 * t400 + t295 * t401) * t292;
t317 = Icges(5,5) * t292 - qJD(4) * t254;
t316 = Icges(5,6) * t292 - qJD(4) * t252;
t308 = t145 * t292 - t285 * t366 + t286 * t368;
t13 = t278 * t482 + t308 * t279;
t307 = t146 * t292 - t285 * t365 + t286 * t367;
t14 = t278 * t483 + t307 * t279;
t15 = t308 * t278 - t279 * t482;
t16 = t307 * t278 - t279 * t483;
t309 = (-t220 * t279 - t148) * t199 - (-t220 * t278 - t147) * t200 + (-t218 + t221) * t292;
t299 = -t285 * t475 + t309 * t286;
t63 = -t145 * t279 - t328;
t64 = -t146 * t279 - t354;
t30 = t199 * t64 - t200 * t63 + t453;
t65 = -t147 * t427 - t419;
t66 = -t148 * t427 + t418;
t97 = -t279 * t340 + t440;
t95 = t97 * t292;
t31 = t199 * t66 - t200 * t65 + t95;
t40 = t285 * t368 + t286 * t366;
t41 = t285 * t367 + t286 * t365;
t306 = t216 * t292 - t285 * t360 + t286 * t359;
t45 = t278 * t479 + t306 * t279;
t46 = t306 * t278 - t279 * t479;
t83 = t148 * t286 + t150 * t285;
t313 = (-t13 * t200 + t131 * t66 + t132 * t65 + t14 * t199 + t289 * t97 + t292 * t45) * t466 + (-t278 * t324 + t279 * t299) * t470 + (t278 * t299 + t279 * t324) * t467 + (t131 * t64 + t132 * t63 - t15 * t200 + t16 * t199 + t289 * t96 + t292 * t46) * t465 + (t309 * t285 + t286 * t475) * t463 + t30 * t376 + t31 * t375 + ((t292 * t66 - t13) * t279 + (t292 * t65 + t14) * t278) * t469 + (t278 * t66 - t279 * t65) * t474 + (t278 * t64 - t279 * t63) * t473 + ((t292 * t64 - t15) * t279 + (t292 * t63 + t16) * t278) * t468 + (t278 * t83 - t279 * t82) * t464 + ((t292 * t83 - t40) * t279 + (t292 * t82 + t41) * t278) * t462;
t105 = t279 * t316 - t347 * t432;
t107 = -t255 * t432 + t279 * t317;
t305 = -qJD(4) * t99 - t105 * t294 + t107 * t295 + t158 * t292;
t106 = t278 * t316 + t292 * t333;
t108 = t278 * t317 + t292 * t335;
t304 = -qJD(4) * t98 - t106 * t294 + t108 * t295 + t157 * t292;
t226 = t347 * qJD(4);
t227 = t255 * qJD(4);
t303 = -qJD(4) * t339 - t226 * t294 + t227 * t295 + t250 * t292;
t302 = (t374 * t455 + (t80 * (-rSges(5,3) - pkin(7)) + t81 * t374) * t278) * t292;
t301 = (t62 * (-pkin(4) * t396 - t222 * t291) + (t61 * (-t280 - t277) - t62 * t296) * t292) * t279;
t117 = -t279 * t338 + t437;
t115 = t117 * t292;
t36 = qJD(4) * t352 + t446;
t37 = qJD(4) * t351 + t115;
t51 = -qJD(4) * t343 + t106 * t295 + t108 * t294;
t52 = -qJD(4) * t342 + t105 * t295 + t107 * t294;
t57 = t278 * t478 + t303 * t279;
t58 = t303 * t278 - t279 * t478;
t300 = (t115 + ((t70 - t137 + (t158 + t443) * t279 + t417) * t279 + t416 * t278) * qJD(4)) * t370 + (t95 + (t64 + (t146 + t445) * t279 + t354 + t419) * t200 + (-t279 * t363 - t495 + t63) * t199) * t467 + (t83 + t97) * t474 + (t82 + t96) * t473 + (t117 + t99) * t472 + (t116 + t98) * t471 + (-t453 + (t66 + t495) * t200 + (t278 * t363 - t122 + t65) * t199 + ((t146 + t345) * t199 + t363 * t200) * t279 + t30) * t470 + (t41 + t45) * t469 + (t36 - t446 + ((t279 * t361 - t416 + t72) * t279 + (t278 * t361 + t364 + t71) * t278) * qJD(4)) * t373 + (t57 + t52) * t372 + (-qJD(4) * t338 + t226 * t295 + t227 * t294 + t285 * t359 + t286 * t360) * t292 + (t40 + t46 + t31) * t468 + (t51 + t58 + t37) * t371 + (t218 * t286 + t220 * t285 + Icges(4,3) + t339) * t289;
t189 = t256 * t279;
t188 = t256 * t278;
t167 = t204 * t292 + t388;
t86 = qJD(4) * t341 + qJD(1);
t44 = t164 * t190 - t165 * t191 + qJDD(1) + (t109 * t279 + t110 * t278) * qJD(4);
t20 = t305 * t278 - t279 * t481;
t19 = t304 * t278 - t279 * t480;
t18 = t278 * t481 + t305 * t279;
t17 = t278 * t480 + t304 * t279;
t1 = [m(5) * t44 + m(6) * t10 + t393 * qJDD(1) + (-m(5) - m(6) - t393) * g(3); Icges(3,3) * qJDD(2) + t300 + (t488 * (t204 + t275) + t489 * (-t203 - t461) + (-t178 - t388 + t167) * t166) * m(4) + (g(1) * t211 - g(2) * t212 + (t211 ^ 2 + t212 ^ 2) * qJDD(2)) * m(3) + (t61 * (t322 - t388) + t301 + (-t389 - t312 + t61 + t477) * t62 + t491 * (t119 + t275) + t492 * (t118 - t461)) * m(6) + (t80 * (t380 - t388) + t302 + (-t389 - t326 + t80 + t484) * t81 + t490 * (t126 + t275) + t493 * (t125 - t461)) * m(5); t300 + (t382 * t454 + t301 + (-t327 + t477) * t62 + (t322 - t412) * t61 + t491 * t119 + t492 * t118) * m(6) + (t302 + (t379 + t484) * t81 + (t380 - t198) * t80 + (t487 + t490) * t126 + t493 * t125) * m(5) + (-t166 * t178 - t167 * t441 + (t166 * t292 + t488) * t204 + (t167 * t292 - t489) * t203) * m(4); ((t292 * t72 - t17) * t279 + (t292 * t71 + t18) * t278) * t372 + t36 * t376 + t37 * t375 + (t117 * t289 + t190 * t72 + t191 * t71 + t292 * t57 + (-t17 * t279 + t18 * t278) * qJD(4)) * t466 + (t116 * t289 + t190 * t70 + t191 * t69 + t292 * t58 + (-t19 * t279 + t20 * t278) * qJD(4)) * t465 + ((-t398 * t437 - t435) * t279 + (t321 + (t323 * t278 + (t436 - t476) * t279) * qJD(4)) * t278) * t370 + ((t294 * t401 + t295 * t400) * t292 + ((t278 * t408 - t279 * t409) * t295 + (t278 * t410 + t279 * t411) * t294) * qJD(4)) * t463 + ((t292 * t70 - t19) * t279 + (t292 * t69 + t20) * t278) * t371 + t313 + ((-t399 * t436 + t435) * t278 + (t321 + (-t476 * t279 + (t437 + t323) * t278) * qJD(4)) * t279) * t373 + t351 * t472 + t352 * t471 + (t278 * t99 - t279 * t98) * t464 + ((t292 * t99 - t51) * t279 + (t292 * t98 + t52) * t278) * t462 + (-(-t62 * t384 + (t353 * t295 + t53 * (-t278 ^ 2 - t279 ^ 2) * t294) * qJD(4)) * pkin(4) + t53 * t392 + (t21 * t325 + t61 * t355 + t10 * t141 + t53 * t111 + (-t53 * t140 + t325 * t62) * t292) * t279 + (t22 * t325 + t62 * t355 - t10 * t140 + t53 * t112 + (t61 * t222 - t415 * t53) * t292) * t278 - g(3) * (t223 + t288) - (g(1) * t279 + g(2) * t278) * t325 + t494) * m(6) + (t44 * t341 + t86 * ((t109 + t151) * t279 + (-t165 * t292 + t110) * t278) + t350 * t235 + ((-t292 * t81 - t47) * t279 + (-t48 + t487) * t278) * t256 - (t188 * t80 - t189 * t81) * t292 - (t86 * (-t188 * t278 - t189 * t279) + t350 * t257) * qJD(4) + g(1) * t189 + g(2) * t188 - g(3) * t257) * m(5); t313 + (g(1) * t176 + g(2) * t175 - g(3) * t223 + t53 * (-t153 * t432 + t392) + t353 * t196 + ((-t292 * t62 - t21) * t279 + (-t22 + t454) * t278) * t222 + t494) * m(6);];
tau = t1;
