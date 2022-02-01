% Calculate vector of inverse dynamics joint torques for
% S5RPRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:48:34
% EndTime: 2022-01-20 09:48:54
% DurationCPUTime: 10.94s
% Computational Cost: add. (18888->602), mult. (11804->779), div. (0->0), fcn. (9175->10), ass. (0->349)
t301 = qJ(1) + pkin(9);
t292 = qJ(3) + t301;
t285 = sin(t292);
t278 = t285 * pkin(7);
t286 = cos(t292);
t300 = qJD(1) + qJD(3);
t303 = sin(qJ(4));
t407 = qJD(4) * t303;
t393 = t285 * t407;
t242 = pkin(4) * t393;
t307 = -pkin(8) - pkin(7);
t442 = t285 * t307;
t416 = -t300 * t442 - t242;
t305 = cos(qJ(4));
t296 = t305 * pkin(4);
t287 = t296 + pkin(3);
t473 = pkin(3) - t287;
t114 = (-t286 * t473 - t278) * t300 + t416;
t408 = qJD(4) * t300;
t243 = t285 * t408;
t405 = qJD(5) * t300;
t132 = t285 * t405 + t243 + (-qJDD(4) - qJDD(5)) * t286;
t208 = t286 * pkin(3) + t278;
t181 = t208 * t300;
t193 = -qJDD(4) * t286 + t243;
t302 = qJ(4) + qJ(5);
t294 = cos(t302);
t284 = t294 * rSges(6,1);
t293 = sin(t302);
t225 = -rSges(6,2) * t293 + t284;
t299 = qJD(4) + qJD(5);
t198 = t225 * t299;
t204 = t286 * t299;
t469 = rSges(6,2) * t294;
t224 = rSges(6,1) * t293 + t469;
t298 = qJDD(1) + qJDD(3);
t290 = sin(t301);
t291 = cos(t301);
t309 = qJD(1) ^ 2;
t304 = sin(qJ(1));
t306 = cos(qJ(1));
t340 = (-qJDD(1) * t304 - t306 * t309) * pkin(1);
t317 = (-qJDD(1) * t290 - t291 * t309) * pkin(2) + t340;
t279 = t286 * pkin(7);
t207 = pkin(3) * t285 - t279;
t264 = t286 * t307;
t415 = -t285 * t287 - t264;
t140 = t207 + t415;
t447 = t285 * t293;
t235 = rSges(6,2) * t447;
t446 = t285 * t294;
t154 = rSges(6,1) * t446 - t286 * rSges(6,3) - t235;
t397 = t140 - t154 - t207;
t432 = t305 * qJD(4) ^ 2;
t402 = t299 * t469;
t440 = t286 * t293;
t403 = rSges(6,2) * t440;
t435 = t293 * t299;
t395 = -t300 * t403 + (-rSges(6,1) * t435 - t402) * t285;
t439 = t286 * t294;
t501 = rSges(6,1) * t439 + t285 * rSges(6,3);
t96 = t300 * t501 + t395;
t17 = t132 * t224 - t198 * t204 + (t193 * t303 - t286 * t432) * pkin(4) + (-t114 - t181 - t96) * t300 + t397 * t298 + t317;
t512 = -g(1) + t17;
t433 = t300 * t303;
t398 = t286 * t433;
t406 = qJD(4) * t305;
t401 = rSges(5,2) * t406;
t394 = rSges(5,1) * t393 + rSges(5,2) * t398 + t285 * t401;
t436 = t286 * t305;
t500 = rSges(5,1) * t436 + t285 * rSges(5,3);
t112 = t300 * t500 - t394;
t472 = rSges(5,1) * t305;
t262 = -rSges(5,2) * t303 + t472;
t237 = t262 * qJD(4);
t260 = rSges(5,1) * t303 + rSges(5,2) * t305;
t409 = qJD(4) * t286;
t444 = t285 * t303;
t414 = rSges(5,2) * t444 + t286 * rSges(5,3);
t443 = t285 * t305;
t166 = rSges(5,1) * t443 - t414;
t419 = -t166 - t207;
t47 = -t237 * t409 + t193 * t260 + (-t112 - t181) * t300 + t419 * t298 + t317;
t511 = -g(1) + t47;
t445 = t285 * t300;
t245 = rSges(4,2) * t445;
t438 = t286 * t300;
t180 = rSges(4,1) * t438 - t245;
t205 = rSges(4,1) * t285 + rSges(4,2) * t286;
t510 = t180 * t300 + t205 * t298 + g(1) - t317;
t248 = pkin(7) * t438;
t391 = t286 * t407;
t371 = pkin(4) * t391;
t113 = -t371 - t248 + (t285 * t473 - t264) * t300;
t192 = qJDD(4) * t285 + t286 * t408;
t131 = qJDD(5) * t285 + t286 * t405 + t192;
t203 = t285 * t299;
t282 = pkin(2) * t291;
t297 = t306 * pkin(1);
t288 = qJDD(1) * t297;
t474 = pkin(1) * t304;
t369 = -pkin(2) * t290 - t474;
t334 = qJDD(1) * t282 + t309 * t369 + t288;
t324 = t300 * (-pkin(3) * t445 + t248) + t298 * t208 + t334;
t372 = t286 * t287 - t442;
t141 = t372 - t208;
t155 = -t403 + t501;
t427 = t141 + t155;
t400 = t294 * t445;
t418 = rSges(6,3) * t438 + t300 * t235;
t95 = -t286 * t402 + (-t286 * t435 - t400) * rSges(6,1) + t418;
t18 = -t131 * t224 - t198 * t203 + (t113 + t95) * t300 + t427 * t298 + (-t192 * t303 - t285 * t432) * pkin(4) + t324;
t509 = -g(2) + t18;
t347 = rSges(5,2) * t285 * t433 + rSges(5,3) * t438 - t286 * t401;
t111 = (-t300 * t443 - t391) * rSges(5,1) + t347;
t437 = t286 * t303;
t167 = -rSges(5,2) * t437 + t500;
t410 = qJD(4) * t285;
t48 = t111 * t300 + t167 * t298 - t192 * t260 - t237 * t410 + t324;
t508 = -g(2) + t48;
t277 = t286 * rSges(4,1);
t206 = -rSges(4,2) * t285 + t277;
t434 = t300 * t205;
t507 = t206 * t298 - t300 * t434 - g(2) + t334;
t283 = Icges(6,4) * t294;
t222 = Icges(6,1) * t293 + t283;
t358 = -Icges(6,2) * t293 + t283;
t506 = t222 + t358;
t232 = Icges(6,4) * t447;
t151 = Icges(6,1) * t446 - Icges(6,5) * t286 - t232;
t149 = Icges(6,4) * t446 - Icges(6,2) * t447 - Icges(6,6) * t286;
t457 = t149 * t293;
t357 = -t151 * t294 + t457;
t339 = t357 * t285;
t219 = Icges(6,5) * t294 - Icges(6,6) * t293;
t341 = t219 * t286;
t148 = Icges(6,3) * t285 + t341;
t462 = Icges(6,4) * t293;
t223 = Icges(6,1) * t294 - t462;
t345 = t223 * t286;
t152 = Icges(6,5) * t285 + t345;
t430 = t285 * t148 + t152 * t439;
t505 = -t339 - t430;
t10 = t131 * t154 - t132 * t155 - t140 * t192 - t141 * t193 + t203 * t96 + t204 * t95 + qJDD(2) + (t113 * t286 + t114 * t285) * qJD(4);
t177 = t224 * t285;
t178 = t224 * t286;
t53 = t154 * t203 + t155 * t204 + qJD(2) + (-t140 * t285 + t141 * t286) * qJD(4);
t338 = -t204 * t224 - t371;
t349 = t369 * qJD(1);
t319 = t349 + t338;
t60 = t300 * t397 + t319;
t499 = t282 + t297;
t348 = t499 * qJD(1);
t396 = t208 + t427;
t424 = t203 * t224 + t242;
t61 = t300 * t396 + t348 - t424;
t504 = -t60 * (t177 * t300 - t204 * t225) - t53 * (-t203 * t177 - t178 * t204) - t61 * (-t300 * t178 - t203 * t225) + t10 * (t285 * t154 + t286 * t155);
t392 = t260 * t409;
t325 = t349 - t392;
t78 = t300 * t419 + t325;
t503 = t300 * t78;
t153 = t300 * t166;
t199 = t300 * t207;
t502 = -t153 - t199;
t214 = t291 * rSges(3,1) - rSges(3,2) * t290;
t201 = t214 + t297;
t295 = Icges(5,4) * t305;
t359 = -Icges(5,2) * t303 + t295;
t258 = Icges(5,1) * t303 + t295;
t142 = t300 * t154;
t497 = t300 * t140 - t142 - t199;
t218 = Icges(6,5) * t293 + Icges(6,6) * t294;
t330 = Icges(6,3) * t300 - t218 * t299;
t343 = t358 * t286;
t150 = Icges(6,6) * t285 + t343;
t456 = t150 * t293;
t496 = -t219 * t445 + t286 * t330 + t300 * (-t152 * t294 + t456);
t495 = t285 * t330 + (t341 + t357) * t300;
t255 = Icges(5,5) * t305 - Icges(5,6) * t303;
t254 = Icges(5,5) * t303 + Icges(5,6) * t305;
t327 = Icges(5,3) * t300 - qJD(4) * t254;
t463 = Icges(5,4) * t303;
t259 = Icges(5,1) * t305 - t463;
t346 = t259 * t286;
t164 = Icges(5,5) * t285 + t346;
t344 = t359 * t286;
t162 = Icges(5,6) * t285 + t344;
t454 = t162 * t303;
t354 = -t164 * t305 + t454;
t494 = -t255 * t445 + t286 * t327 + t300 * t354;
t342 = t255 * t286;
t251 = Icges(5,4) * t444;
t163 = Icges(5,1) * t443 - Icges(5,5) * t286 - t251;
t161 = Icges(5,4) * t443 - Icges(5,2) * t444 - Icges(5,6) * t286;
t455 = t161 * t303;
t355 = -t163 * t305 + t455;
t493 = t285 * t327 + (t342 + t355) * t300;
t220 = Icges(6,2) * t294 + t462;
t352 = t220 * t293 - t222 * t294;
t492 = t219 * t299 + t300 * t352;
t256 = Icges(5,2) * t305 + t463;
t350 = t256 * t303 - t258 * t305;
t491 = t255 * qJD(4) + t300 * t350;
t159 = Icges(5,5) * t443 - Icges(5,6) * t444 - Icges(5,3) * t286;
t69 = -t286 * t159 - t285 * t355;
t421 = -Icges(5,2) * t443 + t163 - t251;
t423 = t258 * t285 + t161;
t490 = -t303 * t421 - t423 * t305;
t489 = t203 * (-t220 * t286 + t152) - t204 * (-Icges(6,2) * t446 + t151 - t232) + t300 * t506;
t488 = m(3) + m(4);
t487 = t131 / 0.2e1;
t486 = t132 / 0.2e1;
t485 = t192 / 0.2e1;
t484 = t193 / 0.2e1;
t483 = -t203 / 0.2e1;
t482 = t203 / 0.2e1;
t481 = -t204 / 0.2e1;
t480 = t204 / 0.2e1;
t479 = t285 / 0.2e1;
t478 = -t286 / 0.2e1;
t477 = t298 / 0.2e1;
t476 = -t300 / 0.2e1;
t475 = t300 / 0.2e1;
t468 = t286 * t78;
t467 = t300 * t60;
t452 = t218 * t286;
t98 = -t285 * t352 - t452;
t466 = t98 * t300;
t449 = t254 * t286;
t116 = -t285 * t350 - t449;
t459 = t116 * t300;
t146 = t206 * t300 + t348;
t458 = t146 * t205;
t453 = t218 * t285;
t451 = t220 * t299;
t450 = t254 * t285;
t448 = t255 * t300;
t147 = Icges(6,5) * t446 - Icges(6,6) * t447 - Icges(6,3) * t286;
t431 = -t285 * t147 - t151 * t439;
t429 = -t285 * t159 - t163 * t436;
t160 = Icges(5,3) * t285 + t342;
t428 = t285 * t160 + t164 * t436;
t422 = -t258 * t286 - t162;
t420 = -t256 * t286 + t164;
t126 = t167 + t208;
t413 = -t256 + t259;
t412 = t258 + t359;
t404 = t285 * t96 + (t142 + t95) * t286;
t390 = t445 / 0.2e1;
t389 = t438 / 0.2e1;
t388 = -pkin(3) - t472;
t387 = -t410 / 0.2e1;
t386 = t410 / 0.2e1;
t385 = -t409 / 0.2e1;
t384 = t409 / 0.2e1;
t337 = -pkin(4) * t303 - t224;
t332 = Icges(6,5) * t300 - t222 * t299;
t382 = -t149 * t299 + t285 * t332 + t300 * t345;
t381 = -t150 * t299 - t223 * t445 + t286 * t332;
t331 = Icges(6,6) * t300 - t451;
t380 = t151 * t299 + t285 * t331 + t300 * t343;
t379 = t152 * t299 + t286 * t331 - t358 * t445;
t137 = t164 * t443;
t378 = t286 * t160 - t137;
t377 = -t147 + t456;
t375 = -t159 + t454;
t374 = t506 * t299;
t373 = t223 * t299 - t451;
t370 = -pkin(4) * t406 - t198;
t122 = t152 * t446;
t367 = t150 * t447 - t122;
t263 = rSges(2,1) * t306 - rSges(2,2) * t304;
t261 = rSges(2,1) * t304 + rSges(2,2) * t306;
t213 = rSges(3,1) * t290 + rSges(3,2) * t291;
t365 = -t285 * t61 - t286 * t60;
t70 = -t162 * t444 - t378;
t364 = t285 * t70 - t286 * t69;
t71 = -t161 * t437 - t429;
t72 = -t162 * t437 + t428;
t363 = t285 * t72 - t286 * t71;
t202 = t260 * t410;
t79 = t126 * t300 - t202 + t348;
t362 = -t285 * t79 - t468;
t82 = t149 * t294 + t151 * t293;
t100 = t161 * t305 + t163 * t303;
t101 = t162 * t305 + t164 * t303;
t353 = t166 * t285 + t167 * t286;
t351 = t256 * t305 + t258 * t303;
t120 = -t154 + t415;
t336 = t203 * t452 - t204 * t453 - t219 * t300;
t335 = -t303 * t420 + t305 * t422;
t125 = t285 * t388 + t279 + t414;
t121 = t155 + t372;
t145 = t349 - t434;
t333 = (-t303 * t412 + t305 * t413) * t300;
t329 = Icges(5,5) * t300 - qJD(4) * t258;
t328 = Icges(5,6) * t300 - qJD(4) * t256;
t321 = t147 * t300 - t293 * t380 + t294 * t382;
t13 = t285 * t495 + t321 * t286;
t320 = t148 * t300 - t293 * t379 + t294 * t381;
t14 = t285 * t496 + t320 * t286;
t15 = t321 * t285 - t286 * t495;
t16 = t320 * t285 - t286 * t496;
t63 = -t147 * t286 - t339;
t64 = -t148 * t286 - t367;
t30 = t203 * t64 - t204 * t63 + t466;
t65 = -t149 * t440 - t431;
t66 = -t150 * t440 + t430;
t99 = -t286 * t352 + t453;
t97 = t99 * t300;
t31 = t203 * t66 - t204 * t65 + t97;
t322 = (-t222 * t286 - t150) * t203 - (-t222 * t285 - t149) * t204 + (-t220 + t223) * t300;
t312 = -t293 * t489 + t322 * t294;
t40 = t293 * t382 + t294 * t380;
t41 = t293 * t381 + t294 * t379;
t318 = t218 * t300 - t293 * t374 + t294 * t373;
t45 = t285 * t492 + t318 * t286;
t46 = t318 * t285 - t286 * t492;
t83 = t150 * t294 + t152 * t293;
t326 = (-t13 * t204 + t131 * t66 + t132 * t65 + t14 * t203 + t298 * t99 + t300 * t45) * t479 + (-t285 * t336 + t286 * t312) * t483 + (t285 * t312 + t286 * t336) * t480 + (t131 * t64 + t132 * t63 - t15 * t204 + t16 * t203 + t298 * t98 + t300 * t46) * t478 + (t322 * t293 + t294 * t489) * t476 + t30 * t390 + t31 * t389 + ((t300 * t66 - t13) * t286 + (t300 * t65 + t14) * t285) * t482 + (t285 * t66 - t286 * t65) * t487 + (t285 * t64 - t286 * t63) * t486 + ((t300 * t64 - t15) * t286 + (t300 * t63 + t16) * t285) * t481 + (t285 * t83 - t286 * t82) * t477 + ((t300 * t83 - t40) * t286 + (t300 * t82 + t41) * t285) * t475;
t107 = t286 * t328 - t359 * t445;
t109 = -t259 * t445 + t286 * t329;
t316 = -qJD(4) * t101 - t107 * t303 + t109 * t305 + t160 * t300;
t108 = t285 * t328 + t300 * t344;
t110 = t285 * t329 + t300 * t346;
t315 = -qJD(4) * t100 - t108 * t303 + t110 * t305 + t159 * t300;
t228 = t359 * qJD(4);
t229 = t259 * qJD(4);
t314 = -qJD(4) * t351 - t228 * t303 + t229 * t305 + t254 * t300;
t117 = -t286 * t350 + t450;
t115 = t117 * t300;
t36 = qJD(4) * t364 + t459;
t37 = qJD(4) * t363 + t115;
t51 = -qJD(4) * t355 + t108 * t305 + t110 * t303;
t52 = -qJD(4) * t354 + t107 * t305 + t109 * t303;
t57 = t285 * t491 + t314 * t286;
t58 = t314 * t285 - t286 * t491;
t313 = (t115 + ((t70 - t137 + (t160 + t455) * t286 + t429) * t286 + t428 * t285) * qJD(4)) * t384 + (t97 + (t64 + (t148 + t457) * t286 + t367 + t431) * t204 + (-t286 * t377 - t505 + t63) * t203) * t480 + (t83 + t99) * t487 + (t82 + t98) * t486 + (t117 + t101) * t485 + (t116 + t100) * t484 + (t30 - t466 + (t66 + t505) * t204 + (t377 * t285 - t122 + t65) * t203 + ((t148 + t357) * t203 + t377 * t204) * t286) * t483 + (t41 + t45) * t482 + (-t459 + ((t286 * t375 - t428 + t72) * t286 + (t285 * t375 + t378 + t71) * t285) * qJD(4) + t36) * t387 + (t52 + t57) * t386 + (-qJD(4) * t350 + t228 * t305 + t229 * t303 + t293 * t373 + t294 * t374) * t300 + (t40 + t46 + t31) * t481 + (t51 + t58 + t37) * t385 + (t220 * t294 + t222 * t293 + Icges(4,3) + t351) * t298;
t311 = t78 * t394 + t79 * (-rSges(5,1) * t391 + t248 + t347) + (t388 * t468 + (t78 * (-rSges(5,3) - pkin(7)) + t79 * t388) * t285) * t300;
t310 = t60 * (-rSges(6,3) * t445 - t395 - t416) + t61 * (-rSges(6,1) * t400 - t287 * t445 + t418) + (t61 * (-pkin(4) * t407 - t224 * t299) + (t60 * (-t287 - t284) - t61 * t307) * t300) * t286;
t191 = t260 * t286;
t190 = t260 * t285;
t86 = qJD(4) * t353 + qJD(2);
t44 = t166 * t192 - t167 * t193 + qJDD(2) + (t111 * t286 + t112 * t285) * qJD(4);
t22 = t316 * t285 - t286 * t494;
t21 = t315 * t285 - t286 * t493;
t20 = t285 * t494 + t316 * t286;
t19 = t285 * t493 + t315 * t286;
t1 = [t313 + (Icges(2,3) + Icges(3,3)) * qJDD(1) + (t145 * t245 + (-t145 * t277 - t458) * t300 + (-t145 * t499 + t146 * t369) * qJD(1) + t507 * (t206 + t499) - t510 * (-t205 + t369)) * m(4) + ((qJDD(1) * t214 - g(2) + t288) * t201 + (-g(1) - qJDD(1) * t213 + t340 + (-0.2e1 * t214 + 0.2e1 * t201 - t297) * t309) * (-t213 - t474)) * m(3) + ((t261 ^ 2 + t263 ^ 2) * qJDD(1) + g(1) * t261 - g(2) * t263) * m(2) + (-(-t60 + t319 + t497) * t61 + (t369 * t61 - t499 * t60) * qJD(1) + t310 + t509 * (t121 + t499) + t512 * (t120 + t369)) * m(6) + (-(-t78 + t325 + t502) * t79 + (t369 * t79 - t499 * t78) * qJD(1) + t311 + t508 * (t126 + t499) + t511 * (t125 + t369)) * m(5); t488 * qJDD(2) + m(5) * t44 + m(6) * t10 + (-m(5) - m(6) - t488) * g(3); t313 + (t310 - t60 * t424 - t61 * (t338 + t497) + t396 * t467 + t509 * t121 + t512 * t120) * m(6) + (t311 - t78 * t202 - t79 * (-t392 + t502) + (t503 + t508) * t126 + t511 * t125) * m(5) + (-t145 * t180 - t146 * t434 + t458 * t300 + (t145 * t300 + t507) * t206 + t510 * t205) * m(4); (t117 * t298 + t192 * t72 + t193 * t71 + t300 * t57 + (-t19 * t286 + t20 * t285) * qJD(4)) * t479 + (t116 * t298 + t192 * t70 + t193 * t69 + t300 * t58 + (-t21 * t286 + t22 * t285) * qJD(4)) * t478 + ((t101 * t300 - t51) * t286 + (t100 * t300 + t52) * t285) * t475 + t326 + ((-t410 * t449 + t448) * t285 + (t333 + (-t490 * t286 + (t450 + t335) * t285) * qJD(4)) * t286) * t387 + ((t303 * t413 + t305 * t412) * t300 + ((t285 * t420 - t286 * t421) * t305 + (t285 * t422 + t286 * t423) * t303) * qJD(4)) * t476 + ((t300 * t70 - t21) * t286 + (t300 * t69 + t22) * t285) * t385 + ((t300 * t72 - t19) * t286 + (t300 * t71 + t20) * t285) * t386 + ((-t409 * t450 - t448) * t286 + (t333 + (t335 * t285 + (t449 - t490) * t286) * qJD(4)) * t285) * t384 + t363 * t485 + t364 * t484 + (-t100 * t286 + t101 * t285) * t477 + t36 * t390 + t37 * t389 + (-g(3) * (t225 + t296) - (g(1) * t286 + g(2) * t285) * t337 - (-t61 * t398 + (t365 * t305 + t53 * (-t285 ^ 2 - t286 ^ 2) * t303) * qJD(4)) * pkin(4) + t53 * t404 + (t17 * t337 + t60 * t370 + t10 * t141 + t53 * t113 + (-t53 * t140 + t337 * t61) * t300) * t286 + (t18 * t337 + t61 * t370 - t10 * t140 + t53 * t114 + (t60 * t224 - t427 * t53) * t300) * t285 + t504) * m(6) + (-(t190 * t78 - t191 * t79) * t300 - (t86 * (-t190 * t285 - t191 * t286) + t362 * t262) * qJD(4) + t44 * t353 + t86 * ((t111 + t153) * t286 + (-t167 * t300 + t112) * t285) + t362 * t237 + ((-t300 * t79 - t47) * t286 + (-t48 + t503) * t285) * t260 + g(1) * t191 + g(2) * t190 - g(3) * t262) * m(5); t326 + (t53 * (-t155 * t445 + t404) + t365 * t198 + ((-t300 * t61 - t17) * t286 + (-t18 + t467) * t285) * t224 + g(1) * t178 + g(2) * t177 - g(3) * t225 + t504) * m(6);];
tau = t1;
