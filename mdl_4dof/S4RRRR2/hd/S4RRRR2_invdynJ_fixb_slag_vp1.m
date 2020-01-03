% Calculate vector of inverse dynamics joint torques for
% S4RRRR2
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRRR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR2_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR2_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR2_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:07
% EndTime: 2019-12-31 17:23:17
% DurationCPUTime: 7.33s
% Computational Cost: add. (13297->575), mult. (11320->748), div. (0->0), fcn. (8923->8), ass. (0->334)
t288 = qJ(3) + qJ(4);
t280 = cos(t288);
t265 = Icges(5,4) * t280;
t278 = sin(t288);
t208 = Icges(5,1) * t278 + t265;
t344 = -Icges(5,2) * t278 + t265;
t492 = t208 + t344;
t291 = sin(qJ(1));
t452 = pkin(1) * qJD(1);
t387 = t291 * t452;
t289 = qJ(1) + qJ(2);
t279 = sin(t289);
t281 = cos(t289);
t211 = rSges(3,1) * t279 + rSges(3,2) * t281;
t287 = qJD(1) + qJD(2);
t433 = t211 * t287;
t166 = -t387 - t433;
t421 = t281 * t287;
t247 = pkin(6) * t421;
t294 = -pkin(7) - pkin(6);
t257 = t281 * t294;
t290 = sin(qJ(3));
t393 = qJD(3) * t290;
t375 = t281 * t393;
t354 = pkin(3) * t375;
t292 = cos(qJ(3));
t283 = t292 * pkin(3);
t275 = t283 + pkin(2);
t457 = pkin(2) - t275;
t108 = -t354 - t247 + (t279 * t457 - t257) * t287;
t272 = t279 * pkin(6);
t376 = t279 * t393;
t235 = pkin(3) * t376;
t423 = t279 * t294;
t401 = -t287 * t423 - t235;
t109 = (-t281 * t457 - t272) * t287 + t401;
t394 = qJD(3) * t287;
t191 = qJDD(3) * t279 + t281 * t394;
t391 = qJD(4) * t287;
t130 = qJDD(4) * t279 + t281 * t391 + t191;
t236 = t279 * t394;
t131 = t279 * t391 + t236 + (-qJDD(3) - qJDD(4)) * t281;
t273 = t281 * pkin(6);
t214 = pkin(2) * t279 - t273;
t400 = -t279 * t275 - t257;
t140 = t214 + t400;
t215 = t281 * pkin(2) + t272;
t355 = t281 * t275 - t423;
t141 = t355 - t215;
t429 = t278 * t279;
t231 = rSges(5,2) * t429;
t427 = t279 * t280;
t152 = rSges(5,1) * t427 - t281 * rSges(5,3) - t231;
t428 = t278 * t281;
t389 = rSges(5,2) * t428;
t422 = t280 * t281;
t483 = rSges(5,1) * t422 + t279 * rSges(5,3);
t153 = -t389 + t483;
t192 = -qJDD(3) * t281 + t236;
t286 = qJD(3) + qJD(4);
t202 = t279 * t286;
t203 = t281 * t286;
t426 = t279 * t287;
t384 = t280 * t426;
t453 = rSges(5,2) * t280;
t388 = t286 * t453;
t403 = rSges(5,3) * t421 + t287 * t231;
t92 = -t281 * t388 + (-t203 * t278 - t384) * rSges(5,1) + t403;
t455 = rSges(5,1) * t278;
t379 = -t202 * t455 - t279 * t388 - t287 * t389;
t93 = t287 * t483 + t379;
t10 = t130 * t152 - t131 * t153 - t140 * t191 - t141 * t192 + t202 * t93 + t203 * t92 + (t108 * t281 + t109 * t279) * qJD(3);
t210 = t453 + t455;
t174 = t210 * t279;
t175 = t210 * t281;
t268 = t280 * rSges(5,1);
t212 = -rSges(5,2) * t278 + t268;
t54 = t152 * t202 + t153 * t203 + (-t140 * t279 + t141 * t281) * qJD(3);
t325 = -t203 * t210 - t354;
t310 = t325 - t387;
t381 = t140 - t152 - t214;
t59 = t287 * t381 + t310;
t412 = t141 + t153;
t380 = t215 + t412;
t293 = cos(qJ(1));
t386 = t293 * t452;
t409 = t202 * t210 + t235;
t60 = t287 * t380 + t386 - t409;
t491 = -(t174 * t287 - t203 * t212) * t59 - t54 * (-t202 * t174 - t175 * t203) - t60 * (-t287 * t175 - t202 * t212) + t10 * (t279 * t152 + t281 * t153);
t179 = t212 * t286;
t182 = t215 * t287;
t285 = qJDD(1) + qJDD(2);
t296 = qJD(1) ^ 2;
t327 = (-qJDD(1) * t291 - t293 * t296) * pkin(1);
t417 = t292 * qJD(3) ^ 2;
t17 = t131 * t210 - t179 * t203 + (t192 * t290 - t281 * t417) * pkin(3) + t327 + (-t109 - t182 - t93) * t287 + t381 * t285;
t490 = t17 - g(1);
t284 = t293 * pkin(1);
t458 = pkin(1) * t291;
t353 = qJDD(1) * t284 - t296 * t458;
t335 = t287 * (-pkin(2) * t426 + t247) + t285 * t215 + t353;
t18 = -t130 * t210 - t179 * t202 + (t108 + t92) * t287 + t412 * t285 + (-t191 * t290 - t279 * t417) * pkin(3) + t335;
t489 = t18 - g(2);
t418 = t287 * t290;
t382 = t281 * t418;
t392 = qJD(3) * t292;
t385 = rSges(4,2) * t392;
t378 = rSges(4,1) * t376 + rSges(4,2) * t382 + t279 * t385;
t419 = t281 * t292;
t482 = rSges(4,1) * t419 + t279 * rSges(4,3);
t107 = t287 * t482 - t378;
t456 = rSges(4,1) * t292;
t250 = -rSges(4,2) * t290 + t456;
t221 = t250 * qJD(3);
t248 = rSges(4,1) * t290 + rSges(4,2) * t292;
t395 = qJD(3) * t281;
t425 = t279 * t290;
t397 = rSges(4,2) * t425 + t281 * rSges(4,3);
t424 = t279 * t292;
t163 = rSges(4,1) * t424 - t397;
t404 = -t163 - t214;
t46 = -t221 * t395 + t192 * t248 + (-t107 - t182) * t287 + t404 * t285 + t327;
t488 = t46 - g(1);
t334 = rSges(4,2) * t279 * t418 + rSges(4,3) * t421 - t281 * t385;
t106 = (-t287 * t424 - t375) * rSges(4,1) + t334;
t420 = t281 * t290;
t164 = -rSges(4,2) * t420 + t482;
t396 = qJD(3) * t279;
t47 = t106 * t287 + t164 * t285 - t191 * t248 - t221 * t396 + t335;
t487 = t47 - g(2);
t181 = rSges(3,1) * t421 - rSges(3,2) * t426;
t486 = -t181 * t287 - t211 * t285 - g(1) + t327;
t213 = t281 * rSges(3,1) - rSges(3,2) * t279;
t485 = t213 * t285 - t287 * t433 - g(2) + t353;
t377 = t248 * t395;
t324 = -t377 - t387;
t90 = t287 * t404 + t324;
t484 = t287 * t90;
t282 = Icges(4,4) * t292;
t345 = -Icges(4,2) * t290 + t282;
t242 = Icges(4,1) * t290 + t282;
t151 = t287 * t163;
t196 = t287 * t214;
t481 = -rSges(4,1) * t375 + t151 + t196 + t247 + t334;
t445 = Icges(5,4) * t278;
t209 = Icges(5,1) * t280 - t445;
t332 = t209 * t281;
t150 = Icges(5,5) * t279 + t332;
t205 = Icges(5,5) * t280 - Icges(5,6) * t278;
t204 = Icges(5,5) * t278 + Icges(5,6) * t280;
t316 = Icges(5,3) * t287 - t204 * t286;
t330 = t344 * t281;
t148 = Icges(5,6) * t279 + t330;
t440 = t148 * t278;
t480 = -t205 * t426 + t281 * t316 + t287 * (-t150 * t280 + t440);
t328 = t205 * t281;
t147 = Icges(5,4) * t427 - Icges(5,2) * t429 - Icges(5,6) * t281;
t228 = Icges(5,4) * t429;
t149 = Icges(5,1) * t427 - Icges(5,5) * t281 - t228;
t343 = t147 * t278 - t149 * t280;
t479 = t279 * t316 + (t328 + t343) * t287;
t239 = Icges(4,5) * t292 - Icges(4,6) * t290;
t238 = Icges(4,5) * t290 + Icges(4,6) * t292;
t313 = Icges(4,3) * t287 - qJD(3) * t238;
t446 = Icges(4,4) * t290;
t243 = Icges(4,1) * t292 - t446;
t333 = t243 * t281;
t162 = Icges(4,5) * t279 + t333;
t331 = t345 * t281;
t160 = Icges(4,6) * t279 + t331;
t437 = t160 * t290;
t340 = -t162 * t292 + t437;
t478 = -t239 * t426 + t281 * t313 + t287 * t340;
t329 = t239 * t281;
t254 = Icges(4,4) * t425;
t161 = Icges(4,1) * t424 - Icges(4,5) * t281 - t254;
t159 = Icges(4,4) * t424 - Icges(4,2) * t425 - Icges(4,6) * t281;
t438 = t159 * t290;
t341 = -t161 * t292 + t438;
t477 = t279 * t313 + (t329 + t341) * t287;
t206 = Icges(5,2) * t280 + t445;
t338 = t206 * t278 - t208 * t280;
t476 = t205 * t286 + t287 * t338;
t240 = Icges(4,2) * t292 + t446;
t336 = t240 * t290 - t242 * t292;
t475 = t239 * qJD(3) + t287 * t336;
t157 = Icges(4,5) * t424 - Icges(4,6) * t425 - Icges(4,3) * t281;
t71 = -t157 * t281 - t279 * t341;
t139 = t287 * t152;
t474 = -rSges(5,1) * t384 - t287 * t140 - t275 * t426 + t139 + t196 + t403;
t406 = -Icges(4,2) * t424 + t161 - t254;
t408 = t242 * t279 + t159;
t473 = -t290 * t406 - t292 * t408;
t472 = t202 * (-t206 * t281 + t150) - t203 * (-Icges(5,2) * t427 + t149 - t228) + t287 * t492;
t471 = t130 / 0.2e1;
t470 = t131 / 0.2e1;
t469 = t191 / 0.2e1;
t468 = t192 / 0.2e1;
t467 = -t202 / 0.2e1;
t466 = t202 / 0.2e1;
t465 = -t203 / 0.2e1;
t464 = t203 / 0.2e1;
t463 = t279 / 0.2e1;
t462 = -t281 / 0.2e1;
t461 = t285 / 0.2e1;
t460 = -t287 / 0.2e1;
t459 = t287 / 0.2e1;
t451 = t281 * t90;
t450 = t287 * t59;
t435 = t204 * t281;
t94 = -t279 * t338 - t435;
t449 = t94 * t287;
t431 = t238 * t281;
t115 = -t279 * t336 - t431;
t442 = t115 * t287;
t145 = Icges(5,5) * t427 - Icges(5,6) * t429 - Icges(5,3) * t281;
t441 = t145 * t281;
t436 = t204 * t279;
t434 = t206 * t286;
t432 = t238 * t279;
t430 = t239 * t287;
t416 = -t279 * t145 - t149 * t422;
t146 = Icges(5,3) * t279 + t328;
t415 = t279 * t146 + t150 * t422;
t414 = -t279 * t157 - t161 * t419;
t158 = Icges(4,3) * t279 + t329;
t413 = t279 * t158 + t162 * t419;
t407 = -t242 * t281 - t160;
t405 = -t240 * t281 + t162;
t127 = t164 + t215;
t399 = -t240 + t243;
t398 = t242 + t345;
t390 = t279 * t93 + (t139 + t92) * t281;
t374 = t426 / 0.2e1;
t373 = t421 / 0.2e1;
t372 = -pkin(2) - t456;
t371 = -t396 / 0.2e1;
t370 = t396 / 0.2e1;
t369 = -t395 / 0.2e1;
t368 = t395 / 0.2e1;
t323 = -pkin(3) * t290 - t210;
t318 = Icges(5,5) * t287 - t208 * t286;
t366 = -t147 * t286 + t279 * t318 + t287 * t332;
t365 = -t148 * t286 - t209 * t426 + t281 * t318;
t317 = Icges(5,6) * t287 - t434;
t364 = t149 * t286 + t279 * t317 + t287 * t330;
t363 = t150 * t286 + t281 * t317 - t344 * t426;
t122 = t150 * t427;
t362 = t146 * t281 - t122;
t136 = t162 * t424;
t361 = t158 * t281 - t136;
t360 = -t145 + t440;
t358 = -t157 + t437;
t357 = t492 * t286;
t356 = t209 * t286 - t434;
t352 = -pkin(3) * t392 - t179;
t251 = rSges(2,1) * t293 - rSges(2,2) * t291;
t249 = rSges(2,1) * t291 + rSges(2,2) * t293;
t351 = -t279 * t60 - t281 * t59;
t72 = -t160 * t425 - t361;
t350 = t279 * t72 - t281 * t71;
t73 = -t159 * t420 - t414;
t74 = -t160 * t420 + t413;
t349 = t279 * t74 - t281 * t73;
t197 = t248 * t396;
t91 = t127 * t287 - t197 + t386;
t348 = -t279 * t91 - t451;
t82 = t147 * t280 + t149 * t278;
t98 = t159 * t292 + t161 * t290;
t99 = t160 * t292 + t162 * t290;
t339 = t163 * t279 + t164 * t281;
t337 = t240 * t292 + t242 * t290;
t326 = t343 * t279;
t117 = -t152 + t400;
t322 = t202 * t435 - t203 * t436 - t205 * t287;
t321 = -t290 * t405 + t292 * t407;
t126 = t279 * t372 + t273 + t397;
t320 = -rSges(5,3) * t426 - t379 - t401;
t118 = t153 + t355;
t319 = (-t290 * t398 + t292 * t399) * t287;
t315 = Icges(4,5) * t287 - qJD(3) * t242;
t314 = Icges(4,6) * t287 - qJD(3) * t240;
t306 = t145 * t287 - t278 * t364 + t280 * t366;
t13 = t279 * t479 + t306 * t281;
t305 = t146 * t287 - t278 * t363 + t280 * t365;
t14 = t279 * t480 + t305 * t281;
t15 = t306 * t279 - t281 * t479;
t16 = t305 * t279 - t281 * t480;
t307 = (-t208 * t281 - t148) * t202 - (-t208 * t279 - t147) * t203 + (-t206 + t209) * t287;
t297 = -t278 * t472 + t307 * t280;
t62 = -t326 - t441;
t63 = -t148 * t429 - t362;
t30 = t202 * t63 - t203 * t62 + t449;
t64 = -t147 * t428 - t416;
t65 = -t148 * t428 + t415;
t95 = -t281 * t338 + t436;
t81 = t95 * t287;
t31 = t202 * t65 - t203 * t64 + t81;
t42 = t278 * t366 + t280 * t364;
t43 = t278 * t365 + t280 * t363;
t304 = t204 * t287 - t278 * t357 + t280 * t356;
t44 = t279 * t476 + t304 * t281;
t45 = t304 * t279 - t281 * t476;
t83 = t148 * t280 + t150 * t278;
t311 = (-t13 * t203 + t130 * t65 + t131 * t64 + t14 * t202 + t285 * t95 + t287 * t44) * t463 + (-t279 * t322 + t281 * t297) * t467 + (t279 * t297 + t281 * t322) * t464 + (t130 * t63 + t131 * t62 - t15 * t203 + t16 * t202 + t285 * t94 + t287 * t45) * t462 + (t307 * t278 + t280 * t472) * t460 + t30 * t374 + t31 * t373 + ((t287 * t65 - t13) * t281 + (t287 * t64 + t14) * t279) * t466 + (t279 * t65 - t281 * t64) * t471 + (t279 * t63 - t281 * t62) * t470 + ((t287 * t63 - t15) * t281 + (t287 * t62 + t16) * t279) * t465 + (t279 * t83 - t281 * t82) * t461 + ((t287 * t83 - t42) * t281 + (t287 * t82 + t43) * t279) * t459;
t102 = t281 * t314 - t345 * t426;
t104 = -t243 * t426 + t281 * t315;
t303 = -qJD(3) * t99 - t102 * t290 + t104 * t292 + t158 * t287;
t103 = t279 * t314 + t287 * t331;
t105 = t279 * t315 + t287 * t333;
t302 = -qJD(3) * t98 - t103 * t290 + t105 * t292 + t157 * t287;
t217 = t345 * qJD(3);
t218 = t243 * qJD(3);
t301 = -qJD(3) * t337 - t217 * t290 + t218 * t292 + t238 * t287;
t300 = (t372 * t451 + (t90 * (-rSges(4,3) - pkin(6)) + t91 * t372) * t279) * t287;
t299 = (t60 * (-pkin(3) * t393 - t210 * t286) + (t59 * (-t275 - t268) - t60 * t294) * t287) * t281;
t116 = -t281 * t336 + t432;
t110 = t116 * t287;
t36 = qJD(3) * t350 + t442;
t37 = qJD(3) * t349 + t110;
t50 = -qJD(3) * t341 + t103 * t292 + t105 * t290;
t51 = -qJD(3) * t340 + t102 * t292 + t104 * t290;
t55 = t279 * t475 + t301 * t281;
t56 = t301 * t279 - t281 * t475;
t298 = (t81 + (t63 + (t147 * t281 + t148 * t279) * t278 + t362 + t416) * t203 + (-t149 * t427 + t441 + t62 + (t147 * t279 - t148 * t281) * t278 + t415) * t202) * t464 + (t110 + ((t72 - t136 + (t158 + t438) * t281 + t414) * t281 + t413 * t279) * qJD(3)) * t368 + (t83 + t95) * t471 + (t82 + t94) * t470 + (t116 + t99) * t469 + (t115 + t98) * t468 + (-t449 + (t65 - t326 - t415) * t203 + (t279 * t360 - t122 + t64) * t202 + ((t146 + t343) * t202 + t360 * t203) * t281 + t30) * t467 + (t43 + t44) * t466 + (-t442 + ((t281 * t358 - t413 + t74) * t281 + (t279 * t358 + t361 + t73) * t279) * qJD(3) + t36) * t371 + (t51 + t55) * t370 + (-qJD(3) * t336 + t217 * t292 + t218 * t290 + t278 * t356 + t280 * t357) * t287 + (t42 + t45 + t31) * t465 + (t50 + t56 + t37) * t369 + (t206 * t280 + t208 * t278 + Icges(3,3) + t337) * t285;
t190 = t248 * t281;
t189 = t248 * t279;
t167 = t213 * t287 + t386;
t97 = t339 * qJD(3);
t22 = t303 * t279 - t281 * t478;
t21 = t302 * t279 - t281 * t477;
t20 = t279 * t478 + t303 * t281;
t19 = t279 * t477 + t302 * t281;
t1 = [Icges(2,3) * qJDD(1) + t298 + (t485 * (t213 + t284) + t486 * (-t211 - t458) + (-t181 - t386 + t167) * t166) * m(3) + ((t249 ^ 2 + t251 ^ 2) * qJDD(1) + g(1) * t249 - g(2) * t251) * m(2) + (t59 * (t320 - t386) + t299 + (-t387 - t310 + t59 + t474) * t60 + t489 * (t118 + t284) + t490 * (t117 - t458)) * m(5) + (t90 * (t378 - t386) + t300 + (-t387 - t324 + t90 + t481) * t91 + t487 * (t127 + t284) + t488 * (t126 - t458)) * m(4); t298 + (t380 * t450 + t299 + (-t325 + t474) * t60 + (t320 - t409) * t59 + t489 * t118 + t490 * t117) * m(5) + (t300 + (t377 + t481) * t91 + (t378 - t197) * t90 + (t484 + t487) * t127 + t488 * t126) * m(4) + (-t166 * t181 - t167 * t433 + (t166 * t287 + t485) * t213 + (t167 * t287 - t486) * t211) * m(3); ((-t396 * t431 + t430) * t279 + (t319 + (-t473 * t281 + (t432 + t321) * t279) * qJD(3)) * t281) * t371 + (t279 * t99 - t281 * t98) * t461 + t349 * t469 + t350 * t468 + t311 + ((t287 * t99 - t50) * t281 + (t287 * t98 + t51) * t279) * t459 + ((t287 * t74 - t19) * t281 + (t287 * t73 + t20) * t279) * t370 + ((t287 * t72 - t21) * t281 + (t287 * t71 + t22) * t279) * t369 + ((-t395 * t432 - t430) * t281 + (t319 + (t321 * t279 + (t431 - t473) * t281) * qJD(3)) * t279) * t368 + (t115 * t285 + t191 * t72 + t192 * t71 + t287 * t56 + (-t21 * t281 + t22 * t279) * qJD(3)) * t462 + (t116 * t285 + t191 * t74 + t192 * t73 + t287 * t55 + (-t19 * t281 + t20 * t279) * qJD(3)) * t463 + t36 * t374 + t37 * t373 + ((t290 * t399 + t292 * t398) * t287 + ((t279 * t405 - t281 * t406) * t292 + (t279 * t407 + t281 * t408) * t290) * qJD(3)) * t460 + (t54 * t390 + (t17 * t323 + t59 * t352 + t10 * t141 + t54 * t108 + (-t54 * t140 + t323 * t60) * t287) * t281 + (t18 * t323 + t60 * t352 - t10 * t140 + t54 * t109 + (t59 * t210 - t412 * t54) * t287) * t279 - (-t60 * t382 + (t351 * t292 + t54 * (-t279 ^ 2 - t281 ^ 2) * t290) * qJD(3)) * pkin(3) - g(3) * (t212 + t283) - (g(1) * t281 + g(2) * t279) * t323 + t491) * m(5) + (-(t189 * t90 - t190 * t91) * t287 - (t97 * (-t189 * t279 - t190 * t281) + t348 * t250) * qJD(3) + (t163 * t191 - t164 * t192 + (t106 * t281 + t107 * t279) * qJD(3)) * t339 + t97 * ((t106 + t151) * t281 + (-t164 * t287 + t107) * t279) + t348 * t221 + ((-t287 * t91 - t46) * t281 + (-t47 + t484) * t279) * t248 + g(1) * t190 + g(2) * t189 - g(3) * t250) * m(4); t311 + (t54 * (-t153 * t426 + t390) + t351 * t179 + ((-t287 * t60 - t17) * t281 + (-t18 + t450) * t279) * t210 + g(1) * t175 + g(2) * t174 - g(3) * t212 + t491) * m(5);];
tau = t1;
