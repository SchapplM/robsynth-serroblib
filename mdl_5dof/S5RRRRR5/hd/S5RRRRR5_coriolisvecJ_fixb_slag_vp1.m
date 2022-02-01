% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 12:01:40
% EndTime: 2022-01-20 12:01:58
% DurationCPUTime: 9.19s
% Computational Cost: add. (22097->564), mult. (13212->715), div. (0->0), fcn. (10256->10), ass. (0->351)
t281 = qJ(1) + qJ(2);
t275 = qJ(3) + t281;
t266 = sin(t275);
t267 = cos(t275);
t193 = rSges(4,1) * t266 + rSges(4,2) * t267;
t279 = qJD(1) + qJD(2);
t269 = qJD(3) + t279;
t178 = t269 * t193;
t283 = sin(qJ(1));
t478 = pkin(1) * qJD(1);
t406 = t283 * t478;
t271 = sin(t281);
t441 = t271 * t279;
t413 = pkin(2) * t441;
t128 = -t406 - t413 - t178;
t280 = qJ(4) + qJ(5);
t272 = cos(t280);
t263 = Icges(6,4) * t272;
t270 = sin(t280);
t207 = Icges(6,1) * t270 + t263;
t359 = -Icges(6,2) * t270 + t263;
t519 = t207 + t359;
t282 = sin(qJ(4));
t451 = t266 * t282;
t422 = rSges(5,2) * t451 + t267 * rSges(5,3);
t284 = cos(qJ(4));
t450 = t266 * t284;
t155 = rSges(5,1) * t450 - t422;
t134 = t269 * t155;
t261 = t267 * pkin(8);
t197 = pkin(3) * t266 - t261;
t181 = t269 * t197;
t448 = t267 * t269;
t221 = pkin(8) * t448;
t416 = qJD(4) * t284;
t404 = rSges(5,2) * t416;
t443 = t269 * t282;
t345 = rSges(5,2) * t266 * t443 + rSges(5,3) * t448 - t267 * t404;
t417 = qJD(4) * t282;
t395 = t267 * t417;
t518 = -rSges(5,1) * t395 + t134 + t181 + t221 + t345;
t286 = -pkin(9) - pkin(8);
t248 = t267 * t286;
t484 = pkin(4) * t284;
t268 = pkin(3) + t484;
t424 = -t266 * t268 - t248;
t130 = t197 + t424;
t118 = t269 * t130;
t453 = t266 * t270;
t227 = rSges(6,2) * t453;
t452 = t266 * t272;
t141 = rSges(6,1) * t452 - t267 * rSges(6,3) - t227;
t127 = t269 * t141;
t403 = t269 * t452;
t427 = rSges(6,3) * t448 + t269 * t227;
t454 = t266 * t269;
t517 = -rSges(6,1) * t403 - t268 * t454 - t118 + t127 + t181 + t427;
t224 = Icges(6,4) * t453;
t139 = Icges(6,1) * t452 - Icges(6,5) * t267 - t224;
t137 = Icges(6,4) * t452 - Icges(6,2) * t453 - Icges(6,6) * t267;
t467 = t137 * t270;
t358 = -t139 * t272 + t467;
t337 = t358 * t266;
t204 = Icges(6,5) * t272 - Icges(6,6) * t270;
t338 = t204 * t267;
t136 = Icges(6,3) * t266 + t338;
t471 = Icges(6,4) * t270;
t208 = Icges(6,1) * t272 - t471;
t342 = t208 * t267;
t140 = Icges(6,5) * t266 + t342;
t446 = t267 * t272;
t438 = t266 * t136 + t140 * t446;
t516 = -t337 - t438;
t273 = cos(t281);
t210 = rSges(3,1) * t271 + rSges(3,2) * t273;
t458 = t210 * t279;
t170 = -t406 - t458;
t278 = qJD(4) + qJD(5);
t479 = rSges(6,2) * t272;
t407 = t278 * t479;
t442 = t270 * t278;
t515 = -rSges(6,1) * t442 - t407;
t514 = 0.2e1 * qJD(4);
t423 = rSges(6,1) * t446 + t266 * rSges(6,3);
t444 = t267 * t284;
t512 = rSges(5,1) * t444 + t266 * rSges(5,3);
t260 = t266 * pkin(8);
t198 = t267 * pkin(3) + t260;
t274 = Icges(5,4) * t284;
t360 = -Icges(5,2) * t282 + t274;
t241 = Icges(5,1) * t282 + t274;
t192 = t267 * t278;
t158 = t269 * t192;
t480 = rSges(6,2) * t270;
t481 = rSges(6,1) * t272;
t211 = -t480 + t481;
t185 = t211 * t278;
t191 = t266 * t278;
t209 = rSges(6,1) * t270 + t479;
t288 = qJD(1) ^ 2;
t487 = pkin(1) * t283;
t415 = t288 * t487;
t485 = pkin(2) * t279 ^ 2;
t347 = -t271 * t485 - t415;
t329 = t269 * (-pkin(3) * t454 + t221) + t347;
t373 = pkin(4) * t395;
t411 = qJD(4) ^ 2 * t484;
t86 = -t267 * t407 + (-t267 * t442 - t403) * rSges(6,1) + t427;
t483 = pkin(3) - t268;
t96 = -t373 - t221 + (t266 * t483 - t248) * t269;
t32 = -t266 * t411 - t158 * t209 - t185 * t191 + (t86 + t96 - t373) * t269 + t329;
t333 = -t192 * t209 - t373;
t311 = t333 - t413;
t298 = t311 - t406;
t57 = (t130 - t141 - t197) * t269 + t298;
t285 = cos(qJ(1));
t405 = t285 * t478;
t440 = t273 * t279;
t412 = pkin(2) * t440;
t335 = t405 + t412;
t449 = t266 * t286;
t374 = t267 * t268 - t449;
t131 = t374 - t198;
t447 = t267 * t270;
t408 = rSges(6,2) * t447;
t142 = -t408 + t423;
t435 = -t131 - t142;
t400 = t198 - t435;
t396 = t266 * t417;
t231 = pkin(4) * t396;
t432 = t191 * t209 + t231;
t58 = t269 * t400 + t335 - t432;
t290 = (-t32 * t480 + t58 * (-pkin(4) * t417 + t515) + (t57 * (-t268 - t481) - t58 * t286) * t269) * t267;
t477 = t269 * t57;
t511 = t400 * t477 + t290;
t509 = -t413 + t518;
t238 = Icges(5,5) * t284 - Icges(5,6) * t282;
t237 = Icges(5,5) * t282 + Icges(5,6) * t284;
t314 = Icges(5,3) * t269 - qJD(4) * t237;
t472 = Icges(5,4) * t282;
t242 = Icges(5,1) * t284 - t472;
t343 = t242 * t267;
t154 = Icges(5,5) * t266 + t343;
t341 = t360 * t267;
t152 = Icges(5,6) * t266 + t341;
t463 = t152 * t282;
t355 = -t154 * t284 + t463;
t508 = -t238 * t454 + t267 * t314 + t269 * t355;
t339 = t238 * t267;
t234 = Icges(5,4) * t451;
t153 = Icges(5,1) * t450 - Icges(5,5) * t267 - t234;
t151 = Icges(5,4) * t450 - Icges(5,2) * t451 - Icges(5,6) * t267;
t464 = t151 * t282;
t356 = -t153 * t284 + t464;
t507 = t266 * t314 + (t339 + t356) * t269;
t203 = Icges(6,5) * t270 + Icges(6,6) * t272;
t317 = Icges(6,3) * t269 - t203 * t278;
t340 = t359 * t267;
t138 = Icges(6,6) * t266 + t340;
t466 = t138 * t270;
t506 = -t204 * t454 + t267 * t317 + t269 * (-t140 * t272 + t466);
t505 = t266 * t317 + (t338 + t358) * t269;
t205 = Icges(6,2) * t272 + t471;
t354 = t205 * t270 - t207 * t272;
t504 = t204 * t278 + t269 * t354;
t244 = rSges(5,1) * t282 + rSges(5,2) * t284;
t419 = qJD(4) * t266;
t188 = t244 * t419;
t445 = t267 * t282;
t156 = -rSges(5,2) * t445 + t512;
t330 = t156 + t198;
t363 = -t269 * t330 + t188;
t239 = Icges(5,2) * t284 + t472;
t353 = t239 * t282 - t241 * t284;
t503 = t238 * qJD(4) + t269 * t353;
t149 = Icges(5,5) * t450 - Icges(5,6) * t451 - Icges(5,3) * t267;
t65 = -t149 * t267 - t266 * t356;
t502 = -t413 + t517;
t157 = t269 * t191;
t399 = t266 * t515 - t269 * t408;
t87 = t423 * t269 + t399;
t426 = -t269 * t449 - t231;
t97 = (-t267 * t483 - t260) * t269 + t426;
t12 = t141 * t158 - t142 * t157 + t191 * t87 + t192 * t86 + ((t96 - t118) * t267 + (-t131 * t269 + t97) * t266) * qJD(4);
t168 = t209 * t266;
t169 = t209 * t267;
t50 = t141 * t191 + t142 * t192 + (-t130 * t266 + t131 * t267) * qJD(4);
t501 = -t57 * (t168 * t269 - t192 * t211) - t50 * (-t191 * t168 - t169 * t192) - t58 * (-t269 * t169 - t191 * t211) + t12 * (t266 * t141 + t267 * t142);
t482 = rSges(5,1) * t284;
t392 = -pkin(3) - t482;
t418 = qJD(4) * t267;
t397 = t244 * t418;
t334 = -t397 - t413;
t310 = t334 - t406;
t74 = (-t155 - t197) * t269 + t310;
t476 = t74 * t267;
t75 = t335 - t363;
t292 = (t392 * t476 + (t74 * (-rSges(5,3) - pkin(8)) + t75 * t392) * t266) * t269;
t401 = t267 * t443;
t398 = rSges(5,1) * t396 + rSges(5,2) * t401 + t266 * t404;
t500 = t292 + (-t363 + t398) * t74;
t429 = -Icges(5,2) * t450 + t153 - t234;
t431 = t241 * t266 + t151;
t499 = -t282 * t429 - t284 * t431;
t498 = t191 * (-t205 * t267 + t140) - t192 * (-Icges(6,2) * t452 + t139 - t224) + t269 * t519;
t497 = t157 / 0.2e1;
t496 = t158 / 0.2e1;
t495 = -t191 / 0.2e1;
t494 = t191 / 0.2e1;
t493 = -t192 / 0.2e1;
t492 = t192 / 0.2e1;
t491 = t266 / 0.2e1;
t490 = -t267 / 0.2e1;
t489 = -t269 / 0.2e1;
t488 = t269 / 0.2e1;
t486 = pkin(2) * t271;
t276 = t285 * pkin(1);
t460 = t203 * t267;
t92 = -t266 * t354 - t460;
t475 = t92 * t269;
t456 = t237 * t267;
t110 = -t266 * t353 - t456;
t468 = t110 * t269;
t194 = t267 * rSges(4,1) - rSges(4,2) * t266;
t462 = t194 * t269;
t461 = t203 * t266;
t459 = t205 * t278;
t457 = t237 * t266;
t455 = t238 * t269;
t135 = Icges(6,5) * t452 - Icges(6,6) * t453 - Icges(6,3) * t267;
t439 = -t266 * t135 - t139 * t446;
t437 = -t266 * t149 - t153 * t444;
t150 = Icges(5,3) * t266 + t339;
t436 = t266 * t150 + t154 * t444;
t430 = -t241 * t267 - t152;
t428 = -t239 * t267 + t154;
t421 = -t239 + t242;
t420 = t241 + t360;
t414 = t288 * t276;
t410 = t266 * t87 + (t127 + t86) * t267;
t394 = t454 / 0.2e1;
t393 = t448 / 0.2e1;
t391 = -t419 / 0.2e1;
t388 = t418 / 0.2e1;
t386 = -pkin(4) * t282 - t209;
t212 = t273 * rSges(3,1) - rSges(3,2) * t271;
t319 = Icges(6,5) * t269 - t207 * t278;
t384 = -t137 * t278 + t266 * t319 + t269 * t342;
t383 = -t138 * t278 - t208 * t454 + t267 * t319;
t318 = Icges(6,6) * t269 - t459;
t382 = t139 * t278 + t266 * t318 + t269 * t340;
t381 = t140 * t278 + t267 * t318 - t359 * t454;
t124 = t154 * t450;
t380 = t267 * t150 - t124;
t379 = -t135 + t466;
t377 = -t149 + t463;
t376 = t519 * t278;
t375 = t208 * t278 - t459;
t187 = rSges(3,1) * t440 - rSges(3,2) * t441;
t160 = rSges(4,1) * t448 - rSges(4,2) * t454;
t370 = -pkin(4) * t416 - t185;
t265 = pkin(2) * t273;
t369 = t194 + t265;
t114 = t140 * t452;
t368 = t138 * t453 - t114;
t366 = -rSges(5,2) * t282 + t482;
t365 = -t58 * t266 - t57 * t267;
t364 = -t75 * t266 - t476;
t88 = t137 * t272 + t139 * t270;
t107 = t151 * t284 + t153 * t282;
t108 = t152 * t284 + t154 * t282;
t351 = t374 + t423;
t66 = -t152 * t451 - t380;
t349 = (t266 * t66 - t267 * t65) * qJD(4);
t67 = -t151 * t445 - t437;
t68 = -t152 * t445 + t436;
t348 = (t266 * t68 - t267 * t67) * qJD(4);
t346 = -t273 * t485 - t414;
t344 = t265 + t351;
t91 = (t155 * t266 + t156 * t267) * qJD(4);
t332 = -t193 - t486;
t331 = -t141 + t424;
t328 = t191 * t460 - t192 * t461 - t204 * t269;
t327 = -t160 - t412;
t326 = -t282 * t428 + t284 * t430;
t325 = t266 * t392 + t261 + t422;
t324 = t265 + t330;
t323 = -rSges(6,3) * t454 - t399 - t426;
t300 = t135 * t269 - t270 * t382 + t272 * t384;
t13 = t505 * t266 + t300 * t267;
t299 = t136 * t269 - t270 * t381 + t272 * t383;
t14 = t506 * t266 + t299 * t267;
t15 = t300 * t266 - t505 * t267;
t16 = t299 * t266 - t506 * t267;
t61 = -t135 * t267 - t337;
t62 = -t136 * t267 - t368;
t28 = t191 * t62 - t192 * t61 + t475;
t304 = (-t207 * t267 - t138) * t191 - (-t207 * t266 - t137) * t192 + (-t205 + t208) * t269;
t289 = -t498 * t270 + t304 * t272;
t63 = -t137 * t447 - t439;
t64 = -t138 * t447 + t438;
t93 = -t267 * t354 + t461;
t90 = t93 * t269;
t29 = t191 * t64 - t192 * t63 + t90;
t40 = t270 * t384 + t272 * t382;
t41 = t270 * t383 + t272 * t381;
t297 = t203 * t269 - t270 * t376 + t272 * t375;
t44 = t504 * t266 + t297 * t267;
t45 = t297 * t266 - t504 * t267;
t89 = t138 * t272 + t140 * t270;
t322 = (-t13 * t192 + t14 * t191 + t157 * t63 + t158 * t64 + t269 * t44) * t491 + (-t266 * t328 + t267 * t289) * t495 + (t266 * t289 + t267 * t328) * t492 + (-t15 * t192 + t157 * t61 + t158 * t62 + t16 * t191 + t269 * t45) * t490 + (t304 * t270 + t498 * t272) * t489 + t28 * t394 + t29 * t393 + ((t269 * t64 - t13) * t267 + (t269 * t63 + t14) * t266) * t494 + (t266 * t62 - t267 * t61) * t497 + (t266 * t64 - t267 * t63) * t496 + ((t269 * t62 - t15) * t267 + (t269 * t61 + t16) * t266) * t493 + ((t269 * t89 - t40) * t267 + (t269 * t88 + t41) * t266) * t488;
t321 = t331 - t486;
t320 = (-t282 * t420 + t284 * t421) * t269;
t316 = Icges(5,5) * t269 - qJD(4) * t241;
t315 = Icges(5,6) * t269 - qJD(4) * t239;
t309 = t325 - t486;
t306 = t323 - t412;
t100 = t267 * t315 - t360 * t454;
t102 = -t242 * t454 + t267 * t316;
t296 = -qJD(4) * t108 - t100 * t282 + t102 * t284 + t150 * t269;
t101 = t266 * t315 + t269 * t341;
t103 = t266 * t316 + t269 * t343;
t295 = -qJD(4) * t107 - t101 * t282 + t103 * t284 + t149 * t269;
t214 = t360 * qJD(4);
t215 = t242 * qJD(4);
t294 = -t214 * t282 + t215 * t284 + t237 * t269 + (-t239 * t284 - t241 * t282) * qJD(4);
t111 = -t267 * t353 + t457;
t109 = t111 * t269;
t34 = t349 + t468;
t35 = t109 + t348;
t48 = -qJD(4) * t356 + t101 * t284 + t103 * t282;
t49 = -qJD(4) * t355 + t100 * t284 + t102 * t282;
t53 = t503 * t266 + t294 * t267;
t54 = t294 * t266 - t503 * t267;
t291 = (t90 + (t62 + (t136 + t467) * t267 + t368 + t439) * t192 + (-t267 * t379 - t516 + t61) * t191) * t492 + (t109 + ((t66 - t124 + (t150 + t464) * t267 + t437) * t267 + t436 * t266) * qJD(4)) * t388 + (t88 + t92) * t497 + (t89 + t93) * t496 + (t28 - t475 + (t64 + t516) * t192 + (t379 * t266 - t114 + t63) * t191 + ((t136 + t358) * t191 + t379 * t192) * t267) * t495 + (t41 + t44) * t494 + (t34 - t468 + ((t267 * t377 - t436 + t68) * t267 + (t266 * t377 + t380 + t67) * t266) * qJD(4)) * t391 + (t49 + t53) * t419 / 0.2e1 + (-qJD(4) * t353 + t214 * t284 + t215 * t282 + t270 * t375 + t272 * t376) * t269 + (t40 + t45 + t29) * t493 - (t35 + t48 + t54) * t418 / 0.2e1 + ((t107 + t110) * t266 + (t108 + t111) * t267) * qJD(4) * t488;
t219 = t366 * qJD(4);
t180 = t244 * t267;
t179 = t244 * t266;
t171 = t212 * t279 + t405;
t161 = t198 * t269;
t144 = -t187 * t279 - t414;
t143 = -t279 * t458 - t415;
t129 = t335 + t462;
t113 = -t160 * t269 + t346;
t112 = -t178 * t269 + t347;
t106 = t512 * t269 - t398;
t105 = (-t269 * t450 - t395) * rSges(5,1) + t345;
t56 = -t219 * t418 + (-t106 - t161 + t188) * t269 + t346;
t55 = t105 * t269 + (-t219 * t266 - t244 * t448) * qJD(4) + t329;
t33 = -t267 * t411 + t157 * t209 - t185 * t192 + (-t161 - t87 - t97 + t231) * t269 + t346;
t1 = [t291 + m(3) * (t144 * (-t210 - t487) + t143 * (t212 + t276) + (-t187 - t405 + t171) * t170) + (t33 * (t321 - t487) + t57 * (t306 - t405) + t32 * (t276 + t344) + t290 + (-t298 + t57 - t406 + t502) * t58) * m(6) + (t56 * (t309 - t487) + t74 * (-t335 + t398) + t55 * (t276 + t324) + t292 + (-t310 + t74 - t406 + t509) * t75) * m(5) + m(4) * (t113 * (t332 - t487) + t112 * (t276 + t369) + (t327 - t405 + t129) * t128); t291 + (t32 * t344 + t321 * t33 + (-t311 + t502) * t58 + (t306 + t412 - t432) * t57 + t511) * m(6) + (t56 * t309 + t55 * t324 + (-t334 + t509) * t75 + t500) * m(5) + (t112 * t369 + t113 * t332 + (t327 + t412 + t462) * t128) * m(4) + (-(-t170 * t212 - t171 * t210) * t279 + t143 * t212 - t144 * t210 - t170 * t187 - t171 * t458) * m(3); t291 + (t32 * t351 + t33 * t331 + (-t333 + t517) * t58 + (t323 - t432) * t57 + t511) * m(6) + (t56 * t325 + t55 * t330 + (t397 + t518) * t75 + t500) * m(5) + (-(-t128 * t194 - t129 * t193) * t269 + t112 * t194 - t113 * t193 - t128 * t160 - t129 * t178) * m(4); ((-t418 * t457 - t455) * t267 + (t320 + (t326 * t266 + (t456 - t499) * t267) * qJD(4)) * t266) * t388 + ((-t419 * t456 + t455) * t266 + (t320 + (-t499 * t267 + (t457 + t326) * t266) * qJD(4)) * t267) * t391 + t322 + ((t282 * t421 + t284 * t420) * t269 + ((t266 * t428 - t429 * t267) * t284 + (t266 * t430 + t267 * t431) * t282) * qJD(4)) * t489 + ((t108 * t269 - t48) * t267 + (t107 * t269 + t49) * t266) * t488 + (t269 * t53 + ((-t507 * t266 - t295 * t267 + t269 * t68) * t267 + (t508 * t266 + t296 * t267 + t269 * t67) * t266) * t514) * t491 + (t269 * t54 + ((-t295 * t266 + t507 * t267 + t269 * t66) * t267 + (t296 * t266 - t508 * t267 + t269 * t65) * t266) * t514) * t490 + (t349 + t34) * t394 + (t348 + t35) * t393 + (-(-t58 * t401 + (t365 * t284 + t50 * (-t266 ^ 2 - t267 ^ 2) * t282) * qJD(4)) * pkin(4) + t50 * t410 + (t33 * t386 + t57 * t370 + t12 * t131 + t50 * t96 + (-t50 * t130 + t386 * t58) * t269) * t267 + (t32 * t386 + t58 * t370 - t12 * t130 + t50 * t97 + (t57 * t209 + t435 * t50) * t269) * t266 + t501) * m(6) + (0.2e1 * t91 * ((t105 + t134) * t267 + (-t156 * t269 + t106) * t266) + t364 * t219 + ((-t269 * t75 - t56) * t267 + (t269 * t74 - t55) * t266) * t244 - (t179 * t74 - t180 * t75) * t269 - (t91 * (-t179 * t266 - t180 * t267) + t364 * t366) * qJD(4)) * m(5); t322 + (t50 * (-t142 * t454 + t410) + t365 * t185 + ((-t269 * t58 - t33) * t267 + (-t32 + t477) * t266) * t209 + t501) * m(6);];
tauc = t1(:);
