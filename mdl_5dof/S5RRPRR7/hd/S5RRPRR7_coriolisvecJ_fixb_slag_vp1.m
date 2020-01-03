% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR7_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR7_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR7_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR7_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR7_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:15:28
% EndTime: 2019-12-31 20:15:38
% DurationCPUTime: 6.97s
% Computational Cost: add. (13452->551), mult. (11596->709), div. (0->0), fcn. (8900->8), ass. (0->331)
t288 = sin(qJ(4));
t405 = qJD(4) * t288;
t525 = rSges(5,2) * t405;
t287 = qJ(1) + qJ(2);
t279 = sin(t287);
t281 = cos(t287);
t290 = cos(qJ(4));
t460 = Icges(5,4) * t288;
t354 = Icges(5,2) * t290 + t460;
t158 = -Icges(5,6) * t279 + t281 * t354;
t459 = Icges(5,4) * t290;
t356 = Icges(5,1) * t288 + t459;
t327 = t356 * t281;
t160 = -Icges(5,5) * t279 + t327;
t103 = t158 * t288 - t160 * t290;
t285 = qJD(1) + qJD(2);
t325 = t354 * t285;
t243 = -Icges(5,2) * t288 + t459;
t493 = -Icges(5,6) * t285 + qJD(4) * t243;
t107 = t279 * t493 + t281 * t325;
t245 = Icges(5,1) * t290 - t460;
t492 = -Icges(5,5) * t285 + qJD(4) * t245;
t109 = t279 * t492 + t285 * t327;
t157 = Icges(5,6) * t281 + t279 * t354;
t433 = t279 * t290;
t254 = Icges(5,4) * t433;
t434 = t279 * t288;
t159 = Icges(5,1) * t434 + Icges(5,5) * t281 + t254;
t348 = t157 * t290 + t159 * t288;
t524 = -qJD(4) * t348 + t103 * t285 - t107 * t288 + t109 * t290;
t412 = t243 + t356;
t413 = -t354 + t245;
t523 = (t288 * t412 - t290 * t413) * t285;
t289 = sin(qJ(1));
t464 = pkin(1) * qJD(1);
t396 = t289 * t464;
t212 = rSges(3,1) * t279 + rSges(3,2) * t281;
t440 = t212 * t285;
t163 = -t396 - t440;
t286 = qJ(4) + qJ(5);
t278 = sin(t286);
t280 = cos(t286);
t358 = rSges(6,1) * t278 + rSges(6,2) * t280;
t436 = t279 * t280;
t229 = Icges(6,4) * t436;
t437 = t278 * t279;
t145 = Icges(6,1) * t437 + Icges(6,5) * t281 + t229;
t457 = Icges(6,4) * t280;
t355 = Icges(6,1) * t278 + t457;
t326 = t355 * t281;
t146 = -Icges(6,5) * t279 + t326;
t284 = qJD(4) + qJD(5);
t201 = t284 * t279;
t202 = t281 * t284;
t206 = -Icges(6,2) * t278 + t457;
t458 = Icges(6,4) * t278;
t208 = Icges(6,1) * t280 - t458;
t353 = Icges(6,2) * t280 + t458;
t302 = t201 * (t206 * t281 + t146) - t202 * (-Icges(6,2) * t437 + t145 + t229) - t285 * (-t353 + t208);
t143 = Icges(6,6) * t281 + t279 * t353;
t144 = -Icges(6,6) * t279 + t281 * t353;
t303 = t201 * (-t208 * t281 + t144) - t202 * (-t208 * t279 + t143) - t285 * (t206 + t355);
t522 = t303 * t278 - t302 * t280;
t521 = 0.2e1 * qJD(4);
t520 = rSges(5,2) * t290;
t466 = rSges(6,2) * t278;
t467 = rSges(6,1) * t280;
t213 = -t466 + t467;
t173 = t213 * t279;
t404 = qJD(4) * t290;
t386 = t281 * t404;
t236 = pkin(4) * t386;
t259 = qJD(3) * t281;
t291 = cos(qJ(1));
t395 = t291 * t464;
t362 = -t259 + t395;
t338 = -t236 + t362;
t277 = t281 * pkin(2);
t214 = t279 * qJ(3) + t277;
t276 = t281 * pkin(7);
t376 = t214 + t276;
t272 = t281 * rSges(6,3);
t148 = rSges(6,1) * t437 + rSges(6,2) * t436 + t272;
t292 = -pkin(8) - pkin(7);
t257 = pkin(4) * t434;
t505 = t257 - t276;
t166 = -t281 * t292 + t505;
t424 = -t148 - t166;
t498 = -t202 * t213 + t285 * (t376 - t424);
t60 = t338 + t498;
t519 = t60 * (t285 * t173 + t202 * t358);
t351 = Icges(6,5) * t278 + Icges(6,6) * t280;
t142 = -Icges(6,3) * t279 + t281 * t351;
t64 = -t281 * t142 - t144 * t436 - t146 * t437;
t345 = t206 * t280 + t208 * t278;
t204 = Icges(6,5) * t280 - Icges(6,6) * t278;
t444 = t204 * t281;
t96 = t279 * t345 + t444;
t518 = t201 * t64 + t96 * t285;
t174 = t213 * t281;
t398 = t284 * t467;
t198 = t281 * t398;
t397 = t284 * t466;
t94 = t281 * t397 - t198 + (t279 * t358 + t272) * t285;
t517 = t173 * t201 + t202 * t174 + t281 * t94;
t516 = t213 * t285;
t432 = t281 * t285;
t349 = t144 * t280 + t146 * t278;
t515 = t281 * t349;
t337 = -rSges(4,2) * t281 + t279 * rSges(4,3) + t214;
t514 = t285 * t337;
t469 = pkin(7) * t279;
t513 = t285 * t469;
t346 = t158 * t290 + t160 * t288;
t510 = t346 * t281;
t322 = t351 * t285;
t352 = Icges(5,5) * t288 + Icges(5,6) * t290;
t509 = t352 * t285;
t273 = t281 * rSges(5,3);
t161 = rSges(5,1) * t434 + rSges(5,2) * t433 + t273;
t508 = t161 + t376;
t429 = t285 * t292;
t228 = t281 * t429;
t507 = t198 + t228;
t388 = t279 * t404;
t234 = pkin(4) * t388;
t506 = t201 * t213 + t234;
t261 = t281 * qJ(3);
t210 = pkin(2) * t279 - t261;
t192 = t285 * t210;
t258 = qJD(3) * t279;
t363 = t258 - t396;
t339 = -t192 + t363;
t435 = t279 * t285;
t401 = pkin(7) * t435;
t504 = -t339 + t401 - t396;
t211 = t279 * rSges(4,2) + t281 * rSges(4,3);
t411 = rSges(4,2) * t435 + rSges(4,3) * t432;
t503 = -t285 * t211 + t411;
t441 = t208 * t284;
t502 = -Icges(6,5) * t285 + t441;
t442 = t206 * t284;
t501 = -Icges(6,6) * t285 + t442;
t359 = rSges(5,1) * t288 + t520;
t162 = -t279 * rSges(5,3) + t281 * t359;
t147 = t285 * t162;
t252 = rSges(5,1) * t290 - rSges(5,2) * t288;
t407 = qJD(4) * t279;
t194 = t252 * t407;
t389 = t279 * t405;
t431 = t281 * t288;
t394 = t285 * t431;
t391 = t432 * t520 + (t388 + t394) * rSges(5,1);
t239 = qJ(3) * t432;
t414 = t239 + t258;
t500 = -rSges(5,2) * t389 - t147 - t194 + t391 + t414;
t268 = t279 * rSges(6,3);
t149 = t281 * t358 - t268;
t134 = t285 * t149;
t402 = pkin(4) * t431;
t165 = t402 + (pkin(7) + t292) * t279;
t152 = t285 * t165;
t390 = pkin(4) * t394 + t279 * t429 + t234;
t392 = t279 * t398 + t358 * t432;
t499 = t390 + t392 + t414 - t134 - t152;
t406 = qJD(4) * t281;
t195 = t252 * t406;
t497 = t285 * t508 - t195;
t496 = -Icges(6,3) * t285 + t204 * t284;
t106 = t279 * t325 - t281 * t493;
t108 = -t281 * t492 + t356 * t435;
t156 = -Icges(5,3) * t279 + t281 * t352;
t495 = qJD(4) * t103 + t106 * t290 + t108 * t288 + t156 * t285;
t241 = Icges(5,5) * t290 - Icges(5,6) * t288;
t494 = -Icges(5,3) * t285 + qJD(4) * t241;
t218 = t354 * qJD(4);
t219 = t356 * qJD(4);
t491 = qJD(4) * (t243 * t288 - t245 * t290) + t218 * t290 + t219 * t288 + t241 * t285;
t155 = Icges(5,3) * t281 + t279 * t352;
t347 = t157 * t288 - t159 * t290;
t490 = qJD(4) * t347 - t107 * t290 - t109 * t288 + t155 * t285;
t419 = t243 * t281 + t160;
t421 = -t245 * t281 + t158;
t488 = t288 * t421 - t290 * t419;
t420 = -Icges(5,2) * t434 + t159 + t254;
t422 = -t245 * t279 + t157;
t487 = t288 * t422 - t290 * t420;
t367 = t355 * t284 + t442;
t368 = -t353 * t284 + t441;
t486 = t204 * t285 + t278 * t367 - t280 * t368;
t324 = t353 * t285;
t372 = t146 * t284 - t279 * t324 + t281 * t501;
t374 = t144 * t284 - t281 * t502 + t355 * t435;
t485 = t142 * t285 + t278 * t374 - t280 * t372;
t141 = Icges(6,3) * t281 + t279 * t351;
t373 = t145 * t284 + t279 * t501 + t281 * t324;
t375 = t143 * t284 - t279 * t502 - t285 * t326;
t484 = t141 * t285 + t278 * t375 - t280 * t373;
t483 = t279 ^ 2;
t482 = -pkin(2) - pkin(7);
t175 = t284 * t435;
t481 = -t175 / 0.2e1;
t176 = t285 * t202;
t480 = t176 / 0.2e1;
t479 = -t201 / 0.2e1;
t478 = t201 / 0.2e1;
t477 = -t202 / 0.2e1;
t476 = t202 / 0.2e1;
t475 = t279 / 0.2e1;
t474 = -t281 / 0.2e1;
t473 = t281 / 0.2e1;
t472 = -t285 / 0.2e1;
t471 = t285 / 0.2e1;
t470 = pkin(1) * t289;
t282 = t291 * pkin(1);
t167 = t279 * t204;
t97 = t281 * t345 - t167;
t462 = t97 * t285;
t119 = t390 + t401;
t334 = -rSges(6,3) * t285 - t397;
t95 = t279 * t334 + t392;
t461 = -t119 - t95;
t183 = t279 * t241;
t344 = t243 * t290 + t245 * t288;
t116 = t281 * t344 - t183;
t449 = t116 * t285;
t447 = t142 * t279;
t438 = t241 * t281;
t423 = t149 + t165;
t341 = -t210 + t211;
t410 = t258 - t192;
t409 = t261 - t268;
t408 = qJD(3) * t285;
t403 = -rSges(5,3) + t482;
t294 = qJD(1) ^ 2;
t400 = t294 * t470;
t399 = t294 * t282;
t63 = t281 * t141 + t143 * t436 + t145 * t437;
t68 = t281 * t155 + t157 * t433 + t159 * t434;
t69 = -t281 * t156 - t158 * t433 - t160 * t434;
t385 = -t435 / 0.2e1;
t384 = t432 / 0.2e1;
t383 = -t407 / 0.2e1;
t381 = -t406 / 0.2e1;
t378 = pkin(4) * t290 + t213;
t377 = -t210 - t469;
t216 = t281 * rSges(3,1) - rSges(3,2) * t279;
t370 = t174 * t285 - t201 * t358;
t365 = t281 * t408 - t399;
t182 = rSges(3,1) * t432 - rSges(3,2) * t435;
t361 = t258 + t506;
t78 = t194 + (t162 + t377) * t285 + t363;
t79 = t362 + t497;
t357 = t279 * t78 - t281 * t79;
t350 = t143 * t280 + t145 * t278;
t85 = t144 * t278 - t146 * t280;
t342 = t285 * (-pkin(2) * t435 + t414) + t279 * t408 - t400;
t336 = -rSges(5,1) * t386 + t281 * t525;
t335 = t63 + t447;
t333 = t257 + t148 + t214;
t330 = (t279 * t69 + t281 * t68) * qJD(4);
t138 = t279 * t155;
t70 = -t281 * t348 + t138;
t71 = -t156 * t279 + t510;
t329 = (t279 * t71 + t281 * t70) * qJD(4);
t328 = t259 - t336;
t99 = (-t161 * t279 - t162 * t281) * qJD(4);
t321 = pkin(4) * t288 + t358;
t283 = t285 ^ 2;
t320 = -t276 * t283 + t365;
t317 = t167 * t202 - t201 * t444 - t322;
t314 = t279 * t322 - t281 * t496 - t285 * t349;
t313 = t279 * t496 + t281 * t322 + t285 * t350;
t312 = t279 * t509 - t281 * t494 - t285 * t346;
t311 = t279 * t494 + t281 * t509 + t285 * t348;
t126 = t279 * t141;
t65 = -t281 * t350 + t126;
t310 = -t351 * t284 + t285 * t345;
t309 = -t352 * qJD(4) + t285 * t344;
t308 = -t283 * t469 + t342;
t12 = t313 * t279 + t281 * t484;
t13 = t314 * t279 - t281 * t485;
t14 = -t279 * t484 + t313 * t281;
t15 = t279 * t485 + t314 * t281;
t27 = t202 * t63 + t518;
t66 = -t447 + t515;
t28 = t201 * t66 + t202 * t65 - t462;
t42 = -t278 * t373 - t280 * t375;
t43 = t278 * t372 + t280 * t374;
t44 = t310 * t279 + t281 * t486;
t45 = -t279 * t486 + t310 * t281;
t84 = -t143 * t278 + t145 * t280;
t306 = (t12 * t202 + t13 * t201 - t175 * t65 + t176 * t66 + t285 * t44) * t475 + (t317 * t279 - t281 * t522) * t479 + (t279 * t522 + t317 * t281) * t477 + (t14 * t202 + t15 * t201 - t175 * t63 + t176 * t64 + t285 * t45) * t473 + (t278 * t302 + t280 * t303) * t472 + t27 * t385 + t28 * t384 + ((t285 * t66 + t12) * t281 + (-t285 * t65 + t13) * t279) * t478 + (t279 * t64 + t281 * t63) * t481 + (t279 * t66 + t281 * t65) * t480 + ((t285 * t64 + t14) * t281 + (-t285 * t63 + t15) * t279) * t476 + ((t285 * t85 + t42) * t281 + (-t285 * t84 + t43) * t279) * t471;
t299 = t279 * t482 + t162 + t261;
t117 = t285 * t341 + t363;
t118 = t362 + t514;
t298 = (-t117 * t277 + (t117 * (-rSges(4,3) - qJ(3)) - t118 * pkin(2)) * t279) * t285;
t297 = (t78 * t403 * t281 + (t78 * (-qJ(3) - t359) + t79 * t403) * t279) * t285;
t115 = t279 * t344 + t438;
t114 = t115 * t285;
t34 = t114 + t330;
t35 = t329 - t449;
t49 = qJD(4) * t346 - t106 * t288 + t108 * t290;
t52 = t309 * t279 + t281 * t491;
t53 = -t279 * t491 + t309 * t281;
t296 = -t176 * t97 / 0.2e1 + t85 * t480 + (t114 + ((-t70 + t138 + t69) * t279 + (t71 - t510 + (t156 - t348) * t279 + t68) * t281) * qJD(4)) * t383 + ((t335 + t66 - t515) * t202 + t518) * t479 + (t84 + t96) * t481 + (t28 + t462 + (t349 * t279 - t126 + t64) * t202 + (t335 - t63) * t201 + ((t142 + t350) * t202 - t349 * t201) * t281) * t477 + (t42 + t45) * t476 + (t35 + t449 + (t156 * t483 + (-t138 + t69 + (t156 + t348) * t281) * t281) * qJD(4)) * t381 + (-qJD(4) * t344 + t218 * t288 - t219 * t290 - t278 * t368 - t280 * t367) * t285 + (t43 + t44 + t27) * t478 + (t49 + t52 + t34) * t407 / 0.2e1 + (t53 + t524) * t406 / 0.2e1 + (t281 * t116 + (-t347 + t115) * t279) * qJD(4) * t472;
t180 = t358 * t284;
t293 = qJD(4) ^ 2;
t32 = t293 * t402 + t175 * t213 + t180 * t202 + (t234 - t461) * t285 + t308;
t120 = t285 * t505 - t228 - t236;
t153 = t214 * t285 - t259;
t33 = -t293 * t257 + t176 * t213 - t180 * t201 + (-t120 - t153 - t94 + t236) * t285 + t320;
t59 = -t396 + (t377 + t423) * t285 + t361;
t295 = (t33 * (-pkin(2) + t292) - t60 * t397 + (t59 * (-qJ(3) - t321) + t60 * (-rSges(6,3) - pkin(2))) * t285) * t279 + (t33 * t321 + t59 * (-pkin(2) * t285 + t334) - t32 * t292) * t281;
t249 = rSges(4,2) * t432;
t222 = t359 * qJD(4);
t191 = t252 * t281;
t190 = t252 * t279;
t164 = t216 * t285 + t395;
t137 = -t182 * t285 - t399;
t136 = -t285 * t440 - t400;
t129 = t281 * t149;
t113 = (-rSges(5,3) * t285 - t525) * t279 + t391;
t112 = (t279 * t359 + t273) * t285 + t336;
t82 = (-rSges(4,3) * t435 - t153 + t249) * t285 + t365;
t81 = t285 * t411 + t342;
t56 = -t148 * t201 - t149 * t202 + (-t165 * t281 - t166 * t279) * qJD(4);
t55 = -t222 * t407 + (-t112 - t153 + t195) * t285 + t320;
t54 = t113 * t285 + (t222 * t281 + t252 * t435) * qJD(4) + t308;
t16 = -t148 * t176 + t149 * t175 - t201 * t95 + t202 * t94 + ((-t166 * t285 + t120) * t281 + (-t119 + t152) * t279) * qJD(4);
t1 = [m(3) * (t137 * (-t212 - t470) + t136 * (t216 + t282) + (-t182 - t395 + t164) * t163) + t296 + (t33 * (t409 - t470) + t59 * (-t338 + t507) + t32 * (t282 + t333) + t295 + (t59 + t499 + t504 - t506) * t60) * m(6) + (t55 * (t299 - t470) + t78 * (t328 - t395) + t54 * (t282 + t508) + t297 + (t78 + t500 + t504) * t79) * m(5) + (t82 * (t341 - t470) + t117 * (t249 - t362) + t81 * (t282 + t337) + t298 + (t239 + t363 + t117 - t339 + t503) * t118) * m(4); t296 + (t32 * t333 + t33 * t409 + t295 + (t192 - t361 + t499 + t513) * t60 + (t498 + t507) * t59) * m(6) + (t55 * t299 + t54 * t508 + t297 + (-t410 + t500 + t513) * t79 + (-t259 + t328 + t497) * t78) * m(5) + (t81 * t337 + t82 * t341 + t298 + (-t410 + t414 + t503) * t118 + (t249 + t514) * t117) * m(4) + (-(-t163 * t216 - t164 * t212) * t285 + t136 * t216 - t137 * t212 - t163 * t182 - t164 * t440) * m(3); 0.2e1 * (t32 * t474 + t33 * t475) * m(6) + 0.2e1 * (t474 * t54 + t475 * t55) * m(5) + 0.2e1 * (t474 * t81 + t475 * t82) * m(4); ((-t288 * t413 - t290 * t412) * t285 + ((t279 * t421 - t281 * t422) * t290 + (t279 * t419 - t281 * t420) * t288) * qJD(4)) * t472 + (t524 * t281 + (t285 * t347 + t49) * t279) * t471 + t306 + ((-t407 * t438 - t509) * t279 + (t523 + (t487 * t281 + (t183 - t488) * t279) * qJD(4)) * t281) * t383 + ((t183 * t406 - t509) * t281 + (-t523 + (t488 * t279 + (-t438 - t487) * t281) * qJD(4)) * t279) * t381 + (t285 * t52 + ((t311 * t279 + t281 * t490 + t285 * t71) * t281 + (t312 * t279 - t281 * t495 - t285 * t70) * t279) * t521) * t475 + (t285 * t53 + ((-t279 * t490 + t311 * t281 + t285 * t69) * t281 + (t279 * t495 + t312 * t281 - t285 * t68) * t279) * t521) * t473 + (t330 + t34) * t385 + (t329 + t35) * t384 + (-t59 * (-pkin(4) * t389 + t370) - t519 - t16 * t129 + (-t16 * t165 + t60 * t180 - t32 * t378 + t59 * t516) * t281 + (t33 * t378 + t59 * (-pkin(4) * t405 - t180) + t16 * t424 + t60 * t516) * t279 + (-(-t281 ^ 2 - t483) * pkin(4) * t404 + (t424 * t285 + t120) * t281 + (t423 * t285 + t461) * t279 + t517) * t56) * m(6) + (-(t190 * t79 + t191 * t78) * t285 - (t99 * (-t190 * t279 - t191 * t281) - t357 * t359) * qJD(4) + 0.2e1 * t99 * ((-t161 * t285 + t112) * t281 + (-t113 + t147) * t279) - t357 * t222 + ((t285 * t78 - t54) * t281 + (t285 * t79 + t55) * t279) * t252) * m(5); t306 + (-t519 - t59 * t370 + t16 * (-t148 * t279 - t129) - (t279 * t59 - t281 * t60) * t180 + ((t285 * t59 - t32) * t281 + (t285 * t60 + t33) * t279) * t213 + (-t148 * t432 + (-t95 + t134) * t279 + t517) * t56) * m(6);];
tauc = t1(:);
