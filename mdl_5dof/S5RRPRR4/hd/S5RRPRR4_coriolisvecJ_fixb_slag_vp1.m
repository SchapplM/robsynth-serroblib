% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:03
% EndTime: 2022-01-20 10:48:19
% DurationCPUTime: 7.58s
% Computational Cost: add. (17736->548), mult. (11158->709), div. (0->0), fcn. (8628->10), ass. (0->347)
t278 = qJ(4) + qJ(5);
t271 = cos(t278);
t262 = Icges(6,4) * t271;
t269 = sin(t278);
t205 = Icges(6,1) * t269 + t262;
t346 = -Icges(6,2) * t269 + t262;
t509 = t205 + t346;
t279 = qJ(1) + qJ(2);
t268 = pkin(9) + t279;
t265 = cos(t268);
t272 = cos(t279);
t266 = pkin(2) * t272;
t508 = -t265 * rSges(4,1) - t266;
t264 = sin(t268);
t437 = t265 * t271;
t498 = rSges(6,1) * t437 + t264 * rSges(6,3);
t507 = -t266 - t498;
t445 = t264 * t269;
t217 = Icges(6,4) * t445;
t444 = t264 * t271;
t139 = Icges(6,1) * t444 - Icges(6,5) * t265 - t217;
t137 = Icges(6,4) * t444 - Icges(6,2) * t445 - Icges(6,6) * t265;
t456 = t137 * t269;
t345 = -t139 * t271 + t456;
t322 = t345 * t264;
t202 = Icges(6,5) * t271 - Icges(6,6) * t269;
t323 = t202 * t265;
t136 = Icges(6,3) * t264 + t323;
t460 = Icges(6,4) * t269;
t206 = Icges(6,1) * t271 - t460;
t327 = t206 * t265;
t140 = Icges(6,5) * t264 + t327;
t429 = t264 * t136 + t140 * t437;
t506 = -t322 - t429;
t400 = rSges(6,1) * t444;
t284 = -pkin(8) - pkin(7);
t247 = t265 * t284;
t282 = cos(qJ(4));
t470 = pkin(4) * t282;
t267 = pkin(3) + t470;
t415 = -t264 * t267 - t247;
t270 = sin(t279);
t472 = pkin(2) * t270;
t505 = -t400 - t472 + t415;
t281 = sin(qJ(1));
t465 = pkin(1) * qJD(1);
t396 = t281 * t465;
t208 = rSges(3,1) * t270 + rSges(3,2) * t272;
t277 = qJD(1) + qJD(2);
t449 = t208 * t277;
t168 = -t396 - t449;
t276 = qJD(4) + qJD(5);
t466 = rSges(6,2) * t271;
t397 = t276 * t466;
t433 = t269 * t276;
t504 = -rSges(6,1) * t433 - t397;
t432 = t270 * t277;
t405 = pkin(2) * t432;
t503 = t405 - t396;
t502 = 2 * qJD(4);
t280 = sin(qJ(4));
t435 = t265 * t280;
t434 = t265 * t282;
t235 = rSges(5,1) * t434;
t497 = t264 * rSges(5,3) + t235;
t156 = -rSges(5,2) * t435 + t497;
t259 = t264 * pkin(7);
t261 = t265 * pkin(3);
t193 = t261 + t259;
t376 = t193 + t266;
t501 = t156 + t376;
t500 = -rSges(4,2) * t264 - t508;
t220 = rSges(6,2) * t445;
t499 = -t265 * rSges(6,3) - t220;
t273 = Icges(5,4) * t282;
t347 = -Icges(5,2) * t280 + t273;
t240 = Icges(5,1) * t280 + t273;
t442 = t264 * t280;
t413 = rSges(5,2) * t442 + t265 * rSges(5,3);
t441 = t264 * t282;
t155 = rSges(5,1) * t441 - t413;
t141 = t277 * t155;
t260 = t265 * pkin(7);
t192 = pkin(3) * t264 - t260;
t186 = t277 * t192;
t436 = t265 * t277;
t230 = pkin(7) * t436;
t407 = qJD(4) * t282;
t394 = rSges(5,2) * t407;
t431 = t277 * t280;
t330 = rSges(5,2) * t264 * t431 + rSges(5,3) * t436 - t265 * t394;
t408 = qJD(4) * t280;
t387 = t265 * t408;
t496 = -rSges(5,1) * t387 + t141 + t186 + t230 + t330;
t201 = Icges(6,5) * t269 + Icges(6,6) * t271;
t309 = Icges(6,3) * t277 - t201 * t276;
t443 = t264 * t277;
t325 = t346 * t265;
t138 = Icges(6,6) * t264 + t325;
t455 = t138 * t269;
t495 = -t202 * t443 + t265 * t309 + t277 * (-t140 * t271 + t455);
t494 = t264 * t309 + (t323 + t345) * t277;
t237 = Icges(5,5) * t282 - Icges(5,6) * t280;
t236 = Icges(5,5) * t280 + Icges(5,6) * t282;
t306 = Icges(5,3) * t277 - qJD(4) * t236;
t461 = Icges(5,4) * t280;
t241 = Icges(5,1) * t282 - t461;
t328 = t241 * t265;
t153 = Icges(5,5) * t264 + t328;
t326 = t347 * t265;
t151 = Icges(5,6) * t264 + t326;
t453 = t151 * t280;
t342 = -t153 * t282 + t453;
t493 = -t237 * t443 + t265 * t306 + t277 * t342;
t324 = t237 * t265;
t233 = Icges(5,4) * t442;
t152 = Icges(5,1) * t441 - Icges(5,5) * t265 - t233;
t150 = Icges(5,4) * t441 - Icges(5,2) * t442 - Icges(5,6) * t265;
t454 = t150 * t280;
t343 = -t152 * t282 + t454;
t492 = t264 * t306 + (t324 + t343) * t277;
t203 = Icges(6,2) * t271 + t460;
t340 = t203 * t269 - t205 * t271;
t491 = t202 * t276 + t277 * t340;
t243 = rSges(5,1) * t280 + rSges(5,2) * t282;
t410 = qJD(4) * t264;
t187 = t243 * t410;
t490 = t277 * t501 - t187;
t238 = Icges(5,2) * t282 + t461;
t339 = t280 * t238 - t240 * t282;
t489 = t237 * qJD(4) + t277 * t339;
t148 = Icges(5,5) * t441 - Icges(5,6) * t442 - Icges(5,3) * t265;
t66 = -t265 * t148 - t264 * t343;
t188 = t264 * t276;
t207 = rSges(6,1) * t269 + t466;
t388 = t264 * t408;
t226 = pkin(4) * t388;
t219 = t265 * t267;
t440 = t264 * t284;
t363 = t219 - t440;
t129 = t363 - t193;
t438 = t265 * t269;
t398 = rSges(6,2) * t438;
t143 = -t398 + t498;
t426 = -t129 - t143;
t488 = t277 * (t376 - t426) - t188 * t207 - t226;
t128 = t192 + t415;
t120 = t277 * t128;
t142 = t400 + t499;
t130 = t277 * t142;
t418 = rSges(6,3) * t436 + t277 * t220;
t487 = t130 + t186 + t418 - t120;
t361 = pkin(4) * t387;
t469 = pkin(3) - t267;
t108 = -t361 - t230 + (t264 * t469 - t247) * t277;
t416 = -t277 * t440 - t226;
t109 = (-t265 * t469 - t259) * t277 + t416;
t362 = t276 * t277;
t166 = t264 * t362;
t167 = t265 * t362;
t189 = t265 * t276;
t90 = -t265 * t397 + (-t265 * t433 - t271 * t443) * rSges(6,1) + t418;
t391 = t504 * t264 - t277 * t398;
t91 = t277 * t498 + t391;
t12 = t142 * t167 - t143 * t166 + t188 * t91 + t189 * t90 + ((t108 - t120) * t265 + (-t129 * t277 + t109) * t264) * qJD(4);
t164 = t207 * t264;
t165 = t207 * t265;
t467 = rSges(6,2) * t269;
t209 = rSges(6,1) * t271 - t467;
t51 = t142 * t188 + t143 * t189 + qJD(3) + (-t128 * t264 + t129 * t265) * qJD(4);
t321 = -t189 * t207 - t361;
t305 = t321 - t396;
t377 = -t192 - t472;
t56 = (t128 - t142 + t377) * t277 + t305;
t283 = cos(qJ(1));
t395 = t283 * t465;
t57 = t395 + t488;
t486 = -t56 * (t164 * t277 - t189 * t209) - t51 * (-t188 * t164 - t165 * t189) - t57 * (-t277 * t165 - t188 * t209) + t12 * (t264 * t142 + t265 * t143);
t420 = -Icges(5,2) * t441 + t152 - t233;
t422 = t240 * t264 + t150;
t485 = -t280 * t420 - t282 * t422;
t484 = t188 * (-t203 * t265 + t140) - t189 * (-Icges(6,2) * t444 + t139 - t217) + t277 * t509;
t275 = t277 ^ 2;
t483 = t166 / 0.2e1;
t482 = t167 / 0.2e1;
t481 = -t188 / 0.2e1;
t480 = t188 / 0.2e1;
t479 = -t189 / 0.2e1;
t478 = t189 / 0.2e1;
t477 = t264 / 0.2e1;
t476 = -t265 / 0.2e1;
t475 = -t277 / 0.2e1;
t474 = t277 / 0.2e1;
t473 = pkin(1) * t281;
t471 = pkin(2) * t275;
t274 = t283 * pkin(1);
t468 = rSges(5,1) * t282;
t263 = t272 * rSges(3,1);
t451 = t201 * t265;
t93 = -t264 * t340 - t451;
t464 = t93 * t277;
t447 = t236 * t265;
t111 = -t264 * t339 - t447;
t457 = t111 * t277;
t452 = t201 * t264;
t450 = t203 * t276;
t448 = t236 * t264;
t446 = t237 * t277;
t135 = Icges(6,5) * t444 - Icges(6,6) * t445 - Icges(6,3) * t265;
t430 = -t264 * t135 - t139 * t437;
t428 = -t264 * t148 - t152 * t434;
t149 = Icges(5,3) * t264 + t324;
t427 = t264 * t149 + t153 * t434;
t421 = -t240 * t265 - t151;
t419 = -t238 * t265 + t153;
t412 = -t238 + t241;
t411 = t240 + t347;
t409 = qJD(4) * t265;
t76 = t395 + t490;
t406 = t76 * t472;
t404 = (qJD(4) ^ 2) * t470;
t286 = qJD(1) ^ 2;
t403 = t286 * t473;
t402 = t286 * t274;
t401 = t264 * t91 + (t130 + t90) * t265;
t392 = t265 * t431;
t390 = rSges(5,1) * t388 + rSges(5,2) * t392 + t264 * t394;
t389 = t243 * t409;
t386 = t443 / 0.2e1;
t385 = t436 / 0.2e1;
t384 = -pkin(3) - t468;
t383 = -t410 / 0.2e1;
t380 = t409 / 0.2e1;
t190 = rSges(4,1) * t264 + rSges(4,2) * t265;
t319 = -t190 - t472;
t375 = -pkin(4) * t280 - t207;
t210 = -rSges(3,2) * t270 + t263;
t311 = Icges(6,5) * t277 - t205 * t276;
t373 = -t137 * t276 + t264 * t311 + t277 * t327;
t372 = -t138 * t276 - t206 * t443 + t265 * t311;
t310 = Icges(6,6) * t277 - t450;
t371 = t139 * t276 + t264 * t310 + t277 * t325;
t370 = t140 * t276 + t265 * t310 - t346 * t443;
t125 = t153 * t441;
t369 = t265 * t149 - t125;
t368 = -t135 + t455;
t366 = -t148 + t453;
t365 = t509 * t276;
t364 = t206 * t276 - t450;
t185 = -rSges(3,2) * t432 + t263 * t277;
t183 = t209 * t276;
t357 = -pkin(4) * t407 - t183;
t113 = t140 * t444;
t355 = t138 * t445 - t113;
t353 = -rSges(5,2) * t280 + t468;
t352 = -t264 * t57 - t265 * t56;
t320 = -t389 - t396;
t75 = (-t155 + t377) * t277 + t320;
t351 = -t264 * t76 - t265 * t75;
t350 = -t391 - t416;
t81 = t137 * t271 + t139 * t269;
t97 = t150 * t282 + t152 * t280;
t98 = t151 * t282 + t153 * t280;
t341 = t155 * t264 + t156 * t265;
t67 = -t151 * t442 - t369;
t334 = (t264 * t67 - t265 * t66) * qJD(4);
t68 = -t150 * t435 - t428;
t69 = -t151 * t435 + t427;
t333 = (t264 * t69 - t265 * t68) * qJD(4);
t332 = -t270 * t471 - t403;
t331 = -t272 * t471 - t402;
t329 = t363 - t507;
t318 = t277 * (-pkin(3) * t443 + t230) + t332;
t317 = t188 * t451 - t189 * t452 - t202 * t277;
t316 = -t280 * t419 + t282 * t421;
t297 = t135 * t277 - t269 * t371 + t271 * t373;
t13 = t264 * t494 + t297 * t265;
t296 = t136 * t277 - t269 * t370 + t271 * t372;
t14 = t264 * t495 + t296 * t265;
t15 = t297 * t264 - t265 * t494;
t16 = t296 * t264 - t265 * t495;
t62 = -t135 * t265 - t322;
t63 = -t136 * t265 - t355;
t28 = t188 * t63 - t189 * t62 + t464;
t300 = (-t205 * t265 - t138) * t188 - (-t205 * t264 - t137) * t189 + (-t203 + t206) * t277;
t288 = -t269 * t484 + t300 * t271;
t64 = -t137 * t438 - t430;
t65 = -t138 * t438 + t429;
t94 = -t265 * t340 + t452;
t92 = t94 * t277;
t29 = t188 * t65 - t189 * t64 + t92;
t40 = t269 * t373 + t271 * t371;
t41 = t269 * t372 + t271 * t370;
t295 = t201 * t277 - t269 * t365 + t271 * t364;
t44 = t264 * t491 + t295 * t265;
t45 = t295 * t264 - t265 * t491;
t82 = t138 * t271 + t140 * t269;
t314 = (-t13 * t189 + t14 * t188 + t166 * t64 + t167 * t65 + t277 * t44) * t477 + (-t264 * t317 + t265 * t288) * t481 + (t264 * t288 + t265 * t317) * t478 + (-t15 * t189 + t16 * t188 + t166 * t62 + t167 * t63 + t277 * t45) * t476 + (t300 * t269 + t271 * t484) * t475 + t28 * t386 + t29 * t385 + ((t277 * t65 - t13) * t265 + (t277 * t64 + t14) * t264) * t480 + (t264 * t63 - t265 * t62) * t483 + (t264 * t65 - t265 * t64) * t482 + ((t277 * t63 - t15) * t265 + (t277 * t62 + t16) * t264) * t479 + ((t277 * t82 - t40) * t265 + (t277 * t81 + t41) * t264) * t474;
t313 = -t499 + t505;
t312 = (-t280 * t411 + t282 * t412) * t277;
t308 = Icges(5,5) * t277 - qJD(4) * t240;
t307 = Icges(5,6) * t277 - qJD(4) * t238;
t304 = t264 * t384 + t260 + t413 - t472;
t106 = (-t277 * t441 - t387) * rSges(5,1) + t330;
t107 = t277 * t497 - t390;
t302 = (t106 + t141) * t265 + (-t156 * t277 + t107) * t264;
t102 = t265 * t307 - t347 * t443;
t104 = -t241 * t443 + t265 * t308;
t294 = -qJD(4) * t98 - t102 * t280 + t104 * t282 + t149 * t277;
t103 = t264 * t307 + t277 * t326;
t105 = t264 * t308 + t277 * t328;
t293 = -qJD(4) * t97 - t103 * t280 + t105 * t282 + t148 * t277;
t213 = t347 * qJD(4);
t214 = t241 * qJD(4);
t292 = -t213 * t280 + t214 * t282 + t236 * t277 + (-t238 * t282 - t240 * t280) * qJD(4);
t133 = t277 * t319 - t396;
t134 = t277 * t500 + t395;
t291 = (t133 * t508 + t134 * t319) * t277;
t112 = -t265 * t339 + t448;
t110 = t112 * t277;
t32 = t334 + t457;
t33 = t110 + t333;
t49 = -qJD(4) * t343 + t103 * t282 + t105 * t280;
t50 = -qJD(4) * t342 + t102 * t282 + t104 * t280;
t54 = t264 * t489 + t292 * t265;
t55 = t292 * t264 - t265 * t489;
t290 = (t92 + (t63 + (t136 + t456) * t265 + t355 + t430) * t189 + (-t265 * t368 - t506 + t62) * t188) * t478 + (t110 + ((t67 - t125 + (t149 + t454) * t265 + t428) * t265 + t427 * t264) * qJD(4)) * t380 + (t81 + t93) * t483 + (t82 + t94) * t482 + (-t464 + (t65 + t506) * t189 + (t368 * t264 - t113 + t64) * t188 + ((t136 + t345) * t188 + t368 * t189) * t265 + t28) * t481 + (t41 + t44) * t480 + (t32 - t457 + ((t265 * t366 - t427 + t69) * t265 + (t264 * t366 + t369 + t68) * t264) * qJD(4)) * t383 + (t50 + t54) * t410 / 0.2e1 + (-qJD(4) * t339 + t213 * t282 + t214 * t280 + t269 * t364 + t271 * t365) * t277 + (t40 + t45 + t29) * t479 - (t49 + t55 + t33) * t409 / 0.2e1 + ((t111 + t97) * t264 + (t112 + t98) * t265) * qJD(4) * t474;
t289 = (t75 * (-t235 - t261 - t266) - t406 + (t75 * (-rSges(5,3) - pkin(7)) + t76 * t384) * t264) * t277;
t36 = -t264 * t404 - t167 * t207 - t183 * t188 + (t108 + t90 - t361) * t277 + t318;
t287 = (t56 * (-t219 + t507) + t57 * t505) * t277 + (-t36 * t467 + t57 * (-pkin(4) * t408 + t504)) * t265;
t227 = rSges(4,2) * t443;
t222 = t353 * qJD(4);
t182 = t277 * t190;
t178 = t243 * t265;
t177 = t243 * t264;
t170 = t193 * t277;
t169 = t210 * t277 + t395;
t147 = -t185 * t277 - t402;
t146 = -t277 * t449 - t403;
t118 = -t277 * (rSges(4,1) * t436 - t227) + t331;
t117 = -t190 * t275 + t332;
t83 = qJD(4) * t341 + qJD(3);
t61 = -t222 * t409 + (-t107 - t170 + t187) * t277 + t331;
t60 = t106 * t277 + (-t222 * t264 - t243 * t436) * qJD(4) + t318;
t46 = t302 * qJD(4);
t37 = -t265 * t404 + t166 * t207 - t183 * t189 + (-t109 - t170 - t91 + t226) * t277 + t331;
t1 = [m(3) * (t147 * (-t208 - t473) + t146 * (t210 + t274) + (-t185 - t395 + t169) * t168) + t290 + (t37 * (t313 - t473) + t56 * (t350 - t395) + t36 * (t274 + t329) + t287 + (-t305 + t56 + t487 + t503) * t57) * m(6) + (t61 * (t304 - t473) + t75 * (t390 - t395) + t60 * (t274 + t501) + t289 + (-t320 + t75 + t496 + t503) * t76) * m(5) + (t118 * (t319 - t473) + t133 * (t227 - t395) + t117 * (t274 + t500) + t291 + (t133 + t182 + t405) * t134) * m(4); t290 + (t37 * t313 + t36 * t329 + t287 + (t277 * t472 - t321 + t487) * t57 + (t350 + t488) * t56) * m(6) + (t406 * t277 + t61 * t304 + t60 * t501 + t289 + (t389 + t496) * t76 + (t390 + t490) * t75) * m(5) + (t117 * t500 + t118 * t319 + t133 * t227 + t291 + t134 * t182 - (-t133 * t500 - t134 * t472) * t277) * m(4) + (t146 * t210 - t147 * t208 - t168 * t185 - t169 * t449 - (-t168 * t210 - t169 * t208) * t277) * m(3); m(5) * t46 + m(6) * t12; t314 + ((t277 * t98 - t49) * t265 + (t277 * t97 + t50) * t264) * t474 + ((t280 * t412 + t282 * t411) * t277 + ((t264 * t419 - t265 * t420) * t282 + (t264 * t421 + t265 * t422) * t280) * qJD(4)) * t475 + ((-t409 * t448 - t446) * t265 + (t312 + (t316 * t264 + (t447 - t485) * t265) * qJD(4)) * t264) * t380 + ((-t410 * t447 + t446) * t264 + (t312 + (-t485 * t265 + (t448 + t316) * t264) * qJD(4)) * t265) * t383 + (t277 * t54 + ((-t264 * t492 - t293 * t265 + t277 * t69) * t265 + (t264 * t493 + t294 * t265 + t277 * t68) * t264) * t502) * t477 + (t277 * t55 + ((-t293 * t264 + t265 * t492 + t277 * t67) * t265 + (t294 * t264 - t265 * t493 + t277 * t66) * t264) * t502) * t476 + (t32 + t334) * t386 + (t33 + t333) * t385 + (t51 * t401 + (t37 * t375 + t56 * t357 + t12 * t129 + t51 * t108 + (-t51 * t128 + t375 * t57) * t277) * t265 + (t36 * t375 + t57 * t357 - t12 * t128 + t51 * t109 + (t56 * t207 + t426 * t51) * t277) * t264 - (-t57 * t392 + (t352 * t282 + t51 * (-t264 ^ 2 - t265 ^ 2) * t280) * qJD(4)) * pkin(4) + t486) * m(6) + (-(t177 * t75 - t178 * t76) * t277 - (t83 * (-t177 * t264 - t178 * t265) + t351 * t353) * qJD(4) + t46 * t341 + t83 * t302 + t351 * t222 + ((-t277 * t76 - t61) * t265 + (t277 * t75 - t60) * t264) * t243) * m(5); t314 + (t51 * (-t143 * t443 + t401) + t352 * t183 + ((-t277 * t57 - t37) * t265 + (t277 * t56 - t36) * t264) * t207 + t486) * m(6);];
tauc = t1(:);
