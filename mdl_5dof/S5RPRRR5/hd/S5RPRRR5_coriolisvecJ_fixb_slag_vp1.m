% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR5_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR5_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 09:48:35
% EndTime: 2022-01-20 09:48:52
% DurationCPUTime: 8.56s
% Computational Cost: add. (17538->515), mult. (11058->689), div. (0->0), fcn. (8580->10), ass. (0->318)
t280 = qJD(1) ^ 2;
t273 = qJ(4) + qJ(5);
t267 = cos(t273);
t259 = Icges(6,4) * t267;
t266 = sin(t273);
t203 = Icges(6,1) * t266 + t259;
t337 = -Icges(6,2) * t266 + t259;
t491 = t203 + t337;
t272 = qJ(1) + pkin(9);
t265 = qJ(3) + t272;
t260 = sin(t265);
t430 = t260 * t266;
t213 = Icges(6,4) * t430;
t261 = cos(t265);
t429 = t260 * t267;
t139 = Icges(6,1) * t429 - Icges(6,5) * t261 - t213;
t137 = Icges(6,4) * t429 - Icges(6,2) * t430 - Icges(6,6) * t261;
t440 = t137 * t266;
t336 = -t139 * t267 + t440;
t314 = t336 * t260;
t200 = Icges(6,5) * t267 - Icges(6,6) * t266;
t315 = t200 * t261;
t136 = Icges(6,3) * t260 + t315;
t445 = Icges(6,4) * t266;
t204 = Icges(6,1) * t267 - t445;
t319 = t204 * t261;
t140 = Icges(6,5) * t260 + t319;
t422 = t261 * t267;
t414 = t260 * t136 + t140 * t422;
t490 = -t314 - t414;
t270 = qJD(4) + qJD(5);
t452 = rSges(6,2) * t267;
t386 = t270 * t452;
t418 = t266 * t270;
t488 = -rSges(6,1) * t418 - t386;
t486 = 2 * qJD(4);
t274 = sin(qJ(4));
t427 = t260 * t274;
t398 = rSges(5,2) * t427 + t261 * rSges(5,3);
t276 = cos(qJ(4));
t426 = t260 * t276;
t153 = rSges(5,1) * t426 - t398;
t271 = qJD(1) + qJD(3);
t141 = t271 * t153;
t255 = t261 * pkin(7);
t189 = pkin(3) * t260 - t255;
t183 = t271 * t189;
t485 = -t141 - t183;
t399 = rSges(6,1) * t422 + t260 * rSges(6,3);
t419 = t261 * t276;
t484 = rSges(5,1) * t419 + t260 * rSges(5,3);
t254 = t260 * pkin(7);
t190 = t261 * pkin(3) + t254;
t483 = pkin(2) * cos(t272) + cos(qJ(1)) * pkin(1);
t268 = Icges(5,4) * t276;
t338 = -Icges(5,2) * t274 + t268;
t236 = Icges(5,1) * t274 + t268;
t278 = -pkin(8) - pkin(7);
t242 = t261 * t278;
t458 = pkin(4) * t276;
t262 = pkin(3) + t458;
t400 = -t260 * t262 - t242;
t128 = t189 + t400;
t118 = t271 * t128;
t216 = rSges(6,2) * t430;
t142 = rSges(6,1) * t429 - t261 * rSges(6,3) - t216;
t130 = t271 * t142;
t481 = t118 - t130 - t183;
t199 = Icges(6,5) * t266 + Icges(6,6) * t267;
t301 = Icges(6,3) * t271 - t199 * t270;
t428 = t260 * t271;
t317 = t337 * t261;
t138 = Icges(6,6) * t260 + t317;
t439 = t138 * t266;
t480 = -t200 * t428 + t261 * t301 + t271 * (-t140 * t267 + t439);
t479 = t260 * t301 + (t315 + t336) * t271;
t233 = Icges(5,5) * t276 - Icges(5,6) * t274;
t232 = Icges(5,5) * t274 + Icges(5,6) * t276;
t298 = Icges(5,3) * t271 - qJD(4) * t232;
t446 = Icges(5,4) * t274;
t237 = Icges(5,1) * t276 - t446;
t320 = t237 * t261;
t151 = Icges(5,5) * t260 + t320;
t318 = t338 * t261;
t149 = Icges(5,6) * t260 + t318;
t437 = t149 * t274;
t333 = -t151 * t276 + t437;
t477 = -t233 * t428 + t261 * t298 + t271 * t333;
t316 = t233 * t261;
t229 = Icges(5,4) * t427;
t150 = Icges(5,1) * t426 - Icges(5,5) * t261 - t229;
t148 = Icges(5,4) * t426 - Icges(5,2) * t427 - Icges(5,6) * t261;
t438 = t148 * t274;
t334 = -t150 * t276 + t438;
t476 = t260 * t298 + (t316 + t334) * t271;
t201 = Icges(6,2) * t267 + t445;
t331 = t201 * t266 - t203 * t267;
t475 = t200 * t270 + t271 * t331;
t238 = rSges(5,1) * t274 + rSges(5,2) * t276;
t394 = qJD(4) * t260;
t184 = t238 * t394;
t420 = t261 * t274;
t154 = -rSges(5,2) * t420 + t484;
t310 = t154 + t190;
t474 = -t271 * t310 + t184;
t234 = Icges(5,2) * t276 + t446;
t330 = t234 * t274 - t236 * t276;
t473 = t233 * qJD(4) + t271 * t330;
t146 = Icges(5,5) * t426 - Icges(5,6) * t427 - Icges(5,3) * t261;
t66 = -t261 * t146 - t260 * t334;
t421 = t261 * t271;
t226 = pkin(7) * t421;
t392 = qJD(4) * t274;
t376 = t261 * t392;
t352 = pkin(4) * t376;
t457 = pkin(3) - t262;
t108 = -t352 - t226 + (t260 * t457 - t242) * t271;
t378 = t260 * t392;
t222 = pkin(4) * t378;
t425 = t260 * t278;
t401 = -t271 * t425 - t222;
t109 = (-t261 * t457 - t254) * t271 + t401;
t354 = t261 * t262 - t425;
t129 = t354 - t190;
t423 = t261 * t266;
t387 = rSges(6,2) * t423;
t143 = -t387 + t399;
t353 = t270 * t271;
t164 = t260 * t353;
t165 = t261 * t353;
t185 = t260 * t270;
t186 = t261 * t270;
t384 = t267 * t428;
t403 = rSges(6,3) * t421 + t271 * t216;
t90 = -t261 * t386 + (-t261 * t418 - t384) * rSges(6,1) + t403;
t380 = t488 * t260 - t271 * t387;
t91 = t271 * t399 + t380;
t12 = t142 * t165 - t143 * t164 + t185 * t91 + t186 * t90 + ((t108 - t118) * t261 + (-t129 * t271 + t109) * t260) * qJD(4);
t205 = rSges(6,1) * t266 + t452;
t162 = t205 * t260;
t163 = t205 * t261;
t453 = rSges(6,2) * t266;
t455 = rSges(6,1) * t267;
t206 = -t453 + t455;
t51 = t142 * t185 + t143 * t186 + qJD(2) + (-t128 * t260 + t129 * t261) * qJD(4);
t313 = -t186 * t205 - t352;
t348 = -pkin(2) * sin(t272) - sin(qJ(1)) * pkin(1);
t326 = t348 * qJD(1);
t289 = t326 + t313;
t57 = (t128 - t142 - t189) * t271 + t289;
t325 = t483 * qJD(1);
t411 = -t129 - t143;
t381 = t190 - t411;
t408 = t185 * t205 + t222;
t58 = t271 * t381 + t325 - t408;
t472 = -t57 * (t162 * t271 - t186 * t206) - t51 * (-t185 * t162 - t163 * t186) - t58 * (-t271 * t163 - t185 * t206) + t12 * (t260 * t142 + t261 * t143);
t405 = -Icges(5,2) * t426 + t150 - t229;
t407 = t236 * t260 + t148;
t471 = -t274 * t405 - t276 * t407;
t470 = t185 * (-t201 * t261 + t140) - t186 * (-Icges(6,2) * t429 + t139 - t213) + t271 * t491;
t469 = t164 / 0.2e1;
t468 = t165 / 0.2e1;
t467 = -t185 / 0.2e1;
t466 = t185 / 0.2e1;
t465 = -t186 / 0.2e1;
t464 = t186 / 0.2e1;
t463 = t260 / 0.2e1;
t462 = -t261 / 0.2e1;
t461 = -t271 / 0.2e1;
t460 = t271 / 0.2e1;
t456 = rSges(5,1) * t276;
t253 = t261 * rSges(4,1);
t393 = qJD(4) * t261;
t377 = t238 * t393;
t297 = t326 - t377;
t75 = (-t153 - t189) * t271 + t297;
t451 = t261 * t75;
t450 = t271 * t57;
t435 = t199 * t261;
t93 = -t260 * t331 - t435;
t449 = t93 * t271;
t432 = t232 * t261;
t111 = -t260 * t330 - t432;
t442 = t111 * t271;
t188 = -rSges(4,2) * t260 + t253;
t134 = t188 * t271 + t325;
t187 = rSges(4,1) * t260 + rSges(4,2) * t261;
t441 = t134 * t187;
t436 = t199 * t260;
t434 = t201 * t270;
t433 = t232 * t260;
t431 = t233 * t271;
t417 = t271 * t187;
t416 = t271 * t274;
t135 = Icges(6,5) * t429 - Icges(6,6) * t430 - Icges(6,3) * t261;
t415 = -t260 * t135 - t139 * t422;
t413 = -t260 * t146 - t150 * t419;
t147 = Icges(5,3) * t260 + t316;
t412 = t260 * t147 + t151 * t419;
t406 = -t236 * t261 - t149;
t404 = -t234 * t261 + t151;
t397 = -t234 + t237;
t396 = t236 + t338;
t391 = qJD(4) * t276;
t390 = (qJD(4) ^ 2) * t458;
t389 = t260 * t91 + (t130 + t90) * t261;
t385 = rSges(5,2) * t391;
t382 = t261 * t416;
t379 = rSges(5,1) * t378 + rSges(5,2) * t382 + t260 * t385;
t375 = t428 / 0.2e1;
t374 = t421 / 0.2e1;
t373 = -pkin(3) - t456;
t372 = -t394 / 0.2e1;
t369 = t393 / 0.2e1;
t367 = -pkin(4) * t274 - t205;
t303 = Icges(6,5) * t271 - t203 * t270;
t364 = -t137 * t270 + t260 * t303 + t271 * t319;
t363 = -t138 * t270 - t204 * t428 + t261 * t303;
t302 = Icges(6,6) * t271 - t434;
t362 = t139 * t270 + t260 * t302 + t271 * t317;
t361 = t140 * t270 + t261 * t302 - t337 * t428;
t125 = t151 * t426;
t360 = t261 * t147 - t125;
t359 = -t135 + t439;
t357 = -t146 + t437;
t356 = t491 * t270;
t355 = t204 * t270 - t434;
t182 = t206 * t270;
t349 = -pkin(4) * t391 - t182;
t113 = t140 * t429;
t346 = t138 * t430 - t113;
t343 = -rSges(5,2) * t274 + t456;
t342 = -t260 * t58 - t261 * t57;
t76 = t325 - t474;
t341 = -t260 * t76 - t451;
t81 = t137 * t267 + t139 * t266;
t97 = t148 * t276 + t150 * t274;
t98 = t149 * t276 + t151 * t274;
t332 = t153 * t260 + t154 * t261;
t329 = t354 + t399;
t328 = t348 * t280;
t327 = t483 * t280;
t67 = -t149 * t427 - t360;
t323 = (t260 * t67 - t261 * t66) * qJD(4);
t68 = -t148 * t420 - t413;
t69 = -t149 * t420 + t412;
t322 = (t260 * t69 - t261 * t68) * qJD(4);
t321 = rSges(5,2) * t260 * t416 + rSges(5,3) * t421 - t261 * t385;
t311 = -t142 + t400;
t309 = t271 * (-pkin(3) * t428 + t226) + t328;
t308 = t185 * t435 - t186 * t436 - t200 * t271;
t307 = -t274 * t404 + t276 * t406;
t306 = t260 * t373 + t255 + t398;
t133 = t326 - t417;
t291 = t135 * t271 - t266 * t362 + t267 * t364;
t13 = t479 * t260 + t291 * t261;
t290 = t136 * t271 - t266 * t361 + t267 * t363;
t14 = t480 * t260 + t290 * t261;
t15 = t291 * t260 - t261 * t479;
t16 = t290 * t260 - t261 * t480;
t62 = -t135 * t261 - t314;
t63 = -t136 * t261 - t346;
t28 = t185 * t63 - t186 * t62 + t449;
t294 = (-t203 * t261 - t138) * t185 - (-t203 * t260 - t137) * t186 + (-t201 + t204) * t271;
t283 = -t266 * t470 + t294 * t267;
t64 = -t137 * t423 - t415;
t65 = -t138 * t423 + t414;
t94 = -t261 * t331 + t436;
t92 = t94 * t271;
t29 = t185 * t65 - t186 * t64 + t92;
t40 = t266 * t364 + t267 * t362;
t41 = t266 * t363 + t267 * t361;
t288 = t199 * t271 - t266 * t356 + t267 * t355;
t44 = t475 * t260 + t288 * t261;
t45 = t288 * t260 - t261 * t475;
t82 = t138 * t267 + t140 * t266;
t305 = (-t13 * t186 + t14 * t185 + t164 * t64 + t165 * t65 + t271 * t44) * t463 + (-t260 * t308 + t261 * t283) * t467 + (t260 * t283 + t261 * t308) * t464 + (-t15 * t186 + t16 * t185 + t164 * t62 + t165 * t63 + t271 * t45) * t462 + (t294 * t266 + t267 * t470) * t461 + t28 * t375 + t29 * t374 + ((t271 * t65 - t13) * t261 + (t271 * t64 + t14) * t260) * t466 + (t260 * t63 - t261 * t62) * t469 + (t260 * t65 - t261 * t64) * t468 + ((t271 * t63 - t15) * t261 + (t271 * t62 + t16) * t260) * t465 + ((t271 * t82 - t40) * t261 + (t271 * t81 + t41) * t260) * t460;
t304 = (-t274 * t396 + t276 * t397) * t271;
t300 = Icges(5,5) * t271 - qJD(4) * t236;
t299 = Icges(5,6) * t271 - qJD(4) * t234;
t106 = (-t271 * t426 - t376) * rSges(5,1) + t321;
t107 = t271 * t484 - t379;
t295 = (t106 + t141) * t261 + (-t154 * t271 + t107) * t260;
t102 = t261 * t299 - t338 * t428;
t104 = -t237 * t428 + t261 * t300;
t287 = -qJD(4) * t98 - t102 * t274 + t104 * t276 + t147 * t271;
t103 = t260 * t299 + t271 * t318;
t105 = t260 * t300 + t271 * t320;
t286 = -qJD(4) * t97 - t103 * t274 + t105 * t276 + t146 * t271;
t209 = t338 * qJD(4);
t210 = t237 * qJD(4);
t285 = -t209 * t274 + t210 * t276 + t232 * t271 + (-t234 * t276 - t236 * t274) * qJD(4);
t112 = -t261 * t330 + t433;
t110 = t112 * t271;
t32 = t323 + t442;
t33 = t110 + t322;
t49 = -qJD(4) * t334 + t103 * t276 + t105 * t274;
t50 = -qJD(4) * t333 + t102 * t276 + t104 * t274;
t54 = t473 * t260 + t285 * t261;
t55 = t285 * t260 - t261 * t473;
t284 = (t92 + (t63 + (t136 + t440) * t261 + t346 + t415) * t186 + (-t261 * t359 - t490 + t62) * t185) * t464 + (t110 + ((t67 - t125 + (t147 + t438) * t261 + t413) * t261 + t412 * t260) * qJD(4)) * t369 + (t81 + t93) * t469 + (t82 + t94) * t468 + (t28 - t449 + (t65 + t490) * t186 + (t359 * t260 - t113 + t64) * t185 + ((t136 + t336) * t185 + t359 * t186) * t261) * t467 + (t41 + t44) * t466 + (t32 - t442 + ((t261 * t357 - t412 + t69) * t261 + (t260 * t357 + t360 + t68) * t260) * qJD(4)) * t372 + (t50 + t54) * t394 / 0.2e1 + (-qJD(4) * t330 + t209 * t276 + t210 * t274 + t266 * t355 + t267 * t356) * t271 + (t40 + t45 + t29) * t465 - (t49 + t55 + t33) * t393 / 0.2e1 + ((t111 + t97) * t260 + (t112 + t98) * t261) * qJD(4) * t460;
t282 = t75 * t379 + t76 * (-rSges(5,1) * t376 + t226 + t321) + (t373 * t451 + (t75 * (-rSges(5,3) - pkin(7)) + t76 * t373) * t260) * t271;
t36 = -t260 * t390 - t165 * t205 - t182 * t185 + (t108 + t90 - t352) * t271 + t309;
t281 = t57 * (-rSges(6,3) * t428 - t380 - t401) + t58 * (-rSges(6,1) * t384 - t262 * t428 + t403) + (-t36 * t453 + t58 * (-pkin(4) * t392 + t488) + (t57 * (-t262 - t455) - t58 * t278) * t271) * t261;
t223 = rSges(4,2) * t428;
t218 = t343 * qJD(4);
t178 = t238 * t261;
t177 = t238 * t260;
t168 = t190 * t271;
t167 = rSges(4,1) * t421 - t223;
t120 = -t167 * t271 - t327;
t119 = -t271 * t417 + t328;
t83 = qJD(4) * t332 + qJD(2);
t61 = -t218 * t393 - t327 + (-t107 - t168 + t184) * t271;
t60 = t106 * t271 + (-t218 * t260 - t238 * t421) * qJD(4) + t309;
t46 = t295 * qJD(4);
t37 = -t261 * t390 + t164 * t205 - t182 * t186 - t327 + (-t109 - t168 - t91 + t222) * t271;
t1 = [m(4) * (t120 * (-t187 + t348) + t133 * t223 + t119 * (t188 + t483) + (-t133 * t253 - t441) * t271 + (-t133 * t483 + t134 * t348) * qJD(1)) + t284 + (t37 * (t311 + t348) + t36 * (t329 + t483) + (t348 * t58 - t483 * t57) * qJD(1) + t281 - (-t57 + t289 + t481) * t58) * m(6) + (t61 * (t306 + t348) + t60 * (t310 + t483) + (t348 * t76 - t483 * t75) * qJD(1) + t282 - (-t75 + t297 + t485) * t76) * m(5); m(5) * t46 + m(6) * t12; t284 + (-t57 * t408 - t58 * (t313 + t481) + t381 * t450 + t37 * t311 + t36 * t329 + t281) * m(6) + (-t75 * t474 - t76 * (-t377 + t485) + t61 * t306 + t60 * t310 + t282) * m(5) + (t119 * t188 - t120 * t187 - t133 * t167 - t134 * t417 - (-t133 * t188 - t441) * t271) * m(4); t305 + ((t271 * t98 - t49) * t261 + (t271 * t97 + t50) * t260) * t460 + ((t274 * t397 + t276 * t396) * t271 + ((t260 * t404 - t261 * t405) * t276 + (t260 * t406 + t261 * t407) * t274) * qJD(4)) * t461 + ((-t394 * t432 + t431) * t260 + (t304 + (-t471 * t261 + (t433 + t307) * t260) * qJD(4)) * t261) * t372 + ((-t393 * t433 - t431) * t261 + (t304 + (t307 * t260 + (t432 - t471) * t261) * qJD(4)) * t260) * t369 + (t271 * t54 + ((-t476 * t260 - t286 * t261 + t271 * t69) * t261 + (t477 * t260 + t287 * t261 + t271 * t68) * t260) * t486) * t463 + (t271 * t55 + ((-t286 * t260 + t261 * t476 + t271 * t67) * t261 + (t287 * t260 - t261 * t477 + t271 * t66) * t260) * t486) * t462 + (t32 + t323) * t375 + (t33 + t322) * t374 + (t51 * t389 + (t37 * t367 + t57 * t349 + t12 * t129 + t51 * t108 + (-t51 * t128 + t367 * t58) * t271) * t261 + (t36 * t367 + t58 * t349 - t12 * t128 + t51 * t109 + (t57 * t205 + t411 * t51) * t271) * t260 - (-t58 * t382 + (t342 * t276 + t51 * (-t260 ^ 2 - t261 ^ 2) * t274) * qJD(4)) * pkin(4) + t472) * m(6) + (-(t177 * t75 - t178 * t76) * t271 - (t83 * (-t177 * t260 - t178 * t261) + t341 * t343) * qJD(4) + t46 * t332 + t83 * t295 + t341 * t218 + ((-t271 * t76 - t61) * t261 + (t271 * t75 - t60) * t260) * t238) * m(5); t305 + (t51 * (-t143 * t428 + t389) + t342 * t182 + ((-t271 * t58 - t37) * t261 + (-t36 + t450) * t260) * t205 + t472) * m(6);];
tauc = t1(:);
