% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:05
% EndTime: 2019-12-05 17:44:43
% DurationCPUTime: 18.31s
% Computational Cost: add. (18911->734), mult. (21699->1010), div. (0->0), fcn. (20721->10), ass. (0->358)
t295 = pkin(9) + qJ(4);
t287 = qJ(5) + t295;
t280 = sin(t287);
t281 = cos(t287);
t303 = cos(qJ(1));
t300 = cos(pkin(8));
t302 = sin(qJ(1));
t454 = t300 * t302;
t211 = t280 * t454 + t281 * t303;
t449 = t303 * t280;
t452 = t302 * t281;
t212 = t300 * t452 - t449;
t298 = sin(pkin(8));
t456 = t298 * t302;
t114 = -Icges(6,5) * t212 + Icges(6,6) * t211 - Icges(6,3) * t456;
t213 = t300 * t449 - t452;
t453 = t300 * t303;
t214 = t280 * t302 + t281 * t453;
t455 = t298 * t303;
t115 = Icges(6,5) * t214 - Icges(6,6) * t213 + Icges(6,3) * t455;
t201 = Icges(6,4) * t214;
t118 = -Icges(6,2) * t213 + Icges(6,6) * t455 + t201;
t200 = Icges(6,4) * t213;
t122 = -Icges(6,1) * t214 - Icges(6,5) * t455 + t200;
t346 = t118 * t213 + t214 * t122;
t469 = Icges(6,4) * t212;
t117 = Icges(6,2) * t211 - Icges(6,6) * t456 - t469;
t199 = Icges(6,4) * t211;
t120 = -Icges(6,1) * t212 - Icges(6,5) * t456 + t199;
t487 = t117 * t211 - t120 * t212;
t524 = t346 + t487 + (-t114 * t302 - t115 * t303) * t298;
t296 = qJD(4) + qJD(5);
t378 = t296 * t298;
t246 = t302 * t378;
t255 = -t296 * t300 + qJD(1);
t46 = t114 * t455 - t117 * t213 + t120 * t214;
t190 = -Icges(6,3) * t300 + (Icges(6,5) * t281 - Icges(6,6) * t280) * t298;
t467 = Icges(6,4) * t281;
t191 = -Icges(6,6) * t300 + (-Icges(6,2) * t280 + t467) * t298;
t468 = Icges(6,4) * t280;
t192 = -Icges(6,5) * t300 + (Icges(6,1) * t281 - t468) * t298;
t65 = t190 * t455 - t191 * t213 + t192 * t214;
t523 = t246 * t46 - t255 * t65;
t285 = sin(t295);
t448 = t303 * t285;
t286 = cos(t295);
t451 = t302 * t286;
t230 = t300 * t448 - t451;
t460 = t285 * t302;
t231 = t286 * t453 + t460;
t134 = Icges(5,5) * t231 - Icges(5,6) * t230 + Icges(5,3) * t455;
t218 = Icges(5,4) * t231;
t137 = -Icges(5,2) * t230 + Icges(5,6) * t455 + t218;
t217 = Icges(5,4) * t230;
t141 = -Icges(5,1) * t231 - Icges(5,5) * t455 + t217;
t60 = -t134 * t300 - (t137 * t285 + t141 * t286) * t298;
t54 = -t115 * t300 - (t118 * t280 + t122 * t281) * t298;
t521 = t302 * qJ(2);
t345 = -t137 * t230 - t231 * t141;
t228 = t285 * t454 + t286 * t303;
t229 = t300 * t451 - t448;
t472 = Icges(5,4) * t229;
t136 = Icges(5,2) * t228 - Icges(5,6) * t456 - t472;
t216 = Icges(5,4) * t228;
t139 = -Icges(5,1) * t229 - Icges(5,5) * t456 + t216;
t444 = -t136 * t228 + t139 * t229;
t520 = -t345 - t444;
t301 = -pkin(6) - qJ(3);
t269 = t301 * t455;
t297 = sin(pkin(9));
t458 = t297 * t302;
t402 = pkin(3) * t458;
t299 = cos(pkin(9));
t282 = pkin(3) * t299 + pkin(2);
t254 = pkin(4) * t286 + t282;
t490 = pkin(4) * t285;
t491 = pkin(3) * t297;
t367 = t490 + t491;
t424 = -t254 * t453 - t302 * t367;
t293 = -pkin(7) + t301;
t459 = t293 * t298;
t461 = t282 * t300;
t101 = t402 - t269 + (t459 + t461) * t303 + t424;
t354 = rSges(6,1) * t214 - rSges(6,2) * t213;
t124 = rSges(6,3) * t455 + t354;
t414 = t293 - t301;
t419 = t254 - t282;
t175 = t298 * t419 + t300 * t414;
t193 = -rSges(6,3) * t300 + (rSges(6,1) * t281 - rSges(6,2) * t280) * t298;
t247 = t296 * t455;
t277 = -qJD(4) * t300 + qJD(1);
t403 = qJD(4) * t303;
t391 = t298 * t403;
t519 = t101 * t277 - t124 * t255 + t175 * t391 + t193 * t247;
t445 = t137 * t228 + t141 * t229;
t45 = -t115 * t456 + t118 * t211 + t122 * t212;
t450 = t302 * t299;
t358 = rSges(4,2) * (t297 * t453 - t450) - rSges(4,1) * (t299 * t453 + t458);
t457 = t297 * t303;
t423 = (-t300 * t450 + t457) * rSges(4,1) + (t297 * t454 + t299 * t303) * rSges(4,2);
t356 = -rSges(5,1) * t231 + rSges(5,2) * t230;
t145 = rSges(5,3) * t455 - t356;
t198 = -rSges(5,3) * t300 + (rSges(5,1) * t286 - rSges(5,2) * t285) * t298;
t511 = -t145 * t277 + t198 * t391;
t64 = -t190 * t456 + t191 * t211 - t192 * t212;
t507 = t247 * t45 + t255 * t64;
t404 = qJD(4) * t302;
t406 = qJD(3) * t303;
t506 = t298 * (t198 * t404 + t406);
t505 = t298 * (t175 * t404 + t406);
t412 = qJD(1) * t302;
t396 = t298 * t412;
t410 = qJD(2) * t302;
t503 = t298 * t406 + t410;
t133 = -Icges(5,5) * t229 + Icges(5,6) * t228 - Icges(5,3) * t456;
t51 = t133 * t455 - t136 * t230 + t139 * t231;
t315 = t302 * (Icges(5,2) * t229 + t139 + t216) - t303 * (-Icges(5,2) * t231 - t141 - t217);
t316 = t302 * (-Icges(5,1) * t228 + t136 - t472) - t303 * (Icges(5,1) * t230 + t137 + t218);
t209 = (-Icges(6,2) * t281 - t468) * t298;
t502 = -t246 * (Icges(6,2) * t212 + t120 + t199) + t247 * (-Icges(6,2) * t214 - t122 - t200) + t255 * (t192 + t209);
t210 = (-Icges(6,1) * t280 - t467) * t298;
t501 = -t246 * (-Icges(6,1) * t211 + t117 - t469) + t247 * (Icges(6,1) * t213 + t118 + t201) + t255 * (t191 - t210);
t348 = qJD(1) * t378;
t237 = t302 * t348;
t500 = -t237 / 0.2e1;
t238 = t303 * t348;
t499 = -t238 / 0.2e1;
t498 = t246 / 0.2e1;
t497 = -t246 / 0.2e1;
t496 = -t247 / 0.2e1;
t495 = t247 / 0.2e1;
t493 = -t300 / 0.2e1;
t492 = pkin(2) * t300;
t489 = t303 * pkin(1);
t411 = qJD(1) * t303;
t395 = t298 * t411;
t130 = qJD(1) * t213 + t212 * t296;
t131 = -qJD(1) * t214 + t211 * t296;
t439 = rSges(6,1) * t131 + rSges(6,2) * t130;
t78 = -rSges(6,3) * t395 + t439;
t260 = t301 * t395;
t383 = t300 * t419;
t482 = pkin(4) * qJD(4);
t376 = t300 * t285 * t482;
t399 = t286 * t482;
t397 = t293 * t395 + t302 * t376 + t303 * t399;
t82 = -t260 + (-pkin(4) * t460 - t303 * t383) * qJD(1) + t397;
t488 = -t78 - t82;
t486 = rSges(3,1) * t300;
t483 = rSges(6,3) * t298;
t72 = Icges(6,5) * t131 + Icges(6,6) * t130 - Icges(6,3) * t395;
t74 = Icges(6,4) * t131 + Icges(6,2) * t130 - Icges(6,6) * t395;
t76 = Icges(6,1) * t131 + Icges(6,4) * t130 - Icges(6,5) * t395;
t27 = -t300 * t72 + ((-t117 * t296 + t76) * t281 + (-t120 * t296 - t74) * t280) * t298;
t479 = t27 * t246;
t195 = -Icges(5,3) * t300 + (Icges(5,5) * t286 - Icges(5,6) * t285) * t298;
t470 = Icges(5,4) * t286;
t196 = -Icges(5,6) * t300 + (-Icges(5,2) * t285 + t470) * t298;
t471 = Icges(5,4) * t285;
t197 = -Icges(5,5) * t300 + (Icges(5,1) * t286 - t471) * t298;
t70 = t195 * t455 - t196 * t230 + t197 * t231;
t478 = t277 * t70;
t128 = qJD(1) * t211 - t214 * t296;
t129 = -qJD(1) * t212 - t213 * t296;
t71 = Icges(6,5) * t129 + Icges(6,6) * t128 - Icges(6,3) * t396;
t73 = Icges(6,4) * t129 + Icges(6,2) * t128 - Icges(6,6) * t396;
t75 = Icges(6,1) * t129 + Icges(6,4) * t128 - Icges(6,5) * t396;
t28 = -t300 * t71 + ((-t118 * t296 + t75) * t281 + (t122 * t296 - t73) * t280) * t298;
t477 = t28 * t247;
t53 = -t114 * t300 + (-t117 * t280 + t120 * t281) * t298;
t476 = t53 * t238;
t475 = t54 * t237;
t474 = -rSges(3,3) - qJ(2);
t215 = (-rSges(6,1) * t280 - rSges(6,2) * t281) * t298;
t188 = t296 * t215;
t355 = -rSges(6,1) * t129 - rSges(6,2) * t128;
t77 = -rSges(6,3) * t396 - t355;
t473 = t188 * t455 + t300 * t77;
t466 = qJ(3) * t298;
t463 = t133 * t302;
t462 = t254 * t300;
t418 = pkin(3) * t457 + t301 * t456;
t420 = t293 * t456 + t303 * t367;
t100 = -t302 * t383 - t418 + t420;
t428 = -rSges(6,1) * t212 + rSges(6,2) * t211;
t123 = -rSges(6,3) * t456 + t428;
t447 = -t100 - t123;
t446 = t101 - t124;
t162 = qJD(1) * t230 + qJD(4) * t229;
t163 = -qJD(1) * t231 + qJD(4) * t228;
t434 = rSges(5,1) * t163 + rSges(5,2) * t162;
t433 = t188 * t456 + t193 * t395;
t227 = (-Icges(5,1) * t285 - t470) * t298;
t430 = t196 - t227;
t226 = (-Icges(5,2) * t286 - t471) * t298;
t429 = t197 + t226;
t427 = -rSges(5,1) * t229 + rSges(5,2) * t228;
t426 = t358 * qJD(1);
t288 = qJD(2) * t303;
t264 = t489 + t521;
t257 = qJD(1) * t264;
t415 = t288 - t257;
t425 = (t288 + t415) * qJD(1);
t290 = t303 * qJ(2);
t262 = -pkin(1) * t302 + t290;
t353 = t466 + t492;
t422 = -t302 * t353 + t262;
t249 = t353 * t303;
t421 = -t249 - t264;
t394 = t300 * t412;
t417 = pkin(2) * t394 + qJ(3) * t396;
t416 = rSges(3,2) * t456 + rSges(3,3) * t303;
t413 = qJD(1) * t298;
t409 = qJD(3) * t298;
t408 = qJD(3) * t300;
t407 = qJD(3) * t302;
t405 = qJD(4) * t298;
t401 = rSges(3,1) * t453;
t393 = t298 * t407;
t398 = qJD(1) * (-t353 * t411 - t393) + t425;
t390 = -t456 / 0.2e1;
t389 = t455 / 0.2e1;
t388 = -pkin(1) - t486;
t387 = -t413 / 0.2e1;
t386 = -t405 / 0.2e1;
t385 = t405 / 0.2e1;
t384 = -qJ(2) - t491;
t154 = rSges(6,1) * t211 + rSges(6,2) * t212;
t155 = -rSges(6,1) * t213 - rSges(6,2) * t214;
t382 = -t154 * t247 - t155 * t246;
t381 = t154 * t255 + t215 * t246;
t380 = -t155 * t255 + t215 * t247;
t379 = t298 ^ 2 * qJD(4) ^ 2 * t490;
t326 = t466 + (pkin(2) - t282) * t300;
t311 = t303 * t326 - t402;
t375 = qJD(1) * (qJD(1) * t311 + t260) + t398;
t374 = t302 * t387;
t373 = t303 * t387;
t372 = t302 * t386;
t371 = t302 * t385;
t370 = t303 * t386;
t369 = t303 * t385;
t368 = qJD(1) * t386;
t366 = t288 - t393;
t284 = pkin(1) * t412;
t243 = qJ(2) * t411 - t284 + t410;
t365 = t417 - t243 - t503;
t360 = -rSges(3,2) * t298 + t486;
t359 = t423 * qJD(1);
t160 = qJD(1) * t228 - qJD(4) * t231;
t319 = t230 * qJD(4);
t161 = -qJD(1) * t229 - t319;
t357 = rSges(5,1) * t161 + rSges(5,2) * t160;
t329 = (-rSges(5,1) * t285 - rSges(5,2) * t286) * t298;
t207 = qJD(4) * t329;
t252 = t282 * t394;
t337 = -qJD(1) * t418 + t252 + t365 - t417;
t89 = -rSges(5,3) * t396 + t357;
t42 = t207 * t391 - t277 * t89 + (t337 - t506) * qJD(1);
t90 = -rSges(5,3) * t395 + t434;
t43 = t277 * t90 + (t207 * t404 + (t198 * t403 - t407) * qJD(1)) * t298 + t375;
t352 = t302 * t43 + t303 * t42;
t49 = -t133 * t456 - t444;
t50 = -t134 * t456 + t445;
t351 = -t50 * t302 - t49 * t303;
t144 = -rSges(5,3) * t456 + t427;
t318 = t410 + (t302 * t326 + t418 + t422) * qJD(1);
t57 = t144 * t277 + t318 + t506;
t173 = t269 + t311;
t312 = (t173 + t421) * qJD(1) + t366;
t58 = t312 + t511;
t350 = t302 * t57 + t303 * t58;
t349 = -t302 * t89 - t303 * t90;
t347 = -t462 - t483;
t344 = -t144 * t303 - t145 * t302;
t343 = (Icges(5,5) * t228 + Icges(5,6) * t229) * t302 - (-Icges(5,5) * t230 - Icges(5,6) * t231) * t303;
t342 = t302 * t368;
t341 = t303 * t368;
t340 = -qJ(2) - t367;
t339 = -rSges(5,3) * t298 - pkin(1) - t461;
t338 = -pkin(1) + t347;
t336 = -rSges(3,3) * t302 - t401;
t328 = (-t302 * t49 + t303 * t50) * t298;
t52 = t134 * t455 + t345;
t327 = (-t302 * t51 + t303 * t52) * t298;
t81 = t252 + (t298 * t414 - t462) * t412 + (t285 * t411 - t319) * pkin(4);
t29 = -t303 * t379 + t188 * t247 - t193 * t237 - t255 * t77 - t277 * t81 + (t337 - t505) * qJD(1);
t39 = -t408 - t123 * t247 - t124 * t246 + (-t100 * t303 + t101 * t302) * t405;
t325 = t29 * (t124 * t300 + t193 * t455) + t39 * t123 * t396;
t225 = (-Icges(5,5) * t285 - Icges(5,6) * t286) * t298;
t208 = (-Icges(6,5) * t280 - Icges(6,6) * t281) * t298;
t322 = -qJD(1) * t249 - t257 + t366;
t321 = -(Icges(6,5) * t211 + Icges(6,6) * t212) * t246 + (-Icges(6,5) * t213 - Icges(6,6) * t214) * t247 + t208 * t255;
t320 = qJD(1) * t173 + t322;
t317 = -rSges(4,3) * t298 - pkin(1) - t353;
t314 = -rSges(4,3) * t455 + t358;
t11 = t117 * t128 + t120 * t129 - t213 * t74 + t214 * t76 + (-t114 * t412 + t303 * t72) * t298;
t12 = t118 * t128 - t122 * t129 - t213 * t73 + t214 * t75 + (-t115 * t412 + t303 * t71) * t298;
t13 = t117 * t130 + t120 * t131 + t211 * t74 - t212 * t76 + (-t114 * t411 - t302 * t72) * t298;
t14 = t118 * t130 - t122 * t131 + t211 * t73 - t212 * t75 + (-t115 * t411 - t302 * t71) * t298;
t44 = -t114 * t456 + t487;
t17 = -t246 * t44 + t507;
t47 = t115 * t455 - t346;
t18 = t247 * t47 - t523;
t185 = t296 * t208;
t186 = t296 * t209;
t187 = t296 * t210;
t34 = t128 * t191 + t129 * t192 - t186 * t213 + t187 * t214 + (t185 * t303 - t190 * t412) * t298;
t35 = t130 * t191 + t131 * t192 + t186 * t211 - t187 * t212 + (-t185 * t302 - t190 * t411) * t298;
t55 = -t185 * t300 + ((-t191 * t296 + t187) * t281 + (-t192 * t296 - t186) * t280) * t298;
t48 = t55 * t255;
t313 = (-t13 * t246 + t14 * t247 - t237 * t45 - t238 * t44 + t255 * t35) * t390 + t18 * t374 + t17 * t373 + (-t300 * t65 + (-t302 * t46 + t303 * t47) * t298) * t500 + (-t11 * t246 + t12 * t247 - t237 * t47 - t238 * t46 + t255 * t34) * t389 + (-t300 * t64 + (-t302 * t44 + t303 * t45) * t298) * t499 + (-t300 * t35 + (-t13 * t302 + t14 * t303 + (-t302 * t45 - t303 * t44) * qJD(1)) * t298) * t497 + (-t300 * t34 + (-t11 * t302 + t12 * t303 + (-t302 * t47 - t303 * t46) * qJD(1)) * t298) * t495 + (-t475 - t476 + t477 + t48 - t479) * t493 + t255 * (-t300 * t55 + (-t27 * t302 + t28 * t303 + (-t302 * t54 - t303 * t53) * qJD(1)) * t298) / 0.2e1 + (t211 * t502 + t212 * t501 - t321 * t456) * t498 + (-t213 * t502 - t214 * t501 + t321 * t455) * t496 - (-t321 * t300 + (-t280 * t502 - t281 * t501) * t298) * t255 / 0.2e1;
t83 = Icges(5,5) * t161 + Icges(5,6) * t160 - Icges(5,3) * t396;
t84 = Icges(5,5) * t163 + Icges(5,6) * t162 - Icges(5,3) * t395;
t85 = Icges(5,4) * t161 + Icges(5,2) * t160 - Icges(5,6) * t396;
t86 = Icges(5,4) * t163 + Icges(5,2) * t162 - Icges(5,6) * t395;
t87 = Icges(5,1) * t161 + Icges(5,4) * t160 - Icges(5,5) * t396;
t88 = Icges(5,1) * t163 + Icges(5,4) * t162 - Icges(5,5) * t395;
t308 = (-(t136 * t160 + t139 * t161 - t230 * t86 + t231 * t88 + (-t133 * t412 + t303 * t84) * t298) * t302 + (t137 * t160 - t141 * t161 - t230 * t85 + t231 * t87 + (-t134 * t412 + t303 * t83) * t298) * t303 + (-t52 * t302 - t51 * t303) * qJD(1)) * t298;
t307 = (qJD(1) * t351 - (t136 * t162 + t139 * t163 + t228 * t86 - t229 * t88 + (-t133 * t411 - t302 * t84) * t298) * t302 + (t137 * t162 - t141 * t163 + t228 * t85 - t229 * t87 + (-t134 * t411 - t302 * t83) * t298) * t303) * t298;
t31 = -t300 * t84 + (-t285 * t86 + t286 * t88 + (-t136 * t286 - t139 * t285) * qJD(4)) * t298;
t32 = -t300 * t83 + (-t285 * t85 + t286 * t87 + (-t137 * t286 + t141 * t285) * qJD(4)) * t298;
t59 = -t133 * t300 + (-t136 * t285 + t139 * t286) * t298;
t306 = (-t302 * t31 + t303 * t32 + (-t302 * t60 - t303 * t59) * qJD(1)) * t298;
t276 = rSges(3,2) * t455;
t271 = rSges(3,2) * t395;
t233 = -t276 - t336;
t220 = t230 * pkin(4);
t219 = t228 * pkin(4);
t206 = qJD(4) * t227;
t205 = qJD(4) * t226;
t204 = qJD(4) * t225;
t182 = t193 * t456;
t178 = t288 + (-t233 - t264) * qJD(1);
t177 = t410 + (-rSges(3,1) * t454 + t262 + t416) * qJD(1);
t171 = -rSges(5,1) * t230 - rSges(5,2) * t231;
t170 = rSges(5,1) * t228 + rSges(5,2) * t229;
t143 = (qJD(1) * t336 + t271) * qJD(1) + t425;
t142 = (-rSges(3,3) * t411 - t243 + (qJD(1) * t360 - qJD(2)) * t302) * qJD(1);
t93 = (t314 + t421) * qJD(1) + t366;
t92 = (-rSges(4,3) * t456 + t422 + t423) * qJD(1) + t503;
t80 = ((-rSges(4,3) * t411 - t407) * t298 + t426) * qJD(1) + t398;
t79 = ((rSges(4,3) * t412 - t406) * t298 - t359 + t365) * qJD(1);
t69 = -t195 * t456 + t196 * t228 - t197 * t229;
t68 = t344 * t405 - t408;
t66 = t69 * t277;
t62 = -t204 * t300 + (-t205 * t285 + t206 * t286 + (-t196 * t286 - t197 * t285) * qJD(4)) * t298;
t61 = t62 * t277;
t41 = t162 * t196 + t163 * t197 + t205 * t228 - t206 * t229 + (-t195 * t411 - t204 * t302) * t298;
t40 = t160 * t196 + t161 * t197 - t205 * t230 + t206 * t231 + (-t195 * t412 + t204 * t303) * t298;
t38 = t312 + t519;
t37 = t100 * t277 + t123 * t255 + t193 * t246 + t318 + t505;
t33 = ((t144 * t302 - t145 * t303) * qJD(1) + t349) * t405;
t30 = -t302 * t379 + t188 * t246 + t193 * t238 + t255 * t78 + t277 * t82 + (t175 * t403 - t407) * t413 + t375;
t26 = qJD(4) * t327 + t478;
t25 = qJD(4) * t328 + t66;
t10 = t123 * t237 - t124 * t238 - t246 * t77 - t247 * t78 + (-t302 * t81 - t303 * t82 + (t100 * t302 + t101 * t303) * qJD(1)) * t405;
t1 = [t48 + (t66 + (t445 * t303 + (t298 * t463 - t52 - t520) * t302) * t405) * t370 - t475 / 0.2e1 - t476 / 0.2e1 - t479 / 0.2e1 + t61 + t477 / 0.2e1 + (-(t47 + t524) * t246 + t507) * t496 + t35 * t497 + t64 * t499 + t65 * t500 + ((-t44 + t524) * t247 + t18 + t523) * t498 + (t34 + t17) * t495 + (-t478 + (-(-t445 - t51) * t302 + t520 * t303 + (-t134 * t302 ^ 2 + (-t134 * t303 - t463) * t303) * t298 + t351) * t405 + t26) * t371 + (t60 + t70) * t342 + (t59 + t69) * t341 + ((t302 * t338 + t290 + t420 + t428) * t30 + (-t354 + t424 + (-pkin(1) + t459 - t483) * t303 - t521) * t29 + (t284 + t355 + (t376 - t409) * t303 + (-qJD(2) - t399) * t302 + (t340 * t303 + (-t347 - t459) * t302) * qJD(1)) * t38 + (-t320 + t38 + t288 + t397 + t439 - t409 * t302 + (t302 * t340 + t303 * t338) * qJD(1) - t519) * t37) * m(6) + (-(t320 + t511 - t58) * t57 + t42 * (t269 + t356) + t58 * (t252 + t284 - t357) + t43 * (t290 + t418 + t427) + t57 * (t260 + t288 + t434) + (t42 * t339 - t58 * t409) * t303 + (-t58 * qJD(2) + t43 * t339 + t42 * t384 - t57 * t409) * t302 + ((t57 * t384 + t58 * (rSges(5,3) - t301) * t298) * t302 + (t339 * t57 + t384 * t58) * t303) * qJD(1)) * m(5) + (t79 * t358 + t93 * (t284 - t359 + t417) + t80 * (t290 + t423) + t92 * (t288 + t426) + (-t79 * qJ(2) + t93 * (rSges(4,3) * t413 - qJD(2)) + t80 * t317 + t92 * (-qJ(2) * qJD(1) - t409)) * t302 + (t79 * (-pkin(1) - t492) + (t79 * (-rSges(4,3) - qJ(3)) - t93 * qJD(3)) * t298 + (-t93 * qJ(2) + t317 * t92) * qJD(1)) * t303 - (qJD(1) * t314 + t322 - t93) * t92) * m(4) + (t142 * (t276 - t401 - t489) + t178 * t284 + t143 * (t290 + t416) + t177 * (t271 + t288) + (-t178 * qJD(2) + t142 * t474 + t143 * t388) * t302 + ((t177 * t388 + t178 * t474) * t303 + (t177 * t474 + t178 * t360) * t302) * qJD(1) - (-qJD(1) * t233 - t178 + t415) * t177) * m(3) + (t31 + t41) * t372 + (t32 + t40 + t25) * t369; m(3) * (t142 * t303 + t143 * t302) + m(4) * (t302 * t80 + t303 * t79) + m(5) * t352 + m(6) * (t29 * t303 + t30 * t302); 0.2e1 * (-m(5) * t33 / 0.2e1 - m(6) * t10 / 0.2e1) * t300 + 0.2e1 * (m(4) * (-t302 * t79 + t303 * t80) / 0.2e1 + m(5) * (-t42 * t302 + t303 * t43) / 0.2e1 + m(6) * (-t29 * t302 + t30 * t303) / 0.2e1) * t298; (qJD(4) * t306 + t61) * t493 + t277 * (-t300 * t62 + t306) / 0.2e1 + (-t300 * t69 + t328) * t341 + (-t300 * t70 + t327) * t342 + t26 * t374 + t25 * t373 + (-t300 * t41 + t307) * t372 + (-t300 * t40 + t308) * t369 + ((t225 * t455 - t230 * t429 - t231 * t430) * t277 + (t230 * t315 + t231 * t316 - t343 * t455) * t405) * t370 + ((-t225 * t456 + t228 * t429 + t229 * t430) * t277 + (-t228 * t315 - t229 * t316 + t343 * t456) * t405) * t371 + (qJD(4) * t307 + t277 * t41) * t390 + (qJD(4) * t308 + t277 * t40) * t389 - t277 * (-t300 * t225 * t277 + ((-t285 * t429 - t286 * t430) * t277 + ((t285 * t315 + t286 * t316) * t298 + t343 * t300) * qJD(4)) * t298) / 0.2e1 + t313 + (-t39 * ((-t219 * t303 + t220 * t302) * t405 + t382) - t38 * (t220 * t277 + t380) - t37 * (t219 * t277 + t381) + t38 * t473 + t30 * t182 + t37 * t433 + (-t101 * t29 + t30 * t447 + t37 * t488 + t38 * t81) * t300 + ((t10 * t447 + t39 * t488 + t29 * t175 + (t175 * t37 + t39 * t446) * qJD(1)) * t303 + (t10 * t446 + t39 * (-t77 - t81) + t30 * t175 + (t39 * t100 + t38 * (-t175 - t193)) * qJD(1)) * t302) * t298 + t325) * m(6) + ((-t144 * t43 + t145 * t42 - t57 * t90 + t58 * t89) * t300 + (t33 * t344 + t68 * (t144 * t412 - t145 * t411 + t349) + t350 * t207 + ((-t302 * t58 + t303 * t57) * qJD(1) + t352) * t198) * t298 - (t170 * t57 - t171 * t58) * t277 - (t68 * (-t170 * t303 - t171 * t302) + t350 * t329) * t405) * m(5); t313 + (t30 * (-t123 * t300 + t182) + (t10 * (-t123 * t303 - t124 * t302) + t39 * (-t124 * t411 - t302 * t77 - t303 * t78)) * t298 + t325 - t382 * t39 + (-t193 * t396 - t380 + t473) * t38 + (-t300 * t78 - t381 + t433) * t37) * m(6);];
tauc = t1(:);
