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
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:16:17
% EndTime: 2022-01-23 09:16:38
% DurationCPUTime: 13.66s
% Computational Cost: add. (19601->701), mult. (21719->967), div. (0->0), fcn. (20747->10), ass. (0->344)
t326 = pkin(9) + qJ(4);
t313 = sin(t326);
t314 = cos(t326);
t334 = cos(qJ(1));
t331 = cos(pkin(8));
t333 = sin(qJ(1));
t463 = t331 * t333;
t246 = t313 * t463 + t314 * t334;
t460 = t334 * t313;
t247 = t314 * t463 - t460;
t329 = sin(pkin(8));
t465 = t329 * t333;
t149 = Icges(5,5) * t247 - Icges(5,6) * t246 + Icges(5,3) * t465;
t231 = Icges(5,4) * t247;
t152 = -Icges(5,2) * t246 + Icges(5,6) * t465 + t231;
t230 = Icges(5,4) * t246;
t156 = -Icges(5,1) * t247 - Icges(5,5) * t465 + t230;
t63 = -t149 * t331 - (t152 * t313 + t156 * t314) * t329;
t315 = qJ(5) + t326;
t304 = sin(t315);
t305 = cos(t315);
t225 = t304 * t463 + t305 * t334;
t461 = t334 * t304;
t226 = t305 * t463 - t461;
t128 = Icges(6,5) * t226 - Icges(6,6) * t225 + Icges(6,3) * t465;
t214 = Icges(6,4) * t226;
t131 = -Icges(6,2) * t225 + Icges(6,6) * t465 + t214;
t213 = Icges(6,4) * t225;
t135 = -Icges(6,1) * t226 - Icges(6,5) * t465 + t213;
t56 = -t128 * t331 - (t131 * t304 + t135 * t305) * t329;
t283 = pkin(2) * t331 + qJ(3) * t329 + pkin(1);
t272 = t283 * t334;
t318 = t333 * qJ(2);
t323 = t334 * pkin(1);
t435 = t323 + t318;
t279 = qJD(1) * t435;
t430 = qJD(3) * t333;
t295 = t329 * t430;
t317 = qJD(2) * t334;
t437 = t317 - t295;
t392 = -qJD(1) * (-t323 + t272) - t279 + t437;
t444 = t272 + t318;
t330 = cos(pkin(9));
t306 = t330 * pkin(3) + pkin(2);
t332 = qJ(3) + pkin(6);
t262 = t306 * t331 + t329 * t332 + pkin(1);
t328 = sin(pkin(9));
t493 = pkin(3) * t328;
t303 = qJ(2) + t493;
t526 = t262 * t334 + t303 * t333;
t537 = qJD(1) * (-t444 + t526) - t392;
t327 = qJD(4) + qJD(5);
t396 = t327 * t329;
t267 = t333 * t396;
t268 = t334 * t396;
t275 = -t327 * t331 + qJD(1);
t227 = t305 * t333 - t331 * t461;
t462 = t331 * t334;
t228 = t304 * t333 + t305 * t462;
t464 = t329 * t334;
t48 = t128 * t464 + t227 * t131 - t135 * t228;
t130 = Icges(6,5) * t228 + Icges(6,6) * t227 + Icges(6,3) * t464;
t469 = Icges(6,4) * t228;
t133 = Icges(6,2) * t227 + Icges(6,6) * t464 + t469;
t215 = Icges(6,4) * t227;
t136 = Icges(6,1) * t228 + Icges(6,5) * t464 + t215;
t49 = t130 * t464 + t227 * t133 + t228 * t136;
t205 = -Icges(6,3) * t331 + (Icges(6,5) * t305 - Icges(6,6) * t304) * t329;
t467 = Icges(6,4) * t305;
t206 = -Icges(6,6) * t331 + (-Icges(6,2) * t304 + t467) * t329;
t468 = Icges(6,4) * t304;
t207 = -Icges(6,5) * t331 + (Icges(6,1) * t305 - t468) * t329;
t69 = t205 * t464 + t206 * t227 + t207 * t228;
t18 = t267 * t48 + t268 * t49 + t69 * t275;
t320 = t334 * qJ(2);
t399 = -t262 * t333 + t303 * t334;
t385 = -t320 + t399;
t518 = pkin(4) * t314;
t274 = t306 + t518;
t466 = t274 * t331;
t405 = pkin(1) + t466;
t492 = pkin(4) * t313;
t277 = t492 + t493;
t324 = -pkin(7) - t332;
t513 = t334 * t277 + t324 * t465;
t102 = t333 * t405 + t385 - t513;
t379 = rSges(6,1) * t226 - rSges(6,2) * t225;
t137 = rSges(6,3) * t465 + t379;
t190 = (t324 + t332) * t331 + (t274 - t306) * t329;
t427 = qJD(4) * t333;
t412 = t329 * t427;
t173 = t190 * t412;
t208 = -rSges(6,3) * t331 + (rSges(6,1) * t305 - rSges(6,2) * t304) * t329;
t428 = qJD(4) * t331;
t302 = qJD(1) - t428;
t536 = -t102 * t302 - t137 * t275 + t208 * t267 + t173;
t534 = -t152 * t246 - t156 * t247;
t248 = t314 * t333 - t331 * t460;
t249 = t313 * t333 + t314 * t462;
t533 = t248 * t152 - t156 * t249;
t46 = t128 * t465 - t131 * t225 - t135 * t226;
t527 = t149 * t464;
t269 = t333 * t277;
t271 = t283 * t333;
t369 = t271 + t385;
t494 = pkin(1) * t333;
t252 = -t271 + t494;
t296 = qJD(3) * t464;
t287 = -t320 + t494;
t316 = qJD(2) * t333;
t438 = -qJD(1) * t287 + t316;
t393 = qJD(1) * t252 + t296 + t438;
t441 = t296 + t316;
t525 = -qJD(1) * t369 - t393 + t441;
t524 = rSges(4,2) * (t328 * t463 + t330 * t334) - rSges(4,1) * (-t328 * t334 + t330 * t463);
t162 = t249 * rSges(5,1) + t248 * rSges(5,2) + rSges(5,3) * t464;
t523 = t162 * t302 + t537;
t522 = -(t328 * t333 + t330 * t462) * rSges(4,1) - (-t328 * t462 + t330 * t333) * rSges(4,2);
t433 = qJD(1) * t333;
t160 = rSges(5,1) * t247 - rSges(5,2) * t246 + rSges(5,3) * t465;
t212 = -rSges(5,3) * t331 + (rSges(5,1) * t314 - rSges(5,2) * t313) * t329;
t521 = -t160 * t302 + t212 * t412;
t139 = t228 * rSges(6,1) + t227 * rSges(6,2) + rSges(6,3) * t464;
t353 = t274 * t462 - t324 * t464 + t269 + t435;
t512 = t353 - t526;
t520 = t139 * t275 - t208 * t268 + t512 * t302 + t537;
t53 = t527 + t533;
t515 = t53 - t527;
t171 = -rSges(6,1) * t225 - rSges(6,2) * t226;
t229 = (-rSges(6,1) * t304 - rSges(6,2) * t305) * t329;
t204 = t327 * t229;
t432 = qJD(1) * t334;
t413 = t329 * t432;
t145 = qJD(1) * t227 - t226 * t327;
t146 = qJD(1) * t228 - t225 * t327;
t380 = rSges(6,1) * t146 + rSges(6,2) * t145;
t83 = rSges(6,3) * t413 + t380;
t511 = t171 * t275 + t204 * t465 + t208 * t413 - t267 * t229 + t331 * t83;
t370 = rSges(3,1) * t462 - rSges(3,2) * t464 + t333 * rSges(3,3);
t509 = qJD(1) * t370;
t47 = t130 * t465 - t225 * t133 + t226 * t136;
t151 = Icges(5,5) * t249 + Icges(5,6) * t248 + Icges(5,3) * t464;
t472 = Icges(5,4) * t249;
t154 = Icges(5,2) * t248 + Icges(5,6) * t464 + t472;
t232 = Icges(5,4) * t248;
t157 = Icges(5,1) * t249 + Icges(5,5) * t464 + t232;
t52 = t151 * t465 - t246 * t154 + t247 * t157;
t447 = t252 - t287;
t354 = (t369 + t447) * qJD(1) + t441;
t61 = t354 + t521;
t479 = t334 * t61;
t426 = qJD(4) * t334;
t411 = t329 * t426;
t62 = -t212 * t411 + t523;
t507 = (t333 * t62 + t479) * qJD(1);
t506 = t333 * (-Icges(5,2) * t247 - t156 - t230) + t334 * (-Icges(5,2) * t249 + t157 + t232);
t223 = (-Icges(6,2) * t305 - t468) * t329;
t505 = t267 * (-Icges(6,2) * t226 - t135 - t213) + t268 * (-Icges(6,2) * t228 + t136 + t215) + t275 * (t207 + t223);
t414 = t329 * t433;
t257 = t327 * t414;
t504 = -t257 / 0.2e1;
t258 = qJD(1) * t268;
t503 = t258 / 0.2e1;
t502 = -t267 / 0.2e1;
t501 = t267 / 0.2e1;
t499 = t268 / 0.2e1;
t497 = -t331 / 0.2e1;
t496 = t333 / 0.2e1;
t495 = -t334 / 0.2e1;
t143 = qJD(1) * t225 - t228 * t327;
t144 = -qJD(1) * t226 + t227 * t327;
t456 = t144 * rSges(6,1) + t143 * rSges(6,2);
t82 = -rSges(6,3) * t414 + t456;
t281 = t303 * t432;
t310 = qJ(2) * t432;
t421 = pkin(4) * t426;
t343 = -t313 * t331 * t421 + t277 * t432 + t324 * t414 + t427 * t518 + t310;
t398 = t262 - t466;
t86 = -t281 + (-pkin(1) + t398) * t433 + t343;
t491 = -t82 - t86;
t490 = rSges(3,1) * t331;
t77 = Icges(6,5) * t146 + Icges(6,6) * t145 + Icges(6,3) * t413;
t79 = Icges(6,4) * t146 + Icges(6,2) * t145 + Icges(6,6) * t413;
t81 = Icges(6,1) * t146 + Icges(6,4) * t145 + Icges(6,5) * t413;
t29 = -t331 * t77 + ((-t131 * t327 + t81) * t305 + (t135 * t327 - t79) * t304) * t329;
t484 = t29 * t267;
t76 = Icges(6,5) * t144 + Icges(6,6) * t143 - Icges(6,3) * t414;
t78 = Icges(6,4) * t144 + Icges(6,2) * t143 - Icges(6,6) * t414;
t80 = Icges(6,1) * t144 + Icges(6,4) * t143 - Icges(6,5) * t414;
t30 = -t331 * t76 + ((-t133 * t327 + t80) * t305 + (-t136 * t327 - t78) * t304) * t329;
t483 = t30 * t268;
t209 = -Icges(5,3) * t331 + (Icges(5,5) * t314 - Icges(5,6) * t313) * t329;
t470 = Icges(5,4) * t314;
t210 = -Icges(5,6) * t331 + (-Icges(5,2) * t313 + t470) * t329;
t471 = Icges(5,4) * t313;
t211 = -Icges(5,5) * t331 + (Icges(5,1) * t314 - t471) * t329;
t74 = t209 * t465 - t210 * t246 + t211 * t247;
t482 = t302 * t74;
t51 = t149 * t465 + t534;
t481 = t333 * t51;
t36 = t354 + t536;
t478 = t36 * t190;
t37 = -t190 * t411 + t520;
t477 = t37 * t204;
t361 = (-rSges(5,1) * t313 - rSges(5,2) * t314) * t329;
t221 = qJD(4) * t361;
t425 = qJD(1) * qJD(2);
t439 = t310 + t316;
t446 = qJD(1) * (-pkin(1) * t433 + t439) + t333 * t425;
t394 = t446 + (0.2e1 * t296 + (pkin(1) - t283) * t433) * qJD(1);
t445 = t262 - t283;
t376 = qJD(1) * (-t433 * t445 + t281 - t310) + t394;
t429 = qJD(4) * t329;
t176 = qJD(1) * t246 - qJD(4) * t249;
t177 = -qJD(1) * t247 + qJD(4) * t248;
t453 = t177 * rSges(5,1) + t176 * rSges(5,2);
t94 = -rSges(5,3) * t414 + t453;
t42 = t302 * t94 + (t212 * t433 - t221 * t334) * t429 + t376;
t476 = t42 * t334;
t308 = t334 * t425;
t282 = t303 * t433;
t309 = qJ(2) * t433;
t443 = t282 - t309;
t312 = pkin(1) * t432;
t436 = t317 - t309;
t263 = t312 - t436;
t452 = -t283 * t432 - t263 - t295 + t312;
t418 = -t432 * t445 - t443 + t452;
t178 = qJD(1) * t248 - qJD(4) * t247;
t351 = t246 * qJD(4);
t179 = qJD(1) * t249 - t351;
t381 = rSges(5,1) * t179 + rSges(5,2) * t178;
t95 = rSges(5,3) * t413 + t381;
t43 = t221 * t412 - t302 * t95 + t308 + ((t212 * t426 - t430) * t329 + t418) * qJD(1);
t475 = t43 * t333;
t474 = t56 * t258;
t57 = -t130 * t331 + (-t133 * t304 + t136 * t305) * t329;
t473 = t57 * t257;
t459 = -t512 - t139;
t245 = (-Icges(5,1) * t313 - t470) * t329;
t450 = -t210 + t245;
t244 = (-Icges(5,2) * t314 - t471) * t329;
t449 = t211 + t244;
t448 = t524 * qJD(1);
t442 = rSges(3,2) * t414 + rSges(3,3) * t432;
t440 = rSges(3,2) * t465 + t334 * rSges(3,3);
t431 = qJD(3) * t331;
t424 = rSges(3,1) * t463;
t54 = t151 * t464 + t248 * t154 + t249 * t157;
t416 = rSges(4,3) * t464 - t522;
t415 = -t295 + t436;
t410 = t465 / 0.2e1;
t409 = t464 / 0.2e1;
t408 = -pkin(1) - t490;
t407 = -t429 / 0.2e1;
t406 = t429 / 0.2e1;
t172 = rSges(6,1) * t227 - rSges(6,2) * t228;
t402 = t268 * t171 - t172 * t267;
t401 = t275 * t172 - t229 * t268;
t397 = t329 ^ 2 * qJD(4) ^ 2 * t492;
t391 = -t414 / 0.2e1;
t390 = qJD(1) * t409;
t389 = t333 * t407;
t388 = t333 * t406;
t387 = t334 * t407;
t386 = t334 * t406;
t382 = t522 * qJD(1);
t378 = t333 * t61 - t334 * t62;
t377 = -t333 * t94 + t334 * t95;
t374 = t160 * t334 - t162 * t333;
t373 = (-Icges(5,5) * t246 - Icges(5,6) * t247) * t333 + (Icges(5,5) * t248 - Icges(5,6) * t249) * t334;
t372 = qJD(1) * t389;
t371 = qJD(1) * t386;
t368 = -rSges(6,3) * t329 - t405;
t360 = (t334 * t52 + t481) * t329;
t359 = (t333 * t53 + t334 * t54) * t329;
t355 = t329 * t373;
t224 = (-Icges(6,1) * t304 - t467) * t329;
t243 = (-Icges(5,5) * t313 - Icges(5,6) * t314) * t329;
t222 = (-Icges(6,5) * t304 - Icges(6,6) * t305) * t329;
t352 = (-Icges(6,5) * t225 - Icges(6,6) * t226) * t267 + (Icges(6,5) * t227 - Icges(6,6) * t228) * t268 + t222 * t275;
t350 = (Icges(5,1) * t248 - t154 - t472) * t334 + (-Icges(5,1) * t246 - t152 - t231) * t333;
t348 = -rSges(4,3) * t465 + t524;
t87 = t312 - pkin(4) * t351 + (t269 + (-t324 * t329 - t398) * t334) * qJD(1) - t443;
t28 = -t333 * t397 + t204 * t267 + t208 * t258 - t275 * t83 - t302 * t87 + t308 + ((t190 * t426 - t430) * t329 + t418) * qJD(1);
t41 = -t431 + t137 * t268 - t139 * t267 + (t102 * t334 - t333 * t512) * t429;
t9 = -t137 * t257 - t139 * t258 - t267 * t82 + t268 * t83 + (-t333 * t86 + t334 * t87 + (-t102 * t333 - t334 * t512) * qJD(1)) * t429;
t347 = t28 * (t331 * t137 + t208 * t465) + (t137 * t9 + t41 * t83) * t464;
t68 = t205 * t465 - t206 * t225 + t207 * t226;
t344 = t329 * t352;
t11 = t131 * t143 - t135 * t144 + t227 * t79 + t228 * t81 + (-t128 * t433 + t334 * t77) * t329;
t12 = t133 * t143 + t136 * t144 + t227 * t78 + t228 * t80 + (-t130 * t433 + t334 * t76) * t329;
t13 = t131 * t145 - t135 * t146 - t225 * t79 + t226 * t81 + (t128 * t432 + t333 * t77) * t329;
t14 = t133 * t145 + t136 * t146 - t225 * t78 + t226 * t80 + (t130 * t432 + t333 * t76) * t329;
t341 = (Icges(6,1) * t227 - t133 - t469) * t268 + (-Icges(6,1) * t225 - t131 - t214) * t267 + (-t206 + t224) * t275;
t201 = t327 * t222;
t202 = t327 * t223;
t203 = t327 * t224;
t39 = t143 * t206 + t144 * t207 + t202 * t227 + t203 * t228 + (t201 * t334 - t205 * t433) * t329;
t40 = t145 * t206 + t146 * t207 - t202 * t225 + t203 * t226 + (t201 * t333 + t205 * t432) * t329;
t60 = -t201 * t331 + ((-t206 * t327 + t203) * t305 + (-t207 * t327 - t202) * t304) * t329;
t50 = t60 * t275;
t342 = (t13 * t267 + t14 * t268 - t257 * t47 + t258 * t46 + t275 * t40) * t410 - (-t352 * t331 + (-t304 * t505 + t305 * t341) * t329) * t275 / 0.2e1 + t18 * t391 + (t267 * t46 + t268 * t47 + t275 * t68) * t390 + (-t331 * t69 + (t333 * t48 + t334 * t49) * t329) * t504 + (t11 * t267 + t12 * t268 - t257 * t49 + t258 * t48 + t275 * t39) * t409 + (-t331 * t68 + (t333 * t46 + t334 * t47) * t329) * t503 + (-t331 * t40 + (t13 * t333 + t14 * t334 + (-t333 * t47 + t334 * t46) * qJD(1)) * t329) * t501 + (-t331 * t39 + (t11 * t333 + t12 * t334 + (-t333 * t49 + t334 * t48) * qJD(1)) * t329) * t499 + (-t473 + t474 + t483 + t50 + t484) * t497 + t275 * (-t331 * t60 + (t29 * t333 + t30 * t334 + (-t333 * t57 + t334 * t56) * qJD(1)) * t329) / 0.2e1 + (-t225 * t505 + t226 * t341 + t333 * t344) * t502 - (t227 * t505 + t341 * t228 + t334 * t344) * t268 / 0.2e1;
t88 = Icges(5,5) * t177 + Icges(5,6) * t176 - Icges(5,3) * t414;
t89 = Icges(5,5) * t179 + Icges(5,6) * t178 + Icges(5,3) * t413;
t90 = Icges(5,4) * t177 + Icges(5,2) * t176 - Icges(5,6) * t414;
t91 = Icges(5,4) * t179 + Icges(5,2) * t178 + Icges(5,6) * t413;
t92 = Icges(5,1) * t177 + Icges(5,4) * t176 - Icges(5,5) * t414;
t93 = Icges(5,1) * t179 + Icges(5,4) * t178 + Icges(5,5) * t413;
t339 = ((t152 * t176 - t156 * t177 + t248 * t91 + t249 * t93 + (-t149 * t433 + t334 * t89) * t329) * t333 + (t154 * t176 + t157 * t177 + t248 * t90 + t249 * t92 + (-t151 * t433 + t334 * t88) * t329) * t334 + (-t54 * t333 + t53 * t334) * qJD(1)) * t329;
t338 = ((t152 * t178 - t156 * t179 - t246 * t91 + t247 * t93 + (t149 * t432 + t333 * t89) * t329) * t333 + (t154 * t178 + t157 * t179 - t246 * t90 + t247 * t92 + (t151 * t432 + t333 * t88) * t329) * t334 + (-t52 * t333 + t51 * t334) * qJD(1)) * t329;
t31 = -t331 * t89 + (-t313 * t91 + t314 * t93 + (-t152 * t314 + t156 * t313) * qJD(4)) * t329;
t32 = -t331 * t88 + (-t313 * t90 + t314 * t92 + (-t154 * t314 - t157 * t313) * qJD(4)) * t329;
t64 = -t151 * t331 + (-t154 * t313 + t157 * t314) * t329;
t337 = (t31 * t333 + t32 * t334 + (-t333 * t64 + t63 * t334) * qJD(1)) * t329;
t251 = t424 - t440;
t234 = t248 * pkin(4);
t233 = t246 * pkin(4);
t220 = qJD(4) * t245;
t219 = qJD(4) * t244;
t218 = qJD(4) * t243;
t195 = t208 * t414;
t193 = t279 - t317 + t509;
t192 = t316 + (-t251 - t287) * qJD(1);
t187 = rSges(5,1) * t248 - rSges(5,2) * t249;
t186 = -rSges(5,1) * t246 - rSges(5,2) * t247;
t159 = t308 + (-t263 - t509) * qJD(1);
t158 = qJD(1) * (-qJD(1) * t424 + t442) + t446;
t100 = qJD(1) * t416 - t392;
t99 = (t348 + t447) * qJD(1) + t441;
t85 = t308 + ((-rSges(4,3) * t432 - t430) * t329 + t382 + t452) * qJD(1);
t84 = qJD(1) * (-rSges(4,3) * t414 + t448) + t394;
t75 = t209 * t464 + t210 * t248 + t211 * t249;
t73 = t374 * t429 - t431;
t70 = t75 * t302;
t66 = -t218 * t331 + (-t219 * t313 + t220 * t314 + (-t210 * t314 - t211 * t313) * qJD(4)) * t329;
t65 = t66 * t302;
t45 = t178 * t210 + t179 * t211 - t219 * t246 + t220 * t247 + (t209 * t432 + t218 * t333) * t329;
t44 = t176 * t210 + t177 * t211 + t219 * t248 + t220 * t249 + (-t209 * t433 + t218 * t334) * t329;
t38 = ((-t160 * t333 - t162 * t334) * qJD(1) + t377) * t429;
t27 = qJD(1) * t173 - t204 * t268 + t208 * t257 + t275 * t82 + t302 * t86 + t334 * t397 + t376;
t26 = qJD(4) * t359 + t70;
t25 = qJD(4) * t360 + t482;
t1 = [t50 + t483 / 0.2e1 + t484 / 0.2e1 - t473 / 0.2e1 + t474 / 0.2e1 + t65 + t39 * t499 + t68 * t503 + t69 * t504 + t18 * t502 + (t70 + ((t51 - t534 + t54) * t334 + t515 * t333) * t429) * t389 + (t40 + t18) * t501 + (t32 + t44) * t386 + (t64 + t75) * t372 + (t63 + t74) * t371 + (-t478 * t334 * t429 + t27 * (t353 + t139) + (t368 * t333 + t320 - t379 + t513) * t28 + (t368 * t433 + t343 + t456 + t525 - t536) * t37 + (t314 * t421 - t312 - t380 + t415 + t428 * t492 * t333 + ((-t466 + (-rSges(6,3) + t324) * t329) * t334 - t269) * qJD(1) + t520) * t36) * m(6) + (-t212 * t479 * t429 + t43 * (-t160 + t399) + t42 * (t162 + t526) + (-rSges(5,3) * t329 - t262) * t507 + (t281 + t453 - t521 + t525) * t62 + (-t282 - t381 + t437 + t523) * t61) * m(5) + (t85 * (-t271 + t320 + t348) + t99 * (t382 + t415) + t84 * (t416 + t444) + t100 * (t296 + t439 + t448) + (t100 * t333 + t334 * t99) * qJD(1) * (-rSges(4,3) * t329 - t283) - (qJD(1) * t348 + t393 - t99) * t100) * m(4) + (-(-qJD(1) * t251 - t192 + t438) * t193 + t159 * (t333 * t408 + t320 + t440) - t192 * t263 + t158 * (t370 + t435) + t193 * (t439 + t442) + (t192 * (rSges(3,2) * t329 - t490) * t334 + (-t192 * rSges(3,3) + t193 * t408) * t333) * qJD(1)) * m(3) + (t31 + t45 + t26) * t388 + (t25 - t482 + ((t515 - t52 - t533) * t334 - t481) * t429) * t387; 0.2e1 * (t27 * t495 + t28 * t496) * m(6) + 0.2e1 * (t475 / 0.2e1 - t476 / 0.2e1) * m(5) + 0.2e1 * (t495 * t84 + t496 * t85) * m(4) + 0.2e1 * (t158 * t495 + t159 * t496) * m(3); 0.2e1 * (-m(5) * t38 / 0.2e1 - m(6) * t9 / 0.2e1) * t331 + 0.2e1 * (m(4) * (t333 * t84 + t334 * t85) / 0.2e1 + m(5) * (t333 * t42 + t334 * t43) / 0.2e1 + m(6) * (t27 * t333 + t28 * t334) / 0.2e1) * t329; ((t243 * t465 - t246 * t449 + t247 * t450) * t302 + (-t246 * t506 + t247 * t350 + t333 * t355) * t429) * t389 + t25 * t390 + (-t331 * t75 + t359) * t372 + (-t331 * t74 + t360) * t371 + ((t243 * t464 + t248 * t449 + t249 * t450) * t302 + (t248 * t506 + t350 * t249 + t334 * t355) * t429) * t387 + t26 * t391 - t302 * (-t331 * t243 * t302 + ((-t313 * t449 + t314 * t450) * t302 + ((-t313 * t506 + t314 * t350) * t329 - t373 * t331) * qJD(4)) * t329) / 0.2e1 + (qJD(4) * t339 + t302 * t44) * t409 + (qJD(4) * t338 + t302 * t45) * t410 + t342 + (-t331 * t45 + t338) * t388 + (-t331 * t44 + t339) * t386 + (qJD(4) * t337 + t65) * t497 + t302 * (-t331 * t66 + t337) / 0.2e1 + (-t37 * (t234 * t302 + t401) - t41 * ((-t233 * t334 - t234 * t333) * t429 + t402) + t37 * t195 + (t28 * t102 + t27 * t459 + t37 * t491) * t331 + ((t27 * (-t190 - t208) - t477 + t9 * t102 + t41 * t87 + (t41 * t459 + t478) * qJD(1)) * t334 + (t28 * t190 + t9 * t459 + t41 * t491 + (t37 * t190 + t41 * (-t102 - t137)) * qJD(1)) * t333) * t329 + t347 + (-t233 * t302 + t331 * t87 + t511) * t36) * m(6) + ((t160 * t43 - t162 * t42 + t61 * t95 - t62 * t94) * t331 + (t38 * t374 + t73 * (-t160 * t433 - t162 * t432 + t377) + t378 * t221 + (t475 - t476 + t507) * t212) * t329 - (-t186 * t61 + t187 * t62) * t302 - (t73 * (t186 * t334 - t187 * t333) + t378 * t361) * t429) * m(5); t342 + (-t27 * t331 * t139 + ((-qJD(1) * t139 * t41 - t208 * t27 - t477) * t334 + (-t9 * t139 + t41 * (-qJD(1) * t137 - t82)) * t333) * t329 + t347 - t41 * t402 + (-t331 * t82 + t195 - t401) * t37 + t511 * t36) * m(6);];
tauc = t1(:);
