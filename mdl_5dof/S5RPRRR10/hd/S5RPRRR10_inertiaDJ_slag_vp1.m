% Calculate time derivative of joint inertia matrix for
% S5RPRRR10
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR10_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR10_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR10_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR10_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:09:44
% EndTime: 2019-12-31 19:10:07
% DurationCPUTime: 11.67s
% Computational Cost: add. (34021->777), mult. (37124->1101), div. (0->0), fcn. (35898->10), ass. (0->397)
t290 = pkin(9) + qJ(3);
t283 = sin(t290);
t498 = -qJD(1) * t283 / 0.2e1;
t299 = cos(qJ(1));
t284 = cos(t290);
t411 = qJD(3) * t284;
t374 = t411 / 0.2e1;
t297 = sin(qJ(1));
t414 = qJD(1) * t297;
t386 = t283 * t414;
t497 = t299 * t374 - t386 / 0.2e1;
t413 = qJD(1) * t299;
t375 = t413 / 0.2e1;
t496 = t283 * t375 + t297 * t374;
t410 = qJD(3) * t297;
t379 = t284 * t410;
t309 = t283 * t413 + t379;
t300 = -pkin(8) - pkin(7);
t471 = pkin(7) + t300;
t495 = t284 * t471;
t465 = Icges(4,4) * t283;
t343 = Icges(4,1) * t284 - t465;
t232 = Icges(4,5) * t297 + t299 * t343;
t443 = t232 * t284;
t464 = Icges(4,4) * t284;
t339 = -Icges(4,2) * t283 + t464;
t230 = Icges(4,6) * t297 + t299 * t339;
t448 = t230 * t283;
t323 = -t443 + t448;
t494 = t297 * t323;
t231 = -Icges(4,5) * t299 + t297 * t343;
t445 = t231 * t284;
t229 = -Icges(4,6) * t299 + t297 * t339;
t450 = t229 * t283;
t324 = -t445 + t450;
t493 = t299 * t324;
t292 = qJ(4) + qJ(5);
t287 = sin(t292);
t288 = cos(t292);
t333 = Icges(6,5) * t288 - Icges(6,6) * t287;
t291 = qJD(4) + qJD(5);
t441 = t283 * t291;
t158 = (-Icges(6,5) * t287 - Icges(6,6) * t288) * t441 + (Icges(6,3) * t283 + t284 * t333) * qJD(3);
t460 = Icges(6,4) * t288;
t336 = -Icges(6,2) * t287 + t460;
t218 = -Icges(6,6) * t284 + t283 * t336;
t453 = t218 * t287;
t492 = -qJD(3) * t453 - t158;
t296 = sin(qJ(4));
t298 = cos(qJ(4));
t334 = Icges(5,5) * t298 - Icges(5,6) * t296;
t407 = qJD(4) * t283;
t175 = (-Icges(5,5) * t296 - Icges(5,6) * t298) * t407 + (Icges(5,3) * t283 + t284 * t334) * qJD(3);
t462 = Icges(5,4) * t298;
t337 = -Icges(5,2) * t296 + t462;
t224 = -Icges(5,6) * t284 + t283 * t337;
t451 = t224 * t296;
t491 = -qJD(3) * t451 - t175;
t289 = t297 * rSges(4,3);
t439 = t283 * t299;
t490 = -rSges(4,2) * t439 + t289;
t295 = -pkin(6) - qJ(2);
t294 = cos(pkin(9));
t280 = pkin(2) * t294 + pkin(1);
t470 = rSges(4,1) * t284;
t357 = -rSges(4,2) * t283 + t470;
t316 = -t280 - t357;
t210 = (rSges(4,3) - t295) * t299 + t316 * t297;
t438 = t284 * t299;
t234 = rSges(4,1) * t438 + t490;
t372 = t299 * t280 - t297 * t295;
t211 = t234 + t372;
t489 = t210 * t299 + t211 * t297;
t335 = Icges(4,5) * t284 - Icges(4,6) * t283;
t227 = -Icges(4,3) * t299 + t297 * t335;
t415 = qJD(1) * t284;
t364 = -qJD(4) + t415;
t408 = qJD(3) * t299;
t381 = t283 * t408;
t488 = t297 * t364 + t381;
t369 = -t291 + t415;
t487 = t297 * t369 + t381;
t382 = t283 * t410;
t486 = t299 * t369 - t382;
t318 = rSges(3,1) * t294 - rSges(3,2) * sin(pkin(9)) + pkin(1);
t468 = rSges(3,3) + qJ(2);
t222 = t297 * t468 + t299 * t318;
t435 = t288 * t297;
t436 = t287 * t299;
t242 = -t284 * t436 + t435;
t434 = t288 * t299;
t437 = t287 * t297;
t243 = t284 * t434 + t437;
t172 = t243 * rSges(6,1) + t242 * rSges(6,2) + rSges(6,3) * t439;
t275 = pkin(3) * t438;
t247 = pkin(7) * t439 + t275;
t432 = t296 * t297;
t279 = pkin(4) * t432;
t281 = pkin(4) * t298 + pkin(3);
t317 = t281 * t438 - t300 * t439 + t279;
t184 = t317 - t247;
t423 = t172 + t184;
t240 = -t284 * t437 - t434;
t241 = t284 * t435 - t436;
t352 = -t241 * rSges(6,1) - t240 * rSges(6,2);
t440 = t283 * t297;
t171 = rSges(6,3) * t440 - t352;
t472 = pkin(3) - t281;
t306 = -t283 * t471 - t284 * t472;
t431 = t296 * t299;
t404 = pkin(4) * t431;
t183 = t297 * t306 - t404;
t424 = t171 + t183;
t485 = -t297 * t424 - t299 * t423;
t484 = 2 * m(4);
t483 = 2 * m(5);
t482 = 2 * m(6);
t481 = t297 ^ 2;
t480 = t299 ^ 2;
t479 = -t284 / 0.2e1;
t478 = t297 / 0.2e1;
t477 = t299 / 0.2e1;
t476 = -rSges(5,3) - pkin(7);
t261 = rSges(4,1) * t283 + rSges(4,2) * t284;
t475 = m(4) * t261;
t474 = pkin(3) * t284;
t473 = qJD(1) / 0.2e1;
t469 = rSges(6,3) * t283;
t467 = -rSges(6,3) + t300;
t370 = -t284 * t291 + qJD(1);
t322 = t297 * t370;
t154 = -t287 * t486 + t288 * t322;
t155 = t287 * t322 + t288 * t486;
t353 = t155 * rSges(6,1) + t154 * rSges(6,2);
t104 = rSges(6,3) * t309 + t353;
t378 = t284 * t408;
t466 = t104 * t439 + t171 * t378;
t463 = Icges(5,4) * t296;
t461 = Icges(6,4) * t287;
t429 = t298 * t299;
t249 = -t284 * t432 - t429;
t430 = t297 * t298;
t250 = t284 * t430 - t431;
t355 = -rSges(5,1) * t250 - rSges(5,2) * t249;
t198 = rSges(5,3) * t440 - t355;
t456 = t198 * t299;
t340 = Icges(6,1) * t288 - t461;
t219 = -Icges(6,5) * t284 + t283 * t340;
t452 = t219 * t288;
t449 = t229 * t284;
t447 = t230 * t284;
t446 = t231 * t283;
t444 = t232 * t283;
t442 = t283 * t281;
t176 = (-Icges(5,2) * t298 - t463) * t407 + (Icges(5,6) * t283 + t284 * t337) * qJD(3);
t433 = t296 * t176;
t428 = t299 * t295;
t321 = t370 * t299;
t152 = t287 * t487 + t288 * t321;
t153 = t287 * t321 - t288 * t487;
t394 = t153 * rSges(6,1) + t152 * rSges(6,2) + rSges(6,3) * t378;
t103 = -rSges(6,3) * t386 + t394;
t269 = pkin(7) * t378;
t405 = qJD(4) * t298;
t401 = pkin(4) * t405;
t388 = qJD(1) * t404 + t297 * t401 + t300 * t386;
t406 = qJD(4) * t296;
t402 = pkin(4) * t406;
t427 = t103 - t269 + (pkin(7) * t414 + t408 * t472) * t283 + ((-qJD(3) * t300 - t402) * t299 + t472 * t414) * t284 + t388;
t351 = rSges(6,1) * t288 - rSges(6,2) * t287;
t161 = (-rSges(6,1) * t287 - rSges(6,2) * t288) * t441 + (t284 * t351 + t469) * qJD(3);
t380 = t283 * t406;
t191 = -pkin(4) * t380 + qJD(3) * t306;
t426 = -t161 - t191;
t220 = -rSges(6,3) * t284 + t283 * t351;
t412 = qJD(3) * t283;
t425 = t172 * t412 + t220 * t386;
t354 = rSges(5,1) * t298 - rSges(5,2) * t296;
t182 = (-rSges(5,1) * t296 - rSges(5,2) * t298) * t407 + (rSges(5,3) * t283 + t284 * t354) * qJD(3);
t359 = pkin(7) * t283 + t474;
t257 = t359 * qJD(3);
t422 = -t182 - t257;
t251 = -t284 * t431 + t430;
t252 = t284 * t429 + t432;
t199 = t252 * rSges(5,1) + t251 * rSges(5,2) + rSges(5,3) * t439;
t421 = -t199 - t247;
t134 = t284 * t171 + t220 * t440;
t216 = -t283 * t472 + t495;
t420 = t216 + t220;
t226 = -rSges(5,3) * t284 + t283 * t354;
t263 = t283 * pkin(3) - t284 * pkin(7);
t419 = -t226 - t263;
t246 = t359 * t297;
t418 = t297 * t246 + t299 * t247;
t286 = qJD(2) * t299;
t417 = t295 * t414 + t286;
t228 = Icges(4,3) * t297 + t299 * t335;
t416 = qJD(1) * t228;
t409 = qJD(3) * t298;
t365 = -qJD(4) * t284 + qJD(1);
t320 = t298 * t365;
t180 = t297 * t320 + (-t299 * t364 + t382) * t296;
t319 = t365 * t296;
t181 = t364 * t429 + (-t283 * t409 + t319) * t297;
t113 = Icges(5,5) * t181 + Icges(5,6) * t180 + Icges(5,3) * t309;
t115 = Icges(5,4) * t181 + Icges(5,2) * t180 + Icges(5,6) * t309;
t117 = Icges(5,1) * t181 + Icges(5,4) * t180 + Icges(5,5) * t309;
t192 = Icges(5,5) * t250 + Icges(5,6) * t249 + Icges(5,3) * t440;
t194 = Icges(5,4) * t250 + Icges(5,2) * t249 + Icges(5,6) * t440;
t196 = Icges(5,1) * t250 + Icges(5,4) * t249 + Icges(5,5) * t440;
t328 = -t194 * t296 + t196 * t298;
t33 = (qJD(3) * t328 - t113) * t284 + (qJD(3) * t192 - t115 * t296 + t117 * t298 + (-t194 * t298 - t196 * t296) * qJD(4)) * t283;
t341 = Icges(5,1) * t298 - t463;
t177 = (-Icges(5,1) * t296 - t462) * t407 + (Icges(5,5) * t283 + t284 * t341) * qJD(3);
t223 = -Icges(5,3) * t284 + t283 * t334;
t225 = -Icges(5,5) * t284 + t283 * t341;
t59 = t175 * t440 + t176 * t249 + t177 * t250 + t180 * t224 + t181 * t225 + t223 * t309;
t400 = t33 / 0.2e1 + t59 / 0.2e1;
t178 = t296 * t488 + t299 * t320;
t179 = -t298 * t488 + t299 * t319;
t308 = t378 - t386;
t112 = Icges(5,5) * t179 + Icges(5,6) * t178 + Icges(5,3) * t308;
t114 = Icges(5,4) * t179 + Icges(5,2) * t178 + Icges(5,6) * t308;
t116 = Icges(5,1) * t179 + Icges(5,4) * t178 + Icges(5,5) * t308;
t193 = Icges(5,5) * t252 + Icges(5,6) * t251 + Icges(5,3) * t439;
t195 = Icges(5,4) * t252 + Icges(5,2) * t251 + Icges(5,6) * t439;
t197 = Icges(5,1) * t252 + Icges(5,4) * t251 + Icges(5,5) * t439;
t327 = -t195 * t296 + t197 * t298;
t34 = (qJD(3) * t327 - t112) * t284 + (qJD(3) * t193 - t114 * t296 + t116 * t298 + (-t195 * t298 - t197 * t296) * qJD(4)) * t283;
t58 = t175 * t439 + t176 * t251 + t177 * t252 + t178 * t224 + t179 * t225 + t223 * t308;
t399 = t34 / 0.2e1 + t58 / 0.2e1;
t398 = t291 * t288 * t218;
t123 = t223 * t440 + t224 * t249 + t225 * t250;
t93 = -t192 * t284 + t283 * t328;
t397 = t93 / 0.2e1 + t123 / 0.2e1;
t124 = t223 * t439 + t224 * t251 + t225 * t252;
t94 = -t193 * t284 + t283 * t327;
t396 = -t94 / 0.2e1 - t124 / 0.2e1;
t160 = (-Icges(6,1) * t287 - t460) * t441 + (Icges(6,5) * t283 + t284 * t340) * qJD(3);
t217 = -Icges(6,3) * t284 + t283 * t333;
t395 = t283 * t288 * t160 + t217 * t412 + t411 * t452;
t393 = -t257 + t426;
t392 = t283 * t298 * t177 + t284 * t225 * t409 + t223 * t412;
t391 = t179 * rSges(5,1) + t178 * rSges(5,2) + rSges(5,3) * t378;
t268 = pkin(3) * t382;
t390 = t297 * (pkin(7) * t309 + qJD(1) * t275 - t268) + t299 * (-pkin(7) * t386 + t269 + (-t284 * t414 - t381) * pkin(3)) + t246 * t413;
t389 = -t263 - t420;
t387 = t226 * t414;
t377 = t440 / 0.2e1;
t376 = t439 / 0.2e1;
t373 = t299 * t420;
t201 = t419 * t299;
t371 = -t281 * t284 - t280;
t368 = t284 * t402;
t367 = t299 * t401;
t366 = t284 * t104 + t161 * t440 + t309 * t220;
t143 = t389 * t299;
t128 = -t217 * t284 + (t452 - t453) * t283;
t165 = Icges(6,5) * t241 + Icges(6,6) * t240 + Icges(6,3) * t440;
t167 = Icges(6,4) * t241 + Icges(6,2) * t240 + Icges(6,6) * t440;
t169 = Icges(6,1) * t241 + Icges(6,4) * t240 + Icges(6,5) * t440;
t332 = -t167 * t287 + t169 * t288;
t85 = -t165 * t284 + t283 * t332;
t166 = Icges(6,5) * t243 + Icges(6,6) * t242 + Icges(6,3) * t439;
t168 = Icges(6,4) * t243 + Icges(6,2) * t242 + Icges(6,6) * t439;
t170 = Icges(6,1) * t243 + Icges(6,4) * t242 + Icges(6,5) * t439;
t331 = -t168 * t287 + t170 * t288;
t86 = -t166 * t284 + t283 * t331;
t347 = t85 * t297 + t86 * t299;
t110 = t217 * t440 + t218 * t240 + t219 * t241;
t75 = t165 * t440 + t167 * t240 + t169 * t241;
t76 = t166 * t440 + t168 * t240 + t170 * t241;
t350 = t297 * t75 + t299 * t76;
t40 = -t110 * t284 + t283 * t350;
t111 = t217 * t439 + t218 * t242 + t219 * t243;
t77 = t165 * t439 + t167 * t242 + t169 * t243;
t78 = t166 * t439 + t168 * t242 + t170 * t243;
t349 = t297 * t77 + t299 * t78;
t41 = -t111 * t284 + t283 * t349;
t100 = Icges(6,4) * t155 + Icges(6,2) * t154 + Icges(6,6) * t309;
t102 = Icges(6,1) * t155 + Icges(6,4) * t154 + Icges(6,5) * t309;
t98 = Icges(6,5) * t155 + Icges(6,6) * t154 + Icges(6,3) * t309;
t19 = t100 * t242 + t102 * t243 + t152 * t167 + t153 * t169 + t165 * t308 + t439 * t98;
t101 = Icges(6,1) * t153 + Icges(6,4) * t152 + Icges(6,5) * t308;
t97 = Icges(6,5) * t153 + Icges(6,6) * t152 + Icges(6,3) * t308;
t99 = Icges(6,4) * t153 + Icges(6,2) * t152 + Icges(6,6) * t308;
t20 = t101 * t243 + t152 * t168 + t153 * t170 + t166 * t308 + t242 * t99 + t439 * t97;
t159 = (-Icges(6,2) * t288 - t461) * t441 + (Icges(6,6) * t283 + t284 * t336) * qJD(3);
t48 = t152 * t218 + t153 * t219 + t158 * t439 + t159 * t242 + t160 * t243 + t217 * t308;
t57 = t297 * t78 - t299 * t77;
t5 = (qJD(3) * t349 - t48) * t284 + (-qJD(1) * t57 + qJD(3) * t111 + t19 * t297 + t20 * t299) * t283;
t21 = t100 * t240 + t102 * t241 + t154 * t167 + t155 * t169 + t165 * t309 + t440 * t98;
t22 = t101 * t241 + t154 * t168 + t155 * t170 + t166 * t309 + t240 * t99 + t440 * t97;
t49 = t154 * t218 + t155 * t219 + t158 * t440 + t159 * t240 + t160 * t241 + t217 * t309;
t56 = t297 * t76 - t299 * t75;
t6 = (qJD(3) * t350 - t49) * t284 + (-qJD(1) * t56 + qJD(3) * t110 + t21 * t297 + t22 * t299) * t283;
t358 = t41 * t378 + t5 * t439 + t6 * t440 + (-t128 * t284 + t283 * t347) * t412 + t309 * t40;
t356 = rSges(5,1) * t181 + rSges(5,2) * t180;
t348 = t297 * t86 - t299 * t85;
t89 = t192 * t440 + t194 * t249 + t196 * t250;
t90 = t193 * t440 + t195 * t249 + t197 * t250;
t62 = t297 * t90 - t299 * t89;
t346 = t297 * t89 + t299 * t90;
t91 = t192 * t439 + t194 * t251 + t196 * t252;
t92 = t193 * t439 + t195 * t251 + t197 * t252;
t63 = t297 * t92 - t299 * t91;
t345 = t297 * t91 + t299 * t92;
t344 = t297 * t93 + t299 * t94;
t338 = Icges(4,2) * t284 + t465;
t326 = -t199 * t297 + t456;
t325 = -t198 * t297 - t199 * t299;
t315 = qJD(3) * t261;
t314 = t183 * t299 - t297 * t423;
t312 = qJD(3) * t338;
t311 = qJD(3) * (-Icges(4,5) * t283 - Icges(4,6) * t284);
t310 = t283 * t476 - t280 - t474;
t307 = t283 * t467 + t371;
t12 = qJD(1) * t349 - t19 * t299 + t20 * t297;
t13 = qJD(1) * t350 - t21 * t299 + t22 * t297;
t25 = (qJD(3) * t332 - t98) * t284 + (qJD(3) * t165 + (-t167 * t291 + t102) * t288 + (-t169 * t291 - t100) * t287) * t283;
t26 = (qJD(3) * t331 - t97) * t284 + (qJD(3) * t166 + (-t168 * t291 + t101) * t288 + (-t170 * t291 - t99) * t287) * t283;
t305 = -t299 * t6 / 0.2e1 + t5 * t478 + t13 * t377 + t12 * t376 + (qJD(1) * t347 - t25 * t299 + t26 * t297) * t479 + t40 * t414 / 0.2e1 + t41 * t375 + t348 * t412 / 0.2e1 + t497 * t57 + t496 * t56;
t304 = rSges(4,2) * t386 + rSges(4,3) * t413 - t299 * t315;
t303 = t297 * t310 - t428;
t221 = -t297 * t318 + t299 * t468;
t125 = t128 * t412;
t60 = t492 * t284 + (-t398 + (-t219 * t291 - t159) * t287) * t283 + t395;
t7 = t125 + (qJD(3) * t347 - t60) * t284 + (-qJD(1) * t348 + t25 * t297 + t26 * t299) * t283;
t302 = -t284 * t7 - t386 * t41 + t358;
t301 = t125 + (t25 + t49) * t377 + (t26 + t48) * t376 + (t111 + t86) * t497 + (t110 + t85) * t496;
t285 = qJD(2) * t297;
t256 = t357 * qJD(3);
t248 = t263 * t414;
t233 = -rSges(4,3) * t299 + t297 * t357;
t206 = -qJD(1) * t222 + t286;
t205 = t221 * qJD(1) + t285;
t200 = t419 * t297;
t186 = t297 * t311 + t416;
t185 = -qJD(1) * t227 + t299 * t311;
t157 = t171 * t439;
t148 = t261 * t410 + (t299 * t316 - t289) * qJD(1) + t417;
t147 = t285 + (-t428 + (-t280 - t470) * t297) * qJD(1) + t304;
t145 = t372 - t421;
t144 = t303 + t355;
t142 = t389 * t297;
t141 = -t199 * t284 - t226 * t439;
t140 = t198 * t284 + t226 * t440;
t139 = t228 * t297 - t323 * t299;
t138 = t227 * t297 - t493;
t137 = -t228 * t299 - t494;
t136 = -t227 * t299 - t297 * t324;
t135 = -t172 * t284 - t220 * t439;
t133 = -t223 * t284 + (t225 * t298 - t451) * t283;
t132 = t317 + t372 + t172;
t131 = (pkin(4) * t296 - t295) * t299 + t307 * t297 + t352;
t130 = t133 * t412;
t129 = t326 * t283;
t127 = qJD(1) * t201 + t297 * t422;
t126 = t299 * t422 + t248 + t387;
t122 = -t367 + t268 + (-t368 + (-t442 - t495) * qJD(3)) * t297 + (t299 * t306 + t279) * qJD(1);
t120 = -t172 * t440 + t157;
t119 = rSges(5,3) * t309 + t356;
t118 = -rSges(5,3) * t386 + t391;
t109 = -t325 + t418;
t88 = t310 * t413 + t379 * t476 + t268 - t356 + t417;
t87 = -pkin(3) * t381 + qJD(1) * t303 + t269 + t285 + t391;
t84 = -t283 * t373 - t284 * t423;
t83 = t183 * t284 + t216 * t440 + t134;
t74 = qJD(1) * t143 + t297 * t393;
t73 = t299 * t393 + t414 * t420 + t248;
t72 = t283 * t314 + t157;
t71 = t367 + (t368 + (t284 * t467 + t442) * qJD(3)) * t297 + (t299 * t307 - t279) * qJD(1) - t353 + t417;
t70 = t285 + (-t368 + (-t284 * t300 - t442) * qJD(3)) * t299 + (-t428 + (t371 - t469) * t297) * qJD(1) + t388 + t394;
t69 = t418 - t485;
t68 = (t226 * t410 + t119) * t284 + (-qJD(3) * t198 + t182 * t297 + t226 * t413) * t283;
t67 = (-t226 * t408 - t118) * t284 + (qJD(3) * t199 - t182 * t299 + t387) * t283;
t66 = t491 * t284 + (-t433 + (-t224 * t298 - t225 * t296) * qJD(4)) * t283 + t392;
t65 = -t171 * t412 + t366;
t64 = -t161 * t439 + (-t220 * t408 - t103) * t284 + t425;
t51 = t326 * t411 + (qJD(1) * t325 - t118 * t297 + t119 * t299) * t283;
t50 = t118 * t299 + t119 * t297 + (t297 * t421 + t456) * qJD(1) + t390;
t45 = -t124 * t284 + t283 * t345;
t44 = -t123 * t284 + t283 * t346;
t42 = -t172 * t379 + (-t103 * t297 + (-t171 * t297 - t172 * t299) * qJD(1)) * t283 + t466;
t32 = (t216 * t410 + t122) * t284 + (-qJD(3) * t424 + t191 * t297 + t216 * t413) * t283 + t366;
t31 = (-qJD(3) * t373 - t427) * t284 + (qJD(3) * t184 + t216 * t414 + t299 * t426) * t283 + t425;
t30 = t112 * t440 + t114 * t249 + t116 * t250 + t180 * t195 + t181 * t197 + t193 * t309;
t29 = t113 * t440 + t115 * t249 + t117 * t250 + t180 * t194 + t181 * t196 + t192 * t309;
t28 = t112 * t439 + t114 * t251 + t116 * t252 + t178 * t195 + t179 * t197 + t193 * t308;
t27 = t113 * t439 + t115 * t251 + t117 * t252 + t178 * t194 + t179 * t196 + t192 * t308;
t18 = t427 * t299 + (t104 + t122) * t297 + (t424 * t299 + (-t247 - t423) * t297) * qJD(1) + t390;
t17 = t314 * t411 + (qJD(1) * t485 + t122 * t299 - t427 * t297) * t283 + t466;
t16 = qJD(1) * t346 - t29 * t299 + t297 * t30;
t15 = qJD(1) * t345 - t27 * t299 + t28 * t297;
t9 = (qJD(3) * t346 - t59) * t284 + (-t62 * qJD(1) + qJD(3) * t123 + t29 * t297 + t299 * t30) * t283;
t8 = (qJD(3) * t345 - t58) * t284 + (-t63 * qJD(1) + qJD(3) * t124 + t27 * t297 + t28 * t299) * t283;
t1 = [-t225 * t380 - t287 * t219 * t441 + 0.2e1 * m(3) * (t205 * t222 + t206 * t221) + (t147 * t211 + t148 * t210) * t484 + (t131 * t71 + t132 * t70) * t482 + (t144 * t88 + t145 * t87) * t483 + t395 + t392 + (t343 - t338) * t412 + (Icges(4,1) * t283 + t339 + t464) * t411 + (t491 + t492) * t284 + (-t287 * t159 - t224 * t405 - t398 - t433) * t283; m(6) * (t297 * t71 - t299 * t70 + (t131 * t299 + t132 * t297) * qJD(1)) + m(5) * (t297 * t88 - t299 * t87 + (t144 * t299 + t145 * t297) * qJD(1)) + m(4) * (t489 * qJD(1) - t147 * t299 + t148 * t297) + m(3) * (-t205 * t299 + t206 * t297 + (t221 * t299 + t222 * t297) * qJD(1)); 0; (-t25 / 0.2e1 - t49 / 0.2e1 + (qJD(1) * t230 - t297 * t312) * t479 + t232 * t498 + (t450 / 0.2e1 - t445 / 0.2e1) * qJD(3) - t400) * t299 + (t26 / 0.2e1 + t48 / 0.2e1 + (-qJD(1) * t229 - t299 * t312) * t284 / 0.2e1 + t231 * t498 + (-t448 / 0.2e1 + t443 / 0.2e1) * qJD(3) + t399) * t297 + m(4) * ((-t147 * t297 - t148 * t299) * t261 - t489 * t256) + m(6) * (t131 * t73 + t132 * t74 + t142 * t70 + t143 * t71) + m(5) * (t126 * t144 + t127 * t145 + t200 * t87 + t201 * t88) + (t481 / 0.2e1 + t480 / 0.2e1) * t335 * qJD(3) + ((-t211 * t475 + t447 / 0.2e1 + t444 / 0.2e1 + t111 / 0.2e1 + t86 / 0.2e1 - t396) * t299 + (t210 * t475 + t110 / 0.2e1 + t85 / 0.2e1 + t449 / 0.2e1 + t446 / 0.2e1 + t397) * t297) * qJD(1); m(5) * (t126 * t297 - t127 * t299 + (t200 * t297 + t201 * t299) * qJD(1)) + m(6) * (t297 * t73 - t299 * t74 + (t142 * t297 + t143 * t299) * qJD(1)); (t142 * t74 + t143 * t73 + t18 * t69) * t482 - t299 * t13 + t297 * t12 + t297 * t15 + (t109 * t50 + t126 * t201 + t127 * t200) * t483 - t299 * t16 - t299 * ((t299 * t186 + (t137 + t493) * qJD(1)) * t299 + (t136 * qJD(1) + (-t230 * t411 - t232 * t412 + t416) * t297 + (-t185 + (t446 + t449) * qJD(3) - t323 * qJD(1)) * t299) * t297) + t297 * ((t297 * t185 + (t138 + t494) * qJD(1)) * t297 + (t139 * qJD(1) + (t229 * t411 + t231 * t412) * t299 + (-t186 + (-t444 - t447) * qJD(3) + (t228 - t324) * qJD(1)) * t297) * t299) + ((t233 * t297 + t234 * t299) * ((qJD(1) * t233 + t304) * t299 + (-t297 * t315 + (-t234 + t490) * qJD(1)) * t297) + (t480 + t481) * t261 * t256) * t484 + (-t136 * t299 + t137 * t297 + t56 + t62) * t414 + (-t138 * t299 + t139 * t297 + t57 + t63) * t413; (-t60 - t66 + (t297 * t397 - t299 * t396) * qJD(3)) * t284 + m(6) * (t131 * t32 + t132 * t31 + t70 * t84 + t71 * t83) + m(5) * (t140 * t88 + t141 * t87 + t144 * t68 + t145 * t67) + t130 + t301 + (t399 * t299 + t400 * t297 + (t297 * t396 + t299 * t397) * qJD(1)) * t283; m(5) * (t297 * t68 - t299 * t67 + (t140 * t299 + t141 * t297) * qJD(1)) + m(6) * (t297 * t32 - t299 * t31 + (t297 * t84 + t299 * t83) * qJD(1)); (qJD(3) * (t297 * t94 - t299 * t93) / 0.2e1 + t15 * t477 + t16 * t478 + (t62 * t477 - t297 * t63 / 0.2e1) * qJD(1)) * t283 + (t63 * t374 + t45 * t473 - t9 / 0.2e1 + (qJD(1) * t94 - t33) * t479) * t299 + (t62 * t374 + t44 * t473 + t8 / 0.2e1 + (qJD(1) * t93 + t34) * t479) * t297 + m(6) * (t142 * t31 + t143 * t32 + t17 * t69 + t18 * t72 + t73 * t83 + t74 * t84) + m(5) * (t109 * t51 + t126 * t140 + t127 * t141 + t129 * t50 + t200 * t67 + t201 * t68) + t305; (t17 * t72 + t31 * t84 + t32 * t83) * t482 + (t129 * t51 + t140 * t68 + t141 * t67) * t483 + (t66 * t284 - t130 - t7 + (-t284 * t344 + t297 * t44 + t299 * t45) * qJD(3)) * t284 + (t299 * t8 + t297 * t9 - t284 * (t297 * t33 + t299 * t34) + (-t133 * t284 + t283 * t344) * qJD(3) + ((-t284 * t93 + t44) * t299 + (t94 * t284 - t41 - t45) * t297) * qJD(1)) * t283 + t358; m(6) * (t131 * t65 + t132 * t64 + t134 * t71 + t135 * t70) + t301 - t60 * t284; m(6) * (t297 * t65 - t299 * t64 + (t134 * t299 + t135 * t297) * qJD(1)); m(6) * (t120 * t18 + t134 * t73 + t135 * t74 + t142 * t64 + t143 * t65 + t42 * t69) + t305; m(6) * (t120 * t17 + t134 * t32 + t135 * t31 + t42 * t72 + t64 * t84 + t65 * t83) + t302; (t120 * t42 + t134 * t65 + t135 * t64) * t482 + t302;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
