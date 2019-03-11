% Calculate time derivative of joint inertia matrix for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR2_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR2_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR2_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:41:16
% EndTime: 2019-03-09 02:41:37
% DurationCPUTime: 13.57s
% Computational Cost: add. (18169->645), mult. (17550->889), div. (0->0), fcn. (15881->10), ass. (0->337)
t245 = qJ(3) + pkin(10);
t240 = sin(t245);
t242 = cos(t245);
t302 = Icges(6,4) * t242 - Icges(6,5) * t240;
t249 = sin(qJ(3));
t252 = cos(qJ(3));
t494 = Icges(4,5) * t252 + Icges(5,5) * t242 - Icges(4,6) * t249 - Icges(5,6) * t240;
t498 = -t302 + t494;
t414 = Icges(6,6) * t240;
t424 = Icges(5,4) * t240;
t497 = -t414 - t424 + (-Icges(5,2) - Icges(6,3)) * t242;
t413 = Icges(6,6) * t242;
t423 = Icges(5,4) * t242;
t496 = -t413 - t423 + (-Icges(5,1) - Icges(6,2)) * t240;
t495 = Icges(6,1) + Icges(4,3) + Icges(5,3);
t246 = qJ(1) + pkin(9);
t241 = sin(t246);
t243 = cos(t246);
t304 = -Icges(5,2) * t240 + t423;
t139 = -Icges(5,6) * t243 + t241 * t304;
t295 = -Icges(6,3) * t240 + t413;
t459 = Icges(6,5) * t243 + t241 * t295;
t474 = t139 + t459;
t140 = Icges(5,6) * t241 + t243 * t304;
t143 = Icges(6,5) * t241 - t243 * t295;
t473 = t140 - t143;
t309 = Icges(5,1) * t242 - t424;
t141 = -Icges(5,5) * t243 + t241 * t309;
t297 = Icges(6,2) * t242 - t414;
t458 = Icges(6,4) * t243 + t241 * t297;
t472 = t141 + t458;
t142 = Icges(5,5) * t241 + t243 * t309;
t145 = Icges(6,4) * t241 - t243 * t297;
t471 = t142 - t145;
t493 = t497 * qJD(3);
t492 = t496 * qJD(3);
t488 = t495 * t241 + t243 * t498;
t426 = Icges(4,4) * t249;
t311 = Icges(4,1) * t252 - t426;
t166 = Icges(4,5) * t241 + t243 * t311;
t400 = t166 * t252;
t425 = Icges(4,4) * t252;
t306 = -Icges(4,2) * t249 + t425;
t163 = Icges(4,6) * t241 + t243 * t306;
t405 = t163 * t249;
t461 = t240 * t473 - t242 * t471 - t400 + t405;
t491 = t241 * t461;
t490 = t461 * t243;
t489 = t241 * t498 - t495 * t243;
t487 = (-Icges(4,5) * t249 - Icges(4,6) * t252 + (Icges(6,5) - Icges(5,6)) * t242 + (Icges(6,4) - Icges(5,5)) * t240) * qJD(3);
t165 = -Icges(4,5) * t243 + t241 * t311;
t402 = t165 * t252;
t162 = -Icges(4,6) * t243 + t241 * t306;
t407 = t162 * t249;
t486 = t240 * t474 - t472 * t242 - t402 + t407;
t373 = qJD(3) * t243;
t348 = t240 * t373;
t379 = qJD(1) * t241;
t485 = t242 * t379 + t348;
t484 = t489 * t241 - t491 + (-t486 - t488) * t243;
t483 = t488 * qJD(1);
t482 = t241 * t471 + t243 * t472;
t481 = t241 * t473 + t243 * t474;
t480 = t241 ^ 2;
t479 = t243 ^ 2;
t478 = qJD(3) / 0.2e1;
t442 = pkin(3) * t252;
t237 = pkin(2) + t442;
t216 = t243 * t237;
t247 = -qJ(4) - pkin(7);
t436 = -pkin(7) - t247;
t135 = -pkin(2) * t243 + t241 * t436 + t216;
t432 = rSges(5,1) * t242;
t331 = -rSges(5,2) * t240 + t432;
t150 = -rSges(5,3) * t243 + t241 * t331;
t395 = t242 * t243;
t232 = t241 * rSges(5,3);
t399 = t240 * t243;
t470 = -rSges(5,2) * t399 + t232;
t151 = rSges(5,1) * t395 + t470;
t198 = rSges(5,1) * t240 + rSges(5,2) * t242;
t275 = qJD(3) * t198;
t236 = t243 * pkin(7);
t394 = t243 * t247;
t437 = pkin(2) - t237;
t475 = t241 * t437;
t134 = t236 + t394 - t475;
t230 = qJD(4) * t241;
t371 = qJD(3) * t249;
t361 = pkin(3) * t371;
t354 = qJD(4) * t243 + t241 * t361 + t247 * t379;
t378 = qJD(1) * t243;
t439 = pkin(7) * t241;
t365 = t134 * t378 + t241 * ((-t243 * t437 - t439) * qJD(1) - t354) + t243 * (-t243 * t361 + t230 + (t243 * t436 + t475) * qJD(1));
t353 = t240 * t379;
t385 = rSges(5,2) * t353 + rSges(5,3) * t378;
t20 = (qJD(1) * t150 - t243 * t275 + t385) * t243 + (-t241 * t275 + (-t135 - t151 + t470) * qJD(1)) * t241 + t365;
t455 = 2 * m(5);
t477 = t20 * t455;
t431 = rSges(4,2) * t249;
t433 = rSges(4,1) * t252;
t332 = -t431 + t433;
t429 = rSges(4,3) * t243;
t167 = t241 * t332 - t429;
t364 = t243 * t431;
t233 = t241 * rSges(4,3);
t384 = t243 * t433 + t233;
t169 = -t364 + t384;
t224 = rSges(4,1) * t249 + rSges(4,2) * t252;
t276 = qJD(3) * t224;
t377 = qJD(1) * t249;
t351 = t241 * t377;
t255 = rSges(4,2) * t351 + rSges(4,3) * t378 - t243 * t276;
t36 = (qJD(1) * t167 + t255) * t243 + (-t241 * t276 + (-t169 - t364 + t233) * qJD(1)) * t241;
t456 = 2 * m(4);
t476 = t36 * t456;
t234 = t241 * rSges(6,1);
t469 = -rSges(6,2) * t395 + t234;
t444 = sin(qJ(1)) * pkin(1);
t468 = t236 - t444;
t467 = t241 * t486 + t243 * t489;
t464 = t241 * t488 - t490;
t463 = -t241 * t487 - t483;
t462 = -qJD(1) * t489 + t243 * t487;
t454 = 2 * m(6);
t453 = 2 * m(7);
t452 = m(6) / 0.2e1;
t451 = m(7) / 0.2e1;
t450 = t240 / 0.2e1;
t449 = t241 / 0.2e1;
t448 = -t243 / 0.2e1;
t447 = t243 / 0.2e1;
t446 = rSges(7,3) + pkin(8);
t445 = m(4) * t224;
t443 = pkin(3) * t249;
t441 = pkin(4) * t240;
t440 = pkin(4) * t242;
t235 = t241 * pkin(5);
t244 = cos(qJ(1)) * pkin(1);
t438 = qJD(1) / 0.2e1;
t248 = sin(qJ(6));
t251 = cos(qJ(6));
t421 = Icges(7,4) * t248;
t301 = Icges(7,2) * t251 + t421;
t367 = qJD(6) * t242;
t420 = Icges(7,4) * t251;
t116 = (Icges(7,2) * t248 - t420) * t367 + (Icges(7,6) * t242 + t240 * t301) * qJD(3);
t307 = Icges(7,1) * t248 + t420;
t164 = Icges(7,5) * t240 - t242 * t307;
t298 = Icges(7,5) * t248 + Icges(7,6) * t251;
t113 = (-Icges(7,5) * t251 + Icges(7,6) * t248) * t367 + (Icges(7,3) * t242 + t240 * t298) * qJD(3);
t158 = Icges(7,3) * t240 - t242 * t298;
t161 = Icges(7,6) * t240 - t242 * t301;
t370 = qJD(3) * t251;
t372 = qJD(3) * t248;
t374 = qJD(3) * t242;
t312 = t248 * t161 * t367 + t158 * t374 + (t161 * t370 + t164 * t372 + t113) * t240;
t119 = (-Icges(7,1) * t251 + t421) * t367 + (Icges(7,5) * t242 + t240 * t307) * qJD(3);
t408 = t119 * t248;
t69 = t158 * t240 + (-t161 * t251 - t164 * t248) * t242;
t435 = ((-t408 + (-qJD(6) * t164 - t116) * t251) * t242 + t312) * t240 + t69 * t374;
t342 = qJD(6) * t240 + qJD(1);
t256 = t242 * t370 - t248 * t342;
t341 = qJD(1) * t240 + qJD(6);
t282 = t241 * t341;
t86 = t243 * t256 - t251 * t282;
t257 = t242 * t372 + t251 * t342;
t87 = t243 * t257 - t248 * t282;
t434 = t87 * rSges(7,1) + t86 * rSges(7,2);
t430 = rSges(6,2) * t240;
t428 = -rSges(6,3) - qJ(5);
t410 = qJ(5) * t240;
t409 = qJ(5) * t242;
t398 = t241 * t242;
t397 = t241 * t248;
t396 = t241 * t251;
t393 = t243 * t248;
t392 = t243 * t251;
t175 = t240 * t392 - t397;
t176 = t240 * t393 + t396;
t108 = t176 * rSges(7,1) + t175 * rSges(7,2) + rSges(7,3) * t395;
t391 = pkin(8) * t395 + t108 + t235;
t177 = t240 * t396 + t393;
t178 = t240 * t397 - t392;
t330 = -rSges(7,1) * t178 - rSges(7,2) * t177;
t109 = rSges(7,3) * t398 - t330;
t390 = -pkin(5) * t243 + pkin(8) * t398 + t109;
t389 = t241 * t134 + t243 * t135;
t172 = pkin(4) * t395 + qJ(5) * t399;
t388 = -t135 - t172;
t196 = -t409 + t441;
t219 = pkin(3) * t351;
t387 = t196 * t379 + t219;
t347 = t242 * t373;
t368 = qJD(5) * t240;
t386 = qJ(5) * t347 + t243 * t368;
t383 = t479 + t480;
t376 = qJD(3) * t240;
t375 = qJD(3) * t241;
t369 = qJD(3) * t252;
t366 = -pkin(4) - t446;
t360 = pkin(3) * t369;
t102 = Icges(7,5) * t176 + Icges(7,6) * t175 + Icges(7,3) * t395;
t104 = Icges(7,4) * t176 + Icges(7,2) * t175 + Icges(7,6) * t395;
t106 = Icges(7,1) * t176 + Icges(7,4) * t175 + Icges(7,5) * t395;
t293 = t104 * t251 + t106 * t248;
t40 = Icges(7,5) * t87 + Icges(7,6) * t86 - Icges(7,3) * t485;
t42 = Icges(7,4) * t87 + Icges(7,2) * t86 - Icges(7,6) * t485;
t44 = Icges(7,1) * t87 + Icges(7,4) * t86 - Icges(7,5) * t485;
t10 = (qJD(3) * t293 + t40) * t240 + (qJD(3) * t102 - t248 * t44 - t251 * t42 + (t104 * t248 - t106 * t251) * qJD(6)) * t242;
t19 = t113 * t395 + t116 * t175 + t119 * t176 - t158 * t485 + t161 * t86 + t164 * t87;
t359 = t10 / 0.2e1 + t19 / 0.2e1;
t103 = Icges(7,5) * t178 + Icges(7,6) * t177 + Icges(7,3) * t398;
t105 = Icges(7,4) * t178 + Icges(7,2) * t177 + Icges(7,6) * t398;
t107 = Icges(7,1) * t178 + Icges(7,4) * t177 + Icges(7,5) * t398;
t292 = t105 * t251 + t107 * t248;
t349 = t240 * t375;
t350 = t242 * t378;
t262 = -t349 + t350;
t281 = t243 * t341;
t84 = t241 * t256 + t251 * t281;
t85 = t241 * t257 + t248 * t281;
t39 = Icges(7,5) * t85 + Icges(7,6) * t84 + Icges(7,3) * t262;
t41 = Icges(7,4) * t85 + Icges(7,2) * t84 + Icges(7,6) * t262;
t43 = Icges(7,1) * t85 + Icges(7,4) * t84 + Icges(7,5) * t262;
t11 = (qJD(3) * t292 + t39) * t240 + (qJD(3) * t103 - t248 * t43 - t251 * t41 + (t105 * t248 - t107 * t251) * qJD(6)) * t242;
t18 = t113 * t398 + t116 * t177 + t119 * t178 + t158 * t262 + t161 * t84 + t164 * t85;
t358 = t11 / 0.2e1 + t18 / 0.2e1;
t32 = t102 * t240 - t242 * t293;
t48 = t158 * t395 + t161 * t175 + t164 * t176;
t357 = -t32 / 0.2e1 - t48 / 0.2e1;
t33 = t103 * t240 - t242 * t292;
t49 = t158 * t398 + t161 * t177 + t164 * t178;
t356 = t33 / 0.2e1 + t49 / 0.2e1;
t355 = t230 + t386;
t346 = -t376 / 0.2e1;
t345 = -t196 - t443;
t344 = -t198 - t443;
t329 = rSges(7,1) * t248 + rSges(7,2) * t251;
t168 = rSges(7,3) * t240 - t242 * t329;
t343 = pkin(8) * t240 + t168;
t325 = t410 + t440;
t171 = t325 * t241;
t340 = t241 * t171 + t243 * t172 + t389;
t339 = rSges(6,1) * t378 + rSges(6,2) * t485 + rSges(6,3) * t347;
t204 = pkin(4) * t349;
t338 = t204 + t354;
t327 = rSges(6,3) * t242 + t430;
t337 = t327 + t345;
t336 = t85 * rSges(7,1) + t84 * rSges(7,2);
t335 = -t241 * t247 + t216 + t244;
t334 = -t302 * qJD(3) / 0.2e1 + t494 * t478;
t155 = t344 * t243;
t328 = -rSges(6,2) * t242 + rSges(6,3) * t240;
t326 = -t394 - t444;
t27 = t102 * t395 + t104 * t175 + t106 * t176;
t28 = t103 * t395 + t105 * t175 + t107 * t176;
t320 = t241 * t28 + t243 * t27;
t16 = t241 * t27 - t243 * t28;
t29 = t102 * t398 + t104 * t177 + t106 * t178;
t30 = t103 * t398 + t105 * t177 + t107 * t178;
t319 = t241 * t30 + t243 * t29;
t17 = t241 * t29 - t243 * t30;
t318 = t241 * t33 + t243 * t32;
t317 = t241 * t32 - t243 * t33;
t260 = t242 * t366 - t237 - t410;
t254 = t241 * t260 - t444;
t50 = (pkin(5) - t247) * t243 + t254 + t330;
t273 = t335 + t172;
t51 = t273 + t391;
t316 = t241 * t51 + t243 * t50;
t65 = -t109 * t240 + t168 * t398;
t66 = t108 * t240 - t168 * t395;
t315 = t241 * t66 + t243 * t65;
t313 = t240 * t428 - t237;
t259 = (rSges(6,2) - pkin(4)) * t242 + t313;
t72 = -t444 + (rSges(6,1) - t247) * t243 + t259 * t241;
t152 = rSges(6,3) * t399 + t469;
t73 = t152 + t273;
t314 = t241 * t73 + t243 * t72;
t291 = t108 * t241 - t109 * t243;
t280 = t171 * t378 + t241 * (pkin(4) * t350 + t241 * t368 - t204 + (t240 * t378 + t241 * t374) * qJ(5)) + t243 * (-pkin(4) * t485 - qJ(5) * t353 + t386) + t365;
t170 = qJD(3) * t325 - qJD(5) * t242;
t279 = -t328 * qJD(3) - t170 - t360;
t278 = -pkin(2) - t332;
t124 = t337 * t243;
t277 = -t237 - t331;
t274 = -t343 + t345;
t77 = t274 * t243;
t122 = (-rSges(7,1) * t251 + rSges(7,2) * t248) * t367 + (rSges(7,3) * t242 + t240 * t329) * qJD(3);
t258 = -t122 - t170 + (-pkin(8) * t242 - t442) * qJD(3);
t229 = pkin(5) * t378;
t208 = t332 * qJD(3);
t189 = t331 * qJD(3);
t154 = t344 * t241;
t153 = -rSges(6,1) * t243 + t241 * t328;
t127 = t439 + t244 + (pkin(2) - t431) * t243 + t384;
t126 = t241 * t278 + t429 + t468;
t123 = t337 * t241;
t112 = t151 + t335;
t111 = -t444 + (rSges(5,3) - t247) * t243 + t277 * t241;
t83 = -t198 * t378 - t189 * t241 + (-t241 * t369 - t243 * t377) * pkin(3);
t82 = t198 * t379 + t219 + (-t189 - t360) * t243;
t76 = t274 * t241;
t75 = t224 * t375 + (-t244 + (-rSges(4,3) - pkin(7)) * t241 + t278 * t243) * qJD(1);
t74 = ((-pkin(2) - t433) * t241 + t468) * qJD(1) + t255;
t64 = t198 * t375 + (t243 * t277 - t232 - t244) * qJD(1) + t354;
t63 = t230 + qJD(3) * t155 + ((-t237 - t432) * t241 + t326) * qJD(1) + t385;
t61 = qJD(1) * t124 + t241 * t279;
t60 = t243 * t279 - t327 * t379 + t387;
t47 = t291 * t242;
t46 = -rSges(7,3) * t485 + t434;
t45 = rSges(7,3) * t262 + t336;
t38 = (-t368 + (t242 * t428 - t430) * qJD(3)) * t241 + (t243 * t259 - t234 - t244) * qJD(1) + t338;
t37 = (-t441 - t443) * t373 + ((t313 - t440) * t241 + t326) * qJD(1) + t339 + t355;
t35 = qJD(1) * t77 + t241 * t258;
t34 = t243 * t258 + t343 * t379 + t387;
t31 = t152 * t243 + t153 * t241 + t340;
t26 = (-t368 + (t240 * t446 - t409) * qJD(3)) * t241 + (t243 * t260 - t235 - t244) * qJD(1) - t336 + t338;
t25 = t229 + (t240 * t366 - t443) * t373 + (t254 - t394) * qJD(1) + t355 + t434;
t24 = t241 * t390 + t243 * t391 + t340;
t22 = (-t168 * t375 - t45) * t240 + (-qJD(3) * t109 + t122 * t241 + t168 * t378) * t242;
t21 = (t168 * t373 + t46) * t240 + (qJD(3) * t108 - t122 * t243 + t168 * t379) * t242;
t15 = (qJD(1) * t153 + t339) * t243 + (t327 * t375 + (-t152 + t388 + t469) * qJD(1)) * t241 + t280;
t14 = t291 * t376 + (-t241 * t46 + t243 * t45 + (-t108 * t243 - t109 * t241) * qJD(1)) * t242;
t13 = t240 * t49 + t242 * t319;
t12 = t48 * t240 + t242 * t320;
t9 = -t103 * t485 + t105 * t86 + t107 * t87 + t175 * t41 + t176 * t43 + t39 * t395;
t8 = -t102 * t485 + t104 * t86 + t106 * t87 + t175 * t42 + t176 * t44 + t395 * t40;
t7 = t103 * t262 + t105 * t84 + t107 * t85 + t177 * t41 + t178 * t43 + t39 * t398;
t6 = t102 * t262 + t104 * t84 + t106 * t85 + t177 * t42 + t178 * t44 + t398 * t40;
t5 = (-pkin(8) * t348 + qJD(1) * t390 + t229 + t46) * t243 + (-pkin(8) * t349 + t45 + (t388 - t391 + t235) * qJD(1)) * t241 + t280;
t4 = qJD(1) * t320 + t241 * t8 - t243 * t9;
t3 = qJD(1) * t319 + t241 * t6 - t243 * t7;
t2 = (-qJD(3) * t320 + t19) * t240 + (-qJD(1) * t16 + qJD(3) * t48 + t241 * t9 + t243 * t8) * t242;
t1 = (-qJD(3) * t319 + t18) * t240 + (-qJD(1) * t17 + qJD(3) * t49 + t241 * t7 + t243 * t6) * t242;
t23 = [(t25 * t51 + t26 * t50) * t453 + (t37 * t73 + t38 * t72) * t454 + (t111 * t64 + t112 * t63) * t455 + (t126 * t75 + t127 * t74) * t456 + t312 - t251 * t164 * t367 + (-Icges(4,2) * t252 + t311 - t426) * t371 + (Icges(4,1) * t249 + t306 + t425) * t369 + (-t116 * t251 - t408) * t242 + (t309 + t297 + t497) * t376 + (t295 + t304 - t496) * t374; 0; 0; m(7) * (t25 * t76 + t26 * t77 + t34 * t50 + t35 * t51) + m(6) * (t123 * t37 + t124 * t38 + t60 * t72 + t61 * t73) + m(5) * (t111 * t82 + t112 * t83 + t154 * t63 + t155 * t64) + (m(4) * (-t126 * t208 - t224 * t75) + t334 * t243 + (t407 / 0.2e1 - t402 / 0.2e1) * qJD(3) - t358) * t243 + (m(4) * (-t127 * t208 - t224 * t74) + t334 * t241 + (-t405 / 0.2e1 + t400 / 0.2e1) * qJD(3) + t359) * t241 + ((-qJD(3) * t473 + t243 * t492) * t449 + (-qJD(3) * t474 + t241 * t492) * t448) * t240 + ((qJD(3) * t471 + t243 * t493) * t449 + (qJD(3) * t472 + t241 * t493) * t448) * t242 + ((t471 * t448 - t472 * t449) * t240 + (t473 * t448 - t474 * t449) * t242 + (-t127 * t445 + (t140 / 0.2e1 - t143 / 0.2e1) * t242 + (t142 / 0.2e1 - t145 / 0.2e1) * t240 - t357) * t243 + (t126 * t445 + (t139 / 0.2e1 + t459 / 0.2e1) * t242 + (t141 / 0.2e1 + t458 / 0.2e1) * t240 + t356) * t241) * qJD(1); m(4) * t36 + m(5) * t20 + m(6) * t15 + m(7) * t5; (t24 * t5 + t34 * t77 + t35 * t76) * t453 + t16 * t378 + t17 * t379 + (t123 * t61 + t124 * t60 + t15 * t31) * t454 + (t154 * t83 + t155 * t82 + t20 * t389) * t455 + t383 * t224 * t208 * t456 + (t151 * t477 + t169 * t476 + t467 * t379 + t463 * t479 - t3 + (-t243 * t486 - t484) * t378) * t243 + (t4 + t150 * t477 + t167 * t476 + t462 * t480 + t464 * t378 + ((t163 * t369 + t166 * t371 + t463 - t483) * t241 + (t162 * t369 + t165 * t371 + t462) * t243 + t482 * t376 + t481 * t374 + ((-t163 * t252 - t166 * t249) * t241 + (-t162 * t252 - t165 * t249) * t243 - t481 * t242 - t482 * t240) * qJD(3) + ((-t486 + t488) * t241 + t490 + t464 + t467) * qJD(1)) * t243 + (t484 + t491) * t379) * t241; m(7) * (qJD(1) * t316 + t241 * t26 - t243 * t25) + m(5) * (t241 * t64 - t243 * t63 + (t111 * t243 + t112 * t241) * qJD(1)) + m(6) * (qJD(1) * t314 + t241 * t38 - t243 * t37); 0; m(7) * (t241 * t34 - t243 * t35 + (t241 * t76 + t243 * t77) * qJD(1)) + m(6) * (t241 * t60 - t243 * t61 + (t123 * t241 + t124 * t243) * qJD(1)) + m(5) * (t241 * t82 - t243 * t83 + (t154 * t241 + t155 * t243) * qJD(1)); 0; 0.2e1 * (t314 * t452 + t316 * t451) * t374 + 0.2e1 * ((t241 * t25 + t243 * t26 + t378 * t51 - t379 * t50) * t451 + (t241 * t37 + t243 * t38 + t378 * t73 - t379 * t72) * t452) * t240; (m(6) + m(7)) * t376; 0.2e1 * ((t373 * t77 + t375 * t76 - t5) * t451 + (t123 * t375 + t124 * t373 - t15) * t452) * t242 + 0.2e1 * ((qJD(3) * t24 + t241 * t35 + t243 * t34 + t378 * t76 - t379 * t77) * t451 + (qJD(3) * t31 + t123 * t378 - t124 * t379 + t241 * t61 + t243 * t60) * t452) * t240; 0; 0.4e1 * (t452 + t451) * (-0.1e1 + t383) * t240 * t374; m(7) * (t21 * t51 + t22 * t50 + t25 * t66 + t26 * t65) + (-t241 * t356 + t243 * t357) * t376 + (t359 * t243 + t358 * t241 + (t241 * t357 + t243 * t356) * qJD(1)) * t242 + t435; m(7) * t14; m(7) * (t14 * t24 + t21 * t76 + t22 * t77 + t34 * t65 + t35 * t66 - t47 * t5) + (t12 * t438 + t16 * t346 - t1 / 0.2e1 + (qJD(1) * t32 - t11) * t450) * t243 + (t2 / 0.2e1 + t17 * t346 + t13 * t438 + (qJD(1) * t33 + t10) * t450) * t241 + (t3 * t449 + t4 * t447 + t317 * t478 + (t17 * t447 - t241 * t16 / 0.2e1) * qJD(1)) * t242; m(7) * (qJD(1) * t315 - t21 * t243 + t22 * t241); m(7) * ((qJD(3) * t315 - t14) * t242 + (-qJD(3) * t47 + t21 * t241 + t22 * t243 + (-t241 * t65 + t243 * t66) * qJD(1)) * t240); (-t14 * t47 + t21 * t66 + t22 * t65) * t453 + ((-t243 * t12 - t241 * t13 - t240 * t318) * qJD(3) + t435) * t240 + (t243 * t2 + t241 * t1 + t240 * (t10 * t243 + t11 * t241) + (t69 * t240 + t242 * t318) * qJD(3) + (-t241 * t12 + t243 * t13 - t240 * t317) * qJD(1)) * t242;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t23(1) t23(2) t23(4) t23(7) t23(11) t23(16); t23(2) t23(3) t23(5) t23(8) t23(12) t23(17); t23(4) t23(5) t23(6) t23(9) t23(13) t23(18); t23(7) t23(8) t23(9) t23(10) t23(14) t23(19); t23(11) t23(12) t23(13) t23(14) t23(15) t23(20); t23(16) t23(17) t23(18) t23(19) t23(20) t23(21);];
Mq  = res;
