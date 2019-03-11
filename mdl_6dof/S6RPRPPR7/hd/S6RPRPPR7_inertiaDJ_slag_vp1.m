% Calculate time derivative of joint inertia matrix for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
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
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR7_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR7_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:56:14
% EndTime: 2019-03-09 02:56:38
% DurationCPUTime: 18.12s
% Computational Cost: add. (11812->754), mult. (17980->1058), div. (0->0), fcn. (16203->8), ass. (0->385)
t259 = qJ(3) + pkin(9);
t248 = cos(t259);
t247 = sin(t259);
t442 = Icges(6,6) * t247;
t452 = Icges(5,4) * t247;
t509 = -t442 - t452 + (Icges(5,1) + Icges(6,2)) * t248;
t441 = Icges(6,6) * t248;
t451 = Icges(5,4) * t248;
t508 = t441 + t451 + (-Icges(5,2) - Icges(6,3)) * t247;
t507 = t508 * qJD(3);
t506 = t509 * qJD(3);
t265 = sin(qJ(1));
t268 = cos(qJ(1));
t318 = Icges(5,2) * t248 + t452;
t141 = Icges(5,6) * t268 + t265 * t318;
t306 = Icges(6,3) * t248 + t442;
t146 = Icges(6,5) * t268 - t265 * t306;
t505 = -t141 + t146;
t308 = Icges(6,2) * t247 + t441;
t147 = Icges(6,4) * t265 + t268 * t308;
t323 = Icges(5,1) * t247 + t451;
t489 = -Icges(5,5) * t265 + t268 * t323;
t504 = t147 + t489;
t143 = Icges(5,5) * t268 + t265 * t323;
t148 = Icges(6,4) * t268 - t265 * t308;
t503 = -t148 + t143;
t145 = Icges(6,5) * t265 + t268 * t306;
t491 = -Icges(5,6) * t265 + t268 * t318;
t502 = t491 + t145;
t264 = sin(qJ(3));
t267 = cos(qJ(3));
t396 = qJD(3) * t267;
t402 = qJD(1) * t268;
t501 = t264 * t402 + t265 * t396;
t398 = qJD(3) * t265;
t366 = t248 * t398;
t274 = t247 * t402 + t366;
t500 = -qJD(3) / 0.2e1;
t263 = sin(qJ(6));
t266 = cos(qJ(6));
t357 = -qJD(1) * t248 - qJD(6);
t293 = t357 * t268;
t358 = qJD(6) * t248 + qJD(1);
t397 = qJD(3) * t266;
t87 = t266 * t293 + (t247 * t397 + t263 * t358) * t265;
t368 = t247 * t398;
t421 = t265 * t266;
t88 = -t358 * t421 + (t293 + t368) * t263;
t43 = t88 * rSges(7,1) + t87 * rSges(7,2) + rSges(7,3) * t274;
t499 = -pkin(8) * t274 - t43;
t453 = Icges(4,4) * t267;
t325 = Icges(4,1) * t264 + t453;
t163 = Icges(4,5) * t268 + t265 * t325;
t432 = t163 * t264;
t454 = Icges(4,4) * t264;
t320 = Icges(4,2) * t267 + t454;
t161 = Icges(4,6) * t268 + t265 * t320;
t435 = t161 * t267;
t296 = t432 + t435;
t498 = t268 * t296;
t297 = t146 * t248 + t148 * t247;
t497 = t268 * t297;
t300 = t141 * t248 + t143 * t247;
t496 = t268 * t300;
t463 = rSges(5,2) * t248;
t349 = rSges(5,1) * t247 + t463;
t286 = t268 * t349;
t350 = rSges(4,1) * t264 + rSges(4,2) * t267;
t287 = t268 * t350;
t424 = t263 * t268;
t179 = -t248 * t421 - t424;
t419 = t266 * t268;
t422 = t265 * t263;
t180 = -t248 * t422 + t419;
t428 = t247 * t265;
t111 = t180 * rSges(7,1) + t179 * rSges(7,2) + rSges(7,3) * t428;
t257 = t268 * pkin(5);
t495 = -pkin(8) * t428 - t111 - t257;
t346 = rSges(6,2) * t248 - rSges(6,3) * t247;
t494 = qJD(3) * t346;
t311 = Icges(5,5) * t247 + Icges(5,6) * t248;
t493 = -Icges(5,3) * t265 + t268 * t311;
t313 = Icges(4,5) * t264 + Icges(4,6) * t267;
t492 = -Icges(4,3) * t265 + t268 * t313;
t316 = Icges(6,4) * t247 + Icges(6,5) * t248;
t149 = Icges(6,1) * t265 + t268 * t316;
t490 = -Icges(4,6) * t265 + t268 * t320;
t488 = -Icges(4,5) * t265 + t268 * t325;
t177 = t248 * t419 - t422;
t178 = t248 * t424 + t421;
t348 = -t178 * rSges(7,1) - t177 * rSges(7,2);
t427 = t247 * t268;
t110 = -rSges(7,3) * t427 - t348;
t347 = rSges(7,1) * t263 + rSges(7,2) * t266;
t136 = rSges(7,3) * t248 + t247 * t347;
t395 = qJD(3) * t268;
t403 = qJD(1) * t265;
t365 = t248 * t395;
t372 = t247 * t403;
t273 = -t365 + t372;
t367 = t247 * t395;
t270 = t265 * t357 - t367;
t89 = t266 * t270 - t358 * t424;
t90 = t263 * t270 + t358 * t419;
t353 = t90 * rSges(7,1) + t89 * rSges(7,2);
t44 = rSges(7,3) * t273 + t353;
t392 = qJD(6) * t247;
t91 = (rSges(7,1) * t266 - rSges(7,2) * t263) * t392 + (-rSges(7,3) * t247 + t248 * t347) * qJD(3);
t21 = (-t136 * t395 - t44) * t248 + (qJD(3) * t110 + t136 * t403 - t268 * t91) * t247;
t22 = (-t136 * t398 + t43) * t248 + (-qJD(3) * t111 - t136 * t402 - t265 * t91) * t247;
t63 = t111 * t248 - t136 * t428;
t64 = -t248 * t110 - t136 * t427;
t487 = qJD(1) * (t265 * t63 + t268 * t64) + t21 * t265 - t22 * t268;
t486 = 2 * m(4);
t485 = 2 * m(5);
t484 = 2 * m(6);
t483 = 2 * m(7);
t260 = t265 ^ 2;
t261 = t268 ^ 2;
t482 = m(6) / 0.2e1;
t481 = m(7) / 0.2e1;
t480 = -pkin(1) - pkin(5);
t479 = t248 / 0.2e1;
t478 = -t264 / 0.2e1;
t477 = t265 / 0.2e1;
t476 = t267 / 0.2e1;
t475 = t268 / 0.2e1;
t474 = -rSges(6,1) - pkin(1);
t473 = rSges(3,2) - pkin(1);
t472 = -rSges(5,3) - pkin(1);
t471 = rSges(7,3) + pkin(8);
t465 = rSges(4,2) * t264;
t222 = rSges(4,1) * t267 - t465;
t470 = m(4) * t222;
t469 = pkin(3) * t264;
t468 = pkin(3) * t267;
t467 = pkin(7) * t268;
t466 = t265 * pkin(5);
t464 = rSges(5,2) * t247;
t462 = rSges(6,2) * t247;
t461 = rSges(6,3) * t248;
t460 = t265 * rSges(6,1);
t459 = t265 * rSges(4,3);
t458 = t265 * rSges(5,3);
t256 = t268 * rSges(6,1);
t255 = t268 * rSges(4,3);
t254 = t268 * rSges(5,3);
t423 = t264 * t265;
t242 = pkin(3) * t423;
t262 = -qJ(4) - pkin(7);
t388 = pkin(3) * t395;
t351 = qJD(4) * t265 - t262 * t402 - t267 * t388;
t123 = t268 * ((t242 - t467) * qJD(1) + t351);
t426 = t248 * t265;
t381 = qJ(5) * t426;
t376 = pkin(4) * t365 + qJ(5) * t367 + qJD(1) * t381;
t394 = qJD(5) * t248;
t457 = t123 + t268 * (pkin(4) * t372 + t268 * t394 - t376);
t354 = pkin(3) * t501 + qJD(4) * t268 + t262 * t403;
t127 = pkin(7) * t403 + t354;
t377 = pkin(4) * t274 + qJ(5) * t368;
t393 = qJD(5) * t265;
t456 = -t127 - (-qJ(5) * t402 - t393) * t248 - t377;
t449 = Icges(7,4) * t263;
t448 = Icges(7,4) * t266;
t438 = qJ(5) * t248;
t322 = Icges(7,1) * t263 + t448;
t135 = Icges(7,5) * t248 + t247 * t322;
t437 = t135 * t263;
t436 = t161 * t264;
t434 = t490 * t264;
t433 = t490 * t267;
t431 = t163 * t267;
t430 = t488 * t264;
t429 = t488 * t267;
t425 = t248 * t268;
t420 = t265 * t267;
t390 = pkin(8) * t427;
t418 = t110 - t390 + t466;
t152 = rSges(5,1) * t428 + rSges(5,2) * t426 + t254;
t176 = t242 + (-pkin(7) - t262) * t268;
t417 = -t152 - t176;
t410 = t265 * t262 + t268 * t469;
t175 = -t265 * pkin(7) - t410;
t157 = t268 * t175;
t228 = pkin(4) * t427;
t382 = qJ(5) * t425;
t166 = -t228 + t382;
t416 = t268 * t166 + t157;
t227 = pkin(4) * t428;
t415 = -t227 + t381 - t176;
t414 = -t166 - t175;
t197 = pkin(4) * t248 + qJ(5) * t247;
t243 = pkin(3) * t420;
t413 = t265 * t197 + t243;
t411 = qJD(1) * t243 + t264 * t388;
t409 = qJ(2) * t402 + qJD(2) * t265;
t408 = t268 * pkin(1) + t265 * qJ(2);
t407 = t260 + t261;
t139 = Icges(5,3) * t268 + t265 * t311;
t406 = qJD(1) * t139;
t150 = Icges(6,1) * t268 - t265 * t316;
t405 = qJD(1) * t150;
t159 = Icges(4,3) * t268 + t265 * t313;
t404 = qJD(1) * t159;
t401 = qJD(3) * t247;
t400 = qJD(3) * t248;
t399 = qJD(3) * t264;
t391 = -rSges(4,3) - pkin(1) - pkin(7);
t389 = pkin(3) * t399;
t105 = Icges(7,5) * t180 + Icges(7,6) * t179 + Icges(7,3) * t428;
t107 = Icges(7,4) * t180 + Icges(7,2) * t179 + Icges(7,6) * t428;
t109 = Icges(7,1) * t180 + Icges(7,4) * t179 + Icges(7,5) * t428;
t304 = t107 * t266 + t109 * t263;
t37 = Icges(7,5) * t88 + Icges(7,6) * t87 + Icges(7,3) * t274;
t39 = Icges(7,4) * t88 + Icges(7,2) * t87 + Icges(7,6) * t274;
t41 = Icges(7,1) * t88 + Icges(7,4) * t87 + Icges(7,5) * t274;
t11 = (qJD(3) * t304 + t37) * t248 + (-qJD(3) * t105 + t263 * t41 + t266 * t39 + (-t107 * t263 + t109 * t266) * qJD(6)) * t247;
t310 = Icges(7,5) * t263 + Icges(7,6) * t266;
t133 = Icges(7,3) * t248 + t247 * t310;
t315 = Icges(7,2) * t266 + t449;
t134 = Icges(7,6) * t248 + t247 * t315;
t84 = (Icges(7,5) * t266 - Icges(7,6) * t263) * t392 + (-Icges(7,3) * t247 + t248 * t310) * qJD(3);
t85 = (-Icges(7,2) * t263 + t448) * t392 + (-Icges(7,6) * t247 + t248 * t315) * qJD(3);
t86 = (Icges(7,1) * t266 - t449) * t392 + (-Icges(7,5) * t247 + t248 * t322) * qJD(3);
t16 = t133 * t274 + t87 * t134 + t88 * t135 + t179 * t85 + t180 * t86 + t428 * t84;
t387 = t11 / 0.2e1 + t16 / 0.2e1;
t104 = Icges(7,5) * t178 + Icges(7,6) * t177 - Icges(7,3) * t427;
t106 = Icges(7,4) * t178 + Icges(7,2) * t177 - Icges(7,6) * t427;
t108 = Icges(7,1) * t178 + Icges(7,4) * t177 - Icges(7,5) * t427;
t305 = t106 * t266 + t108 * t263;
t38 = Icges(7,5) * t90 + Icges(7,6) * t89 + Icges(7,3) * t273;
t40 = Icges(7,4) * t90 + Icges(7,2) * t89 + Icges(7,6) * t273;
t42 = Icges(7,1) * t90 + Icges(7,4) * t89 + Icges(7,5) * t273;
t10 = (qJD(3) * t305 + t38) * t248 + (-qJD(3) * t104 + t263 * t42 + t266 * t40 + (-t106 * t263 + t108 * t266) * qJD(6)) * t247;
t17 = t133 * t273 + t89 * t134 + t90 * t135 + t177 * t85 + t178 * t86 - t427 * t84;
t386 = -t17 / 0.2e1 - t10 / 0.2e1;
t31 = t105 * t248 + t247 * t304;
t46 = t133 * t428 + t134 * t179 + t135 * t180;
t385 = t31 / 0.2e1 + t46 / 0.2e1;
t30 = t104 * t248 + t247 * t305;
t45 = -t133 * t427 + t177 * t134 + t178 * t135;
t384 = t45 / 0.2e1 + t30 / 0.2e1;
t383 = t480 * t265;
t151 = qJD(5) * t247 + (-pkin(4) * t247 + t438) * qJD(3);
t369 = t267 * t402;
t237 = pkin(3) * t369;
t380 = t265 * t151 + t197 * t402 + t237;
t345 = t461 + t462;
t379 = t265 * t345 - t256 + t415;
t378 = t197 * t403 + t411;
t375 = -rSges(5,1) * t274 - t402 * t463;
t374 = rSges(4,1) * t501 + rSges(4,2) * t369;
t169 = rSges(4,1) * t423 + rSges(4,2) * t420 + t255;
t253 = t268 * qJ(2);
t373 = t253 + t410;
t363 = qJD(6) * t134 * t263;
t362 = -qJ(2) - t469;
t361 = -t197 - t468;
t360 = pkin(8) * t248 + t136;
t204 = t350 * qJD(3);
t359 = t204 * t407;
t356 = t415 + t495;
t355 = t228 + t373;
t352 = t316 * qJD(3) / 0.2e1 + (t311 + t313) * t500;
t199 = rSges(5,1) * t248 - t464;
t289 = t354 + t409;
t271 = t289 + t377;
t23 = (-t382 + t383) * qJD(1) + t271 - t248 * t393 - t499;
t251 = qJD(2) * t268;
t290 = t251 - t351;
t272 = t290 + t376;
t24 = (qJD(3) * t471 - qJD(5)) * t425 + (t480 * t268 + ((-pkin(4) - t471) * t247 + t362) * t265) * qJD(1) + t272 - t353;
t343 = -t23 * t268 + t265 * t24;
t26 = -t104 * t427 + t177 * t106 + t178 * t108;
t27 = -t105 * t427 + t177 * t107 + t178 * t109;
t340 = t26 * t268 - t265 * t27;
t18 = t26 * t265 + t268 * t27;
t28 = t104 * t428 + t106 * t179 + t108 * t180;
t29 = t105 * t428 + t107 * t179 + t109 * t180;
t339 = t265 * t29 - t268 * t28;
t19 = t28 * t265 + t268 * t29;
t338 = t265 * t31 - t268 * t30;
t337 = t30 * t265 + t268 * t31;
t329 = qJD(1) * t360;
t32 = t265 * t329 + (pkin(8) * t401 - t151 - t91) * t268 + t378;
t33 = t268 * t329 + (t91 + (-pkin(8) * t247 - t469) * qJD(3)) * t265 + t380;
t336 = t33 * t265 - t268 * t32;
t208 = rSges(6,3) * t368;
t285 = -t462 + (-rSges(6,3) - qJ(5)) * t248;
t269 = t265 * t474 + t268 * t285;
t34 = t208 + (-rSges(6,2) * qJD(3) - qJD(5)) * t426 + t269 * qJD(1) + t271;
t35 = (-t394 - t494) * t268 + (t474 * t268 + (t461 + (rSges(6,2) - pkin(4)) * t247 + t362) * t265) * qJD(1) + t272;
t335 = t265 * t35 - t268 * t34;
t48 = t383 + (t247 * t471 - t438) * t268 + t348 + t355;
t291 = -t262 * t268 + t242 + t408;
t284 = t227 + t291;
t49 = t284 - t381 - t495;
t334 = t265 * t48 - t268 * t49;
t188 = t345 * qJD(3);
t51 = -t346 * t403 + (-t151 - t188) * t268 + t378;
t52 = -t346 * t402 + (t188 - t389) * t265 + t380;
t333 = t52 * t265 - t268 * t51;
t332 = t265 * t64 - t268 * t63;
t69 = t269 + t355;
t70 = t265 * t285 + t256 + t284;
t330 = t265 * t69 - t268 * t70;
t326 = Icges(4,1) * t267 - t454;
t321 = -Icges(4,2) * t264 + t453;
t317 = Icges(6,4) * t248 - Icges(6,5) * t247;
t314 = Icges(4,5) * t267 - Icges(4,6) * t264;
t312 = Icges(5,5) * t248 - Icges(5,6) * t247;
t303 = t110 * t265 + t111 * t268;
t299 = -t247 * t489 - t248 * t491;
t298 = t145 * t248 + t147 * t247;
t295 = -t430 - t433;
t294 = (t482 + t481) * t401;
t292 = t266 * t135 * t392 + t400 * t437 + (t134 * t397 + t84) * t248 + (t263 * t86 + t266 * t85) * t247;
t288 = rSges(3,3) * t268 + t265 * t473;
t283 = t299 * t265;
t282 = t298 * t265;
t281 = t295 * t265;
t280 = qJD(3) * t326;
t278 = qJD(3) * t321;
t189 = t349 * qJD(3);
t174 = -rSges(3,2) * t268 + t265 * rSges(3,3) + t408;
t173 = t253 + t288;
t170 = t459 - t287;
t154 = t268 * t345 + t460;
t153 = t458 - t286;
t138 = (-t199 - t468) * t268;
t137 = t199 * t265 + t243;
t131 = t251 + (t473 * t268 + (-rSges(3,3) - qJ(2)) * t265) * qJD(1);
t130 = qJD(1) * t288 + t409;
t129 = t169 + t408 + t467;
t128 = t265 * t391 + t253 + t287;
t118 = qJD(1) * t492 + t314 * t398;
t117 = -t314 * t395 + t404;
t115 = (t346 + t361) * t268;
t114 = -t265 * t346 + t413;
t113 = t291 + t152;
t112 = t265 * t472 + t286 + t373;
t103 = t317 * t395 + t405;
t102 = -qJD(1) * t149 - t317 * t398;
t93 = qJD(1) * t493 + t312 * t398;
t92 = -t312 * t395 + t406;
t80 = t199 * t402 + t237 + (-t189 - t389) * t265;
t79 = t189 * t268 + t199 * t403 + t411;
t74 = t251 + t222 * t395 + (t391 * t268 + (-qJ(2) - t350) * t265) * qJD(1);
t73 = (-rSges(4,2) * t399 + qJD(1) * t391) * t265 + t374 + t409;
t72 = (-t360 + t361) * t268;
t71 = t265 * t360 + t413;
t68 = -t265 * t492 - t268 * t295;
t67 = t265 * t159 - t498;
t66 = -t268 * t492 + t281;
t65 = t159 * t268 + t265 * t296;
t62 = t150 * t268 - t265 * t297;
t61 = t149 * t268 - t282;
t60 = t265 * t150 + t497;
t59 = t265 * t149 + t268 * t298;
t58 = -t265 * t493 - t268 * t299;
t57 = t265 * t139 - t496;
t56 = -t268 * t493 + t283;
t55 = t139 * t268 + t265 * t300;
t54 = t199 * t395 + (t472 * t268 + (-t349 + t362) * t265) * qJD(1) + t290;
t53 = (-rSges(5,2) * t401 + qJD(1) * t472) * t265 + t289 - t375;
t50 = t133 * t248 + (t134 * t266 + t437) * t247;
t47 = t303 * t247;
t36 = t154 * t268 + t265 * t379 + t416;
t25 = t265 * t356 + t268 * t418 + t416;
t20 = ((-qJD(3) * t133 - t363) * t247 + t292) * t248;
t15 = t261 * t494 + (rSges(6,2) * t366 - t208 + t456) * t265 + ((t379 + t256) * t268 + (-t154 + t414 + t460) * t265) * qJD(1) + t457;
t14 = t303 * t400 + (t265 * t44 + t268 * t43 + (t110 * t268 - t265 * t111) * qJD(1)) * t247;
t13 = t247 * t339 + t46 * t248;
t12 = -t247 * t340 + t45 * t248;
t9 = (-pkin(8) * t365 + t44) * t268 + (t456 + t499) * t265 + ((t356 + t257) * t268 + (t390 + t414 - t418 + t466) * t265) * qJD(1) + t457;
t8 = t105 * t273 + t89 * t107 + t90 * t109 + t177 * t39 + t178 * t41 - t37 * t427;
t7 = t104 * t273 + t89 * t106 + t90 * t108 + t177 * t40 + t178 * t42 - t38 * t427;
t6 = t105 * t274 + t87 * t107 + t88 * t109 + t179 * t39 + t180 * t41 + t37 * t428;
t5 = t104 * t274 + t87 * t106 + t88 * t108 + t179 * t40 + t180 * t42 + t38 * t428;
t4 = qJD(1) * t340 + t7 * t265 + t268 * t8;
t3 = -qJD(1) * t339 + t5 * t265 + t268 * t6;
t2 = (-qJD(3) * t340 + t17) * t248 + (qJD(1) * t18 - qJD(3) * t45 + t265 * t8 - t268 * t7) * t247;
t1 = (qJD(3) * t339 + t16) * t248 + (qJD(1) * t19 - qJD(3) * t46 + t265 * t6 - t268 * t5) * t247;
t75 = [0.2e1 * m(3) * (t130 * t174 + t131 * t173) + (t128 * t74 + t129 * t73) * t486 + (t112 * t54 + t113 * t53) * t485 + (t34 * t70 + t35 * t69) * t484 + (t23 * t49 + t24 * t48) * t483 - t325 * t396 + t320 * t399 + t292 - t264 * t280 - t267 * t278 + (-t323 - t308) * t400 - t507 * t248 + (t318 + t306 - t133) * t401 + (-t363 - t506) * t247; m(7) * ((t265 * t49 + t268 * t48) * qJD(1) + t343) + m(5) * (t265 * t54 - t268 * t53 + (t112 * t268 + t113 * t265) * qJD(1)) + m(6) * ((t265 * t70 + t268 * t69) * qJD(1) + t335) + m(4) * (t265 * t74 - t268 * t73 + (t128 * t268 + t129 * t265) * qJD(1)) + m(3) * (-t130 * t268 + t265 * t131 + (t173 * t268 + t174 * t265) * qJD(1)); 0; m(5) * (t112 * t80 + t113 * t79 + t137 * t54 + t138 * t53) + m(6) * (t114 * t35 + t115 * t34 + t51 * t70 + t52 * t69) + m(7) * (t23 * t72 + t24 * t71 + t32 * t49 + t33 * t48) + (m(4) * (t129 * t204 - t222 * t73) + (qJD(1) * t490 + t265 * t278) * t478 + (qJD(1) * t488 + t265 * t280) * t476 + t352 * t268 + (-t435 / 0.2e1 - t432 / 0.2e1) * qJD(3) + t387) * t268 + (m(4) * (-t128 * t204 + t222 * t74) + (qJD(1) * t161 - t321 * t395) * t478 + (qJD(1) * t163 - t326 * t395) * t476 + t352 * t265 + (t433 / 0.2e1 + t430 / 0.2e1) * qJD(3) - t386) * t265 + ((qJD(3) * t502 - t395 * t509) * t477 + (qJD(3) * t505 + t265 * t506) * t475 + (t475 * t504 + t477 * t503) * qJD(1)) * t248 + ((qJD(3) * t504 + t395 * t508) * t477 + (-qJD(3) * t503 - t265 * t507) * t475 + (-t502 * t475 + t477 * t505) * qJD(1)) * t247 + ((t129 * t470 + t436 / 0.2e1 - t431 / 0.2e1 + (t148 / 0.2e1 - t143 / 0.2e1) * t248 + (-t146 / 0.2e1 + t141 / 0.2e1) * t247 - t385) * t265 + (t434 / 0.2e1 - t429 / 0.2e1 + t128 * t470 + (-t147 / 0.2e1 - t489 / 0.2e1) * t248 + (t145 / 0.2e1 + t491 / 0.2e1) * t247 + t384) * t268) * qJD(1); m(5) * (t80 * t265 - t268 * t79 + (t137 * t268 + t138 * t265) * qJD(1)) + m(6) * ((t114 * t268 + t115 * t265) * qJD(1) + t333) + m(7) * ((t265 * t72 + t268 * t71) * qJD(1) + t336) - m(4) * t359; (t25 * t9 + t32 * t72 + t33 * t71) * t483 + t265 * t4 + t268 * t3 + (t114 * t52 + t115 * t51 + t15 * t36) * t484 + (t137 * t80 + t138 * t79 + (t153 * t268 + t265 * t417 + t157) * (t123 + (-t127 + t375) * t265 + (-t199 * t261 + t260 * t464) * qJD(3) + ((t417 + t254) * t268 + (-t153 - t175 + t286 + t458) * t265) * qJD(1))) * t485 + t265 * ((t265 * t117 + (-t67 + t281) * qJD(1)) * t265 + (t68 * qJD(1) + (t161 * t399 - t163 * t396 + t404) * t268 + (t118 + (t429 - t434) * qJD(3) + t296 * qJD(1)) * t265) * t268) + t265 * ((t265 * t92 + (-t57 + t283) * qJD(1)) * t265 + (t58 * qJD(1) + (t141 * t401 - t143 * t400 + t406) * t268 + (t93 + (-t247 * t491 + t248 * t489) * qJD(3) + t300 * qJD(1)) * t265) * t268) + t265 * ((t265 * t103 + (-t60 - t282) * qJD(1)) * t265 + (t59 * qJD(1) + (-t146 * t401 + t148 * t400 + t405) * t268 + (t102 + (-t145 * t247 + t147 * t248) * qJD(3) - t297 * qJD(1)) * t265) * t268) + t268 * ((t118 * t268 + (t66 + t498) * qJD(1)) * t268 + (-t65 * qJD(1) + (-t396 * t488 + t399 * t490) * t265 + (t117 + (t431 - t436) * qJD(3) + (-t159 + t295) * qJD(1)) * t268) * t265) + ((-t265 * t169 + t170 * t268) * (-t265 * t374 + (-t222 * t261 + t260 * t465) * qJD(3) + ((-t169 + t255) * t268 + (-t170 + t287 + t459) * t265) * qJD(1)) - t222 * t359) * t486 + t268 * ((t268 * t93 + (t56 + t496) * qJD(1)) * t268 + (-t55 * qJD(1) + (-t400 * t489 + t401 * t491) * t265 + (t92 + (-t141 * t247 + t143 * t248) * qJD(3) + (-t139 + t299) * qJD(1)) * t268) * t265) + t268 * ((t268 * t102 + (t61 - t497) * qJD(1)) * t268 + (-t62 * qJD(1) + (t145 * t401 - t147 * t400) * t265 + (t103 + (t146 * t247 - t148 * t248) * qJD(3) + (-t150 - t298) * qJD(1)) * t268) * t265) + (-t19 + (-t55 - t62 - t65) * t268 + (-t56 - t61 - t66) * t265) * t403 + (t18 + (t57 + t60 + t67) * t268 + (t58 + t59 + t68) * t265) * t402; m(7) * (-qJD(1) * t334 + t265 * t23 + t24 * t268) + m(5) * (t265 * t53 + t268 * t54 + (-t112 * t265 + t113 * t268) * qJD(1)) + m(6) * (-qJD(1) * t330 + t265 * t34 + t268 * t35); 0; m(7) * (t265 * t32 + t268 * t33 + (-t265 * t71 + t268 * t72) * qJD(1)) + m(6) * (t265 * t51 + t268 * t52 + (-t114 * t265 + t115 * t268) * qJD(1)) + m(5) * (t265 * t79 + t268 * t80 + (-t137 * t265 + t138 * t268) * qJD(1)); 0; 0.2e1 * (t330 * t482 + t334 * t481) * t401 + 0.2e1 * ((-t402 * t48 - t403 * t49 - t343) * t481 + (-t402 * t69 - t403 * t70 - t335) * t482) * t248; 0.2e1 * t407 * t294; 0.2e1 * ((-t395 * t72 + t398 * t71 + t9) * t481 + (t114 * t398 - t115 * t395 + t15) * t482) * t247 + 0.2e1 * ((qJD(3) * t25 - t402 * t71 - t403 * t72 - t336) * t481 + (qJD(3) * t36 - t114 * t402 - t115 * t403 - t333) * t482) * t248; 0; 0.4e1 * (0.1e1 - t407) * t248 * t294; m(7) * (t21 * t48 + t22 * t49 + t23 * t63 + t24 * t64) + t20 + (t265 * t385 - t268 * t384) * t400 + (-qJD(3) * t50 + t386 * t268 + t387 * t265 + (t265 * t384 + t268 * t385) * qJD(1)) * t247; m(7) * t487; m(7) * (t14 * t25 + t21 * t71 + t22 * t72 + t32 * t63 + t33 * t64 + t47 * t9) + (qJD(1) * t12 / 0.2e1 - t18 * t400 / 0.2e1 + t1 / 0.2e1 + (qJD(1) * t30 + t11) * t479) * t268 + (t2 / 0.2e1 + t19 * t400 / 0.2e1 - qJD(1) * t13 / 0.2e1 + (-qJD(1) * t31 + t10) * t479) * t265 + (-t268 * t4 / 0.2e1 + t3 * t477 + t337 * t500 + (t18 * t477 + t19 * t475) * qJD(1)) * t247; m(7) * (-qJD(1) * t332 + t21 * t268 + t22 * t265); m(7) * ((qJD(3) * t332 + t14) * t247 + (qJD(3) * t47 - t487) * t248); (t14 * t47 + t21 * t64 + t22 * t63) * t483 + (t20 + (-t268 * t12 + t265 * t13 + t248 * t338) * qJD(3)) * t248 + (t265 * t1 - t268 * t2 + t248 * (-t10 * t268 + t11 * t265) + (-t247 * t338 - 0.2e1 * t50 * t248) * qJD(3) + (t265 * t12 + t268 * t13 + t248 * t337) * qJD(1)) * t247;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t75(1) t75(2) t75(4) t75(7) t75(11) t75(16); t75(2) t75(3) t75(5) t75(8) t75(12) t75(17); t75(4) t75(5) t75(6) t75(9) t75(13) t75(18); t75(7) t75(8) t75(9) t75(10) t75(14) t75(19); t75(11) t75(12) t75(13) t75(14) t75(15) t75(20); t75(16) t75(17) t75(18) t75(19) t75(20) t75(21);];
Mq  = res;
