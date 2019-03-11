% Calculate time derivative of joint inertia matrix for
% S6RPRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 02:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_inertiaDJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR6_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR6_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR6_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:53:18
% EndTime: 2019-03-09 02:53:43
% DurationCPUTime: 15.30s
% Computational Cost: add. (16714->823), mult. (19424->1156), div. (0->0), fcn. (17779->10), ass. (0->399)
t522 = Icges(6,5) / 0.2e1;
t521 = Icges(6,6) / 0.2e1;
t520 = Icges(6,3) / 0.2e1;
t283 = sin(qJ(1));
t491 = -t283 / 0.2e1;
t519 = -qJD(1) / 0.2e1;
t282 = sin(qJ(3));
t284 = cos(qJ(3));
t408 = qJD(3) * t284;
t285 = cos(qJ(1));
t413 = qJD(1) * t285;
t518 = t282 * t413 + t283 * t408;
t275 = qJ(3) + pkin(9);
t263 = sin(t275);
t265 = cos(t275);
t409 = qJD(3) * t283;
t374 = t265 * t409;
t517 = t263 * t413 + t374;
t407 = qJD(3) * t285;
t375 = t263 * t407;
t414 = qJD(1) * t283;
t516 = t265 * t414 + t375;
t376 = t263 * t409;
t379 = t265 * t413;
t288 = t376 - t379;
t274 = pkin(10) + qJ(6);
t262 = sin(t274);
t264 = cos(t274);
t464 = Icges(7,4) * t264;
t326 = -Icges(7,2) * t262 + t464;
t151 = Icges(7,6) * t263 + t265 * t326;
t465 = Icges(7,4) * t262;
t332 = Icges(7,1) * t264 - t465;
t152 = Icges(7,5) * t263 + t265 * t332;
t515 = -t151 * t262 + t152 * t264;
t278 = sin(pkin(10));
t279 = cos(pkin(10));
t333 = Icges(6,1) * t279 - Icges(6,4) * t278;
t159 = Icges(6,5) * t263 + t265 * t333;
t514 = -t159 / 0.2e1;
t259 = pkin(5) * t279 + pkin(4);
t281 = -pkin(8) - qJ(5);
t365 = qJD(6) * t263 + qJD(1);
t308 = t365 * t283;
t364 = qJD(1) * t263 + qJD(6);
t500 = t285 * t364 + t374;
t86 = -t262 * t500 - t264 * t308;
t87 = -t262 * t308 + t264 * t500;
t48 = t87 * rSges(7,1) + t86 * rSges(7,2) + rSges(7,3) * t288;
t513 = t259 * t517 + t281 * t379 + t48;
t358 = rSges(5,1) * t263 + rSges(5,2) * t265;
t298 = t285 * t358;
t359 = rSges(4,1) * t282 + rSges(4,2) * t284;
t299 = t285 * t359;
t469 = Icges(4,4) * t282;
t330 = Icges(4,2) * t284 + t469;
t186 = Icges(4,6) * t285 + t283 * t330;
t468 = Icges(4,4) * t284;
t336 = Icges(4,1) * t282 + t468;
t188 = Icges(4,5) * t285 + t283 * t336;
t312 = t186 * t284 + t188 * t282;
t295 = t312 * t285;
t467 = Icges(5,4) * t263;
t328 = Icges(5,2) * t265 + t467;
t166 = Icges(5,6) * t285 + t283 * t328;
t466 = Icges(5,4) * t265;
t334 = Icges(5,1) * t263 + t466;
t168 = Icges(5,5) * t285 + t283 * t334;
t314 = t166 * t265 + t168 * t263;
t297 = t314 * t285;
t481 = pkin(4) - t259;
t512 = t481 * t263;
t437 = t283 * t262;
t443 = t264 * t285;
t176 = -t263 * t437 + t443;
t436 = t283 * t264;
t177 = t262 * t285 + t263 * t436;
t442 = t265 * t283;
t100 = t177 * rSges(7,1) + t176 * rSges(7,2) - rSges(7,3) * t442;
t440 = t278 * t285;
t253 = pkin(5) * t440;
t445 = t263 * t283;
t511 = t259 * t445 + t281 * t442 + t100 + t253;
t510 = -t151 * t264 - t152 * t262;
t377 = t284 * t413;
t384 = rSges(4,1) * t518 + rSges(4,2) * t377;
t403 = -rSges(4,3) - pkin(1) - pkin(7);
t410 = qJD(3) * t282;
t419 = qJ(2) * t413 + qJD(2) * t283;
t88 = (-rSges(4,2) * t410 + qJD(1) * t403) * t283 + t384 + t419;
t479 = rSges(4,2) * t282;
t237 = rSges(4,1) * t284 - t479;
t268 = qJD(2) * t285;
t89 = t268 + t237 * t407 + (t403 * t285 + (-qJ(2) - t359) * t283) * qJD(1);
t509 = t283 * t89 - t285 * t88;
t449 = t259 * t263;
t471 = rSges(7,3) - t281;
t508 = -t265 * t471 + t449;
t373 = t265 * t407;
t507 = t283 * t364 - t373;
t322 = Icges(5,5) * t263 + Icges(5,6) * t265;
t506 = -Icges(5,3) * t283 + t285 * t322;
t324 = Icges(4,5) * t282 + Icges(4,6) * t284;
t505 = -Icges(4,3) * t283 + t285 * t324;
t504 = -Icges(5,6) * t283 + t285 * t328;
t503 = -Icges(4,6) * t283 + t285 * t330;
t502 = -Icges(5,5) * t283 + t285 * t334;
t501 = -Icges(4,5) * t283 + t285 * t336;
t444 = t263 * t285;
t178 = t262 * t444 + t436;
t179 = -t263 * t443 + t437;
t354 = -t179 * rSges(7,1) - t178 * rSges(7,2);
t441 = t265 * t285;
t101 = rSges(7,3) * t441 - t354;
t353 = rSges(7,1) * t264 - rSges(7,2) * t262;
t154 = rSges(7,3) * t263 + t265 * t353;
t309 = t285 * t365;
t84 = -t262 * t507 + t264 * t309;
t85 = t262 * t309 + t264 * t507;
t361 = t85 * rSges(7,1) + t84 * rSges(7,2);
t47 = -rSges(7,3) * t516 + t361;
t404 = qJD(6) * t265;
t93 = (-rSges(7,1) * t262 - rSges(7,2) * t264) * t404 + (rSges(7,3) * t265 - t263 * t353) * qJD(3);
t21 = (-t154 * t407 - t47) * t263 + (-qJD(3) * t101 - t154 * t414 + t285 * t93) * t265;
t22 = (-t154 * t409 + t48) * t263 + (qJD(3) * t100 + t154 * t413 + t283 * t93) * t265;
t58 = t100 * t263 + t154 * t442;
t59 = -t263 * t101 + t154 * t441;
t499 = qJD(1) * (t283 * t58 + t285 * t59) + t21 * t283 - t22 * t285;
t498 = 2 * m(4);
t497 = 2 * m(5);
t496 = 2 * m(6);
t495 = 2 * m(7);
t276 = t283 ^ 2;
t277 = t285 ^ 2;
t494 = m(6) / 0.2e1;
t493 = m(7) / 0.2e1;
t492 = t263 / 0.2e1;
t489 = t285 / 0.2e1;
t488 = rSges(3,2) - pkin(1);
t487 = -rSges(5,3) - pkin(1);
t486 = m(4) * t237;
t485 = pkin(3) * t282;
t484 = pkin(3) * t284;
t483 = pkin(4) * t263;
t482 = pkin(7) * t285;
t273 = t285 * pkin(1);
t320 = Icges(7,5) * t264 - Icges(7,6) * t262;
t150 = Icges(7,3) * t263 + t265 * t320;
t411 = qJD(3) * t265;
t412 = qJD(3) * t263;
t90 = (-Icges(7,5) * t262 - Icges(7,6) * t264) * t404 + (Icges(7,3) * t265 - t263 * t320) * qJD(3);
t92 = (-Icges(7,1) * t262 - t464) * t404 + (Icges(7,5) * t265 - t263 * t332) * qJD(3);
t286 = t265 * t264 * t92 + t150 * t411 + t263 * t90 - t412 * t515;
t91 = (-Icges(7,2) * t264 - t465) * t404 + (Icges(7,6) * t265 - t263 * t326) * qJD(3);
t476 = t262 * t91;
t55 = t150 * t263 + t265 * t515;
t480 = ((qJD(6) * t510 - t476) * t265 + t286) * t263 + t55 * t411;
t478 = rSges(5,2) * t263;
t477 = rSges(6,3) * t265;
t475 = t283 * rSges(4,3);
t474 = t283 * rSges(5,3);
t272 = t285 * rSges(4,3);
t271 = t285 * rSges(5,3);
t432 = -qJ(5) - t281;
t367 = t432 * t265;
t470 = -(t367 + t512) * qJD(3) - t93;
t453 = t186 * t282;
t452 = t503 * t282;
t451 = t188 * t284;
t450 = t501 * t284;
t448 = t259 * t265;
t447 = t263 * t278;
t446 = t263 * t279;
t439 = t279 * t285;
t438 = t282 * t283;
t435 = t283 * t278;
t434 = t283 * t279;
t433 = t283 * t284;
t241 = pkin(4) * t444;
t402 = pkin(5) * t435;
t431 = t101 + t402 + t241 + (t367 - t449) * t285;
t256 = pkin(3) * t438;
t280 = -qJ(4) - pkin(7);
t400 = pkin(3) * t407;
t360 = qJD(4) * t283 - t280 * t413 - t284 * t400;
t135 = t285 * ((t256 - t482) * qJD(1) + t360);
t240 = pkin(4) * t445;
t386 = pkin(4) * t373 + qJ(5) * t516;
t406 = qJD(5) * t265;
t430 = t285 * (qJD(1) * t240 + t285 * t406 - t386) + t135;
t362 = pkin(3) * t518 + qJD(4) * t285 + t280 * t414;
t139 = pkin(7) * t414 + t362;
t387 = pkin(4) * t517 + qJ(5) * t376;
t405 = qJD(5) * t283;
t429 = -(-qJ(5) * t413 - t405) * t265 - t387 - t139;
t428 = t263 * t432 - t265 * t481 + t154;
t173 = rSges(5,1) * t445 + rSges(5,2) * t442 + t271;
t203 = t256 + (-pkin(7) - t280) * t285;
t427 = -t173 - t203;
t420 = t283 * t280 + t285 * t485;
t202 = -t283 * pkin(7) - t420;
t180 = t285 * t202;
t191 = qJ(5) * t441 - t241;
t426 = t285 * t191 + t180;
t198 = -t263 * t435 + t439;
t199 = t263 * t434 + t440;
t425 = t199 * rSges(6,1) + t198 * rSges(6,2);
t190 = -qJ(5) * t442 + t240;
t424 = -t190 - t203;
t423 = -t191 - t202;
t215 = pkin(4) * t265 + qJ(5) * t263;
t257 = pkin(3) * t433;
t422 = t283 * t215 + t257;
t421 = qJD(1) * t257 + t282 * t400;
t418 = t283 * qJ(2) + t273;
t417 = t276 + t277;
t164 = Icges(5,3) * t285 + t283 * t322;
t416 = qJD(1) * t164;
t184 = Icges(4,3) * t285 + t283 * t324;
t415 = qJD(1) * t184;
t401 = pkin(3) * t410;
t96 = Icges(7,4) * t177 + Icges(7,2) * t176 - Icges(7,6) * t442;
t98 = Icges(7,1) * t177 + Icges(7,4) * t176 - Icges(7,5) * t442;
t349 = t262 * t96 - t264 * t98;
t94 = Icges(7,5) * t177 + Icges(7,6) * t176 - Icges(7,3) * t442;
t32 = t263 * t94 - t265 * t349;
t49 = -t150 * t442 + t151 * t176 + t152 * t177;
t399 = t32 / 0.2e1 + t49 / 0.2e1;
t97 = Icges(7,4) * t179 + Icges(7,2) * t178 + Icges(7,6) * t441;
t99 = Icges(7,1) * t179 + Icges(7,4) * t178 + Icges(7,5) * t441;
t348 = t262 * t97 - t264 * t99;
t95 = Icges(7,5) * t179 + Icges(7,6) * t178 + Icges(7,3) * t441;
t33 = t263 * t95 - t265 * t348;
t50 = t150 * t441 + t178 * t151 + t179 * t152;
t398 = -t33 / 0.2e1 - t50 / 0.2e1;
t110 = Icges(6,5) * t199 + Icges(6,6) * t198 - Icges(6,3) * t442;
t397 = t110 * t442;
t396 = t110 * t441;
t200 = t263 * t440 + t434;
t300 = t263 * t439 - t435;
t111 = -Icges(6,5) * t300 + Icges(6,6) * t200 + Icges(6,3) * t441;
t395 = t111 * t442;
t394 = t111 * t441;
t393 = rSges(6,3) * t442 + t424 - t425;
t143 = -qJD(1) * t200 - t278 * t374;
t144 = qJD(1) * t300 + t279 * t374;
t392 = -t144 * rSges(6,1) - t143 * rSges(6,2) - rSges(6,3) * t376;
t170 = qJD(5) * t263 + (qJ(5) * t265 - t483) * qJD(3);
t250 = pkin(3) * t377;
t391 = t283 * t170 + t215 * t413 + t250;
t390 = t215 * t414 + t421;
t385 = -rSges(5,1) * t517 - rSges(5,2) * t379;
t192 = rSges(4,1) * t438 + rSges(4,2) * t433 + t272;
t270 = t285 * qJ(2);
t383 = t270 + t420;
t382 = -pkin(5) * t278 - pkin(1);
t371 = -qJ(2) - t485;
t370 = -t215 - t484;
t369 = t265 * (-rSges(6,3) - qJ(5));
t221 = t359 * qJD(3);
t368 = t221 * t417;
t366 = qJD(1) * t428;
t363 = t190 + t424 - t511;
t216 = rSges(5,1) * t265 - t478;
t141 = qJD(1) * t198 + t278 * t373;
t142 = qJD(1) * t199 - t279 * t373;
t357 = -t142 * rSges(6,1) - t141 * rSges(6,2);
t356 = rSges(6,1) * t300 - t200 * rSges(6,2);
t355 = rSges(6,1) * t279 - rSges(6,2) * t278;
t304 = t362 + t419;
t24 = (qJD(1) * t382 - t281 * t412 - t406) * t283 + t304 + t513;
t306 = t268 - t360;
t25 = (-t406 + (t263 * t471 + t448) * qJD(3)) * t285 + (t382 * t285 + (t371 - t508) * t283) * qJD(1) + t306 - t361;
t351 = -t24 * t285 + t283 * t25;
t26 = t283 * t366 + (-t170 + t470) * t285 + t390;
t27 = t285 * t366 + (-t401 - t470) * t283 + t391;
t350 = -t26 * t285 + t27 * t283;
t28 = t176 * t96 + t177 * t98 - t442 * t94;
t29 = t176 * t97 + t177 * t99 - t442 * t95;
t18 = t28 * t285 + t29 * t283;
t347 = t28 * t283 - t285 * t29;
t30 = t178 * t96 + t179 * t98 + t441 * t94;
t31 = t178 * t97 + t179 * t99 + t441 * t95;
t19 = t31 * t283 + t285 * t30;
t346 = t283 * t30 - t285 * t31;
t345 = t33 * t283 + t285 * t32;
t344 = t283 * t32 - t285 * t33;
t289 = -t283 * pkin(1) + t285 * t369;
t39 = qJD(1) * t289 - t265 * t405 + t304 + t387 - t392;
t40 = (rSges(6,3) * t412 - t406) * t285 + (-t273 + (t371 + t477 - t483) * t283) * qJD(1) + t306 + t357 + t386;
t343 = t283 * t40 - t285 * t39;
t153 = (-t263 * t355 + t477) * qJD(3);
t161 = rSges(6,3) * t263 + t265 * t355;
t52 = t161 * t414 + (-t153 - t170) * t285 + t390;
t53 = t161 * t413 + (t153 - t401) * t283 + t391;
t342 = t53 * t283 - t285 * t52;
t56 = t382 * t283 + t285 * t508 + t354 + t383;
t307 = -t280 * t285 + t256 + t418;
t57 = t307 + t511;
t341 = t283 * t56 - t285 * t57;
t340 = t283 * t59 - t285 * t58;
t74 = t241 + t289 + t356 + t383;
t75 = t283 * t369 + t240 + t307 + t425;
t338 = t283 * t74 - t285 * t75;
t337 = Icges(4,1) * t284 - t469;
t335 = Icges(5,1) * t265 - t467;
t331 = -Icges(4,2) * t282 + t468;
t329 = -Icges(5,2) * t263 + t466;
t327 = Icges(6,4) * t279 - Icges(6,2) * t278;
t325 = Icges(4,5) * t284 - Icges(4,6) * t282;
t323 = Icges(5,5) * t265 - Icges(5,6) * t263;
t321 = Icges(6,5) * t279 - Icges(6,6) * t278;
t319 = t100 * t285 + t101 * t283;
t313 = -t263 * t502 - t265 * t504;
t311 = -t282 * t501 - t284 * t503;
t310 = (t494 + t493) * t412;
t303 = rSges(3,3) * t285 + t283 * t488;
t296 = t313 * t283;
t294 = t311 * t283;
t293 = qJD(3) * t337;
t292 = qJD(3) * t335;
t291 = qJD(3) * t331;
t290 = qJD(3) * t329;
t208 = t358 * qJD(3);
t197 = -rSges(3,2) * t285 + t283 * rSges(3,3) + t418;
t196 = t270 + t303;
t193 = t475 - t299;
t174 = t474 - t298;
t163 = (-t216 - t484) * t285;
t162 = t216 * t283 + t257;
t158 = Icges(6,6) * t263 + t265 * t327;
t156 = t268 + (t488 * t285 + (-rSges(3,3) - qJ(2)) * t283) * qJD(1);
t155 = qJD(1) * t303 + t419;
t149 = (Icges(6,5) * t265 - t263 * t333) * qJD(3);
t148 = (Icges(6,6) * t265 - t263 * t327) * qJD(3);
t146 = t192 + t418 + t482;
t145 = t283 * t403 + t270 + t299;
t129 = qJD(1) * t505 + t325 * t409;
t128 = -t325 * t407 + t415;
t125 = t307 + t173;
t124 = t283 * t487 + t298 + t383;
t123 = rSges(6,3) * t441 - t356;
t117 = qJD(1) * t506 + t323 * t409;
t116 = -t323 * t407 + t416;
t115 = -Icges(6,1) * t300 + Icges(6,4) * t200 + Icges(6,5) * t441;
t114 = Icges(6,1) * t199 + Icges(6,4) * t198 - Icges(6,5) * t442;
t113 = -Icges(6,4) * t300 + Icges(6,2) * t200 + Icges(6,6) * t441;
t112 = Icges(6,4) * t199 + Icges(6,2) * t198 - Icges(6,6) * t442;
t106 = (-t161 + t370) * t285;
t105 = t161 * t283 + t422;
t104 = t216 * t413 + t250 + (-t208 - t401) * t283;
t103 = t208 * t285 + t216 * t414 + t421;
t79 = -t283 * t505 - t285 * t311;
t78 = t283 * t184 - t295;
t77 = -t285 * t505 + t294;
t76 = t184 * t285 + t283 * t312;
t73 = -t283 * t506 - t285 * t313;
t72 = t283 * t164 - t297;
t71 = -t285 * t506 + t296;
t70 = t164 * t285 + t283 * t314;
t69 = Icges(6,1) * t144 + Icges(6,4) * t143 + Icges(6,5) * t288;
t68 = Icges(6,1) * t142 + Icges(6,4) * t141 - Icges(6,5) * t516;
t67 = Icges(6,4) * t144 + Icges(6,2) * t143 + Icges(6,6) * t288;
t66 = Icges(6,4) * t142 + Icges(6,2) * t141 - Icges(6,6) * t516;
t63 = t216 * t407 + (t487 * t285 + (-t358 + t371) * t283) * qJD(1) + t306;
t62 = (-rSges(5,2) * t412 + qJD(1) * t487) * t283 + t304 - t385;
t61 = (t370 - t428) * t285;
t60 = t283 * t428 + t422;
t51 = t319 * t265;
t46 = Icges(7,1) * t87 + Icges(7,4) * t86 + Icges(7,5) * t288;
t45 = Icges(7,1) * t85 + Icges(7,4) * t84 - Icges(7,5) * t516;
t44 = Icges(7,4) * t87 + Icges(7,2) * t86 + Icges(7,6) * t288;
t43 = Icges(7,4) * t85 + Icges(7,2) * t84 - Icges(7,6) * t516;
t42 = Icges(7,5) * t87 + Icges(7,6) * t86 + Icges(7,3) * t288;
t41 = Icges(7,5) * t85 + Icges(7,6) * t84 - Icges(7,3) * t516;
t38 = t123 * t285 + t283 * t393 + t426;
t37 = t200 * t113 - t115 * t300 + t394;
t36 = t200 * t112 - t114 * t300 + t396;
t35 = t113 * t198 + t115 * t199 - t395;
t34 = t112 * t198 + t114 * t199 - t397;
t23 = t283 * t363 + t285 * t431 + t426;
t17 = t150 * t288 + t86 * t151 + t87 * t152 + t176 * t91 + t177 * t92 - t442 * t90;
t16 = -t150 * t516 + t84 * t151 + t85 * t152 + t178 * t91 + t179 * t92 + t441 * t90;
t15 = t285 * (-rSges(6,3) * t375 - t357) + (t392 + t429) * t283 + (t393 * t285 + (-t123 + t423) * t283) * qJD(1) + t430;
t14 = t319 * t412 + (-t283 * t47 - t285 * t48 + (t283 * t100 - t101 * t285) * qJD(1)) * t265;
t13 = t50 * t263 - t265 * t346;
t12 = t49 * t263 - t265 * t347;
t11 = (qJD(3) * t348 + t41) * t263 + (qJD(3) * t95 - t262 * t43 + t264 * t45 + (-t262 * t99 - t264 * t97) * qJD(6)) * t265;
t10 = (qJD(3) * t349 + t42) * t263 + (qJD(3) * t94 - t262 * t44 + t264 * t46 + (-t262 * t98 - t264 * t96) * qJD(6)) * t265;
t9 = t95 * t376 + t176 * t43 + t177 * t45 + t86 * t97 + t87 * t99 + (-t283 * t41 - t413 * t95) * t265;
t8 = t94 * t376 + t176 * t44 + t177 * t46 + t86 * t96 + t87 * t98 + (-t283 * t42 - t413 * t94) * t265;
t7 = -t95 * t375 + t178 * t43 + t179 * t45 + t84 * t97 + t85 * t99 + (t285 * t41 - t414 * t95) * t265;
t6 = -t94 * t375 + t178 * t44 + t179 * t46 + t84 * t96 + t85 * t98 + (t285 * t42 - t414 * t94) * t265;
t5 = (t47 + (t263 * t281 - t448) * t407 + t386) * t285 + (t281 * t376 + t387 + t429 - t513) * t283 + ((t363 + t253) * t285 + (t402 + ((-qJ(5) + t281) * t265 - t512) * t285 + t423 - t431) * t283) * qJD(1) + t430;
t4 = -qJD(1) * t347 + t9 * t283 + t285 * t8;
t3 = -qJD(1) * t346 + t7 * t283 + t285 * t6;
t2 = (qJD(3) * t347 + t17) * t263 + (-qJD(1) * t18 + qJD(3) * t49 - t283 * t8 + t285 * t9) * t265;
t1 = (qJD(3) * t346 + t16) * t263 + (-qJD(1) * t19 + qJD(3) * t50 - t283 * t6 + t285 * t7) * t265;
t20 = [(t24 * t57 + t25 * t56) * t495 + (t39 * t75 + t40 * t74) * t496 + (t124 * t63 + t125 * t62) * t497 + (t145 * t89 + t146 * t88) * t498 + 0.2e1 * m(3) * (t155 * t197 + t156 * t196) - t336 * t408 - t334 * t411 + t286 - t282 * t293 - t284 * t291 + t330 * t410 + t510 * t404 + (Icges(6,3) * t411 - t292) * t263 + (t158 * t278 - t159 * t279 - t263 * t321 + t328) * t412 + (Icges(6,3) * t412 - t148 * t278 + t149 * t279 + t321 * t411 - t290 - t476) * t265; m(7) * ((t283 * t57 + t285 * t56) * qJD(1) + t351) + m(6) * ((t283 * t75 + t285 * t74) * qJD(1) + t343) + m(5) * (t283 * t63 - t285 * t62 + (t124 * t285 + t125 * t283) * qJD(1)) + m(4) * ((t145 * t285 + t146 * t283) * qJD(1) + t509) + m(3) * (-t155 * t285 + t283 * t156 + (t196 * t285 + t197 * t283) * qJD(1)); 0; m(4) * (t509 * t237 - (t145 * t283 - t146 * t285) * t221) + m(7) * (t24 * t61 + t25 * t60 + t26 * t57 + t27 * t56) + m(6) * (t105 * t40 + t106 * t39 + t52 * t75 + t53 * t74) + m(5) * (t103 * t125 + t104 * t124 + t162 * t63 + t163 * t62) + ((t143 * t521 + t144 * t522 + t288 * t520 + t290 * t491 + t504 * t519) * t285 + (t142 * t522 + t141 * t521 - t516 * t520 + t166 * t519 + t329 * t407 / 0.2e1) * t283) * t263 + (-t295 / 0.2e1 - t297 / 0.2e1 + (t313 + t311) * t491) * qJD(3) + ((t146 * t486 + t453 / 0.2e1 - t451 / 0.2e1 - t158 * t198 / 0.2e1 + t199 * t514 + (t166 / 0.2e1 - t110 / 0.2e1) * t263 + (t112 * t278 / 0.2e1 - t114 * t279 / 0.2e1) * t265 - t399) * t283 + (t145 * t486 + t200 * t158 / 0.2e1 + t300 * t514 + t452 / 0.2e1 - t450 / 0.2e1 + (t504 / 0.2e1 + t111 / 0.2e1) * t263 + (-t502 / 0.2e1 - t113 * t278 / 0.2e1 + t115 * t279 / 0.2e1) * t265 - t398) * t285) * qJD(1) + (-t322 - t324) * qJD(3) * (t277 / 0.2e1 + t276 / 0.2e1) + (t141 * t158 + t142 * t159 + t200 * t148 - t300 * t149 + t16 - (qJD(1) * t186 - t331 * t407) * t282 + (qJD(1) * t188 - t337 * t407) * t284 + t11 + (-t278 * t66 + t279 * t68 - t335 * t407) * t265 + (t111 * t265 + t113 * t447 - t115 * t446) * qJD(3)) * t283 / 0.2e1 + (t143 * t158 + t144 * t159 + t198 * t148 + t199 * t149 - (qJD(1) * t503 + t283 * t291) * t282 + (qJD(1) * t501 + t283 * t293) * t284 + t17 + t10 + (qJD(1) * t502 - t278 * t67 + t279 * t69 + t283 * t292) * t265 + (t110 * t265 + t112 * t447 - t114 * t446) * qJD(3)) * t489; m(5) * (-t103 * t285 + t104 * t283 + (t162 * t285 + t163 * t283) * qJD(1)) + m(6) * ((t105 * t285 + t106 * t283) * qJD(1) + t342) + m(7) * ((t283 * t61 + t285 * t60) * qJD(1) + t350) - m(4) * t368; t285 * t4 + t283 * t3 + (t23 * t5 + t26 * t61 + t27 * t60) * t495 + (t105 * t53 + t106 * t52 + t38 * t15) * t496 + (t162 * t104 + t163 * t103 + (t174 * t285 + t283 * t427 + t180) * (t135 + (-t139 + t385) * t283 + (-t216 * t277 + t276 * t478) * qJD(3) + ((t427 + t271) * t285 + (-t174 - t202 + t298 + t474) * t283) * qJD(1))) * t497 + ((-t283 * t192 + t193 * t285) * (-t283 * t384 + (-t237 * t277 + t276 * t479) * qJD(3) + ((-t192 + t272) * t285 + (-t193 + t299 + t475) * t283) * qJD(1)) - t237 * t368) * t498 + t285 * ((t285 * t129 + (t77 + t295) * qJD(1)) * t285 + (-t76 * qJD(1) + (-t408 * t501 + t410 * t503) * t283 + (t128 + (t451 - t453) * qJD(3) + (-t184 + t311) * qJD(1)) * t285) * t283) + t285 * ((t143 * t112 + t144 * t114 + t198 * t67 + t199 * t69 + (t35 - t396) * qJD(1)) * t285 + (t143 * t113 + t144 * t115 + t198 * t66 + t199 * t68 + (-t34 - t394) * qJD(1)) * t283) + t283 * ((t141 * t113 + t142 * t115 + t200 * t66 - t300 * t68 + (-t36 - t395) * qJD(1)) * t283 + (t141 * t112 + t142 * t114 + t200 * t67 - t300 * t69 + (t37 - t397) * qJD(1)) * t285) + t285 * ((t285 * t117 + (t71 + t297) * qJD(1)) * t285 + (-t70 * qJD(1) + (-t411 * t502 + t412 * t504) * t283 + (t116 + (-t166 * t263 + t168 * t265) * qJD(3) + (-t164 + t313) * qJD(1)) * t285) * t283) + t283 * ((t283 * t128 + (-t78 + t294) * qJD(1)) * t283 + (t79 * qJD(1) + (t186 * t410 - t188 * t408 + t415) * t285 + (t129 + (t450 - t452) * qJD(3) + t312 * qJD(1)) * t283) * t285) + t283 * ((t283 * t116 + (-t72 + t296) * qJD(1)) * t283 + (t73 * qJD(1) + (t166 * t412 - t168 * t411 + t416) * t285 + (t117 + (-t263 * t504 + t265 * t502) * qJD(3) + t314 * qJD(1)) * t283) * t285) + (-t18 + (-t34 - t70 - t76) * t285 + (-t35 - t71 - t77) * t283) * t414 + (t19 + (t36 + t72 + t78) * t285 + (t37 + t73 + t79) * t283) * t413; m(7) * (-qJD(1) * t341 + t283 * t24 + t25 * t285) + m(6) * (-qJD(1) * t338 + t283 * t39 + t285 * t40) + m(5) * (t283 * t62 + t285 * t63 + (-t124 * t283 + t125 * t285) * qJD(1)); 0; m(7) * (t283 * t26 + t27 * t285 + (-t283 * t60 + t285 * t61) * qJD(1)) + m(6) * (t283 * t52 + t285 * t53 + (-t105 * t283 + t106 * t285) * qJD(1)) + m(5) * (t283 * t103 + t104 * t285 + (-t162 * t283 + t163 * t285) * qJD(1)); 0; 0.2e1 * (t338 * t494 + t341 * t493) * t412 + 0.2e1 * ((-t413 * t56 - t414 * t57 - t351) * t493 + (-t413 * t74 - t414 * t75 - t343) * t494) * t265; 0.2e1 * t417 * t310; 0.2e1 * ((-t407 * t61 + t409 * t60 + t5) * t493 + (t105 * t409 - t106 * t407 + t15) * t494) * t263 + 0.2e1 * ((qJD(3) * t23 - t413 * t60 - t414 * t61 - t350) * t493 + (qJD(3) * t38 - t105 * t413 - t106 * t414 - t342) * t494) * t265; 0; 0.4e1 * (0.1e1 - t417) * t265 * t310; m(7) * (t21 * t56 + t22 * t57 + t24 * t58 + t25 * t59) + (t283 * t399 + t285 * t398) * t412 + ((t11 / 0.2e1 + t16 / 0.2e1) * t285 + (-t10 / 0.2e1 - t17 / 0.2e1) * t283 + (t283 * t398 - t285 * t399) * qJD(1)) * t265 + t480; m(7) * t499; m(7) * (t14 * t23 + t21 * t60 + t22 * t61 + t26 * t58 + t27 * t59 - t5 * t51) + (t2 / 0.2e1 + qJD(1) * t13 / 0.2e1 - t19 * t412 / 0.2e1 + (qJD(1) * t33 + t10) * t492) * t285 + (t12 * t519 + t1 / 0.2e1 + t18 * t412 / 0.2e1 + (-qJD(1) * t32 + t11) * t492) * t283 + (t3 * t489 + t4 * t491 + qJD(3) * t345 / 0.2e1 + (t19 * t491 - t285 * t18 / 0.2e1) * qJD(1)) * t265; m(7) * (-qJD(1) * t340 + t21 * t285 + t22 * t283); m(7) * ((qJD(3) * t340 + t14) * t263 + (-qJD(3) * t51 - t499) * t265); (-t14 * t51 + t21 * t59 + t22 * t58) * t495 + ((t283 * t12 - t285 * t13 + t263 * t344) * qJD(3) + t480) * t263 + (-t283 * t2 + t285 * t1 + t263 * (-t10 * t283 + t11 * t285) + (t55 * t263 - t265 * t344) * qJD(3) + (-t285 * t12 - t283 * t13 - t263 * t345) * qJD(1)) * t265;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t20(1) t20(2) t20(4) t20(7) t20(11) t20(16); t20(2) t20(3) t20(5) t20(8) t20(12) t20(17); t20(4) t20(5) t20(6) t20(9) t20(13) t20(18); t20(7) t20(8) t20(9) t20(10) t20(14) t20(19); t20(11) t20(12) t20(13) t20(14) t20(15) t20(20); t20(16) t20(17) t20(18) t20(19) t20(20) t20(21);];
Mq  = res;
