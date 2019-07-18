% Calculate time derivative of joint inertia matrix for
% S5RRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
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
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(4,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_inertiaDJ_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR1_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR1_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR1_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:20:46
% EndTime: 2019-07-18 17:21:24
% DurationCPUTime: 15.86s
% Computational Cost: add. (16870->704), mult. (24544->1000), div. (0->0), fcn. (22905->8), ass. (0->376)
t296 = cos(qJ(2));
t293 = sin(qJ(2));
t478 = Icges(4,4) * t293;
t480 = Icges(3,4) * t293;
t520 = -t478 - t480 + (-Icges(3,2) - Icges(4,2)) * t296;
t477 = Icges(4,4) * t296;
t479 = Icges(3,4) * t296;
t519 = -t477 - t479 + (-Icges(3,1) - Icges(4,1)) * t293;
t518 = (t519 * t293 + t296 * t520) * qJD(2);
t294 = sin(qJ(1));
t288 = t294 ^ 2;
t297 = cos(qJ(1));
t289 = t297 ^ 2;
t428 = t288 + t289;
t290 = qJ(2) + qJ(4);
t280 = sin(t290);
t287 = qJD(2) + qJD(4);
t295 = cos(qJ(5));
t292 = sin(qJ(5));
t473 = Icges(6,4) * t295;
t351 = -Icges(6,2) * t292 + t473;
t281 = cos(t290);
t449 = t281 * t287;
t474 = Icges(6,4) * t292;
t116 = t351 * t449 + (Icges(6,6) * t287 + (-Icges(6,2) * t295 - t474) * qJD(5)) * t280;
t185 = -Icges(6,6) * t281 + t280 * t351;
t357 = Icges(6,1) * t295 - t474;
t186 = -Icges(6,5) * t281 + t280 * t357;
t517 = -t116 * t292 + (-t185 * t295 - t186 * t292) * qJD(5);
t475 = Icges(5,4) * t281;
t352 = -Icges(5,2) * t280 + t475;
t191 = Icges(5,6) * t294 + t297 * t352;
t476 = Icges(5,4) * t280;
t358 = Icges(5,1) * t281 - t476;
t193 = Icges(5,5) * t294 + t297 * t358;
t337 = t191 * t280 - t193 * t281;
t516 = t294 * t337;
t190 = -Icges(5,6) * t297 + t294 * t352;
t192 = -Icges(5,5) * t297 + t294 * t358;
t338 = t190 * t280 - t192 * t281;
t515 = t297 * t338;
t356 = -Icges(3,2) * t293 + t479;
t209 = Icges(3,6) * t294 + t297 * t356;
t362 = Icges(3,1) * t296 - t480;
t213 = Icges(3,5) * t294 + t297 * t362;
t333 = t209 * t293 - t213 * t296;
t315 = t333 * t294;
t208 = -Icges(3,6) * t297 + t294 * t356;
t212 = -Icges(3,5) * t297 + t294 * t362;
t334 = t208 * t293 - t212 * t296;
t316 = t334 * t297;
t354 = -Icges(4,2) * t293 + t477;
t207 = Icges(4,6) * t294 + t297 * t354;
t360 = Icges(4,1) * t296 - t478;
t211 = Icges(4,5) * t294 + t297 * t360;
t335 = t207 * t293 - t211 * t296;
t317 = t335 * t294;
t206 = -Icges(4,6) * t297 + t294 * t354;
t210 = -Icges(4,5) * t297 + t294 * t360;
t336 = t206 * t293 - t210 * t296;
t318 = t336 * t297;
t347 = Icges(6,5) * t295 - Icges(6,6) * t292;
t115 = t347 * t449 + (Icges(6,3) * t287 + (-Icges(6,5) * t292 - Icges(6,6) * t295) * qJD(5)) * t280;
t459 = t185 * t292;
t514 = -t287 * t459 - t115;
t440 = t295 * t297;
t444 = t292 * t294;
t220 = -t281 * t444 - t440;
t442 = t294 * t295;
t443 = t292 * t297;
t221 = t281 * t442 - t443;
t370 = -t221 * rSges(6,1) - t220 * rSges(6,2);
t451 = t280 * t294;
t145 = rSges(6,3) * t451 - t370;
t222 = -t281 * t443 + t442;
t223 = t281 * t440 + t444;
t450 = t280 * t297;
t146 = t223 * rSges(6,1) + t222 * rSges(6,2) + rSges(6,3) * t450;
t513 = -t294 * t145 - t297 * t146;
t369 = rSges(6,1) * t295 - rSges(6,2) * t292;
t118 = t369 * t449 + (rSges(6,3) * t287 + (-rSges(6,1) * t292 - rSges(6,2) * t295) * qJD(5)) * t280;
t187 = -rSges(6,3) * t281 + t280 * t369;
t424 = qJD(1) * t294;
t178 = t187 * t424;
t512 = -t118 * t297 + t178;
t439 = t296 * t297;
t488 = rSges(4,2) * t293;
t511 = t294 * rSges(4,3) - t297 * t488;
t216 = rSges(4,1) * t439 + t511;
t282 = t294 * qJ(3);
t245 = pkin(1) * t439 + t282;
t175 = t216 + t245;
t420 = qJD(2) * t296;
t423 = qJD(1) * t297;
t510 = -t293 * t423 - t294 * t420;
t441 = t294 * t296;
t481 = t297 * rSges(3,3);
t489 = rSges(3,2) * t293;
t215 = rSges(3,1) * t441 - t294 * t489 - t481;
t374 = rSges(3,1) * t296 - t489;
t485 = rSges(3,3) * t294;
t217 = t297 * t374 + t485;
t509 = t215 * t297 - t217 * t294;
t232 = rSges(5,1) * t280 + rSges(5,2) * t281;
t298 = pkin(2) + pkin(1);
t491 = pkin(1) - t298;
t246 = t491 * t293;
t492 = t293 * pkin(1);
t389 = t246 - t492;
t377 = -t232 + t389;
t163 = t377 * t294;
t164 = t377 * t297;
t508 = t163 * t294 + t164 * t297;
t348 = Icges(5,5) * t281 - Icges(5,6) * t280;
t188 = -Icges(5,3) * t297 + t294 * t348;
t507 = qJD(1) * t188;
t349 = Icges(4,5) * t296 - Icges(4,6) * t293;
t202 = -Icges(4,3) * t297 + t294 * t349;
t350 = Icges(3,5) * t296 - Icges(3,6) * t293;
t204 = -Icges(3,3) * t297 + t294 * t350;
t381 = qJD(1) * t281 - qJD(5);
t445 = t287 * t297;
t411 = t280 * t445;
t506 = t294 * t381 + t411;
t229 = Icges(5,2) * t281 + t476;
t230 = Icges(5,1) * t280 + t475;
t332 = t229 * t280 - t230 * t281;
t504 = qJD(1) * t332 + t348 * t287;
t503 = 2 * m(3);
t502 = 2 * m(4);
t501 = 2 * m(5);
t500 = 2 * m(6);
t499 = t294 / 0.2e1;
t498 = -t297 / 0.2e1;
t497 = -rSges(4,1) - pkin(1);
t496 = -rSges(6,3) - pkin(4);
t257 = rSges(3,1) * t293 + rSges(3,2) * t296;
t495 = m(3) * t257;
t494 = m(5) * t232;
t490 = rSges(5,1) * t281;
t487 = rSges(4,2) * t296;
t486 = rSges(5,2) * t280;
t484 = rSges(4,3) * t297;
t139 = Icges(6,5) * t221 + Icges(6,6) * t220 + Icges(6,3) * t451;
t141 = Icges(6,4) * t221 + Icges(6,2) * t220 + Icges(6,6) * t451;
t143 = Icges(6,1) * t221 + Icges(6,4) * t220 + Icges(6,5) * t451;
t346 = -t141 * t292 + t143 * t295;
t382 = -qJD(5) * t281 + qJD(1);
t331 = t295 * t382;
t447 = t287 * t294;
t127 = t294 * t331 + (t280 * t447 - t297 * t381) * t292;
t330 = t382 * t292;
t446 = t287 * t295;
t128 = t381 * t440 + (-t280 * t446 + t330) * t294;
t410 = t281 * t447;
t307 = t280 * t423 + t410;
t69 = Icges(6,5) * t128 + Icges(6,6) * t127 + Icges(6,3) * t307;
t71 = Icges(6,4) * t128 + Icges(6,2) * t127 + Icges(6,6) * t307;
t73 = Icges(6,1) * t128 + Icges(6,4) * t127 + Icges(6,5) * t307;
t19 = (t287 * t346 - t69) * t281 + (t139 * t287 - t292 * t71 + t295 * t73 + (-t141 * t295 - t143 * t292) * qJD(5)) * t280;
t483 = t19 * t297;
t140 = Icges(6,5) * t223 + Icges(6,6) * t222 + Icges(6,3) * t450;
t142 = Icges(6,4) * t223 + Icges(6,2) * t222 + Icges(6,6) * t450;
t144 = Icges(6,1) * t223 + Icges(6,4) * t222 + Icges(6,5) * t450;
t345 = -t142 * t292 + t144 * t295;
t125 = t292 * t506 + t297 * t331;
t126 = -t295 * t506 + t297 * t330;
t402 = t280 * t424;
t409 = t281 * t445;
t306 = -t402 + t409;
t68 = Icges(6,5) * t126 + Icges(6,6) * t125 + Icges(6,3) * t306;
t70 = Icges(6,4) * t126 + Icges(6,2) * t125 + Icges(6,6) * t306;
t72 = Icges(6,1) * t126 + Icges(6,4) * t125 + Icges(6,5) * t306;
t20 = (t287 * t345 - t68) * t281 + (t140 * t287 - t292 * t70 + t295 * t72 + (-t142 * t295 - t144 * t292) * qJD(5)) * t280;
t482 = t20 * t294;
t284 = t294 * rSges(5,3);
t454 = t229 * t287;
t453 = t230 * t287;
t452 = t280 * t287;
t448 = t281 * t297;
t438 = t296 * t298;
t291 = -pkin(3) - qJ(3);
t437 = t297 * t291;
t383 = -t291 * t294 + t297 * t438;
t181 = t383 - t245;
t436 = -t181 - t245;
t372 = -t486 + t490;
t194 = -t297 * rSges(5,3) + t294 * t372;
t195 = rSges(5,1) * t448 - rSges(5,2) * t450 + t284;
t129 = t294 * t194 + t297 * t195;
t400 = t293 * t424;
t270 = pkin(1) * t400;
t435 = -t246 * t424 + t270;
t283 = t297 * qJ(3);
t244 = pkin(1) * t441 - t283;
t434 = t294 * t244 + t297 * t245;
t433 = rSges(5,2) * t402 + rSges(5,3) * t423;
t422 = qJD(2) * t293;
t395 = t298 * t422;
t432 = -t291 * t424 - t294 * t395;
t431 = rSges(4,2) * t400 + rSges(4,3) * t423;
t421 = qJD(2) * t294;
t397 = t293 * t421;
t269 = pkin(1) * t397;
t279 = qJD(3) * t297;
t430 = -t269 - t279;
t275 = qJ(3) * t423;
t278 = qJD(3) * t294;
t429 = t275 + t278;
t189 = Icges(5,3) * t294 + t297 * t348;
t427 = qJD(1) * t189;
t203 = Icges(4,3) * t294 + t297 * t349;
t426 = qJD(1) * t203;
t205 = Icges(3,3) * t294 + t297 * t350;
t425 = qJD(1) * t205;
t419 = qJD(2) * t297;
t417 = pkin(4) * t449;
t263 = pkin(4) * t448;
t415 = pkin(1) * t420;
t61 = -t139 * t281 + t280 * t346;
t184 = -Icges(6,3) * t281 + t280 * t347;
t78 = t184 * t451 + t185 * t220 + t186 * t221;
t414 = t61 / 0.2e1 + t78 / 0.2e1;
t62 = -t140 * t281 + t280 * t345;
t79 = t184 * t450 + t185 * t222 + t186 * t223;
t413 = t62 / 0.2e1 + t79 / 0.2e1;
t117 = t357 * t449 + (Icges(6,5) * t287 + (-Icges(6,1) * t292 - t473) * qJD(5)) * t280;
t408 = t280 * t295 * t117 + t281 * t186 * t446 + t184 * t452;
t401 = t281 * t424;
t308 = -t401 - t411;
t323 = t232 * t287;
t407 = t294 * (-t294 * t323 + (t297 * t372 + t284) * qJD(1)) + t297 * (rSges(5,1) * t308 - rSges(5,2) * t409 + t433) + t194 * t423;
t406 = t126 * rSges(6,1) + t125 * rSges(6,2) + rSges(6,3) * t409;
t396 = t293 * t419;
t398 = t296 * t424;
t405 = t294 * (qJD(1) * t245 + t430) + t297 * ((-t396 - t398) * pkin(1) + t429) + t244 * t423;
t404 = t279 - t432;
t403 = t187 * t423;
t393 = t449 / 0.2e1;
t392 = t296 * t491;
t391 = t424 / 0.2e1;
t390 = t423 / 0.2e1;
t256 = rSges(4,1) * t293 + t487;
t388 = -t256 - t492;
t132 = -qJD(1) * t190 - t297 * t454;
t387 = t193 * t287 + t132;
t133 = qJD(1) * t191 - t294 * t454;
t386 = t192 * t287 + t133;
t134 = -qJD(1) * t192 - t297 * t453;
t385 = -t191 * t287 + t134;
t135 = qJD(1) * t193 - t294 * t453;
t384 = t190 * t287 - t135;
t180 = -t294 * t392 + t283 + t437;
t380 = t294 * t180 + t297 * t181 + t434;
t80 = pkin(4) * t280 * t428 - t513;
t379 = -pkin(4) * t452 - t118;
t378 = -t187 + t389;
t376 = qJD(2) * t392 - t415;
t74 = -rSges(6,3) * t402 + t406;
t371 = t128 * rSges(6,1) + t127 * rSges(6,2);
t75 = rSges(6,3) * t307 + t371;
t375 = t145 * t423 + t294 * t75 + t297 * t74 + t417 * t428;
t373 = rSges(4,1) * t296 - t488;
t50 = t139 * t451 + t141 * t220 + t143 * t221;
t51 = t140 * t451 + t142 * t220 + t144 * t221;
t39 = t294 * t51 - t297 * t50;
t368 = t294 * t50 + t297 * t51;
t52 = t139 * t450 + t141 * t222 + t143 * t223;
t53 = t140 * t450 + t142 * t222 + t144 * t223;
t40 = t294 * t53 - t297 * t52;
t367 = t294 * t52 + t297 * t53;
t366 = t294 * t62 - t297 * t61;
t365 = t294 * t61 + t297 * t62;
t228 = Icges(5,5) * t280 + Icges(5,6) * t281;
t319 = t287 * t228;
t130 = -t297 * t319 - t507;
t131 = -t294 * t319 + t427;
t13 = t125 * t141 + t126 * t143 + t139 * t306 + t222 * t71 + t223 * t73 + t450 * t69;
t14 = t125 * t142 + t126 * t144 + t140 * t306 + t222 * t70 + t223 * t72 + t450 * t68;
t8 = qJD(1) * t367 - t13 * t297 + t14 * t294;
t88 = -t188 * t297 - t294 * t338;
t89 = -t189 * t297 - t516;
t90 = t188 * t294 - t515;
t91 = t189 * t294 - t297 * t337;
t364 = t39 * t424 + t40 * t423 + (-t90 * t423 - t88 * t424) * t297 + (t8 + (t91 * qJD(1) + (t133 * t280 - t135 * t281 + t190 * t449 + t192 * t452 - t507) * t297) * t297 + t89 * t424 + t91 * t423 + ((t90 + t516) * qJD(1) + (-t131 + t385 * t281 - t387 * t280 + (t189 - t338) * qJD(1)) * t297 + t130 * t294) * t294) * t294;
t363 = -t438 - t490;
t344 = t145 * t297 - t146 * t294;
t305 = t363 + t486;
t167 = (rSges(5,3) - t291) * t297 + t305 * t294;
t168 = t195 + t383;
t339 = t167 * t297 + t168 * t294;
t201 = t372 * t287;
t329 = -t201 + t376;
t328 = t378 * t297;
t326 = t294 * (t269 + (-t297 * t392 - t282) * qJD(1) + t432) + t297 * (-t275 + t491 * t396 + (t441 * t491 - t437) * qJD(1)) + t180 * t423 + t405;
t325 = t296 * t497 + t488;
t324 = t280 * t496 - t438;
t322 = qJD(2) * t256;
t310 = qJD(2) * (-Icges(3,5) * t293 - Icges(3,6) * t296);
t309 = qJD(2) * (-Icges(4,5) * t293 - Icges(4,6) * t296);
t12 = (t131 * t297 + (t89 + t515) * qJD(1)) * t297 + (t88 * qJD(1) + (-t132 * t280 + t134 * t281 - t191 * t449 - t193 * t452 + t427) * t294 + (-t130 + t384 * t281 + t386 * t280 + (-t188 - t337) * qJD(1)) * t297) * t294;
t15 = t127 * t141 + t128 * t143 + t139 * t307 + t220 * t71 + t221 * t73 + t451 * t69;
t16 = t127 * t142 + t128 * t144 + t140 * t307 + t220 * t70 + t221 * t72 + t451 * t68;
t9 = qJD(1) * t368 - t15 * t297 + t16 * t294;
t304 = (-t12 - t9) * t297 + t364;
t303 = t376 + t379;
t24 = t280 * t368 - t281 * t78;
t25 = t280 * t367 - t281 * t79;
t29 = t115 * t450 + t116 * t222 + t117 * t223 + t125 * t185 + t126 * t186 + t184 * t306;
t3 = (t287 * t367 - t29) * t281 + (-qJD(1) * t40 + t13 * t294 + t14 * t297 + t287 * t79) * t280;
t30 = t115 * t451 + t116 * t220 + t117 * t221 + t127 * t185 + t128 * t186 + t184 * t307;
t4 = (t287 * t368 - t30) * t281 + (-qJD(1) * t39 + t15 * t294 + t16 * t297 + t287 * t78) * t280;
t302 = t3 * t499 + t4 * t498 - t281 * (qJD(1) * t365 + t482 - t483) / 0.2e1 + t24 * t391 + t25 * t390 + t366 * t452 / 0.2e1 + t9 * t451 / 0.2e1 + t8 * t450 / 0.2e1 + (t297 * t393 - t402 / 0.2e1) * t40 + (t280 * t390 + t294 * t393) * t39;
t301 = t294 * t324 - t437;
t197 = t352 * t287;
t198 = t358 * t287;
t300 = qJD(1) * t228 + (t198 - t454) * t281 + (-t197 - t453) * t280;
t299 = -t483 / 0.2e1 + t482 / 0.2e1 + (t280 * t385 + t281 * t387 + t294 * t504 + t300 * t297 + t29) * t499 + (-t280 * t384 + t281 * t386 + t300 * t294 - t297 * t504 + t30) * t498 + (t190 * t281 + t192 * t280 - t228 * t297 - t294 * t332 + t61 + t78) * t391 + (t191 * t281 + t193 * t280 + t228 * t294 - t297 * t332 + t62 + t79) * t390;
t261 = t294 * t281 * pkin(4);
t248 = qJD(1) * t263;
t241 = t374 * qJD(2);
t240 = t373 * qJD(2);
t214 = t294 * t373 - t484;
t200 = t388 * t297;
t199 = t388 * t294;
t174 = t294 * t325 + t283 + t484;
t166 = -t187 * t297 + t263;
t165 = -t187 * t294 + t261;
t162 = -rSges(3,1) * t397 + (rSges(3,1) * t439 + t485) * qJD(1) + t510 * rSges(3,2);
t161 = -t257 * t419 + (-t294 * t374 + t481) * qJD(1);
t152 = t294 * t310 + t425;
t151 = -qJD(1) * t204 + t297 * t310;
t150 = t294 * t309 + t426;
t149 = -qJD(1) * t202 + t297 * t309;
t148 = pkin(1) * t510 - t240 * t294 - t256 * t423;
t147 = t256 * t424 + t270 + (-t240 - t415) * t297;
t114 = t263 + t328;
t113 = t294 * t378 + t261;
t107 = t256 * t421 + ((-rSges(4,3) - qJ(3)) * t294 + t325 * t297) * qJD(1) - t430;
t106 = t497 * t398 + (t293 * t497 - t487) * t419 + t429 + t431;
t105 = pkin(4) * t450 + t146 + t383;
t104 = t301 + t370;
t103 = t205 * t294 - t297 * t333;
t102 = t204 * t294 - t316;
t101 = t203 * t294 - t297 * t335;
t100 = t202 * t294 - t318;
t99 = -t205 * t297 - t315;
t98 = -t204 * t297 - t294 * t334;
t97 = -t203 * t297 - t317;
t96 = -t202 * t297 - t294 * t336;
t95 = t232 * t447 + (t297 * t305 - t284) * qJD(1) + t404;
t94 = t278 + t363 * t424 + (-qJD(1) * t291 - t323 - t395) * t297 + t433;
t93 = -t146 * t281 - t187 * t450;
t92 = t145 * t281 + t187 * t451;
t87 = qJD(1) * t164 + t294 * t329;
t86 = t232 * t424 + t297 * t329 + t435;
t85 = -t184 * t281 + (t186 * t295 - t459) * t280;
t84 = t85 * t452;
t83 = t344 * t280;
t82 = t294 * t379 + t248 - t403;
t81 = pkin(4) * t308 + t512;
t65 = t380 + t129;
t64 = qJD(1) * t328 + t294 * t303 + t248;
t63 = -pkin(4) * t401 + t297 * t303 + t178 + t435;
t60 = -t195 * t424 + t407;
t55 = t324 * t423 + t410 * t496 - t371 + t404;
t54 = t278 + (-t395 + t417) * t297 + t301 * qJD(1) + t406;
t45 = t80 + t380;
t44 = (t187 * t447 + t75) * t281 + (t118 * t294 - t145 * t287 + t403) * t280;
t43 = (-t187 * t445 - t74) * t281 + (t146 * t287 + t512) * t280;
t41 = t280 * t517 + t514 * t281 + t408;
t32 = -t146 * t424 + t375;
t31 = (-t195 + t436) * t424 + t326 + t407;
t26 = t344 * t449 + (qJD(1) * t513 - t294 * t74 + t297 * t75) * t280;
t21 = (-t146 + t436) * t424 + t326 + t375;
t1 = [(t104 * t55 + t105 * t54) * t500 + (t167 * t95 + t168 * t94) * t501 + (t106 * t175 + t107 * t174) * t502 + (t161 * t217 + t162 * t215) * t503 + t408 + t230 * t449 - t229 * t452 + (t197 + t514) * t281 + (t360 + t362 + t520) * t422 + (t354 + t356 - t519) * t420 + (t198 + t517) * t280; m(3) * ((-t161 * t294 + t162 * t297) * t257 + t509 * t241) + t299 + m(6) * (t104 * t63 + t105 * t64 + t113 * t54 + t114 * t55) + m(5) * (t163 * t94 + t164 * t95 + t167 * t86 + t168 * t87) + m(4) * (t106 * t199 + t107 * t200 + t147 * t174 + t148 * t175) + ((-t217 * t495 + (t207 / 0.2e1 + t209 / 0.2e1) * t296 + (t211 / 0.2e1 + t213 / 0.2e1) * t293) * t297 + (-t215 * t495 + (t206 / 0.2e1 + t208 / 0.2e1) * t296 + (t210 / 0.2e1 + t212 / 0.2e1) * t293) * t294) * qJD(1) + (t518 * t297 + ((-t206 - t208) * t296 + (-t210 - t212) * t293) * qJD(1)) * t499 + (t518 * t294 + ((t207 + t209) * t296 + (t211 + t213) * t293) * qJD(1)) * t498 + (t318 / 0.2e1 + t316 / 0.2e1 - t315 / 0.2e1 - t317 / 0.2e1 + (t350 + t349) * (t289 / 0.2e1 + t288 / 0.2e1)) * qJD(2); (t200 * t147 + t199 * t148 + (t214 * t294 + t216 * t297 + t434) * ((qJD(1) * t214 - t297 * t322 + t431) * t297 + (-t294 * t322 + (t511 - t175) * qJD(1)) * t294 + t405)) * t502 + t364 + t294 * ((t294 * t151 + (t102 + t315) * qJD(1)) * t294 + (t103 * qJD(1) + (t208 * t420 + t212 * t422) * t297 + (-t152 + (-t209 * t296 - t213 * t293) * qJD(2) + (t205 - t334) * qJD(1)) * t294) * t297) - t297 * ((t297 * t152 + (t99 + t316) * qJD(1)) * t297 + (t98 * qJD(1) + (-t209 * t420 - t213 * t422 + t425) * t294 + (-t151 + (t208 * t296 + t212 * t293) * qJD(2) - t333 * qJD(1)) * t297) * t294) + t294 * ((t294 * t149 + (t100 + t317) * qJD(1)) * t294 + (t101 * qJD(1) + (t206 * t420 + t210 * t422) * t297 + (-t150 + (-t207 * t296 - t211 * t293) * qJD(2) + (t203 - t336) * qJD(1)) * t294) * t297) - t297 * ((t297 * t150 + (t97 + t318) * qJD(1)) * t297 + (t96 * qJD(1) + (-t207 * t420 - t211 * t422 + t426) * t294 + (-t149 + (t206 * t296 + t210 * t293) * qJD(2) - t335 * qJD(1)) * t297) * t294) + (t113 * t64 + t114 * t63 + t21 * t45) * t500 + (t163 * t87 + t164 * t86 + t31 * t65) * t501 + ((t215 * t294 + t217 * t297) * (qJD(1) * t509 + t297 * t161 + t294 * t162) + t428 * t257 * t241) * t503 - t297 * t12 - t297 * t9 + ((-t96 - t98) * t297 + (t97 + t99) * t294) * t424 + ((-t100 - t102) * t297 + (t101 + t103) * t294) * t423; m(6) * (t294 * t55 - t297 * t54 + (t104 * t297 + t105 * t294) * qJD(1)) + m(5) * (qJD(1) * t339 + t294 * t95 - t297 * t94) + m(4) * (-t106 * t297 + t107 * t294 + (t174 * t297 + t175 * t294) * qJD(1)); m(6) * (t294 * t63 - t297 * t64 + (t113 * t294 + t114 * t297) * qJD(1)) + m(5) * (qJD(1) * t508 + t294 * t86 - t297 * t87) + m(4) * (t147 * t294 - t148 * t297 + (t199 * t294 + t200 * t297) * qJD(1)); 0; t299 + m(6) * (t104 * t81 + t105 * t82 + t165 * t54 + t166 * t55) - m(5) * t339 * t201 + (-t294 * t94 - t297 * t95 + (t167 * t294 - t168 * t297) * qJD(1)) * t494; m(6) * (t113 * t82 + t114 * t81 + t165 * t64 + t166 * t63 + t21 * t80 + t32 * t45) + m(5) * (t129 * t31 - t201 * t508 + t60 * t65) + (-t294 * t87 - t297 * t86 + (-t163 * t297 + t164 * t294) * qJD(1)) * t494 + t304; m(6) * (t294 * t81 - t297 * t82 + (t165 * t294 + t166 * t297) * qJD(1)); (t201 * t232 * t428 + t129 * t60) * t501 + (t165 * t82 + t166 * t81 + t32 * t80) * t500 + t304; m(6) * (t104 * t44 + t105 * t43 + t54 * t93 + t55 * t92) + t84 + (-t41 + (t294 * t414 + t297 * t413) * t287) * t281 + ((t20 / 0.2e1 + t29 / 0.2e1) * t297 + (t19 / 0.2e1 + t30 / 0.2e1) * t294 + (-t294 * t413 + t297 * t414) * qJD(1)) * t280; t302 + m(6) * (t113 * t43 + t114 * t44 + t21 * t83 + t26 * t45 + t63 * t92 + t64 * t93); m(6) * (t294 * t44 - t297 * t43 + (t294 * t93 + t297 * t92) * qJD(1)); t302 + m(6) * (t165 * t43 + t166 * t44 + t26 * t80 + t32 * t83 + t81 * t92 + t82 * t93); (t26 * t83 + t43 * t93 + t44 * t92) * t500 + (t41 * t281 - t84 + (t294 * t24 + t297 * t25 - t281 * t365) * t287) * t281 + (t297 * t3 + t294 * t4 + t365 * t452 + (-t19 * t294 - t20 * t297 - t287 * t85) * t281 + (t297 * t24 - t294 * t25 + t281 * t366) * qJD(1)) * t280;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq  = res;
